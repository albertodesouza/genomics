"""AlphaGenome prediction access, coordinate remapping and binning.

Three coordinate systems are supported for every request:

``reference``  reference-window offsets (= genomic positions); each haplotype's predictions are
               remapped through its own indels, deleted reference bases become NaN.
``haplotype``  raw prediction-array index (the haplotype's own consensus coordinates).
``aligned``    the bcftools_chain expanded axis used to build CNN training tensors (needs the
               alignment service; insertion columns are shared across the cohort).

Binned outputs (``resolution`` > 1 bp per row, e.g. 128 bp ChIP-seq) are expanded to bases on the
fly: every base reads the row of the bin its haplotype position falls in.

The pseudo-sample ``@reference`` is AlphaGenome's prediction of the reference window itself
(``references/windows/<gene>/predictions_ref``); its columns are matched to the dataset's track
order by ontology term / strand / assay, and it is identical in genomic and haplotype coordinates.
"""
from __future__ import annotations

import io
import os
import struct
import sys
import threading
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.visualizer.cache import DiskArrayCache, LRUCache, stable_key
from genomics.visualizer.coords import (
    HaplotypeEvents,
    ReferenceMap,
    SampleVariants,
    parse_events_batch,
    parse_window_vcf,
    reference_map,
    unpack_events,
)
from genomics.visualizer.datasets import Dataset, read_fasta_bytes
from genomics.visualizer.jobs import JobCancelled, is_cancelled

COORDINATE_SYSTEMS = ("reference", "haplotype", "aligned")
DIPLOID = "H1+H2"
MAX_BINS = 8192
MAX_SERIES = 32
MAX_GROUP_TRACKS = 16
REFERENCE_SAMPLE = "@reference"
ProgressFn = Callable[[float, str], None]


def _read_npz_member_fast(path: Path, member: str = "values.npy") -> Optional[np.ndarray]:
    """Inflate one ``.npz`` member with a single GIL-free ``zlib`` call (read-only result).

    ``np.load`` streams members through ``zipfile`` in small chunks, which serialises badly across
    threads; bulk cohort jobs load thousands of arrays, so this path matters. Returns ``None``
    when the member is absent or uses an unsupported layout, letting the caller fall back.
    """
    with open(path, "rb") as f:
        try:
            info = zipfile.ZipFile(f).getinfo(member)
        except (KeyError, zipfile.BadZipFile):
            return None
        if info.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED) or info.flag_bits & 0x1:
            return None
        f.seek(info.header_offset)
        header = f.read(30)
        name_len, extra_len = struct.unpack_from("<HH", header, 26)
        f.seek(info.header_offset + 30 + name_len + extra_len)
        raw = f.read(info.compress_size)
    buf = zlib.decompress(raw, -15, max(info.file_size, 1)) if info.compress_type == zipfile.ZIP_DEFLATED else raw
    stream = io.BytesIO(buf[:65536])
    version = np.lib.format.read_magic(stream)
    if version == (1, 0):
        shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
    elif version == (2, 0):
        shape, fortran, dtype = np.lib.format.read_array_header_2_0(stream)
    else:
        return None
    if dtype.hasobject:
        return None
    return np.frombuffer(buf, dtype=dtype, offset=stream.tell()).reshape(shape, order="F" if fortran else "C")


def npz_header(path: Path) -> Dict[str, object]:
    """Members, ``values`` shape and stored ``resolution`` of a prediction ``.npz`` without inflating ``values``."""
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        members = [n[:-4] if n.endswith(".npy") else n for n in names]
        shape: Optional[Tuple[int, ...]] = None
        if "values.npy" in names:
            with archive.open("values.npy") as handle:
                version = np.lib.format.read_magic(handle)
                reader = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
                shape = tuple(int(n) for n in reader(handle)[0])
        resolution = None
        if "resolution.npy" in names:
            with archive.open("resolution.npy") as handle:
                resolution = int(np.lib.format.read_array(handle))
    if shape is None:  # legacy layouts (track_<i> members, unnamed arrays)
        shape = tuple(load_prediction_matrix(path).shape)
    return {"members": members, "shape": shape, "resolution": resolution}


def load_prediction_matrix(path: Path) -> np.ndarray:
    """Load a prediction ``.npz`` as a float32 ``(length, tracks)`` matrix."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Prediction file not found: {path}")
    array = _read_npz_member_fast(path)
    if array is not None:
        return _as_matrix(array)
    with np.load(path) as data:
        keys = list(data.files)
        if "values" in keys:
            array = np.asarray(data["values"])
        else:
            track_keys = sorted((k for k in keys if k.startswith("track_") and k[6:].isdigit()), key=lambda k: int(k[6:]))
            if track_keys:
                array = np.column_stack([np.asarray(data[k]).reshape(-1) for k in track_keys])
            else:
                array_keys = [k for k in keys if not k.endswith("_metadata")]
                if not array_keys:
                    raise ValueError(f"No array payload in {path}")
                array = np.asarray(data[array_keys[0]])
    return _as_matrix(array)


def _as_matrix(array: np.ndarray) -> np.ndarray:
    if array.ndim == 0:
        array = array.reshape(1, 1)
    elif array.ndim == 1:
        array = array.reshape(-1, 1)
    elif array.ndim > 2:
        array = array.reshape(-1, array.shape[-1])
    return np.ascontiguousarray(array, dtype=np.float32)


def bin_matrix(values: np.ndarray, bins: int) -> Dict[str, np.ndarray]:
    """Bin ``(n, tracks)`` values into ``bins`` columns: NaN-aware mean, min and max.

    Returns arrays shaped ``(tracks, bins)`` plus ``edges`` (bin start offsets, length bins+1).
    When ``n <= bins`` every position is its own bin (no smoothing).
    """
    n = values.shape[0]
    tracks = values.shape[1] if values.ndim == 2 else 1
    values = values.reshape(n, tracks)
    if n == 0:
        empty = np.zeros((tracks, 0), np.float32)
        return {"mean": empty, "min": empty, "max": empty, "edges": np.zeros(1, np.int64)}
    if n <= bins:
        t = np.ascontiguousarray(values.T, dtype=np.float32)
        return {"mean": t, "min": t, "max": t, "edges": np.arange(n + 1, dtype=np.int64)}
    edges = np.unique(np.linspace(0, n, bins + 1).astype(np.int64))
    starts = edges[:-1]
    finite = np.isfinite(values)
    with np.errstate(invalid="ignore", divide="ignore"):
        sums = np.add.reduceat(np.where(finite, values, 0.0).astype(np.float64), starts, axis=0)
        counts = np.add.reduceat(finite.astype(np.int32), starts, axis=0)
        mean = (sums / counts).astype(np.float32)
        vmin = np.minimum.reduceat(np.where(finite, values, np.inf), starts, axis=0)
        vmax = np.maximum.reduceat(np.where(finite, values, -np.inf), starts, axis=0)
    empty = counts == 0
    vmin[empty] = np.nan
    vmax[empty] = np.nan
    return {
        "mean": np.ascontiguousarray(mean.T),
        "min": np.ascontiguousarray(vmin.T.astype(np.float32)),
        "max": np.ascontiguousarray(vmax.T.astype(np.float32)),
        "edges": edges,
    }


def _track_key(record: Dict[str, object], with_name: bool) -> Tuple:
    base = (str(record.get("ontology_curie") or record.get("biosample_name")), str(record.get("strand")))
    return base + ((str(record.get("name") or record.get("Assay title") or ""),) if with_name else ())


def match_track_columns(stored: List[Dict[str, object]], canonical: List[Dict[str, object]]) -> List[int]:
    """Column of ``stored`` holding each ``canonical`` track (same ontology / strand / assay), or -1."""
    has_name = lambda rows: all(r.get("name") or r.get("Assay title") for r in rows)  # noqa: E731
    with_name = has_name(stored) and has_name(canonical)
    lookup: Dict[Tuple, List[int]] = {}
    for i, record in enumerate(stored):
        lookup.setdefault(_track_key(record, with_name), []).append(i)
    used: Dict[Tuple, int] = {}
    columns = []
    for record in canonical:
        key = _track_key(record, with_name)
        options = lookup.get(key) or []
        k = used.get(key, 0)
        columns.append(options[k] if k < len(options) else -1)
        used[key] = k + 1
    return columns


@dataclass(frozen=True)
class SeriesSpec:
    sample: str
    haplotype: str  # H1, H2 or H1+H2 (ignored for the reference)

    @property
    def label(self) -> str:
        return "Reference genome" if self.sample == REFERENCE_SAMPLE else f"{self.sample} {self.haplotype}"


def parse_series(text: str) -> List[SeriesSpec]:
    specs = []
    for raw in text.split(","):
        raw = raw.strip()
        if not raw:
            continue
        sample, _, hap = raw.partition(":")
        specs.append(SeriesSpec(sample.strip(), (hap.strip() or "H1").upper()))
    return specs


class SignalService:
    def __init__(self, cache_bytes: int, disk_cache: DiskArrayCache, alignment=None, workers: int = 8):
        self.arrays = LRUCache(max(cache_bytes * 3 // 4, 64 << 20))
        self.small = LRUCache(max(cache_bytes // 8, 32 << 20))  # variants, maps, sequences
        self.results = LRUCache(max(cache_bytes // 8, 32 << 20))
        self.disk = disk_cache
        self.alignment = alignment
        self.workers = max(1, workers)
        # Bulk jobs are dominated by zlib inflation (GIL-free), so they use more threads.
        self.bulk_workers = max(self.workers, min(16, os.cpu_count() or 4))

    # -- primitives ---------------------------------------------------------------
    def prediction(self, dataset: Dataset, sample: str, gene: str, haplotype: str, output: str, cache: bool = True) -> np.ndarray:
        path = dataset.prediction_path(sample, gene, haplotype, output)
        if not cache:
            hit = self.arrays.get(("pred", str(path)))
            return hit if hit is not None else load_prediction_matrix(path)
        return self.arrays.get_or_load(("pred", str(path)), lambda: load_prediction_matrix(path))

    def reference_prediction(self, dataset: Dataset, gene: str, output: str) -> np.ndarray:
        """Reference-window prediction with columns in the dataset's track order (NaN when absent)."""
        path = dataset.reference_prediction_path(gene, output)
        if not path.exists():
            raise FileNotFoundError(f"No AlphaGenome {output} prediction of the reference window for {gene} (predict haplotype 'ref')")

        def load() -> np.ndarray:
            matrix = load_prediction_matrix(path)
            stored = dataset.reference_outputs(gene).get(output) or []
            canonical = [t.get("metadata") or {} for t in (dataset.gene_info(gene)["outputs"].get(output) or {}).get("tracks", [])]
            if not canonical or not stored:
                return matrix
            columns = match_track_columns(stored, canonical)
            out = np.full((matrix.shape[0], len(columns)), np.nan, dtype=np.float32)
            for j, c in enumerate(columns):
                if 0 <= c < matrix.shape[1]:
                    out[:, j] = matrix[:, c]
            return out

        return self.arrays.get_or_load(("refpred", str(path), path.stat().st_mtime_ns), load)

    def reference_indexed_window(self, dataset: Dataset, gene: str, matrix: np.ndarray, res: int, coords: str, start: int, end: int) -> np.ndarray:
        """``(end-start, tracks)`` of a matrix indexed by reference offset (rows of ``res`` bases) in ``coords``."""
        out = np.full((end - start, matrix.shape[1]), np.nan, dtype=np.float32)
        positions = matrix.shape[0] * res
        if coords in ("reference", "haplotype"):
            lo, hi = max(start, 0), min(end, positions)
            if hi > lo:
                out[lo - start:hi - start] = matrix[lo:hi] if res == 1 else matrix[np.arange(lo, hi) // res]
            return out
        if coords == "aligned":
            axis = self._alignment().axis(dataset, gene)
            slots = np.asarray(axis["insertion_slots"], dtype=np.int64)
            pos = np.arange(start, end, dtype=np.int64)
            k = np.searchsorted(slots, pos)
            is_slot = (k < slots.size) & (slots[np.minimum(k, max(slots.size - 1, 0))] == pos) if slots.size else np.zeros(pos.size, bool)
            ref = int(axis["ref_start_offset"]) + pos - k
            valid = ~is_slot & (ref >= 0) & (ref < positions)
            out[np.nonzero(valid)[0]] = matrix[ref[valid] // res]
            return out
        raise ValueError(f"Unknown coordinate system: {coords}")

    def _parse_vcf(self, dataset: Dataset, sample: str, gene: str) -> SampleVariants:
        path = dataset.sample_vcf_path(sample, gene)
        window = dataset.window(gene)
        if path is None or window.start is None:
            raise FileNotFoundError(f"No window VCF for {sample}/{gene}; reference coordinates unavailable")
        parsed = parse_window_vcf(path, window.start)
        self.small.put(("events", dataset.id, sample, gene), parsed.haplotypes)
        return parsed

    def variants(self, dataset: Dataset, sample: str, gene: str) -> SampleVariants:
        return self.small.get_or_load(("vcf", dataset.id, sample, gene), lambda: self._parse_vcf(dataset, sample, gene))

    def events(self, dataset: Dataset, sample: str, gene: str) -> Dict[str, HaplotypeEvents]:
        """Per-haplotype applied indel events (compact; cached separately from full records)."""
        return self.small.get_or_load(("events", dataset.id, sample, gene), lambda: self._parse_vcf(dataset, sample, gene).haplotypes)

    def ref_map(self, dataset: Dataset, sample: str, gene: str, haplotype: str, cache: bool = True) -> ReferenceMap:
        def load() -> ReferenceMap:
            events = self.events(dataset, sample, gene).get(haplotype) or HaplotypeEvents.empty()
            return reference_map(events, self.reference_length(dataset, gene))

        if not cache:
            hit = self.small.get(("refmap", dataset.id, sample, gene, haplotype))
            return hit if hit is not None else load()
        return self.small.get_or_load(("refmap", dataset.id, sample, gene, haplotype), load)

    def reference_sequence(self, dataset: Dataset, gene: str) -> bytes:
        path = dataset.reference_dir(gene) / "ref.window.fa"
        return self.small.get_or_load(("refseq", str(path)), lambda: read_fasta_bytes(path))

    def haplotype_sequence(self, dataset: Dataset, sample: str, gene: str, haplotype: str) -> bytes:
        path = dataset.haplotype_fasta(sample, gene, haplotype)
        return self.small.get_or_load(("hapseq", str(path)), lambda: read_fasta_bytes(path))

    def reference_length(self, dataset: Dataset, gene: str) -> int:
        window = dataset.window(gene)
        if window.length:
            return int(window.length)
        return len(self.reference_sequence(dataset, gene))

    @staticmethod
    def resolution(dataset: Dataset, gene: str, output: str) -> int:
        """Bases per stored row of ``output`` (1 for per-base tracks)."""
        return int((dataset.gene_info(gene)["outputs"].get(output) or {}).get("resolution") or 1)

    def domain_length(self, dataset: Dataset, gene: str, output: str, coords: str) -> int:
        if coords == "reference":
            return self.reference_length(dataset, gene)
        if coords == "aligned":
            return int(self._alignment().axis(dataset, gene)["expanded_length"])
        info = dataset.gene_info(gene)
        length = (info["outputs"].get(output) or {}).get("length")
        return int(length or self.reference_length(dataset, gene))

    def _alignment(self):
        if self.alignment is None:
            raise RuntimeError("Aligned coordinates are not available in this server")
        return self.alignment

    # -- coordinate mapping -------------------------------------------------------
    def haplotype_window(
        self,
        dataset: Dataset,
        sample: str,
        gene: str,
        haplotype: str,
        output: str,
        coords: str,
        start: int,
        end: int,
        tracks: Optional[Sequence[int]] = None,
        cache: bool = True,
    ) -> np.ndarray:
        """``(end-start, tracks)`` values for one haplotype in the requested coordinates."""
        if sample == REFERENCE_SAMPLE:
            matrix = self.reference_prediction(dataset, gene, output)
            if tracks is not None and list(tracks) != list(range(matrix.shape[1])):
                matrix = matrix[:, list(tracks)]
            return self.reference_indexed_window(dataset, gene, matrix, self.resolution(dataset, gene, output), coords, start, end)
        if haplotype == DIPLOID:
            h1 = self.haplotype_window(dataset, sample, gene, "H1", output, coords, start, end, tracks, cache)
            h2 = self.haplotype_window(dataset, sample, gene, "H2", output, coords, start, end, tracks, cache)
            f1, f2 = np.isfinite(h1), np.isfinite(h2)
            total = np.where(f1, h1, 0.0) + np.where(f2, h2, 0.0)
            count = f1.astype(np.float32) + f2
            with np.errstate(invalid="ignore", divide="ignore"):
                return (total / count).astype(np.float32)
        matrix = self.prediction(dataset, sample, gene, haplotype, output, cache=cache)
        if tracks is not None and list(tracks) != list(range(matrix.shape[1])):
            matrix = matrix[:, list(tracks)]
        res = self.resolution(dataset, gene, output)
        positions = matrix.shape[0] * res  # haplotype bases covered by the stored rows
        if coords == "haplotype":
            out = np.full((end - start, matrix.shape[1]), np.nan, dtype=np.float32)
            lo, hi = max(start, 0), min(end, positions)
            if hi > lo:
                out[lo - start:hi - start] = matrix[lo:hi] if res == 1 else matrix[np.arange(lo, hi) // res]
            return out
        if coords == "reference":
            ref_map = self.ref_map(dataset, sample, gene, haplotype, cache=cache)
            local = ref_map.local[max(start, 0):min(end, ref_map.local.size)]
            out = np.full((end - start, matrix.shape[1]), np.nan, dtype=np.float32)
            valid = (local >= 0) & (local < positions)
            offset = max(start, 0) - start
            target = np.nonzero(valid)[0] + offset
            out[target] = matrix[local[valid] // res]
            return out
        if coords == "aligned":
            expanded, source = self._alignment().entry_arrays(dataset, gene, sample, haplotype)
            out = np.full((end - start, matrix.shape[1]), np.nan, dtype=np.float32)
            mask = (expanded >= start) & (expanded < end) & (source >= 0) & (source < positions)
            out[expanded[mask] - start] = matrix[source[mask] // res]
            return out
        raise ValueError(f"Unknown coordinate system: {coords}")

    # -- endpoints ------------------------------------------------------------------
    def series_payload(
        self,
        dataset: Dataset,
        gene: str,
        output: str,
        series: List[SeriesSpec],
        tracks: List[int],
        coords: str,
        start: int,
        end: int,
        bins: int,
    ) -> Dict[str, object]:
        if not series:
            raise ValueError("Select at least one sample")
        if len(series) > MAX_SERIES:
            raise ValueError(f"At most {MAX_SERIES} series per request")
        domain = self.domain_length(dataset, gene, output, coords)
        start, end = _clamp_range(start, end, domain)
        bins = max(1, min(int(bins), MAX_BINS))
        results = []

        def one(spec: SeriesSpec):
            try:
                values = self.haplotype_window(dataset, spec.sample, gene, spec.haplotype, output, coords, start, end, tracks)
                return spec, bin_matrix(values, bins), None
            except Exception as exc:
                return spec, None, str(exc)

        with ThreadPoolExecutor(max_workers=min(self.workers, len(series))) as pool:
            for spec, binned, error in pool.map(one, series):
                item: Dict[str, object] = {"sample": spec.sample, "haplotype": spec.haplotype, "label": spec.label}
                if error:
                    item["error"] = error
                else:
                    item.update({"mean": binned["mean"], "min": binned["min"], "max": binned["max"]})
                results.append(item)
        ok = [r for r in results if "error" not in r]
        if not ok:
            raise RuntimeError("; ".join(str(r["error"]) for r in results[:3]))
        edges = bin_matrix(np.zeros((end - start, 1), np.float32), bins)["edges"]
        return {
            "gene": gene,
            "output": output,
            "coords": coords,
            "tracks": tracks,
            "start": start,
            "end": end,
            "domain": domain,
            "edges": (edges + start).astype(np.int64),
            "series": results,
        }

    # -- cohort aggregates --------------------------------------------------------
    def group_aggregate_key(
        self, dataset: Dataset, gene: str, output: str, tracks: List[int], haplotypes: List[str], coords: str, samples: List[str]
    ) -> str:
        return stable_key(
            {
                "v": 2,
                "dataset": dataset.fingerprint,
                "gene": gene,
                "output": output,
                "tracks": tracks,
                "haplotypes": haplotypes,
                "coords": coords,
                "samples": stable_key(sorted(samples)),
                "axis": self._alignment().axis_key(dataset, gene) if coords == "aligned" else None,
            }
        )

    def cached_group_aggregate(self, key: str) -> Optional[Dict[str, np.ndarray]]:
        hit = self.results.get(("group", key))
        if hit is not None:
            return hit
        hit = self.disk.load("group_aggregates", key)
        if hit is not None:
            self.results.put(("group", key), hit)
        return hit

    def compute_group_aggregate(
        self,
        dataset: Dataset,
        gene: str,
        output: str,
        tracks: List[int],
        haplotypes: List[str],
        coords: str,
        samples: List[str],
        progress: ProgressFn,
    ) -> Dict[str, np.ndarray]:
        """Full-domain per-position mean/std/count over ``samples`` x ``haplotypes``."""
        key = self.group_aggregate_key(dataset, gene, output, tracks, haplotypes, coords, samples)
        cached = self.cached_group_aggregate(key)
        if cached is not None:
            return cached
        if not samples:
            raise ValueError("Group has no samples")
        domain = self.domain_length(dataset, gene, output, coords)
        if coords == "reference":
            self.ensure_variants(dataset, gene, samples, progress)
        sums = np.zeros((domain, len(tracks)), np.float64)
        sumsq = np.zeros((domain, len(tracks)), np.float64)
        counts = np.zeros((domain, len(tracks)), np.int32)
        # Lock striping by row chunk lets workers accumulate different chunks concurrently.
        bounds = np.linspace(0, domain, 17).astype(np.int64)
        chunk_locks = [threading.Lock() for _ in range(len(bounds) - 1)]
        lock = threading.Lock()
        jobs = [(s, h) for s in samples for h in haplotypes]
        done = [0]
        failures: List[str] = []

        def one(job_index: int) -> None:
            if is_cancelled(progress):
                return
            sample, hap = jobs[job_index]
            try:
                values = self.haplotype_window(dataset, sample, gene, hap, output, coords, 0, domain, tracks, cache=False)
            except Exception as exc:
                with lock:
                    failures.append(f"{sample} {hap}: {exc}")
                    done[0] += 1
                return
            finite = np.isfinite(values)
            clean = np.where(finite, values, np.float32(0.0))
            n_chunks = len(chunk_locks)
            for step in range(n_chunks):
                c = (job_index + step) % n_chunks
                a, b = int(bounds[c]), int(bounds[c + 1])
                part = clean[a:b]
                with chunk_locks[c]:
                    np.add(sums[a:b], part, out=sums[a:b])
                    np.add(sumsq[a:b], np.square(part, dtype=np.float64), out=sumsq[a:b])
                    np.add(counts[a:b], finite[a:b], out=counts[a:b])
            with lock:
                done[0] += 1
                if done[0] % 8 == 0 or done[0] == len(jobs):
                    progress(0.1 + 0.9 * done[0] / len(jobs), f"Aggregated {done[0]}/{len(jobs)} haplotypes")

        progress(0.1, f"Aggregating {len(jobs)} haplotypes")
        with ThreadPoolExecutor(max_workers=self.bulk_workers) as pool:
            list(pool.map(one, range(len(jobs))))
        if is_cancelled(progress):
            raise JobCancelled()
        if len(failures) == len(jobs):
            raise RuntimeError(f"No haplotype could be loaded: {failures[0]}")
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = sums / counts
            var = np.maximum(sumsq / counts - mean ** 2, 0.0)
        result = {
            "mean": mean.astype(np.float32),
            "std": np.sqrt(var).astype(np.float32),
            "count": counts.astype(np.int32),
            "n_haplotypes": np.asarray([len(jobs) - len(failures)]),
            "n_failed": np.asarray([len(failures)]),
        }
        self.results.put(("group", key), result)
        self.disk.save("group_aggregates", key, result)
        return result

    def ensure_variants(self, dataset: Dataset, gene: str, samples: List[str], progress: ProgressFn) -> None:
        """Load indel events for many samples: disk index first, then parse the rest in parallel.

        Parsed events are persisted per (dataset, gene) and validated against each VCF's mtime, so
        reading a whole cohort's VCFs is a one-time cost per gene.
        """
        missing = [s for s in samples if self.small.get(("events", dataset.id, s, gene)) is None]
        if len(missing) < 16:
            return
        window = dataset.window(gene)
        if window.start is None:
            return
        disk_key = stable_key({"v": 1, "dataset": str(dataset.path), "gene": gene})
        stored = self.disk.load("indel_events", disk_key)
        index: Dict[str, Tuple[int, np.ndarray, np.ndarray]] = {}
        if stored is not None:
            bounds = np.concatenate([[0], np.cumsum(stored["counts"])])
            for i, sample in enumerate(stored["samples"].tolist()):
                lo, hi = int(bounds[i]), int(bounds[i + 1])
                index[sample] = (int(stored["mtimes"][i]), stored["events"][lo:hi], stored["haps"][lo:hi])
        todo: List[Tuple[str, str, int, int]] = []
        for sample in missing:
            path = dataset.sample_vcf_path(sample, gene)
            if path is None:
                continue
            mtime = path.stat().st_mtime_ns
            known = index.get(sample)
            if known is not None and known[0] == mtime:
                self.small.put(("events", dataset.id, sample, gene), unpack_events(known[1], known[2]))
            else:
                todo.append((sample, str(path), int(window.start), mtime))
        if not todo:
            return
        progress(0.0, f"Reading indels for {len(todo)} samples")
        mtimes = {sample: mtime for sample, _p, _s, mtime in todo}
        for sample, events, haps in self._parse_events_parallel([(s, p, w) for s, p, w, _m in todo], progress):
            if haps.size == 1 and haps[0] < 0:
                continue  # unreadable VCF; the per-sample path will report it
            self.small.put(("events", dataset.id, sample, gene), unpack_events(events, haps))
            index[sample] = (mtimes[sample], events, haps)
        names = sorted(index)
        self.disk.save(
            "indel_events",
            disk_key,
            {
                "samples": np.asarray(names),
                "mtimes": np.asarray([index[s][0] for s in names], dtype=np.int64),
                "counts": np.asarray([len(index[s][1]) for s in names], dtype=np.int64),
                "events": np.concatenate([index[s][1] for s in names]).astype(np.int64) if names else np.zeros((0, 3), np.int64),
                "haps": np.concatenate([index[s][2] for s in names]).astype(np.int8) if names else np.zeros(0, np.int8),
            },
        )

    def _parse_events_parallel(self, jobs: List[Tuple[str, str, int]], progress: ProgressFn):
        """Parse VCFs in a forkserver process pool (BGZF inflation thrashes the GIL in threads)."""
        chunks = [jobs[i:i + 24] for i in range(0, len(jobs), 24)]
        results = []
        try:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor, as_completed

            context = multiprocessing.get_context("forkserver")
            context.set_forkserver_preload([])  # do not re-import the caller's __main__
            workers = max(2, min(16, (os.cpu_count() or 4) - 2, len(chunks)))
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
                futures = [pool.submit(parse_events_batch, chunk) for chunk in chunks]
                for i, future in enumerate(as_completed(futures), start=1):
                    results.extend(future.result())
                    progress(0.1 * i / len(chunks), f"Read indels for {min(i * 24, len(jobs))}/{len(jobs)} samples")
            return results
        except (ValueError, OSError, ImportError, RuntimeError) as exc:
            # No forkserver on this platform (or pool failure): parse in threads instead.
            print(f"visualizer: process pool unavailable ({exc}); parsing VCFs in threads", file=sys.stderr)
        with ThreadPoolExecutor(max_workers=self.bulk_workers) as pool:
            for i, batch in enumerate(pool.map(parse_events_batch, chunks), start=1):
                results.extend(batch)
                progress(0.1 * i / len(chunks), f"Read indels for {min(i * 24, len(jobs))}/{len(jobs)} samples")
        return results

    def group_payload(
        self,
        aggregates: Dict[str, Dict[str, np.ndarray]],
        sizes: Dict[str, int],
        tracks: List[int],
        domain: int,
        start: int,
        end: int,
        bins: int,
    ) -> Dict[str, object]:
        start, end = _clamp_range(start, end, domain)
        bins = max(1, min(int(bins), MAX_BINS))
        groups = []
        edges = None
        for name, agg in aggregates.items():
            mean = agg["mean"][start:end]
            std = agg["std"][start:end]
            binned = bin_matrix(mean, bins)
            upper = bin_matrix(mean + std, bins)["mean"]
            lower = bin_matrix(mean - std, bins)["mean"]
            edges = binned["edges"]
            groups.append(
                {
                    "group": name,
                    "samples": sizes.get(name, 0),
                    "haplotypes": int(agg["n_haplotypes"][0]),
                    "failed": int(agg["n_failed"][0]),
                    "mean": binned["mean"],
                    "min": binned["min"],
                    "max": binned["max"],
                    "upper": upper,
                    "lower": lower,
                }
            )
        if edges is None:
            edges = np.arange(1)
        return {"start": start, "end": end, "domain": domain, "tracks": tracks, "edges": (edges + start).astype(np.int64), "groups": groups}

    # -- population matrix -------------------------------------------------------
    def population_matrix(
        self,
        dataset: Dataset,
        gene: str,
        output: str,
        track: int,
        haplotype: str,
        coords: str,
        samples: List[str],
        start: int,
        end: int,
        bins: int,
        progress: ProgressFn,
    ) -> Dict[str, object]:
        domain = self.domain_length(dataset, gene, output, coords)
        start, end = _clamp_range(start, end, domain)
        bins = max(1, min(int(bins), 2048))
        key = ("popmatrix", dataset.fingerprint, gene, output, track, haplotype, coords, stable_key(samples), start, end, bins)
        hit = self.results.get(key)
        if hit is not None:
            return hit
        if coords == "reference":
            self.ensure_variants(dataset, gene, samples, progress)
        n_cols = min(bins, end - start)
        matrix = np.full((len(samples), n_cols), np.nan, dtype=np.float32)
        done = [0]
        lock = threading.Lock()
        failed: List[str] = []

        def one(i: int) -> None:
            if is_cancelled(progress):
                return
            sample = samples[i]
            try:
                values = self.haplotype_window(dataset, sample, gene, haplotype, output, coords, start, end, [track], cache=False)
                row = bin_matrix(values, bins)["mean"][0]
                matrix[i, : row.size] = row[:n_cols]
            except Exception:
                failed.append(sample)
            with lock:
                done[0] += 1
                if done[0] % 16 == 0:
                    progress(0.1 + 0.9 * done[0] / len(samples), f"Loaded {done[0]}/{len(samples)} samples")

        progress(0.1, f"Loading {len(samples)} samples")
        with ThreadPoolExecutor(max_workers=self.bulk_workers) as pool:
            list(pool.map(one, range(len(samples))))
        if is_cancelled(progress):
            raise JobCancelled()
        edges = bin_matrix(np.zeros((end - start, 1), np.float32), bins)["edges"]
        payload = {
            "start": start,
            "end": end,
            "domain": domain,
            "edges": (edges + start).astype(np.int64),
            "samples": samples,
            "matrix": matrix,
            "failed": failed[:50],
            "failed_count": len(failed),
        }
        self.results.put(key, payload)
        return payload


def _clamp_range(start: int, end: int, domain: int) -> Tuple[int, int]:
    start = max(0, min(int(start), max(domain - 1, 0)))
    end = max(start + 1, min(int(end), domain))
    return start, end
