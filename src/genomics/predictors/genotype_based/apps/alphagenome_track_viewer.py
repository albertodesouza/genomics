from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse

import numpy as np

from genomics.predictors.genotype_based.alignment.bcftools_chain_mapper import BcftoolsChainMapper
from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner
from genomics.workspace import DEFAULT_CONSENSUS_DATASET_DIR


DEFAULT_POINTS = 1000
MAX_POINTS = 20_000
MAX_SERIES = 24
DEFAULT_HEATMAP_TRACKS = 128
MAX_HEATMAP_TRACKS = 512
MAX_HEATMAP_SAMPLES = 120
DEFAULT_VIEW_LENGTH = 32768

PIGMENTATION_BY_POPULATION = {
    "YRI": "strong pigmentation",
    "ESN": "strong pigmentation",
    "LWK": "strong pigmentation",
    "MSL": "strong pigmentation",
    "GWD": "strong pigmentation",
    "FIN": "weak pigmentation",
    "CEU": "weak pigmentation",
    "GBR": "weak pigmentation",
}


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _safe_int(value: str, default: int, minimum: int, maximum: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return max(minimum, min(maximum, parsed))


def _json_scalar(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if not np.isfinite(value):
            return None
        return float(value)
    return value


def _sanitize_json(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): _sanitize_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_sanitize_json(v) for v in value]
    if isinstance(value, tuple):
        return [_sanitize_json(v) for v in value]
    return value


def _load_prediction_array(path: Path) -> Tuple[np.ndarray, str, List[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Prediction file not found: {path}")
    with np.load(path) as data:
        keys = list(data.files)
        if "values" in data:
            return np.asarray(data["values"]), "values", keys
        if "track_0" in data:
            track_keys = sorted(k for k in keys if k.startswith("track_") and k[6:].isdigit())
            if track_keys:
                tracks = [np.asarray(data[k]) for k in track_keys]
                if len(tracks) == 1:
                    return tracks[0], track_keys[0], keys
                return np.column_stack(tracks), "+".join(track_keys), keys
        array_keys = [k for k in keys if not k.endswith("_metadata")]
        if not array_keys:
            raise ValueError(f"No array payload found in {path}")
        key = array_keys[0]
        return np.asarray(data[key]), key, keys


def _extract_signal(array: np.ndarray, track: int, start: int, length: int) -> np.ndarray:
    if array.ndim == 0:
        signal = array.reshape(1)
    elif array.ndim == 1:
        signal = array
    else:
        tracks = array.shape[-1]
        if track < 0 or track >= tracks:
            raise ValueError(f"track={track} out of range for array shape {tuple(array.shape)}")
        signal = array.reshape((-1, tracks))[:, track]

    start_index = max(0, start - 1)
    end_index = min(signal.shape[0], start_index + length)
    if start_index >= signal.shape[0] or end_index <= start_index:
        return np.asarray([], dtype=np.float32)
    return np.asarray(signal[start_index:end_index], dtype=np.float32)


def _extract_full_signal(array: np.ndarray, track: int) -> np.ndarray:
    if array.ndim == 0:
        return array.reshape(1).astype(np.float32)
    if array.ndim == 1:
        return np.asarray(array, dtype=np.float32)
    tracks = array.shape[-1]
    if track < 0 or track >= tracks:
        raise ValueError(f"track={track} out of range for array shape {tuple(array.shape)}")
    return np.asarray(array.reshape((-1, tracks))[:, track], dtype=np.float32)


def _extract_matrix(array: np.ndarray, start: int, length: int, track_start: int, track_count: int) -> np.ndarray:
    if array.ndim == 0:
        matrix = array.reshape((1, 1))
    elif array.ndim == 1:
        matrix = array.reshape((-1, 1))
    else:
        tracks = array.shape[-1]
        matrix = array.reshape((-1, tracks))

    left_track = max(0, track_start)
    right_track = min(matrix.shape[1], left_track + max(1, track_count))
    if left_track >= matrix.shape[1] or right_track <= left_track:
        raise ValueError(f"track_start={track_start} out of range for array shape {tuple(array.shape)}")

    start_index = max(0, start - 1)
    end_index = min(matrix.shape[0], start_index + length)
    if start_index >= matrix.shape[0] or end_index <= start_index:
        return np.empty((0, right_track - left_track), dtype=np.float32)
    return np.asarray(matrix[start_index:end_index, left_track:right_track], dtype=np.float32)


def _extract_full_matrix(array: np.ndarray, track_start: int, track_count: int) -> np.ndarray:
    if array.ndim == 0:
        matrix = array.reshape((1, 1))
    elif array.ndim == 1:
        matrix = array.reshape((-1, 1))
    else:
        tracks = array.shape[-1]
        matrix = array.reshape((-1, tracks))

    left_track = max(0, track_start)
    right_track = min(matrix.shape[1], left_track + max(1, track_count))
    if left_track >= matrix.shape[1] or right_track <= left_track:
        raise ValueError(f"track_start={track_start} out of range for array shape {tuple(array.shape)}")
    return np.asarray(matrix[:, left_track:right_track], dtype=np.float32)


def _align_signal_window(
    signal: np.ndarray,
    entry: Dict[str, Any],
    expanded_length: int,
    start: int,
    length: int,
) -> np.ndarray:
    end = min(expanded_length, start + length - 1)
    if end < start:
        return np.asarray([], dtype=np.float32)
    window = np.full(end - start + 1, np.nan, dtype=np.float32)
    copy_from = entry.get("copy_from_indices", [])
    copy_to = entry.get("expanded_indices", [])
    source_start = int(entry.get("source_start_idx", 0))
    absolute_source = str(entry.get("mapping_method")) == "bcftools_chain"
    for source_raw, target_raw in zip(copy_from, copy_to):
        source = int(source_raw) if absolute_source else source_start + int(source_raw)
        target = int(target_raw)
        if start <= target + 1 <= end and 0 <= source < signal.size:
            window[target + 1 - start] = signal[source]
    return window


def _align_matrix_window(
    matrix: np.ndarray,
    entry: Dict[str, Any],
    expanded_length: int,
    start: int,
    length: int,
) -> np.ndarray:
    end = min(expanded_length, start + length - 1)
    if end < start:
        return np.empty((0, matrix.shape[1]), dtype=np.float32)
    window = np.full((end - start + 1, matrix.shape[1]), np.nan, dtype=np.float32)
    copy_from = entry.get("copy_from_indices", [])
    copy_to = entry.get("expanded_indices", [])
    source_start = int(entry.get("source_start_idx", 0))
    absolute_source = str(entry.get("mapping_method")) == "bcftools_chain"
    for source_raw, target_raw in zip(copy_from, copy_to):
        source = int(source_raw) if absolute_source else source_start + int(source_raw)
        target = int(target_raw)
        if start <= target + 1 <= end and 0 <= source < matrix.shape[0]:
            window[target + 1 - start, :] = matrix[source, :]
    return window


def _downsample(signal: np.ndarray, start: int, points: int) -> Tuple[List[int], List[Optional[float]], str]:
    if signal.size == 0:
        return [], [], "empty"
    if signal.size <= points:
        x = list(range(start, start + int(signal.size)))
        y = [_json_scalar(v) for v in signal]
        return x, y, "none"

    bucket_count = max(1, points)
    edges = np.linspace(0, signal.size, num=bucket_count + 1, dtype=np.int64)
    xs: List[int] = []
    ys: List[Optional[float]] = []
    for i in range(bucket_count):
        left = int(edges[i])
        right = int(edges[i + 1])
        if right <= left:
            continue
        bucket = signal[left:right]
        if bucket.size == 0:
            continue
        finite = bucket[np.isfinite(bucket)]
        xs.append(start + (left + right - 1) // 2)
        ys.append(None if finite.size == 0 else _json_scalar(np.mean(finite)))
    return xs, ys, "mean"


def _downsample_matrix(matrix: np.ndarray, start: int, points: int) -> Tuple[List[int], np.ndarray, str]:
    if matrix.size == 0 or matrix.shape[0] == 0:
        return [], np.empty((matrix.shape[1] if matrix.ndim == 2 else 0, 0), dtype=np.float32), "empty"
    if matrix.shape[0] <= points:
        x = list(range(start, start + int(matrix.shape[0])))
        return x, matrix.T.astype(np.float32, copy=False), "none"

    bucket_count = max(1, points)
    edges = np.linspace(0, matrix.shape[0], num=bucket_count + 1, dtype=np.int64)
    xs: List[int] = []
    values = np.full((matrix.shape[1], bucket_count), np.nan, dtype=np.float32)
    used = 0
    for i in range(bucket_count):
        left = int(edges[i])
        right = int(edges[i + 1])
        if right <= left:
            continue
        bucket = matrix[left:right, :]
        if bucket.size == 0:
            continue
        xs.append(start + (left + right - 1) // 2)
        with np.errstate(invalid="ignore"):
            finite_count = np.sum(np.isfinite(bucket), axis=0)
            sums = np.nansum(bucket, axis=0)
            means = np.divide(sums, finite_count, out=np.full(bucket.shape[1], np.nan, dtype=np.float32), where=finite_count > 0)
        values[:, used] = means.astype(np.float32)
        used += 1
    return xs, values[:, :used], "mean"


class AlphaGenomeRepository:
    def __init__(self, dataset_dir: Path, consensus_dataset_dir: Optional[Path] = None):
        self.dataset_dir = Path(dataset_dir).resolve()
        self.consensus_dataset_dir = Path(consensus_dataset_dir).resolve() if consensus_dataset_dir else DEFAULT_CONSENSUS_DATASET_DIR
        self.cache_dir = self.dataset_dir / ".alphagenome_track_viewer_cache"
        metadata_path = self.dataset_dir / "dataset_metadata.json"
        if not metadata_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json not found: {metadata_path}")
        self.metadata = _load_json(metadata_path)
        self.individuals = self._load_individuals()
        self.individual_rows = self._load_individual_rows()
        self.genes = self._load_genes()
        self.outputs = self._load_outputs()
        self._aligner: Optional[DynamicIndelAligner] = None
        self._chain_mapper: Optional[BcftoolsChainMapper] = None

    def _cache_path(self, namespace: str, payload: Dict[str, Any]) -> Path:
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        digest = hashlib.sha256(encoded).hexdigest()
        return self.cache_dir / namespace / f"{digest}.json.gz"

    def _load_cache(self, path: Path) -> Optional[Dict[str, Any]]:
        if not path.exists():
            return None
        try:
            with gzip.open(path, "rt", encoding="utf-8") as f:
                payload = json.load(f)
            if isinstance(payload, dict):
                payload["cache"] = {"hit": True, "path": str(path)}
                return payload
        except Exception:
            return None
        return None

    def _save_cache(self, path: Path, payload: Dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = path.with_suffix(path.suffix + ".tmp")
        cache_payload = dict(payload)
        cache_payload["cache"] = {"hit": False, "path": str(path)}
        with gzip.open(tmp_path, "wt", encoding="utf-8") as f:
            json.dump(cache_payload, f, allow_nan=False, separators=(",", ":"))
        tmp_path.replace(path)

    def _load_individuals(self) -> List[str]:
        raw = self.metadata.get("individuals", [])
        if raw and isinstance(raw[0], dict):
            values = [str(item.get("sample_id") or item.get("id") or "") for item in raw]
        else:
            values = [str(item) for item in raw]
        return sorted(v for v in values if v)

    def _load_individual_rows(self) -> List[Dict[str, Any]]:
        pedigree = self.metadata.get("individuals_pedigree", {}) or {}
        if not isinstance(pedigree, dict):
            pedigree = {}
        rows = []
        for sample_id in self.individuals:
            meta = pedigree.get(sample_id, {}) or {}
            if not isinstance(meta, dict):
                meta = {}
            rows.append({
                "sample_id": sample_id,
                "population": str(meta.get("population", "")),
                "superpopulation": str(meta.get("superpopulation", "")),
                "pigmentation": str(meta.get("pigmentation", "")),
                "sex": meta.get("sex_label", meta.get("sex", "")),
                "family_id": str(meta.get("family_id", "")),
            })
        return rows

    def _load_genes(self) -> List[str]:
        raw = self.metadata.get("genes", [])
        if raw and isinstance(raw[0], dict):
            values = [str(item.get("gene") or item.get("name") or item.get("id") or "") for item in raw]
        else:
            values = [str(item) for item in raw]
        return sorted(v for v in values if v)

    def _load_outputs(self) -> List[str]:
        outputs = [str(v).split(".")[-1].lower() for v in self.metadata.get("alphagenome_outputs", []) if str(v)]
        if outputs:
            return sorted(set(outputs))

        first_sample = self.individuals[0] if self.individuals else None
        first_gene = self.genes[0] if self.genes else None
        if not first_sample or not first_gene:
            return []
        pred_dir = self.dataset_dir / "individuals" / first_sample / "windows" / first_gene / "predictions_H1"
        if not pred_dir.exists():
            return []
        return sorted(path.stem for path in pred_dir.glob("*.npz"))

    def _default_view_window(self) -> Dict[str, int]:
        window_size = int(self.metadata.get("window_size") or 524288)
        length = min(DEFAULT_VIEW_LENGTH, window_size)
        start = max(1, (window_size - length) // 2 + 1)
        return {"start": start, "length": length}

    def _distribution(self, field: str) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for row in self.individual_rows:
            key = str(row.get(field) or "")
            if not key:
                continue
            counts[key] = counts.get(key, 0) + 1
        return dict(sorted(counts.items()))

    def _class_label_for_row(self, row: Dict[str, Any], class_field: str) -> str:
        if class_field == "pigmentation":
            explicit = str(row.get("pigmentation") or "")
            if explicit:
                return explicit
            return PIGMENTATION_BY_POPULATION.get(str(row.get("population") or ""), "")
        return str(row.get(class_field) or "")

    def _class_distribution(self, class_field: str) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for row in self.individual_rows:
            key = self._class_label_for_row(row, class_field)
            if not key:
                continue
            counts[key] = counts.get(key, 0) + 1
        return dict(sorted(counts.items()))

    def _samples_for_class(self, class_field: str, class_value: str) -> List[str]:
        return [
            str(row["sample_id"])
            for row in self.individual_rows
            if self._class_label_for_row(row, class_field) == class_value
        ]

    def options_payload(self) -> Dict[str, Any]:
        default_window = self._default_view_window()
        return {
            "dataset_dir": str(self.dataset_dir),
            "dataset_name": self.metadata.get("dataset_name"),
            "samples": self.individuals,
            "individuals": self.individual_rows,
            "populations": self._distribution("population"),
            "superpopulations": self._distribution("superpopulation"),
            "pigmentations": self._class_distribution("pigmentation"),
            "genes": self.genes,
            "haplotypes": ["H1", "H2"],
            "outputs": self.outputs,
            "defaults": {
                "sample": self.individuals[0] if self.individuals else None,
                "gene": self.genes[0] if self.genes else None,
                "haplotype": "H1",
                "output": self.outputs[0] if self.outputs else "rna_seq",
                "track": 0,
                "start": default_window["start"],
                "length": default_window["length"],
                "points": DEFAULT_POINTS,
            },
        }

    def _prediction_dir(self, sample: str, gene: str, haplotype: str) -> Path:
        if haplotype not in {"H1", "H2"}:
            raise ValueError("haplotype must be H1 or H2")
        return self.dataset_dir / "individuals" / sample / "windows" / gene / f"predictions_{haplotype}"

    def _load_track_metadata(self, pred_dir: Path, output: str) -> Optional[List[Dict[str, Any]]]:
        path = pred_dir / f"{output}_metadata.json"
        if not path.exists():
            return None
        payload = _load_json(path)
        metadata = payload.get("metadata")
        return metadata if isinstance(metadata, list) else None

    def track_options_payload(self, sample: str, gene: str, haplotype: str, output: str) -> Dict[str, Any]:
        pred_dir = self._prediction_dir(sample, gene, haplotype)
        npz_path = pred_dir / f"{output}.npz"
        array, _array_key, _npz_keys = _load_prediction_array(npz_path)
        metadata = self._load_track_metadata(pred_dir, output) or []
        track_count = 1 if array.ndim <= 1 else int(array.shape[-1])
        tracks = []
        for idx in range(track_count):
            meta = metadata[idx] if idx < len(metadata) else {}
            label = self._track_label(idx, meta)
            tracks.append({"index": idx, "label": label, "metadata": meta})
        return {"track_count": track_count, "tracks": tracks, "array_shape": list(array.shape), "array_dtype": str(array.dtype)}

    def _track_label(self, idx: int, meta: Optional[Dict[str, Any]]) -> str:
        if not meta:
            return f"Track {idx}"
        bits = []
        if meta.get("biosample_name"):
            bits.append(str(meta["biosample_name"]))
        elif meta.get("name"):
            bits.append(str(meta["name"]))
        if meta.get("ontology_curie"):
            bits.append(str(meta["ontology_curie"]))
        if meta.get("strand"):
            bits.append(f"strand {meta['strand']}")
        if meta.get("Assay title"):
            bits.append(str(meta["Assay title"]))
        return f"{idx}: " + " | ".join(bits)

    def _load_individual_metadata_summary(self, sample_id: str) -> Dict[str, Any]:
        base = next((row for row in self.individual_rows if row["sample_id"] == sample_id), {"sample_id": sample_id})
        metadata_path = self.dataset_dir / "individuals" / sample_id / "individual_metadata.json"
        summary = dict(base)
        if metadata_path.exists():
            try:
                meta = _load_json(metadata_path)
                for key in ("longevity", "family_id", "created_at", "last_updated"):
                    if key in meta:
                        summary[key] = meta[key]
                if isinstance(meta.get("frog_likelihood"), list) and meta["frog_likelihood"]:
                    vals = [float(v) for v in meta["frog_likelihood"] if isinstance(v, (int, float)) and math.isfinite(float(v))]
                    if vals:
                        summary["frog_likelihood_max"] = max(vals)
                        summary["frog_likelihood_argmax"] = vals.index(max(vals))
                if isinstance(meta.get("windows"), list):
                    summary["window_count"] = len(meta["windows"])
            except Exception as exc:
                summary["metadata_error"] = str(exc)
        return summary

    def _get_aligner(self) -> DynamicIndelAligner:
        if self._aligner is None:
            self._aligner = DynamicIndelAligner(self.dataset_dir, selected_sample_ids=self.individuals, center_window_size=DEFAULT_VIEW_LENGTH)
        return self._aligner

    def _get_chain_mapper(self) -> BcftoolsChainMapper:
        aligner = self._get_aligner()
        if self._chain_mapper is None:
            self._chain_mapper = BcftoolsChainMapper(
                dataset_dir=self.dataset_dir,
                consensus_dataset_dir=self.consensus_dataset_dir,
                aligner=aligner,
            )
        return self._chain_mapper

    def model_window_payload(self, gene: str, align: bool) -> Dict[str, Any]:
        if not align:
            default_window = self._default_view_window()
            return {
                "gene": gene,
                "align": False,
                "start": default_window["start"],
                "length": default_window["length"],
                "x_axis": "sample_sequence",
                "model_window": None,
            }
        aligner = self._get_aligner()
        aligner.build_alignment_axis_for_gene(gene, self.individuals)
        center_slice = aligner.get_reference_centered_expanded_slice(gene, DEFAULT_VIEW_LENGTH)
        return {
            "gene": gene,
            "align": True,
            "start": int(center_slice["expanded_start"]) + 1,
            "length": int(center_slice["expanded_end"]) - int(center_slice["expanded_start"]),
            "x_axis": "bcftools_chain_expanded_alignment",
            "model_window": center_slice,
        }

    def _effective_window(
        self,
        gene: str,
        start: int,
        length: int,
        aligner: Optional[DynamicIndelAligner],
        expanded_length: Optional[int],
    ) -> Tuple[int, int, Optional[Dict[str, int]]]:
        if aligner is None:
            return start, length, None
        center_slice = aligner.get_reference_centered_expanded_slice(gene, DEFAULT_VIEW_LENGTH)
        effective_start = int(center_slice["expanded_start"]) + 1
        effective_length = int(center_slice["expanded_end"]) - int(center_slice["expanded_start"])
        if expanded_length is not None:
            effective_start = max(1, min(effective_start, int(expanded_length or 1)))
            effective_length = max(1, min(effective_length, int(expanded_length or 1) - effective_start + 1))
        return effective_start, effective_length, center_slice

    def track_payload(
        self,
        sample: str,
        gene: str,
        haplotype: str,
        output: str,
        track: int,
        start: int,
        length: int,
        points: int,
    ) -> Dict[str, Any]:
        pred_dir = self._prediction_dir(sample, gene, haplotype)
        npz_path = pred_dir / f"{output}.npz"
        array, array_key, npz_keys = _load_prediction_array(npz_path)
        metadata = self._load_track_metadata(pred_dir, output)
        signal = _extract_signal(array, track=track, start=start, length=length)
        finite = signal[np.isfinite(signal)]
        x, y, method = _downsample(signal, start=max(1, start), points=points)

        track_count = 1 if array.ndim <= 1 else int(array.shape[-1])
        track_meta = metadata[track] if metadata and 0 <= track < len(metadata) else None
        return {
            "sample": sample,
            "gene": gene,
            "haplotype": haplotype,
            "output": output,
            "track": track,
            "track_count": track_count,
            "track_metadata": track_meta,
            "array_shape": list(array.shape),
            "array_dtype": str(array.dtype),
            "array_key": array_key,
            "npz_keys": npz_keys,
            "source_file": str(npz_path),
            "requested": {"start": start, "length": length, "points": points},
            "returned_points": len(y),
            "downsample": method,
            "stats": {
                "min": None if finite.size == 0 else _json_scalar(np.min(finite)),
                "max": None if finite.size == 0 else _json_scalar(np.max(finite)),
                "mean": None if finite.size == 0 else _json_scalar(np.mean(finite)),
            },
            "x": x,
            "y": y,
        }

    def compare_payload(
        self,
        samples: List[str],
        gene: str,
        haplotypes: List[str],
        output: str,
        track: int,
        start: int,
        length: int,
        points: int,
        align: bool,
    ) -> Dict[str, Any]:
        samples = [sample for sample in samples if sample]
        haplotypes = [hap for hap in haplotypes if hap in {"H1", "H2"}]
        if not samples:
            raise ValueError("Select at least one sample")
        if not haplotypes:
            raise ValueError("Select at least one haplotype")
        if len(samples) * len(haplotypes) > MAX_SERIES:
            raise ValueError(f"Too many series: max {MAX_SERIES} sample/haplotype combinations")

        aligner = self._get_aligner() if align else None
        chain_mapper = self._get_chain_mapper() if align else None
        if aligner is not None:
            aligner.build_alignment_axis_for_gene(gene, self.individuals)
        expanded_length = aligner.get_expanded_length(gene) if aligner else None
        effective_start, effective_length, model_window = self._effective_window(gene, start, length, aligner, expanded_length)
        series = []
        array_shape: Optional[List[int]] = None
        array_dtype: Optional[str] = None
        track_count: Optional[int] = None
        track_meta = None
        individual_summaries = [self._load_individual_metadata_summary(sample) for sample in samples]

        for sample in samples:
            for haplotype in haplotypes:
                pred_dir = self._prediction_dir(sample, gene, haplotype)
                npz_path = pred_dir / f"{output}.npz"
                array, array_key, npz_keys = _load_prediction_array(npz_path)
                metadata = self._load_track_metadata(pred_dir, output)
                full_signal = _extract_full_signal(array, track=track)
                if aligner is not None:
                    entry = chain_mapper.get_haplotype_entry(gene, sample, haplotype) if chain_mapper is not None else None
                    if entry is None:
                        raise RuntimeError(f"Alignment entry missing for {sample} {haplotype}")
                    signal = _align_signal_window(full_signal, entry, int(expanded_length), effective_start, effective_length)
                    x_start = effective_start
                else:
                    signal = _extract_signal(array, track=track, start=start, length=length)
                    x_start = max(1, start)
                finite = signal[np.isfinite(signal)]
                x, y, method = _downsample(signal, start=x_start, points=points)
                series.append({
                    "label": f"{sample}_{haplotype}",
                    "sample": sample,
                    "haplotype": haplotype,
                    "x": x,
                    "y": y,
                    "returned_points": len(y),
                    "downsample": method,
                    "stats": {
                        "min": None if finite.size == 0 else _json_scalar(np.min(finite)),
                        "max": None if finite.size == 0 else _json_scalar(np.max(finite)),
                        "mean": None if finite.size == 0 else _json_scalar(np.mean(finite)),
                    },
                    "source_file": str(npz_path),
                })
                array_shape = list(array.shape)
                array_dtype = str(array.dtype)
                track_count = 1 if array.ndim <= 1 else int(array.shape[-1])
                if metadata and 0 <= track < len(metadata):
                    track_meta = metadata[track]

        return {
            "gene": gene,
            "output": output,
            "track": track,
            "track_count": track_count,
            "track_metadata": track_meta,
            "array_shape": array_shape or [],
            "array_dtype": array_dtype,
            "requested": {"start": start, "length": length, "points": points, "align": align},
            "effective_window": {"start": effective_start, "length": effective_length},
            "model_window": model_window,
            "x_axis": "bcftools_chain_expanded_alignment" if align else "sample_sequence",
            "expanded_length": expanded_length,
            "individuals": individual_summaries,
            "series": series,
        }

    def heatmap_payload(
        self,
        mode: str,
        sample: str,
        class_field: str,
        class_value: str,
        gene: str,
        haplotypes: List[str],
        output: str,
        start: int,
        length: int,
        points: int,
        track_start: int,
        track_count: int,
        align: bool,
        max_samples: int,
    ) -> Dict[str, Any]:
        mode = mode if mode in {"sample", "class_mean"} else "sample"
        haplotypes = [hap for hap in haplotypes if hap in {"H1", "H2"}]
        if not haplotypes:
            raise ValueError("Select at least one haplotype")
        if track_count < 1:
            raise ValueError("track_count must be positive")

        samples = [sample] if mode == "sample" else self._samples_for_class(class_field, class_value)
        samples = [s for s in samples if s]
        if not samples:
            raise ValueError("No samples found for heatmap selection")
        truncated_samples = False
        if mode != "class_mean":
            max_samples = max(1, min(MAX_HEATMAP_SAMPLES, max_samples))
            truncated_samples = len(samples) > max_samples
            samples = samples[:max_samples]

        aligner = self._get_aligner() if align else None
        chain_mapper = self._get_chain_mapper() if align else None
        if aligner is not None:
            aligner.build_alignment_axis_for_gene(gene, self.individuals)
        expanded_length = aligner.get_expanded_length(gene) if aligner else None
        effective_start, effective_length, model_window = self._effective_window(gene, start, length, aligner, expanded_length)

        cache_path: Optional[Path] = None
        if mode == "class_mean":
            cache_key = {
                "version": 3,
                "dataset_dir": str(self.dataset_dir),
                "consensus_dataset_dir": str(self.consensus_dataset_dir),
                "mode": mode,
                "class_field": class_field,
                "class_value": class_value,
                "gene": gene,
                "haplotypes": haplotypes,
                "output": output,
                "effective_start": effective_start,
                "effective_length": effective_length,
                "model_window": model_window,
                "points": points,
                "track_start": track_start,
                "track_count": track_count,
                "align": align,
                "samples": samples,
            }
            cache_path = self._cache_path("heatmap_class_mean", cache_key)
            cached = self._load_cache(cache_path)
            if cached is not None:
                return cached

        matrices: List[np.ndarray] = []
        x_values: List[int] = []
        array_shape: Optional[List[int]] = None
        array_dtype: Optional[str] = None
        track_count_total: Optional[int] = None
        track_meta: List[Dict[str, Any]] = []
        missing: List[str] = []

        for current_sample in samples:
            for haplotype in haplotypes:
                pred_dir = self._prediction_dir(current_sample, gene, haplotype)
                npz_path = pred_dir / f"{output}.npz"
                try:
                    array, _array_key, _npz_keys = _load_prediction_array(npz_path)
                    metadata = self._load_track_metadata(pred_dir, output)
                    full_matrix = _extract_full_matrix(array, track_start=track_start, track_count=track_count)
                    if aligner is not None:
                        entry = chain_mapper.get_haplotype_entry(gene, current_sample, haplotype) if chain_mapper is not None else None
                        if entry is None:
                            missing.append(f"{current_sample}_{haplotype}: missing alignment")
                            continue
                        matrix = _align_matrix_window(full_matrix, entry, int(expanded_length), effective_start, effective_length)
                        x_start = effective_start
                    else:
                        matrix = _extract_matrix(array, start=start, length=length, track_start=track_start, track_count=track_count)
                        x_start = max(1, start)
                    current_x, sampled, _method = _downsample_matrix(matrix, start=x_start, points=points)
                    if not x_values:
                        x_values = current_x
                    matrices.append(sampled)
                    array_shape = list(array.shape)
                    array_dtype = str(array.dtype)
                    track_count_total = 1 if array.ndim <= 1 else int(array.shape[-1])
                    if metadata and not track_meta:
                        end_track = min(len(metadata), track_start + sampled.shape[0])
                        track_meta = metadata[track_start:end_track]
                except Exception as exc:
                    missing.append(f"{current_sample}_{haplotype}: {exc}")

        if not matrices:
            raise RuntimeError("No heatmap matrix could be loaded")

        min_rows = min(matrix.shape[0] for matrix in matrices)
        min_cols = min(matrix.shape[1] for matrix in matrices)
        stack = np.stack([matrix[:min_rows, :min_cols] for matrix in matrices], axis=0)
        with np.errstate(invalid="ignore"):
            finite_count = np.sum(np.isfinite(stack), axis=0)
            sums = np.nansum(stack, axis=0)
            heatmap = np.divide(sums, finite_count, out=np.full((min_rows, min_cols), np.nan, dtype=np.float32), where=finite_count > 0)
        finite = heatmap[np.isfinite(heatmap)]
        if finite.size:
            p99 = float(np.nanpercentile(finite, 99))
            vmax = p99 if p99 > 0 else float(np.nanmax(finite))
        else:
            vmax = 1.0
        if not math.isfinite(vmax) or vmax <= 0:
            vmax = 1.0

        track_labels = []
        for offset in range(min_rows):
            meta = track_meta[offset] if offset < len(track_meta) else {}
            track_labels.append(self._track_label(track_start + offset, meta))

        result = {
            "mode": mode,
            "class_field": class_field,
            "class_value": class_value,
            "sample": sample,
            "samples_used": samples,
            "samples_requested": len(self._samples_for_class(class_field, class_value)) if mode == "class_mean" else len(samples),
            "samples_truncated": truncated_samples,
            "haplotypes": haplotypes,
            "gene": gene,
            "output": output,
            "track_start": track_start,
            "track_rows": min_rows,
            "track_count_total": track_count_total,
            "track_labels": track_labels,
            "array_shape": array_shape or [],
            "array_dtype": array_dtype,
            "requested": {"start": start, "length": length, "points": points, "align": align, "track_count": track_count, "max_samples": None if mode == "class_mean" else max_samples},
            "effective_window": {"start": effective_start, "length": effective_length},
            "model_window": model_window,
            "x_axis": "bcftools_chain_expanded_alignment" if align else "sample_sequence",
            "expanded_length": expanded_length,
            "x": x_values[:min_cols],
            "matrix": [[_json_scalar(v) for v in row] for row in heatmap],
            "stats": {
                "min": None if finite.size == 0 else _json_scalar(np.nanmin(finite)),
                "max": None if finite.size == 0 else _json_scalar(np.nanmax(finite)),
                "mean": None if finite.size == 0 else _json_scalar(np.nanmean(finite)),
                "vmax": _json_scalar(vmax),
            },
            "missing": missing[:50],
            "missing_count": len(missing),
            "cache": {"hit": False, "path": str(cache_path) if cache_path else None},
        }
        if cache_path is not None:
            self._save_cache(cache_path, result)
        return result


class AlphaGenomeTrackViewerHandler(BaseHTTPRequestHandler):
    repository: AlphaGenomeRepository

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(_sanitize_json(payload), allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, text: str, content_type: str = "text/html; charset=utf-8", status: int = 200) -> None:
        body = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/":
                self._send_text(INDEX_HTML)
            elif parsed.path == "/api/options":
                self._send_json(self.repository.options_payload())
            elif parsed.path == "/api/tracks":
                self._handle_tracks(parsed.query)
            elif parsed.path == "/api/track-options":
                self._handle_track_options(parsed.query)
            elif parsed.path == "/api/heatmap":
                self._handle_heatmap(parsed.query)
            elif parsed.path == "/api/model-window":
                self._handle_model_window(parsed.query)
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _handle_tracks(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.options_payload()["defaults"]
        sample = qs.get("sample", [defaults["sample"] or ""])[0]
        samples = [item for raw in qs.get("samples", []) for item in raw.split(",") if item] or [sample]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        haplotype = qs.get("haplotype", ["H1"])[0]
        haplotypes = [item for raw in qs.get("haplotypes", []) for item in raw.split(",") if item] or [haplotype]
        output = qs.get("output", [defaults["output"] or "rna_seq"])[0]
        track = _safe_int(qs.get("track", ["0"])[0], 0, 0, 1_000_000)
        start = _safe_int(qs.get("start", ["1"])[0], 1, 1, 1_000_000_000)
        length = _safe_int(qs.get("length", ["2000"])[0], 2000, 1, 50_000_000)
        points = _safe_int(qs.get("points", [str(DEFAULT_POINTS)])[0], DEFAULT_POINTS, 1, MAX_POINTS)
        align = qs.get("align", ["false"])[0].lower() in {"1", "true", "yes", "on"}
        self._send_json(self.repository.compare_payload(samples, gene, haplotypes, output, track, start, length, points, align))

    def _handle_track_options(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.options_payload()["defaults"]
        sample = qs.get("sample", [defaults["sample"] or ""])[0]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        haplotype = qs.get("haplotype", ["H1"])[0]
        output = qs.get("output", [defaults["output"] or "rna_seq"])[0]
        self._send_json(self.repository.track_options_payload(sample, gene, haplotype, output))

    def _handle_model_window(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.options_payload()["defaults"]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        align = qs.get("align", ["false"])[0].lower() in {"1", "true", "yes", "on"}
        self._send_json(self.repository.model_window_payload(gene, align))

    def _handle_heatmap(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.options_payload()["defaults"]
        sample = qs.get("sample", [defaults["sample"] or ""])[0]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        haplotype = qs.get("haplotype", ["H1"])[0]
        haplotypes = [item for raw in qs.get("haplotypes", []) for item in raw.split(",") if item] or [haplotype]
        output = qs.get("output", [defaults["output"] or "rna_seq"])[0]
        mode = qs.get("mode", ["sample"])[0]
        class_field = qs.get("class_field", ["superpopulation"])[0]
        class_value = qs.get("class_value", [""])[0]
        start = _safe_int(qs.get("start", ["1"])[0], 1, 1, 1_000_000_000)
        length = _safe_int(qs.get("length", ["2000"])[0], 2000, 1, 50_000_000)
        points = _safe_int(qs.get("points", [str(DEFAULT_POINTS)])[0], DEFAULT_POINTS, 1, MAX_POINTS)
        track_start = _safe_int(qs.get("track_start", ["0"])[0], 0, 0, 1_000_000)
        track_count = _safe_int(qs.get("track_count", [str(DEFAULT_HEATMAP_TRACKS)])[0], DEFAULT_HEATMAP_TRACKS, 1, MAX_HEATMAP_TRACKS)
        max_samples = _safe_int(qs.get("max_samples", ["60"])[0], 60, 1, MAX_HEATMAP_SAMPLES)
        align = qs.get("align", ["false"])[0].lower() in {"1", "true", "yes", "on"}
        self._send_json(self.repository.heatmap_payload(
            mode=mode,
            sample=sample,
            class_field=class_field,
            class_value=class_value,
            gene=gene,
            haplotypes=haplotypes,
            output=output,
            start=start,
            length=length,
            points=points,
            track_start=track_start,
            track_count=track_count,
            align=align,
            max_samples=max_samples,
        ))


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>AlphaGenome Track Viewer</title>
  <style>
    :root { color-scheme: dark; --bg: #10141f; --panel: #172033; --line: #67d2ff; --muted: #9ca9bf; --text: #edf4ff; --warn: #ffd166; --bad: #ff6b6b; --ok: #4dd187; }
    * { box-sizing: border-box; }
    body { margin: 0; background: radial-gradient(circle at top left, #1d3150, var(--bg) 42rem); color: var(--text); font-family: system-ui, -apple-system, Segoe UI, sans-serif; }
    main { width: min(1920px, calc(100vw - 24px)); margin: 14px auto 30px; }
    h1 { margin: 0 0 4px; font-size: clamp(1.6rem, 3vw, 2.5rem); letter-spacing: -0.04em; }
    .sub { color: var(--muted); margin-bottom: 18px; }
    .panel { background: color-mix(in srgb, var(--panel), transparent 4%); border: 1px solid #2b3b59; border-radius: 18px; padding: 16px; box-shadow: 0 20px 50px #0007; }
    form { display: grid; grid-template-columns: repeat(8, minmax(0, 1fr)); gap: 12px; align-items: end; }
    label { display: grid; gap: 5px; color: var(--muted); font-size: .82rem; }
    select, input, button { width: 100%; border: 1px solid #344869; border-radius: 10px; background: #0c1220; color: var(--text); padding: 9px 10px; font: inherit; }
    button { cursor: pointer; background: linear-gradient(135deg, #1878ff, #21b6d7); border: 0; font-weight: 700; }
    .wide { grid-column: span 2; }
    .full { grid-column: 1 / -1; }
    .filter-grid { display:grid; grid-template-columns: repeat(3, minmax(0,1fr)); gap:12px; }
    .filter-panel { border:1px solid #344869; border-radius:14px; background:#0c1220; padding:10px; min-height:0; }
    .filter-panel h3 { margin:0 0 8px; font-size:.86rem; color:#dce8fb; display:flex; justify-content:space-between; }
    .check-list { max-height:180px; overflow:auto; margin-top:8px; }
    .check-list label { display:flex; grid-template-columns:none; align-items:center; gap:8px; color:#dce8fb; font-size:.82rem; padding:3px 0; }
    .check-list input { width:auto; }
    .button-row { display:grid; grid-template-columns: repeat(3, 1fr); gap:8px; margin-top:8px; }
    .button-row.two { grid-template-columns: repeat(2, 1fr); }
    .secondary { background:#22324a; border:1px solid #344869; }
    .sample-note { color: var(--muted); font-size: .78rem; }
    #status { min-height: 1.5rem; margin: 14px 0; color: var(--muted); }
    #status.error { color: var(--bad); }
    .checks { display:flex; gap:14px; align-items:center; flex-wrap:wrap; color:var(--text); }
    .checks label { display:flex; grid-template-columns:none; align-items:center; gap:6px; color:var(--text); }
    .checks input { width:auto; }
    .stats { display: grid; grid-template-columns: repeat(5, 1fr); gap: 10px; margin-bottom: 12px; }
    .stat { background: #0c1220; border: 1px solid #273854; border-radius: 14px; padding: 12px; min-width: 0; }
    .stat b { display: block; font-size: .78rem; color: var(--muted); font-weight: 600; }
    .stat span { display: block; margin-top: 4px; overflow-wrap: anywhere; }
    .viz-grid { display:grid; grid-template-columns:minmax(900px,1fr) 420px; gap:12px; align-items:stretch; }
    .plot-wrap { height: min(68vh, 760px); min-height: 460px; background: #080d17; border: 1px solid #273854; border-radius: 16px; overflow: hidden; }
    .heatmap-wrap { margin-top:12px; background:#080d17; border:1px solid #273854; border-radius:16px; padding:12px; }
    .heatmap-controls { display:grid; grid-template-columns: repeat(8, minmax(0, 1fr)); gap:10px; align-items:end; margin-bottom:10px; }
    .heatmap-canvas-wrap { height:min(62vh, 720px); min-height:420px; overflow:auto; border-radius:12px; border:1px solid #1d2c45; background:#02050a; }
    .heatmap-panels { display:grid; grid-template-columns: repeat(auto-fit, minmax(520px, 1fr)); gap:12px; margin-top:10px; }
    .heatmap-panel { background:#050914; border:1px solid #1d2c45; border-radius:14px; padding:10px; min-width:0; }
    .heatmap-panel h3 { margin:0 0 8px; font-size:.92rem; color:#dce8fb; overflow-wrap:anywhere; }
    canvas { display:block; image-rendering:pixelated; }
    .tabs { display:flex; gap:8px; margin:14px 0 10px; }
    .tab { width:auto; padding:9px 14px; background:#22324a; border:1px solid #344869; }
    .tab.active { background:linear-gradient(135deg, #1878ff, #21b6d7); border:0; }
    .hidden { display:none; }
    .side-panel { background:#0c1220; border:1px solid #273854; border-radius:16px; padding:12px; overflow:auto; max-height:min(68vh,760px); }
    .person { border-bottom:1px solid #273854; padding:8px 0; }
    .person b { color:#ffd166; }
    .kv { display:grid; grid-template-columns:110px 1fr; gap:5px; font-size:.82rem; color:#dbe7f3; }
    .kv span:first-child { color:var(--muted); }
    svg { display: block; width: 100%; height: 100%; }
    .meta { color: var(--muted); font-size: .88rem; margin-top: 12px; white-space: pre-wrap; overflow-wrap: anywhere; }
    @media (max-width: 1000px) { .viz-grid { grid-template-columns:1fr; } }
    @media (max-width: 900px) { form, .heatmap-controls { grid-template-columns: repeat(2, 1fr); } .wide { grid-column: span 2; } .stats { grid-template-columns: repeat(2, 1fr); } }
  </style>
</head>
<body>
<main>
  <h1>AlphaGenome Track Viewer</h1>
  <div class="sub">Compare AlphaGenome tracks across samples and haplotypes, optionally aligned with the bcftools_chain mapping used for training.</div>
  <section class="panel">
    <form id="controls">
      <div class="full filter-grid">
        <div class="filter-panel">
          <h3><span>Superpopulations</span><span id="superpopCount"></span></h3>
          <input id="superpopSearch" placeholder="Search superpopulation">
          <div class="button-row two"><button type="button" class="secondary" id="superpopAll">All</button><button type="button" class="secondary" id="superpopNone">None</button></div>
          <div class="check-list" id="superpopChecks"></div>
        </div>
        <div class="filter-panel">
          <h3><span>Populations</span><span id="populationCount"></span></h3>
          <input id="populationSearch" placeholder="Search population">
          <div class="button-row two"><button type="button" class="secondary" id="populationAll">All</button><button type="button" class="secondary" id="populationNone">None</button></div>
          <div class="check-list" id="populationChecks"></div>
        </div>
        <div class="filter-panel">
          <h3><span>Individuals</span><span id="sampleCount"></span></h3>
          <input id="sampleFilter" placeholder="Search sample, e.g. HG00096">
          <div class="button-row"><button type="button" class="secondary" id="sampleFirst5">First 5</button><button type="button" class="secondary" id="sampleVisible">Visible</button><button type="button" class="secondary" id="sampleNone">None</button></div>
          <div class="check-list" id="sampleChecks"></div>
          <span class="sample-note" id="sampleNote">Loading samples...</span>
        </div>
      </div>
      <label class="wide">Gene<select id="gene"></select></label>
      <label>Haplotypes<div class="checks"><label><input id="hapH1" type="checkbox" checked> H1</label><label><input id="hapH2" type="checkbox"> H2</label></div></label>
      <label>Output<select id="output"></select></label>
      <label class="wide">Track<select id="track"></select></label>
      <label>Start<input id="start" type="number" value="1" min="1" title="When aligned, this is set automatically to the model's global aligned window"></label>
      <label>Length<input id="length" type="number" value="2000" min="1" title="When aligned, this is the model window size"></label>
      <label>Points<input id="points" type="number" value="1000" min="1" max="20000"></label>
      <label>Alignment<div class="checks"><label><input id="align" type="checkbox"> align with bcftools_chain</label></div></label>
      <button type="submit">Refresh now</button>
      <button type="button" id="prevWindow">Previous window</button>
      <button type="button" id="nextWindow">Next window</button>
      <button type="button" id="zoomIn">Zoom in</button>
      <button type="button" id="zoomOut">Zoom out</button>
      <button type="button" id="centerWindow">Center 32768</button>
    </form>
    <div id="status">Loading options...</div>
    <div class="stats">
      <div class="stat"><b>Series</b><span id="seriesCount">-</span></div>
      <div class="stat"><b>Axis</b><span id="axis">-</span></div>
      <div class="stat"><b>Track count</b><span id="trackCount">-</span></div>
      <div class="stat"><b>Shape</b><span id="shape">-</span></div>
      <div class="stat"><b>Returned</b><span id="returned">-</span></div>
    </div>
    <div class="tabs">
      <button type="button" class="tab active" id="lineTab">Line plot</button>
      <button type="button" class="tab" id="heatmapTab">RNA-seq grayscale heatmap</button>
    </div>
    <div class="viz-grid" id="lineSection">
      <div class="plot-wrap"><svg id="plot" role="img" aria-label="Signal plot"></svg></div>
      <aside class="side-panel">
        <h2 style="margin-top:0">Selected individuals</h2>
        <div id="individualPanel">No data loaded.</div>
      </aside>
    </div>
    <div class="heatmap-wrap hidden" id="heatmapSection">
      <div class="heatmap-controls">
        <label>Source<select id="heatmapMode"><option value="sample">Selected individual</option><option value="class_mean">Class mean</option></select></label>
        <label>Class field<select id="heatmapClassField"><option value="superpopulation">Superpopulation</option><option value="pigmentation">Pigmentation</option></select></label>
        <label class="wide">Classes<select id="heatmapClassValue" multiple size="4"></select></label>
        <label>First track<input id="heatmapTrackStart" type="number" value="0" min="0"></label>
        <label>Track rows<input id="heatmapTrackCount" type="number" value="128" min="1" max="512"></label>
        <label>Max samples<input id="heatmapMaxSamples" type="number" value="60" min="1" max="120"></label>
        <button type="button" id="loadHeatmap">Render heatmap</button>
      </div>
      <div class="sample-note" id="heatmapInfo">Use output rna_seq to visualize intensity across tracks. In class mean mode, selected individuals are ignored. The plot scrolls horizontally at pixel level.</div>
      <div class="heatmap-panels" id="heatmapPanels">
        <div class="heatmap-panel">
          <h3>Heatmap</h3>
          <div class="heatmap-canvas-wrap"><canvas id="heatmap" role="img" aria-label="RNA-seq grayscale heatmap"></canvas></div>
        </div>
      </div>
    </div>
    <div class="meta" id="meta"></div>
  </section>
</main>
<script>
const $ = (id) => document.getElementById(id);
const fmt = (v) => v === null || v === undefined || Number.isNaN(v) ? '-' : Number(v).toPrecision(6);
let allSamples = [];
let allIndividuals = [];
let populationRows = [];
let superpopulationRows = [];
let pigmentationRows = [];
let defaultWindow = {start: 1, length: 32768};
let autoLoadTimer = null;

function apiUrl(url) {
  if (!url.startsWith('/api')) return url;
  const parts = window.location.pathname.split('/').filter(Boolean);
  if (parts[0] === 'apps' && parts[1]) return `/apps/${parts[1]}${url}`;
  return url;
}

function fillSelect(select, values) {
  select.innerHTML = '';
  for (const value of values || []) {
    const option = document.createElement('option');
    option.value = value;
    option.textContent = value;
    select.appendChild(option);
  }
}

function selectedSelectValues(select) {
  return Array.from(select.selectedOptions || []).map(option => option.value).filter(Boolean);
}

function fillTrackSelect(tracks) {
  const previous = $('track').value;
  $('track').innerHTML = '';
  for (const track of tracks || []) {
    const option = document.createElement('option');
    option.value = track.index;
    option.textContent = track.label;
    option.title = JSON.stringify(track.metadata || {}, null, 2);
    $('track').appendChild(option);
  }
  if (previous) $('track').value = previous;
  if (!$('track').value && $('track').options.length) $('track').selectedIndex = 0;
}

function fillClassValues() {
  const field = $('heatmapClassField').value;
  const rows = field === 'pigmentation' ? pigmentationRows : superpopulationRows;
  fillSelect($('heatmapClassValue'), rows.map(row => row.id));
  Array.from($('heatmapClassValue').options).slice(0, 2).forEach(option => { option.selected = true; });
}

function updateHeatmapModeControls() {
  const classMean = $('heatmapMode').value === 'class_mean';
  $('heatmapMaxSamples').disabled = classMean;
  $('heatmapMaxSamples').title = classMean ? 'Class mean uses all individuals in the selected class' : '';
}

function updateWindowControlsState() {
  const aligned = $('align').checked;
  $('start').disabled = aligned;
  $('length').disabled = aligned;
  ['prevWindow','nextWindow','zoomIn','zoomOut','centerWindow'].forEach(id => { $(id).disabled = aligned; });
  $('start').title = aligned ? 'Automatically set to the training model global aligned window for the selected gene' : '';
  $('length').title = aligned ? 'Automatically set to the training model window size (32768 or shorter at boundaries)' : '';
}

async function syncModelWindow() {
  updateWindowControlsState();
  if (!$('gene').value) return;
  const params = new URLSearchParams({gene: $('gene').value, align: $('align').checked ? 'true' : 'false'});
  const response = await fetch(apiUrl('/api/model-window?' + params.toString()));
  const data = await response.json();
  if (!response.ok || data.error) throw new Error(data.error || response.statusText);
  $('start').value = data.start;
  $('length').value = data.length;
  defaultWindow = {start: data.start || 1, length: data.length || 32768};
  return data;
}

function checkedValues(id) { return Array.from(document.querySelectorAll(`#${id} input:checked`)).map(el => el.value); }

function setChecks(id, checked) {
  document.querySelectorAll(`#${id} input`).forEach(el => { el.checked = checked; });
  updateCounts();
  if (id === 'superpopChecks') renderPopulations();
  else if (id !== 'sampleChecks') renderSamples();
  scheduleLoad();
}

function updateCounts() {
  $('superpopCount').textContent = `${checkedValues('superpopChecks').length}/${document.querySelectorAll('#superpopChecks input').length}`;
  $('populationCount').textContent = `${checkedValues('populationChecks').length}/${document.querySelectorAll('#populationChecks input').length}`;
  $('sampleCount').textContent = `${checkedValues('sampleChecks').length}/${document.querySelectorAll('#sampleChecks input').length}`;
}

function renderCheckList(id, rows, searchId, countId, formatter) {
  const previous = new Set(checkedValues(id));
  const hadPrevious = document.querySelectorAll(`#${id} input`).length > 0;
  const needle = $(searchId).value.trim().toLowerCase();
  const filtered = rows.filter(row => !needle || String(row.id).toLowerCase().includes(needle));
  $(id).innerHTML = filtered.map(row => {
    const checked = (!hadPrevious || previous.has(row.id)) ? 'checked' : '';
    return `<label><input type="checkbox" value="${row.id}" ${checked}>${formatter(row)}</label>`;
  }).join('');
  $(countId).textContent = `${checkedValues(id).length}/${filtered.length}`;
  document.querySelectorAll(`#${id} input`).forEach(el => el.addEventListener('change', () => {
    updateCounts();
    if (id === 'superpopChecks') renderPopulations();
    else renderSamples();
    scheduleLoad();
  }));
}

function renderPopulations() {
  const previous = new Set(checkedValues('populationChecks'));
  const hadPrevious = document.querySelectorAll('#populationChecks input').length > 0;
  const superpops = new Set(checkedValues('superpopChecks'));
  const counts = new Map();
  for (const row of allIndividuals) {
    if (superpops.size && !superpops.has(row.superpopulation || '')) continue;
    const pop = row.population || '';
    if (!pop) continue;
    counts.set(pop, (counts.get(pop) || 0) + 1);
  }
  const rows = Array.from(counts.entries()).sort((a, b) => a[0].localeCompare(b[0])).map(([id, count]) => ({id, count}));
  const needle = $('populationSearch').value.trim().toLowerCase();
  const filtered = rows.filter(row => !needle || row.id.toLowerCase().includes(needle));
  $('populationChecks').innerHTML = filtered.map(row => {
    const checked = (!hadPrevious || previous.has(row.id)) ? 'checked' : '';
    return `<label><input type="checkbox" value="${row.id}" ${checked}>${row.id} (${row.count})</label>`;
  }).join('');
  $('populationCount').textContent = `${checkedValues('populationChecks').length}/${filtered.length}`;
  document.querySelectorAll('#populationChecks input').forEach(el => el.addEventListener('change', () => { updateCounts(); renderSamples(); scheduleLoad(); }));
  renderSamples();
}

function renderSamples() {
  const previous = new Set(checkedValues('sampleChecks'));
  const hadPrevious = document.querySelectorAll('#sampleChecks input').length > 0;
  const needle = $('sampleFilter').value.trim().toLowerCase();
  const superpops = new Set(checkedValues('superpopChecks'));
  const pops = new Set(checkedValues('populationChecks'));
  let rows = allIndividuals;
  if (superpops.size) rows = rows.filter(row => superpops.has(row.superpopulation));
  if (pops.size) rows = rows.filter(row => pops.has(row.population));
  if (needle) rows = rows.filter(row => row.sample_id.toLowerCase().includes(needle));
  const visible = rows.slice(0, 500);
  $('sampleChecks').innerHTML = visible.map((row, idx) => {
    const checked = (!hadPrevious ? idx < 3 : previous.has(row.sample_id)) ? 'checked' : '';
    return `<label><input type="checkbox" value="${row.sample_id}" ${checked}>${row.sample_id} | ${row.superpopulation || '-'} / ${row.population || '-'}</label>`;
  }).join('');
  $('sampleNote').textContent = `Showing ${visible.length} of ${rows.length} matching samples. Max 24 sample/haplotype series.`;
  document.querySelectorAll('#sampleChecks input').forEach(el => el.addEventListener('change', () => { updateCounts(); scheduleLoad(); }));
  updateCounts();
}

const colors = ['#67d2ff', '#ff7b72', '#7ee787', '#d2a8ff', '#ffd166', '#ffa657', '#79c0ff', '#f778ba'];
const formatAxis = (v) => Math.round(Number(v)).toLocaleString('en-US', {useGrouping: false});

function scheduleLoad(delay = 350) {
  clearTimeout(autoLoadTimer);
  autoLoadTimer = setTimeout(() => {
    if (!$('heatmapSection').classList.contains('hidden')) loadHeatmap().catch(showError);
    else loadTrack().catch(showError);
  }, delay);
}

function drawPlot(series) {
  const svg = $('plot');
  const width = svg.clientWidth || 800;
  const height = svg.clientHeight || 360;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.innerHTML = '';
  const pad = {l: 54, r: 18, t: 18, b: 34};
  const all = [];
  for (const s of series) for (let i = 0; i < s.y.length; i++) if (Number.isFinite(s.y[i])) all.push({x: s.x[i], y: s.y[i]});
  const finite = all;
  if (!finite.length) {
    svg.innerHTML = `<text x="${width / 2}" y="${height / 2}" text-anchor="middle" fill="#9ca9bf">No finite values</text>`;
    return;
  }
  const xmin = Math.min(...finite.map(p => p.x));
  const xmax = Math.max(...finite.map(p => p.x));
  let ymin = Math.min(...finite.map(p => p.y));
  let ymax = Math.max(...finite.map(p => p.y));
  if (ymin === ymax) { ymin -= 1; ymax += 1; }
  const sx = (v) => pad.l + ((v - xmin) / Math.max(1, xmax - xmin)) * (width - pad.l - pad.r);
  const sy = (v) => pad.t + (1 - ((v - ymin) / (ymax - ymin))) * (height - pad.t - pad.b);
  const grid = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  grid.setAttribute('stroke', '#25344e');
  grid.setAttribute('stroke-width', '1');
  for (let i = 0; i <= 4; i++) {
    const yy = pad.t + i * (height - pad.t - pad.b) / 4;
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('x1', pad.l); line.setAttribute('x2', width - pad.r); line.setAttribute('y1', yy); line.setAttribute('y2', yy);
    grid.appendChild(line);
  }
  svg.appendChild(grid);
  const centerX = Math.round((xmin + xmax) / 2);
  const centerLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
  centerLine.setAttribute('x1', sx(centerX)); centerLine.setAttribute('x2', sx(centerX));
  centerLine.setAttribute('y1', pad.t); centerLine.setAttribute('y2', height - pad.b);
  centerLine.setAttribute('stroke', '#ffd166'); centerLine.setAttribute('stroke-width', '1.5'); centerLine.setAttribute('stroke-dasharray', '4 4');
  svg.appendChild(centerLine);
  const centerLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
  centerLabel.setAttribute('x', sx(centerX) + 4); centerLabel.setAttribute('y', pad.t + 14); centerLabel.setAttribute('fill', '#ffd166'); centerLabel.setAttribute('font-size', '12');
  centerLabel.textContent = `center ${formatAxis(centerX)}`;
  svg.appendChild(centerLabel);
  series.forEach((s, si) => {
    const pts = s.y.map((v, i) => ({x: s.x[i], y: v})).filter(p => Number.isFinite(p.y));
    const d = pts.map((p, i) => `${i ? 'L' : 'M'}${sx(p.x).toFixed(2)},${sy(p.y).toFixed(2)}`).join(' ');
    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    path.setAttribute('d', d);
    path.setAttribute('fill', 'none');
    path.setAttribute('stroke', colors[si % colors.length]);
    path.setAttribute('stroke-width', '2');
    path.setAttribute('opacity', '0.9');
    path.dataset.seriesIndex = String(si);
    path.classList.add('series-path');
    svg.appendChild(path);
  });
  const labels = [
    [pad.l, height - 10, xmin], [width - pad.r, height - 10, xmax],
    [8, sy(ymin), ymin], [8, sy(ymax), ymax]
  ];
  for (const [lx, ly, text] of labels) {
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', lx); t.setAttribute('y', ly); t.setAttribute('fill', '#9ca9bf'); t.setAttribute('font-size', '12');
    if (lx > width / 2) t.setAttribute('text-anchor', 'end');
    t.textContent = lx === 8 ? Number(text).toPrecision(5) : formatAxis(text);
    svg.appendChild(t);
  }
  series.slice(0, 12).forEach((s, i) => {
    const y = 18 + i * 16;
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('x1', width - 190); line.setAttribute('x2', width - 170); line.setAttribute('y1', y); line.setAttribute('y2', y);
    line.setAttribute('stroke', colors[i % colors.length]); line.setAttribute('stroke-width', '3');
    line.dataset.seriesIndex = String(i);
    line.style.cursor = 'pointer';
    line.addEventListener('mouseenter', () => highlightSeries(i));
    line.addEventListener('mouseleave', clearHighlight);
    svg.appendChild(line);
    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', width - 165); label.setAttribute('y', y + 4); label.setAttribute('fill', '#dbe7f3'); label.setAttribute('font-size', '12');
    label.textContent = s.label;
    label.style.cursor = 'pointer';
    label.dataset.seriesIndex = String(i);
    label.addEventListener('mouseenter', () => highlightSeries(i));
    label.addEventListener('mouseleave', clearHighlight);
    svg.appendChild(label);
  });
}

function setTab(tab) {
  const heatmap = tab === 'heatmap';
  $('lineTab').classList.toggle('active', !heatmap);
  $('heatmapTab').classList.toggle('active', heatmap);
  $('lineSection').classList.toggle('hidden', heatmap);
  $('heatmapSection').classList.toggle('hidden', !heatmap);
  if (heatmap) loadHeatmap().catch(showError);
  else loadTrack().catch(showError);
}

function drawHeatmap(data, canvas) {
  const matrix = data.matrix || [];
  const rows = matrix.length;
  const cols = rows ? matrix[0].length : 0;
  const wrap = canvas.parentElement;
  if (!rows || !cols) {
    wrap.innerHTML = '<div style="display:grid;place-items:center;height:100%;color:#9ca9bf">No heatmap values</div>';
    return;
  }
  const left = 170;
  const top = 18;
  const bottom = 34;
  const right = 120;
  const cellW = 4;
  const cellH = Math.max(8, Math.floor((wrap.clientHeight - top - bottom) / Math.max(rows, 1)));
  const cssWidth = left + cols * cellW + right;
  const cssHeight = Math.max(220, top + rows * cellH + bottom);
  const vmin = 0;
  const vmax = Math.max(Number(data.stats && data.stats.max) || 0, 1e-12);

  wrap.innerHTML = '';
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', `0 0 ${cssWidth} ${cssHeight}`);
  svg.style.width = `${cssWidth}px`;
  svg.style.height = `${cssHeight}px`;
  svg.style.display = 'block';
  svg.style.background = '#02050a';
  const cells = document.createElementNS(ns, 'g');
  cells.setAttribute('shape-rendering', 'crispEdges');
  const fragment = document.createDocumentFragment();
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const v = Number(matrix[r][c]);
      const rect = document.createElementNS(ns, 'rect');
      rect.setAttribute('x', left + c * cellW);
      rect.setAttribute('y', top + r * cellH);
      rect.setAttribute('width', cellW);
      rect.setAttribute('height', cellH);
      if (!Number.isFinite(v)) {
        rect.setAttribute('fill', '#122033');
      } else {
        const g = Math.max(0, Math.min(255, Math.round(255 * (v / vmax))));
        rect.setAttribute('fill', `rgb(${g},${g},${g})`);
        rect.setAttribute('data-value', String(v));
      }
      rect.setAttribute('data-track-row', String(r));
      rect.setAttribute('data-position-index', String(c));
      rect.appendChild(document.createElementNS(ns, 'title')).textContent = `${(data.track_labels && data.track_labels[r]) || 'track ' + r}\nposition=${(data.x || [])[c] ?? c}\nvalue=${Number.isFinite(v) ? v : 'NA'}`;
      fragment.appendChild(rect);
    }
  }
  cells.appendChild(fragment);
  svg.appendChild(cells);

  const step = Math.max(1, Math.ceil(rows / 28));
  for (let r = 0; r < rows; r += step) {
    const label = (data.track_labels && data.track_labels[r]) || `Track ${Number(data.track_start || 0) + r}`;
    const text = document.createElementNS(ns, 'text');
    text.setAttribute('x', 8);
    text.setAttribute('y', top + r * cellH + cellH / 2 + 4);
    text.setAttribute('fill', '#9ca9bf');
    text.setAttribute('font-size', '12');
    text.textContent = label.length > 24 ? label.slice(0, 24) + '...' : label;
    svg.appendChild(text);
  }
  const x = data.x || [];
  if (x.length) {
    const startText = document.createElementNS(ns, 'text');
    startText.setAttribute('x', left);
    startText.setAttribute('y', cssHeight - 10);
    startText.setAttribute('fill', '#9ca9bf');
    startText.setAttribute('font-size', '12');
    startText.textContent = formatAxis(x[0]);
    svg.appendChild(startText);
    const endText = document.createElementNS(ns, 'text');
    endText.setAttribute('x', cssWidth - 18);
    endText.setAttribute('y', cssHeight - 10);
    endText.setAttribute('text-anchor', 'end');
    endText.setAttribute('fill', '#9ca9bf');
    endText.setAttribute('font-size', '12');
    endText.textContent = formatAxis(x[x.length - 1]);
    svg.appendChild(endText);
  }
  const barX = cssWidth - 86;
  const barY = 12;
  const barW = 70;
  const barH = 8;
  const defs = document.createElementNS(ns, 'defs');
  const grad = document.createElementNS(ns, 'linearGradient');
  grad.setAttribute('id', `gray_${Math.random().toString(36).slice(2)}`);
  grad.setAttribute('x1', '0%');
  grad.setAttribute('x2', '100%');
  const stop0 = document.createElementNS(ns, 'stop');
  stop0.setAttribute('offset', '0%');
  stop0.setAttribute('stop-color', '#000');
  const stop1 = document.createElementNS(ns, 'stop');
  stop1.setAttribute('offset', '100%');
  stop1.setAttribute('stop-color', '#fff');
  grad.appendChild(stop0);
  grad.appendChild(stop1);
  defs.appendChild(grad);
  svg.appendChild(defs);
  const bar = document.createElementNS(ns, 'rect');
  bar.setAttribute('x', barX);
  bar.setAttribute('y', barY);
  bar.setAttribute('width', barW);
  bar.setAttribute('height', barH);
  bar.setAttribute('fill', `url(#${grad.getAttribute('id')})`);
  svg.appendChild(bar);
  const scaleText = document.createElementNS(ns, 'text');
  scaleText.setAttribute('x', barX);
  scaleText.setAttribute('y', barY + 24);
  scaleText.setAttribute('fill', '#9ca9bf');
  scaleText.setAttribute('font-size', '12');
  scaleText.textContent = `${fmt(vmin)}..${fmt(vmax)}`;
  svg.appendChild(scaleText);
  wrap.appendChild(svg);
}

function highlightSeries(index) {
  document.querySelectorAll('#plot .series-path').forEach(path => {
    const active = Number(path.dataset.seriesIndex) === index;
    path.setAttribute('opacity', active ? '1' : '0.12');
    path.setAttribute('stroke-width', active ? '4' : '1.5');
  });
}

function clearHighlight() {
  document.querySelectorAll('#plot .series-path').forEach(path => {
    path.setAttribute('opacity', '0.9');
    path.setAttribute('stroke-width', '2');
  });
}

function renderIndividualsPanel(rows) {
  if (!rows || !rows.length) {
    $('individualPanel').textContent = 'No selected individual metadata.';
    return;
  }
  $('individualPanel').innerHTML = rows.map(row => `
    <div class="person">
      <b>${row.sample_id}</b>
      <div class="kv"><span>Superpop</span><span>${row.superpopulation || '-'}</span></div>
      <div class="kv"><span>Population</span><span>${row.population || '-'}</span></div>
      <div class="kv"><span>Sex</span><span>${row.sex || '-'}</span></div>
      <div class="kv"><span>Family</span><span>${row.family_id || '-'}</span></div>
      <div class="kv"><span>Longevity</span><span>${row.longevity ?? '-'}</span></div>
      <div class="kv"><span>Windows</span><span>${row.window_count ?? '-'}</span></div>
      <div class="kv"><span>Frog max</span><span>${row.frog_likelihood_max == null ? '-' : Number(row.frog_likelihood_max).toPrecision(4)}</span></div>
    </div>
  `).join('');
}

async function loadTrackOptions() {
  const sample = selectedSamples()[0] || (allIndividuals[0] && allIndividuals[0].sample_id) || '';
  if (!sample || !$('gene').value || !$('output').value) return;
  const haplotype = $('hapH1').checked ? 'H1' : 'H2';
  const params = new URLSearchParams({sample, gene: $('gene').value, haplotype, output: $('output').value});
  const response = await fetch(apiUrl('/api/track-options?' + params.toString()));
  const data = await response.json();
  if (!response.ok || data.error) throw new Error(data.error || response.statusText);
  fillTrackSelect(data.tracks || []);
  $('shape').textContent = `[${data.array_shape.join(', ')}] ${data.array_dtype}`;
}

function selectedSamples() { return checkedValues('sampleChecks'); }
function selectedHaplotypes() {
  const h = [];
  if ($('hapH1').checked) h.push('H1');
  if ($('hapH2').checked) h.push('H2');
  return h;
}

function moveWindow(direction) {
  const start = Math.max(1, Number($('start').value || 1));
  const length = Math.max(1, Number($('length').value || 1));
  $('start').value = Math.max(1, start + direction * length);
  if (!$('heatmapSection').classList.contains('hidden')) loadHeatmap().catch(showError);
  else loadTrack().catch(showError);
}

function zoom(factor) {
  const start = Math.max(1, Number($('start').value || 1));
  const length = Math.max(1, Number($('length').value || 1));
  const center = start + Math.floor(length / 2);
  const nextLength = Math.max(10, Math.round(length * factor));
  $('length').value = nextLength;
  $('start').value = Math.max(1, center - Math.floor(nextLength / 2));
  if (!$('heatmapSection').classList.contains('hidden')) loadHeatmap().catch(showError);
  else loadTrack().catch(showError);
}

function centerDefaultWindow() {
  $('start').value = defaultWindow.start;
  $('length').value = defaultWindow.length;
  if (!$('heatmapSection').classList.contains('hidden')) loadHeatmap().catch(showError);
  else loadTrack().catch(showError);
}

async function loadTrack() {
  $('status').className = '';
  $('status').textContent = 'Loading track...';
  const samples = selectedSamples();
  const haplotypes = selectedHaplotypes();
  if (!samples.length) throw new Error('Select at least one sample');
  if (!haplotypes.length) throw new Error('Select at least one haplotype');
  const params = new URLSearchParams({
    samples: samples.join(','),
    gene: $('gene').value,
    haplotypes: haplotypes.join(','),
    output: $('output').value,
    track: $('track').value,
    start: $('start').value,
    length: $('length').value,
    points: $('points').value,
    align: $('align').checked ? 'true' : 'false',
  });
  const response = await fetch(apiUrl('/api/tracks?' + params.toString()));
  const data = await response.json();
  if (!response.ok || data.error) throw new Error(data.error || response.statusText);
  $('seriesCount').textContent = data.series.length;
  $('axis').textContent = data.x_axis;
  $('trackCount').textContent = data.track_count;
  $('shape').textContent = `[${data.array_shape.join(', ')}] ${data.array_dtype}`;
  $('returned').textContent = data.series.map(s => `${s.label}:${s.returned_points}`).join(', ');
  $('meta').textContent = `axis=${data.x_axis}\nexpanded_length=${data.expanded_length || '-'}\ntrack_metadata=${JSON.stringify(data.track_metadata || {}, null, 2)}\nseries=${JSON.stringify(data.series.map(s => ({label:s.label, stats:s.stats, source:s.source_file})), null, 2)}`;
  renderIndividualsPanel(data.individuals || []);
  drawPlot(data.series);
  $('status').textContent = 'Loaded.';
}

function heatmapParamsFor(classValue = '') {
  const samples = selectedSamples();
  const haplotypes = selectedHaplotypes();
  return new URLSearchParams({
    mode: $('heatmapMode').value,
    sample: samples[0] || '',
    class_field: $('heatmapClassField').value,
    class_value: classValue,
    gene: $('gene').value,
    haplotypes: haplotypes.join(','),
    output: $('output').value,
    start: $('start').value,
    length: $('length').value,
    points: $('points').value,
    track_start: $('heatmapTrackStart').value,
    track_count: $('heatmapTrackCount').value,
    max_samples: $('heatmapMaxSamples').value,
    align: $('align').checked ? 'true' : 'false',
  });
}

function renderHeatmapPanels(count, titles) {
  $('heatmapPanels').innerHTML = '';
  const canvases = [];
  for (let i = 0; i < count; i++) {
    const panel = document.createElement('div');
    panel.className = 'heatmap-panel';
    const title = document.createElement('h3');
    title.textContent = titles[i] || `Heatmap ${i + 1}`;
    const wrap = document.createElement('div');
    wrap.className = 'heatmap-canvas-wrap';
    const canvas = document.createElement('div');
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', title.textContent);
    wrap.appendChild(canvas);
    panel.appendChild(title);
    panel.appendChild(wrap);
    $('heatmapPanels').appendChild(panel);
    canvases.push(canvas);
  }
  return canvases;
}

function interleaveHeatmaps(payloads) {
  if (!payloads.length) return null;
  const minRows = Math.min(...payloads.map(item => (item.matrix || []).length));
  const minCols = Math.min(...payloads.map(item => (item.x || []).length));
  const matrix = [];
  const trackLabels = [];
  const values = [];
  for (let row = 0; row < minRows; row++) {
    for (const item of payloads) {
      const sourceRow = (item.matrix && item.matrix[row]) || [];
      const clipped = sourceRow.slice(0, minCols);
      matrix.push(clipped);
      const baseLabel = (item.track_labels && item.track_labels[row]) || `Track ${Number(item.track_start || 0) + row}`;
      trackLabels.push(`${item.class_value || item.sample}: ${baseLabel}`);
      for (const v of clipped) if (Number.isFinite(Number(v))) values.push(Number(v));
    }
  }
  const sorted = values.slice().sort((a, b) => a - b);
  const p99Index = sorted.length ? Math.min(sorted.length - 1, Math.floor(sorted.length * 0.99)) : 0;
  const vmax = sorted.length ? Math.max(sorted[p99Index], 1e-12) : 1;
  const mean = values.length ? values.reduce((acc, v) => acc + v, 0) / values.length : null;
  return {
    ...payloads[0],
    mode: 'class_mean_interleaved',
    class_value: payloads.map(item => item.class_value).join(' + '),
    samples_used: payloads.flatMap(item => item.samples_used || []),
    samples_requested: payloads.reduce((acc, item) => acc + Number(item.samples_requested || 0), 0),
    samples_truncated: payloads.some(item => item.samples_truncated),
    track_rows: matrix.length,
    track_labels: trackLabels,
    x: (payloads[0].x || []).slice(0, minCols),
    matrix,
    stats: {
      min: sorted.length ? sorted[0] : null,
      max: sorted.length ? sorted[sorted.length - 1] : null,
      mean,
      vmax,
    },
    missing_count: payloads.reduce((acc, item) => acc + Number(item.missing_count || 0), 0),
    cache: null,
  };
}

async function fetchHeatmap(classValue = '') {
  const response = await fetch(apiUrl('/api/heatmap?' + heatmapParamsFor(classValue).toString()));
  const data = await response.json();
  if (!response.ok || data.error) throw new Error(data.error || response.statusText);
  return data;
}

async function loadHeatmap() {
  $('status').className = '';
  $('status').textContent = 'Loading RNA-seq heatmap...';
  const samples = selectedSamples();
  const haplotypes = selectedHaplotypes();
  if (!haplotypes.length) throw new Error('Select at least one haplotype');
  const mode = $('heatmapMode').value;
  updateHeatmapModeControls();
  if (mode === 'sample' && !samples.length) throw new Error('Select one sample for sample heatmap');
  const classValues = mode === 'class_mean' ? selectedSelectValues($('heatmapClassValue')) : [''];
  if (mode === 'class_mean' && !classValues.length) throw new Error('Select at least one class for class mean heatmap');
  const titles = mode === 'class_mean' ? classValues : [`${samples[0]}, ${haplotypes.join('+')}`];
  const canvases = renderHeatmapPanels(mode === 'class_mean' ? 1 : titles.length, mode === 'class_mean' ? [`Interleaved: ${titles.join(' | ')}`] : titles);
  const payloads = [];
  for (let idx = 0; idx < classValues.length; idx++) {
    $('status').textContent = `Loading RNA-seq heatmap ${idx + 1}/${classValues.length}...`;
    const data = await fetchHeatmap(classValues[idx]);
    payloads.push(data);
  }
  const data = mode === 'class_mean' ? interleaveHeatmaps(payloads) : payloads[0];
  if (mode === 'class_mean') drawHeatmap(data, canvases[0]);
  else payloads.forEach((item, idx) => drawHeatmap(item, canvases[idx]));
  const source = mode === 'class_mean'
    ? payloads.map(item => `${item.class_value}: ${item.samples_used.length}/${item.samples_requested}${item.cache && item.cache.path ? ` ${item.cache.hit ? 'hit' : 'miss'}` : ''}`).join(' | ')
    : `${data.sample}, haplotypes ${data.haplotypes.join('+')}`;
  const cacheInfo = data.cache && data.cache.path ? ` | cache=${data.cache.hit ? 'hit' : 'miss'}` : '';
  const win = data.effective_window ? `window=${data.effective_window.start}+${data.effective_window.length - 1}` : 'window=-';
  $('heatmapInfo').textContent = `Heatmap ${data.track_rows} tracks x ${data.x.length} positions | ${source} | axis=${data.x_axis} | ${win} | grayscale=absolute 0..max | SVG rect plot, horizontal scroll, 4px/position | stats min=${fmt(data.stats.min)} max=${fmt(data.stats.max)} mean=${fmt(data.stats.mean)} | missing=${payloads.map(item => item.missing_count).join(',')}${mode === 'class_mean' ? '' : cacheInfo}`;
  $('returned').textContent = `${data.track_rows}x${data.x.length}`;
  $('axis').textContent = data.x_axis;
  $('trackCount').textContent = data.track_count_total;
  $('shape').textContent = `[${data.array_shape.join(', ')}] ${data.array_dtype}`;
  $('status').textContent = 'Heatmap loaded.';
}

async function init() {
  const response = await fetch(apiUrl('/api/options'));
  const options = await response.json();
  if (!response.ok || options.error) throw new Error(options.error || response.statusText);
  allSamples = options.samples || [];
  renderSamples('');
  fillSelect($('gene'), options.genes);
  fillSelect($('output'), options.outputs && options.outputs.length ? options.outputs : ['rna_seq']);
  const d = options.defaults || {};
  if (d.gene) $('gene').value = d.gene;
  if (d.output) $('output').value = d.output;
  if (d.start) $('start').value = d.start;
  if (d.length) $('length').value = d.length;
  defaultWindow = {start: d.start || 1, length: d.length || 32768};
  $('status').textContent = `Loaded options from ${options.dataset_dir}`;
  allIndividuals = options.individuals || (options.samples || []).map(sample => ({sample_id: sample, population: '', superpopulation: ''}));
  populationRows = Object.entries(options.populations || {}).map(([id, count]) => ({id, count}));
  superpopulationRows = Object.entries(options.superpopulations || {}).map(([id, count]) => ({id, count}));
  pigmentationRows = Object.entries(options.pigmentations || {}).map(([id, count]) => ({id, count}));
  renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`);
  renderPopulations();
  fillClassValues();
  updateHeatmapModeControls();
  updateWindowControlsState();
  await syncModelWindow();
  await loadTrackOptions();
  if (selectedSamples().length && $('gene').value) await loadTrack();
}

$('controls').addEventListener('submit', (event) => { event.preventDefault(); loadTrack().catch(showError); });
$('lineTab').addEventListener('click', () => setTab('line'));
$('heatmapTab').addEventListener('click', () => setTab('heatmap'));
$('heatmapClassField').addEventListener('change', () => { fillClassValues(); loadHeatmap().catch(showError); });
$('heatmapMode').addEventListener('change', () => { updateHeatmapModeControls(); loadHeatmap().catch(showError); });
$('heatmapClassValue').addEventListener('change', () => loadHeatmap().catch(showError));
$('loadHeatmap').addEventListener('click', () => loadHeatmap().catch(showError));
['heatmapTrackStart','heatmapTrackCount','heatmapMaxSamples'].forEach(id => $(id).addEventListener('change', () => loadHeatmap().catch(showError)));
$('superpopSearch').addEventListener('input', () => { renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`); renderPopulations(); });
$('populationSearch').addEventListener('input', renderPopulations);
$('sampleFilter').addEventListener('input', renderSamples);
$('superpopAll').addEventListener('click', () => setChecks('superpopChecks', true));
$('superpopNone').addEventListener('click', () => setChecks('superpopChecks', false));
$('populationAll').addEventListener('click', () => setChecks('populationChecks', true));
$('populationNone').addEventListener('click', () => setChecks('populationChecks', false));
$('sampleFirst5').addEventListener('click', () => { document.querySelectorAll('#sampleChecks input').forEach((el, i) => { el.checked = i < 5; }); updateCounts(); scheduleLoad(); });
$('sampleVisible').addEventListener('click', () => setChecks('sampleChecks', true));
$('sampleNone').addEventListener('click', () => setChecks('sampleChecks', false));
$('prevWindow').addEventListener('click', () => moveWindow(-1));
$('nextWindow').addEventListener('click', () => moveWindow(1));
$('zoomIn').addEventListener('click', () => zoom(0.5));
$('zoomOut').addEventListener('click', () => zoom(2));
$('centerWindow').addEventListener('click', centerDefaultWindow);
$('gene').addEventListener('change', () => syncModelWindow().then(() => loadTrackOptions()).then(() => scheduleLoad()).catch(showError));
$('align').addEventListener('change', () => syncModelWindow().then(() => scheduleLoad()).catch(showError));
['output','hapH1','hapH2'].forEach(id => {
  $(id).addEventListener('change', () => loadTrackOptions().then(() => scheduleLoad()).catch(showError));
});
['track','start','length','points'].forEach(id => {
  $(id).addEventListener('change', () => scheduleLoad());
});
function showError(err) { $('status').className = 'error'; $('status').textContent = err.message || String(err); }
init().catch(showError);
</script>
</body>
</html>
"""


def serve(dataset_dir: Path, host: str, port: int, consensus_dataset_dir: Optional[Path] = None) -> None:
    repository = AlphaGenomeRepository(dataset_dir, consensus_dataset_dir=consensus_dataset_dir)
    AlphaGenomeTrackViewerHandler.repository = repository
    server = ThreadingHTTPServer((host, port), AlphaGenomeTrackViewerHandler)
    print(f"AlphaGenome Track Viewer: http://{host}:{port}")
    print(f"Dataset: {repository.dataset_dir}")
    print(f"Consensus dataset: {repository.consensus_dataset_dir}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Local AlphaGenome track viewer")
    parser.add_argument("dataset_dir", type=Path, help="Dataset directory containing dataset_metadata.json")
    parser.add_argument("--host", default="127.0.0.1", help="Bind host")
    parser.add_argument("--port", type=int, default=8774, help="Bind port")
    parser.add_argument("--consensus-dataset-dir", type=Path, default=DEFAULT_CONSENSUS_DATASET_DIR)
    args = parser.parse_args(argv)
    serve(args.dataset_dir, args.host, args.port, args.consensus_dataset_dir)


if __name__ == "__main__":
    main()
