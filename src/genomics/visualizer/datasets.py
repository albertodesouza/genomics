"""Dataset catalog for the visualizer.

A dataset is any directory following the canonical layout::

    dataset_metadata.json
    references/windows/<gene>/ref.window.fa (+ window_metadata.json)
    individuals/<sample>/windows/<gene>/<sample>.<H>.window.fixed.fa
    individuals/<sample>/windows/<gene>/<sample>.window[.consensus_ready].vcf.gz
    individuals/<sample>/windows/<gene>/predictions_<H>/<output>.npz (+ <output>_metadata.json)

Track outputs may be binned (``resolution`` bp per row, e.g. 128 for ChIP-seq); contact maps and
splice junctions are listed separately (``other_outputs``) since they are not 1-D tracks.

Nothing here is specific to 1000 Genomes: sample facets are taken from whatever fields the
pedigree records carry (plus an optional annotation table), genes are discovered from the
reference windows on disk, and outputs/haplotypes/tracks from the prediction files.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

# Population-derived labels used by the pigmentation experiments. Applied only to datasets whose
# samples carry a ``population`` field with at least one of these codes.
DERIVED_FIELDS: Dict[str, Tuple[str, Dict[str, str]]] = {
    "pigmentation": (
        "population",
        {
            "YRI": "strong pigmentation",
            "ESN": "strong pigmentation",
            "LWK": "strong pigmentation",
            "MSL": "strong pigmentation",
            "GWD": "strong pigmentation",
            "FIN": "weak pigmentation",
            "CEU": "weak pigmentation",
            "GBR": "weak pigmentation",
        },
    ),
}
MAX_FACET_VALUES = 200
HIDDEN_FIELDS = {"sex"}  # numeric duplicate of sex_label
FASTA_HEADER_RE = re.compile(r"^>?(?P<chrom>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)")


def load_json(path: Path) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _scalar(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def slugify(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "-", text).strip("-").lower()
    return slug or "dataset"


@dataclass
class GeneWindow:
    gene: str
    chromosome: Optional[str]
    start: Optional[int]  # genomic 1-based position of reference offset 0
    end: Optional[int]
    length: Optional[int]
    strand: Optional[str]
    metadata: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "gene": self.gene,
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "strand": self.strand,
        }


def read_fasta_header(path: Path) -> Optional[str]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            line = f.readline().strip()
    except OSError:
        return None
    return line if line.startswith(">") else None


def read_fasta_bytes(path: Path) -> bytes:
    with open(path, "rb") as f:
        raw = f.read()
    lines = raw.split(b"\n")
    if lines and lines[0].startswith(b">"):
        lines = lines[1:]
    return b"".join(line.strip() for line in lines).upper()


def load_annotation_table(path: Path) -> Dict[str, Dict[str, Any]]:
    """Sample annotation CSV/TSV keyed by the first column (or ``sample_id``/``sample``)."""
    text = Path(path).read_text(encoding="utf-8")
    dialect = "excel-tab" if "\t" in text.splitlines()[0] else "excel"
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    if not reader.fieldnames:
        return {}
    key = next((name for name in ("sample_id", "sample", "Sample", "IID", "id") if name in reader.fieldnames), reader.fieldnames[0])
    rows: Dict[str, Dict[str, Any]] = {}
    for row in reader:
        sample = str(row.get(key) or "").strip()
        if sample:
            rows[sample] = {k: v for k, v in row.items() if k != key and k}
    return rows


class Dataset:
    def __init__(
        self,
        path: Path,
        dataset_id: Optional[str] = None,
        name: Optional[str] = None,
        annotations: Optional[Path] = None,
    ):
        self.path = Path(path).expanduser().resolve()
        self.metadata_path = self.path / "dataset_metadata.json"
        if not self.metadata_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json not found in {self.path}")
        self.metadata: Dict[str, Any] = load_json(self.metadata_path)
        if not isinstance(self.metadata, dict):
            raise ValueError(f"dataset_metadata.json must contain an object: {self.metadata_path}")
        self.name = name or str(self.metadata.get("dataset_name") or self.path.name)
        self.id = dataset_id or slugify(self.path.name)
        self.annotations_path = Path(annotations).resolve() if annotations else None
        self._lock = threading.Lock()
        self._gene_cache: Dict[str, Dict[str, Any]] = {}
        self._windows: Dict[str, GeneWindow] = {}
        self.field_kinds: Dict[str, str] = {}  # forced kinds of fields added with set_field
        self.samples, self.fields = self._build_samples()
        self.sample_index = {row["sample_id"]: i for i, row in enumerate(self.samples)}
        self.genes = self._discover_genes()

    # -- identity -------------------------------------------------------------
    @property
    def fingerprint(self) -> str:
        stat = self.metadata_path.stat()
        payload = f"{self.path}|{stat.st_mtime_ns}|{stat.st_size}|{self.annotations_path or ''}"
        return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:16]

    # -- samples ---------------------------------------------------------------
    def _raw_sample_ids(self) -> List[str]:
        raw = self.metadata.get("individuals", []) or []
        ids: List[str] = []
        for item in raw:
            if isinstance(item, dict):
                value = item.get("sample_id") or item.get("id") or ""
            else:
                value = item
            if str(value):
                ids.append(str(value))
        if not ids:
            individuals_dir = self.path / "individuals"
            if individuals_dir.is_dir():
                ids = sorted(p.name for p in individuals_dir.iterdir() if p.is_dir())
        return ids

    def _build_samples(self) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        pedigree = self.metadata.get("individuals_pedigree", {}) or {}
        if not isinstance(pedigree, dict):
            pedigree = {}
        annotations = load_annotation_table(self.annotations_path) if self.annotations_path else {}
        rows: List[Dict[str, Any]] = []
        for sample_id in self._raw_sample_ids():
            row: Dict[str, Any] = {"sample_id": sample_id}
            meta = pedigree.get(sample_id, {})
            if isinstance(meta, dict):
                for key, value in meta.items():
                    if key not in HIDDEN_FIELDS or "sex_label" not in meta:
                        row[str(key)] = _scalar(value)
            for key, value in annotations.get(sample_id, {}).items():
                row[str(key)] = _scalar(value)
            rows.append(row)

        for derived, (source, mapping) in DERIVED_FIELDS.items():
            if any(derived in row for row in rows):
                continue
            if not any(str(row.get(source, "")) in mapping for row in rows):
                continue
            for row in rows:
                label = mapping.get(str(row.get(source, "")))
                if label:
                    row[derived] = label

        return rows, self._describe_fields(rows)

    @staticmethod
    def _describe_fields(rows: List[Dict[str, Any]], kinds: Optional[Dict[str, str]] = None) -> List[Dict[str, Any]]:
        names: List[str] = []
        for row in rows:
            for key in row:
                if key not in names:
                    names.append(key)
        preferred = ["sample_id", "superpopulation", "population", "sex_label", "pigmentation"]
        names.sort(key=lambda n: (preferred.index(n) if n in preferred else len(preferred), n))
        fields = []
        for name in names:
            values = [row.get(name) for row in rows if row.get(name) not in (None, "")]
            counts: Dict[str, int] = {}
            numeric = bool(values) and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values)
            for value in values:
                counts[str(value)] = counts.get(str(value), 0) + 1
            categorical = name != "sample_id" and 1 < len(counts) <= MAX_FACET_VALUES and not (numeric and len(counts) > 12)
            if name == "family_id" and len(counts) > 0.5 * max(len(rows), 1):
                categorical = False
            forced = (kinds or {}).get(name)
            if forced:
                categorical = forced == "categorical"
                numeric = numeric and forced == "numeric"
            fields.append(
                {
                    "name": name,
                    "label": name.replace("_label", "").replace("_", " ").strip().capitalize(),
                    "kind": "categorical" if categorical else ("numeric" if numeric else "text"),
                    "distinct": len(counts),
                    "missing": len(rows) - len(values),
                    "counts": sorted(({"value": k, "count": v} for k, v in counts.items()), key=lambda x: (-x["count"], x["value"]))
                    if categorical
                    else [],
                }
            )
        return fields

    def set_field(self, name: str, values: Dict[str, Any], kind: Optional[str] = None) -> None:
        """Add or replace a sample field (e.g. a model's derived target) and refresh the facets.
        ``kind`` ("categorical", "numeric" or "text") overrides the kind guessed from the values."""
        with self._lock:
            for row in self.samples:
                value = values.get(row["sample_id"])
                if value is None:
                    row.pop(name, None)
                else:
                    row[name] = _scalar(value)
            if kind and values:
                self.field_kinds[name] = kind
            else:
                self.field_kinds.pop(name, None)
            self.fields = self._describe_fields(self.samples, self.field_kinds)

    def reset_gene_cache(self) -> None:
        """Forget discovered outputs/tracks (after new predictions were written)."""
        with self._lock:
            self._gene_cache.clear()

    def facet_fields(self) -> List[str]:
        return [f["name"] for f in self.fields if f["kind"] == "categorical"]

    def filter_samples(self, filters: Optional[Dict[str, List[str]]] = None, ids: Optional[Iterable[str]] = None) -> List[str]:
        selected = list(ids) if ids is not None else [row["sample_id"] for row in self.samples]
        if not filters:
            return [s for s in selected if s in self.sample_index]
        wanted = {key: {str(v) for v in values} for key, values in filters.items() if values}
        out = []
        for sample_id in selected:
            idx = self.sample_index.get(sample_id)
            if idx is None:
                continue
            row = self.samples[idx]
            if all(str(row.get(key, "")) in values for key, values in wanted.items()):
                out.append(sample_id)
        return out

    def group_samples(self, field_name: str, values: List[str], sample_ids: List[str]) -> Dict[str, List[str]]:
        groups: Dict[str, List[str]] = {value: [] for value in values}
        for sample_id in sample_ids:
            row = self.samples[self.sample_index[sample_id]]
            key = str(row.get(field_name, ""))
            if key in groups:
                groups[key].append(sample_id)
        return groups

    def sample_detail(self, sample_id: str) -> Dict[str, Any]:
        idx = self.sample_index.get(sample_id)
        if idx is None:
            raise KeyError(f"Unknown sample: {sample_id}")
        detail: Dict[str, Any] = {"sample": self.samples[idx]}
        meta_path = self.path / "individuals" / sample_id / "individual_metadata.json"
        if meta_path.exists():
            try:
                detail["individual_metadata"] = load_json(meta_path)
            except Exception as exc:  # keep the detail view useful with one malformed file
                detail["individual_metadata_error"] = str(exc)
        windows_dir = self.path / "individuals" / sample_id / "windows"
        windows = []
        if windows_dir.is_dir():
            for gene_dir in sorted(p for p in windows_dir.iterdir() if p.is_dir()):
                haps = sorted(p.name[len("predictions_"):] for p in gene_dir.glob("predictions_*") if p.is_dir())
                windows.append(
                    {
                        "gene": gene_dir.name,
                        "haplotypes": haps,
                        "outputs": sorted({npz.stem for hap in haps for npz in (gene_dir / f"predictions_{hap}").glob("*.npz")}),
                        "has_vcf": self.sample_vcf_path(sample_id, gene_dir.name) is not None,
                    }
                )
        detail["windows"] = windows
        return detail

    # -- genes -------------------------------------------------------------------
    def _discover_genes(self) -> List[str]:
        names = []
        for raw in self.metadata.get("genes", []) or []:
            value = raw.get("gene") or raw.get("name") if isinstance(raw, dict) else raw
            if value:
                names.append(str(value))
        refs = self.path / "references" / "windows"
        if refs.is_dir():
            names.extend(p.name for p in refs.iterdir() if p.is_dir())
        return sorted(set(names))

    def window(self, gene: str) -> GeneWindow:
        cached = self._windows.get(gene)
        if cached is not None:
            return cached
        if gene not in self.genes:
            raise KeyError(f"Unknown gene: {gene}")
        meta: Dict[str, Any] = {}
        meta_path = self.reference_dir(gene) / "window_metadata.json"
        if meta_path.exists():
            try:
                meta = load_json(meta_path)
            except Exception:
                meta = {}
        catalog = (self.metadata.get("window_catalog") or {}).get(gene) or {}
        chrom = meta.get("chromosome") or catalog.get("chromosome")
        start = meta.get("start") or catalog.get("start")
        end = meta.get("end") or catalog.get("end")
        # The FASTA header is authoritative for where reference offset 0 sits.
        header = read_fasta_header(self.reference_dir(gene) / "ref.window.fa")
        match = FASTA_HEADER_RE.match(header or "")
        if match:
            chrom = match.group("chrom")
            start = int(match.group("start"))
            end = int(match.group("end"))
        start = int(start) if start is not None else None
        end = int(end) if end is not None else None
        strand = (self.metadata.get("gene_strands") or {}).get(gene)
        window = GeneWindow(
            gene=gene,
            chromosome=str(chrom) if chrom else None,
            start=start,
            end=end,
            length=(end - start + 1) if start is not None and end is not None else None,
            strand=strand,
            metadata=meta or catalog,
        )
        self._windows[gene] = window
        return window

    def reference_dir(self, gene: str) -> Path:
        return self.path / "references" / "windows" / gene

    def gene_dir(self, sample_id: str, gene: str) -> Path:
        return self.path / "individuals" / sample_id / "windows" / gene

    def prediction_path(self, sample_id: str, gene: str, haplotype: str, output: str) -> Path:
        return self.gene_dir(sample_id, gene) / f"predictions_{haplotype}" / f"{output}.npz"

    def reference_prediction_path(self, gene: str, output: str) -> Path:
        """AlphaGenome prediction of the reference window (``predict-dataset --haplotypes ref``)."""
        return self.reference_dir(gene) / "predictions_ref" / f"{output}.npz"

    def reference_outputs(self, gene: str) -> Dict[str, List[Dict[str, Any]]]:
        """{output: track metadata records} of the stored reference predictions of a window."""
        folder = self.reference_dir(gene) / "predictions_ref"
        out: Dict[str, List[Dict[str, Any]]] = {}
        if not folder.is_dir():
            return out
        for npz in sorted(folder.glob("*.npz")):
            meta_path = npz.with_name(f"{npz.stem}_metadata.json")
            records: Any = []
            if meta_path.exists():
                try:
                    payload = load_json(meta_path)
                    records = payload.get("metadata", payload) if isinstance(payload, dict) else payload
                except Exception:
                    records = []
            out[npz.stem] = [r for r in records if isinstance(r, dict)] if isinstance(records, list) else []
        return out

    def haplotype_fasta(self, sample_id: str, gene: str, haplotype: str) -> Path:
        return self.gene_dir(sample_id, gene) / f"{sample_id}.{haplotype}.window.fixed.fa"

    def sample_vcf_path(self, sample_id: str, gene: str) -> Optional[Path]:
        gene_dir = self.gene_dir(sample_id, gene)
        for name in (f"{sample_id}.window.consensus_ready.vcf.gz", f"{sample_id}.window.vcf.gz", f"{sample_id}.window.vcf"):
            candidate = gene_dir / name
            if candidate.exists():
                return candidate
        return None

    def representative_sample(self, gene: str) -> Optional[str]:
        for row in self.samples:
            if (self.gene_dir(row["sample_id"], gene)).is_dir():
                return row["sample_id"]
        return None

    def gene_info(self, gene: str) -> Dict[str, Any]:
        """Outputs, haplotypes and per-output track metadata for a gene (cached)."""
        with self._lock:
            cached = self._gene_cache.get(gene)
        if cached is not None:
            return cached
        window = self.window(gene)
        sample = self.representative_sample(gene)
        haplotypes: List[str] = []
        outputs: Dict[str, Dict[str, Any]] = {}
        others: Dict[str, Dict[str, Any]] = {}
        if sample:
            gene_dir = self.gene_dir(sample, gene)
            haplotypes = sorted(p.name[len("predictions_"):] for p in gene_dir.glob("predictions_*") if p.is_dir())
            for hap in haplotypes[:1]:
                for npz in sorted((gene_dir / f"predictions_{hap}").glob("*.npz")):
                    described = self._describe_output(npz, window.length)
                    if described.get("kind") == "tracks" and described["tracks"]:
                        outputs[npz.stem] = described
                    else:
                        others[npz.stem] = described
        info = {
            **window.as_dict(),
            "representative_sample": sample,
            "haplotypes": haplotypes,
            "outputs": outputs,
            "other_outputs": others,
            "has_reference": (self.reference_dir(gene) / "ref.window.fa").exists(),
            "reference_outputs": sorted(self.reference_outputs(gene)),
            "model_window": self.model_window(gene),
        }
        with self._lock:
            self._gene_cache[gene] = info
        return info

    def _describe_output(self, npz_path: Path, window_length: Optional[int] = None) -> Dict[str, Any]:
        """Tracks, kind (tracks / contact_map / junctions), resolution and length in bp of one output."""
        from genomics.visualizer.signals import npz_header

        meta_path = npz_path.with_name(f"{npz_path.stem}_metadata.json")
        payload: Any = {}
        if meta_path.exists():
            try:
                payload = load_json(meta_path)
            except Exception:
                payload = {}
        records = payload.get("metadata", payload) if isinstance(payload, dict) else payload
        records = records if isinstance(records, list) else []
        extra = payload if isinstance(payload, dict) else {}
        try:
            header = npz_header(npz_path)
        except Exception as exc:
            return {"error": str(exc), "tracks": [], "length": None, "kind": "tracks"}
        shape = list(header["shape"])
        kind = extra.get("kind") or ("junctions" if "starts" in header["members"] else "contact_map" if len(shape) == 3 else "tracks")
        resolution = extra.get("resolution") or header.get("resolution")
        if kind == "tracks" and not resolution:
            # Older files carry no resolution: binned when the window length is a multiple of the rows.
            rows = shape[0] if shape else 0
            resolution = window_length // rows if window_length and rows and rows < window_length and window_length % rows == 0 else 1
        n_tracks = shape[-1] if len(shape) >= 2 else (1 if shape else 0)
        tracks: List[Dict[str, Any]] = []
        for idx in range(n_tracks):
            meta = records[idx] if idx < len(records) and isinstance(records[idx], dict) else {}
            tracks.append({"index": idx, "label": track_label(idx, meta), "short": track_short_label(idx, meta), "metadata": {k: _scalar(v) for k, v in meta.items()}})
        out: Dict[str, Any] = {"tracks": tracks, "kind": kind, "resolution": int(resolution) if resolution else None, "dtype": "float32", "shape": shape}
        if kind == "junctions":
            out.update(length=window_length, junctions=shape[0] if shape else 0)
        else:
            out.update(bins=shape[0] if shape else 0, length=(shape[0] if shape else 0) * int(resolution or 1))
        return out

    def model_window(self, gene: str, size: int = 32768) -> Optional[Dict[str, int]]:
        """Reference-centred training window, matching DynamicIndelAligner's centring."""
        window = self.window(gene)
        meta = window.metadata or {}
        try:
            meta_length = int(meta["end"]) - int(meta["start"]) + 1
        except (KeyError, TypeError, ValueError):
            meta_length = window.length
        if not meta_length:
            return None
        size = min(size, meta_length)
        start = max(0, meta_length // 2 - size // 2)
        end = min(meta_length, start + size)
        start = max(0, end - size)
        return {"start": start, "end": end, "size": end - start}

    # -- summary -------------------------------------------------------------------
    def summary(self) -> Dict[str, Any]:
        genes = []
        for gene in self.genes:
            try:
                window = self.window(gene).as_dict()
            except Exception:
                window = {"gene": gene}
            window["in_metadata"] = gene in {str(g) for g in self.metadata.get("genes", []) or []}
            genes.append(window)
        return {
            "id": self.id,
            "name": self.name,
            "path": str(self.path),
            "created_at": self.metadata.get("creation_date"),
            "last_updated": self.metadata.get("last_updated"),
            "layout_version": self.metadata.get("layout_version"),
            "window_size": self.metadata.get("window_size"),
            "outputs": [str(v).lower() for v in self.metadata.get("alphagenome_outputs", []) or []],
            "ontologies": self.metadata.get("ontologies", []) or [],
            "ontology_details": self.metadata.get("ontology_details", {}) or {},
            "sample_count": len(self.samples),
            "gene_count": len(self.genes),
            "genes": genes,
            "fields": self.fields,
            "annotations_path": str(self.annotations_path) if self.annotations_path else None,
        }

    def listing(self) -> Dict[str, Any]:
        return {"id": self.id, "name": self.name, "path": str(self.path), "sample_count": len(self.samples), "gene_count": len(self.genes)}


def track_label(idx: int, meta: Optional[Dict[str, Any]]) -> str:
    if not meta:
        return f"Track {idx}"
    bits = []
    mark = meta.get("transcription_factor") or meta.get("histone_mark")
    if mark:
        bits.append(str(mark))
    name = meta.get("biosample_name") or meta.get("name")
    if name:
        bits.append(str(name))
    if meta.get("ontology_curie"):
        bits.append(str(meta["ontology_curie"]))
    if meta.get("Assay title"):
        bits.append(str(meta["Assay title"]))
    if meta.get("strand") and meta.get("strand") != ".":
        bits.append(f"strand {meta['strand']}")
    return " · ".join(bits) or f"Track {idx}"


def track_short_label(idx: int, meta: Optional[Dict[str, Any]]) -> str:
    if not meta:
        return f"T{idx}"
    name = str(meta.get("biosample_name") or meta.get("ontology_curie") or meta.get("name") or f"T{idx}")
    mark = meta.get("transcription_factor") or meta.get("histone_mark")
    if mark:
        name = f"{mark} · {name}"
    strand = meta.get("strand")
    return f"{name} ({strand})" if strand and strand != "." else name


class DatasetCatalog:
    def __init__(self) -> None:
        self._datasets: Dict[str, Dataset] = {}
        self._lock = threading.Lock()

    def add(self, path: Path, dataset_id: Optional[str] = None, annotations: Optional[Path] = None) -> Dataset:
        resolved = Path(path).expanduser().resolve()
        with self._lock:
            for existing in self._datasets.values():
                if existing.path == resolved:
                    return existing
        dataset = Dataset(resolved, dataset_id=dataset_id, annotations=annotations)
        with self._lock:
            base = dataset.id
            suffix = 2
            while dataset.id in self._datasets:
                dataset.id = f"{base}-{suffix}"
                suffix += 1
            self._datasets[dataset.id] = dataset
        return dataset

    def reload(self, dataset_id: str) -> Dataset:
        """Re-read a dataset from disk (new samples, genes or predictions), keeping its id."""
        old = self.get(dataset_id)
        fresh = Dataset(old.path, dataset_id=old.id, annotations=old.annotations_path)
        with self._lock:
            self._datasets[old.id] = fresh
        return fresh

    def remove(self, dataset_id: str) -> None:
        with self._lock:
            self._datasets.pop(dataset_id, None)

    def get(self, dataset_id: str) -> Dataset:
        dataset = self._datasets.get(dataset_id)
        if dataset is None:
            raise KeyError(f"Unknown dataset: {dataset_id}")
        return dataset

    def all(self) -> List[Dataset]:
        return list(self._datasets.values())


class DatasetMemory:
    """Datasets opened or imported from the UI, reopened when the visualizer starts again.

    Stored in ``~/.config/genomics/visualizer_datasets.json`` (``$XDG_CONFIG_HOME`` honoured).
    """

    def __init__(self, path: Optional[Path] = None):
        if path is None:
            import os

            base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
            path = Path(base) / "genomics" / "visualizer_datasets.json"
        self.path = Path(path)
        self._lock = threading.Lock()

    def entries(self) -> List[Dict[str, Optional[str]]]:
        try:
            data = load_json(self.path)
        except (OSError, ValueError):
            return []
        items = data.get("datasets") if isinstance(data, dict) else None
        return [i for i in items or [] if isinstance(i, dict) and i.get("path")]

    def _write(self, items: List[Dict[str, Optional[str]]]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps({"datasets": items}, indent=2) + "\n", encoding="utf-8")

    def remember(self, path: Path, annotations: Optional[Path]) -> None:
        with self._lock:
            items = [i for i in self.entries() if Path(str(i["path"])) != Path(path)]
            items.append({"path": str(path), "annotations": str(annotations) if annotations else None})
            self._write(items)

    def missing(self) -> List[str]:
        """Remembered datasets whose directory (or its dataset_metadata.json) is gone."""
        return [str(i["path"]) for i in self.entries() if not (Path(str(i["path"])) / "dataset_metadata.json").exists()]

    def forget(self, path: Path) -> None:
        with self._lock:
            self._write([i for i in self.entries() if Path(str(i["path"])) != Path(path)])


def to_jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value
