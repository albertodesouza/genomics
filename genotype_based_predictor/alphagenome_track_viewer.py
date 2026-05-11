from __future__ import annotations

import argparse
import json
import math
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse

import numpy as np

from genotype_based_predictor.dynamic_indel_alignment import DynamicIndelAligner


DEFAULT_POINTS = 1000
MAX_POINTS = 20_000
MAX_SERIES = 24
DEFAULT_VIEW_LENGTH = 32768


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
    for source_raw, target_raw in zip(copy_from, copy_to):
        source = int(source_raw)
        target = int(target_raw)
        if start <= target + 1 <= end and 0 <= source < signal.size:
            window[target + 1 - start] = signal[source]
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
        xs.append(start + (left + right - 1) // 2)
        ys.append(_json_scalar(np.nanmean(bucket)))
    return xs, ys, "mean"


class AlphaGenomeRepository:
    def __init__(self, dataset_dir: Path):
        self.dataset_dir = Path(dataset_dir).resolve()
        metadata_path = self.dataset_dir / "dataset_metadata.json"
        if not metadata_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json not found: {metadata_path}")
        self.metadata = _load_json(metadata_path)
        self.individuals = self._load_individuals()
        self.individual_rows = self._load_individual_rows()
        self.genes = self._load_genes()
        self.outputs = self._load_outputs()
        self._aligner: Optional[DynamicIndelAligner] = None

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

    def options_payload(self) -> Dict[str, Any]:
        default_window = self._default_view_window()
        return {
            "dataset_dir": str(self.dataset_dir),
            "dataset_name": self.metadata.get("dataset_name"),
            "samples": self.individuals,
            "individuals": self.individual_rows,
            "populations": self._distribution("population"),
            "superpopulations": self._distribution("superpopulation"),
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
            self._aligner = DynamicIndelAligner(self.dataset_dir, selected_sample_ids=self.individuals, center_window_size=None)
        return self._aligner

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
        if aligner is not None:
            aligner.build_alignment_axis_for_gene(gene, self.individuals)
        expanded_length = aligner.get_expanded_length(gene) if aligner else None
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
                    entry = aligner.get_haplotype_entry(gene, sample, haplotype)
                    if entry is None:
                        raise RuntimeError(f"Alignment entry missing for {sample} {haplotype}")
                    signal = _align_signal_window(full_signal, entry, int(expanded_length), start, length)
                    x_start = start
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
            "x_axis": "expanded_alignment" if align else "sample_sequence",
            "expanded_length": expanded_length,
            "individuals": individual_summaries,
            "series": series,
        }


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
    .side-panel { background:#0c1220; border:1px solid #273854; border-radius:16px; padding:12px; overflow:auto; max-height:min(68vh,760px); }
    .person { border-bottom:1px solid #273854; padding:8px 0; }
    .person b { color:#ffd166; }
    .kv { display:grid; grid-template-columns:110px 1fr; gap:5px; font-size:.82rem; color:#dbe7f3; }
    .kv span:first-child { color:var(--muted); }
    svg { display: block; width: 100%; height: 100%; }
    .meta { color: var(--muted); font-size: .88rem; margin-top: 12px; white-space: pre-wrap; overflow-wrap: anywhere; }
    @media (max-width: 1000px) { .viz-grid { grid-template-columns:1fr; } }
    @media (max-width: 900px) { form { grid-template-columns: repeat(2, 1fr); } .wide { grid-column: span 2; } .stats { grid-template-columns: repeat(2, 1fr); } }
  </style>
</head>
<body>
<main>
  <h1>AlphaGenome Track Viewer</h1>
  <div class="sub">Compare AlphaGenome tracks across samples and haplotypes, optionally aligned by DynamicIndelAligner.</div>
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
      <label>Start<input id="start" type="number" value="1" min="1"></label>
      <label>Length<input id="length" type="number" value="2000" min="1"></label>
      <label>Points<input id="points" type="number" value="1000" min="1" max="20000"></label>
      <label>Alignment<div class="checks"><label><input id="align" type="checkbox"> align with DynamicIndelAligner</label></div></label>
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
    <div class="viz-grid">
      <div class="plot-wrap"><svg id="plot" role="img" aria-label="Signal plot"></svg></div>
      <aside class="side-panel">
        <h2 style="margin-top:0">Selected individuals</h2>
        <div id="individualPanel">No data loaded.</div>
      </aside>
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
  autoLoadTimer = setTimeout(() => loadTrack().catch(showError), delay);
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
  loadTrack().catch(showError);
}

function zoom(factor) {
  const start = Math.max(1, Number($('start').value || 1));
  const length = Math.max(1, Number($('length').value || 1));
  const center = start + Math.floor(length / 2);
  const nextLength = Math.max(10, Math.round(length * factor));
  $('length').value = nextLength;
  $('start').value = Math.max(1, center - Math.floor(nextLength / 2));
  loadTrack().catch(showError);
}

function centerDefaultWindow() {
  $('start').value = defaultWindow.start;
  $('length').value = defaultWindow.length;
  loadTrack().catch(showError);
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
  renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`);
  renderPopulations();
  await loadTrackOptions();
  if (selectedSamples().length && $('gene').value) await loadTrack();
}

$('controls').addEventListener('submit', (event) => { event.preventDefault(); loadTrack().catch(showError); });
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
['gene','output','hapH1','hapH2'].forEach(id => {
  $(id).addEventListener('change', () => loadTrackOptions().then(() => scheduleLoad()).catch(showError));
});
['track','start','length','points','align'].forEach(id => {
  $(id).addEventListener('change', () => scheduleLoad());
});
function showError(err) { $('status').className = 'error'; $('status').textContent = err.message || String(err); }
init().catch(showError);
</script>
</body>
</html>
"""


def serve(dataset_dir: Path, host: str, port: int) -> None:
    repository = AlphaGenomeRepository(dataset_dir)
    AlphaGenomeTrackViewerHandler.repository = repository
    server = ThreadingHTTPServer((host, port), AlphaGenomeTrackViewerHandler)
    print(f"AlphaGenome Track Viewer: http://{host}:{port}")
    print(f"Dataset: {repository.dataset_dir}")
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
    args = parser.parse_args(argv)
    serve(args.dataset_dir, args.host, args.port)


if __name__ == "__main__":
    main()
