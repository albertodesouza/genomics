"""Region scalars: one number per sample from an AlphaGenome track over a region, as a sample field.

A scalar is defined by a window, an output, a track and a region (``[start, end)`` in reference offsets
of the window), e.g. "CAGE (+) of melanocyte summed over TSS ± 500 bp". Its value for a sample is the
track's mean over the region per haplotype (the per-haplotype region means of the Variant page, so the
two share one disk cache), averaged over H1 and H2, then summed over the region's length when asked,
and optionally log2(x + 1) transformed.

The value becomes a numeric sample field (a Samples column, exported with the CSV) and, when ``bins`` is
set, a categorical companion field ``<name>_bin`` with quantile bins ``Q1`` (lowest) ... ``Qk`` over all
samples of the dataset, so the scalar can filter the cohort and split group means on Tracks.

Definitions are kept per dataset path in ``~/.config/genomics/visualizer_scalars.json`` (or only in
memory, without a path) and reapplied from the cache whenever the dataset's samples are listed.
"""
from __future__ import annotations

import json
import math
import re
import threading
import weakref
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from genomics.visualizer.datasets import Dataset, load_json

STATS = ("mean", "sum")
TRANSFORMS = ("none", "log2p1")
MAX_BINS = 9
MAX_REGION = 2_000_000
NAME_RE = re.compile(r"^[a-z][a-z0-9_]{0,47}$")


class ScalarError(ValueError):
    pass


def bin_field(name: str) -> str:
    return f"{name}_bin"


class ScalarStore:
    """Scalar definitions per dataset path, persisted to a JSON file (or in memory without a path)."""

    def __init__(self, path: Optional[Path] = None):
        self.path = Path(path) if path is not None else None
        self._memory: Dict[str, List[Dict[str, Any]]] = {}
        self._lock = threading.Lock()

    @classmethod
    def default(cls) -> "ScalarStore":
        import os

        base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
        return cls(Path(base) / "genomics" / "visualizer_scalars.json")

    def _read(self) -> Dict[str, List[Dict[str, Any]]]:
        if self.path is None:
            return self._memory
        try:
            data = load_json(self.path)
        except (OSError, ValueError):
            return {}
        items = data.get("datasets") if isinstance(data, dict) else None
        return {str(k): [s for s in v if isinstance(s, dict)] for k, v in (items or {}).items() if isinstance(v, list)}

    def _write(self, data: Dict[str, List[Dict[str, Any]]]) -> None:
        if self.path is None:
            self._memory = data
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"datasets": data}, indent=2) + "\n", encoding="utf-8")
        tmp.replace(self.path)

    def list(self, dataset: Dataset) -> List[Dict[str, Any]]:
        with self._lock:
            return [dict(s) for s in self._read().get(str(dataset.path), [])]

    def put(self, dataset: Dataset, spec: Dict[str, Any]) -> None:
        with self._lock:
            data = self._read()
            items = [s for s in data.get(str(dataset.path), []) if s.get("name") != spec["name"]]
            data[str(dataset.path)] = items + [spec]
            self._write(data)

    def delete(self, dataset: Dataset, name: str) -> bool:
        with self._lock:
            data = self._read()
            items = data.get(str(dataset.path), [])
            kept = [s for s in items if s.get("name") != name]
            if len(kept) == len(items):
                return False
            data[str(dataset.path)] = kept
            self._write(data)
            return True


def scalar_values(spec: Dict[str, Any], region: Dict[str, np.ndarray]) -> Dict[str, float]:
    """Per-sample value of ``spec`` from per-haplotype region means (see ``GenotypeService.region_means``)."""
    track = int(spec["track"])
    means = region["means"]
    if track >= means.shape[2]:
        raise ScalarError(f"{spec['output']} has {means.shape[2]} tracks")
    haps = means[:, :, track].astype(np.float64)  # (n, 2)
    ok = np.isfinite(haps)
    with np.errstate(invalid="ignore", divide="ignore"):
        values = np.where(ok, haps, 0.0).sum(axis=1) / ok.sum(axis=1)
    if spec.get("stat") == "sum":
        values = values * (int(spec["end"]) - int(spec["start"]))
    if spec.get("transform") == "log2p1":
        with np.errstate(invalid="ignore"):
            values = np.log2(np.maximum(values, 0.0) + 1.0)
    out = {}
    for sample, value in zip(region["samples"].tolist(), values.tolist()):
        if math.isfinite(value):
            out[str(sample)] = float(f"{value:.6g}")
    return out


def quantile_bins(values: Dict[str, float], bins: int) -> Dict[str, Any]:
    """``{"labels": {sample: "Q1".."Qk"}, "edges": [...]}``; ties share a bin, so fewer bins may be used."""
    if bins < 2 or not values:
        return {"labels": {}, "edges": []}
    arr = np.asarray(list(values.values()), dtype=np.float64)
    edges = np.unique(np.quantile(arr, np.linspace(0.0, 1.0, bins + 1)[1:-1]))
    labels = {}
    for sample, value in values.items():
        labels[sample] = f"Q{int(np.searchsorted(edges, value, side='right')) + 1}"
    return {"labels": labels, "edges": [float(e) for e in edges]}


def histogram(values: List[float], n_bins: int = 30) -> Dict[str, Any]:
    arr = np.asarray([v for v in values if math.isfinite(v)], dtype=np.float64)
    if not arr.size:
        return {"n": 0, "edges": [], "counts": []}
    lo, hi = float(arr.min()), float(arr.max())
    if hi <= lo:
        hi = lo + (abs(lo) or 1.0) * 1e-6
    counts, edges = np.histogram(arr, bins=n_bins, range=(lo, hi))
    q = np.quantile(arr, [0.25, 0.5, 0.75])
    return {"n": int(arr.size), "edges": edges.tolist(), "counts": counts.tolist(), "mean": float(arr.mean()),
            "q1": float(q[0]), "median": float(q[1]), "q3": float(q[2]), "min": lo, "max": float(arr.max())}


class ScalarService:
    """Derived sample fields kept per dataset: region scalars (``kind`` "region", the default), genotype
    PCs (``kind`` "pca", fields ``pc1``..``pcK``, values from the ancestry PCA cache) and matched groups
    (``kind`` "match", one categorical field whose values are stored with the definition)."""

    def __init__(self, genotypes, store: Optional[ScalarStore] = None, ancestry=None):
        self.genotypes = genotypes
        self.ancestry = ancestry
        self.store = store or ScalarStore()
        # dataset -> {spec name: (spec JSON, field names set)}
        self._applied: "weakref.WeakKeyDictionary[Dataset, Dict[str, Any]]" = weakref.WeakKeyDictionary()
        self._lock = threading.Lock()

    # -- definitions -------------------------------------------------------------------
    def _own_fields(self, dataset: Dataset) -> set:
        with self._lock:
            return {f for _, fields in self._applied.get(dataset, {}).values() for f in fields}

    def check_name(self, dataset: Dataset, name: str, fields: List[str], replace: bool, what: str = "region scalar") -> None:
        """``name`` is a valid, free spec name and ``fields`` do not shadow the dataset's own fields."""
        if not NAME_RE.match(name):
            raise ScalarError("Name must start with a letter and use only a-z, 0-9 and _ (at most 48 characters)")
        specs = {sp["name"]: sp for sp in self.store.list(dataset)}
        taken = {f["name"] for f in dataset.fields} - self._own_fields(dataset)
        clash = [f for f in fields if f in taken]
        if clash:
            raise ScalarError(f"{clash[0]!r} is already a sample field")
        if name in specs and not replace:
            raise ScalarError(f"A {what} named {name!r} already exists")
        for other in specs.values():
            if other["name"] != name and set(fields) & set(self.spec_fields(other)):
                raise ScalarError(f"{other['name']!r} already defines {sorted(set(fields) & set(self.spec_fields(other)))[0]!r}")

    @staticmethod
    def spec_fields(spec: Dict[str, Any]) -> List[str]:
        kind = spec.get("kind", "region")
        if kind == "pca":
            return [f"pc{i + 1}" for i in range(int(spec.get("count") or 0))]
        if kind == "match":
            return [spec["name"]]
        return [spec["name"], bin_field(spec["name"])]

    def validate(self, dataset: Dataset, body: Dict[str, Any], gene_info: Dict[str, Any], replace: bool = False) -> Dict[str, Any]:
        if not isinstance(body, dict):
            raise ScalarError("Expected a JSON object")
        name = str(body.get("name") or "").strip().lower()
        self.check_name(dataset, name, [name, bin_field(name)], replace)
        output = str(body.get("output") or "")
        outputs = gene_info.get("outputs") or {}
        if output not in outputs:
            raise ScalarError(f"{gene_info.get('gene')} has no {output!r} predictions")
        try:
            track = int(body.get("track"))
            start, end = int(body.get("start")), int(body.get("end"))
            bins = int(body.get("bins") or 0)
        except (TypeError, ValueError):
            raise ScalarError("track, start, end and bins must be integers")
        n_tracks = len(outputs[output].get("tracks") or [])
        if not 0 <= track < n_tracks:
            raise ScalarError(f"{output} has {n_tracks} tracks")
        length = int(gene_info.get("length") or 0)
        if start < 0 or end <= start or (length and end > length):
            raise ScalarError(f"Region must lie within the window (0-{length})")
        if end - start > MAX_REGION:
            raise ScalarError("Region too long")
        stat = str(body.get("stat") or "mean")
        transform = str(body.get("transform") or "none")
        if stat not in STATS or transform not in TRANSFORMS:
            raise ScalarError(f"stat must be one of {STATS}, transform one of {TRANSFORMS}")
        if bins and not 2 <= bins <= MAX_BINS:
            raise ScalarError(f"bins must be 0 or 2-{MAX_BINS}")
        track_meta = outputs[output]["tracks"][track]
        return {
            "kind": "region", "name": name, "gene": gene_info["gene"], "output": output, "track": track,
            "track_label": str(body.get("track_label") or track_meta.get("label") or f"track {track}"),
            "start": start, "end": end, "region": str(body.get("region") or "")[:200],
            "stat": stat, "transform": transform, "bins": bins,
            "description": str(body.get("description") or "")[:300],
        }

    def specs(self, dataset: Dataset, kind: Optional[str] = None) -> List[Dict[str, Any]]:
        return [sp for sp in self.store.list(dataset) if kind is None or sp.get("kind", "region") == kind]

    def spec(self, dataset: Dataset, name: str) -> Dict[str, Any]:
        for spec in self.specs(dataset):
            if spec["name"] == name:
                return spec
        raise KeyError(f"No derived field named {name!r}")

    # -- values ------------------------------------------------------------------------
    def _region(self, dataset: Dataset, spec: Dict[str, Any], progress=None):
        args = (dataset, spec["gene"], spec["output"], int(spec["start"]), int(spec["end"]))
        if progress is None:
            return self.genotypes.cached_region_means(*args)
        return self.genotypes.region_means(*args, progress=progress)

    def region_key(self, dataset: Dataset, spec: Dict[str, Any]) -> str:
        return self.genotypes.region_key(dataset, spec["gene"], spec["output"], int(spec["start"]), int(spec["end"]))

    def is_cached(self, dataset: Dataset, spec: Dict[str, Any]) -> bool:
        return self._region(dataset, spec) is not None

    def fields_for(self, dataset: Dataset, spec: Dict[str, Any]) -> Optional[Dict[str, Tuple[Dict[str, Any], str]]]:
        """{field: (values by sample, kind)} from caches only; None while a cache is missing."""
        kind = spec.get("kind", "region")
        if kind == "match":
            return {spec["name"]: ({str(k): str(v) for k, v in (spec.get("values") or {}).items()}, "categorical")}
        if kind == "pca":
            if self.ancestry is None:
                return None
            result = self.ancestry.cached(dataset, spec["params"])
            if result is None:
                return None
            samples = result["samples"].tolist()
            scores = result["scores"]
            return {f"pc{i + 1}": ({s: float(f"{v:.6g}") for s, v in zip(samples, scores[:, i].tolist())}, "numeric")
                    for i in range(min(int(spec["count"]), scores.shape[1]))}
        region = self._region(dataset, spec)
        if region is None:
            return None
        values = scalar_values(spec, region)
        return {spec["name"]: (values, "numeric"),
                bin_field(spec["name"]): (quantile_bins(values, int(spec.get("bins") or 0))["labels"], "categorical")}

    def compute(self, dataset: Dataset, spec: Dict[str, Any], progress) -> None:
        """Apply a region scalar, computing its region means first when ``progress`` is given."""
        if progress is not None:
            self._region(dataset, spec, progress)
        fields = self.fields_for(dataset, spec)
        if fields is None:
            raise ScalarError(f"No cached values for {spec['name']}")
        self.apply(dataset, spec, fields)

    def apply(self, dataset: Dataset, spec: Dict[str, Any], fields: Dict[str, Tuple[Dict[str, Any], str]]) -> None:
        with self._lock:
            previous = self._applied.get(dataset, {}).get(spec["name"])
        for stale in (previous[1] if previous else set()) - set(fields):
            dataset.set_field(stale, {})
        for name, (values, kind) in fields.items():
            dataset.set_field(name, values, kind=kind)
        with self._lock:
            self._applied.setdefault(dataset, {})[spec["name"]] = (json.dumps(spec, sort_keys=True), {f for f, (v, _) in fields.items() if v})

    def unapply(self, dataset: Dataset, name: str) -> None:
        with self._lock:
            previous = self._applied.get(dataset, {}).pop(name, None)
        for field in previous[1] if previous else ():
            dataset.set_field(field, {})

    def save(self, dataset: Dataset, spec: Dict[str, Any]) -> None:
        """Store a definition whose values are available and apply it."""
        fields = self.fields_for(dataset, spec)
        if fields is None:
            raise ScalarError(f"No cached values for {spec['name']}")
        self.apply(dataset, spec, fields)
        self.store.put(dataset, spec)

    def sync(self, dataset: Dataset) -> None:
        """Apply stored definitions whose values are cached, and drop fields of deleted ones (cheap)."""
        specs = {sp["name"]: sp for sp in self.specs(dataset)}
        with self._lock:
            applied = {k: v[0] for k, v in self._applied.get(dataset, {}).items()}
        for name in set(applied) - set(specs):
            self.unapply(dataset, name)
        for name, spec in specs.items():
            if applied.get(name) == json.dumps(spec, sort_keys=True):
                continue
            try:
                fields = self.fields_for(dataset, spec)
            except Exception:  # a definition whose window or output disappeared stays listed, without values
                fields = None
            if fields is not None:
                self.apply(dataset, spec, fields)

    def is_applied(self, dataset: Dataset, spec: Dict[str, Any]) -> bool:
        with self._lock:
            entry = self._applied.get(dataset, {}).get(spec["name"])
        return bool(entry) and entry[0] == json.dumps(spec, sort_keys=True)

    def describe(self, dataset: Dataset, spec: Dict[str, Any]) -> Dict[str, Any]:
        ready = self.is_applied(dataset, spec)
        out = {k: v for k, v in spec.items() if k != "values"}
        out.update(ready=ready, kind=spec.get("kind", "region"), fields=self.spec_fields(spec))
        if out["kind"] == "region":
            out.update(field=spec["name"], bin_field=bin_field(spec["name"]) if spec.get("bins") else None)
            if ready:
                column = [row.get(spec["name"]) for row in dataset.samples]
                values = [float(v) for v in column if isinstance(v, (int, float))]
                out["histogram"] = histogram(values)
                out["edges"] = quantile_bins(dict(enumerate(values)), int(spec.get("bins") or 0))["edges"]
        elif out["kind"] == "match":
            counts: Dict[str, int] = {}
            for v in (spec.get("values") or {}).values():
                counts[str(v)] = counts.get(str(v), 0) + 1
            out["counts"] = counts
        return out
