"""Perturbation Lab: edit a haplotype in silico, re-predict it with AlphaGenome and re-score it with
a trained model, inside the visualizer process.

Any torch run under the runs roots can be loaded (its config, processed cache and checkpoint;
see :class:`~genomics.predictors.genotype_based.analysis.model_context.GenotypeModelContext`).
Edits are given in reference (genomic) window offsets and applied to each chosen haplotype's own
sequence through its indel map, so a selection made on the shared genomic axis edits the right
bases of every haplotype. Every operation keeps the haplotype length, so the edited prediction is
remapped onto genomic coordinates with the original haplotype's map:

``overwrite``  every base of the region becomes ``base``
``scramble``   the region's bases are shuffled (``seed``)
``reference``  bases aligned to the reference are reverted to it (removes SNVs; indels stay)
``sequence``   the region is replaced by ``sequence`` (same length as the haplotype region)

AlphaGenome is called once per edited haplotype for the model's input output plus any display
outputs, with the ontology terms of the tracks on disk, and the result is reordered into the
on-disk track order so the model sees exactly the columns it was trained on. Results are cached
on disk by sequence hash.
"""
from __future__ import annotations

import hashlib
import json
import threading
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeout
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.visualizer.cache import stable_key
from genomics.visualizer.datasets import Dataset, load_json
from genomics.visualizer.sequences import MAX_LETTER_SPAN
from genomics.visualizer.signals import _clamp_range, bin_matrix

OPS = ("overwrite", "scramble", "reference", "sequence")
HAPLOTYPES = ("H1", "H2")
MAX_EDITS = 64
PREDICT_TIMEOUT = 900.0
ProgressFn = Callable[[float, str], None]


class PerturbError(ValueError):
    """Bad request (unknown sample, invalid edit, no model loaded)."""


def _run_config(path: Path, cache: Dict[str, Tuple[int, Dict[str, Any]]]) -> Dict[str, Any]:
    config_path = path / "config.yaml"
    stat = config_path.stat()
    hit = cache.get(str(config_path))
    if hit and hit[0] == stat.st_mtime_ns:
        return hit[1]
    import yaml

    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except Exception:
        raw = {}
    cache[str(config_path)] = (stat.st_mtime_ns, raw)
    return raw


def _track_key(record: Dict[str, Any], with_name: bool) -> Tuple:
    base = (str(record.get("ontology_curie")), str(record.get("strand")))
    return base + ((str(record.get("name") or record.get("Assay title") or ""),) if with_name else ())


def reorder_columns(values: np.ndarray, predicted: List[Dict[str, Any]], canonical: List[Dict[str, Any]]) -> np.ndarray:
    """Reorder AlphaGenome columns into the on-disk track order (same ontology/strand/assay)."""
    with_name = all(r.get("name") or r.get("Assay title") for r in canonical) and all(r.get("name") or r.get("Assay title") for r in predicted)
    lookup: Dict[Tuple, List[int]] = {}
    for i, record in enumerate(predicted):
        lookup.setdefault(_track_key(record, with_name), []).append(i)
    columns = []
    used: Dict[Tuple, int] = {}
    for record in canonical:
        key = _track_key(record, with_name)
        options = lookup.get(key) or []
        k = used.get(key, 0)
        if k >= len(options):
            raise PerturbError(f"AlphaGenome returned no track for {key} (ontology/strand of the stored predictions)")
        columns.append(options[k])
        used[key] = k + 1
    return np.ascontiguousarray(values[:, columns], dtype=np.float32)


class PerturbService:
    def __init__(self, app, default_config: Optional[Path] = None):
        self.app = app
        self.default_config = Path(default_config) if default_config else None
        self.context = None
        self.context_id: Optional[Tuple[str, str]] = None
        self._lock = threading.RLock()
        self._configs: Dict[str, Tuple[int, Dict[str, Any]]] = {}
        self._edited: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
        self._results: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
        self._labels: Dict[str, Optional[str]] = {}

    # -- models ---------------------------------------------------------------------------------
    def models(self) -> Dict[str, Any]:
        runs = []
        for run_id, path in self.app.experiments.run_dirs().items():
            if not (path / "config.yaml").exists():
                continue
            checkpoints = sorted(p.name for p in (path / "models").glob("*.pt")) if (path / "models").is_dir() else []
            raw = _run_config(path, self._configs)
            di = raw.get("dataset_input") or {}
            model = raw.get("model") or {}
            output = raw.get("output") or {}
            model_type = str(model.get("type") or "")
            reasons = []
            if not checkpoints:
                reasons.append("no checkpoint")
            if model_type.upper() not in ("NN", "CNN", "CNN2"):
                reasons.append(f"{model_type or 'unknown'} model (the lab needs NN, CNN or CNN2)")
            if output.get("prediction_target") == "frog_likelihood":
                reasons.append("regression target")
            runs.append({
                "id": run_id,
                "name": path.name,
                "path": str(path),
                "model": model_type,
                "target": output.get("prediction_target"),
                "classes": output.get("known_classes"),
                "genes": di.get("genes_to_use") or [],
                "outputs": di.get("alphagenome_outputs") or [],
                "layout": di.get("tensor_layout"),
                "dataset": di.get("dataset_dir") or di.get("dataset_id"),
                "checkpoints": checkpoints,
                "default_checkpoint": "best_accuracy.pt" if "best_accuracy.pt" in checkpoints else (checkpoints[0] if checkpoints else None),
                "compatible": not reasons,
                "reasons": reasons,
                "mtime": path.stat().st_mtime,
            })
        runs.sort(key=lambda r: (not r["compatible"], -r["mtime"]))
        return {"runs": runs, "default": self._default_run(runs), "loaded": self.describe() if self.context else None}

    def _default_run(self, runs: List[Dict[str, Any]]) -> Optional[str]:
        compatible = [r for r in runs if r["compatible"]]
        if self.default_config and self.default_config.exists():
            try:
                from genomics.predictors.genotype_based.config import generate_experiment_name, load_config

                name = generate_experiment_name(load_config(self.default_config))
                for run in compatible:
                    if run["name"] == name:
                        return run["id"]
            except Exception:
                pass
        for run in compatible:
            if "pigmentation" in run["name"] and run["layout"] == "haplotype_channels":
                return run["id"]
        return compatible[0]["id"] if compatible else None

    def load(self, run_id: str, checkpoint: Optional[str], progress: ProgressFn) -> Dict[str, Any]:
        from genomics.predictors.genotype_based.analysis.model_context import GenotypeModelContext

        runs = self.app.experiments.run_dirs()
        path = runs.get(run_id)
        if path is None:
            raise PerturbError(f"Unknown run: {run_id}")
        checkpoint = checkpoint or "best_accuracy.pt"
        with self._lock:
            if self.context is not None and self.context_id == (run_id, checkpoint):
                return self.describe()
        context = GenotypeModelContext(path / "config.yaml", checkpoint, experiment_dir=path, progress=progress)
        with self._lock:
            self.context = context
            self.context_id = (run_id, checkpoint)
            self._labels = context.labels()
            self._edited.clear()
            self._results.clear()
        dataset = self.dataset()
        label_field = self._label_field(dataset)
        progress(1.0, "Model ready")
        return self.describe(label_field)

    def _require(self):
        if self.context is None:
            raise PerturbError("Load a model first")
        return self.context

    def dataset(self) -> Dataset:
        context = self._require()
        return self.app.catalog.add(Path(context.dataset_dir))

    def _label_field(self, dataset: Dataset) -> str:
        """The model's target as a sample field (derived targets are added to the dataset)."""
        context = self._require()
        target = context.target
        names = {f["name"] for f in dataset.fields}
        derived = (getattr(context.config.output, "derived_targets", None) or {})
        if target not in names or target in derived:
            dataset.set_field(target, {s: v for s, v in self._labels.items() if v is not None})
        return target

    def describe(self, label_field: Optional[str] = None) -> Dict[str, Any]:
        context = self._require()
        dataset = self.dataset()
        counts: Dict[str, int] = {}
        for value in self._labels.values():
            if value is not None:
                counts[value] = counts.get(value, 0) + 1
        samples = [
            {"id": s, "label": self._labels.get(s), "split": context.split_of(s)}
            for s in dataset.sample_index
            if s in self._labels
        ]
        return {
            "run": self.context_id[0],
            "checkpoint": self.context_id[1],
            "config": str(context.config_path),
            "model": context.config.model.type,
            "target": context.target,
            "label_field": label_field or context.target,
            "classes": context.class_names,
            "class_counts": counts,
            "genes": [g for g in context.genes if g in dataset.genes],
            "outputs": context.outputs,
            "ontology_terms": list(context.config.dataset_input.ontology_terms or []),
            "window_center_size": context.window_center_size,
            "dataset_id": dataset.id,
            "dataset_path": str(dataset.path),
            "samples": samples,
            "device": str(context.device),
        }

    # -- edits ----------------------------------------------------------------------------------
    @staticmethod
    def normalize_edits(edits: Any) -> List[Dict[str, Any]]:
        if not isinstance(edits, list):
            raise PerturbError("edits must be a list")
        if len(edits) > MAX_EDITS:
            raise PerturbError(f"At most {MAX_EDITS} edits")
        out = []
        for raw in edits:
            if not isinstance(raw, dict):
                raise PerturbError("each edit must be an object")
            op = raw.get("op")
            if op not in OPS:
                raise PerturbError(f"op must be one of {', '.join(OPS)}")
            start, end = int(raw.get("start", -1)), int(raw.get("end", -1))
            if start < 0 or end <= start:
                raise PerturbError("each edit needs 0 <= start < end (reference window offsets)")
            haps = [h for h in (raw.get("haplotypes") or list(HAPLOTYPES)) if h in HAPLOTYPES]
            if not haps:
                raise PerturbError("each edit needs haplotypes H1 and/or H2")
            edit = {"op": op, "start": start, "end": end, "haplotypes": haps}
            if op == "overwrite":
                base = str(raw.get("base") or "").upper()
                if base not in ("A", "C", "G", "T", "N"):
                    raise PerturbError("overwrite needs base A, C, G, T or N")
                edit["base"] = base
            elif op == "scramble":
                edit["seed"] = int(raw.get("seed") or 0)
            elif op == "sequence":
                seq = "".join(str(raw.get("sequence") or "").split()).upper()
                if not seq or set(seq) - set("ACGTN"):
                    raise PerturbError("sequence edits need A/C/G/T/N bases")
                edit["sequence"] = seq
            out.append(edit)
        return out

    def _haplotype(self, dataset: Dataset, sample: str, gene: str, hap: str) -> np.ndarray:
        return np.frombuffer(self.app.signals.haplotype_sequence(dataset, sample, gene, hap), dtype=np.uint8)

    def edited_haplotype(self, dataset: Dataset, sample: str, gene: str, hap: str, edits: List[Dict[str, Any]]) -> Tuple[np.ndarray, List[Dict[str, Any]]]:
        """Apply ``edits`` (reference offsets) to one haplotype; returns the sequence and what changed."""
        original = self._haplotype(dataset, sample, gene, hap)
        seq = original.copy()
        ref = np.frombuffer(self.app.signals.reference_sequence(dataset, gene), dtype=np.uint8)
        ref_map = self.app.signals.ref_map(dataset, sample, gene, hap)
        resolved = []
        for index, edit in enumerate(edits):
            if hap not in edit["haplotypes"]:
                continue
            start, end = _clamp_range(edit["start"], edit["end"], ref_map.local.size)
            local = ref_map.local[start:end]
            mapped = local[(local >= 0) & (local < seq.size)]
            if mapped.size == 0:
                raise PerturbError(f"Edit {index + 1}: the region is deleted on {sample} {hap}")
            lo, hi = int(mapped.min()), int(mapped.max()) + 1
            before = seq[lo:hi].copy()
            if edit["op"] == "overwrite":
                seq[lo:hi] = ord(edit["base"])
            elif edit["op"] == "scramble":
                seq[lo:hi] = np.random.default_rng(edit["seed"]).permutation(seq[lo:hi])
            elif edit["op"] == "reference":
                positions = np.arange(start, end)
                ok = (local >= 0) & (local < seq.size) & (positions < ref.size)
                seq[local[ok]] = ref[positions[ok]]
            else:
                replacement = np.frombuffer(edit["sequence"].encode("ascii"), dtype=np.uint8)
                if replacement.size != hi - lo:
                    raise PerturbError(f"Edit {index + 1}: the sequence has {replacement.size} bases but the {hap} region has {hi - lo}")
                seq[lo:hi] = replacement
            resolved.append({
                "index": index,
                "haplotype": hap,
                "op": edit["op"],
                "start": start,
                "end": end,
                "local_start": lo,
                "local_end": hi,
                "changed": int(np.count_nonzero(before != seq[lo:hi])),
            })
        return seq, resolved

    def sequence_view(self, sample: str, gene: str, haplotypes: Sequence[str], edits: Any, start: int, end: int, bins: int) -> Dict[str, Any]:
        dataset = self.dataset()
        if gene not in dataset.genes:
            raise PerturbError(f"Unknown gene: {gene}")
        edits = self.normalize_edits(edits or [])
        domain = self.app.signals.reference_length(dataset, gene)
        start, end = _clamp_range(start, end, domain)
        span = end - start
        letters = span <= MAX_LETTER_SPAN
        ref = np.frombuffer(self.app.signals.reference_sequence(dataset, gene), dtype=np.uint8)
        ref_row = np.full(span, ord("N"), dtype=np.uint8)
        hi_ref = min(end, ref.size)
        if hi_ref > start:
            ref_row[: hi_ref - start] = ref[start:hi_ref]
        rows = []
        for hap in [h for h in haplotypes if h in HAPLOTYPES]:
            ref_map = self.app.signals.ref_map(dataset, sample, gene, hap)
            original = self._haplotype(dataset, sample, gene, hap)
            edited, resolved = self.edited_haplotype(dataset, sample, gene, hap, edits)
            local = ref_map.local[start:end]
            valid = (local >= 0) & (local < original.size)
            orig_row = np.full(span, ord("-"), dtype=np.uint8)
            edit_row = orig_row.copy()
            orig_row[valid] = original[local[valid]]
            edit_row[valid] = edited[local[valid]]
            changed = orig_row != edit_row
            variant = (orig_row != ref_row) & valid & (orig_row != ord("N"))
            item: Dict[str, Any] = {"haplotype": hap, "changed": int(changed.sum()), "edits": resolved}
            if letters:
                item["original"] = orig_row.tobytes().decode("ascii", "replace")
                item["edited"] = edit_row.tobytes().decode("ascii", "replace")
            else:
                edges = np.unique(np.linspace(0, span, min(bins, span) + 1).astype(np.int64))
                widths = np.diff(edges).astype(np.float32)
                item["changed_density"] = (np.add.reduceat(changed.astype(np.float32), edges[:-1]) / widths).astype(np.float32)
                item["variant_density"] = (np.add.reduceat(variant.astype(np.float32), edges[:-1]) / widths).astype(np.float32)
            rows.append(item)
        payload: Dict[str, Any] = {"gene": gene, "sample": sample, "start": start, "end": end, "domain": domain, "mode": "letters" if letters else "density", "rows": rows}
        if letters:
            payload["reference"] = ref_row.tobytes().decode("ascii", "replace")
        else:
            payload["edges"] = (np.unique(np.linspace(0, span, min(bins, span) + 1).astype(np.int64)) + start).astype(np.int64)
        return payload

    # -- prediction -----------------------------------------------------------------------------
    def score(self, sample: str) -> Dict[str, Any]:
        context = self._require()
        if sample not in self._labels:
            raise PerturbError(f"Unknown sample for this model's dataset: {sample}")
        with self._lock:
            probs = context.baseline(sample)
        return {"sample": sample, "label": self._labels.get(sample), "split": context.split_of(sample), "classes": context.class_names, "probabilities": probs}

    def _canonical(self, dataset: Dataset, sample: str, gene: str, hap: str, output: str) -> List[Dict[str, Any]]:
        path = dataset.prediction_path(sample, gene, hap, output).with_name(f"{output}_metadata.json")
        if not path.exists():
            raise PerturbError(f"No track metadata for {sample}/{gene}/{hap}/{output}; cannot match AlphaGenome tracks")
        payload = load_json(path)
        records = payload.get("metadata", payload) if isinstance(payload, dict) else payload
        if not isinstance(records, list) or not records:
            raise PerturbError(f"Empty track metadata: {path}")
        return records

    def _predict(self, sequence: bytes, outputs: List[str], terms: List[str]) -> Dict[str, Tuple[np.ndarray, List[Dict[str, Any]]]]:
        """AlphaGenome prediction for one sequence (disk-cached by content)."""
        key = hashlib.sha1(sequence + json.dumps([sorted(outputs), sorted(terms)]).encode()).hexdigest()
        cache_dir = Path(self.app.cache_dir) / "perturb" if self.app.cache_dir else None
        if cache_dir is not None and (cache_dir / f"{key}.npz").exists():
            try:
                with np.load(cache_dir / f"{key}.npz") as data:
                    meta = json.loads((cache_dir / f"{key}.json").read_text(encoding="utf-8"))
                    return {o: (np.asarray(data[o]), meta[o]) for o in outputs}
            except Exception:
                pass
        backend = getattr(self.app, "alphagenome", None)
        if backend is None:
            raise PerturbError("No AlphaGenome backend configured")
        reasons = backend.reasons()
        if reasons:
            raise PerturbError(f"AlphaGenome backend not ready: {'; '.join(reasons)}")
        from alphagenome.models import dna_client

        client = backend.create_client()
        requested = [getattr(dna_client.OutputType, o.upper()) for o in outputs]
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(client.predict_sequence, sequence.decode("ascii"), requested_outputs=requested, ontology_terms=terms or None)
            try:
                prediction = future.result(timeout=PREDICT_TIMEOUT)
            except FutureTimeout:
                raise RuntimeError(f"AlphaGenome did not answer within {PREDICT_TIMEOUT:.0f}s")
        from genomics.workflows.alphagenome.predict_dataset import metadata_records

        result = {}
        for output in outputs:
            track = getattr(prediction, output, None)
            if track is None:
                raise RuntimeError(f"AlphaGenome returned no {output}")
            result[output] = (np.asarray(track.values, dtype=np.float32), metadata_records(track))
        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
            tmp = cache_dir / f".{key}.tmp.npz"
            np.savez_compressed(tmp, **{o: v for o, (v, _m) in result.items()})
            tmp.replace(cache_dir / f"{key}.npz")
            (cache_dir / f"{key}.json").write_text(json.dumps({o: m for o, (_v, m) in result.items()}), encoding="utf-8")
        return result

    def apply_key(self, sample: str, gene: str, edits: List[Dict[str, Any]], outputs: List[str]) -> str:
        dataset = self.dataset()
        return stable_key({"v": 1, "run": self.context_id, "dataset": dataset.fingerprint, "sample": sample, "gene": gene, "edits": edits, "outputs": sorted(outputs)})

    def apply(self, sample: str, gene: str, edits: List[Dict[str, Any]], outputs: List[str], progress: ProgressFn) -> Dict[str, Any]:
        context = self._require()
        dataset = self.dataset()
        if sample not in self._labels:
            raise PerturbError(f"Unknown sample for this model's dataset: {sample}")
        if gene not in context.genes:
            raise PerturbError(f"{gene} is not one of the model's genes")
        if not edits:
            raise PerturbError("Add at least one edit")
        key = self.apply_key(sample, gene, edits, outputs)
        started = time.time()
        haps = sorted({h for e in edits for h in e["haplotypes"]})
        wanted = list(dict.fromkeys([*context.outputs, *outputs]))
        on_disk = set(dataset.gene_info(gene)["outputs"])
        wanted = [o for o in wanted if o in on_disk or o in context.outputs]
        overrides: Dict[Tuple[str, str], Dict[str, Tuple[np.ndarray, Optional[list]]]] = {}
        edited_arrays: Dict[str, Dict[str, np.ndarray]] = {}
        resolved_all: List[Dict[str, Any]] = []
        for i, hap in enumerate(haps):
            progress(0.05 + 0.85 * i / len(haps), f"Re-predicting {sample} {gene} {hap} with AlphaGenome")
            seq, resolved = self.edited_haplotype(dataset, sample, gene, hap, edits)
            resolved_all.extend(resolved)
            canonical = {o: self._canonical(dataset, sample, gene, hap, o) for o in wanted}
            terms = sorted({str(r.get("ontology_curie")) for records in canonical.values() for r in records if r.get("ontology_curie")})
            predicted = self._predict(seq.tobytes(), wanted, terms)
            edited_arrays[hap] = {}
            for output in wanted:
                values, records = predicted[output]
                matrix = reorder_columns(values, records, canonical[output])
                edited_arrays[hap][output] = matrix
                if output in context.outputs:
                    overrides.setdefault((gene, hap), {})[output] = (matrix, canonical[output])
        progress(0.92, "Scoring with the model")
        with self._lock:
            baseline = context.baseline(sample)
            edited = context.score(sample, overrides)
        window = dataset.model_window(gene, context.window_center_size) or {}
        for item in resolved_all:
            item["inside_model_window"] = bool(window) and item["end"] > window.get("start", 0) and item["start"] < window.get("end", 0)
        classes = context.class_names
        result = {
            "key": key,
            "sample": sample,
            "gene": gene,
            "haplotypes": haps,
            "label": self._labels.get(sample),
            "split": context.split_of(sample),
            "classes": classes,
            "baseline": baseline,
            "edited": edited,
            "delta": edited - baseline,
            "predicted_baseline": classes[int(np.argmax(baseline))],
            "predicted_edited": classes[int(np.argmax(edited))],
            "edits": resolved_all,
            "outputs": wanted,
            "model_window": window,
            "elapsed": round(time.time() - started, 2),
        }
        with self._lock:
            self._edited[key] = {"sample": sample, "gene": gene, "arrays": edited_arrays}
            self._edited.move_to_end(key)
            while len(self._edited) > 6:
                self._edited.popitem(last=False)
            self._results[key] = result
            while len(self._results) > 64:
                self._results.popitem(last=False)
        return result

    def result(self, key: str) -> Optional[Dict[str, Any]]:
        return self._results.get(key)

    def signal(self, key: str, output: str, tracks: List[int], start: int, end: int, bins: int) -> Dict[str, Any]:
        entry = self._edited.get(key)
        if entry is None:
            raise PerturbError("This edit is no longer in memory; apply it again")
        dataset = self.dataset()
        sample, gene = entry["sample"], entry["gene"]
        domain = self.app.signals.reference_length(dataset, gene)
        start, end = _clamp_range(start, end, domain)
        bins = max(1, min(int(bins), 8192))
        series = []
        for hap, arrays in entry["arrays"].items():
            matrix = arrays.get(output)
            if matrix is None:
                continue
            baseline = self.app.signals.haplotype_window(dataset, sample, gene, hap, output, "reference", start, end, tracks)
            ref_map = self.app.signals.ref_map(dataset, sample, gene, hap)
            local = ref_map.local[start:end]
            cols = matrix[:, tracks] if tracks else matrix
            res = self.app.signals.resolution(dataset, gene, output)
            edited = np.full((end - start, cols.shape[1]), np.nan, dtype=np.float32)
            valid = (local >= 0) & (local < cols.shape[0] * res)
            edited[np.nonzero(valid)[0]] = cols[local[valid] // res]
            for kind, values in (("baseline", baseline), ("edited", edited)):
                binned = bin_matrix(values, bins)
                series.append({"haplotype": hap, "kind": kind, "mean": binned["mean"], "min": binned["min"], "max": binned["max"]})
        if not series:
            raise PerturbError(f"No edited {output} prediction for this edit")
        edges = bin_matrix(np.zeros((end - start, 1), np.float32), bins)["edges"]
        return {"key": key, "output": output, "tracks": tracks, "start": start, "end": end, "domain": domain, "edges": (edges + start).astype(np.int64), "series": series}
