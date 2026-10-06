"""Builds the detached tasks the visualizer can start: dataset import, AlphaGenome predictions,
training and evaluation (see :mod:`genomics.visualizer.tasks`).

Each builder validates the form values, writes the inputs the step needs (import spec, sample
list, training config) into the task directory and returns the commands to run. The commands are
the same ``python -m`` modules the ``genomics`` CLI runs, so every task is reproducible from its
``task.json`` without the visualizer.
"""
from __future__ import annotations

import copy
import json
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

from genomics.visualizer.datasets import Dataset, slugify
from genomics.workflows.alphagenome.outputs import ALL_OUTPUTS, OUTPUT_SPECS, normalize_outputs

IMPORT_MODULE = "genomics.workflows.dataset_builders.vcf_import"
PREDICT_MODULE = "genomics.workflows.alphagenome.predict_dataset"
TRAIN_MODULE = "genomics.predictors.genotype_based.experiments.train"
EVALUATE_MODULE = "genomics.predictors.genotype_based.experiments.evaluate_checkpoint"
PREDICT_OUTPUTS = ALL_OUTPUTS
WINDOW_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
WINDOW_SIZES = (16384, 131072, 524288, 1048576)
MODEL_TYPES = ("CNN2", "CNN", "NN", "LOGREG", "RF", "SVM", "XGBOOST")
DEFAULT_BASE_CONFIG = "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml"
SPLITS = ("train", "val", "test")


class LaunchError(ValueError):
    """Invalid task parameters (reported to the user)."""


def _python() -> str:
    return sys.executable


def _repo_root() -> Path:
    from genomics.workspace import repo_root

    return repo_root()


def _strings(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        value = re.split(r"[\s,;]+", value)
    return [str(v).strip() for v in value if str(v).strip()]


# ---------------------------------------------------------------------------------------- import
def default_reference() -> Optional[str]:
    try:
        from genomics.core.reference_registry import resolve_reference

        ref = resolve_reference()
        return str(ref.path) if Path(ref.path).exists() else None
    except Exception:
        return None


def default_gtf() -> Optional[str]:
    try:
        from genomics.workspace import DEFAULT_DATASET_DIR

        candidate = Path(DEFAULT_DATASET_DIR) / "gtf_cache.feather"
        return str(candidate) if candidate.exists() else None
    except Exception:
        return None


def default_import_root() -> str:
    try:
        from genomics.workspace import data_root

        return str(Path(data_root()) / "imported")
    except Exception:
        return str(Path.home() / "genomics_datasets")


def import_defaults() -> Dict[str, Any]:
    return {
        "reference_fasta": default_reference(),
        "gtf": default_gtf(),
        "output_root": default_import_root(),
        "window_sizes": list(WINDOW_SIZES),
        "window_size": 524288,
        "outputs": list(PREDICT_OUTPUTS),
        "output_specs": import_defaults_outputs(),
        "region_presets": region_presets(),
    }


def import_defaults_outputs() -> Dict[str, Dict[str, Any]]:
    return {name: spec.as_dict() for name, spec in OUTPUT_SPECS.items()}


def region_presets() -> List[Dict[str, Any]]:
    from genomics.workflows.alphagenome.regulatory_regions import preset_regions

    return preset_regions()


def import_task(body: Dict[str, Any], alphagenome_ready: bool) -> Tuple[str, List[Dict[str, Any]], Dict[str, Any], Dict[str, str]]:
    """Returns (title, steps, params, files) for a VCF import (optionally followed by predictions)."""
    name = str(body.get("name") or "").strip()
    if not name:
        raise LaunchError("Give the dataset a name")
    output_dir = str(body.get("output_dir") or "").strip() or str(Path(default_import_root()) / slugify(name))
    out = Path(output_dir).expanduser()
    if out.exists() and any(out.iterdir()) and not (out / "import_spec.json").exists():
        raise LaunchError(f"{out} exists and is not an imported dataset; choose an empty or new directory")
    vcf = str(body.get("vcf") or "").strip()
    reference = str(body.get("reference_fasta") or "").strip()
    if not vcf or not reference:
        raise LaunchError("The VCF and the reference FASTA are required")
    from genomics.workflows.dataset_builders.vcf_import.builder import reference_available

    if not reference_available(reference):
        raise LaunchError(f"Reference FASTA not found: {reference}")
    window_size = int(body.get("window_size") or 524288)
    if window_size not in WINDOW_SIZES:
        raise LaunchError(f"Window size must be an AlphaGenome input length: {', '.join(map(str, WINDOW_SIZES))}")
    genes = _strings(body.get("genes"))
    regions = parse_regions(body.get("regions"))
    if not genes and not regions:
        raise LaunchError("Add at least one gene or region")
    samples = _strings(body.get("samples"))
    metadata = body.get("metadata") or {}
    if not isinstance(metadata, dict):
        raise LaunchError("metadata must map sample ids to fields")
    spec = {
        "name": name,
        "output_dir": str(out),
        "vcf": vcf,
        "vcf_overrides": {str(k): str(v) for k, v in (body.get("vcf_overrides") or {}).items()},
        "reference_fasta": reference,
        "window_size": window_size,
        "genes": genes,
        "regions": regions,
        "gtf": str(body.get("gtf") or "").strip() or None,
        "samples": samples,
        "metadata": {str(s): {str(k): v for k, v in (fields or {}).items() if v not in (None, "")} for s, fields in metadata.items()},
        "variant_filter": "snps" if body.get("variant_filter") == "snps" else "all",
        "workers": int(body.get("workers") or 8),
        "keep_raw": True,
    }
    steps = [{"title": "Build windows from the VCF", "command": [_python(), "-m", IMPORT_MODULE, "--spec", "{task_dir}/import_spec.json"], "weight": 1.0}]
    params = {"name": name, "output_dir": str(out), "vcf": vcf, "samples": len(samples) or None, "genes": genes + [r["name"] for r in regions], "window_size": window_size}
    predict = body.get("predict") or {}
    if predict.get("enabled"):
        outputs, terms = _predict_outputs(predict)
        command = [_python(), "-m", PREDICT_MODULE, str(out), "--outputs", ",".join(outputs), *_ontology_args(terms)]
        steps.append({"title": "AlphaGenome predictions", "command": command, "weight": 3.0})
        params["predict"] = {"outputs": outputs, "ontology_terms": terms or "all", "backend_ready": alphagenome_ready}
    return f"Import {name}", steps, params, {"import_spec.json": json.dumps(spec, indent=2)}


def parse_regions(raw_regions: Any) -> List[Dict[str, Any]]:
    """``NAME=chr:start[-end]`` strings or ``{name, chrom, start, end}`` dicts (1-based, inclusive)."""
    regions = []
    for raw in raw_regions or []:
        match = re.match(r"^\s*([\w.\-]+)\s*=\s*([\w.]+):(\d[\d,]*)(?:-(\d[\d,]*))?\s*$", raw) if isinstance(raw, str) else None
        if isinstance(raw, dict):
            try:
                region = {"name": str(raw.get("name") or "").strip(), "chrom": str(raw.get("chrom") or "").strip(), "start": int(raw["start"]), "end": int(raw.get("end") or raw["start"])}
            except (KeyError, TypeError, ValueError):
                raise LaunchError(f"Region needs name, chrom and start: {raw!r}")
        elif match:
            start = int(match.group(3).replace(",", ""))
            end = int((match.group(4) or match.group(3)).replace(",", ""))
            region = {"name": match.group(1), "chrom": match.group(2), "start": start, "end": end}
        elif str(raw).strip():
            raise LaunchError(f"Region must look like NAME=chr:start-end, got {raw!r}")
        else:
            continue
        if not WINDOW_NAME_RE.match(region["name"]) or not region["chrom"]:
            raise LaunchError(f"Region names may use letters, digits, '_', '.' and '-': {region['name']!r}")
        if region["end"] < region["start"] or region["start"] < 1:
            raise LaunchError(f"Region {region['name']}: end must be >= start >= 1")
        regions.append(region)
    return regions


def _predict_outputs(body: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    """(outputs, ontology CURIEs); no CURIEs means every tissue / cell type."""
    try:
        outputs = normalize_outputs(_strings(body.get("outputs")))
    except ValueError as exc:
        raise LaunchError(str(exc))
    if not outputs:
        raise LaunchError(f"Choose outputs among {', '.join(PREDICT_OUTPUTS)}")
    if body.get("all_tissues"):
        return outputs, []
    terms = _strings(body.get("ontology_terms"))
    if not terms:
        if all(not OUTPUT_SPECS[o].tissue_specific for o in outputs):
            return outputs, []
        raise LaunchError("Choose at least one tissue / cell type (ontology CURIE, e.g. CL:1000458)")
    bad = [t for t in terms if not re.match(r"^[A-Za-z]+:\d+$", t)]
    if bad:
        raise LaunchError(f"Ontology terms look like CL:1000458 or UBERON:0002107; got {', '.join(bad)}")
    return outputs, list(dict.fromkeys(terms))


def _ontology_args(terms: List[str]) -> List[str]:
    return ["--ontology", ",".join(terms)] if terms else ["--all-tissues"]


# ---------------------------------------------------------------------------- adding windows
def extension_source(dataset: Dataset) -> Dict[str, Any]:
    """Where new windows for ``dataset`` come from: its VCF, reference and gene table.

    ``problems`` lists what is missing (then only existing windows can be predicted).
    """
    spec: Dict[str, Any] = {}
    spec_path = dataset.path / "import_spec.json"
    if spec_path.exists():
        try:
            spec = json.loads(spec_path.read_text(encoding="utf-8"))
        except ValueError:
            spec = {}
    raw = dataset.metadata.get("raw_variant_sources") or {}
    vcf = str(spec.get("vcf") or raw.get("vcf_pattern") or raw.get("vcf_path") or "")
    reference = str(spec.get("reference_fasta") or raw.get("reference_fasta") or default_reference() or "")
    gtf = spec.get("gtf") or next((str(dataset.path / n) for n in ("gtf_cache.feather", "gtf_cache.parquet") if (dataset.path / n).exists()), None) or default_gtf()
    window_size = int(dataset.metadata.get("window_size") or spec.get("window_size") or 0)
    problems = []
    if not vcf:
        problems.append("the dataset does not record its source VCF")
    else:
        from genomics.workflows.dataset_builders.vcf_import.builder import resolve_vcf_files

        if not resolve_vcf_files(vcf):
            problems.append(f"source VCF not found: {vcf}")
    from genomics.workflows.dataset_builders.vcf_import.builder import reference_available

    if not reference or not reference_available(reference):
        problems.append("reference FASTA not found" + (f": {reference}" if reference else ""))
    if window_size not in WINDOW_SIZES:
        problems.append(f"window size {window_size or 'unknown'} is not an AlphaGenome input length")
    from genomics.workflows.dataset_builders.vcf_import.builder import missing_tools

    missing = missing_tools()
    if missing:
        problems.append(f"missing tools: {', '.join(missing)}")
    return {
        "vcf": vcf or None,
        "vcf_overrides": spec.get("vcf_overrides") or raw.get("vcf_overrides") or {},
        "reference_fasta": reference or None,
        "gtf": str(gtf) if gtf else None,
        "window_size": window_size or None,
        "variant_filter": spec.get("variant_filter") or "all",
        "problems": problems,
        "available": not problems,
    }


def plan_new_windows(dataset: Dataset, genes: List[str], regions: List[Dict[str, Any]]) -> Tuple[List[str], List[Dict[str, Any]], List[str]]:
    """(genes to build, regions to build, names already in the dataset), validated against the gene table."""
    existing = set(dataset.genes)
    known = [g for g in genes if g in existing] + [r["name"] for r in regions if r["name"] in existing]
    genes = [g for g in dict.fromkeys(genes) if g not in existing]
    regions = [r for r in regions if r["name"] not in existing]
    if not genes and not regions:
        return [], [], known
    source = extension_source(dataset)
    if source["problems"]:
        raise LaunchError(f"Cannot build new windows for this dataset: {'; '.join(source['problems'])}")
    from genomics.workflows.dataset_builders.vcf_import.builder import ImportSpecError, fasta_chrom_names, reference_fai, resolve_targets

    try:
        fai = reference_fai(source["reference_fasta"])
    except ImportSpecError as exc:
        raise LaunchError(str(exc))
    if fai.exists():
        gtf = Path(source["gtf"]).expanduser() if source["gtf"] else None
        try:
            targets = resolve_targets({"window_size": source["window_size"], "genes": genes, "regions": regions}, fasta_chrom_names(fai), gtf)
        except ImportSpecError as exc:
            raise LaunchError(str(exc))
        clashes = [t.name for t in targets if t.name in existing]
        if clashes:
            raise LaunchError(f"Already windows in this dataset: {', '.join(clashes)}")
    return genes, regions, known


def extension_spec(dataset: Dataset, genes: List[str], regions: List[Dict[str, Any]]) -> Dict[str, Any]:
    source = extension_source(dataset)
    return {
        "extend": True,
        "output_dir": str(dataset.path),
        "vcf": source["vcf"],
        "vcf_overrides": source["vcf_overrides"],
        "reference_fasta": source["reference_fasta"],
        "gtf": source["gtf"],
        "window_size": source["window_size"],
        "genes": genes,
        "regions": [{k: r[k] for k in ("name", "chrom", "start", "end")} for r in regions],
        "samples": [row["sample_id"] for row in dataset.samples],
        "variant_filter": source["variant_filter"],
        "workers": 8,
        "keep_raw": True,
    }


# ------------------------------------------------------------------------------------ predictions
def predict_task(dataset: Dataset, body: Dict[str, Any]) -> Tuple[str, List[Dict[str, Any]], Dict[str, Any], Dict[str, str]]:
    outputs, terms = _predict_outputs(body)
    new_genes, new_regions, known = plan_new_windows(dataset, _strings(body.get("new_genes")), parse_regions(body.get("new_regions")))
    chosen = [g for g in _strings(body.get("genes")) if g in dataset.genes] + known
    new_names = new_genes + [r["name"] for r in new_regions]
    if not chosen and not new_names:
        chosen = list(dataset.genes)
    genes = list(dict.fromkeys(chosen + new_names))
    if not genes:
        raise LaunchError("Choose at least one window, or add genes / regions")
    samples = [s for s in _strings(body.get("samples")) if s in dataset.sample_index]
    haplotypes = [h for h in _strings(body.get("haplotypes")) if h in ("H1", "H2", "ref")] or ["H1", "H2"]
    reference_only = haplotypes == ["ref"]
    command = [_python(), "-m", PREDICT_MODULE, str(dataset.path), "--outputs", ",".join(outputs), *_ontology_args(terms),
               "--genes", ",".join(genes), "--haplotypes", ",".join(haplotypes)]
    files: Dict[str, str] = {}
    steps: List[Dict[str, Any]] = []
    if new_names:
        files["extend_spec.json"] = json.dumps(extension_spec(dataset, new_genes, new_regions), indent=2)
        steps.append({"title": f"Build {len(new_names)} new window(s) for every sample", "command": [_python(), "-m", IMPORT_MODULE, "--spec", "{task_dir}/extend_spec.json"], "weight": 1.0})
    if samples and len(samples) < len(dataset.samples) and not reference_only:
        files["samples.txt"] = "\n".join(samples) + "\n"
        command += ["--samples-file", "{task_dir}/samples.txt"]
    if body.get("overwrite"):
        command.append("--overwrite")
    steps.append({"title": "AlphaGenome predictions", "command": command, "weight": 3.0})
    n_samples = 0 if reference_only else len(samples) or len(dataset.samples)
    sample_haps = [h for h in haplotypes if h != "ref"]
    params = {
        "dataset_id": dataset.id,
        "dataset_path": str(dataset.path),
        "outputs": outputs,
        "ontology_terms": terms or "all",
        "genes": genes,
        "new_windows": new_names,
        "samples": n_samples,
        "haplotypes": haplotypes,
        "windows": n_samples * len(genes) * len(sample_haps) + (len(genes) if "ref" in haplotypes else 0),
        "reference": "ref" in haplotypes,
        "overwrite": bool(body.get("overwrite")),
    }
    what = ", ".join(o.lower() for o in outputs) if len(outputs) < 5 else f"{len(outputs)} outputs"
    who = "reference genome" if reference_only else f"{n_samples} samples{' + reference' if 'ref' in haplotypes else ''}"
    title = f"AlphaGenome {what} for {dataset.name} ({who} × {len(genes)} windows{f', {len(new_names)} new' if new_names else ''})"
    return title, steps, params, files


# ------------------------------------------------------------------------------------- training
def base_configs() -> List[Dict[str, Any]]:
    """Training config templates under configs/predictors/genotype_based (path, model, target)."""
    root = _repo_root() / "configs" / "predictors" / "genotype_based"
    items = []
    for path in sorted(root.rglob("*.yaml")):
        rel = path.relative_to(_repo_root()).as_posix()
        items.append({"path": rel, "name": path.relative_to(root).as_posix()})
    return items


def _load_base(body: Dict[str, Any], runs: Dict[str, Path]) -> Tuple[Dict[str, Any], str]:
    run_id = body.get("base_run")
    if run_id:
        path = runs.get(run_id)
        if path is None or not (path / "config.yaml").exists():
            raise LaunchError(f"Unknown run: {run_id}")
        return yaml.safe_load((path / "config.yaml").read_text(encoding="utf-8")) or {}, f"run {path.name}"
    rel = str(body.get("base_config") or DEFAULT_BASE_CONFIG)
    path = (_repo_root() / rel).resolve()
    root = (_repo_root() / "configs").resolve()
    if root not in path.parents or not path.exists():
        raise LaunchError(f"Config not found under configs/: {rel}")
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}, rel


def train_options(dataset: Dataset) -> Dict[str, Any]:
    outputs: Dict[str, Any] = {}
    for gene in dataset.genes:
        info = dataset.gene_info(gene)
        for name, out in (info.get("outputs") or {}).items():
            if name in outputs or (out.get("resolution") or 1) != 1:
                continue
            terms: Dict[str, Dict[str, Any]] = {}
            for track in out.get("tracks") or []:
                meta = track.get("metadata") or {}
                curie = meta.get("ontology_curie")
                if not curie:
                    continue
                entry = terms.setdefault(curie, {"curie": curie, "name": meta.get("biosample_name") or "", "strands": []})
                if meta.get("strand") not in entry["strands"]:
                    entry["strands"].append(meta.get("strand"))
            outputs[name] = {"tracks": len(out.get("tracks") or []), "ontologies": list(terms.values())}
        if len(outputs) >= 8:
            break
    return {
        "base_configs": base_configs(),
        "default_base": DEFAULT_BASE_CONFIG,
        "model_types": list(MODEL_TYPES),
        "outputs": outputs,
        "genes": list(dataset.genes),
        "fields": [f for f in dataset.fields if f["kind"] == "categorical"],
    }


def _rows_per_gene(n_ontologies: int, strands: int, feature_mode: str, di: Dict[str, Any]) -> int:
    masks = 2 + (1 if di.get("indel_include_valid_mask", True) else 0) + (1 if di.get("indel_include_snp_mask") else 0)
    signals = n_ontologies * strands
    if feature_mode == "masks_only":
        return masks
    if feature_mode == "signals_only":
        return signals
    return signals + masks


def build_train_config(dataset: Dataset, body: Dict[str, Any], runs: Dict[str, Path], runs_root: Path, cache_root: Path) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Training config for ``dataset`` from a template plus the form's choices."""
    base, base_label = _load_base(body, runs)
    config = copy.deepcopy(base)
    di = config.setdefault("dataset_input", {})
    for key in ("dataset_id", "view_path", "sample_ids_path", "superpopulations_to_use", "populations_to_use", "gene_window_metadata", "gene_order", "gene_track_strands", "normalization_params_path", "reference_predictions_dataset_dir", "reference_predictions_sample_id"):
        di.pop(key, None)
    di["dataset_dir"] = str(dataset.path)
    di["consensus_dataset_dir"] = str(dataset.path)
    genes = [g for g in _strings(body.get("genes")) if g in dataset.genes]
    if not genes:
        raise LaunchError("Choose at least one gene window")
    di["genes_to_use"] = genes
    output = str(body.get("output") or (di.get("alphagenome_outputs") or ["rna_seq"])[0]).lower()
    info_outputs = {}
    for gene in genes:
        info_outputs.update(dataset.gene_info(gene).get("outputs") or {})
    if output not in info_outputs:
        raise LaunchError(f"No {output} predictions in this dataset for the chosen genes (run AlphaGenome predictions first)")
    if (info_outputs[output].get("resolution") or 1) != 1:
        raise LaunchError(f"{output} is stored at {info_outputs[output]['resolution']} bp resolution; training uses per-base outputs")
    di["alphagenome_outputs"] = [output]
    tracks = [t.get("metadata") or {} for t in info_outputs[output].get("tracks") or []]
    available = {m.get("ontology_curie") for m in tracks if m.get("ontology_curie")}
    terms = [t for t in _strings(body.get("ontology_terms")) if t in available] or sorted(available)
    if not terms:
        raise LaunchError(f"The {output} tracks have no ontology metadata")
    di["ontology_terms"] = terms
    strands = sorted({str(m.get("strand")) for m in tracks if m.get("ontology_curie") in terms})
    if strands == ["."]:
        di["track_strands"] = ["."]
        n_strands = 1
    else:
        di.pop("track_strands", None)
        n_strands = 2
    wcs = int(body.get("window_center_size") or di.get("window_center_size") or 32768)
    lengths = [dataset.window(g).length or 0 for g in genes]
    if wcs <= 0 or (lengths and wcs > min(lengths)):
        raise LaunchError(f"Window center size must be between 1 and {min(lengths)} bp")
    di["window_center_size"] = wcs
    feature_mode = str(body.get("feature_mode") or di.get("feature_mode") or "signals_only")
    if feature_mode not in ("signals_only", "signals_and_masks", "masks_only"):
        raise LaunchError("feature_mode must be signals_only, signals_and_masks or masks_only")
    di["feature_mode"] = feature_mode
    cohort = [s for s in _strings(body.get("sample_ids")) if s in dataset.sample_index]
    di["sample_ids"] = cohort if cohort and len(cohort) < len(dataset.samples) else None
    di["processed_cache_dir"] = str(cache_root / slugify(dataset.id))
    group = slugify(str(body.get("run_name") or "")) or f"visualizer-{time.strftime('%Y%m%d-%H%M%S')}"
    di["results_dir"] = str(runs_root / group)
    di["cache_processed_tensors"] = True

    out = config.setdefault("output", {})
    field = str(body.get("target_field") or "")
    if field not in {f["name"] for f in dataset.fields}:
        raise LaunchError("Choose the sample field to predict")
    class_map: Dict[str, List[str]] = {}
    for value, label in (body.get("class_map") or {}).items():
        label = str(label or "").strip()
        if label:
            class_map.setdefault(label, []).append(str(value))
    if len(class_map) < 2:
        raise LaunchError("Choose at least two classes")
    identity = all(len(v) == 1 and v[0] == k for k, v in class_map.items())
    if identity:
        out["prediction_target"] = field
        out["derived_targets"] = {}
    else:
        target = slugify(str(body.get("target_name") or f"{field}_classes")).replace("-", "_").replace(".", "_")
        out["prediction_target"] = target
        out["derived_targets"] = {target: {"source_field": field, "exclude_unmapped": True, "class_map": class_map}}
    out["known_classes"] = sorted(class_map)

    model = config.setdefault("model", {})
    model_type = str(body.get("model_type") or model.get("type") or "CNN2").upper()
    if model_type not in MODEL_TYPES:
        raise LaunchError(f"Model type must be one of {', '.join(MODEL_TYPES)}")
    model["type"] = model_type
    rows = _rows_per_gene(len(terms), n_strands, feature_mode, di)
    cnn2 = model.setdefault("cnn2", {})
    for key in ("kernel_stage1", "stride_stage1"):
        width = (cnn2.get(key) or [rows, 32])[1]
        cnn2[key] = [rows, width]
    cnn = model.setdefault("cnn", {})
    if isinstance(cnn.get("kernel_size"), list):
        cnn["kernel_size"] = [rows, cnn["kernel_size"][1]]
    if isinstance(cnn.get("stride"), list):
        cnn["stride"] = [rows, cnn["stride"][1]]

    training = config.setdefault("training", {})
    for key, cast in (("num_epochs", int), ("batch_size", int), ("learning_rate", float)):
        if body.get(key) not in (None, ""):
            training[key] = cast(body[key])
    split = config.setdefault("data_split", {})
    fractions = [float(body.get(f"{s}_split") or split.get(f"{s}_split") or 0) for s in SPLITS]
    if any(f <= 0 for f in fractions) or abs(sum(fractions) - 1.0) > 1e-6:
        raise LaunchError("Train/val/test fractions must be positive and add up to 1")
    for name, value in zip(SPLITS, fractions):
        split[f"{name}_split"] = value
    if body.get("seed") not in (None, ""):
        split["random_seed"] = int(body["seed"])
    has_families = any(row.get("family_id") not in (None, "", row["sample_id"]) for row in dataset.samples)
    split.setdefault("family_split_mode", "family_aware")
    if not has_families:
        split["family_split_mode"] = "ignore"
    config.setdefault("wandb", {})["use_wandb"] = False
    config.setdefault("checkpointing", {})["load_checkpoint"] = None
    config["mode"] = "train"
    config.setdefault("metadata", {})["note"] = f"Started from the visualizer ({base_label}) on {dataset.name}"
    summary = {
        "dataset_id": dataset.id,
        "dataset_path": str(dataset.path),
        "base": base_label,
        "target": out["prediction_target"],
        "classes": out["known_classes"],
        "genes": genes,
        "output": output,
        "ontology_terms": terms,
        "model": model_type,
        "rows_per_gene": rows,
        "input_shape": [2 * rows * len(genes), wcs],
        "results_dir": di["results_dir"],
        "samples": len(cohort) if di["sample_ids"] else len(dataset.samples),
        "epochs": training.get("num_epochs"),
    }
    return config, summary


def validate_config(config: Dict[str, Any], scratch: Path) -> None:
    """Typed validation with the genotype config schema (no torch needed)."""
    from genomics.predictors.genotype_based.config import load_config

    scratch.mkdir(parents=True, exist_ok=True)
    path = scratch / "validate_config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    try:
        load_config(path)
    except Exception as exc:
        raise LaunchError(f"Invalid training config: {exc}")
    finally:
        path.unlink(missing_ok=True)


def train_task(config: Dict[str, Any], summary: Dict[str, Any], evaluate_test: bool) -> Tuple[str, List[Dict[str, Any]], Dict[str, Any], Dict[str, str]]:
    epochs = int(summary.get("epochs") or (config.get("training") or {}).get("num_epochs") or 0)
    steps = [{
        "title": "Training",
        "command": [_python(), "-m", TRAIN_MODULE, "{task_dir}/config.yaml"],
        "weight": 9.0,
        "history_glob": f"{summary['results_dir']}/*/models/training_history.json",
        "epochs": epochs,
    }]
    if evaluate_test:
        steps.append({
            "title": "Evaluation on the test split",
            "command": [_python(), "-m", EVALUATE_MODULE, "{task_dir}/config.yaml", "--checkpoint", "best_accuracy", "--split", "test", "--output-name", "test_best_accuracy"],
            "weight": 1.0,
        })
    title = f"Train {summary['model']} on {summary['target']} ({len(summary['genes'])} genes, {summary['samples']} samples)"
    return title, steps, {**summary, "evaluate_test": evaluate_test}, {"config.yaml": yaml.safe_dump(config, sort_keys=False, allow_unicode=True)}


def evaluate_task(run_id: str, run_dir: Path, body: Dict[str, Any]) -> Tuple[str, List[Dict[str, Any]], Dict[str, Any], Dict[str, str]]:
    split = str(body.get("split") or "test")
    if split not in SPLITS:
        raise LaunchError("split must be train, val or test")
    checkpoints = sorted(p.name for p in (run_dir / "models").glob("*.pt")) if (run_dir / "models").is_dir() else []
    is_sklearn = (run_dir / "models" / "sklearn_baseline.joblib").exists()
    checkpoint = str(body.get("checkpoint") or ("best_accuracy.pt" if "best_accuracy.pt" in checkpoints else (checkpoints[0] if checkpoints else "best_accuracy")))
    if checkpoints and checkpoint not in checkpoints and not is_sklearn:
        raise LaunchError(f"Checkpoint {checkpoint} not in this run ({', '.join(checkpoints)})")
    if not (run_dir / "config.yaml").exists():
        raise LaunchError("This run has no config.yaml")
    stem = Path(checkpoint).stem
    output_name = f"{split}_{stem}"
    command = [_python(), "-m", EVALUATE_MODULE, str(run_dir / "config.yaml"), "--checkpoint", stem, "--split", split, "--experiment-dir", str(run_dir), "--output-name", output_name]
    params = {"run": run_id, "run_dir": str(run_dir), "checkpoint": checkpoint, "split": split, "results_file": f"{output_name}_results.json"}
    return f"Evaluate {run_dir.name} ({stem}) on {split}", [{"title": f"Evaluation on {split}", "command": command}], params, {}
