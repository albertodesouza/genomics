from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


def _load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"YAML inválido: {path}")
    return payload


def _write_yaml(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _append_jsonl(path: Path, row: Dict[str, Any]) -> None:
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(row, ensure_ascii=True) + "\n")


def _sanitize(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value)


def _result_key(gene: str, ontology: str) -> str:
    return f"{gene}::{ontology}"


def _discover_genes(dataset_dir: Path) -> List[str]:
    metadata = _load_json(dataset_dir / "dataset_metadata.json")
    catalog = metadata.get("window_catalog")
    if isinstance(catalog, dict) and catalog:
        return sorted(catalog.keys())
    individuals = metadata.get("individuals", [])
    if not individuals:
        raise ValueError(f"Nenhum indivíduo em {dataset_dir / 'dataset_metadata.json'}")
    sample_meta = _load_json(dataset_dir / "individuals" / str(individuals[0]) / "individual_metadata.json")
    windows = sample_meta.get("windows", [])
    if isinstance(windows, dict):
        return sorted(windows.keys())
    return sorted(str(item) for item in windows)


def _discover_ontologies(dataset_dir: Path) -> List[str]:
    metadata = _load_json(dataset_dir / "dataset_metadata.json")
    ontology_details = metadata.get("ontology_details", {})
    if isinstance(ontology_details, dict) and ontology_details:
        return sorted(ontology_details.keys())
    ontologies = metadata.get("ontologies", [])
    cleaned = sorted({str(item).strip() for item in ontologies if str(item).strip()})
    if cleaned:
        return cleaned
    raise ValueError(f"Nenhuma ontologia encontrada em {dataset_dir}")


def _load_gene_list(path: Path) -> List[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.strip().startswith("#")]


def _read_metrics(experiment_dir: Path) -> Dict[str, Any]:
    history_path = experiment_dir / "models" / "training_history.json"
    test_path = experiment_dir / "test_results.json"
    val_path = experiment_dir / "val_results.json"
    train_path = experiment_dir / "train_results.json"
    history = _load_json(history_path) if history_path.exists() else {}
    test = _load_json(test_path) if test_path.exists() else {}
    val = _load_json(val_path) if val_path.exists() else {}
    train = _load_json(train_path) if train_path.exists() else {}
    return {
        "train_accuracy": train.get("accuracy"),
        "val_accuracy": val.get("accuracy"),
        "test_accuracy": test.get("accuracy"),
        "train_f1": train.get("f1"),
        "val_f1": val.get("f1"),
        "test_f1": test.get("f1"),
        "best_val_accuracy": max(history.get("val_accuracy", []) or [None]),
        "epochs_recorded": len(history.get("epoch", [])),
    }


def _load_existing(path: Path) -> Dict[str, Dict[str, Any]]:
    if not path.exists():
        return {}
    rows: Dict[str, Dict[str, Any]] = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            rows[_result_key(row["gene"], row.get("ontology") or "all_ontologies")] = row
    return rows


def _write_results(path: Path, rows: List[Dict[str, Any]]) -> None:
    fields = [
        "gene", "ontology", "status", "returncode", "experiment_dir", "config_path", "log_path",
        "train_accuracy", "val_accuracy", "test_accuracy", "train_f1", "val_f1", "test_f1",
        "best_val_accuracy", "epochs_recorded", "started_at", "finished_at",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _build_config(base: Dict[str, Any], gene: str, ontology: Optional[str], output_root: Path, prefix: str) -> Dict[str, Any]:
    config = deepcopy(base)
    dataset_input = config.setdefault("dataset_input", {})
    gene_slug = _sanitize(gene)
    ontology_slug = _sanitize(ontology or "all_ontologies")
    dataset_input["genes_to_use"] = [gene]
    dataset_input["processed_cache_dir"] = str(output_root / "cache" / ontology_slug / gene_slug)
    dataset_input["results_dir"] = str(output_root / "runs")
    if ontology:
        dataset_input["ontology_terms"] = [ontology]
    else:
        dataset_input.pop("ontology_terms", None)

    model = config.setdefault("model", {})
    tracks_per_gene = 2 if ontology else None
    if tracks_per_gene is not None:
        cnn = model.get("cnn") or {}
        if isinstance(cnn.get("kernel_size"), list) and cnn["kernel_size"][0] > tracks_per_gene:
            cnn["kernel_size"][0] = tracks_per_gene
        if isinstance(cnn.get("stride"), list) and cnn["stride"][0] > tracks_per_gene:
            cnn["stride"][0] = tracks_per_gene
        cnn2 = model.setdefault("cnn2", {})
        kernel = list(cnn2.get("kernel_stage1", [tracks_per_gene, 32]))
        stride = list(cnn2.get("stride_stage1", [tracks_per_gene, 32]))
        if kernel[0] > tracks_per_gene:
            kernel[0] = tracks_per_gene
        if stride[0] > tracks_per_gene:
            stride[0] = tracks_per_gene
        cnn2["kernel_stage1"] = kernel
        cnn2["stride_stage1"] = stride

    wandb = config.setdefault("wandb", {})
    wandb["run_name"] = f"{prefix}_{gene_slug}_{ontology_slug}"
    return config


def _run_train(config_path: Path, log_path: Path) -> tuple[int, Optional[Path]]:
    command = [sys.executable, "-m", "genomics.predictors.genotype_based.experiments.train", str(config_path)]
    experiment_dir: Optional[Path] = None
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log.write(line)
            if "Diretório do experimento:" in line:
                experiment_dir = Path(line.split("Diretório do experimento:", 1)[1].strip())
        return proc.wait(), experiment_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Run genotype_based_predictor single-gene screen")
    parser.add_argument("config", type=Path, help="Base YAML config")
    parser.add_argument("--dataset-dir", type=Path, default=None, help="Dataset dir for gene discovery; defaults to config dataset_dir")
    parser.add_argument("--genes-file", type=Path, default=None)
    parser.add_argument("--genes", nargs="*", default=None)
    parser.add_argument("--output-root", type=Path, default=Path("results/genotype_based_predictor/single_gene_screen"))
    parser.add_argument("--run-name-prefix", default="single_gene")
    parser.add_argument("--split-by-ontology", action="store_true")
    parser.add_argument("--ontology", action="append", default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    base_config = _load_yaml(args.config.resolve())
    dataset_input = base_config.get("dataset_input", {}) if isinstance(base_config.get("dataset_input"), dict) else {}
    dataset_dir = args.dataset_dir or Path(dataset_input.get("dataset_dir", "/dados/GENOMICS_DATA/v1/1kG_high_coverage"))
    dataset_dir = dataset_dir.resolve()
    output_root = args.output_root.resolve()
    temp_dir = output_root / "temp_configs"
    logs_dir = output_root / "logs"
    output_root.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)

    if args.genes:
        genes = args.genes
    elif args.genes_file:
        genes = _load_gene_list(args.genes_file.resolve())
    else:
        genes = list(dataset_input.get("genes_to_use") or _discover_genes(dataset_dir))

    available = set(_discover_genes(dataset_dir))
    missing = [gene for gene in genes if gene not in available]
    if missing:
        raise ValueError(f"Genes ausentes no dataset: {', '.join(missing)}")

    if args.split_by_ontology:
        ontologies = args.ontology or _discover_ontologies(dataset_dir)
    else:
        if args.ontology:
            raise ValueError("--ontology requer --split-by-ontology")
        ontologies = [None]

    manifest = []
    for ontology in ontologies:
        for gene in genes:
            cfg = _build_config(base_config, gene, ontology, output_root, args.run_name_prefix)
            ontology_slug = _sanitize(ontology or "all_ontologies")
            config_path = temp_dir / f"{ontology_slug}__{_sanitize(gene)}.yaml"
            _write_yaml(config_path, cfg)
            manifest.append({"gene": gene, "ontology": ontology or "all_ontologies", "config_path": str(config_path)})

    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if args.dry_run:
        print(f"Manifest: {manifest_path}")
        print(f"Configs: {temp_dir}")
        return

    results_csv = output_root / "results.csv"
    results_jsonl = output_root / "results.jsonl"
    rows = _load_existing(results_csv)
    for item in manifest:
        gene = item["gene"]
        ontology = item["ontology"]
        key = _result_key(gene, ontology)
        previous = rows.get(key)
        if previous and previous.get("status") == "completed" and not args.force:
            print(f"[skip] {gene} / {ontology}")
            continue
        started = datetime.utcnow().isoformat()
        config_path = Path(item["config_path"])
        log_path = logs_dir / f"{_sanitize(ontology)}__{_sanitize(gene)}.log"
        print(f"[run] {gene} / {ontology}")
        returncode, experiment_dir = _run_train(config_path, log_path)
        finished = datetime.utcnow().isoformat()
        row: Dict[str, Any] = {
            "gene": gene,
            "ontology": ontology,
            "status": "completed" if returncode == 0 else "failed",
            "returncode": returncode,
            "experiment_dir": str(experiment_dir) if experiment_dir else "",
            "config_path": str(config_path),
            "log_path": str(log_path),
            "train_accuracy": "",
            "val_accuracy": "",
            "test_accuracy": "",
            "train_f1": "",
            "val_f1": "",
            "test_f1": "",
            "best_val_accuracy": "",
            "epochs_recorded": "",
            "started_at": started,
            "finished_at": finished,
        }
        if returncode == 0 and experiment_dir:
            row.update({k: v if v is not None else "" for k, v in _read_metrics(experiment_dir).items()})
        rows[key] = row
        ordered = [rows[key] for key in sorted(rows)]
        _write_results(results_csv, ordered)
        _append_jsonl(results_jsonl, row)
        if returncode != 0:
            print(f"[fail] {gene} / {ontology}: {log_path}")
        else:
            print(f"[done] {gene} / {ontology}: test_acc={row['test_accuracy']}")

    print(f"Manifest: {manifest_path}")
    print(f"Results CSV: {results_csv}")


if __name__ == "__main__":
    main()
