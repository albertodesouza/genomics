#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import yaml
from datasets import load_from_disk


def load_yaml(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_yaml(path: Path, payload: Dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)


def load_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def append_jsonl(path: Path, row: Dict) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, ensure_ascii=True) + "\n")


def discover_genes_from_hf_dataset(dataset_path: Path) -> List[str]:
    gene_windows = load_from_disk(str(dataset_path / "gene_windows"))
    return sorted(set(gene_windows["gene"]))


def sanitize_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value)


def extract_group(gene: str, curated_genes: set[str]) -> str:
    return "curated11" if gene in curated_genes else "random"


def build_temp_config(base_config: Dict, gene: str, run_name_prefix: str) -> Dict:
    config = deepcopy(base_config)
    gene_slug = sanitize_name(gene)
    run_slug = f"{run_name_prefix}_{gene_slug}"
    config["dataset_input"]["genes_to_use"] = [gene]
    config.setdefault("wandb", {})
    config["wandb"]["run_name"] = run_slug
    return config


def run_gene(script_module: str, temp_config_path: Path, output_dir: Path, log_path: Path) -> int:
    command = [sys.executable, "-m", script_module, "--config", str(temp_config_path), "--output-dir", str(output_dir)]
    with log_path.open("w", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(temp_config_path.parent.parent),
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_file.write(line)
            log_file.flush()
        return process.wait()


def read_metrics(experiment_dir: Path) -> Dict:
    results_path = experiment_dir / "results.json"
    if not results_path.exists():
        return {}
    results = load_json(results_path)
    return {
        "train_accuracy": _last_value(results.get("train_accuracy")),
        "val_accuracy": _last_value(results.get("val_accuracy")),
        "test_accuracy": results.get("test_accuracy"),
        "train_f1": _last_value(results.get("train_f1")),
        "val_f1": _last_value(results.get("val_f1")),
        "test_f1": results.get("test_f1"),
        "test_macro_f1": results.get("test_macro_f1"),
        "best_val_accuracy": _max_value(results.get("val_accuracy")),
        "epochs_recorded": len(results.get("train_loss", [])),
    }


def _last_value(values):
    if isinstance(values, list) and values:
        return values[-1]
    return None


def _max_value(values):
    if isinstance(values, list) and values:
        return max(values)
    return None


def load_existing_results(results_csv: Path) -> Dict[str, Dict]:
    if not results_csv.exists():
        return {}
    rows: Dict[str, Dict] = {}
    with results_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows[row["gene"]] = row
    return rows


def write_results_csv(results_csv: Path, rows: List[Dict]) -> None:
    fieldnames = [
        "gene",
        "group",
        "status",
        "returncode",
        "experiment_dir",
        "log_path",
        "train_accuracy",
        "val_accuracy",
        "test_accuracy",
        "train_f1",
        "val_f1",
        "test_f1",
        "test_macro_f1",
        "best_val_accuracy",
        "epochs_recorded",
        "started_at",
        "finished_at",
    ]
    with results_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def row_needs_rerun(row: Optional[Dict]) -> bool:
    if not row:
        return True
    if row.get("status") != "completed":
        return True
    if row.get("test_accuracy") in ("", None):
        return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Run resumable single-gene screening with neural_prediction")
    parser.add_argument("--base-config", required=True, help="Base YAML config for neural_prediction")
    parser.add_argument("--output-root", required=True, help="Directory for manifests, temp configs, logs, and runs")
    parser.add_argument("--run-name-prefix", default="neural_prediction_single", help="Prefix for run names")
    parser.add_argument("--force", action="store_true", help="Rerun genes even if they already completed")
    parser.add_argument("--genes-file", help="Optional newline-delimited gene list")
    args = parser.parse_args()

    base_config_path = Path(args.base_config).resolve()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifests_dir = output_root / "manifests"
    temp_configs_dir = output_root / "temp_configs"
    logs_dir = output_root / "logs"
    runs_dir = output_root / "runs"
    manifests_dir.mkdir(exist_ok=True)
    temp_configs_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    runs_dir.mkdir(exist_ok=True)

    base_config = load_yaml(base_config_path)
    dataset_path = Path(base_config["dataset_input"]["hf_dataset_path"]).resolve()
    all_genes = discover_genes_from_hf_dataset(dataset_path)

    if args.genes_file:
        requested_genes = [line.strip() for line in Path(args.genes_file).read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    else:
        requested_genes = all_genes

    missing = sorted(set(requested_genes) - set(all_genes))
    if missing:
        raise ValueError(f"Genes not found in HF dataset: {', '.join(missing)}")

    curated_genes = set(base_config.get("dataset_input", {}).get("genes_to_use", []))
    manifest_rows = [{"gene": gene, "group": extract_group(gene, curated_genes)} for gene in requested_genes]
    with (manifests_dir / "gene_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest_rows, handle, indent=2)
        handle.write("\n")

    results_csv = output_root / "results.csv"
    results_jsonl = output_root / "results.jsonl"
    existing = load_existing_results(results_csv)
    all_rows = list(existing.values())
    row_by_gene = {row["gene"]: row for row in all_rows}

    for entry in manifest_rows:
        gene = entry["gene"]
        current = row_by_gene.get(gene)
        if not args.force and not row_needs_rerun(current):
            continue

        temp_config = build_temp_config(base_config, gene, args.run_name_prefix)
        temp_config_path = temp_configs_dir / f"{sanitize_name(gene)}.yaml"
        write_yaml(temp_config_path, temp_config)

        experiment_dir = runs_dir / sanitize_name(gene)
        log_path = logs_dir / f"{sanitize_name(gene)}.log"
        started_at = datetime.now().isoformat()
        returncode = run_gene("neural_prediction.main", temp_config_path, experiment_dir, log_path)
        finished_at = datetime.now().isoformat()
        metrics = read_metrics(experiment_dir)
        row = {
            "gene": gene,
            "group": entry["group"],
            "status": "completed" if returncode == 0 else "failed",
            "returncode": returncode,
            "experiment_dir": str(experiment_dir),
            "log_path": str(log_path),
            "train_accuracy": metrics.get("train_accuracy", ""),
            "val_accuracy": metrics.get("val_accuracy", ""),
            "test_accuracy": metrics.get("test_accuracy", ""),
            "train_f1": metrics.get("train_f1", ""),
            "val_f1": metrics.get("val_f1", ""),
            "test_f1": metrics.get("test_f1", ""),
            "test_macro_f1": metrics.get("test_macro_f1", ""),
            "best_val_accuracy": metrics.get("best_val_accuracy", ""),
            "epochs_recorded": metrics.get("epochs_recorded", ""),
            "started_at": started_at,
            "finished_at": finished_at,
        }
        row_by_gene[gene] = row
        all_rows = [row_by_gene[key] for key in sorted(row_by_gene)]
        write_results_csv(results_csv, all_rows)
        append_jsonl(results_jsonl, row)


if __name__ == "__main__":
    main()
