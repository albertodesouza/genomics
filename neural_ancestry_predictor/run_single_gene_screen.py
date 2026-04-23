#!/usr/bin/env python3

import argparse
import csv
import json
import subprocess
import sys
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml


def load_yaml(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_yaml(path: Path, data: Dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def load_json(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def append_jsonl(path: Path, row: Dict) -> None:
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(row, ensure_ascii=True) + "\n")


def discover_genes_from_dataset(dataset_dir: Path) -> List[str]:
    metadata_path = dataset_dir / "dataset_metadata.json"
    metadata = load_json(metadata_path)
    individuals = metadata.get("individuals", [])
    if not individuals:
        raise ValueError(f"No individuals found in {metadata_path}")

    sample_id = individuals[0]
    sample_metadata_path = dataset_dir / "individuals" / sample_id / "individual_metadata.json"
    sample_metadata = load_json(sample_metadata_path)
    genes = sample_metadata.get("windows", [])
    if not genes:
        raise ValueError(f"No windows found in {sample_metadata_path}")
    return genes


def discover_ontologies_from_dataset(dataset_dir: Path) -> List[str]:
    metadata_path = dataset_dir / "dataset_metadata.json"
    metadata = load_json(metadata_path)
    ontology_details = metadata.get("ontology_details", {})
    if ontology_details:
        return sorted(ontology_details.keys())

    ontologies = metadata.get("ontologies", [])
    cleaned = sorted({ontology.strip() for ontology in ontologies if ontology and ontology.strip()})
    if cleaned:
        return cleaned

    raise ValueError(f"No ontologies found in {metadata_path}")


def build_gene_dataset_map(dataset_dirs: List[Path]) -> Dict[str, Path]:
    gene_to_dataset: Dict[str, Path] = {}
    for dataset_dir in dataset_dirs:
        genes = discover_genes_from_dataset(dataset_dir)
        for gene in genes:
            if gene in gene_to_dataset:
                raise ValueError(
                    f"Gene {gene} appears in multiple datasets: "
                    f"{gene_to_dataset[gene]} and {dataset_dir}"
                )
            gene_to_dataset[gene] = dataset_dir
    return gene_to_dataset


def build_shared_ontology_list(dataset_dirs: List[Path]) -> List[str]:
    ontology_sets = [set(discover_ontologies_from_dataset(dataset_dir)) for dataset_dir in dataset_dirs]
    if not ontology_sets:
        raise ValueError("No datasets provided")

    shared = set.intersection(*ontology_sets)
    if not shared:
        raise ValueError("Provided datasets do not share any ontology CURIEs")

    return sorted(shared)


def load_gene_list(path: Path) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.append(gene)
    return genes


def extract_macro_f1(classification_report: str) -> Optional[float]:
    for line in classification_report.splitlines():
        stripped = line.strip()
        if stripped.startswith("macro avg"):
            parts = stripped.split()
            if len(parts) >= 5:
                try:
                    return float(parts[4])
                except ValueError:
                    return None
    return None


def read_metrics(experiment_dir: Path) -> Dict:
    training_history_path = experiment_dir / "models" / "training_history.json"
    test_results_path = experiment_dir / "test_results.json"
    if not test_results_path.exists():
        test_results_path = experiment_dir / "test_final_results.json"
    val_results_path = experiment_dir / "val_results.json"
    train_results_path = experiment_dir / "train_results.json"

    history = load_json(training_history_path) if training_history_path.exists() else {}
    test_results = load_json(test_results_path) if test_results_path.exists() else {}
    val_results = load_json(val_results_path) if val_results_path.exists() else {}
    train_results = load_json(train_results_path) if train_results_path.exists() else {}

    best_val_accuracy = None
    if history.get("val_accuracy"):
        best_val_accuracy = max(history["val_accuracy"])

    final_val_accuracy = None
    if history.get("val_accuracy"):
        final_val_accuracy = history["val_accuracy"][-1]

    final_train_accuracy = None
    if history.get("train_accuracy"):
        final_train_accuracy = history["train_accuracy"][-1]

    return {
        "train_accuracy": train_results.get("accuracy", final_train_accuracy),
        "val_accuracy": val_results.get("accuracy", final_val_accuracy),
        "test_accuracy": test_results.get("accuracy"),
        "train_f1": train_results.get("f1"),
        "val_f1": val_results.get("f1"),
        "test_f1": test_results.get("f1"),
        "test_macro_f1": extract_macro_f1(test_results.get("classification_report", "")),
        "best_val_accuracy": best_val_accuracy,
        "epochs_recorded": len(history.get("epoch", [])),
    }


def recover_metrics_from_existing_run(row: Dict) -> Optional[Dict]:
    experiment_dir_value = row.get("experiment_dir", "")
    if not experiment_dir_value:
        return None

    experiment_dir = Path(experiment_dir_value)
    if not experiment_dir.exists():
        return None

    test_results_path = experiment_dir / "test_results.json"
    if not test_results_path.exists():
        test_results_path = experiment_dir / "test_final_results.json"
    if not test_results_path.exists():
        return None

    try:
        return read_metrics(experiment_dir)
    except Exception:
        return None


def row_needs_rerun(row: Optional[Dict]) -> bool:
    if not row:
        return True
    if row.get("status") != "completed":
        return True
    if row.get("test_accuracy") in ("", None):
        return True
    return False


def sanitize_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value)


def result_key(gene: str, ontology: str) -> str:
    return f"{gene}::{ontology}"


def build_temp_config(
    base_config: Dict,
    dataset_dir: Path,
    gene: str,
    ontology: Optional[str],
    tracks_per_gene: Optional[int],
    output_root: Path,
    run_name_prefix: str,
) -> Dict:
    config = deepcopy(base_config)
    gene_slug = sanitize_name(gene)
    ontology_slug = sanitize_name(ontology or "all_ontologies")
    dataset_slug = sanitize_name(dataset_dir.name)
    run_slug = f"{run_name_prefix}_{gene_slug}_{ontology_slug}"

    config["dataset_input"]["dataset_dir"] = str(dataset_dir)
    config["dataset_input"]["genes_to_use"] = [gene]
    if ontology is None:
        config["dataset_input"].pop("ontology_terms", None)
    else:
        config["dataset_input"]["ontology_terms"] = [ontology]
    config["dataset_input"]["processed_cache_dir"] = str(
        output_root / "runs" / ontology_slug / dataset_slug / gene_slug
    )

    config.setdefault("wandb", {})
    config["wandb"]["run_name"] = run_slug

    config.setdefault("checkpointing", {})
    config["checkpointing"]["checkpoint_dir"] = "models"
    # Keep only the final checkpoint so automatic test evaluation still works,
    # but avoid accumulating best/periodic checkpoints across dozens of genes.
    config["checkpointing"]["save_during_training"] = False
    config["checkpointing"]["save_frequency"] = 10_000_000

    if tracks_per_gene is not None:
        cnn_config = config.get("model", {}).get("cnn")
        if cnn_config:
            kernel_size = list(cnn_config.get("kernel_size", []))
            if len(kernel_size) == 2 and kernel_size[0] > tracks_per_gene:
                kernel_size[0] = tracks_per_gene
                cnn_config["kernel_size"] = kernel_size

            stride = cnn_config.get("stride")
            if isinstance(stride, list) and len(stride) == 2 and stride[0] > tracks_per_gene:
                stride[0] = tracks_per_gene
                cnn_config["stride"] = stride

        cnn2_config = config.get("model", {}).setdefault("cnn2", {})
        kernel_stage1 = list(cnn2_config.get("kernel_stage1", [6, 32]))
        if len(kernel_stage1) == 2 and kernel_stage1[0] > tracks_per_gene:
            kernel_stage1[0] = tracks_per_gene
            cnn2_config["kernel_stage1"] = kernel_stage1

        stride_stage1 = list(cnn2_config.get("stride_stage1", [6, 32]))
        if len(stride_stage1) == 2 and stride_stage1[0] > tracks_per_gene:
            stride_stage1[0] = tracks_per_gene
            cnn2_config["stride_stage1"] = stride_stage1

    return config


def run_gene(
    script_path: Path,
    temp_config_path: Path,
    log_path: Path,
) -> Tuple[int, Optional[Path]]:
    command = [sys.executable, str(script_path), "--config", str(temp_config_path)]

    experiment_dir = None
    with open(log_path, "w", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(script_path.parent),
        )

        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_file.write(line)
            log_file.flush()

            marker = "Diret\u00f3rio do experimento:"
            if marker in line:
                experiment_dir = Path(line.split(marker, 1)[1].strip())

        returncode = process.wait()

    return returncode, experiment_dir


def load_existing_results(results_csv: Path) -> Dict[str, Dict]:
    if not results_csv.exists():
        return {}

    existing: Dict[str, Dict] = {}
    with open(results_csv, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ontology = row.get("ontology") or "unknown"
            existing[result_key(row["gene"], ontology)] = row
    return existing


def write_results_csv(results_csv: Path, rows: List[Dict]) -> None:
    fieldnames = [
        "gene",
        "ontology",
        "group",
        "dataset_dir",
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
    with open(results_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run resumable single-gene pigmentation screening")
    parser.add_argument("--base-config", required=True, help="Base pigmentation YAML config")
    parser.add_argument("--dataset", action="append", required=True, help="Dataset directory; repeatable")
    parser.add_argument("--genes-file", help="Optional newline-delimited gene list")
    parser.add_argument("--output-root", required=True, help="Directory for manifests, temp configs, logs, and runs")
    parser.add_argument("--run-name-prefix", default="pigment_single", help="Prefix for W&B run names")
    parser.add_argument(
        "--split-by-ontology",
        action="store_true",
        help="Run each gene once per ontology CURIE instead of using all ontologies together",
    )
    parser.add_argument("--force", action="store_true", help="Rerun genes even if they already completed")
    args = parser.parse_args()

    script_path = Path(__file__).with_name("neural_ancestry_predictor.py")
    base_config_path = Path(args.base_config).resolve()
    output_root = Path(args.output_root).resolve()
    dataset_dirs = [Path(p).resolve() for p in args.dataset]

    output_root.mkdir(parents=True, exist_ok=True)
    manifests_dir = output_root / "manifests"
    temp_configs_dir = output_root / "temp_configs"
    logs_dir = output_root / "logs"
    manifests_dir.mkdir(exist_ok=True)
    temp_configs_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)

    base_config = load_yaml(base_config_path)
    gene_to_dataset = build_gene_dataset_map(dataset_dirs)

    if args.genes_file:
        requested_genes = load_gene_list(Path(args.genes_file).resolve())
    else:
        requested_genes = sorted(gene_to_dataset.keys())

    missing = [gene for gene in requested_genes if gene not in gene_to_dataset]
    if missing:
        raise ValueError(f"Genes not found in provided datasets: {', '.join(missing)}")

    curated_genes = set(base_config.get("dataset_input", {}).get("genes_to_use", []))
    ontologies: List[Optional[str]]
    if args.split_by_ontology:
        ontologies = build_shared_ontology_list(dataset_dirs)
        if not ontologies:
            raise ValueError("Could not determine ontology CURIEs from dataset metadata")
    else:
        ontologies = [None]

    manifest_rows = []
    for ontology in ontologies:
        for gene in requested_genes:
            dataset_dir = gene_to_dataset[gene]
            manifest_rows.append(
                {
                    "gene": gene,
                    "ontology": ontology or "all_ontologies",
                    "dataset_dir": str(dataset_dir),
                    "group": "curated11" if gene in curated_genes else "random",
                }
            )

    manifest_path = manifests_dir / "gene_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest_rows, f, indent=2)

    results_csv = output_root / "results.csv"
    results_jsonl = output_root / "results.jsonl"
    existing = load_existing_results(results_csv)
    all_rows = list(existing.values())
    row_by_key = {
        result_key(row["gene"], row.get("ontology") or "unknown"): row for row in all_rows
    }

    recovered_any = False
    for key, row in list(row_by_key.items()):
        if row.get("status") != "completed":
            continue
        if row.get("test_accuracy") not in ("", None):
            continue

        recovered = recover_metrics_from_existing_run(row)
        if recovered is None:
            continue

        row.update({key: recovered.get(key, "") for key in recovered})
        recovered_any = True

    if recovered_any:
        all_rows = [row_by_key[key] for key in sorted(row_by_key)]
        write_results_csv(results_csv, all_rows)

    for entry in manifest_rows:
        gene = entry["gene"]
        ontology = entry["ontology"]
        dataset_dir = Path(entry["dataset_dir"])
        group = entry["group"]
        key = result_key(gene, ontology)

        if not args.force:
            previous = row_by_key.get(key)
            if previous and not row_needs_rerun(previous):
                print(f"[skip] {gene} / {ontology}: already completed")
                continue

        started_at = datetime.utcnow().isoformat()
        temp_config = build_temp_config(
            base_config,
            dataset_dir,
            gene,
            ontology,
            2 if ontology is not None else None,
            output_root,
            args.run_name_prefix,
        )
        temp_config_path = temp_configs_dir / f"{sanitize_name(ontology)}__{sanitize_name(gene)}.yaml"
        log_path = logs_dir / f"{sanitize_name(ontology)}__{sanitize_name(gene)}.log"
        write_yaml(temp_config_path, temp_config)

        print(f"[run] {gene} / {ontology} ({group}) -> {dataset_dir.name}")
        print(f"      log: {log_path}")
        returncode, experiment_dir = run_gene(script_path, temp_config_path, log_path)
        finished_at = datetime.utcnow().isoformat()

        row = {
            "gene": gene,
            "ontology": ontology,
            "group": group,
            "dataset_dir": str(dataset_dir),
            "status": "failed" if returncode != 0 else "completed",
            "returncode": returncode,
            "experiment_dir": str(experiment_dir) if experiment_dir else "",
            "log_path": str(log_path),
            "train_accuracy": "",
            "val_accuracy": "",
            "test_accuracy": "",
            "train_f1": "",
            "val_f1": "",
            "test_f1": "",
            "test_macro_f1": "",
            "best_val_accuracy": "",
            "epochs_recorded": "",
            "started_at": started_at,
            "finished_at": finished_at,
        }

        if returncode == 0 and experiment_dir:
            metrics = read_metrics(experiment_dir)
            row.update({key: metrics.get(key, "") for key in metrics})

        row_by_key[key] = row
        all_rows = [row_by_key[key] for key in sorted(row_by_key)]
        write_results_csv(results_csv, all_rows)
        append_jsonl(results_jsonl, row)

        if returncode != 0:
            print(f"[fail] {gene} / {ontology}: see {log_path}")
        else:
            print(
                f"[done] {gene} / {ontology}: "
                f"test_acc={row['test_accuracy']} macro_f1={row['test_macro_f1']}"
            )

    print(f"Manifest: {manifest_path}")
    print(f"Results CSV: {results_csv}")
    print(f"Results JSONL: {results_jsonl}")


if __name__ == "__main__":
    main()
