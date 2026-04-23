from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from datasets import load_from_disk

from .fileio import load_json


def validate_hf_dataset_store(dataset_dir: Path) -> Dict:
    manifest_path = dataset_dir / "manifest.json"
    samples_path = dataset_dir / "samples"
    gene_windows_path = dataset_dir / "gene_windows"

    errors: List[str] = []
    warnings: List[str] = []

    if not manifest_path.exists():
        errors.append(f"Missing manifest: {manifest_path}")
    if not samples_path.exists():
        errors.append(f"Missing samples dataset: {samples_path}")
    if not gene_windows_path.exists():
        errors.append(f"Missing gene_windows dataset: {gene_windows_path}")

    if errors:
        return {
            "valid": False,
            "errors": errors,
            "warnings": warnings,
            "summary": {},
        }

    manifest = load_json(manifest_path)
    samples = load_from_disk(str(samples_path))
    gene_windows = load_from_disk(str(gene_windows_path))

    sample_ids = set(samples["sample_id"])
    gene_sample_ids = set(gene_windows["sample_id"])
    missing_sample_refs = sorted(gene_sample_ids - sample_ids)
    if missing_sample_refs:
        errors.append(
            f"Found {len(missing_sample_refs)} gene_windows rows referencing unknown sample_ids"
        )

    manifest_sample_count = manifest.get("summary", {}).get("total_samples")
    manifest_gene_window_count = manifest.get("summary", {}).get("total_gene_windows")
    if manifest_sample_count != len(samples):
        errors.append(
            f"Manifest sample count mismatch: manifest={manifest_sample_count} actual={len(samples)}"
        )
    if manifest_gene_window_count != len(gene_windows):
        errors.append(
            "Manifest gene_window count mismatch: "
            f"manifest={manifest_gene_window_count} actual={len(gene_windows)}"
        )

    genes_in_rows = sorted(set(gene_windows["gene"]))
    genes_in_manifest = manifest.get("summary", {}).get("genes", [])
    if genes_in_rows != genes_in_manifest:
        errors.append("Manifest gene list mismatch")

    if len(sample_ids) != len(samples):
        errors.append("Duplicate sample_id rows found in samples dataset")

    gene_keys = list(zip(gene_windows["sample_id"], gene_windows["gene"]))
    if len(set(gene_keys)) != len(gene_keys):
        errors.append("Duplicate (sample_id, gene) rows found in gene_windows dataset")

    required_sample_columns = ["source_dataset", "source_path", "source_sample_path"]
    for column in required_sample_columns:
        if column not in samples.column_names:
            errors.append(f"Missing required samples provenance column: {column}")

    required_gene_columns = [
        "source_dataset",
        "source_path",
        "source_sample_path",
        "source_window_path",
        "source_metadata_path",
        "ref_sequence",
        "h1_sequence",
        "h2_sequence",
    ]
    for column in required_gene_columns:
        if column not in gene_windows.column_names:
            errors.append(f"Missing required gene_windows column: {column}")

    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "samples": len(samples),
            "gene_windows": len(gene_windows),
            "genes": len(genes_in_rows),
        },
    }


def inspect_hf_dataset_store(dataset_dir: Path) -> Dict:
    manifest = load_json(dataset_dir / "manifest.json")
    samples = load_from_disk(str(dataset_dir / "samples"))
    gene_windows = load_from_disk(str(dataset_dir / "gene_windows"))

    genes = sorted(set(gene_windows["gene"]))
    source_datasets = sorted(set(samples["source_dataset"]))
    populations = sorted({value for value in samples["population"] if value})
    superpopulations = sorted({value for value in samples["superpopulation"] if value})
    rows_with_sequences = sum(1 for value in gene_windows["ref_sequence"] if value)

    gene_counts: Dict[str, int] = {}
    for gene in gene_windows["gene"]:
        gene_counts[gene] = gene_counts.get(gene, 0) + 1

    top_genes = sorted(gene_counts.items(), key=lambda item: (-item[1], item[0]))[:20]

    return {
        "manifest": manifest,
        "summary": {
            "samples": len(samples),
            "gene_windows": len(gene_windows),
            "genes": len(genes),
            "source_datasets": source_datasets,
            "populations": populations,
            "superpopulations": superpopulations,
            "rows_with_sequences": rows_with_sequences,
            "top_genes": top_genes,
        },
    }
