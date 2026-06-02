"""Helpers to replace duplicated datasets with lightweight subset manifests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

from rich.console import Console

console = Console()


MANIFEST_VERSION = 1


def _load_json(path: Path) -> Dict:
    with open(path) as f:
        return json.load(f)


def _write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def create_subset_manifest(dataset_dir: Path, canonical_dataset_dir: Path, output_dir: Path | None = None) -> Path:
    """Creates a lightweight manifest preserving a subset dataset definition."""
    dataset_dir = Path(dataset_dir)
    canonical_dataset_dir = Path(canonical_dataset_dir)
    output_dir = Path(output_dir) if output_dir else dataset_dir.parent / f"{dataset_dir.name}_manifest"

    meta = _load_json(dataset_dir / "dataset_metadata.json")
    canonical_meta = _load_json(canonical_dataset_dir / "dataset_metadata.json")

    individuals: List[str] = list(meta.get("individuals", []))
    canonical_individuals = set(canonical_meta.get("individuals", []))
    missing = sorted(set(individuals) - canonical_individuals)
    if missing:
        raise ValueError(
            f"{dataset_dir} não é subconjunto de {canonical_dataset_dir}. Exemplos ausentes: {missing[:5]}"
        )

    manifest = {
        "manifest_type": "subset_dataset",
        "manifest_version": MANIFEST_VERSION,
        "source_dataset_dir": str(dataset_dir.resolve()),
        "canonical_dataset_dir": str(canonical_dataset_dir.resolve()),
        "dataset_name": meta.get("dataset_name", dataset_dir.name),
        "total_individuals": len(individuals),
        "individuals": individuals,
        "window_size": meta.get("window_size"),
        "genes": meta.get("genes", []),
        "alphagenome_outputs": meta.get("alphagenome_outputs", []),
        "population_distribution": meta.get("population_distribution", {}),
        "superpopulation_distribution": meta.get("superpopulation_distribution", {}),
        "individuals_pedigree": meta.get("individuals_pedigree", {}),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_json(output_dir / "subset_manifest.json", manifest)
    _write_json(output_dir / "dataset_metadata.json", {
        "dataset_layout": "subset_manifest",
        "manifest_version": MANIFEST_VERSION,
        "canonical_dataset_dir": str(canonical_dataset_dir.resolve()),
        "dataset_name": manifest["dataset_name"],
        "total_individuals": manifest["total_individuals"],
        "individuals": individuals,
        "window_size": manifest["window_size"],
        "genes": manifest["genes"],
        "alphagenome_outputs": manifest["alphagenome_outputs"],
        "population_distribution": manifest["population_distribution"],
        "superpopulation_distribution": manifest["superpopulation_distribution"],
        "individuals_pedigree": manifest["individuals_pedigree"],
    })
    console.print(f"[green]Manifesto criado:[/green] {output_dir}")
    return output_dir
