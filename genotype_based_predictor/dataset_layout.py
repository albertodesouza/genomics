"""Utilities to validate and materialize the on-disk dataset layout."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Dict, Iterable

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

console = Console()


LAYOUT_VERSION = 1


def is_dataset_dir(dataset_dir: Path) -> bool:
    marker = Path(dataset_dir) / "layout_metadata.json"
    if not marker.exists():
        return False
    try:
        with open(marker) as f:
            meta = json.load(f)
        return meta.get("layout_version") == LAYOUT_VERSION
    except Exception:
        return False


def _load_json(path: Path) -> Dict:
    with open(path) as f:
        return json.load(f)


def _write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def _copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists():
        shutil.copy2(src, dst)


def _iter_individual_dirs(dataset_dir: Path) -> Iterable[Path]:
    individuals_dir = dataset_dir / "individuals"
    if not individuals_dir.exists():
        return []
    return sorted([p for p in individuals_dir.iterdir() if p.is_dir()])


def materialize_dataset(source_dir: Path, target_dir: Path) -> Path:
    """Creates the normalized dataset view consumed by this package."""
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)

    layout_marker = target_dir / "layout_metadata.json"
    if layout_marker.exists() and is_dataset_dir(target_dir):
        console.print(f"[green]Dataset já existe:[/green] {target_dir}")
        return target_dir

    dataset_meta = _load_json(source_dir / "dataset_metadata.json")
    target_dir.mkdir(parents=True, exist_ok=True)

    ref_windows_dir = target_dir / "references" / "windows"
    individuals_out_dir = target_dir / "individuals"

    global_window_metadata: Dict[str, Dict] = {}
    gene_order = list(dataset_meta.get("genes", []))

    individual_dirs = list(_iter_individual_dirs(source_dir))
    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Materializando individuos...", total=len(individual_dirs))
        for ind_dir in individual_dirs:
            sample_id = ind_dir.name
            ind_meta_path = ind_dir / "individual_metadata.json"
            if not ind_meta_path.exists():
                progress.update(task, advance=1)
                continue

            ind_meta = _load_json(ind_meta_path)
            windows = list(ind_meta.get("windows", []))
            window_metadata = ind_meta.get("window_metadata", {})

            new_ind_meta = {
                "sample_id": ind_meta.get("sample_id", sample_id),
                "family_id": ind_meta.get("family_id", "0"),
                "sex": ind_meta.get("sex", 0),
                "population": ind_meta.get("population", ""),
                "superpopulation": ind_meta.get("superpopulation", ""),
                "longevity": ind_meta.get("longevity", False),
                "windows": windows,
                "created_at": ind_meta.get("created_at"),
                "last_updated": ind_meta.get("last_updated"),
            }
            if "frog_likelihood" in ind_meta:
                new_ind_meta["frog_likelihood"] = ind_meta["frog_likelihood"]
            if "frog_population_names" in ind_meta:
                new_ind_meta["frog_population_names"] = ind_meta["frog_population_names"]

            _write_json(individuals_out_dir / sample_id / "individual_metadata.json", new_ind_meta)

            for window_name in windows:
                src_window_dir = ind_dir / "windows" / window_name
                if not src_window_dir.exists():
                    continue

                if window_name in window_metadata and window_name not in global_window_metadata:
                    global_window_metadata[window_name] = dict(window_metadata[window_name])
                    _write_json(ref_windows_dir / window_name / "window_metadata.json", global_window_metadata[window_name])

                ref_src = src_window_dir / "ref.window.fa"
                if ref_src.exists():
                    _copy_file(ref_src, ref_windows_dir / window_name / "ref.window.fa")

                dst_window_dir = individuals_out_dir / sample_id / "windows" / window_name
                for file_name in (f"{sample_id}.H1.window.fixed.fa", f"{sample_id}.H2.window.fixed.fa"):
                    src = src_window_dir / file_name
                    if src.exists():
                        _copy_file(src, dst_window_dir / file_name)

            for pred_dir_name in ("predictions_H1", "predictions_H2"):
                src_pred_dir = src_window_dir / pred_dir_name
                dst_pred_dir = dst_window_dir / pred_dir_name
                if not src_pred_dir.exists():
                    continue
                dst_pred_dir.mkdir(parents=True, exist_ok=True)
                for pred_file in src_pred_dir.iterdir():
                    if pred_file.is_file():
                        _copy_file(pred_file, dst_pred_dir / pred_file.name)

            # Preserve raw per-window VCF artefacts when available.
            for suffix in (
                f"{sample_id}.window.vcf.gz",
                f"{sample_id}.window.vcf.gz.tbi",
                f"{sample_id}.window.consensus_ready.vcf.gz",
                f"{sample_id}.window.consensus_ready.vcf.gz.tbi",
            ):
                src_variant = src_window_dir / suffix
                if src_variant.exists():
                    _copy_file(src_variant, dst_window_dir / suffix)

            progress.update(task, advance=1)

    out_meta = dict(dataset_meta)
    out_meta["layout_version"] = LAYOUT_VERSION
    out_meta["source_dataset_dir"] = str(source_dir.resolve())
    out_meta["window_catalog"] = global_window_metadata
    out_meta.setdefault("raw_variant_sources", {
        "vcf_pattern": os.environ.get("KG1000_VCF_PATTERN"),
        "vcf_root_dir": os.environ.get("KG1000_VCF_ROOT_DIR"),
        "notes": "Set KG1000_VCF_PATTERN or KG1000_VCF_ROOT_DIR to preserve the chromosome-level 1000G VCF source used to generate consensus FASTAs.",
    })
    if gene_order:
        out_meta["genes"] = gene_order

    for window_name, meta in global_window_metadata.items():
        chromosome = str(meta.get("chromosome", ""))
        chrom_token = chromosome.replace("chr", "") if chromosome else None
        if chrom_token:
            meta.setdefault("raw_variant_source", {
                "chromosome": chromosome,
                "vcf_pattern": out_meta.get("raw_variant_sources", {}).get("vcf_pattern"),
                "vcf_root_dir": out_meta.get("raw_variant_sources", {}).get("vcf_root_dir"),
                "vcf_pattern_example": None if not out_meta.get("raw_variant_sources", {}).get("vcf_pattern") else out_meta["raw_variant_sources"]["vcf_pattern"].replace("{chrom}", chrom_token),
            })
    _write_json(target_dir / "dataset_metadata.json", out_meta)
    _write_json(layout_marker, {
        "layout_version": LAYOUT_VERSION,
        "source_dataset_dir": str(source_dir.resolve()),
        "target_dataset_dir": str(target_dir.resolve()),
    })

    console.print(f"[green]Dataset materializado em[/green] {target_dir}")
    return target_dir
