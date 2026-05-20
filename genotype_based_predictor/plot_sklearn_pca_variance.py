#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build/reuse sklearn PCA cache and plot cumulative explained variance."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import joblib
import numpy as np
from rich.console import Console

from genotype_based_predictor.config import get_dataset_cache_dir, load_config
from genotype_based_predictor.data_pipeline import prepare_data

NEURAL_PREDICTOR_DIR = ROOT / "neural_ancestry_predictor"
if str(NEURAL_PREDICTOR_DIR) not in sys.path:
    sys.path.insert(0, str(NEURAL_PREDICTOR_DIR))

from sklearn_pca_cache import (  # noqa: E402
    METADATA_FILENAME,
    SCALER_PCA_FILENAME,
    ensure_sklearn_pca_cache,
)

console = Console()


def _first_k_at_threshold(cumulative: np.ndarray, thresholds: List[float]) -> Dict[float, Optional[int]]:
    out: Dict[float, Optional[int]] = {}
    for threshold in thresholds:
        idx = np.searchsorted(cumulative, threshold, side="left")
        out[threshold] = None if idx >= len(cumulative) else int(idx + 1)
    return out


def _save_variance_json(
    path: Path,
    ratios: np.ndarray,
    cumulative: np.ndarray,
    thresholds: List[float],
    crosses: Dict[float, Optional[int]],
    metadata: Dict,
) -> None:
    payload = {
        "source": "sklearn_pca_cache",
        "cache_dir": metadata.get("cache_dir"),
        "effective_k": int(metadata.get("pca_n_components_effective", len(cumulative))),
        "pca_components_requested": int(metadata.get("pca_components_requested", len(cumulative))),
        "n_train": int(metadata.get("n_train", 0)),
        "n_val": int(metadata.get("n_val", 0)),
        "n_test": int(metadata.get("n_test", 0)),
        "n_features": int(metadata.get("n_features_original", 0)),
        "total_explained_variance_ratio": float(cumulative[-1]) if len(cumulative) else 0.0,
        "thresholds": {
            f"{threshold:.2f}": {
                "first_k": None if crosses[threshold] is None else int(crosses[threshold]),
                "cumulative_at_first_k": None
                if crosses[threshold] is None
                else float(cumulative[crosses[threshold] - 1]),
            }
            for threshold in thresholds
        },
        "explained_variance_ratio": [float(x) for x in ratios],
        "cumulative_explained_variance_ratio": [float(x) for x in cumulative],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    console.print(f"[green]✓ JSON salvo em: {path}[/green]")


def _save_plot(path: Path, ratios: np.ndarray, cumulative: np.ndarray, metadata: Dict) -> None:
    import os

    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    k_axis = np.arange(1, len(cumulative) + 1)
    crosses = _first_k_at_threshold(cumulative, [0.8, 0.9, 0.95])
    backend = metadata.get("pca_backend", "incremental")
    title_backend = "IncrementalPCA" if backend == "incremental" else str(backend)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axes[0].plot(k_axis, cumulative, color="C0", linewidth=1.5)
    for threshold, kk in crosses.items():
        axes[0].axhline(threshold, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        if kk is not None:
            axes[0].axvline(kk, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
    axes[0].set_xlabel("Number of PCA components")
    axes[0].set_ylabel("Cumulative explained variance ratio")
    axes[0].set_title("Cumulative variance")
    axes[0].set_ylim(0.0, min(1.05, max(0.05, float(cumulative[-1]) * 1.05)))
    axes[0].grid(True, alpha=0.3)

    axes[1].bar(k_axis, ratios, width=1.0, color="C1", alpha=0.85)
    axes[1].set_xlabel("PCA component")
    axes[1].set_ylabel("Explained variance ratio")
    axes[1].set_title("Per-component variance")
    axes[1].grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        f"{title_backend} variance "
        f"(k={metadata.get('pca_n_components_effective')}, "
        f"n_train={metadata.get('n_train')}, D={metadata.get('n_features_original')})",
        fontsize=11,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    console.print(f"[green]✓ Figura salva em: {path}[/green]")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build/reuse sklearn PCA cache and plot explained variance.")
    parser.add_argument("config", type=Path, help="Pipeline YAML config.")
    parser.add_argument("--max-components", type=int, default=None, help="Override model.sklearn.pca_components.")
    parser.add_argument("--force", action="store_true", help="Rebuild PCA cache even if it already exists.")
    parser.add_argument("--output", type=Path, required=True, help="PNG output path.")
    parser.add_argument("--json-output", type=Path, required=True, help="JSON output path.")
    args = parser.parse_args()

    config = load_config(args.config)
    if args.max_components is not None:
        config.model.sklearn.pca_components = int(args.max_components)
    if args.force:
        config.debug.force_pca_cache_rebuild = True

    experiment_dir = Path(config.dataset_input.processed_cache_dir) / "_pca_variance_plot"
    experiment_dir.mkdir(parents=True, exist_ok=True)
    _full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)

    dataset_cache_dir = get_dataset_cache_dir(config)
    pca_dir = ensure_sklearn_pca_cache(
        config.model_dump(),
        dataset_cache_dir,
        train_loader,
        val_loader,
        test_loader,
        force=config.debug.force_pca_cache_rebuild,
        log=console.print,
        rich_console=console,
    )
    bundle = joblib.load(pca_dir / SCALER_PCA_FILENAME)
    pca = bundle["pca"]
    with open(pca_dir / METADATA_FILENAME, "r", encoding="utf-8") as f:
        metadata = json.load(f)
    metadata["cache_dir"] = str(pca_dir.resolve())

    ratios = np.asarray(pca.explained_variance_ratio_, dtype=np.float64)
    cumulative = np.cumsum(ratios)
    thresholds = [0.8, 0.85, 0.9, 0.95, 0.99]
    crosses = _first_k_at_threshold(cumulative, thresholds)

    console.print("[cyan]Primeiro k com variância acumulada >= limiar:[/cyan]")
    for threshold in thresholds:
        kk = crosses[threshold]
        if kk is None:
            console.print(f"  {threshold:.0%}: não atingido com k<={len(cumulative)}")
        else:
            console.print(f"  {threshold:.0%}: k >= {kk} (cum={cumulative[kk - 1]:.4f})")

    _save_variance_json(args.json_output, ratios, cumulative, thresholds, crosses, metadata)
    _save_plot(args.output, ratios, cumulative, metadata)
    console.print(f"[green]✓ Cache PCA: {pca_dir}[/green]")


if __name__ == "__main__":
    main()
