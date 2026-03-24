#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit StandardScaler + IncrementalPCA on the training split (same pipeline as PCA cache)
and plot cumulative explained variance ratio vs number of components.

IncrementalPCA only exposes variance for the k components you fit: set --max-components
to the largest k you want on the x-axis (capped by n_train and LAPACK-safe limit for D).

Example:
    python3 plot_sklearn_pca_variance.py --config configs/genes_1000.yaml
    python3 plot_sklearn_pca_variance.py --config configs/genes_1000.yaml --max-components 495 \\
        --output runs/pca_var_curve.png

Note: explained_variance_ratio_ from IncrementalPCA is a streaming approximation;
use it as a guide for choosing k, not as exact population PCA.
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from sklearn_pca_cache import (
    compute_sklearn_pca_effective_k,
    fit_incremental_pca_on_train,
    fit_standard_scaler_incremental,
    load_neural_ancestry_predictor_for_cli,
)


def _first_k_at_threshold(cumulative: np.ndarray, thresholds: List[float]) -> Dict[float, Optional[int]]:
    """1-based component index where cumulative variance first reaches threshold."""
    out: Dict[float, Optional[int]] = {}
    for t in thresholds:
        idx = np.searchsorted(cumulative, t, side="left")
        if idx >= len(cumulative):
            out[t] = None
        else:
            out[t] = int(idx + 1)
    return out


def run_plot(
    config_path: Path,
    output_path: Path,
    max_components: Optional[int],
    show: bool,
) -> None:
    import os

    if not show:
        os.environ.setdefault("MPLBACKEND", "Agg")
    nap = load_neural_ancestry_predictor_for_cli()
    if not show:
        import matplotlib

        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    console = nap.console
    load_config = nap.load_config
    get_dataset_cache_dir = nap.get_dataset_cache_dir
    prepare_data = nap.prepare_data
    validate_cache = nap.validate_cache

    config = load_config(config_path)
    dataset_cache_dir = Path(get_dataset_cache_dir(config))
    if not dataset_cache_dir.exists() or not validate_cache(dataset_cache_dir, config):
        console.print(
            f"[red]Cache do dataset inválido ou ausente: {dataset_cache_dir}[/red]\n"
            "[yellow]Gere o cache do dataset antes (modo train do neural_ancestry_predictor).[/yellow]"
        )
        sys.exit(1)

    exp_dir = Path(config["dataset_input"]["processed_cache_dir"]) / "_pca_variance_plot"
    exp_dir.mkdir(parents=True, exist_ok=True)
    _full, train_loader, _val_loader, _test_loader = prepare_data(config, exp_dir)

    n_train = len(train_loader.dataset)
    first = next(iter(train_loader))
    n_features = int(np.prod(first[0].shape[1:]))

    cfg = copy.deepcopy(config)
    if max_components is not None:
        cfg.setdefault("model", {}).setdefault("sklearn", {})["pca_components"] = int(max_components)

    align = bool(cfg.get("model", {}).get("sklearn", {}).get("pca_align_n_train", False))
    effective_k, pca_req = compute_sklearn_pca_effective_k(
        cfg, n_train=n_train, n_features=n_features, log=console.print
    )

    console.print(
        f"[cyan]Ajustando PCA para curva de variância: k={effective_k} "
        f"(pedido no YAML/CLI={pca_req}, n_train={n_train}, D={n_features})[/cyan]"
    )

    scaler = fit_standard_scaler_incremental(train_loader, rich_console=console, progress_desc="Variância PCA: StandardScaler")
    pca = fit_incremental_pca_on_train(
        train_loader,
        scaler,
        effective_k,
        log=console.print,
        forbid_tail_padding=align,
        rich_console=console,
        progress_desc="Variância PCA: IncrementalPCA",
    )

    ratios = np.asarray(pca.explained_variance_ratio_, dtype=np.float64)
    cumulative = np.cumsum(ratios)
    k_axis = np.arange(1, len(cumulative) + 1)

    thresholds = [0.8, 0.85, 0.9, 0.95, 0.99]
    crosses = _first_k_at_threshold(cumulative, thresholds)
    console.print("[cyan]Primeiro k com variância acumulada ≥ limiar (aprox. IncrementalPCA):[/cyan]")
    for t in thresholds:
        kk = crosses[t]
        if kk is None:
            console.print(f"  • {t:.0%}: não atingido com k≤{effective_k}")
        else:
            console.print(f"  • {t:.0%}: k ≥ {kk}  (cum={cumulative[kk - 1]:.4f})")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)

    ax0 = axes[0]
    ax0.plot(k_axis, cumulative, color="C0", linewidth=1.5)
    for t in (0.8, 0.9, 0.95):
        ax0.axhline(t, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax0.set_xlabel("Number of components")
    ax0.set_ylabel("Cumulative explained variance ratio")
    ax0.set_title("IncrementalPCA (train) — cumulative variance")
    ax0.set_ylim(0.0, min(1.05, float(cumulative[-1]) * 1.05 + 0.02))
    ax0.grid(True, alpha=0.3)

    ax1 = axes[1]
    ax1.bar(k_axis, ratios, width=1.0, color="C1", alpha=0.85)
    ax1.set_xlabel("Component index")
    ax1.set_ylabel("Explained variance ratio (per component)")
    ax1.set_title("Per-component ratio (scree-style)")
    ax1.grid(True, axis="y", alpha=0.3)

    fig.suptitle(f"PCA variance diagnostic (k={effective_k}, n_train={n_train}, D={n_features})", fontsize=11)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    console.print(f"[green]Figura guardada em {output_path.resolve()}[/green]")
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot cumulative explained variance for IncrementalPCA on the training set."
    )
    parser.add_argument("--config", type=str, required=True, help="YAML de configuração")
    parser.add_argument(
        "--max-components",
        type=int,
        default=None,
        help="Override model.sklearn.pca_components for this fit (curve length = effective k after caps/align)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="pca_cumulative_variance.png",
        help="PNG path for the figure",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive matplotlib window after saving",
    )
    args = parser.parse_args()
    run_plot(Path(args.config), Path(args.output), args.max_components, args.show)


if __name__ == "__main__":
    main()
