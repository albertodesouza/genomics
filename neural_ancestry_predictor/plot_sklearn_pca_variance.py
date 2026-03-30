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

    # Replot sem refazer o fit (usa cache existente em pca_cache/):
    python3 plot_sklearn_pca_variance.py \\
        --from-cache /dados/.../pca_cache/rna_seq_..._pca300 \\
        --output runs/pca_var_curve_v2.png

Note: explained_variance_ratio_ from IncrementalPCA is a streaming approximation;
use it as a guide for choosing k, not as exact population PCA.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import joblib
import numpy as np

from sklearn_pca_cache import (
    compute_sklearn_pca_effective_k,
    fit_incremental_pca_on_train,
    fit_standard_scaler_incremental,
    load_neural_ancestry_predictor_for_cli,
    SCALER_PCA_FILENAME,
    METADATA_FILENAME,
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
    *,
    paper_clean: bool = False,
    paper_k: int = 300,
    paper_threshold: float = 0.95,
    paper_no_title: bool = False,
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

    _save_variance_plot(
        cumulative,
        ratios,
        k_axis,
        effective_k,
        n_train,
        n_features,
        output_path,
        show,
        paper_clean=paper_clean,
        paper_k=paper_k,
        paper_threshold=paper_threshold,
        paper_no_title=paper_no_title,
    )


def _save_variance_plot(
    cumulative: np.ndarray,
    ratios: np.ndarray,
    k_axis: np.ndarray,
    effective_k: int,
    n_train: int,
    n_features: int,
    output_path: Path,
    show: bool,
    *,
    paper_clean: bool = False,
    paper_k: int = 300,
    paper_threshold: float = 0.95,
    paper_no_title: bool = False,
) -> None:
    """Gera e salva a figura de variância PCA (cumulative + scree)."""
    import os
    if not show:
        os.environ.setdefault("MPLBACKEND", "Agg")
        import matplotlib
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    thresholds = [0.8, 0.85, 0.9, 0.95, 0.99]
    crosses = _first_k_at_threshold(cumulative, thresholds)

    if paper_clean:
        fig, ax = plt.subplots(1, 1, figsize=(7.2, 4.2), constrained_layout=True)
        ax.plot(k_axis, cumulative, color="black", linewidth=1.8, label="Cumulative variance")
        ax.axhline(
            paper_threshold,
            color="dimgray",
            linestyle="--",
            linewidth=1.2,
            label=f"{paper_threshold:.0%} threshold",
        )
        k_mark = int(max(1, min(paper_k, len(cumulative))))
        y_mark = float(cumulative[k_mark - 1])
        ax.scatter([k_mark], [y_mark], color="tab:red", s=36, zorder=3)
        ax.axvline(k_mark, color="tab:red", linestyle=":", linewidth=1.0, alpha=0.9)
        ax.annotate(
            f"k={k_mark}, cum={y_mark:.3f}",
            xy=(k_mark, y_mark),
            xytext=(8, 10),
            textcoords="offset points",
            fontsize=9,
            color="tab:red",
        )
        ax.set_xlabel("Number of components (k)")
        ax.set_ylabel("Cumulative explained variance ratio")
        ax.set_ylim(0.0, min(1.02, float(max(cumulative[-1], paper_threshold)) + 0.03))
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, loc="lower right")
        if not paper_no_title:
            ax.set_title("IncrementalPCA cumulative variance")
    else:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)

        ax0 = axes[0]
        ax0.plot(k_axis, cumulative, color="C0", linewidth=1.5)
        for t in (0.8, 0.9, 0.95):
            kk = crosses.get(t)
            if kk is not None:
                ax0.axvline(kk, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
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
    print(f"Figura guardada em {output_path.resolve()}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def run_plot_from_cache(
    cache_dir: Path,
    output_path: Path,
    show: bool,
    *,
    paper_clean: bool = False,
    paper_k: int = 300,
    paper_threshold: float = 0.95,
    paper_no_title: bool = False,
) -> None:
    """
    Replot a curva de variância PCA carregando o objeto PCA já treinado
    do scaler_pca.joblib em cache_dir, sem refazer o fit.
    """
    cache_dir = Path(cache_dir)
    joblib_path = cache_dir / SCALER_PCA_FILENAME
    meta_path = cache_dir / METADATA_FILENAME

    if not joblib_path.exists():
        print(f"Erro: {joblib_path} não encontrado.")
        sys.exit(1)

    bundle = joblib.load(joblib_path)
    pca = bundle["pca"]

    # Ler metadados para exibir n_train e n_features no título
    n_train, n_features, effective_k = 0, 0, 0
    if meta_path.exists():
        with open(meta_path, "r") as f:
            meta = json.load(f)
        n_train = int(meta.get("n_train", 0))
        n_features = int(meta.get("n_features_original", 0))
        effective_k = int(meta.get("pca_n_components_effective", len(pca.explained_variance_ratio_)))
    else:
        effective_k = len(pca.explained_variance_ratio_)
        print(f"Aviso: {meta_path} não encontrado; título do gráfico ficará incompleto.")

    print(f"Carregando PCA do cache: {joblib_path}")
    print(f"  k={effective_k}, n_train={n_train}, D={n_features}")

    ratios = np.asarray(pca.explained_variance_ratio_, dtype=np.float64)
    cumulative = np.cumsum(ratios)
    k_axis = np.arange(1, len(cumulative) + 1)

    thresholds = [0.8, 0.85, 0.9, 0.95, 0.99]
    crosses = _first_k_at_threshold(cumulative, thresholds)
    print("Primeiro k com variância acumulada >= limiar:")
    for t in thresholds:
        kk = crosses[t]
        if kk is None:
            print(f"  {t:.0%}: não atingido com k<={effective_k}")
        else:
            print(f"  {t:.0%}: k >= {kk}  (cum={cumulative[kk - 1]:.4f})")

    _save_variance_plot(
        cumulative,
        ratios,
        k_axis,
        effective_k,
        n_train,
        n_features,
        output_path,
        show,
        paper_clean=paper_clean,
        paper_k=paper_k,
        paper_threshold=paper_threshold,
        paper_no_title=paper_no_title,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot cumulative explained variance for IncrementalPCA on the training set."
    )
    # Modo 1: fit novo (requer --config)
    parser.add_argument("--config", type=str, default=None, help="YAML de configuração (modo fit)")
    parser.add_argument(
        "--max-components",
        type=int,
        default=None,
        help="Override model.sklearn.pca_components for this fit (curve length = effective k after caps/align)",
    )
    # Modo 2: replot a partir do cache (não faz fit)
    parser.add_argument(
        "--from-cache",
        type=str,
        default=None,
        metavar="CACHE_DIR",
        help=(
            "Caminho para um diretório de cache PCA existente (ex: pca_cache/rna_seq_..._pca300). "
            "Carrega scaler_pca.joblib e replota sem refazer o fit."
        ),
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
    parser.add_argument(
        "--paper-clean",
        action="store_true",
        help="Generate a clean single-panel plot suitable for papers.",
    )
    parser.add_argument(
        "--paper-k",
        type=int,
        default=300,
        help="k value to highlight with a marker in --paper-clean mode.",
    )
    parser.add_argument(
        "--paper-threshold",
        type=float,
        default=0.95,
        help="Horizontal dashed threshold in --paper-clean mode (default 0.95).",
    )
    parser.add_argument(
        "--paper-no-title",
        action="store_true",
        help="Hide title in --paper-clean mode.",
    )
    args = parser.parse_args()

    if args.from_cache:
        run_plot_from_cache(
            Path(args.from_cache),
            Path(args.output),
            args.show,
            paper_clean=args.paper_clean,
            paper_k=args.paper_k,
            paper_threshold=args.paper_threshold,
            paper_no_title=args.paper_no_title,
        )
    elif args.config:
        run_plot(
            Path(args.config),
            Path(args.output),
            args.max_components,
            args.show,
            paper_clean=args.paper_clean,
            paper_k=args.paper_k,
            paper_threshold=args.paper_threshold,
            paper_no_title=args.paper_no_title,
        )
    else:
        parser.error("Forneça --config (fit novo) ou --from-cache (replot do cache).")


if __name__ == "__main__":
    main()
