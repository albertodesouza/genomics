#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hyperparameter search for sklearn baselines (SVM, RF, XGBoost) on cached PCA features.

Uses the same PCA cache as normal training (StandardScaler + IncrementalPCA matrices on disk).
Each trial: fit classifier on X_train, evaluate on train and validation; save metrics and plots.

Example:
    python3 tune_sklearn_baselines.py --config configs/genes_1000.yaml --output-dir runs/sklearn_tune_001

Custom grids (JSON):
    python3 tune_sklearn_baselines.py --config configs/genes_1000.yaml \\
        --grid-json configs/sklearn_tune_grid.example.json --output-dir runs/sklearn_tune_custom
"""

from __future__ import annotations

import os

# Prefer non-interactive backend before neural_ancestry_predictor runs (DISPLAY can make it pick TkAgg).
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import copy
import csv
import itertools
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from sklearn_pca_cache import METADATA_FILENAME, ensure_sklearn_pca_cache, load_neural_ancestry_predictor_for_cli


def _use_agg_backend() -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)


# Default search spaces (override with --grid-json)
DEFAULT_GRIDS: Dict[str, Dict[str, List[Any]]] = {
    "SVM": {
        "C": [0.01, 0.1, 1.0, 10.0],
        "calibrate_probabilities": [False, True],
    },
    "RF": {
        "n_estimators": [100, 200, 400],
        "max_depth": [None, 10, 20],
    },
    "XGBOOST": {
        "n_estimators": [100, 200],
        "max_depth": [4, 6, 8],
        "learning_rate": [0.05, 0.1],
    },
}


def _grid_product(grid: Dict[str, List[Any]]) -> Iterable[Dict[str, Any]]:
    keys = list(grid.keys())
    values = [grid[k] for k in keys]
    for combo in itertools.product(*values):
        yield dict(zip(keys, combo))


def _merge_config(
    base: Dict[str, Any],
    model_type: str,
    svm: Optional[Dict[str, Any]] = None,
    rf: Optional[Dict[str, Any]] = None,
    xgb: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    cfg = copy.deepcopy(base)
    cfg.setdefault("model", {})["type"] = model_type
    sk = cfg.setdefault("model", {}).setdefault("sklearn", {})
    if svm:
        sk.setdefault("svm", {}).update(svm)
    if rf:
        sk.setdefault("random_forest", {}).update(rf)
    if xgb:
        sk.setdefault("xgboost", {}).update(xgb)
    return cfg


def _short_label(model: str, svm: Dict, rf: Dict, xgb: Dict) -> str:
    if model == "SVM":
        cal = svm.get("calibrate_probabilities", False)
        return f"C={svm.get('C')} cal={cal}"
    if model == "RF":
        md = rf.get("max_depth")
        md_s = "None" if md is None else str(md)
        return f"n={rf.get('n_estimators')} md={md_s}"
    if model == "XGBOOST":
        return f"nt={xgb.get('n_estimators')} d={xgb.get('max_depth')} lr={xgb.get('learning_rate')}"
    return model


def _load_grid_json(path: Path) -> Dict[str, Dict[str, List[Any]]]:
    with open(path, "r") as f:
        raw = json.load(f)
    out: Dict[str, Dict[str, List[Any]]] = {}
    for model, spec in raw.items():
        m = model.upper()
        if m not in ("SVM", "RF", "XGBOOST"):
            continue
        out[m] = {k: (v if isinstance(v, list) else [v]) for k, v in spec.items()}
    return out


def _run_trials(
    nap: Any,
    config: Dict[str, Any],
    dataset_cache_dir: Path,
    train_loader: Any,
    val_loader: Any,
    test_loader: Any,
    full_dataset: Any,
    models: List[str],
    grids: Dict[str, Dict[str, List[Any]]],
    force_pca: bool,
) -> Tuple[List[Dict[str, Any]], Path]:
    _ensure = nap._ensure_sklearn_classification_target
    _seed = nap._sklearn_effective_random_seed
    build_clf = nap.build_sklearn_classifier
    metrics_dict = nap.sklearn_metrics_dict

    _ensure(config)
    random_seed = _seed(config)

    pca_dir = ensure_sklearn_pca_cache(
        config,
        dataset_cache_dir,
        train_loader,
        val_loader,
        test_loader,
        force=force_pca,
        log=nap.console.print,
        rich_console=nap.console,
    )
    X_train = np.load(pca_dir / "X_train.npy")
    y_train = np.load(pca_dir / "y_train.npy")
    X_val = np.load(pca_dir / "X_val.npy")
    y_val = np.load(pca_dir / "y_val.npy")

    with open(pca_dir / METADATA_FILENAME, "r") as f:
        pca_meta = json.load(f)
    effective_k = int(pca_meta["pca_n_components_effective"])

    results: List[Dict[str, Any]] = []
    trial_id = 0

    for model in models:
        grid = grids.get(model)
        if not grid:
            nap.console.print(f"[yellow]Sem grid para {model}, a ignorar.[/yellow]")
            continue
        if model == "XGBOOST":
            try:
                import xgboost  # noqa: F401
            except ImportError:
                nap.console.print("[yellow]XGBOOST: pacote xgboost não instalado — a saltar.[/yellow]")
                continue

        for params in _grid_product(grid):
            svm_u: Dict[str, Any] = {}
            rf_u: Dict[str, Any] = {}
            xgb_u: Dict[str, Any] = {}
            if model == "SVM":
                svm_u = dict(params)
            elif model == "RF":
                rf_u = dict(params)
            elif model == "XGBOOST":
                xgb_u = dict(params)

            cfg_t = _merge_config(config, model, svm=svm_u or None, rf=rf_u or None, xgb=xgb_u or None)
            t0 = time.perf_counter()
            clf = build_clf(cfg_t, model, random_seed, n_train_samples=len(y_train))
            clf.fit(X_train, y_train)
            fit_s = time.perf_counter() - t0

            y_pred_tr = clf.predict(X_train)
            y_pred_va = clf.predict(X_val)
            m_tr = metrics_dict(y_train, y_pred_tr, full_dataset)
            m_va = metrics_dict(y_val, y_pred_va, full_dataset)

            trial_id += 1
            row = {
                "trial_id": trial_id,
                "model": model,
                "pca_k": effective_k,
                "pca_cache": pca_dir.name,
                "params": params,
                "svm_params": svm_u,
                "rf_params": rf_u,
                "xgb_params": xgb_u,
                "label": _short_label(model, svm_u, rf_u, xgb_u),
                "train_accuracy": m_tr["accuracy"],
                "train_precision": m_tr["precision"],
                "train_recall": m_tr["recall"],
                "train_f1": m_tr["f1"],
                "val_accuracy": m_va["accuracy"],
                "val_precision": m_va["precision"],
                "val_recall": m_va["recall"],
                "val_f1": m_va["f1"],
                "fit_seconds": round(fit_s, 4),
            }
            results.append(row)
            nap.console.print(
                f"[cyan]trial {trial_id}[/cyan] {model} {_short_label(model, svm_u, rf_u, xgb_u)} "
                f"| val_acc={m_va['accuracy']:.4f} val_f1={m_va['f1']:.4f} ({fit_s:.2f}s)"
            )

    return results, pca_dir


def _save_results_json_csv(results: List[Dict[str, Any]], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    serializable = []
    for r in results:
        d = {k: v for k, v in r.items() if k not in ("svm_params", "rf_params", "xgb_params")}
        d["params_json"] = json.dumps(r["params"], sort_keys=True)
        serializable.append(d)
    with open(out_dir / "tuning_results.json", "w") as f:
        json.dump(serializable, f, indent=2)

    if not results:
        return
    fieldnames = [
        "trial_id",
        "model",
        "pca_k",
        "pca_cache",
        "params_json",
        "label",
        "train_accuracy",
        "train_precision",
        "train_recall",
        "train_f1",
        "val_accuracy",
        "val_precision",
        "val_recall",
        "val_f1",
        "fit_seconds",
    ]
    with open(out_dir / "tuning_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for d in serializable:
            w.writerow(d)


def _plot_validation_summary(results: List[Dict[str, Any]], out_dir: Path) -> None:
    _use_agg_backend()
    import matplotlib.pyplot as plt

    if not results:
        return
    models = sorted({r["model"] for r in results})
    metric_keys = [
        ("val_accuracy", "Validation accuracy"),
        ("val_precision", "Validation precision (weighted)"),
        ("val_recall", "Validation recall (weighted)"),
        ("val_f1", "Validation F1 (weighted)"),
    ]

    best_per_model: Dict[str, Dict[str, Any]] = {}
    for m in models:
        rows = [r for r in results if r["model"] == m]
        best = max(rows, key=lambda x: x["val_accuracy"])
        best_per_model[m] = best

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    for ax, (key, title) in zip(axes.flat, metric_keys):
        xs = np.arange(len(models))
        vals = [best_per_model[m][key] for m in models]
        ax.bar(xs, vals, color=["C0", "C1", "C2"][: len(models)])
        ax.set_xticks(xs)
        ax.set_xticklabels(models)
        ax.set_ylim(0.0, 1.05)
        ax.set_title(title)
        ax.set_ylabel("Score")
        for i, v in enumerate(vals):
            ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=9)
    fig.suptitle("Best trial per model (by validation accuracy)", fontsize=12)
    fig.savefig(out_dir / "tuning_best_per_model_val_metrics.png", dpi=160)
    plt.close(fig)


def _plot_top_trials_per_model(
    results: List[Dict[str, Any]],
    out_dir: Path,
    top_k: int = 12,
    label_max_len: int = 44,
) -> None:
    _use_agg_backend()
    import matplotlib.pyplot as plt

    if not results:
        return
    models = sorted({r["model"] for r in results})
    n = len(models)
    fig_h = max(4.5, 1.8 * n)
    fig, axes = plt.subplots(n, 1, figsize=(11, fig_h), constrained_layout=True)
    if n == 1:
        axes = [axes]
    for ax, model in zip(axes, models):
        rows = [r for r in results if r["model"] == model]
        rows = sorted(rows, key=lambda x: -x["val_accuracy"])[:top_k]
        if not rows:
            continue
        labels = [r["label"][:label_max_len] for r in rows]
        vals = [r["val_accuracy"] for r in rows]
        y_pos = np.arange(len(labels))
        ax.barh(y_pos, vals, color="steelblue")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Validation accuracy")
        ax.set_xlim(0.0, 1.02)
        ax.set_title(f"{model}: top {len(rows)} trials by validation accuracy")
        for i, v in enumerate(vals):
            ax.text(v + 0.008, i, f"{v:.3f}", va="center", fontsize=7)
    fig.savefig(out_dir / "tuning_top_trials_val_accuracy.png", dpi=160)
    plt.close(fig)


def _plot_ranked_val_accuracy_curves(results: List[Dict[str, Any]], out_dir: Path) -> None:
    _use_agg_backend()
    import matplotlib.pyplot as plt

    if not results:
        return
    models = sorted({r["model"] for r in results})
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for i, model in enumerate(models):
        rows = [r for r in results if r["model"] == model]
        vals = sorted((r["val_accuracy"] for r in rows), reverse=True)
        ranks = np.arange(1, len(vals) + 1)
        ax.plot(ranks, vals, "o-", markersize=3, linewidth=1.2, label=f"{model} (n={len(vals)})", c=f"C{i}")
    ax.set_xlabel("Trial rank (1 = best)")
    ax.set_ylabel("Validation accuracy")
    ax.set_title("Ranked validation accuracy curves (less cluttered)")
    ax.set_ylim(0.0, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(out_dir / "tuning_ranked_val_accuracy_curves.png", dpi=160)
    plt.close(fig)


def _plot_pair_heatmaps_if_two_params(results: List[Dict[str, Any]], out_dir: Path) -> None:
    """When exactly 2 params vary for a model, create a validation-accuracy heatmap."""
    _use_agg_backend()
    import matplotlib.pyplot as plt

    for model in sorted({r["model"] for r in results}):
        rows = [r for r in results if r["model"] == model]
        if len(rows) < 3:
            continue
        pnames = sorted({k for r in rows for k in r["params"].keys()})
        varying = []
        for pn in pnames:
            vals = {json.dumps(r["params"].get(pn), sort_keys=True) for r in rows}
            if len(vals) > 1:
                varying.append(pn)
        if len(varying) != 2:
            continue
        p0, p1 = varying
        xvals = sorted({str(r["params"][p0]) for r in rows})
        yvals = sorted({str(r["params"][p1]) for r in rows})
        grid = np.full((len(yvals), len(xvals)), np.nan, dtype=np.float64)
        x_idx = {v: i for i, v in enumerate(xvals)}
        y_idx = {v: i for i, v in enumerate(yvals)}
        for r in rows:
            xi = x_idx[str(r["params"][p0])]
            yi = y_idx[str(r["params"][p1])]
            grid[yi, xi] = max(grid[yi, xi], float(r["val_accuracy"])) if not np.isnan(grid[yi, xi]) else float(r["val_accuracy"])

        fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
        im = ax.imshow(grid, cmap="viridis", vmin=np.nanmin(grid), vmax=np.nanmax(grid))
        fig.colorbar(im, ax=ax, label="Validation accuracy")
        ax.set_xticks(np.arange(len(xvals)))
        ax.set_yticks(np.arange(len(yvals)))
        ax.set_xticklabels(xvals, rotation=45, ha="right")
        ax.set_yticklabels(yvals)
        ax.set_xlabel(p0)
        ax.set_ylabel(p1)
        ax.set_title(f"{model}: validation accuracy heatmap")
        for yi in range(len(yvals)):
            for xi in range(len(xvals)):
                if np.isnan(grid[yi, xi]):
                    continue
                ax.text(xi, yi, f"{grid[yi, xi]:.3f}", ha="center", va="center", color="white", fontsize=7)
        fig.savefig(out_dir / f"tuning_{model.lower()}_val_accuracy_heatmap.png", dpi=160)
        plt.close(fig)


def _plot_train_vs_val_accuracy(results: List[Dict[str, Any]], out_dir: Path) -> None:
    _use_agg_backend()
    import matplotlib.pyplot as plt

    if not results:
        return
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    models = sorted({r["model"] for r in results})
    for i, m in enumerate(models):
        pts = [r for r in results if r["model"] == m]
        ax.scatter(
            [p["train_accuracy"] for p in pts],
            [p["val_accuracy"] for p in pts],
            label=m,
            alpha=0.75,
            s=36,
            c=f"C{i}",
        )
    lims = [0.0, 1.05]
    ax.plot(lims, lims, "k--", alpha=0.35, label="y = x")
    ax.set_xlabel("Train accuracy")
    ax.set_ylabel("Validation accuracy")
    ax.set_title("Train vs validation accuracy (all trials)")
    ax.legend()
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.grid(True, alpha=0.3)
    fig.savefig(out_dir / "tuning_train_vs_val_accuracy.png", dpi=160)
    plt.close(fig)


def _plot_single_param_sweeps(results: List[Dict[str, Any]], out_dir: Path) -> None:
    """If exactly one numeric param varies (others fixed), draw a line plot vs validation accuracy."""
    _use_agg_backend()
    import matplotlib.pyplot as plt

    for model in sorted({r["model"] for r in results}):
        rows = [r for r in results if r["model"] == model]
        if len(rows) < 2:
            continue
        param_names = set()
        for r in rows:
            param_names.update(r["params"].keys())
        varying: List[str] = []
        for pname in sorted(param_names):
            vals = {json.dumps(r["params"].get(pname), sort_keys=True) for r in rows}
            if len(vals) > 1:
                varying.append(pname)
        if len(varying) != 1:
            continue
        pname = varying[0]
        numeric_ok = all(
            isinstance(r["params"].get(pname), (int, float)) and not isinstance(r["params"].get(pname), bool)
            for r in rows
        )
        if not numeric_ok:
            continue
        rows_s = sorted(rows, key=lambda r: float(r["params"][pname]))
        x = [float(r["params"][pname]) for r in rows_s]
        y = [r["val_accuracy"] for r in rows_s]
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
        ax.plot(x, y, "o-", color="C0")
        ax.set_xlabel(pname)
        ax.set_ylabel("Validation accuracy")
        ax.set_title(f"{model}: validation accuracy vs {pname}")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0.0, 1.02)
        safe = pname.replace("/", "_")
        fig.savefig(out_dir / f"tuning_{model.lower()}_{safe}_vs_val_accuracy.png", dpi=160)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Sklearn baseline hyperparameter tuning on cached PCA features.")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument(
        "--models",
        type=str,
        default="SVM,RF,XGBOOST",
        help="Comma-separated: SVM,RF,XGBOOST",
    )
    parser.add_argument("--grid-json", type=str, default=None, help="Optional JSON overriding default grids.")
    parser.add_argument(
        "--force-pca-rebuild",
        action="store_true",
        help="Pass force=True to ensure_sklearn_pca_cache (slow).",
    )
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--top-k-bars", type=int, default=12, help="Top trials per model in bar chart.")
    parser.add_argument("--label-max-len", type=int, default=44, help="Max chars in trial labels for bar charts.")
    args = parser.parse_args()

    nap = load_neural_ancestry_predictor_for_cli()
    load_config = nap.load_config
    get_dataset_cache_dir = nap.get_dataset_cache_dir
    prepare_data = nap.prepare_data
    validate_cache = nap.validate_cache

    config_path = Path(args.config)
    config = load_config(config_path)
    dataset_cache_dir = Path(get_dataset_cache_dir(config))
    if not dataset_cache_dir.exists() or not validate_cache(dataset_cache_dir, config):
        nap.console.print(
            f"[red]Cache do dataset inválido ou ausente: {dataset_cache_dir}[/red]\n"
            "[yellow]Gere o cache do dataset antes (treino neural_ancestry_predictor).[/yellow]"
        )
        sys.exit(1)

    models = [m.strip().upper() for m in args.models.split(",") if m.strip()]
    grids = copy.deepcopy(DEFAULT_GRIDS)
    if args.grid_json:
        custom = _load_grid_json(Path(args.grid_json))
        for m, g in custom.items():
            grids[m] = g

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    exp_dir = Path(config["dataset_input"]["processed_cache_dir"]) / "_sklearn_tune_workspace"
    exp_dir.mkdir(parents=True, exist_ok=True)
    _full, train_loader, val_loader, test_loader = prepare_data(config, exp_dir)

    nap.console.print(f"[cyan]Sklearn tuning | models={models} | output={out_dir.resolve()}[/cyan]")

    results, pca_dir = _run_trials(
        nap,
        config,
        dataset_cache_dir,
        train_loader,
        val_loader,
        test_loader,
        _full,
        models,
        grids,
        force_pca=args.force_pca_rebuild,
    )

    meta = {
        "config": str(config_path.resolve()),
        "pca_cache_dir": str(pca_dir.resolve()),
        "models": models,
        "n_trials": len(results),
    }
    with open(out_dir / "tuning_run_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    _save_results_json_csv(results, out_dir)

    if not args.no_plots and results:
        _plot_validation_summary(results, out_dir)
        _plot_top_trials_per_model(
            results,
            out_dir,
            top_k=args.top_k_bars,
            label_max_len=args.label_max_len,
        )
        _plot_ranked_val_accuracy_curves(results, out_dir)
        _plot_pair_heatmaps_if_two_params(results, out_dir)
        _plot_train_vs_val_accuracy(results, out_dir)
        _plot_single_param_sweeps(results, out_dir)
        nap.console.print(f"[green]Figuras guardadas em {out_dir.resolve()}[/green]")

    if results:
        best = max(results, key=lambda r: r["val_accuracy"])
        nap.console.print(
            f"[green]Melhor trial (validation accuracy): #{best['trial_id']} {best['model']} {best['label']} "
            f"val_acc={best['val_accuracy']:.4f} val_f1={best['val_f1']:.4f}[/green]"
        )
    nap.console.print("[green]Concluído.[/green]")


if __name__ == "__main__":
    main()
