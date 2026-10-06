from __future__ import annotations

import argparse
import csv
import itertools
import json
import shutil
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np
from rich.console import Console
from rich.table import Table

from genomics.core import update_manifest
from genomics.core.metrics import save_classification_plots, save_results_json
from genomics.predictors.genotype_based.config import PipelineConfig, generate_dataset_name, load_config
from genomics.predictors.genotype_based.config import get_dataset_cache_dir, get_experiment_runs_dir


console = Console()


def _metric_value(results: Dict[str, Any], short_key: str) -> float:
    weighted_key = {
        "accuracy": "weighted_accuracy",
        "precision": "weighted_precision",
        "recall": "weighted_recall",
        "f1": "weighted_f1_score",
    }.get(short_key)
    value = results.get(weighted_key, results.get(short_key, 0.0)) if weighted_key else results.get(short_key, 0.0)
    return float(value)


def _search_dir(config: PipelineConfig) -> Path:
    search = config.hyperparameter_search
    if search.output_dir:
        return Path(search.output_dir)
    run_name = search.run_name or f"search_{generate_dataset_name(config)}"
    return get_experiment_runs_dir(config) / run_name


def _rf_candidates(config: PipelineConfig) -> Iterable[Dict[str, Any]]:
    rf = config.hyperparameter_search.random_forest
    if not rf.enabled:
        return []
    return (
        {
            "model_type": "RF",
            "params": {
                "n_estimators": n_estimators,
                "max_depth": max_depth,
                "class_weight": class_weight,
                "random_state": rf.random_state,
                "n_jobs": rf.n_jobs,
            },
        }
        for n_estimators, max_depth, class_weight in itertools.product(rf.n_estimators, rf.max_depth, rf.class_weight)
    )


def _xgboost_candidates(config: PipelineConfig) -> Iterable[Dict[str, Any]]:
    xgb = config.hyperparameter_search.xgboost
    if not xgb.enabled:
        return []
    return (
        {
            "model_type": "XGBOOST",
            "params": {
                "n_estimators": n_estimators,
                "max_depth": max_depth,
                "learning_rate": learning_rate,
                "subsample": subsample,
                "colsample_bytree": colsample_bytree,
                "tree_method": tree_method,
                "random_state": xgb.random_state,
                "n_jobs": xgb.n_jobs,
                "eval_metric": xgb.eval_metric,
            },
        }
        for n_estimators, max_depth, learning_rate, subsample, colsample_bytree, tree_method in itertools.product(
            xgb.n_estimators,
            xgb.max_depth,
            xgb.learning_rate,
            xgb.subsample,
            xgb.colsample_bytree,
            xgb.tree_method,
        )
    )


def _candidate_config(config: PipelineConfig, candidate: Dict[str, Any]) -> PipelineConfig:
    candidate_config = config.model_copy(deep=True)
    model_type = candidate["model_type"]
    params = candidate["params"]
    candidate_config.model.type = model_type
    if model_type == "RF":
        rf = candidate_config.model.sklearn.random_forest
        rf.n_estimators = int(params["n_estimators"])
        rf.max_depth = params["max_depth"]
        rf.class_weight = params["class_weight"]
        rf.random_state = int(params["random_state"])
        rf.n_jobs = int(params["n_jobs"])
    elif model_type == "XGBOOST":
        xgb = candidate_config.model.sklearn.xgboost
        xgb.n_estimators = int(params["n_estimators"])
        xgb.max_depth = int(params["max_depth"])
        xgb.learning_rate = float(params["learning_rate"])
        xgb.subsample = float(params["subsample"])
        xgb.colsample_bytree = float(params["colsample_bytree"])
        xgb.tree_method = str(params["tree_method"])
        xgb.random_state = int(params["random_state"])
        xgb.n_jobs = int(params["n_jobs"])
        xgb.eval_metric = str(params["eval_metric"])
    return candidate_config


def _train_candidate(
    *,
    config: PipelineConfig,
    candidate: Dict[str, Any],
    pca_dir: Path,
    val_dataset: Any,
    full_dataset: Any,
    search_dir: Path,
    rank_index: int,
) -> Dict[str, Any]:
    import joblib

    from genomics.predictors.genotype_based.models.sklearn_models import (
        SKLEARN_ARTIFACT_FILENAME,
        build_sklearn_classifier,
        run_sklearn_eval_and_save,
        sklearn_metrics_dict,
    )

    candidate_config = _candidate_config(config, candidate)
    model_type = candidate["model_type"]
    params = candidate["params"]
    X_train = np.load(pca_dir / "X_train.npy")
    y_train = np.load(pca_dir / "y_train.npy")
    valid = y_train >= 0
    X_train, y_train = X_train[valid], y_train[valid]
    clf = build_sklearn_classifier(candidate_config, model_type, int(config.data_split.random_seed or 42), n_train=len(y_train))
    clf.fit(X_train, y_train)

    X_val = np.load(pca_dir / "X_val.npy")
    y_val = np.load(pca_dir / "y_val.npy")
    results = sklearn_metrics_dict(y_val, clf.predict(X_val), full_dataset)
    candidate_name = f"{rank_index:03d}_{model_type.lower()}"
    candidate_dir = search_dir / "candidates" / candidate_name
    models_dir = candidate_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = models_dir / SKLEARN_ARTIFACT_FILENAME
    joblib.dump(
        {
            "classifier": clf,
            "model_type": model_type,
            "pca_n_components": int(_read_pca_k(pca_dir)),
            "pca_cache_dir": str(pca_dir.resolve()),
            "dataset_cache_dir": str(get_dataset_cache_dir(config).resolve()),
            "search_params": params,
        },
        artifact_path,
    )
    run_sklearn_eval_and_save(results, candidate_dir, "val", wandb_run=None, split_name="val")
    row = {
        "candidate": candidate_name,
        "model_type": model_type,
        "params": params,
        "val_weighted_accuracy": _metric_value(results, "accuracy"),
        "val_weighted_precision": _metric_value(results, "precision"),
        "val_weighted_recall": _metric_value(results, "recall"),
        "val_weighted_f1_score": _metric_value(results, "f1"),
        "artifact_path": str(artifact_path),
        "candidate_dir": str(candidate_dir),
    }
    return row


def _read_pca_k(pca_dir: Path) -> int:
    from genomics.core.sklearn_pca_cache import METADATA_FILENAME as SKLEARN_PCA_METADATA_FILENAME

    with open(pca_dir / SKLEARN_PCA_METADATA_FILENAME, "r", encoding="utf-8") as f:
        return int(json.load(f)["pca_n_components_effective"])


def _save_search_outputs(rows: List[Dict[str, Any]], best: Dict[str, Any], search_dir: Path, metric: str) -> None:
    search_dir.mkdir(parents=True, exist_ok=True)
    with open(search_dir / "search_results.json", "w", encoding="utf-8") as f:
        json.dump({"selection_metric": metric, "best": best, "results": rows}, f, indent=2, ensure_ascii=False)
    with open(search_dir / "search_results.csv", "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "rank",
            "candidate",
            "model_type",
            "val_weighted_accuracy",
            "val_weighted_precision",
            "val_weighted_recall",
            "val_weighted_f1_score",
            "params",
            "artifact_path",
        ]
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            csv_row = {key: row.get(key) for key in fieldnames}
            csv_row["rank"] = rank
            csv_row["params"] = json.dumps(row["params"], sort_keys=True)
            writer.writerow(csv_row)
    save_results_json(best, search_dir / "best_summary.json", console)
    _save_search_plot(rows, search_dir / "plots")


def _search_plot_labels(rows: List[Dict[str, Any]]) -> List[str]:
    return [f"{row['candidate']}\n{row['model_type']}" for row in rows]


def _save_search_plot(rows: List[Dict[str, Any]], plots_dir: Path) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError:
        console.print("[yellow]matplotlib indisponível; plot da busca não foi gerado.[/yellow]")
        return
    labels = _search_plot_labels(rows)
    values = [row["val_weighted_accuracy"] for row in rows]
    fig, ax = plt.subplots(figsize=(max(10, len(rows) * 0.35), 5))
    bars = ax.bar(labels, values)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Validation accuracy")
    ax.set_title("Hyperparameter search validation accuracy")
    ax.tick_params(axis="x", rotation=90)
    for bar, value in zip(bars, values):
        ax.annotate(f"{value:.3f}", xy=(bar.get_x() + bar.get_width() / 2, value), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    fig.savefig(plots_dir / "validation_accuracy_by_candidate.png", dpi=160)
    plt.close(fig)


def _copy_best_artifact(best: Dict[str, Any], search_dir: Path) -> Path:
    from genomics.predictors.genotype_based.models.sklearn_models import SKLEARN_ARTIFACT_FILENAME

    best_dir = search_dir / "best"
    models_dir = best_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(best["artifact_path"], models_dir / SKLEARN_ARTIFACT_FILENAME)
    shutil.copytree(Path(best["candidate_dir"]) / "plots", best_dir / "plots", dirs_exist_ok=True)
    for name in ("val_results.json",):
        src = Path(best["candidate_dir"]) / name
        if src.exists():
            shutil.copy2(src, best_dir / name)
    return best_dir


def run_search(config_path: Path) -> Path:
    from genomics.core.sklearn_pca_cache import ensure_sklearn_pca_cache
    from genomics.predictors.genotype_based.data.pipeline import prepare_data

    config = load_config(config_path)
    if not config.hyperparameter_search.enabled:
        raise ValueError("hyperparameter_search.enabled deve ser true para genotype search")
    search_dir = _search_dir(config).resolve()
    search_dir.mkdir(parents=True, exist_ok=True)
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, search_dir)
    pca_dir = ensure_sklearn_pca_cache(
        config.model_dump(),
        get_dataset_cache_dir(config),
        train_loader,
        val_loader,
        test_loader,
        force=config.debug.force_pca_cache_rebuild,
        log=console.print,
        rich_console=console,
    )
    candidates = list(_rf_candidates(config)) + list(_xgboost_candidates(config))
    if not candidates:
        raise ValueError("Nenhum candidato habilitado em hyperparameter_search")
    rows = []
    for idx, candidate in enumerate(candidates, start=1):
        console.print(f"[cyan]Treinando candidato {idx}/{len(candidates)}:[/cyan] {candidate['model_type']} {candidate['params']}")
        rows.append(
            _train_candidate(
                config=config,
                candidate=candidate,
                pca_dir=pca_dir,
                val_dataset=val_loader.dataset,
                full_dataset=full_ds,
                search_dir=search_dir,
                rank_index=idx,
            )
        )
    metric_key = f"val_weighted_{config.hyperparameter_search.selection_metric}"
    if config.hyperparameter_search.selection_metric == "f1":
        metric_key = "val_weighted_f1_score"
    rows.sort(key=lambda row: row.get(metric_key, 0.0), reverse=True)
    best = rows[0]
    best_dir = _copy_best_artifact(best, search_dir)
    _save_search_outputs(rows, {**best, "best_dir": str(best_dir)}, search_dir, metric_key)
    update_manifest(search_dir, status="completed", best_dir=str(best_dir), selection_metric=metric_key)

    table = Table(title="Hyperparameter search results", show_header=True)
    table.add_column("Rank", justify="right")
    table.add_column("Model")
    table.add_column("Val Acc", justify="right")
    table.add_column("Val F1", justify="right")
    table.add_column("Params")
    for rank, row in enumerate(rows, start=1):
        table.add_row(str(rank), row["model_type"], f"{row['val_weighted_accuracy']:.4f}", f"{row['val_weighted_f1_score']:.4f}", json.dumps(row["params"], sort_keys=True))
    console.print(table)
    console.print(f"[bold green]Melhor modelo:[/bold green] {best['model_type']} → {best_dir}")
    return search_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Grid search for genotype sklearn baselines using validation accuracy")
    parser.add_argument("config_path", type=Path)
    args = parser.parse_args()
    run_search(args.config_path.resolve())


if __name__ == "__main__":
    main()
