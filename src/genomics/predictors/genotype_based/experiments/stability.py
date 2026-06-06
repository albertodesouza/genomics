from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from rich.console import Console
from torch.utils.data import ConcatDataset, DataLoader, Subset

from genomics.core.metrics import save_results_json
from genomics.core.sklearn_pca_cache import (
    compute_sklearn_pca_effective_k,
    fit_incremental_pca_on_train,
    fit_standard_scaler_incremental,
    fit_streaming_randomized_pca_on_train,
    stack_scaled_pca_batches,
)
from genomics.predictors.genotype_based.config import (
    PipelineConfig,
    get_experiment_runs_dir,
    load_config,
)
from genomics.predictors.genotype_based.data.pipeline import _make_data_loaders, prepare_data


console = Console()
SKLEARN_BASELINE_TYPES = frozenset({"SVM", "LOGREG", "RF", "XGBOOST"})


class _IndexedDataset:
    """Wraps a dataset and returns a stable local index for PCA routines."""

    def __init__(self, dataset: Any):
        self.dataset = dataset
        for name in ("target_to_idx", "idx_to_target", "config"):
            if hasattr(dataset, name):
                setattr(self, name, getattr(dataset, name))

    def __len__(self) -> int:
        return len(self.dataset)

    def __getitem__(self, idx: int):
        features, target, _old_idx = self.dataset[idx]
        return features, target, int(idx)

    def get_num_classes(self) -> int:
        if hasattr(self.dataset, "get_num_classes"):
            return self.dataset.get_num_classes()
        return len(getattr(self, "idx_to_target", {}))

    def get_input_shape(self):
        if hasattr(self.dataset, "get_input_shape"):
            return self.dataset.get_input_shape()
        features, _target, _idx = self[0]
        shape = features.shape
        if len(shape) == 3:
            return (shape[0] * shape[1], shape[2])
        if len(shape) == 2:
            return (shape[0], shape[1])
        return (1, shape[0])


def _output_dir(config: PipelineConfig) -> Path:
    st = config.stability_analysis
    if st.output_dir:
        return Path(st.output_dir)
    run_name = st.run_name or f"stability_{config.model.type.lower()}"
    return get_experiment_runs_dir(config) / run_name


def _targets_for_dataset(dataset: Any) -> np.ndarray:
    targets = []
    for idx in range(len(dataset)):
        _features, target, _item_idx = dataset[idx]
        targets.append(int(target.item() if hasattr(target, "item") else target))
    return np.asarray(targets, dtype=np.int64)


def _split_plan(config: PipelineConfig, y_dev: np.ndarray) -> List[Dict[str, Any]]:
    st = config.stability_analysis
    seed = int(st.random_seed)
    indices = np.arange(len(y_dev))
    val_fraction = st.val_split if st.val_split is not None else config.data_split.val_split / max(config.data_split.train_split + config.data_split.val_split, 1e-12)
    strategy = "repeated_random_split" if st.strategy == "randomized_split" else st.strategy

    if strategy == "repeated_random_split":
        if st.stratify:
            from sklearn.model_selection import StratifiedShuffleSplit

            splitter = StratifiedShuffleSplit(n_splits=st.n_repeats, test_size=val_fraction, random_state=seed)
            pairs = splitter.split(indices, y_dev)
        else:
            from sklearn.model_selection import ShuffleSplit

            splitter = ShuffleSplit(n_splits=st.n_repeats, test_size=val_fraction, random_state=seed)
            pairs = splitter.split(indices)
        return [
            {"name": f"repeat_{i:03d}", "kind": strategy, "seed": seed + i - 1, "train_idx": train_idx.tolist(), "val_idx": val_idx.tolist()}
            for i, (train_idx, val_idx) in enumerate(pairs, start=1)
        ]

    if strategy == "cross_validation":
        if st.stratify:
            from sklearn.model_selection import StratifiedKFold

            splitter = StratifiedKFold(n_splits=st.n_splits, shuffle=True, random_state=seed)
            pairs = splitter.split(indices, y_dev)
        else:
            from sklearn.model_selection import KFold

            splitter = KFold(n_splits=st.n_splits, shuffle=True, random_state=seed)
            pairs = splitter.split(indices)
        return [
            {"name": f"fold_{i:03d}", "kind": strategy, "seed": seed, "train_idx": train_idx.tolist(), "val_idx": val_idx.tolist()}
            for i, (train_idx, val_idx) in enumerate(pairs, start=1)
        ]

    raise ValueError(f"stability_analysis.strategy nao suportada: {st.strategy}")


def _loader_for(dataset: Any, indices: Sequence[int], config: PipelineConfig, *, train: bool) -> DataLoader:
    train_loader, val_loader, _test_loader = _make_data_loaders(
        Subset(dataset, list(indices)),
        Subset(dataset, list(indices)),
        Subset(dataset, []),
        config,
    )
    return train_loader if train else val_loader


def _fit_transform_split(
    config: PipelineConfig,
    train_loader: DataLoader,
    val_loader: DataLoader,
    output_dir: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    first = next(iter(train_loader))
    n_features = int(np.prod(first[0].shape[1:]))
    n_train = len(train_loader.dataset)
    sk = config.model.sklearn
    k, _pca_req = compute_sklearn_pca_effective_k(config.model_dump(), n_train=n_train, n_features=n_features, log=console.print)
    scaler = fit_standard_scaler_incremental(train_loader, rich_console=console)
    output_dir.mkdir(parents=True, exist_ok=True)
    if sk.pca_backend == "randomized_streaming":
        pca = fit_streaming_randomized_pca_on_train(
            train_loader,
            scaler,
            k,
            output_dir,
            oversampling=sk.randomized_pca_oversampling,
            n_iter=sk.randomized_pca_n_iter,
            feature_chunk_size=sk.randomized_pca_feature_chunk_size,
            dtype=sk.randomized_pca_dtype,
            random_state=config.data_split.random_seed or config.stability_analysis.random_seed,
            log=console.print,
            rich_console=console,
        )
    elif sk.pca_backend == "incremental":
        pca = fit_incremental_pca_on_train(
            train_loader,
            scaler,
            k,
            log=console.print,
            forbid_tail_padding=sk.pca_align_n_train,
            rich_console=console,
        )
    else:
        raise ValueError(f"pca_backend nao suportado: {sk.pca_backend}")
    X_train, y_train = stack_scaled_pca_batches(train_loader, scaler, pca, rich_console=console, progress_desc="Stability PCA: transform train")
    X_val, y_val = stack_scaled_pca_batches(val_loader, scaler, pca, rich_console=console, progress_desc="Stability PCA: transform val")
    return X_train, y_train, X_val, y_val, int(k)


def _summarize(rows: List[Dict[str, Any]], metric_names: Sequence[str]) -> Dict[str, Any]:
    summary = {}
    for metric in metric_names:
        values = np.asarray([float(row[f"val_{metric}"]) for row in rows], dtype=np.float64)
        summary[metric] = {
            "mean": float(values.mean()) if values.size else 0.0,
            "std": float(values.std(ddof=1)) if values.size > 1 else 0.0,
            "min": float(values.min()) if values.size else 0.0,
            "max": float(values.max()) if values.size else 0.0,
        }
    return summary


def _save_outputs(output_dir: Path, config: PipelineConfig, rows: List[Dict[str, Any]], plan: List[Dict[str, Any]]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    metric_names = ("accuracy", "precision", "recall", "f1")
    summary = {
        "strategy": config.stability_analysis.strategy,
        "test_split_policy": "fixed_from_base_split_not_used_for_resampling",
        "selection_metric": config.stability_analysis.selection_metric,
        "num_runs": len(rows),
        "metrics": _summarize(rows, metric_names),
        "runs": rows,
        "split_plan": [
            {key: value for key, value in item.items() if key not in ("train_idx", "val_idx")}
            for item in plan
        ],
    }
    save_results_json(summary, output_dir / "stability_results.json", console)
    with open(output_dir / "stability_results.csv", "w", newline="", encoding="utf-8") as f:
        fieldnames = ["run", "kind", "seed", "n_train", "n_val", "pca_components", "val_accuracy", "val_precision", "val_recall", "val_f1"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})
    with open(output_dir / "split_plan.json", "w", encoding="utf-8") as f:
        json.dump(plan, f, indent=2)


def run_stability(config_path: Path) -> Path:
    from genomics.predictors.genotype_based.models.sklearn_models import (
        build_sklearn_classifier,
        sklearn_metrics_dict,
    )

    config = load_config(config_path)
    if not config.stability_analysis.enabled:
        raise ValueError("stability_analysis.enabled deve ser true para genotype stability")
    if config.model.type.upper() not in SKLEARN_BASELINE_TYPES:
        raise ValueError("genotype stability suporta inicialmente apenas modelos sklearn: SVM, LOGREG, RF, XGBOOST")
    config.evaluation.confidence_intervals.enabled = False
    output_dir = _output_dir(config).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    prepare_dir = output_dir / "prepare_cache"
    prepare_dir.mkdir(parents=True, exist_ok=True)
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, prepare_dir)
    train_ds = train_loader.dataset
    val_ds = val_loader.dataset
    test_ds = test_loader.dataset
    dev_ds = _IndexedDataset(ConcatDataset([train_ds, val_ds]))
    y_dev = _targets_for_dataset(dev_ds)
    plan = _split_plan(config, y_dev)
    console.print(
        f"[green]Stability dev={len(dev_ds)} test_fixo={len(test_ds)} runs={len(plan)} "
        f"strategy={config.stability_analysis.strategy}[/green]"
    )

    rows: List[Dict[str, Any]] = []
    for idx, split in enumerate(plan, start=1):
        name = split["name"]
        console.print(f"[cyan]Run {idx}/{len(plan)}: {name}[/cyan]")
        train_loader = _loader_for(dev_ds, split["train_idx"], config, train=True)
        val_loader = _loader_for(dev_ds, split["val_idx"], config, train=False)
        X_train, y_train, X_val, y_val, k = _fit_transform_split(config, train_loader, val_loader, output_dir / "pca" / name)
        valid = y_train >= 0
        X_train, y_train = X_train[valid], y_train[valid]
        clf = build_sklearn_classifier(config, config.model.type, int(split.get("seed", config.stability_analysis.random_seed)), n_train=len(y_train))
        clf.fit(X_train, y_train)
        results = sklearn_metrics_dict(y_val, clf.predict(X_val), full_ds, config=None)
        row = {
            "run": name,
            "kind": split["kind"],
            "seed": split.get("seed"),
            "n_train": len(split["train_idx"]),
            "n_val": len(split["val_idx"]),
            "pca_components": k,
            "val_accuracy": float(results.get("accuracy", 0.0)),
            "val_precision": float(results.get("precision", 0.0)),
            "val_recall": float(results.get("recall", 0.0)),
            "val_f1": float(results.get("f1", 0.0)),
        }
        rows.append(row)
        save_results_json(results, output_dir / "runs" / name / "val_results.json", console)
    _save_outputs(output_dir, config, rows, plan)
    console.print(f"[bold green]✓ Stability results salvos em: {output_dir}[/bold green]")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Genotype stability analysis with fixed test split")
    parser.add_argument("config_path", type=Path)
    args = parser.parse_args()
    run_stability(args.config_path)


if __name__ == "__main__":
    main()
