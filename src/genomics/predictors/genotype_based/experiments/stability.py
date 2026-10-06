from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

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
    _stratified_fit_indices,
)
from genomics.predictors.genotype_based.config import (
    PipelineConfig,
    get_experiment_runs_dir,
    load_config,
)
from genomics.predictors.genotype_based.data.pipeline import _collate_fn, prepare_data


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

    def get_sample_id(self, idx: int) -> str:
        dataset = self.dataset
        if hasattr(dataset, "datasets") and hasattr(dataset, "cumulative_sizes"):
            previous = 0
            for child, cumulative in zip(dataset.datasets, dataset.cumulative_sizes):
                if idx < cumulative:
                    local_idx = idx - previous
                    if hasattr(child, "get_sample_id"):
                        return str(child.get_sample_id(local_idx))
                    break
                previous = cumulative
        if hasattr(dataset, "get_sample_id"):
            return str(dataset.get_sample_id(idx))
        return f"#{idx + 1}"


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


def _sample_ids_for_dataset(dataset: Any) -> List[str]:
    sample_ids = []
    for idx in range(len(dataset)):
        if hasattr(dataset, "get_sample_id"):
            sample_ids.append(str(dataset.get_sample_id(idx)))
            continue
        metadata = dataset.get_sample_metadata(idx) if hasattr(dataset, "get_sample_metadata") else None
        if metadata and metadata.get("sample_id"):
            sample_ids.append(str(metadata["sample_id"]))
        else:
            sample_ids.append(f"#{idx + 1}")
    return sample_ids


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


def _split_plan_path(config: PipelineConfig, output_dir: Path) -> Path:
    value = config.stability_analysis.split_plan_path
    if value:
        return Path(value)
    return output_dir / "split_plan.json"


def _plan_with_sample_ids(
    config: PipelineConfig,
    plan: List[Dict[str, Any]],
    dev_sample_ids: Sequence[str],
    test_sample_ids: Sequence[str],
) -> Dict[str, Any]:
    return {
        "format_version": 1,
        "strategy": config.stability_analysis.strategy,
        "random_seed": config.stability_analysis.random_seed,
        "stratify": config.stability_analysis.stratify,
        "n_repeats": config.stability_analysis.n_repeats,
        "n_splits": config.stability_analysis.n_splits,
        "val_split": config.stability_analysis.val_split,
        "test_split_policy": "fixed_from_base_split_not_used_for_resampling",
        "dev_sample_ids": list(dev_sample_ids),
        "test_sample_ids": list(test_sample_ids),
        "runs": [
            {
                "name": item["name"],
                "kind": item["kind"],
                "seed": item.get("seed"),
                "train_sample_ids": [dev_sample_ids[idx] for idx in item["train_idx"]],
                "val_sample_ids": [dev_sample_ids[idx] for idx in item["val_idx"]],
            }
            for item in plan
        ],
    }


def _plan_from_sample_ids(payload: Dict[str, Any], dev_sample_ids: Sequence[str]) -> List[Dict[str, Any]]:
    if int(payload.get("format_version", 0)) != 1:
        raise ValueError("split_plan_path tem format_version incompatível")
    sample_to_idx = {sample_id: idx for idx, sample_id in enumerate(dev_sample_ids)}
    plan = []
    for item in payload.get("runs", []):
        missing = [sid for sid in item.get("train_sample_ids", []) + item.get("val_sample_ids", []) if sid not in sample_to_idx]
        if missing:
            raise ValueError(f"split_plan_path contém sample_id ausente no dev atual: {missing[:5]}")
        plan.append(
            {
                "name": item["name"],
                "kind": item.get("kind", payload.get("strategy", "repeated_random_split")),
                "seed": item.get("seed"),
                "train_idx": [sample_to_idx[sid] for sid in item.get("train_sample_ids", [])],
                "val_idx": [sample_to_idx[sid] for sid in item.get("val_sample_ids", [])],
            }
        )
    if not plan:
        raise ValueError("split_plan_path não contém runs")
    return plan


def _load_or_create_split_plan(
    config: PipelineConfig,
    output_dir: Path,
    y_dev: np.ndarray,
    dev_sample_ids: Sequence[str],
    test_sample_ids: Sequence[str],
) -> List[Dict[str, Any]]:
    path = _split_plan_path(config, output_dir)
    if config.stability_analysis.reuse_split_plan and path.exists():
        with open(path, "r", encoding="utf-8") as f:
            payload = json.load(f)
        console.print(f"[green]Split plan compartilhado carregado:[/green] {path}")
        return _plan_from_sample_ids(payload, dev_sample_ids)

    plan = _split_plan(config, y_dev)
    payload = _plan_with_sample_ids(config, plan, dev_sample_ids, test_sample_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    console.print(f"[green]Split plan compartilhado salvo:[/green] {path}")
    return plan


def _loader_for(dataset: Any, indices: Sequence[int], config: PipelineConfig, *, train: bool) -> DataLoader:
    # Stability repeatedly scans cached tensors for PCA/scaler fitting. Sequential
    # access avoids random shard hopping and keeps lazy shard memory bounded.
    return DataLoader(
        Subset(dataset, list(indices)),
        batch_size=config.training.batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=_collate_fn,
    )


def _pca_fit_loader_from_labels(
    train_loader: DataLoader,
    labels: Optional[np.ndarray],
    config: PipelineConfig,
) -> DataLoader:
    sk = config.model.sklearn
    fraction = float(sk.pca_fit_sample_fraction)
    min_samples = int(sk.pca_fit_min_samples or 0)
    min_per_class = int(sk.pca_fit_min_samples_per_class or 0)
    if labels is None or (fraction >= 1.0 and min_samples <= 0 and min_per_class <= 0):
        return train_loader

    fit_indices = _stratified_fit_indices(
        np.asarray(labels, dtype=np.int64),
        fraction=fraction,
        min_samples=min_samples,
        min_samples_per_class=min_per_class,
        random_seed=sk.pca_fit_random_seed,
        stratify=sk.pca_fit_stratify,
    )
    if len(fit_indices) >= len(train_loader.dataset):
        return train_loader
    console.print(
        f"[cyan]PCA fit subset: {len(fit_indices)}/{len(train_loader.dataset)} "
        f"amostras (fraction={fraction}, stratify={sk.pca_fit_stratify})[/cyan]"
    )
    return DataLoader(
        Subset(train_loader.dataset, fit_indices),
        batch_size=train_loader.batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=train_loader.collate_fn,
    )


def _fit_transform_split(
    config: PipelineConfig,
    train_loader: DataLoader,
    val_loader: DataLoader,
    output_dir: Path,
    train_labels: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    t0 = time.perf_counter()
    console.print(
        f"[cyan]Stability PCA: preparando fit loader "
        f"(train={len(train_loader.dataset)} val={len(val_loader.dataset)})[/cyan]"
    )
    fit_loader = _pca_fit_loader_from_labels(train_loader, train_labels, config)
    console.print(f"[cyan]Stability PCA: lendo primeiro batch de {len(fit_loader.dataset)} amostras de fit...[/cyan]")
    first = next(iter(fit_loader))
    n_features = int(np.prod(first[0].shape[1:]))
    n_fit = len(fit_loader.dataset)
    sk = config.model.sklearn
    k, _pca_req = compute_sklearn_pca_effective_k(config.model_dump(), n_train=n_fit, n_features=n_features, log=console.print)
    console.print(f"[cyan]Stability PCA: StandardScaler n_fit={n_fit} n_features={n_features} k={k}[/cyan]")
    scaler = fit_standard_scaler_incremental(fit_loader, rich_console=console)
    output_dir.mkdir(parents=True, exist_ok=True)
    if sk.pca_backend == "randomized_streaming":
        console.print(f"[cyan]Stability PCA: randomized_streaming em {output_dir}[/cyan]")
        pca = fit_streaming_randomized_pca_on_train(
            fit_loader,
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
        console.print(f"[cyan]Stability PCA: incremental k={k}[/cyan]")
        pca = fit_incremental_pca_on_train(
            fit_loader,
            scaler,
            k,
            log=console.print,
            forbid_tail_padding=sk.pca_align_n_train,
            rich_console=console,
        )
    else:
        raise ValueError(f"pca_backend nao suportado: {sk.pca_backend}")
    console.print("[cyan]Stability PCA: transformando train completo...[/cyan]")
    X_train, y_train = stack_scaled_pca_batches(train_loader, scaler, pca, rich_console=console, progress_desc="Stability PCA: transform train")
    console.print("[cyan]Stability PCA: transformando val...[/cyan]")
    X_val, y_val = stack_scaled_pca_batches(val_loader, scaler, pca, rich_console=console, progress_desc="Stability PCA: transform val")
    console.print(
        f"[green]✓ Stability PCA pronto em {time.perf_counter() - t0:.1f}s: "
        f"X_train={X_train.shape} X_val={X_val.shape}[/green]"
    )
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
        "split_plan_path": str(_split_plan_path(config, output_dir)),
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
    with open(output_dir / "split_plan_indices.json", "w", encoding="utf-8") as f:
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
    dev_sample_ids = _sample_ids_for_dataset(dev_ds)
    test_sample_ids = _sample_ids_for_dataset(test_ds)
    plan = _load_or_create_split_plan(config, output_dir, y_dev, dev_sample_ids, test_sample_ids)
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
        train_labels = y_dev[np.asarray(split["train_idx"], dtype=np.int64)]
        X_train, y_train, X_val, y_val, k = _fit_transform_split(
            config,
            train_loader,
            val_loader,
            output_dir / "pca" / name,
            train_labels=train_labels,
        )
        valid = y_train >= 0
        X_train, y_train = X_train[valid], y_train[valid]
        clf = build_sklearn_classifier(config, config.model.type, int(split.get("seed", config.stability_analysis.random_seed)), n_train=len(y_train))
        fit_t0 = time.perf_counter()
        console.print(
            f"[cyan]Treinando {config.model.type.upper()} stability {name}: "
            f"{len(y_train)} amostras x {X_train.shape[1]} componentes PCA...[/cyan]"
        )
        clf.fit(X_train, y_train)
        console.print(f"[green]✓ {config.model.type.upper()} stability {name} treinado em {time.perf_counter() - fit_t0:.1f}s[/green]")
        console.print(f"[cyan]Avaliando validation {name} ({len(y_val)} amostras)...[/cyan]")
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
        console.print(
            f"[green]✓ {name}: "
            f"acc={row['val_accuracy']:.4f} precision={row['val_precision']:.4f} "
            f"recall={row['val_recall']:.4f} f1={row['val_f1']:.4f}[/green]"
        )
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
