from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.predictors.snp_ancestry.markers import load_statistics
from genomics.predictors.snp_ancestry.train_ml import (
    SUPPORTED_MODELS,
    _make_classifier,
    build_feature_matrix,
    evaluate_predictions,
    impute_with_train_mean,
    load_marker_ids,
    split_sample_ids,
)


DEFAULT_REMOVE_TOP = (0, 1, 5, 10, 50, 100)


def ablated_marker_ids(marker_ids: Sequence[str], remove_top: int) -> List[str]:
    if remove_top < 0:
        raise ValueError("remove_top values must be >= 0")
    return list(marker_ids[remove_top:])


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def _write_tsv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "model",
        "remove_top",
        "n_markers",
        "metric_split",
        "accuracy",
        "balanced_accuracy",
        "macro_f1",
        "weighted_f1",
        "n_samples",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _choose_metric_split(split_metrics: Dict[str, dict], preferred: str) -> Tuple[str, dict]:
    if preferred != "auto":
        if preferred not in split_metrics:
            raise ValueError(f"Requested metric split '{preferred}' has no samples")
        return preferred, split_metrics[preferred]
    for split_name in ("test", "val", "train"):
        if split_name in split_metrics:
            return split_name, split_metrics[split_name]
    raise ValueError("No split metrics were produced")


def _prepare_arrays(config: dict, stats: dict, marker_ids: Sequence[str]) -> Tuple[dict, List[str]]:
    by_split, meta, level = split_sample_ids(config)
    if not by_split["train"]:
        raise ValueError("No train samples found in configured splits")

    arrays = {}
    for split_name in ("train", "val", "test"):
        sample_ids = by_split[split_name]
        if sample_ids:
            X, y, kept = build_feature_matrix(sample_ids, marker_ids, config, stats, meta, level)
        else:
            X, y, kept = np.empty((0, len(marker_ids))), np.array([], dtype=object), []
        arrays[split_name] = {"X": X, "y": y, "sample_ids": kept}
    return arrays, level


def run_ablation(args: argparse.Namespace) -> int:
    try:
        from sklearn.preprocessing import LabelEncoder
    except ImportError as exc:
        raise SystemExit(
            "Missing scikit-learn dependency. Install with `python3 -m pip install -e \".[snp-ancestry]\"` or `python3 -m pip install scikit-learn`."
        ) from exc

    from genomics.predictors.snp_ancestry.pipeline import load_config

    config = load_config(args.config)
    stats, statistics_path = load_statistics(args.config, args.statistics)
    marker_ids = load_marker_ids(args.markers, args.top_markers)
    available = set(stats.get("allele_frequencies", {}))
    marker_ids = [rsid for rsid in marker_ids if rsid in available]
    if not marker_ids:
        raise SystemExit("No marker IDs overlap the statistics file")

    remove_top_values = sorted(set(args.remove_top))
    arrays, level = _prepare_arrays(config, stats, marker_ids)
    label_encoder = LabelEncoder()
    y_train = arrays["train"]["y"]
    y_train_encoded = label_encoder.fit_transform(y_train)
    labels = list(label_encoder.classes_)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    summary = {
        "metadata": {
            "config": str(args.config),
            "statistics": str(statistics_path),
            "markers": str(args.markers),
            "level": level,
            "base_marker_count": len(marker_ids),
            "models": list(args.models),
            "remove_top": remove_top_values,
            "labels": labels,
            "created_at": datetime.now().isoformat(),
        },
        "runs": [],
    }

    for remove_top in remove_top_values:
        kept_marker_ids = ablated_marker_ids(marker_ids, remove_top)
        if not kept_marker_ids:
            print(f"remove_top={remove_top}: skipped because no markers remain")
            continue
        cols = np.arange(remove_top, len(marker_ids))
        X_train_raw = arrays["train"]["X"][:, cols]
        X_val_raw = arrays["val"]["X"][:, cols]
        X_test_raw = arrays["test"]["X"][:, cols]
        X_train, (X_val, X_test), _ = impute_with_train_mean(X_train_raw, X_val_raw, X_test_raw)

        for model_name in args.models:
            model = _make_classifier(model_name, args.random_seed)
            model.fit(X_train, y_train_encoded)
            split_metrics = {}
            for split_name, X_split in (("train", X_train), ("val", X_val), ("test", X_test)):
                y_split = arrays[split_name]["y"]
                if len(y_split) == 0:
                    continue
                pred_encoded = model.predict(X_split)
                y_pred = label_encoder.inverse_transform(pred_encoded)
                metrics = evaluate_predictions(y_split, y_pred, labels)
                metrics["n_samples"] = int(len(y_split))
                split_metrics[split_name] = metrics

            metric_split, metrics = _choose_metric_split(split_metrics, args.metric_split)
            row = {
                "model": model_name,
                "remove_top": int(remove_top),
                "n_markers": len(kept_marker_ids),
                "metric_split": metric_split,
                "accuracy": float(metrics["accuracy"]),
                "balanced_accuracy": float(metrics["balanced_accuracy"]),
                "macro_f1": float(metrics["macro_f1"]),
                "weighted_f1": float(metrics["weighted_f1"]),
                "n_samples": int(metrics["n_samples"]),
            }
            rows.append(row)
            summary["runs"].append({**row, "metrics": split_metrics})
            print(
                f"{model_name} remove_top={remove_top} n_markers={len(kept_marker_ids)} "
                f"{metric_split}_accuracy={metrics['accuracy']:.4f} "
                f"{metric_split}_balanced_accuracy={metrics['balanced_accuracy']:.4f}"
            )

    _write_tsv(args.output_dir / "ablation.tsv", rows)
    _write_json(args.output_dir / "summary.json", summary)
    print(f"Output: {args.output_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ablate top SNP ancestry AIMs and retrain sklearn baselines")
    parser.add_argument("--config", type=Path, required=True, help="SNP ancestry YAML config used for splits and genotype files")
    parser.add_argument("--markers", type=Path, required=True, help="AIM TSV generated by `genomics snp-ancestry markers`")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--statistics", type=Path, default=None, help="Explicit statistics JSON path; defaults to the config-derived path")
    parser.add_argument("--models", nargs="+", choices=SUPPORTED_MODELS, default=["logistic"])
    parser.add_argument("--remove-top", nargs="+", type=int, default=list(DEFAULT_REMOVE_TOP), help="Top marker counts to remove before retraining")
    parser.add_argument("--top-markers", type=int, default=None, help="Use only the first N rows from the marker TSV before ablation")
    parser.add_argument("--metric-split", choices=["auto", "train", "val", "test"], default="auto")
    parser.add_argument("--random-seed", type=int, default=13)
    parser.set_defaults(func=run_ablation)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
