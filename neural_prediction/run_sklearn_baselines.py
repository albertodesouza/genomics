#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
from sklearn.svm import SVC

from .config import load_config
from .sklearn_pca_cache import METADATA_FILENAME, ensure_sklearn_pca_cache


def metrics_dict(y_true, y_pred) -> Dict[str, float]:
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, average="weighted", zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, average="weighted", zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, average="weighted", zero_division=0)),
    }


def build_classifier(model_name: str, params: Dict[str, Any]):
    if model_name == "SVM":
        return SVC(C=float(params.get("C", 1.0)), probability=bool(params.get("calibrate_probabilities", False)), kernel="rbf")
    if model_name == "RF":
        return RandomForestClassifier(
            n_estimators=int(params.get("n_estimators", 200)),
            max_depth=params.get("max_depth", None),
            random_state=13,
            n_jobs=-1,
        )
    if model_name == "XGBOOST":
        import xgboost as xgb

        return xgb.XGBClassifier(
            n_estimators=int(params.get("n_estimators", 200)),
            max_depth=int(params.get("max_depth", 6)),
            learning_rate=float(params.get("learning_rate", 0.1)),
            random_state=13,
            n_jobs=-1,
            eval_metric="mlogloss",
        )
    raise ValueError(f"Unsupported sklearn baseline: {model_name}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run sklearn baselines on neural_prediction PCA cache")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--force-pca", action="store_true")
    parser.add_argument("--models", nargs="*", default=["SVM", "RF", "XGBOOST"])
    args = parser.parse_args()

    config = load_config(Path(args.config))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pca_dir = ensure_sklearn_pca_cache(config, force=args.force_pca)
    X_train = np.load(pca_dir / "X_train.npy")
    y_train = np.load(pca_dir / "y_train.npy")
    X_val = np.load(pca_dir / "X_val.npy")
    y_val = np.load(pca_dir / "y_val.npy")
    metadata = json.loads((pca_dir / METADATA_FILENAME).read_text(encoding="utf-8"))

    grids = {
        "SVM": [{"C": c, "calibrate_probabilities": cal} for c in [0.1, 1.0, 10.0] for cal in [False, True]],
        "RF": [{"n_estimators": n, "max_depth": d} for n in [100, 200] for d in [None, 10, 20]],
        "XGBOOST": [{"n_estimators": n, "max_depth": d, "learning_rate": lr} for n in [100, 200] for d in [4, 6] for lr in [0.05, 0.1]],
    }

    results: List[Dict[str, Any]] = []
    trial_id = 0
    for model_name in args.models:
        if model_name == "XGBOOST":
            try:
                import xgboost  # noqa: F401
            except ImportError:
                continue
        for params in grids[model_name]:
            trial_id += 1
            classifier = build_classifier(model_name, params)
            classifier.fit(X_train, y_train)
            y_pred_train = classifier.predict(X_train)
            y_pred_val = classifier.predict(X_val)
            train_metrics = metrics_dict(y_train, y_pred_train)
            val_metrics = metrics_dict(y_val, y_pred_val)
            results.append(
                {
                    "trial_id": trial_id,
                    "model": model_name,
                    "pca_k": metadata["pca_components_effective"],
                    "params_json": json.dumps(params, sort_keys=True),
                    "train_accuracy": train_metrics["accuracy"],
                    "train_precision": train_metrics["precision"],
                    "train_recall": train_metrics["recall"],
                    "train_f1": train_metrics["f1"],
                    "val_accuracy": val_metrics["accuracy"],
                    "val_precision": val_metrics["precision"],
                    "val_recall": val_metrics["recall"],
                    "val_f1": val_metrics["f1"],
                }
            )

    with (output_dir / "tuning_results.json").open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
        handle.write("\n")

    if results:
        with (output_dir / "tuning_results.csv").open("w", encoding="utf-8", newline="") as handle:
            fieldnames = list(results[0].keys())
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
