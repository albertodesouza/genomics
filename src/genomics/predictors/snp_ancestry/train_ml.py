from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from genomics.predictors.snp_ancestry.markers import load_statistics
from genomics.predictors.snp_ancestry.pipeline import (
    _genotype_dose,
    _load_23andme_genotypes,
    load_config,
    load_configured_splits,
    sample_metadata,
)


SUPPORTED_MODELS = ("logistic", "random_forest")


def load_marker_ids(path: Path, top: Optional[int] = None) -> List[str]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "rsid" not in reader.fieldnames:
            raise ValueError(f"Marker TSV must contain an 'rsid' column: {path}")
        marker_ids = [row["rsid"] for row in reader if row.get("rsid")]
    if top is not None:
        marker_ids = marker_ids[:top]
    if not marker_ids:
        raise ValueError(f"No markers found in {path}")
    return marker_ids


def split_sample_ids(config: dict) -> Tuple[Dict[str, List[str]], Dict[str, dict], str]:
    splits = load_configured_splits(config)
    meta = sample_metadata(splits, config)
    level = config.get("prediction", {}).get("level", "superpopulation")
    by_split: Dict[str, List[str]] = {"train": [], "val": [], "test": []}
    for sid, item in meta.items():
        split_name = item.get("split")
        if split_name in by_split and item.get(level) is not None:
            by_split[split_name].append(sid)
    return by_split, meta, level


def build_feature_matrix(
    sample_ids: Sequence[str],
    marker_ids: Sequence[str],
    config: dict,
    stats: dict,
    meta: Dict[str, dict],
    level: str,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    inp = config["input"]
    pred_cfg = config.get("prediction", {})
    ind_dir = Path(inp["individuals_dir"])
    fname_tpl = config.get("conversion", {}).get("output_filename", "{sample_id}_23andme.txt")
    haplotype_mode = pred_cfg.get("haplotype_mode", "H1+H2")
    ref_alleles = stats.get("ref_alleles", {})
    needed = set(marker_ids)
    X = np.full((len(sample_ids), len(marker_ids)), np.nan, dtype=np.float32)
    y: List[str] = []
    kept_samples: List[str] = []

    kept_idx = 0
    for sid in sample_ids:
        path = ind_dir / sid / fname_tpl.format(sample_id=sid)
        if not path.exists():
            continue
        genotypes = _load_23andme_genotypes(str(path), needed)
        for col_idx, rsid in enumerate(marker_ids):
            genotype = genotypes.get(rsid)
            if genotype is None:
                continue
            dose = _genotype_dose(genotype, ref_alleles.get(rsid, ""), haplotype_mode)
            if dose is not None:
                X[kept_idx, col_idx] = float(dose)
        y.append(str(meta[sid][level]))
        kept_samples.append(sid)
        kept_idx += 1

    if not kept_samples:
        raise ValueError("No samples with existing 23andMe files were found")
    kept_rows = len(kept_samples)
    return X[:kept_rows], np.array(y, dtype=object), kept_samples


def impute_with_train_mean(X_train: np.ndarray, *others: np.ndarray) -> Tuple[np.ndarray, List[np.ndarray], np.ndarray]:
    means = np.nanmean(X_train, axis=0)
    means = np.where(np.isnan(means), 0.0, means).astype(np.float32)
    return _impute(X_train, means), [_impute(arr, means) for arr in others], means


def _impute(arr: np.ndarray, means: np.ndarray) -> np.ndarray:
    out = np.array(arr, dtype=np.float32, copy=True)
    missing = np.isnan(out)
    if missing.any():
        out[missing] = np.take(means, np.where(missing)[1])
    return out


def _make_classifier(model_name: str, random_seed: int):
    if model_name == "logistic":
        from sklearn.linear_model import LogisticRegression

        return LogisticRegression(max_iter=2000, multi_class="auto", class_weight="balanced", random_state=random_seed)
    if model_name == "random_forest":
        from sklearn.ensemble import RandomForestClassifier

        return RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced",
            random_state=random_seed,
            n_jobs=-1,
        )
    raise ValueError(f"Unsupported model '{model_name}'. Choose one of: {', '.join(SUPPORTED_MODELS)}")


def evaluate_predictions(y_true: Sequence[str], y_pred: Sequence[str], labels: Sequence[str]) -> dict:
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, confusion_matrix, f1_score

    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, labels=list(labels), average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_true, y_pred, labels=list(labels), average="weighted", zero_division=0)),
        "confusion_matrix": confusion_matrix(y_true, y_pred, labels=list(labels)).tolist(),
        "classification_report": classification_report(y_true, y_pred, labels=list(labels), zero_division=0),
    }


def _feature_importance(model, model_name: str, marker_ids: Sequence[str]) -> List[dict]:
    if model_name == "random_forest" and hasattr(model, "feature_importances_"):
        scores = np.asarray(model.feature_importances_, dtype=float)
    elif model_name == "logistic" and hasattr(model, "coef_"):
        scores = np.mean(np.abs(np.asarray(model.coef_, dtype=float)), axis=0)
    else:
        scores = np.zeros(len(marker_ids), dtype=float)
    order = np.argsort(-scores)
    return [
        {"rank": int(rank), "rsid": marker_ids[int(idx)], "importance": float(scores[int(idx)])}
        for rank, idx in enumerate(order, start=1)
    ]


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def _write_predictions(path: Path, sample_ids: Sequence[str], y_true: Sequence[str], y_pred: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "true", "predicted", "correct"], delimiter="\t")
        writer.writeheader()
        for sid, true, pred in zip(sample_ids, y_true, y_pred):
            writer.writerow({"sample_id": sid, "true": true, "predicted": pred, "correct": str(true == pred).lower()})


def _write_importance(path: Path, rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "rsid", "importance"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_train_ml(args: argparse.Namespace) -> int:
    try:
        import joblib
        from sklearn.preprocessing import LabelEncoder
    except ImportError as exc:
        raise SystemExit(
            "Missing scikit-learn/joblib dependencies. Install with `python3 -m pip install -e \".[snp-ancestry]\"` or `python3 -m pip install scikit-learn joblib`."
        ) from exc

    config = load_config(args.config)
    stats, statistics_path = load_statistics(args.config, args.statistics)
    marker_ids = load_marker_ids(args.markers, args.top_markers)
    available = set(stats.get("allele_frequencies", {}))
    marker_ids = [rsid for rsid in marker_ids if rsid in available]
    if not marker_ids:
        raise SystemExit("No marker IDs overlap the statistics file")

    by_split, meta, level = split_sample_ids(config)
    train_ids = by_split["train"]
    val_ids = by_split["val"]
    test_ids = by_split["test"]
    if not train_ids:
        raise SystemExit("No train samples found in configured splits")

    X_train_raw, y_train, train_kept = build_feature_matrix(train_ids, marker_ids, config, stats, meta, level)
    X_val_raw, y_val, val_kept = build_feature_matrix(val_ids, marker_ids, config, stats, meta, level) if val_ids else (np.empty((0, len(marker_ids))), np.array([], dtype=object), [])
    X_test_raw, y_test, test_kept = build_feature_matrix(test_ids, marker_ids, config, stats, meta, level) if test_ids else (np.empty((0, len(marker_ids))), np.array([], dtype=object), [])
    X_train, (X_val, X_test), impute_means = impute_with_train_mean(X_train_raw, X_val_raw, X_test_raw)

    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)
    labels = list(label_encoder.classes_)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "metadata": {
            "config": str(args.config),
            "statistics": str(statistics_path),
            "markers": str(args.markers),
            "level": level,
            "n_markers": len(marker_ids),
            "labels": labels,
            "created_at": datetime.now().isoformat(),
        },
        "models": {},
    }

    for model_name in args.models:
        model = _make_classifier(model_name, args.random_seed)
        model.fit(X_train, y_train_encoded)
        model_dir = output_dir / model_name
        model_dir.mkdir(parents=True, exist_ok=True)

        split_payload = {}
        for split_name, X_split, y_split, kept in (
            ("train", X_train, y_train, train_kept),
            ("val", X_val, y_val, val_kept),
            ("test", X_test, y_test, test_kept),
        ):
            if len(y_split) == 0:
                continue
            pred_encoded = model.predict(X_split)
            y_pred = label_encoder.inverse_transform(pred_encoded)
            metrics = evaluate_predictions(y_split, y_pred, labels)
            metrics["n_samples"] = int(len(y_split))
            split_payload[split_name] = metrics
            _write_predictions(model_dir / f"predictions_{split_name}.tsv", kept, y_split, y_pred)

        importance = _feature_importance(model, model_name, marker_ids)
        _write_importance(model_dir / "feature_importance.tsv", importance)
        joblib.dump(
            {
                "model": model,
                "label_encoder": label_encoder,
                "marker_ids": marker_ids,
                "impute_means": impute_means,
                "level": level,
                "config": str(args.config),
                "statistics": str(statistics_path),
            },
            model_dir / "model.joblib",
        )
        model_summary = {
            "metrics": split_payload,
            "artifact": str(model_dir / "model.joblib"),
            "feature_importance": str(model_dir / "feature_importance.tsv"),
        }
        _write_json(model_dir / "metrics.json", model_summary)
        summary["models"][model_name] = model_summary
        test_metrics = split_payload.get("test") or split_payload.get("val") or split_payload.get("train")
        print(
            f"{model_name}: accuracy={test_metrics['accuracy']:.4f} "
            f"balanced_accuracy={test_metrics['balanced_accuracy']:.4f} "
            f"macro_f1={test_metrics['macro_f1']:.4f}"
        )

    _write_json(output_dir / "summary.json", summary)
    print(f"Output: {output_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train sklearn baselines on SNP ancestry AIMs")
    parser.add_argument("--config", type=Path, required=True, help="SNP ancestry YAML config used for splits and genotype files")
    parser.add_argument("--markers", type=Path, required=True, help="AIM TSV generated by `genomics snp-ancestry markers`")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--statistics", type=Path, default=None, help="Explicit statistics JSON path; defaults to the config-derived path")
    parser.add_argument("--models", nargs="+", choices=SUPPORTED_MODELS, default=["logistic", "random_forest"])
    parser.add_argument("--top-markers", type=int, default=None, help="Use only the first N rows from the marker TSV")
    parser.add_argument("--random-seed", type=int, default=13)
    parser.set_defaults(func=run_train_ml)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
