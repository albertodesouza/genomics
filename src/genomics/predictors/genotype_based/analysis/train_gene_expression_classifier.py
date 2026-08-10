"""Train a linear classifier on scalar "total gene expression" features.

Each individual is reduced to one number per (gene, ontology_term, strand): the mean
AlphaGenome RNA-seq signal across that gene's window, combined across both haplotypes. This is
the "poor man's PrediXcan" TWAS arm -- classical TWAS methods associate one predicted-expression
scalar per gene with phenotype, rather than a full positional signal -- used to check how much of
CNN2AncestryPredictor's performance survives when the model only ever sees gene-level expression
totals instead of positional detail. See
notebooks/pigmentation_alphagenome_cnn_variant_scoring.ipynb Section 9 for how this is compared
against the full-tensor CNN/NN/sklearn arms and the snp_ancestry GWAS-style baselines.

Reuses the exact same dataset_input/output config and processed-tensor cache as the other
pigmentation model arms (only `dataset_input`/`output` are read; `model`/`training` are ignored
since this script does its own lightweight LogisticRegression fit), so results are directly
comparable.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

import numpy as np

from genomics.core.metrics import auc_score, classification_metrics, save_results_json
from genomics.predictors.genotype_based.analysis.gene_expression_scalar import (
    materialize_scalar_split,
    scalar_feature_names,
)
from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
from genomics.predictors.genotype_based.data.pipeline import (
    _make_runtime_processed_datasets,
    _resolve_runtime_dataset_dir,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train a LogisticRegression on scalar mean-RNA-seq gene-expression features"
    )
    parser.add_argument("config_path", type=str, help="genotype_based YAML config (dataset_input/output are used; model/training are ignored)")
    parser.add_argument("--reduction", choices=["mean", "sum"], default="mean", help="How to combine per-position signal within a gene window (default: mean, i.e. average RNA-seq signal)")
    parser.add_argument("--haplotype-mode", choices=["separate", "sum"], default="separate", help="'separate' keeps H1/H2 as two features per track (default, matching the full-tensor CNN/NN arms); 'sum' merges them into one diploid-total feature per track (assumes each haplotype's contribution adds independently)")
    parser.add_argument("--C", type=float, default=1.0, help="LogisticRegression inverse regularization strength")
    parser.add_argument("--out-dir", type=str, default=None, help="Defaults to <results_dir>/scalar_gene_expression_logreg_<reduction>_<haplotype_mode>")
    args = parser.parse_args()

    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)

    runtime_dataset_dir = _resolve_runtime_dataset_dir(config)
    cache_dir = get_dataset_cache_dir(config)
    print(f"Runtime dataset dir: {runtime_dataset_dir}")
    print(f"Cache dir          : {cache_dir}")

    full_ds, train_ds, val_ds, test_ds = _make_runtime_processed_datasets(runtime_dataset_dir, cache_dir, config)

    feature_names = scalar_feature_names(config, haplotype_mode=args.haplotype_mode)
    print(f"Reducing to {len(feature_names)} scalar features per individual (reduction={args.reduction}, haplotype_mode={args.haplotype_mode})...")

    X_train, y_train, train_ids = materialize_scalar_split(train_ds, config, reduction=args.reduction, haplotype_mode=args.haplotype_mode)
    X_val, y_val, val_ids = materialize_scalar_split(val_ds, config, reduction=args.reduction, haplotype_mode=args.haplotype_mode)
    X_test, y_test, test_ids = materialize_scalar_split(test_ds, config, reduction=args.reduction, haplotype_mode=args.haplotype_mode)
    print(f"train={len(y_train)} val={len(y_val)} test={len(y_test)}")

    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)

    seed = config.data_split.random_seed if config.data_split.random_seed not in (None, -1) else 42
    model = LogisticRegression(C=args.C, max_iter=5000, class_weight="balanced", random_state=seed)
    model.fit(X_train_s, y_train)

    idx_to_target = full_ds.idx_to_target
    class_list: List[int] = model.classes_.tolist()
    class_names = [idx_to_target[c] for c in class_list]
    strong_label = next((idx for idx, name in idx_to_target.items() if name == "strong pigmentation"), None)

    out_dir = Path(args.out_dir) if args.out_dir else Path(config.dataset_input.results_dir) / f"scalar_gene_expression_logreg_{args.reduction}_{args.haplotype_mode}"
    out_dir.mkdir(parents=True, exist_ok=True)

    for split_name, X_raw, y, sample_ids in (
        ("train", X_train, y_train, train_ids),
        ("val", X_val, y_val, val_ids),
        ("test", X_test, y_test, test_ids),
    ):
        if len(y) == 0:
            continue
        X_s = scaler.transform(X_raw)
        y_pred = model.predict(X_s)
        metrics = classification_metrics(y, y_pred, class_names)
        if strong_label is not None and strong_label in class_list:
            strong_col = class_list.index(strong_label)
            y_score = model.predict_proba(X_s)[:, strong_col]
            y_bin = (np.asarray(y) == strong_label).astype(int)
            metrics["auc_vs_strong_pigmentation"] = auc_score(y_bin, y_score)
        save_results_json(metrics, out_dir / f"{split_name}_results.json")
        auc_note = ""
        if "auc_vs_strong_pigmentation" in metrics:
            auc_note = f" auc_vs_strong={metrics['auc_vs_strong_pigmentation']['auc']:.4f}"
        print(f"{split_name}: accuracy={metrics['accuracy']:.4f}{auc_note} n={len(y)}")

    coef = model.coef_[0] if model.coef_.shape[0] == 1 else np.mean(np.abs(model.coef_), axis=0)
    ranked = sorted(zip(feature_names, coef.tolist()), key=lambda kv: -abs(kv[1]))
    with open(out_dir / "feature_importance.json", "w", encoding="utf-8") as f:
        json.dump([{"feature": name, "coefficient": value} for name, value in ranked], f, indent=2)

    import joblib

    joblib.dump(
        {
            "scaler": scaler,
            "model": model,
            "feature_names": feature_names,
            "idx_to_target": idx_to_target,
            "strong_label": strong_label,
            "reduction": args.reduction,
            "haplotype_mode": args.haplotype_mode,
        },
        out_dir / "model.joblib",
    )

    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
