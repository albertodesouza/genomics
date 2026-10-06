#!/usr/bin/env python3
"""Single-gene arms vs the eleven-gene panel: information carried vs reliance taken.

Re-analysis only. Reads the two single-gene knockdown replays, the two single-gene run
directories' test metrics, and the published eleven-gene raw_center_crop knockdown.

The comparison the arms are for:

    single-gene test accuracy      what the window alone can do -- the INFORMATION it carries
    single-gene |Delta|            the same 100 bp scramble with reliance at its maximum
    panel |Delta|                  what the eleven-gene classifier actually took from it
    panel |Delta| / single |Delta| the fraction of the available response the panel model kept

A gene whose single-gene accuracy is high but whose panel |Delta| is a small fraction of its
single-gene |Delta| is informative but REDUNDANT -- the panel model routed around it. That is a
different statement from "the classifier does not use this gene", and the panel measurement alone
cannot make it.

Usage:
  python3 scripts/experiments/single_gene_report.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home/breno/I2CA/genomics")
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
PUBLISHED = KO_DIR / "pigmentation_test_split_knockout_no_alignment.csv"
NULL = KO_DIR / "pigmentation_test_split_null_scramble.csv"
OUT = REPO_ROOT / "results/genotype_based_predictor/single_gene_arms.json"
RUN_NAME = "cnn2_pigmentation_rna_seq_H1+H2_raw_center_crop_32768_log_s1k6x32f16_s2f32_s3f64_gpavg_fc256_L100-40_relu_0.5_adam"
ARMS = {"slc24a5": "SLC24A5", "tyr": "TYR"}
METHOD = "biology_tss"


def p_weak(strong, weak):
    s = np.vstack([np.asarray(strong, float), np.asarray(weak, float)])
    s = s - s.max(axis=0, keepdims=True)
    e = np.exp(s)
    return e[1] / e.sum(axis=0)


def load_delta(path):
    df = pd.read_csv(path)
    df["p_weak_base"] = p_weak(df["baseline_strong_logit"], df["baseline_weak_logit"])
    df["p_weak_pert"] = p_weak(df["perturbed_strong_logit"], df["perturbed_weak_logit"])
    df["delta"] = df["p_weak_pert"] - df["p_weak_base"]
    return df


def main():
    pub = load_delta(PUBLISHED)
    pub = pub[pub["method"] == METHOD]
    pub_g = pub.groupby("gene")["delta"].agg(["mean", "std"])
    null = load_delta(NULL)
    pooled_p99_11 = float(null["delta"].abs().quantile(0.99))

    rows = []
    print("=" * 100)
    print("SINGLE-GENE ARMS -- one gene in the input, nothing to route around")
    print("=" * 100)
    for arm, gene in ARMS.items():
        ko_path = KO_DIR / f"pigmentation_test_split_knockout_single_{arm}.csv"
        run_dir = REPO_ROOT / f"results/genotype_based_predictor/runs_single_gene_{arm}" / RUN_NAME
        test_json = run_dir / "test_best_accuracy_results.json"
        if not ko_path.exists():
            print(f"  {arm}: no knockdown replay yet ({ko_path.name}) -- skipping")
            continue
        ko = load_delta(ko_path)
        ko_m = ko[ko["method"] == METHOD]
        acc = None
        cm = None
        if test_json.exists():
            tj = json.loads(test_json.read_text())
            acc = float(tj["weighted_accuracy"])
            cm = tj["confusion_matrix"]

        single_delta = float(ko_m["delta"].mean())
        single_abs = abs(single_delta)
        single_in = float(ko_m["delta_in"].mean())
        flip = float(ko_m["flipped"].mean())
        acc_base = float(ko_m["correct_baseline"].mean())
        acc_pert = float(ko_m["correct_perturbed"].mean())
        panel_delta = float(pub_g.loc[gene, "mean"])
        panel_abs = abs(panel_delta)

        rec = {
            "arm": arm, "gene": gene,
            "single_gene_test_accuracy": acc, "single_gene_confusion_matrix": cm,
            "single_gene_baseline_accuracy_on_replayed_rows": acc_base,
            "single_gene_accuracy_after_knockdown": acc_pert,
            "single_gene_delta": single_delta, "single_gene_abs_delta": single_abs,
            "single_gene_delta_in": single_in,
            "single_gene_ratio": single_abs / single_in if single_in else None,
            "single_gene_flip_rate": flip,
            "panel_delta": panel_delta, "panel_abs_delta": panel_abs,
            "panel_flip_rate": float(pub[pub["gene"] == gene]["flipped"].mean()),
            "retained_fraction": panel_abs / single_abs if single_abs else None,
            "panel_clears_null_p99": panel_abs > pooled_p99_11,
        }
        # per-method breakdown, since all three anchorings were replayed
        rec["by_method"] = {
            m: {"delta": float(sub["delta"].mean()), "abs_delta": float(abs(sub["delta"].mean())),
                "delta_in": float(sub["delta_in"].mean()), "flip_rate": float(sub["flipped"].mean()),
                "n": int(len(sub))}
            for m, sub in ko.groupby("method")}
        rows.append(rec)

        print()
        print(f"--- {gene} (arm {arm}) ---")
        print(f"  single-gene test accuracy      : {acc:.4f}" if acc is not None else "  single-gene test accuracy      : n/a")
        if cm:
            print(f"  confusion matrix               : {cm}")
        print(f"  accuracy on the 162 replayed   : {acc_base:.4f} baseline -> {acc_pert:.4f} after knockdown")
        print(f"  single-gene Delta (P(weak))    : {single_delta:+.4f}   |Delta| = {single_abs:.4f}")
        print(f"  single-gene |Delta_in|         : {single_in:.1f}    ratio = {single_abs/single_in:.3e}")
        print(f"  single-gene flip rate          : {flip:.3f}")
        print(f"  panel Delta (11-gene model)    : {panel_delta:+.4f}   |Delta| = {panel_abs:.4f}  "
              f"flip = {rec['panel_flip_rate']:.3f}")
        print(f"  retained by the panel model    : {rec['retained_fraction']:.3f} of the single-gene response")
        print("  by anchoring method:")
        for m, v in rec["by_method"].items():
            print(f"      {m:<18} Delta = {v['delta']:+.4f}  |Delta_in| = {v['delta_in']:>9.1f}  flip = {v['flip_rate']:.3f}")

    if len(rows) == 2:
        a, b = rows
        print()
        print("=" * 100)
        print(f"  {a['gene']} vs {b['gene']}")
        print(f"    information (single-gene accuracy) : {a['single_gene_test_accuracy']:.4f} vs {b['single_gene_test_accuracy']:.4f}")
        print(f"    reliance at maximum (single |Delta|): {a['single_gene_abs_delta']:.4f} vs {b['single_gene_abs_delta']:.4f}")
        print(f"    reliance in the panel (panel |Delta|): {a['panel_abs_delta']:.4f} vs {b['panel_abs_delta']:.4f}")
        print(f"    retained fraction                   : {a['retained_fraction']:.3f} vs {b['retained_fraction']:.3f}")

    OUT.write_text(json.dumps({"method": METHOD, "pooled_null_p99_11gene": pooled_p99_11,
                               "arms": rows}, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
