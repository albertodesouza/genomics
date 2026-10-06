#!/usr/bin/env python3
"""Rank the fourteen genes -- eleven pigmentation, three irrelevant -- on one classifier.

Re-analysis only: reads the CSVs written by specificity_14gene_knockdown.py and
specificity_14gene_null.py. No AlphaGenome calls, no GPU, no training.

Statistics are defined exactly as in scripts/experiments/readout_normalisation.py, so the numbers
are comparable to the published eleven-gene ones in kind (they are NOT comparable in value: this
is a different classifier, and the paper says so):

    Delta        mean change in P(weak) over the 162 test individuals, signed
    |Delta|      its absolute value -- the per-gene reliance statistic
    |Delta_in|   delivered perturbation reaching the CNN input, summed over haplotypes
    ratio        |Delta| / |Delta_in| -- response per unit delivered
    clears_null  |Delta| above the 99th percentile of the pooled null-scramble |delta|,
                 measured on THIS checkpoint

The pre-registered reading, fixed in the config before the retrain:
    all three controls flat with |Delta_in| in the responding range -> specificity established;
    any control responding at panel-gene magnitude -> the probe responds to delivered
    perturbation at any well-expressed locus; mixed -> report per gene, never average.

Usage:
  python3 scripts/experiments/specificity_14gene_report.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home/breno/I2CA/genomics")
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
KNOCKOUT = KO_DIR / "pigmentation_test_split_knockout_14gene.csv"
NULL = KO_DIR / "pigmentation_test_split_null_14gene.csv"
PUBLISHED = KO_DIR / "pigmentation_test_split_knockout_no_alignment.csv"
OUT = REPO_ROOT / "results/genotype_based_predictor/specificity_14gene_ranking.json"
PANEL = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"]
CONTROLS = ["TPM2", "SMCR8", "PSMC4"]


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
    ko = load_delta(KNOCKOUT)
    expected_rows = ko["sample_id"].nunique() * ko["gene"].nunique()
    if len(ko) != expected_rows:
        raise SystemExit(f"ABORT: knockdown CSV has {len(ko)} rows, expected {expected_rows}. "
                         "The run is incomplete -- do not read a partial ranking.")
    # The null is written incrementally by a long run. Reading a partial file silently produces a
    # p99 computed from whichever genes happened to finish first, and every gene then "clears" it.
    have_null = False
    null = None
    if NULL.exists():
        null = load_delta(NULL)
        n_null = len(null)
        if n_null < expected_rows:
            print(f"WARNING: null CSV has {n_null}/{expected_rows} rows and is still being written. "
                  "clears_null is SUPPRESSED until it completes.\n")
            null = None
        else:
            have_null = True

    g = ko.groupby("gene")
    per = g["delta"].agg(["mean", "std", "count"]).rename(
        columns={"mean": "delta", "std": "delta_sd", "count": "n"})
    per["abs_delta"] = per["delta"].abs()
    per["mean_abs_delta"] = g["delta"].apply(lambda s: s.abs().mean())
    per["delta_in"] = g["delta_in"].mean()
    per["flip_rate"] = g["flipped"].mean()
    per["gene_class"] = g["gene_class"].first()
    per["ratio"] = per["abs_delta"] / per["delta_in"]

    if have_null:
        pooled_p99 = float(null["delta"].abs().quantile(0.99))
        per = per.join(null.groupby("gene")["delta"].apply(lambda s: s.abs().mean()).rename("null_mean_abs"))
        per["clears_null"] = per["abs_delta"] > pooled_p99
        per["null_delta_in"] = null.groupby("gene")["delta_in"].mean()
    else:
        pooled_p99 = None
        per["null_mean_abs"] = np.nan
        per["clears_null"] = np.nan

    per = per.sort_values("abs_delta", ascending=False)
    per["rank_abs_delta"] = np.arange(1, len(per) + 1)
    per["rank_ratio"] = per["ratio"].rank(ascending=False).astype(int)

    print("=" * 92)
    print("FOURTEEN-GENE RANKING   (11 pigmentation panel genes + 3 pigmentation-irrelevant controls)")
    print(f"classifier: runs_specificity_14gene, method biology_tss, n = {int(per['n'].max())} test individuals")
    if have_null:
        print(f"pooled null |delta| p99 on this checkpoint = {pooled_p99:.5f}")
    else:
        print("NO NULL AVAILABLE YET -- clears_null is undefined.")
    print("=" * 92)
    hdr = f"{'#':>2} {'gene':<9}{'class':<9}{'Delta':>10}{'|Delta|':>10}{'|Delta_in|':>12}{'ratio':>11}{'flip':>7}{'null':>9}{'clears':>7}"
    print(hdr)
    for i, (gene, r) in enumerate(per.iterrows(), 1):
        cn = "" if not have_null else ("yes" if r["clears_null"] else "no")
        nm = "" if not have_null else f"{r['null_mean_abs']:.5f}"
        print(f"{i:>2} {gene:<9}{r['gene_class']:<9}{r['delta']:>10.4f}{r['abs_delta']:>10.4f}"
              f"{r['delta_in']:>12.1f}{r['ratio']:>11.3e}{r['flip_rate']:>7.3f}{nm:>9}{cn:>7}")

    # ---- does the panel/control split show up in the ranking? --------------------
    ctrl = per[per["gene_class"] == "control"]
    panel = per[per["gene_class"] == "panel"]
    best_ctrl_rank = int(ctrl["rank_abs_delta"].min())
    n_panel_below = int((panel["rank_abs_delta"] > best_ctrl_rank).sum())
    u, pu = stats.mannwhitneyu(panel["abs_delta"], ctrl["abs_delta"], alternative="greater")
    print()
    print(f"best control by |Delta|: {ctrl['rank_abs_delta'].idxmin()} at rank {best_ctrl_rank}/14; "
          f"{n_panel_below} panel genes rank below it")
    print(f"Mann-Whitney |Delta| panel > control: U = {u:.1f}, p = {pu:.4f}  (n = 11 vs 3)")
    ur, pur = stats.mannwhitneyu(panel["ratio"], ctrl["ratio"], alternative="greater")
    print(f"Mann-Whitney ratio  panel > control: U = {ur:.1f}, p = {pur:.4f}")
    print()
    print("delivered magnitude check -- the controls were selected to bracket TYR/MFSD12/MC1R:")
    for c in CONTROLS:
        if c not in per.index:
            continue
        d = per.loc[c, "delta_in"]
        lo = panel[panel["delta_in"] <= d].sort_values("delta_in").tail(1)
        hi = panel[panel["delta_in"] >= d].sort_values("delta_in").head(1)
        span = (f"{lo.index[0]} {lo['delta_in'].iloc[0]:.0f}" if len(lo) else "--") + "  <  " + \
               (f"{hi.index[0]} {hi['delta_in'].iloc[0]:.0f}" if len(hi) else "--")
        print(f"  {c:<7} |Delta_in| = {d:>9.0f}   between panel genes: {span}")

    # ---- comparison with the published 11-gene ranking (rank agreement only) -----
    comp = None
    if PUBLISHED.exists():
        pub = load_delta(PUBLISHED)
        pub = pub[pub["method"] == "biology_tss"]
        pub_g = pub.groupby("gene")["delta"].mean().abs()
        common = [x for x in per.index if x in pub_g.index]
        rho, prho = stats.spearmanr(per.loc[common, "abs_delta"], pub_g.loc[common])
        signs = int((np.sign(per.loc[common, "delta"]) == np.sign(
            pub.groupby("gene")["delta"].mean().loc[common])).sum())
        comp = {"spearman_rho": float(rho), "p": float(prho), "n_genes": len(common),
                "sign_agreement": f"{signs}/{len(common)}"}
        print()
        print(f"vs the published 11-gene classifier (different model, ranking comparison only): "
              f"Spearman rho = {rho:.3f} (p = {prho:.4f}), signs agree {signs}/{len(common)}")

    payload = {
        "classifier": "runs_specificity_14gene",
        "method": "biology_tss",
        "n_individuals": int(per["n"].max()),
        "pooled_null_p99": pooled_p99,
        "per_gene": [
            {"gene": gene, "gene_class": r["gene_class"], "rank": int(r["rank_abs_delta"]),
             "delta": float(r["delta"]), "delta_sd": float(r["delta_sd"]),
             "abs_delta": float(r["abs_delta"]), "mean_abs_delta": float(r["mean_abs_delta"]),
             "delta_in": float(r["delta_in"]), "ratio": float(r["ratio"]),
             "rank_ratio": int(r["rank_ratio"]), "flip_rate": float(r["flip_rate"]),
             "null_mean_abs": (None if not have_null else float(r["null_mean_abs"])),
             "clears_null": (None if not have_null else bool(r["clears_null"]))}
            for gene, r in per.iterrows()],
        "panel_vs_control": {
            "best_control": ctrl["rank_abs_delta"].idxmin(), "best_control_rank": best_ctrl_rank,
            "panel_genes_ranked_below_best_control": n_panel_below,
            "mannwhitney_abs_delta": {"U": float(u), "p": float(pu)},
            "mannwhitney_ratio": {"U": float(ur), "p": float(pur)}},
        "vs_published_11gene": comp,
    }
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
