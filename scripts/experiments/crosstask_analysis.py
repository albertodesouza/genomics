#!/usr/bin/env python3
"""E16 analysis -- compare the pigmentation and superpopulation reliance rankings.

Reads the two replays and asks whether the probe returns the same gene ordering when the only
thing that changed is the label the probed classifier was fitted to.

Reported for the full 162 and, separately, for the 52 pigmentation test individuals that are NOT
in the superpopulation model's training split -- the superpopulation model saw 110 of them during
training, and although a perturbation readout is not an accuracy measurement, the held-out subset
is the version of the comparison that cannot be objected to.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home/breno/I2CA/genomics")
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
PIG = KO_DIR / "pigmentation_test_split_knockout.csv"
XT = KO_DIR / "pigmentation_test_split_knockout_superpopulation.csv"
OUT = REPO_ROOT / "results/genotype_based_predictor/crosstask_transfer.json"
METHOD = "biology_tss"


def p_weak(df):
    s = np.vstack([df["baseline_strong_logit"], df["baseline_weak_logit"]])
    s = s - s.max(0, keepdims=True); e = np.exp(s); b = e[1] / e.sum(0)
    s = np.vstack([df["perturbed_strong_logit"], df["perturbed_weak_logit"]])
    s = s - s.max(0, keepdims=True); e = np.exp(s); p = e[1] / e.sum(0)
    return p - b


def report(pig, xt, label):
    pg = pig.groupby("gene")["delta"].mean().abs().rename("pig_absdelta")
    # class-agnostic magnitude for the 5-class model, plus the drop in its own baseline class
    xg = xt.groupby("gene").agg(xt_tv=("tv_distance", "mean"),
                                xt_dbase=("delta_baseline_class", "mean"),
                                xt_flip=("flipped", "mean"))
    j = pg.to_frame().join(xg).sort_values("pig_absdelta", ascending=False)
    rho, p = stats.spearmanr(j["pig_absdelta"], j["xt_tv"])
    top_pig = list(j.index[:3])
    top_xt = list(j.sort_values("xt_tv", ascending=False).index[:3])
    print(f"\n=== {label}  (n_individuals = {pig['sample_id'].nunique()}) ===")
    with pd.option_context("display.float_format", lambda v: f"{v:11.6g}"):
        print(j.to_string())
    print(f"Spearman(pigmentation |Delta|, superpopulation TV) = {rho:.3f}  (p = {p:.4f})")
    print(f"top-3 pigmentation : {top_pig}")
    print(f"top-3 superpop     : {top_xt}")
    print(f"superpopulation flip rate, all genes pooled: {xt['flipped'].mean():.4f}")
    return {"label": label, "n_individuals": int(pig["sample_id"].nunique()),
            "spearman_rho": float(rho), "spearman_p": float(p),
            "top3_pigmentation": top_pig, "top3_superpopulation": top_xt,
            "pooled_flip_rate": float(xt["flipped"].mean()),
            "per_gene": json.loads(j.reset_index().to_json(orient="records"))}


def main():
    pig = pd.read_csv(PIG); pig = pig[pig["method"] == METHOD].copy()
    pig["delta"] = p_weak(pig)
    xt = pd.read_csv(XT); xt = xt[xt["method"] == METHOD].copy()

    payload = {"method": METHOD, "comparisons": []}
    payload["comparisons"].append(report(pig, xt, "all pigmentation test individuals"))

    if "in_superpop_train" in xt.columns and xt["in_superpop_train"].notna().any():
        held = xt[xt["in_superpop_train"] == 0]
        if len(held):
            ids = set(held["sample_id"])
            payload["comparisons"].append(
                report(pig[pig["sample_id"].isin(ids)], held,
                       "held out of the superpopulation training split"))
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
