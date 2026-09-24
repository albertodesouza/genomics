#!/usr/bin/env python3
"""Knockdown readout: the old 3-ontology/2-strand arms against the new 1-ontology arms.

Both sides are one classifier per gene, replayed on the same 162 held-out individuals
with the same 100 bp promoter scramble at the same MANE Select TSS, and both CSVs carry
delta_log_odds, so the response is compared on one scale with no conversion.

What differs between the two sides is THREE things at once, not one:
  channels      6 per haplotype (CL:1000458 + CL:0000346 + CL:2000092, both strands)
                -> 1 (CL:1000458, the gene's own strand)
  alignment     raw_center_crop (none) -> bcftools_chain
  pooling       global average, 400 epochs, scheduler on
                -> global max, 600 epochs, no scheduler
so a difference below is attributable to the configuration change as a whole and not to
the channel count by itself.

Only method=biology_tss is used on both sides: the old CSVs carry three anchorings and
the new one, so an unfiltered average would compare a different number of scrambles.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RES = Path("/home/breno/I2CA/genomics/results/genotype_based_predictor")
OLD = RES / "knockout_bulk"
NEW = RES / "knockout_bulk" / "poolmax"
METHOD = "biology_tss"

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
DRAW1 = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
         "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
NOT_PIGM = ["EDAR", "TCHH"]        # in the old draw, excluded from the new analysis


def summarise(path: Path, gene: str) -> dict | None:
    if not path.exists():
        return None
    d = pd.read_csv(path)
    d = d[d["method"] == METHOD]
    if d.empty:
        return None
    assert set(d["gene"]) == {gene}, (gene, set(d["gene"]))
    # Balanced accuracy straight from the replayed rows: every row carries the true
    # label and whether the unperturbed prediction was right, so the two recalls are
    # available without the checkpoint. This is the baseline accuracy on exactly the
    # individuals the knockdown was replayed on, on both sides.
    rec = d.groupby("true_label")["correct_baseline"].mean()
    # Probability-scale Delta, the statistic the previous paper reported, recomputed from
    # the logits so the two scales are compared on identical rows.
    def pw(a, b):
        return 1.0 / (1.0 + np.exp(-(b - a)))
    dprob = (pw(d["perturbed_strong_logit"], d["perturbed_weak_logit"])
             - pw(d["baseline_strong_logit"], d["baseline_weak_logit"]))
    return {
        "n": len(d),
        "bal": float(rec.mean()),
        "delta_prob": float(dprob.mean()),
        "delta": float(d["delta_log_odds"].mean()),
        "sd": float(d["delta_log_odds"].std(ddof=1)),
        "flip": float(d["flipped"].mean()),
        "delta_in": float(d["delta_in"].mean()),
        "acc_base": float(d["correct_baseline"].mean()),
        "acc_pert": float(d["correct_perturbed"].mean()),
    }


def bal_acc_from_cm(cm) -> float:
    cm = np.asarray(cm, dtype=float)
    return float(np.mean(np.diag(cm) / cm.sum(axis=1)))


def main():
    arms = {a["gene"]: a for a in
            json.loads((RES / "single_gene_arms.json").read_text())["arms"]}
    new_tbl = pd.read_csv(RES / "poolmax_final_table.csv").set_index("gene")

    rows = []
    for g in PANEL + DRAW1 + NOT_PIGM:
        o = summarise(OLD / f"pigmentation_test_split_knockout_single_{g.lower()}.csv", g)
        n = summarise(NEW / f"single_{g.lower()}.csv", g)
        if o is None or n is None:
            print(f"  (missing) {g}: old={o is not None} new={n is not None}")
            continue
        rows.append({
            "gene": g,
            "class": "panel" if g in PANEL else ("excluded" if g in NOT_PIGM else "control"),
            "old_bal": o["bal"], "new_bal": n["bal"],
            "new_bal_ckpt": float(new_tbl.loc[g, "bal_acc"]) if g in new_tbl.index else np.nan,
            "old_dprob": o["delta_prob"], "new_dprob": n["delta_prob"],
            "old_delta": o["delta"], "new_delta": n["delta"],
            "old_sd": o["sd"], "new_sd": n["sd"],
            "old_flip": o["flip"], "new_flip": n["flip"],
            "old_din": o["delta_in"], "new_din": n["delta_in"],
        })
    df = pd.DataFrame(rows)
    df["old_abs"] = df["old_delta"].abs()
    df["new_abs"] = df["new_delta"].abs()

    pd.set_option("display.width", 200)
    print("\n=== per gene, method=biology_tss, Delta on the log-odds scale "
          "(9 panel + 11 control draw 1; EDAR/TCHH shown but excluded from the tests)")
    show = df[["gene", "class", "old_bal", "new_bal", "old_delta", "new_delta",
               "old_dprob", "new_dprob", "old_flip", "new_flip"]].copy()
    print(show.to_string(index=False, float_format=lambda v: f"{v:8.3f}"))

    core = df[df["class"] != "excluded"]
    print("\n=== panel vs control, Mann-Whitney one-sided (panel > control) on |Delta|")
    for tag, col in (("old 6-channel", "old_abs"), ("new 1-channel", "new_abs")):
        p = core[core["class"] == "panel"][col]
        c = core[core["class"] == "control"][col]
        u = stats.mannwhitneyu(p, c, alternative="greater")
        print(f"  {tag:14s} panel median {p.median():7.3f}  control median {c.median():7.3f}"
              f"  p = {u.pvalue:.4f}  AUC = {u.statistic/(len(p)*len(c)):.4f}")

    print("\n=== the previous paper's statistic: |Delta| on the PROBABILITY scale")
    for tag, col in (("old 6-channel", "old_dprob"), ("new 1-channel", "new_dprob")):
        p = core[core["class"] == "panel"][col].abs()
        c = core[core["class"] == "control"][col].abs()
        u = stats.mannwhitneyu(p, c, alternative="greater")
        print(f"  {tag:14s} panel median {p.median():7.4f}  control median {c.median():7.4f}"
              f"  p = {u.pvalue:.4f}  AUC = {u.statistic/(len(p)*len(c)):.4f}")
    # With EDAR and TCHH back in, which is how the previous paper ran it (22 arms).
    full = df.copy()
    full["cls2"] = np.where(full["class"] == "control", "control", "panel")
    p = full[full["cls2"] == "panel"]["old_dprob"].abs()
    c = full[full["cls2"] == "control"]["old_dprob"].abs()
    u = stats.mannwhitneyu(p, c, alternative="greater")
    print(f"  old, 11 vs 11 (EDAR+TCHH back in, as the previous paper ran it): "
          f"panel median {p.median():.4f} control median {c.median():.4f} "
          f"p = {u.pvalue:.4f} AUC = {u.statistic/(len(p)*len(c)):.4f}")

    print("\n=== the same on balanced accuracy (both computed on the replayed rows)")
    for tag, col in (("old 6-channel", "old_bal"), ("new 1-channel", "new_bal")):
        p = core[core["class"] == "panel"][col].dropna()
        c = core[core["class"] == "control"][col].dropna()
        u = stats.mannwhitneyu(p, c, alternative="greater")
        print(f"  {tag:14s} panel median {p.median():7.4f}  control median {c.median():7.4f}"
              f"  p = {u.pvalue:.4f}  AUC = {u.statistic/(len(p)*len(c)):.4f}")

    print("\n=== do the two setups agree on which gene responds?")
    for tag, a, b in (("|Delta|", "old_abs", "new_abs"),
                      ("signed Delta", "old_delta", "new_delta"),
                      ("balanced accuracy", "old_bal", "new_bal")):
        s = stats.spearmanr(core[a], core[b], nan_policy="omit")
        print(f"  {tag:18s} Spearman rho = {s.statistic:+.3f}  p = {s.pvalue:.4f}  "
              f"(n = {core[[a, b]].dropna().shape[0]})")
    sign_agree = int((np.sign(core["old_delta"]) == np.sign(core["new_delta"])).sum())
    print(f"  sign of Delta agrees on {sign_agree}/{len(core)} arms")

    print("\n=== magnitude of the response, both scales")
    for tag, col in (("old 6-channel", "old_abs"), ("new 1-channel", "new_abs")):
        v = core[col]
        print(f"  {tag:14s} median {v.median():8.3f}  IQR [{v.quantile(.25):.3f}, "
              f"{v.quantile(.75):.3f}]  max {v.max():8.3f} ({core.loc[v.idxmax(),'gene']})")
    print("\n=== flip rate")
    for tag, col in (("old 6-channel", "old_flip"), ("new 1-channel", "new_flip")):
        v = core[col]
        print(f"  {tag:14s} median {v.median():.3f}  max {v.max():.3f} "
              f"({core.loc[v.idxmax(),'gene']})")

    print("\n=== raw Delta_in is NOT comparable across setups "
          "(summed over 6 channels vs 1), reported only to show the factor")
    r = (core["old_din"] / core["new_din"]).replace([np.inf, -np.inf], np.nan).dropna()
    print(f"  old/new ratio: median {r.median():.2f}  range [{r.min():.2f}, {r.max():.2f}]")

    out = Path("/tmp/claude-1001/-home-breno-I2CA/"
               "bca70c00-f8fc-4512-b818-f91b2702a636/scratchpad/setup_comparison.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
