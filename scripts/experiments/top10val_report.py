#!/usr/bin/env python3
"""What the val-selected top-10 arm shows, against its own null and against the old arm.

Three questions, in the order they have to be answered:

  1. Did selecting on validation instead of test unsaturate the classifier? The old arm
     reached bal_acc 1.0000 on the split it was selected on, with a median baseline margin
     of 10.9 log-odds, which left the knockdown nothing to move (8 of its 10 genes flipped
     no prediction). If the new arm is equally saturated, the fix did not work and the
     knockdown is still being read on a model with no decision to perturb.
  2. Does each gene's response clear its OWN null band? The null is the same 100 bp
     scramble placed at a random window position instead of the promoter, one draw per
     (individual, gene), replayed through this same model. This is the test the backup
     paper passed and the new paper did not have.
  3. Does |Delta| separate panel from control on this arm, now that the composition is
     5 and 5 rather than 7 and 3?

Reads finished CSVs only. No API calls, no GPU.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit

RES = Path("/home/breno/I2CA/genomics/results/genotype_based_predictor")
KD = RES / "knockout_bulk"
PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]


def load(p: Path) -> pd.DataFrame:
    d = pd.read_csv(p)
    d["lo_b"] = d["baseline_strong_logit"] - d["baseline_weak_logit"]
    d["lo_p"] = d["perturbed_strong_logit"] - d["perturbed_weak_logit"]
    d["dprob"] = expit(d["lo_p"]) - expit(d["lo_b"])
    return d


def auc_p(panel_vals, ctrl_vals):
    """AUC = P(a panel gene outranks a control gene), with the one-sided Mann-Whitney."""
    u = stats.mannwhitneyu(panel_vals, ctrl_vals, alternative="greater")
    return u.statistic / (len(panel_vals) * len(ctrl_vals)), u.pvalue


def main() -> int:
    kd = load(KD / "top10val" / "top10val_knockdown.csv")
    null = load(KD / "top10val" / "top10val_null.csv")
    genes = sorted(kd["gene"].unique())

    print("=== 1. saturation: is there a decision left to move? ===")
    for name, p in (("old top10 (test-selected)", KD / "top10" / "top10_knockdown.csv"),
                    ("new top10val (val-selected)", KD / "top10val" / "top10val_knockdown.csv")):
        d = load(p)
        b = d.drop_duplicates("sample_id")["lo_b"]
        pc = expit(b.abs())
        print(f"  {name:28s} |log-odds| median {b.abs().median():6.2f}   "
              f"P(correct) median {pc.median():.6f}   "
              f"slope {(pc*(1-pc)).median():.2e}   "
              f"in [0.01,0.99] {((expit(b)>0.01)&(expit(b)<0.99)).mean()*100:5.1f}%   "
              f"flips {int(d['flipped'].sum()):4d}")

    print("\n=== 2. each gene against its own null band ===")
    print(f"  {'gene':9s} {'class':8s} {'|D|':>9s} {'null p50':>9s} {'null p95':>9s} "
          f"{'null max':>9s} {'ratio':>8s} {'clears':>7s} {'flip':>5s} {'nullflip':>8s}")
    rows = []
    for g in genes:
        k = kd[kd["gene"] == g]
        n = null[null["gene"] == g]
        obs = abs(float(k["delta_log_odds"].mean()))
        na = n["delta_log_odds"].abs()
        clears = obs > na.max()
        rows.append({"gene": g, "panel": g in PANEL, "absd": obs,
                     "null_p95": float(na.quantile(.95)), "null_max": float(na.max()),
                     "clears": clears})
        print(f"  {g:9s} {'panel' if g in PANEL else 'control':8s} {obs:9.4f} "
              f"{na.median():9.4f} {na.quantile(.95):9.4f} {na.max():9.4f} "
              f"{obs/max(na.median(),1e-12):8.1f} {'yes' if clears else 'NO':>7s} "
              f"{int(k['flipped'].sum()):5d} {int(n['flipped'].sum()):8d}")
    r = pd.DataFrame(rows)
    print(f"\n  pooled null: median |D| {null['delta_log_odds'].abs().median():.4f}, "
          f"p95 {null['delta_log_odds'].abs().quantile(.95):.4f}, "
          f"max {null['delta_log_odds'].abs().max():.4f}, "
          f"flips {int(null['flipped'].sum())}/{len(null)}")
    print(f"  genes clearing their own null: {int(r['clears'].sum())}/{len(r)} "
          f"(panel {int(r[r.panel]['clears'].sum())}/{int(r.panel.sum())}, "
          f"control {int(r[~r.panel]['clears'].sum())}/{int((~r.panel).sum())})")

    print("\n=== 3. panel vs control on this arm ===")
    pv, cv = r[r.panel]["absd"].to_numpy(), r[~r.panel]["absd"].to_numpy()
    a, p = auc_p(pv, cv)
    print(f"  |Delta|  panel median {np.median(pv):.4f} (n={len(pv)})  "
          f"control median {np.median(cv):.4f} (n={len(cv)})  AUC {a:.4f}  p {p:.4f}")
    print("\n  per gene, ordered by |Delta|:")
    for _, x in r.sort_values("absd", ascending=False).iterrows():
        print(f"    {x['gene']:9s} {'panel' if x['panel'] else 'control':8s} "
              f"{x['absd']:9.4f}  {'clears null' if x['clears'] else 'within null'}")
    out = KD / "top10val" / "top10val_report.json"
    out.write_text(json.dumps(r.to_dict("records"), indent=1))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
