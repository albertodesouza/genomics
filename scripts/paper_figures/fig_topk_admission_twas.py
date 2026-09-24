#!/usr/bin/env python3
"""Same pre-declared pAUC head-to-head as fig_topk_admission.py, CNN vs the TWAS arm.

Reviewer ask (Q5): a TWAS-style baseline. scripts/experiments/twas_scalar_sweep.py
produces twas_final_table.csv (the "poor man's PrediXcan" scalar-expression
LogisticRegression, one arm per gene, reusing the exact tensor cache / gene set as the
CNN). This script answers two questions with the identical machinery already used for the
CNN-vs-GWAS comparison (imported, not reimplemented, so the two are the same statistic):

  1. Does the TWAS arm's own bal_acc separate panel from control at all (the Section 4.2
     analogue, auc_and_p)?
  2. At the threshold that recovers k panel genes, how many controls has the CNN admitted
     against how many the TWAS arm has admitted (the Section 4.3 analogue, admitted/pauc/
     nulls/pvals from fig_topk_admission.py)?

Writes a separate JSON/tables so nothing the paper already reads is touched.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from _common import CONTROL, PANEL, RES, TABDIR, auc_and_p, load_final_table
from fig_topk_admission import CAP, admitted, nulls, pauc, pvals

TWAS_TABLE = RES / "twas_final_table.csv"


def load_twas() -> pd.DataFrame:
    df = pd.read_csv(TWAS_TABLE)
    missing = (set(PANEL) | set(CONTROL)) - set(df["gene"])
    if missing:
        raise SystemExit(f"ABORT: twas_final_table.csv is missing {sorted(missing)}")
    return df


def main(args):
    twas = load_twas()
    cnn_df = load_final_table()
    genes = sorted(set(twas["gene"]) & set(cnn_df["gene"]))
    mask = np.array([g in PANEL for g in genes], bool)
    n, n1 = len(genes), int(mask.sum())
    n0 = n - n1

    twas_score = dict(zip(twas["gene"], twas["bal_acc"]))
    cnn_score = dict(zip(cnn_df["gene"], cnn_df["bal_acc"]))

    # (1) TWAS's own panel-vs-control separation, same statistic Section 4.2 reports for
    # the CNN and the paper's random-label control uses throughout.
    auc_twas, p_twas, mode = auc_and_p(twas_score, genes, list(mask))
    print(f"TWAS arm alone: AUC = {auc_twas:.4f}, p = {p_twas:.4f} ({mode}, n={n})")

    # (2) Head-to-head against the CNN, identical pAUC/nulls/jackknife machinery as the
    # GWAS comparison.
    r_cnn = pd.Series([cnn_score[g] for g in genes]).rank().to_numpy()
    r_twas = pd.Series([twas_score[g] for g in genes]).rank().to_numpy()
    Cc = admitted(r_cnn, mask, n1, n0)[0]
    Ct = admitted(r_twas, mask, n1, n0)[0]
    lab, flip = nulls(r_cnn, r_twas, mask, args.b, args.seed)

    sweep = []
    for cap in [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00]:
        pc, pt = pauc(Cc, cap, n1, n0), pauc(Ct, cap, n1, n0)
        d = float(pc - pt)
        p1a, p2a, ca = pvals(d, *lab, cap, n1, n0)
        p1b, p2b, cb = pvals(d, *flip, cap, n1, n0)
        sweep.append({"cap": cap, "pauc_cnn": float(pc), "pauc_twas": float(pt), "d": d,
                      "ceiling": float(1.0 - pt),
                      "p_label_2s": p2a, "p_flip_2s": p2b})
        print(f"  cap {cap:.2f}  pAUC_cnn {pc:.4f}  pAUC_twas {pt:.4f}  d {d:+.4f}  "
              f"p2(label) {p2a:.4f}  p2(flip) {p2b:.4f}")

    jack = []
    for g in genes:
        keep = [x for x in genes if x != g]
        mk = np.array([x in PANEL for x in keep], bool)
        k1, k0 = int(mk.sum()), len(keep) - int(mk.sum())
        rc = pd.Series([cnn_score[x] for x in keep]).rank().to_numpy()
        rt = pd.Series([twas_score[x] for x in keep]).rank().to_numpy()
        dj = float(pauc(admitted(rc, mk, k1, k0)[0], CAP, k1, k0)
                   - pauc(admitted(rt, mk, k1, k0)[0], CAP, k1, k0))
        lj, _ = nulls(rc, rt, mk, args.jb, args.seed)
        jack.append({"dropped": g, "panel": bool(g in PANEL), "d": dj,
                     "p_label_2s": pvals(dj, *lj, CAP, k1, k0)[1]})
    worst = max(jack, key=lambda r: r["p_label_2s"])
    print(f"jackknife worst: drop {worst['dropped']} d {worst['d']:+.4f} "
          f"p2 {worst['p_label_2s']:.4f}")

    out = {"n": n, "n_panel": n1, "n_control": n0, "cap_primary": CAP,
           "twas_alone": {"auc": auc_twas, "p": p_twas, "mode": mode},
           "sweep": sweep, "jackknife": jack, "jackknife_worst": worst}
    dst = RES / "gwas_ranking" / "topk_admission_twas.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dst}")

    lines = [r"\begin{tabular}{rrrrrr}", r"\toprule",
             r"$c$ & pAUC classifier & pAUC TWAS & $d$ & $p$ (label) & $p$ (flip) \\",
             r"\midrule"]
    for s in sweep:
        mark = r"\ \textbf{*}" if s["cap"] == CAP else ""
        lines.append(f"{s['cap']:.2f}{mark} & {s['pauc_cnn']:.4f} & {s['pauc_twas']:.4f} & "
                     f"{s['d']:+.4f} & {s['p_label_2s']:.4f} & {s['p_flip_2s']:.4f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    p = TABDIR / "tab-pauc-sweep-twas.tex"
    p.write_text("\n".join(lines) + "\n")
    print(f"wrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", type=int, default=200_000)
    ap.add_argument("--jb", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=13)
    main(ap.parse_args())
