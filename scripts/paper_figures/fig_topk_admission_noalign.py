#!/usr/bin/env python3
"""Same pre-declared pAUC head-to-head as fig_topk_admission.py / _twas.py / _rawseq.py,
CNN vs the no-alignment (naive-crop) arm.

Reviewer ask (Q1, scoped to "no alignment" only): read each individual's own predicted
AlphaGenome track at a fixed, un-realigned position (the central 32,768bp of the raw
524,288bp window, no Fit/Apply coordinate mapping) instead of the aligned tensor Section
3.3 describes. scripts/experiments/noalign_sweep.py produces noalign_final_table.csv.
This script answers the same two questions fig_topk_admission_twas.py / _rawseq.py ask of
their arms, with the identical machinery (imported, not reimplemented):

  1. Does the no-alignment arm's own bal_acc separate panel from control at all?
  2. At the threshold that recovers k panel genes, how many controls has the CNN admitted
     against how many the no-alignment arm has admitted?

Writes a separate JSON/table so nothing the paper already reads is touched.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from _common import CONTROL, PANEL, RES, TABDIR, auc_and_p, load_final_table
from fig_topk_admission import CAP, admitted, nulls, pauc, pvals

NOALIGN_TABLE = RES / "noalign_final_table.csv"


def load_noalign() -> pd.DataFrame:
    df = pd.read_csv(NOALIGN_TABLE)
    missing = (set(PANEL) | set(CONTROL)) - set(df["gene"])
    if missing:
        raise SystemExit(f"ABORT: noalign_final_table.csv is missing {sorted(missing)}")
    return df


def main(args):
    noalign = load_noalign()
    cnn_df = load_final_table()
    genes = sorted(set(noalign["gene"]) & set(cnn_df["gene"]))
    mask = np.array([g in PANEL for g in genes], bool)
    n, n1 = len(genes), int(mask.sum())
    n0 = n - n1

    noalign_score = dict(zip(noalign["gene"], noalign["bal_acc"]))
    cnn_score = dict(zip(cnn_df["gene"], cnn_df["bal_acc"]))

    auc_rs, p_rs, mode = auc_and_p(noalign_score, genes, list(mask))
    print(f"no-alignment arm alone: AUC = {auc_rs:.4f}, p = {p_rs:.4f} ({mode}, n={n})")

    delta_panel = float(np.mean([noalign_score[g] - cnn_score[g] for g in genes if g in PANEL]))
    delta_control = float(np.mean([noalign_score[g] - cnn_score[g] for g in genes if g not in PANEL]))
    print(f"mean bal_acc delta (no_align - published): panel {delta_panel:+.4f}, "
          f"control {delta_control:+.4f}")

    r_cnn = pd.Series([cnn_score[g] for g in genes]).rank().to_numpy()
    r_rs = pd.Series([noalign_score[g] for g in genes]).rank().to_numpy()
    Cc = admitted(r_cnn, mask, n1, n0)[0]
    Cr = admitted(r_rs, mask, n1, n0)[0]
    lab, flip = nulls(r_cnn, r_rs, mask, args.b, args.seed)

    sweep = []
    for cap in [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00]:
        pc, pr = pauc(Cc, cap, n1, n0), pauc(Cr, cap, n1, n0)
        d = float(pc - pr)
        p1a, p2a, ca = pvals(d, *lab, cap, n1, n0)
        p1b, p2b, cb = pvals(d, *flip, cap, n1, n0)
        sweep.append({"cap": cap, "pauc_cnn": float(pc), "pauc_noalign": float(pr), "d": d,
                      "ceiling": float(1.0 - pr),
                      "p_label_2s": p2a, "p_flip_2s": p2b})
        print(f"  cap {cap:.2f}  pAUC_cnn {pc:.4f}  pAUC_noalign {pr:.4f}  d {d:+.4f}  "
              f"p2(label) {p2a:.4f}  p2(flip) {p2b:.4f}")

    jack = []
    for g in genes:
        keep = [x for x in genes if x != g]
        mk = np.array([x in PANEL for x in keep], bool)
        k1, k0 = int(mk.sum()), len(keep) - int(mk.sum())
        rc = pd.Series([cnn_score[x] for x in keep]).rank().to_numpy()
        rr = pd.Series([noalign_score[x] for x in keep]).rank().to_numpy()
        dj = float(pauc(admitted(rc, mk, k1, k0)[0], CAP, k1, k0)
                   - pauc(admitted(rr, mk, k1, k0)[0], CAP, k1, k0))
        lj, _ = nulls(rc, rr, mk, args.jb, args.seed)
        jack.append({"dropped": g, "panel": bool(g in PANEL), "d": dj,
                     "p_label_2s": pvals(dj, *lj, CAP, k1, k0)[1]})
    worst = max(jack, key=lambda r: r["p_label_2s"])
    print(f"jackknife worst: drop {worst['dropped']} d {worst['d']:+.4f} "
          f"p2 {worst['p_label_2s']:.4f}")

    out = {"n": n, "n_panel": n1, "n_control": n0, "cap_primary": CAP,
           "noalign_alone": {"auc": auc_rs, "p": p_rs, "mode": mode},
           "delta_panel": delta_panel, "delta_control": delta_control,
           "sweep": sweep, "jackknife": jack, "jackknife_worst": worst}
    dst = RES / "gwas_ranking" / "topk_admission_noalign.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dst}")

    lines = [r"\begin{tabular}{rrrrrr}", r"\toprule",
             r"$c$ & pAUC classifier & pAUC no-align & $d$ & $p$ (label) & $p$ (flip) \\",
             r"\midrule"]
    for s in sweep:
        mark = r"\ \textbf{*}" if s["cap"] == CAP else ""
        lines.append(f"{s['cap']:.2f}{mark} & {s['pauc_cnn']:.4f} & {s['pauc_noalign']:.4f} & "
                     f"{s['d']:+.4f} & {s['p_label_2s']:.4f} & {s['p_flip_2s']:.4f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    p = TABDIR / "tab-pauc-sweep-noalign.tex"
    p.write_text("\n".join(lines) + "\n")
    print(f"wrote {p}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", type=int, default=200_000)
    ap.add_argument("--jb", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=13)
    main(ap.parse_args())
