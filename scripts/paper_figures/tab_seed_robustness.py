#!/usr/bin/env python3
"""Reviewer ask (Q3): does the per-gene classifier's accuracy, and the panel-vs-control
separation it supports, depend on the arbitrary choice of training seed?

Reads results/genotype_based_predictor/seed_robustness_table.csv (published seed 13 plus
retrained arms under seeds 7, 21, 99; scripts/experiments/seed_robustness_sweep.py) and
writes the coverage/AUC-per-seed table plus the per-gene stability numbers the paper's
Discussion cites. Regenerates as the backfill (some (gene, seed) arms still pending at the
time of writing) completes -- rerun this script to refresh the numbers before submission.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _common import CONTROL, PANEL, RES, TABDIR, auc_and_p

SEED_TABLE = RES / "seed_robustness_table.csv"
PUBLISHED_SEED = 13
SEEDS = [7, PUBLISHED_SEED, 21, 99]


def main():
    df = pd.read_csv(SEED_TABLE)
    genes_all = PANEL + CONTROL
    n_total = len(genes_all)

    rows = []
    pivot = df.pivot(index="gene", columns="seed", values="bal_acc")
    pub_rank = pivot[PUBLISHED_SEED]

    for s in SEEDS:
        sub = df[df["seed"] == s].set_index("gene")["bal_acc"].to_dict()
        genes_s = [g for g in genes_all if g in sub]
        mask_s = [g in PANEL for g in genes_s]
        auc, p, _mode = auc_and_p(sub, genes_s, mask_s)
        if s == PUBLISHED_SEED:
            rho, rho_p = 1.0, 0.0
        else:
            common = pivot[[PUBLISHED_SEED, s]].dropna()
            rho, rho_p = spearmanr(common[PUBLISHED_SEED], common[s])
        rows.append({"seed": s, "n": len(genes_s), "n_total": n_total,
                      "auc": auc, "p": p, "rho": rho, "rho_p": rho_p})

    complete = pivot.dropna()
    stds = complete.std(axis=1)
    n_complete = len(complete)

    lines = []
    lines.append(r"\begin{tabular}{rrrrrl}")
    lines.append(r"\toprule")
    lines.append(r"seed & $n$ genes & AUC & $p$ & $\rho$ (vs.\ seed 13) & \\")
    lines.append(r"\midrule")
    for r in rows:
        tag = r"\textit{published}" if r["seed"] == PUBLISHED_SEED else ""
        rho_str = "---" if r["seed"] == PUBLISHED_SEED else f"{r['rho']:.3f}"
        lines.append(f"{r['seed']} & {r['n']}/{r['n_total']} & {r['auc']:.4f} & "
                      f"{r['p']:.4f} & {rho_str} & {tag} \\\\")
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    out = TABDIR / "tab-seed-robustness.tex"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")

    print()
    print(f"n genes with all {len(SEEDS)} seeds complete: {n_complete}/{n_total}")
    print(f"mean per-gene std(bal_acc) across seeds: {stds.mean():.4f}")
    print(f"max per-gene std(bal_acc): {stds.max():.4f} ({stds.idxmax()})")
    print(f"total (gene, seed) arms complete: {len(df)}/{n_total * len(SEEDS)}")
    for r in rows:
        if r["seed"] != PUBLISHED_SEED:
            print(f"Spearman(seed13, seed{r['seed']}) n={r['n']}: "
                  f"rho={r['rho']:.3f} p={r['rho_p']:.2e}")


if __name__ == "__main__":
    main()
