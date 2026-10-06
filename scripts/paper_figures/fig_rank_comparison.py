#!/usr/bin/env python3
"""Figure 3: the same genes ranked by variant-level association and by the classifier.

A slope chart, because the question is not what either score is but whether the panel
sits higher in one ordering than in the other. Rank 1 is the strongest gene under that
method. Each method is annotated with its own panel-against-control AUC and the
permutation p-value of that AUC; the head-to-head difference is annotated between them,
calibrated by a paired permutation that reassigns the panel/control labels while keeping
each gene's two scores together.
"""
from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import (C_CTRL, C_DRAW2, C_DRAW3, C_PANEL, DRAW2, DRAW3, PANEL, auc_and_p,
                     load_final_table, load_gwas, savefig)


def paired_headtohead(a: dict, b: dict, genes, is_panel, b_max=2_000_000, seed=13):
    """Two-sided p for AUC(a) - AUC(b) under reassignment of the panel/control labels.

    The permutation is over label assignments rather than over scores, so each gene keeps
    its pair (a_g, b_g) together and the correlation between the two methods is preserved
    -- which is what makes this a test of the difference and not of the two margins.
    """
    from itertools import combinations
    from math import comb

    ra = pd.Series([a[g] for g in genes]).rank().to_numpy()
    rb = pd.Series([b[g] for g in genes]).rank().to_numpy()
    m = np.asarray(is_panel, dtype=bool)
    n, n1 = len(genes), int(m.sum())
    den = n1 * (n - n1)
    off = n1 * (n1 + 1) / 2

    def d(mask_idx):
        return ((ra[mask_idx].sum(axis=-1) - off) - (rb[mask_idx].sum(axis=-1) - off)) / den

    obs = float(d(np.flatnonzero(m)))
    total = comb(n, n1)
    if total <= b_max:
        idx = np.array(list(combinations(range(n), n1)), dtype=np.int64)
        mode = f"exact ({total:,})"
    else:
        rng = np.random.default_rng(seed)
        idx = np.concatenate([np.argsort(rng.random((min(200_000, b_max - lo), n)),
                                        axis=1)[:, :n1]
                              for lo in range(0, b_max, 200_000)], axis=0)
        mode = f"Monte Carlo (B = {len(idx):,} of {total:,})"
    null = d(idx)
    p = float((np.abs(null) >= abs(obs) - 1e-12).mean())
    crit = float(np.quantile(np.abs(null), 0.95))
    return obs, p, crit, mode


def main(args):
    df = load_final_table()
    col = {"test": "bal_acc", "val": "bal_acc_val"}[args.split]
    other_col = {"test": "bal_acc_val", "val": "bal_acc"}[args.split]
    other_label = {"test": "val", "val": "test"}[args.split]
    if other_col not in df.columns:
        raise SystemExit("ABORT: final table has no bal_acc_val column; rerun "
                         "scripts/experiments/poolmax_final_table.py")
    cnn = dict(zip(df["gene"], df[col]))
    other = dict(zip(df["gene"], df[other_col]))
    gw = load_gwas(args.arm, str(args.n))["gene_scores"][args.res]

    genes = sorted(set(cnn) & set(gw))
    is_panel = [g in PANEL for g in genes]
    gwas = {g: gw[g][args.est] for g in genes}

    a_g, p_g, mode_g = auc_and_p(gwas, genes, is_panel)
    a_c, p_c, mode_c = auc_and_p(cnn, genes, is_panel)
    # The CNN ranking on the other split, reported in the figure rather than left to the
    # reader: the two splits rank these arms differently (Spearman +0.895), so a single
    # number invites the question of whether the comparison turns on the choice.
    a_o, p_o, _ = auc_and_p(other, genes, is_panel)
    dauc, p_h, crit, mode_h = paired_headtohead(cnn, gwas, genes, is_panel)

    # Rank 1 = strongest under that method.
    rg = pd.Series(gwas).rank(ascending=False, method="min")
    rc = pd.Series(cnn).rank(ascending=False, method="min")

    # Ties share a rank under method="min", so their labels would be drawn on top of
    # each other; nudge the label (not the point) when a rank is occupied twice.
    def nudges(r):
        out, seen = {}, {}
        for g in sorted(genes, key=lambda x: (r[x], x)):
            k = r[g]
            seen[k] = seen.get(k, 0) + 1
            out[g] = 0.0
        for g in genes:
            k = r[g]
            if seen[k] > 1:
                sib = sorted([x for x in genes if r[x] == k])
                out[g] = (sib.index(g) - (len(sib) - 1) / 2) * 0.34
        return out

    ng, nc = nudges(rg), nudges(rc)

    fig, ax = plt.subplots(figsize=(7.6, 9.0))
    for g in genes:
        col = (C_PANEL if g in PANEL else
               C_DRAW2 if g in DRAW2 else C_DRAW3 if g in DRAW3 else C_CTRL)
        lw = 1.6 if g in PANEL else 1.0
        ax.plot([0, 1], [rg[g], rc[g]], color=col, lw=lw,
                alpha=0.95 if g in PANEL else 0.7, zorder=2 if g in PANEL else 1)
        ax.scatter([0, 1], [rg[g], rc[g]], s=26, c=col, zorder=3, lw=0.5,
                   edgecolors="white")
        ax.annotate(g, (-0.035, rg[g] + ng[g]), ha="right", va="center", fontsize=8,
                    color=col)
        ax.annotate(g, (1.035, rc[g] + nc[g]), ha="left", va="center", fontsize=8,
                    color=col)

    ax.set_xlim(-0.30, 1.30)
    ax.set_ylim(len(genes) + 0.6, 0.4)
    ax.set_yticks(range(1, len(genes) + 1))
    ax.set_ylabel("rank (1 = strongest under that method)")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"GWAS\n{EST_LABEL[args.est]}, {args.res}",
                        f"AlphaGenome + CNN\nbalanced accuracy ({args.split} split)"],
                       fontsize=9.5)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.grid(axis="y", lw=0.25, alpha=0.35)

    h = [plt.Line2D([], [], color=c, lw=2) for c in (C_PANEL, C_CTRL, C_DRAW2, C_DRAW3)]
    ax.legend(h, [f"panel ($n = {sum(is_panel)}$)", "control, draw 1", "control, draw 2",
                 "control, draw 3"],
              fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.085), ncol=4,
              frameon=False)

    note = (f"panel-vs-control AUC:   GWAS {a_g:.4f} (p = {p_g:.4f})"
            f"    CNN {a_c:.4f} (p = {p_c:.4f})\n"
            f"head to head:   dAUC = {dauc:+.4f},  two-sided p = {p_h:.4f},"
            f"  detectable at |dAUC| \u2265 {crit:.4f}\n"
            f"cross-split check:   CNN on {other_label} {a_o:.4f} (p = {p_o:.4f})")
    fig.text(0.5, 0.012, note, ha="center", va="bottom", fontsize=8.5,
             family="monospace",
             bbox=dict(fc="white", ec="0.7", lw=0.6, pad=5))
    ax.set_title(f"{len(genes)} genes ranked two ways "
                 f"({args.arm} GWAS arm, {args.res} resolution)", fontsize=10)
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    stem = f"rank-comparison-{args.est}-{args.res}"
    savefig(fig, f"{stem}.png" if args.split == "test" else f"{stem}-{args.split}.png")

    print(f"n = {len(genes)} genes, {sum(is_panel)} panel")
    print(f"GWAS {args.est}/{args.res}: AUC {a_g:.4f} p {p_g:.4f}  [{mode_g}]")
    print(f"CNN  bal_acc ({args.split:4s}) : AUC {a_c:.4f} p {p_c:.4f}  [{mode_c}]")
    print(f"CNN  bal_acc ({other_label:4s}) : AUC {a_o:.4f} p {p_o:.4f}  (cross-split check)")
    print(f"dAUC (CNN - GWAS) {dauc:+.4f} two-sided p {p_h:.4f} "
          f"crit {crit:.4f}  [{mode_h}]")
    print("Spearman(GWAS, CNN) = "
          f"{pd.Series(gwas).corr(pd.Series(cnn), method='spearman'):+.3f}")


# Literal unicode: mathtext is off for these labels, so "x" would read as a variable and
# "--" as two hyphens.
EST_LABEL = {"S_sidak": "\u0160id\u00e1k \u00d7 Li & Ji $m_{eff}$",
             "S_min": "raw min $p$",
             "S_gates": "GATES",
             "S_meanchi2": "mean $\\chi^2$"}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=33)
    ap.add_argument("--arm", default="uncorrected")
    ap.add_argument("--res", default="window", choices=["window", "crop"])
    ap.add_argument("--est", default="S_sidak", choices=list(EST_LABEL))
    ap.add_argument("--split", default="test", choices=["test", "val"],
                    help="which split's balanced accuracy is the CNN ranking")
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    main(ap.parse_args())
