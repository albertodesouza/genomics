#!/usr/bin/env python3
"""Figure 2: balanced accuracy of the AlphaGenome+CNN arms, panel against control.

Two distributions over the same axis, with every arm drawn as a point because at 9
against 33 the violin is an interpolation and the points are the data. The
majority-class baseline of the split is drawn, since that is what each accuracy is read
against.

BOTH SPLITS ARE SHOWN, side by side, because the two rankings disagree (Spearman on the
two splits is well short of 1) and a reader should be able to see that the
panel-against-control separation survives the disagreement rather than being an artefact
of whichever split is plotted. The left panel is validation, the right panel is the
held-out test split the accuracies are reported on elsewhere in the paper.

No selection is drawn on top of the points. An earlier version of this figure ringed the
ten arms carried into a multi-gene probe selected on validation; that probe and its
selection belong to a version of the paper that no longer exists, and ringing points by a
criterion the text does not otherwise mention was confusing on its own terms, since the
ring was a top ten of validation shown on a test axis that orders arms differently.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _common import (C_CTRL, C_DRAW2, C_DRAW3, C_PANEL, CONTROL, DRAW2, DRAW3, PANEL,
                     auc_and_p, load_final_table, savefig)

# Majority-class accuracy of each split: 115/162 on test, 108/149 on validation. Read off
# the class supports rather than reused across splits, which differ.
BASELINE = {"bal_acc": 0.7099, "bal_acc_val": 108 / 149}
SPLIT_LABEL = {"bal_acc_val": "validation", "bal_acc": "test"}


def draw_split(ax, df, col, rng):
    """One panel: panel-vs-control violins of `col`, every arm a point, top-10 ringed."""
    panel = df[df["gene"].isin(PANEL)]
    ctrl = df[df["gene"].isin(CONTROL)]
    pv, cv = panel[col].to_numpy(), ctrl[col].to_numpy()

    # The same estimator and the same null as Fig. rank-comparison, taken from _common
    # rather than recomputed here: this is the identical comparison on the identical
    # genes, and scipy's Mann-Whitney p and the label-assignment permutation p differ in
    # the third decimal, which would put two numbers for one test in one paper.
    genes = list(df["gene"])
    auc, pval, mode = auc_and_p(dict(zip(df["gene"], df[col])), genes,
                                [g in PANEL for g in genes])

    parts = ax.violinplot([pv, cv], positions=[0, 1], widths=0.72,
                          showextrema=False, showmedians=False)
    for body, c in zip(parts["bodies"], (C_PANEL, C_CTRL)):
        body.set_facecolor(c)
        body.set_alpha(0.22)
        body.set_edgecolor(c)
        body.set_linewidth(1.0)

    for pos, sub in ((0, panel), (1, ctrl)):
        for _, r in sub.iterrows():
            # Controls are drawn in shades by draw so each later draw stays visible as a
            # subset rather than merging into one anonymous cloud.
            c = (C_PANEL if r["gene"] in PANEL else
                 C_DRAW2 if r["gene"] in DRAW2 else C_DRAW3 if r["gene"] in DRAW3 else C_CTRL)
            ax.scatter(pos + rng.uniform(-0.13, 0.13), r[col], s=34, c=c,
                       lw=0.6, edgecolors="white", zorder=3)
        med = float(np.median(sub[col]))
        c = C_PANEL if pos == 0 else C_CTRL
        ax.hlines(med, pos - 0.30, pos + 0.30, color=c, lw=2.0, zorder=4)
        ax.annotate(f"median {med:.3f}", (pos + 0.33, med), fontsize=7.5, va="center",
                    color=c)

    base = BASELINE[col]
    ax.axhline(base, color="tab:red", ls="--", lw=1.0)
    ax.annotate(f"majority class {base:.4f}", (0.01, base),
                xycoords=("axes fraction", "data"), ha="left", va="bottom",
                fontsize=7, color="tab:red",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"panel\n($n = {len(pv)}$)",
                        f"control\n($n = {len(cv)}$, three draws)"])
    ax.set_xlim(-0.6, 1.72)
    # En-dash as a literal character: mathtext is off, so "--" would render as two
    # hyphens rather than as a dash.
    ax.set_title(f"{SPLIT_LABEL[col]} split\n$\\mathrm{{AUC}} = {auc:.4f}$ "
                 f"($p = {pval:.4f}$, permutation)", fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", lw=0.3, alpha=0.4)
    return pval, auc, float(np.median(pv)), float(np.median(cv)), mode


def main():
    df = load_final_table()
    if "bal_acc_val" not in df.columns:
        raise SystemExit("ABORT: final table has no bal_acc_val column; rerun "
                         "scripts/experiments/poolmax_final_table.py")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.4), sharey=True)
    stats_out = {}
    for ax, col in zip(axes, ("bal_acc_val", "bal_acc")):
        # Same jitter draw in both panels, so a given arm sits at the same horizontal
        # offset left and right and the eye can follow it across the splits.
        stats_out[col] = draw_split(ax, df, col, np.random.default_rng(13))
    axes[0].set_ylabel("balanced accuracy")

    fig.suptitle("One classifier per gene: balanced accuracy by class", fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 0.99))
    savefig(fig, "balacc-violin.png")

    for col, (p, auc, mp, mc, mode) in stats_out.items():
        print(f"{SPLIT_LABEL[col]:11s} panel median {mp:.4f}  control median {mc:.4f}  "
              f"AUC = {auc:.4f}  p = {p:.4g}  [{mode}]")
    ctrl = df[df["gene"].isin(CONTROL)]
    pmin = df[df["gene"].isin(PANEL)]["bal_acc"].min()
    print(f"controls above the worst panel arm (test): "
          f"{int((ctrl['bal_acc'] > pmin).sum())}/{len(ctrl)}")
    best = ctrl.loc[ctrl["bal_acc"].idxmax()]
    print(f"best control (test) {best['gene']} {best['bal_acc']:.4f} vs best panel "
          f"{df[df['gene'].isin(PANEL)]['bal_acc'].max():.4f}")


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    main()
