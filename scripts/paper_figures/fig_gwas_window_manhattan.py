#!/usr/bin/env python3
"""Manhattan of the selected windows only, laid end to end -- panel against control.

The backup paper carried this as panel (b) of its Manhattan figure, over the eleven
panel windows the published scan covered. Here it covers every window the gene ranking
scores, so for the first time the panel windows and the control windows appear on the
same axis: that comparison is the point, because the claim being read against this
figure is that the variant-level estimator returns no usable ordering over these genes.

Data source. The per-variant -log10 p comes from the plink2 run of the gene-ranking arm
(scripts/experiments/gwas_gene_ranking.py), extracted from its work directory into
results/genotype_based_predictor/gwas_ranking/union_uncorrected.neglogp.tsv.gz so this
figure does not depend on a /tmp directory surviving. The extraction is verified against
the ranking JSON: the variant count inside each of the 44 windows equals that window's
`m`, for all 44.

Windows are ordered by the primary gene score within each class, matching
fig_gwas.py, so the figure reads left to right as the ranking the estimator induces --
and the reader can see that the ranking is not reflected in the point clouds.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from _common import (C_CTRL, C_DRAW2, C_DRAW3, C_PANEL, DRAW1, DRAW2, DRAW3, GENOMEWIDE,
                     PANEL, RES, gene_class, load_gwas, savefig, window_interval)

VARIANTS = RES / "gwas_ranking" / "union_uncorrected.neglogp.tsv.gz"
THRESH = -np.log10(GENOMEWIDE)
GROUPS = [("panel", PANEL, C_PANEL), ("control (draw 1)", DRAW1, C_CTRL),
          ("control (draw 2)", DRAW2, C_DRAW2), ("control (draw 3)", DRAW3, C_DRAW3)]
EST = "S_sidak"


def load_variants():
    v = pd.read_csv(VARIANTS, sep="\t")
    v["CHROM"] = v["CHROM"].astype(str).str.replace("chr", "", regex=False)
    return v


def per_window(v, genes):
    """{gene: array of -log10 p inside its 524,288 bp window}, with the count checked
    against the ranking JSON so a silently truncated extraction cannot be plotted."""
    gw = load_gwas("uncorrected", "44")["gene_scores"]["window"]
    out = {}
    for g in genes:
        c, s, e = window_interval(g)
        c = str(c).replace("chr", "")
        sel = v[(v["CHROM"] == c) & (v["POS"] >= s) & (v["POS"] <= e)]
        assert len(sel) == gw[g]["m"], f"{g}: {len(sel)} variants, JSON says {gw[g]['m']}"
        out[g] = sel["NEG_LOG10_P"].to_numpy()
    return out, gw


def main() -> int:
    v = load_variants()
    genes_all = [g for _, gs, _ in GROUPS for g in gs]
    vals, gw = per_window(v, genes_all)

    # Ordered by the primary estimator inside each class, as in fig_gwas.py.
    ordered, colours, bounds = [], [], []
    for name, gs, colour in GROUPS:
        for g in sorted(gs, key=lambda x: -gw[x][EST]):
            ordered.append(g)
            colours.append(colour)
        bounds.append((name, len(ordered)))

    fig, ax = plt.subplots(figsize=(6.9, 3.9))
    rng = np.random.default_rng(13)
    fracs = {}
    for i, (g, colour) in enumerate(zip(ordered, colours)):
        y = vals[g]
        # Jitter inside the slot: position within a window is not comparable between
        # windows, so spreading the points is honest and the alternative (true position)
        # would suggest a within-window structure the figure is not about.
        x = i + 0.5 + rng.uniform(-0.38, 0.38, size=len(y))
        if i % 2 == 0:
            ax.axvspan(i, i + 1, color="0.93", lw=0, zorder=0)
        below = y < THRESH
        ax.scatter(x[below], y[below], s=0.7, color=colour, alpha=0.20, lw=0, zorder=2)
        ax.scatter(x[~below], y[~below], s=0.7, color=colour, alpha=0.80, lw=0, zorder=3)
        fracs[g] = float((~below).mean())

    ax.axhline(THRESH, color="#c0392b", lw=0.9, ls="--", zorder=4)
    ax.set_xlim(0, len(ordered))
    ax.set_ylim(-1.5, max(v["NEG_LOG10_P"].max() * 1.02, THRESH + 4))
    ax.set_xticks(np.arange(len(ordered)) + 0.5)
    ax.set_xticklabels(ordered, rotation=90, fontsize=6.2)
    ax.set_ylabel("$-\\log_{10} p$", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    # Against the right edge with a white patch behind it: every window is dense right
    # down to the line, so there is no interior gap this label could sit in.
    ax.text(len(ordered) - 0.25, THRESH + 1.2, "$p = 5\\times10^{-8}$", fontsize=7,
            color="#c0392b", ha="right", va="bottom",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))

    # Class brackets above the axis, so the three groups are separable without colour.
    ytop = ax.get_ylim()[1]
    prev = 0
    for name, end in bounds:
        sub = ordered[prev:end]
        lo, hi = min(fracs[g] for g in sub), max(fracs[g] for g in sub)
        ax.text((prev + end) / 2, ytop * 1.008,
                f"{name} ($n = {len(sub)}$)\n{lo*100:.0f}\u2013{hi*100:.0f}% above the line",
                fontsize=7.2, ha="center", va="bottom", linespacing=1.35)
        if end < len(ordered):
            ax.axvline(end, color="0.35", lw=0.8, zorder=5)
        prev = end

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    savefig(fig, "gwas-window-manhattan.png")
    plt.close(fig)

    n_tot = sum(len(vals[g]) for g in ordered)
    n_hit = sum(int((vals[g] >= THRESH).sum()) for g in ordered)
    print(f"\n{len(ordered)} windows, {n_tot:,} variants, "
          f"{n_hit:,} above 5e-8 ({n_hit/n_tot*100:.1f}%)")
    for name, gs, _ in GROUPS:
        sub = [g for g in ordered if g in gs]
        t = sum(len(vals[g]) for g in sub)
        h = sum(int((vals[g] >= THRESH).sum()) for g in sub)
        print(f"  {name:18s} {len(sub):2d} windows, {t:6,} variants, "
              f"{h/t*100:5.1f}% above; per-window "
              f"{min(fracs[g] for g in sub)*100:.1f}-{max(fracs[g] for g in sub)*100:.1f}%")
    worst = min(fracs, key=fracs.get); best = max(fracs, key=fracs.get)
    print(f"  lowest  {worst} ({gene_class(worst)}) {fracs[worst]*100:.1f}%")
    print(f"  highest {best} ({gene_class(best)}) {fracs[best]*100:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
