#!/usr/bin/env python3
"""Figure 1a/1b: the variant-level association arm.

1a  chromosome 15 Manhattan, uncorrected arm, with the three panel windows it carries
    (OCA2, HERC2, SLC24A5) shaded. Shown for one chromosome rather than genome-wide
    because the question is not whether the canonical loci are detected but what the
    rest of the chromosome looks like while they are.
1b  gene score S_g = -log10 p_gene for every panel and control window, one marker per
    estimator. The best-SNP family and the aggregating family are drawn on separate
    axes because their scores are two orders of magnitude apart and a shared axis would
    crush the family that ranks the panel higher.

Usage: python fig_gwas.py [--n 33] [--arm uncorrected]
"""
from __future__ import annotations

import argparse
import gzip

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _common import (C_CTRL, C_DRAW2, C_DRAW3, C_PANEL, DRAW1, DRAW2, DRAW3, GENOMEWIDE,
                     PANEL, RES, crop_interval, load_gwas, savefig, window_interval)

CHR15_PANEL = ["OCA2", "HERC2", "SLC24A5"]


def read_chr15(path, keep_below=0.02, seed=13):
    """-log10 p per position. Points below `keep_below` of the axis are thinned, since a
    Manhattan plot of 4.5e5 variants is 95% an unreadable band at the bottom."""
    pos, val = [], []
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_pos, i_p = head.index("POS"), head.index("P")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p = f[i_p]
            if p in (".", "NA", ""):
                continue
            pos.append(int(f[i_pos]))
            val.append(float(p))
    pos = np.asarray(pos, dtype=np.int64)
    v = np.asarray(val, dtype=float)
    v = np.maximum(v, 1e-320)
    y = -np.log10(v)
    rng = np.random.default_rng(seed)
    lo = y < keep_below * y.max()
    thin = lo & (rng.random(len(y)) > 0.05)
    return pos[~thin], y[~thin], len(y)


def panel_a(args):
    src = RES / "gwas" / f"chr15_{args.arm}.sumstats.tsv.gz"
    pos, y, n_total = read_chr15(src)

    fig, ax = plt.subplots(figsize=(11, 3.4))
    # OCA2 and HERC2 are adjacent and their windows overlap by 234,773 bp, so their
    # labels are staggered rather than centred: the collision in the figure is the same
    # collision that makes assigning rs12913832 to a gene hard.
    lab_y = {"OCA2": 1.105, "HERC2": 1.015, "SLC24A5": 1.015}
    for g in CHR15_PANEL:
        _, ws, we = window_interval(g)
        _, cs, ce = crop_interval(g)
        ax.axvspan(ws / 1e6, we / 1e6, color="tab:blue", alpha=0.12, lw=0, zorder=0)
        ax.axvspan(cs / 1e6, ce / 1e6, color="tab:blue", alpha=0.35, lw=0, zorder=0)
        ax.annotate(g, ((ws + we) / 2e6, lab_y[g]), xycoords=("data", "axes fraction"),
                    ha="center", va="bottom", fontsize=8, color="tab:blue")
    ax.scatter(pos / 1e6, y, s=1.2, c="#555555", alpha=0.45, lw=0, rasterized=True)
    ax.axhline(-np.log10(GENOMEWIDE), color="tab:red", lw=1.0, ls="--")
    ax.annotate(r"genome-wide $5\times10^{-8}$", (0.004, -np.log10(GENOMEWIDE)),
                xycoords=("axes fraction", "data"), ha="left", va="bottom", fontsize=7.5,
                color="tab:red",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    ax.set_xlabel("Position on chromosome 15 (Mb)")
    ax.set_ylabel(r"$-\log_{10} p$")
    ax.set_xlim(pos.min() / 1e6, pos.max() / 1e6)
    ax.set_ylim(0, y.max() * 1.04)
    ax.set_title(f"Chromosome 15, {args.arm} arm: {n_total:,} variants tested. "
                 f"Light shading = 524,288 bp window, dark = 32,768 bp crop.",
                 fontsize=9, pad=24)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    savefig(fig, "gwas-chr15-manhattan.png")
    plt.close(fig)


def panel_b(args):
    d = load_gwas(args.arm, str(args.n))
    for res in ("window", "crop"):
        sc = d["gene_scores"][res]
        groups = [("panel", [g for g in PANEL if g in sc], C_PANEL),
                  ("control (draw 1)", [g for g in DRAW1 if g in sc], C_CTRL),
                  ("control (draw 2)", [g for g in DRAW2 if g in sc], C_DRAW2),
                  ("control (draw 3)", [g for g in DRAW3 if g in sc], C_DRAW3)]
        # Within each group genes are ordered by the primary estimator, so the figure is
        # read left to right as a ranking and the groups stay visually separate.
        order, xs, colors, boundaries = [], [], [], []
        x = 0.0
        for name, gs, col in groups:
            gs = sorted(gs, key=lambda g: -sc[g]["S_sidak"])
            for g in gs:
                order.append(g)
                xs.append(x)
                colors.append(col)
                x += 1.0
            boundaries.append((name, x - len(gs) / 2 - 0.5, col))
            x += 1.0
        xs = np.asarray(xs)

        fig, (ax, ax2) = plt.subplots(
            2, 1, figsize=(12, 5.6), sharex=True,
            gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.08})

        # Small x offsets: the three best-SNP estimators agree to within a fraction of a
        # marker on most windows, so drawn at the same x the later ones hide the earlier.
        ax.scatter(xs - 0.24, [sc[g]["S_min"] for g in order], s=46, facecolors="none",
                   edgecolors=colors, lw=1.0, zorder=2)
        ax.scatter(xs, [sc[g]["S_sidak"] for g in order], s=40, c=colors, marker="o",
                   zorder=3)
        ax.scatter(xs + 0.24, [sc[g]["S_gates"] for g in order], s=34, c=colors, marker="v",
                   zorder=3)
        ax.axhline(-np.log10(GENOMEWIDE), color="tab:red", lw=0.9, ls="--")
        ax.annotate(r"$5\times10^{-8}$", (0.998, -np.log10(GENOMEWIDE)),
                    xycoords=("axes fraction", "data"), ha="right", va="bottom",
                    fontsize=7, color="tab:red")
        ax.set_ylabel(r"$S_g = -\log_{10} p_{\rm gene}$" "\nbest-SNP family", fontsize=9)
        h = [plt.Line2D([], [], marker=m, ls="none", mfc=f, mec="k", color="k", ms=6.5)
             for m, f in (("o", "none"), ("o", "k"), ("v", "k"))]
        ax.legend(h, ["raw min $p$",
                      # Literal unicode: mathtext is off, so "x" reads as a variable.
                      "\u0160id\u00e1k \u00d7 Li & Ji $m_{eff}$ (primary)",
                      "GATES (extended Simes)"], fontsize=8, loc="center left", ncol=1,
                  frameon=False)

        ax2.scatter(xs, [sc[g]["S_meanchi2"] for g in order], s=44, c=colors, marker="s")
        ax2.set_ylabel("$S_g$\nmean $\\chi^2$", fontsize=9)
        ax2.set_xticks(xs)
        ax2.set_xticklabels(order, rotation=90, fontsize=7.5)
        for t, c in zip(ax2.get_xticklabels(), colors):
            t.set_color(c)
        for a in (ax, ax2):
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
            a.grid(axis="y", lw=0.3, alpha=0.4)
        for name, xc, col in boundaries:
            ax.annotate(name, (xc, 1.02), xycoords=("data", "axes fraction"),
                        ha="center", va="bottom", fontsize=9, color=col)
        ax2.set_xlim(-1, xs.max() + 1)
        lam = d["lambda_gc_windows"]
        fig.suptitle(f"Gene score over the {'524,288 bp window' if res == 'window' else '32,768 bp crop'}, "
                     f"{args.arm} arm ($\\lambda_{{GC}} = {lam:.1f}$), "
                     f"{len(order)} genes", fontsize=9.5, y=0.995)
        savefig(fig, f"gwas-gene-scores-{res}.png")
        plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=33)
    ap.add_argument("--arm", default="uncorrected")
    ap.add_argument("--only", choices=["a", "b"], default=None)
    a = ap.parse_args()
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.facecolor": "white",
                         "text.usetex": False, "mathtext.default": "regular"})
    if a.only in (None, "a"):
        panel_a(a)
    if a.only in (None, "b"):
        panel_b(a)
