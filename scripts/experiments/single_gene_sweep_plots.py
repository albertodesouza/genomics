#!/usr/bin/env python3
"""Panel-vs-control figures for the single-gene sweep.

Two axes, because they disagree. The eleven control genes are a random draw (the
`random_11_1` panel, screened only to exclude any GO annotation under pigmentation,
melanin biosynthesis or melanocyte differentiation), so both comparisons are tests on a
randomly drawn control panel rather than on a hand-picked one.

  (a) accuracy  -- what one gene's window alone supports. Separates (p = 0.033).
  (b) |Delta|   -- the same 100 bp promoter scramble through that arm's own checkpoint.
                   Does NOT separate (p = 0.106).
  (c) the two against each other, so a reader can see that the gene topping the response
      axis is a control and that the two axes rank the genes differently.

Reads results/genotype_based_predictor/single_gene_sweep_report.json only.

Usage:
  python3 scripts/experiments/single_gene_sweep_plots.py \
      --outdir /home/breno/I2CA/paper-knockdown-probe/figures
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO_ROOT = Path("/home/breno/I2CA/genomics")
SRC = REPO_ROOT / "results" / "genotype_based_predictor" / "single_gene_sweep_report.json"
WIDTH = 5.5

# Categorical slots 1-3 of the validated reference palette; the documented all-pairs-safe
# subset, which is what a two-class scatter needs.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8c8b85"
GRID = "#e3e3df"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.labelcolor": INK, "text.color": INK,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

CLS = {"panel": BLUE, "control": ORANGE}


def strip(ax, df, col, ylabel, pval, baseline=None, logy=False):
    """One quantity, both classes, every gene drawn. Jitter is deterministic: the reader
    should get the same picture on a redraw, and a moved point is a different claim."""
    rng = np.random.default_rng(13)
    for i, cls in enumerate(("panel", "control")):
        g = df[df["class"] == cls]
        x = i + rng.uniform(-0.14, 0.14, len(g))
        ax.scatter(x, g[col], s=26, facecolor=CLS[cls], edgecolor="white",
                   linewidth=0.5, zorder=3, clip_on=False)
        med = g[col].median()
        ax.plot([i - 0.28, i + 0.28], [med, med], color=INK, lw=1.4, zorder=4)
        ax.annotate(f"{med:.3f}", (i + 0.30, med), fontsize=6.5, color=INK2,
                    va="center", ha="left")
    if baseline is not None:
        ax.axhline(baseline, color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=1)
        ax.annotate("majority class", (-0.46, baseline), fontsize=6.5, color=MUTED,
                    va="top", ha="left",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9))
    if logy:
        ax.set_yscale("log")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"panel\n(n={(df['class']=='panel').sum()})",
                        f"control\n(n={(df['class']=='control').sum()})"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Mann–Whitney $p$ = {pval:.3f}", color=INK2, fontsize=7.5, pad=5)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path,
                    default=Path("/home/breno/I2CA/paper-knockdown-probe/figures"))
    args = ap.parse_args()
    if not SRC.exists():
        raise SystemExit(f"ABORT: {SRC} missing; run single_gene_sweep_report.py first.")
    rep = json.loads(SRC.read_text())
    df = pd.DataFrame(rep["arms"])
    df = df[df["class"].isin(("panel", "control"))].copy()
    if df.empty:
        raise SystemExit("ABORT: no classified arms in the report.")
    baseline = float(rep["baseline"])
    p_acc = float(rep["mannwhitney"]["acc"]["p"])
    p_del = float(rep["mannwhitney"]["abs_delta_single"]["p"])

    fig = plt.figure(figsize=(WIDTH, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], wspace=0.42, hspace=0.42)
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    strip(ax_a, df, "acc", "test accuracy", p_acc, baseline=baseline)
    strip(ax_b, df, "abs_delta_single", r"$|\Delta|$, promoter scramble", p_del)
    ax_b.set_ylabel(r"$|\Delta|$")

    # (c) the two axes against each other, every gene labelled.
    for cls in ("panel", "control"):
        g = df[df["class"] == cls]
        ax_c.scatter(g.acc, g.abs_delta_single, s=30, facecolor=CLS[cls],
                     edgecolor="white", linewidth=0.5, zorder=3, label=cls)
    ax_c.axvline(baseline, color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=1)
    # Labels are placed greedily against already-placed boxes: with 22 genes and two
    # near-exact accuracy ties (TYRP1/PSMC4 at .9136, TCHH/TPM2 at .7160) a fixed offset
    # overplots names, and an unreadable name is worse than a leader line.
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    placed = []
    order = df.sort_values("abs_delta_single", ascending=False)
    cands = [(4, 3), (4, -7), (-4, 3), (-4, -7), (4, 9), (-4, 9), (0, 11), (0, -13)]
    for _, r in order.iterrows():
        for dx, dy in cands:
            t = ax_c.annotate(r.gene, (r.acc, r.abs_delta_single), fontsize=5.6, color=INK2,
                              xytext=(dx, dy), textcoords="offset points",
                              ha="left" if dx >= 0 else "right", va="center")
            bb = t.get_window_extent(renderer=rend).expanded(1.06, 1.24)
            if not any(bb.overlaps(q) for q in placed):
                placed.append(bb)
                if abs(dx) > 4 or abs(dy) > 7:      # far enough to need a leader
                    ax_c.annotate("", xy=(r.acc, r.abs_delta_single),
                                  xytext=(dx * 0.55, dy * 0.55), textcoords="offset points",
                                  arrowprops=dict(arrowstyle="-", lw=0.35, color=MUTED,
                                                  shrinkA=0, shrinkB=0))
                break
            t.remove()
        else:
            t = ax_c.annotate(r.gene, (r.acc, r.abs_delta_single), fontsize=5.6,
                              color=MUTED, xytext=(4, 3), textcoords="offset points")
            placed.append(t.get_window_extent(renderer=rend))
    ax_c.set_xlabel("test accuracy, one gene alone")
    ax_c.margins(x=0.06, y=0.09)
    ax_c.set_ylabel(r"$|\Delta|$")
    ax_c.grid(color=GRID, lw=0.5, zorder=0)
    ax_c.set_axisbelow(True)
    for side in ("top", "right"):
        ax_c.spines[side].set_visible(False)
    ax_c.legend(handles=[Line2D([], [], marker="o", ls="", markersize=4.5,
                                markerfacecolor=CLS[c], markeredgecolor="white",
                                label={"panel": "pigmentation panel",
                                       "control": "random control"}[c])
                         for c in ("panel", "control")],
                loc="upper left", frameon=False, handletextpad=0.3, borderpad=0.1)

    fig.canvas.draw()
    for ax, lab in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
        bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
        fig.text(bb.x0, bb.y1 + 0.014, lab, fontsize=9, fontweight="bold", color=INK,
                 ha="left", va="bottom")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / "single-gene-panel-vs-control.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")

    # What the figure is asserting, printed so the caption cannot drift from the data.
    for cls in ("panel", "control"):
        g = df[df["class"] == cls]
        print(f"  {cls:8} n={len(g):2}  acc median={g.acc.median():.4f} "
              f"[{g.acc.min():.4f},{g.acc.max():.4f}]  "
              f"|D| median={g.abs_delta_single.median():.4f} "
              f"[{g.abs_delta_single.min():.4f},{g.abs_delta_single.max():.4f}]")
    top = df.loc[df.abs_delta_single.idxmax()]
    print(f"  highest |D| overall: {top.gene} ({top['class']}) at {top.abs_delta_single:.3f}")
    below = df[df.acc < baseline]
    if len(below):
        print(f"  below the majority baseline: {', '.join(below.gene)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
