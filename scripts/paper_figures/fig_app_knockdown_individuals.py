#!/usr/bin/env python3
"""Per-individual knockdown arrows, one panel per gene of the top-10 model.

The backup paper carried this figure for a single gene (SLC24A5) and for the three
target-location strategies; here there is one anchoring method, so the axis freed up is
spent on the gene, giving one panel for each of the ten genes the probe is applied to.

Every panel reads the SAME classifier -- the top-10 multi-gene model -- so the baseline
log-odds of an individual is a property of the individual and not of the panel. That is
what makes the figure readable across panels: individuals are sorted once, by baseline
log-odds, and that order is reused in all ten panels, so the vertical position of a row
means the same thing everywhere and the panels differ only in the displacement field.
The x-axis is shared for the same reason.

An arrow is one individual: dot at the baseline log-odds, head at the value after the
100bp promoter scramble of that panel's gene, coloured by TRUE class. Colour here means
class, not gene class -- the gene's panel/control membership is in the panel title, so
the two uses of the scheme cannot be confused inside one figure.

Split into two floats of five panels because ten panels at a legible arrow density
overrun the page; the split follows the balanced-accuracy rank that defined the top-10,
so each float is a contiguous block of that ranking rather than an arbitrary half.

Reads the replay CSV only. No API calls, no GPU.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from _common import C_CTRL, C_PANEL, TOP10_GENES, gene_class, load_top10, savefig

# The old paper's convention, kept: dark for strong pigmentation, amber for weak.
C_STRONG = C_PANEL
C_WEAK = C_CTRL
LAB = {"strong pigmentation": "strong", "weak pigmentation": "weak"}


def prepare():
    d, _ = load_top10()
    d = d.copy()
    d["pre"] = d["baseline_strong_logit"] - d["baseline_weak_logit"]
    d["post"] = d["perturbed_strong_logit"] - d["perturbed_weak_logit"]
    # delta_log_odds is written by the replay; recomputing it from the logits is a cheap
    # check that the two columns describe the same quantity in the same direction.
    err = float(np.abs((d["post"] - d["pre"]) - d["delta_log_odds"]).max())
    assert err < 1e-4, f"post-pre disagrees with delta_log_odds by {err}"

    base = d.drop_duplicates("sample_id").set_index("sample_id")["pre"]
    assert d.groupby("sample_id")["pre"].nunique().max() == 1, "baseline varies within an individual"
    order = base.sort_values().index
    rank = {s: i for i, s in enumerate(order)}
    d["y"] = d["sample_id"].map(rank)
    return d


def one_figure(d, genes, name, xlim):
    fig, axes = plt.subplots(3, 2, figsize=(6.6, 7.1), sharex=True, sharey=True)
    flat = axes.ravel()
    n_ind = d["sample_id"].nunique()

    for ax, gene in zip(flat, genes):
        g = d[d["gene"] == gene]
        ax.axvline(0.0, color="0.55", lw=0.7, ls="--", zorder=1)
        for lab, colour in ((("strong pigmentation"), C_STRONG), (("weak pigmentation"), C_WEAK)):
            s = g[g["true_label"] == lab]
            if s.empty:
                continue
            ax.hlines(s["y"], np.minimum(s["pre"], s["post"]), np.maximum(s["pre"], s["post"]),
                      color=colour, lw=0.5, alpha=0.85, zorder=2)
            ax.scatter(s["pre"], s["y"], s=1.1, color=colour, lw=0, zorder=3)
            for sel, marker in ((s["post"] > s["pre"], ">"), (s["post"] < s["pre"], "<")):
                h = s[sel]
                if not h.empty:
                    ax.scatter(h["post"], h["y"], s=3.0, marker=marker, color=colour,
                               lw=0, zorder=4)

        flips = int(g["flipped"].sum())
        m = float(g["delta_log_odds"].mean())
        # HERC2 and MC1R sit at 1e-3; two decimals would print them as +0.00 and read as
        # a rounding artefact rather than as the near-exact nothing they are.
        mtxt = f"{m:+.3f}" if abs(m) < 0.01 else f"{m:+.2f}"
        ax.set_title(f"{gene} ({gene_class(gene).replace('control (draw 1)', 'control')})\n"
                     f"mean $\\Delta$ = {mtxt}, "
                     f"{flips}/{n_ind} flipped", fontsize=8.5, pad=4)
        ax.set_yticks([])
        ax.tick_params(labelsize=7.5)
        ax.set_xlim(*xlim)
        ax.set_ylim(-3, n_ind + 2)

    # The sixth slot carries the legend, so no panel has to give up space to it.
    leg = flat[len(genes)]
    leg.axis("off")
    leg.legend(handles=[
        Line2D([], [], color=C_STRONG, lw=1.4, marker="<", ms=4, label="strong pigmentation"),
        Line2D([], [], color=C_WEAK, lw=1.4, marker="<", ms=4, label="weak pigmentation"),
        Line2D([], [], color="0.55", lw=0.7, ls="--", label="decision boundary"),
    ], loc="center", frameon=False, fontsize=8.5,
        title="one arrow per test individual\n(dot: before, head: after)", title_fontsize=8.5)
    for ax in flat[len(genes) + 1:]:
        ax.axis("off")

    # sharex hides the tick labels on every row but the last, and the last row is only
    # half full -- so the lowest USED panel of each column gets the axis back, rather
    # than the panel that happens to sit in the bottom row.
    nrow = axes.shape[0]
    for col in range(axes.shape[1]):
        used = [r for r in range(nrow) if r * axes.shape[1] + col < len(genes)]
        if not used:
            continue
        ax = axes[max(used)][col]
        ax.tick_params(labelbottom=True)
        ax.set_xlabel("strong/weak log-odds", fontsize=8.5)
    axes[nrow // 2][0].set_ylabel(
        f"individuals, sorted by baseline log-odds ($n = {n_ind}$)", fontsize=8.5)
    fig.tight_layout(h_pad=1.4, w_pad=1.0)
    savefig(fig, name)
    plt.close(fig)


def main() -> int:
    d = prepare()
    lo = float(min(d["pre"].min(), d["post"].min()))
    hi = float(max(d["pre"].max(), d["post"].max()))
    pad = 0.04 * (hi - lo)
    xlim = (lo - pad, hi + pad)

    one_figure(d, TOP10_GENES[:5], "app-knockdown-individuals-a.png", xlim)
    one_figure(d, TOP10_GENES[5:], "app-knockdown-individuals-b.png", xlim)

    print("\nshared x-range: [%.2f, %.2f]" % xlim)
    for gene in TOP10_GENES:
        g = d[d["gene"] == gene]
        print(f"  {gene:9s} {gene_class(gene):17s} mean {g['delta_log_odds'].mean():+8.3f}  "
              f"median {g['delta_log_odds'].median():+8.3f}  "
              f"flipped {int(g['flipped'].sum()):3d}  "
              f"toward strong {int((g['delta_log_odds'] > 0).sum()):3d}/{len(g)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
