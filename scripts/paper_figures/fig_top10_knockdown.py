#!/usr/bin/env python3
"""The knockdown arm: one classifier over the top-10 genes, one gene scrambled at a time.

Produces
  tables/tab-top10-perturbation.tex   per-gene readouts, with the AFR/EUR split
  figures/top10-delta-violin.png      the response by class, panel against control

The unit is the gene, as it must be: Delta is a mean over the 162 held-out individuals
and the question is which gene the model relies on. Seven of the ten are panel genes and
three are controls, by accuracy and not by class, which is what makes a specificity test
possible on this arm at all -- and what limits it, since three points do not support a
density.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from _common import (C_CTRL, C_PANEL, PANEL, TABDIR, load_top10, load_top10_null,
                     savefig, top10_selection_drift)

MIN_FOR_DENSITY = 5     # below this a violin body is a kernel over almost nothing


def thousands(x, n=0):
    return f"{x:,.{n}f}".replace(",", r"\,")


def table(reg, rows, null):
    """Per-gene readouts. Ordered by |Delta|, which is the readout of reliance."""
    lfc = rows.groupby("gene")["expr_log2fc"].mean()
    # The null band is per gene, not pooled: the old paper showed the null's width varies
    # by orders of magnitude between genes, so a pooled threshold would call a gene with a
    # tight null significant and one with a wide null not, at the same |Delta|.
    nb = null.groupby("gene")["delta_log_odds"].apply(lambda v: v.abs().max())
    r = reg.reindex(reg["delta"].abs().sort_values(ascending=False).index)
    lines = [r"\begin{tabular}{llrrrrrrrr}", r"\toprule",
             r"gene & class & $\Delta_{\mathrm{in}}$ & $\Delta_{\mathrm{expr}}$ & "
             r"$\log_2$FC & $\Delta$ & $\Delta_{\mathrm{AFR}}$ & "
             r"$\Delta_{\mathrm{EUR}}$ & unif. & null max \\",
             r"\midrule"]
    for _, x in r.iterrows():
        g = x["gene"]
        f = lfc.get(g, float("nan"))
        lines.append(
            f"{g} & {x['class']} & {thousands(x['delta_in'])} & "
            f"{thousands(x['delta_expr'])} & "
            + (r"\textcolor{red}{" + f"{f:+.2f}" + "}" if f > 0 else f"{f:+.2f}") + " & "
            f"{x['delta']:+.3f} $\\pm$ {x['delta_sd']:.2f} & "
            f"{x['delta_afr']:+.3f} & {x['delta_eur']:+.3f} & "
            f"{x['uniformidade']:.3f} & "
            # Bold where the response fails to clear its own null: the reader should be
            # able to find those three rows without doing the comparison by eye.
            + (f"{nb[g]:.3f}" if abs(x["delta"]) > nb[g]
               else r"\textbf{" + f"{nb[g]:.3f}" + "}")
            + " \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    out = TABDIR / "tab-top10-perturbation.tex"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


def violin(reg):
    pan = reg[reg["class"] == "panel"]["delta"]
    ctl = reg[reg["class"] == "control"]["delta"]
    u = stats.mannwhitneyu(pan.abs(), ctl.abs(), alternative="greater")
    auc = u.statistic / (len(pan) * len(ctl))

    meds: dict[str, float] = {}
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for pos, v, col, lab in ((0, pan, C_PANEL, "panel"), (1, ctl, C_CTRL, "control")):
        y = np.log10(v.abs().to_numpy())
        if len(v) >= MIN_FOR_DENSITY:
            parts = ax.violinplot([y], positions=[pos], widths=0.72,
                                  showextrema=False, showmedians=False)
            b = parts["bodies"][0]
            b.set_facecolor(col); b.set_alpha(0.22)
            b.set_edgecolor(col); b.set_linewidth(1.0)
        else:
            # Three points. A kernel density over them would be an artefact of the
            # bandwidth, so the points and the median are drawn and nothing else.
            ax.annotate(f"$n = {len(v)}$: no density drawn", (pos, 0.02),
                        xycoords=("data", "axes fraction"), ha="center", va="bottom",
                        fontsize=7.5, color=col, style="italic")
        med = float(v.abs().median())
        ax.hlines(np.log10(med), pos - 0.24, pos + 0.24, color=col, lw=2.0, zorder=4)
        meds[lab] = med

    # Labels are placed in a fixed column to the right of each group and pushed apart
    # to a minimum spacing, with a leader to the point: at ten genes over four decades a
    # jittered label collides with its neighbour and with the median annotation.
    ylo, yhi = np.log10(0.0008), np.log10(14.0)
    ax.set_ylim(ylo, yhi)
    min_gap = 0.055 * (yhi - ylo)
    for pos, sub, col in ((0, reg[reg["class"] == "panel"], C_PANEL),
                          (1, reg[reg["class"] == "control"], C_CTRL)):
        pts = sorted(((np.log10(abs(x["delta"])), x["gene"]) for _, x in sub.iterrows()))
        ys = [y for y, _ in pts]
        placed = []
        for y in ys:                        # single upward pass is enough for n <= 7
            placed.append(y if not placed else max(y, placed[-1] + min_gap))
        shift = max(0.0, placed[-1] - yhi + min_gap) if placed else 0.0
        placed = [y - shift for y in placed]
        for (y, gene), ly in zip(pts, placed):
            ax.scatter(pos, y, s=40, c=col, lw=0.7, edgecolors="white", zorder=3)
            ax.annotate(gene, xy=(pos + 0.06, y), xytext=(pos + 0.36, ly),
                        fontsize=7.5, va="center", ha="left", color=col,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color=col,
                                        shrinkA=2, shrinkB=2))

    t = [0.001, 0.01, 0.1, 1, 10]
    ax.set_yticks([np.log10(v) for v in t])
    ax.set_yticklabels([f"{v:g}" for v in t])
    ax.set_ylabel(r"$|\Delta|$, shift in log odds of strong over weak")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"panel\n($n = {len(pan)}$)", f"control\n($n = {len(ctl)}$)"])
    ax.set_xlim(-0.35, 1.95)
    ax.set_title("Knockdown response on the top-10 model, one gene at a time\n"
                 f"median $|\\Delta|$: panel {meds['panel']:.3f}, "
                 f"control {meds['control']:.3f}"
                 "    |    "
                 f"Mann--Whitney one-sided $p = {u.pvalue:.4f}$, AUC $= {auc:.4f}$"
                 + (" (controls respond more)" if auc < 0.5 else ""), fontsize=9.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", lw=0.3, alpha=0.4)
    fig.tight_layout()
    savefig(fig, "top10-delta-violin.png")

    print(f"panel  n={len(pan)} median |delta| {pan.abs().median():.4f} "
          f"range [{pan.abs().min():.4f}, {pan.abs().max():.4f}]")
    print(f"ctrl   n={len(ctl)} median |delta| {ctl.abs().median():.4f} "
          f"range [{ctl.abs().min():.4f}, {ctl.abs().max():.4f}]")
    print(f"MWU one-sided p = {u.pvalue:.4f}  U = {u.statistic:.1f}  AUC = {auc:.4f}")
    order = reg.reindex(reg["delta"].abs().sort_values(ascending=False).index)
    print("ranking by |Delta|: " + ", ".join(
        f"{g}{'*' if c == 'control' else ''}"
        for g, c in zip(order["gene"], order["class"])) + "   (* = control)")


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    rows, reg = load_top10()
    d = top10_selection_drift()
    print(f"selection: trained on the top 10 of the arms available then; over all "
          f"{d['n_arms_now']} arms now it would be {d['enters']} in, {d['leaves']} out")
    null = load_top10_null()
    table(reg, rows, null)
    violin(reg)
