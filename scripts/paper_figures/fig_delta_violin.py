#!/usr/bin/env python3
"""Knockdown figure: the distribution of Delta, panel against control.

One point per arm -- Delta is already a mean over the held-out individuals, and the
claim this figure is read against is about genes, so the gene is the unit.

Both panels are drawn on a compressed vertical scale, because Delta spans four orders of
magnitude (TRHR at +40 against LACTB2 at -0.0005) and a linear axis would collapse every
body into a sliver at zero. The compression is applied to the values BEFORE the kernel
density is estimated, and the ticks are labelled in the original units, so the width of a
body is a density in the space it is drawn in rather than a linear density stretched by
the axis. Because Mann-Whitney is rank-based and both transforms are monotone, the
annotated test is identical on the transformed and the untransformed values: the
transform is presentational and load-bearing for nothing.

  (a) signed Delta, on an inverse-hyperbolic-sine scale, which is linear near zero and
      logarithmic in both tails and is therefore defined for the negative arms;
  (b) |Delta|, the readout of reliance, on a log scale.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from _common import (C_CTRL, C_DRAW2, C_PANEL, CONTROL, DRAW2, PANEL, load_final_table,
                     savefig)

# Arms named on the figure: the ones whose position is the result a reader would ask
# about, rather than an arbitrary top-k.
LABEL_ABOVE = 5.0


def asinh_fwd(x):
    return np.arcsinh(np.asarray(x, dtype=float))


def asinh_inv(y):
    return np.sinh(np.asarray(y, dtype=float))


def violin(ax, data_t, pos, col, width=0.70):
    parts = ax.violinplot([data_t], positions=[pos], widths=width,
                          showextrema=False, showmedians=False)
    b = parts["bodies"][0]
    b.set_facecolor(col)
    b.set_alpha(0.22)
    b.set_edgecolor(col)
    b.set_linewidth(1.0)


def strip(ax, df, pos, fwd, rng, label_above=None):
    """One marker per arm, jittered, with the two control draws in different shades."""
    for _, r in df.iterrows():
        c = C_PANEL if r["gene"] in PANEL else (C_DRAW2 if r["gene"] in DRAW2 else C_CTRL)
        x = pos + rng.uniform(-0.13, 0.13)
        yt = float(fwd(r["delta"] if label_above is None else abs(r["delta"])))
        ax.scatter(x, yt, s=34, c=c, lw=0.6, edgecolors="white", zorder=3)
        if abs(r["delta"]) >= LABEL_ABOVE:
            ax.annotate(r["gene"], (x + 0.16, yt), fontsize=7, va="center", color=c)


def ticks(ax, values, fwd, fmt):
    ax.set_yticks([float(fwd(v)) for v in values])
    ax.set_yticklabels([fmt(v) for v in values])


def main():
    df = load_final_table()
    pan = df[df["gene"].isin(PANEL)]
    ctl = df[df["gene"].isin(CONTROL)]
    rng = np.random.default_rng(13)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.4, 4.9))

    # ---- (a) signed Delta
    for pos, sub, col in ((0, pan, C_PANEL), (1, ctl, C_CTRL)):
        violin(axA, asinh_fwd(sub["delta"].to_numpy()), pos, col)
        strip(axA, sub, pos, asinh_fwd, rng)
        med = float(sub["delta"].median())
        axA.hlines(asinh_fwd(med), pos - 0.32, pos + 0.32, color=col, lw=2.0, zorder=4)
        # Above the bar and centred: to the left it would run into the y label at pos 0
        # and into the zero line's annotation at pos 1.
        axA.annotate(f"median {med:+.3f}", (pos, asinh_fwd(med)), fontsize=8,
                     va="bottom", ha="center", color=col, zorder=5,
                     xytext=(0, 3), textcoords="offset points",
                     bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0))
    axA.axhline(0.0, color="tab:red", ls="--", lw=1.0)
    axA.annotate("no displacement", (0.008, 0.0), xycoords=("axes fraction", "data"),
                 fontsize=7.5, color="tab:red", ha="left", va="bottom",
                 bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0))
    ticks(axA, [-3, -1, -0.3, 0, 0.3, 1, 3, 10, 30], asinh_fwd,
          lambda v: f"{v:+g}" if v else "0")
    u = stats.mannwhitneyu(pan["delta"], ctl["delta"], alternative="greater")
    axA.set_ylabel(r"$\Delta$, shift in log odds of strong over weak")
    axA.set_title(r"(a) signed $\Delta$" "\n"
                  f"Mann\u2013Whitney one-sided $p = {u.pvalue:.3f}$, "
                  f"AUC $= {u.statistic / (len(pan) * len(ctl)):.4f}$", fontsize=9)

    # ---- (b) |Delta|
    for pos, sub, col in ((0, pan, C_PANEL), (1, ctl, C_CTRL)):
        violin(axB, np.log10(sub["delta"].abs().to_numpy()), pos, col)
        strip(axB, sub, pos, lambda v: np.log10(np.abs(v)), rng, label_above=True)
        med = float(sub["delta"].abs().median())
        axB.hlines(np.log10(med), pos - 0.32, pos + 0.32, color=col, lw=2.0, zorder=4)
        axB.annotate(f"median {med:.3f}", (pos, np.log10(med)), fontsize=8,
                     va="bottom", ha="center", color=col, zorder=5,
                     xytext=(0, 3), textcoords="offset points",
                     bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0))
    ticks(axB, [0.001, 0.01, 0.1, 1, 10, 40], lambda v: np.log10(v),
          lambda v: f"{v:g}")
    u2 = stats.mannwhitneyu(pan["delta"].abs(), ctl["delta"].abs(), alternative="greater")
    axB.set_ylabel(r"$|\Delta|$")
    a2 = u2.statistic / (len(pan) * len(ctl))
    axB.set_title(r"(b) $|\Delta|$, the readout of reliance" "\n"
                  f"Mann\u2013Whitney one-sided $p = {u2.pvalue:.3f}$, "
                  f"AUC $= {a2:.4f}$"
                  + (" (controls respond more)" if a2 < 0.5 else ""), fontsize=9)

    for a in (axA, axB):
        a.set_xticks([0, 1])
        a.set_xticklabels([f"panel\n($n = {len(pan)}$)",
                           f"control\n($n = {len(ctl)}$, two random draws)"])
        a.set_xlim(-0.72, 1.62)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.grid(axis="y", lw=0.3, alpha=0.4)
        a.tick_params(labelsize=8)

    h = [plt.Line2D([], [], marker="o", ls="none", color=c, ms=6.5)
         for c in (C_PANEL, C_CTRL, C_DRAW2)]
    fig.legend(h, ["panel gene", "control, draw 1", "control, draw 2"],
               loc="lower center", ncol=3, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    savefig(fig, "delta-violin.png")

    print(f"panel  median delta {pan['delta'].median():+.4f}  "
          f"median |delta| {pan['delta'].abs().median():.4f}")
    print(f"ctrl   median delta {ctl['delta'].median():+.4f}  "
          f"median |delta| {ctl['delta'].abs().median():.4f}")
    print(f"signed: p = {u.pvalue:.4f}  AUC = {u.statistic/(len(pan)*len(ctl)):.4f}")
    print(f"|delta|: p = {u2.pvalue:.4f}  AUC = {u2.statistic/(len(pan)*len(ctl)):.4f}"
          "   (AUC < 0.5 means the controls respond MORE)")
    top = df.reindex(df["delta"].abs().sort_values(ascending=False).index).head(5)
    print("largest |delta|:")
    print(top[["gene", "classe", "delta", "delta_sd", "delta_in_norm"]].to_string(index=False))


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    main()
