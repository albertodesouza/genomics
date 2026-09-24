#!/usr/bin/env python3
"""Response against the two things that could be driving it, on the top-10 arm.

The backup paper carried one panel of this family: response against DELIVERED
perturbation, one point per gene, which is panel (b) here regenerated on the multi-gene
model. Panel (a) is the one the biology actually rests on -- response against the change
in the gene's own predicted expression -- and the figure exists to put the two side by
side, because the null in (a) only means something next to the correlation in (b).

WHY THE X-AXIS OF (a) IS A FOLD CHANGE AND NOT RAW Delta_expr. Raw Delta_expr is a sum
of RNA-seq units over the gene's exons, so it is dominated by how much the gene is
expressed at baseline rather than by how hard it was hit: DDB1 records 41,823 against
MC1R's 1,091, but their baselines are 172,399 and 2,227, and in relative terms the
ordering reverses (-0.40 against -0.96 in log2). It is the same comparability defect that
raw Delta_in has. The log2 fold change is the cross-gene comparable form, and unlike
`delta_expr` (a per-position absolute sum, so non-negative by construction) it keeps the
sign that says whether the intervention was a knockdown at all.

Reads the replay CSV only. No API calls, no GPU.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import stats

from _common import C_CTRL, C_PANEL, gene_class, load_top10, savefig


def per_gene():
    d, _ = load_top10()
    g = d.groupby("gene").agg(
        delta=("delta_log_odds", "mean"),
        l2fc=("expr_log2fc", "mean"),
        dexpr=("delta_expr", "mean"),
        eb=("expr_baseline", "mean"),
        dinn=("delta_in_norm", "mean"),
    )
    g["absd"] = g["delta"].abs()
    g["panel"] = [gene_class(i) == "panel" for i in g.index]
    return g


def _p(v):
    """Two decimals unless that would print a significant p as 0.00."""
    return f"{v:.2f}" if v >= 0.01 else f"{v:.4f}"


def rho(x, y):
    r = stats.spearmanr(x, y)
    return r.statistic, r.pvalue


# Panel (b) is the panel carrying the delivery comparison, so its geometry is kept clean
# and the labels are placed directly rather than in a leader column.
def label_offsets(ax, xs, ys, names, pad=3.0):
    """Place each label at the first of eight compass offsets that collides with nothing.

    Replaces a hand-tuned offset table, which was correct for one gene set and raised a
    KeyError the moment the composition changed. Collision is tested in display space
    against the already-placed labels and against every data point, so the result adapts
    to whatever ten genes the selection returns.
    """
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    placed = []

    def clashes(bb):
        grown = bb.expanded(1.0 + pad / max(bb.width, 1), 1.0 + pad / max(bb.height, 1))
        if any(grown.overlaps(q) for q in placed):
            return True
        return any(grown.contains(px, py) for px, py in pts)

    cands = [(6, 5, "left", "bottom"), (6, -5, "left", "top"),
             (-6, 5, "right", "bottom"), (-6, -5, "right", "top"),
             (0, 9, "center", "bottom"), (0, -9, "center", "top"),
             (0, 19, "center", "bottom"), (0, -19, "center", "top"),
             (0, 29, "center", "bottom"), (0, -29, "center", "top")]
    for x, y, n in zip(xs, ys, names):
        chosen = None
        for dx, dy, ha, va in cands:
            t = ax.annotate(n, (x, y), textcoords="offset points", xytext=(dx, dy),
                            ha=ha, va=va, fontsize=7.5)
            bb = t.get_window_extent(rend)
            if not clashes(bb):
                chosen = (t, bb)
                break
            t.remove()
        if chosen is None:   # every slot taken: keep the last try rather than drop a label
            dx, dy, ha, va = cands[-1]
            t = ax.annotate(n, (x, y), textcoords="offset points", xytext=(dx, dy),
                            ha=ha, va=va, fontsize=7.5)
            chosen = (t, t.get_window_extent(rend))
        placed.append(chosen[1])


def leader_column(ax, xs, ys, names, x_col, min_gap_frac=0.062, logy=False):
    """Labels in a fixed vertical column at x_col, de-collided upward, each joined to its
    point by a thin leader. Needed because seven of the ten genes sit inside one small
    cluster, where any direct offset overlaps a neighbour. xy is the arrow target and
    xytext is where the text goes.
    """
    y0, y1 = ax.get_ylim()
    tf = (lambda v: np.log10(v)) if logy else (lambda v: np.asarray(v, dtype=float))
    inv = (lambda v: 10.0 ** v) if logy else (lambda v: v)
    y0, y1 = float(tf(y0)), float(tf(y1))
    gap = min_gap_frac * (y1 - y0)
    order = np.argsort(np.asarray(ys))
    xs, ys, names = np.asarray(xs), np.asarray(ys), np.asarray(names)
    slots = tf(ys[order]).astype(float).copy()
    for i in range(1, len(slots)):
        if slots[i] - slots[i - 1] < gap:
            slots[i] = slots[i - 1] + gap
    shift = max(0.0, slots[-1] - (y1 - 0.4 * gap))
    slots -= shift
    for k, i in enumerate(order):
        ax.annotate(names[i], xy=(xs[i], ys[i]), xytext=(x_col, inv(slots[k])),
                    fontsize=7.5, ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.45, color="0.6",
                                    shrinkA=1.5, shrinkB=2.5))


def main() -> int:
    g = per_gene()
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(6.9, 3.3))
    colours = [C_PANEL if p else C_CTRL for p in g["panel"]]

    # ---- (a) the expression readout -------------------------------------------------
    # Signed against signed: the quadrants are the biologically meaningful thing. A
    # pigmentation gene whose loss of function lightens skin should sit lower-left --
    # expression down, decision pushed toward weak.
    axa.axhline(0, color="0.6", lw=0.7, ls="--", zorder=1)
    axa.axvline(0, color="0.6", lw=0.7, ls="--", zorder=1)
    axa.scatter(g["l2fc"], g["delta"], s=34, c=colours, zorder=3, lw=0)
    # The two genes the scramble did not knock down get a ring rather than a shaded
    # half-axis: it names them instead of colouring a region no other gene occupies,
    # and it leaves the right-hand side free for the label column.
    failed = g[g["l2fc"] > 0]
    axa.scatter(failed["l2fc"], failed["delta"], s=104, facecolors="none",
                edgecolors="#c0392b", lw=1.1, zorder=4)
    # Computed from the data, never hardcoded: a fixed window silently dropped two genes
    # off the axis when the top-10 composition changed, and a clipped point in a
    # correlation figure is the kind of error that reads as a result.
    xlo, xhi = float(g["l2fc"].min()), float(g["l2fc"].max())
    ylo, yhi = float(g["delta"].min()), float(g["delta"].max())
    xpad, ypad = 0.12 * (xhi - xlo), 0.10 * (yhi - ylo)
    # Right margin widened to hold the label column; left/bottom just padded.
    axa.set_xlim(xlo - xpad, xhi + xpad + 0.62 * (xhi - xlo))
    axa.set_ylim(ylo - ypad, yhi + ypad)
    assert g["l2fc"].between(*axa.get_xlim()).all() and g["delta"].between(*axa.get_ylim()).all()
    leader_column(axa, g["l2fc"], g["delta"], g.index,
                  x_col=xhi + 0.35 * xpad + 0.10 * (xhi - xlo))
    rs, ps = rho(g["l2fc"], g["delta"])
    ra, pa = rho(g["l2fc"], g["absd"])
    axa.set_xlabel("predicted expression change, $\\log_2$ fold", fontsize=9)
    axa.set_ylabel("response $\\Delta$ (log-odds)", fontsize=9)
    axa.set_title(f"(a) against the gene's own expression\n"
                  f"signed $\\rho = {rs:+.3f}$ ($p = {_p(ps)}$); "
                  f"$|\\Delta|$: $\\rho = {ra:+.3f}$ ($p = {_p(pa)}$)", fontsize=8.5)

    # ---- (b) the delivery confound, the backup paper's figure ------------------------
    axb.scatter(g["dinn"], g["absd"], s=34, c=colours, zorder=3, lw=0)
    axb.set_xscale("log")
    axb.set_yscale("log")
    axb.set_xlim(100, 3.2e4)
    axb.set_ylim(2.2e-4, 60)
    label_offsets(axb, g["dinn"], g["absd"], g.index)
    rb, pb = rho(g["dinn"], g["absd"])
    lx, ly = np.log10(g["dinn"]), np.log10(g["absd"])
    sl, ic = np.polyfit(lx, ly, 1)
    xx = np.array([g["dinn"].min() * 0.55, g["dinn"].max() * 1.8])
    axb.plot(xx, 10 ** (ic + sl * np.log10(xx)), color="0.45", lw=0.9, ls="--", zorder=2)
    axb.set_xlabel("delivered perturbation $\\Delta_{\\mathrm{in}}$ (normalised)", fontsize=9)
    axb.set_ylabel("$|\\Delta|$ (log-odds)", fontsize=9)
    axb.set_title(f"(b) against how much was delivered\n"
                  f"$\\rho = {rb:+.3f}$ ($p = {pb:.4f}$), slope ${sl:.2f}$", fontsize=8.5)

    for ax in (axa, axb):
        ax.tick_params(labelsize=8)
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=5.5, color=C_PANEL, label="panel gene"),
        Line2D([], [], marker="o", ls="", ms=5.5, color=C_CTRL, label="control gene"),
        Line2D([], [], marker="o", ls="", ms=7.5, markerfacecolor="none",
               markeredgecolor="#c0392b", label="expression rose: intervention failed"),
    ], loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    savefig(fig, "app-expr-vs-response-top10.png")
    plt.close(fig)

    print(f"\n(a) log2FC vs signed Delta  rho={rs:+.3f} p={ps:.4f}")
    print(f"    log2FC vs |Delta|       rho={ra:+.3f} p={pa:.4f}")
    print(f"(b) norm Delta_in vs |Delta| rho={rb:+.3f} p={pb:.4f} slope={sl:.3f}")
    r, p = rho(g["dexpr"], g["absd"])
    print(f"    raw Delta_expr vs |Delta| rho={r:+.3f} p={p:.4f}  (not plotted: see header)")
    r, p = rho(g["eb"], g["dexpr"])
    print(f"    raw Delta_expr vs baseline expression rho={r:+.3f} p={p:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
