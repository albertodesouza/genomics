#!/usr/bin/env python3
"""Appendix figure A2: the knockdown before and after, on both readouts.

(a) Predicted exonic expression of the gene, summed over its MANE Select exons and over
    both haplotypes, before the scramble and after it. This is the panel that says
    whether the intervention was a knockdown at all: an arrow pointing down is a loss of
    predicted transcription, an arrow pointing up is an intervention that failed on its
    own terms, and the frozen model decides which of the two happened.
(b) The classifier's decision on the same individuals, as mean log-odds of strong over
    weak before and after, drawn separately for the strong and the weak individuals
    because a displacement that moves both classes the same way changes no decision.

Reads the knockdown replay CSVs only.
"""
from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import (C_CTRL, C_PANEL, CONTROL, KD_DIR, PANEL, TOP10_GENES, savefig,
                     load_top10)

C_STRONG, C_WEAK = "#2b2118", "#d9a441"


def load(genes, method="biology_tss"):
    frames = []
    for g in genes:
        p = KD_DIR / f"single_{g.lower()}.csv"
        if not p.exists():
            print(f"skip {g}: no replay CSV")
            continue
        d = pd.read_csv(p)
        d = d[d["method"] == method]
        if d.empty:
            print(f"skip {g}: no {method} rows")
            continue
        frames.append(d)
    if not frames:
        raise SystemExit("no replay CSVs")
    return pd.concat(frames, ignore_index=True)


def main(args):
    if args.source == "top10":
        # The arm the knockdown is reported from: one classifier over the top-10 genes,
        # one gene scrambled at a time.
        df, _ = load_top10(args.method)
        genes_in_order = [g for g in TOP10_GENES]
    else:
        df = load(PANEL + CONTROL, args.method)
        genes_in_order = PANEL + CONTROL
    signed = "expr_mod" in df.columns and df["expr_mod"].notna().any()
    if not signed:
        print("NOTE: expr_mod/expr_log2fc absent from these CSVs (written by an earlier "
              "version of the replay); panel (a) falls back to the magnitude of the "
              "exonic change, which cannot show direction. Re-run the replay to fix.")

    genes = [g for g in genes_in_order if g in set(df["gene"])]
    order = (df.groupby("gene")["delta_log_odds"].mean().abs()
             .reindex(genes).sort_values(ascending=False).index.tolist())
    y = {g: i for i, g in enumerate(order)}

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.4, 0.30 * len(order) + 2.2),
        gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.30})

    # ---- (a) predicted exonic expression, before -> after
    # Panel (a) carries no class contrast of its own -- panel/control is already on the
    # shared y tick labels -- so its arrows are neutral grey and the two class colours are
    # reserved for panel (b), where they mean strong/weak individuals.
    A_COL = "#4d4d4d"
    lim_a, lim_b = [], []
    for g in order:
        d = df[df["gene"] == g]
        col = A_COL
        pre = d["expr_baseline"].mean()
        post = d["expr_mod"].mean() if signed else pre - d["delta_expr"].mean()
        axA.annotate("", xy=(post, y[g]), xytext=(pre, y[g]),
                     arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3,
                                     shrinkA=0, shrinkB=0, mutation_scale=9))
        axA.scatter([pre], [y[g]], s=22, facecolors="white", edgecolors=col, lw=1.1,
                    zorder=3)
        lim_a += [pre, post]
        lfc = np.log2(max(post, 1e-9) / max(pre, 1e-9))
        axA.annotate(f"{lfc:+.2f}", (1.006, y[g]), xycoords=("axes fraction", "data"),
                     fontsize=7, va="center",
                     color="tab:red" if lfc > 0 else "#4d4d4d")
    axA.set_xscale("log")
    # Explicit limits: only the "before" circles are drawn with scatter, so autoscaling
    # would leave every arrowhead outside the pre range invisible.
    lo, hi = min(v for v in lim_a if v > 0), max(lim_a)
    axA.set_xlim(lo / 3.0, hi * 3.0)
    axA.set_yticks(range(len(order)))
    axA.set_yticklabels(order, fontsize=8)
    for t, g in zip(axA.get_yticklabels(), order):
        t.set_color(C_PANEL if g in PANEL else C_CTRL)
    axA.set_ylim(len(order) - 0.5, -0.5)
    axA.set_xlabel("predicted exonic RNA-seq, summed over exons and haplotypes\n"
                   "(open circle = before, arrowhead = after; log scale)", fontsize=8)
    if signed:
        axA.set_title("(a) Did the scramble reduce predicted expression?\n"
                      r"annotation = $\log_2$ fold change; red = increased", fontsize=9)
    else:
        # Without expr_mod the only available quantity is the MAGNITUDE of the exonic
        # change, so "after" is forced below "before" by construction and every fold
        # change is negative whether or not the intervention actually reduced anything.
        # Say so on the figure, so a stale render cannot be read as a result.
        axA.set_title("(a) Exonic change, DIRECTION NOT AVAILABLE\n"
                      "these CSVs carry only the magnitude; re-run the replay",
                      fontsize=9, color="tab:red")

    # ---- (b) decision, before -> after, by true class
    for g in order:
        d = df[df["gene"] == g]
        for cls, col, dy in (("strong pigmentation", C_STRONG, -0.17),
                             ("weak pigmentation", C_WEAK, +0.17)):
            s = d[d["true_label"] == cls]
            if s.empty:
                continue
            pre = (s["baseline_strong_logit"] - s["baseline_weak_logit"]).mean()
            post = (s["perturbed_strong_logit"] - s["perturbed_weak_logit"]).mean()
            axB.annotate("", xy=(post, y[g] + dy), xytext=(pre, y[g] + dy),
                         arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3,
                                         shrinkA=0, shrinkB=0, mutation_scale=9))
            axB.scatter([pre], [y[g] + dy], s=20, facecolors="white", edgecolors=col,
                        lw=1.1, zorder=3)
            lim_b += [pre, post]
    pad = 0.06 * (max(lim_b) - min(lim_b))
    axB.set_xlim(min(lim_b) - pad, max(lim_b) + pad)
    axB.axvline(0, color="tab:red", ls="--", lw=1.0)
    axB.annotate("decision boundary", (0, -0.9), fontsize=7.5, color="tab:red",
                 ha="center", va="bottom")
    axB.set_yticks(range(len(order)))
    axB.set_yticklabels([])
    axB.set_ylim(len(order) - 0.5, -0.5)
    axB.set_xlabel(r"$\log$ odds of strong over weak"
                   "\n(open circle = before, arrowhead = after)", fontsize=8)
    axB.set_title("(b) Where did the decision move?\n"
                  "one arrow per true class, mean over its individuals", fontsize=9)

    for a in (axA, axB):
        a.grid(axis="y", lw=0.25, alpha=0.35)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.tick_params(labelsize=8)

    h = [plt.Line2D([], [], color="#4d4d4d", lw=2),
         plt.Line2D([], [], color=C_STRONG, lw=2), plt.Line2D([], [], color=C_WEAK, lw=2),
         plt.Line2D([], [], color=C_PANEL, lw=0, marker="s", ms=6),
         plt.Line2D([], [], color=C_CTRL, lw=0, marker="s", ms=6)]
    fig.legend(h, ["expression, all arms (a)",
                   "strong individuals (b)", "weak individuals (b)",
                   "gene label: panel", "gene label: control"],
               loc="lower center", ncol=5, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    n_ind = df.groupby("gene")["sample_id"].nunique().max()
    what = ("top-10 multi-gene model, one gene scrambled at a time"
            if args.source == "top10" else "single-gene models")
    fig.suptitle(f"Promoter knockdown on the {what}, {args.method} anchoring: "
                 f"{len(order)} genes, {n_ind} held-out individuals each, "
                 f"ordered by $|\\Delta|$", fontsize=10.5)
    fig.tight_layout(rect=(0, 0.035, 1, 0.985))
    savefig(fig, f"app-knockdown-prepost-{args.source}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", default="biology_tss")
    ap.add_argument("--source", default="top10", choices=["top10", "single"])
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    main(ap.parse_args())
