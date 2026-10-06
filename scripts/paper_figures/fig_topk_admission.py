#!/usr/bin/env python3
"""Controls admitted at the thresholds that recover panel genes, and the test on it.

The head-to-head statistic the paper used to report was Delta AUC, which is unusable at
this panel size for a reason that is arithmetic rather than statistical: with the
association arm already at AUC 0.724, the largest attainable Delta AUC is 0.276, and the
paired permutation's own 95th percentile sits at 0.271 -- 98% of the way to the ceiling.
Only a near-perfect ranking (AUC >= 0.995) could ever have cleared it, so a
non-significant result there carries no information either way.

This module tests the claim that is actually localised where the two methods differ: at
the threshold that recovers k panel genes, how many control genes has the method already
admitted?

    C_M(k) = # controls method M ranks at or above its k-th best panel gene
             (strict count + 1/2 per tie, the midrank convention auc_and_p uses)

C(k)/n0 is the false-positive rate at the k-th panel gene, so {C(k)} is the ROC read
horizontally, and the integral of the claim over the top of the ranking is the partial
AUC over FPR <= c:

    pAUC(c) = (1 / (c * n1)) * sum_k max(0, c - C(k)/n0)

the average TPR over FPR in [0, c]. Two properties make this the right summary rather
than a hand-picked one: at c = 1 it reduces exactly to the Mann-Whitney AUC (so the cap
sweep interpolates between "the top only" and the AUC the paper already reports), and its
own ceiling is 1 regardless of how the association arm performs, which is what removes
the Delta AUC ceiling.

Two nulls, because they answer different questions and disagree here:

  (a) label reassignment -- resample which 9 of the 42 genes are the panel. The null the
      rest of the paper already uses for every panel-against-control AUC. Asks whether the
      panel's top-of-ranking advantage under the classifier exceeds what an arbitrary
      gene set would show. Primary.
  (b) method flip -- keep panel membership fixed (it is literature curation, not a draw)
      and exchange the two methods' ranks per gene, re-ranking afterwards so both
      pseudo-rankings stay proper permutations. Stricter: it treats the two methods'
      disagreements as noise that could have gone either way. Reported as sensitivity.

Writes the JSON, the two appendix tables and the figure. Reads finished artefacts only.
"""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import (C_CTRL, C_DRAW3, C_PANEL, PANEL, RES, TABDIR, load_final_table,
                     load_gwas, savefig)

# Pre-declared primary cap: 0.15 of 33 controls is exactly 5 admitted false positives,
# the scale a gene screen is actually read at. It is deliberately not the cap that
# maximises the result (0.10 does, see the sweep) -- picking that one is what makes a cap
# choice unfalsifiable.
CAP = 0.15
ESTS = ["S_sidak", "S_meanchi2"]
EST_LABEL = {"S_sidak": "Šidák", "S_meanchi2": "mean $\\chi^2$"}


def admitted(S, mask, n1, n0):
    """C(k) for k = 1..n1. S: (B,n) scores, mask: (B,n) or (n,) panel mask -> (B,n1)."""
    S = np.atleast_2d(S)
    M = np.broadcast_to(mask, S.shape)
    out = np.empty((S.shape[0], n1))
    for b0 in range(0, S.shape[0], 20_000):
        sl = slice(b0, min(b0 + 20_000, S.shape[0]))
        s, mm = S[sl], M[sl]
        # k-th best panel score, descending; controls sorted so the comparison broadcasts
        p = np.sort(np.where(mm, s, -np.inf), axis=1)[:, ::-1][:, :n1]
        c = np.sort(np.where(~mm, s, np.inf), axis=1)[:, :n0]
        out[sl] = ((c[:, :, None] > p[:, None, :]).sum(axis=1)
                   + 0.5 * (c[:, :, None] == p[:, None, :]).sum(axis=1))
    return out


def pauc(C, cap, n1, n0):
    """Average TPR over FPR in [0, cap], from the C(k) curve. cap = 1 gives the AUC."""
    return np.maximum(0.0, cap - C / n0).sum(axis=-1) / (cap * n1)


def nulls(r_cnn, r_gw, mask, B, seed=13):
    """(label-reassignment, method-flip) replicate draws, each as (C_cnn, C_gwas)."""
    n, n1 = len(mask), int(mask.sum())
    n0 = n - n1
    rng = np.random.default_rng(seed)
    idx = np.argsort(rng.random((B, n)), axis=1)[:, :n1]
    MA = np.zeros((B, n), bool)
    np.put_along_axis(MA, idx, True, axis=1)
    lab = (admitted(np.broadcast_to(r_cnn, (B, n)), MA, n1, n0),
           admitted(np.broadcast_to(r_gw, (B, n)), MA, n1, n0))
    F = rng.random((B, n)) < 0.5
    SC = pd.DataFrame(np.where(F, r_gw, r_cnn)).rank(axis=1).to_numpy()
    SG = pd.DataFrame(np.where(F, r_cnn, r_gw)).rank(axis=1).to_numpy()
    flip = (admitted(SC, mask, n1, n0), admitted(SG, mask, n1, n0))
    return lab, flip


def pvals(d_obs, Cc, Cg, cap, n1, n0):
    """(one-sided, two-sided) p and the 95% critical value, for one null's draws."""
    dn = pauc(Cc, cap, n1, n0) - pauc(Cg, cap, n1, n0)
    return (float((dn >= d_obs - 1e-12).mean()),
            float((np.abs(dn) >= abs(d_obs) - 1e-12).mean()),
            float(np.quantile(np.abs(dn), 0.95)))


def main(args):
    df = load_final_table()
    cnn = dict(zip(df["gene"], df["bal_acc"]))
    gw = load_gwas(args.arm, str(args.n))["gene_scores"][args.res]
    genes = sorted(set(cnn) & set(gw))
    mask = np.array([g in PANEL for g in genes], bool)
    n, n1 = len(genes), int(mask.sum())
    n0 = n - n1
    r_cnn = pd.Series([cnn[g] for g in genes]).rank().to_numpy()
    print(f"n = {n} ({n1} panel, {n0} control); B = {args.b:,}; primary cap = {CAP}")

    out = {"n": n, "n_panel": n1, "n_control": n0, "cap_primary": CAP, "B": args.b,
           "arm": args.arm, "resolution": args.res, "estimators": {}}
    curves = {}
    for est in ESTS:
        r_gw = pd.Series([gw[g][est] for g in genes]).rank().to_numpy()
        Cc = admitted(r_cnn, mask, n1, n0)[0]
        Cg = admitted(r_gw, mask, n1, n0)[0]
        curves[est] = (Cc, Cg)
        lab, flip = nulls(r_cnn, r_gw, mask, args.b, args.seed)

        sweep = []
        for cap in [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00]:
            pc, pg = pauc(Cc, cap, n1, n0), pauc(Cg, cap, n1, n0)
            d = float(pc - pg)
            p1a, p2a, ca = pvals(d, *lab, cap, n1, n0)
            p1b, p2b, cb = pvals(d, *flip, cap, n1, n0)
            sweep.append({"cap": cap, "pauc_cnn": float(pc), "pauc_gwas": float(pg),
                          "d": d, "ceiling": float(1.0 - pg),
                          "p_label_1s": p1a, "p_label_2s": p2a, "crit_label": ca,
                          "p_flip_1s": p1b, "p_flip_2s": p2b, "crit_flip": cb})
            print(f"  {est:11s} cap {cap:.2f}  d {d:+.4f}  ceil {1-pg:.4f}  "
                  f"p2(label) {p2a:.4f}  p1(flip) {p1b:.4f}")

        # Leave-one-gene-out at the primary cap: no single gene may carry the result.
        jack = []
        for g in genes:
            keep = [x for x in genes if x != g]
            mk = np.array([x in PANEL for x in keep], bool)
            k1, k0 = int(mk.sum()), len(keep) - int(mk.sum())
            rc = pd.Series([cnn[x] for x in keep]).rank().to_numpy()
            rgw = pd.Series([gw[x][est] for x in keep]).rank().to_numpy()
            dj = float(pauc(admitted(rc, mk, k1, k0)[0], CAP, k1, k0)
                       - pauc(admitted(rgw, mk, k1, k0)[0], CAP, k1, k0))
            lj, _ = nulls(rc, rgw, mk, args.jb, args.seed)
            jack.append({"dropped": g, "panel": bool(g in PANEL), "d": dj,
                         "p_label_2s": pvals(dj, *lj, CAP, k1, k0)[1]})
        worst = max(jack, key=lambda r: r["p_label_2s"])
        print(f"  {est:11s} jackknife worst: drop {worst['dropped']} "
              f"d {worst['d']:+.4f} p2 {worst['p_label_2s']:.4f}")

        out["estimators"][est] = {
            "C_cnn": Cc.tolist(), "C_gwas": Cg.tolist(),
            "auc_cnn": float(pauc(Cc, 1.0, n1, n0)), "auc_gwas": float(pauc(Cg, 1.0, n1, n0)),
            "sweep": sweep, "jackknife": jack,
            "jackknife_worst": worst,
            "primary": next(s for s in sweep if s["cap"] == CAP)}

    if args.legacy31:
        out["legacy_22_control"] = legacy31(args)

    dst = RES / "gwas_ranking" / "topk_admission.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dst}")
    tables(out)
    figure(out, curves, genes, mask, n0)
    return out


def legacy31(args):
    """The same pre-declared test on the superseded 22-control build, for the record.

    The discussion cites these to show that growing the control set sharpened the
    comparison that can be tested and left the one that cannot unmoved, so the numbers
    have to be regenerable rather than remembered.
    """
    df = pd.read_csv(RES / "poolmax_final_table_31arms_backup.csv")
    cnn = dict(zip(df["gene"], df["bal_acc"]))
    gw = json.loads((RES / "gwas_ranking" / "gene_ranking_uncorrected_33.json")
                    .read_text())["gene_scores"][args.res]
    genes = sorted(set(cnn) & set(gw))
    mask = np.array([g in PANEL for g in genes], bool)
    n1, n0 = int(mask.sum()), len(genes) - int(mask.sum())
    rc = pd.Series([cnn[g] for g in genes]).rank().to_numpy()
    res = {"n": len(genes), "n_control": n0}
    for est in ESTS:
        rg = pd.Series([gw[g][est] for g in genes]).rank().to_numpy()
        Cc = admitted(rc, mask, n1, n0)[0]
        Cg = admitted(rg, mask, n1, n0)[0]
        d = float(pauc(Cc, CAP, n1, n0) - pauc(Cg, CAP, n1, n0))
        lab, flip = nulls(rc, rg, mask, args.b, args.seed)
        p1a, p2a, ca = pvals(d, *lab, CAP, n1, n0)
        p1b, p2b, _ = pvals(d, *flip, CAP, n1, n0)
        res[est] = {"d": d, "crit_label": ca, "p_label_2s": p2a, "p_flip_1s": p1b,
                    "dauc": float(pauc(Cc, 1.0, n1, n0) - pauc(Cg, 1.0, n1, n0))}
        print(f"  [22 ctrl] {est:11s} d {d:+.4f}  p2(label) {p2a:.4f}  "
              f"p1(flip) {p1b:.4f}  dAUC {res[est]['dauc']:+.4f}")
    return res


def tables(out):
    n0 = out["n_control"]
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"\multicolumn{2}{l}{} & \multicolumn{2}{c}{$\mathrm{pAUC}(c)$} & & "
             r"\multicolumn{2}{c}{two-sided $p$} \\",
             r"\cmidrule(lr){3-4}\cmidrule(lr){6-7}",
             r"estimator & $c$ & classifier & assoc.\ & $d$ & label & method \\",
             r"\midrule"]
    for est in ESTS:
        lines.append(rf"\multicolumn{{7}}{{l}}{{\textit{{{EST_LABEL[est]}}}}} \\")
        for s in out["estimators"][est]["sweep"]:
            mark = r"\ \textbf{*}" if s["cap"] == out["cap_primary"] else ""
            lines.append(f"& {s['cap']:.2f}{mark} & {s['pauc_cnn']:.4f} & "
                         f"{s['pauc_gwas']:.4f} & {s['d']:+.4f} & "
                         f"{s['p_label_2s']:.4f} & {s['p_flip_2s']:.4f} \\\\")
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}"]
    p = TABDIR / "tab-pauc-sweep.tex"
    p.write_text("\n".join(lines) + "\n")
    print(f"wrote {p}")

    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
             r"$k$ & $C_{\mathrm{cls}}(k)$ & $C_{\text{\v S}}(k)$ & "
             r"$C_{\bar{\chi}^2}(k)$ & $\Delta C$ vs {\v S}id{\'a}k \\", r"\midrule"]
    a, b = out["estimators"]["S_sidak"], out["estimators"]["S_meanchi2"]
    for k in range(out["n_panel"]):
        lines.append(f"{k+1} & {a['C_cnn'][k]:.1f} & {a['C_gwas'][k]:.1f} & "
                     f"{b['C_gwas'][k]:.1f} & {a['C_gwas'][k]-a['C_cnn'][k]:+.1f} \\\\")
    lines += [r"\midrule",
              rf"\multicolumn{{5}}{{l}}{{out of {n0} control genes}} \\",
              r"\bottomrule", r"\end{tabular}"]
    p = TABDIR / "tab-admitted.tex"
    p.write_text("\n".join(lines) + "\n")
    print(f"wrote {p}")


def figure(out, curves, genes, mask, n0):
    n1 = out["n_panel"]
    cap = out["cap_primary"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    tp = np.arange(1, n1 + 1) / n1
    style = {"S_sidak": (C_CTRL, "-"), "S_meanchi2": (C_DRAW3, "--")}

    def step(ax, C, colour, ls, label):
        x, y = [0.0], [0.0]
        for f, t in zip(C / n0, tp):
            x += [f, f]
            y += [y[-1], t]
        x += [1.0]
        y += [y[-1]]
        ax.plot(x, y, color=colour, ls=ls, lw=1.8, label=label, zorder=3)

    ax1.axvspan(0, cap, color="0.92", zorder=0)
    ax1.plot([0, 1], [0, 1], color="0.75", lw=0.8, ls=":", zorder=1)
    step(ax1, curves["S_sidak"][0], C_PANEL, "-", "classifier (bal. acc.)")
    for est in ESTS:
        step(ax1, curves[est][1], *style[est], f"GWAS, {EST_LABEL[est]}")
    ax1.set_xlim(0, 0.45)
    ax1.set_ylim(0, 1.02)
    ax1.set_xlabel("fraction of control genes admitted (FPR)")
    ax1.set_ylabel("fraction of panel genes recovered")
    ax1.set_title(f"ROC, top region ($c = {cap}$ shaded)", fontsize=10)
    ax1.legend(fontsize=8, loc="lower right", frameon=False)
    ax1.grid(lw=0.25, alpha=0.35)

    ks = np.arange(1, n1 + 1)
    ax2.plot(ks, curves["S_sidak"][0], "o-", color=C_PANEL, lw=1.8, ms=4,
             label="classifier")
    for est in ESTS:
        c, ls = style[est]
        ax2.plot(ks, curves[est][1], "s" + ls, color=c, lw=1.8, ms=4,
                 label=f"GWAS, {EST_LABEL[est]}")
    ax2.axhline(cap * n0, color="0.6", lw=0.8, ls=":")
    ax2.text(1.05, cap * n0 + 1.0, f"$c = {cap}$ ({cap*n0:.0f} of {n0})",
             fontsize=7.5, color="0.4")
    ax2.set_xticks(ks)
    ax2.set_xlabel("panel genes recovered, $k$")
    ax2.set_ylabel("control genes admitted, $C(k)$")
    ax2.set_title("the same curve, read as a screen", fontsize=10)
    ax2.legend(fontsize=8, loc="upper left", frameon=False)
    ax2.grid(lw=0.25, alpha=0.35)

    s = out["estimators"]["S_sidak"]["primary"]
    fig.text(0.5, -0.04,
             f"pAUC($c = {cap}$): classifier {s['pauc_cnn']:.4f} against "
             f"{s['pauc_gwas']:.4f}, $d = {s['d']:+.4f}$ "
             f"(ceiling {s['ceiling']:.4f}), two-sided $p = {s['p_label_2s']:.4f}$",
             ha="center", va="top", fontsize=8.5, family="monospace",
             bbox=dict(fc="white", ec="0.7", lw=0.6, pad=5))
    fig.tight_layout()
    savefig(fig, "topk-admission.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="uncorrected")
    ap.add_argument("--n", type=int, default=44)
    ap.add_argument("--res", default="window")
    ap.add_argument("--b", type=int, default=200_000)
    ap.add_argument("--jb", type=int, default=50_000, help="B for each jackknife refit")
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--legacy31", action="store_true", default=True,
                    help="also run the pre-declared test on the 22-control build")
    main(ap.parse_args())
