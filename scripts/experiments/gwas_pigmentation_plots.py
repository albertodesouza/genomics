#!/usr/bin/env python3
"""Manhattan and Q-Q figures for the pigmentation GWAS (gwas_pigmentation.py).

Reads the exported summary statistics under
results/genotype_based_predictor/gwas/ and writes two figures sized for the
5.5 in text block of the ICLR style file.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/gwas_pigmentation_plots.py \
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
GWAS = REPO_ROOT / "results" / "genotype_based_predictor" / "gwas"
GW = 5e-8
WIDTH = 5.5

# Categorical slots 1-3 of the validated reference palette; these three are the
# documented all-pairs-safe subset, which is what a scatter needs.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8c8b85"
# Alternating window tones carry grouping, not identity: the axis labels carry that.
TONE_A, TONE_B = "#8fa3b8", "#41556b"
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


def load(name: str) -> pd.DataFrame:
    d = pd.read_csv(GWAS / f"{name}.sumstats.tsv.gz", sep="\t", dtype={"#CHROM": str})
    d = d[(d.ERRCODE.astype(str) == ".") & d.P.notna()].copy()
    d["nlp"] = -np.log10(np.clip(d.P.values, 1e-300, 1.0))
    return d


def tidy(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


# ------------------------------------------------------------------ manhattan


def manhattan(summary: dict, outdir: Path) -> Path:
    chr15 = load("chr15_uncorrected")
    panel = load("panel_uncorrected")
    sites = pd.read_csv(GWAS / "panel_sites.tsv", sep="\t")
    genes = list(summary["build"]["panel"]["windows"])
    thr = -np.log10(GW)
    # Points reach ~96; the band above 100 is left clear for the variant labels.
    ymax = 124.0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(WIDTH, 4.9),
                                   gridspec_kw={"hspace": 0.62})

    def cloud(ax, x, y):
        """Below the threshold in the recessive tone, above it in the dark one."""
        sig = y >= thr
        ax.scatter(x[~sig], y[~sig], s=0.9, c=TONE_A, lw=0, rasterized=True, zorder=2)
        ax.scatter(x[sig], y[sig], s=0.9, c=TONE_B, lw=0, rasterized=True, zorder=3)
        ax.axhline(thr, color=ORANGE, lw=1.0, ls="--", zorder=4)

    # (a) chromosome 15 --------------------------------------------------
    wins = [v for g, v in summary["build"]["panel"]["windows"].items() if v["chrom"] == "chr15"]
    for v in wins:
        ax1.axvspan(v["start"] / 1e6, v["end"] / 1e6, color="#f0c987", alpha=0.85, lw=0, zorder=1)
    cloud(ax1, chr15.POS.values / 1e6, chr15.nlp.values)
    a = summary["arms"]["chr15"]["uncorrected"]
    ax1.text(0.5, 1.03, f"{a['n_genomewide']:,} of {a['n_tests']:,} SNPs "
             f"({a['frac_genomewide']*100:.0f}%) above the line", transform=ax1.transAxes,
             ha="center", va="bottom", fontsize=7, color=INK2)
    ax1.text(101, thr + 2, r"$p = 5\times10^{-8}$", ha="right", va="bottom",
             color=ORANGE, fontsize=6.8, zorder=5)
    ax1.text(28.07, 120, "OCA2 / HERC2", ha="center", va="top", fontsize=6.2,
             color="#8a6520", zorder=5)
    ax1.text(48.13, 120, "SLC24A5", ha="center", va="top", fontsize=6.2,
             color="#8a6520", zorder=5)

    label_x = {"rs12913832": 33.0, "rs1426654": 61.0}
    for rsid, k in summary["chr15_known_variants"].items():
        if not k.get("present"):
            continue
        x, y = k["pos"] / 1e6, -np.log10(k["p"])
        ax1.scatter([x], [y], s=30, facecolor="none", edgecolor=INK, lw=1.1, zorder=6)
        ax1.annotate(f"{rsid} — rank {k['rank']:,}",
                     xy=(x, y), xytext=(label_x[rsid], 106),
                     fontsize=6.6, color=INK, ha="left", va="center", zorder=7,
                     arrowprops=dict(arrowstyle="-", lw=0.6, color=INK, shrinkA=1, shrinkB=3))

    ax1.set_title("(a) all of chromosome 15", loc="left", pad=14)
    ax1.set_xlabel("chromosome 15 position (Mb)", labelpad=2)

    # (b) the eleven input windows --------------------------------------
    merged = panel.merge(sites, on=["POS"], how="left")
    gap = 160000.0
    offsets, ticks, cursor = {}, [], 0.0
    for i, g in enumerate(genes):
        w = summary["build"]["panel"]["windows"][g]
        span = w["end"] - w["start"]
        offsets[g] = cursor
        if i % 2 == 1:
            ax2.axvspan(cursor / 1e6, (cursor + span) / 1e6, color="#f2f2ef", lw=0, zorder=0)
        ticks.append((cursor + span / 2) / 1e6)
        cursor += span + gap
    x_all, y_all = [], []
    for g in genes:
        sub = merged[merged.GENE == g]
        w = summary["build"]["panel"]["windows"][g]
        x_all.append((offsets[g] + (sub.POS.values - w["start"])) / 1e6)
        y_all.append(sub.nlp.values)
    cloud(ax2, np.concatenate(x_all), np.concatenate(y_all))
    b = summary["arms"]["panel"]["uncorrected"]
    ax2.text(0.5, 1.03, f"{b['n_genomewide']:,} of {b['n_tests']:,} SNPs "
             f"({b['frac_genomewide']*100:.0f}%) above the line",
             transform=ax2.transAxes, ha="center", va="bottom",
             fontsize=7, color=INK2)
    ax2.set_title("(b) the eleven 524 kbp windows the classifier is shown", loc="left", pad=14)
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(genes, rotation=45, ha="right", fontsize=6.8)
    ax2.set_xlim(-0.1, (cursor - gap) / 1e6 + 0.1)
    ax2.tick_params(axis="x", length=0)

    for ax, xlim in ((ax1, None), (ax2, True)):
        ax.set_ylabel(r"$-\log_{10} p$", labelpad=2)
        ax.set_ylim(0, ymax)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.grid(axis="y", color=GRID, lw=0.5)
        tidy(ax)
    ax1.set_xlim(chr15.POS.min() / 1e6 - 1, chr15.POS.max() / 1e6 + 1)

    out = outdir / "gwas-manhattan.png"
    fig.savefig(out)
    plt.close(fig)
    return out


# ------------------------------------------------------------------ q-q


def qq_points(p: np.ndarray, keep_head: int = 4000, thin: int = 40):
    p = np.sort(np.clip(np.asarray(p, float), 1e-300, 1.0))
    n = p.size
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p)
    idx = np.concatenate([np.arange(min(keep_head, n)),
                          np.arange(min(keep_head, n), n, thin)])
    return exp[idx], obs[idx], n


def qq(summary: dict, outdir: Path) -> Path:
    from scipy import stats

    fig, ax = plt.subplots(figsize=(WIDTH, 3.3))

    # 95% concentration band under the null, from the beta order statistics.
    n_ref = summary["arms"]["chr15"]["uncorrected"]["n_tests"]
    k = np.unique(np.round(np.geomspace(1, n_ref, 500)).astype(int))
    lo = stats.beta.ppf(0.025, k, n_ref - k + 1)
    hi = stats.beta.ppf(0.975, k, n_ref - k + 1)
    ax.fill_between(-np.log10((k - 0.5) / n_ref), -np.log10(hi), -np.log10(lo),
                    color=GRID, lw=0, zorder=1, label="95% null band")

    series = [
        ("chr15_uncorrected", "chr15", "uncorrected", BLUE, "-",
         "pigmentation label, uncorrected"),
        ("chr15_founders", "chr15", "founders", ORANGE, "-",
         "pigmentation label, 840 pedigree founders"),
        ("chr15_randomlabel", "chr15", "randomlabel", AQUA, "-",
         "ancestry-independent label (control)"),
    ]
    for name, arm, tag, colour, ls, label in series:
        d = load(name)
        exp, obs, n = qq_points(d.P.values)
        lam = summary["arms"][arm][tag]["lambda_gc"]
        ax.plot(exp, obs, ls, color=colour, lw=1.4, zorder=3,
                label=rf"{label}   $\lambda_{{\mathrm{{GC}}}} = {lam:.2f}$")

    lim = max(ax.get_xlim()[1], 1.0)
    top = ax.get_ylim()[1]
    ax.plot([0, lim], [0, lim], color=MUTED, lw=0.8, ls=":", zorder=2,
            label="expected under the null")
    ax.axhline(-np.log10(GW), color=MUTED, lw=0.7, ls="--", zorder=2)
    ax.text(0.02, -np.log10(GW) + 1.5, r"$p = 5\times10^{-8}$", transform=ax.get_yaxis_transform(),
            fontsize=6.8, color=INK2, va="bottom")

    pc = summary["arms"]["chr15"]["pc10"]
    ax.text(0.975, 0.40,
            f"with ten genotype PCs as covariates: 0 of {pc['n_tests']:,} Firth fits converged\n"
            r"(PC1 separates the label completely, $|r| = "
            f"{summary['arms']['chr15']['structure']['pc1_point_biserial_r']:.3f}$)",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.8, color=INK2)

    ax.set_xlabel(r"expected $-\log_{10} p$ under the null")
    ax.set_ylabel(r"observed $-\log_{10} p$")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, top)
    ax.grid(color=GRID, lw=0.5)
    tidy(ax)
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none",
              framealpha=0.92, borderaxespad=0.3, handlelength=1.6, labelspacing=0.35)

    out = outdir / "gwas-qq.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path,
                    default=Path("/home/breno/I2CA/paper-knockdown-probe/figures"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads((GWAS / "gwas_pigmentation.json").read_text())
    for f in (manhattan(summary, args.outdir), qq(summary, args.outdir)):
        print("wrote", f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
