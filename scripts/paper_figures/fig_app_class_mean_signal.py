#!/usr/bin/env python3
"""Appendix figure A1: class-mean predicted RNA-seq per gene, with exons marked.

For each arm, the mean normalised RNA-seq the classifier reads, averaged over the strong
individuals and over the weak individuals separately, with a +/-1 SD band and the exons
of the gene's MANE Select transcript shaded. This is the quantity the classifier
discriminates on, so it is the direct answer to "is the signal separable at all", and the
exon overlay says whether the separable part of the window is the transcript or something
else in the 32,768 bp crop.

Reads the on-disk processed tensor shards only. Nothing is trained, predicted or
downloaded beyond the GENCODE annotation, which is cached under notebooks/.cache.

Usage: python fig_app_class_mean_signal.py [--split test] [--genes TYR OCA2 ...]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import yaml

from _common import C_PANEL, DATASET, DRAW1, DRAW2, DRAW3, PANEL, REPO, savefig

sys.path.insert(0, str(REPO / "notebooks"))
from genotype_cnn_alignment_deeplift_summary.annotations import (  # noqa: E402
    build_transcript_extractors, gene_boundary_positions, gene_exon_local_boxes,
    load_gtf, window_genomic_axis,
)

CLASS_NAMES = ["strong pigmentation", "weak pigmentation"]   # index order of the label
C_STRONG, C_WEAK = C_PANEL, "#d9a441"
GTF_CACHE = REPO / "notebooks" / ".cache"


def arm_cache_dir(gene: str) -> Path:
    """The processed-tensor cache of this gene's own poolmax arm.

    Resolved through the config the sweep wrote next to the checkpoint, the same way the
    knockdown replay resolves it, so a figure can never pick up another arm's tensors.
    """
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if not hits:
        raise FileNotFoundError(f"no poolmax config for {gene}")
    if len(hits) > 1:
        raise RuntimeError(f"{gene}: {len(hits)} poolmax configs, ambiguous")
    cfg = yaml.safe_load(hits[0].read_text())
    genes = cfg["dataset_input"]["genes_to_use"]
    if genes != [gene]:
        raise RuntimeError(f"{gene}: config gene list is {genes}")
    d = REPO / cfg["dataset_input"]["processed_cache_dir"] / "datasets"
    views = [p for p in d.iterdir() if p.is_dir() and (p / ".cache_complete").exists()]
    if len(views) != 1:
        raise RuntimeError(f"{gene}: {len(views)} complete views under {d}")
    return views[0]


def class_curves(view: Path, split: str):
    """(mean, sd, n) of the normalised signal per class, pooled over both haplotypes.

    The two haplotype channels share one fitted divisor and are exchangeable by
    construction (Section: Representation), so pooling them is the same reduction the
    normalisation already assumes rather than an extra averaging step.
    """
    sums, sqs, ns = {}, {}, {}
    for shard in sorted(view.glob(f"{split}_data_shard_*.pt")):
        for x, y in torch.load(shard, map_location="cpu", weights_only=False):
            k = int(y)
            v = x.reshape(x.shape[0], -1).double().numpy()     # (haplotype, position)
            sums[k] = sums.get(k, 0.0) + v.sum(axis=0)
            sqs[k] = sqs.get(k, 0.0) + (v ** 2).sum(axis=0)
            ns[k] = ns.get(k, 0) + v.shape[0]
    out = {}
    for k in sorted(sums):
        m = sums[k] / ns[k]
        var = np.maximum(sqs[k] / ns[k] - m ** 2, 0.0)
        out[k] = (m, np.sqrt(var), ns[k] // 2)                 # n individuals
    return out


def smooth(v, w=64):
    if w <= 1:
        return v
    k = np.ones(w) / w
    return np.convolve(v, k, mode="same")


def main(args):
    gtf = load_gtf(GTF_CACHE)
    _, _, tx_mane, tx_coding, gene_id_map = build_transcript_extractors(gtf)

    if args.genes:
        groups = [("selected", args.genes)]
    else:
        # One figure per group: 29+ arms at three columns is ten rows, which overflows a
        # single page and LaTeX then refuses to place the float.
        groups = [("panel", PANEL), ("control-draw1", DRAW1), ("control-draw2", DRAW2),
                 ("control-draw3", DRAW3)]
    for tag, genes in groups:
        one_figure(tag, genes, args, tx_mane, tx_coding, gene_id_map)


def one_figure(tag, genes, args, tx_mane, tx_coding, gene_id_map):
    ncol = 3
    avail = []
    for g in genes:
        try:
            avail.append((g, arm_cache_dir(g)))
        except (FileNotFoundError, RuntimeError) as e:
            print(f"skip {g}: {e}")
    if not avail:
        print(f"no arms available for {tag}")
        return

    nrow = -(-len(avail) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 2.05 * nrow),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)

    for i, (gene, view) in enumerate(avail):
        ax = axes[i // ncol][i % ncol]
        ax.set_visible(True)
        curves = class_curves(view, args.split)
        chrom, w0, L = window_genomic_axis(DATASET, gene, 32768)
        boxes, strand = gene_exon_local_boxes(gene, chrom, w0, L, gene_id_map,
                                              tx_mane, tx_coding)
        tss, tes = gene_boundary_positions(boxes, strand)

        for s, e in boxes:
            ax.axvspan(s, e, color="tab:green", alpha=0.14, lw=0, zorder=0)
        if tss is not None:
            ax.axvline(tss, color="tab:green", lw=0.9, ls=":", alpha=0.9, zorder=1)
            ax.axvline(tes, color="tab:green", lw=0.9, ls="-.", alpha=0.9, zorder=1)

        n_txt = []
        for k, col in ((0, C_STRONG), (1, C_WEAK)):
            if k not in curves:
                continue
            m, sd, n = curves[k]
            ms, ss = smooth(m), smooth(sd)
            x = np.arange(len(ms))
            ax.fill_between(x, ms - ss, ms + ss, color=col, alpha=0.28, lw=0, zorder=2)
            ax.plot(x, ms, color=col, lw=1.0, zorder=3)
            n_txt.append(f"{CLASS_NAMES[k].split()[0]} $n={n}$")

        ax.set_xlim(0, L)
        ax.set_title(f"{gene} ({'panel' if gene in PANEL else 'control'}, "
                     f"{chrom} {strand}) -- " + ", ".join(n_txt), fontsize=8)
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if i % ncol == 0:
            ax.set_ylabel("normalised\nRNA-seq", fontsize=7.5)
        if i // ncol == nrow - 1:
            ax.set_xlabel("position in 32,768 bp crop", fontsize=7.5)

    # Legend and titles are placed in figure fractions, so their reserved band has to
    # shrink as rows are added or a 1-row figure loses its x labels to the legend.
    pad = min(0.24, 0.62 / nrow)
    h = [plt.Line2D([], [], color=C_STRONG, lw=2),
         plt.Line2D([], [], color=C_WEAK, lw=2),
         plt.Rectangle((0, 0), 1, 1, fc="tab:green", alpha=0.25)]
    fig.legend(h, ["strong pigmentation (AFR)", "weak pigmentation (EUR)",
                   "exon, MANE Select transcript"],
               loc="lower center", ncol=3, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, 0.002))
    fig.suptitle(f"Class-mean predicted RNA-seq of CL:1000458 on the gene's own strand, "
                 f"{args.split} split, {tag} arms "
                 f"(shaded band = $\\pm 1$ SD across individuals)",
                 fontsize=10.5, y=1.0 - 0.02 * pad)
    fig.tight_layout(rect=(0, pad * 0.30, 1, 1.0 - pad * 0.10))
    savefig(fig, f"app-class-mean-signal-{tag}-{args.split}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="test")
    ap.add_argument("--genes", nargs="*", default=None)
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    main(ap.parse_args())
