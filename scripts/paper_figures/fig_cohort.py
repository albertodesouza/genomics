#!/usr/bin/env python3
"""Method figure: which individuals the label retains, and how they fall across splits.

(a) the eight retained populations, one bar each, coloured by the class the population
    is mapped to -- the figure that makes plain that the label is a relabelling of
    continental ancestry rather than a noisy correlate of it;
(b) the *family-aware split*: the train/validation/test partition used to train and
    evaluate the CNN classifier throughout the rest of the paper. It keeps every
    individual and assigns whole families to a single split, so a parent's genome is
    never in training while the child's is scored.

The founders subset used by the variant-association test of Sec. 3.5 is a separate,
unrelated procedure and is not shown here -- it is reported as a number in the text
where that test is described.

Population membership comes from the cohort's own sample table; the split counts are
read off the processed tensor shards of one arm, since every arm shares one split.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from _common import RES, savefig
from fig_app_class_mean_signal import arm_cache_dir

STRONG_POPS = ["YRI", "ESN", "LWK", "MSL", "GWD"]
WEAK_POPS = ["FIN", "CEU", "GBR"]
C_STRONG, C_WEAK = "#2b2118", "#d9a441"
REF_ARM = "SLC24A5"


def split_counts(gene=REF_ARM):
    view = arm_cache_dir(gene)
    out = {}
    for split in ("train", "val", "test"):
        c = {}
        for shard in sorted(view.glob(f"{split}_data_shard_*.pt")):
            for _, y in torch.load(shard, map_location="cpu", weights_only=False):
                c[int(y)] = c.get(int(y), 0) + 1
        out[split] = c
    return out


def main():
    info = pd.read_csv(RES / "gwas" / "sample_info.tsv", sep="\t")
    pops = STRONG_POPS + WEAK_POPS
    n = info.groupby("POP").size()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 3.5),
                                   gridspec_kw={"width_ratios": [2.0, 1.0]})

    x = np.arange(len(pops))
    cols = [C_STRONG] * len(STRONG_POPS) + [C_WEAK] * len(WEAK_POPS)
    ax1.bar(x, [n.get(p, 0) for p in pops], color=cols, width=0.74)
    for i, p in enumerate(pops):
        ax1.annotate(f"{n.get(p, 0)}", (i, n.get(p, 0)), ha="center", va="bottom",
                     fontsize=7.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(pops)
    for t, c in zip(ax1.get_xticklabels(), cols):
        t.set_color(c)
    ax1.set_ylabel("individuals")
    ax1.set_title(f"(a) Population $\\rightarrow$ class label ($n = {len(info)}$)",
                  fontsize=9)

    sc = split_counts()
    splits = ["train", "val", "test"]
    xs = np.arange(len(splits))
    w = 0.36
    for k, (col, lab) in enumerate(((C_STRONG, "strong (AFR)"), (C_WEAK, "weak (EUR)"))):
        v = [sc[s].get(k, 0) for s in splits]
        ax2.bar(xs + (k - 0.5) * w, v, width=w, color=col, label=lab)
        for i, y in enumerate(v):
            ax2.annotate(f"{y}", (xs[i] + (k - 0.5) * w, y), ha="center", va="bottom",
                         fontsize=7.5)
    ax2.set_xticks(xs)
    ax2.set_xticklabels([f"{s}\n({sum(sc[s].values())})" for s in splits])
    ax2.legend(fontsize=8, frameon=False)
    ax2.set_title("(b) Family-aware train/val/test split\n(used to train and evaluate the classifier)",
                  fontsize=9)

    for a in (ax1, ax2):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.grid(axis="y", lw=0.3, alpha=0.4)
        a.tick_params(labelsize=8)
    fig.tight_layout()
    savefig(fig, "cohort-splits.png")
    print({s: sc[s] for s in splits})
    print(f"founders total {int(info['FOUNDER'].sum())}")


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 9, "savefig.facecolor": "white",
                         "mathtext.default": "regular"})
    main()
