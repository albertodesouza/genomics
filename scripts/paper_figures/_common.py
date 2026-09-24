"""Shared loading for the paper figures.

Every figure reads finished artefacts only -- the GWAS ranking JSON, the poolmax final
table and the knockdown replay CSVs -- so a figure can be regenerated without re-running
any arm. Nothing here trains, predicts or bills the AlphaGenome API.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/home/breno/I2CA/genomics")
PAPER = Path("/home/breno/I2CA/phenotype-gene-discovery-alphagenome")
FIGDIR = PAPER / "figures"
TABDIR = PAPER / "tables"
RES = REPO / "results" / "genotype_based_predictor"
DATASET = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")

FINAL_TABLE = RES / "poolmax_final_table.csv"
KD_DIR = RES / "knockout_bulk" / "poolmax"

# The top-10 multi-gene arm: the genes with the highest balanced accuracy among the
# single-gene arms, gathered into one classifier, with the knockdown applied one gene at
# a time to that model. This is the arm the knockdown results are reported from.
#
# Selected on the VALIDATION split, over all 31 single-gene arms. The superseded arm
# (results/.../knockout_bulk/top10/) ranked on the TEST split and over only the 20 arms
# that existed before the second control draw; its files are kept for the comparison in
# the text but nothing in the paper reads them.
TOP10_CSV = RES / "knockout_bulk" / "top10val" / "top10val_knockdown.csv"
TOP10_NULL_CSV = RES / "knockout_bulk" / "top10val" / "top10val_null.csv"
TOP10_REGIME = RES / "top10val_regime.json"
TOP10_GENES = ["SLC24A5", "SLC45A2", "SPRED2", "HERC2", "DDB1",
               "LRRC36", "FYB1", "MFSD12", "FRA10AC1", "OR51S1"]

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
DRAW1 = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
         "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]
DRAW3 = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6",
         "KIAA0319", "RIDA", "SEM1", "SFMBT2", "SUMF2"]
CONTROL = DRAW1 + DRAW2 + DRAW3

# Colour scheme, consistent across every figure: panel dark, controls in the three draws
# distinguished but plainly one family.
C_PANEL = "#2b2118"
C_CTRL = "#d9a441"
C_DRAW2 = "#b8863b"
C_DRAW3 = "#8f6423"
GENOMEWIDE = 5e-8

FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)


def gene_class(g: str) -> str:
    if g in PANEL:
        return "panel"
    if g in DRAW1:
        return "control (draw 1)"
    if g in DRAW2:
        return "control (draw 2)"
    if g in DRAW3:
        return "control (draw 3)"
    raise KeyError(g)


def load_final_table() -> pd.DataFrame:
    """The per-arm accuracy + knockdown table the sweep drivers write."""
    df = pd.read_csv(FINAL_TABLE)
    df["classe"] = df["gene"].map(lambda g: "painel" if g in PANEL else "controle")
    df["draw"] = df["gene"].map(gene_class)
    return df


def load_gwas(tag: str = "uncorrected", n: str = "33") -> dict:
    """The gene-ranking JSON. `tag` is the GWAS arm, `n` the number of windows built."""
    for cand in (f"gene_ranking_{tag}_{n}.json", f"gene_ranking{n}_{tag}.json",
                 f"gene_ranking_{tag}.json"):
        p = RES / "gwas_ranking" / cand
        if p.exists():
            return json.loads(p.read_text())
    raise FileNotFoundError(f"no gene_ranking JSON for tag={tag} n={n}")


def load_top10(method: str = "biology_tss") -> tuple[pd.DataFrame, pd.DataFrame]:
    """(per-individual rows, per-gene summary) for the top-10 multi-gene knockdown.

    The per-gene summary is the regime report written alongside the replay, which already
    carries the AFR/EUR split and the uniformity and slope diagnostics; the per-individual
    rows are returned as well because the violin and the before/after figure need them.
    """
    d = pd.read_csv(TOP10_CSV)
    if "method" in d.columns:
        d = d[d["method"] == method]
    assert set(d["gene"]) == set(TOP10_GENES), sorted(set(d["gene"]))
    reg = pd.DataFrame(json.loads(TOP10_REGIME.read_text()))
    reg["class"] = reg["classe"].map({"painel": "panel", "controle": "control"})
    return d, reg


def load_top10_null() -> pd.DataFrame:
    """The null band: the same 100 bp scramble at a RANDOM window position instead of the
    promoter, one draw per (individual, gene), replayed through the same top-10 model.

    Panel and control genes were written under different method spellings
    (random_windowd0 and random_windowd1 -- different draws of the same procedure), so the
    method column is not filtered here; the gene set is asserted instead.
    """
    d = pd.read_csv(TOP10_NULL_CSV)
    assert set(d["gene"]) == set(TOP10_GENES), sorted(set(d["gene"]))
    return d


def top10_selection_drift() -> dict:
    """Whether the trained composition still matches a fresh selection on the same rule.

    The rule is: the ten highest balanced accuracies on the VALIDATION split, over every
    single-gene arm. Kept as a check rather than a caveat -- it returned a non-empty drift
    for the superseded test-selected arm, and should return an empty one here.
    """
    t = load_final_table().sort_values("bal_acc", ascending=False)
    now = t.head(10)["gene"].tolist()
    return {
        "trained": list(TOP10_GENES),
        "now": now,
        "enters": sorted(set(now) - set(TOP10_GENES)),
        "leaves": sorted(set(TOP10_GENES) - set(now)),
        "n_arms_now": len(t),
        "note": "load_final_table() carries TEST bal_acc; the selection used VAL, so a "
                "difference here is expected and is not drift in the selection rule",
    }


def window_interval(gene: str) -> tuple[str, int, int]:
    """(chrom, start_1based, end_1based) of the 524,288 bp window as built on disk."""
    meta = json.loads(
        (DATASET / "references" / "windows" / gene / "window_metadata.json").read_text())
    return meta["chromosome"], int(meta["start"]), int(meta["end"])


def crop_interval(gene: str, size: int = 32768) -> tuple[str, int, int]:
    """The central `size` bp actually read by the classifier, in reference coordinates."""
    chrom, s, e = window_interval(gene)
    full = e - s + 1
    size = min(size, full)
    cs = max(0, full // 2 - size // 2)
    return chrom, s + cs, s + cs + size - 1


def auc_and_p(scores: dict[str, float], genes: list[str], is_panel: list[bool],
              b_max: int = 2_000_000, seed: int = 13) -> tuple[float, float, str]:
    """AUC = P(panel outranks control), with its own panel/control permutation null.

    Exact when the number of label assignments is enumerable, Monte Carlo otherwise. Ties
    contribute 1/2, which is what makes this the Mann-Whitney statistic rather than a
    strict-inequality count.

    The p is ONE-SIDED (Pr[null >= obs]), which is what every per-method AUC p-value in
    the paper reports. fig_topk_admission.py reports its own head-to-head p two-sided, so
    the two are not directly comparable without halving or doubling one of them.
    """
    from itertools import combinations
    from math import comb

    v = np.array([scores[g] for g in genes], dtype=float)
    mask = np.array(is_panel, dtype=bool)
    n, n1 = len(v), int(mask.sum())
    ranks = pd.Series(v).rank().to_numpy()          # average ranks: ties at 1/2
    obs = (ranks[mask].sum() - n1 * (n1 + 1) / 2) / (n1 * (n - n1))

    total = comb(n, n1)
    if total <= b_max:
        idx = np.array(list(combinations(range(n), n1)), dtype=np.int64)
        mode = f"exact ({total:,} assignments)"
    else:
        rng = np.random.default_rng(seed)
        idx = np.concatenate([np.argsort(rng.random((min(200_000, b_max - lo), n)),
                                        axis=1)[:, :n1]
                              for lo in range(0, b_max, 200_000)], axis=0)
        mode = f"Monte Carlo (B = {len(idx):,} of {total:,})"
    null = (ranks[idx].sum(axis=1) - n1 * (n1 + 1) / 2) / (n1 * (n - n1))
    p = float((null >= obs - 1e-12).mean())
    return float(obs), p, mode


def savefig(fig, name: str) -> Path:
    out = FIGDIR / name
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")
    return out
