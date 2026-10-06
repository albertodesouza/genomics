#!/usr/bin/env python3
"""Gene scores for the 22 windows, and the head-to-head test against the CNN ranking.

Consumes the work directory of gwas_gene_ranking.py. Two halves:

SCORING. For each 524,288 bp window (and, separately, for the 32,768 bp crop the CNN
actually reads), four gene scores over the window's MAF-filtered variants -- the raw
minimum, the Sidak-corrected minimum with Li & Ji effective test count, the MAGMA-style
mean chi-squared, and GATES -- all as S_g = -log10 p_gene in log space. The variant
correlation matrix R is rank-deficient (1072 individuals, thousands of variants), so
its non-zero eigenvalues are taken from the 1072x1072 Gram matrix of the standardised
genotypes, which is exact and cheap where an m x m decomposition would not be.

COMPARISON. Both estimators score the same genes -- every gene that has both a window
score and a trained classifier, which is nine pigmentation genes and the controls of both
random draws -- so the ranking each induces can be tested against the other rather than
each against its own null. The gene list is derived from the CNN table rather than fixed
here, so adding control arms widens the comparison instead of silently leaving them out. The statistic is the AUC of panel against
control, i.e. the probability that a randomly drawn panel gene outranks a randomly
drawn control, which is the Mann-Whitney U divided by n1*n2 and is on the same scale
for both methods. The head-to-head null is exact: the panel/control assignment is
exchangeable under it, and there are only C(20,9) = 167,960 assignments, so every one
is enumerated and both AUCs are recomputed from the same fixed score vectors on each.
Enumerating rather than sampling matters here because the two scores are observed on
the same genes and are therefore correlated; permuting the labels preserves that
correlation, where two independent one-sample tests would not. Above C(n, n1) = MAX_EXACT
the enumeration falls back to sampling and the mode is recorded in the output. DeLong's asymptotic test
for two correlated AUCs is reported alongside as a cross-check.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.special import gammaln, gammaincc

REPO_ROOT = Path("/home/breno/I2CA/genomics")
RESULTS = REPO_ROOT / "results" / "genotype_based_predictor" / "gwas_ranking"
CNN_TABLE = REPO_ROOT / "results" / "genotype_based_predictor" / "poolmax_final_table.csv"
LN10 = math.log(10.0)
GW_L = -math.log10(5e-8)
NOT_PIGMENTATION = {"EDAR", "TCHH"}


def _log(m: str) -> None:
    print(m, flush=True)


# ---------------------------------------------------------------- gene scores


def eig_nonzero(Z: np.ndarray) -> np.ndarray:
    """Non-zero eigenvalues of R = Z'Z/(n-1) via the n x n Gram matrix.

    Z is n x m standardised genotypes. R is m x m of rank <= n-1, and the non-zero
    spectrum of Z'Z equals that of ZZ', so the small side is decomposed. Returns the
    eigenvalues of R, zeros omitted; sum(lambda) = trace(R) = m up to the rank deficit.
    """
    n = Z.shape[0]
    G = (Z @ Z.T) / (n - 1)
    w = np.linalg.eigvalsh(G.astype(np.float64))
    return np.clip(w, 0.0, None)


def m_eff_li_ji(lam: np.ndarray, m: int) -> float:
    """Li & Ji (2005): each eigenvalue contributes 1 if >= 1, plus its fractional part."""
    v = float(np.sum((lam >= 1.0).astype(float) + (lam - np.floor(lam))))
    return float(min(max(v, 1.0), m))


def m_eff_nyholt(sum_lam_sq: float, m: int) -> float:
    """Cheverud/Nyholt: 1 + (m-1)(1 - Var(lambda)/m), Var from sum of squares."""
    var_lam = sum_lam_sq / m - 1.0
    v = 1.0 + (m - 1.0) * (1.0 - var_lam / m)
    return float(min(max(v, 1.0), m))


def s_sidak(L_min: float, m_eff: float) -> float:
    """-log10 of 1 - (1 - p_min)^m_eff, exact for large p_min, asymptotic for small."""
    p_min = 10.0 ** (-L_min)
    if m_eff * p_min < 1e-8:
        return L_min - math.log10(m_eff)
    p = -math.expm1(m_eff * math.log1p(-p_min))
    return -math.log10(min(max(p, 1e-300), 1.0))


def log_chi2_sf(x: float, d: float) -> float:
    """Natural log of P(chi2_d > x), valid where the survival function underflows.

    scipy's logsf is the log of an already-underflowed tail, so it returns -inf once
    the aggregated statistic runs past a double. Beyond that point the upper incomplete
    gamma is taken from its asymptotic expansion, which needs y >> a and is exactly the
    regime that underflows. A continental contrast puts these windows there.
    """
    a, y = d / 2.0, x / 2.0
    v = gammaincc(a, y)
    if v > 1e-300:
        return math.log(v)
    if y <= 4.0 * a:                              # series unreliable; Wilson-Hilferty
        z = ((x / d) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * d))) / math.sqrt(2.0 / (9.0 * d))
        return float(stats.norm.logsf(z))
    s, term = 1.0, 1.0
    for k in range(1, 200):
        term *= (a - k) / y
        if abs(term) < 1e-18 * abs(s):
            break
        s += term
    return (a - 1.0) * math.log(y) - y + math.log(max(s, 1e-300)) - gammaln(a)


def s_meanchi2(chi2: np.ndarray, sum_lam_sq: float, m: int) -> float:
    """Satterthwaite tail of the window's chi-squared sum, a quadratic form in N(0, R).

    Q = z'z with z ~ N(0, R) is distributed as sum_i lambda_i * chi2(1), whose first two
    moments are m and 2*trace(R^2); matching them to a scaled chi2 gives the tail below.
    At this label the result runs to thousands of -log10 units, where it is a ranking
    device and not a calibrated probability -- the same caveat S_min carries at 88.
    """
    Q = float(np.sum(chi2))
    c = sum_lam_sq / m
    d = m * m / sum_lam_sq
    if not (c > 0 and d > 0):
        return float("nan")
    return float(-log_chi2_sf(Q / c, d) / LN10)


def s_gates(L_sorted: np.ndarray, Z_sorted: np.ndarray, m_eff_tot: float) -> tuple[float, int]:
    """GATES: -log10 min_k m_eff * p_(k) / m_eff(k), nested counts by Nyholt.

    m_eff(k) needs only the sum of squared off-diagonal correlations among the k
    smallest p-values, which is a cumulative sum, so every k is covered by one gemm.
    """
    n, m = Z_sorted.shape
    R = (Z_sorted.T @ Z_sorted) / (n - 1)
    np.square(R, out=R)
    rowsum = np.tril(R, -1).sum(axis=1)
    S = np.cumsum(rowsum, dtype=np.float64)          # S[k-1] = sum_{j<i<=k} r_ij^2
    k = np.arange(1, m + 1, dtype=np.float64)
    var_lam = 2.0 * S / k
    m_eff_k = np.where(k >= 2, 1.0 + (k - 1.0) * (1.0 - var_lam / k), 1.0)
    m_eff_k = np.clip(m_eff_k, 1.0, k)
    # S_g = max_k ( L_(k) + log10 m_eff(k) - log10 m_eff_tot )
    obj = L_sorted + np.log10(m_eff_k) - math.log10(m_eff_tot)
    j = int(np.argmax(obj))
    return float(obj[j]), j + 1


def score_window(L: np.ndarray, chi2: np.ndarray, Z: np.ndarray) -> dict:
    m = int(L.size)
    if m == 0:
        return {"m": 0}
    lam = eig_nonzero(Z)
    sum_lam_sq = float(np.sum(lam ** 2))
    me_lj = m_eff_li_ji(lam, m)
    me_ny = m_eff_nyholt(sum_lam_sq, m)
    o = np.argsort(-L)                                # ascending p == descending -log10 p
    s_g, k_star = s_gates(L[o], np.ascontiguousarray(Z[:, o]), me_ny)
    return {
        "m": m,
        "m_eff_li_ji": me_lj,
        "m_eff_nyholt": me_ny,
        "n_genomewide": int(np.sum(L > GW_L)),
        "frac_genomewide": float(np.mean(L > GW_L)),
        "S_min": float(np.max(L)),
        "S_sidak": s_sidak(float(np.max(L)), me_lj),
        "S_sidak_nyholt": s_sidak(float(np.max(L)), me_ny),
        "S_meanchi2": s_meanchi2(chi2, sum_lam_sq, m),
        "S_gates": s_g,
        "gates_k": k_star,
    }


# ---------------------------------------------------------------- comparison


def midranks(x: np.ndarray) -> np.ndarray:
    return stats.rankdata(x, method="average")


def auc_from_ranks(r: np.ndarray, idx: tuple, n1: int, n0: int) -> float:
    return (float(r[list(idx)].sum()) - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def delong_cov(scores: dict, panel: np.ndarray, ctrl: np.ndarray) -> tuple:
    """DeLong structural components for two correlated AUCs on shared case/control sets."""
    keys = list(scores)
    n1, n0 = panel.sum(), ctrl.sum()
    V10, V01, aucs = {}, {}, {}
    for k in keys:
        x = scores[k][panel]
        y = scores[k][ctrl]
        psi = lambda a, b: (a > b) * 1.0 + (a == b) * 0.5
        V10[k] = np.array([psi(xi, y).mean() for xi in x])
        V01[k] = np.array([psi(x, yj).mean() for yj in y])
        aucs[k] = float(V10[k].mean())
    S10 = np.cov(np.vstack([V10[k] for k in keys]), ddof=1)
    S01 = np.cov(np.vstack([V01[k] for k in keys]), ddof=1)
    S = S10 / n1 + S01 / n0
    return aucs, np.atleast_2d(S)


MAX_EXACT = 2_000_000


def _label_sums(ranks: dict, n: int, n1: int, rng_seed: int = 13) -> tuple[dict, str, int]:
    """Rank sums of the panel group under every (or many) panel/control assignments.

    The assignment is the exchangeable object, so permuting it and recomputing both
    AUCs from the same fixed score vectors preserves the correlation between the two
    methods -- which is the whole reason a paired test is available here. Enumerated
    exactly while C(n, n1) is tractable; above that, sampled, with the count reported
    so the p-value's resolution is on record.
    """
    total = math.comb(n, n1)
    if total <= MAX_EXACT:
        subsets = np.array(list(combinations(range(n), n1)), dtype=np.int64)
        mode, draws = "exact", total
    else:
        rng = np.random.default_rng(rng_seed)
        draws = MAX_EXACT
        chunks = []
        for lo in range(0, draws, 200_000):
            b = min(200_000, draws - lo)
            chunks.append(np.argsort(rng.random((b, n)), axis=1)[:, :n1])
        subsets = np.concatenate(chunks, axis=0)
        mode = "montecarlo"
    return {k: np.asarray(v)[subsets].sum(axis=1) for k, v in ranks.items()}, mode, draws


def head_to_head(names: list[str], score_map: dict, is_panel: np.ndarray) -> dict:
    n = len(names)
    n1 = int(is_panel.sum())
    n0 = n - n1
    ranks = {k: midranks(v) for k, v in score_map.items()}
    obs_idx = tuple(np.flatnonzero(is_panel))
    obs = {k: auc_from_ranks(ranks[k], obs_idx, n1, n0) for k in score_map}
    sums, mode, draws = _label_sums(ranks, n, n1)
    _log(f"  {mode}: {draws} panel/control assignments over {n} genes")
    null = {k: (v - n1 * (n1 + 1) / 2.0) / (n1 * n0) for k, v in sums.items()}
    out = {"n_panel": n1, "n_control": n0, "n_assignments": int(draws),
           "null_mode": mode, "per_method": {}}
    for k in score_map:
        out["per_method"][k] = {
            "auc": obs[k],
            "mean_rank_panel": float(np.mean([n - ranks[k][i] + 1 for i in obs_idx])),
            "mean_rank_control": float(np.mean([n - ranks[k][i] + 1
                                                for i in range(n) if not is_panel[i]])),
            "p_exact_one_sided": float(np.mean(null[k] >= obs[k] - 1e-12)),
        }
    return out, obs, null, ranks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path,
                    default=Path(os.environ.get("GWAS_RANK_WORK", "/tmp/gwas_gene_ranking")))
    ap.add_argument("--arm", default="uncorrected",
                    choices=("uncorrected", "founders", "randomlabel", "pc10"))
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    os.chdir(REPO_ROOT)
    RESULTS.mkdir(parents=True, exist_ok=True)
    work = args.work
    out_path = args.out or (RESULTS / f"gene_ranking_{args.arm}.json")

    info = json.loads((work / "build_info.json").read_text())
    windows = info["windows"]
    geno = np.load(work / "geno.npy")                          # individuals x union sites
    col = {}
    with open(work / "union_sites.tsv") as f:
        next(f)
        for line in f:
            c, p, r, a, j = line.rstrip("\n").split("\t")
            col[(c, int(p), r, a)] = int(j)

    name = "RANDOM" if args.arm == "randomlabel" else "PIGM"
    sf = work / f"union_{args.arm}.{name}.glm.logistic.hybrid"
    hdr = None
    rec = []
    with open(sf) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = {k: i for i, k in enumerate(p)}
                continue
            if p[hdr["TEST"]] != "ADD" or p[hdr["ERRCODE"]] != ".":
                continue
            try:
                L = float(p[hdr["NEG_LOG10_P"]])
                z = float(p[hdr["Z_STAT"]] if "Z_STAT" in hdr else p[hdr["T_STAT"]])
            except (ValueError, KeyError):
                continue
            if not (math.isfinite(L) and math.isfinite(z)):
                continue
            rec.append((f"chr{p[hdr['#CHROM']]}", int(p[hdr["POS"]]),
                        p[hdr["REF"]] if "REF" in hdr else "", p[hdr["ALT"]] if "ALT" in hdr else "",
                        L, z, p[hdr["ID"]]))
    _log(f"{args.arm}: {len(rec)} converged ADD tests")
    by_chrom: dict = {}
    for chrom, pos, ref, alt, L, z, vid in rec:
        parts = vid.split(":")
        key = (chrom, pos, parts[2], parts[3]) if len(parts) == 4 else None
        j = col.get(key) if key else None
        if j is None:
            continue
        by_chrom.setdefault(chrom, []).append((pos, L, z, j))
    for c in by_chrom:
        by_chrom[c].sort()
    matched = sum(len(v) for v in by_chrom.values())
    _log(f"matched to union columns: {matched}")

    lam_gc = float(np.median(np.array([z for _, _, z, _ in
                                       sum(by_chrom.values(), [])]) ** 2) / stats.chi2.ppf(0.5, 1))
    _log(f"lambda_GC over the {len(windows)} windows: {lam_gc:.2f}")

    scores = {"window": {}, "crop": {}}
    for g, w in windows.items():
        arr = by_chrom.get(w["chrom"], [])
        for tag, lo, hi in (("window", w["start"], w["end"]),
                            ("crop", w["crop_start"], w["crop_end"])):
            sel = [(L, z, j) for pos, L, z, j in arr if lo <= pos <= hi]
            if not sel:
                scores[tag][g] = {"m": 0}
                continue
            L = np.array([s[0] for s in sel])          # -log10 p, straight from plink2
            chi2 = np.array([s[1] for s in sel]) ** 2  # Z^2, exact where chi2.isf underflows
            idx = np.array([s[2] for s in sel])
            G = geno[:, idx].astype(np.float32)
            mu = G.mean(axis=0)
            sd = G.std(axis=0)
            keep = sd > 0
            if not keep.all():
                G, L, chi2, keep_n = G[:, keep], L[keep], chi2[keep], int(keep.sum())
                mu, sd = mu[keep], sd[keep]
            Z = (G - mu) / sd
            scores[tag][g] = score_window(L, chi2, Z)
            _log(f"  {g:9s} {tag:7s} m={scores[tag][g]['m']:6d} "
                 f"m_eff(LiJi)={scores[tag][g]['m_eff_li_ji']:8.1f} "
                 f"S_min={scores[tag][g]['S_min']:7.1f} S_sidak={scores[tag][g]['S_sidak']:7.1f} "
                 f"S_chi2={scores[tag][g]['S_meanchi2']:8.1f} S_gates={scores[tag][g]['S_gates']:7.1f}")

    # Both splits are read. The head-to-head ranks on test, which is the split the
    # classifier is reported on; the validation ranking is carried alongside because it is
    # the ranking that selected the probe's ten genes, and the two orderings differ
    # (Spearman +0.895), so a downstream reader must be able to tell which is which.
    rows = list(csv.DictReader(open(CNN_TABLE)))
    if "bal_acc_val" not in rows[0]:
        raise SystemExit("ABORT: poolmax_final_table.csv has no bal_acc_val column; "
                         "rerun scripts/experiments/poolmax_final_table.py")
    cnn = {r["gene"]: (float(r["bal_acc"]), r["classe"]) for r in rows}
    cnn_val = {r["gene"]: float(r["bal_acc_val"]) for r in rows}
    genes = [g for g in windows if g in cnn and g not in NOT_PIGMENTATION]
    genes.sort()
    is_panel = np.array([cnn[g][1] == "painel" for g in genes])
    _log(f"\nhead-to-head over {len(genes)} genes: "
         f"{int(is_panel.sum())} panel, {int((~is_panel).sum())} control")

    result = {"arm": args.arm, "lambda_gc_windows": lam_gc, "gene_scores": scores,
              "cnn": {g: cnn[g][0] for g in genes},
              "cnn_val": {g: cnn_val[g] for g in genes}, "genes": genes,
              "is_panel": [bool(b) for b in is_panel], "comparisons": {}}

    for tag in ("window", "crop"):
        smap = {"CNN_bal_acc": np.array([cnn[g][0] for g in genes]),
                "CNN_bal_acc_val": np.array([cnn_val[g] for g in genes])}
        for k in ("S_min", "S_sidak", "S_meanchi2", "S_gates"):
            smap[f"GWAS_{k}"] = np.array([scores[tag][g].get(k, np.nan) for g in genes])
        h2h, obs, null, ranks = head_to_head(genes, smap, is_panel)
        pairs = {}
        for k in smap:
            # Only GWAS estimators are comparators here; the second CNN column is a
            # cross-split check and belongs in per_method, not in cnn_vs_gwas.
            if k.startswith("CNN_"):
                continue
            d = obs["CNN_bal_acc"] - obs[k]
            nd = null["CNN_bal_acc"] - null[k]
            pairs[k] = {
                "delta_auc_cnn_minus_gwas": float(d),
                "p_exact_two_sided": float(np.mean(np.abs(nd) >= abs(d) - 1e-12)),
                "p_exact_one_sided_cnn_higher": float(np.mean(nd >= d - 1e-12)),
                "spearman_rho_vs_cnn": float(stats.spearmanr(smap["CNN_bal_acc"], smap[k]).statistic),
            }
        aucs, S = delong_cov(smap, is_panel, ~is_panel)
        keys = list(smap)
        ci = keys.index("CNN_bal_acc")
        for k in pairs:
            j = keys.index(k)
            var = S[ci, ci] + S[j, j] - 2 * S[ci, j]
            z = (aucs["CNN_bal_acc"] - aucs[k]) / math.sqrt(var) if var > 0 else float("nan")
            pairs[k]["delong_z"] = float(z)
            pairs[k]["delong_p_two_sided"] = float(2 * stats.norm.sf(abs(z))) if math.isfinite(z) else float("nan")
        result["comparisons"][tag] = {"per_method": h2h["per_method"],
                                      "n_panel": h2h["n_panel"], "n_control": h2h["n_control"],
                                      "n_assignments": h2h["n_assignments"],
                                      "null_mode": h2h["null_mode"],
                                      "cnn_vs_gwas": pairs}

    out_path.write_text(json.dumps(result, indent=1))
    _log(f"\nwrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
