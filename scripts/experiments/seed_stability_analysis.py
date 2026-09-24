#!/usr/bin/env python3
"""E22 analysis -- how repeatable is the knockdown readout across scramble draws?

Draw 0 is the headline run (seed 0); draws 1-4 come from seed_stability_knockout.py with the
target position held fixed, so the permutation seed is the only thing that differs. Draw 0 is
restricted to the same 60-individual subsample the later draws used, so all five draws are
estimates of the same quantity on the same individuals.

Reports the between-draw spread of the per-gene cohort mean (the statistic the paper quotes), the
stability of the induced gene ordering, and -- for contrast -- the much larger per-individual
spread, which is what makes averaging necessary.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home/breno/I2CA/genomics")
KO = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
HEAD = KO / "pigmentation_test_split_knockout.csv"
SEEDS = KO / "pigmentation_test_split_knockout_seeds.csv"
OUT = REPO_ROOT / "results/genotype_based_predictor/seed_stability.json"
METHOD = "biology_tss"


def delta(df: pd.DataFrame) -> np.ndarray:
    s = np.vstack([df["baseline_strong_logit"], df["baseline_weak_logit"]])
    s = s - s.max(0, keepdims=True); e = np.exp(s); b = e[1] / e.sum(0)
    s = np.vstack([df["perturbed_strong_logit"], df["perturbed_weak_logit"]])
    s = s - s.max(0, keepdims=True); e = np.exp(s); p = e[1] / e.sum(0)
    return p - b


def main() -> None:
    sd = pd.read_csv(SEEDS); sd["delta"] = delta(sd)
    genes = sorted(sd["gene"].unique())
    ids = sorted(sd["sample_id"].unique())

    hd = pd.read_csv(HEAD)
    hd = hd[(hd["method"] == METHOD) & (hd["gene"].isin(genes)) & (hd["sample_id"].isin(ids))].copy()
    hd["delta"] = delta(hd); hd["draw"] = 0
    both = pd.concat([hd[["sample_id", "gene", "draw", "delta"]],
                      sd[["sample_id", "gene", "draw", "delta"]]], ignore_index=True)

    print(f"{both['draw'].nunique()} draws x {len(genes)} genes x {len(ids)} individuals "
          f"= {len(both)} interventions")

    # ---- per-gene cohort mean, per draw ----
    piv = both.pivot_table(index="gene", columns="draw", values="delta", aggfunc="mean")
    piv = piv.reindex(piv.abs().max(axis=1).sort_values(ascending=False).index)
    summ = pd.DataFrame({
        "mean": piv.mean(axis=1), "sd": piv.std(axis=1, ddof=1),
        "min": piv.min(axis=1), "max": piv.max(axis=1)})
    summ["cv_pct"] = 100 * summ["sd"] / summ["mean"].abs()
    summ["range"] = summ["max"] - summ["min"]

    print("\n=== per-gene cohort-mean Delta, one column per scramble draw ===")
    with pd.option_context("display.width", 200, "display.float_format", lambda v: f"{v:9.5f}"):
        print(piv.to_string())
        print("\n=== between-draw spread of that cohort mean ===")
        print(summ.to_string())

    max_sd = float(summ["sd"].max()); max_sd_gene = summ["sd"].idxmax()
    max_cv = float(summ["cv_pct"].max()); max_cv_gene = summ["cv_pct"].idxmax()
    print(f"\nlargest between-draw SD  : {max_sd:.4f}  ({max_sd_gene})")
    print(f"largest between-draw CV  : {max_cv:.1f}%  ({max_cv_gene})")

    # ---- ordering stability: pairwise Spearman over |Delta| ----
    rhos = [stats.spearmanr(piv[a].abs(), piv[b].abs()).statistic
            for a, b in itertools.combinations(piv.columns, 2)]
    print(f"\npairwise Spearman of the |Delta| ordering across draws: "
          f"min {min(rhos):.3f}, median {float(np.median(rhos)):.3f}, max {max(rhos):.3f} "
          f"({len(rhos)} pairs)")
    ident = sum(1 for a, b in itertools.combinations(piv.columns, 2)
                if list(piv[a].abs().sort_values(ascending=False).index)
                == list(piv[b].abs().sort_values(ascending=False).index))
    print(f"draw pairs giving an identical gene ordering: {ident}/{len(rhos)}")

    # ---- ordering stability, properly: concordance over all draws at once ----
    # Pairwise Spearman shares draws between pairs, so its median/min summarise nothing
    # with a distribution. Kendall's W scores all m draws jointly against a known null.
    ranks = piv.abs().rank(axis=0, ascending=False)          # genes x draws
    n_obj, m_raters = ranks.shape
    S = ((ranks.sum(axis=1) - m_raters * (n_obj + 1) / 2) ** 2).sum()
    kendall_w = float(12 * S / (m_raters ** 2 * (n_obj ** 3 - n_obj)))
    chi2 = m_raters * (n_obj - 1) * kendall_w
    chi2_p = float(stats.chi2.sf(chi2, n_obj - 1))
    # The chi2 null needs n >> 7 to be accurate; at n=7, m=5 it is badly conservative
    # (1.7e-4 against 5e-6), so the paper quotes the permutation null instead.
    rng = np.random.default_rng(0)
    n_perm, hits = 200_000, 0
    for _ in range(n_perm):
        P = np.array([rng.permutation(n_obj) + 1 for _ in range(m_raters)]).T
        s_ = ((P.sum(axis=1) - m_raters * (n_obj + 1) / 2) ** 2).sum()
        if 12 * s_ / (m_raters ** 2 * (n_obj ** 3 - n_obj)) >= kendall_w:
            hits += 1
    w_p = (hits + 1) / (n_perm + 1)
    print(f"\nKendall's W over all {m_raters} draws x {n_obj} genes: W = {kendall_w:.3f}")
    print(f"  permutation null ({n_perm} draws): p = {w_p:.2g}")
    print(f"  chi2 approximation             : p = {chi2_p:.2g} (conservative at this n)")
    # W is a monotone rescaling of the mean pairwise Spearman: W = ((m-1)*rho_bar + 1)/m.
    # It carries the same information, but as one statistic with a null rather than a
    # summary over dependent pairs.
    print(f"  identity check, ((m-1)*mean_rho+1)/m = "
          f"{((m_raters - 1) * float(np.mean(rhos)) + 1) / m_raters:.6f}")

    # Ordering churn is dominated by genes the instrument does not resolve. Restrict to the
    # genes whose mean magnitude clears their own between-draw spread by >= RESOLVED_RATIO.
    # The natural gap in this panel is between MC1R (9.8) and DDB1 (7.6), which is also the
    # boundary the paper draws between "resolved" and "adequate".
    RESOLVED_RATIO = 8.0
    ratio = summ["mean"].abs() / summ["sd"]
    resolved = list(ratio[ratio >= RESOLVED_RATIO].index)
    sub = piv.loc[resolved].abs()
    orders = {tuple(sub[c].sort_values(ascending=False).index) for c in sub.columns}
    print(f"resolved genes (|mean|/SD >= {RESOLVED_RATIO:g}): {resolved}")
    print(f"  distinct orderings across the {m_raters} draws: {len(orders)}")
    if len(orders) == 1:
        print(f"  identical in every draw: {' > '.join(next(iter(orders)))}")

    # ---- per-individual spread, for contrast ----
    per_ind = both.groupby(["gene", "sample_id"])["delta"].std(ddof=1)
    print("\n=== per-individual SD across draws (mean over individuals) ===")
    pi = per_ind.groupby("gene").mean().reindex(piv.index)
    with pd.option_context("display.float_format", lambda v: f"{v:9.5f}"):
        print(pd.DataFrame({"per_individual_sd": pi, "cohort_mean_sd": summ["sd"],
                            "ratio": pi / summ["sd"]}).to_string())

    OUT.write_text(json.dumps({
        "n_draws": int(both["draw"].nunique()), "n_genes": len(genes), "n_individuals": len(ids),
        "max_between_draw_sd": max_sd, "max_between_draw_sd_gene": max_sd_gene,
        "max_between_draw_cv_pct": max_cv, "max_between_draw_cv_gene": max_cv_gene,
        "ordering_spearman": {"min": float(min(rhos)), "median": float(np.median(rhos)),
                              "max": float(max(rhos)), "n_pairs": len(rhos),
                              "identical_orderings": int(ident)},
        "ordering_concordance": {"kendall_w": kendall_w, "p_permutation": w_p,
                                 "n_permutations": n_perm,
                                 "chi2": float(chi2), "df": int(n_obj - 1), "p_chi2": chi2_p,
                                 "mean_pairwise_spearman": float(np.mean(rhos)),
                                 "resolved_ratio_threshold": RESOLVED_RATIO,
                                 "resolved_genes": resolved,
                                 "n_distinct_orderings_resolved": len(orders),
                                 "resolved_ordering": list(next(iter(orders))) if len(orders) == 1 else None},
        "per_gene": json.loads(summ.reset_index().to_json(orient="records")),
        "per_draw_means": json.loads(piv.reset_index().to_json(orient="records")),
    }, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
