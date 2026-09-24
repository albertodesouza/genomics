#!/usr/bin/env python3
"""E21 -- normalised readout and ancestry-stratification of the knockdown response.

Two re-analyses of runs already on disk. No AlphaGenome calls, no GPU, no training.

(3) NORMALISED READOUT. The paper reports |Delta| (mean signed change in P(weak)) as the
    per-gene reliance statistic, and separately notes that |Delta| tracks |Delta_in|, the total
    absolute change arriving at the classifier's input. If the ranking is partly a ranking of how
    much perturbation each gene happens to deliver, then the defensible statistic is the response
    PER UNIT DELIVERED. We report the raw ratio, and -- because the ratio is unstable where
    |Delta_in| is near zero -- the residual from a log-log fit of |Delta| on |Delta_in| over the
    genes that clear the null.

(4a) ANCESTRY STRATIFICATION. Section IV-A of the paper argues the residual occupancy imprint is
    bit-identical before and after the edit and therefore "cancels exactly" in Delta. It does not:
    the classifier is nonlinear, so its sensitivity at a point depends on the whole input. If the
    imprint modulates the response, Delta should vary with ancestry beyond what the label explains.
    The pigmentation split is AFR-vs-EUR by construction, so the informative contrast is BETWEEN
    POPULATIONS WITHIN a superpopulation (YRI/ESN/LWK/MSL/GWD; FIN/CEU/GBR), where the label is
    constant. Kruskal-Wallis within each superpopulation, per gene.

Usage:
  python3 scripts/experiments/readout_normalisation.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home/breno/I2CA/genomics")
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
KNOCKOUT = KO_DIR / "pigmentation_test_split_knockout.csv"
EXPRDELTA = KO_DIR / "knockdown_expression_delta.csv"
NULL = KO_DIR / "pigmentation_test_split_null_scramble.csv"
OUT = REPO_ROOT / "results/genotype_based_predictor/readout_normalisation.json"

METHOD = "biology_tss"  # the annotation-based strategy the paper's Table II reports


def p_weak(strong_logit: pd.Series, weak_logit: pd.Series) -> np.ndarray:
    """Softmax over the two class logits, returning P(weak)."""
    s = np.vstack([strong_logit.to_numpy(float), weak_logit.to_numpy(float)])
    s = s - s.max(axis=0, keepdims=True)
    e = np.exp(s)
    return (e[1] / e.sum(axis=0))


def load_delta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["p_weak_base"] = p_weak(df["baseline_strong_logit"], df["baseline_weak_logit"])
    df["p_weak_pert"] = p_weak(df["perturbed_strong_logit"], df["perturbed_weak_logit"])
    df["delta"] = df["p_weak_pert"] - df["p_weak_base"]
    return df


def main() -> None:
    ko = load_delta(KNOCKOUT)
    ko = ko[ko["method"] == METHOD].copy()
    null = load_delta(NULL)

    # ---- per-gene Delta and null ------------------------------------------------
    per_gene = ko.groupby("gene")["delta"].agg(["mean", "std", "count"]).rename(
        columns={"mean": "delta", "std": "delta_sd", "count": "n"})
    per_gene["abs_delta"] = per_gene["delta"].abs()
    null_g = null.groupby("gene")["delta"].apply(lambda s: s.abs().mean()).rename("null_mean_abs")
    pooled_p99 = null["delta"].abs().quantile(0.99)
    per_gene = per_gene.join(null_g)
    per_gene["clears_null"] = per_gene["abs_delta"] > pooled_p99

    # ---- delivered perturbation: |Delta_in| summed over haplotypes ---------------
    ed = pd.read_csv(EXPRDELTA)
    ed = ed[ed["method"] == METHOD]
    din = ed.groupby(["sample_id", "gene"])["crop_abs_delta"].sum().rename("delta_in").reset_index()
    din_g = din.groupby("gene")["delta_in"].mean().rename("delta_in")
    per_gene = per_gene.join(din_g)

    # ---- (3) normalised readout --------------------------------------------------
    per_gene["ratio"] = per_gene["abs_delta"] / per_gene["delta_in"]
    est = per_gene[per_gene["clears_null"]].copy()
    lx = np.log10(est["delta_in"].to_numpy(float))
    ly = np.log10(est["abs_delta"].to_numpy(float))
    slope, intercept, r, p, se = stats.linregress(lx, ly)
    per_gene["logfit_resid"] = np.log10(per_gene["abs_delta"]) - (
        intercept + slope * np.log10(per_gene["delta_in"]))
    rho, rho_p = stats.spearmanr(per_gene["abs_delta"], per_gene["delta_in"])

    per_gene = per_gene.sort_values("abs_delta", ascending=False)

    print("=" * 78)
    print("(3) NORMALISED READOUT   method =", METHOD)
    print("=" * 78)
    print(f"Spearman |Delta| vs |Delta_in| over all 11 genes: rho = {rho:.3f} (p = {rho_p:.4f})")
    print(f"log-log fit on the {len(est)} genes clearing the pooled null p99 = {pooled_p99:.4f}:")
    print(f"    log10|Delta| = {intercept:.3f} + {slope:.3f} * log10|Delta_in|   (r = {r:.3f}, p = {p:.4f})")
    print()
    cols = ["abs_delta", "delta", "delta_in", "ratio", "logfit_resid", "clears_null"]
    with pd.option_context("display.width", 200, "display.float_format", lambda v: f"{v:11.6g}"):
        print(per_gene[cols].to_string())
    print()
    print("ranked by response per unit delivered (genes clearing the null only):")
    r_est = per_gene[per_gene["clears_null"]].sort_values("ratio", ascending=False)
    for g, row in r_est.iterrows():
        print(f"    {g:<9} ratio = {row['ratio']:.3e}   resid = {row['logfit_resid']:+.3f}")

    # ---- (4a) ancestry stratification within superpopulation ---------------------
    print()
    print("=" * 78)
    print("(4a) DELTA STRATIFIED BY POPULATION, WITHIN SUPERPOPULATION")
    print("=" * 78)
    strat = {}
    for sup, sub in ko.groupby("superpopulation"):
        pops = sorted(sub["population"].unique())
        print(f"\n{sup}  populations = {pops}")
        rows = []
        for gene, gsub in sub.groupby("gene"):
            groups = [g["delta"].to_numpy(float) for _, g in gsub.groupby("population")]
            groups = [g for g in groups if len(g) >= 3]
            if len(groups) < 2:
                continue
            H, pv = stats.kruskal(*groups)
            means = gsub.groupby("population")["delta"].mean()
            spread = float(means.max() - means.min())
            rows.append({"gene": gene, "H": H, "p": pv, "spread": spread,
                         "mean_delta": float(gsub["delta"].mean())})
        rdf = pd.DataFrame(rows).sort_values("p")
        # Holm-Bonferroni over the genes tested in this superpopulation
        m = len(rdf)
        rdf["p_holm"] = np.minimum.accumulate(
            (rdf["p"].to_numpy() * (m - np.arange(m)))[::-1])[::-1].clip(max=1.0)
        strat[sup] = rdf.to_dict("records")
        with pd.option_context("display.float_format", lambda v: f"{v:11.6g}"):
            print(rdf.to_string(index=False))

    payload = {
        "method": METHOD,
        "pooled_null_p99": float(pooled_p99),
        "spearman_absdelta_vs_deltain": {"rho": float(rho), "p": float(rho_p)},
        "loglog_fit": {"slope": float(slope), "intercept": float(intercept),
                       "r": float(r), "p": float(p), "stderr": float(se),
                       "n_genes": int(len(est))},
        "per_gene": json.loads(per_gene.reset_index().to_json(orient="records")),
        "stratification": strat,
    }
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
