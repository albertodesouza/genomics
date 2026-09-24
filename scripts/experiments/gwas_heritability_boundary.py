#!/usr/bin/env python3
"""Whether a GRM-based mixed model degenerates on this cohort the way PC-correction does.

`gwas_gene_ranking.py`'s `union_pc10` arm already shows the fixed-effect remedy fails:
every one of 58,270 Firth fits with ten genotype PCs as covariates returns ERRCODE
UNFINISHED (see `results/.../gwas_ranking/pc_correction_summary.json`, written by this
script's sibling check). The paper's argument extends that failure to mixed models on
mechanistic grounds -- the GRM's leading eigenvectors are the same PCs, so a random
polygenic effect is collinear with the label the same way PC1 is -- but a mechanism is
not a measurement. This script is the measurement: it fits

    y = mu * 1 + u + e,   u ~ N(0, sigma_g^2 G),   e ~ N(0, sigma_e^2 I)

by REML, G the standardised-genotype GRM (GCTA definition) over the same 58,270 MAF-
filtered variants `union_pc10` is run on, y the 0/1 pigmentation label treated as
quantitative -- the same approximation a first-pass GREML heritability check always
makes, not a liability-threshold model. The question is not the point estimate of h^2,
it is whether the REML profile likelihood has an interior maximum or is monotone
increasing to the h^2 = 1 boundary of the parameter space; the latter is the mixed-model
analogue of "0 of 58,270 Firth fits converged" and is what the paper's argument predicts.

Reads `/tmp/gwas_gene_ranking/{geno.npy,union_sites.tsv,union_maf.pvar,pheno.tsv}`,
written by `gwas_gene_ranking.py --stage build` on 2026-09-11 for the current panel
(9 genes) + both control draws (11 + 11) union. Writes
`results/genotype_based_predictor/gwas_ranking/heritability_boundary.json`.
Zero API cost, no plink2 call, no new genotype build.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

WORK = Path("/tmp/gwas_gene_ranking")
OUT = (Path("/home/breno/I2CA/genomics/results/genotype_based_predictor")
       / "gwas_ranking" / "heritability_boundary.json")


def load_maf_columns() -> np.ndarray:
    """Column indices into geno.npy for the 58,270 sites union_maf.pvar kept."""
    col_of = {}
    with open(WORK / "union_sites.tsv") as f:
        next(f)
        for line in f:
            chrom, pos, ref, alt, col = line.rstrip("\n").split("\t")
            col_of[(chrom[3:], int(pos), ref, alt)] = int(col)
    cols = []
    with open(WORK / "union_maf.pvar") as f:
        next(f)
        for line in f:
            chrom, pos, _id, ref, alt = line.rstrip("\n").split("\t")
            cols.append(col_of[(chrom, int(pos), ref, alt)])
    return np.asarray(cols, dtype=np.int64)


def load_label() -> np.ndarray:
    ids, pigm = [], []
    with open(WORK / "pheno.tsv") as f:
        next(f)
        for line in f:
            iid, p, _r = line.rstrip("\n").split("\t")
            ids.append(iid)
            pigm.append(int(p))
    y = np.asarray(pigm, dtype=float) - 1.0  # plink 1/2 -> 0/1
    return y


def grm(geno_cols: np.ndarray) -> np.ndarray:
    """GCTA-style GRM: G = Z Z' / M, Z the per-site-standardised genotype matrix."""
    x = geno_cols.astype(np.float64)
    p = x.mean(axis=0) / 2.0
    denom = np.sqrt(2.0 * p * (1.0 - p))
    keep = denom > 1e-9
    z = (x[:, keep] - 2.0 * p[keep]) / denom[keep]
    m = int(keep.sum())
    return (z @ z.T) / m, m


def reml_profile(y: np.ndarray, evals: np.ndarray, u_ty: np.ndarray, u_t1: np.ndarray,
                  h2_grid: np.ndarray) -> np.ndarray:
    """REML log-likelihood (up to an additive constant) for an intercept-only mixed
    model, profiled over mu and the total variance at each h2, via the eigendecomposition
    trick (Kang et al. 2008): with G = U diag(evals) U', d_i(h2) = h2*evals_i + (1-h2),
    the GLS estimate of mu and the residual sum of squares are weighted sums over the
    rotated data, and

      logL_REML(h2) = -0.5 * [ (n-1) log(2 pi) + (n-1) + (n-1) log(RSS(h2)/(n-1))
                                + sum_i log(d_i(h2)) + log(sum_i u_t1_i^2 / d_i(h2)) ]

    n - 1 rather than n because X = [1] has rank 1; the term is dropped consistently
    across the grid so only relative values matter.
    """
    n = len(y)
    out = np.empty_like(h2_grid)
    for j, h2 in enumerate(h2_grid):
        d = h2 * evals + (1.0 - h2)
        d = np.maximum(d, 1e-10)
        w1 = (u_t1 * u_t1 / d).sum()
        wy1 = (u_t1 * u_ty / d).sum()
        mu = wy1 / w1
        resid = u_ty - u_t1 * mu
        rss = (resid * resid / d).sum()
        sigma2 = rss / (n - 1)
        out[j] = -0.5 * (np.log(d).sum() + (n - 1) * np.log(sigma2) + np.log(w1))
    return out


def main() -> int:
    cols = load_maf_columns()
    y = load_label()
    geno = np.load(WORK / "geno.npy", mmap_mode="r")
    geno_cols = np.asarray(geno[:, cols])
    g, m_eff = grm(geno_cols)
    n = g.shape[0]

    evals, evecs = np.linalg.eigh(g)
    evals = np.clip(evals, 1e-8, None)
    u_ty = evecs.T @ (y - y.mean())
    u_t1 = evecs.T @ np.ones(n)

    h2_grid = np.linspace(0.0, 1.0, 4001)
    logl = reml_profile(y, evals, u_ty, u_t1, h2_grid)
    i_max = int(np.argmax(logl))
    h2_hat = float(h2_grid[i_max])

    result = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "n_individuals": n,
        "n_variants_grm": m_eff,
        "model": "y = mu*1 + u + e, u~N(0,sigma_g^2 G), e~N(0,sigma_e^2 I); "
                 "y is the 0/1 pigmentation label treated as quantitative "
                 "(linear GREML, not a liability-threshold model)",
        "h2_grid_step": float(h2_grid[1] - h2_grid[0]),
        "h2_hat_argmax": h2_hat,
        "logL_at_h2": {
            "0.00": float(logl[0]),
            "0.50": float(logl[np.searchsorted(h2_grid, 0.50)]),
            "0.90": float(logl[np.searchsorted(h2_grid, 0.90)]),
            "0.99": float(logl[np.searchsorted(h2_grid, 0.99)]),
            "1.00": float(logl[-1]),
        },
        "logl_monotone_increasing": bool(np.all(np.diff(logl) >= -1e-9)),
        "logl_curve": {"h2": h2_grid.tolist(), "logL": logl.tolist()},
        "grm_top_eigenvalues": evals[-5:][::-1].tolist(),
        "grm_eigenvalue_sum": float(evals.sum()),
        "note": "logl_monotone_increasing=True means REML places the maximum at the "
                "h2=1 boundary of the parameter space rather than at an interior point "
                "-- the mixed-model analogue of the 0-of-58270 Firth non-convergence "
                "under fixed-effect PC correction, not a calibrated heritability.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({k: v for k, v in result.items() if k != "logl_curve"},
                               indent=2))
    (OUT.parent / "heritability_boundary_curve.json").write_text(
        json.dumps({"h2": result["logl_curve"]["h2"], "logL": result["logl_curve"]["logL"]}))
    print(json.dumps({k: v for k, v in result.items()
                       if k not in ("logl_curve",)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
