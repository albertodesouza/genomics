#!/usr/bin/env python3
"""Deterministic, magnitude-matched gene ablation -- the knockdown without the scramble.

WHY THIS EXISTS
---------------
The 100 bp promoter scramble has three documented pathologies, all of which the
paper reports and none of which it can remove:

  1. Permutation noise is common-mode across individuals, so cohort averaging
     does not reduce it. TYRP1 clears the random-window null at 0.251 and then
     reverses sign between draws (between-draw SD 0.180, |mean|/SD = 0.4).
  2. Delivered magnitude spans three orders across the panel (TYR 54,297 units,
     TCHH 53), and the classifier responds to *absolute* change, so the raw
     ranking is partly a ranking of delivered magnitude rather than of reliance.
  3. It perturbs the promoter, which for five of eleven genes falls entirely
     outside the 32 kbp crop the classifier reads.

This script perturbs the gene's channel block *directly in the classifier's input
tensor*, which addresses all three at once:

  - deterministic: no seed, so pathology 1 cannot occur and the effective sample
    size that repeatability analysis takes away is recovered;
  - magnitude is a free parameter, so pathology 2 is removed by construction
    rather than corrected post hoc by a six-point log-log fit;
  - it reaches the whole block, so pathology 3 does not apply.

The cost is that a scaled-down block is out of distribution for the classifier,
which is declarable and is why this is reported *alongside* the scramble rather
than instead of it. The comparison between the two is itself the result: the
scramble measures reliance mediated by the frozen model's response, the ablation
measures reliance on the gene's channel block, and where they disagree the
disagreement localises whether the bottleneck is the frozen model or the CNN.

THE DESIGN, AND WHY A SWEEP RATHER THAN ONE ABLATION
----------------------------------------------------
Zeroing a gene's block entirely delivers a perturbation equal to that gene's own
total signal, so a full ablation inherits exactly the magnitude confound it was
meant to remove -- a highly expressed gene is perturbed harder.

Instead we scale gene g's block by lambda in [0,1] and record both

    delivered(g, lambda) = (1 - lambda) * ||block_g||_1        (an identity)
    Delta(g, lambda)     = mean over individuals of
                           P(weak | scaled) - P(weak | original)

over a grid of lambda. That traces each gene's response as a function of
delivered magnitude. A matched comparison at any budget B is then a *read-off*
from these curves at delivered = B -- no fit, no extrapolation, and no gene
measured against a line it helped determine. Genes whose full block is smaller
than B simply cannot reach that budget, and the script says so rather than
extrapolating; that inability is itself the finding for the dead-zone genes.

PRE-REGISTERED READING
----------------------
The claim under test is that SLC24A5 is read more efficiently per unit delivered
perturbation than anything else in the panel -- currently supported by a
leave-one-out log-log fit over six genes that is significant (p = 0.0022) but
anchored at its low end by a single gene, MC1R, without which the relationship
vanishes.

  - If the magnitude-matched ablation puts SLC24A5 on top at a shared budget,
    that claim survives without depending on any fit.
  - If it does not, the claim was an artefact of the normalisation and should be
    withdrawn.

Either outcome is worth having, and this measurement can produce it today: no
AlphaGenome calls, no retraining, forward passes only.

Usage:
  python3 scripts/experiments/deterministic_gene_ablation.py
  ... --config <path>     # default: the probed DITA tracks-only classifier
  ... --limit 20          # smoke test on a few individuals
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _p in (REPO_ROOT / "src", REPO_ROOT / "notebooks", REPO_ROOT / "scripts" / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
os.chdir(REPO_ROOT)

import numpy as np
import torch

# The probed classifier: tracks-only, no genotype channels, shared axis.
DEFAULT_CONFIG = "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml"
OUT_PATH = REPO_ROOT / "results" / "genotype_based_predictor" / "deterministic_gene_ablation.json"

# Default lambda grid: dense near 1 (small perturbations, where the response is
# most likely to be locally linear) and reaching 0 (full ablation).
LAMBDA_GRID = [1.0, 0.9999, 0.999, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1, 0.0]

# Published scramble readout on this same classifier, for the comparison that is
# the point of the experiment (Table 3, MANE-TSS strategy, |Delta|).
SCRAMBLE_ABS_DELTA = {
    "SLC24A5": 0.503, "TYR": 0.337, "SLC45A2": 0.259, "TYRP1": 0.251,
    "DDB1": 0.240, "MFSD12": 0.240, "MC1R": 0.055, "HERC2": 0.018,
    "OCA2": 0.002, "EDAR": 0.002, "TCHH": 0.000,
}
SCRAMBLE_DELTA_IN = {
    "TYR": 54296.6, "DDB1": 34150.3, "TYRP1": 26589.5, "SLC45A2": 23730.6,
    "MFSD12": 11868.3, "SLC24A5": 6537.0, "OCA2": 2109.0, "MC1R": 1391.0,
    "HERC2": 798.6, "EDAR": 59.0, "TCHH": 53.0,
}


def _log(m: str) -> None:
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def build_ctx(config_path: str, device):
    from dotenv import load_dotenv
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import build_model, load_checkpoint
    from genomics.predictors.genotype_based.config import (
        generate_experiment_name, get_dataset_cache_dir, get_experiment_runs_dir, load_config)
    from genomics.predictors.genotype_based.data.pipeline import (
        _make_runtime_processed_datasets, _resolve_runtime_dataset_dir)

    load_dotenv(Path.home() / ".env")
    cfg = load_config(REPO_ROOT / config_path)
    cache_dir = get_dataset_cache_dir(cfg)
    ds_dir = _resolve_runtime_dataset_dir(cfg)
    _log(f"Loading dataset ({ds_dir}, cache={cache_dir})...")
    with quiet_pipeline_logs():
        full_ds, train_ds, val_ds, test_ds = _make_runtime_processed_datasets(ds_dir, cache_dir, cfg)
    exp = get_experiment_runs_dir(cfg) / generate_experiment_name(cfg)
    ck = exp / "models" / "best_accuracy.pt"
    if not ck.exists():
        raise FileNotFoundError(f"checkpoint not found: {ck}")
    _log(f"Loading checkpoint {ck}")
    model = build_model(cfg, full_ds, device)
    model = load_checkpoint(model, ck, device)
    model.eval()

    idx_to_target = full_ds.idx_to_target
    strong = [i for i, n in idx_to_target.items() if n == "strong pigmentation"][0]
    weak = next(i for i in range(2) if i != strong)
    return cfg, full_ds, test_ds, model, strong, weak


def gene_block_slicer(sample_shape, genes, layout):
    """Return (rows_per_gene, slicer), where slicer(batched_tensor, k) is a *view*
    of gene k's block.

    Two layouts are in use and they stack differently, but in both the gene axis
    is the second-to-last one of a per-individual tensor, so the slicer indexes
    `dim -2` and works on a batched tensor unchanged:

      haplotype_channels (DITA)  (2, n_genes*C, L)   -- haplotype-major; within a
                                 haplotype, genes are concatenated along channels.
      raw_center_crop            (n_genes*2*C, L)    -- a flat vstack in
                                 gene-major, haplotype-minor order.

    Both orderings follow `genes_to_use`, which is what the dataset iterates over
    when it builds the tensor, so the mapping is exact rather than inferred.
    """
    n = len(genes)
    expected_rank = {"haplotype_channels": 3, "raw_center_crop": 2}
    if layout not in expected_rank:
        raise ValueError(f"unsupported tensor_layout {layout!r}")
    if len(sample_shape) != expected_rank[layout]:
        raise ValueError(f"unexpected {layout} per-individual shape {sample_shape}")
    gene_axis = sample_shape[-2]
    if gene_axis % n:
        raise ValueError(f"gene axis {gene_axis} not divisible by {n} genes ({layout})")
    per = gene_axis // n
    return per, (lambda t, k: t[..., k * per:(k + 1) * per, :])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--limit", type=int, default=0, help="use only the first N test individuals")
    ap.add_argument("--out", default=str(OUT_PATH))
    ap.add_argument("--lambdas", help="comma-separated lambda grid (default: see LAMBDA_GRID)")
    ap.add_argument("--device", choices=("auto", "cuda", "cpu"), default="auto",
                    help="auto picks cuda when available. Use cpu when a training "
                         "run holds the GPU: this job is forward-pass only, so it "
                         "runs on CPU in tens of minutes rather than failing.")
    ap.add_argument("--batch", type=int, default=64)
    args = ap.parse_args()

    lambdas = ([float(x) for x in args.lambdas.split(",")] if args.lambdas else list(LAMBDA_GRID))
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    _log(f"device = {device}")

    cfg, full_ds, test_ds, model, strong_idx, weak_idx = build_ctx(args.config, device)
    genes = list(cfg.dataset_input.genes_to_use)
    layout = cfg.dataset_input.tensor_layout
    _log(f"{len(genes)} genes, layout={layout}, weak class index = {weak_idx}")

    # Materialise the test split once. These are the *final* normalised tensors
    # the classifier consumes, so the ablation is applied at exactly the point
    # the paper's Delta is defined.
    n = len(test_ds) if not args.limit else min(args.limit, len(test_ds))
    _log(f"Materialising {n} test individuals...")
    X = torch.stack([test_ds[i][0] for i in range(n)]).to(device)
    _log(f"Test tensor: {tuple(X.shape)}")

    per, block_of = gene_block_slicer(tuple(X.shape[1:]), genes, layout)
    _log(f"{per} rows/channels per gene block")

    # Per-gene block L1, per individual and in aggregate. This is the denominator
    # that makes 'delivered' an identity rather than an estimate.
    block_l1 = {}
    for k, g in enumerate(genes):
        b = block_of(X, k)
        block_l1[g] = float(b.abs().sum().item()) / n          # mean per individual
    order = sorted(genes, key=lambda g: -block_l1[g])
    _log("block L1 per individual, descending: " +
         ", ".join(f"{g}={block_l1[g]:.1f}" for g in order))

    # Correctness guard on the gene->block mapping. A wrong axis silently gives
    # one gene the whole tensor and every other gene nothing, which is exactly
    # what happened on the first run of this script, so check rather than trust.
    dead = [g for g in genes if block_l1[g] <= 0.0]
    if dead:
        raise RuntimeError(
            f"{len(dead)} gene block(s) have zero L1 ({', '.join(dead)}). The "
            f"gene->block mapping for layout {layout!r} is wrong; refusing to "
            f"produce numbers from it.")
    ref = {g: SCRAMBLE_DELTA_IN[g] for g in genes if g in SCRAMBLE_DELTA_IN}
    if len(ref) >= 5:
        from scipy.stats import spearmanr
        common = [g for g in genes if g in ref]
        rho, p = spearmanr([block_l1[g] for g in common], [ref[g] for g in common])
        _log(f"guard: block L1 vs published |Delta_in| over {len(common)} genes: "
             f"rho={rho:.3f} (p={p:.4f}). A near-zero or negative rho would mean "
             f"the gene->block mapping is scrambled.")

    @torch.no_grad()
    def _forward(t):
        out = []
        for i in range(0, t.shape[0], args.batch):
            out.append(model(t[i:i + args.batch]))
        return torch.cat(out)

    def p_weak(t):
        return torch.softmax(_forward(t), dim=1)[:, weak_idx]

    def readouts(t):
        """P(weak) and the log-odds of weak-vs-strong, from one forward pass.

        The probability scale compresses non-uniformly: dP/dlogit = p(1-p), so the same shift in
        the decision function registers as a large Delta for an individual near the boundary and
        as nothing for a saturated one. Baseline confidence is ancestry-structured in this cohort,
        which makes a P(weak) mean an ancestry-weighted mean. The log-odds is a linear functional
        of the penultimate representation and carries no such weighting, so both are recorded and
        the paper can choose."""
        logits = _forward(t)
        strong_idx = 1 - weak_idx
        return (torch.softmax(logits, dim=1)[:, weak_idx],
                logits[:, weak_idx] - logits[:, strong_idx])

    base, base_lo = readouts(X)
    _log(f"baseline mean P(weak) = {base.mean().item():.4f}; "
         f"predicted weak: {(base > 0.5).sum().item()}/{n}")

    def ablate(k, lam):
        Xa = X.clone()
        block_of(Xa, k).mul_(lam)
        pw, lo = readouts(Xa)
        d = pw - base
        dl = lo - base_lo
        del Xa
        return {
            "lambda": float(lam),
            "delivered": float((1.0 - lam) * block_l1[genes[k]]),
            "delta": float(d.mean().item()),
            "mean_abs_delta": float(d.abs().mean().item()),
            "delta_logodds": float(dl.mean().item()),
            "mean_abs_delta_logodds": float(dl.abs().mean().item()),
            "median_abs_delta_logodds": float(dl.abs().median().item()),
            "flip_rate": float(((pw > 0.5) != (base > 0.5)).float().mean().item()),
        }

    # ---- (a) per-gene response curve, its own scale ----------------------
    curves = {}
    for k, g in enumerate(genes):
        pts = [{"lambda": 1.0, "delivered": 0.0, "delta": 0.0,
                "mean_abs_delta": 0.0, "delta_logodds": 0.0,
                "mean_abs_delta_logodds": 0.0, "median_abs_delta_logodds": 0.0,
                "flip_rate": 0.0}]
        pts += [ablate(k, lam) for lam in lambdas if lam != 1.0]
        curves[g] = pts
        _log(f"{g:<8} full ablation (lambda=0): delta={pts[-1]['delta']:+.4f}, "
             f"delivered={pts[-1]['delivered']:.1f}")

    # ---- (b) EXACT magnitude matching, no interpolation -------------------
    # For a shared budget B, scale gene g by lambda_g = 1 - B/||block_g||_1, so
    # every gene receives exactly B units of absolute change. This is what the
    # post-hoc log-log normalisation of the scramble readout approximates by
    # fitting; here it is arranged by construction, which is the whole point.
    #
    # B is capped at the smallest block L1, since a gene cannot deliver more
    # than it has. That cap is set by the dead-zone genes (EDAR, TCHH), which is
    # itself the finding: no budget both they and the responding genes can reach
    # is large enough to move the classifier.
    b_all = min(block_l1.values())
    smallest = min(block_l1, key=block_l1.get)
    _log("")
    _log(f"largest budget EVERY gene can deliver: {b_all:.1f} units, capped by {smallest}")

    # The all-gene cap is set by the dead-zone genes and is so small that the
    # comparison at it is essentially a read of the initial slope. That is a
    # legitimate quantity, but it is not the one a reader wants for the genes
    # that actually respond, so also report the cap computed over the responding
    # subset -- the genes that clear the random-window null screen. Which genes
    # are excluded, and why, is recorded in the output.
    RESPONDING = [g for g in genes if SCRAMBLE_ABS_DELTA.get(g, 0.0) >= 0.05]
    b_resp = min(block_l1[g] for g in RESPONDING) if RESPONDING else b_all
    excluded = [g for g in genes if g not in RESPONDING]
    _log(f"largest budget every RESPONDING gene can deliver: {b_resp:.1f} units "
         f"(excludes {', '.join(excluded)})")

    budget_sets = {
        "all_genes": ([b_all, b_all / 2, b_all / 10], genes),
        "responding_only": ([b_resp, b_resp / 2, b_resp / 10], RESPONDING),
    }
    matched = {}
    for label, (budgets, subset) in budget_sets.items():
        matched[label] = {}
        for B in budgets:
            row = {}
            for g in subset:
                k = genes.index(g)
                lam = 1.0 - B / block_l1[g]
                r = ablate(k, lam)
                r["budget"] = float(B)
                row[g] = r
            matched[label][f"{B:.4f}"] = row

    for label in budget_sets:
        _log("")
        _log(f"=== magnitude-matched ranking [{label}] "
             f"(every gene delivered the SAME budget) ===")
        for B, row in matched[label].items():
            rank = sorted(row.items(), key=lambda kv: -abs(kv[1]["delta"]))
            _log(f"budget {B} [P(weak)] : " + ", ".join(f"{g}={r['delta']:+.5f}" for g, r in rank))
            rank_lo = sorted(row.items(), key=lambda kv: -abs(kv[1]["delta_logodds"]))
            _log(f"budget {B} [log-odds]: " + ", ".join(f"{g}={r['delta_logodds']:+.5f}" for g, r in rank_lo))

    top_scr = max(SCRAMBLE_ABS_DELTA, key=SCRAMBLE_ABS_DELTA.get)
    top_row = matched["responding_only"][f"{b_resp:.4f}"]
    top_abl = max(top_row, key=lambda g: abs(top_row[g]["delta"]))

    # The full ranking comparison, not just the top gene: agreement at the top
    # with disagreement below is a different finding from agreement throughout.
    from scipy.stats import spearmanr
    common = [g for g in top_row if g in SCRAMBLE_ABS_DELTA]
    abl = [abs(top_row[g]["delta"]) for g in common]

    # The matched ablation delivers a common budget to every gene; the raw scramble |Delta|
    # does not, and readout_normalisation.json measures its dependence on delivered
    # magnitude at rho = 0.773. Correlating the two is therefore confounded by construction
    # and is NOT reported: the scramble has to be put on a per-unit-delivered footing first.
    rho_rank, p_rank = spearmanr(abl, [SCRAMBLE_ABS_DELTA[g] for g in common])

    ratio = [SCRAMBLE_ABS_DELTA[g] / SCRAMBLE_DELTA_IN[g] for g in common]
    rho_ratio, p_ratio = spearmanr(abl, ratio)
    loo = json.loads((OUT_PATH.parent / "readout_normalisation_loo.json").read_text())
    resid = [loo["loo"][g]["held_out_residual_log10"] for g in common]
    rho_resid, p_resid = spearmanr(abl, resid)
    _log("")
    _log(f"matched ablation vs scramble, ranking over the {len(common)} responding genes:")
    _log(f"  per unit delivered : rho = {rho_ratio:+.3f} (p = {p_ratio:.3f})")
    _log(f"  leave-one-out resid: rho = {rho_resid:+.3f} (p = {p_resid:.3f})")
    _log(f"  [not reported, confounded] vs raw |Delta|: rho = {rho_rank:+.3f} (p = {p_rank:.3f})")
    _log("  two defensible normalisations disagree at n = 7: no ranking claim is made.")
    _log("")
    _log(f"scramble top gene: {top_scr}; magnitude-matched ablation top gene: {top_abl}")
    _log("The pre-registered reading: agreement means the efficiency claim for "
         "SLC24A5 survives without depending on the six-point log-log fit; "
         "disagreement means it was an artefact of that normalisation.")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "generated": datetime.now(timezone.utc).isoformat(),
        "config": args.config,
        "layout": layout,
        "genes": genes,
        "n_individuals": n,
        "rows_per_gene_block": per,
        "lambda_grid": lambdas,
        "baseline_mean_p_weak": float(base.mean().item()),
        "baseline_mean_logodds": float(base_lo.mean().item()),
        "baseline_saturation": {
            "fraction_p1mp_below_0p01": float(((base * (1 - base)) < 0.01).float().mean().item()),
            "fraction_p1mp_below_0p001": float(((base * (1 - base)) < 0.001).float().mean().item()),
        },
        "block_l1_per_individual": block_l1,
        "curves": curves,
        "matched_exact": matched,
        "matched_budget_cap_all_genes": float(b_all),
        "matched_budget_cap_responding": float(b_resp),
        "responding_genes": RESPONDING,
        "excluded_from_responding_cap": excluded,
        "scramble_vs_ablation_rank_spearman": {
            "note": "vs raw |Delta|; confounded by delivered magnitude, not reported in the paper",
            "rho_vs_raw": float(rho_rank), "p_vs_raw": float(p_rank),
            "rho_vs_per_unit_delivered": float(rho_ratio), "p_vs_per_unit_delivered": float(p_ratio),
            "rho_vs_loo_residual": float(rho_resid), "p_vs_loo_residual": float(p_resid)},
        "scramble_abs_delta": SCRAMBLE_ABS_DELTA,
        "scramble_delta_in": SCRAMBLE_DELTA_IN,
        "note": "delivered = (1-lambda) * mean per-individual block L1, an identity "
                "rather than a measurement. matched_exact sets lambda_g = 1 - B/L1_g "
                "per gene so every gene receives exactly budget B; there is no "
                "interpolation and no fit anywhere in this file.",
    }, indent=2), encoding="utf-8")
    _log(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
