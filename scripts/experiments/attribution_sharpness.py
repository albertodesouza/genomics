#!/usr/bin/env python3
"""E18.3 -- attribution sharpness, DITA vs. no_alignment, both remapped to reference coordinates.

The direct test of the coordinate-coherence mechanism DITA is supposed to buy: if projecting onto
the shared axis removes indel-induced blurring, the aligned model's attribution should be
measurably *more concentrated* (sharper) than the unaligned model's, once both are compared on the
same footing. E18.1 (scripts/experiments/exon_promoter_enrichment.py) already established the
per-gene attribution "mass curve" (sum of |class-mean log-odds DeepLIFT attribution|, pooled over
every ontology/strand track) in each model's own local axis; this script reuses that construction
and adds the piece E18.1 didn't need: mapping DITA's curve *back* to reference coordinates so it is
directly comparable to `no_alignment`'s, which is already in reference coordinates by construction.

Why this is a fair comparison rather than an apples-to-oranges one. The DynamicIndelAligner's
reference sub-crop (verified empirically in exon_promoter_enrichment.py's investigation, gene=TYR:
`ref_start_offset=245760`) is centered by the exact same formula as `no_alignment`'s naive crop
(`window_genomic_axis`: `center_offset - size//2`, `center_offset=full_ref_length//2=262144`,
`262144-16384=245760`) -- both models' windows cover the *identical* 32,768 bp reference span for
every gene (cross-checked in E18.1's output: e.g. EDAR annotated exon+promoter footprint is exactly
374 bp under both models, HERC2 exactly 3041 bp under both). So remapping DITA's curve to this same
reference span, folding each insertion slot's attribution into the reference position it was
inserted after (the natural inverse of how `expanded_index_map`/`insertion_slots_by_ref` were built
-- see dita_box_to_local in exon_promoter_enrichment.py for the forward direction), gives two
length-32768 curves over the *same* chromosome, same start, same strand, comparable index for
index. The fold is mass-conserving: every local position in [0, window_size) is either an anchor's
own slot or an insertion slot immediately following exactly one anchor, so summing the two never
double-counts or drops mass.

Sharpness metrics (both computed on the reference-frame curve, so higher/lower has a shared meaning
across models):
  - top-k mass fraction, k in {1%, 5%, 10%} of the 32,768-position window: fraction of total |mass|
    contained in the k*L positions with the largest values. Higher = sharper.
  - participation ratio, expressed as a percentage of window length: (sum m)^2 / sum(m^2) / L * 100,
    the standard inverse-participation-ratio "effective number of positions carrying the mass",
    which needs no arbitrary top-k cutoff. Lower = sharper (fewer effective positions).

No AlphaGenome call anywhere in this script; no new checkpoints. Reuses the same two pigmentation
checkpoints (`dita_no_masks`, `no_alignment`) E18.1 already runs DeepLIFT on.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/attribution_sharpness.py
  ... --genes TYR,SLC24A5   # smoke test
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _path in (REPO_ROOT / "src", REPO_ROOT / "notebooks"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
os.chdir(REPO_ROOT)

import numpy as np

CONFIG_PATHS = {
    "no_alignment": REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment.yaml",
    "dita_no_masks": REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml",
}
TOP_K_FRACTIONS = (0.01, 0.05, 0.10)


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


def gene_mass_curve(gene, config, attr_by_target, target_idx, other_idx):
    """Identical construction to exon_promoter_enrichment.gene_mass_curve -- duplicated (not
    imported) so this script stays runnable standalone, matching the other scripts/experiments/
    convention. See that module's docstring for why this particular pooling was chosen."""
    from genotype_cnn_alignment_deeplift_summary.attribution import gene_track_curves

    logodds_target = attr_by_target[(target_idx, target_idx)] - attr_by_target[(target_idx, other_idx)]
    logodds_other = attr_by_target[(other_idx, other_idx)] - attr_by_target[(other_idx, target_idx)]
    curves_target = gene_track_curves(logodds_target, config, gene)
    curves_other = gene_track_curves(logodds_other, config, gene)
    keys = set(curves_target) | set(curves_other)
    length = next(iter(curves_target.values()), next(iter(curves_other.values()))).shape[0]
    mass = np.zeros(length, dtype=np.float64)
    for k in keys:
        if k in curves_target:
            mass += np.abs(curves_target[k])
        if k in curves_other:
            mass += np.abs(curves_other[k])
    return mass


def fold_dita_curve_to_reference(local_mass, axis, slice_info):
    """local_mass: 1D array of length slice_info['window_size'], indexed by DITA's own expanded
    local axis (as produced by gene_track_curves on a haplotype_channels tensor). Returns a 1D
    array of length axis['ref_length'], indexed by reference position (ref_idx, 0-based, relative
    to axis['alignment_start_1based']) -- the same reference span `no_alignment`'s curve already
    uses. Every insertion slot's mass is folded into the reference anchor it was inserted after."""
    ref_length = axis["ref_length"]
    expanded_index_map = axis["expanded_index_map"]  # {ref_idx: expanded_idx}
    insertion_slots_by_ref = axis["insertion_slots_by_ref"]  # {ref_idx: [expanded_idx, ...]}
    expanded_start = slice_info["expanded_start"]
    window_size = slice_info["window_size"]

    ref_curve = np.zeros(ref_length, dtype=np.float64)
    for ref_idx in range(ref_length):
        local_i = expanded_index_map[ref_idx] - expanded_start
        if 0 <= local_i < window_size:
            ref_curve[ref_idx] += local_mass[local_i]
        for slot_expanded in insertion_slots_by_ref.get(ref_idx, ()):
            slot_local = slot_expanded - expanded_start
            if 0 <= slot_local < window_size:
                ref_curve[ref_idx] += local_mass[slot_local]
    return ref_curve


def sharpness_metrics(curve, top_k_fractions=TOP_K_FRACTIONS):
    total = curve.sum()
    L = len(curve)
    if total <= 0:
        return {**{f"top{int(k*100)}pct_frac": float("nan") for k in top_k_fractions},
                "participation_ratio_pct": float("nan"), "total_mass": float(total)}
    sorted_desc = np.sort(curve)[::-1]
    cum = np.cumsum(sorted_desc)
    out = {}
    for k in top_k_fractions:
        n = max(1, int(round(k * L)))
        out[f"top{int(k*100)}pct_frac"] = float(cum[n - 1] / total)
    participation_ratio = (curve.sum() ** 2) / np.sum(curve ** 2)  # positions, not %
    out["participation_ratio_pct"] = float(100.0 * participation_ratio / L)
    out["total_mass"] = float(total)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=str, default=None, help="comma-separated subset, for smoke tests")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    import torch

    from genotype_cnn_alignment_deeplift_summary.attribution import (
        compute_attributions_by_target_class, compute_class_mean_deeplift,
    )
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import load_predictor_bundle

    from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
    from genomics.predictors.genotype_based.data.pipeline import _make_runtime_processed_datasets, _resolve_runtime_dataset_dir

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"device={device}")

    _log("Loading both pigmentation checkpoints (dita_no_masks, no_alignment)...")
    bundle = load_predictor_bundle(CONFIG_PATHS, device, class_names_from="dita_no_masks")
    configs, loaders, models, class_names = bundle.configs, bundle.loaders, bundle.models, bundle.class_names
    assert len(class_names) == 2
    target_idx = class_names.index("strong pigmentation") if "strong pigmentation" in class_names else 0
    other_idx = next(i for i in range(2) if i != target_idx)

    genes = list(configs["dita_no_masks"].dataset_input.genes_to_use)
    if args.genes:
        wanted = set(args.genes.split(","))
        genes = [g for g in genes if g in wanted]
    window_center_size = configs["dita_no_masks"].dataset_input.window_center_size

    _log("Rebuilding the DITA alignment axis (same config/cache as the trained run)...")
    dita_cfg = load_config(CONFIG_PATHS["dita_no_masks"])
    cache_dir = get_dataset_cache_dir(dita_cfg)
    runtime_dataset_dir = _resolve_runtime_dataset_dir(dita_cfg)
    with quiet_pipeline_logs():
        dita_full_ds, _, _, _ = _make_runtime_processed_datasets(runtime_dataset_dir, cache_dir, dita_cfg)
    aligner = dita_full_ds.dynamic_indel_aligner

    per_model_curves = {}
    for model_name in ("dita_no_masks", "no_alignment"):
        _log(f"=== {model_name}: DeepLIFT class-mean attribution over the test split ===")
        test_dataset = loaders[model_name]["test"].dataset
        deeplift_results = compute_class_mean_deeplift(
            models[model_name], test_dataset, (target_idx, other_idx), class_names, model_name=model_name,
        )
        attr_by_target = compute_attributions_by_target_class(
            models[model_name], test_dataset, target_idx, other_idx, deeplift_results, class_names,
        )
        curves = {}
        for gene in genes:
            local_mass = gene_mass_curve(gene, configs[model_name], attr_by_target, target_idx, other_idx)
            if model_name == "dita_no_masks":
                axis = aligner.get_alignment_axis(gene)
                slice_info = aligner.get_reference_centered_expanded_slice(gene, window_center_size)
                curves[gene] = fold_dita_curve_to_reference(local_mass, axis, slice_info)
            else:
                curves[gene] = local_mass
        per_model_curves[model_name] = curves

    gene_rows = []
    for gene in genes:
        dita_curve = per_model_curves["dita_no_masks"][gene]
        noalign_curve = per_model_curves["no_alignment"][gene]
        assert len(dita_curve) == len(noalign_curve), \
            f"{gene}: reference-frame length mismatch ({len(dita_curve)} vs {len(noalign_curve)})"
        dita_metrics = sharpness_metrics(dita_curve)
        noalign_metrics = sharpness_metrics(noalign_curve)
        row = {"gene": gene, "window_bp": len(dita_curve), "dita_no_masks": dita_metrics, "no_alignment": noalign_metrics}
        gene_rows.append(row)
        _log(f"  {gene}: DITA top5%={dita_metrics['top5pct_frac']:.4f} PR%={dita_metrics['participation_ratio_pct']:.2f}  |  "
             f"no_align top5%={noalign_metrics['top5pct_frac']:.4f} PR%={noalign_metrics['participation_ratio_pct']:.2f}")

    # -- summary across genes --
    def paired(metric):
        d = np.array([r["dita_no_masks"][metric] for r in gene_rows])
        n = np.array([r["no_alignment"][metric] for r in gene_rows])
        valid = ~(np.isnan(d) | np.isnan(n))
        return d[valid], n[valid]

    summary = {}
    for metric in ("top1pct_frac", "top5pct_frac", "top10pct_frac", "participation_ratio_pct"):
        d, n = paired(metric)
        if len(d) == 0:
            continue
        sharper_higher = metric != "participation_ratio_pct"
        n_dita_sharper = int((d > n).sum()) if sharper_higher else int((d < n).sum())
        try:
            from scipy.stats import wilcoxon
            stat, p = wilcoxon(d, n)
            p = float(p)
        except Exception as e:
            p = None
        summary[metric] = {
            "n_genes": int(len(d)), "n_dita_sharper": n_dita_sharper,
            "dita_mean": float(d.mean()), "no_alignment_mean": float(n.mean()),
            "wilcoxon_p": p,
        }
        _log(f"SUMMARY {metric}: DITA sharper in {n_dita_sharper}/{len(d)} genes "
             f"(mean DITA={d.mean():.4f} vs no_align={n.mean():.4f}, Wilcoxon p={p})")

    out = args.out or REPO_ROOT / "results/genotype_based_predictor/attribution_sharpness.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"genes": genes, "per_gene": gene_rows, "summary": summary}, indent=2))
    _log(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
