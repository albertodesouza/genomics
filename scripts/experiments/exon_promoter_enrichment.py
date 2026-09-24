#!/usr/bin/env python3
"""E18.1 -- exon/promoter enrichment of DeepLIFT attribution mass against a permutation null.

For each of the two tracks-only pigmentation models (`dita_no_masks`, aligned; `no_alignment`,
unaligned) and each of the eleven genes, this measures whether the class-mean DeepLIFT attribution
toward the strong/weak pigmentation log-odds concentrates inside annotated exons and the promoter,
compared to a null built from length-matched intervals placed at uniformly random positions in the
same window (1000 draws per gene per model). This is the test the paper's Section
"Attribution along the shared axis" explicitly declines to report until it exists: a per-position
attribution curve is only interpretable as a claim about genomic position if it is measurably
non-uniform with respect to functional annotation, and until now that had not been checked for
either model.

Per-gene attribution curve. `compute_attributions_by_target_class` (attribution.py) gives four
DeepLIFT tensors per gene, keyed by (population, toward-logit); a population's own log-odds
attribution is attr[(p,p)] - attr[(p,other)]. This script sums the *absolute* per-track curves
(via `gene_track_curves`) of both populations' own log-odds attribution, pooled over every
(ontology, strand) track for that gene, into one non-negative "attribution mass" curve of length
`window_center_size` (32768). This is a deliberate pooling choice (documented, not hidden): it
scores where the model's attention concentrates, independent of which class it favours, which is
the property the exon/promoter question is actually about.

Coordinate mapping (the part that would silently invalidate the result if wrong -- see
in-line comments):
  - `no_alignment` (raw_center_crop): axis index i is exactly reference position
    `window_start_0based + i` (annotations.window_genomic_axis, confirmed exact for this layout).
    Exon boxes come directly from `annotations.gene_exon_local_boxes` in this same coordinate
    system; the promoter box is `TSS +/- 2000 bp` computed the same way.
  - `dita_no_masks` (haplotype_channels): the DynamicIndelAligner builds its expanded axis over
    only the central `window_center_size` reference bp (NOT the full 524,288 bp window -- verified
    empirically: for TYR, `ref_length=32768`, `alignment_start_1based=89220434`, which is exactly
    the gene's 524,288 bp window start plus `ref_start_offset=245760`). A genomic position maps to
    the model's local axis via `ref_idx = pos_1based - alignment_start_1based`, then
    `expanded_idx = expanded_index_map[ref_idx]` (insertions add extra slots after an anchor;
    deletions/unfilled slots are zero-filled, not removed -- the map is still every anchor's own
    fixed slot), then `local_i = expanded_idx - expanded_start` from
    `get_reference_centered_expanded_slice(gene, window_center_size)`. Exon/promoter boxes are
    obtained by calling `annotations.gene_exon_local_boxes`/TSS lookup against the gene's FULL
    524,288 bp window (so the returned "local" coordinates are 0-based offsets from that window's
    own start), converting to absolute genomic position, then applying the mapping above. A box
    entirely outside the aligner's 32,768 bp reference sub-crop (this happens for genes where the
    annotated TSS lies outside the naive center crop -- the paper already documents this for EDAR,
    HERC2, OCA2, SLC45A2, TYR) contributes zero width, which is reported per-gene, not hidden.
  - Simplification kept for both layouts: a box's mapped span is taken as [first mapped position,
    last mapped position], not per-base -- for DITA this ignores at most a handful of insertion
    slots trailing the box's last reference anchor (network-wide across TYR's whole window, indels
    add only 44 of 32812 expanded slots), a second-order effect on an already coarse ~kbp-scale box.

Permutation null. Per gene per model: take the merged (non-overlapping) annotated intervals'
individual widths, and on each of 1000 draws place each width at an independent uniformly random
start in [0, L - width), take the union of the drawn intervals (so overlaps between drawn intervals
are not double-counted, matching how the real annotated fraction is computed), and record the
fraction of total |mass| falling inside that union. The reported p-value is one-sided:
(1 + #{null draws >= observed}) / (1000 + 1).

No AlphaGenome call anywhere in this script -- everything after the two DeepLIFT passes (one per
model, each a full pass over the 162-individual test split) is array arithmetic on precomputed
attribution tensors. Runs on the existing `dita_no_masks` and `no_alignment` pigmentation
checkpoints.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/exon_promoter_enrichment.py
  ... --n-perm 1000 --seed 13
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

PROMOTER_HALF_WIDTH = 2000  # symmetric window around the annotated TSS; see module docstring

CONFIG_PATHS = {
    "no_alignment": REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment.yaml",
    "dita_no_masks": REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml",
}
GTF_CACHE_DIR = REPO_ROOT / "notebooks" / ".cache"


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


# --------------------------------------------------------------------------------------------
# Interval arithmetic
# --------------------------------------------------------------------------------------------

def merge_intervals(boxes):
    """[(s, e), ...] (possibly overlapping) -> sorted, merged, non-overlapping [(s, e), ...]."""
    boxes = sorted(b for b in boxes if b[1] > b[0])
    if not boxes:
        return []
    merged = [list(boxes[0])]
    for s, e in boxes[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def mass_in_intervals(prefix_sum, intervals):
    """prefix_sum[i] = sum(curve[:i]); intervals must already be merged/non-overlapping."""
    return sum(prefix_sum[e] - prefix_sum[s] for s, e in intervals)


def permutation_null(curve, widths, n_perm, rng):
    """Fraction of total |curve| mass falling in a union of len(widths) length-matched intervals
    placed at independent uniformly random positions, repeated n_perm times."""
    L = len(curve)
    prefix_sum = np.concatenate([[0.0], np.cumsum(curve)])
    total = prefix_sum[-1]
    fractions = np.empty(n_perm)
    for p in range(n_perm):
        boxes = []
        for w in widths:
            if w >= L:
                boxes.append((0, L))
                continue
            s = int(rng.integers(0, L - w + 1))
            boxes.append((s, s + w))
        fractions[p] = mass_in_intervals(prefix_sum, merge_intervals(boxes)) / total
    return fractions


# --------------------------------------------------------------------------------------------
# Coordinate mapping
# --------------------------------------------------------------------------------------------

def dita_box_to_local(axis, slice_info, abs_start_0based, abs_end_0based):
    """Absolute 0-based half-open genomic [abs_start_0based, abs_end_0based) -> local DITA axis
    [lo, hi) in [0, window_size), or None if the box has no overlap with the aligner's reference
    sub-crop. See module docstring for the derivation."""
    ref_lo = (abs_start_0based + 1) - axis["alignment_start_1based"]
    ref_hi_incl = abs_end_0based - axis["alignment_start_1based"]
    ref_lo = max(ref_lo, 0)
    ref_hi_incl = min(ref_hi_incl, axis["ref_length"] - 1)
    if ref_hi_incl < ref_lo:
        return None
    expanded_index_map = axis["expanded_index_map"]
    expanded_lo = expanded_index_map[ref_lo]
    expanded_hi_excl = expanded_index_map[ref_hi_incl] + 1
    local_lo = expanded_lo - slice_info["expanded_start"]
    local_hi = expanded_hi_excl - slice_info["expanded_start"]
    local_lo = max(local_lo, 0)
    local_hi = min(local_hi, slice_info["window_size"])
    if local_hi <= local_lo:
        return None
    return (local_lo, local_hi)


def no_alignment_box_to_local(window_start_0based, local_length, abs_start_0based, abs_end_0based):
    lo = max(abs_start_0based - window_start_0based, 0)
    hi = min(abs_end_0based - window_start_0based, local_length)
    if hi <= lo:
        return None
    return (lo, hi)


# --------------------------------------------------------------------------------------------
# Attribution mass curves
# --------------------------------------------------------------------------------------------

def gene_mass_curve(gene, config, deeplift_results, attr_by_target, target_idx, other_idx):
    """Sum of |per-track log-odds attribution curve|, pooled over every (ontology, strand) track,
    for both populations' own log-odds attribution -- see module docstring."""
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


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--genes", type=str, default=None, help="comma-separated subset, for smoke tests")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    import torch

    from genotype_cnn_alignment_deeplift_summary.annotations import (
        build_tss_tables, build_transcript_extractors, gene_exon_local_boxes, get_gene_tss, load_gtf,
        window_genomic_axis,
    )
    from genotype_cnn_alignment_deeplift_summary.attribution import (
        compute_attributions_by_target_class, compute_class_mean_deeplift,
    )
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import load_predictor_bundle

    from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
    from genomics.predictors.genotype_based.data.pipeline import _make_runtime_processed_datasets, _resolve_runtime_dataset_dir

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"device={device}")

    _log("Loading GENCODE v46 GTF (cached)...")
    gtf = load_gtf(GTF_CACHE_DIR)
    gtf_mane, gtf_protein_coding, transcript_extractor_mane, transcript_extractor_coding, gene_id_map = \
        build_transcript_extractors(gtf)
    tss_df_mane, tss_df_coding = build_tss_tables(gtf_mane, gtf_protein_coding)

    _log("Loading both pigmentation checkpoints (dita_no_masks, no_alignment)...")
    bundle = load_predictor_bundle(CONFIG_PATHS, device, class_names_from="dita_no_masks")
    configs, loaders, models, class_names = bundle.configs, bundle.loaders, bundle.models, bundle.class_names
    assert len(class_names) == 2
    target_idx = class_names.index("strong pigmentation") if "strong pigmentation" in class_names else 0
    other_idx = next(i for i in range(2) if i != target_idx)
    _log(f"class_names={class_names} target_idx={target_idx} other_idx={other_idx}")

    genes = list(configs["dita_no_masks"].dataset_input.genes_to_use)
    if args.genes:
        wanted = set(args.genes.split(","))
        genes = [g for g in genes if g in wanted]
    window_center_size = configs["dita_no_masks"].dataset_input.window_center_size
    dataset_dir = Path(configs["dita_no_masks"].dataset_input.dataset_dir)

    _log("Rebuilding the DITA alignment axis (same config/cache as the trained run)...")
    dita_cfg = load_config(CONFIG_PATHS["dita_no_masks"])
    cache_dir = get_dataset_cache_dir(dita_cfg)
    runtime_dataset_dir = _resolve_runtime_dataset_dir(dita_cfg)
    with quiet_pipeline_logs():
        dita_full_ds, _, _, _ = _make_runtime_processed_datasets(runtime_dataset_dir, cache_dir, dita_cfg)
    aligner = dita_full_ds.dynamic_indel_aligner

    rng = np.random.default_rng(args.seed)
    per_model_results = {}

    for model_name in ("dita_no_masks", "no_alignment"):
        _log(f"=== {model_name}: DeepLIFT class-mean attribution over the test split ===")
        test_dataset = loaders[model_name]["test"].dataset
        deeplift_results = compute_class_mean_deeplift(
            models[model_name], test_dataset, (target_idx, other_idx), class_names, model_name=model_name,
        )
        # Two more full DeepLIFT passes (cross-attribution), computed once per model -- NOT per
        # gene, since it scores the whole input tensor; gene_track_curves does the per-gene split.
        attr_by_target = compute_attributions_by_target_class(
            models[model_name], test_dataset, target_idx, other_idx, deeplift_results, class_names,
        )

        gene_rows = []
        for gene in genes:
            mass = gene_mass_curve(gene, configs[model_name], deeplift_results, attr_by_target, target_idx, other_idx)
            L = len(mass)

            # -- annotated boxes, in this model's own local axis coordinates --
            boxes = []
            note = ""
            if model_name == "no_alignment":
                chrom, window_start_0based, local_length = window_genomic_axis(dataset_dir, gene, window_center_size)
                exon_boxes, strand = gene_exon_local_boxes(
                    gene, chrom, window_start_0based, local_length, gene_id_map,
                    transcript_extractor_mane, transcript_extractor_coding,
                )
                boxes.extend(exon_boxes)
                tss_chrom, tss_pos_0based, tss_strand = get_gene_tss(gene, tss_df_mane, tss_df_coding)
                promoter = no_alignment_box_to_local(
                    window_start_0based, local_length,
                    tss_pos_0based - PROMOTER_HALF_WIDTH, tss_pos_0based + PROMOTER_HALF_WIDTH,
                )
                if promoter:
                    boxes.append(promoter)
                else:
                    note = "promoter (TSS+/-2kbp) entirely outside the 32768bp crop"
            else:
                axis = aligner.get_alignment_axis(gene)
                slice_info = aligner.get_reference_centered_expanded_slice(gene, window_center_size)
                # full 524,288bp window -> local boxes are absolute-0-based-genomic minus the
                # window's own start, i.e. we recover absolute genomic coords by re-adding it.
                meta = json.loads((dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text())
                chrom = meta["chromosome"]
                full_window_start_0based = int(meta["start"]) - 1
                full_ref_length = int(meta["end"]) - int(meta["start"]) + 1
                exon_boxes_full, strand = gene_exon_local_boxes(
                    gene, chrom, full_window_start_0based, full_ref_length, gene_id_map,
                    transcript_extractor_mane, transcript_extractor_coding,
                )
                n_dropped = 0
                for s, e in exon_boxes_full:
                    mapped = dita_box_to_local(axis, slice_info,
                                                full_window_start_0based + s, full_window_start_0based + e)
                    if mapped:
                        boxes.append(mapped)
                    else:
                        n_dropped += 1
                tss_chrom, tss_pos_0based, tss_strand = get_gene_tss(gene, tss_df_mane, tss_df_coding)
                promoter = dita_box_to_local(axis, slice_info,
                                              tss_pos_0based - PROMOTER_HALF_WIDTH, tss_pos_0based + PROMOTER_HALF_WIDTH)
                if promoter:
                    boxes.append(promoter)
                else:
                    note = "promoter (TSS+/-2kbp) entirely outside the 32768bp reference sub-crop"
                if n_dropped:
                    note = (note + "; " if note else "") + f"{n_dropped}/{len(exon_boxes_full)} exon boxes fell entirely outside the crop"

            merged = merge_intervals(boxes)
            annotated_bp = sum(e - s for s, e in merged)
            prefix_sum = np.concatenate([[0.0], np.cumsum(mass)])
            total_mass = prefix_sum[-1]
            observed_frac = mass_in_intervals(prefix_sum, merged) / total_mass if merged and total_mass > 0 else 0.0

            widths = [e - s for s, e in merged]
            if widths and total_mass > 0:
                null_fracs = permutation_null(mass, widths, args.n_perm, rng)
                p_value = (1 + int((null_fracs >= observed_frac).sum())) / (args.n_perm + 1)
                null_mean, null_p95 = float(null_fracs.mean()), float(np.quantile(null_fracs, 0.95))
            else:
                null_fracs = np.array([])
                p_value, null_mean, null_p95 = float("nan"), float("nan"), float("nan")

            row = {
                "gene": gene, "window_bp": L, "annotated_bp": annotated_bp, "n_boxes": len(merged),
                "observed_frac": observed_frac, "null_mean_frac": null_mean, "null_p95_frac": null_p95,
                "p_value": p_value, "note": note,
            }
            gene_rows.append(row)
            _log(f"  {gene}: annotated={annotated_bp}/{L}bp observed_frac={observed_frac:.4f} "
                 f"null_mean={null_mean:.4f} p={p_value:.4f}" + (f"  [{note}]" if note else ""))

        per_model_results[model_name] = gene_rows

    out = args.out or REPO_ROOT / "results/genotype_based_predictor/exon_promoter_enrichment.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "promoter_half_width_bp": PROMOTER_HALF_WIDTH, "n_perm": args.n_perm, "seed": args.seed,
        "genes": genes, "results": per_model_results,
    }, indent=2))
    _log(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
