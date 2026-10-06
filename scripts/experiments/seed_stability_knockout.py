#!/usr/bin/env python3
"""E22 -- is the knockdown readout repeatable across scramble draws?

Every cell of the paper's knockdown table rests on ONE random permutation of the target 100 bp.
The permutation is the instrument's only stochastic component, so if a different permutation of
the same bases gives a materially different Delta, the readout is not a measurement and no effect
size in the paper carries meaning. Variance ACROSS INDIVIDUALS is already available (n=162);
variance ACROSS DRAWS is not measured anywhere.

Design. The target position is held FIXED at exactly the value the headline run used -- read from
`pigmentation_test_split_knockout.csv` rather than recomputed -- so the permutation seed is the
only thing that differs between draws. The headline run called `apply_scramble` with its default
seed=0 and wrote a cache key with no seed in it, so the existing data IS draw 0; this script adds
draws 1..N with the seed present in the cache key (without which draw k would silently read back
draw 0's cached prediction).

Scope. The genes that clear the pooled null, and a deterministic subsample of individuals: the
statistic of interest is the between-draw spread of the per-gene COHORT MEAN, and a 60-individual
mean estimates that well enough to see whether it moves. 4 extra draws x 7 genes x 60 individuals
= 1680 rows, ~3360 AlphaGenome calls, ~4 h at the observed 0.12 rows/s, ~7 GB of cache.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/seed_stability_knockout.py
  ... --limit-individuals 2 --genes TYR --draws 1     # smoke test
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _path in (REPO_ROOT / "src", REPO_ROOT / "notebooks", REPO_ROOT / "scripts" / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd

KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
DEFAULT_IN = KO_DIR / "pigmentation_test_split_knockout.csv"
DEFAULT_OUT = KO_DIR / "pigmentation_test_split_knockout_seeds.csv"
METHOD = "biology_tss"
ESTABLISHED = ["SLC24A5", "TYR", "SLC45A2", "TYRP1", "DDB1", "MFSD12", "MC1R"]
MIN_FREE_GB = 15.0
SUBSAMPLE_SEED = 13

CSV_FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "method",
    "draw", "scramble_seed", "scramble_window", "chrom",
    "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit",
    "perturbed_strong_logit", "perturbed_weak_logit",
    "baseline_pred", "perturbed_pred", "flipped",
]


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)
def _free_gb(p: Path) -> float: return shutil.disk_usage(p).free / 2**30


def seeded_knockout(ctx, sample_id, gene, targets, seed, scramble_window=100):
    """Identical to knockout_gene_experiment for method=biology_tss, except the scramble seed is
    explicit and enters the prediction cache key."""
    from alphagenome.data import genome
    from genotype_cnn_alignment_deeplift_summary.haplotype import load_haplotype_fasta
    from genotype_cnn_alignment_deeplift_summary.knockout import (
        apply_scramble, build_individual_tensor, load_raw_prediction, predict_modified_sequence,
        reorder_to_canonical, run_cnn_logits,
    )

    window_meta = json.loads(
        (ctx.dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text())
    chrom = window_meta["chromosome"]
    start_1based = int(window_meta["start"])
    full_ref_length = int(window_meta["end"]) - start_1based + 1
    interval = genome.Interval(chromosome=chrom, start=start_1based - 1,
                               end=start_1based - 1 + full_ref_length)

    hap_tracks = {}
    for haplotype in ("H1", "H2"):
        target = int(targets[haplotype])
        seq = load_haplotype_fasta(ctx.dataset_dir, sample_id, gene, haplotype)
        modified_seq, _o, _s, _st, _en = apply_scramble(
            seq, target, window_size=scramble_window, seed=seed)
        orig_array, orig_meta = load_raw_prediction(ctx.dataset_dir, sample_id, gene, haplotype)
        canonical_order = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
        # seed IS part of the key: without it draw k reads back draw 0's cached prediction.
        cache_key = (f"{sample_id}_{gene}_{haplotype}_{METHOD}_{target}"
                     f"_scramble{scramble_window}_seed{seed}")
        mod_values, mod_meta = predict_modified_sequence(ctx, cache_key, modified_seq, interval)
        hap_tracks[haplotype] = (
            reorder_to_canonical(mod_values, pd.DataFrame(mod_meta), canonical_order), orig_meta)

    baseline_tensor = build_individual_tensor(ctx, sample_id)
    perturbed_tensor = build_individual_tensor(
        ctx, sample_id, overrides={(gene, h): hap_tracks[h] for h in ("H1", "H2")})
    return {
        "chrom": chrom,
        "baseline_logits": run_cnn_logits(ctx, baseline_tensor).cpu().numpy(),
        "perturbed_logits": run_cnn_logits(ctx, perturbed_tensor).cpu().numpy(),
        "scramble_window": scramble_window,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_path", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--draws", type=int, default=4, help="ADDITIONAL draws; the headline run is draw 0")
    ap.add_argument("--limit-individuals", type=int, default=60)
    ap.add_argument("--genes", type=str, default=",".join(ESTABLISHED))
    args = ap.parse_args()

    import torch
    from bulk_knockout_pigmentation import _build_context

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"Using device: {device}")
    ctx, ko_genes, test_sample_ids, pedigree, class_names, strong_idx, weak_idx = _build_context(device)

    genes = [g for g in args.genes.split(",") if g]
    unknown = [g for g in genes if g not in ko_genes]
    if unknown:
        raise SystemExit(f"ABORT: unknown genes {unknown}")

    src = pd.read_csv(args.in_path)
    src = src[src["method"] == METHOD]
    targets = {(r["sample_id"], r["gene"]): {"H1": r["h1_target_local_idx"],
                                             "H2": r["h2_target_local_idx"]}
               for _, r in src.iterrows()}

    ids = [s for s in test_sample_ids if (s, genes[0]) in targets]
    if args.limit_individuals and args.limit_individuals < len(ids):
        rng = np.random.default_rng(SUBSAMPLE_SEED)
        ids = sorted(rng.choice(ids, size=args.limit_individuals, replace=False).tolist())

    total = len(ids) * len(genes) * args.draws
    _log(f"Plan: {len(ids)} individuals x {len(genes)} genes x {args.draws} extra draws = {total} rows "
         f"(~{2*total} AlphaGenome calls). Free disk: {_free_gb(REPO_ROOT):.1f} GB")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if args.out.exists():
        prev = pd.read_csv(args.out)
        done = set(zip(prev["sample_id"], prev["gene"], prev["draw"]))
        _log(f"Resuming: {len(done)} rows present.")
    write_header = not args.out.exists()

    n_written = n_skipped = n_failed = 0
    t0 = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        if write_header:
            w.writeheader(); fh.flush()
        for sample_id in ids:
            for gene in genes:
                for draw in range(1, args.draws + 1):
                    if (sample_id, gene, draw) in done:
                        n_skipped += 1
                        continue
                    if (sample_id, gene) not in targets:
                        n_skipped += 1
                        continue
                    if _free_gb(REPO_ROOT) < MIN_FREE_GB:
                        _log(f"ABORT: free disk below {MIN_FREE_GB} GB. {n_written} written; rerun to resume.")
                        return 2
                    try:
                        r = seeded_knockout(ctx, sample_id, gene, targets[(sample_id, gene)], seed=draw)
                        bl, pl = r["baseline_logits"], r["perturbed_logits"]
                        bp, pp = class_names[int(bl.argmax())], class_names[int(pl.argmax())]
                        ped_row = pedigree.get(sample_id, {})
                        w.writerow({
                            "sample_id": sample_id,
                            "population": ped_row.get("population"),
                            "superpopulation": ped_row.get("superpopulation"),
                            "true_label": ctx.full_ds._get_target_value(ped_row),
                            "gene": gene, "method": METHOD, "draw": draw, "scramble_seed": draw,
                            "scramble_window": r["scramble_window"], "chrom": r["chrom"],
                            "h1_target_local_idx": int(targets[(sample_id, gene)]["H1"]),
                            "h2_target_local_idx": int(targets[(sample_id, gene)]["H2"]),
                            "baseline_strong_logit": float(bl[ctx.strong_idx]),
                            "baseline_weak_logit": float(bl[ctx.weak_idx]),
                            "perturbed_strong_logit": float(pl[ctx.strong_idx]),
                            "perturbed_weak_logit": float(pl[ctx.weak_idx]),
                            "baseline_pred": bp, "perturbed_pred": pp, "flipped": bp != pp,
                        })
                        fh.flush()
                        n_written += 1
                        if n_written % 25 == 0:
                            el = time.monotonic() - t0
                            _log(f"progress: {n_written} written, {n_skipped} skipped, {n_failed} failed "
                                 f"({n_written/el:.2f} rows/s, free={_free_gb(REPO_ROOT):.1f} GB, "
                                 f"last={sample_id}/{gene}/d{draw})")
                    except Exception as exc:  # noqa: BLE001
                        n_failed += 1
                        _log(f"FAILED {sample_id}/{gene}/d{draw}: {exc!r}")
    _log(f"Done. {n_written} written, {n_skipped} skipped, {n_failed} failed. Output: {args.out}")
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
