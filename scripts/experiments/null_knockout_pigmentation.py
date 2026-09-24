#!/usr/bin/env python3
"""E5 -- null-scramble baseline for the promoter-knockout experiment.

Identical to scripts/experiments/bulk_knockout_pigmentation.py in every respect except WHERE the
100 bp window is scrambled: instead of a promoter location, the window centre is drawn uniformly
at random from the same gene window, excluding a radius around the annotated MANE Select TSS so a
"null" draw can never land on the promoter it is meant to be a null for. Same scramble procedure,
same window size, same shuffle seed, same CNN, same individuals -- only the location differs, so
the resulting |delta| distribution is the reference the real per-gene effects must be read against.

Deliberately implemented as a separate script rather than a new method inside
`knockout_gene_experiment`: that function produced the paper's headline result and is left
byte-identical. The shared primitives (apply_scramble, predict_modified_sequence,
build_individual_tensor, run_cnn_logits) are imported, not copied, so the intervention itself is
the same code.

One draw per (individual, gene) by default -> 162 x 11 = 1782 rows, 3564 AlphaGenome calls,
matching the per-gene n of each real method. Resumable per row; predictions are disk-cached and
reused, so an interrupted run is never re-billed.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/null_knockout_pigmentation.py
  ... --limit-individuals 2 --genes MC1R,TYR --out /tmp/null_smoke.csv    # smoke test
"""
from __future__ import annotations

import argparse
import csv
import hashlib
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

CSV_FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "method",
    "scramble_window", "chrom", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit", "perturbed_strong_logit", "perturbed_weak_logit",
    "baseline_pred", "perturbed_pred", "flipped",
    "draw", "exclude_radius", "h1_dist_to_tss", "h2_dist_to_tss",
]
METHOD = "random_window"
DEFAULT_OUT = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_null_scramble.csv"
MIN_FREE_GB = 15.0


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


def _free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 2**30


def pick_random_target(seq_len, tss_local_idx, exclude_radius, window_size, seed):
    """Uniform centre in [half, seq_len-half) at least `exclude_radius` from the TSS.

    Rejection sampling with a hard cap, then a deterministic fallback that walks away from the
    TSS -- so a gene whose TSS sits near the middle of its window can never spin forever.
    """
    rng = np.random.default_rng(seed)
    half = window_size // 2
    lo, hi = half, seq_len - half
    if hi <= lo:
        raise ValueError(f"window too short: seq_len={seq_len}")
    for _ in range(1000):
        c = int(rng.integers(lo, hi))
        if abs(c - tss_local_idx) >= exclude_radius:
            return c
    left_ok = tss_local_idx - exclude_radius - 1
    right_ok = tss_local_idx + exclude_radius + 1
    if left_ok >= lo:
        return int(left_ok)
    if right_ok < hi:
        return int(right_ok)
    raise ValueError(f"no admissible window for seq_len={seq_len}, tss={tss_local_idx}, r={exclude_radius}")


def _seed_for(sample_id, gene, haplotype, draw):
    h = hashlib.sha256(f"{sample_id}|{gene}|{haplotype}|{draw}".encode()).digest()
    return int.from_bytes(h[:8], "big")


def null_gene_experiment(ctx, sample_id, gene, draw, exclude_radius, scramble_window=100):
    from alphagenome.data import genome
    from genotype_cnn_alignment_deeplift_summary.annotations import get_gene_tss
    from genotype_cnn_alignment_deeplift_summary.haplotype import haplotype_local_idx, load_haplotype_fasta
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

    _, tss_pos_0based, _ = get_gene_tss(gene, ctx.tss_df_mane, ctx.tss_df_coding)

    hap_tracks, targets, dists = {}, {}, {}
    for haplotype in ("H1", "H2"):
        seq = load_haplotype_fasta(ctx.dataset_dir, sample_id, gene, haplotype)
        tss_hap_local_idx = haplotype_local_idx(
            ctx.dataset_dir, sample_id, gene, haplotype, start_1based, tss_pos_0based + 1)
        target = pick_random_target(
            len(seq), tss_hap_local_idx, exclude_radius, scramble_window,
            _seed_for(sample_id, gene, haplotype, draw))
        targets[haplotype] = target
        dists[haplotype] = int(target - tss_hap_local_idx)

        modified_seq, _orig_seg, _scr_seg, _s, _e = apply_scramble(
            seq, target, window_size=scramble_window)
        orig_array, orig_meta = load_raw_prediction(ctx.dataset_dir, sample_id, gene, haplotype)
        canonical_order = [(m["ontology_curie"], m["strand"]) for m in orig_meta]

        cache_key = f"{sample_id}_{gene}_{haplotype}_{METHOD}d{draw}_{target}_scramble{scramble_window}"
        mod_values, mod_meta_records = predict_modified_sequence(ctx, cache_key, modified_seq, interval)
        mod_reordered = reorder_to_canonical(mod_values, pd.DataFrame(mod_meta_records), canonical_order)
        hap_tracks[haplotype] = (mod_reordered, orig_meta)

    baseline_tensor = build_individual_tensor(ctx, sample_id)
    overrides = {(gene, hap): hap_tracks[hap] for hap in ("H1", "H2")}
    perturbed_tensor = build_individual_tensor(ctx, sample_id, overrides=overrides)
    baseline_logits = run_cnn_logits(ctx, baseline_tensor)
    perturbed_logits = run_cnn_logits(ctx, perturbed_tensor)
    return {
        "chrom": chrom, "targets": targets, "dists": dists,
        "baseline_logits": baseline_logits.cpu().numpy(),
        "perturbed_logits": perturbed_logits.cpu().numpy(),
        "scramble_window": scramble_window,
    }


def _already_done(out_path: Path) -> set:
    if not out_path.exists():
        return set()
    existing = pd.read_csv(out_path)
    if "draw" not in existing.columns:
        return set()
    return set(zip(existing["sample_id"], existing["gene"], existing["draw"]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--draws", type=int, default=1, help="random windows per (individual, gene)")
    ap.add_argument("--exclude-radius", type=int, default=5000,
                    help="bp around the annotated TSS that a null draw may not land in")
    ap.add_argument("--limit-individuals", type=int, default=None)
    ap.add_argument("--genes", type=str, default=None)
    args = ap.parse_args()

    import torch
    from bulk_knockout_pigmentation import _build_context

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"Using device: {device}")
    ctx, ko_genes, test_sample_ids, pedigree, class_names, strong_idx, weak_idx = _build_context(device)

    genes = [g for g in ko_genes if g in args.genes.split(",")] if args.genes else ko_genes
    sample_ids = test_sample_ids[: args.limit_individuals] if args.limit_individuals else test_sample_ids
    total = len(sample_ids) * len(genes) * args.draws
    _log(f"Plan: {len(sample_ids)} individuals x {len(genes)} genes x {args.draws} draws = {total} rows "
         f"(~{2*total} AlphaGenome calls) -> {args.out}")
    _log(f"Excluding +/-{args.exclude_radius} bp around the annotated TSS. Free disk: {_free_gb(REPO_ROOT):.1f} GB")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = _already_done(args.out)
    if done:
        _log(f"Resuming: {len(done)} rows already present.")
    write_header = not args.out.exists()

    n_written = n_skipped = n_failed = 0
    t0 = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        if write_header:
            w.writeheader(); fh.flush()
        for sample_id in sample_ids:
            for gene in genes:
                for draw in range(args.draws):
                    if (sample_id, gene, draw) in done:
                        n_skipped += 1
                        continue
                    if _free_gb(REPO_ROOT) < MIN_FREE_GB:
                        _log(f"ABORT: free disk below {MIN_FREE_GB} GB. {n_written} rows written; rerun to resume.")
                        return 2
                    try:
                        r = null_gene_experiment(ctx, sample_id, gene, draw, args.exclude_radius)
                        bl, pl = r["baseline_logits"], r["perturbed_logits"]
                        bp, pp = class_names[int(bl.argmax())], class_names[int(pl.argmax())]
                        ped_row = pedigree.get(sample_id, {})
                        w.writerow({
                            "sample_id": sample_id,
                            "population": ped_row.get("population"),
                            "superpopulation": ped_row.get("superpopulation"),
                            "true_label": ctx.full_ds._get_target_value(ped_row),
                            "gene": gene, "method": METHOD,
                            "scramble_window": r["scramble_window"], "chrom": r["chrom"],
                            "h1_target_local_idx": r["targets"]["H1"],
                            "h2_target_local_idx": r["targets"]["H2"],
                            "baseline_strong_logit": float(bl[ctx.strong_idx]),
                            "baseline_weak_logit": float(bl[ctx.weak_idx]),
                            "perturbed_strong_logit": float(pl[ctx.strong_idx]),
                            "perturbed_weak_logit": float(pl[ctx.weak_idx]),
                            "baseline_pred": bp, "perturbed_pred": pp, "flipped": bp != pp,
                            "draw": draw, "exclude_radius": args.exclude_radius,
                            "h1_dist_to_tss": r["dists"]["H1"], "h2_dist_to_tss": r["dists"]["H2"],
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
