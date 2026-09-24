#!/usr/bin/env python3
"""E6 -- intervention validity: did the promoter scramble actually change anything?

For every (individual, gene, method, haplotype) row of the bulk knockout CSV this measures two
distinct things, which the paper currently conflates in every 0.00 cell of Table II:

  1. Did AlphaGenome's prediction move at all?  Measured on the RAW per-haplotype prediction,
     integrated over the gene's exons (haplotype-local coordinates) and over the whole window.
  2. Did that movement REACH THE CNN?  Measured on the aligned+cropped signal rows that the CNN
     actually consumes, obtained by calling the real `_process_window_haplotype_channels`, so no
     assumption about where the crop lands can be wrong.

A gene can score zero in Table II for two very different reasons -- the classifier ignores it, or
the edit never entered its input at all -- and (2) separates them.

Read-only w.r.t. the repo and the AlphaGenome API: every modified prediction is taken from the
existing on-disk cache. Rows whose cache entry is missing are reported, never re-billed.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/knockdown_expression_delta.py
  ... --limit 50    # smoke test
"""
from __future__ import annotations

import argparse
import csv
import json
import os
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

DEFAULT_IN = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout.csv"
DEFAULT_OUT = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/knockdown_expression_delta.csv"
AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"

FIELDS = [
    "sample_id", "gene", "method", "haplotype", "target_local_idx",
    # (1) raw AlphaGenome movement
    "raw_abs_delta_total", "raw_abs_delta_exons", "raw_sum_orig_exons", "raw_sum_mod_exons",
    "raw_rel_delta_exons", "raw_n_changed_positions", "raw_first_changed", "raw_last_changed",
    # (2) movement that survives alignment + crop into the CNN input
    "crop_len", "crop_abs_delta", "crop_n_changed_units", "crop_frac_changed",
    "target_inside_crop",
]


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_path", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--limit", type=int, default=None, help="only the first N CSV rows (smoke test)")
    args = ap.parse_args()

    import torch
    from bulk_knockout_pigmentation import _build_context
    from genotype_cnn_alignment_deeplift_summary.annotations import gene_exon_haplotype_local_boxes
    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction, reorder_to_canonical

    device = torch.device("cpu")
    _log("Building context (no AlphaGenome calls will be made; cache only)...")
    ctx, ko_genes, _test_ids, _ped, _cn, _si, _wi = _build_context(device)

    df = pd.read_csv(args.in_path)
    if args.limit:
        df = df.head(args.limit)
    _log(f"{len(df)} knockout rows x 2 haplotypes = {2*len(df)} measurements -> {args.out}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if args.out.exists():
        prev = pd.read_csv(args.out)
        done = set(zip(prev["sample_id"], prev["gene"], prev["method"], prev["haplotype"]))
        _log(f"Resuming: {len(done)} measurements already present.")
    write_header = not args.out.exists()

    n_ok = n_missing = n_failed = 0
    t0 = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            w.writeheader(); fh.flush()

        for _, row in df.iterrows():
            sample_id, gene, method = row["sample_id"], row["gene"], row["method"]
            for haplotype, idx_col in (("H1", "h1_target_local_idx"), ("H2", "h2_target_local_idx")):
                if (sample_id, gene, method, haplotype) in done:
                    continue
                target_idx = int(row[idx_col])
                key = f"{sample_id}_{gene}_{haplotype}_{method}_{target_idx}_scramble{int(row['scramble_window'])}"
                npz = AG_CACHE / f"seq_{key}.npz"
                meta_json = AG_CACHE / f"seq_{key}_meta.json"
                if not (npz.exists() and meta_json.exists()):
                    n_missing += 1
                    continue
                try:
                    orig_array, orig_meta = load_raw_prediction(ctx.dataset_dir, sample_id, gene, haplotype)
                    canonical = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
                    mod_values = np.load(npz)["values"]
                    mod_meta = json.loads(meta_json.read_text())
                    mod_array = reorder_to_canonical(mod_values, pd.DataFrame(mod_meta), canonical)

                    # ---- (1) raw movement, haplotype-local ----
                    diff = np.abs(mod_array.astype(np.float64) - orig_array.astype(np.float64))
                    per_pos = diff.sum(axis=1)
                    changed = np.flatnonzero(per_pos > 1e-6)
                    boxes, _strand = gene_exon_haplotype_local_boxes(
                        ctx.dataset_dir, sample_id, gene, haplotype, ctx.gene_id_map,
                        ctx.transcript_extractor_mane, ctx.transcript_extractor_coding,
                    )
                    exon_mask = np.zeros(orig_array.shape[0], dtype=bool)
                    for a, b in boxes:
                        a = max(0, int(a)); b = min(orig_array.shape[0], int(b))
                        if b > a:
                            exon_mask[a:b] = True
                    sum_orig_ex = float(orig_array[exon_mask].sum()) if exon_mask.any() else float("nan")
                    sum_mod_ex = float(mod_array[exon_mask].sum()) if exon_mask.any() else float("nan")

                    # ---- (2) movement reaching the CNN, via the real code path ----
                    base = ctx.full_ds._process_window_haplotype_channels(
                        sample_id, gene, haplotype, {"rna_seq": orig_array}, {"rna_seq": orig_meta})
                    pert = ctx.full_ds._process_window_haplotype_channels(
                        sample_id, gene, haplotype, {"rna_seq": mod_array}, {"rna_seq": orig_meta})
                    if base is None or pert is None:
                        raise RuntimeError("processed window returned None")
                    sb, sp = base[0], pert[0]
                    cdiff = np.abs(np.asarray(sp, dtype=np.float64) - np.asarray(sb, dtype=np.float64))
                    cper = cdiff.sum(axis=0) if cdiff.ndim > 1 else cdiff
                    cchanged = int((cper > 1e-6).sum())
                    crop_len = int(cper.shape[0])

                    w.writerow({
                        "sample_id": sample_id, "gene": gene, "method": method,
                        "haplotype": haplotype, "target_local_idx": target_idx,
                        "raw_abs_delta_total": float(diff.sum()),
                        "raw_abs_delta_exons": float(diff[exon_mask].sum()) if exon_mask.any() else float("nan"),
                        "raw_sum_orig_exons": sum_orig_ex,
                        "raw_sum_mod_exons": sum_mod_ex,
                        "raw_rel_delta_exons": (sum_mod_ex - sum_orig_ex) / sum_orig_ex if sum_orig_ex else float("nan"),
                        "raw_n_changed_positions": int(changed.size),
                        "raw_first_changed": int(changed[0]) if changed.size else -1,
                        "raw_last_changed": int(changed[-1]) if changed.size else -1,
                        "crop_len": crop_len,
                        "crop_abs_delta": float(cdiff.sum()),
                        "crop_n_changed_units": cchanged,
                        "crop_frac_changed": cchanged / crop_len if crop_len else float("nan"),
                        "target_inside_crop": bool(cchanged > 0),
                    })
                    fh.flush()
                    n_ok += 1
                    if n_ok % 100 == 0:
                        el = time.monotonic() - t0
                        _log(f"progress: {n_ok} ok, {n_missing} cache-miss, {n_failed} failed "
                             f"({n_ok/el:.2f}/s, last={sample_id}/{gene}/{method}/{haplotype})")
                except Exception as exc:  # noqa: BLE001
                    n_failed += 1
                    _log(f"FAILED {sample_id}/{gene}/{method}/{haplotype}: {exc!r}")

    _log(f"Done. {n_ok} measured, {n_missing} cache-miss, {n_failed} failed. Output: {args.out}")
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
