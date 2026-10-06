#!/usr/bin/env python3
"""E18.2 -- does the UNALIGNED model respond to the same knockdown as the aligned one?

The paper argues that DITA's value is interpretive: a shared axis is a precondition for reading
attribution as a genomic statement. The empirical half of that claim is whether the aligned model
is actually the better probe. This replays the *already-computed* promoter scrambles through the
`no_alignment` (raw_center_crop) classifier and asks whether it recovers the same genes.

Costs nothing at the AlphaGenome API: every scrambled prediction is read from the on-disk cache
written by the original bulk run. Only the downstream tensor build + CNN forward pass are redone,
against a different model and a different tensor layout.

Correctness guard: before replaying anything, the script rebuilds one individual's UNPERTURBED
tensor by hand and asserts it matches what the dataset itself produces for that individual. If the
row order or normalisation of the raw_center_crop layout were wrong, every downstream number would
be quietly wrong, so this refuses to run rather than emit plausible garbage.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/unaligned_knockout_replay.py
  ... --limit 20    # smoke test
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
import torch

CONFIG = "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment.yaml"
DEFAULT_IN = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout.csv"
DEFAULT_OUT = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout_no_alignment.csv"
AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"

FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "method",
    "scramble_window", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit", "perturbed_strong_logit", "perturbed_weak_logit",
    "baseline_pred", "perturbed_pred", "flipped", "model",
]


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


def build_ctx(device):
    from dotenv import load_dotenv
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import build_model, load_checkpoint
    from genomics.predictors.genotype_based.config import (
        generate_experiment_name, get_dataset_cache_dir, get_experiment_runs_dir, load_config)
    from genomics.predictors.genotype_based.data.pipeline import (
        _make_runtime_processed_datasets, _resolve_runtime_dataset_dir)

    load_dotenv(Path.home() / ".env")
    cfg = load_config(REPO_ROOT / CONFIG)
    cache_dir = get_dataset_cache_dir(cfg)
    ds_dir = _resolve_runtime_dataset_dir(cfg)
    _log(f"Loading no_alignment dataset ({ds_dir}, cache={cache_dir})...")
    with quiet_pipeline_logs():
        full_ds, _a, _b, _c = _make_runtime_processed_datasets(ds_dir, cache_dir, cfg)
    exp = get_experiment_runs_dir(cfg) / generate_experiment_name(cfg)
    ck = exp / "models" / "best_accuracy.pt"
    _log(f"Loading checkpoint {ck}...")
    model = build_model(cfg, full_ds, device)
    model = load_checkpoint(model, ck, device)
    class_names = full_ds.get_class_names()
    idx_to_target = full_ds.idx_to_target
    strong = [i for i, n in idx_to_target.items() if n == "strong pigmentation"][0]
    weak = next(i for i in range(2) if i != strong)
    return cfg, full_ds, model, class_names, strong, weak


def build_tensor(cfg, full_ds, dataset_dir, sample_id, overrides=None):
    """raw_center_crop layout: for each gene, H1 rows then H2 rows, vstacked (no shared axis)."""
    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction
    overrides = overrides or {}
    rows = []
    for gene in cfg.dataset_input.genes_to_use:
        for hap in ("H1", "H2"):
            if (gene, hap) in overrides:
                array, meta = overrides[(gene, hap)]
            else:
                array, meta = load_raw_prediction(dataset_dir, sample_id, gene, hap)
            r = full_ds._process_haplotype_raw_center_crop({"rna_seq": array}, {"rna_seq": meta})
            if r is None:
                raise RuntimeError(f"no rows for {sample_id}/{gene}/{hap}")
            rows.append(r)
    stacked = np.vstack(rows).astype(np.float32)
    # Use the dataset's own normaliser, not the bare apply_normalization: for this layout it is
    # what __getitem__ calls, and the equality guard in main() fails loudly if that ever diverges.
    return full_ds._normalize_features_tensor(torch.FloatTensor(stacked))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_path", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction, reorder_to_canonical

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"device={device}")
    cfg, full_ds, model, class_names, strong_idx, weak_idx = build_ctx(device)
    dataset_dir = Path(cfg.dataset_input.dataset_dir)

    df = pd.read_csv(args.in_path)
    if args.limit:
        df = df.head(args.limit)

    # ---- correctness guard ----
    probe_id = df.iloc[0]["sample_id"]
    mine = build_tensor(cfg, full_ds, dataset_dir, probe_id)
    sid_list = [full_ds._sample_id_for_base_index(b) for b in full_ds.valid_sample_indices]
    if probe_id not in sid_list:
        _log(f"WARNING: {probe_id} not in this dataset's index; skipping the equality guard.")
    else:
        theirs, _t = full_ds[sid_list.index(probe_id)]
        theirs = theirs.squeeze()
        mine_s = mine.squeeze()
        if mine_s.shape != theirs.shape:
            _log(f"ABORT: shape mismatch, mine={tuple(mine_s.shape)} dataset={tuple(theirs.shape)}")
            return 2
        md = float((mine_s - theirs).abs().max())
        if md > 1e-4:
            _log(f"ABORT: reconstructed baseline tensor differs from the dataset's own (max|diff|={md:.3e}).")
            return 2
        _log(f"Guard OK: reconstructed tensor matches dataset output for {probe_id} (max|diff|={md:.2e}, shape={tuple(mine_s.shape)}).")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if args.out.exists():
        prev = pd.read_csv(args.out)
        done = set(zip(prev["sample_id"], prev["gene"], prev["method"]))
        _log(f"Resuming: {len(done)} rows present.")
    write_header = not args.out.exists()

    ped = full_ds.dataset_metadata.get("individuals_pedigree", {})
    baseline_cache = {}
    n_ok = n_missing = n_failed = 0
    t0 = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            w.writeheader(); fh.flush()
        for _, row in df.iterrows():
            sid, gene, method = row["sample_id"], row["gene"], row["method"]
            if (sid, gene, method) in done:
                continue
            try:
                overrides = {}
                miss = False
                for hap, col in (("H1", "h1_target_local_idx"), ("H2", "h2_target_local_idx")):
                    tgt = int(row[col])
                    key = f"{sid}_{gene}_{hap}_{method}_{tgt}_scramble{int(row['scramble_window'])}"
                    npz, mj = AG_CACHE / f"seq_{key}.npz", AG_CACHE / f"seq_{key}_meta.json"
                    if not (npz.exists() and mj.exists()):
                        miss = True; break
                    orig_array, orig_meta = load_raw_prediction(dataset_dir, sid, gene, hap)
                    canonical = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
                    mod = reorder_to_canonical(np.load(npz)["values"],
                                               pd.DataFrame(json.loads(mj.read_text())), canonical)
                    overrides[(gene, hap)] = (mod, orig_meta)
                if miss:
                    n_missing += 1
                    continue
                if sid not in baseline_cache:
                    with torch.no_grad():
                        bt = build_tensor(cfg, full_ds, dataset_dir, sid)
                        baseline_cache[sid] = model(bt.unsqueeze(0).float().to(device))[0].cpu().numpy()
                bl = baseline_cache[sid]
                with torch.no_grad():
                    pt = build_tensor(cfg, full_ds, dataset_dir, sid, overrides=overrides)
                    pl = model(pt.unsqueeze(0).float().to(device))[0].cpu().numpy()
                bp, pp = class_names[int(bl.argmax())], class_names[int(pl.argmax())]
                pr = ped.get(sid, {})
                w.writerow({
                    "sample_id": sid, "population": pr.get("population"),
                    "superpopulation": pr.get("superpopulation"),
                    "true_label": full_ds._get_target_value(pr),
                    "gene": gene, "method": method,
                    "scramble_window": int(row["scramble_window"]),
                    "h1_target_local_idx": int(row["h1_target_local_idx"]),
                    "h2_target_local_idx": int(row["h2_target_local_idx"]),
                    "baseline_strong_logit": float(bl[strong_idx]),
                    "baseline_weak_logit": float(bl[weak_idx]),
                    "perturbed_strong_logit": float(pl[strong_idx]),
                    "perturbed_weak_logit": float(pl[weak_idx]),
                    "baseline_pred": bp, "perturbed_pred": pp, "flipped": bp != pp,
                    "model": "no_alignment",
                })
                fh.flush()
                n_ok += 1
                if n_ok % 100 == 0:
                    el = time.monotonic() - t0
                    _log(f"progress: {n_ok} ok, {n_missing} cache-miss, {n_failed} failed "
                         f"({n_ok/el:.2f}/s, last={sid}/{gene}/{method})")
            except Exception as exc:  # noqa: BLE001
                n_failed += 1
                _log(f"FAILED {sid}/{gene}/{method}: {exc!r}")
    _log(f"Done. {n_ok} replayed, {n_missing} cache-miss, {n_failed} failed. Output: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
