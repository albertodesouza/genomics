#!/usr/bin/env python3
"""E16 -- cross-task transfer: does a SUPERPOPULATION-trained model rely on the same genes?

The pigmentation label is confounded with continental ancestry by construction, so a gene ranking
recovered from the pigmentation classifier admits two readings: the classifier uses these loci
because they carry pigmentation-relevant predicted function, or because they are the loci whose
predicted function best separates African from European haplotype backgrounds. A perturbation that
is not a variant removes the confound from the INTERVENTION, but not from the FUNCTION BEING
PROBED -- that was fixed at training time by the label.

This replays the already-computed promoter scrambles through the superpopulation-trained
tracks-only classifier, which was never told anything about pigmentation. Same DITA tensor layout,
same genes, same crop, same cached AlphaGenome predictions; only the model and the label set
differ (5 classes instead of 2).

Reading the outcome (either is publishable under a method framing; both are informative):
  rankings AGREE   -> the readout tracks ancestry-discriminative predicted function, and the
                      pigmentation ranking cannot be read as pigmentation-specific.
  rankings DIFFER  -> the readout is label-specific, i.e. it measures what THIS classifier uses
                      rather than a fixed property of the representation. Materially strengthens
                      the claim that the probe measures model reliance.

Costs nothing at the AlphaGenome API: every scrambled prediction is read from the on-disk cache
written by the original bulk run. Only the tensor build + CNN forward pass are redone.

Correctness guard: before replaying anything the script rebuilds one individual's UNPERTURBED
tensor by hand and asserts it matches what the dataset itself produces. It refuses to run rather
than emit plausible garbage.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/crosstask_knockout_replay.py
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

CONFIG = "configs/predictors/genotype_based/superpopulation/superpopulation.yaml"
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
DEFAULT_IN = KO_DIR / "pigmentation_test_split_knockout.csv"
DEFAULT_OUT = KO_DIR / "pigmentation_test_split_knockout_superpopulation.csv"
AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"

FIELDS = [
    "sample_id", "population", "superpopulation", "gene", "method",
    "scramble_window", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_pred", "perturbed_pred", "flipped",
    "baseline_probs", "perturbed_probs",
    "delta_baseline_class", "tv_distance", "in_superpop_train",
]


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


def softmax(x: np.ndarray) -> np.ndarray:
    e = np.exp(x - x.max())
    return e / e.sum()


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
    _log(f"Loading superpopulation DITA dataset ({ds_dir}, cache={cache_dir})...")
    with quiet_pipeline_logs():
        full_ds, train_ds, val_ds, test_ds = _make_runtime_processed_datasets(ds_dir, cache_dir, cfg)
    exp = get_experiment_runs_dir(cfg) / generate_experiment_name(cfg)
    ck = exp / "models" / "best_accuracy.pt"
    if not ck.exists():
        raise SystemExit(f"ABORT: no checkpoint at {ck}")
    _log(f"Loading checkpoint {ck}...")
    model = build_model(cfg, full_ds, device)
    model = load_checkpoint(model, ck, device)
    model.eval()
    return cfg, full_ds, train_ds, model, full_ds.get_class_names()


def build_tensor(cfg, full_ds, dataset_dir, sample_id, overrides=None):
    """haplotype_channels (DITA) layout: H1 rows and H2 rows as two channels."""
    from genomics.predictors.genotype_based.data.normalization import apply_normalization
    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction
    overrides = overrides or {}
    h1_rows, h2_rows = [], []
    for gene in cfg.dataset_input.genes_to_use:
        for hap, rows in (("H1", h1_rows), ("H2", h2_rows)):
            if (gene, hap) in overrides:
                array, meta = overrides[(gene, hap)]
            else:
                array, meta = load_raw_prediction(dataset_dir, sample_id, gene, hap)
            result = full_ds._process_window_haplotype_channels(
                sample_id, gene, hap, {"rna_seq": array}, {"rna_seq": meta})
            if result is None:
                raise RuntimeError(f"no rows for {sample_id}/{gene}/{hap}")
            signals, _masks = result
            rows.append(signals)
    stacked = np.stack([np.concatenate(h1_rows, axis=0), np.concatenate(h2_rows, axis=0)], axis=0)
    return apply_normalization(torch.from_numpy(stacked), full_ds.normalization_params)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_path", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction, reorder_to_canonical

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"device={device}")
    cfg, full_ds, train_ds, model, class_names = build_ctx(device)
    dataset_dir = Path(cfg.dataset_input.dataset_dir)
    _log(f"classes = {class_names}")

    # Which of the pigmentation test individuals the superpopulation model was TRAINED on.
    # Recorded per row so the replay can be read on the held-out subset alone if needed.
    try:
        if hasattr(train_ds, "valid_sample_indices"):
            base = list(train_ds.valid_sample_indices)
        else:  # torch Subset over full_ds: .indices are positions into valid_sample_indices
            base = [full_ds.valid_sample_indices[i] for i in train_ds.indices]
        train_ids = {full_ds._sample_id_for_base_index(b) for b in base}
    except Exception as exc:  # pragma: no cover - diagnostic only
        _log(f"WARNING: could not resolve superpopulation train ids ({exc}); marking unknown.")
        train_ids = None

    df = pd.read_csv(args.in_path)
    if args.limit:
        df = df.head(args.limit)
    if train_ids is not None:
        overlap = sorted(set(df["sample_id"]) & train_ids)
        _log(f"Pigmentation test individuals also in the superpopulation TRAIN split: "
             f"{len(overlap)}/{df['sample_id'].nunique()}")

    # ---- correctness guard ----
    probe_id = df.iloc[0]["sample_id"]
    sid_list = [full_ds._sample_id_for_base_index(b) for b in full_ds.valid_sample_indices]
    if probe_id not in sid_list:
        _log(f"WARNING: {probe_id} not in this dataset's index; skipping the equality guard.")
    else:
        mine = build_tensor(cfg, full_ds, dataset_dir, probe_id).squeeze()
        theirs, _t = full_ds[sid_list.index(probe_id)]
        theirs = theirs.squeeze()
        if mine.shape != theirs.shape:
            _log(f"ABORT: shape mismatch, mine={tuple(mine.shape)} dataset={tuple(theirs.shape)}")
            return 2
        md = float((mine - theirs).abs().max())
        if md > 1e-4:
            _log(f"ABORT: reconstructed baseline differs from the dataset's own (max|diff|={md:.3e}).")
            return 2
        _log(f"Guard OK: max|diff|={md:.2e}, shape={tuple(mine.shape)}")

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
                bp_i, pp_i = int(bl.argmax()), int(pl.argmax())
                b_prob, p_prob = softmax(bl), softmax(pl)
                pr = ped.get(sid, {})
                w.writerow({
                    "sample_id": sid, "population": pr.get("population"),
                    "superpopulation": pr.get("superpopulation"),
                    "gene": gene, "method": method,
                    "scramble_window": int(row["scramble_window"]),
                    "h1_target_local_idx": int(row["h1_target_local_idx"]),
                    "h2_target_local_idx": int(row["h2_target_local_idx"]),
                    "baseline_pred": class_names[bp_i], "perturbed_pred": class_names[pp_i],
                    "flipped": int(bp_i != pp_i),
                    "baseline_probs": ";".join(f"{v:.6f}" for v in b_prob),
                    "perturbed_probs": ";".join(f"{v:.6f}" for v in p_prob),
                    # signed drop in confidence in the model's own baseline call
                    "delta_baseline_class": float(p_prob[bp_i] - b_prob[bp_i]),
                    # class-agnostic magnitude of the induced shift
                    "tv_distance": float(0.5 * np.abs(p_prob - b_prob).sum()),
                    "in_superpop_train": (int(sid in train_ids) if train_ids is not None else ""),
                })
                n_ok += 1
                if n_ok % 200 == 0:
                    fh.flush()
                    _log(f"{n_ok} replayed ({time.monotonic()-t0:.0f}s), missing={n_missing}, failed={n_failed}")
            except Exception as exc:
                n_failed += 1
                _log(f"FAILED {sid}/{gene}/{method}: {exc}")
    _log(f"DONE: ok={n_ok} missing={n_missing} failed={n_failed} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
