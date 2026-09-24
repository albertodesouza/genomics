#!/usr/bin/env python3
"""Knock a gene down on a classifier that was trained on THAT GENE ALONE.

In the eleven-gene panel the readout measures RELIANCE: how much this classifier's decision moves
when a gene's promoter is scrambled. Reliance is a property of the trained model, not of the gene:
a gene can be flat because its window carries nothing about the label, or because it carries
something a cheaper gene carries too and the classifier routed around it. The panel measurement
cannot tell those apart.

A single-gene classifier removes the choice. With one gene in the input there is nothing to route
to, so:
    single-gene test accuracy   = the INFORMATION the window carries about the label
    single-gene |Delta|         = the same 100 bp promoter scramble, with reliance at its maximum
    panel |Delta| / single |Delta| = how much of the available response the panel model gave up

Two arms, chosen for opposite efficiency in the published panel (readout_normalisation.json):
    SLC24A5  |Delta| 0.503 at |Delta_in|  6,537   ratio 7.7e-5  (most efficient in the panel)
    TYR      |Delta| 0.337 at |Delta_in| 54,297   ratio 6.2e-6  (12x less efficient)

Costs nothing at the AlphaGenome API: every scrambled prediction is read from the on-disk cache
written by the published bulk run. All three anchoring methods are replayed, since they are all
already cached. Refuses to run rather than emit numbers if a cache entry is missing.

Correctness guard: rebuilds one individual's UNPERTURBED tensor by hand and asserts it matches
what the dataset itself produces, before replaying anything.

Usage:
  python3 scripts/experiments/single_gene_knockdown_replay.py --arm slc24a5
  ... --limit 10     # smoke test
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
for _p in (REPO_ROOT / "src", REPO_ROOT / "notebooks"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd
import torch

AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
REFERENCE_KO = KO_DIR / "pigmentation_test_split_knockout_no_alignment.csv"
# The three specificity control genes were introduced by the 14-gene run and never appear in
# the published panel CSV. That run anchored on biology_tss only -- the paper's primary
# strategy and the only one whose scrambles were ever generated for them -- so a control arm
# replays one method where a panel arm replays three. Verified before use: on the 1782 rows
# the two files share, h1/h2_target_local_idx and scramble_window agree exactly, so the cache
# keys built from either file are the same and no arm re-bills AlphaGenome.
CONTROL_KO = KO_DIR / "pigmentation_test_split_knockout_14gene.csv"
# Third fallback, for a gene built after the last knockdown run: it appears in neither CSV,
# so scramble_prefetch.py records the loci it derived when it populated the cache. Same
# deterministic haplotype_local_idx call, so the keys rebuilt from it are the same keys.
MANIFEST_KO = KO_DIR / "scramble_manifest.csv"
# The first two arms were written by hand; the sweep generates the rest under single_gene/.
LEGACY_ARMS = {
    "slc24a5": "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment_single_slc24a5.yaml",
    "tyr": "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment_single_tyr.yaml",
}
SWEEP_CONFIG_DIR = Path("configs/predictors/genotype_based/pigmentation/single_gene")


def resolve_arm(arm: str):
    """arm -> (gene, config path). Accepts the gene name in any case."""
    key = arm.lower()
    if key in LEGACY_ARMS:
        return arm.upper(), LEGACY_ARMS[key]
    cand = SWEEP_CONFIG_DIR / f"pigmentation_binary_single_{key}.yaml"
    if (REPO_ROOT / cand).exists():
        return arm.upper(), str(cand)
    raise SystemExit(f"ABORT: no config for arm {arm!r}. Looked at "
                     f"{LEGACY_ARMS.get(key, '(no legacy entry)')} and {cand}. "
                     "Generate it with scripts/experiments/single_gene_sweep.py --write-configs.")
FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "arm", "method",
    "scramble_window", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit", "perturbed_strong_logit", "perturbed_weak_logit",
    "delta_log_odds", "baseline_pred", "perturbed_pred", "flipped", "correct_baseline", "correct_perturbed",
    "h1_crop_abs_delta", "h2_crop_abs_delta", "delta_in",
    "h1_expr_abs_delta", "h2_expr_abs_delta", "delta_expr", "delta_expr_rel", "expr_baseline",
    "delta_in_norm", "n_genes_in_model", "expr_mod", "delta_expr_signed", "expr_log2fc",
]


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", required=True, help="gene name, e.g. SLC24A5 / tpm2")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--config", type=Path, default=None,
                    help="explicit config of an already-trained isolated run (e.g. the "
                         "config.yaml a poolmax run wrote next to its checkpoint). The run "
                         "directory it resolves to must still be gene-isolated.")
    ap.add_argument("--allow-multigene", action="store_true",
                    help="permit a config whose gene list contains the arm plus others. The "
                         "arm is then the gene whose promoter is scrambled while the remaining "
                         "genes stay at baseline -- the panel-style measurement, which reports "
                         "RELIANCE (how much this model leans on the gene) rather than the "
                         "information the window carries. The model's gene count is recorded "
                         "per row as n_genes_in_model so the provenance travels with the data.")
    ap.add_argument("--no-expr", action="store_true",
                    help="skip the exon-integrated expression delta (delta_expr). It needs the "
                         "local GENCODE cache but no API call; on by default because a knockdown "
                         "table cannot separate 'the CNN ignored it' from 'expression never "
                         "moved' without it.")
    ap.add_argument("--method", default=None,
                    help="replay only this anchoring method. Needed to compare a panel gene "
                         "against a control: the panel CSV carries three methods and the "
                         "control sources one or two, so an unfiltered replay averages "
                         "different numbers of scrambles per gene.")
    args = ap.parse_args()
    gene, config_path = resolve_arm(args.arm)
    if args.config is not None:
        config_path = args.config
    arm = args.arm.lower()
    out = args.out or KO_DIR / f"pigmentation_test_split_knockout_single_{arm}.csv"

    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction, reorder_to_canonical
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import build_model, load_checkpoint
    from genomics.predictors.genotype_based.config import (
        generate_experiment_name, get_dataset_cache_dir, get_experiment_runs_dir, load_config)
    from genomics.predictors.genotype_based.data.pipeline import _make_runtime_processed_datasets, _resolve_runtime_dataset_dir

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"device={device}  arm={arm}  gene={gene}  config={config_path}")

    cfg = load_config(REPO_ROOT / config_path)
    genes_all = list(cfg.dataset_input.genes_to_use)
    if args.allow_multigene:
        if gene not in genes_all:
            raise SystemExit(f"ABORT: arm {gene!r} is not in the config's gene list {genes_all}")
    elif genes_all != [gene]:
        raise SystemExit(f"ABORT: config gene list is {genes_all}, expected [{gene!r}]. "
                         f"Pass --allow-multigene for the panel-style measurement.")
    # Checked BEFORE the dataset is materialised: building a 1.6 GB tensor view only to discover
    # the arm was never trained wastes minutes and disk on a machine that has little of either.
    exp = get_experiment_runs_dir(cfg) / generate_experiment_name(cfg)
    # The guard exists so a replay can never silently read a multi-gene panel checkpoint and
    # report it as a single-gene number. What matters is isolation, not the directory's name:
    # require the arm to appear in the path, and the config to name exactly this one gene
    # (already asserted above). Without --config the historical layout is still required.
    if args.config is None:
        if f"runs_single_gene_{arm}" not in str(exp):
            raise SystemExit(f"ABORT: checkpoint dir is not the isolated single-gene run: {exp}")
    elif args.allow_multigene:
        pass  # a multi-gene run directory is not named after any one arm
    elif arm not in str(exp).lower():
        raise SystemExit(f"ABORT: --config resolves to {exp}, which does not name arm {arm!r}; "
                         f"refusing to replay a checkpoint that may not be this gene's.")
    ck = exp / "models" / "best_accuracy.pt"
    if not ck.exists():
        raise SystemExit(f"ABORT: no checkpoint at {ck}. Train this arm first: "
                         f"python3 scripts/experiments/single_gene_sweep.py --genes {gene}")

    cache_dir = get_dataset_cache_dir(cfg)
    ds_dir = _resolve_runtime_dataset_dir(cfg)
    with quiet_pipeline_logs():
        full_ds, _tr, _va, _te = _make_runtime_processed_datasets(ds_dir, cache_dir, cfg)
    _log(f"checkpoint {ck}")
    model = build_model(cfg, full_ds, device)
    model = load_checkpoint(model, ck, device)
    model.eval()

    class_names = full_ds.get_class_names()
    strong_idx = [i for i, n in full_ds.idx_to_target.items() if n == "strong pigmentation"][0]
    weak_idx = next(i for i in range(2) if i != strong_idx)
    ped = full_ds.dataset_metadata.get("individuals_pedigree", {})
    dataset_dir = Path(cfg.dataset_input.dataset_dir)

    def crop_rows(array, meta):
        r = full_ds._process_haplotype_raw_center_crop({"rna_seq": array}, {"rna_seq": meta})
        if r is None:
            raise RuntimeError("raw_center_crop returned None")
        return r

    def dita_signal(sid, g, hap, arr, meta):
        """The signal rows the CNN actually consumes, via the dataset's own DITA path.

        Under tensor_layout='haplotype_channels' the tensor is (haplotype, channel, L) and the
        alignment to haplotype coordinates happens inside this call. The older
        `_process_haplotype_raw_center_crop` above produced the flat pre-DITA layout, which
        vstack-ed to (2, L) and no longer matches what the model was trained on.
        """
        r = full_ds._process_window_haplotype_channels(sid, g, hap, {"rna_seq": arr}, {"rna_seq": meta})
        if r is None:
            raise RuntimeError(f"haplotype_channels returned None for {sid}/{g}/{hap}")
        signals, masks = r
        if full_ds.feature_mode == "signals_only":
            return signals
        if full_ds.feature_mode == "masks_only":
            return masks
        return np.concatenate([signals, masks], axis=0)

    def build_tensor(sid, overrides=None):
        overrides = overrides or {}
        per_hap_rows = {"H1": [], "H2": []}
        for g in genes_all:
            for hap in ("H1", "H2"):
                arr, meta = overrides.get((g, hap)) or load_raw_prediction(dataset_dir, sid, g, hap)
                per_hap_rows[hap].append(dita_signal(sid, g, hap, arr, meta))
        stacked = np.stack([np.concatenate(per_hap_rows["H1"], axis=0),
                            np.concatenate(per_hap_rows["H2"], axis=0)])
        return full_ds._normalize_features_tensor(torch.FloatTensor(stacked.astype(np.float32)))

    src = REFERENCE_KO
    df = pd.read_csv(src)
    df = df[df["gene"] == gene].copy()
    if df.empty:
        src = CONTROL_KO
        df = pd.read_csv(src)
        df = df[df["gene"] == gene].copy()
    if df.empty and MANIFEST_KO.exists():
        src = MANIFEST_KO
        df = pd.read_csv(src)
        df = df[df["gene"] == gene].copy()
    if df.empty:
        raise SystemExit(
            f"ABORT: no scramble rows for {gene!r} in {REFERENCE_KO.name}, {CONTROL_KO.name} "
            f"or {MANIFEST_KO.name}. The scrambles must exist before an arm can be replayed.")
    if args.method:
        before = len(df)
        df = df[df["method"] == args.method].copy()
        if df.empty:
            raise SystemExit(f"ABORT: {gene} has no rows for method {args.method!r} in "
                             f"{src.name} (had {before} rows across "
                             f"{sorted(pd.read_csv(src).query('gene == @gene')['method'].unique())})")
    if args.limit:
        df = df.head(args.limit)
    _log(f"{len(df)} reference rows from {src.name} "
         f"({df['method'].nunique()} methods, {df['sample_id'].nunique()} individuals)")

    # ---- exon-integrated expression delta ----
    # delta_expr answers a different question from delta_in: delta_in is how much perturbation
    # entered the CNN's input, delta_expr is how much AlphaGenome's predicted expression of the
    # gene actually moved. A gene can have a large delta_in and a flat delta_expr (the scramble
    # landed in the window but not on anything transcribed), and the paper's 0.00 cells need
    # the two separated. Integrated over exons, since the expression of a gene is its transcript
    # signal, not the whole 32 kb window. Read from the local GENCODE cache: no API call.
    exon_ctx = None
    if not args.no_expr:
        from genotype_cnn_alignment_deeplift_summary.annotations import (
            build_transcript_extractors, gene_exon_haplotype_local_boxes, load_gtf)
        _log("loading GENCODE annotations for delta_expr ...")
        _gtf = load_gtf(REPO_ROOT / "notebooks" / ".cache" / "annotations")
        _mane, _pc, _ex_mane, _ex_cod, _gid = build_transcript_extractors(_gtf)
        exon_ctx = (_gid, _ex_mane, _ex_cod, gene_exon_haplotype_local_boxes)
        _log("annotations ready")

    _exon_mask_cache = {}

    def exon_mask_for(sid, hap, n_pos):
        """Boolean mask over this haplotype's local positions that are exonic for this gene."""
        ck = (sid, hap, n_pos)
        if ck not in _exon_mask_cache:
            gid, ex_mane, ex_cod, boxes_fn = exon_ctx
            boxes, _strand = boxes_fn(dataset_dir, sid, gene, hap, gid, ex_mane, ex_cod)
            m = np.zeros(n_pos, dtype=bool)
            for a, b in boxes:
                a, b = max(0, int(a)), min(n_pos, int(b))
                if b > a:
                    m[a:b] = True
            if not m.any():
                raise RuntimeError(f"no exonic positions for {gene}/{sid}/{hap}; delta_expr "
                                   f"would be undefined rather than zero")
            _exon_mask_cache[ck] = m
        return _exon_mask_cache[ck]

    # ---- correctness guard ----
    sid_list = [full_ds._sample_id_for_base_index(b) for b in full_ds.valid_sample_indices]
    probe = df.iloc[0]["sample_id"]
    if probe not in sid_list:
        raise SystemExit(f"ABORT: {probe} not in this dataset index.")
    mine = build_tensor(probe).squeeze()
    theirs = full_ds[sid_list.index(probe)][0].squeeze()
    if mine.shape != theirs.shape:
        raise SystemExit(f"ABORT: shape {tuple(mine.shape)} vs {tuple(theirs.shape)}")
    md = float((mine - theirs).abs().max())
    if md > 1e-4:
        raise SystemExit(f"ABORT: reconstructed baseline differs (max|diff|={md:.3e})")
    _log(f"guard OK: max|diff|={md:.2e}, tensor {tuple(mine.shape)}")

    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        prev = pd.read_csv(out)
        done = set(zip(prev["sample_id"], prev["method"]))
    write_header = not out.exists()

    baseline_cache = {}
    n_ok = n_missing = n_fail = 0
    t0 = time.monotonic()
    with open(out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            w.writeheader(); fh.flush()
        for _, row in df.iterrows():
            sid, method = row["sample_id"], row["method"]
            if (sid, method) in done:
                continue
            try:
                overrides, per_hap = {}, {}
                per_hap_expr, per_hap_base, per_hap_mod = {}, {}, {}
                miss = False
                for hap, col in (("H1", "h1_target_local_idx"), ("H2", "h2_target_local_idx")):
                    tgt = int(row[col])
                    key = f"{sid}_{gene}_{hap}_{method}_{tgt}_scramble{int(row['scramble_window'])}"
                    npz, mj = AG_CACHE / f"seq_{key}.npz", AG_CACHE / f"seq_{key}_meta.json"
                    if not (npz.exists() and mj.exists()):
                        miss = True
                        break
                    orig_arr, orig_meta = load_raw_prediction(dataset_dir, sid, gene, hap)
                    canonical = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
                    mod = reorder_to_canonical(np.load(npz)["values"],
                                               pd.DataFrame(json.loads(mj.read_text())), canonical)
                    overrides[(gene, hap)] = (mod, orig_meta)
                    base_r = dita_signal(sid, gene, hap, orig_arr, orig_meta)
                    pert_r = dita_signal(sid, gene, hap, mod, orig_meta)
                    per_hap[hap] = float(np.abs(np.asarray(pert_r, float) - np.asarray(base_r, float)).sum())
                    if exon_ctx is not None:
                        o64 = np.asarray(orig_arr, dtype=np.float64)
                        m64 = np.asarray(mod, dtype=np.float64)
                        em = exon_mask_for(sid, hap, o64.shape[0])
                        # delta_expr is a MAGNITUDE (sum of |change|); on its own it cannot
                        # say whether the scramble lowered predicted expression, which is the
                        # only direction that makes the intervention a knockdown. The signed
                        # exon sums below settle that.
                        per_hap_expr[hap] = float(np.abs(m64[em] - o64[em]).sum())
                        per_hap_base[hap] = float(o64[em].sum())
                        per_hap_mod[hap] = float(m64[em].sum())
                if miss:
                    n_missing += 1
                    continue
                # delta_in above is in raw AlphaGenome units, which is the published
                # convention but is NOT comparable across genes: normalization divides by
                # this gene's own log_max, and for a near-silent melanocyte window that
                # divisor is tiny (PRSS55 0.0014, LRRC36 0.0151). A raw delta_in of 0.9 in
                # PRSS55 reaches the model as ~612. delta_in_norm is the same perturbation
                # measured on the normalized tensor the model actually consumes, so it is
                # the only one of the two that can be compared between two genes.
                if sid not in baseline_cache:
                    base_t = build_tensor(sid)
                    with torch.no_grad():
                        baseline_cache[sid] = (
                            model(base_t.unsqueeze(0).float().to(device))[0].cpu().numpy(), base_t)
                bl, base_t = baseline_cache[sid]
                pert_t = build_tensor(sid, overrides=overrides)
                with torch.no_grad():
                    pl = model(pert_t.unsqueeze(0).float().to(device))[0].cpu().numpy()
                delta_in_norm = float((pert_t - base_t).abs().sum())
                pr = ped.get(sid, {})
                true_label = full_ds._get_target_value(pr)
                bpred, ppred = class_names[int(bl.argmax())], class_names[int(pl.argmax())]
                w.writerow({
                    "sample_id": sid, "population": pr.get("population"),
                    "superpopulation": pr.get("superpopulation"), "true_label": true_label,
                    "gene": gene, "arm": arm, "method": method,
                    "scramble_window": int(row["scramble_window"]),
                    "h1_target_local_idx": int(row["h1_target_local_idx"]),
                    "h2_target_local_idx": int(row["h2_target_local_idx"]),
                    "baseline_strong_logit": float(bl[strong_idx]), "baseline_weak_logit": float(bl[weak_idx]),
                    "perturbed_strong_logit": float(pl[strong_idx]), "perturbed_weak_logit": float(pl[weak_idx]),
                    "delta_log_odds": float((pl[strong_idx] - pl[weak_idx]) - (bl[strong_idx] - bl[weak_idx])),
                    "baseline_pred": bpred, "perturbed_pred": ppred,
                    "flipped": int(bpred != ppred),
                    "correct_baseline": int(bpred == true_label), "correct_perturbed": int(ppred == true_label),
                    "h1_crop_abs_delta": per_hap["H1"], "h2_crop_abs_delta": per_hap["H2"],
                    "delta_in": per_hap["H1"] + per_hap["H2"],
                    "delta_in_norm": delta_in_norm,
                    "n_genes_in_model": len(genes_all),
                    "h1_expr_abs_delta": per_hap_expr.get("H1", float("nan")),
                    "h2_expr_abs_delta": per_hap_expr.get("H2", float("nan")),
                    "delta_expr": (per_hap_expr["H1"] + per_hap_expr["H2"]
                                   if per_hap_expr else float("nan")),
                    "expr_baseline": (per_hap_base["H1"] + per_hap_base["H2"]
                                      if per_hap_base else float("nan")),
                    "delta_expr_rel": ((per_hap_expr["H1"] + per_hap_expr["H2"])
                                       / (per_hap_base["H1"] + per_hap_base["H2"])
                                       if per_hap_base and (per_hap_base["H1"] + per_hap_base["H2"]) > 0
                                       else float("nan")),
                    "expr_mod": (per_hap_mod["H1"] + per_hap_mod["H2"]
                                 if per_hap_mod else float("nan")),
                    "delta_expr_signed": ((per_hap_mod["H1"] + per_hap_mod["H2"])
                                          - (per_hap_base["H1"] + per_hap_base["H2"])
                                          if per_hap_mod else float("nan")),
                    "expr_log2fc": (float(np.log2((per_hap_mod["H1"] + per_hap_mod["H2"])
                                                  / (per_hap_base["H1"] + per_hap_base["H2"])))
                                    if per_hap_mod and (per_hap_base["H1"] + per_hap_base["H2"]) > 0
                                    and (per_hap_mod["H1"] + per_hap_mod["H2"]) > 0
                                    else float("nan")),
                })
                fh.flush()
                n_ok += 1
                if n_ok % 100 == 0:
                    _log(f"{n_ok} rows ({n_ok/(time.monotonic()-t0):.1f}/s), missing={n_missing}, failed={n_fail}")
            except Exception as exc:  # noqa: BLE001
                n_fail += 1
                _log(f"FAILED {sid}/{method}: {exc!r}")
    _log(f"DONE: ok={n_ok} missing={n_missing} failed={n_fail} -> {out}")
    if n_missing:
        _log("WARNING: rows were skipped for a missing cache entry; the replay is incomplete.")
    return 0 if (n_fail == 0 and n_missing == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
