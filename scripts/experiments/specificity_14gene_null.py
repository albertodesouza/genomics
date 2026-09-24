#!/usr/bin/env python3
"""Null band for the fourteen-gene ranking: the same scramble, at a location that is NOT a promoter.

specificity_14gene_knockdown.py measures |Delta| for eleven pigmentation genes and three
pigmentation-irrelevant controls on the fourteen-gene checkpoint. Reading "the controls are flat"
requires a noise floor measured on THAT checkpoint: a 100 bp scramble somewhere in the same window
that is not the promoter. Identical in every respect -- same window size, same shuffle seed, same
individuals, same classifier -- except where the centre lands.

The draw is the published null's, reproduced exactly: a uniform centre at least `exclude_radius`
bp from the annotated MANE TSS, seeded by sha256(sample|gene|haplotype|draw), so for the eleven
panel genes every scrambled prediction is already on disk from the published null run and is
replayed at zero API cost. Only the three control windows are new: 162 x 3 x 2 = 972 calls.

The draw logic is imported from null_knockout_pigmentation.py rather than copied, so the two nulls
cannot silently diverge.

Usage:
  python3 scripts/experiments/specificity_14gene_null.py
  ... --limit-individuals 2 --genes TPM2      # smoke test
  ... --no-api                                # replay the eleven cached panel genes only
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
for _p in (REPO_ROOT / "src", REPO_ROOT / "notebooks", REPO_ROOT / "scripts" / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd
import torch

CONFIG = REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment_14gene.yaml"
AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
DEFAULT_OUT = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_null_14gene.csv"
METHOD = "random_window"
SCRAMBLE = 100
PANEL = {"MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"}

FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "gene_class", "method",
    "scramble_window", "chrom", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit", "perturbed_strong_logit", "perturbed_weak_logit",
    "delta_log_odds", "baseline_pred", "perturbed_pred", "flipped",
    "h1_crop_abs_delta", "h2_crop_abs_delta", "delta_in",
    "h1_raw_abs_delta", "h2_raw_abs_delta", "from_cache",
    "draw", "exclude_radius", "h1_dist_to_tss", "h2_dist_to_tss",
]


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--limit-individuals", type=int, default=None)
    ap.add_argument("--genes", default=None)
    ap.add_argument("--no-api", action="store_true",
                    help="replay cached scrambles only; never call AlphaGenome")
    ap.add_argument("--draw", type=int, default=0, help="null draw index (matches the published null)")
    ap.add_argument("--exclude-radius", type=int, default=5000,
                    help="bp around the annotated TSS a null draw may not land in")
    args = ap.parse_args()

    from dotenv import load_dotenv
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    from genotype_cnn_alignment_deeplift_summary.annotations import (
        build_tss_tables, build_transcript_extractors, load_gtf, get_gene_tss)
    from genotype_cnn_alignment_deeplift_summary.haplotype import haplotype_local_idx, load_haplotype_fasta
    from genotype_cnn_alignment_deeplift_summary.knockout import (
        apply_scramble, load_raw_prediction, reorder_to_canonical)
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import build_model, load_checkpoint
    from null_knockout_pigmentation import pick_random_target, _seed_for
    from genomics.predictors.genotype_based.config import (
        generate_experiment_name, get_dataset_cache_dir, get_experiment_runs_dir, load_config)
    from genomics.predictors.genotype_based.data.pipeline import (
        _load_split_index, _make_runtime_processed_datasets, _resolve_runtime_dataset_dir)

    load_dotenv(Path.home() / ".env")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"device={device}")

    cfg = load_config(CONFIG)
    genes_all = list(cfg.dataset_input.genes_to_use)
    cache_dir = get_dataset_cache_dir(cfg)
    ds_dir = _resolve_runtime_dataset_dir(cfg)
    _log(f"loading 14-gene dataset ({ds_dir})...")
    with quiet_pipeline_logs():
        full_ds, _tr, _va, _te = _make_runtime_processed_datasets(ds_dir, cache_dir, cfg)

    exp = get_experiment_runs_dir(cfg) / generate_experiment_name(cfg)
    ck = exp / "models" / "best_accuracy.pt"
    if not ck.exists():
        raise SystemExit(f"ABORT: no checkpoint at {ck}")
    if "runs_specificity_14gene" not in str(exp):
        raise SystemExit(f"ABORT: checkpoint is not the isolated 14-gene run: {exp}")
    _log(f"checkpoint {ck}")
    model = build_model(cfg, full_ds, device)
    model = load_checkpoint(model, ck, device)
    model.eval()

    class_names = full_ds.get_class_names()
    idx_to_target = full_ds.idx_to_target
    strong_idx = [i for i, n in idx_to_target.items() if n == "strong pigmentation"][0]
    weak_idx = next(i for i in range(2) if i != strong_idx)
    split = _load_split_index(cache_dir)
    test_ids = split["test"]
    ped = full_ds.dataset_metadata.get("individuals_pedigree", {})
    dataset_dir = Path(cfg.dataset_input.dataset_dir)

    gtf = load_gtf(REPO_ROOT / "notebooks" / ".cache" / "annotations")
    gtf_mane, gtf_pc, _, _, _ = build_transcript_extractors(gtf)
    tss_mane, tss_pc = build_tss_tables(gtf_mane, gtf_pc)

    ag_client = None
    if not args.no_api:
        key = os.environ.get("ALPHAGENOME_API_KEY")
        if not key:
            raise SystemExit("ALPHAGENOME_API_KEY not found (use --no-api to replay cache only).")
        ag_client = dna_client.create(api_key=key)
    organism = dna_client.Organism.HOMO_SAPIENS
    ontology_terms = list(cfg.dataset_input.ontology_terms)

    genes = [g for g in genes_all if g in args.genes.split(",")] if args.genes else genes_all
    sample_ids = test_ids[: args.limit_individuals] if args.limit_individuals else test_ids

    # window geometry, per gene
    # Window length is taken from the reference FASTA, NOT from window_metadata's end-start+1:
    # the three control windows were written by a later builder whose `end` is 2 bp long
    # (524,290 against a 524,288 bp FASTA), and AlphaGenome rejects an interval whose width does
    # not match the sequence. The panel windows agree either way, so this is a strict repair.
    geom = {}
    for g in genes_all:
        wdir = dataset_dir / "references" / "windows" / g
        wm = json.loads((wdir / "window_metadata.json").read_text())
        chrom, s1 = wm["chromosome"], int(wm["start"])
        L = sum(len(line.strip()) for line in (wdir / "ref.window.fa").read_text().splitlines()
                if not line.startswith(">"))
        meta_L = int(wm["end"]) - s1 + 1
        if meta_L != L:
            _log(f"NOTE: {g} window_metadata length {meta_L} != reference FASTA length {L}; using the FASTA.")
        _c, tss0, _strand = get_gene_tss(g, tss_mane, tss_pc)
        geom[g] = {"chrom": chrom, "start_1based": s1, "len": L, "tss_0based": tss0,
                   "interval": genome.Interval(chromosome=chrom, start=s1 - 1, end=s1 - 1 + L)}

    def raw_pred(sid, gene, hap):
        return load_raw_prediction(dataset_dir, sid, gene, hap)

    def crop_rows(sid, gene, hap, array, meta):
        """The raw_center_crop rows the CNN actually consumes for this gene/haplotype."""
        r = full_ds._process_haplotype_raw_center_crop({"rna_seq": array}, {"rna_seq": meta})
        if r is None:
            raise RuntimeError(f"no rows for {sid}/{gene}/{hap}")
        return r

    def build_tensor(sid, overrides=None):
        """raw_center_crop layout: per gene, H1 rows then H2 rows, vstacked (no shared axis).

        Uses the dataset's own normaliser rather than the bare apply_normalization -- it is what
        __getitem__ calls for this layout, and the equality guard below fails loudly if the two
        ever diverge."""
        overrides = overrides or {}
        rows = []
        for g in genes_all:
            for hap in ("H1", "H2"):
                arr, meta = overrides.get((g, hap)) or raw_pred(sid, g, hap)
                rows.append(crop_rows(sid, g, hap, arr, meta))
        stacked = np.vstack(rows).astype(np.float32)
        return full_ds._normalize_features_tensor(torch.FloatTensor(stacked))

    def predict_scrambled(cache_key, sequence, interval):
        npz, mj = AG_CACHE / f"seq_{cache_key}.npz", AG_CACHE / f"seq_{cache_key}_meta.json"
        if npz.exists() and mj.exists():
            return np.load(npz)["values"], json.loads(mj.read_text()), True
        if ag_client is None:
            raise RuntimeError(f"cache miss for {cache_key} and --no-api is set")
        out = ag_client.predict_sequence(
            sequence, organism=organism, requested_outputs=[dna_client.OutputType.RNA_SEQ],
            ontology_terms=ontology_terms, interval=interval)
        vals = out.rna_seq.values
        recs = out.rna_seq.metadata[["ontology_curie", "strand"]].to_dict("records")
        AG_CACHE.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(npz, values=vals)
        mj.write_text(json.dumps(recs))
        return vals, recs, False

    # ---- correctness guard: rebuild one baseline tensor by hand, compare to the dataset's own
    sid_list = [full_ds._sample_id_for_base_index(b) for b in full_ds.valid_sample_indices]
    probe = sample_ids[0]
    if probe in sid_list:
        mine = build_tensor(probe).squeeze()
        theirs = full_ds[sid_list.index(probe)][0].squeeze()
        if mine.shape != theirs.shape:
            raise SystemExit(f"ABORT: shape {tuple(mine.shape)} vs {tuple(theirs.shape)}")
        md = float((mine - theirs).abs().max())
        if md > 1e-4:
            raise SystemExit(f"ABORT: reconstructed baseline differs (max|diff|={md:.3e})")
        _log(f"guard OK: max|diff|={md:.2e}, tensor {tuple(mine.shape)}")
    else:
        _log(f"WARNING: {probe} not in this dataset index; guard skipped.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if args.out.exists():
        prev = pd.read_csv(args.out)
        done = set(zip(prev["sample_id"], prev["gene"]))
        _log(f"resuming: {len(done)} rows present")
    write_header = not args.out.exists()

    total = len(sample_ids) * len(genes)
    _log(f"plan: {len(sample_ids)} individuals x {len(genes)} genes = {total} rows -> {args.out}")

    baseline_cache = {}
    n_ok = n_fail = n_api = 0
    t0 = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            w.writeheader(); fh.flush()
        for sid in sample_ids:
            for gene in genes:
                if (sid, gene) in done:
                    continue
                try:
                    gm = geom[gene]
                    overrides, per_hap = {}, {}
                    cached_all = True
                    for hap in ("H1", "H2"):
                        tss_loc = haplotype_local_idx(dataset_dir, sid, gene, hap,
                                                      gm["start_1based"], gm["tss_0based"] + 1)
                        seq = load_haplotype_fasta(dataset_dir, sid, gene, hap)
                        loc = pick_random_target(len(seq), tss_loc, args.exclude_radius, SCRAMBLE,
                                                 _seed_for(sid, gene, hap, args.draw))
                        mod_seq, _o, _s, _a, _b = apply_scramble(seq, loc, window_size=SCRAMBLE)
                        orig_arr, orig_meta = raw_pred(sid, gene, hap)
                        canonical = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
                        key = f"{sid}_{gene}_{hap}_{METHOD}d{args.draw}_{loc}_scramble{SCRAMBLE}"
                        vals, recs, hit = predict_scrambled(key, mod_seq, gm["interval"])
                        cached_all = cached_all and hit
                        if not hit:
                            n_api += 1
                        mod = reorder_to_canonical(vals, pd.DataFrame(recs), canonical)
                        overrides[(gene, hap)] = (mod, orig_meta)

                        # delivered perturbation, raw and as it reaches the CNN input
                        base = crop_rows(sid, gene, hap, orig_arr, orig_meta)
                        pert = crop_rows(sid, gene, hap, mod, orig_meta)
                        per_hap[hap] = {
                            "loc": loc, "dist_to_tss": loc - tss_loc,
                            "crop_abs_delta": float(np.abs(np.asarray(pert, float) - np.asarray(base, float)).sum()),
                            "raw_abs_delta": float(np.abs(mod.astype(float) - orig_arr.astype(float)).sum()),
                        }

                    if sid not in baseline_cache:
                        with torch.no_grad():
                            bt = build_tensor(sid)
                            baseline_cache[sid] = model(bt.unsqueeze(0).float().to(device))[0].cpu().numpy()
                    bl = baseline_cache[sid]
                    with torch.no_grad():
                        pt = build_tensor(sid, overrides=overrides)
                        pl = model(pt.unsqueeze(0).float().to(device))[0].cpu().numpy()

                    b_lo = float(bl[strong_idx] - bl[weak_idx])
                    p_lo = float(pl[strong_idx] - pl[weak_idx])
                    pr = ped.get(sid, {})
                    w.writerow({
                        "sample_id": sid, "population": pr.get("population"),
                        "superpopulation": pr.get("superpopulation"),
                        "true_label": full_ds._get_target_value(pr),
                        "gene": gene, "gene_class": "panel" if gene in PANEL else "control",
                        "method": METHOD, "scramble_window": SCRAMBLE, "chrom": gm["chrom"],
                        "h1_target_local_idx": per_hap["H1"]["loc"],
                        "h2_target_local_idx": per_hap["H2"]["loc"],
                        "baseline_strong_logit": float(bl[strong_idx]),
                        "baseline_weak_logit": float(bl[weak_idx]),
                        "perturbed_strong_logit": float(pl[strong_idx]),
                        "perturbed_weak_logit": float(pl[weak_idx]),
                        "delta_log_odds": p_lo - b_lo,
                        "baseline_pred": class_names[int(bl.argmax())],
                        "perturbed_pred": class_names[int(pl.argmax())],
                        "flipped": int(bl.argmax() != pl.argmax()),
                        "h1_crop_abs_delta": per_hap["H1"]["crop_abs_delta"],
                        "h2_crop_abs_delta": per_hap["H2"]["crop_abs_delta"],
                        "delta_in": per_hap["H1"]["crop_abs_delta"] + per_hap["H2"]["crop_abs_delta"],
                        "h1_raw_abs_delta": per_hap["H1"]["raw_abs_delta"],
                        "h2_raw_abs_delta": per_hap["H2"]["raw_abs_delta"],
                        "from_cache": int(cached_all),
                        "draw": args.draw, "exclude_radius": args.exclude_radius,
                        "h1_dist_to_tss": per_hap["H1"]["dist_to_tss"],
                        "h2_dist_to_tss": per_hap["H2"]["dist_to_tss"],
                    })
                    fh.flush()
                    n_ok += 1
                    if n_ok % 25 == 0:
                        el = time.monotonic() - t0
                        _log(f"{n_ok}/{total} rows ({n_ok/el:.2f}/s), api_calls={n_api}, failed={n_fail}, last={sid}/{gene}")
                except Exception as exc:  # noqa: BLE001
                    n_fail += 1
                    _log(f"FAILED {sid}/{gene}: {exc!r}")
    _log(f"DONE: ok={n_ok} failed={n_fail} api_calls={n_api} -> {args.out}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
