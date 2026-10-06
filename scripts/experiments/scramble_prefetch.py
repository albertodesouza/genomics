#!/usr/bin/env python3
"""Populate the promoter-scramble cache for genes that have no classifier yet.

`single_gene_knockdown_replay.py` refuses to emit numbers on a cache miss, and rightly
so -- a partially cached replay silently drops individuals. But a newly built control
gene has no scrambles cached and no arm trained, so nothing can generate them: the
existing generators (`specificity_14gene_knockdown.py`, `..._null.py`) each need a
checkpoint to forward through, and for a new gene that checkpoint does not exist until
after this cache is populated. This script breaks that circle by doing only the half
that costs money and none of the half that needs a model.

It reproduces the two cache keys bit for bit, so the replays afterwards are pure cache
hits:

    biology_tss    {sid}_{gene}_{hap}_biology_tss_{loc}_scramble100
    random_window  {sid}_{gene}_{hap}_random_windowd{draw}_{loc}_scramble100

Both the scramble (`apply_scramble`, fixed seed, composition-preserving) and the null
draw (`pick_random_target` under `_seed_for`) are imported from the same modules the
published runs used, never reimplemented, so a gene that *is* already cached comes back
100% hits -- which is what `--verify-cached` asserts.

Usage:
  python3 scripts/experiments/scramble_prefetch.py --genes PPP1R3E,ECHDC3
  ... --verify-cached TYR          # no-op self-test: must report 324/324 hits, 0 calls
  ... --no-api                     # audit what is missing without spending anything
"""
from __future__ import annotations

import argparse
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

AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
CONFIG = "configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment_14gene.yaml"
# The 162 test individuals are read from a persisted knockdown CSV rather than from a
# dataset view: the views are deleted per arm by the sweep, and the split is identical
# across every gene subset (seed 13, family_aware, strict_determinism).
SPLIT_SOURCE = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout_14gene.csv"
EXPECTED_TEST = 162
# The replay needs the (sample, method, target index) triples to rebuild the cache keys, and
# it reads them from a knockdown CSV. A gene built after the last knockdown appears in no such
# CSV, so the loci computed here -- the same ones the keys were built from -- are written out
# as a manifest for the replay to fall back on.
MANIFEST = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/scramble_manifest.csv"
SCRAMBLE = 100
EXCLUDE_RADIUS = 5000


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="comma-separated gene symbols")
    ap.add_argument("--draw", type=int, default=1, help="null draw index (default 1)")
    ap.add_argument("--exclude-radius", type=int, default=EXCLUDE_RADIUS)
    ap.add_argument("--no-api", action="store_true", help="report misses, spend nothing")
    ap.add_argument("--skip-null", action="store_true")
    ap.add_argument("--verify-cached", action="store_true",
                    help="assert every key is already cached (self-test on a known gene)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from dotenv import load_dotenv
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    from genotype_cnn_alignment_deeplift_summary.annotations import (
        build_tss_tables, build_transcript_extractors, load_gtf, get_gene_tss)
    from genotype_cnn_alignment_deeplift_summary.haplotype import (
        haplotype_local_idx, load_haplotype_fasta)
    from genotype_cnn_alignment_deeplift_summary.knockout import apply_scramble
    from genomics.predictors.genotype_based.config import load_config
    sys.path.insert(0, str(REPO_ROOT / "scripts" / "experiments"))
    from null_knockout_pigmentation import pick_random_target, _seed_for

    load_dotenv(Path.home() / ".env")
    cfg = load_config(CONFIG)
    dataset_dir = Path(cfg.dataset_input.dataset_dir)
    ontology_terms = list(cfg.dataset_input.ontology_terms)

    if not SPLIT_SOURCE.exists():
        raise SystemExit(f"ABORT: {SPLIT_SOURCE} missing; cannot resolve the test split.")
    test_ids = list(pd.read_csv(SPLIT_SOURCE)["sample_id"].unique())
    if len(test_ids) != EXPECTED_TEST:
        raise SystemExit(f"ABORT: split source has {len(test_ids)} individuals, expected {EXPECTED_TEST}")
    _log(f"test split: {len(test_ids)} individuals (from {SPLIT_SOURCE.name})")

    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    gtf = load_gtf(REPO_ROOT / "notebooks" / ".cache" / "annotations")
    gtf_mane, gtf_pc, _, _, _ = build_transcript_extractors(gtf)
    tss_mane, tss_pc = build_tss_tables(gtf_mane, gtf_pc)

    geom = {}
    for g in genes:
        wdir = dataset_dir / "references" / "windows" / g
        wm_path = wdir / "window_metadata.json"
        if not wm_path.exists():
            raise SystemExit(f"ABORT: no window built for {g} ({wm_path} missing). "
                             f"Run the window build first.")
        wm = json.loads(wm_path.read_text())
        chrom, s1 = wm["chromosome"], int(wm["start"])
        # Window length from the reference FASTA, not from window_metadata's end-start+1:
        # the control windows are written by a builder whose `end` is 2 bp long, and
        # AlphaGenome rejects an interval whose width does not match the sequence.
        L = sum(len(line.strip()) for line in (wdir / "ref.window.fa").read_text().splitlines()
                if not line.startswith(">"))
        if int(wm["end"]) - s1 + 1 != L:
            _log(f"NOTE: {g} window_metadata length {int(wm['end']) - s1 + 1} != FASTA {L}; using the FASTA.")
        _c, tss0, _strand = get_gene_tss(g, tss_mane, tss_pc)
        geom[g] = {"chrom": chrom, "start_1based": s1, "len": L, "tss_0based": tss0,
                   "interval": genome.Interval(chromosome=chrom, start=s1 - 1, end=s1 - 1 + L)}
        _log(f"{g}: {chrom}:{s1} len={L} tss={tss0}")

    ag_client = None
    if not (args.no_api or args.verify_cached):
        key = os.environ.get("ALPHAGENOME_API_KEY")
        if not key:
            raise SystemExit("ALPHAGENOME_API_KEY not found.")
        ag_client = dna_client.create(api_key=key)
    organism = dna_client.Organism.HOMO_SAPIENS
    AG_CACHE.mkdir(parents=True, exist_ok=True)

    methods = ["biology_tss"] if args.skip_null else ["biology_tss", "random_window"]
    manifest = {}          # (gene, sid, method) -> {h1, h2}
    n_hit = n_call = n_miss = n_fail = 0
    total = len(test_ids) * len(genes) * 2 * len(methods)
    _log(f"plan: {len(test_ids)} x {len(genes)} genes x 2 hap x {len(methods)} methods = {total} keys")
    t0 = time.monotonic()

    for gene in genes:
        gm = geom[gene]
        for sid in test_ids:
            for hap in ("H1", "H2"):
                try:
                    tss_loc = haplotype_local_idx(dataset_dir, sid, gene, hap,
                                                  gm["start_1based"], gm["tss_0based"] + 1)
                    seq = load_haplotype_fasta(dataset_dir, sid, gene, hap)
                except Exception as e:                                    # noqa: BLE001
                    _log(f"FAIL {sid}/{gene}/{hap}: {type(e).__name__}: {e}")
                    n_fail += 1
                    continue
                for method in methods:
                    if method == "biology_tss":
                        loc = tss_loc
                        ck = f"{sid}_{gene}_{hap}_biology_tss_{loc}_scramble{SCRAMBLE}"
                    else:
                        loc = pick_random_target(len(seq), tss_loc, args.exclude_radius, SCRAMBLE,
                                                 _seed_for(sid, gene, hap, args.draw))
                        ck = f"{sid}_{gene}_{hap}_random_windowd{args.draw}_{loc}_scramble{SCRAMBLE}"
                    key_method = ck.split(f"{sid}_{gene}_{hap}_")[1].rsplit(f"_{loc}_", 1)[0]
                    manifest.setdefault((gene, sid, key_method), {})[hap] = loc
                    npz, mj = AG_CACHE / f"seq_{ck}.npz", AG_CACHE / f"seq_{ck}_meta.json"
                    if npz.exists() and mj.exists():
                        n_hit += 1
                        continue
                    if args.verify_cached:
                        raise SystemExit(f"ABORT: --verify-cached but {ck} is missing")
                    if ag_client is None:
                        n_miss += 1
                        continue
                    try:
                        mod_seq, *_ = apply_scramble(seq, loc, window_size=SCRAMBLE)
                        out = ag_client.predict_sequence(
                            mod_seq, organism=organism,
                            requested_outputs=[dna_client.OutputType.RNA_SEQ],
                            ontology_terms=ontology_terms, interval=gm["interval"])
                        np.savez_compressed(npz, values=out.rna_seq.values)
                        mj.write_text(json.dumps(
                            out.rna_seq.metadata[["ontology_curie", "strand"]].to_dict("records")))
                        n_call += 1
                    except Exception as e:                                # noqa: BLE001
                        _log(f"FAIL {ck}: {type(e).__name__}: {e}")
                        n_fail += 1
                        # A half-written pair would read as a hit next run and poison the replay.
                        for p in (npz, mj):
                            if p.exists():
                                p.unlink()
            done = n_hit + n_call + n_miss + n_fail
            if done and done % 100 == 0:
                el = time.monotonic() - t0
                _log(f"{done}/{total} ({done/el:.2f}/s) hits={n_hit} calls={n_call} "
                     f"missing={n_miss} failed={n_fail} last={sid}/{gene}")

    # Written whole rather than appended: a row is only meaningful if both haplotypes
    # resolved, and a rerun recomputes the same deterministic loci anyway.
    import csv as _csv
    rows = [{"gene": g, "sample_id": s_, "method": m, "scramble_window": SCRAMBLE,
             "h1_target_local_idx": d["H1"], "h2_target_local_idx": d["H2"]}
            for (g, s_, m), d in sorted(manifest.items()) if "H1" in d and "H2" in d]
    if rows:
        prior = []
        if MANIFEST.exists():
            with open(MANIFEST, newline="", encoding="utf-8") as fh:
                keep = {(r["gene"], r["sample_id"], r["method"]) for r in rows}
                prior = [r for r in _csv.DictReader(fh)
                         if (r["gene"], r["sample_id"], r["method"]) not in keep]
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
            w = _csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for r in prior + rows:
                w.writerow(r)
        _log(f"manifest: {len(rows)} rows for this run, {len(prior)} kept -> {MANIFEST.name}")

    _log(f"DONE: hits={n_hit} calls={n_call} missing={n_miss} failed={n_fail} of {total}")
    if args.verify_cached:
        _log("verify-cached passed: every key already on disk.")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
