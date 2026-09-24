#!/usr/bin/env python3
"""Offline bulk promoter-location knockout experiment for the pigmentation `dita_no_masks` CNN2.

For every (individual, gene, method) in the test split, scrambles a 100 bp window at the
candidate promoter location identified by `method` on both haplotypes of that individual's real
consensus sequence, re-predicts with AlphaGenome, and records the shift in the CNN's
log-odds(strong/weak pigmentation). Three methods, each haplotype-specific (each haplotype's own
consensus sequence is re-predicted with AlphaGenome to locate the target, not a shared
reference-genome location):
  - biology_tss: the GENCODE MANE Select transcript's literature TSS.
  - cage_melanocyte: the CAGE (melanocyte) summit within +/-5kb of that TSS.
  - cage_gene_start: the CAGE (melanocyte) summit within +/-5kb of the gene's strand-aware
    genomic 5' edge (no TSS annotation used to anchor the search).

This mirrors notebooks/genotype_cnn_alignment_deeplift_summary.ipynb's Section 9
(`knockout_gene_experiment`) exactly -- see that notebook for the single-individual walkthrough
this bulk run aggregates. Results replace
results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout.csv (the prior
CSV used a stale reference-genome-only CAGE location and an `alphagenome_atac` method that was
later disabled for lacking ontology-matched signal; both are gone here).

Meant to run unattended (nohup/setsid) over hours: ~162 individuals x 11 genes x 3 methods, with
per-row CSV writes and skip-on-resume, so a crash or interruption loses at most the in-flight row.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/bulk_knockout_pigmentation.py
  ... --limit-individuals 2 --genes MC1R,TYR                      # smoke test
  ... --out results/genotype_based_predictor/knockout_bulk/smoke.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _path in (REPO_ROOT / "src", REPO_ROOT / "notebooks"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
os.chdir(REPO_ROOT)

CSV_FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "method",
    "scramble_window", "chrom", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit", "perturbed_strong_logit", "perturbed_weak_logit",
    "baseline_pred", "perturbed_pred", "flipped",
]
METHODS = ["biology_tss", "cage_melanocyte", "cage_gene_start"]
DEFAULT_OUT = REPO_ROOT / "results" / "genotype_based_predictor" / "knockout_bulk" / "pigmentation_test_split_knockout.csv"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


def _build_context(device):
    from dotenv import load_dotenv
    from alphagenome.models import dna_client

    from genotype_cnn_alignment_deeplift_summary.annotations import build_tss_tables, build_transcript_extractors, load_gtf
    from genotype_cnn_alignment_deeplift_summary.context import KnockoutContext
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import build_model, load_checkpoint

    from genomics.predictors.genotype_based.config import (
        generate_experiment_name, get_dataset_cache_dir, get_experiment_runs_dir, load_config,
    )
    from genomics.predictors.genotype_based.data.pipeline import (
        _load_split_index, _make_runtime_processed_datasets, _resolve_runtime_dataset_dir,
    )

    load_dotenv(Path.home() / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY not found in the environment or ~/.env.")
    ag_client = dna_client.create(api_key=api_key)
    organism = dna_client.Organism.HOMO_SAPIENS

    ko_config = load_config(REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml")
    ko_genes = list(ko_config.dataset_input.genes_to_use)
    ko_ontology_terms = list(ko_config.dataset_input.ontology_terms)

    cache_dir = get_dataset_cache_dir(ko_config)
    runtime_dataset_dir = _resolve_runtime_dataset_dir(ko_config)
    _log(f"Loading dataset from {runtime_dataset_dir} (cache={cache_dir})...")
    with quiet_pipeline_logs():
        ko_full_ds, _train_ds, _val_ds, _test_ds = _make_runtime_processed_datasets(runtime_dataset_dir, cache_dir, ko_config)

    experiment_dir = get_experiment_runs_dir(ko_config) / generate_experiment_name(ko_config)
    checkpoint_path = experiment_dir / "models" / "best_accuracy.pt"
    _log(f"Loading checkpoint {checkpoint_path}...")
    model = build_model(ko_config, ko_full_ds, device)
    model = load_checkpoint(model, checkpoint_path, device)

    class_names = ko_full_ds.get_class_names()
    idx_to_target = ko_full_ds.idx_to_target
    strong_idx = [i for i, name in idx_to_target.items() if name == "strong pigmentation"][0]
    weak_idx = next(i for i in range(2) if i != strong_idx)

    split_index = _load_split_index(cache_dir)
    test_sample_ids = split_index["test"]
    pedigree = ko_full_ds.dataset_metadata.get("individuals_pedigree", {})

    gtf_cache_dir = REPO_ROOT / "notebooks" / ".cache" / "annotations"
    _log("Loading GTF annotations...")
    gtf = load_gtf(gtf_cache_dir)
    gtf_mane, gtf_protein_coding, transcript_extractor_mane, transcript_extractor_coding, gene_id_map = (
        build_transcript_extractors(gtf)
    )
    tss_df_mane, tss_df_coding = build_tss_tables(gtf_mane, gtf_protein_coding)
    gtf_gene_rows = gtf[gtf["Feature"] == "gene"]

    # Same cache dir the interactive notebook uses for CAGE curves / scrambled-sequence
    # predictions, so any AlphaGenome calls already made there (or by a prior run of this script)
    # are reused instead of re-billed.
    ag_cache_dir = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
    ag_cache_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = REPO_ROOT / "notebooks" / ".cache" / "genotype_cnn_alignment_deeplift"
    fig_dir.mkdir(parents=True, exist_ok=True)

    class_colors = {strong_idx: "#2b2118", weak_idx: "#d9a441"}
    ctx = KnockoutContext(
        ag_client=ag_client, organism=organism, config=ko_config, genes=ko_genes,
        ontology_terms=ko_ontology_terms, model=model, device=device, dataset_dir=Path(ko_config.dataset_input.dataset_dir),
        full_ds=ko_full_ds, normalization_params=ko_full_ds.normalization_params,
        strong_idx=strong_idx, weak_idx=weak_idx, class_names=class_names, class_colors=class_colors,
        target_idx=strong_idx, other_idx=weak_idx, cache_dir=ag_cache_dir, fig_dir=fig_dir,
        tss_df_mane=tss_df_mane, tss_df_coding=tss_df_coding, gtf_gene_rows=gtf_gene_rows,
        gene_id_map=gene_id_map, transcript_extractor_mane=transcript_extractor_mane,
        transcript_extractor_coding=transcript_extractor_coding,
    )
    return ctx, ko_genes, test_sample_ids, pedigree, class_names, strong_idx, weak_idx


def _compute_row(ctx, sample_id, gene, method, pedigree, class_names):
    from genotype_cnn_alignment_deeplift_summary.knockout import knockout_gene_experiment

    result = knockout_gene_experiment(ctx, sample_id, gene, method)
    baseline_logits, perturbed_logits = result["baseline_logits"], result["perturbed_logits"]
    baseline_pred = class_names[int(baseline_logits.argmax())]
    perturbed_pred = class_names[int(perturbed_logits.argmax())]
    pedigree_row = pedigree.get(sample_id, {})
    true_label = ctx.full_ds._get_target_value(pedigree_row)

    return {
        "sample_id": sample_id,
        "population": pedigree_row.get("population"),
        "superpopulation": pedigree_row.get("superpopulation"),
        "true_label": true_label,
        "gene": gene,
        "method": method,
        "scramble_window": result["scramble_window"],
        "chrom": result["chrom"],
        "h1_target_local_idx": result["hap_tracks"]["H1"]["marker_local_idx"][method],
        "h2_target_local_idx": result["hap_tracks"]["H2"]["marker_local_idx"][method],
        "baseline_strong_logit": float(baseline_logits[ctx.strong_idx]),
        "baseline_weak_logit": float(baseline_logits[ctx.weak_idx]),
        "perturbed_strong_logit": float(perturbed_logits[ctx.strong_idx]),
        "perturbed_weak_logit": float(perturbed_logits[ctx.weak_idx]),
        "baseline_pred": baseline_pred,
        "perturbed_pred": perturbed_pred,
        "flipped": baseline_pred != perturbed_pred,
    }


def _already_done(out_path: Path) -> set:
    if not out_path.exists():
        return set()
    import pandas as pd
    existing = pd.read_csv(out_path)
    return set(zip(existing["sample_id"], existing["gene"], existing["method"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit-individuals", type=int, default=None)
    parser.add_argument("--genes", type=str, default=None, help="comma-separated subset of genes")
    parser.add_argument("--methods", type=str, default=None, help="comma-separated subset of methods")
    args = parser.parse_args()

    import torch
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    _log(f"Using device: {device}")

    ctx, ko_genes, test_sample_ids, pedigree, class_names, strong_idx, weak_idx = _build_context(device)

    genes = [g for g in ko_genes if g in args.genes.split(",")] if args.genes else ko_genes
    methods = [m for m in METHODS if m in args.methods.split(",")] if args.methods else METHODS
    sample_ids = test_sample_ids[: args.limit_individuals] if args.limit_individuals else test_sample_ids

    total = len(sample_ids) * len(genes) * len(methods)
    _log(f"Plan: {len(sample_ids)} individuals x {len(genes)} genes x {len(methods)} methods = {total} rows -> {args.out}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = _already_done(args.out)
    if done:
        _log(f"Resuming: {len(done)} rows already present in {args.out}, skipping those.")
    write_header = not args.out.exists()

    n_written, n_skipped, n_failed = 0, 0, 0
    t_start = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if write_header:
            writer.writeheader()
            f.flush()

        for sample_id in sample_ids:
            for gene in genes:
                for method in methods:
                    if (sample_id, gene, method) in done:
                        n_skipped += 1
                        continue
                    try:
                        row = _compute_row(ctx, sample_id, gene, method, pedigree, class_names)
                    except Exception as exc:  # noqa: BLE001 - keep the run alive, log and move on
                        n_failed += 1
                        _log(f"FAILED {sample_id}/{gene}/{method}: {exc!r}")
                        continue
                    writer.writerow(row)
                    f.flush()
                    n_written += 1
                    if n_written % 25 == 0:
                        elapsed = time.monotonic() - t_start
                        rate = n_written / elapsed if elapsed > 0 else 0
                        _log(f"progress: {n_written} written, {n_skipped} skipped, {n_failed} failed "
                             f"({rate:.2f} rows/s, last={sample_id}/{gene}/{method})")

    _log(f"Done. {n_written} written, {n_skipped} skipped (resumed), {n_failed} failed. Output: {args.out}")
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
