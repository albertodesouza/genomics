#!/usr/bin/env python3
"""One-off: stage new candidate pigmentation genes into `notebooks/data/` (reference
window + 1000 Genomes variant slice) and generate their AlphaGenome RNA-seq predictions
for every individual in the notebook series' dataset.

This only registers the genes in `notebooks/data/genes.json` (the general window/variant
registry) -- it deliberately does NOT touch `notebooks/data/experiment.json`'s active
`"genes"` list, so the live 01-05 notebook pipeline is unaffected while this runs. Once
generation finishes, adding a gene to the active list and re-running notebook 3/4/5 is a
separate, fast step (no new API calls needed at that point).

Safe to interrupt and resume: `generate_predictions` skips any (gene, sample, haplotype)
that's already on disk, so re-running this script picks up where it left off.

Usage (inside tmux, from the repo root):

    conda activate genomics
    cd notebooks
    python3 ../scripts/maintenance/add_pigmentation_genes.py 2>&1 | tee -a /tmp/add_genes.log

Estimated cost for the current gene list (measured empirically over the prior 10,716-call
run: 585.1 min / 10,716 calls =~ 3.3s/call):
1 gene x 1072 individuals x 2 haplotypes = 2,144 AlphaGenome calls -> roughly 2 hours.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

NOTEBOOK_DIR = Path(__file__).resolve().parent.parent.parent / "notebooks"
sys.path.insert(0, str(NOTEBOOK_DIR))

import pandas as pd  # noqa: E402

from utils import prepare_data  # noqa: E402
from utils import setup, predictions as predictions_utils, data as data_utils  # noqa: E402

# Literature-curated addition to the active panel (MC1R, TYRP1, TYR, SLC45A2, MFSD12,
# OCA2, HERC2, SLC24A5): ASIP, the canonical MC1R antagonist and one of the most
# consistently replicated pigmentation GWAS hits (Kanetsky et al. 2002) -- see the
# literature table at the top of notebooks/01_overview_and_vep.ipynb for the full panel
# and per-gene citations. Previously generated candidates (MITF, TPCN2, IRF4, KITLG,
# BNC2 -- from an earlier, since-abandoned automated ranking pass) are left in place on
# disk but not reprocessed here.
NEW_GENES = ["ASIP"]


def stage_gene_data(genes: list, gtf: pd.DataFrame, chr_prefix: str, sample_ids: list) -> dict:
    """Extracts each new gene's reference window FASTA + sliced variant VCF (skips any
    gene already present in genes.json), returns the updated genes.json dict."""
    data_dir = NOTEBOOK_DIR / "data"
    genes_json = json.loads((data_dir / "genes.json").read_text())

    for gene in genes:
        if gene in genes_json:
            print(f"[INFO] {gene}: already registered in genes.json, skipping data prep")
            continue
        chrom, start, end = prepare_data.gene_window(gtf, gene, chr_prefix)
        prepare_data.extract_ref_window_fasta(
            prepare_data.DEFAULT_REF_FASTA, chrom, start, end, data_dir / "references" / gene / "ref.window.fa"
        )
        print(f"[INFO] {gene}: window={chrom}:{start}-{end}")
        vcf_path = prepare_data.DEFAULT_VCF_PATTERN.format(chrom=chrom)
        prepare_data.slice_gene_variants(gene, chrom, start, end, vcf_path, sample_ids)
        genes_json[gene] = {"chrom": chrom, "start": start, "end": end}

    (data_dir / "genes.json").write_text(json.dumps(genes_json, indent=2))
    print(f"[INFO] genes.json updated ({len(genes_json)} total registered genes)")
    return genes_json


def main() -> None:
    prepare_data.require_tools()

    data_dir = NOTEBOOK_DIR / "data"
    individuals = json.loads((data_dir / "individuals.json").read_text())
    sample_ids = sorted(individuals)

    print(f"[1/2] Staging reference/variant data for {len(NEW_GENES)} gene(s): {NEW_GENES}")
    chr_prefix = prepare_data.detect_chr_prefix(prepare_data.DEFAULT_REF_FASTA)
    print("[INFO] Loading GTF (cached after first download)...")
    gtf = pd.read_feather(prepare_data.GTF_URL)
    gene_rows = stage_gene_data(NEW_GENES, gtf, chr_prefix, sample_ids)

    total_calls = len(NEW_GENES) * len(sample_ids) * 2
    print(f"\n[2/2] Generating predictions: {len(NEW_GENES)} genes x {len(sample_ids)} individuals x 2 "
          f"haplotypes = {total_calls} AlphaGenome calls (already-cached ones are skipped automatically)")

    ctx = setup.load_experiment(NOTEBOOK_DIR)
    client, _ = setup.create_client()

    t0 = time.time()

    def on_progress(done, total, gene, sample_id, haplotype, status):
        if status == "generated" and (done % 20 == 0 or done == total):
            elapsed = time.time() - t0
            rate = elapsed / done if done else 0
            eta_h = (total - done) * rate / 3600
            print(f"[{done}/{total}] {gene} {sample_id} {haplotype}: {status} "
                  f"({elapsed/60:.1f} min elapsed, ~{eta_h:.1f}h remaining at current rate)")

    results = predictions_utils.generate_predictions(
        client,
        genes=NEW_GENES,
        gene_rows=gene_rows,
        gene_paths_fn=data_utils.gene_paths,
        sample_ids=sample_ids,
        ontology_terms=ctx.ONTOLOGY_TERMS,  # all 3 terms, matching the existing prediction cache
        output_dir=ctx.PREDICTIONS_CACHE_DIR,
        consensus_cache_dir=ctx.CONSENSUS_CACHE_DIR,
        on_progress=on_progress,
    )

    n_generated = sum(1 for r in results if r["status"] == "generated")
    n_cached = sum(1 for r in results if r["status"] == "cached")
    print(f"\n[DONE] {n_generated} generated, {n_cached} already cached, {len(results)} total "
          f"in {(time.time() - t0) / 60:.1f} min.")


if __name__ == "__main__":
    main()
