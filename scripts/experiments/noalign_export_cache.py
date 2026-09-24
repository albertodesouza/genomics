#!/usr/bin/env python3
"""No-alignment ablation (Q1, scoped to "no alignment" only, per user request): does the
DynamicIndelAligner coordinate-sharing step (Section 3.3 of the paper) matter, or would a
naive crop -- reading the same nominal offset into each individual's own AlphaGenome
output, with no indel-aware remapping -- rank the genes just as well?

Reuses the already-cached raw per-individual AlphaGenome predictions at
individuals/<sample>/windows/<gene>/predictions_H{1,2}/rna_seq.npz, shape (524288, 6):
one row per position in THAT INDIVIDUAL'S OWN personal coordinates (the DynamicIndelAligner
consumes this same array and remaps it; here we skip that remapping entirely). No new
AlphaGenome API calls.

Column selection mirrors the published single-strand config exactly: melanocyte
(CL:1000458) at the gene's own strand, read from that gene's published config
(track_strands), via the metadata JSON sitting next to the .npz (column order is stable
across genes: 0/3 = CL:0000346 +/-, 1/4 = CL:1000458 +/-, 2/5 = CL:2000092 +/-, but we
read it rather than assume it).

Crop: the central 32768 positions of the raw 524288-length array (indices
[262144-16384 : 262144+16384]) -- the SAME nominal window-relative offset for every
individual, which is exactly what "no alignment" means: an insertion upstream of this
offset in one individual and not another silently shifts what each of them has at this
index, uncorrected.

Output: results/cache/genotype_based_predictor/noalign_signal/<gene>.npz, one array per
sample_id x haplotype ("HG00096__H1", "HG00096__H2"), float32 (32768,), raw (pre-log1p)
signal -- normalization is fit and applied at train time, on the train split only, exactly
like the aligned arm (Eq. 2), so this cache holds the one thing that must NOT be shared
with the aligned arm's cache: the un-remapped signal itself.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path("/home/breno/I2CA/genomics")
DATASET_DIR = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
CACHE_DIR = REPO / "results/cache/genotype_based_predictor/noalign_signal"
SELECTED_SAMPLES_CSV = DATASET_DIR / "selected_samples.csv"
RAW_LEN = 524288
CROP_SIZE = 32768
CENTER = RAW_LEN // 2
LO, HI = CENTER - CROP_SIZE // 2, CENTER + CROP_SIZE // 2


def _cohort_sample_ids() -> list[str]:
    ids = []
    with open(SELECTED_SAMPLES_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.append(row["SampleID"])
    return ids


def _gene_strand(gene: str) -> str:
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {gene} resolves to {len(hits)} published configs, expected 1")
    cfg = yaml.safe_load(hits[0].read_text())
    strands = cfg["dataset_input"]["track_strands"]
    if len(strands) != 1:
        raise SystemExit(f"ABORT: {gene} published config is not single-strand: {strands}")
    return strands[0]


def _melanocyte_column(meta_path: Path, strand: str) -> int:
    meta = json.loads(meta_path.read_text())["metadata"]
    hits = [i for i, m in enumerate(meta)
            if m["ontology_curie"] == "CL:1000458" and m["strand"] == strand]
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {meta_path} has {len(hits)} melanocyte/{strand} columns")
    return hits[0]


def export_gene(gene: str, force: bool = False) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = CACHE_DIR / f"{gene.upper()}.npz"
    if out_path.exists() and not force:
        print(f"{gene}: cache exists, skipping ({out_path})")
        return

    strand = _gene_strand(gene)
    sample_ids = _cohort_sample_ids()
    arrays: dict[str, np.ndarray] = {}
    col = None
    for sid in sample_ids:
        for h in ("H1", "H2"):
            base = DATASET_DIR / "individuals" / sid / "windows" / gene / f"predictions_{h}"
            npz_path = base / "rna_seq.npz"
            if col is None:
                col = _melanocyte_column(base / "rna_seq_metadata.json", strand)
            d = np.load(npz_path)
            v = d["values"]
            if v.shape[0] != RAW_LEN:
                raise SystemExit(f"ABORT: {npz_path} has {v.shape[0]} rows, expected {RAW_LEN}")
            arrays[f"{sid}__{h}"] = v[LO:HI, col].astype(np.float32)
    np.savez_compressed(out_path, **arrays)
    size_mb = out_path.stat().st_size / 1e6
    print(f"{gene}: strand={strand} col={col} DONE: {len(sample_ids)} individuals -> "
          f"{out_path} ({size_mb:.1f} MB)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    export_gene(args.gene, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
