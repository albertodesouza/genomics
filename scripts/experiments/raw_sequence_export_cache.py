#!/usr/bin/env python3
"""Build a one-hot raw-sequence cache per gene, aligned to the same coordinate axis
AlphaGenome's RNA-seq tensor uses -- the sequence-only control for reviewer Q2c ("feed
raw sequence directly, no AlphaGenome, with an otherwise identical CNN").

Reuses genomics.predictors.genotype_based.alignment.export_aligned_dna.export_aligned_dna
verbatim -- the same DynamicIndelAligner that builds the RNA-seq arm's shared coordinate
system, called with center_window_size=32768 (its own default, matching the crop the CNN
reads) -- so this control needs no new alignment code, only a one-hot encoding step on
top of an already-existing, already-tested export path.

Pure CPU/IO (bcftools + numpy), no CUDA, so it can run independently of and concurrently
with whatever GPU training is queued.

Output: results/cache/genotype_based_predictor/raw_sequence_onehot/<gene>.npz, holding
one array per sample_id x haplotype, e.g. keys "HG00096__H1", "HG00096__H2", each
float32 (4, 32768) in A/C/G/T channel order. Any character outside ACGT (gap padding
from an insertion the individual does not carry, N, etc.) encodes as an all-zero column
-- the same "no information here" convention the RNA-seq tensor uses for its own gap
positions.

Usage:
  python3 scripts/experiments/raw_sequence_export_cache.py <config.yaml> --gene OCA2
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import yaml

from genomics.predictors.genotype_based.alignment.export_aligned_dna import export_aligned_dna

REPO = Path("/home/breno/I2CA/genomics")
CACHE_DIR = REPO / "results/cache/genotype_based_predictor/raw_sequence_onehot"
SELECTED_SAMPLES_CSV = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage/selected_samples.csv")
CROP_SIZE = 32768

BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def one_hot(seq: str) -> np.ndarray:
    arr = np.zeros((4, len(seq)), dtype=np.float32)
    for i, ch in enumerate(seq):
        idx = BASE_IDX.get(ch)
        if idx is not None:
            arr[idx, i] = 1.0
    return arr


def center_crop(seq: str, size: int = CROP_SIZE) -> str:
    """Same central-crop convention as scripts/paper_figures/_common.py:crop_interval."""
    full = len(seq)
    if full < size:
        raise ValueError(f"sequence length {full} shorter than crop size {size}")
    cs = max(0, full // 2 - size // 2)
    return seq[cs:cs + size]


def _cohort_sample_ids() -> list[str]:
    ids = []
    with open(SELECTED_SAMPLES_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.append(row["SampleID"])
    return ids


def _config_restricted_to_cohort(config_path: Path) -> Path:
    """Write a temp config copy with sample_ids restricted to the 1072-individual cohort,
    so export_aligned_dna does not waste time aligning the other ~2130 1000G individuals
    that are never used by this dataset."""
    cfg = yaml.safe_load(config_path.read_text())
    cfg["dataset_input"]["sample_ids"] = _cohort_sample_ids()
    tmp_path = CACHE_DIR / f"_cohort_restricted_{config_path.stem}.yaml"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    tmp_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return tmp_path


def build_cache(config_path: Path, gene: str, force: bool = False) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = CACHE_DIR / f"{gene.upper()}.npz"
    if out_path.exists() and not force:
        print(f"{gene:<9} SKIP: {out_path} already present")
        return out_path

    restricted_config_path = _config_restricted_to_cohort(config_path)
    tsv_path = CACHE_DIR / f"_{gene.upper()}_aligned.tsv"
    try:
        export_aligned_dna(
            config_path=restricted_config_path,
            output_path=tsv_path,
            gene_name=gene,
            sample_limit=None,
            all_samples=True,
            center_window_size=32768,
        )
    except RuntimeError as exc:
        if "not found in the header" not in str(exc):
            raise
        # This gene's chromosome VCF drops one of the 1072 cohort individuals (a
        # per-chromosome QC exclusion, not a bug in the restriction list) -- bcftools
        # refuses a sample_ids list it can't fully satisfy rather than skipping the
        # missing one, so fall back to the unrestricted (all-3202) export just for this
        # gene instead of guessing which id to drop.
        print(f"{gene:<9} WARN: cohort-restricted sample list not fully present in this "
              f"gene's VCF ({exc}); retrying without the restriction (slower)")
        export_aligned_dna(
            config_path=config_path,
            output_path=tsv_path,
            gene_name=gene,
            sample_limit=None,
            all_samples=True,
            center_window_size=32768,
        )
    finally:
        restricted_config_path.unlink(missing_ok=True)

    arrays = {}
    raw_lengths = set()
    with open(tsv_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "sample_id":
                continue
            sample_id, h1_seq, h2_seq = row[0], row[1], row[2]
            if sample_id == "REF":
                continue
            raw_lengths.add(len(h1_seq))
            raw_lengths.add(len(h2_seq))
            arrays[f"{sample_id}__H1"] = one_hot(center_crop(h1_seq))
            arrays[f"{sample_id}__H2"] = one_hot(center_crop(h2_seq))

    print(f"{gene:<9} raw expanded alignment length(s) {raw_lengths} -> cropped to {CROP_SIZE}")
    if not arrays:
        raise SystemExit(f"ABORT: {gene} produced no sample alignments from {tsv_path}")

    np.savez_compressed(out_path, **arrays)
    tsv_path.unlink()
    n_samples = len(arrays) // 2
    print(f"{gene:<9} DONE: {n_samples} individuals -> {out_path} "
          f"({out_path.stat().st_size / 1024**2:.1f} MB)")
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config_path", type=str)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    build_cache(Path(args.config_path).resolve(), args.gene, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
