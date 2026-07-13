"""Imports precomputed AlphaGenome RNA-seq predictions from the pre-built
`1kG_high_coverage` dataset into this notebook's own prediction cache
(`notebooks/.cache/predictions/`), realigning them onto reference coordinates the same
way `predictions.py`'s own API-based pipeline does -- so imported and freshly-generated
results are interchangeable, and `predictions.generate_predictions(..., resume=True)`
treats an imported (gene, sample, haplotype) exactly like one it generated itself
(skips it).

This is the one piece of the "predictions" pipeline allowed to depend on the external,
pre-built dataset -- meant to be run once from `notebooks/prepare_data.ipynb`, not from
`simplified_notebook.ipynb`. It reuses that dataset's raw (unrealigned) per-haplotype
`predictions_H{1,2}/rna_seq.npz` and its `<sample>.window.consensus_ready.vcf.gz`
(the same VCF `bcftools consensus` built the predicted sequence from), so the results
are realigned with this repo's current (indel-aware, `<DEL>`-aware) logic regardless of
what logic -- if any -- was used when that dataset was originally built.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import numpy as np

from . import realignment
from .predictions import prediction_paths, save_rna_seq


def load_raw_prediction(dataset_dir: Path, sample_id: str, gene: str, haplotype: str) -> np.ndarray:
    pred_dir = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
    npz_path = pred_dir / "rna_seq.npz"
    if not npz_path.exists():
        raise FileNotFoundError(npz_path)
    return np.load(npz_path)["values"]


def import_haplotype(
    dataset_dir: Path,
    *,
    gene: str,
    gene_row: dict,
    sample_id: str,
    haplotype: str,
):
    """Loads one (gene, sample, haplotype)'s precomputed raw prediction + its
    consensus-ready VCF from `dataset_dir`, and returns a realigned `TrackData` --
    the metadata (track names/ontology terms) is carried over unchanged, only the
    per-base values are realigned."""
    import pandas as pd
    from alphagenome.data import genome, track_data

    window_width = gene_row["end"] - gene_row["start"]
    interval = genome.Interval(gene_row["chrom"], gene_row["start"], gene_row["end"])

    raw_values = load_raw_prediction(dataset_dir, sample_id, gene, haplotype)
    if raw_values.shape[0] != window_width:
        raise ValueError(
            f"{sample_id} {gene} {haplotype}: precomputed prediction length {raw_values.shape[0]} "
            f"!= window width {window_width} -- notebooks/data/genes.json's window for this gene "
            "doesn't match the dataset this prediction was computed from (see "
            "prepare_data.load_dataset_window)."
        )

    meta_path = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}" / "rna_seq_metadata.json"
    with open(meta_path) as f:
        metadata_records = json.load(f).get("metadata")
    metadata_df = pd.DataFrame(metadata_records) if metadata_records is not None else None

    vcf_cons_path = (
        dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.window.consensus_ready.vcf.gz"
    )
    if not vcf_cons_path.exists():
        raise FileNotFoundError(vcf_cons_path)
    variants = realignment.load_haplotype_variants(vcf_cons_path, haplotype)
    realigned_values = realignment.realign_consensus_values(raw_values, gene_row["start"], window_width, variants)

    return track_data.TrackData(values=realigned_values, metadata=metadata_df, resolution=1, interval=interval)


def import_precomputed_predictions(
    dataset_dir: Path,
    *,
    genes: Iterable[str],
    gene_rows: Dict[str, dict],
    sample_ids: Iterable[str],
    output_dir: Path,
    haplotypes: Iterable[str] = ("H1", "H2"),
    resume: bool = True,
    on_progress: Optional[Callable[[int, int, str, str, str, str], None]] = None,
) -> List[dict]:
    """Imports every (gene, sample, haplotype) combination it can find in
    `dataset_dir`, writing to the same cache layout `predictions.generate_predictions`
    uses. Returns result records like `generate_predictions` does, with status one of
    "cached" (already present), "imported", or "missing_source" (not found in
    `dataset_dir` -- e.g. that sample/gene wasn't part of the pre-built dataset).

    Resumable: like `generate_predictions`, existing outputs are skipped when
    `resume=True` (the default), so this is safe to re-run or interrupt.
    """
    output_dir = Path(output_dir)
    genes = list(genes)
    sample_ids = list(sample_ids)
    haplotypes = list(haplotypes)
    total = len(genes) * len(sample_ids) * len(haplotypes)

    results = []
    done = 0
    for gene in genes:
        gene_row = gene_rows[gene]
        for sample_id in sample_ids:
            for haplotype in haplotypes:
                npz_path, meta_path = prediction_paths(output_dir, gene, sample_id, haplotype)
                done += 1
                if resume and npz_path.exists() and meta_path.exists():
                    status = "cached"
                else:
                    try:
                        track = import_haplotype(
                            dataset_dir, gene=gene, gene_row=gene_row, sample_id=sample_id, haplotype=haplotype
                        )
                        save_rna_seq(track, npz_path, meta_path)
                        status = "imported"
                    except FileNotFoundError:
                        status = "missing_source"
                results.append({
                    "gene": gene, "sample_id": sample_id, "haplotype": haplotype,
                    "npz_path": npz_path, "meta_path": meta_path, "status": status,
                })
                if on_progress:
                    on_progress(done, total, gene, sample_id, haplotype, status)
    return results
