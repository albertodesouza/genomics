"""GTF-derived gene geometry and per-individual RNA-seq feature extraction, shared across the
`0N_*.ipynb` notebooks. Everything here is parameterized explicitly (cache dir, ontology terms,
gtf, ...) rather than closing over notebook globals, matching `predictions.py`/`consensus.py`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def gene_strand(gtf: pd.DataFrame, gene: str) -> str:
    from alphagenome.data import gene_annotation

    gtf_mane = gene_annotation.filter_to_mane_select_transcript(gene_annotation.filter_protein_coding(gtf))
    return gtf_mane.loc[
        (gtf_mane["gene_name"] == gene) & (gtf_mane["Feature"] == "transcript"), "Strand"
    ].iloc[0]


def compute_gene_zoom_slices(gtf: pd.DataFrame, genes: List[str], full_window: int) -> Dict[str, Tuple[int, int]]:
    """Window-relative (lo, hi) slice for each gene's own GTF **gene** feature span (every
    transcript/isoform), centered within a `full_window`-wide prediction window -- e.g. the
    500kb AlphaGenome prediction window these genes' windows were built with."""
    center_idx = full_window // 2
    gene_zoom_slices = {}
    for gene in genes:
        gene_span = gtf.loc[(gtf["Feature"] == "gene") & (gtf["gene_name"] == gene), ["Start", "End"]].iloc[0]
        zoom_width = int(gene_span["End"] - gene_span["Start"])
        half = zoom_width // 2
        gene_zoom_slices[gene] = (center_idx - half, center_idx + half)
    return gene_zoom_slices


def _mane_exons(gtf: pd.DataFrame, gene: str) -> pd.DataFrame:
    from alphagenome.data import gene_annotation

    gtf_protein_coding = gene_annotation.filter_protein_coding(gtf)
    gtf_mane_transcripts = gene_annotation.filter_to_mane_select_transcript(gtf_protein_coding)
    mane_tx_id = gtf_mane_transcripts.loc[
        (gtf_mane_transcripts["gene_name"] == gene) & (gtf_mane_transcripts["Feature"] == "transcript"),
        "transcript_id",
    ].iloc[0]
    return gtf_protein_coding.loc[
        (gtf_protein_coding["transcript_id"] == mane_tx_id) & (gtf_protein_coding["Feature"] == "exon")
    ]


def mane_exon_mask(gtf: pd.DataFrame, gene: str, window_start: int, window_size: int) -> np.ndarray:
    """Boolean mask, True at positions (relative to the gene's prediction window) that fall
    inside an exon of the gene's MANE Select transcript -- everything else (introns,
    flanking sequence) is False."""
    exons = _mane_exons(gtf, gene)
    mask = np.zeros(window_size, dtype=bool)
    for _, exon in exons.iterrows():
        lo = max(0, int(exon["Start"] - window_start))
        hi = min(window_size, int(exon["End"] - window_start))
        mask[lo:hi] = True
    return mask


def mane_exon_indices(gtf: pd.DataFrame, gene: str, window_start: int, strand: str) -> Tuple[np.ndarray, List[int]]:
    """Window-relative indices of every base inside a MANE Select transcript exon, ordered
    5'->3' along the spliced mRNA -- so consecutive entries are adjacent on the mature
    transcript even though they jump across (removed) introns in genomic coordinates. Also
    returns the cumulative exon-boundary positions in that concatenated coordinate, for
    marking splice junctions on a plot."""
    exons = _mane_exons(gtf, gene).sort_values("Start", ascending=(strand == "+"))

    chunks, boundaries, cum = [], [], 0
    for _, exon in exons.iterrows():
        lo = int(exon["Start"] - window_start)
        hi = int(exon["End"] - window_start)
        chunk = np.arange(lo, hi)
        if strand == "-":
            chunk = chunk[::-1]
        chunks.append(chunk)
        cum += len(chunk)
        boundaries.append(cum)
    return np.concatenate(chunks), boundaries[:-1]


def filter_track_ontologies(track, ontology_terms: List[str]):
    """Restricts a TrackData to just the given ontology_curie terms -- a post-hoc column
    filter, since the on-disk cache already has all 3 terms baked in from generation time."""
    mask = track.metadata["ontology_curie"].isin(ontology_terms).to_numpy()
    return track.filter_tracks(mask)


def class_mean_signal(predictions_cache_dir: Path, gene: str, sample_ids: List[str], strand: str) -> np.ndarray:
    """Mean, across every individual's cached (H1+H2, all-ontology-terms-summed) RNA-seq
    track, restricted to `strand`. Reads only from `predictions_cache_dir`."""
    from .predictions import load_combined_rna_seq

    total = None
    for sample_id in sample_ids:
        track = load_combined_rna_seq(predictions_cache_dir, gene, sample_id)
        track = track.filter_to_positive_strand() if strand == "+" else track.filter_to_negative_strand()
        signal = track.values.sum(axis=1)
        total = signal if total is None else total + signal
    return total / len(sample_ids)


def individual_exon_means(
    predictions_cache_dir: Path,
    gene: str,
    sample_ids: List[str],
    strand: str,
    mask: np.ndarray,
    ontology_terms: List[str],
) -> np.ndarray:
    """Per-individual (H1+H2, restricted to `ontology_terms`) mean RNA-seq signal, restricted
    to the gene's MANE Select transcript exons -- one scalar per individual."""
    from .predictions import load_combined_rna_seq

    means = np.empty(len(sample_ids))
    for i, sample_id in enumerate(sample_ids):
        track = load_combined_rna_seq(predictions_cache_dir, gene, sample_id)
        track = track.filter_to_positive_strand() if strand == "+" else track.filter_to_negative_strand()
        track = filter_track_ontologies(track, ontology_terms)
        signal = track.values.sum(axis=1)
        means[i] = signal[mask].mean()
    return means

def individual_exon_log_means(
    predictions_cache_dir: Path,
    gene: str,
    sample_ids: List[str],
    strand: str,
    mask: np.ndarray,
    ontology_terms: List[str],
) -> np.ndarray:
    """Per-individual (H1+H2, restricted to `ontology_terms`) log(mean + 0.001) RNA-seq signal, restricted
    to the gene's MANE Select transcript exons -- one scalar per individual."""
    means = individual_exon_means(
        predictions_cache_dir,
        gene,
        sample_ids,
        strand,
        mask,
        ontology_terms 
    )
    return np.log(means + 0.001)

def individual_gene_means(
    predictions_cache_dir: Path,
    gene: str,
    sample_ids: List[str],
    strand: str,
    lo: int,
    hi: int,
    ontology_terms: List[str],
) -> np.ndarray:
    """Per-individual (H1+H2, restricted to `ontology_terms`) mean RNA-seq signal, over the
    gene's entire GTF span (`compute_gene_zoom_slices` -- every transcript/isoform, not just
    the MANE Select transcript's exons) -- one scalar per individual."""
    means = np.empty(len(sample_ids))
    for i, sample_id in enumerate(sample_ids):
        signal = individual_gene_span_signal(predictions_cache_dir, gene, sample_id, strand, lo, hi, ontology_terms)
        means[i] = signal.mean()
    return means


def individual_gene_log_means(
    predictions_cache_dir: Path,
    gene: str,
    sample_ids: List[str],
    strand: str,
    lo: int,
    hi: int,
    ontology_terms: List[str],
) -> np.ndarray:
    """Per-individual (H1+H2, restricted to `ontology_terms`) log(mean + 0.001) RNA-seq signal,
    over the gene's entire GTF span -- one scalar per individual."""
    means = individual_gene_means(predictions_cache_dir, gene, sample_ids, strand, lo, hi, ontology_terms)
    return np.log(means + 0.001)


def individual_gene_span_signal(
    predictions_cache_dir: Path,
    gene: str,
    sample_id: str,
    strand: str,
    lo: int,
    hi: int,
    ontology_terms: List[str],
) -> np.ndarray:
    """Combined (H1+H2, restricted to `ontology_terms`) RNA-seq signal for one individual,
    cropped to the gene's own GTF span (`compute_gene_zoom_slices` -- every transcript/isoform,
    not just the MANE Select transcript's exons)."""
    from .predictions import load_combined_rna_seq

    track = load_combined_rna_seq(predictions_cache_dir, gene, sample_id)
    track = track.filter_to_positive_strand() if strand == "+" else track.filter_to_negative_strand()
    track = filter_track_ontologies(track, ontology_terms)
    return track.values.sum(axis=1)[lo:hi].astype(np.float32)
