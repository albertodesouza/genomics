"""Generate (and cache) AlphaGenome RNA-seq predictions for individuals in this
notebook's bundled dataset, on their own diploid consensus sequence per gene,
realigned back onto reference coordinates.

Resumable: each (gene, sample, haplotype) result is checked for an existing,
already-saved file before doing any API work, so a full run can be safely
re-launched after being interrupted (rate limits, crashes, etc.) without redoing
completed work or re-spending API quota on it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

import numpy as np

from . import consensus as consensus_mod
from . import realignment


def predict_haplotype(
    client,
    *,
    gene: str,
    gene_row: dict,
    sample_id: str,
    haplotype: str,
    requested_outputs: list,
    ontology_terms: list,
    ref_window_fa: Path,
    gene_vcf_path: Path,
    consensus_cache_dir: Path,
):
    """Builds the consensus sequence, predicts on it, and realigns the RNA-seq track
    back onto reference coordinates (max-collapse over insertions, zero-fill over
    deletions -- see `realignment.realign_consensus_values`). Returns an AlphaGenome
    output object shaped like a direct reference prediction."""
    import dataclasses

    from alphagenome.data import genome, track_data

    window_width = gene_row["end"] - gene_row["start"]
    interval = genome.Interval(gene_row["chrom"], gene_row["start"], gene_row["end"])

    sequence, vcf_cons_path = consensus_mod.build_haplotype_fasta(
        sample_id=sample_id,
        gene=gene,
        window_width=window_width,
        ref_window_fa=ref_window_fa,
        gene_vcf_path=gene_vcf_path,
        haplotype=haplotype,
        cache_dir=consensus_cache_dir,
    )
    if len(sequence) != window_width:
        raise ValueError(
            f"{sample_id} {gene} {haplotype}: sequence length {len(sequence)} != "
            f"window width {window_width}"
        )

    output = client.predict_sequence(
        sequence=sequence,
        requested_outputs=requested_outputs,
        ontology_terms=ontology_terms,
        interval=interval,
    )
    if output.rna_seq is not None:
        variants = realignment.load_haplotype_variants(vcf_cons_path, haplotype)
        realigned_values = realignment.realign_consensus_values(
            output.rna_seq.values, gene_row["start"], window_width, variants
        )
        realigned_rna_seq = track_data.TrackData(
            values=realigned_values,
            metadata=output.rna_seq.metadata,
            resolution=output.rna_seq.resolution,
            interval=interval,
            uns=output.rna_seq.uns,
        )
        output = dataclasses.replace(output, rna_seq=realigned_rna_seq)
    return output


def prediction_paths(output_dir: Path, gene: str, sample_id: str, haplotype: str) -> tuple:
    """(npz_path, meta_path) for one (gene, sample, haplotype)'s cached prediction --
    shared naming convention with `import_precomputed.py`, so imported and
    API-generated results land in and are found at the same paths."""
    gene_dir = output_dir / gene
    gene_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{sample_id}.{haplotype}"
    return gene_dir / f"{stem}.npz", gene_dir / f"{stem}.json"


def save_rna_seq(track, npz_path: Path, meta_path: Path) -> None:
    np.savez_compressed(npz_path, values=track.values.astype(np.float32))
    metadata_records = track.metadata.to_dict(orient="records") if hasattr(track.metadata, "to_dict") else None
    with open(meta_path, "w") as f:
        json.dump(
            {
                "chrom": track.interval.chromosome,
                "start": track.interval.start,
                "end": track.interval.end,
                "resolution": track.resolution,
                "metadata": metadata_records,
            },
            f,
        )


def load_rna_seq(npz_path: Path, meta_path: Path):
    """Inverse of `save_rna_seq` -- reconstructs an AlphaGenome `TrackData`."""
    import pandas as pd
    from alphagenome.data import genome, track_data

    values = np.load(npz_path)["values"]
    with open(meta_path) as f:
        meta = json.load(f)
    interval = genome.Interval(meta["chrom"], meta["start"], meta["end"])
    metadata = pd.DataFrame(meta["metadata"]) if meta["metadata"] is not None else None
    return track_data.TrackData(values=values, metadata=metadata, resolution=meta["resolution"], interval=interval)


def predict_reference_interval(client, *, gene: str, interval, requested_outputs: list, ontology_terms: list, cache_dir: Path):
    """Cached `client.predict_interval(...)` on the unmodified reference sequence --
    used by the notebook's reference-only overview plots (one call per gene, no
    per-individual consensus). Keyed on gene + interval bounds + ontology terms, so a
    top-to-bottom notebook re-run reuses the cache instead of re-spending API calls on
    predictions that can't change."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = f"{gene}.{interval.start}-{interval.end}.{'_'.join(ontology_terms)}"
    npz_path, meta_path = cache_dir / f"{key}.npz", cache_dir / f"{key}.json"

    if npz_path.exists() and meta_path.exists():
        return load_rna_seq(npz_path, meta_path)

    output = client.predict_interval(
        interval=interval,
        requested_outputs=requested_outputs,
        ontology_terms=ontology_terms,
    )
    save_rna_seq(output.rna_seq, npz_path, meta_path)
    return output.rna_seq

def predict_variant_effect(client, *, gene: str, interval, variant, requested_outputs: list, ontology_terms: list):
    output = client.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=requested_outputs,
        ontology_terms=ontology_terms,
    )
    return output.reference.rna_seq, output.alternate.rna_seq

def generate_predictions(
    client,
    *,
    genes: Iterable[str],
    gene_rows: Dict[str, dict],
    gene_paths_fn: Callable[[str], tuple],
    sample_ids: Iterable[str],
    ontology_terms: list,
    output_dir: Path,
    consensus_cache_dir: Path,
    requested_outputs: Optional[list] = None,
    haplotypes: Iterable[str] = ("H1", "H2"),
    resume: bool = True,
    on_progress: Optional[Callable[[int, int, str, str, str, str], None]] = None,
) -> List[dict]:
    """Generates and caches realigned RNA-seq predictions for every
    (gene, sample_id, haplotype) combination. Returns a list of result records:
    {"gene", "sample_id", "haplotype", "npz_path", "meta_path", "status"}, where
    status is "cached" (skipped, already on disk) or "generated" (fresh API call).

    Safe to interrupt and re-run: existing (gene, sample, haplotype) outputs are
    skipped when `resume=True` (the default), so this can be called repeatedly to
    incrementally grow coverage (e.g. widen `sample_ids` or `genes` over time)
    without redoing completed work.
    """
    from alphagenome.models import dna_client

    if requested_outputs is None:
        requested_outputs = [dna_client.OutputType.RNA_SEQ]

    output_dir = Path(output_dir)
    genes = list(genes)
    sample_ids = list(sample_ids)
    haplotypes = list(haplotypes)
    total = len(genes) * len(sample_ids) * len(haplotypes)

    results = []
    done = 0
    for gene in genes:
        gene_row = gene_rows[gene]
        ref_window_fa, gene_vcf_path = gene_paths_fn(gene)
        for sample_id in sample_ids:
            for haplotype in haplotypes:
                npz_path, meta_path = prediction_paths(output_dir, gene, sample_id, haplotype)
                done += 1
                if resume and npz_path.exists() and meta_path.exists():
                    status = "cached"
                else:
                    output = predict_haplotype(
                        client,
                        gene=gene,
                        gene_row=gene_row,
                        sample_id=sample_id,
                        haplotype=haplotype,
                        requested_outputs=requested_outputs,
                        ontology_terms=ontology_terms,
                        ref_window_fa=ref_window_fa,
                        gene_vcf_path=gene_vcf_path,
                        consensus_cache_dir=consensus_cache_dir,
                    )
                    save_rna_seq(output.rna_seq, npz_path, meta_path)
                    status = "generated"
                results.append({
                    "gene": gene, "sample_id": sample_id, "haplotype": haplotype,
                    "npz_path": npz_path, "meta_path": meta_path, "status": status,
                })
                if on_progress:
                    on_progress(done, total, gene, sample_id, haplotype, status)
    return results
