"""GENCODE GTF loading and gene/exon/TSS coordinate helpers."""

import json

import pandas as pd
import requests
from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils

from .haplotype import haplotype_local_idx

GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)


def load_gtf(cache_dir):
    # Same cache location/URL as notebooks/utils/annotations.py's load_gtf, so this reuses
    # that notebook's ~tens-of-MB download if it already ran -- reimplemented here (rather than
    # importing notebooks.utils) to keep this notebook's dependency on the rest of the repo
    # scoped to genomics.* plus this notebook's own isolated package.
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / "gencode.v46.annotation.gtf.gz.feather"
    if not cache_path.exists():
        r = requests.get(GTF_URL, timeout=120)
        r.raise_for_status()
        cache_path.write_bytes(r.content)
    return pd.read_feather(cache_path)


def build_transcript_extractors(gtf):
    """Returns (gtf_mane, gtf_protein_coding, transcript_extractor_mane,
    transcript_extractor_coding, gene_id_map) derived from a loaded GENCODE GTF dataframe."""
    gtf_mane = gene_annotation.filter_to_mane_select_transcript(gene_annotation.filter_protein_coding(gtf))
    gtf_protein_coding = gene_annotation.filter_protein_coding(gtf)
    transcript_extractor_mane = transcript_utils.TranscriptExtractor(gtf_mane)
    transcript_extractor_coding = transcript_utils.TranscriptExtractor(gtf_protein_coding)
    gene_id_map = gtf.drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
    return gtf_mane, gtf_protein_coding, transcript_extractor_mane, transcript_extractor_coding, gene_id_map


def build_tss_tables(gtf_mane, gtf_protein_coding):
    return gene_annotation.extract_tss(gtf_mane), gene_annotation.extract_tss(gtf_protein_coding)


def window_genomic_axis(dataset_dir, gene, window_center_size):
    '''Reference-coordinate window start/length, using the same centering formula
    `no_alignment`'s `_extract_center_raw` applies when building its training tensors -- exact
    for that model and approximate for the DITA models (bcftools_chain-based expanded axis).'''
    meta_path = dataset_dir / "references" / "windows" / gene / "window_metadata.json"
    meta = json.loads(meta_path.read_text())
    chrom = meta["chromosome"]
    start_1based = int(meta["start"])
    end_1based = int(meta["end"])
    full_ref_length = end_1based - start_1based + 1

    size = min(window_center_size, full_ref_length)
    center_offset = full_ref_length // 2
    crop_start = max(0, center_offset - size // 2)
    crop_end = min(full_ref_length, crop_start + size)
    if crop_end - crop_start < size:
        crop_start = max(0, crop_end - size)

    window_start_0based = (start_1based - 1) + crop_start  # alphagenome genome.Interval is 0-based half-open
    local_length = crop_end - crop_start
    return chrom, window_start_0based, local_length


def gene_exon_local_boxes(gene, chrom, window_start_0based, local_length, gene_id_map,
                           transcript_extractor_mane, transcript_extractor_coding):
    interval = genome.Interval(chrom, window_start_0based, window_start_0based + local_length)
    target_gene_id = gene_id_map.get(gene)

    transcripts = transcript_extractor_mane.extract(interval)
    matched = [t for t in transcripts if t.gene_id == target_gene_id] if target_gene_id else transcripts
    if not matched:
        # Fallback for genes without a MANE Select tag in this GENCODE release: take the
        # protein-coding transcript with the most exons for this gene within the window.
        transcripts = [t for t in transcript_extractor_coding.extract(interval) if t.gene_id == target_gene_id]
        matched = [max(transcripts, key=lambda t: len(t.exons))] if transcripts else []

    boxes, strand = [], None
    for t in matched:
        strand = t.strand
        for exon in t.exons:
            s = max(0, exon.start - window_start_0based)
            e = min(local_length, exon.end - window_start_0based)
            if e > s:
                boxes.append((s, e))
    return sorted(boxes), strand


def gene_exon_haplotype_local_boxes(dataset_dir, sample_id, gene, haplotype, gene_id_map,
                                     transcript_extractor_mane, transcript_extractor_coding):
    """Same exon boxes as `gene_exon_local_boxes`, but each boundary converted through
    `haplotype_local_idx` into this haplotype's own indel-corrected local coordinates -- for
    overlaying onto individual-specific (CAGE/knockout) plots, which are drawn in per-haplotype
    consensus-sequence coordinates rather than the raw reference window axis."""
    meta_path = dataset_dir / "references" / "windows" / gene / "window_metadata.json"
    meta = json.loads(meta_path.read_text())
    chrom = meta["chromosome"]
    start_1based = int(meta["start"])
    full_ref_length = int(meta["end"]) - start_1based + 1
    window_start_0based = start_1based - 1

    ref_boxes, strand = gene_exon_local_boxes(
        gene, chrom, window_start_0based, full_ref_length, gene_id_map,
        transcript_extractor_mane, transcript_extractor_coding,
    )
    hap_boxes = []
    for s, e in ref_boxes:
        start_pos_1based = window_start_0based + s + 1
        end_pos_1based = window_start_0based + e  # exclusive edge -> 1-based pos of last included base
        hap_s = haplotype_local_idx(dataset_dir, sample_id, gene, haplotype, start_1based, start_pos_1based)
        hap_e = haplotype_local_idx(dataset_dir, sample_id, gene, haplotype, start_1based, end_pos_1based) + 1
        hap_boxes.append((hap_s, hap_e))
    return sorted(hap_boxes), strand


def gene_boundary_positions(exon_boxes, strand):
    """(tss_local, tes_local) derived from the already-extracted exon boxes. Transcription runs
    left-to-right on the `+` strand (TSS at the leftmost exon start) and right-to-left on the `-`
    strand (TSS at the rightmost exon end), so which physical edge is "start" vs. "end" flips with
    strand even though the (s, e) boxes themselves are strand-agnostic."""
    if not exon_boxes:
        return None, None
    left = min(s for s, _ in exon_boxes)
    right = max(e for _, e in exon_boxes)
    return (right, left) if strand == "-" else (left, right)


def get_gene_tss(gene, tss_df_mane, tss_df_coding):
    '''Returns (chrom, tss_pos_0based, strand) from the MANE Select transcript, falling back to
    any protein-coding transcript's TSS if no MANE Select tag exists for this gene.'''
    rows = tss_df_mane[tss_df_mane["gene_name"] == gene]
    if rows.empty:
        rows = tss_df_coding[tss_df_coding["gene_name"] == gene]
    if rows.empty:
        raise ValueError(f"No TSS found for gene {gene!r}")
    row = rows.iloc[0]
    return row["Chromosome"], int(row["Start"]), row["Strand"]


def get_gene_start(gene, gtf_gene_rows):
    '''Returns (chrom, gene_start_pos_0based) -- the raw "gene" GTF feature's lower genomic
    coordinate, taken as-is regardless of strand. Deliberately NOT the same as get_gene_tss: for
    a minus-strand gene this point is the transcript's 3' end, not its TSS.'''
    rows = gtf_gene_rows[gtf_gene_rows["gene_name"] == gene]
    if rows.empty:
        raise ValueError(f"No gene feature found for {gene!r}")
    row = rows.iloc[0]
    return row["Chromosome"], int(row["Start"])
