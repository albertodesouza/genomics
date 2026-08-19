"""Individual-specific melanocyte CAGE signal, re-predicted from real consensus haplotypes."""

import json

import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client

from .annotations import get_gene_start, get_gene_tss
from .haplotype import haplotype_local_idx, load_haplotype_fasta


def find_cage_curve_individual(ctx, sample_id, gene, haplotype):
    '''Full-window (524,288 bp) melanocyte CAGE signal from AlphaGenome's own re-prediction of
    this individual's real consensus sequence for this haplotype (predict_sequence on
    <sample>.H{1,2}.window.fixed.fa) -- not the reference genome. Cached per sample/gene/haplotype
    since, unlike a reference-based version, the signal genuinely differs by individual and
    haplotype. Returns the dark/light curves over the *entire* window so both the TSS-anchored and
    gene-start-anchored searches can search their own sub-window without any extra AlphaGenome
    calls.'''
    json_cache_path = ctx.cache_dir / f"cage_curve_{sample_id}_{gene}_{haplotype}.json"
    curve_cache_path = ctx.cache_dir / f"cage_curve_{sample_id}_{gene}_{haplotype}.npz"
    if json_cache_path.exists() and curve_cache_path.exists():
        result = json.loads(json_cache_path.read_text())
        curves = np.load(curve_cache_path)
        result["dark_curve"] = curves["dark"]
        result["light_curve"] = curves["light"]
        return result

    window_meta = json.loads(
        (ctx.dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text()
    )
    chrom = window_meta["chromosome"]
    start_1based = int(window_meta["start"])
    full_ref_length = int(window_meta["end"]) - start_1based + 1
    assert full_ref_length == dna_client.SEQUENCE_LENGTH_500KB

    interval = genome.Interval(chromosome=chrom, start=start_1based - 1, end=start_1based - 1 + full_ref_length)
    seq = load_haplotype_fasta(ctx.dataset_dir, sample_id, gene, haplotype)
    cage_output = ctx.ag_client.predict_sequence(
        seq, organism=ctx.organism, requested_outputs=[dna_client.OutputType.CAGE],
        ontology_terms=["CL:0002566", "CL:0002567"], interval=interval,
    )
    cage_values = cage_output.cage.values
    cage_meta = cage_output.cage.metadata
    dark_cols = [i for i, row in cage_meta.iterrows() if row["ontology_curie"] == "CL:0002566"]
    light_cols = [i for i, row in cage_meta.iterrows() if row["ontology_curie"] == "CL:0002567"]
    # Max across the 2 strand tracks within each subtype -- CAGE is stranded, and either strand's
    # initiation signal is informative for "is this an active promoter here."
    dark_curve = cage_values[:, dark_cols].max(axis=1)
    light_curve = cage_values[:, light_cols].max(axis=1)

    result = {"sample_id": sample_id, "gene": gene, "haplotype": haplotype, "chrom": chrom, "start_1based": start_1based}
    json_cache_path.write_text(json.dumps(result))
    np.savez_compressed(curve_cache_path, dark=dark_curve, light=light_curve)
    result["dark_curve"] = dark_curve
    result["light_curve"] = light_curve
    return result


def _cage_summit_in_curve(curve, anchor_local_idx, search_radius):
    combined = np.maximum(curve["dark_curve"], curve["light_curve"])
    lo = max(0, anchor_local_idx - search_radius)
    hi = min(len(combined), anchor_local_idx + search_radius)
    summit_local_idx = lo + int(np.argmax(combined[lo:hi]))
    return {
        "anchor_local_idx": anchor_local_idx, "summit_local_idx": summit_local_idx,
        "search_lo": lo, "search_hi": hi, "peak_value": float(combined[lo:hi].max()),
    }


def find_cage_promoter_individual(ctx, sample_id, gene, haplotype, search_radius=5000):
    '''CAGE summit within +/-5kb of the GENCODE TSS, searched on this individual's own
    haplotype-specific signal.'''
    curve = find_cage_curve_individual(ctx, sample_id, gene, haplotype)
    window_meta = json.loads(
        (ctx.dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text()
    )
    start_1based = int(window_meta["start"])
    _, tss_pos_0based, _ = get_gene_tss(gene, ctx.tss_df_mane, ctx.tss_df_coding)
    tss_hap_local_idx = haplotype_local_idx(ctx.dataset_dir, sample_id, gene, haplotype, start_1based, tss_pos_0based + 1)
    summit = _cage_summit_in_curve(curve, tss_hap_local_idx, search_radius)
    return {
        "sample_id": sample_id, "gene": gene, "haplotype": haplotype,
        "dark_curve": curve["dark_curve"], "light_curve": curve["light_curve"], **summit,
    }


def find_cage_promoter_no_tss_individual(ctx, sample_id, gene, haplotype, search_radius=5000):
    '''CAGE summit within +/-5kb of the gene's strand-naive genomic start (no biological TSS
    prior), searched on this individual's own haplotype-specific signal. Reuses the same cached
    full-window curve `find_cage_promoter_individual` already fetched -- no extra AlphaGenome
    calls, only the search window/anchor differs.'''
    curve = find_cage_curve_individual(ctx, sample_id, gene, haplotype)
    window_meta = json.loads(
        (ctx.dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text()
    )
    start_1based = int(window_meta["start"])
    _, gene_start_pos_0based = get_gene_start(gene, ctx.gtf_gene_rows)
    gene_start_hap_local_idx = haplotype_local_idx(
        ctx.dataset_dir, sample_id, gene, haplotype, start_1based, gene_start_pos_0based + 1
    )
    summit = _cage_summit_in_curve(curve, gene_start_hap_local_idx, search_radius)
    return {
        "sample_id": sample_id, "gene": gene, "haplotype": haplotype,
        "dark_curve": curve["dark_curve"], "light_curve": curve["light_curve"], **summit,
    }
