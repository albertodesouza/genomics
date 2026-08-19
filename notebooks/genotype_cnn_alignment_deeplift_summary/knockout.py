"""Promoter-location knockout experiments: scramble a candidate promoter site on both haplotypes
of a real individual's consensus sequence, re-predict with AlphaGenome, and measure the shift in
the CNN's log-odds(strong/weak pigmentation)."""

import json

import numpy as np
import pandas as pd
import torch
from alphagenome.data import genome
from alphagenome.models import dna_client

from genomics.predictors.genotype_based.data.normalization import apply_normalization

from .annotations import get_gene_tss
from .cage import find_cage_promoter_individual
from .haplotype import haplotype_local_idx, load_haplotype_fasta


def load_raw_prediction(dataset_dir, sample_id, gene, haplotype):
    d = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
    values = np.load(d / "rna_seq.npz")["values"]
    meta = json.loads((d / "rna_seq_metadata.json").read_text())["metadata"]
    return values, meta


def apply_scramble(sequence, center_idx, window_size=100, seed=0):
    '''Shuffles the window_size-bp window centered on center_idx (clamped to stay in-bounds and
    always exactly window_size long). A fixed seed keeps results reproducible across reruns --
    the shuffle *pattern* is the same every call, but the bases it's applied to differ per
    gene/haplotype/method, so the actual scrambled sequence still differs each time.'''
    half = window_size // 2
    start = max(0, center_idx - half)
    end = min(len(sequence), start + window_size)
    if end - start < window_size:
        start = max(0, end - window_size)
    original_segment = list(sequence[start:end])
    rng = np.random.default_rng(seed)
    scrambled = original_segment.copy()
    rng.shuffle(scrambled)
    scrambled_segment = "".join(scrambled)
    modified = sequence[:start] + scrambled_segment + sequence[end:]
    return modified, "".join(original_segment), scrambled_segment, start, end


def reorder_to_canonical(values, metadata_df, canonical_order):
    lookup = {(row["ontology_curie"], row["strand"]): i for i, row in metadata_df.iterrows()}
    col_idx = [lookup[key] for key in canonical_order]
    return values[:, col_idx]


def predict_modified_sequence(ctx, cache_key, sequence, interval):
    cache_path = ctx.cache_dir / f"seq_{cache_key}.npz"
    meta_cache_path = ctx.cache_dir / f"seq_{cache_key}_meta.json"
    if cache_path.exists() and meta_cache_path.exists():
        return np.load(cache_path)["values"], json.loads(meta_cache_path.read_text())

    output = ctx.ag_client.predict_sequence(
        sequence, organism=ctx.organism, requested_outputs=[dna_client.OutputType.RNA_SEQ],
        ontology_terms=ctx.ontology_terms, interval=interval,
    )
    values = output.rna_seq.values
    meta_records = output.rna_seq.metadata[["ontology_curie", "strand"]].to_dict("records")
    np.savez_compressed(cache_path, values=values)
    meta_cache_path.write_text(json.dumps(meta_records))
    return values, meta_records


def build_individual_tensor(ctx, sample_id, overrides=None):
    '''overrides: {(gene, haplotype): (raw_array, track_meta)} substituted in place of the
    on-disk cached prediction for that gene/haplotype; every other gene/haplotype uses the
    individual's real, unmodified prediction.'''
    overrides = overrides or {}
    h1_rows, h2_rows = [], []
    for gene in ctx.genes:
        for haplotype, rows in (("H1", h1_rows), ("H2", h2_rows)):
            if (gene, haplotype) in overrides:
                array, meta = overrides[(gene, haplotype)]
            else:
                array, meta = load_raw_prediction(ctx.dataset_dir, sample_id, gene, haplotype)
            result = ctx.full_ds._process_window_haplotype_channels(
                sample_id, gene, haplotype, {"rna_seq": array}, {"rna_seq": meta},
            )
            if result is None:
                raise RuntimeError(f"could not build tensor for {sample_id}/{gene}/{haplotype}")
            signals, _masks = result
            rows.append(signals)
    stacked = np.stack([np.concatenate(h1_rows, axis=0), np.concatenate(h2_rows, axis=0)], axis=0)
    return apply_normalization(torch.from_numpy(stacked), ctx.normalization_params)


def run_cnn_logits(ctx, tensor):
    with torch.no_grad():
        return ctx.model(tensor.unsqueeze(0).float().to(ctx.device))[0]


def cnn_crop_bounds(full_ref_length, window_center_size):
    # Identical centering formula to ReferenceRealignMapper._center_crop_bounds /
    # annotations.window_genomic_axis: the actual window_center_size-bp slice the CNN receives
    # as input, taken from the middle of the full 524,288 bp AlphaGenome window.
    size = min(window_center_size, full_ref_length)
    center_offset = full_ref_length // 2
    crop_start = max(0, center_offset - size // 2)
    crop_end = min(full_ref_length, crop_start + size)
    if crop_end - crop_start < size:
        crop_start = max(0, crop_end - size)
    return crop_start, crop_end


def knockout_gene_experiment(ctx, sample_id, gene, method, scramble_window=100):
    window_meta = json.loads(
        (ctx.dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text()
    )
    chrom = window_meta["chromosome"]
    start_1based = int(window_meta["start"])
    full_ref_length = int(window_meta["end"]) - start_1based + 1
    interval = genome.Interval(chromosome=chrom, start=start_1based - 1, end=start_1based - 1 + full_ref_length)

    _, tss_pos_0based, _ = get_gene_tss(gene, ctx.tss_df_mane, ctx.tss_df_coding)
    tss_local_idx = tss_pos_0based - (start_1based - 1)

    location_labels = {
        "biology_tss": "TSS (GENCODE MANE Select)",
        "cage_melanocyte": "CAGE summit (dark/light melanocyte, individual-specific)",
    }
    if method not in location_labels:
        raise ValueError(f"unknown method: {method!r}")
    location_label = location_labels[method]

    hap_tracks = {}
    for haplotype in ("H1", "H2"):
        # biology_tss: a single reference position (GTF-derived, individual-invariant) converted
        # to this haplotype's own indel-corrected local coordinate. cage_melanocyte: found
        # directly on this haplotype's own re-predicted signal (find_cage_promoter_individual)
        # -- already haplotype-local, no reference round-trip needed or meaningful here, since
        # the summit itself can genuinely differ between H1 and H2, not just drift-shift.
        tss_hap_local_idx = haplotype_local_idx(
            ctx.dataset_dir, sample_id, gene, haplotype, start_1based, tss_pos_0based + 1
        )
        cage_hap = find_cage_promoter_individual(ctx, sample_id, gene, haplotype)
        marker_local_idx = {"biology_tss": tss_hap_local_idx, "cage_melanocyte": cage_hap["summit_local_idx"]}
        hap_target_local_idx = marker_local_idx[method]

        seq = load_haplotype_fasta(ctx.dataset_dir, sample_id, gene, haplotype)
        modified_seq, orig_segment, scrambled_segment, seg_start, seg_end = apply_scramble(
            seq, hap_target_local_idx, window_size=scramble_window,
        )
        orig_array, orig_meta = load_raw_prediction(ctx.dataset_dir, sample_id, gene, haplotype)
        canonical_order = [(m["ontology_curie"], m["strand"]) for m in orig_meta]

        cache_key = f"{sample_id}_{gene}_{haplotype}_{method}_{hap_target_local_idx}_scramble{scramble_window}"
        mod_values, mod_meta_records = predict_modified_sequence(ctx, cache_key, modified_seq, interval)
        mod_reordered = reorder_to_canonical(mod_values, pd.DataFrame(mod_meta_records), canonical_order)

        hap_tracks[haplotype] = {
            "orig_array": orig_array, "mod_array": mod_reordered, "meta": orig_meta,
            "scramble_start": seg_start, "scramble_end": seg_end,
            "orig_segment": orig_segment, "scrambled_segment": scrambled_segment,
            "marker_local_idx": marker_local_idx, "cage_peak_value": cage_hap["peak_value"],
        }

    baseline_tensor = build_individual_tensor(ctx, sample_id)
    overrides = {
        (gene, hap): (hap_tracks[hap]["mod_array"], hap_tracks[hap]["meta"]) for hap in ("H1", "H2")
    }
    perturbed_tensor = build_individual_tensor(ctx, sample_id, overrides=overrides)

    baseline_logits = run_cnn_logits(ctx, baseline_tensor)
    perturbed_logits = run_cnn_logits(ctx, perturbed_tensor)
    baseline_log_odds = (baseline_logits[ctx.strong_idx] - baseline_logits[ctx.weak_idx]).item()
    perturbed_log_odds = (perturbed_logits[ctx.strong_idx] - perturbed_logits[ctx.weak_idx]).item()

    return {
        "sample_id": sample_id, "gene": gene, "method": method, "location_label": location_label,
        "chrom": chrom, "tss_local_idx": tss_local_idx, "scramble_window": scramble_window,
        "hap_tracks": hap_tracks,
        "baseline_logits": baseline_logits.cpu().numpy(),
        "perturbed_logits": perturbed_logits.cpu().numpy(),
        "baseline_log_odds": baseline_log_odds,
        "perturbed_log_odds": perturbed_log_odds,
        "delta_log_odds": perturbed_log_odds - baseline_log_odds,
    }
