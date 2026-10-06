"""Realign an AlphaGenome track predicted on an individual's consensus sequence
back onto the reference genome's coordinate axis.

AlphaGenome's client provides no multi-variant / consensus-sequence prediction API
(only single-variant `predict_variant`) and no coordinate-lifting utility for
indels -- the caller is responsible for reconciling positions whenever the
sequence fed to `predict_sequence` differs in length from the reference at any
point. This module follows AlphaGenome's own documented recipe for reconciling
ALT-length variants back onto the reference: an inserted stretch (extra bases in
the individual's sequence, absent from the reference) collapses to the maximum
track value observed across the inserted span, assigned to the single
reference slot the insertion occurred at; a deleted stretch (reference bases
absent from the individual's sequence) has no predicted signal and is zero-filled.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import List, Tuple

import numpy as np

HAPLOTYPE_ALLELE_INDEX = {"H1": 0, "H2": 1}


def load_haplotype_variants(vcf_path: Path, haplotype: str) -> List[Tuple[int, str, str]]:
    """Returns (pos_1based, ref, alt) for variants this haplotype carries the ALT of.

    `vcf_path` is expected to be a single-sample, biallelic-only VCF (the
    per-sample, per-window `<sample>.window.consensus_ready.vcf.gz` files this
    repo's dataset builder already produces as `bcftools consensus`'s own input),
    sorted ascending by position.
    """
    allele_index = HAPLOTYPE_ALLELE_INDEX[haplotype]
    variants: List[Tuple[int, str, str]] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos, ref, alt, info, gt = int(fields[1]), fields[3], fields[4], fields[7], fields[-1]
            if "," in alt:
                continue  # multiallelic sites are not expected in this per-sample VCF
            alleles = gt.replace("/", "|").split("|")
            if alleles[allele_index] != "1":
                continue
            if alt.startswith("<"):
                # Symbolic ALT (e.g. "<DEL>", left in by the dataset builder's consensus-ready
                # filter): bcftools consensus resolves its true extent from INFO/END, not from
                # the literal ALT string -- using len(alt) here (e.g. len("<DEL>") == 5) would
                # misread a real multi-kb deletion as a few-base insertion and desynchronize
                # ref_idx/cons_idx for the remainder of the window. Reproduce bcftools's own
                # anchor-base-plus-deleted-span representation instead.
                end = next(
                    (int(kv[len("END="):]) for kv in info.split(";") if kv.startswith("END=")),
                    None,
                )
                if end is None or end < pos:
                    continue  # can't resolve a safe span; skip rather than misalign
                ref = "N" * (end - pos + 1)
                alt = "N"
            variants.append((pos, ref, alt))
    variants.sort(key=lambda v: v[0])
    return variants


def realign_consensus_values(
    consensus_values: np.ndarray,
    window_start: int,
    window_width: int,
    variants: List[Tuple[int, str, str]],
) -> np.ndarray:
    """Realigns a `(positions, num_tracks)` array from consensus-sequence space
    to reference-genome space.

    `window_start` must use the same 1-based convention as VCF `POS` and as
    `window_metadata.json`'s own `start` field (so index 0 of the output lines
    up with whatever `genome.Interval` the rest of the notebook already builds
    from that same `start` value).
    """
    if consensus_values.ndim != 2:
        raise ValueError(f"expected a (positions, num_tracks) array, got shape {consensus_values.shape}")

    num_tracks = consensus_values.shape[1]
    realigned = np.zeros((window_width, num_tracks), dtype=consensus_values.dtype)
    consensus_length = consensus_values.shape[0]

    ref_idx = 0  # next reference-relative index (0-based) to fill
    cons_idx = 0  # next consensus-sequence index (0-based) to read from

    for vcf_pos, ref, alt in variants:
        ref_local = vcf_pos - window_start
        if ref_local < ref_idx or ref_local >= window_width:
            continue  # out-of-order/overlapping variant, or beyond the window

        # Copy the untouched stretch before this variant 1:1.
        gap = min(ref_local - ref_idx, consensus_length - cons_idx)
        if gap > 0:
            realigned[ref_idx:ref_idx + gap] = consensus_values[cons_idx:cons_idx + gap]
            ref_idx += gap
            cons_idx += gap
        if cons_idx >= consensus_length:
            break

        ref_len, alt_len = len(ref), len(alt)

        if alt_len >= ref_len:
            # SNV (ref_len == alt_len) or insertion (alt_len > ref_len): the anchor plus any
            # inserted bases collapse to `ref_len` reference slot(s) via max over the ALT span.
            span_end = min(cons_idx + alt_len, consensus_length)
            collapsed = consensus_values[cons_idx:span_end].max(axis=0)
            fill_len = min(ref_len, window_width - ref_idx)
            realigned[ref_idx:ref_idx + fill_len] = collapsed
        else:
            # Deletion: the first `alt_len` ref slot(s) (the anchor) keep the real consensus
            # signal; the remaining `ref_len - alt_len` reference-only positions have no
            # consensus counterpart and are zero-filled.
            anchor_values = consensus_values[cons_idx:min(cons_idx + alt_len, consensus_length)]
            anchor_fill_len = min(alt_len, window_width - ref_idx)
            realigned[ref_idx:ref_idx + anchor_fill_len] = anchor_values[:anchor_fill_len]
            # realigned[ref_idx + alt_len : ref_idx + ref_len] is already 0 from initialization.

        ref_idx += ref_len
        cons_idx += alt_len

    # Copy the remaining untouched tail 1:1.
    remaining = min(window_width - ref_idx, consensus_length - cons_idx)
    if remaining > 0:
        realigned[ref_idx:ref_idx + remaining] = consensus_values[cons_idx:cons_idx + remaining]

    return realigned
