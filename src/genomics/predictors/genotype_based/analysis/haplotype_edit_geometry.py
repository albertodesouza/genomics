"""Pure-Python helpers for interactively editing (overwrite/scramble) a haplotype's consensus
sequence at an arbitrary user-selected region, and for translating between reference-genome
coordinates and a haplotype's own indel-shifted FASTA coordinates.

Promoted from ad-hoc cells in notebooks/genotype_cnn_alignment_deeplift.ipynb (Section 9), which
only ever scrambled a fixed-size window centered on one of three hardcoded knockout locations.
Here the region is arbitrary (a user's mouse-drag selection) and "overwrite with a nucleotide" is
a new operation the notebook didn't need.

No torch/pandas/AlphaGenome dependency -- this module only touches plain strings and VCF text, so
it stays unit-testable without a dataset fixture or a trained model.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Union

VALID_OVERWRITE_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class EditResult:
    modified_sequence: str
    original_segment: str
    edited_segment: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive


def haplotype_indel_events(vcf_path: Union[str, Path], haplotype: str) -> List[Tuple[int, int]]:
    """(pos_1based, len(alt) - len(ref)) for every indel this haplotype's phased genotype carries
    in the given per-gene VCF, sorted by position. Only positions where this haplotype's own
    allele changes length are returned (a single-base SNV contributes nothing).

    ``haplotype`` must be "H1" or "H2" -- these map to the first/second allele of the phased GT.
    """
    if haplotype not in ("H1", "H2"):
        raise ValueError(f"haplotype deve ser H1 ou H2, recebeu {haplotype!r}")
    hap_index = 0 if haplotype == "H1" else 1

    events: List[Tuple[int, int]] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            pos = int(cols[1])
            ref = cols[3]
            alt_alleles = cols[4].split(",")
            gt = cols[9].split(":")[0]
            sep = "|" if "|" in gt else "/"
            alleles = gt.split(sep)
            if len(alleles) != 2:
                continue
            allele_str = alleles[hap_index]
            if allele_str in ("0", "."):
                continue
            try:
                allele_num = int(allele_str)
            except ValueError:
                continue
            if not (1 <= allele_num <= len(alt_alleles)):
                continue
            delta = len(alt_alleles[allele_num - 1]) - len(ref)
            if delta != 0:
                events.append((pos, delta))
    events.sort()
    return events


def haplotype_local_idx(events: List[Tuple[int, int]], start_1based: int, target_pos_1based: int) -> int:
    """Reference position ``target_pos_1based`` -> 0-based local index into this haplotype's own
    fixed FASTA, correcting for every indel this haplotype carries between the window start and
    the target position."""
    drift = sum(delta for pos, delta in events if pos < target_pos_1based)
    return (target_pos_1based - start_1based) + drift


def genomic_pos_from_haplotype_local_idx(events: List[Tuple[int, int]], start_1based: int, local_idx: int) -> int:
    """Approximate inverse of :func:`haplotype_local_idx`: given a 0-based local FASTA index,
    return the reference genomic position (1-based) it corresponds to.

    Away from any indel span this is an exact inverse. Inside a deletion span the forward mapping
    isn't injective (no local base exists for those reference positions at all), so this returns
    the reference position reached by walking forward from the window start along the haplotype's
    own coordinate drift, i.e. the smallest genomic position whose forward-mapped local index
    would be >= ``local_idx``.
    """
    drift = 0
    prev_pos = start_1based
    for event_pos, delta in events:
        span_local_end = (event_pos - start_1based) + drift
        if local_idx < span_local_end:
            return start_1based + (local_idx - drift)
        drift += delta
        prev_pos = event_pos
    return start_1based + (local_idx - drift)


def apply_overwrite(sequence: str, start: int, end: int, base: str) -> EditResult:
    """Replace ``sequence[start:end]`` with ``base`` repeated ``end - start`` times.

    ``base`` must be a single uppercase nucleotide in {A, C, G, T} -- lowercase, ``N``, and ``X``
    are rejected explicitly since they are not real edits a user should be able to author here.
    """
    if base not in VALID_OVERWRITE_BASES:
        raise ValueError(f"base deve ser um de {sorted(VALID_OVERWRITE_BASES)}, recebeu {base!r}")
    if not (0 <= start < end <= len(sequence)):
        raise ValueError(f"regiao invalida start={start} end={end} para sequencia de tamanho {len(sequence)}")

    original_segment = sequence[start:end]
    edited_segment = base * (end - start)
    modified = sequence[:start] + edited_segment + sequence[end:]
    return EditResult(
        modified_sequence=modified,
        original_segment=original_segment,
        edited_segment=edited_segment,
        start=start,
        end=end,
    )


def apply_scramble_region(sequence: str, start: int, end: int, seed: int = 0) -> EditResult:
    """Shuffle ``sequence[start:end]`` in place. A fixed ``seed`` keeps the shuffle pattern
    reproducible across calls with the same region length."""
    if not (0 <= start < end <= len(sequence)):
        raise ValueError(f"regiao invalida start={start} end={end} para sequencia de tamanho {len(sequence)}")

    import numpy as np

    original_segment = list(sequence[start:end])
    rng = np.random.default_rng(seed)
    scrambled = original_segment.copy()
    rng.shuffle(scrambled)
    edited_segment = "".join(scrambled)
    modified = sequence[:start] + edited_segment + sequence[end:]
    return EditResult(
        modified_sequence=modified,
        original_segment="".join(original_segment),
        edited_segment=edited_segment,
        start=start,
        end=end,
    )


def apply_scramble(sequence: str, center_idx: int, window_size: int = 100, seed: int = 0) -> EditResult:
    """Back-compat wrapper matching the notebook's original ``apply_scramble`` signature: scrambles
    a ``window_size``-bp window centered on ``center_idx`` (clamped to stay in-bounds and always
    exactly ``window_size`` long, unless the sequence itself is shorter)."""
    half = window_size // 2
    start = max(0, center_idx - half)
    end = min(len(sequence), start + window_size)
    if end - start < window_size:
        start = max(0, end - window_size)
    return apply_scramble_region(sequence, start, end, seed=seed)
