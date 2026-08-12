"""bcftools_chain-aligned per-base sequence data for the interactive sequence-logo/letter plots.

Two views over the same underlying alignment:

- ``individual_aligned_bases``: one individual's actual bases at each position of a gene's
  bcftools_chain expanded axis, with ``X`` for positions where that individual has no aligned
  base (e.g. inside a deletion, or if their consensus data is unavailable at all).
- ``class_base_frequency``: per-position A/C/G/T/gap frequency across every individual in a
  pigmentation class, aligned the same way.

Both are built directly on ``BcftoolsChainMapper.get_haplotype_entry`` (the same primitive
``apps/aligned_dna_viewer.py``'s ``_read_chain_aligned_window_with_masks`` uses), rather than going
through ``alignment/export_aligned_dna.py``'s TSV-export path -- that path is a batch/offline CLI
convenience, not an on-demand per-request API, and would mean writing/parsing a TSV file for
every viewer request.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

BASES = ("A", "C", "G", "T")
GAP = "gap"


def _read_fasta_sequence(path: Path) -> str:
    lines: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines).upper()


def _aligned_bases_with_source_or_none(
    mapper, gene: str, sample_id: str, haplotype: str, start: int, end: int
) -> Optional[Tuple[str, List[Optional[int]]]]:
    """1-based inclusive [start, end] window of ``sample_id``'s haplotype on the gene's
    bcftools_chain expanded axis. Returns ``None`` (not a string of ``X``) when the individual has
    no consensus data at all for this gene/haplotype, so callers can distinguish "no data" from
    "real gap within otherwise-present data". Alongside the aligned bases, also returns each
    position's 0-based index into the haplotype's own (gap-free) FASTA -- ``None`` at any position
    with no aligned base -- so a UI selection made in this expanded/aligned coordinate space can be
    translated back into the raw FASTA coordinates an edit actually needs."""
    entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
    if entry is None:
        return None
    fasta_path = entry.get("fasta_path") or (
        mapper.consensus_dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.{haplotype}.window.fixed.fa"
    )
    fasta = _read_fasta_sequence(Path(fasta_path))
    expanded_to_source = {
        int(target) + 1: int(source)
        for source, target in zip(entry.get("copy_from_indices", []), entry.get("expanded_indices", []))
    }
    bases = []
    source_indices: List[Optional[int]] = []
    for pos in range(start, end + 1):
        source_idx = expanded_to_source.get(pos)
        if source_idx is None or source_idx < 0 or source_idx >= len(fasta):
            bases.append("X")
            source_indices.append(None)
        else:
            bases.append(fasta[source_idx])
            source_indices.append(source_idx)
    return "".join(bases), source_indices


def _aligned_bases_or_none(mapper, gene: str, sample_id: str, haplotype: str, start: int, end: int) -> Optional[str]:
    result = _aligned_bases_with_source_or_none(mapper, gene, sample_id, haplotype, start, end)
    return None if result is None else result[0]


def individual_aligned_bases(mapper, gene: str, sample_id: str, haplotype: str, start: int, end: int) -> str:
    """Point C: one individual's actual aligned bases over [start, end] (1-based, inclusive), with
    ``X`` for any position with no aligned base -- including the case where the individual has no
    consensus data at all for this gene/haplotype."""
    aligned = _aligned_bases_or_none(mapper, gene, sample_id, haplotype, start, end)
    if aligned is None:
        return "X" * max(end - start + 1, 0)
    return aligned


def individual_aligned_bases_with_source(
    mapper, gene: str, sample_id: str, haplotype: str, start: int, end: int
) -> Tuple[str, List[Optional[int]]]:
    """Same as :func:`individual_aligned_bases`, but also returns each rendered position's 0-based
    index into the haplotype's own raw FASTA (``None`` at gaps), so a client-side drag-selection
    over the rendered (expanded-axis) positions can be translated into the FASTA-local
    ``(start, end)`` range :func:`haplotype_edit_geometry.apply_overwrite`/``apply_scramble_region``
    expect."""
    result = _aligned_bases_with_source_or_none(mapper, gene, sample_id, haplotype, start, end)
    if result is None:
        length = max(end - start + 1, 0)
        return "X" * length, [None] * length
    return result


def class_base_frequency(
    mapper,
    gene: str,
    sample_ids: Sequence[str],
    haplotypes: Sequence[str],
    start: int,
    end: int,
) -> Dict[str, Union[np.ndarray, int]]:
    """Point B: per-position base frequency across every (sample_id, haplotype) combination,
    aligned via bcftools_chain.

    Individuals with no consensus data at all for this gene/haplotype are excluded from the
    denominator entirely (``n_used`` < ``n_requested``); individuals that do have data but carry a
    real gap (e.g. a deletion) at a given position count towards the "gap" frequency there.

    Returns a dict with keys "A", "C", "G", "T", "gap" (each a float array of length
    ``end - start + 1``, values summing to 1.0 per position when ``n_used > 0``), plus
    ``n_used`` and ``n_requested`` (= ``len(sample_ids) * len(haplotypes)``).
    """
    length = max(end - start + 1, 0)
    counts = {symbol: np.zeros(length, dtype=np.float64) for symbol in (*BASES, GAP)}
    pairs = [(sample_id, haplotype) for sample_id in sample_ids for haplotype in haplotypes]
    n_requested = len(pairs)
    n_used = 0

    # Fetching each (sample, haplotype)'s aligned bases is the expensive step (a first-time
    # bcftools_chain entry build shells out to `bcftools consensus`); a class spans many
    # individuals, so fetch concurrently and only do the (cheap) tallying serially.
    with ThreadPoolExecutor(max_workers=min(16, len(pairs)) or 1) as pool:
        aligned_results = pool.map(
            lambda pair: _aligned_bases_or_none(mapper, gene, pair[0], pair[1], start, end), pairs
        )
        for aligned in aligned_results:
            if aligned is None:
                continue
            n_used += 1
            for idx, base in enumerate(aligned):
                symbol = base if base in BASES else GAP
                counts[symbol][idx] += 1.0

    frequencies: Dict[str, Union[np.ndarray, int]] = {}
    if n_used > 0:
        for symbol, arr in counts.items():
            frequencies[symbol] = arr / n_used
    else:
        for symbol in counts:
            frequencies[symbol] = counts[symbol]
    frequencies["n_used"] = n_used
    frequencies["n_requested"] = n_requested
    return frequencies
