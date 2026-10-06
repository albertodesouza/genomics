"""Unit tests for alignment_frequency.py against a fake BcftoolsChainMapper-like collaborator.

BcftoolsChainMapper.get_haplotype_entry itself shells out to the ``bcftools`` binary (not
available in this environment), so these tests exercise alignment_frequency's own position-mapping
and frequency-tallying logic against a lightweight stand-in that returns entries shaped exactly
like the real mapper's (``copy_from_indices``/``expanded_indices`` pairs, ``fasta_path``), the same
contract ``apps/aligned_dna_viewer.py``'s ``_read_chain_aligned_window_with_masks`` and this
module's ``_aligned_bases_or_none`` both consume.
"""
from pathlib import Path

from genomics.predictors.genotype_based.analysis.alignment_frequency import (
    class_base_frequency,
    individual_aligned_bases,
    individual_aligned_bases_with_source,
)

GENE = "MC1R"


class FakeMapper:
    def __init__(self, entries):
        self.consensus_dataset_dir = Path("/unused")
        self._entries = entries

    def get_haplotype_entry(self, gene, sample_id, haplotype):
        return self._entries.get((gene, sample_id, haplotype))


def _write_fasta(path: Path, sequence: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">seq\n{sequence}\n")
    return path


def _full_match_entry(fasta_path: Path, length: int) -> dict:
    # 0-based source == 0-based target for every position (no indels).
    return {
        "fasta_path": str(fasta_path),
        "copy_from_indices": list(range(length)),
        "expanded_indices": list(range(length)),
    }


def test_individual_aligned_bases_full_match(tmp_path):
    fasta_path = _write_fasta(tmp_path / "ind1.fa", "ACGTA")
    mapper = FakeMapper({(GENE, "IND1", "H1"): _full_match_entry(fasta_path, 5)})

    assert individual_aligned_bases(mapper, GENE, "IND1", "H1", 1, 5) == "ACGTA"


def test_individual_aligned_bases_gap_at_deletion(tmp_path):
    fasta_path = _write_fasta(tmp_path / "ind2.fa", "ACTA")
    # Expanded (1-based) positions 1,2,4,5 map to fasta source indices 0,1,2,3; position 3 is a
    # deletion with no corresponding source base at all.
    entry = {
        "fasta_path": str(fasta_path),
        "copy_from_indices": [0, 1, 2, 3],
        "expanded_indices": [0, 1, 3, 4],  # 0-based targets -> 1-based positions 1,2,4,5
    }
    mapper = FakeMapper({(GENE, "IND2", "H1"): entry})
    assert individual_aligned_bases(mapper, GENE, "IND2", "H1", 1, 5) == "ACXTA"


def test_individual_aligned_bases_missing_individual_is_all_x(tmp_path):
    mapper = FakeMapper({})
    assert individual_aligned_bases(mapper, GENE, "GHOST", "H1", 1, 5) == "XXXXX"


def test_class_base_frequency_excludes_missing_individual_from_denominator(tmp_path):
    ind1_fasta = _write_fasta(tmp_path / "ind1.fa", "ACGTA")
    ind2_fasta = _write_fasta(tmp_path / "ind2.fa", "ACTA")

    entries = {
        (GENE, "IND1", "H1"): _full_match_entry(ind1_fasta, 5),
        (GENE, "IND2", "H1"): {
            "fasta_path": str(ind2_fasta),
            "copy_from_indices": [0, 1, 2, 3],
            "expanded_indices": [0, 1, 3, 4],
        },
        # IND3/H1 deliberately absent -> get_haplotype_entry returns None.
    }
    mapper = FakeMapper(entries)

    freq = class_base_frequency(mapper, GENE, ["IND1", "IND2", "IND3"], ["H1"], 1, 5)

    assert freq["n_requested"] == 3
    assert freq["n_used"] == 2

    # position 1: IND1='A', IND2='A' -> A frequency 1.0
    assert freq["A"][0] == 1.0
    # position 2: IND1='C', IND2='C' -> C frequency 1.0
    assert freq["C"][1] == 1.0
    # position 3: IND1='G', IND2=gap -> G=0.5, gap=0.5
    assert freq["G"][2] == 0.5
    assert freq["gap"][2] == 0.5
    # position 4: IND1='T', IND2='T' -> T frequency 1.0
    assert freq["T"][3] == 1.0
    # position 5: IND1='A', IND2='A' -> A frequency 1.0
    assert freq["A"][4] == 1.0

    for idx in range(5):
        total = sum(freq[symbol][idx] for symbol in ("A", "C", "G", "T", "gap"))
        assert abs(total - 1.0) < 1e-9


def test_individual_aligned_bases_with_source_maps_gaps_to_none(tmp_path):
    fasta_path = _write_fasta(tmp_path / "ind2.fa", "ACTA")
    entry = {
        "fasta_path": str(fasta_path),
        "copy_from_indices": [0, 1, 2, 3],
        "expanded_indices": [0, 1, 3, 4],
    }
    mapper = FakeMapper({(GENE, "IND2", "H1"): entry})

    bases, source_indices = individual_aligned_bases_with_source(mapper, GENE, "IND2", "H1", 1, 5)

    assert bases == "ACXTA"
    assert source_indices == [0, 1, None, 2, 3]


def test_individual_aligned_bases_with_source_missing_individual(tmp_path):
    mapper = FakeMapper({})
    bases, source_indices = individual_aligned_bases_with_source(mapper, GENE, "GHOST", "H1", 1, 4)
    assert bases == "XXXX"
    assert source_indices == [None, None, None, None]


def test_class_base_frequency_all_missing_returns_zero_used():
    mapper = FakeMapper({})
    freq = class_base_frequency(mapper, GENE, ["GHOST1", "GHOST2"], ["H1", "H2"], 1, 3)
    assert freq["n_used"] == 0
    assert freq["n_requested"] == 4
    for symbol in ("A", "C", "G", "T", "gap"):
        assert list(freq[symbol]) == [0.0, 0.0, 0.0]
