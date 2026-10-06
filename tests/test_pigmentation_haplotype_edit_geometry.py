import gzip

import pytest

from genomics.predictors.genotype_based.analysis.haplotype_edit_geometry import (
    apply_overwrite,
    apply_scramble,
    apply_scramble_region,
    genomic_pos_from_haplotype_local_idx,
    haplotype_indel_events,
    haplotype_local_idx,
)

VCF_HEADER = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1",
    ]
)


def _write_vcf(path, rows):
    lines = [VCF_HEADER]
    for pos, ref, alt, gt in rows:
        lines.append(f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}")
    with gzip.open(path, "wt") as f:
        f.write("\n".join(lines) + "\n")


def test_haplotype_indel_events_only_this_haplotype(tmp_path):
    vcf_path = tmp_path / "gene.window.vcf.gz"
    _write_vcf(
        vcf_path,
        [
            (100, "A", "AGG", "1|0"),  # H1 insertion (+2), H2 ref
            (200, "ACGT", "A", "0|1"),  # H1 ref, H2 deletion (-3)
            (300, "A", "T", "1|1"),  # SNV both haplotypes, no length change
        ],
    )

    h1_events = haplotype_indel_events(vcf_path, "H1")
    assert h1_events == [(100, 2)]

    h2_events = haplotype_indel_events(vcf_path, "H2")
    assert h2_events == [(200, -3)]


def test_haplotype_indel_events_rejects_bad_haplotype(tmp_path):
    vcf_path = tmp_path / "gene.window.vcf.gz"
    _write_vcf(vcf_path, [])
    with pytest.raises(ValueError):
        haplotype_indel_events(vcf_path, "H3")


def test_haplotype_local_idx_no_events_is_identity():
    events = []
    start_1based = 1000
    assert haplotype_local_idx(events, start_1based, 1000) == 0
    assert haplotype_local_idx(events, start_1based, 1050) == 50


def test_haplotype_local_idx_after_insertion_shifts_forward():
    events = [(1010, 5)]  # +5bp insertion at ref pos 1010
    start_1based = 1000
    assert haplotype_local_idx(events, start_1based, 1005) == 5  # before the insertion: unaffected
    assert haplotype_local_idx(events, start_1based, 1020) == 20 + 5  # after: shifted by +5


def test_haplotype_local_idx_after_deletion_shifts_backward():
    events = [(1010, -4)]  # 4bp deletion at ref pos 1010
    start_1based = 1000
    assert haplotype_local_idx(events, start_1based, 1005) == 5
    assert haplotype_local_idx(events, start_1based, 1020) == 20 - 4


def test_genomic_pos_inverts_local_idx_away_from_indels():
    start_1based = 1000
    for events in ([], [(1010, 5)], [(1010, -4)]):
        for target_pos in (1000, 1002, 1005, 1030, 1050):
            local_idx = haplotype_local_idx(events, start_1based, target_pos)
            recovered = genomic_pos_from_haplotype_local_idx(events, start_1based, local_idx)
            assert recovered == target_pos, (events, target_pos)


def test_apply_overwrite_replaces_region_with_repeated_base():
    result = apply_overwrite("ACGTACGT", 2, 5, "A")
    assert result.modified_sequence == "ACAAACGT"
    assert result.original_segment == "GTA"
    assert result.edited_segment == "AAA"
    assert result.start == 2
    assert result.end == 5


@pytest.mark.parametrize("bad_base", ["N", "X", "a", "AC", ""])
def test_apply_overwrite_rejects_invalid_base(bad_base):
    with pytest.raises(ValueError):
        apply_overwrite("ACGTACGT", 0, 2, bad_base)


def test_apply_overwrite_rejects_out_of_range_region():
    with pytest.raises(ValueError):
        apply_overwrite("ACGT", 3, 2, "A")
    with pytest.raises(ValueError):
        apply_overwrite("ACGT", 0, 10, "A")


def test_apply_scramble_region_preserves_length_and_only_touches_span():
    sequence = "ACGTACGTACGT"
    result = apply_scramble_region(sequence, 3, 9, seed=0)
    assert len(result.modified_sequence) == len(sequence)
    assert result.modified_sequence[:3] == sequence[:3]
    assert result.modified_sequence[9:] == sequence[9:]
    assert sorted(result.edited_segment) == sorted(sequence[3:9])


def test_apply_scramble_region_same_seed_is_reproducible():
    sequence = "ACGTACGTACGTACGT"
    r1 = apply_scramble_region(sequence, 0, 16, seed=42)
    r2 = apply_scramble_region(sequence, 0, 16, seed=42)
    assert r1.modified_sequence == r2.modified_sequence


def test_apply_scramble_wrapper_matches_center_window_semantics():
    sequence = "A" * 200
    result = apply_scramble("C" * 50 + sequence[50:], 100, window_size=20, seed=1)
    assert result.end - result.start == 20
    assert result.start == 90 and result.end == 110
