import gzip
import json

import numpy as np

from genomics.predictors.genotype_based.analysis.indel_track_realignment import (
    load_haplotype_variants,
    realign_consensus_values,
)
from genomics.predictors.genotype_based.alignment.reference_realign_mapper import ReferenceRealignMapper


def _vcf_lines(records):
    header = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    body = [
        f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\tGT\t{gt}"
        for pos, ref, alt, info, gt in records
    ]
    return "\n".join(header + body) + "\n"


def _write_vcf_gz(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write(_vcf_lines(records))


# ---------------------------------------------------------------------------
# realign_consensus_values / load_haplotype_variants (zero coverage before this file)
# ---------------------------------------------------------------------------


def test_realign_no_variants_is_identity_copy():
    consensus = np.arange(10, dtype=np.float32).reshape(-1, 1)
    out = realign_consensus_values(consensus, window_start=100, window_width=10, variants=[])
    assert out.shape == (10, 1)
    assert np.array_equal(out[:, 0], consensus[:, 0])


def test_realign_snp_does_not_shift_downstream_positions():
    consensus = np.arange(10, dtype=np.float32).reshape(-1, 1)
    # SNP at window_start + 2 (ref_local index 2): ref_len == alt_len == 1, no length change.
    variants = [(102, "A", "T")]
    out = realign_consensus_values(consensus, window_start=100, window_width=10, variants=variants)
    assert np.array_equal(out[:, 0], consensus[:, 0])


def test_realign_insertion_collapses_via_max_onto_anchor():
    # Reference window is 6 bases wide. The individual's consensus carries a 1bp->3bp insertion
    # at the 2nd reference base, so consensus is 8 bases long for this same window.
    consensus = np.array([0.0, 1.0, 5.0, 9.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32).reshape(-1, 1)
    variants = [(101, "A", "TGG")]  # anchor at window_start+1 (ref_local index 1)
    out = realign_consensus_values(consensus, window_start=100, window_width=6, variants=variants)
    assert out.shape == (6, 1)
    # position 0: untouched copy
    assert out[0, 0] == 0.0
    # position 1 (the anchor): max over the inserted span [1,2,3] in consensus -> max(1,5,9) = 9
    assert out[1, 0] == 9.0
    # remaining positions continue reading from consensus_idx=4 onward, 1:1
    assert np.array_equal(out[2:, 0], consensus[4:8, 0])


def test_realign_deletion_zero_fills_missing_reference_positions():
    # Reference window is 6 bases wide; a 3bp->1bp deletion removes 2 reference-only bases.
    consensus = np.array([0.0, 1.0, 2.0, 3.0], dtype=np.float32).reshape(-1, 1)
    variants = [(101, "ACG", "A")]  # anchor at ref_local index 1, deletes ref_local 2,3
    out = realign_consensus_values(consensus, window_start=100, window_width=6, variants=variants)
    assert out.shape == (6, 1)
    assert out[0, 0] == 0.0
    assert out[1, 0] == 1.0  # anchor keeps real consensus signal
    assert out[2, 0] == 0.0  # deleted reference-only position
    assert out[3, 0] == 0.0  # deleted reference-only position
    assert np.array_equal(out[4:, 0], consensus[2:4, 0])


def test_load_haplotype_variants_filters_by_haplotype_and_genotype(tmp_path):
    vcf_path = tmp_path / "sample.window.consensus_ready.vcf.gz"
    _write_vcf_gz(vcf_path, [
        (100, "A", "T", "AC=1", "0|0"),   # homref, never carried by either haplotype
        (200, "C", "G", "AC=1", "1|0"),   # H1 carries ALT, H2 does not
        (300, "G", "A", "AC=1", "0|1"),   # H2 carries ALT, H1 does not
        (400, "AT", "A", "AC=1", "1|1"),  # both haplotypes carry the deletion
    ])

    h1_variants = load_haplotype_variants(vcf_path, "H1")
    h2_variants = load_haplotype_variants(vcf_path, "H2")

    assert h1_variants == [(200, "C", "G"), (400, "AT", "A")]
    assert h2_variants == [(300, "G", "A"), (400, "AT", "A")]


def test_load_haplotype_variants_resolves_symbolic_alt_from_info_end(tmp_path):
    vcf_path = tmp_path / "sample.window.consensus_ready.vcf.gz"
    _write_vcf_gz(vcf_path, [
        (500, "N", "<DEL>", "END=505", "1|0"),  # symbolic 6bp deletion, H1 only
    ])

    h1_variants = load_haplotype_variants(vcf_path, "H1")
    h2_variants = load_haplotype_variants(vcf_path, "H2")

    assert h1_variants == [(500, "N" * 6, "N")]
    assert h2_variants == []


def _make_dataset_dir(tmp_path, gene, chrom, start, end):
    window_dir = tmp_path / "references" / "windows" / gene
    window_dir.mkdir(parents=True, exist_ok=True)
    (window_dir / "window_metadata.json").write_text(json.dumps({
        "type": "gene",
        "chromosome": chrom,
        "start": start,
        "end": end,
        "window_size": end - start,  # deliberately the off-by-one convention seen in real datasets
    }))
    return tmp_path


# ---------------------------------------------------------------------------
# ReferenceRealignMapper end-to-end (synthetic dataset dir, no external files needed)
# ---------------------------------------------------------------------------


def test_reference_realign_mapper_no_variants_center_crops_identity(tmp_path):
    gene = "TESTGENE"
    start, end = 1000, 1019  # 1-based inclusive -> full_ref_length = 20
    dataset_dir = _make_dataset_dir(tmp_path, gene, "chr1", start, end)

    sample_id = "SAMPLE1"
    vcf_path = (
        dataset_dir / "individuals" / sample_id / "windows" / gene
        / f"{sample_id}.window.consensus_ready.vcf.gz"
    )
    _write_vcf_gz(vcf_path, [])  # no variants at all -> both haplotypes are pure reference

    mapper = ReferenceRealignMapper(
        dataset_dir=dataset_dir,
        consensus_dataset_dir=dataset_dir,
        window_center_size=10,
    )
    row = np.arange(20, dtype=np.float32)
    out = mapper.get_haplotype_track(gene, sample_id, "H1", row)

    assert out is not None
    assert out.shape == (10,)
    # centered crop of a 20-long identity row to width 10 -> positions [5:15)
    assert np.array_equal(out, row[5:15])


def test_reference_realign_mapper_missing_vcf_returns_none(tmp_path):
    gene = "TESTGENE"
    dataset_dir = _make_dataset_dir(tmp_path, gene, "chr1", 1000, 1019)
    mapper = ReferenceRealignMapper(
        dataset_dir=dataset_dir,
        consensus_dataset_dir=dataset_dir,
        window_center_size=10,
    )
    out = mapper.get_haplotype_track(gene, "NOSUCHSAMPLE", "H1", np.arange(20, dtype=np.float32))
    assert out is None
    assert mapper.profile_stats["missing_vcf"] == 1


def test_reference_realign_mapper_insertion_matches_manual_realign(tmp_path):
    gene = "TESTGENE"
    start, end = 1000, 1005  # full_ref_length = 6
    dataset_dir = _make_dataset_dir(tmp_path, gene, "chr1", start, end)

    sample_id = "SAMPLE1"
    vcf_path = (
        dataset_dir / "individuals" / sample_id / "windows" / gene
        / f"{sample_id}.window.consensus_ready.vcf.gz"
    )
    # Insertion anchored at position 1001 (ref_local index 1), carried on H1 only.
    _write_vcf_gz(vcf_path, [(1001, "A", "TGG", "AC=1", "1|0")])

    mapper = ReferenceRealignMapper(
        dataset_dir=dataset_dir,
        consensus_dataset_dir=dataset_dir,
        window_center_size=6,
    )
    # Raw AlphaGenome prediction is always exactly `full_ref_length` long (the dataset builder
    # pads/truncates the consensus sequence to the reference window size before ever submitting it
    # to AlphaGenome -- see `individual_consensus._adjust_to_target_size`), even for a haplotype
    # that carries an insertion: the extra inserted bases get folded in via the max-collapse, and
    # the truncation this implies means the trailing reference-only positions have no consensus
    # counterpart left and legitimately realign to zero.
    row = np.array([0.0, 1.0, 5.0, 9.0, 2.0, 3.0], dtype=np.float32)
    out = mapper.get_haplotype_track(gene, sample_id, "H1", row)

    expected = realign_consensus_values(
        row.reshape(-1, 1), window_start=start, window_width=6,
        variants=[(1001, "A", "TGG")],
    )[:, 0]
    assert out is not None
    assert np.array_equal(out, expected)
    assert np.array_equal(out, np.array([0.0, 9.0, 2.0, 3.0, 0.0, 0.0], dtype=np.float32))

    # H2 does not carry the insertion -> straight reference-length pass-through, no manual variant
    # list needed for this assertion since H2's genotype token is "0".
    row_h2 = np.arange(6, dtype=np.float32)
    out_h2 = mapper.get_haplotype_track(gene, sample_id, "H2", row_h2)
    assert out_h2 is not None
    assert np.array_equal(out_h2, row_h2)
