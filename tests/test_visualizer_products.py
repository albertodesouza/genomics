"""Tests for the gene-products view (genomics.visualizer.products).

Small hand-built transcripts check splicing on each haplotype through its indel map, translation
from the annotated start, phase-aware amino-acid changes, NMD, and the AlphaGenome junction
evidence (PSI, candidate isoforms); the HTTP route is exercised on the synthetic dataset of
``test_visualizer`` with a GTF cache table written next to it.
"""
import json
import time

import numpy as np
import pytest

from genomics.visualizer.coords import HaplotypeEvents, reference_map
from genomics.visualizer.products import (
    REFERENCE,
    HaplotypeView,
    Junctions,
    build_product,
    describe_change,
    fasta,
    junction_kind,
    nmd_status,
    retain_intron,
    reverse_complement,
    splice_in,
    splicing,
    translate,
    variant_region,
)
from test_visualizer import L, WINDOW_START, dataset_dir  # noqa: F401 (fixture)

# 5' flank | exon 1 (5'UTR CC + ATG GCT AAA) | intron | exon 2 (GGC TGG TAA + 3'UTR) | 3' flank
UP, EXON1, INTRON, EXON2, DOWN = "TTTTT", "CCATGGCTAAA", "GTAAGTTTTTTTTTTCAG", "GGCTGGTAATTT", "TTTTT"
SEQ = UP + EXON1 + INTRON + EXON2 + DOWN
E1 = (len(UP), len(UP) + len(EXON1))  # [5, 16)
E2 = (E1[1] + len(INTRON), E1[1] + len(INTRON) + len(EXON2))  # [34, 46)
EXONS = [list(E1), list(E2)]
START = E1[0] + 2  # the A of ATG
STOP = E2[0] + 6  # first base of TAA
PROTEIN = "MAKGW"


def view(seq=SEQ, name=REFERENCE, events=(), length=len(SEQ)):
    """A haplotype view of a ``length``-bp reference window; ``events`` are (offset, ref_len,
    alt_len) indels already applied to ``seq``."""
    if name == REFERENCE:
        return HaplotypeView.reference(seq.encode())
    ev = HaplotypeEvents(np.array([e[0] for e in events], np.int64), np.array([e[1] for e in events], np.int32), np.array([e[2] for e in events], np.int32))
    return HaplotypeView.from_map(name, seq.encode(), reference_map(ev, length))


def edit(pos, ref, alt, seq=SEQ):
    assert seq[pos:pos + len(ref)] == ref
    return seq[:pos] + alt + seq[pos + len(ref):]


def product(seq=SEQ, events=(), name="H1", exons=EXONS, strand="+", start=START, stop=STOP):
    return build_product(view(seq, name, events), exons, strand, start, stop, with_sequence=True)


def change(seq, events=()):
    return describe_change(product(SEQ, name=REFERENCE), product(seq, events))


def test_translate_and_reverse_complement():
    assert translate("ATGGCTAAATAAGGG") == ("MAK", True)
    assert translate("ATGGC") == ("M", False)
    assert reverse_complement("ATGCN") == "NGCAT"


def test_reference_transcript_is_spliced_and_translated():
    p = product(name=REFERENCE)
    assert p["mrna"] == EXON1 + EXON2
    assert p["protein"] == PROTEIN and p["stop_found"] and not p["start_lost"]
    assert (p["utr5"], p["utr3"]) == (2, 3)
    assert p["frame_at_ref_stop"] == 0 and p["nmd"]["predicted"] is False


def test_minus_strand_gives_the_same_product():
    n = len(SEQ)
    rc = reverse_complement(SEQ)
    exons = [[n - e, n - s] for s, e in EXONS][::-1]
    p = build_product(view(rc), exons, "-", n - 1 - START, n - 1 - STOP, with_sequence=True)
    assert p["mrna"] == EXON1 + EXON2 and p["protein"] == PROTEIN


def test_two_snvs_in_one_codon_on_one_haplotype_give_one_amino_acid_change():
    # GCT (Ala) -> GAA (Glu) needs both SNVs on the same haplotype; separately they give Asp / Ala.
    seq = edit(START + 4, "CT", "AA")
    c = change(seq)
    assert c == {"class": "missense", "hgvs": "p.Ala2Glu", "count": 1, "position": 2}


@pytest.mark.parametrize("pos,ref,alt,cls,hgvs", [
    (START + 5, "T", "C", "synonymous", "p.(=)"),  # GCT -> GCC
    (E1[0], "C", "G", "utr_change", "p.(=)"),
    (START + 6, "A", "T", "stop_gained", "p.Lys3Ter"),
    (START + 2, "G", "A", "start_lost", "p.Met1?"),
    (STOP, "T", "C", "stop_lost", "p.Ter6GlnextTer?"),  # TAA -> CAA, no stop before the transcript ends
])
def test_snv_consequences(pos, ref, alt, cls, hgvs):
    c = change(edit(pos, ref, alt))
    assert c["class"] == cls
    if hgvs:
        assert c["hgvs"] == hgvs


def test_frameshift_and_inframe_insertion_through_the_indel_map():
    # delete one A of AAA (VCF style: anchor at START+5, ref "TA" -> "T")
    deleted = edit(START + 6, "A", "")
    c = change(deleted, events=[(START + 5, 2, 1)])
    # ATG GCT AAG GCT GGT AAT TT: Lys3 survives (AAA -> AAG), the frame shifts from Gly4 on
    assert c["class"] == "frameshift" and c["hgvs"] == "p.Gly4AlafsTer?"
    # insert GCT after the T of GCT: one extra Ala, frame kept
    inserted = edit(START + 5, "T", "TGCT")
    p = product(inserted, events=[(START + 5, 1, 4)])
    assert p["protein"] == "MAAKGW" and p["frame_at_ref_stop"] == 0
    assert describe_change(product(name=REFERENCE), p) == {"class": "inframe_indel", "hgvs": "p.Ala2_Lys3insAla", "position": 3, "length_change": 1}


def test_nmd_rule_and_escapes():
    assert nmd_status([200, 200, 200], 10, 100, True)["predicted"] is False  # start-proximal
    assert nmd_status([200, 200, 200], 10, 250, True)["predicted"] is True
    assert nmd_status([200, 200, 200], 10, 360, True)["predicted"] is False  # within 50 nt of the last junction
    assert nmd_status([200, 900, 200], 10, 400, True)["predicted"] is False  # long exon
    assert nmd_status([600], 10, 300, True)["reason"] == "single exon"


def test_splice_in_and_retain_intron():
    exons = [[0, 10], [20, 30], [40, 50]]
    assert splice_in(exons, (10, 40)) == [[0, 10], [40, 50]]  # exon skipping
    assert splice_in(exons, (10, 25)) == [[0, 10], [25, 30], [40, 50]]  # acceptor inside exon 2
    assert splice_in(exons, (14, 20)) == [[0, 14], [20, 30], [40, 50]]  # donor inside the intron
    assert retain_intron(exons, (30, 40)) == [[0, 10], [20, 50]]
    introns = [(10, 20), (30, 40)]
    assert junction_kind((10, 40), introns, "+") == "exon skipping (1 exon)"
    assert junction_kind((10, 25), introns, "+") == "alternative acceptor"
    assert junction_kind((10, 25), introns, "-") == "alternative donor"


def test_variant_region():
    tx = {"strand": "+", "exons": EXONS, "cds": [[START, E1[1]], [E2[0], STOP]]}
    assert variant_region(tx, START, START + 1) == "coding"
    assert variant_region(tx, E1[0], E1[0] + 1) == "5' UTR"
    assert variant_region(tx, E2[1] - 1, E2[1]) == "3' UTR"
    assert variant_region(tx, E1[1], E1[1] + 1) == "splice donor"
    assert variant_region(tx, E2[0] - 1, E2[0]) == "splice acceptor"
    assert variant_region(tx, E1[1] + 5, E1[1] + 6) == "splice region (intron)"
    assert variant_region(tx, E1[1] + 9, E1[1] + 10) == "intron"
    assert variant_region(tx, 0, 1) == "upstream"
    assert variant_region({**tx, "strand": "-"}, E1[1], E1[1] + 1) == "splice acceptor"
    assert variant_region(tx, 5000, 5001) is None


def _junction_file(path, junctions, tracks=1):
    starts = np.array([j[0] for j in junctions], np.int64)
    ends = np.array([j[1] for j in junctions], np.int64)
    values = np.array([j[2] for j in junctions], np.float32).reshape(len(junctions), tracks)
    np.savez_compressed(path, starts=starts, ends=ends, strands=np.array(["+"] * len(junctions)), values=values)
    path.with_name(f"{path.stem}_metadata.json").write_text(json.dumps({"metadata": [{"biosample_name": f"cell {i}"} for i in range(tracks)]}), encoding="utf-8")
    return path


def test_junctions_are_mapped_back_through_insertions(tmp_path):
    seq = edit(3, "T", "TAAA")[:len(SEQ)]  # 3 bases inserted upstream, trimmed to the window length
    hap = view(seq, "H1", events=[(3, 1, 4)])
    path = _junction_file(tmp_path / "splice_junctions.npz", [(E1[1] + 3, E2[0] + 3, 2.0), (E1[1] + 3, E2[0] + 13, 1.0)])
    j = Junctions.load(path, hap, "+")
    assert set(j.values) == {(E1[1], E2[0]), (E1[1], E2[0] + 10)}
    psi = j.psi("+", np.array([0.0]))
    assert psi[(E1[1], E2[0])][0][0] == pytest.approx(2 / 3)  # share of the donor
    assert psi[(E1[1], E2[0])][1][0] == pytest.approx(1.0)  # only junction at its acceptor
    # the floor damps weak sites: 0.1 of signal against a floor of 1 is a 0.1 share, not 1.0
    weak = Junctions([], {(1, 5): np.array([0.1])}).psi("+", np.array([1.0]))
    assert weak[(1, 5)][0][0] == pytest.approx(0.1)


def test_splicing_finds_exon_skipping_and_intron_retention(tmp_path):
    seq = "".join("ACGT"[i % 4] for i in range(60))
    exons = [[0, 10], [20, 30], [40, 50]]
    tx = {"id": "T1", "name": "G-201", "strand": "+", "exons": exons, "cds": [], "mane": True}
    views = {REFERENCE: view(seq), "H1": view(seq, "H1", length=len(seq)), "H2": view(seq, "H2", length=len(seq))}
    paths = {
        REFERENCE: _junction_file(tmp_path / "ref.npz", [(10, 20, 2.0), (30, 40, 2.0)]),
        "H1": _junction_file(tmp_path / "h1.npz", [(10, 20, 0.5), (30, 40, 0.5), (10, 40, 1.5)]),
        "H2": tmp_path / "missing.npz",
    }
    ref_products = {"T1": build_product(views[REFERENCE], exons, "+", None, with_sequence=True)}
    sp = splicing({"name": "G"}, tx, [tx], views, paths, "+", ref_products, with_sequence=False)
    assert sp["available"] and sp["missing"] == ["H2"] and sp["expressed"] == [True]
    kinds = {c["kind"]: c for c in sp["candidates"]}
    skip = kinds["exon skipping (1 exon)"]
    assert skip["exons"] == [[0, 10], [40, 50]]
    assert skip["junction"]["H1"]["usage"][0] == pytest.approx(0.75) and skip["junction"]["ref"]["usage"][0] == 0
    assert skip["products"]["H1"]["mrna_length"] == 20
    # both introns lose the same share (to the skipping junction): no intron stands out as retained
    assert "intron retention (junction lost)" not in kinds
    assert [r["H1"]["delta"][0] for r in sp["base_introns"]] == [pytest.approx(-0.75), pytest.approx(-0.75)]


def test_splicing_flags_an_intron_whose_junction_alone_collapses(tmp_path):
    seq = "".join("ACGT"[i % 4] for i in range(80))
    exons = [[0, 10], [20, 30], [40, 50], [60, 70]]
    tx = {"id": "T1", "name": "G-201", "strand": "+", "exons": exons, "cds": [], "mane": True}
    views = {REFERENCE: view(seq), "H1": view(seq, "H1", length=len(seq)), "H2": view(seq, "H2", length=len(seq))}
    paths = {
        REFERENCE: _junction_file(tmp_path / "ref.npz", [(10, 20, 2.0), (30, 40, 2.0), (50, 60, 2.0)]),
        # half the gene's expression on H1, and the middle junction down to a tenth
        "H1": _junction_file(tmp_path / "h1.npz", [(10, 20, 1.0), (30, 40, 0.1), (50, 60, 1.0)]),
        "H2": _junction_file(tmp_path / "h2.npz", [(10, 20, 1.0), (30, 40, 1.0), (50, 60, 1.0)]),
    }
    ref_products = {"T1": build_product(views[REFERENCE], exons, "+", None, with_sequence=True)}
    sp = splicing({"name": "G"}, tx, [tx], views, paths, "+", ref_products, with_sequence=False)
    (retained,) = [c for c in sp["candidates"] if c["kind"].startswith("intron retention")]
    assert (retained["junction"]["start"], retained["junction"]["end"]) == (30, 40)
    assert retained["exons"] == [[0, 10], [20, 50], [60, 70]] and retained["max_delta"] == pytest.approx(0.9)


# ----------------------------------------------------------------------------- HTTP API
def _write_gtf_cache(dataset_dir, gene_id=None):  # noqa: F811
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")
    g = lambda off: off + WINDOW_START - 1  # window offset -> 0-based genomic Start (pyranges)
    rows = [
        ("gene", 5, 140, None, None, None),
        ("transcript", 5, 140, "ENST1.1", "basic,MANE_Select,Ensembl_canonical,cds_start_NF,mRNA_start_NF", "."),
        ("exon", 5, 40, "ENST1.1", None, "."), ("exon", 60, 140, "ENST1.1", None, "."),
        ("CDS", 12, 40, "ENST1.1", None, "1"), ("CDS", 60, 101, "ENST1.1", None, "0"),
    ]
    df = pd.DataFrame({
        "Chromosome": "chr1", "Feature": [r[0] for r in rows], "Start": [g(r[1]) for r in rows], "End": [g(r[2]) for r in rows], "Strand": "+",
        "gene_name": "GENE1", "gene_type": "protein_coding", "transcript_id": [r[3] for r in rows], "transcript_name": [("GENE1-201" if r[3] else None) for r in rows],
        "transcript_type": [("protein_coding" if r[3] else None) for r in rows], "tag": [r[4] for r in rows], "exon_number": None, "Frame": [r[5] or "." for r in rows],
    })
    if gene_id:
        df["gene_id"] = gene_id
    df.to_feather(dataset_dir / "gtf_cache.feather")


def test_http_products_and_fasta(dataset_dir, tmp_path):  # noqa: F811
    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.server import HttpError, Query, RawResponse, VisualizerApp

    _write_gtf_cache(dataset_dir)
    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=128 << 20, workers=2, runs_roots=[], remote=False)

    def get(path, **params):
        from urllib.parse import urlencode

        for _ in range(400):
            out = app.dispatch("GET", path, Query(urlencode(params)), None)
            if not (isinstance(out, dict) and out.get("pending")):
                return out
            time.sleep(0.02)
        raise AssertionError("job did not finish")

    try:
        base = f"/api/d/{ds.id}"
        data = get(f"{base}/products", gene="GENE1", sample="S1")
        assert data["target"] == "GENE1" and data["base"] == "ENST1.1" and data["strand"] == "+"
        (tx,) = data["transcripts"]
        assert tx["incomplete"] == {"start": True, "end": False} and tx["mane"]
        ref = tx["products"]["ref"]
        assert ref["coding"] and ref["mrna_length"] == 35 + 80 and "mrna" not in ref
        assert ref["cds_start"] == 7 + 1  # 7 nt of 5' UTR, then the phase-1 base skipped
        for hap in ("H1", "H2"):
            assert tx["products"][hap]["change"]["class"] in {"missense", "frameshift", "inframe_indel", "stop_gained", "synonymous", "no_change", "utr_change", "stop_lost"}
        assert {v["region"] for v in data["variants"]} <= {"coding", "5' UTR", "3' UTR", "intron", "upstream", "downstream", "splice region (intron)",
                                                             "splice donor", "splice acceptor", "coding, splice region", "5' UTR, splice region", "3' UTR, splice region"}
        assert data["splicing"]["available"] is False  # the fixture has no splice_junctions predictions
        with_seq = get(f"{base}/products", gene="GENE1", sample="S1", sequences="mrna")
        assert len(with_seq["transcripts"][0]["products"]["H1"]["mrna"]) == with_seq["transcripts"][0]["products"]["H1"]["mrna_length"]
        raw = get(f"{base}/products/fasta", gene="GENE1", sample="S1", kind="protein")
        assert isinstance(raw, RawResponse) and raw.data.decode().startswith(">GENE1|GENE1-201|ref")
        assert fasta(with_seq, "mrna").count(">") == 3
        with pytest.raises(HttpError):
            get(f"{base}/products", gene="GENE1", sample="../S1")
        missing = get(f"{base}/products", gene="GENE1", sample="NOPE")
        assert missing["transcripts"] == [] and "no GENE1 window" in missing["message"]
    finally:
        app.shutdown()
        (dataset_dir / "gtf_cache.feather").unlink()
