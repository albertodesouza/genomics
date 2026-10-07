"""Tests for the expression (relative vs absolute) and structure views of the gene-products page.

Pure helpers are checked directly; the HTTP routes run on the synthetic dataset of
``test_visualizer`` with constant RNA-seq predictions (so every fold is known), a faked HPA
download and GTEx, and remote lookups otherwise disabled.
"""
import hashlib
import io
import json
import time
import zipfile
from urllib.parse import urlencode

import numpy as np
import pytest

from genomics.visualizer.alphagenome import AlphaGenomeUnavailable, predict_cached
from genomics.visualizer.expression import ReferenceLevels, exon_offsets, fold_label, mean_over, productive_share
from genomics.visualizer.products import HaplotypeView
from genomics.visualizer.structures import compare_proteins, fold_window, mean_b_factor
from test_visualizer import L, _write_reference_prediction, dataset_dir  # noqa: F401 (fixture)
from test_visualizer_products import _write_gtf_cache
from test_visualizer_variants import FakeRemote

# ----------------------------------------------------------------------------- helpers


def test_compare_proteins_modes():
    assert compare_proteins("MAKGW", "MAKGW") == {"mode": "identical"}
    assert compare_proteins("MAKGW", "MAEGW") == {"mode": "substitution", "substitutions": [{"pos": 3, "ref": "K", "alt": "E"}]}
    assert compare_proteins("MAKGW", "MAK") == {"mode": "truncated", "length": 3, "lost_from": 4}
    assert compare_proteins("MAKGW", "MAKRSTV") == {"mode": "altered", "first_change": 4}
    assert compare_proteins("MAKGW", "MAKGWQQ") == {"mode": "altered", "first_change": 6}


def test_fold_window_keeps_context_and_limit():
    assert fold_window(300, 250) == [1, 300]
    assert fold_window(1000, 600) == [450, 849]  # 150 residues of context before the change
    assert fold_window(1000, 990) == [601, 1000]  # clamped to the end, still 400 long
    assert fold_window(1000, 20) == [1, 400]


def test_mean_b_factor_reads_ca_atoms():
    pdb = "\n".join([
        "ATOM      1  N   MET A   1      13.733   9.331  31.688  1.00 10.00           N",
        "ATOM      2  CA  MET A   1      12.801   8.455  30.982  1.00 80.00           C",
        "ATOM      3  CA  ALA A   2      12.801   8.455  30.982  1.00 60.00           C",
    ])
    assert mean_b_factor(pdb) == 70.0


def test_exon_mean_fold_and_labels():
    transcripts = [{"exons": [[2, 5], [8, 10]]}, {"exons": [[4, 6]]}]
    offsets = exon_offsets(transcripts, 12)
    assert offsets.tolist() == [2, 3, 4, 5, 8, 9]
    values = np.arange(24, dtype=np.float32).reshape(12, 2)
    ref = HaplotypeView.reference(b"A" * 12)
    assert mean_over(values, ref.local, offsets, [1]) == pytest.approx(np.mean([5, 7, 9, 11, 17, 19]))
    shifted = ref.local + 1  # a haplotype with one extra base upstream
    assert mean_over(values, shifted, offsets, [0]) == pytest.approx(np.mean([6, 8, 10, 12, 18, 20]))
    assert [fold_label(v) for v in (1.2, 0.85, 1.05, None)] == ["higher", "lower", "similar", "unknown"]


def test_productive_share_removes_nmd_and_missing_proteins():
    product = {"coding": True, "protein": "MAK", "nmd": {"predicted": False}}
    cands = [
        {"junction": {"H1": {"usage": [0.3]}}, "products": {"H1": {"protein": "MA", "nmd": {"predicted": True}}}},
        {"junction": {"H1": {"usage": [0.2]}}, "products": {"H1": {"protein": "MAKK", "nmd": {"predicted": False}}}},
    ]
    assert productive_share(product, cands, "H1", 0) == pytest.approx(0.7)
    assert productive_share(product, cands, "H1", None) == 1.0
    assert productive_share({**product, "nmd": {"predicted": True}}, cands, "H1", 0) == 0.0


def test_predict_cached_serves_the_cache_without_a_backend(tmp_path):
    seq = b"ACGT" * 4
    key = hashlib.sha1(seq + json.dumps([["rna_seq"], ["CL:1"]]).encode()).hexdigest()
    np.savez_compressed(tmp_path / f"{key}.npz", rna_seq=np.ones((16, 2), np.float32))
    (tmp_path / f"{key}.json").write_text(json.dumps({"rna_seq": [{"ontology_curie": "CL:1"}]}), encoding="utf-8")
    values, meta = predict_cached(None, tmp_path, seq, ["rna_seq"], ["CL:1"])["rna_seq"]
    assert values.shape == (16, 2) and meta[0]["ontology_curie"] == "CL:1"
    with pytest.raises(AlphaGenomeUnavailable):
        predict_cached(None, tmp_path, b"TTTT", ["rna_seq"], ["CL:1"])


class BytesRemote:
    def __init__(self, data):
        self.data = data
        self.offline = False

    def get_bytes(self, url, ttl=None, headers=None, params=None):
        return self.data


def _hpa_zip(rows):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("rna_single_cell_type.tsv", "Gene\tGene name\tCell type\tnCPM\n" + "".join(f"{g}\t{n}\t{c}\t{v}\n" for g, n, c, v in rows))
    return buf.getvalue()


def test_reference_levels_from_hpa_and_gtex(tmp_path):
    from genomics.visualizer.expression import GTEX_ANCHOR, HPA_ANCHOR
    from genomics.visualizer.gtex import GtexClient

    remote = BytesRemote(_hpa_zip([("ENSG0001", "GENE1", "melanocytes", "120.5"), ("ENSG0001", "GENE1", "b-cells", "3")]))
    gtex = GtexClient(FakeRemote({"reference/gene": {"data": [{"geneSymbol": "GENE1", "gencodeId": "ENSG0001.5"}]},
                                  "expression/medianGeneExpression": {"data": [{"median": 4.25, "unit": "TPM"}]}}))
    levels = ReferenceLevels(remote, gtex, tmp_path)
    hpa = levels.level(HPA_ANCHOR, "GENE1", "ENSG0001.7")
    assert hpa["value"] == 120.5 and hpa["unit"] == "nCPM"
    assert (tmp_path / "expression" / "hpa_single_cell_melanocytes.json").exists()  # parsed once, cached
    assert levels.level(HPA_ANCHOR, "GENE1", None)["value"] is None
    assert levels.level(GTEX_ANCHOR, "GENE1", "ENSG0001.7")["value"] == 4.25
    missing = levels.level(GTEX_ANCHOR, "NOPE", None)
    assert missing["value"] is None and "error" in missing


# ----------------------------------------------------------------------------- HTTP API
def _constant_rna_seq(dataset_dir, values):  # noqa: F811
    """RNA-seq of melanocyte of skin (+ strand) = a constant per haplotype: ref, H1, H2 of S1."""
    meta = {"metadata": [{"ontology_curie": "CL:1000458", "biosample_name": "melanocyte of skin", "strand": "+"}]}
    folders = {"ref": dataset_dir / "references" / "windows" / "GENE1" / "predictions_ref"}
    for hap in ("H1", "H2"):
        folders[hap] = dataset_dir / "individuals" / "S1" / "windows" / "GENE1" / f"predictions_{hap}"
    for hap, folder in folders.items():
        folder.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(folder / "rna_seq.npz", values=np.full((L, 1), values[hap], np.float32))
        (folder / "rna_seq_metadata.json").write_text(json.dumps(meta), encoding="utf-8")


@pytest.fixture
def api(dataset_dir, tmp_path):  # noqa: F811
    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.gtex import GtexClient
    from genomics.visualizer.server import Query, VisualizerApp

    _write_reference_prediction(dataset_dir)
    _constant_rna_seq(dataset_dir, {"ref": 2.0, "H1": 3.0, "H2": 1.0})
    _write_gtf_cache(dataset_dir, gene_id="ENSG0001.3")
    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=128 << 20, workers=2, runs_roots=[], remote=False)
    app.expression.references.remote = BytesRemote(_hpa_zip([("ENSG0001", "GENE1", "melanocytes", "100")]))
    app.expression.references.gtex = GtexClient(FakeRemote({}))

    def get(path, method="GET", _gtex=None, **params):
        if _gtex is not None:
            app.gtex = _gtex
        for _ in range(500):
            out = app.dispatch(method, f"/api/d/{ds.id}{path}", Query(urlencode(params)), None)
            if not (isinstance(out, dict) and out.get("pending")):
                return out
            time.sleep(0.02)
        raise AssertionError("job did not finish")

    def tissues(catalog, save):
        save(catalog, tmp_path / "cache" / "alphagenome_catalog.json")
        app._ag_catalog = None
        return app.dispatch("GET", "/api/products/tissues", Query(""), None)

    get.tissues = tissues
    try:
        yield get
    finally:
        app.shutdown()
        (dataset_dir / "gtf_cache.feather").unlink()


def test_http_expression_relative_and_absolute(api):
    data = api("/products/expression", gene="GENE1", sample="S1")
    assert data["available"] and data["target"] == "GENE1"
    tissues = {t["ontology"]: t for t in data["tissues"]}
    mel = tissues["CL:1000458"]
    assert mel["available"] and mel["source"] == "stored"
    fold = mel["mrna"]["fold"]
    assert (fold["H1"], fold["H2"], fold["individual"]) == (pytest.approx(1.5), pytest.approx(0.5), pytest.approx(1.0))
    assert mel["mrna"]["verdict"] == "similar"
    absolute = mel["mrna"]["absolute"]
    assert absolute["ref"] == 100 and absolute["H1"] == pytest.approx(75) and absolute["H2"] == pytest.approx(25) and absolute["individual"] == pytest.approx(100)
    assert mel["expressed"] is True and mel["anchor"]["unit"] == "nCPM"
    assert mel["protein"]["absolute"] is None and "proteomics" in data["notes"]["protein_absolute"]
    share = mel["protein"]["productive_share"]
    if share["ref"] > 0:
        assert mel["protein"]["fold"]["H1"] == pytest.approx(1.5 * share["H1"] / share["ref"])
    # the other tissues are not stored and there is no AlphaGenome backend: reported, not failed
    assert not tissues["EFO:0000572"]["available"] and data["problems"]
    assert tissues["EFO:0000572"]["anchor"]["value"] is None


def test_http_structure_offline_and_rna(api):
    payload = api("/products", gene="GENE1", sample="S1")
    tx = payload["transcripts"][0]["id"]
    data = api("/products/structure", gene="GENE1", sample="S1", tx=tx)
    assert data["available"] and "error" in data["model"]  # remote lookups are disabled in tests
    assert set(data["products"]) == {"ref", "H1", "H2"} and data["products"]["ref"]["mode"] == "identical"
    assert data["esmfold"]["remote"] is False
    from genomics.visualizer.server import HttpError

    with pytest.raises(HttpError, match="no-remote"):
        api("/products/fold", method="POST", gene="GENE1", sample="S1", tx=tx, hap="H1")
    pytest.importorskip("RNA")
    rna = api("/products/rna", gene="GENE1", sample="S1", tx=tx, hap="H1", center="start", width="60")
    assert rna["available"] and rna["hap"] == "H1" and rna["width"] == 60
    ref, hap = rna["reference"], rna["haplotype"]
    assert len(ref["sequence"]) == len(ref["structure"]) == len(ref["xy"]) == 60
    assert set(ref["sequence"]) <= set("ACGUN") and isinstance(rna["ddg"], float)
    assert all(0 <= i < len(hap["sequence"]) for i in hap["changed"])


# ----------------------------------------------------------------------------- similarity
def test_global_alignment_and_protein_similarity():
    from genomics.visualizer.structures import global_alignment, protein_similarity

    ref = "MKTAYIAKQRQISFVKSHFSRQ"
    assert protein_similarity(ref, ref)["similarity"] == 1.0
    one = protein_similarity(ref, ref[:5] + "W" + ref[6:])  # I6W
    assert one["substitutions"] == 1 and one["identity"] == pytest.approx(21 / 22, abs=1e-4) and 0.8 < one["similarity"] < 1
    trunc = protein_similarity(ref, ref[:11])
    assert trunc["coverage"] == pytest.approx(0.5) and trunc["similarity"] < 0.6
    # a deletion inside the protein: the alignment skips it and keeps the rest aligned
    pairs = global_alignment(ref, ref[:8] + ref[11:])
    assert len(pairs) == len(ref) - 3 and pairs[8] == (11, 8)
    assert protein_similarity(ref, "")["coverage"] == 0


def test_tm_score_is_rotation_invariant_and_low_for_unrelated_coordinates():
    from genomics.visualizer.structures import tm_score

    rng = np.random.default_rng(0)
    x = np.cumsum(rng.normal(size=(80, 3)) * 2.5, axis=0)  # a random-walk "chain"
    th = 1.1
    rot = np.array([[np.cos(th), 0, np.sin(th)], [0, 1, 0], [-np.sin(th), 0, np.cos(th)]])
    same = tm_score(x, x @ rot + np.array([3.0, -7.0, 1.0]), 80)
    assert same["tm"] == pytest.approx(1.0, abs=1e-6) and same["rmsd"] < 1e-6
    other = tm_score(x, np.cumsum(rng.normal(size=(80, 3)) * 2.5, axis=0), 80)
    assert other["tm"] < 0.4
    half = tm_score(x[:40], x[:40] @ rot, 80)  # half the reference aligned
    assert half["tm"] == pytest.approx(0.5, abs=1e-6)


def test_residue_context_reads_plddt_and_burial():
    from genomics.visualizer.structures import residue_context

    lines = []
    for i in range(1, 31):  # residues 1-20 packed in a 3 Å cube, 21-30 far away
        x, y, z = (i % 3) * 3.0, ((i // 3) % 3) * 3.0, (i // 9) * 3.0
        if i > 20:
            x, y, z = 100.0 + 10 * i, 0.0, 0.0
        lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{50 + i:6.2f}           C")
    ctx = residue_context("\n".join(lines), [5, 25, 99])
    assert ctx[5]["burial"] == "buried" and ctx[5]["plddt"] == 55.0
    assert ctx[25]["burial"] == "surface" and 99 not in ctx


# ----------------------------------------------------------------------------- per-variant effects
def test_variant_effects_classify_each_variant_alone():
    from types import SimpleNamespace

    from genomics.visualizer.products import build_product, variant_effects
    from test_visualizer_products import E1, E2, EXONS, SEQ, START, STOP, view

    tx = {"id": "T1", "name": "G-201", "strand": "+", "exons": EXONS, "cds": [[START, E1[1]], [E2[0], STOP]], "mane": True}
    ref_view = view()
    ref_products = {"T1": build_product(ref_view, EXONS, "+", START, STOP, with_sequence=True)}
    records = [  # (offset, ref, alt): synonymous GCT>GCC, missense GCT>GAT, nonsense AAA>TAA, frameshift, UTR, intron
        (START + 5, "T", "C"), (START + 4, "C", "A"), (START + 6, "A", "T"), (START + 6, "AA", "A"), (E1[0], "C", "G"), (E1[1] + 9, "T", "A")]
    variants = SimpleNamespace(
        positions=np.array([o + 1000 for o, _r, _a in records]), ids=[f"v{i}" for i in range(len(records))], refs=[r for _o, r, _a in records],
        alt_alleles=[[a] for _o, _r, a in records], genotypes=["0|1"] * len(records), carried=np.array([[0, 1]] * len(records)))
    rows = variant_effects(ref_view, [tx], tx, "+", ref_products, variants, window_start=1000)
    by = {r["id"]: r for r in rows}
    assert [by[f"v{i}"]["category"] for i in range(6)] == ["synonymous", "missense", "nonsense", "frameshift", "utr", "intron"]
    assert by["v1"]["base_change"]["hgvs"] == "p.Ala2Asp" and by["v2"]["base_change"]["hgvs"] == "p.Lys3Ter"
    assert all(r["haplotypes"] == ["H2"] for r in rows)
    assert rows[0]["category"] == "frameshift"  # most severe first


def test_one_letter_and_gtex_intron_chain_matching():
    from genomics.visualizer.gtex import GtexClient
    from genomics.visualizer.report import ReportService, _one_letter

    assert _one_letter("p.Arg151Cys") == "R151C" and _one_letter("p.[Met708Ile;Glu920Val]") is None
    exons_v26 = [{"transcriptId": "ENST9.2", "start": s, "end": e} for s, e in [(1005, 1040), (1061, 1140)]]
    gtex = GtexClient(FakeRemote({
        "reference/gene": {"data": [{"geneSymbol": "GENE1", "gencodeId": "ENSG0001.5"}]},
        "expression/medianTranscriptExpression": {"data": [{"transcriptId": "ENST9.2", "median": 8.0}, {"transcriptId": "ENST2.1", "median": 2.0}]},
        "reference/exon": {"data": exons_v26},
    }))
    service = ReportService(type("App", (), {"gtex": gtex})())
    # window offset 0 = position 1001: exons [4, 40) and [60, 140) are 1005-1040 and 1061-1140
    current = [{"id": "ENST1.3", "name": "G-201", "exons": [[4, 40], [60, 140]], "mane": True},  # new id, same intron chain as ENST9
               {"id": "ENST2.1", "name": "G-202", "exons": [[10, 50]]}]  # same id
    out = service._transcript_tpms("GENE1", {"tissue": "Cells_EBV-transformed_lymphocytes", "label": "GTEx LCL", "proxy": False}, current, 1001)
    assert out["values"] == {"ENST1.3": 8.0, "ENST2.1": 2.0} and out["matched_by"] == {"ENST1.3": "ENST9"}


def test_http_report_baseline_and_haplotypes(api):
    from genomics.visualizer.gtex import GtexClient

    gtex = GtexClient(FakeRemote({
        "reference/gene": {"data": [{"geneSymbol": "GENE1", "gencodeId": "ENSG0001.5"}]},
        "expression/medianTranscriptExpression": {"data": [{"transcriptId": "ENST1.1", "median": 5.0}]},
    }))
    report = api("/products/report", gene="GENE1", sample="S1", tissue="CL:1000458", _gtex=gtex)
    assert report["available"] and report["tissue"]["ontology"] == "CL:1000458"
    (base,) = report["baseline"]
    assert base["share"] == 1.0 and base["level"] == 100 and report["level"]["unit"] == "nCPM"
    h1, h2 = report["haplotypes"]["H1"], report["haplotypes"]["H2"]
    assert h1["gene_fold"] == pytest.approx(1.5) and h2["gene_fold"] == pytest.approx(0.5)
    t1 = h1["transcripts"][0]
    nmd = 0.2 if t1["nmd_gained"] else 1.0
    assert t1["fold"] == pytest.approx(1.5 * nmd) and t1["level"] == pytest.approx(50 * 1.5 * nmd)
    assert set(h1["flags"]) >= {"unstable", "premature_stop", "missense", "base_similarity"}
    assert isinstance(h1["variants"], list) and sum(h1["counts"].values()) == len(h1["variants"])


# ----------------------------------------------------------------------------- any AlphaGenome tissue
def test_tissue_spec_sources():
    from genomics.visualizer.expression import POOLED_SHARES, cell_line_key, cell_type_key, tissue_spec

    assert cell_type_key("melanocyte of skin") == cell_type_key("Melanocytes") == "melanocyte"
    assert cell_type_key("T-cell") == cell_type_key("t-cells") and cell_line_key("Hep-G2") == cell_line_key("HepG2")
    assert tissue_spec("CL:1000458", None)["anchor"]["cell_type"] == "melanocytes"  # curated default
    gtex = tissue_spec("UBERON:0004264", {"name": "lower leg skin", "type": "tissue", "gtex_tissue": "Skin_Sun_Exposed_Lower_leg"})
    assert gtex["anchor"]["source"] == "gtex" and gtex["transcripts"]["tissue"] == "Skin_Sun_Exposed_Lower_leg"
    by_term = tissue_spec("UBERON:1", {"name": "x", "type": "tissue"}, gtex_by_ontology={"UBERON:1": "Liver"})
    assert by_term["anchor"]["tissue"] == "Liver" and not by_term["anchor"]["matched_by_name"]
    by_name = tissue_spec("UBERON:2", {"name": "thyroid gland", "type": "tissue"}, gtex_by_name={"thyroid": "Thyroid"})
    assert by_name["anchor"]["tissue"] == "Thyroid" and by_name["anchor"]["matched_by_name"]

    class FakeLevels:
        def match_hpa(self, name, biosample_type):
            return {"source": "hpa_single_cell", "cell_type": "t-cells", "unit": "nCPM", "label": "HPA t-cells", "matched_by_name": True} if name == "T-cell" else None

    cell = tissue_spec("CL:0000084", {"name": "T-cell", "type": "primary_cell"}, FakeLevels())
    assert cell["anchor"]["cell_type"] == "t-cells" and cell["transcripts"] is POOLED_SHARES
    nothing = tissue_spec("CL:9", {"name": "odd cell", "type": "primary_cell"}, FakeLevels())
    assert nothing["anchor"] is None and nothing["label"] == "odd cell"


def test_hpa_match_by_name_scans_once(tmp_path):
    remote = BytesRemote(_hpa_zip([("ENSG0001", "GENE1", "t-cells", "42"), ("ENSG0001", "GENE1", "melanocytes", "7")]))
    levels = ReferenceLevels(remote, None, tmp_path)
    anchor = levels.match_hpa("T-cell", "primary_cell")
    assert anchor["cell_type"] == "t-cells" and anchor["matched_by_name"]
    assert levels.level(anchor, "GENE1", "ENSG0001.2")["value"] == 42.0
    assert levels.match_hpa("odd cell", "primary_cell") is None and levels.match_hpa("liver", "tissue") is None
    assert levels.level(None, "GENE1", "ENSG0001")["value"] is None


def test_http_products_tissues_and_report_for_any_tissue(api):
    from genomics.workflows.alphagenome.catalog import CATALOG_VERSION, save_catalog

    catalog = {"version": CATALOG_VERSION, "ontologies": [{"curie": "EFO:0001187", "name": "HepG2"}], "outputs": {},
               "tracks": {"RNA_SEQ": [{"ontology_curie": "EFO:0001187", "biosample_name": "HepG2", "biosample_type": "cell_line", "strand": "+"},
                                      {"ontology_curie": "UBERON:0004264", "biosample_name": "lower leg skin", "biosample_type": "tissue", "gtex_tissue": "Skin_Sun_Exposed_Lower_leg"}],
                          "SPLICE_JUNCTIONS": [{"ontology_curie": "EFO:0001187"}]}}
    rows = api.tissues(catalog, save_catalog)
    by = {r["curie"]: r for r in rows["tissues"]}
    assert rows["catalog"] and by["EFO:0001187"]["baseline"] == "hpa_cell_line" and by["EFO:0001187"]["junctions"] == 1
    assert by["UBERON:0004264"]["baseline"] == "gtex" and by["CL:1000458"]["default"]
    report = api("/products/report", gene="GENE1", sample="S1", tissue="EFO:0001187")
    assert report["available"] and report["tissue"]["ontology"] == "EFO:0001187" and report["tissue"]["label"] == "HepG2"
    assert report["haplotypes"]["H1"]["gene_fold"] is None and report["expression_problems"]  # not stored, no AlphaGenome in tests
