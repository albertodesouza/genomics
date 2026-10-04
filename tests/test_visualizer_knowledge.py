"""Tests for the visualizer's links to public data: the bigWig reader, observed-data lookups
(FANTOM5 / ENCODE), the HGNC / Gene Ontology knowledge base and the offline remote cache.

No network: bigWig files are written locally by a minimal writer, FANTOM5 / ENCODE / HGNC / QuickGO
responses are small fixtures served by a fake :class:`RemoteCache`.
"""
import http.server
import json
import struct
import threading
import zlib
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

import numpy as np
import pytest

from genomics.visualizer.bigwig import BigWig
from genomics.visualizer.cache import DiskArrayCache
from genomics.visualizer.knowledge import HgncTable, KnowledgeBase
from genomics.visualizer.observed import Fantom5Index, ObservedService, ObservedUnavailable, pick_encode_file
from genomics.visualizer.remote import RemoteCache, RemoteError


# ------------------------------------------------------------------------------------ bigWig
def write_bigwig(path: Path, chroms, sections, two_level: bool = False) -> None:
    """Minimal bigWig writer: ``sections`` = [(chrom, kind, start, end, step, span, items)], one
    zlib-compressed data block each; kind 1 items are (start, end, value), 2 (start, value), 3 values."""
    names = list(chroms)
    key_size = max(len(n) for n in names)
    header_size, tree_header = 64, 32
    chrom_tree = struct.pack("<IIIIQQ", 0x78CA8C91, len(names), key_size, 8, len(names), 0)
    chrom_tree += struct.pack("<BBH", 1, 0, len(names))
    for i, name in enumerate(names):
        chrom_tree += name.encode().ljust(key_size, b"\0") + struct.pack("<II", i, chroms[name])
    data_offset = header_size + len(chrom_tree)
    blocks, leaves = b"", []
    offset = data_offset + 8
    for chrom, kind, start, end, step, span, items in sections:
        cid = names.index(chrom)
        body = struct.pack("<IIIIIBBH", cid, start, end, step, span, kind, 0, len(items))
        for item in items:
            body += struct.pack({1: "<IIf", 2: "<If", 3: "<f"}[kind], *(item if isinstance(item, tuple) else (item,)))
        block = zlib.compress(body)
        leaves.append((cid, start, cid, end, offset, len(block)))
        blocks += block
        offset += len(block)
    index_offset = offset
    rtree = struct.pack("<IIQIIIIQII", 0x2468ACE0, 256, len(leaves), leaves[0][0], leaves[0][1], leaves[-1][2], leaves[-1][3], index_offset, 1, 0)
    leaf_node = struct.pack("<BBH", 1, 0, len(leaves)) + b"".join(struct.pack("<IIIIQQ", *leaf) for leaf in leaves)
    if two_level:
        child = index_offset + 48 + 4 + 24
        rtree += struct.pack("<BBH", 0, 0, 1) + struct.pack("<IIIIQ", leaves[0][0], leaves[0][1], leaves[-1][2], leaves[-1][3], child)
    rtree += leaf_node
    header = struct.pack("<IHHQQQHHQQIQ", 0x888FFC26, 4, 0, header_size, data_offset, index_offset, 0, 0, 0, 0, 1 << 15, 0)
    path.write_bytes(header + chrom_tree + struct.pack("<Q", len(leaves)) + blocks + rtree)


@pytest.fixture
def bigwig_file(tmp_path):
    path = tmp_path / "x.bw"
    write_bigwig(path, {"chr1": 1000, "chr2": 500}, [
        ("chr1", 1, 10, 30, 0, 0, [(10, 12, 1.5), (20, 30, -2.0)]),
        ("chr1", 2, 100, 106, 0, 3, [(100, 4.0), (103, 5.0)]),
        ("chr1", 3, 200, 204, 2, 1, [7.0, 8.0]),
        ("chr2", 1, 0, 5, 0, 0, [(0, 5, 9.0)]),
    ], two_level=True)
    return path


def test_bigwig_reads_bedgraph_varstep_fixedstep_sections(bigwig_file):
    bw = BigWig(bigwig_file)
    assert bw.chroms == {"chr1": (0, 1000), "chr2": (1, 500)}
    values = bw.values("chr1", 0, 300)
    expected = np.zeros(300, np.float32)
    expected[10:12] = 1.5
    expected[20:30] = -2.0
    expected[100:103] = 4.0
    expected[103:106] = 5.0
    expected[200] = 7.0
    expected[202] = 8.0
    np.testing.assert_array_equal(values, expected)
    np.testing.assert_array_equal(bw.values("chr1", 25, 105), expected[25:105])  # partial overlaps are clipped
    np.testing.assert_array_equal(bw.values("2", 3, 8), [9, 9, 0, 0, 0])  # 2 -> chr2
    assert np.isnan(bw.values("chrX", 0, 3, missing=np.nan)).all()


def test_bigwig_over_http_without_range_support(bigwig_file):
    """SimpleHTTPRequestHandler ignores Range: the reader slices the full response."""
    handler = type("H", (http.server.SimpleHTTPRequestHandler,), {"log_message": lambda *a: None})
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), lambda *a, **k: handler(*a, directory=str(bigwig_file.parent), **k))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        bw = BigWig(f"http://127.0.0.1:{server.server_address[1]}/{bigwig_file.name}")
        np.testing.assert_array_equal(BigWig(bigwig_file).values("chr1", 0, 300), bw.values("chr1", 0, 300))
    finally:
        server.shutdown()
        server.server_close()


# ----------------------------------------------------------------------------- fake remote
class FakeRemote:
    """Serves canned responses by URL prefix; records requests."""

    def __init__(self, responses):
        self.responses = responses
        self.calls = []
        self.offline = False

    def _find(self, url, params):
        self.calls.append((url, params))
        for prefix, value in self.responses.items():
            if url.startswith(prefix):
                return value(params) if callable(value) else value
        raise RemoteError(f"no fixture for {url}")

    def get_text(self, url, ttl=0, params=None, headers=None):
        return self._find(url, params)

    def get_json(self, url, ttl=0, params=None):
        return self._find(url, params)


OBO = """format-version: 1.2

[Term]
id: FF:0000090
name: human light melanocyte sample
relationship: derives_from CL:0002567 ! light melanocyte

[Term]
id: FF:0000402
name: human K562 cell line sample
relationship: derives_from CLO:0007050

[Term]
id: FF:11274-116H5
name: Melanocyte - light, donor1
is_a: FF:0000090 ! human light melanocyte sample

[Term]
id: FF:11351-117H1
name: Melanocyte - light, donor2
is_a: FF:0000090 ! human light melanocyte sample

[Term]
id: FF:10824-111C5
name: K562 response to hemin, 3hr
is_a: FF:0000402

[Term]
id: FF:10454-106G9
name: Chronic myelogenous leukemia cell line:K562
is_a: FF:0000402
"""


def _trackdb(entries):
    blocks = []
    for ff, name, tech, files in entries:
        for scale, strand, path in files:
            blocks.append(f"""    track {name}_{scale}_{strand}
    longLabel {name}
    bigDataUrl {path}
    type bigWig
    metadata ontology_id={ff} sequence_tech={tech}
    subGroups sequenceTech={tech} category=primaryCell strand={'forward' if strand == '+' else 'reverse'}
""")
    return "\n".join(blocks)


def test_fantom5_index_matches_ontology_terms_and_cell_line_names():
    url = lambda name, scale, s: f"http://fantom.gsc.riken.jp/5/datahub/hg38/{'ctss' if scale == 'counts' else 'tpm'}/x/{quote(name)}.CNhs1{len(name)}.hg38.{scale}.{'fwd' if s == '+' else 'rev'}.bw"  # noqa: E731
    entries = [
        ("11274-116H5", "Melanocyte - light, donor1", "hCAGE", [(sc, s, url("Melanocyte - light, donor1", sc, s)) for sc in ("counts", "tpm") for s in "+-"]),
        ("11351-117H1", "Melanocyte - light, donor2", "hCAGE", [("tpm", "+", url("Melanocyte - light, donor2", "tpm", "+"))]),
        ("10824-111C5", "K562 response to hemin, 3hr", "LQhCAGE", [("tpm", "+", url("K562 response to hemin, 3hr", "tpm", "+"))]),
        ("10454-106G9", "Chronic myelogenous leukemia cell line:K562", "hCAGE", [("tpm", "+", url("Chronic myelogenous leukemia cell line:K562", "tpm", "+"))]),
    ]
    index = Fantom5Index(OBO, _trackdb(entries))
    libs, how = index.match("CL:0002567")
    assert how == "ontology" and [lib["ff"] for lib in libs] == ["11274-116H5", "11351-117H1"]
    assert libs[0]["urls"]["tpm"]["-"].startswith("https://")  # upgraded to https
    assert index.match("EFO:0002067", "K562", cell_line=False) == ([], "")
    libs, how = index.match("EFO:0002067", "K562", cell_line=True)
    assert how == "name" and libs[0]["name"].endswith("cell line:K562")  # untreated before time courses
    assert libs[0]["name"] == "Chronic myelogenous leukemia cell line:K562" and len(libs) == 2


def test_observed_fantom5_end_to_end_with_local_bigwigs(tmp_path):
    window = SimpleNamespace(chromosome="chr1", start=101, length=50)  # 1-based start: offset 0 = 0-based 100
    files = {}
    for donor, value in (("d1", 2.0), ("d2", 4.0)):
        path = tmp_path / "tpm" / f"{donor}.CNhs1{value:.0f}.hg38.tpm.fwd.bw"  # data hub layout: scale folder, strand suffix
        path.parent.mkdir(exist_ok=True)
        write_bigwig(path, {"chr1": 1000}, [("chr1", 1, 110, 112, 0, 0, [(110, 112, value)])])
        files[donor] = str(path)
    trackdb = _trackdb([
        ("11274-116H5", "Melanocyte - light, donor1", "hCAGE", [("tpm", "+", files["d1"])]),
        ("11351-117H1", "Melanocyte - light, donor2", "hCAGE", [("tpm", "+", files["d2"])]),
    ])
    remote = FakeRemote({"https://fantom.gsc.riken.jp/5/datafiles": OBO, "https://fantom.gsc.riken.jp/5/datahub": trackdb})
    service = ObservedService(remote, DiskArrayCache(tmp_path / "cache"))
    dataset = SimpleNamespace(id="ds", window=lambda gene: window)
    meta = {"ontology_curie": "CL:0002567", "strand": "+", "Assay title": "hCAGE", "data_source": "fantom"}
    result = service.load(dataset, "G", "cage", 1, meta, "tpm")
    assert result["info"]["provider"] == "FANTOM5" and result["info"]["sources_loaded"] == 2
    expected = np.zeros(50, np.float32)
    expected[10:12] = 3.0  # mean of the two donors at 0-based 110-111
    np.testing.assert_array_equal(result["values"], expected)
    assert service.cached(dataset, "G", "cage", 1, "tpm") is result
    assert list((tmp_path / "cache" / "observed").glob("*.npz"))  # per-file windows cached on disk
    missing = service.load(dataset, "G", "cage", 0, {"ontology_curie": "CL:9999999", "strand": "+"}, "tpm")
    assert missing["values"] is None and "No FANTOM5" in missing["info"]["reason"]


def test_observed_sources_for_encode_and_unsupported_tracks():
    def search(params):
        assert params["biosample_ontology.term_id"] == "CL:1000458" and params["assay_title"] == "total RNA-seq"
        files = [
            {"accession": "F1", "href": "/files/F1/@@download/F1.bigWig", "output_type": "plus strand signal of all reads", "file_format": "bigWig", "assembly": "GRCh38", "status": "released"},
            {"accession": "F2", "href": "/files/F2/@@download/F2.bigWig", "output_type": "plus strand signal of unique reads", "file_format": "bigWig", "assembly": "GRCh38", "status": "released", "preferred_default": True},
            {"accession": "F3", "href": "/files/F3/@@download/F3.bigWig", "output_type": "minus strand signal of unique reads", "file_format": "bigWig", "assembly": "GRCh38", "status": "released"},
            {"accession": "F4", "href": "/files/F4/@@download/F4.bigWig", "output_type": "plus strand signal of unique reads", "file_format": "bigWig", "assembly": "hg19", "status": "released"},
        ]
        return {"@graph": [{"accession": "ENCSR1", "biosample_summary": "melanocyte", "date_released": "2020", "files": files}]}

    service = ObservedService(FakeRemote({"https://www.encodeproject.org/search/": search}), DiskArrayCache(None))
    meta = {"ontology_curie": "CL:1000458", "strand": "+", "Assay title": "total RNA-seq", "data_source": "encode"}
    info = service.sources("rna_seq", meta)
    assert info["provider"] == "ENCODE" and [i["file"] for i in info["items"]] == ["https://www.encodeproject.org/files/F2/@@download/F2.bigWig"]
    assert "biosample_ontology.term_id=CL%3A1000458" in info["query"]
    with pytest.raises(ObservedUnavailable, match="GTEx"):
        service.sources("rna_seq", {**meta, "data_source": "gtex"})
    with pytest.raises(ObservedUnavailable):
        service.sources("splice_sites", {"strand": "+"})


def test_pick_encode_file_prefers_strand_output_type_and_pooled_files():
    files = [
        {"output_type": "fold change over control", "file_format": "bigWig", "assembly": "GRCh38", "status": "released", "href": "/a", "biological_replicates": [1]},
        {"output_type": "fold change over control", "file_format": "bigWig", "assembly": "GRCh38", "status": "released", "href": "/b", "biological_replicates": [1, 2]},
        {"output_type": "signal p-value", "file_format": "bigWig", "assembly": "GRCh38", "status": "released", "href": "/c", "preferred_default": True},
    ]
    assert pick_encode_file(files, None)["href"] == "/b"
    assert pick_encode_file(files, "+")["href"] == "/b"  # unstranded library for a stranded track
    assert pick_encode_file([{**files[0], "status": "revoked"}], None) is None


# ------------------------------------------------------------------------------ knowledge
HGNC = "\t".join(["hgnc_id", "symbol", "name", "locus_group", "locus_type", "status", "location", "alias_symbol", "prev_symbol", "gene_group", "gene_group_id", "entrez_id", "ensembl_gene_id", "uniprot_ids", "omim_id"]) + "\n" + "\n".join([
    "\t".join(["HGNC:12450", "TYRP1", "tyrosinase related protein 1", "protein-coding gene", "gene with protein product", "Approved", "9p23", '"GP75|OCA3"', '"TYRP|CAS2"', '"Tyrosinase family"', "3506", "7306", "ENSG00000107165", "P17643", "115501"]),
    "\t".join(["HGNC:12442", "TYR", "tyrosinase", "protein-coding gene", "gene with protein product", "Approved", "11q14.3", "OCA1A", "", '"Tyrosinase family"', "3506", "7299", "ENSG00000077498", "P14679", "606933"]),
    "\t".join(["HGNC:99999", "OLD1", "withdrawn gene", "other", "unknown", "Entry Withdrawn", "", "", "", "", "", "", "", "", ""]),
]) + "\n"


def test_hgnc_table_resolves_symbols_aliases_ids_and_groups():
    table = HgncTable(HGNC)
    assert len(table) == 2  # withdrawn entries are dropped
    for text in ("TYRP1", "tyrp1", "OCA3", "TYRP", "HGNC:12450", "ENSG00000107165.16"):
        assert table.resolve(text) == "TYRP1", text
    hits = table.search("OCA3")
    assert hits[0][0]["symbol"] == "TYRP1" and hits[0][1] == "alias OCA3"
    assert [rec["symbol"] for rec, _ in table.search("tyrosinase")] == ["TYR", "TYRP1"]
    assert table.search_groups("tyrosin") == [{"name": "Tyrosinase family", "id": "3506", "size": 2}]
    assert table.get("OCA3")["uniprot_ids"] == ["P17643"]


def test_knowledge_go_genes_and_gene_annotations_from_quickgo():
    go_tsv = "SYMBOL\tQUALIFIER\tGO TERM\tGO NAME\nTYRP1\tinvolved_in\tGO:0042438\tmelanin biosynthetic process\nP99999\tinvolved_in\tGO:0042438\tx\nTYR\tNOT|involved_in\tGO:0042438\tx\nOCA3\tinvolved_in\tGO:0042438\tmelanin biosynthetic process\n"
    term = {"results": [{"id": "GO:0042438", "name": "melanin biosynthetic process", "aspect": "biological_process", "definition": {"text": "def"}}]}
    annotations = {"pageInfo": {"total": 1}, "results": [
        {"goId": "GO:0042438", "goName": "melanin biosynthetic process", "goAspect": "biological_process", "qualifier": "involved_in", "goEvidence": "TAS"},
        {"goId": "GO:0042438", "goName": "melanin biosynthetic process", "goAspect": "biological_process", "qualifier": "involved_in", "goEvidence": "IEA"},
        {"goId": "GO:0033162", "goName": "melanosome membrane", "goAspect": "cellular_component", "qualifier": "NOT|located_in", "goEvidence": "IBA"},
    ]}
    remote = FakeRemote({
        "https://storage.googleapis.com/": HGNC,
        "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch": go_tsv,
        "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/": term,
        "https://www.ebi.ac.uk/QuickGO/services/annotation/search": annotations,
    })
    kb = KnowledgeBase(remote)
    result = kb.go_genes("go:0042438")
    assert [g["symbol"] for g in result["genes"]] == ["TYRP1"]  # NOT annotations, unknown symbols dropped; aliases merged
    assert "regulates" not in result["relations"] and "regulates" in kb.go_genes("GO:0042438", include_regulation=True)["relations"]
    go = kb.gene_go("OCA3")
    bp = next(a for a in go["aspects"] if a["aspect"] == "biological_process")["terms"][0]
    assert bp["id"] == "GO:0042438" and bp["evidence"] == ["IEA", "TAS"] and not bp["negated"]
    assert next(a for a in go["aspects"] if a["aspect"] == "cellular_component")["terms"][0]["negated"]
    with pytest.raises(ValueError):
        kb.go_genes("melanin")
    with pytest.raises(KeyError):
        kb.gene_go("NOT_A_GENE")


def test_remote_cache_serves_stale_copies_and_respects_offline(tmp_path, monkeypatch):
    cache = RemoteCache(tmp_path)
    calls = []

    def fake_download(url, headers):
        calls.append(url)
        if len(calls) > 1:
            raise RemoteError("down")
        return b'{"ok": 1}'

    monkeypatch.setattr(cache, "_download", fake_download)
    assert cache.get_json("https://example.org/x", params={"a": 1}) == {"ok": 1}
    assert cache.get_json("https://example.org/x", ttl=0, params={"a": 1}) == {"ok": 1}  # expired + unreachable: stale copy
    offline = RemoteCache(tmp_path, offline=True)
    assert offline.get_json("https://example.org/x", params={"a": 1}) == {"ok": 1}
    with pytest.raises(RemoteError, match="--no-remote"):
        offline.get_json("https://example.org/other")
