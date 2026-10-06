"""Variant page backend: cohort genotypes, AlphaGenome effect by genotype, GTEx client, HTTP API."""
import http.client
import json
import threading
import time
from urllib.parse import urlencode

import numpy as np
import pytest

from test_visualizer import L, REF, _free_port, _write_reference_prediction, dataset_dir  # noqa: F401 (fixture)

from genomics.visualizer.cache import DiskArrayCache
from genomics.visualizer.datasets import Dataset
from genomics.visualizer.genotypes import GenotypeService, regress_on_dosage, t_test_p
from genomics.visualizer.gtex import GtexClient, GtexError
from genomics.visualizer.remote import RemoteError
from genomics.visualizer.signals import SignalService

INS = (1031, REF[30], REF[30] + "AAA")  # S1 1|1, S2 0|0, S3 1|0 in the fixture VCFs
SNV = (1011, REF[10], "T")  # S1 1|0, S2 0|1, S3 0|0


def _service(tmp_path):
    signals = SignalService(64 << 20, DiskArrayCache(tmp_path / "cache"))
    return GenotypeService(signals, DiskArrayCache(tmp_path / "cache"))


def test_cohort_genotypes_from_window_vcfs(dataset_dir, tmp_path):  # noqa: F811
    dataset = Dataset(dataset_dir)
    geno = _service(tmp_path).genotypes(dataset, "GENE1")
    assert geno.samples == ["S1", "S2", "S3"]
    site = geno.find(*INS)
    assert geno.dosage(site).tolist() == [2, 0, 1]
    assert geno.haplotype_alleles(site).tolist() == [[1, 1], [0, 0], [1, 0]]
    snv = geno.find(*SNV)
    assert geno.dosage(snv).tolist() == [1, 1, 0]
    assert geno.haplotype_alleles(snv).tolist() == [[1, 0], [0, 1], [0, 0]]
    with pytest.raises(KeyError):
        geno.find(1011, REF[10], "G")
    assert list(geno.positions) == sorted(geno.positions)


def test_cohort_genotypes_are_cached_on_disk(dataset_dir, tmp_path, monkeypatch):  # noqa: F811
    dataset = Dataset(dataset_dir)
    first = _service(tmp_path).genotypes(dataset, "GENE1")
    from genomics.visualizer import genotypes as module

    monkeypatch.setattr(module, "parse_window_vcf", lambda *a, **k: (_ for _ in ()).throw(AssertionError("VCF re-read")))
    again = _service(tmp_path).genotypes(dataset, "GENE1")  # fresh service: memory empty, disk hit
    assert again.n_sites == first.n_sites and again.dosage(again.find(*INS)).tolist() == [2, 0, 1]


def test_regression_and_t_distribution():
    # Reference values from scipy.stats.t.sf (two-sided).
    assert t_test_p(2.0, 10) == pytest.approx(0.0733880347707, rel=1e-9)
    assert t_test_p(-3.5, 600) == pytest.approx(0.000499765454419, rel=1e-9)
    assert t_test_p(0.5, 3) == pytest.approx(0.651447964848151, rel=1e-9)
    r = regress_on_dosage(np.array([1.0, 2.0, 3.0, np.nan]), np.array([2, 0, 1, 1]))
    assert r["n"] == 3 and r["slope"] == pytest.approx(-0.5) and r["intercept"] == pytest.approx(2.5)
    assert np.isnan(regress_on_dosage(np.array([1.0, 2.0]), np.array([0, 1]))["p"])  # too few samples
    assert np.isnan(regress_on_dosage(np.array([1.0, 2.0, 3.0]), np.array([1, 1, 1]))["slope"])  # monomorphic


def test_effect_by_genotype_over_a_region(dataset_dir, tmp_path):  # noqa: F811
    _write_reference_prediction(dataset_dir)
    dataset = Dataset(dataset_dir)
    service = _service(tmp_path)
    geno = service.genotypes(dataset, "GENE1")
    region = service.region_means(dataset, "GENE1", "rna_seq", 0, L)
    assert region["means"].shape == (3, 2, 2)
    np.testing.assert_allclose(region["means"][:, :, 1], [[1, 1], [2, 2], [3, 3]])  # track 1 = per-sample constant
    # The reference prediction has track 1 ("cell B -", constant 7) but not track 0.
    assert region["reference"][1] == pytest.approx(7.0) and np.isnan(region["reference"][0])
    effect = service.effect(geno, geno.find(*INS), region, 1, ["S1", "S2", "S3"])
    assert effect["regression"]["slope"] == pytest.approx(-0.5)  # dosages 2, 0, 1 -> values 1, 2, 3
    assert [g["summary"]["n"] for g in effect["groups"]] == [1, 1, 1]
    assert effect["groups"][0]["samples"] == ["S2"] and effect["groups"][2]["values"] == [1.0]
    # ALT haplotypes: S1 H1, S1 H2, S3 H1 -> (1 + 1 + 3) / 3; REF: S2 H1, S2 H2, S3 H2 -> (2 + 2 + 3) / 3
    assert effect["haplotypes"]["log2fc"] == pytest.approx(np.log2((5 / 3) / (7 / 3)))
    assert effect["direction"] == -1
    # A cohort subset only uses its samples.
    assert service.effect(geno, geno.find(*INS), region, 1, ["S1", "S3"])["n_samples"] == 2


# ----------------------------------------------------------------------------- GTEx client
class FakeRemote:
    """Answers GtexClient's requests from a dict of path -> JSON (or an exception)."""

    def __init__(self, answers):
        self.answers = answers
        self.calls = []

    def get_json(self, url, ttl=None, params=None):
        path = url.split("/api/v2/", 1)[1]
        self.calls.append((path, params))
        key = path if path in self.answers else f"{path}?{params.get('tissueSiteDetailId', '')}"
        answer = self.answers.get(key, self.answers.get(path))
        if answer is None:  # like GTEx: unknown variant / gene / gene x tissue pair
            raise RemoteError("gtexportal.org unreachable: HTTP Error 400: Bad Request")
        if isinstance(answer, Exception):
            raise answer
        return answer


GTEX_ANSWERS = {
    "dataset/variant": {"data": [{"snpId": "rs1", "variantId": "chr1_1031_A_AAAA_b38", "chromosome": "chr1", "pos": 1031, "ref": "A", "alt": "AAAA"}]},
    "association/singleTissueEqtl": {"data": [
        {"geneSymbol": "GENE1", "gencodeId": "ENSG1.1", "tissueSiteDetailId": "Skin_Sun_Exposed_Lower_leg", "nes": -0.2, "pValue": 1e-5},
        {"geneSymbol": "OTHER", "gencodeId": "ENSG2.1", "tissueSiteDetailId": "Whole_Blood", "nes": 0.4, "pValue": 1e-9}]},
    "reference/gene": {"data": [{"geneSymbol": "GENE1", "gencodeId": "ENSG1.1"}]},
    "dataset/tissueSiteDetail": {"data": [{"tissueSiteDetailId": "Skin_Sun_Exposed_Lower_leg", "tissueSiteDetail": "Skin - Sun Exposed (Lower leg)", "colorHex": "0000FF"},
                                          {"tissueSiteDetailId": "Whole_Blood", "tissueSiteDetail": "Whole Blood"}]},
    "association/dyneqtl?Skin_Sun_Exposed_Lower_leg": {"nes": -0.2, "pValue": 1e-5, "tStatistic": -4.4, "maf": 0.3, "homoRefCount": 2, "hetCount": 2, "homoAltCount": 1,
                                                       "data": [0.1, 0.2, -0.1, 0.0, -0.5], "genotypes": [0, 0, 1, 1, 2]},
    "association/dyneqtl?Whole_Blood": RemoteError("gtexportal.org unreachable: HTTP Error 400: Bad Request"),
}


def test_gtex_report_parses_and_flags_untested_tissues():
    client = GtexClient(FakeRemote(GTEX_ANSWERS))
    report = client.report("chr1_1031_A_AAAA_b38", "GENE1", ["Skin_Sun_Exposed_Lower_leg", "Whole_Blood"])
    assert report["in_gtex"] and report["rsid"] == "rs1"
    assert [r["gene"] for r in report["significant"]] == ["OTHER", "GENE1"]  # sorted by p
    skin, blood = report["dynamic"]
    assert skin["nes"] == -0.2 and skin["counts"] == {"0/0": 2, "0/1": 2, "1/1": 1}
    assert skin["values"] == {"0/0": [0.1, 0.2], "0/1": [-0.1, 0.0], "1/1": [-0.5]}
    assert "not tested" in blood["error"]


def test_gtex_distinguishes_not_found_from_unreachable():
    missing = GtexClient(FakeRemote({"dataset/variant": RemoteError("gtexportal.org unreachable: HTTP Error 400: Bad Request")}))
    assert missing.report("chr1_5_A_G_b38", "GENE1", []) == {"variant_id": "chr1_5_A_G_b38", "in_gtex": False, "rsid": None, "dataset": "gtex_v8"}
    offline = GtexClient(FakeRemote({"dataset/variant": RemoteError("gtexportal.org: not cached and remote lookups are disabled (--no-remote)")}))
    with pytest.raises(GtexError):
        offline.report("chr1_5_A_G_b38", "GENE1", [])


# ----------------------------------------------------------------------------- HTTP API
@pytest.fixture
def server(dataset_dir, tmp_path):  # noqa: F811
    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.server import Handler, Server, VisualizerApp

    _write_reference_prediction(dataset_dir)
    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=256 << 20, workers=2, runs_roots=[], remote=False)
    app.gtex = GtexClient(FakeRemote(GTEX_ANSWERS))
    handler = type("VariantHandler", (Handler,), {"app": app})
    port = _free_port()
    srv = Server(("127.0.0.1", port), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()

    def get(path, params=None):
        for _ in range(200):  # poll background jobs
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
            conn.request("GET", f"{path}?{urlencode(params or {})}")
            res = conn.getresponse()
            body = json.loads(res.read() or b"null")
            conn.close()
            if res.status != 200 or not (isinstance(body, dict) and body.get("pending")):
                return res.status, body
            time.sleep(0.05)
        raise AssertionError("job did not finish")

    try:
        yield get, f"/api/d/{ds.id}"
    finally:
        srv.shutdown()
        srv.server_close()
        app.shutdown()


def test_http_variant_site_effect_groups_and_gtex(server):
    get, base = server
    pos, ref, alt = INS
    status, sites = get(f"{base}/variant/sites", {"gene": "GENE1"})
    assert status == 200 and any(s["pos"] == pos and s["alt"] == alt for s in sites["sites"])
    assert {"variant_id", "af", "type"} <= set(sites["sites"][0])

    status, site = get(f"{base}/variant/site", {"gene": "GENE1", "pos": pos, "ref": ref, "alt": alt, "field": "superpopulation"})
    assert status == 200 and site["variant_id"] == f"chr1_{pos}_{ref}_{alt}_b38" and site["type"] == "insertion"
    assert site["cohort"]["counts"] == {"0/0": 1, "0/1": 1, "1/1": 1} and site["cohort"]["af"] == pytest.approx(0.5)
    assert {g["group"]: g["af"] for g in site["cohort"]["by_group"]} == {"AFR": pytest.approx(0.75), "EUR": 0.0}
    status, site = get(f"{base}/variant/site", {"gene": "GENE1", "pos": pos, "ref": ref, "alt": alt, "filters": json.dumps({"superpopulation": ["AFR"]})})
    assert site["cohort"]["cohort"] == 2

    status, effect = get(f"{base}/variant/effect", {"gene": "GENE1", "pos": pos, "ref": ref, "alt": alt, "output": "rna_seq", "track": 1, "start": 0, "end": L})
    assert status == 200 and effect["regression"]["slope"] == pytest.approx(-0.5) and effect["labels"] == ["0/0", "0/1", "1/1"]

    status, groups = get(f"{base}/groups", {"gene": "GENE1", "output": "rna_seq", "tracks": "1", "coords": "reference", "start": 0, "end": L, "bins": 10, "field": f"variant:{SNV[0]}:{SNV[1]}:{SNV[2]}"})
    assert status == 200
    by_name = {g["group"]: g for g in groups["groups"]}
    assert set(by_name) == {f"{SNV[1]}/{SNV[1]}", f"{SNV[1]}/T"} and by_name[f"{SNV[1]}/T"]["samples"] == 2

    assert get(f"{base}/variant/site", {"gene": "GENE1", "pos": pos, "ref": ref, "alt": "C"})[0] == 404
    assert get(f"{base}/variant/site", {"gene": "NOPE", "pos": pos, "ref": ref, "alt": alt})[0] == 404

    status, gtex = get("/api/gtex/variant", {"variant_id": f"chr1_{pos}_{ref}_{alt}_b38", "gene": "GENE1", "tissues": "Skin_Sun_Exposed_Lower_leg"})
    assert status == 200 and gtex["dynamic"][0]["nes"] == -0.2
    assert get("/api/gtex/variant", {"variant_id": "not-a-variant"})[0] == 400
    status, resolved = get("/api/gtex/resolve", {"rsid": "rs1"})
    assert status == 200 and resolved["pos"] == 1031
