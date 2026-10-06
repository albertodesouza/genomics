"""Genotype PCA and PC matching: site selection, PCA, greedy matching, the HTTP API and PC fields."""
import json

import numpy as np
import pytest

from test_visualizer import dataset_dir  # noqa: F401 (fixture)
from test_visualizer_scalars import server  # noqa: F401 (fixture)

from genomics.visualizer.ancestry import dosage_matrix, match_groups, pca_scores, select_sites
from genomics.visualizer.datasets import Dataset
from genomics.visualizer.genotypes import CohortGenotypes


def _geno(sites, n_samples):
    """sites: (pos, ref, alt, carrier haplotype indices)."""
    carriers = [np.asarray(sorted(c), np.int32) for *_, c in sites]
    indptr = np.concatenate([[0], np.cumsum([c.size for c in carriers])]).astype(np.int64)
    return CohortGenotypes(
        samples=[f"S{i}" for i in range(n_samples)], positions=np.asarray([s[0] for s in sites], np.int64),
        refs=[s[1] for s in sites], alts=[s[2] for s in sites], ids=["."] * len(sites), indptr=indptr,
        carriers=np.concatenate(carriers) if carriers else np.zeros(0, np.int32),
    )


def test_select_sites_filters_maf_type_multiallelic_and_spacing():
    geno = _geno([
        (100, "A", "G", [0, 1, 2]),        # MAF 3/10: kept
        (150, "C", "T", [3, 4, 5]),        # too close to 100 with spacing 100
        (300, "A", "AT", [0, 1, 2]),       # indel
        (400, "G", "A", [0]),              # MAF 0.1 < 0.2
        (500, "T", "C", [1, 2]),           # multi-allelic position
        (500, "T", "G", [3]),
        (700, "C", "G", [0, 1, 2, 3, 4, 5, 6, 7, 8]),  # MAF 1/10 (ALT is the major allele)
        (900, "A", "C", [0, 2, 4, 6]),     # kept
    ], n_samples=5)
    assert select_sites(geno, 0.2, 100) == [0, 7]
    assert select_sites(geno, 0.2, 0) == [0, 1, 7]
    np.testing.assert_array_equal(dosage_matrix(geno, [0, 7]), [[2, 1], [1, 1], [0, 1], [0, 1], [0, 0]])


def test_pca_separates_two_populations_and_matching_balances():
    rng = np.random.default_rng(0)
    # Two populations with different allele frequencies at 400 sites, plus a cline inside population 1.
    p1, p2 = rng.uniform(0.05, 0.5, 400), rng.uniform(0.5, 0.95, 400)
    pop = np.r_[np.zeros(60, int), np.ones(40, int)]
    freqs = np.where(pop[:, None] == 0, p1, p2)
    dosages = rng.binomial(2, freqs).astype(np.float32)
    result = pca_scores(dosages, 5)
    pc1 = result["scores"][:, 0]
    assert abs(pc1[pop == 0].mean() - pc1[pop == 1].mean()) > 4 * pc1.std() / 2
    assert result["explained"][0] > result["explained"][1] > 0 and result["n_sites"] == 400
    assert pca_scores(dosages, 5)["scores"][0, 0] == pc1[0]  # deterministic signs

    scores = np.array([[0.0], [1.0], [5.0], [0.1], [0.9], [9.0]], np.float32)
    matched = match_groups(scores, [0, 1, 2], [3, 4, 5], k=1, caliper=0.0)
    assert sorted((a, b) for a, b, _ in matched["pairs"]) == [(0, 3), (1, 4), (2, 5)]  # no caliper
    tight = match_groups(scores, [0, 1, 2], [3, 4, 5], k=1, caliper=0.1)  # caliper 0.1 x SD(=3.3) = 0.33
    assert sorted((a, b) for a, b, _ in tight["pairs"]) == [(0, 3), (1, 4)]
    assert abs(tight["balance_after"][0]) < abs(tight["balance_before"][0])
    with pytest.raises(ValueError):
        match_groups(scores, [], [3], k=1, caliper=0)


def _post(call, path, body):
    return call("POST", path, body=body)


def test_http_pca_fields_and_matching(server, tmp_path):  # noqa: F811
    call, base, app, ds = server
    status, windows = call("GET", f"{base}/ancestry/windows")
    coverage = {w["gene"]: w["coverage"] for w in windows["windows"]}
    assert status == 200 and windows["default"] == ["GENE1"] and coverage == {"EXTRA": 0.0, "GENE1": 1.0}

    # The fixture window has a single SNV: PCA reports that instead of crashing.
    status, error = call("GET", f"{base}/ancestry/pca", {"genes": "GENE1", "min_maf": 0.1, "spacing": 0})
    assert status >= 400 and "sites" in json.dumps(error)

    # Inject a PCA result (as the job would cache it) and use it for PC fields and matching.
    dataset = app.dataset(ds.id)
    params = app.ancestry.params(dataset, ["GENE1"], 0.1, 0, 2)
    scores = np.array([[1.0, 0.0], [-1.0, 0.5], [0.9, -0.5]], np.float32)
    app.ancestry.memory.put(app.ancestry.key(dataset, params), {
        "samples": np.array(["S1", "S2", "S3"]), "scores": scores, "explained": np.array([0.6, 0.2], np.float32),
        "n_sites": np.int64(10), "sites_per_window": np.array([10]),
    })
    query = {"genes": "GENE1", "min_maf": 0.1, "spacing": 0, "components": 2}
    status, pca = call("GET", f"{base}/ancestry/pca", query)
    assert status == 200 and pca["samples"] == ["S1", "S2", "S3"] and pca["scores"][1] == [-1.0, 0.5] and pca["n_sites"] == 10

    status, fields = _post(call, f"{base}/ancestry/fields", {"params": query, "count": 2})
    assert status == 200 and fields["fields"] == ["pc1", "pc2"]
    status, samples = call("GET", f"{base}/samples")
    col = {name: i for i, name in enumerate(samples["columns"])}
    assert [row[col["pc1"]] for row in samples["rows"]] == [1.0, -1.0, pytest.approx(0.9)]
    assert {f["name"]: f["kind"] for f in samples["fields"]}["pc2"] == "numeric"

    body = {"params": query, "field": "superpopulation", "group_a": ["EUR"], "group_b": ["AFR"], "k": 1, "caliper": 0, "name": "eur_afr"}
    status, match = _post(call, f"{base}/ancestry/match", body)
    assert status == 200 and match["pairs"] == 1 and match["label_a"] == "EUR" and match["group_b"] == 2
    status, samples = call("GET", f"{base}/samples")
    col = {name: i for i, name in enumerate(samples["columns"])}
    assert {row[0]: row[col["eur_afr"]] for row in samples["rows"]} == {"S1": None, "S2": "EUR", "S3": "AFR"}  # S3 is nearer S2 on PC1
    assert _post(call, f"{base}/ancestry/match", body)[0] == 400  # name taken
    assert _post(call, f"{base}/ancestry/match", dict(body, name="pc1"))[0] == 400  # a field of the PC set
    assert _post(call, f"{base}/ancestry/match", dict(body, name="x", group_b=["EUR"]))[0] == 400  # overlap

    # Definitions survive a fresh dataset object and listing; removing the PCs drops the fields.
    app.catalog._datasets[ds.id] = Dataset(ds.path)
    status, windows = call("GET", f"{base}/ancestry/windows")
    assert {d["name"]: d["kind"] for d in windows["derived"]} == {"pcs": "pca", "eur_afr": "match"}
    status, samples = call("GET", f"{base}/samples")
    assert {"pc1", "pc2", "eur_afr"} <= set(samples["columns"])
    assert _post(call, f"{base}/ancestry/fields", {"params": query, "count": 0})[0] == 200
    status, samples = call("GET", f"{base}/samples")
    assert "pc1" not in samples["columns"] and "eur_afr" in samples["columns"]
    assert _post(call, f"{base}/scalars/delete", {"name": "eur_afr"})[0] == 200
