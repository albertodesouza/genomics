"""Region scalars: per-sample values from region means, quantile bins, persistence and the HTTP API."""
import http.client
import json
import threading
import time
from urllib.parse import urlencode

import numpy as np
import pytest

from test_visualizer import L, _free_port, _write_reference_prediction, dataset_dir  # noqa: F401 (fixture)

from genomics.visualizer.datasets import Dataset
from genomics.visualizer.scalars import ScalarStore, quantile_bins, scalar_values


def _region():
    means = np.array([[[1.0, 4.0], [3.0, 4.0]], [[2.0, 0.0], [np.nan, 0.0]], [[np.nan, 1.0], [np.nan, 1.0]]], np.float32)
    return {"samples": np.array(["S1", "S2", "S3"]), "means": means, "reference": np.zeros(2, np.float32)}


def test_scalar_values_average_haplotypes_and_transform():
    spec = {"output": "rna_seq", "track": 0, "start": 10, "end": 20, "stat": "mean", "transform": "none"}
    assert scalar_values(spec, _region()) == {"S1": 2.0, "S2": 2.0}  # S2 has one haplotype; S3 none
    assert scalar_values(dict(spec, stat="sum"), _region()) == {"S1": 20.0, "S2": 20.0}  # mean x 10 bp
    assert scalar_values(dict(spec, track=1, transform="log2p1"), _region()) == pytest.approx({"S1": np.log2(5), "S2": 0.0, "S3": 1.0})


def test_quantile_bins():
    values = {f"s{i}": float(i) for i in range(9)}
    labels = quantile_bins(values, 3)["labels"]
    assert [labels[f"s{i}"] for i in range(9)] == ["Q1"] * 3 + ["Q2"] * 3 + ["Q3"] * 3
    tied = quantile_bins({"a": 1.0, "b": 1.0, "c": 1.0, "d": 5.0}, 4)  # ties collapse bins
    assert tied["labels"]["a"] == tied["labels"]["b"] == tied["labels"]["c"] != tied["labels"]["d"]
    assert quantile_bins(values, 0) == {"labels": {}, "edges": []}


def test_store_persists_per_dataset(dataset_dir, tmp_path):  # noqa: F811
    dataset = Dataset(dataset_dir)
    path = tmp_path / "config" / "scalars.json"
    ScalarStore(path).put(dataset, {"name": "a", "track": 0})
    ScalarStore(path).put(dataset, {"name": "a", "track": 1})  # replaces by name
    ScalarStore(path).put(dataset, {"name": "b", "track": 0})
    assert [(s["name"], s["track"]) for s in ScalarStore(path).list(dataset)] == [("a", 1), ("b", 0)]
    assert ScalarStore(path).delete(dataset, "a") and not ScalarStore(path).delete(dataset, "a")
    assert [s["name"] for s in ScalarStore(path).list(dataset)] == ["b"]


@pytest.fixture
def server(dataset_dir, tmp_path):  # noqa: F811
    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.server import Handler, Server, VisualizerApp

    _write_reference_prediction(dataset_dir)
    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    store = ScalarStore(tmp_path / "scalars.json")
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=256 << 20, workers=2, runs_roots=[], remote=False, scalar_store=store)
    handler = type("ScalarHandler", (Handler,), {"app": app})
    port = _free_port()
    srv = Server(("127.0.0.1", port), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()

    def call(method, path, params=None, body=None):
        for _ in range(200):  # poll background jobs
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
            payload = json.dumps(body).encode("utf-8") if body is not None else None
            conn.request(method, f"{path}?{urlencode(params or {})}", body=payload, headers={"Content-Type": "application/json"} if payload else {})
            res = conn.getresponse()
            data = json.loads(res.read() or b"null")
            conn.close()
            if res.status != 200 or not (isinstance(data, dict) and data.get("pending")):
                return res.status, data
            time.sleep(0.05)
        raise AssertionError("job did not finish")

    try:
        yield call, f"/api/d/{ds.id}", app, ds
    finally:
        srv.shutdown()
        srv.server_close()
        app.shutdown()


def test_http_scalar_becomes_sample_fields(server):
    call, base, app, ds = server
    spec = {"name": "rna_b", "gene": "GENE1", "output": "rna_seq", "track": 1, "start": 0, "end": L, "stat": "mean", "bins": 3}
    status, saved = call("POST", f"{base}/scalars", body=spec)
    assert status == 200 and saved["scalar"]["ready"] and saved["scalar"]["histogram"]["n"] == 3

    status, samples = call("GET", f"{base}/samples")
    col = {name: i for i, name in enumerate(samples["columns"])}
    values = {row[col["sample_id"]]: (row[col["rna_b"]], row[col["rna_b_bin"]]) for row in samples["rows"]}
    assert values == {"S1": (1.0, "Q1"), "S2": (2.0, "Q2"), "S3": (3.0, "Q3")}  # track 1 is a per-sample constant
    kinds = {f["name"]: f["kind"] for f in samples["fields"]}
    assert kinds["rna_b"] == "numeric" and kinds["rna_b_bin"] == "categorical"

    # The bin field filters the cohort and splits group means like any categorical field.
    status, groups = call("GET", f"{base}/groups", {"gene": "GENE1", "output": "rna_seq", "tracks": "1", "coords": "reference", "start": 0, "end": L, "bins": 10, "field": "rna_b_bin"})
    assert status == 200 and {g["group"] for g in groups["groups"]} == {"Q1", "Q2", "Q3"}

    assert call("POST", f"{base}/scalars", body=spec)[0] == 400  # duplicate name
    assert call("POST", f"{base}/scalars", body=dict(spec, name="superpopulation"))[0] == 400  # existing field
    assert call("POST", f"{base}/scalars", body=dict(spec, name="x", end=L + 5))[0] == 400  # outside the window
    status, replaced = call("POST", f"{base}/scalars", body=dict(spec, stat="sum", bins=0, replace=True))
    assert status == 200 and replaced["scalar"]["histogram"]["max"] == pytest.approx(3.0 * L)
    status, samples = call("GET", f"{base}/samples")
    assert "rna_b_bin" not in samples["columns"]

    # A fresh dataset object (e.g. after a refresh) gets the stored scalar back from the cache.
    app.catalog._datasets[ds.id] = Dataset(ds.path)
    status, listed = call("GET", f"{base}/scalars")
    assert [s["name"] for s in listed["scalars"]] == ["rna_b"] and listed["scalars"][0]["ready"]

    assert call("POST", f"{base}/scalars/delete", body={"name": "rna_b"})[0] == 200
    status, samples = call("GET", f"{base}/samples")
    assert "rna_b" not in samples["columns"]
    assert call("POST", f"{base}/scalars/delete", body={"name": "rna_b"})[0] == 404
