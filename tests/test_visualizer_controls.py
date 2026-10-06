"""Matched control windows (crop statistic, 1:1 assignment, HTTP) and the notebook client."""
import json
import math

import numpy as np
import pytest

from test_visualizer import dataset_dir  # noqa: F401 (fixture)

from genomics.visualizer.cache import DiskArrayCache
from genomics.visualizer.controls import ControlError, ControlService, crop_signal, match_windows
from genomics.visualizer.datasets import Dataset
from genomics.visualizer.signals import SignalService

L = 64
CROP = 32  # the centred window the "model" reads: offsets 16..48
TERMS = ["CL:1", "CL:2"]


def _window(root, gene, per_track):
    """A window whose reference prediction holds a constant per track, plus one unwanted term."""
    ref = root / "references" / "windows" / gene
    (ref / "predictions_ref").mkdir(parents=True)
    (ref / "ref.window.fa").write_text(f">chr1:1-{L}\n{'ACGT' * (L // 4)}\n")
    values = np.stack([np.full(L, v, np.float32) for v in (*per_track, 1000.0)], axis=1)
    np.savez_compressed(ref / "predictions_ref" / "rna_seq.npz", values=values)
    meta = [{"ontology_curie": "CL:1", "biosample_name": "a", "strand": "+"},
            {"ontology_curie": "CL:2", "biosample_name": "b", "strand": "-"},
            {"ontology_curie": "CL:9", "biosample_name": "other", "strand": "+"}]
    (ref / "predictions_ref" / "rna_seq_metadata.json").write_text(json.dumps({"metadata": meta}))
    # The dataset's own track order is the same, so columns need no remapping.
    sample_dir = root / "individuals" / "S1" / "windows" / gene / "predictions_H1"
    sample_dir.mkdir(parents=True)
    np.savez_compressed(sample_dir / "rna_seq.npz", values=values)
    (sample_dir / "rna_seq_metadata.json").write_text(json.dumps({"metadata": meta}))


@pytest.fixture
def control_dataset(tmp_path):
    """Two panel windows and four candidates, with known signal totals."""
    root = tmp_path / "ds"
    root.mkdir()
    (root / "dataset_metadata.json").write_text(json.dumps({
        "dataset_name": "controls", "individuals": ["S1"],
        "individuals_pedigree": {"S1": {"superpopulation": "AFR"}},
        "genes": ["PANEL_A", "PANEL_B"], "alphagenome_outputs": ["RNA_SEQ"], "window_size": L,
    }))
    # Per-window totals over the crop: 32 rows x the summed track constants.
    for gene, tracks in [("PANEL_A", (1.0, 1.0)), ("PANEL_B", (50.0, 50.0)),
                         ("CTRL_NEAR_A", (1.0, 1.1)), ("CTRL_NEAR_B", (49.0, 49.0)),
                         ("CTRL_TINY", (0.01, 0.01)), ("CTRL_HUGE", (900.0, 900.0))]:
        _window(root, gene, tracks)
    return Dataset(root)


@pytest.fixture
def live_server(dataset_dir, tmp_path):  # noqa: F811
    """The real HTTP server on a port, for the notebook client to talk to."""
    import threading

    from test_visualizer import _free_port, _write_reference_prediction

    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.server import Handler, Server, VisualizerApp

    _write_reference_prediction(dataset_dir)
    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=64 << 20, workers=2, runs_roots=[], remote=False)
    srv = Server(("127.0.0.1", _free_port()), type("ClientHandler", (Handler,), {"app": app}))
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        yield f"http://127.0.0.1:{srv.server_address[1]}", ds.id
    finally:
        srv.shutdown()
        srv.server_close()
        app.shutdown()


@pytest.fixture
def service(tmp_path):
    return ControlService(SignalService(32 << 20, DiskArrayCache(tmp_path / "cache"), workers=1))


def test_crop_signal_sums_only_the_chosen_terms_inside_the_crop(control_dataset, service):
    # 32 crop rows x (1.0 + 1.0): the CL:9 track of 1000.0 and the rows outside the crop are left out.
    assert crop_signal(service.signals, control_dataset, "PANEL_A", "rna_seq", TERMS, CROP) == pytest.approx(64.0)
    assert crop_signal(service.signals, control_dataset, "PANEL_A", "rna_seq", ["CL:1"], CROP) == pytest.approx(32.0)
    assert crop_signal(service.signals, control_dataset, "PANEL_A", "rna_seq", TERMS, L) == pytest.approx(128.0)
    with pytest.raises(ControlError, match="ontology"):
        crop_signal(service.signals, control_dataset, "PANEL_A", "rna_seq", ["CL:404"], CROP)
    with pytest.raises(FileNotFoundError):
        crop_signal(service.signals, control_dataset, "PANEL_A", "cage", TERMS, CROP)


def test_match_windows_minimises_the_total_log_ratio():
    panel = {"P1": 10.0, "P2": 1000.0}
    pairs = match_windows(panel, {"C_LOW": 9.0, "C_MID": 100.0, "C_HIGH": 1100.0})
    assert [(p["gene"], p["control"]) for p in pairs] == [("P1", "C_LOW"), ("P2", "C_HIGH")]
    assert pairs[0]["ratio"] == pytest.approx(0.9) and pairs[1]["ratio"] == pytest.approx(1.1)

    # Order-preserving: the cheaper total wins even when a candidate is nearer one panel gene alone.
    crossed = match_windows({"P1": 10.0, "P2": 20.0}, {"A": 11.0, "B": 21.0})
    assert [(p["gene"], p["control"]) for p in crossed] == [("P1", "A"), ("P2", "B")]
    assert len({p["control"] for p in match_windows({"P1": 1.0, "P2": 1.0}, {"A": 1.0, "B": 1.0})}) == 2

    with pytest.raises(ControlError, match="at least as many"):
        match_windows(panel, {"only": 1.0})
    with pytest.raises(ControlError, match="No panel"):
        match_windows({}, {"a": 1.0})
    assert match_windows({"P": 0.0}, {"C": 0.0})[0]["ratio"] is None  # a zero total is not a crash


def test_matched_panel_picks_unlisted_windows_and_reports_the_pairing(control_dataset, service):
    assert service.candidates(control_dataset, ["PANEL_A"]) == ["CTRL_HUGE", "CTRL_NEAR_A", "CTRL_NEAR_B", "CTRL_TINY"]

    result = service.matched_panel(control_dataset, ["PANEL_A", "PANEL_B"], "rna_seq", TERMS, CROP)
    assert [(p["gene"], p["control"]) for p in result["pairs"]] == [("PANEL_A", "CTRL_NEAR_A"), ("PANEL_B", "CTRL_NEAR_B")]
    assert result["controls"] == ["CTRL_NEAR_A", "CTRL_NEAR_B"]
    assert result["panel_total"] == pytest.approx(64.0 + 3200.0)
    assert result["control_total"] == pytest.approx(67.2 + 3136.0)
    assert result["total_ratio"] == pytest.approx(3203.2 / 3264.0)
    assert result["mean_abs_log_ratio"] == pytest.approx((abs(math.log10(67.2 / 64)) + abs(math.log10(3136 / 3200))) / 2)
    assert result["candidates"] == 4 and result["skipped"] == {}

    # A panel window with no reference prediction is an error naming it, not a silent drop.
    (control_dataset.path / "references" / "windows" / "PANEL_A" / "predictions_ref" / "rna_seq.npz").unlink()
    fresh = Dataset(control_dataset.path)
    with pytest.raises(ControlError, match="PANEL_A"):
        service.matched_panel(fresh, ["PANEL_A", "PANEL_B"], "rna_seq", TERMS, CROP)
    with pytest.raises(ControlError, match="Unknown window"):
        service.matched_panel(fresh, ["NOPE"], "rna_seq", TERMS, CROP)
    with pytest.raises(ControlError, match="at least one"):
        service.matched_panel(fresh, [], "rna_seq", TERMS, CROP)


def test_matched_panel_needs_enough_candidates_with_predictions(control_dataset, service):
    for gene in ("CTRL_TINY", "CTRL_HUGE", "CTRL_NEAR_B"):
        (control_dataset.path / "references" / "windows" / gene / "predictions_ref" / "rna_seq.npz").unlink()
    with pytest.raises(ControlError, match="Only 1 control window"):
        service.matched_panel(Dataset(control_dataset.path), ["PANEL_A", "PANEL_B"], "rna_seq", TERMS, CROP)

    # Every window listed in the metadata is panel, so a panel covering them all has no candidates.
    with pytest.raises(ControlError, match="no windows outside"):
        service.matched_panel(control_dataset, [g for g in control_dataset.genes], "rna_seq", TERMS, CROP)


def test_client_reads_the_same_arrays_as_the_pages(live_server):
    """The notebook client against a live server: arrays, the cohort, and a clear error."""
    import numpy as np

    from genomics.visualizer.client import Visualizer, VisualizerError, _arrays

    url, dataset_id = live_server
    v = Visualizer(url, dataset=dataset_id)
    assert v.dataset == dataset_id and dataset_id in repr(v)
    assert "GENE1" in v.genes()

    samples = v.samples(as_frame=False)
    assert samples["columns"][0] == "sample_id" and len(samples["rows"]) == 3
    assert v.cohort({"superpopulation": ["AFR"]}) == ["S1", "S3"]
    assert v.cohort() == ["S1", "S2", "S3"]

    sig = v.signal("GENE1", "rna_seq", series=["S1:H1"], tracks=[0, 1], start=0, end=40, bins=20)
    assert sig["edges"].shape == (21,) and sig["series"][0]["mean"].shape == (2, 20)
    assert sig["series"][0]["label"] == "S1 H1" and sig["series"][0]["mean"].dtype == np.float32

    with pytest.raises(VisualizerError, match="Unknown gene"):
        v.signal("NOPE", "rna_seq", series=["S1:H1"], end=10)
    with pytest.raises(VisualizerError, match="no dataset|Unknown dataset"):
        Visualizer(url, dataset="not-a-dataset")
    with pytest.raises(VisualizerError, match="Cannot reach"):
        Visualizer("http://127.0.0.1:1", timeout=2)

    # float32 blobs decode to the right shape; integer lists named in ARRAY_KEYS become arrays.
    import base64
    blob = {"$f32": base64.b64encode(np.arange(6, dtype="<f4").tobytes()).decode(), "shape": [2, 3]}
    np.testing.assert_array_equal(_arrays({"mean": blob})["mean"], np.arange(6, dtype=np.float32).reshape(2, 3))
    np.testing.assert_array_equal(_arrays({"edges": [1, 2, None]})["edges"], [1.0, 2.0, np.nan])
    assert _arrays({"label": [1, 2]})["label"] == [1, 2]  # not an array key: left alone
