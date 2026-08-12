"""Tests LiveAlphaGenomePredictor/reorder_to_canonical against a fake alphagenome client, injected
into sys.modules so these tests don't require the real `alphagenome` package (not installed in
every environment that runs the test suite) or a real API key/network access.
"""
import sys
import types

import numpy as np
import pytest

from genomics.predictors.genotype_based.analysis.live_alphagenome_prediction import reorder_to_canonical


def test_reorder_to_canonical_permutes_columns():
    # 3 positions x 3 tracks, tracks currently ordered B, A, C.
    values = np.array(
        [
            [10.0, 20.0, 30.0],
            [11.0, 21.0, 31.0],
            [12.0, 22.0, 32.0],
        ]
    )
    metadata_records = [
        {"ontology_curie": "CL:B", "strand": "+"},
        {"ontology_curie": "CL:A", "strand": "+"},
        {"ontology_curie": "CL:C", "strand": "+"},
    ]
    canonical_order = [("CL:A", "+"), ("CL:B", "+"), ("CL:C", "+")]

    reordered = reorder_to_canonical(values, metadata_records, canonical_order)

    np.testing.assert_array_equal(reordered[:, 0], values[:, 1])  # A
    np.testing.assert_array_equal(reordered[:, 1], values[:, 0])  # B
    np.testing.assert_array_equal(reordered[:, 2], values[:, 2])  # C


def test_reorder_to_canonical_missing_track_raises():
    values = np.zeros((2, 1))
    metadata_records = [{"ontology_curie": "CL:A", "strand": "+"}]
    with pytest.raises(KeyError):
        reorder_to_canonical(values, metadata_records, [("CL:MISSING", "+")])


class _FakeMetadataFrame:
    def __init__(self, records):
        self._records = records

    def __getitem__(self, columns):
        return self

    def to_dict(self, orient):
        assert orient == "records"
        return self._records


class _FakeTrackOutput:
    def __init__(self, values, records):
        self.values = values
        self.metadata = _FakeMetadataFrame(records)


class _FakeSequenceOutput:
    def __init__(self, values, records):
        self.rna_seq = _FakeTrackOutput(values, records)


class _FakeClient:
    def __init__(self):
        self.calls = []

    def predict_sequence(self, sequence, *, organism, requested_outputs, ontology_terms, interval):
        self.calls.append(sequence)
        values = np.array([[1.0, 2.0]])
        records = [{"ontology_curie": "CL:0001", "strand": "+"}, {"ontology_curie": "CL:0002", "strand": "+"}]
        return _FakeSequenceOutput(values, records)


@pytest.fixture
def fake_alphagenome_module(monkeypatch):
    fake_dna_client = types.SimpleNamespace(
        create=lambda api_key: _FakeClient(),
        Organism=types.SimpleNamespace(HOMO_SAPIENS="homo_sapiens"),
        OutputType=types.SimpleNamespace(RNA_SEQ="RNA_SEQ", CAGE="CAGE"),
    )
    fake_models_pkg = types.ModuleType("alphagenome.models")
    fake_models_pkg.dna_client = fake_dna_client
    fake_pkg = types.ModuleType("alphagenome")
    fake_pkg.models = fake_models_pkg

    monkeypatch.setitem(sys.modules, "alphagenome", fake_pkg)
    monkeypatch.setitem(sys.modules, "alphagenome.models", fake_models_pkg)
    monkeypatch.setitem(sys.modules, "alphagenome.models.dna_client", fake_dna_client)
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "fake-key-for-tests")
    return fake_dna_client


def test_predict_sequence_caches_to_disk_and_dedupes_calls(tmp_path, fake_alphagenome_module):
    from genomics.predictors.genotype_based.analysis.live_alphagenome_prediction import LiveAlphaGenomePredictor

    predictor = LiveAlphaGenomePredictor(cache_dir=tmp_path)

    values1, meta1 = predictor.predict_sequence(
        cache_key="sample1_MC1R_H1", sequence="ACGT", interval=object(),
        output_type="RNA_SEQ", ontology_terms=["CL:0001", "CL:0002"],
    )
    assert (tmp_path / "sample1_MC1R_H1.npz").exists()
    assert (tmp_path / "sample1_MC1R_H1_meta.json").exists()
    assert values1.shape == (1, 2)
    assert meta1 == [{"ontology_curie": "CL:0001", "strand": "+"}, {"ontology_curie": "CL:0002", "strand": "+"}]

    # Second call with the same cache_key must hit the on-disk cache, not the fake client again.
    values2, meta2 = predictor.predict_sequence(
        cache_key="sample1_MC1R_H1", sequence="ACGT", interval=object(),
        output_type="RNA_SEQ", ontology_terms=["CL:0001", "CL:0002"],
    )
    np.testing.assert_array_equal(values1, values2)
    assert meta1 == meta2
    assert len(predictor.client.calls) == 1


def test_predict_sequence_rejects_unknown_output_type(tmp_path, fake_alphagenome_module):
    from genomics.predictors.genotype_based.analysis.live_alphagenome_prediction import LiveAlphaGenomePredictor

    predictor = LiveAlphaGenomePredictor(cache_dir=tmp_path)
    with pytest.raises(ValueError):
        predictor.predict_sequence(
            cache_key="k", sequence="ACGT", interval=object(), output_type="NOT_A_TYPE", ontology_terms=None,
        )
