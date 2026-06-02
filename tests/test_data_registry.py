from pathlib import Path

import pytest

from genomics_pipeline.data_registry import (
    classify_dataset_path,
    registered_dataset_ids,
    resolve_dataset,
    resolve_derived_dataset,
)


def test_resolve_registered_datasets_uses_data_root_env(monkeypatch, tmp_path):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path))

    canonical = resolve_dataset("1kg_high_coverage")
    assert canonical.dataset_id == "1kg_high_coverage"
    assert canonical.role == "canonical"
    assert canonical.path == tmp_path / "v1" / "1kG_high_coverage"
    assert canonical.metadata_path == canonical.path / "dataset_metadata.json"

    derived = resolve_derived_dataset("superpopulation")
    assert derived.dataset_id == "variant_transformer_superpopulation"
    assert derived.role == "derived"
    assert derived.source_dataset_id == "1kg_high_coverage"


def test_classify_dataset_paths(monkeypatch, tmp_path):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path))

    assert classify_dataset_path(tmp_path / "v1" / "1kG_high_coverage") == "canonical-ok"
    assert classify_dataset_path(tmp_path / "v1" / "1kG_high_coverage" / "individuals" / "HG00096") == "canonical-ok"
    assert classify_dataset_path(tmp_path / "top3" / "non_longevous_results_genes_1000_all") == "legacy-top3"
    assert classify_dataset_path(tmp_path / "variant_transformer" / "superpopulation") == "derived-ok"
    assert classify_dataset_path(tmp_path / "unregistered" / "dataset") == "dados-unregistered"
    assert classify_dataset_path(Path("/tmp/outside_genomics_data")) == "external-data-path"
    assert classify_dataset_path(None) == "missing-data-path"


def test_unknown_dataset_id_reports_known_ids(monkeypatch, tmp_path):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path))

    with pytest.raises(KeyError) as exc_info:
        resolve_dataset("missing")

    message = str(exc_info.value)
    assert "Dataset id desconhecido" in message
    for dataset_id in registered_dataset_ids():
        assert dataset_id in message
