import csv
import json
from pathlib import Path

import pytest

from genomics.predictors.snp_ancestry.markers import build_marker_table, load_statistics, write_marker_table


def _stats() -> dict:
    return {
        "metadata": {"populations": ["AFR", "EUR", "EAS"]},
        "ref_alleles": {"rs1": "A", "rs2": "G", "rs3": "T"},
        "snp_info": {"rs1": ["chr15", "100"], "rs2": ["chr15", "200"], "rs3": ["chr15", "300"]},
        "allele_frequencies": {
            "rs1": [0.95, 0.05, 0.05],
            "rs2": [0.45, 0.50, 0.55],
            "rs3": [0.80, 0.20, 0.20],
        },
    }


def test_build_marker_table_ranks_by_fst():
    rows = build_marker_table(_stats(), score="fst", top=2)

    assert [row["rsid"] for row in rows] == ["rs1", "rs3"]
    assert rows[0]["rank"] == 1
    assert rows[0]["max_population"] == "AFR"
    assert rows[0]["min_population"] == "EUR"
    assert rows[0]["freq_AFR"] == pytest.approx(0.95)


def test_write_marker_table_outputs_tsv(tmp_path):
    output = tmp_path / "aims.tsv"
    rows = build_marker_table(_stats(), score="max_delta_frequency", top=1)

    write_marker_table(rows, output, _stats()["metadata"]["populations"])

    with open(output, newline="", encoding="utf-8") as f:
        records = list(csv.DictReader(f, delimiter="\t"))

    assert records[0]["rsid"] == "rs1"
    assert records[0]["score_name"] == "max_delta_frequency"
    assert records[0]["freq_AFR"] == "0.95"


def test_load_statistics_accepts_explicit_path(tmp_path):
    config_path = tmp_path / "config.yaml"
    stats_path = tmp_path / "stats.json"
    config_path.write_text("{}\n", encoding="utf-8")
    stats_path.write_text(json.dumps(_stats()), encoding="utf-8")

    stats, path = load_statistics(config_path, stats_path)

    assert path == stats_path
    assert stats["metadata"]["populations"] == ["AFR", "EUR", "EAS"]
