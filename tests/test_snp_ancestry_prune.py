import csv

import pytest

from genomics.predictors.snp_ancestry.prune import load_marker_rows, prune_marker_rows, write_marker_rows


def _rows():
    return [
        {"rank": "1", "rsid": "rs1", "chrom": "chr15", "position": "100", "fst": "0.9"},
        {"rank": "2", "rsid": "rs2", "chrom": "chr15", "position": "120", "fst": "0.8"},
        {"rank": "3", "rsid": "rs3", "chrom": "chr15", "position": "200000", "fst": "0.7"},
        {"rank": "4", "rsid": "rs4", "chrom": "chr16", "position": "110", "fst": "0.6"},
    ]


def test_prune_marker_rows_keeps_highest_rank_per_window():
    kept, removed = prune_marker_rows(_rows(), window_bp=50000)

    assert [row["rsid"] for row in kept] == ["rs1", "rs3", "rs4"]
    assert [row["rank"] for row in kept] == ["1", "2", "3"]
    assert removed[0]["rsid"] == "rs2"
    assert removed[0]["pruned_by"] == "rs1"


def test_prune_marker_rows_rejects_negative_window():
    with pytest.raises(ValueError, match="window_bp"):
        prune_marker_rows(_rows(), window_bp=-1)


def test_load_and_write_marker_rows_roundtrip(tmp_path):
    path = tmp_path / "aims.tsv"
    fieldnames = ["rank", "rsid", "chrom", "position", "fst"]
    write_marker_rows(_rows(), fieldnames, path)

    loaded, loaded_fields = load_marker_rows(path)

    assert loaded_fields == fieldnames
    assert loaded[0]["rsid"] == "rs1"
    with open(path, "r", encoding="utf-8", newline="") as f:
        assert list(csv.DictReader(f, delimiter="\t"))[1]["rsid"] == "rs2"
