import csv
from pathlib import Path

import numpy as np

from genomics.predictors.snp_ancestry.train_ml import build_feature_matrix, impute_with_train_mean, load_marker_ids
from genomics.predictors.snp_ancestry.ablate import ablated_marker_ids


def test_load_marker_ids_reads_top_rows(tmp_path):
    markers = tmp_path / "markers.tsv"
    with open(markers, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "rsid"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"rank": 1, "rsid": "rs1"})
        writer.writerow({"rank": 2, "rsid": "rs2"})

    assert load_marker_ids(markers, top=1) == ["rs1"]


def test_impute_with_train_mean_uses_train_only():
    X_train = np.array([[0.0, np.nan], [2.0, 1.0]], dtype=np.float32)
    X_test = np.array([[np.nan, np.nan]], dtype=np.float32)

    train, [test], means = impute_with_train_mean(X_train, X_test)

    assert means.tolist() == [1.0, 1.0]
    assert train.tolist() == [[0.0, 1.0], [2.0, 1.0]]
    assert test.tolist() == [[1.0, 1.0]]


def test_build_feature_matrix_uses_pipeline_dose_encoding(tmp_path):
    individuals = tmp_path / "individuals"
    sample_dir = individuals / "S1"
    sample_dir.mkdir(parents=True)
    (sample_dir / "S1_23andme.txt").write_text(
        "# header\nrs1\tchr15\t100\tAG\nrs2\tchr15\t200\tGG\n",
        encoding="utf-8",
    )
    config = {
        "input": {"individuals_dir": str(individuals)},
        "conversion": {"output_filename": "{sample_id}_23andme.txt"},
        "prediction": {"haplotype_mode": "H1+H2"},
    }
    stats = {"ref_alleles": {"rs1": "A", "rs2": "G"}}
    meta = {"S1": {"superpopulation": "EUR"}}

    X, y, kept = build_feature_matrix(["S1"], ["rs1", "rs2"], config, stats, meta, "superpopulation")

    assert kept == ["S1"]
    assert y.tolist() == ["EUR"]
    assert X.tolist() == [[1.0, 2.0]]


def test_ablated_marker_ids_removes_ranked_prefix():
    assert ablated_marker_ids(["rs1", "rs2", "rs3"], 0) == ["rs1", "rs2", "rs3"]
    assert ablated_marker_ids(["rs1", "rs2", "rs3"], 2) == ["rs3"]
