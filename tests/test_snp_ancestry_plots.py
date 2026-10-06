import csv
import json

from genomics.predictors.snp_ancestry.plots import discover_model_dirs, metric_labels


def test_metric_labels_reads_classification_report_labels():
    metrics = {
        "confusion_matrix": [[2, 0], [1, 3]],
        "classification_report": """
              precision    recall  f1-score   support

         AFR       1.00      0.67      0.80         3
         EUR       0.75      1.00      0.86         3

    accuracy                           0.83         6
   macro avg       0.88      0.83      0.83         6
weighted avg       0.88      0.83      0.83         6
""",
    }

    assert metric_labels(metrics) == ["AFR", "EUR"]


def test_discover_model_dirs_requires_metrics_json(tmp_path):
    logistic = tmp_path / "logistic"
    logistic.mkdir()
    (logistic / "metrics.json").write_text(json.dumps({"metrics": {}}), encoding="utf-8")
    ignored = tmp_path / "random_forest"
    ignored.mkdir()
    with open(ignored / "feature_importance.tsv", "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "rsid", "importance"], delimiter="\t")
        writer.writeheader()

    assert discover_model_dirs(tmp_path) == [logistic]
