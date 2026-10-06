from genomics.core.metrics import classification_metrics, public_classification_results, with_weighted_metric_aliases


def test_classification_metrics_include_explicit_weighted_aliases():
    results = classification_metrics([0, 0, 1, 1], [0, 1, 1, 1], ["A", "B"])

    assert results["weighted_accuracy"] == results["accuracy"]
    assert results["weighted_precision"] == results["precision"]
    assert results["weighted_recall"] == results["recall"]
    assert results["weighted_f1_score"] == results["f1"]


def test_with_weighted_metric_aliases_backfills_existing_payload():
    payload = {"accuracy": 0.5, "precision": 0.6, "recall": 0.7, "f1": 0.8}

    with_weighted_metric_aliases(payload)

    assert payload["weighted_accuracy"] == 0.5
    assert payload["weighted_precision"] == 0.6
    assert payload["weighted_recall"] == 0.7
    assert payload["weighted_f1_score"] == 0.8


def test_with_weighted_metric_aliases_recomputes_legacy_payload():
    payload = {
        "confusion_matrix": [[8, 2], [1, 9]],
        "per_class_metrics": {
            "A": {"precision": 8 / 9, "recall": 0.8, "f1": 16 / 19, "support": 10},
            "B": {"precision": 9 / 11, "recall": 0.9, "f1": 18 / 20, "support": 10},
        },
    }

    with_weighted_metric_aliases(payload)

    assert payload["weighted_accuracy"] == 0.85
    assert payload["weighted_precision"] > 0.85
    assert payload["weighted_recall"] == 0.85
    assert payload["weighted_f1_score"] > 0.87
    assert payload["accuracy"] == payload["weighted_accuracy"]


def test_public_classification_results_uses_only_weighted_aggregate_names():
    payload = {
        "accuracy": 0.5,
        "precision": 0.6,
        "recall": 0.7,
        "f1": 0.8,
        "confusion_matrix": [[1, 1], [0, 2]],
        "per_class_metrics": {},
    }

    public = public_classification_results(payload)

    assert "accuracy" not in public
    assert "precision" not in public
    assert "recall" not in public
    assert "f1" not in public
    assert public["weighted_accuracy"] == 0.5
    assert public["weighted_precision"] == 0.6
    assert public["weighted_recall"] == 0.7
    assert public["weighted_f1_score"] == 0.8
