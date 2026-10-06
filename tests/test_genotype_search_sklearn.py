def test_search_plot_labels_use_candidate_names():
    from genomics.predictors.genotype_based.experiments import search_sklearn

    rows = [
        {"candidate": "005_rf", "model_type": "RF", "val_accuracy": 0.82},
        {"candidate": "002_xgboost", "model_type": "XGBOOST", "val_accuracy": 0.79},
    ]

    assert search_sklearn._search_plot_labels(rows) == ["005_rf\nRF", "002_xgboost\nXGBOOST"]
