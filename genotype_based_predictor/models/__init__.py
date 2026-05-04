# -*- coding: utf-8 -*-
"""
models/__init__.py - Exporta todos os modelos disponíveis.
"""
from genotype_based_predictor.models.nn_model import NNAncestryPredictor
from genotype_based_predictor.models.cnn_model import CNNAncestryPredictor
from genotype_based_predictor.models.cnn2_model import CNN2AncestryPredictor

SKLEARN_BASELINE_TYPES = frozenset({"SVM", "RF", "XGBOOST"})

__all__ = [
    "NNAncestryPredictor",
    "CNNAncestryPredictor",
    "CNN2AncestryPredictor",
    "SKLEARN_BASELINE_TYPES",
]
