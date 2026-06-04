# -*- coding: utf-8 -*-
"""Model exports for the genotype-based predictor."""

from .nn_model import NNAncestryPredictor
from .cnn_model import CNNAncestryPredictor
from .cnn2_model import CNN2AncestryPredictor

SKLEARN_BASELINE_TYPES = frozenset({"SVM", "RF", "XGBOOST"})

__all__ = [
    "NNAncestryPredictor",
    "CNNAncestryPredictor",
    "CNN2AncestryPredictor",
    "SKLEARN_BASELINE_TYPES",
]
