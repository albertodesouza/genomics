# -*- coding: utf-8 -*-
"""Lightweight package init for genotype_based_predictor."""

from __future__ import annotations

import importlib

__all__ = [
    "config",
    "normalization",
    "utils",
    "dataset",
    "data_splitting",
    "data_pipeline",
    "dataset_layout",
    "genomic_dataset",
    "indel_tensor_builder",
    "dataset_subset_manifest",
    "training",
    "evaluation",
    "interpretability",
    "experiment",
    "models",
]


def __getattr__(name: str):
    if name in __all__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
