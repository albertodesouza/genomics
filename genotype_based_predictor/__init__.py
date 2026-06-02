# -*- coding: utf-8 -*-
"""Lightweight package init for genotype_based_predictor."""

from __future__ import annotations

import importlib

__all__ = [
    "alignment",
    "analysis",
    "apps",
    "config",
    "data",
    "utils",
    "experiments",
    "models",
    "tools",
]


def __getattr__(name: str):
    if name in __all__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
