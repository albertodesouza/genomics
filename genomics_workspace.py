from __future__ import annotations

import os
from pathlib import Path

from genomics_pipeline.data_registry import legacy_top3_ref, resolve_dataset


def repo_root() -> Path:
    return Path(__file__).resolve().parent


def data_root() -> Path:
    return Path(os.environ.get("GENOMICS_DATA_ROOT", "/dados/GENOMICS_DATA")).resolve()


def results_root() -> Path:
    return Path(os.environ.get("GENOMICS_RESULTS_ROOT", repo_root() / "results")).resolve()


def data_path(*parts: str) -> Path:
    return data_root().joinpath(*parts)


def results_path(*parts: str) -> Path:
    return results_root().joinpath(*parts)


def cache_path(*parts: str) -> Path:
    return results_path("cache", *parts)


DEFAULT_DATASET_ID = "1kg_high_coverage"
DEFAULT_DATASET_DIR = resolve_dataset(DEFAULT_DATASET_ID).path
CANONICAL_1KG_HIGH_COVERAGE_DIR = DEFAULT_DATASET_DIR
LEGACY_TOP3_DATASET_DIR = legacy_top3_ref().path
DEFAULT_CONSENSUS_DATASET_DIR = CANONICAL_1KG_HIGH_COVERAGE_DIR
DEFAULT_GENOTYPE_RUNS_ROOT = results_path("genotype_based_predictor", "runs")
DEFAULT_GENOTYPE_CACHE_ROOT = cache_path("genotype_based_predictor")
DEFAULT_VARIANT_TRANSFORMER_DATASET_DIR = resolve_dataset("variant_transformer_superpopulation").path
DEFAULT_VARIANT_TRANSFORMER_RUNS_ROOT = results_path("variant_transformer_predictor", "runs")
