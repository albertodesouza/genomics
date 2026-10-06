from __future__ import annotations

import argparse
from pathlib import Path

from genomics.predictors.genotype_based.config import load_config


def run_search(config_path: Path) -> Path:
    config = load_config(config_path)
    if config.hyperparameter_search.pytorch.enabled:
        from genomics.predictors.genotype_based.experiments.search_pytorch import run_search as run_pytorch_search

        return run_pytorch_search(config_path)

    from genomics.predictors.genotype_based.experiments.search_sklearn import run_search as run_sklearn_search

    return run_sklearn_search(config_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Hyperparameter search for genotype predictors")
    parser.add_argument("config_path", type=Path)
    args = parser.parse_args()
    run_search(args.config_path.resolve())


if __name__ == "__main__":
    main()
