#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config
from .data import create_dataloaders, get_dataset_cache_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify neural_prediction processed dataset and cache")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(Path(args.config))
    data_bundle = create_dataloaders(config)
    cache_dir = get_dataset_cache_dir(config)

    print(f"Dataset length: {len(data_bundle.dataset)}")
    print(f"Input shape: {data_bundle.dataset.get_input_shape()}")
    print(f"Num classes/targets: {data_bundle.dataset.get_num_classes()}")
    print(f"Train batches: {len(data_bundle.train_loader)}")
    print(f"Val batches: {len(data_bundle.val_loader)}")
    print(f"Test batches: {len(data_bundle.test_loader)}")
    print(f"Cache dir: {cache_dir if cache_dir is not None else 'disabled'}")
    if cache_dir is not None:
        print(f"Cache exists: {(cache_dir / 'cache.json').exists()}")


if __name__ == "__main__":
    main()
