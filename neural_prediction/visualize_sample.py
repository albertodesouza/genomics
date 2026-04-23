#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from .config import load_config
from .data import create_dataloaders


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize one processed sample from neural_prediction")
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-index", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(Path(args.config))
    data_bundle = create_dataloaders(config)
    features, target = data_bundle.dataset[args.sample_index]

    plt.figure(figsize=(10, 6))
    plt.imshow(features.numpy(), aspect="auto", cmap="viridis")
    plt.colorbar()
    plt.title(f"Sample {args.sample_index} target={target}")
    plt.xlabel("Position")
    plt.ylabel("Track row")
    plt.tight_layout()
    plt.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
