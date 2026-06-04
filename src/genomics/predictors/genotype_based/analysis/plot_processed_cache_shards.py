"""Plot processed tensor shards for quick visual alignment checks.

This is a versioned copy of the small inspection script used during cache
materialization. It loads a few processed `.pt` shards, plots one channel from
both haplotypes across multiple samples, and writes a PNG for visual review.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import torch


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot processed cache shard signals")
    parser.add_argument("cache_dir", nargs="?", default=".", help="Directory containing train_data_shard_*.pt")
    parser.add_argument("--split", default="train", choices=["train", "val", "test"], help="Dataset split to inspect")
    parser.add_argument("--max-shards", type=int, default=5, help="Maximum number of shards to load")
    parser.add_argument("--max-samples", type=int, default=16, help="Maximum number of samples to plot")
    parser.add_argument("--channel", type=int, default=0, help="Channel index to plot within each haplotype")
    parser.add_argument("--output", default="plot.png", help="Output PNG path")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    shard_paths = sorted(cache_dir.glob(f"{args.split}_data_shard_*.pt"))[: args.max_shards]
    if not shard_paths:
        raise FileNotFoundError(f"No {args.split}_data_shard_*.pt files found in {cache_dir}")

    all_data = []
    for path in shard_paths:
        print(f"Loading {path}...")
        shard = torch.load(path, map_location="cpu")
        all_data.extend(shard)

    print(f"Total samples loaded: {len(all_data)}")

    plt.figure(figsize=(18, 8))
    plotted = 0
    for sample_idx in range(min(args.max_samples, len(all_data))):
        features = all_data[sample_idx][0]
        for hap_idx in range(features.shape[0]):
            if args.channel >= features.shape[1]:
                raise ValueError(f"channel={args.channel} out of range for tensor shape {tuple(features.shape)}")
            y = features[hap_idx, args.channel].cpu().numpy()
            plt.plot(y, alpha=0.5)
            plotted += 1

    plt.title(f"Processed cache signals: split={args.split}, channel={args.channel}, lines={plotted}")
    plt.xlabel("Position")
    plt.ylabel("Value")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
