"""Generate UMAP/PCA visualizations from processed tensor shards.

The script avoids flattening full tensors directly. Each sample is summarized by
per-channel statistics, then embedded in 2D and colored by target class.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def _load_class_names(cache_dir: Path) -> Dict[int, str]:
    meta_path = cache_dir / "metadata.json"
    if not meta_path.exists():
        return {}
    with open(meta_path) as f:
        meta = json.load(f)
    raw = meta.get("class_names", {}) or {}
    return {int(k): str(v) for k, v in raw.items()}


def _load_split_sample_ids(cache_dir: Path, split: str) -> List[str]:
    path = cache_dir / "split_index.json"
    if not path.exists():
        return []
    with open(path) as f:
        payload = json.load(f)
    return [str(x) for x in payload.get(split, [])]


def _iter_split_items(cache_dir: Path, split: str, max_samples: Optional[int]) -> Iterable[Tuple[int, torch.Tensor, torch.Tensor]]:
    loaded = 0
    shard_paths = sorted(cache_dir.glob(f"{split}_data_shard_*.pt"))
    if shard_paths:
        for shard_path in shard_paths:
            shard = torch.load(shard_path, map_location="cpu")
            for features, target in shard:
                yield loaded, features, target
                loaded += 1
                if max_samples is not None and loaded >= max_samples:
                    return
        return

    data_path = cache_dir / f"{split}_data.pt"
    if not data_path.exists():
        raise FileNotFoundError(f"No {split}_data_shard_*.pt or {split}_data.pt found in {cache_dir}")
    data = torch.load(data_path, map_location="cpu")
    for features, target in data:
        yield loaded, features, target
        loaded += 1
        if max_samples is not None and loaded >= max_samples:
            return


def _summarize_tensor(features: torch.Tensor) -> np.ndarray:
    x = features.detach().cpu().float().numpy()
    flat = x.reshape(-1, x.shape[-1])
    means = flat.mean(axis=1)
    stds = flat.std(axis=1)
    mins = flat.min(axis=1)
    maxs = flat.max(axis=1)
    nonzero = (flat != 0).mean(axis=1)
    return np.concatenate([means, stds, mins, maxs, nonzero]).astype(np.float32, copy=False)


def _embed(features: np.ndarray, method: str, seed: int) -> Tuple[np.ndarray, str]:
    scaled = StandardScaler().fit_transform(features)
    if method == "umap":
        try:
            import umap  # type: ignore
        except ImportError:
            print("umap-learn not installed; falling back to PCA")
        else:
            reducer = umap.UMAP(n_components=2, random_state=seed, n_neighbors=15, min_dist=0.1)
            return reducer.fit_transform(scaled), "umap"
    return PCA(n_components=2, random_state=seed).fit_transform(scaled), "pca"


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot UMAP/PCA of processed cache shards")
    parser.add_argument("cache_dir", help="Processed cache directory containing shards")
    parser.add_argument("--split", default="train", choices=["train", "val", "test"], help="Split to visualize")
    parser.add_argument("--max-samples", type=int, default=1000, help="Maximum samples to load")
    parser.add_argument("--method", choices=["umap", "pca"], default="umap", help="Embedding method")
    parser.add_argument("--seed", type=int, default=13, help="Random seed")
    parser.add_argument("--output", default="umap.png", help="Output plot path")
    parser.add_argument("--csv-output", default=None, help="Optional CSV output with coordinates")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    class_names = _load_class_names(cache_dir)
    sample_ids = _load_split_sample_ids(cache_dir, args.split)

    rows = []
    labels = []
    row_sample_ids = []
    for idx, features, target in _iter_split_items(cache_dir, args.split, args.max_samples):
        rows.append(_summarize_tensor(features))
        label = int(target.item() if target.ndim == 0 else target.reshape(-1)[0].item())
        labels.append(label)
        row_sample_ids.append(sample_ids[idx] if idx < len(sample_ids) else str(idx))

    if not rows:
        raise RuntimeError(f"No samples loaded from {cache_dir} split={args.split}")

    X = np.vstack(rows)
    y = np.asarray(labels, dtype=np.int64)
    coords, used_method = _embed(X, args.method, args.seed)

    plt.figure(figsize=(10, 8))
    for label in sorted(set(labels)):
        mask = y == label
        name = class_names.get(label, str(label))
        plt.scatter(coords[mask, 0], coords[mask, 1], s=22, alpha=0.8, label=name)
    plt.title(f"{used_method.upper()} of processed cache ({args.split}, n={len(labels)})")
    plt.xlabel(f"{used_method.upper()} 1")
    plt.ylabel(f"{used_method.upper()} 2")
    plt.legend(title="class", markerscale=1.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved {args.output}")

    if args.csv_output:
        with open(args.csv_output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["sample_id", "label", "class_name", "x", "y"])
            for sid, label, xy in zip(row_sample_ids, labels, coords):
                writer.writerow([sid, label, class_names.get(label, str(label)), float(xy[0]), float(xy[1])])
        print(f"Saved {args.csv_output}")


if __name__ == "__main__":
    main()
