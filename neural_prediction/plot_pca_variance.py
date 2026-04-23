#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot explained variance from a neural_prediction PCA cache")
    parser.add_argument("--pca-cache-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    bundle = joblib.load(Path(args.pca_cache_dir) / "scaler_pca.joblib")
    pca = bundle["pca"]
    variance = np.asarray(pca.explained_variance_ratio_)
    cumulative = np.cumsum(variance)

    plt.figure(figsize=(8, 5))
    plt.plot(np.arange(1, len(cumulative) + 1), cumulative)
    plt.xlabel("Components")
    plt.ylabel("Cumulative explained variance")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
