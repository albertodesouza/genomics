#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot only the cumulative explained variance ratio from a PCA variance JSON file.

The input JSON is produced by plot_sklearn_pca_variance.py --json-output.
By default, the plot highlights the 0.95 cumulative-variance threshold without a
title, adding dashed guide lines and explicit axis tick values at the marked
point.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Optional, Tuple

import numpy as np


os.environ.setdefault("MPLBACKEND", "Agg")


def _first_k_from_json_or_curve(data: dict, cumulative: np.ndarray, threshold: float) -> Optional[int]:
    threshold_key = f"{threshold:.2f}"
    threshold_data = data.get("thresholds", {}).get(threshold_key)
    if isinstance(threshold_data, dict) and threshold_data.get("first_k") is not None:
        return int(threshold_data["first_k"])

    idx = np.searchsorted(cumulative, threshold, side="left")
    if idx >= len(cumulative):
        return None
    return int(idx + 1)


def _load_cumulative_curve(input_json: Path, threshold: float) -> Tuple[np.ndarray, Optional[int]]:
    with open(input_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    if "cumulative_explained_variance_ratio" not in data:
        raise ValueError("JSON does not contain 'cumulative_explained_variance_ratio'.")

    cumulative = np.asarray(data["cumulative_explained_variance_ratio"], dtype=np.float64)
    if cumulative.ndim != 1 or len(cumulative) == 0:
        raise ValueError("'cumulative_explained_variance_ratio' must be a non-empty 1D array.")

    return cumulative, _first_k_from_json_or_curve(data, cumulative, threshold)


def plot_cumulative_from_json(input_json: Path, output_path: Path, threshold: float, show: bool) -> None:
    if not show:
        import matplotlib

        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    cumulative, first_k = _load_cumulative_curve(input_json, threshold)
    k_axis = np.arange(1, len(cumulative) + 1)

    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    ax.plot(k_axis, cumulative, color="black", linewidth=1.8)
    ax.set_xlim(1, int(k_axis[-1]))
    ax.set_ylim(0.0, min(1.02, max(float(cumulative[-1]), threshold) + 0.03))

    if first_k is not None:
        ax.scatter([first_k], [threshold], color="tab:red", s=42, zorder=3)
        ax.plot([first_k, first_k], [0.0, threshold], color="tab:red", linestyle="--", linewidth=1.0)
        ax.plot([1, first_k], [threshold, threshold], color="tab:red", linestyle="--", linewidth=1.0)

        xticks = sorted(set(int(x) for x in ax.get_xticks() if 1 <= x <= k_axis[-1]) | {int(first_k)})
        yticks = sorted(set(float(y) for y in ax.get_yticks() if 0.0 <= y <= 1.0) | {float(threshold)})
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels([str(x) for x in xticks])
        ax.set_yticklabels([f"{y:.2f}" for y in yticks])
    else:
        ax.axhline(threshold, color="tab:red", linestyle="--", linewidth=1.0)
        print(f"Aviso: threshold {threshold:.2f} não foi atingido pela curva.")

    ax.set_xlabel("Number of components (k)")
    ax.set_ylabel("Cumulative explained variance ratio")
    ax.grid(True, alpha=0.25, linewidth=0.6)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    print(f"Figura guardada em {output_path.resolve()}")
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot only cumulative explained variance ratio from PCA variance JSON."
    )
    parser.add_argument(
        "--input-json",
        type=str,
        required=True,
        help="JSON produced by plot_sklearn_pca_variance.py --json-output.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="pca_cumulative_variance_clean.png",
        help="PNG path for the clean cumulative variance figure.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.95,
        help="Cumulative variance threshold to mark (default: 0.95).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive matplotlib window after saving.",
    )
    args = parser.parse_args()

    plot_cumulative_from_json(
        Path(args.input_json),
        Path(args.output),
        args.threshold,
        args.show,
    )


if __name__ == "__main__":
    main()
