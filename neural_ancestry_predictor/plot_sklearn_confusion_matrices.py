#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot confusion matrices from sklearn baseline result JSON files.

Input expected per experiment directory:
  - train_results.json
  - val_results.json
  - test_results.json

Each file should contain:
  - "confusion_matrix": 2D list
  - optional "classification_report": string (used to infer class labels)
"""

from __future__ import annotations

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import json
import re
from pathlib import Path
from typing import List, Optional

import numpy as np


def _use_agg_backend() -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)


def _infer_labels_from_report(report: str, n_classes: int) -> List[str]:
    labels: List[str] = []
    for raw in report.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(("accuracy", "macro avg", "weighted avg")):
            continue
        # Expected rows like:
        # AFR       0.91      0.98      0.95        54
        if re.search(r"\d+\s*$", line):
            tokens = line.split()
            if len(tokens) >= 5:
                label = tokens[0]
                # Skip header row: "precision recall f1-score support"
                if label.lower() != "precision":
                    labels.append(label)
    if len(labels) == n_classes:
        return labels
    return [str(i) for i in range(n_classes)]


def _plot_matrix(cm: np.ndarray, labels: List[str], title: str, output_png: Path) -> None:
    _use_agg_backend()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    im = ax.imshow(cm, cmap="Blues")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title(title)

    max_v = float(np.max(cm)) if cm.size > 0 else 0.0
    threshold = max_v / 2.0 if max_v > 0 else 0.0
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            v = int(cm[i, j])
            color = "white" if v > threshold else "black"
            ax.text(j, i, str(v), ha="center", va="center", color=color, fontsize=8)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=170)
    plt.close(fig)


def _plot_split(exp_dir: Path, split: str, out_dir: Path) -> Optional[Path]:
    p = exp_dir / f"{split}_results.json"
    if not p.exists():
        return None
    with open(p, "r") as f:
        data = json.load(f)

    cm = np.asarray(data.get("confusion_matrix", []), dtype=np.int64)
    if cm.ndim != 2 or cm.shape[0] != cm.shape[1] or cm.shape[0] == 0:
        raise ValueError(f"Confusion matrix inválida em {p}")

    labels = _infer_labels_from_report(data.get("classification_report", ""), cm.shape[0])
    out = out_dir / f"{split}_confusion_matrix.png"
    _plot_matrix(cm, labels, f"{exp_dir.name} - {split}", out)
    return out


def run(experiment_dirs: List[Path], out_subdir: str) -> None:
    for exp_dir in experiment_dirs:
        if not exp_dir.exists():
            print(f"[skip] Diretório não encontrado: {exp_dir}")
            continue
        out_dir = exp_dir / out_subdir
        print(f"\n[exp] {exp_dir}")
        any_done = False
        for split in ("train", "val", "test"):
            out = _plot_split(exp_dir, split, out_dir)
            if out is None:
                print(f"  - {split}: sem arquivo de resultados")
                continue
            print(f"  - {split}: {out}")
            any_done = True
        if not any_done:
            print("  - Nenhum split encontrado (train/val/test_results.json).")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plota matrizes de confusão a partir de *_results.json dos baselines sklearn."
    )
    parser.add_argument(
        "--experiment-dirs",
        nargs="+",
        required=True,
        help="Um ou mais diretórios de experimento (svm/rf/xgboost).",
    )
    parser.add_argument(
        "--out-subdir",
        type=str,
        default="confusion_matrices",
        help="Subpasta dentro de cada experimento para salvar PNGs.",
    )
    args = parser.parse_args()
    run([Path(p) for p in args.experiment_dirs], args.out_subdir)


if __name__ == "__main__":
    main()
