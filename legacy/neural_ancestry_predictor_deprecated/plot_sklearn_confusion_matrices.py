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
from typing import Dict, List, Optional, Tuple

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


def _format_confusion_matrix_md(cm: np.ndarray, labels: List[str]) -> str:
    """
    Return a Markdown table representing the confusion matrix with row=true and column=predicted.
    """
    cm_int = np.asarray(cm, dtype=np.int64)
    header = ["True \\ Predicted"] + list(labels)
    md_lines: List[str] = []
    md_lines.append("| " + " | ".join(header) + " |")
    md_lines.append("|" + "|".join(["---"] * len(header)) + "|")
    for i in range(cm_int.shape[0]):
        row = [labels[i]] + [str(int(x)) for x in cm_int[i, :]]
        md_lines.append("| " + " | ".join(row) + " |")
    return "\n".join(md_lines)

def _safe_div(n: float, d: float) -> float:
    return float(n / d) if d != 0 else 0.0


def _metrics_from_confusion_matrix(cm: np.ndarray, labels: List[str]) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    """
    Compute per-class + macro/weighted + accuracy rows from confusion matrix.
    """
    rows: List[Dict[str, object]] = []
    n = cm.shape[0]
    total_support = int(np.sum(cm))

    precisions: List[float] = []
    recalls: List[float] = []
    f1s: List[float] = []
    supports: List[int] = []

    for i in range(n):
        tp = float(cm[i, i])
        fp = float(np.sum(cm[:, i]) - tp)
        fn = float(np.sum(cm[i, :]) - tp)
        support = int(np.sum(cm[i, :]))

        p = _safe_div(tp, tp + fp)
        r = _safe_div(tp, tp + fn)
        f1 = _safe_div(2.0 * p * r, p + r) if (p + r) > 0 else 0.0

        rows.append(
            {
                "Class": labels[i],
                "Precision": p,
                "Recall": r,
                "F1-score": f1,
                "Support": support,
            }
        )
        precisions.append(p)
        recalls.append(r)
        f1s.append(f1)
        supports.append(support)

    macro = {
        "Class": "Macro average",
        "Precision": float(np.mean(precisions)) if precisions else 0.0,
        "Recall": float(np.mean(recalls)) if recalls else 0.0,
        "F1-score": float(np.mean(f1s)) if f1s else 0.0,
        "Support": total_support,
    }
    weights = np.asarray(supports, dtype=np.float64)
    wsum = float(np.sum(weights))
    weighted = {
        "Class": "Weighted average",
        "Precision": float(np.sum(np.asarray(precisions) * weights) / wsum) if wsum > 0 else 0.0,
        "Recall": float(np.sum(np.asarray(recalls) * weights) / wsum) if wsum > 0 else 0.0,
        "F1-score": float(np.sum(np.asarray(f1s) * weights) / wsum) if wsum > 0 else 0.0,
        "Support": total_support,
    }
    accuracy = {
        "Class": "Accuracy",
        "Precision": float(np.trace(cm) / total_support) if total_support > 0 else 0.0,
        "Recall": "",
        "F1-score": "",
        "Support": total_support,
    }

    return rows + [macro, weighted, accuracy], accuracy


def _write_table_files(
    table_rows: List[Dict[str, object]],
    title: str,
    out_md: Path,
    out_csv: Path,
) -> str:
    out_md.parent.mkdir(parents=True, exist_ok=True)

    md_lines: List[str] = []
    md_lines.append(f"## {title}\n")
    md_lines.append("| Class | Precision | Recall | F1-score | Support |\n")
    md_lines.append("|---|---:|---:|---:|---:|\n")
    for r in table_rows:
        p = r["Precision"]
        rc = r["Recall"]
        f1 = r["F1-score"]
        support = r["Support"]
        p_s = f"{p:.2f}" if isinstance(p, (float, int)) else str(p)
        r_s = f"{rc:.2f}" if isinstance(rc, (float, int)) else str(rc)
        f1_s = f"{f1:.2f}" if isinstance(f1, (float, int)) else str(f1)
        md_lines.append(f"| {r['Class']} | {p_s} | {r_s} | {f1_s} | {support} |\n")

    # Markdown
    md_content = "".join(md_lines)
    with open(out_md, "w") as f:
        f.write(md_content)

    # CSV
    with open(out_csv, "w") as f:
        f.write("Class,Precision,Recall,F1-score,Support\n")
        for r in table_rows:
            p = r["Precision"]
            rc = r["Recall"]
            f1 = r["F1-score"]
            support = r["Support"]
            p_s = f"{p:.6f}" if isinstance(p, (float, int)) else str(p)
            r_s = f"{rc:.6f}" if isinstance(rc, (float, int)) else str(rc)
            f1_s = f"{f1:.6f}" if isinstance(f1, (float, int)) else str(f1)
            f.write(f"{r['Class']},{p_s},{r_s},{f1_s},{support}\n")

    return md_content


def _plot_split(exp_dir: Path, split: str, out_dir: Path) -> Optional[Tuple[Path, str]]:
    p = exp_dir / f"{split}_results.json"
    if not p.exists():
        return None
    with open(p, "r") as f:
        data = json.load(f)

    cm = np.asarray(data.get("confusion_matrix", []), dtype=np.int64)
    if cm.ndim != 2 or cm.shape[0] != cm.shape[1] or cm.shape[0] == 0:
        raise ValueError(f"Confusion matrix inválida em {p}")

    labels = _infer_labels_from_report(data.get("classification_report", ""), cm.shape[0])
    md_matrix = _format_confusion_matrix_md(cm, labels)
    out = out_dir / f"{split}_confusion_matrix.png"
    _plot_matrix(cm, labels, f"{exp_dir.name} - {split}", out)
    return out, md_matrix


def _table_for_split(exp_dir: Path, split: str, out_dir: Path) -> Optional[Tuple[Path, Path]]:
    p = exp_dir / f"{split}_results.json"
    if not p.exists():
        return None
    with open(p, "r") as f:
        data = json.load(f)
    cm = np.asarray(data.get("confusion_matrix", []), dtype=np.int64)
    if cm.ndim != 2 or cm.shape[0] != cm.shape[1] or cm.shape[0] == 0:
        raise ValueError(f"Confusion matrix inválida em {p}")
    labels = _infer_labels_from_report(data.get("classification_report", ""), cm.shape[0])
    table_rows, _acc = _metrics_from_confusion_matrix(cm, labels)
    out_md = out_dir / f"{split}_classification_table.md"
    out_csv = out_dir / f"{split}_classification_table.csv"
    title = f"Classification performance on {split} set ({exp_dir.name})"
    md_content = _write_table_files(table_rows, title, out_md, out_csv)
    # Print the table to terminal (Markdown) so you can copy/paste quickly.
    print(f"\n{title}\n")
    print(md_content)
    return out_md, out_csv


def run(
    experiment_dirs: List[Path],
    out_subdir: str,
    table_subdir: str,
    table_splits: List[str],
    matrix_splits: List[str],
) -> None:
    for exp_dir in experiment_dirs:
        if not exp_dir.exists():
            print(f"[skip] Diretório não encontrado: {exp_dir}")
            continue
        out_dir = exp_dir / out_subdir
        table_dir = exp_dir / table_subdir
        print(f"\n[exp] {exp_dir}")
        any_done = False
        for split in ("train", "val", "test"):
            out_pair = _plot_split(exp_dir, split, out_dir)
            if out_pair is None:
                print(f"  - {split}: sem arquivo de resultados")
                continue
            out_png, md_matrix = out_pair
            print(f"  - {split}: {out_png}")
            any_done = True
            if split in matrix_splits:
                print(f"\nConfusion matrix on {split} set ({exp_dir.name})\n")
                print(md_matrix)
        if not any_done:
            print("  - Nenhum split encontrado (train/val/test_results.json).")
        for split in table_splits:
            out_pair = _table_for_split(exp_dir, split, table_dir)
            if out_pair is None:
                print(f"  - tabela {split}: sem arquivo de resultados")
                continue
            out_md, out_csv = out_pair
            print(f"  - tabela {split}: {out_md}")
            print(f"  - tabela {split} (csv): {out_csv}")


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
    parser.add_argument(
        "--table-subdir",
        type=str,
        default="classification_tables",
        help="Subpasta dentro de cada experimento para salvar tabelas (md/csv).",
    )
    parser.add_argument(
        "--table-splits",
        type=str,
        default="val,test",
        help="Splits para gerar tabela, separados por vírgula (ex: val,test).",
    )
    parser.add_argument(
        "--matrix-splits",
        type=str,
        default="val,test",
        help="Splits para imprimir confusion matrices no terminal (ex: val,test).",
    )
    args = parser.parse_args()
    table_splits = [x.strip() for x in args.table_splits.split(",") if x.strip()]
    matrix_splits = [x.strip() for x in args.matrix_splits.split(",") if x.strip()]
    run(
        [Path(p) for p in args.experiment_dirs],
        args.out_subdir,
        args.table_subdir,
        table_splits,
        matrix_splits,
    )


if __name__ == "__main__":
    main()
