from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
from rich.console import Console
from rich.table import Table


def _fallback_classification_metrics(y_true: Sequence[int], y_pred: Sequence[int], class_names: Sequence[str]) -> Dict[str, Any]:
    labels = list(range(len(class_names)))
    y_true_arr = np.asarray(list(y_true), dtype=int)
    y_pred_arr = np.asarray(list(y_pred), dtype=int)
    cm = np.zeros((len(labels), len(labels)), dtype=int)
    for target, pred in zip(y_true_arr, y_pred_arr):
        if 0 <= target < len(labels) and 0 <= pred < len(labels):
            cm[target, pred] += 1

    total = int(cm.sum())
    acc = float(np.trace(cm) / total) if total else 0.0
    rows = ["              precision    recall  f1-score   support"]
    weighted_p = weighted_r = weighted_f1 = 0.0
    for idx, name in enumerate(class_names):
        tp = float(cm[idx, idx])
        fp = float(cm[:, idx].sum() - cm[idx, idx])
        fn = float(cm[idx, :].sum() - cm[idx, idx])
        support = float(cm[idx, :].sum())
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        weighted_p += precision * support
        weighted_r += recall * support
        weighted_f1 += f1 * support
        rows.append(f"{name:>12} {precision:10.2f} {recall:9.2f} {f1:9.2f} {int(support):9d}")
    if total:
        weighted_p /= total
        weighted_r /= total
        weighted_f1 /= total
    rows.append("")
    rows.append(f"{'accuracy':>12} {'':>10} {'':>9} {acc:9.2f} {total:9d}")
    rows.append(f"{'weighted avg':>12} {weighted_p:10.2f} {weighted_r:9.2f} {weighted_f1:9.2f} {total:9d}")
    return {
        "accuracy": acc,
        "precision": float(weighted_p),
        "recall": float(weighted_r),
        "f1": float(weighted_f1),
        "confusion_matrix": cm.tolist(),
        "classification_report": "\n".join(rows),
        "num_samples": total,
    }


def classification_metrics(y_true: Sequence[int], y_pred: Sequence[int], class_names: Sequence[str]) -> Dict[str, Any]:
    try:
        from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_fscore_support

        labels = list(range(len(class_names)))
        p, r, f1, _ = precision_recall_fscore_support(y_true, y_pred, average="weighted", zero_division=0)
        acc = accuracy_score(y_true, y_pred)
        return {
            "accuracy": float(acc),
            "precision": float(p),
            "recall": float(r),
            "f1": float(f1),
            "confusion_matrix": confusion_matrix(y_true, y_pred, labels=labels).tolist(),
            "classification_report": classification_report(y_true, y_pred, labels=labels, target_names=list(class_names), zero_division=0),
            "num_samples": len(y_true),
        }
    except ImportError:
        return _fallback_classification_metrics(y_true, y_pred, class_names)


def print_classification_metrics(results: Dict[str, Any], title: str, console: Console) -> None:
    table = Table(title=title, show_header=True)
    table.add_column("Métrica")
    table.add_column("Valor", justify="right")
    table.add_row("Accuracy", f"{float(results['accuracy']):.4f}")
    table.add_row("Precision (weighted)", f"{float(results['precision']):.4f}")
    table.add_row("Recall (weighted)", f"{float(results['recall']):.4f}")
    table.add_row("F1 (weighted)", f"{float(results['f1']):.4f}")
    table.add_row("Amostras", str(results.get("num_samples", "")))
    console.print(table)
    report = results.get("classification_report")
    if report:
        console.print(report)


def save_results_json(results: Dict[str, Any], output_path: Path, console: Console | None = None) -> None:
    serializable: Dict[str, Any] = {}
    for key, value in results.items():
        if isinstance(value, np.ndarray):
            serializable[key] = value.tolist()
        else:
            serializable[key] = value
    with open(output_path, "w") as f:
        json.dump(serializable, f, indent=2)
    if console is not None:
        console.print(f"[green]Resultados salvos:[/green] {output_path}")
