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
    per_class = {}
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
        per_class[name] = {
            "precision": float(precision),
            "recall": float(recall),
            "f1": float(f1),
            "support": int(support),
        }
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
        "per_class_metrics": per_class,
        "classification_report": "\n".join(rows),
        "num_samples": total,
    }


def classification_metrics(y_true: Sequence[int], y_pred: Sequence[int], class_names: Sequence[str]) -> Dict[str, Any]:
    try:
        from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_fscore_support

        labels = list(range(len(class_names)))
        p, r, f1, _ = precision_recall_fscore_support(y_true, y_pred, average="weighted", zero_division=0)
        acc = accuracy_score(y_true, y_pred)
        per_p, per_r, per_f1, per_support = precision_recall_fscore_support(y_true, y_pred, labels=labels, average=None, zero_division=0)
        per_class = {
            str(name): {
                "precision": float(per_p[idx]),
                "recall": float(per_r[idx]),
                "f1": float(per_f1[idx]),
                "support": int(per_support[idx]),
            }
            for idx, name in enumerate(class_names)
        }
        return {
            "accuracy": float(acc),
            "precision": float(p),
            "recall": float(r),
            "f1": float(f1),
            "confusion_matrix": confusion_matrix(y_true, y_pred, labels=labels).tolist(),
            "per_class_metrics": per_class,
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
    per_class = results.get("per_class_metrics") or {}
    if per_class:
        pct = Table(title="Per-class metrics", show_header=True)
        pct.add_column("Classe")
        pct.add_column("Precision", justify="right")
        pct.add_column("Recall", justify="right")
        pct.add_column("F1", justify="right")
        pct.add_column("Support", justify="right")
        for class_name, row in per_class.items():
            pct.add_row(
                str(class_name),
                f"{float(row.get('precision', 0.0)):.4f}",
                f"{float(row.get('recall', 0.0)):.4f}",
                f"{float(row.get('f1', 0.0)):.4f}",
                str(row.get("support", 0)),
            )
        console.print(pct)
    confusion = results.get("confusion_matrix") or []
    if confusion:
        labels = list(per_class) if per_class else [str(i) for i in range(len(confusion))]
        cmt = Table(title="Confusion matrix", show_header=True)
        cmt.add_column("true\\pred")
        for label in labels:
            cmt.add_column(str(label), justify="right")
        for idx, row in enumerate(confusion):
            cmt.add_row(str(labels[idx] if idx < len(labels) else idx), *(str(value) for value in row))
        console.print(cmt)


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


def save_classification_plots(results: Dict[str, Any], output_dir: Path, prefix: str, console: Console | None = None) -> None:
    """Save ready-to-use PNG plots for classification metrics when matplotlib is available."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        if console is not None:
            console.print("[yellow]matplotlib indisponível; plots PNG não foram gerados.[/yellow]")
        return

    per_class = results.get("per_class_metrics") or {}
    if per_class:
        labels = list(per_class)
        metrics = ["precision", "recall", "f1"]
        x = np.arange(len(labels))
        width = 0.25
        fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.2), 5))
        for offset, metric in enumerate(metrics):
            values = [float(per_class[label].get(metric, 0.0)) for label in labels]
            ax.bar(x + (offset - 1) * width, values, width, label=metric)
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel("Score")
        ax.set_title(f"{prefix} per-class metrics")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / f"{prefix}_per_class_metrics.png", dpi=160)
        plt.close(fig)

    confusion = results.get("confusion_matrix") or []
    if confusion:
        labels = list(per_class) if per_class else [str(i) for i in range(len(confusion))]
        matrix = np.asarray(confusion, dtype=float)
        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.0), max(5, len(labels) * 0.9)))
        im = ax.imshow(matrix, cmap="Blues")
        ax.set_title(f"{prefix} confusion matrix")
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels)
        threshold = matrix.max() / 2.0 if matrix.size else 0.0
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                ax.text(j, i, str(int(matrix[i, j])), ha="center", va="center", color="white" if matrix[i, j] > threshold else "black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        fig.savefig(output_dir / f"{prefix}_confusion_matrix.png", dpi=160)
        plt.close(fig)
    if console is not None:
        console.print(f"[green]Plots salvos em:[/green] {output_dir}")
