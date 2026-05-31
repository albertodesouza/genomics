from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import torch
from rich.console import Console
from rich.table import Table
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_fscore_support

console = Console()


def move_batch_to_device(batch: Dict, device: torch.device) -> Dict:
    moved = {}
    for key, value in batch.items():
        moved[key] = value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value
    return moved


@torch.no_grad()
def evaluate(model, loader, class_names: List[str], device: torch.device, description: str = "test") -> Dict:
    model.eval()
    preds: List[int] = []
    targets: List[int] = []
    for batch in loader:
        batch = move_batch_to_device(batch, device)
        logits = model(batch)
        preds.extend(logits.argmax(dim=1).detach().cpu().tolist())
        targets.extend(batch["targets"].detach().cpu().tolist())
    if not targets:
        return {}
    labels = list(range(len(class_names)))
    p, r, f1, _ = precision_recall_fscore_support(targets, preds, average="weighted", zero_division=0)
    acc = accuracy_score(targets, preds)
    results = {
        "accuracy": float(acc),
        "precision": float(p),
        "recall": float(r),
        "f1": float(f1),
        "confusion_matrix": confusion_matrix(targets, preds, labels=labels).tolist(),
        "classification_report": classification_report(targets, preds, labels=labels, target_names=class_names, zero_division=0),
        "num_samples": len(targets),
    }
    table = Table(title=f"{description} metrics")
    table.add_column("metric")
    table.add_column("value", justify="right")
    table.add_row("accuracy", f"{acc:.4f}")
    table.add_row("precision", f"{float(p):.4f}")
    table.add_row("recall", f"{float(r):.4f}")
    table.add_row("f1", f"{float(f1):.4f}")
    table.add_row("samples", str(len(targets)))
    console.print(table)
    console.print(results["classification_report"])
    return results


def save_results(results: Dict, output_path: Path) -> None:
    serializable = {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in results.items()}
    with open(output_path, "w") as f:
        json.dump(serializable, f, indent=2)
    console.print(f"[green]Resultados salvos:[/green] {output_path}")
