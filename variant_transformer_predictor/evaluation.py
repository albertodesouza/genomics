from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import torch
from rich.console import Console

from genomics_pipeline.metrics import classification_metrics, print_classification_metrics, save_results_json
from genomics_pipeline.torch_utils import move_to_device

console = Console()


def move_batch_to_device(batch: Dict, device: torch.device) -> Dict:
    return move_to_device(batch, device)


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
    results = classification_metrics(targets, preds, class_names)
    print_classification_metrics(results, f"{description} metrics", console)
    return results


def save_results(results: Dict, output_path: Path) -> None:
    save_results_json(results, output_path, console)
