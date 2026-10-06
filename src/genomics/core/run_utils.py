from __future__ import annotations

from typing import Any, Dict, Tuple

import torch


def select_device(prefer_cuda: bool = True) -> torch.device:
    return torch.device("cuda" if prefer_cuda and torch.cuda.is_available() else "cpu")


def select_split_loader(choice: str, train_loader: Any, val_loader: Any, test_loader: Any) -> Tuple[Any, str]:
    normalized = str(choice or "test").lower()
    if normalized == "train":
        return train_loader, "train"
    if normalized in {"val", "validation"}:
        return val_loader, "val"
    return test_loader, "test"


def training_manifest_fields(history: Dict[str, Any]) -> Dict[str, Any]:
    epochs = history.get("epoch") or []
    val_accuracy = history.get("val_accuracy") or []
    val_loss = history.get("val_loss") or []
    return {
        "training_history": "models/training_history.json",
        "last_epoch": history.get("last_epoch") if history.get("last_epoch") is not None else (epochs[-1] if epochs else None),
        "best_val_accuracy": max(val_accuracy or [0.0]),
        "best_val_loss": min(val_loss or [float("inf")]),
    }
