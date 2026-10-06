from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import torch
import torch.nn as nn
import torch.optim as optim


def resolve_checkpoint_path(models_dir: Path, checkpoint: str) -> Path:
    path = Path(checkpoint)
    if path.is_absolute() or path.exists():
        return path
    if path.suffix != ".pt":
        path = path.with_suffix(".pt")
    return Path(models_dir) / path


def save_checkpoint(
    path: Path,
    model: nn.Module,
    optimizer: Optional[optim.Optimizer],
    epoch: int,
    metrics: Dict[str, float],
    scheduler: Any = None,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    payload: Dict[str, Any] = {
        "epoch": epoch,
        "model_state_dict": model.state_dict(),
        **metrics,
    }
    if optimizer is not None:
        payload["optimizer_state_dict"] = optimizer.state_dict()
    if scheduler is not None:
        payload["scheduler_state_dict"] = scheduler.state_dict()
    if extra:
        payload.update(extra)
    torch.save(payload, path)


def load_checkpoint(
    path: Path,
    model: nn.Module,
    optimizer: Optional[optim.Optimizer] = None,
    scheduler: Any = None,
    device: torch.device | str = "cpu",
    weights_only: bool = False,
) -> Dict[str, Any]:
    state = torch.load(path, map_location=device, weights_only=weights_only)
    model.load_state_dict(state.get("model_state_dict", state))
    if optimizer is not None and state.get("optimizer_state_dict") is not None:
        optimizer.load_state_dict(state["optimizer_state_dict"])
    if scheduler is not None and state.get("scheduler_state_dict") is not None:
        scheduler.load_state_dict(state["scheduler_state_dict"])
    return state
