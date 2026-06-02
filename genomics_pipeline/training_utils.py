from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, MutableMapping, Optional

import torch.optim as optim


def _get_config_value(config: Any, name: str, default: Any = None) -> Any:
    if isinstance(config, dict):
        return config.get(name, default)
    return getattr(config, name, default)


def make_lr_scheduler(optimizer: optim.Optimizer, scheduler_config: Any, *, default_t_max: Optional[int] = None):
    if not _get_config_value(scheduler_config, "enabled", False):
        return None
    scheduler_type = str(_get_config_value(scheduler_config, "type", "plateau")).lower()
    if scheduler_type == "plateau":
        return optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode=_get_config_value(scheduler_config, "mode", "min"),
            factor=_get_config_value(scheduler_config, "factor", 0.5),
            patience=_get_config_value(scheduler_config, "patience", 10),
            min_lr=_get_config_value(scheduler_config, "min_lr", 1e-6),
        )
    if scheduler_type == "cosine":
        return optim.lr_scheduler.CosineAnnealingLR(
            optimizer,
            T_max=_get_config_value(scheduler_config, "T_max", default_t_max or 100),
            eta_min=_get_config_value(scheduler_config, "eta_min", 1e-6),
        )
    if scheduler_type == "step":
        return optim.lr_scheduler.StepLR(
            optimizer,
            step_size=_get_config_value(scheduler_config, "step_size", 30),
            gamma=_get_config_value(scheduler_config, "gamma", 0.1),
        )
    if scheduler_type == "exponential":
        return optim.lr_scheduler.ExponentialLR(
            optimizer,
            gamma=_get_config_value(scheduler_config, "gamma", 0.95),
        )
    if scheduler_type == "multistep":
        return optim.lr_scheduler.MultiStepLR(
            optimizer,
            milestones=_get_config_value(scheduler_config, "milestones", [30, 60, 90]),
            gamma=_get_config_value(scheduler_config, "gamma", 0.1),
        )
    if scheduler_type == "cosine_warm_restarts":
        return optim.lr_scheduler.CosineAnnealingWarmRestarts(
            optimizer,
            T_0=_get_config_value(scheduler_config, "T_0", 50),
            T_mult=_get_config_value(scheduler_config, "T_mult", 1),
            eta_min=_get_config_value(scheduler_config, "eta_min", 1e-6),
        )
    raise ValueError(f"Scheduler nao suportado: {scheduler_type}")


def step_lr_scheduler(scheduler: Any, val_loss: float, *, metric_available: bool = True) -> None:
    if scheduler is None:
        return
    if isinstance(scheduler, optim.lr_scheduler.ReduceLROnPlateau):
        if metric_available:
            scheduler.step(val_loss)
    else:
        scheduler.step()


def new_training_history() -> Dict[str, list]:
    return {"epoch": [], "train_loss": [], "train_accuracy": [], "val_loss": [], "val_accuracy": []}


def append_history_epoch(
    history: MutableMapping[str, list],
    epoch: int,
    train_metrics: Dict[str, float],
    val_metrics: Dict[str, float],
) -> None:
    history.setdefault("epoch", []).append(epoch)
    history.setdefault("train_loss", []).append(train_metrics["loss"])
    history.setdefault("train_accuracy", []).append(train_metrics["accuracy"])
    history.setdefault("val_loss", []).append(val_metrics["loss"])
    history.setdefault("val_accuracy", []).append(val_metrics["accuracy"])


def write_training_history(path: Path, history: Dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(history, f, indent=2)
