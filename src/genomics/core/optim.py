from __future__ import annotations

from typing import Any, Optional

import torch.nn as nn
import torch.optim as optim


def make_optimizer(
    model: nn.Module,
    optimizer: str,
    learning_rate: float,
    weight_decay: float,
    momentum: float = 0.9,
    foreach: Optional[bool] = None,
) -> optim.Optimizer:
    kwargs: dict[str, Any] = {"lr": learning_rate, "weight_decay": weight_decay}
    if foreach is not None and optimizer.lower() in {"adam", "adamw", "sgd"}:
        kwargs["foreach"] = foreach
    opt = optimizer.lower()
    if opt == "adam":
        return optim.Adam(model.parameters(), **kwargs)
    if opt == "adamw":
        return optim.AdamW(model.parameters(), **kwargs)
    if opt == "sgd":
        return optim.SGD(model.parameters(), momentum=momentum, **kwargs)
    raise ValueError(f"Otimizador não suportado: {optimizer}")


def make_optimizer_from_config(model: nn.Module, training_config: Any, foreach: Optional[bool] = None) -> optim.Optimizer:
    return make_optimizer(
        model=model,
        optimizer=training_config.optimizer,
        learning_rate=training_config.learning_rate,
        weight_decay=training_config.weight_decay,
        foreach=foreach,
    )
