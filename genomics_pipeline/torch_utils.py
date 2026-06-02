from __future__ import annotations

from typing import Any

import torch


def move_to_device(value: Any, device: torch.device, non_blocking: bool = True) -> Any:
    if isinstance(value, torch.Tensor):
        return value.to(device, non_blocking=non_blocking)
    if isinstance(value, dict):
        return {key: move_to_device(item, device, non_blocking=non_blocking) for key, item in value.items()}
    if isinstance(value, list):
        return [move_to_device(item, device, non_blocking=non_blocking) for item in value]
    if isinstance(value, tuple):
        return tuple(move_to_device(item, device, non_blocking=non_blocking) for item in value)
    return value
