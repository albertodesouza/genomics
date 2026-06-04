from __future__ import annotations

from typing import List

import torch


def pad_1d(values: List[torch.Tensor], pad_value: int = 0) -> torch.Tensor:
    max_len = max((value.shape[0] for value in values), default=0)
    dtype = values[0].dtype if values else torch.long
    out = torch.full((len(values), max_len), pad_value, dtype=dtype)
    for idx, value in enumerate(values):
        out[idx, : value.shape[0]] = value
    return out


def pad_2d(values: List[torch.Tensor], pad_value: int | float = 0) -> torch.Tensor:
    max_len = max((value.shape[0] for value in values), default=0)
    width = values[0].shape[1] if values else 0
    dtype = values[0].dtype if values else torch.long
    out = torch.full((len(values), max_len, width), pad_value, dtype=dtype)
    for idx, value in enumerate(values):
        out[idx, : value.shape[0], :] = value
    return out
