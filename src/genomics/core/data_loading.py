from __future__ import annotations

from typing import Any, Callable, Optional

import torch


def make_data_loader_generator(seed: Optional[int]) -> Optional[torch.Generator]:
    if seed is None:
        return None
    generator = torch.Generator()
    generator.manual_seed(int(seed))
    return generator


def dataloader_kwargs(
    *,
    batch_size: int,
    shuffle: bool,
    num_workers: int = 0,
    pin_memory: bool = True,
    collate_fn: Optional[Callable[..., Any]] = None,
    generator: Optional[torch.Generator] = None,
    worker_init_fn: Optional[Callable[[int], None]] = None,
    persistent_workers: bool = False,
    prefetch_factor: Optional[int] = None,
) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "batch_size": batch_size,
        "shuffle": shuffle,
        "num_workers": int(num_workers),
        "pin_memory": bool(pin_memory),
    }
    if collate_fn is not None:
        kwargs["collate_fn"] = collate_fn
    if generator is not None:
        kwargs["generator"] = generator
    if worker_init_fn is not None and num_workers > 0:
        kwargs["worker_init_fn"] = worker_init_fn
    if num_workers > 0:
        kwargs["persistent_workers"] = bool(persistent_workers)
        if prefetch_factor is not None:
            kwargs["prefetch_factor"] = prefetch_factor
    return kwargs
