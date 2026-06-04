from __future__ import annotations

import os
import random
from typing import Optional

import numpy as np
import torch
from rich.console import Console


console = Console()


def set_random_seeds(seed: int, strict_determinism: bool = True) -> None:
    """Configure Python, NumPy and PyTorch random seeds."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.cuda.empty_cache()

    if strict_determinism:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        try:
            torch.use_deterministic_algorithms(True)
        except AttributeError:
            torch.set_deterministic(True)
        except Exception as exc:
            console.print(
                "[yellow]Nao foi possivel ativar determinismo estrito no PyTorch. "
                f"Continuando com determinismo parcial. Detalhe: {exc}[/yellow]"
            )
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        console.print(
            f"[green]Seed configurada: {seed} "
            f"(determinismo ESTRITO)[/green]"
        )
    else:
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True
        try:
            torch.use_deterministic_algorithms(False)
        except AttributeError:
            try:
                torch.set_deterministic(False)
            except AttributeError:
                pass
        except Exception:
            pass
        console.print(
            f"[green]Seed configurada: {seed} "
            f"(determinismo PARCIAL)[/green]"
        )


def worker_init_fn(worker_id: int) -> None:
    """Initialize DataLoader worker RNGs from the PyTorch worker seed."""
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def make_torch_generator(seed: Optional[int]) -> Optional[torch.Generator]:
    if seed is None:
        return None
    generator = torch.Generator()
    generator.manual_seed(int(seed))
    return generator
