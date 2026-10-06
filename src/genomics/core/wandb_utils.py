from __future__ import annotations

from pathlib import Path
from typing import Any, Optional


def init_wandb_if_enabled(
    wandb_config: Any,
    resolved_config: Any,
    experiment_dir: Path,
    config_path: Optional[Path] = None,
    console: Any = None,
) -> Any:
    if not getattr(wandb_config, "use_wandb", False):
        return None
    try:
        import wandb
    except ImportError as exc:
        raise ImportError("W&B habilitado, mas o pacote 'wandb' nao esta instalado") from exc

    run_name = getattr(wandb_config, "run_name", None) or Path(experiment_dir).name
    project_name = getattr(wandb_config, "project_name")
    if hasattr(resolved_config, "model_dump"):
        config_payload = resolved_config.model_dump(mode="python")
    else:
        config_payload = resolved_config
    run = wandb.init(project=project_name, name=run_name, config=config_payload, dir=str(experiment_dir))
    if config_path is not None:
        run.config.update({"config_path": str(config_path), "experiment_dir": str(experiment_dir)}, allow_val_change=True)
    if console is not None:
        console.print(f"[green]✓ W&B habilitado:[/green] {project_name}/{run_name}")
    return run


def finish_wandb(run: Any) -> None:
    if run is not None:
        run.finish()
