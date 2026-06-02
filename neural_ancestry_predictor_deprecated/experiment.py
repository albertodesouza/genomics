from __future__ import annotations

from pathlib import Path
from typing import Dict

from rich.console import Console

from genomics_pipeline import setup_experiment_run
from neural_ancestry_predictor_deprecated.config import generate_experiment_name, get_results_dir


console = Console()


def setup_experiment_dir(config: Dict, config_path: str) -> Path:
    """Cria e configura diretório do experimento neural."""
    experiment_name = generate_experiment_name(config)

    run = setup_experiment_run(
        pipeline="neural_ancestry_predictor_deprecated",
        run_name=experiment_name,
        runs_root=get_results_dir(config),
        config_path=Path(config_path),
        resolved_config=config,
        extra={"dataset_cache_root": config["dataset_input"].get("processed_cache_dir")},
    )
    experiment_dir = run.run_dir

    console.print(f"[green]📁 Diretório do experimento:[/green] {experiment_dir}")
    console.print(f"[green]   Nome:[/green] {experiment_name}")

    return experiment_dir
