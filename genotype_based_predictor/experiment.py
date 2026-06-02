# -*- coding: utf-8 -*-
"""
experiment.py — Setup de diretório de experimento e estado de interrupção.
"""

from pathlib import Path

from rich.console import Console

from genotype_based_predictor.config import PipelineConfig, generate_experiment_name, get_experiment_runs_dir
from genomics_pipeline import setup_experiment_run

console = Console()


class InterruptState:
    """
    Flag global para capturar interrupção graciosamente (CTRL+C).

    O Trainer verifica ``interrupted`` a cada época para encerrar
    salvando o checkpoint antes de sair.
    """
    interrupted: bool = False


interrupt_state = InterruptState()


def setup_experiment_dir(config: PipelineConfig, config_path: str) -> Path:
    """
    Cria e configura o diretório do experimento.

    O nome do diretório é gerado automaticamente a partir dos parâmetros
    de configuração para que experimentos distintos nunca se sobrescrevam.

    Estrutura criada
    ----------------
    ::

        <results_dir>/
        └── <experiment_name>/
            ├── config.yaml   ← cópia do YAML usado
            └── models/       ← checkpoints e artefatos

    Parameters
    ----------
    config : PipelineConfig
        Configuração validada do experimento.
    config_path : str
        Caminho do arquivo YAML (para cópia de referência).

    Returns
    -------
    Path
        Caminho do diretório do experimento criado.
    """
    experiment_name = generate_experiment_name(config)
    run = setup_experiment_run(
        pipeline="genotype_based_predictor",
        run_name=experiment_name,
        runs_root=get_experiment_runs_dir(config),
        config_path=Path(config_path),
        resolved_config=config.model_dump(mode="python"),
        extra={"dataset_cache_root": str(Path(config.dataset_input.processed_cache_dir))},
    )
    experiment_dir = run.run_dir

    console.print(f"[green]📁 Experimento:[/green] {experiment_dir}")
    console.print(f"[green]   Nome:[/green] {experiment_name}")

    return experiment_dir
