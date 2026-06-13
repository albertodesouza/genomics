from __future__ import annotations

import argparse
import signal
import subprocess
import sys
from pathlib import Path

import torch
from rich.console import Console

from genomics.predictors.genotype_based.config import load_config
from genomics.predictors.genotype_based.data.pipeline import prepare_data
from .evaluation import run_test_and_save
from .experiment import interrupt_state, setup_experiment_dir
from genomics.predictors.genotype_based.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor, SKLEARN_BASELINE_TYPES
from .training import Trainer
from genomics.predictors.genotype_based.utils import set_random_seeds
from genomics.core import update_manifest
from genomics.core.run_utils import select_device, select_split_loader, training_manifest_fields
from genomics.core.wandb_utils import finish_wandb, init_wandb_if_enabled

console = Console()


def _restore_terminal() -> None:
    if not sys.stdout.isatty():
        return
    try:
        sys.stdout.write("\033[0m\033[?25h")
        sys.stdout.flush()
    except Exception:
        pass
    try:
        with open("/dev/tty") as tty:
            subprocess.run(["stty", "sane"], stdin=tty, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass


def _install_signal_handlers() -> None:
    def _handle_sigint(_signum, _frame):
        interrupt_state.interrupted = True
        _restore_terminal()
        console.print("[yellow]\n⚠ Interrupção solicitada. Encerrando época atual e avaliando teste...[/yellow]")

    signal.signal(signal.SIGINT, _handle_sigint)
    signal.signal(signal.SIGTERM, _handle_sigint)


def _build_model(config, dataset):
    input_shape = dataset.get_input_shape()
    num_classes = dataset.get_num_classes()
    model_type = config.model.type.upper()

    if model_type == "NN":
        return NNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN":
        return CNNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN2":
        return CNN2AncestryPredictor(config, input_shape, num_classes)
    raise ValueError(f"Modelo nao suportado no entrypoint de treino: {config.model.type}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Train genotype_based_predictor models")
    parser.add_argument("config_path", type=str, help="Path to YAML config")
    args = parser.parse_args()

    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    interrupt_state.interrupted = False
    _install_signal_handlers()

    training_seed = config.training.random_seed
    if training_seed is None:
        training_seed = config.data_split.random_seed
    if training_seed is not None and training_seed != -1:
        set_random_seeds(training_seed, config.data_split.strict_determinism)

    device = select_device()
    console.print(f"[green]Device:[/green] {device}")

    wandb_run = None
    experiment_dir = setup_experiment_dir(config, str(config_path))
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    wandb_run = init_wandb_if_enabled(config.wandb, config, experiment_dir, config_path, console)

    try:
        if config.model.type.upper() in SKLEARN_BASELINE_TYPES:
            from genomics.predictors.genotype_based.models.sklearn_models import train_sklearn_baseline

            train_sklearn_baseline(
                config=config,
                model_type=config.model.type,
                train_loader=train_loader,
                val_loader=val_loader,
                test_loader=test_loader,
                full_dataset=full_ds,
                experiment_dir=experiment_dir,
                wandb_run=wandb_run,
            )
            return

        model = _build_model(config, full_ds).to(device)
        trainer = Trainer(
            model=model,
            train_loader=train_loader,
            val_loader=val_loader,
            config=config,
            device=device,
            experiment_dir=experiment_dir,
            wandb_run=wandb_run,
        )
        history = trainer.train()
        update_manifest(
            experiment_dir,
            status="interrupted" if history.get("interrupted") else "completed",
            training_random_seed=training_seed,
            split_random_seed=config.data_split.random_seed,
            **training_manifest_fields(history),
        )

        selected_loader, split_name = select_split_loader("val", train_loader, val_loader, test_loader)
        if history.get("interrupted"):
            console.print(f"[cyan]Executando avaliação pós-interrupção no split '{split_name}'...[/cyan]")
        else:
            best_accuracy_path = experiment_dir / "models" / "best_accuracy.pt"
            if best_accuracy_path.exists():
                checkpoint = torch.load(best_accuracy_path, map_location=device)
                state_dict = checkpoint.get("model_state_dict", checkpoint)
                model.load_state_dict(state_dict)
                console.print(f"[green]✓ Melhor checkpoint carregado para validação final:[/green] {best_accuracy_path}")
            else:
                console.print("[yellow]best_accuracy.pt não encontrado; avaliando pesos finais.[/yellow]")
            console.print(f"[cyan]Executando avaliação final do melhor modelo no split '{split_name}'...[/cyan]")
        output_name = "val_best_accuracy" if not history.get("interrupted") else split_name
        run_test_and_save(model, selected_loader, full_ds, config, device, output_name, experiment_dir, wandb_run)
    finally:
        _restore_terminal()
        if "experiment_dir" in locals():
            update_manifest(experiment_dir, wandb_enabled=wandb_run is not None)
        finish_wandb(wandb_run)


if __name__ == "__main__":
    main()
