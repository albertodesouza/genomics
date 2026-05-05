from __future__ import annotations

import argparse
import signal
from pathlib import Path

import torch
from rich.console import Console

from genotype_based_predictor.config import load_config
from genotype_based_predictor.data_pipeline import prepare_data
from genotype_based_predictor.evaluation import run_test_and_save
from genotype_based_predictor.experiment import interrupt_state, setup_experiment_dir
from genotype_based_predictor.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor
from genotype_based_predictor.training import Trainer
from genotype_based_predictor.utils import set_random_seeds

console = Console()


def _install_signal_handlers() -> None:
    def _handle_sigint(_signum, _frame):
        interrupt_state.interrupted = True
        console.print("[yellow]\n⚠ Interrupção solicitada. Encerrando época atual e avaliando teste...[/yellow]")

    signal.signal(signal.SIGINT, _handle_sigint)


def _select_test_loader(config, train_loader, val_loader, test_loader):
    choice = config.test_dataset.lower()
    if choice == "train":
        return train_loader, "train"
    if choice == "val":
        return val_loader, "val"
    return test_loader, "test"


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

    if config.data_split.random_seed is not None and config.data_split.random_seed != -1:
        set_random_seeds(config.data_split.random_seed, config.data_split.strict_determinism)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    console.print(f"[green]Device:[/green] {device}")

    experiment_dir = setup_experiment_dir(config, str(config_path))
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    _install_signal_handlers()

    model = _build_model(config, full_ds).to(device)
    trainer = Trainer(
        model=model,
        train_loader=train_loader,
        val_loader=val_loader,
        config=config,
        device=device,
        experiment_dir=experiment_dir,
        wandb_run=None,
    )
    history = trainer.train()

    if history.get("interrupted"):
        selected_loader, split_name = _select_test_loader(config, train_loader, val_loader, test_loader)
        console.print(f"[cyan]Executando avaliação pós-interrupção no split '{split_name}'...[/cyan]")
        run_test_and_save(model, selected_loader, full_ds, config, device, split_name, experiment_dir)


if __name__ == "__main__":
    main()
