from __future__ import annotations

import argparse
from pathlib import Path

import torch
from rich.console import Console

from genomics.predictors.genotype_based.config import generate_experiment_name, get_experiment_runs_dir, load_config
from genomics.predictors.genotype_based.data.pipeline import prepare_data
from .evaluation import run_test_and_save
from genomics.predictors.genotype_based.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor, SKLEARN_BASELINE_TYPES

console = Console()


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
    raise ValueError(f"Modelo PyTorch nao suportado: {config.model.type}")


def _select_loader(config, split_name: str, train_loader, val_loader, test_loader):
    split = split_name.lower()
    if split == "train":
        return train_loader
    if split == "val":
        return val_loader
    if split == "test":
        return test_loader
    raise ValueError("--split deve ser train, val ou test")


def _resolve_checkpoint(experiment_dir: Path, checkpoint: str) -> Path:
    path = Path(checkpoint)
    if path.exists():
        return path.resolve()
    if path.suffix != ".pt":
        path = path.with_suffix(".pt")
    return experiment_dir / "models" / path.name


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate a saved genotype_based_predictor checkpoint")
    parser.add_argument("config_path", type=Path, help="Path to the YAML config used for training")
    parser.add_argument("--checkpoint", default="best_accuracy", help="Checkpoint alias/path: best_accuracy, best_loss, final or .pt path")
    parser.add_argument("--split", default="test", choices=["train", "val", "test"], help="Split to evaluate")
    parser.add_argument("--experiment-dir", type=Path, default=None, help="Experiment directory; inferred from config if omitted")
    parser.add_argument("--output-name", default=None, help="Result JSON prefix; default is <split>_<checkpoint_stem>")
    args = parser.parse_args()

    config_path = args.config_path.resolve()
    config = load_config(config_path)
    experiment_dir = args.experiment_dir or (get_experiment_runs_dir(config) / generate_experiment_name(config))
    experiment_dir = experiment_dir.resolve()
    checkpoint_path = _resolve_checkpoint(experiment_dir, args.checkpoint)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    console.print(f"[green]Experimento:[/green] {experiment_dir}")
    if config.model.type.upper() in SKLEARN_BASELINE_TYPES:
        console.print(f"[green]Artefato sklearn:[/green] {experiment_dir / 'models' / 'sklearn_baseline.joblib'}")
    else:
        console.print(f"[green]Checkpoint:[/green] {checkpoint_path}")
    console.print(f"[green]Split:[/green] {args.split}")
    console.print(f"[green]Device:[/green] {device}")

    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    if config.model.type.upper() in SKLEARN_BASELINE_TYPES:
        from genomics.predictors.genotype_based.models.sklearn_models import run_sklearn_test_mode

        config.test_dataset = args.split
        run_sklearn_test_mode(
            config,
            train_loader,
            val_loader,
            test_loader,
            full_ds,
            experiment_dir,
            wandb_run=None,
            output_name=args.output_name,
        )
        return

    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint nao encontrado: {checkpoint_path}")

    loader = _select_loader(config, args.split, train_loader, val_loader, test_loader)
    model = _build_model(config, full_ds).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device)
    state_dict = checkpoint.get("model_state_dict", checkpoint)
    model.load_state_dict(state_dict)

    checkpoint_stem = checkpoint_path.stem
    output_name = args.output_name or f"{args.split}_{checkpoint_stem}"
    run_test_and_save(model, loader, full_ds, config, device, output_name, experiment_dir, wandb_run=None)


if __name__ == "__main__":
    main()
