#!/usr/bin/env python3

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
from pathlib import Path

import torch

try:
    import wandb
except ImportError:
    wandb = None

from .config import ExperimentPaths, load_config
from .data import create_dataloaders
from .model import build_model
from .train import evaluate_model, load_checkpoint_if_available, train_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Refactored neural prediction pipeline using HF genomics datasets")
    parser.add_argument("--config", required=True, help="Path to YAML config")
    parser.add_argument("--output-dir", required=True, help="Directory to save results")
    return parser.parse_args()


def create_experiment_paths(output_dir: Path) -> ExperimentPaths:
    models_dir = output_dir / "models"
    reports_dir = output_dir / "reports"
    models_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    return ExperimentPaths(
        root=output_dir,
        models_dir=models_dir,
        reports_dir=reports_dir,
        config_copy_path=output_dir / "config.yaml",
    )


def main() -> None:
    args = parse_args()
    config = load_config(Path(args.config))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    experiment = create_experiment_paths(output_dir)
    shutil.copyfile(Path(args.config), experiment.config_copy_path)

    data_bundle = create_dataloaders(config)
    model = build_model(
        config=config,
        input_shape=data_bundle.dataset.get_input_shape(),
        num_classes=data_bundle.dataset.get_num_classes(),
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    wandb_run = None
    if config.wandb.use_wandb:
        if wandb is None:
            raise ImportError("wandb is enabled in config but the package is not installed")
        wandb_run = wandb.init(
            project=config.wandb.project_name,
            name=config.wandb.run_name,
            config={
                "prediction_target": config.output.prediction_target,
                "genes_to_use": config.dataset_input.genes_to_use,
                "model_type": config.model.type,
                "mode": config.mode,
            },
        )

    if config.mode == "train":
        history = train_model(model, data_bundle, config, device, experiment.models_dir, wandb_run=wandb_run)
        torch.save(model.state_dict(), experiment.models_dir / "final.pt")
    elif config.mode == "test":
        load_checkpoint_if_available(model, config.checkpointing.load_checkpoint, device)
        history = evaluate_model(model, data_bundle, device)
    else:
        raise ValueError(f"Unsupported mode: {config.mode}")

    with (experiment.root / "results.json").open("w", encoding="utf-8") as handle:
        json.dump(history, handle, indent=2)
        handle.write("\n")

    with (experiment.reports_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump({"mode": config.mode, "timestamp": dt.datetime.now().isoformat()}, handle, indent=2)
        handle.write("\n")

    if wandb_run is not None:
        wandb_run.finish()


if __name__ == "__main__":
    main()
