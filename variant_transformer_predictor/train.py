from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import torch
from rich.console import Console
from torch.utils.data import DataLoader

from variant_transformer_predictor.config import get_experiment_dir, load_config
from variant_transformer_predictor.dataset import VariantTokenDataset, collate_variant_tokens
from variant_transformer_predictor.evaluation import evaluate, save_results
from variant_transformer_predictor.model import VariantTransformerClassifier
from variant_transformer_predictor.training import Trainer

console = Console()


def _make_loaders(config):
    kwargs = dict(
        processed_dir=config.dataset.processed_dir,
        max_sequence_length=config.dataset.max_sequence_length,
        truncate_policy=config.dataset.truncate_policy,
        loading_strategy=config.dataset.loading_strategy,
    )
    train_ds = VariantTokenDataset(split="train", **kwargs)
    val_ds = VariantTokenDataset(split="val", **kwargs)
    test_ds = VariantTokenDataset(split="test", **kwargs)
    train_loader = DataLoader(train_ds, batch_size=config.training.batch_size, shuffle=True, num_workers=config.training.num_workers, collate_fn=collate_variant_tokens)
    val_loader = DataLoader(val_ds, batch_size=config.training.batch_size, shuffle=False, num_workers=config.training.num_workers, collate_fn=collate_variant_tokens)
    test_loader = DataLoader(test_ds, batch_size=config.training.batch_size, shuffle=False, num_workers=config.training.num_workers, collate_fn=collate_variant_tokens)
    return train_ds, train_loader, val_loader, test_loader


def main() -> int:
    parser = argparse.ArgumentParser(description="Treina Transformer esparso baseado em variantes")
    parser.add_argument("config_path")
    args = parser.parse_args()
    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    if config.data_split.random_seed is not None:
        torch.manual_seed(int(config.data_split.random_seed))
    experiment_dir = get_experiment_dir(config)
    experiment_dir.mkdir(parents=True, exist_ok=True)
    (experiment_dir / "models").mkdir(exist_ok=True)
    shutil.copyfile(config_path, experiment_dir / "config.yaml")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    full_ds, train_loader, val_loader, test_loader = _make_loaders(config)
    num_genes = int(full_ds.metadata["num_genes"])
    num_classes = full_ds.get_num_classes()
    model = VariantTransformerClassifier(config.model, num_genes, num_classes, int(full_ds.metadata["l_max"])).to(device)
    console.print(f"[green]Experimento:[/green] {experiment_dir}")
    console.print(f"[green]Device:[/green] {device}")
    console.print(f"[green]Genes/classes:[/green] {num_genes}/{num_classes}")
    wandb_run = None
    if config.wandb.use_wandb:
        import wandb
        wandb_run = wandb.init(project=config.wandb.project_name, name=config.wandb.run_name or experiment_dir.name, config=config.model_dump(mode="python"))
    try:
        trainer = Trainer(model, train_loader, val_loader, config, device, experiment_dir, full_ds.get_class_names(), wandb_run)
        history = trainer.train()
        with open(experiment_dir / "models" / "training_history.json", "w") as f:
            json.dump(history, f, indent=2)
        ckpt_path = experiment_dir / "models" / "best_accuracy.pt"
        if ckpt_path.exists():
            state = torch.load(ckpt_path, map_location=device, weights_only=True)
            model.load_state_dict(state["model_state_dict"])
        results = evaluate(model, test_loader, full_ds.get_class_names(), device, "test_best_accuracy")
        if results:
            save_results(results, experiment_dir / "test_results.json")
    finally:
        if wandb_run is not None:
            wandb_run.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
