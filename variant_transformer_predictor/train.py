from __future__ import annotations

import argparse
import json
from pathlib import Path

from rich.console import Console
from torch.utils.data import DataLoader

from variant_transformer_predictor.config import get_experiment_dir, load_config
from variant_transformer_predictor.dataset import VariantTokenDataset, collate_variant_tokens
from variant_transformer_predictor.evaluation import evaluate, save_results
from variant_transformer_predictor.model import VariantTransformerClassifier
from variant_transformer_predictor.training import Trainer
from genomics_pipeline import set_random_seeds, setup_experiment_run, update_manifest
from genomics_pipeline.checkpointing import load_checkpoint
from genomics_pipeline.data_loading import dataloader_kwargs
from genomics_pipeline.run_utils import select_device, training_manifest_fields
from genomics_pipeline.wandb_utils import finish_wandb, init_wandb_if_enabled

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
    common = dict(
        batch_size=config.training.batch_size,
        num_workers=config.training.num_workers,
        collate_fn=collate_variant_tokens,
    )
    train_loader = DataLoader(train_ds, **dataloader_kwargs(shuffle=True, **common))
    val_loader = DataLoader(val_ds, **dataloader_kwargs(shuffle=False, **common))
    test_loader = DataLoader(test_ds, **dataloader_kwargs(shuffle=False, **common))
    return train_ds, train_loader, val_loader, test_loader


def main() -> int:
    parser = argparse.ArgumentParser(description="Treina Transformer esparso baseado em variantes")
    parser.add_argument("config_path")
    args = parser.parse_args()
    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    if config.data_split.random_seed is not None:
        set_random_seeds(int(config.data_split.random_seed), strict_determinism=False)
    experiment_dir = get_experiment_dir(config)
    run = setup_experiment_run(
        pipeline="variant_transformer_predictor",
        run_name=experiment_dir.name,
        runs_root=experiment_dir.parent,
        config_path=config_path,
        resolved_config=config.model_dump(mode="python"),
        extra={"processed_dir": config.dataset.processed_dir},
    )
    experiment_dir = run.run_dir
    device = select_device()
    full_ds, train_loader, val_loader, test_loader = _make_loaders(config)
    num_genes = int(full_ds.metadata["num_genes"])
    num_classes = full_ds.get_num_classes()
    model = VariantTransformerClassifier(config.model, num_genes, num_classes, int(full_ds.metadata["l_max"])).to(device)
    console.print(f"[green]Experimento:[/green] {experiment_dir}")
    console.print(f"[green]Device:[/green] {device}")
    console.print(f"[green]Genes/classes:[/green] {num_genes}/{num_classes}")
    wandb_run = init_wandb_if_enabled(config.wandb, config, experiment_dir, config_path, console)
    try:
        trainer = Trainer(model, train_loader, val_loader, config, device, experiment_dir, full_ds.get_class_names(), wandb_run)
        history = trainer.train()
        with open(experiment_dir / "models" / "training_history.json", "w") as f:
            json.dump(history, f, indent=2)
        ckpt_path = experiment_dir / "models" / "best_accuracy.pt"
        if ckpt_path.exists():
            load_checkpoint(ckpt_path, model, device=device, weights_only=True)
        results = evaluate(model, test_loader, full_ds.get_class_names(), device, "test_best_accuracy")
        if results:
            save_results(results, experiment_dir / "test_results.json")
        update_manifest(
            experiment_dir,
            status="completed",
            **training_manifest_fields(history),
            test_results="test_results.json" if results else None,
        )
    finally:
        update_manifest(experiment_dir, wandb_enabled=wandb_run is not None)
        finish_wandb(wandb_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
