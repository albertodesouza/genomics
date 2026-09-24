#!/usr/bin/env python3
"""Train+evaluate one gene's arm on raw one-hot DNA sequence, no AlphaGenome.

Reviewer ask (Q2c): feed raw sequence directly (no AlphaGenome) with an otherwise
identical CNN, to isolate the added value of the predicted-functional gate.

Uses the same aligned coordinate axis the RNA-seq arm reads (DynamicIndelAligner, via
scripts/experiments/raw_sequence_export_cache.py -- run that first for this gene) and the
same train/val/test split and labels as the published arm (via `prepare_data(config)` on
the SAME production config), so the only thing that differs from the published arm is the
input: one-hot A/C/G/T (4 channels) per haplotype instead of 1-channel predicted RNA-seq.

The CNN architecture is unchanged except stage 1's kernel/stride height, widened from 1
to 4 so it still spans exactly one haplotype's channels (same convention
single_gene_sweep_melanocyte_strand*.py uses for a restricted track count) -- the config
passed in must already have kernel_stage1/stride_stage1 = [4, ...] (see
raw_sequence_sweep.py's config writer).

Usage:
  python3 scripts/experiments/train_raw_sequence.py <config.yaml> --gene OCA2
"""
from __future__ import annotations

import argparse
import signal
from pathlib import Path

import numpy as np
import torch
from rich.console import Console
from torch.utils.data import DataLoader, Dataset

from genomics.core import update_manifest
from genomics.core.run_utils import select_device, training_manifest_fields
from genomics.predictors.genotype_based.config import load_config
from genomics.predictors.genotype_based.data.pipeline import _collate_fn, prepare_data
from genomics.predictors.genotype_based.experiments.evaluation import run_test_and_save
from genomics.predictors.genotype_based.experiments.experiment import interrupt_state, setup_experiment_dir
from genomics.predictors.genotype_based.experiments.training import Trainer
from genomics.predictors.genotype_based.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor
from genomics.predictors.genotype_based.utils import set_random_seeds

console = Console()
REPO = Path("/home/breno/I2CA/genomics")
ONEHOT_CACHE_DIR = REPO / "results/cache/genotype_based_predictor/raw_sequence_onehot"


class RawSequenceDataset(Dataset):
    """Replaces AlphaGenome-derived features with the gene's cached one-hot DNA tensor,
    keeping the base dataset's labels and (sample, split) membership untouched."""

    def __init__(self, base, onehot: dict):
        self.base = base
        self.onehot = onehot
        for name in ("target_to_idx", "idx_to_target", "config"):
            if hasattr(base, name):
                setattr(self, name, getattr(base, name))

    def __len__(self):
        return len(self.base)

    def __getitem__(self, idx):
        _features, target, orig_idx = self.base[idx]
        sample_id = self.base.get_sample_id(idx)
        h1 = self.onehot[f"{sample_id}__H1"]
        h2 = self.onehot[f"{sample_id}__H2"]
        stacked = torch.from_numpy(np.stack([h1, h2], axis=0))  # (2 hap, 4 channel, 32768)
        return stacked, target, orig_idx


def _loader_with_onehot(loader: DataLoader, onehot: dict, shuffle: bool) -> DataLoader:
    ds = RawSequenceDataset(loader.dataset, onehot)
    return DataLoader(ds, batch_size=loader.batch_size, shuffle=shuffle,
                       collate_fn=_collate_fn, num_workers=loader.num_workers)


def _build_model_with_shape(config, input_shape, num_classes):
    model_type = config.model.type.upper()
    if model_type == "NN":
        return NNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN":
        return CNNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN2":
        return CNN2AncestryPredictor(config, input_shape, num_classes)
    raise ValueError(f"Unsupported model type for the raw-sequence control: {config.model.type}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_path", type=str)
    parser.add_argument("--gene", required=True)
    args = parser.parse_args()

    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    gene = args.gene.upper()
    interrupt_state.interrupted = False

    def _handle(_s, _f):
        interrupt_state.interrupted = True
    signal.signal(signal.SIGINT, _handle)
    signal.signal(signal.SIGTERM, _handle)

    onehot_path = ONEHOT_CACHE_DIR / f"{gene}.npz"
    if not onehot_path.exists():
        raise SystemExit(f"ABORT: {onehot_path} missing; run raw_sequence_export_cache.py for {gene} first")
    onehot = dict(np.load(onehot_path))
    console.print(f"[cyan]One-hot cache:[/cyan] {onehot_path} ({len(onehot)//2} individuals)")

    kernel_h = config.model.cnn2.kernel_stage1[0] if config.model.type.upper() == "CNN2" else None
    if kernel_h != 4:
        raise SystemExit(f"ABORT: expected model.cnn2.kernel_stage1 height 4 for one-hot DNA, got {kernel_h}")

    training_seed = config.training.random_seed
    if training_seed is None:
        training_seed = config.data_split.random_seed
    if training_seed is not None and training_seed != -1:
        set_random_seeds(training_seed, config.data_split.strict_determinism)

    device = select_device()
    console.print(f"[green]Device:[/green] {device}  [green]Gene:[/green] {gene}")

    experiment_dir = setup_experiment_dir(config, str(config_path))
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    num_classes = full_ds.get_num_classes()
    input_shape = (8, 32768)  # 2 haplotypes x 4 one-hot channels, flattened row axis

    train_loader = _loader_with_onehot(train_loader, onehot, shuffle=True)
    val_loader = _loader_with_onehot(val_loader, onehot, shuffle=False)
    test_loader = _loader_with_onehot(test_loader, onehot, shuffle=False)

    model = _build_model_with_shape(config, input_shape, num_classes).to(device)
    trainer = Trainer(model=model, train_loader=train_loader, val_loader=val_loader,
                       config=config, device=device, experiment_dir=experiment_dir, wandb_run=None)
    history = trainer.train()
    update_manifest(experiment_dir,
                     status="interrupted" if history.get("interrupted") else "completed",
                     training_random_seed=training_seed,
                     split_random_seed=config.data_split.random_seed,
                     raw_sequence_control=True,
                     **training_manifest_fields(history))

    best_accuracy_path = experiment_dir / "models" / "best_accuracy.pt"
    if best_accuracy_path.exists():
        checkpoint = torch.load(best_accuracy_path, map_location=device)
        model.load_state_dict(checkpoint.get("model_state_dict", checkpoint))
        console.print(f"[green]Loaded best checkpoint:[/green] {best_accuracy_path}")
    else:
        console.print("[yellow]best_accuracy.pt not found; evaluating final weights.[/yellow]")

    run_test_and_save(model, val_loader, full_ds, config, device, "val_best_accuracy", experiment_dir, None)
    run_test_and_save(model, test_loader, full_ds, config, device, "test_best_accuracy", experiment_dir, None)
    console.print(f"[green]Done:[/green] {experiment_dir}")


if __name__ == "__main__":
    main()
