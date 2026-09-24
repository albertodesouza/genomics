#!/usr/bin/env python3
"""Train+evaluate one gene's arm on the un-aligned (no DynamicIndelAligner) RNA-seq
signal -- reviewer Q1, scoped to "no alignment" only, per user request.

Same published config, same split/labels (via `prepare_data(config)`, only used for
those -- not for features), same architecture and hyperparameters, same log1p
normalisation procedure (Eq. 2 of the paper). The ONLY thing that differs from the
published arm: the (2, 32768) tensor is built from
results/cache/genotype_based_predictor/noalign_signal/<GENE>.npz (a naive crop of each
individual's own raw AlphaGenome output, see noalign_export_cache.py) instead of from the
DynamicIndelAligner-aligned cache. The published config's tensor_layout, kernel/stride and
input shape are all already exactly (2, 32768) for a single-ontology single-strand arm, so
none of that needs to change here -- only the feature source.

The log1p divisor is refit here (max over the TRAIN split only, pooled over haplotype and
position, matching processed_dataset.py's own fit procedure) rather than reused from the
published run, because it must be a property of the un-aligned signal's own train-split
distribution, not of the aligned one.

Usage:
  python3 scripts/experiments/train_noalign.py <config.yaml> --gene OCA2
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
from genomics.predictors.genotype_based.models import CNN2AncestryPredictor
from genomics.predictors.genotype_based.utils import set_random_seeds

console = Console()
REPO = Path("/home/breno/I2CA/genomics")
NOALIGN_CACHE_DIR = REPO / "results/cache/genotype_based_predictor/noalign_signal"


class NoAlignDataset(Dataset):
    """Replaces the aligned RNA-seq tensor with the gene's naive-crop (no indel
    remapping) cache, log1p-normalised with a divisor fit on the train split alone."""

    def __init__(self, base, cache: dict, log_max: float):
        self.base = base
        self.cache = cache
        self.log_max = log_max
        for name in ("target_to_idx", "idx_to_target", "config"):
            if hasattr(base, name):
                setattr(self, name, getattr(base, name))

    def __len__(self):
        return len(self.base)

    def __getitem__(self, idx):
        _features, target, orig_idx = self.base[idx]
        sample_id = self.base.get_sample_id(idx)
        h1 = self.cache[f"{sample_id}__H1"]
        h2 = self.cache[f"{sample_id}__H2"]
        stacked = np.log1p(np.stack([h1, h2], axis=0)) / self.log_max  # (2, 32768)
        return torch.from_numpy(stacked.astype(np.float32)), target, orig_idx


def _fit_log_max(train_loader: DataLoader, cache: dict) -> float:
    """Max over the train split only, pooled over haplotype and position -- same
    reduction processed_dataset.py's own log-normalization fit uses."""
    base = train_loader.dataset
    xmax = 0.0
    for idx in range(len(base)):
        sample_id = base.get_sample_id(idx)
        for h in ("H1", "H2"):
            v = cache[f"{sample_id}__{h}"]
            m = float(v.max())
            if m > xmax:
                xmax = m
    return float(np.log1p(xmax)) if xmax > 0 else 1.0


def _loader_with_noalign(loader: DataLoader, cache: dict, log_max: float, shuffle: bool) -> DataLoader:
    ds = NoAlignDataset(loader.dataset, cache, log_max)
    return DataLoader(ds, batch_size=loader.batch_size, shuffle=shuffle,
                       collate_fn=_collate_fn, num_workers=loader.num_workers)


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

    if config.model.type.upper() != "CNN2":
        raise SystemExit(f"ABORT: no-alignment control only implemented for CNN2, got {config.model.type}")

    cache_path = NOALIGN_CACHE_DIR / f"{gene}.npz"
    if not cache_path.exists():
        raise SystemExit(f"ABORT: {cache_path} missing; run noalign_export_cache.py for {gene} first")
    cache = dict(np.load(cache_path))
    console.print(f"[cyan]No-alignment cache:[/cyan] {cache_path} ({len(cache)//2} individuals)")

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
    input_shape = (2, 32768)  # H1+H2, single channel each -- identical to the published arm

    log_max = _fit_log_max(train_loader, cache)
    console.print(f"[cyan]log1p divisor (train-split max):[/cyan] log_max={log_max:.4f}")

    train_loader = _loader_with_noalign(train_loader, cache, log_max, shuffle=True)
    val_loader = _loader_with_noalign(val_loader, cache, log_max, shuffle=False)
    test_loader = _loader_with_noalign(test_loader, cache, log_max, shuffle=False)

    model = CNN2AncestryPredictor(config, input_shape, num_classes).to(device)
    trainer = Trainer(model=model, train_loader=train_loader, val_loader=val_loader,
                       config=config, device=device, experiment_dir=experiment_dir, wandb_run=None)
    history = trainer.train()
    update_manifest(experiment_dir,
                     status="interrupted" if history.get("interrupted") else "completed",
                     training_random_seed=training_seed,
                     split_random_seed=config.data_split.random_seed,
                     no_alignment_control=True,
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
