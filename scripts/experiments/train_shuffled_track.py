#!/usr/bin/env python3
"""Train+evaluate one gene's CNN arm with its RNA-seq track permuted along position.

Reviewer ask (Q2b): does the classifier need the track's real spatial arrangement, or
would any arbitrary-but-consistent remapping of positions do just as well because the
model is really reading an aggregate, ancestry-correlated signal level rather than
position-specific structure?

ONE FIXED PERMUTATION PER GENE, derived deterministically from the gene name, applied
identically to every individual and every split (train/val/test) and to both haplotypes.
This is the key design choice: the permutation must be the same for everyone, or the
control would also destroy cross-individual comparability and confound "no positional
structure" with "no coherent signal at all". Since it's a fixed relabelling of the 32,768
position axis, every individual's own value distribution is untouched -- only which
position index a value sits at changes, identically for all.

Reuses the EXACT same cached, aligned tensors as the published CNN arm (via
`prepare_data`, the same function `genomics genotype train` calls): no new AlphaGenome
calls, no re-alignment, no new cache written to disk. The permutation is applied at
`__getitem__` time by a thin Dataset wrapper, so the on-disk cache for the gene is
untouched and can still be read by the published (unshuffled) arm concurrently.

The training/eval loop below is a deliberately narrow copy of
experiments/train.py's main() plus the val/test evaluation experiments/
evaluate_checkpoint.py does separately -- both are needed here in one process because the
standard `genomics genotype test` CLI would reload UNSHUFFLED data from the config,
which would evaluate a model trained on shuffled input against unshuffled input.

Usage:
  python3 scripts/experiments/train_shuffled_track.py <config.yaml>
"""
from __future__ import annotations

import argparse
import hashlib
import signal
import sys
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


class PermutedDataset(Dataset):
    """Wraps a (features, target, idx) dataset, permuting features along the last axis."""

    def __init__(self, base, perm: np.ndarray):
        self.base = base
        self.perm = perm
        for name in ("target_to_idx", "idx_to_target", "config"):
            if hasattr(base, name):
                setattr(self, name, getattr(base, name))

    def __len__(self):
        return len(self.base)

    def __getitem__(self, idx):
        features, target, orig_idx = self.base[idx]
        return features[..., self.perm], target, orig_idx

    def get_num_classes(self):
        return self.base.get_num_classes() if hasattr(self.base, "get_num_classes") else len(self.idx_to_target)

    def get_input_shape(self):
        return self.base.get_input_shape() if hasattr(self.base, "get_input_shape") else self[0][0].shape


def _gene_seed(gene: str) -> int:
    return int(hashlib.sha256(gene.encode()).hexdigest()[:8], 16) % (2**31)


def _permuted_loader(loader: DataLoader, perm: np.ndarray, shuffle: bool) -> DataLoader:
    ds = PermutedDataset(loader.dataset, perm)
    return DataLoader(ds, batch_size=loader.batch_size, shuffle=shuffle,
                       collate_fn=_collate_fn, num_workers=loader.num_workers)


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
    raise ValueError(f"Unsupported model type for the shuffled-track control: {config.model.type}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_path", type=str)
    parser.add_argument("--gene", default=None, help="Defaults to config.dataset_input.genes_to_use[0]")
    args = parser.parse_args()

    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    gene = args.gene or config.dataset_input.genes_to_use[0]
    interrupt_state.interrupted = False

    def _handle(_s, _f):
        interrupt_state.interrupted = True
    signal.signal(signal.SIGINT, _handle)
    signal.signal(signal.SIGTERM, _handle)

    training_seed = config.training.random_seed
    if training_seed is None:
        training_seed = config.data_split.random_seed
    if training_seed is not None and training_seed != -1:
        set_random_seeds(training_seed, config.data_split.strict_determinism)

    device = select_device()
    console.print(f"[green]Device:[/green] {device}  [green]Gene:[/green] {gene}")

    experiment_dir = setup_experiment_dir(config, str(config_path))
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)

    sample_features, _, _ = train_loader.dataset[0]
    length = sample_features.shape[-1]
    perm = np.random.default_rng(_gene_seed(gene)).permutation(length)
    console.print(f"[cyan]Position permutation:[/cyan] length={length}, seed={_gene_seed(gene)}, "
                  f"identity={bool(np.array_equal(perm, np.arange(length)))}")

    train_loader = _permuted_loader(train_loader, perm, shuffle=True)
    val_loader = _permuted_loader(val_loader, perm, shuffle=False)
    test_loader = _permuted_loader(test_loader, perm, shuffle=False)

    model = _build_model(config, full_ds).to(device)
    trainer = Trainer(model=model, train_loader=train_loader, val_loader=val_loader,
                       config=config, device=device, experiment_dir=experiment_dir, wandb_run=None)
    history = trainer.train()
    update_manifest(experiment_dir,
                     status="interrupted" if history.get("interrupted") else "completed",
                     training_random_seed=training_seed,
                     split_random_seed=config.data_split.random_seed,
                     shuffled_track_perm_seed=_gene_seed(gene),
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
