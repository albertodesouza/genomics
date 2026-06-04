# -*- coding: utf-8 -*-
"""
training.py — Trainer: loop de treinamento, validação, checkpoints e W&B.
"""

import time
from pathlib import Path
from typing import Any, Dict, Optional

import torch
import torch.nn as nn
import torch.optim as optim
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from torch.utils.data import DataLoader

from genomics.predictors.genotype_based.config import PipelineConfig
from genomics.core.checkpointing import load_checkpoint as load_torch_checkpoint
from genomics.core.checkpointing import resolve_checkpoint_path, save_checkpoint as save_torch_checkpoint
from genomics.core.optim import make_optimizer_from_config
from genomics.core.torch_utils import move_to_device
from genomics.core.training_utils import EpochTrainer, make_lr_scheduler

console = Console()


def _dataset_profile_stats(dataset: Any) -> Optional[Dict[str, float]]:
    seen = set()
    current = dataset
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        stats = getattr(current, "profile_stats", None)
        if isinstance(stats, dict):
            return stats
        current = getattr(current, "dataset", None)
    return None


def _make_optimizer(model: nn.Module, config: PipelineConfig) -> optim.Optimizer:
    return make_optimizer_from_config(model, config.training, foreach=False)


def _make_scheduler(optimizer: optim.Optimizer, config: PipelineConfig):
    try:
        return make_lr_scheduler(optimizer, config.training.lr_scheduler, default_t_max=config.training.num_epochs)
    except ValueError as exc:
        console.print(f"[yellow]{exc}. Sem scheduler.[/yellow]")
        return None


class Trainer(EpochTrainer):
    """
    Gerencia o loop completo de treinamento com validação por época.

    Funcionalidades
    ---------------
    - Salva checkpoints ``best_accuracy.pt`` e ``best_loss.pt``
    - Suporte a todos os LR schedulers (plateau, cosine, step, etc.)
    - Logging no W&B (opcional)
    - Encerramento gracioso via ``InterruptState`` (CTRL+C)

    Parameters
    ----------
    model : nn.Module
    train_loader : DataLoader
    val_loader : DataLoader
    config : PipelineConfig
    device : torch.device
    experiment_dir : Path
    wandb_run : Any, optional
    """

    def __init__(
        self,
        model: nn.Module,
        train_loader: DataLoader,
        val_loader: DataLoader,
        config: PipelineConfig,
        device: torch.device,
        experiment_dir: Path,
        wandb_run: Optional[Any] = None,
    ):
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.config = config
        self.device = device
        self.experiment_dir = experiment_dir
        self.wandb_run = wandb_run

        self.optimizer = _make_optimizer(model, config)
        self.scheduler = _make_scheduler(self.optimizer, config)

        lf = config.training.loss_function
        if lf == "cross_entropy":
            self.criterion = nn.CrossEntropyLoss()
        elif lf == "mse":
            self.criterion = nn.MSELoss()
        else:
            raise ValueError(f"Loss não suportada: {lf}")

        self.is_classification = config.output.prediction_target != "frog_likelihood"
        self.max_samples_per_epoch = config.debug.max_samples_per_epoch
        self.best_val_loss = float("inf")
        self.best_val_accuracy = 0.0
        self.start_epoch = 1

        models_dir = experiment_dir / "models"
        models_dir.mkdir(parents=True, exist_ok=True)
        self.models_dir = models_dir

        if config.training.weight_decay > 0:
            console.print(f"[green]✓ Weight decay (L2): {config.training.weight_decay}[/green]")
        if config.training.lr_scheduler.enabled:
            console.print(f"[green]✓ LR Scheduler: {config.training.lr_scheduler.type}[/green]")
        from .experiment import interrupt_state

        super().__init__(
            num_epochs=config.training.num_epochs,
            start_epoch=self.start_epoch,
            validation_frequency=config.training.validation_frequency,
            early_stopping_patience=config.training.early_stopping_patience,
            save_frequency=config.checkpointing.save_frequency,
            save_during_training=config.checkpointing.save_during_training,
            scheduler=self.scheduler,
            optimizer=self.optimizer,
            wandb_run=wandb_run,
            history_path=self.models_dir / "training_history.json",
            console=console,
            interrupt_check=lambda: bool(interrupt_state.interrupted),
            no_validation_uses_train_metrics=True,
            log_learning_rate=True,
        )
        if config.checkpointing.load_checkpoint:
            self._load_checkpoint(config.checkpointing.load_checkpoint)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _run_epoch(self, loader: DataLoader, train: bool) -> Dict[str, float]:
        self.model.train(train)
        total_loss = total_correct = total_samples = 0
        phase = "train" if train else "val"
        data_wait_s = 0.0
        compute_s = 0.0

        with Progress(
            SpinnerColumn(),
            TextColumn(f"[{phase}]"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("loss={task.fields[loss]:.4f}"),
            TextColumn("acc={task.fields[acc]:.4f}"),
            TextColumn("samples={task.fields[samples]}"),
            TimeElapsedColumn(),
            console=console,
            transient=False,
        ) as progress:
            task = progress.add_task(
                f"{phase}",
                total=len(loader),
                loss=0.0,
                acc=0.0,
                samples=0,
            )

            with torch.set_grad_enabled(train):
                data_t0 = time.perf_counter()
                for batch in loader:
                    data_wait_s += time.perf_counter() - data_t0
                    compute_t0 = time.perf_counter()
                    if len(batch) == 3:
                        features, targets, _idx = batch
                    else:
                        features, targets = batch
                    features = move_to_device(features, self.device)
                    targets = move_to_device(targets, self.device)

                    outputs = self.model(features)

                    if self.is_classification:
                        valid = targets >= 0
                        if valid.sum() == 0:
                            progress.advance(task)
                            compute_s += time.perf_counter() - compute_t0
                            data_t0 = time.perf_counter()
                            continue
                        valid_count = valid.sum().item()
                        loss = self.criterion(outputs[valid], targets[valid])
                        preds = outputs[valid].argmax(dim=1)
                        total_correct += (preds == targets[valid]).sum().item()
                        total_samples += valid_count
                        loss_weight = valid_count
                    else:
                        loss = self.criterion(outputs, targets.float())
                        loss_weight = features.size(0)
                        total_samples += loss_weight

                    if train:
                        self.optimizer.zero_grad(set_to_none=True)
                        loss.backward()
                        self.optimizer.step()

                    total_loss += loss.item() * loss_weight
                    avg_loss = total_loss / max(total_samples, 1)
                    accuracy = total_correct / max(total_samples, 1) if self.is_classification else 0.0
                    progress.update(task, advance=1, loss=avg_loss, acc=accuracy, samples=total_samples)
                    compute_s += time.perf_counter() - compute_t0
                    if self.max_samples_per_epoch is not None and total_samples >= self.max_samples_per_epoch:
                        break
                    data_t0 = time.perf_counter()

        avg_loss = total_loss / max(total_samples, 1)
        accuracy = total_correct / max(total_samples, 1) if self.is_classification else 0.0
        data_fraction = data_wait_s / max(data_wait_s + compute_s, 1e-12)
        console.print(
            f"[cyan][profile {phase}][/cyan] data_wait={data_wait_s:.2f}s "
            f"compute={compute_s:.2f}s data_fraction={data_fraction:.2%}"
        )
        ds_stats = _dataset_profile_stats(getattr(loader, "dataset", None))
        if isinstance(ds_stats, dict) and ds_stats:
            console.print(
                f"[cyan][profile {phase} dataset][/cyan] "
                f"getitem={ds_stats.get('getitem_calls', 0)} "
                f"item_hit={ds_stats.get('item_cache_hits', 0)} "
                f"item_miss={ds_stats.get('item_cache_misses', 0)} "
                f"shard_hit={ds_stats.get('shard_cache_hits', 0)} "
                f"shard_miss={ds_stats.get('shard_cache_misses', 0)} "
                f"torch_load={ds_stats.get('shard_load_calls', 0)} "
                f"load_s={ds_stats.get('shard_load_s', 0.0):.2f}"
            )
        return {
            "loss": avg_loss,
            "accuracy": accuracy,
            "samples": total_samples,
            "data_wait_s": data_wait_s,
            "compute_s": compute_s,
            "data_fraction": data_fraction,
        }

    def _resolve_checkpoint_path(self, checkpoint: str) -> Path:
        return resolve_checkpoint_path(self.models_dir, checkpoint)

    def _load_checkpoint(self, checkpoint: str) -> None:
        checkpoint_path = self._resolve_checkpoint_path(checkpoint)
        if not checkpoint_path.exists():
            raise FileNotFoundError(f"Checkpoint nao encontrado: {checkpoint_path}")

        state = load_torch_checkpoint(
            checkpoint_path,
            self.model,
            optimizer=self.optimizer,
            scheduler=self.scheduler,
            device=self.device,
        )

        self.start_epoch = int(state.get("epoch", 0)) + 1
        self.best_val_loss = float(state.get("best_val_loss", state.get("loss", self.best_val_loss)))
        self.best_val_accuracy = float(state.get("best_val_accuracy", state.get("accuracy", self.best_val_accuracy)))
        if "best_val_loss" not in state:
            best_loss_path = self.models_dir / "best_loss.pt"
            if best_loss_path.exists():
                best_loss_state = torch.load(best_loss_path, map_location="cpu")
                self.best_val_loss = float(best_loss_state.get("loss", self.best_val_loss))
        if "best_val_accuracy" not in state:
            best_accuracy_path = self.models_dir / "best_accuracy.pt"
            if best_accuracy_path.exists():
                best_accuracy_state = torch.load(best_accuracy_path, map_location="cpu")
                self.best_val_accuracy = float(best_accuracy_state.get("accuracy", self.best_val_accuracy))
        console.print(f"[green]✓ Checkpoint carregado:[/green] {checkpoint_path} (retomando em E{self.start_epoch:03d})")

    def _save_checkpoint(self, filename: str, epoch: int, metrics: Dict[str, float]):
        save_torch_checkpoint(
            self.models_dir / filename,
            self.model,
            self.optimizer,
            epoch,
            metrics,
            scheduler=self.scheduler,
            extra={"best_val_loss": self.best_val_loss, "best_val_accuracy": self.best_val_accuracy},
        )

    def save_checkpoint(self, name: str, epoch: int, metrics: Dict[str, float]) -> None:
        self._save_checkpoint(name, epoch, metrics)

    def run_train_epoch(self) -> Dict[str, float]:
        return self._run_epoch(self.train_loader, train=True)

    def run_val_epoch(self) -> Dict[str, float]:
        return self._run_epoch(self.val_loader, train=False)

    def has_validation_data(self) -> bool:
        return len(self.val_loader.dataset) > 0

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def train(self) -> Dict[str, Any]:
        """
        Executa o loop de treinamento completo.

        Returns
        -------
        Dict
            Histórico de métricas por época.
        """
        console.print(
            f"\n[bold cyan]🚀 Iniciando treinamento "
            f"({self.start_epoch}-{self.num_epochs} épocas | "
            f"lr={self.config.training.learning_rate} | "
            f"opt={self.config.training.optimizer})[/bold cyan]"
        )
        history = super().train()
        console.print(
            f"\n[bold green]✓ Treinamento concluído! "
            f"best_val_acc={self.best_val_accuracy:.4f} | "
            f"best_val_loss={self.best_val_loss:.4f}[/bold green]"
        )
        return history
