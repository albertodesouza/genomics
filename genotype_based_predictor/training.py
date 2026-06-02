# -*- coding: utf-8 -*-
"""
training.py — Trainer: loop de treinamento, validação, checkpoints e W&B.
"""

import time
import json
from pathlib import Path
from typing import Any, Dict, Optional

import torch
import torch.nn as nn
import torch.optim as optim
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from torch.utils.data import DataLoader

from genotype_based_predictor.config import PipelineConfig
from genomics_pipeline.checkpointing import load_checkpoint as load_torch_checkpoint
from genomics_pipeline.checkpointing import resolve_checkpoint_path, save_checkpoint as save_torch_checkpoint
from genomics_pipeline.optim import make_optimizer_from_config
from genomics_pipeline.torch_utils import move_to_device

console = Console()


def _make_optimizer(model: nn.Module, config: PipelineConfig) -> optim.Optimizer:
    return make_optimizer_from_config(model, config.training, foreach=False)


def _make_scheduler(optimizer: optim.Optimizer, config: PipelineConfig):
    sc = config.training.lr_scheduler
    if not sc.enabled:
        return None
    t = sc.type.lower()
    if t == "plateau":
        return optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode=sc.mode, factor=sc.factor,
            patience=sc.patience, min_lr=sc.min_lr,
        )
    elif t == "cosine":
        return optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=sc.T_max, eta_min=sc.eta_min,
        )
    elif t == "step":
        return optim.lr_scheduler.StepLR(optimizer, step_size=sc.step_size, gamma=sc.gamma)
    elif t == "exponential":
        return optim.lr_scheduler.ExponentialLR(optimizer, gamma=sc.gamma)
    elif t == "multistep":
        return optim.lr_scheduler.MultiStepLR(optimizer, milestones=sc.milestones, gamma=sc.gamma)
    elif t == "cosine_warm_restarts":
        return optim.lr_scheduler.CosineAnnealingWarmRestarts(
            optimizer, T_0=sc.T_0, T_mult=sc.T_mult, eta_min=sc.eta_min,
        )
    console.print(f"[yellow]Scheduler '{t}' desconhecido. Sem scheduler.[/yellow]")
    return None


class Trainer:
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
        self.num_epochs = config.training.num_epochs
        self.early_stop_patience = config.training.early_stopping_patience
        self.max_samples_per_epoch = config.debug.max_samples_per_epoch

        self.history: Dict[str, list] = {
            "train_loss": [], "train_accuracy": [],
            "val_loss": [], "val_accuracy": [], "epoch": [],
        }
        self.best_val_loss = float("inf")
        self.best_val_accuracy = 0.0
        self.start_epoch = 1

        models_dir = experiment_dir / "models"
        models_dir.mkdir(parents=True, exist_ok=True)
        self.models_dir = models_dir

        if config.checkpointing.load_checkpoint:
            self._load_checkpoint(config.checkpointing.load_checkpoint)

        if config.training.weight_decay > 0:
            console.print(f"[green]✓ Weight decay (L2): {config.training.weight_decay}[/green]")
        if config.training.lr_scheduler.enabled:
            console.print(f"[green]✓ LR Scheduler: {config.training.lr_scheduler.type}[/green]")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _run_epoch(self, loader: DataLoader, train: bool) -> Dict[str, float]:
        self.model.train(train)
        total_loss = total_correct = total_samples = 0
        phase = "train" if train else "val"

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
                for batch in loader:
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
                            continue
                        loss = self.criterion(outputs[valid], targets[valid])
                        preds = outputs[valid].argmax(dim=1)
                        total_correct += (preds == targets[valid]).sum().item()
                        total_samples += valid.sum().item()
                    else:
                        loss = self.criterion(outputs, targets.float())
                        total_samples += features.size(0)

                    if train:
                        self.optimizer.zero_grad()
                        loss.backward()
                        self.optimizer.step()

                    total_loss += loss.item() * features.size(0)
                    avg_loss = total_loss / max(total_samples, 1)
                    accuracy = total_correct / max(total_samples, 1) if self.is_classification else 0.0
                    progress.update(task, advance=1, loss=avg_loss, acc=accuracy, samples=total_samples)
                    if self.max_samples_per_epoch is not None and total_samples >= self.max_samples_per_epoch:
                        break

        avg_loss = total_loss / max(total_samples, 1)
        accuracy = total_correct / max(total_samples, 1) if self.is_classification else 0.0
        return {"loss": avg_loss, "accuracy": accuracy, "samples": total_samples}

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

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def train(self) -> Dict[str, list]:
        """
        Executa o loop de treinamento completo.

        Returns
        -------
        Dict
            Histórico de métricas por época.
        """
        from genotype_based_predictor.experiment import interrupt_state

        console.print(
            f"\n[bold cyan]🚀 Iniciando treinamento "
            f"({self.start_epoch}-{self.num_epochs} épocas | "
            f"lr={self.config.training.learning_rate} | "
            f"opt={self.config.training.optimizer})[/bold cyan]"
        )

        val_metrics: Dict[str, float] = {"loss": float("inf"), "accuracy": 0.0}
        no_improve = 0
        epoch = 0
        val_frequency = max(1, int(self.config.training.validation_frequency))

        for epoch in range(self.start_epoch, self.num_epochs + 1):
            if interrupt_state.interrupted:
                console.print("[yellow]⚠ Treinamento interrompido (CTRL+C)[/yellow]")
                break

            train_metrics = self._run_epoch(self.train_loader, train=True)
            should_validate = (
                len(self.val_loader.dataset) > 0
                and (epoch % val_frequency == 0 or epoch == self.num_epochs)
            )
            if should_validate:
                val_metrics = self._run_epoch(self.val_loader, train=False)
            elif len(self.val_loader.dataset) == 0:
                val_metrics = {"loss": train_metrics["loss"], "accuracy": train_metrics["accuracy"], "samples": 0}

            # Scheduler step
            if self.scheduler is not None:
                if isinstance(self.scheduler, optim.lr_scheduler.ReduceLROnPlateau):
                    if should_validate or len(self.val_loader.dataset) == 0:
                        self.scheduler.step(val_metrics["loss"])
                else:
                    self.scheduler.step()

            # Histórico
            self.history["epoch"].append(epoch)
            self.history["train_loss"].append(train_metrics["loss"])
            self.history["train_accuracy"].append(train_metrics["accuracy"])
            self.history["val_loss"].append(val_metrics["loss"])
            self.history["val_accuracy"].append(val_metrics["accuracy"])

            can_update_best = should_validate or len(self.val_loader.dataset) == 0
            improved_acc = can_update_best and val_metrics["accuracy"] > self.best_val_accuracy
            improved_loss = can_update_best and val_metrics["loss"] < self.best_val_loss

            if improved_acc:
                self.best_val_accuracy = val_metrics["accuracy"]
                if self.config.checkpointing.save_during_training:
                    self._save_checkpoint("best_accuracy.pt", epoch, val_metrics)

            if improved_loss:
                self.best_val_loss = val_metrics["loss"]
                if self.config.checkpointing.save_during_training:
                    self._save_checkpoint("best_loss.pt", epoch, val_metrics)
                no_improve = 0
            elif can_update_best:
                no_improve += 1

            if (
                self.config.checkpointing.save_during_training
                and self.config.checkpointing.save_frequency > 0
                and epoch % self.config.checkpointing.save_frequency == 0
            ):
                self._save_checkpoint(f"epoch_{epoch}.pt", epoch, val_metrics)

            lr_now = self.optimizer.param_groups[0]["lr"]
            acc_tag = " [green]✓best_acc[/green]" if improved_acc else ""
            loss_tag = " [cyan]✓best_loss[/cyan]" if improved_loss else ""
            val_text = (
                f"val loss={val_metrics['loss']:.4f} acc={val_metrics['accuracy']:.4f}"
                if can_update_best else f"val skipped (freq={val_frequency})"
            )
            console.print(
                f"[E{epoch:03d}] "
                f"train loss={train_metrics['loss']:.4f} acc={train_metrics['accuracy']:.4f} | "
                f"{val_text} | "
                f"lr={lr_now:.2e}{acc_tag}{loss_tag}"
            )

            if self.wandb_run:
                log_payload = {
                    "epoch": epoch,
                    "train/loss": train_metrics["loss"],
                    "train/accuracy": train_metrics["accuracy"],
                    "lr": lr_now,
                }
                if can_update_best:
                    log_payload.update({
                        "val/loss": val_metrics["loss"],
                        "val/accuracy": val_metrics["accuracy"],
                    })
                self.wandb_run.log(log_payload)

            if self.early_stop_patience and no_improve >= self.early_stop_patience:
                console.print(f"[yellow]Early stopping após {no_improve} épocas sem melhora[/yellow]")
                break

        self._save_checkpoint("final.pt", epoch, val_metrics)
        console.print(
            f"\n[bold green]✓ Treinamento concluído! "
            f"best_val_acc={self.best_val_accuracy:.4f} | "
            f"best_val_loss={self.best_val_loss:.4f}[/bold green]"
        )
        self.history["interrupted"] = interrupt_state.interrupted
        self.history["last_epoch"] = epoch
        with open(self.models_dir / "training_history.json", "w", encoding="utf-8") as f:
            json.dump(self.history, f, indent=2)
        return self.history
