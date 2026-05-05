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

from genotype_based_predictor.config import PipelineConfig

console = Console()


def _make_optimizer(model: nn.Module, config: PipelineConfig) -> optim.Optimizer:
    t = config.training.optimizer.lower()
    lr = config.training.learning_rate
    wd = config.training.weight_decay
    if t == "adam":
        return optim.Adam(model.parameters(), lr=lr, weight_decay=wd, foreach=False)
    elif t == "adamw":
        return optim.AdamW(model.parameters(), lr=lr, weight_decay=wd, foreach=False)
    elif t == "sgd":
        return optim.SGD(model.parameters(), lr=lr, momentum=0.9, weight_decay=wd, foreach=False)
    raise ValueError(f"Otimizador não suportado: {t}")


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

        self.history: Dict[str, list] = {
            "train_loss": [], "train_accuracy": [],
            "val_loss": [], "val_accuracy": [], "epoch": [],
        }
        self.best_val_loss = float("inf")
        self.best_val_accuracy = 0.0

        models_dir = experiment_dir / "models"
        models_dir.mkdir(parents=True, exist_ok=True)
        self.models_dir = models_dir

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
                    features = features.to(self.device, non_blocking=True)
                    targets = targets.to(self.device, non_blocking=True)

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

        avg_loss = total_loss / max(total_samples, 1)
        accuracy = total_correct / max(total_samples, 1) if self.is_classification else 0.0
        return {"loss": avg_loss, "accuracy": accuracy, "samples": total_samples}

    def _save_checkpoint(self, filename: str, epoch: int, metrics: Dict[str, float]):
        torch.save({
            "epoch": epoch,
            "model_state_dict": self.model.state_dict(),
            "optimizer_state_dict": self.optimizer.state_dict(),
            **metrics,
        }, self.models_dir / filename)

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
            f"({self.num_epochs} épocas | "
            f"lr={self.config.training.learning_rate} | "
            f"opt={self.config.training.optimizer})[/bold cyan]"
        )

        val_metrics: Dict[str, float] = {"loss": float("inf"), "accuracy": 0.0}
        no_improve = 0
        epoch = 0

        for epoch in range(1, self.num_epochs + 1):
            if interrupt_state.interrupted:
                console.print("[yellow]⚠ Treinamento interrompido (CTRL+C)[/yellow]")
                break

            train_metrics = self._run_epoch(self.train_loader, train=True)
            val_metrics = self._run_epoch(self.val_loader, train=False)

            # Scheduler step
            if self.scheduler is not None:
                if isinstance(self.scheduler, optim.lr_scheduler.ReduceLROnPlateau):
                    self.scheduler.step(val_metrics["loss"])
                else:
                    self.scheduler.step()

            # Histórico
            self.history["epoch"].append(epoch)
            self.history["train_loss"].append(train_metrics["loss"])
            self.history["train_accuracy"].append(train_metrics["accuracy"])
            self.history["val_loss"].append(val_metrics["loss"])
            self.history["val_accuracy"].append(val_metrics["accuracy"])

            improved_acc = val_metrics["accuracy"] > self.best_val_accuracy
            improved_loss = val_metrics["loss"] < self.best_val_loss

            if improved_acc:
                self.best_val_accuracy = val_metrics["accuracy"]
                self._save_checkpoint("best_accuracy.pt", epoch, val_metrics)

            if improved_loss:
                self.best_val_loss = val_metrics["loss"]
                self._save_checkpoint("best_loss.pt", epoch, val_metrics)
                no_improve = 0
            else:
                no_improve += 1

            lr_now = self.optimizer.param_groups[0]["lr"]
            acc_tag = " [green]✓best_acc[/green]" if improved_acc else ""
            loss_tag = " [cyan]✓best_loss[/cyan]" if improved_loss else ""
            console.print(
                f"[E{epoch:03d}] "
                f"train loss={train_metrics['loss']:.4f} acc={train_metrics['accuracy']:.4f} | "
                f"val loss={val_metrics['loss']:.4f} acc={val_metrics['accuracy']:.4f} | "
                f"lr={lr_now:.2e}{acc_tag}{loss_tag}"
            )

            if self.wandb_run:
                self.wandb_run.log({
                    "epoch": epoch,
                    "train/loss": train_metrics["loss"],
                    "train/accuracy": train_metrics["accuracy"],
                    "val/loss": val_metrics["loss"],
                    "val/accuracy": val_metrics["accuracy"],
                    "lr": lr_now,
                })

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
        return self.history
