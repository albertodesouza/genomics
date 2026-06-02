from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import torch
import torch.nn as nn
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from genomics_pipeline.checkpointing import load_checkpoint as load_torch_checkpoint
from genomics_pipeline.checkpointing import save_checkpoint as save_torch_checkpoint
from genomics_pipeline.optim import make_optimizer_from_config
from genomics_pipeline.training_utils import append_history_epoch, new_training_history
from variant_transformer_predictor.config import PipelineConfig
from variant_transformer_predictor.evaluation import move_batch_to_device

console = Console()


def make_optimizer(model: nn.Module, config: PipelineConfig):
    return make_optimizer_from_config(model, config.training)


class Trainer:
    def __init__(self, model, train_loader, val_loader, config: PipelineConfig, device: torch.device, experiment_dir: Path, class_names, wandb_run=None):
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.config = config
        self.device = device
        self.experiment_dir = Path(experiment_dir)
        self.class_names = class_names
        self.wandb_run = wandb_run
        self.optimizer = make_optimizer(model, config)
        self.criterion = nn.CrossEntropyLoss()
        self.models_dir = self.experiment_dir / "models"
        self.models_dir.mkdir(parents=True, exist_ok=True)
        self.best_val_accuracy = 0.0
        self.best_val_loss = float("inf")
        self.start_epoch = 1
        if config.checkpointing.load_checkpoint:
            self.load_checkpoint(config.checkpointing.load_checkpoint)

    def load_checkpoint(self, checkpoint: str) -> None:
        path = Path(checkpoint)
        if not path.exists():
            path = self.models_dir / (checkpoint if checkpoint.endswith(".pt") else f"{checkpoint}.pt")
        state = load_torch_checkpoint(path, self.model, optimizer=self.optimizer, device=self.device, weights_only=True)
        self.start_epoch = int(state.get("epoch", 0)) + 1
        self.best_val_accuracy = float(state.get("best_val_accuracy", 0.0))
        self.best_val_loss = float(state.get("best_val_loss", float("inf")))
        console.print(f"[green]Checkpoint carregado:[/green] {path}")

    def save_checkpoint(self, name: str, epoch: int, metrics: Dict[str, float]) -> None:
        save_torch_checkpoint(
            self.models_dir / name,
            self.model,
            self.optimizer,
            epoch,
            metrics,
            extra={"best_val_accuracy": self.best_val_accuracy, "best_val_loss": self.best_val_loss},
        )

    def run_train_epoch(self) -> Dict[str, float]:
        self.model.train()
        total_loss = 0.0
        total_correct = 0
        total = 0
        with Progress(SpinnerColumn(), TextColumn("train"), BarColumn(), TextColumn("{task.completed}/{task.total}"), TimeElapsedColumn(), console=console, transient=True) as progress:
            task = progress.add_task("train", total=len(self.train_loader))
            for batch in self.train_loader:
                batch = move_batch_to_device(batch, self.device)
                logits = self.model(batch)
                loss = self.criterion(logits, batch["targets"])
                self.optimizer.zero_grad(set_to_none=True)
                loss.backward()
                self.optimizer.step()
                bs = int(batch["targets"].shape[0])
                total += bs
                total_loss += float(loss.item()) * bs
                total_correct += int((logits.argmax(dim=1) == batch["targets"]).sum().item())
                progress.advance(task)
        return {"loss": total_loss / max(total, 1), "accuracy": total_correct / max(total, 1), "samples": total}

    @torch.no_grad()
    def run_val_loss(self) -> Dict[str, float]:
        self.model.eval()
        total_loss = 0.0
        total_correct = 0
        total = 0
        for batch in self.val_loader:
            batch = move_batch_to_device(batch, self.device)
            logits = self.model(batch)
            loss = self.criterion(logits, batch["targets"])
            bs = int(batch["targets"].shape[0])
            total += bs
            total_loss += float(loss.item()) * bs
            total_correct += int((logits.argmax(dim=1) == batch["targets"]).sum().item())
        return {"loss": total_loss / max(total, 1), "accuracy": total_correct / max(total, 1), "samples": total}

    def train(self) -> Dict:
        history = new_training_history()
        no_improve = 0
        for epoch in range(self.start_epoch, self.config.training.num_epochs + 1):
            train_metrics = self.run_train_epoch()
            val_metrics = self.run_val_loss() if len(self.val_loader.dataset) else train_metrics
            append_history_epoch(history, epoch, train_metrics, val_metrics)
            improved_acc = val_metrics["accuracy"] > self.best_val_accuracy
            improved_loss = val_metrics["loss"] < self.best_val_loss
            if improved_acc:
                self.best_val_accuracy = val_metrics["accuracy"]
                if self.config.checkpointing.save_during_training:
                    self.save_checkpoint("best_accuracy.pt", epoch, val_metrics)
            if improved_loss:
                self.best_val_loss = val_metrics["loss"]
                no_improve = 0
                if self.config.checkpointing.save_during_training:
                    self.save_checkpoint("best_loss.pt", epoch, val_metrics)
            else:
                no_improve += 1
            if self.config.checkpointing.save_frequency and epoch % self.config.checkpointing.save_frequency == 0:
                self.save_checkpoint(f"epoch_{epoch}.pt", epoch, val_metrics)
            console.print(
                f"[E{epoch:03d}] train loss={train_metrics['loss']:.4f} acc={train_metrics['accuracy']:.4f} | "
                f"val loss={val_metrics['loss']:.4f} acc={val_metrics['accuracy']:.4f}"
            )
            if self.wandb_run:
                self.wandb_run.log({"epoch": epoch, "train/loss": train_metrics["loss"], "train/accuracy": train_metrics["accuracy"], "val/loss": val_metrics["loss"], "val/accuracy": val_metrics["accuracy"]})
            patience = self.config.training.early_stopping_patience
            if patience and no_improve >= patience:
                console.print(f"[yellow]Early stopping apos {no_improve} epocas sem melhora[/yellow]")
                break
        self.save_checkpoint("final.pt", epoch, val_metrics)
        return history
