from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Dict, MutableMapping, Optional

import torch.optim as optim
from rich.console import Console


def _get_config_value(config: Any, name: str, default: Any = None) -> Any:
    if isinstance(config, dict):
        return config.get(name, default)
    return getattr(config, name, default)


def make_lr_scheduler(optimizer: optim.Optimizer, scheduler_config: Any, *, default_t_max: Optional[int] = None):
    if not _get_config_value(scheduler_config, "enabled", False):
        return None
    scheduler_type = str(_get_config_value(scheduler_config, "type", "plateau")).lower()
    if scheduler_type == "plateau":
        return optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode=_get_config_value(scheduler_config, "mode", "min"),
            factor=_get_config_value(scheduler_config, "factor", 0.5),
            patience=_get_config_value(scheduler_config, "patience", 10),
            min_lr=_get_config_value(scheduler_config, "min_lr", 1e-6),
        )
    if scheduler_type == "cosine":
        return optim.lr_scheduler.CosineAnnealingLR(
            optimizer,
            T_max=_get_config_value(scheduler_config, "T_max", default_t_max or 100),
            eta_min=_get_config_value(scheduler_config, "eta_min", 1e-6),
        )
    if scheduler_type == "step":
        return optim.lr_scheduler.StepLR(
            optimizer,
            step_size=_get_config_value(scheduler_config, "step_size", 30),
            gamma=_get_config_value(scheduler_config, "gamma", 0.1),
        )
    if scheduler_type == "exponential":
        return optim.lr_scheduler.ExponentialLR(
            optimizer,
            gamma=_get_config_value(scheduler_config, "gamma", 0.95),
        )
    if scheduler_type == "multistep":
        return optim.lr_scheduler.MultiStepLR(
            optimizer,
            milestones=_get_config_value(scheduler_config, "milestones", [30, 60, 90]),
            gamma=_get_config_value(scheduler_config, "gamma", 0.1),
        )
    if scheduler_type == "cosine_warm_restarts":
        return optim.lr_scheduler.CosineAnnealingWarmRestarts(
            optimizer,
            T_0=_get_config_value(scheduler_config, "T_0", 50),
            T_mult=_get_config_value(scheduler_config, "T_mult", 1),
            eta_min=_get_config_value(scheduler_config, "eta_min", 1e-6),
        )
    raise ValueError(f"Scheduler nao suportado: {scheduler_type}")


def step_lr_scheduler(scheduler: Any, val_loss: float, *, metric_available: bool = True) -> None:
    if scheduler is None:
        return
    if isinstance(scheduler, optim.lr_scheduler.ReduceLROnPlateau):
        if metric_available:
            scheduler.step(val_loss)
    else:
        scheduler.step()


def new_training_history() -> Dict[str, list]:
    return {"epoch": [], "train_loss": [], "train_accuracy": [], "val_loss": [], "val_accuracy": []}


def append_history_epoch(
    history: MutableMapping[str, list],
    epoch: int,
    train_metrics: Dict[str, float],
    val_metrics: Dict[str, float],
) -> None:
    history.setdefault("epoch", []).append(epoch)
    history.setdefault("train_loss", []).append(train_metrics["loss"])
    history.setdefault("train_accuracy", []).append(train_metrics["accuracy"])
    history.setdefault("val_loss", []).append(val_metrics["loss"])
    history.setdefault("val_accuracy", []).append(val_metrics["accuracy"])
    for key in ("data_wait_s", "compute_s", "data_fraction"):
        if key in train_metrics:
            history.setdefault(f"train_{key}", []).append(train_metrics[key])
        if key in val_metrics:
            history.setdefault(f"val_{key}", []).append(val_metrics[key])


def write_training_history(path: Path, history: Dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(history, f, indent=2)


class EpochTrainer:
    """Generic epoch loop for pipeline-specific trainers.

    Subclasses provide data/model-specific epoch implementations and checkpoint
    persistence. This class owns shared bookkeeping: validation cadence, best
    metrics, periodic/final checkpoints, history and W&B logging.
    """

    def __init__(
        self,
        *,
        num_epochs: int,
        start_epoch: int = 1,
        validation_frequency: int = 1,
        early_stopping_patience: Optional[int] = None,
        save_frequency: int = 0,
        save_during_training: bool = True,
        scheduler: Any = None,
        optimizer: Any = None,
        wandb_run: Any = None,
        history_path: Optional[Path] = None,
        console: Optional[Console] = None,
        interrupt_check: Optional[Callable[[], bool]] = None,
        no_validation_uses_train_metrics: bool = True,
        log_learning_rate: bool = True,
    ):
        self.num_epochs = int(num_epochs)
        self.start_epoch = int(start_epoch)
        self.validation_frequency = max(1, int(validation_frequency))
        self.early_stopping_patience = early_stopping_patience
        self.save_frequency = int(save_frequency or 0)
        self.save_during_training = bool(save_during_training)
        self.scheduler = scheduler
        self.optimizer = optimizer
        self.wandb_run = wandb_run
        self.history_path = Path(history_path) if history_path else None
        self.console = console or Console()
        self.interrupt_check = interrupt_check or (lambda: False)
        self.no_validation_uses_train_metrics = no_validation_uses_train_metrics
        self.log_learning_rate = log_learning_rate
        self.history: Dict[str, Any] = new_training_history()
        self.best_val_loss = float("inf")
        self.best_val_accuracy = 0.0

    def run_train_epoch(self) -> Dict[str, float]:
        raise NotImplementedError

    def run_val_epoch(self) -> Dict[str, float]:
        raise NotImplementedError

    def has_validation_data(self) -> bool:
        return True

    def save_checkpoint(self, name: str, epoch: int, metrics: Dict[str, float]) -> None:
        raise NotImplementedError

    def _current_lr(self) -> Optional[float]:
        if self.optimizer is None or not getattr(self.optimizer, "param_groups", None):
            return None
        return float(self.optimizer.param_groups[0]["lr"])

    def _log_epoch(self, epoch: int, train_metrics: Dict[str, float], val_metrics: Dict[str, float], *, can_update_best: bool, improved_acc: bool, improved_loss: bool) -> None:
        lr_now = self._current_lr()
        acc_tag = " [green]✓best_acc[/green]" if improved_acc else ""
        loss_tag = " [cyan]✓best_loss[/cyan]" if improved_loss else ""
        val_text = (
            f"val loss={val_metrics['loss']:.4f} acc={val_metrics['accuracy']:.4f}"
            if can_update_best else f"val skipped (freq={self.validation_frequency})"
        )
        lr_text = f" | lr={lr_now:.2e}" if self.log_learning_rate and lr_now is not None else ""
        self.console.print(
            f"[E{epoch:03d}] "
            f"train loss={train_metrics['loss']:.4f} acc={train_metrics['accuracy']:.4f} | "
            f"{val_text}{lr_text}{acc_tag}{loss_tag}"
        )
        if self.wandb_run:
            payload = {
                "epoch": epoch,
                "train/loss": train_metrics["loss"],
                "train/accuracy": train_metrics["accuracy"],
            }
            if lr_now is not None:
                payload["lr"] = lr_now
            for key in ("data_wait_s", "compute_s", "data_fraction"):
                if key in train_metrics:
                    payload[f"train/{key}"] = train_metrics[key]
            if can_update_best:
                payload.update({"val/loss": val_metrics["loss"], "val/accuracy": val_metrics["accuracy"]})
                for key in ("data_wait_s", "compute_s", "data_fraction"):
                    if key in val_metrics:
                        payload[f"val/{key}"] = val_metrics[key]
            self.wandb_run.log(payload)

    def train(self) -> Dict[str, Any]:
        val_metrics: Dict[str, float] = {"loss": float("inf"), "accuracy": 0.0, "samples": 0}
        no_improve = 0
        epoch = self.start_epoch - 1
        interrupted = False
        for epoch in range(self.start_epoch, self.num_epochs + 1):
            if self.interrupt_check():
                interrupted = True
                self.console.print("[yellow]Treinamento interrompido[/yellow]")
                break

            train_metrics = self.run_train_epoch()
            has_val = self.has_validation_data()
            should_validate = has_val and (epoch % self.validation_frequency == 0 or epoch == self.num_epochs)
            if should_validate:
                val_metrics = self.run_val_epoch()
            elif not has_val and self.no_validation_uses_train_metrics:
                val_metrics = dict(train_metrics)
                val_metrics["samples"] = 0

            can_update_best = should_validate or (not has_val and self.no_validation_uses_train_metrics)
            step_lr_scheduler(self.scheduler, val_metrics["loss"], metric_available=can_update_best)
            append_history_epoch(self.history, epoch, train_metrics, val_metrics)

            improved_acc = can_update_best and val_metrics["accuracy"] > self.best_val_accuracy
            improved_loss = can_update_best and val_metrics["loss"] < self.best_val_loss
            if improved_acc:
                self.best_val_accuracy = val_metrics["accuracy"]
                if self.save_during_training:
                    self.save_checkpoint("best_accuracy.pt", epoch, val_metrics)
            if improved_loss:
                self.best_val_loss = val_metrics["loss"]
                no_improve = 0
                if self.save_during_training:
                    self.save_checkpoint("best_loss.pt", epoch, val_metrics)
            elif can_update_best:
                no_improve += 1

            if self.save_during_training and self.save_frequency > 0 and epoch % self.save_frequency == 0:
                self.save_checkpoint(f"epoch_{epoch}.pt", epoch, val_metrics)

            self._log_epoch(epoch, train_metrics, val_metrics, can_update_best=can_update_best, improved_acc=improved_acc, improved_loss=improved_loss)

            if self.early_stopping_patience and no_improve >= self.early_stopping_patience:
                self.console.print(f"[yellow]Early stopping apos {no_improve} epocas sem melhora[/yellow]")
                break

        self.save_checkpoint("final.pt", epoch, val_metrics)
        self.history["interrupted"] = interrupted
        self.history["last_epoch"] = epoch
        if self.history_path:
            write_training_history(self.history_path, self.history)
        return self.history
