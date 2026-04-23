from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import torch
from torch import nn
from sklearn.metrics import classification_report, confusion_matrix, f1_score

from .config import PredictionConfig


@dataclass
class EpochMetrics:
    loss: float
    accuracy: float
    f1: float
    macro_f1: float
    classification_report: str = ""
    confusion_matrix: List[List[int]] | None = None
    mse: float = 0.0


def evaluate(model: nn.Module, loader, device: torch.device, criterion: nn.Module, is_regression: bool = False) -> EpochMetrics:
    model.eval()
    total_loss = 0.0
    total_correct = 0
    total_examples = 0
    all_targets: List[int] = []
    all_predictions: List[int] = []
    total_squared_error = 0.0
    with torch.no_grad():
        for features, targets in loader:
            features = features.to(device)
            targets = targets.to(device)
            logits = model(features)
            loss = criterion(logits, targets)
            total_loss += loss.item() * targets.size(0)
            if is_regression:
                predictions = logits
                total_squared_error += torch.sum((predictions - targets) ** 2).item()
            else:
                predictions = logits.argmax(dim=1)
                total_correct += (predictions == targets).sum().item()
            total_examples += targets.size(0)
            if not is_regression:
                all_targets.extend(targets.cpu().tolist())
                all_predictions.extend(predictions.cpu().tolist())

    if is_regression:
        return EpochMetrics(
            loss=total_loss / max(total_examples, 1),
            accuracy=0.0,
            f1=0.0,
            macro_f1=0.0,
            mse=total_squared_error / max(total_examples, 1),
        )

    if all_targets:
        f1 = float(f1_score(all_targets, all_predictions, average="weighted", zero_division=0))
        macro_f1 = float(f1_score(all_targets, all_predictions, average="macro", zero_division=0))
        report = classification_report(all_targets, all_predictions, zero_division=0)
        matrix = confusion_matrix(all_targets, all_predictions).tolist()
    else:
        f1 = 0.0
        macro_f1 = 0.0
        report = ""
        matrix = []

    return EpochMetrics(
        loss=total_loss / max(total_examples, 1),
        accuracy=total_correct / max(total_examples, 1),
        f1=f1,
        macro_f1=macro_f1,
        classification_report=report,
        confusion_matrix=matrix,
    )


def evaluate_model(model: nn.Module, data_bundle, device: torch.device) -> Dict:
    is_regression = data_bundle.dataset.prediction_target == "frog_likelihood"
    criterion = nn.MSELoss() if is_regression else nn.CrossEntropyLoss()
    val_metrics = evaluate(model, data_bundle.val_loader, device, criterion, is_regression=is_regression)
    test_metrics = evaluate(model, data_bundle.test_loader, device, criterion, is_regression=is_regression)
    return {
        "val_loss": val_metrics.loss,
        "val_accuracy": val_metrics.accuracy,
        "val_f1": val_metrics.f1,
        "val_macro_f1": val_metrics.macro_f1,
        "val_classification_report": val_metrics.classification_report,
        "val_confusion_matrix": val_metrics.confusion_matrix,
        "val_mse": val_metrics.mse,
        "test_loss": test_metrics.loss,
        "test_accuracy": test_metrics.accuracy,
        "test_f1": test_metrics.f1,
        "test_macro_f1": test_metrics.macro_f1,
        "test_classification_report": test_metrics.classification_report,
        "test_confusion_matrix": test_metrics.confusion_matrix,
        "test_mse": test_metrics.mse,
    }


def load_checkpoint_if_available(model: nn.Module, checkpoint_path: str | None, device: torch.device) -> None:
    if not checkpoint_path:
        return
    state = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(state["model_state_dict"])


def save_checkpoint(model: nn.Module, output_dir: Path, name: str, epoch: int, metrics: Dict) -> Path:
    checkpoint_path = output_dir / f"{name}.pt"
    torch.save(
        {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "metrics": metrics,
        },
        checkpoint_path,
    )
    return checkpoint_path


def train_model(model: nn.Module, data_bundle, config: PredictionConfig, device: torch.device, output_dir: Path, wandb_run=None) -> Dict:
    is_regression = data_bundle.dataset.prediction_target == "frog_likelihood"
    criterion = nn.MSELoss() if is_regression else nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=float(config.training.learning_rate),
        weight_decay=float(config.training.weight_decay),
    )
    load_checkpoint_if_available(model, config.checkpointing.load_checkpoint, device)

    history = {
        "train_loss": [],
        "train_accuracy": [],
        "train_f1": [],
        "train_macro_f1": [],
        "val_loss": [],
        "val_accuracy": [],
        "val_f1": [],
        "val_macro_f1": [],
        "val_mse": [],
    }
    best_val_accuracy = float("-inf")
    best_val_loss = float("inf")
    try:
        for epoch in range(int(config.training.num_epochs)):
            model.train()
            total_loss = 0.0
            total_correct = 0
            total_examples = 0
            train_targets: List[int] = []
            train_predictions: List[int] = []
            train_squared_error = 0.0
            for features, targets in data_bundle.train_loader:
                features = features.to(device)
                targets = targets.to(device)
                optimizer.zero_grad()
                logits = model(features)
                loss = criterion(logits, targets)
                loss.backward()
                optimizer.step()
                total_loss += loss.item() * targets.size(0)
                if is_regression:
                    predictions = logits
                    train_squared_error += torch.sum((predictions - targets) ** 2).item()
                else:
                    predictions = logits.argmax(dim=1)
                    total_correct += (predictions == targets).sum().item()
                total_examples += targets.size(0)
                if not is_regression:
                    train_targets.extend(targets.cpu().tolist())
                    train_predictions.extend(predictions.cpu().tolist())

            train_metrics = EpochMetrics(
                loss=total_loss / max(total_examples, 1),
                accuracy=0.0 if is_regression else total_correct / max(total_examples, 1),
                f1=float(f1_score(train_targets, train_predictions, average="weighted", zero_division=0)) if train_targets else 0.0,
                macro_f1=float(f1_score(train_targets, train_predictions, average="macro", zero_division=0)) if train_targets else 0.0,
                mse=train_squared_error / max(total_examples, 1) if is_regression else 0.0,
            )
            val_metrics = evaluate(model, data_bundle.val_loader, device, criterion, is_regression=is_regression)
            history["train_loss"].append(train_metrics.loss)
            history["train_accuracy"].append(train_metrics.accuracy)
            history["train_f1"].append(train_metrics.f1)
            history["train_macro_f1"].append(train_metrics.macro_f1)
            history["val_loss"].append(val_metrics.loss)
            history["val_accuracy"].append(val_metrics.accuracy)
            history["val_f1"].append(val_metrics.f1)
            history["val_macro_f1"].append(val_metrics.macro_f1)
            history["val_mse"].append(val_metrics.mse)

            if wandb_run is not None:
                wandb_run.log(
                    {
                        "epoch": epoch + 1,
                        "train/loss": train_metrics.loss,
                        "train/accuracy": train_metrics.accuracy,
                        "train/f1": train_metrics.f1,
                        "train/macro_f1": train_metrics.macro_f1,
                        "val/loss": val_metrics.loss,
                        "val/accuracy": val_metrics.accuracy,
                        "val/f1": val_metrics.f1,
                        "val/macro_f1": val_metrics.macro_f1,
                        "val/mse": val_metrics.mse,
                    }
                )

            if config.checkpointing.save_best and not is_regression and val_metrics.accuracy > best_val_accuracy:
                best_val_accuracy = val_metrics.accuracy
                save_checkpoint(
                    model,
                    output_dir,
                    "best_accuracy",
                    epoch,
                    {"val_accuracy": val_metrics.accuracy, "val_f1": val_metrics.f1, "val_macro_f1": val_metrics.macro_f1},
                )

            if config.checkpointing.save_best_loss and val_metrics.loss < best_val_loss:
                best_val_loss = val_metrics.loss
                save_checkpoint(
                    model,
                    output_dir,
                    "best_loss",
                    epoch,
                    {"val_loss": val_metrics.loss, "val_mse": val_metrics.mse},
                )

            if config.checkpointing.save_frequency > 0 and (epoch + 1) % config.checkpointing.save_frequency == 0:
                save_checkpoint(
                    model,
                    output_dir,
                    f"epoch_{epoch + 1}",
                    epoch,
                    {"val_accuracy": val_metrics.accuracy},
                )
    except KeyboardInterrupt:
        save_checkpoint(model, output_dir, "interrupted", max(len(history["train_loss"]) - 1, 0), {"status": "interrupted"})
        raise

    test_metrics = evaluate(model, data_bundle.test_loader, device, criterion, is_regression=is_regression)
    history["test_loss"] = test_metrics.loss
    history["test_accuracy"] = test_metrics.accuracy
    history["test_f1"] = test_metrics.f1
    history["test_macro_f1"] = test_metrics.macro_f1
    history["test_classification_report"] = test_metrics.classification_report
    history["test_confusion_matrix"] = test_metrics.confusion_matrix
    history["test_mse"] = test_metrics.mse

    if wandb_run is not None:
        wandb_run.log(
            {
                "test/loss": test_metrics.loss,
                "test/accuracy": test_metrics.accuracy,
                "test/f1": test_metrics.f1,
                "test/macro_f1": test_metrics.macro_f1,
                "test/mse": test_metrics.mse,
            }
        )
    return history
