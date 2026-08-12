from types import SimpleNamespace

import torch

from genomics.core.training_utils import (
    EpochTrainer,
    append_history_epoch,
    make_lr_scheduler,
    new_training_history,
    step_lr_scheduler,
    write_training_history,
)


class DummyEpochTrainer(EpochTrainer):
    def __init__(self, train_metrics, val_metrics, **kwargs):
        super().__init__(**kwargs)
        self.train_metrics = list(train_metrics)
        self.val_metrics = list(val_metrics)
        self.saved = []

    def run_train_epoch(self):
        return self.train_metrics.pop(0)

    def run_val_epoch(self):
        return self.val_metrics.pop(0)

    def save_checkpoint(self, name, epoch, metrics):
        self.saved.append((name, epoch, dict(metrics)))


def test_training_history_helpers(tmp_path):
    history = new_training_history()
    append_history_epoch(history, 1, {"loss": 0.5, "accuracy": 0.75}, {"loss": 0.4, "accuracy": 0.8})

    assert history == {
        "epoch": [1],
        "train_loss": [0.5],
        "train_accuracy": [0.75],
        "val_loss": [0.4],
        "val_accuracy": [0.8],
    }

    output_path = tmp_path / "models" / "training_history.json"
    write_training_history(output_path, history)
    assert output_path.exists()


def test_make_and_step_plateau_scheduler():
    model = torch.nn.Linear(2, 1)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
    scheduler = make_lr_scheduler(
        optimizer,
        SimpleNamespace(enabled=True, type="plateau", mode="min", factor=0.5, patience=0, min_lr=0.0),
    )

    step_lr_scheduler(scheduler, 1.0)
    step_lr_scheduler(scheduler, 2.0)

    assert optimizer.param_groups[0]["lr"] == 0.05


def test_make_lr_scheduler_disabled_returns_none():
    model = torch.nn.Linear(2, 1)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1)
    assert make_lr_scheduler(optimizer, {"enabled": False}) is None


def test_epoch_trainer_bookkeeping_and_checkpoints(tmp_path):
    trainer = DummyEpochTrainer(
        train_metrics=[
            {"loss": 0.9, "accuracy": 0.1, "samples": 2},
            {"loss": 0.8, "accuracy": 0.2, "samples": 2},
        ],
        val_metrics=[
            {"loss": 0.7, "accuracy": 0.3, "samples": 2},
            {"loss": 0.6, "accuracy": 0.4, "samples": 2},
        ],
        num_epochs=2,
        save_frequency=1,
        history_path=tmp_path / "history.json",
    )

    history = trainer.train()

    assert history["epoch"] == [1, 2]
    assert history["train_loss"] == [0.9, 0.8]
    assert history["val_accuracy"] == [0.3, 0.4]
    assert history["last_epoch"] == 2
    assert trainer.best_val_accuracy == 0.4
    assert trainer.best_val_loss == 0.6
    assert (tmp_path / "history.json").exists()
    saved_names = [item[0] for item in trainer.saved]
    assert saved_names == ["best_accuracy.pt", "best_loss.pt", "epoch_1.pt", "best_accuracy.pt", "best_loss.pt", "epoch_2.pt", "final.pt"]


def test_epoch_trainer_validation_frequency_skips_best_update():
    trainer = DummyEpochTrainer(
        train_metrics=[
            {"loss": 0.9, "accuracy": 0.1, "samples": 2},
            {"loss": 0.8, "accuracy": 0.2, "samples": 2},
        ],
        val_metrics=[{"loss": 0.6, "accuracy": 0.4, "samples": 2}],
        num_epochs=2,
        validation_frequency=2,
        save_frequency=0,
    )

    history = trainer.train()

    assert history["val_loss"] == [float("inf"), 0.6]
    assert trainer.best_val_accuracy == 0.4
    assert [item[0] for item in trainer.saved] == ["best_accuracy.pt", "best_loss.pt", "final.pt"]
