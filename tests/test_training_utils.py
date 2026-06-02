from types import SimpleNamespace

import torch

from genomics_pipeline.training_utils import (
    append_history_epoch,
    make_lr_scheduler,
    new_training_history,
    step_lr_scheduler,
    write_training_history,
)


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
