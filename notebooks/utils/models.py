"""Shared model architecture + persistence helpers for notebook 4 (fits and saves the
pigmentation classifiers) and notebook 5 (loads and evaluates them). Keeping `TinyMLP`'s
definition here -- instead of duplicated inline in both notebooks -- guarantees the
architecture used to save a state_dict is identical to the one used to reload it.
"""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import torch
import torch.nn as nn


class TinyMLP(nn.Module):
    """Deliberately tiny: one hidden layer, since the point is to sanity-check whether any
    small net can exploit the full concatenated signal, not to tune an architecture."""

    def __init__(self, input_dim: int, hidden_dim: int = 32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


def save_sklearn_model(model, path: Path) -> None:
    joblib.dump(model, path)


def load_sklearn_model(path: Path):
    return joblib.load(path)


def save_nn_model(model: "TinyMLP", path: Path) -> None:
    torch.save(model.state_dict(), path)


def load_nn_model(path: Path, input_dim: int, hidden_dim: int = 32, device: str = "cpu") -> "TinyMLP":
    model = TinyMLP(input_dim, hidden_dim).to(device)
    model.load_state_dict(torch.load(path, map_location=device))
    model.eval()
    return model


def save_metadata(metadata: dict, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


def load_metadata(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)
