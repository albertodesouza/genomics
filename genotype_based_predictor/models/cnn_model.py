# -*- coding: utf-8 -*-
"""
cnn_model.py
============

CNN com uma camada convolucional 2D para predição de ancestralidade.

Arquitetura
-----------
Conv2D → ReLU → [MaxPool] → Flatten → MLP → Logits
"""
from typing import List, Optional, Tuple

import torch
import torch.nn as nn
from rich.console import Console

from genotype_based_predictor.config import PipelineConfig
from genotype_based_predictor.models.nn_model import _resolve_gene_rows, _make_activation

console = Console()


class CNNAncestryPredictor(nn.Module):
    """
    CNN com Conv2D para dados genômicos 2D (tracks × posições).

    Parameters
    ----------
    config : PipelineConfig
        Configuração completa.
    input_shape : Tuple[int, int]
        Shape de entrada (num_rows, effective_size).
    num_classes : int
        Número de classes de saída.
    """

    def __init__(self, config: PipelineConfig, input_shape: Tuple[int, int], num_classes: int):
        super().__init__()
        self.config = config
        self.num_classes = num_classes

        self.genes_selected, self.rows_to_use = _resolve_gene_rows(config, input_shape)
        num_rows = len(self.rows_to_use)
        effective_size = input_shape[1]
        self.input_shape = (num_rows, effective_size)

        cnn_cfg = config.model.cnn
        kernel_size = tuple(cnn_cfg.kernel_size)
        num_filters = cnn_cfg.num_filters
        stride = tuple(cnn_cfg.stride) if isinstance(cnn_cfg.stride, list) else (cnn_cfg.stride, cnn_cfg.stride)
        padding = tuple(cnn_cfg.padding) if isinstance(cnn_cfg.padding, list) else (cnn_cfg.padding, cnn_cfg.padding)
        pool_size = tuple(cnn_cfg.pool_size) if cnn_cfg.pool_size else None

        self.conv = nn.Conv2d(1, num_filters, kernel_size=kernel_size, stride=stride, padding=padding)
        self.relu = nn.ReLU()
        self.pool = nn.MaxPool2d(kernel_size=pool_size) if pool_size else None

        with torch.no_grad():
            dummy = torch.zeros(1, 1, num_rows, effective_size)
            out = self.relu(self.conv(dummy))
            if self.pool is not None:
                out = self.pool(out)
            flat_size = out.view(1, -1).shape[1]

        activation_type = config.model.activation
        hidden_layers = config.model.hidden_layers
        dropout_rate = config.model.dropout_rate
        activation = _make_activation(activation_type)

        layers = []
        prev = flat_size
        for h in hidden_layers:
            layers += [nn.Linear(prev, h), activation]
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            prev = h
        layers.append(nn.Linear(prev, num_classes))
        self.fc = nn.Sequential(*layers)

        self._initialize_weights(activation_type)
        console.print(f"[green]Modelo CNN criado:[/green] {num_rows}×{effective_size} → conv{kernel_size}/{num_filters}f → flatten({flat_size}) → {num_classes}")
        console.print(f"  • Params: {self.count_parameters():,}")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 4:
            raise ValueError(f"CNNAncestryPredictor espera entrada (B,2,4,L), recebeu shape={tuple(x.shape)}")
        x = x.view(x.size(0), x.size(1) * x.size(2), x.size(3))
        x = x[:, self.rows_to_use, :]
        x = x.unsqueeze(1)  # (B, 1, rows, cols)
        x = self.relu(self.conv(x))
        if self.pool is not None:
            x = self.pool(x)
        x = x.view(x.size(0), -1)
        return self.fc(x)

    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        return torch.softmax(self.forward(x), dim=1)

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def _initialize_weights(self, activation_type: str):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out", nonlinearity="relu")
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
