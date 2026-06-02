# -*- coding: utf-8 -*-
"""
cnn2_model.py
=============

CNN2 multi-estágio para dados organizados como genes x tracks x posições.

Arquitetura equivalente à CNN2 histórica de ``neural_ancestry_predictor``:
Stage1 agrupa todas as tracks de um gene em uma janela local, preservando a
altura como eixo de genes. Stages 2/3 processam apenas a dimensão genômica.
"""
from typing import Tuple

import torch
import torch.nn as nn
from rich.console import Console

from genotype_based_predictor.config import PipelineConfig
from genotype_based_predictor.models.nn_model import _resolve_gene_rows

console = Console()


class CNN2AncestryPredictor(nn.Module):
    """CNN2 com Stage1 k=(tracks_per_gene, width) para ver o gene completo."""

    def __init__(self, config: PipelineConfig, input_shape: Tuple[int, int], num_classes: int):
        super().__init__()
        self.config = config
        self.num_classes = num_classes

        self.genes_selected, self.rows_to_use = _resolve_gene_rows(config, input_shape)
        num_rows = len(self.rows_to_use)
        effective_size = input_shape[1]
        self.input_shape = (num_rows, effective_size)

        cnn2 = config.model.cnn2
        f1 = cnn2.num_filters_stage1
        f2 = cnn2.num_filters_stage2
        f3 = cnn2.num_filters_stage3
        kernel_s1 = tuple(cnn2.kernel_stage1)
        stride_s1 = tuple(cnn2.stride_stage1)
        kernel_s23 = tuple(cnn2.kernel_stages23)
        stride_s23 = tuple(cnn2.stride_stages23)
        padding_s23 = tuple(cnn2.padding_stages23)
        pool_type = cnn2.global_pool_type
        fc_hidden = cnn2.fc_hidden_size
        dropout_rate = config.model.dropout_rate

        if len(kernel_s1) != 2 or len(stride_s1) != 2:
            raise ValueError("cnn2.kernel_stage1 e cnn2.stride_stage1 devem ter 2 dimensoes")
        if len(kernel_s23) != 2 or len(stride_s23) != 2 or len(padding_s23) != 2:
            raise ValueError("cnn2.kernel/stride/padding_stages23 devem ter 2 dimensoes")
        if num_rows % kernel_s1[0] != 0:
            raise ValueError(
                f"CNN2 Stage1 requer num_rows multiplo de kernel vertical ({kernel_s1[0]}), "
                f"recebeu num_rows={num_rows}"
            )

        conv1_h = (num_rows - kernel_s1[0]) // stride_s1[0] + 1
        conv1_w = (effective_size - kernel_s1[1]) // stride_s1[1] + 1
        conv2_h = (conv1_h + 2 * padding_s23[0] - kernel_s23[0]) // stride_s23[0] + 1
        conv2_w = (conv1_w + 2 * padding_s23[1] - kernel_s23[1]) // stride_s23[1] + 1
        conv3_h = (conv2_h + 2 * padding_s23[0] - kernel_s23[0]) // stride_s23[0] + 1
        conv3_w = (conv2_w + 2 * padding_s23[1] - kernel_s23[1]) // stride_s23[1] + 1
        if min(conv1_h, conv1_w, conv2_h, conv2_w, conv3_h, conv3_w) <= 0:
            raise ValueError(
                "CNN2 produziu dimensao invalida: "
                f"input=({num_rows}, {effective_size}), "
                f"stage1=({conv1_h}, {conv1_w}), stage2=({conv2_h}, {conv2_w}), "
                f"stage3=({conv3_h}, {conv3_w})"
            )

        self.features = nn.Sequential(
            nn.Conv2d(1, f1, kernel_size=kernel_s1, stride=stride_s1),
            nn.ReLU(),
            nn.Conv2d(f1, f2, kernel_size=kernel_s23, stride=stride_s23, padding=padding_s23),
            nn.ReLU(),
            nn.Conv2d(f2, f3, kernel_size=kernel_s23, stride=stride_s23, padding=padding_s23),
            nn.ReLU(),
        )
        if pool_type == "avg":
            self.global_pool = nn.AvgPool2d(kernel_size=(1, conv3_w))
        else:
            self.global_pool = nn.MaxPool2d(kernel_size=(1, conv3_w))
        self.pool_type = pool_type

        flattened_size = f3 * conv3_h
        self.classifier = nn.Sequential(
            nn.Linear(flattened_size, fc_hidden),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(fc_hidden, num_classes),
        )

        self._initialize_weights()
        console.print("[green]Modelo CNN2 criado:[/green]")
        console.print(f"  • Input shape efetivo: {num_rows} x {effective_size} (1 canal)")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {num_rows} linhas")
        console.print(f"  • Stage 1: {f1} filters, kernel={kernel_s1}, stride={stride_s1} → {f1} x {conv1_h} x {conv1_w}")
        console.print(f"  • Stage 2: {f2} filters, kernel={kernel_s23}, stride={stride_s23}, padding={padding_s23} → {f2} x {conv2_h} x {conv2_w}")
        console.print(f"  • Stage 3: {f3} filters, kernel={kernel_s23}, stride={stride_s23}, padding={padding_s23} → {f3} x {conv3_h} x {conv3_w}")
        console.print(f"  • Global Pool ({pool_type}): kernel=(1, {conv3_w}) → {f3} x {conv3_h} x 1")
        console.print(f"  • FC: {flattened_size} → {fc_hidden} → {num_classes}")
        console.print(f"  • Params: {self.count_parameters():,}")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim == 4:
            x = x.view(x.size(0), x.size(1) * x.size(2), x.size(3))
        elif x.ndim != 3:
            raise ValueError(f"CNN2AncestryPredictor espera entrada (B,rows,L) ou (B,2,rows,L), recebeu shape={tuple(x.shape)}")
        x = x[:, self.rows_to_use, :]
        x = x.unsqueeze(1)
        x = self.features(x)
        x = self.global_pool(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)

    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        return torch.softmax(self.forward(x), dim=1)

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out", nonlinearity="relu")
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
