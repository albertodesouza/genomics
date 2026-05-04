# -*- coding: utf-8 -*-
"""
cnn2_model.py
=============

CNN multi-estágio com Global Average/Max Pooling para predição de
ancestralidade a partir de dados genômicos 2D.

Arquitetura
-----------
Stage1: Conv2D(kernel_stage1) → BN → ReLU × n1
Stage2: Conv2D(1×kernel_stage2) → BN → ReLU × n2
Stage3: Conv2D(kernel_stage3×1) → BN → ReLU × n3
→ Global Pool → Flatten → FC → Dropout → Logits
"""
from typing import List, Tuple

import torch
import torch.nn as nn
from rich.console import Console

from genotype_based_predictor.config import PipelineConfig
from genotype_based_predictor.models.nn_model import _resolve_gene_rows

console = Console()


class CNN2AncestryPredictor(nn.Module):
    """
    CNN2 multi-estágio com separação de dimensões horizontais/verticais.

    Stage1 captura padrões locais 2D; Stage2 integra ao longo da dimensão
    genômica (colunas); Stage3 integra ao longo das tracks (linhas).
    Global pooling remove dependência do tamanho de entrada.

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

        cnn2 = config.model.cnn2
        f1 = cnn2.num_filters_stage1
        ks1 = cnn2.kernel_stage1
        ks2 = cnn2.kernel_stage2
        f2 = cnn2.num_filters_stage2
        ks3 = cnn2.kernel_stage3
        f3 = cnn2.num_filters_stage3
        pool_type = cnn2.global_pool_type
        fc_hidden = cnn2.fc_hidden_size
        dropout_rate = config.model.dropout_rate
        num_classes_fc = num_classes

        # Stage 1: 2D conv para capturar padrões locais
        self.stage1 = nn.Sequential(
            nn.Conv2d(1, f1, kernel_size=tuple(ks1), padding=(ks1[0] // 2, ks1[1] // 2), bias=False),
            nn.BatchNorm2d(f1),
            nn.ReLU(),
        )
        # Stage 2: integra na dimensão horizontal
        self.stage2 = nn.Sequential(
            nn.Conv2d(f1, f2, kernel_size=(1, ks2), padding=(0, ks2 // 2), bias=False),
            nn.BatchNorm2d(f2),
            nn.ReLU(),
        )
        # Stage 3: integra na dimensão vertical (tracks)
        self.stage3 = nn.Sequential(
            nn.Conv2d(f2, f3, kernel_size=(ks3, 1), padding=(ks3 // 2, 0), bias=False),
            nn.BatchNorm2d(f3),
            nn.ReLU(),
        )

        # Global pooling
        if pool_type == "avg":
            self.global_pool = nn.AdaptiveAvgPool2d((1, 1))
        else:
            self.global_pool = nn.AdaptiveMaxPool2d((1, 1))

        # FC head
        self.fc = nn.Sequential(
            nn.Linear(f3, fc_hidden),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(fc_hidden, num_classes_fc),
        )

        self._initialize_weights()
        console.print(f"[green]Modelo CNN2 criado:[/green] {num_rows}×{effective_size} → s1f{f1}/s2f{f2}/s3f{f3} → gpool → fc{fc_hidden} → {num_classes}")
        console.print(f"  • Params: {self.count_parameters():,}")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 4:
            raise ValueError(f"CNN2AncestryPredictor espera entrada (B,2,4,L), recebeu shape={tuple(x.shape)}")
        x = x.view(x.size(0), x.size(1) * x.size(2), x.size(3))
        x = x[:, self.rows_to_use, :]  # selecionar genes
        x = x.unsqueeze(1)             # (B, 1, rows, cols)
        x = self.stage1(x)
        x = self.stage2(x)
        x = self.stage3(x)
        x = self.global_pool(x)        # (B, f3, 1, 1)
        x = x.view(x.size(0), -1)     # (B, f3)
        return self.fc(x)

    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        return torch.softmax(self.forward(x), dim=1)

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out", nonlinearity="relu")
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
