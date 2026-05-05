# -*- coding: utf-8 -*-
"""
nn_model.py
===========

Rede neural totalmente conectada (MLP) para predição de ancestralidade.

Arquitetura
-----------
Flatten (2D → 1D) → Camadas ocultas configuráveis → Saída (logits)

Suporte a seleção de genes específicos via `genes_to_use` no config.
"""
from typing import List, Optional, Tuple

import torch
import torch.nn as nn
from rich.console import Console

from genotype_based_predictor.config import PipelineConfig

console = Console()


def _resolve_gene_rows(config: PipelineConfig, input_shape: Tuple[int, int]) -> Tuple[List[str], List[int]]:
    """Determina quais linhas (tracks) usar com base nos genes configurados."""
    if config.dataset_input.tensor_layout == "haplotype_channels":
        genes_to_use = list(config.dataset_input.genes_to_use or [])
        if not genes_to_use:
            raise ValueError("tensor_layout='haplotype_channels' requer ao menos um gene em genes_to_use")
        num_ontologies = len(config.dataset_input.ontology_terms or []) or 1
        expected_rows = (2 * num_ontologies + 6) * len(genes_to_use)
        if input_shape[0] != expected_rows:
            raise ValueError(
                f"tensor_layout='haplotype_channels' espera {expected_rows} linhas derivadas de {len(genes_to_use)} gene(s) e {num_ontologies} ontologia(s), recebeu {input_shape[0]}"
            )
        return genes_to_use, list(range(input_shape[0]))

    GENE_ORDER = config.dataset_input.gene_order or [
        "MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"
    ]
    if not GENE_ORDER or input_shape[0] % len(GENE_ORDER) != 0:
        raise ValueError(
            f"Input shape incompatível com gene_order: {input_shape[0]} linhas para {len(GENE_ORDER)} genes"
        )
    tracks_per_gene = input_shape[0] // len(GENE_ORDER)
    genes_to_use = config.dataset_input.genes_to_use

    if genes_to_use:
        for g in genes_to_use:
            if g not in GENE_ORDER:
                raise ValueError(f"Gene inválido: {g}. Opções: {GENE_ORDER}")
        gene_indices = [i for i, g in enumerate(GENE_ORDER) if g in genes_to_use]
        genes_selected = [GENE_ORDER[i] for i in gene_indices]
        rows_to_use = [r for gi in gene_indices for r in range(gi * tracks_per_gene, (gi + 1) * tracks_per_gene)]
    else:
        genes_selected = list(GENE_ORDER)
        rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))

    return genes_selected, rows_to_use


def _make_activation(activation_type: str) -> nn.Module:
    """Retorna módulo de ativação pelo nome."""
    if activation_type == "relu":
        return nn.ReLU()
    elif activation_type == "tanh":
        return nn.Tanh()
    elif activation_type == "sigmoid":
        return nn.Sigmoid()
    raise ValueError(f"Ativação não suportada: {activation_type}")


class NNAncestryPredictor(nn.Module):
    """
    Rede neural totalmente conectada (MLP) para predição de ancestralidade.

    Recebe features 2D [batch, num_rows, effective_size], faz flatten e
    passa por camadas ocultas configuráveis até a saída (logits).

    Parameters
    ----------
    config : PipelineConfig
        Configuração completa (dataset_input, model, output).
    input_shape : Tuple[int, int]
        Shape de entrada (num_rows, effective_size).
    num_classes : int
        Número de classes (ou tamanho da saída para regressão).
    """

    def __init__(self, config: PipelineConfig, input_shape: Tuple[int, int], num_classes: int):
        super().__init__()
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config.output.prediction_target != "frog_likelihood"

        self.genes_selected, self.rows_to_use = _resolve_gene_rows(config, input_shape)

        num_rows = len(self.rows_to_use)
        effective_size = input_shape[1]
        self.input_shape = (num_rows, effective_size)
        self.input_size = num_rows * effective_size

        hidden_layers = config.model.hidden_layers
        activation_type = config.model.activation
        dropout_rate = config.model.dropout_rate
        activation = _make_activation(activation_type)

        layers = []
        prev_size = self.input_size
        for hidden_size in hidden_layers:
            layers += [nn.Linear(prev_size, hidden_size), activation]
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        layers.append(nn.Linear(prev_size, num_classes))
        self.network = nn.Sequential(*layers)

        self._initialize_weights(activation_type)

        console.print(f"[green]Modelo NN criado:[/green] {num_rows}×{effective_size} → flatten({self.input_size}) → {hidden_layers} → {num_classes}")
        console.print(f"  • Genes: {', '.join(self.genes_selected)} | Params: {self.count_parameters():,}")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 4:
            raise ValueError(f"NNAncestryPredictor espera entrada (B,2,4,L), recebeu shape={tuple(x.shape)}")
        x = x.view(x.size(0), x.size(1) * x.size(2), x.size(3))
        x = x[:, self.rows_to_use, :]
        x = x.view(x.size(0), -1)
        return self.network(x)

    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """Retorna probabilidades (softmax sobre logits)."""
        return torch.softmax(self.forward(x), dim=1)

    def count_parameters(self) -> int:
        """Número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def _initialize_weights(self, activation_type: str):
        total = sum(1 for m in self.modules() if isinstance(m, nn.Linear))
        count = 0
        for m in self.modules():
            if isinstance(m, nn.Linear):
                count += 1
                if count == total:
                    nn.init.xavier_normal_(m.weight)
                elif activation_type == "relu":
                    nn.init.kaiming_normal_(m.weight, mode="fan_in", nonlinearity="relu")
                else:
                    nn.init.xavier_normal_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
