#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
neural_ancestry_predictor.py

Rede Neural para Predição de Ancestralidade a partir de Dados AlphaGenome
==========================================================================

Este módulo implementa uma rede neural configurável via YAML que prediz
ancestralidade (superpopulation, population ou FROG likelihood) a partir
de predições AlphaGenome armazenadas em um dataset PyTorch.

Uso:
    # Treino
    python3 neural_ancestry_predictor.py --config configs/default.yaml

    # Teste
    python3 neural_ancestry_predictor.py --config configs/default.yaml --mode test

Author: ChatGPT (for Alberto)
Created: 2025-11-14
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings
from datetime import datetime

import numpy as np
import yaml
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, Subset
from sklearn.metrics import (
    accuracy_score, precision_recall_fscore_support, 
    confusion_matrix, classification_report
)
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.table import Table
from rich.panel import Panel
import matplotlib
matplotlib.use('TkAgg')  # Backend interativo
import matplotlib.pyplot as plt

# Importar dataset genômico
sys.path.insert(0, str(Path(__file__).parent.parent / "build_non_longevous_dataset"))
from genomic_dataset import GenomicLongevityDataset

console = Console()

# ==============================================================================
# DATA PROCESSING FUNCTIONS
# ==============================================================================

def minmax_keep_zero(x: torch.Tensor, xmax: float) -> torch.Tensor:
    """
    Normalização Min-Max onde:
       - zeros permanecem zeros
       - valores não-zero são divididos pelo máximo não-zero
    
    Args:
        x: tensor PyTorch (qualquer formato)
        xmax: máximo dos valores não-zero (pré-computado)
        
    Returns:
        tensor normalizado no mesmo device (CPU ou GPU)
    """
    if xmax == 0:
        return x
    return x / xmax


def log_normalize(x: torch.Tensor, log_max: float) -> torch.Tensor:
    """
    Normalização logarítmica com log1p para tensores PyTorch.
    Mantém zeros como zeros automaticamente.
    
    Args:
        x: tensor PyTorch
        log_max: log1p do máximo (pré-computado)
        
    Returns:
        tensor normalizado
    """
    if log_max == 0:
        return x
    return torch.log1p(x) / log_max


def zscore_normalize(x: torch.Tensor, mean: float, std: float) -> torch.Tensor:
    """
    Normalização Z-score padrão.
    
    Args:
        x: tensor PyTorch
        mean: média global
        std: desvio padrão global
        
    Returns:
        tensor normalizado
    """
    if std == 0:
        return x
    return (x - mean) / std


def taint_sample(features: torch.Tensor, target_class: int, num_classes: int) -> torch.Tensor:
    """
    Marca amostra com valor sentinela para debug.
    
    Coloca -2 na posição específica baseada na classe:
    posição = classe * (input_vector_size / num_classes)
    
    Exemplo:
        Para 5 classes e vetor de 5500 elementos:
        - Classe 0 -> posição 0 * (5500/5) = 0
        - Classe 1 -> posição 1 * (5500/5) = 1100
        - Classe 2 -> posição 2 * (5500/5) = 2200
        - etc.
    
    Args:
        features: Tensor de entrada
        target_class: Classe de saída (0 a num_classes-1)
        num_classes: Número total de classes
        
    Returns:
        Tensor com tainting aplicado
    """
    features = features.clone()
    input_size = features.shape[0] if features.ndim == 1 else features.numel()
    
    # Calcular posição: classe * (input_size / num_classes)
    position = int(target_class * (input_size / num_classes))
    
    # Garantir que posição está no range válido
    if position < input_size:
        if features.ndim == 1:
            features[position] = -2.0
        else:
            features.view(-1)[position] = -2.0
    
    return features


class ProcessedGenomicDataset(Dataset):
    """
    Dataset wrapper que processa dados genômicos conforme configuração.
    
    Aplica:
    - Extração de trecho central das janelas
    - Downsampling
    - Combinação de haplótipos
    - Normalização
    """
    
    def __init__(
        self,
        base_dataset: GenomicLongevityDataset,
        config: Dict,
        normalization_params: Optional[Dict] = None,
        compute_normalization: bool = True
    ):
        """
        Inicializa dataset processado.
        
        Args:
            base_dataset: Dataset genômico base
            config: Configuração do YAML
            normalization_params: Parâmetros de normalização pré-computados
            compute_normalization: Se True, computa parâmetros de normalização
        """
        self.base_dataset = base_dataset
        self.config = config
        
        # Parâmetros de processamento
        self.alphagenome_outputs = config['dataset_input']['alphagenome_outputs']
        self.haplotype_mode = config['dataset_input']['haplotype_mode']
        self.window_center_size = config['dataset_input']['window_center_size']
        self.downsample_factor = config['dataset_input']['downsample_factor']
        self.prediction_target = config['output']['prediction_target']
        self.normalization_method = config['dataset_input'].get('normalization_method', 'zscore')
        
        # Mapeamentos para targets
        self.target_to_idx = {}
        self.idx_to_target = {}
        
        # Computar ou carregar parâmetros de normalização
        if normalization_params is not None:
            self.normalization_params = normalization_params
        elif compute_normalization:
            console.print("[yellow]Computando parâmetros de normalização...[/yellow]")
            self.normalization_params = self._compute_normalization_params()
        else:
            self.normalization_params = {'mean': 0.0, 'std': 1.0}
        
        # Criar mapeamentos de targets
        self._create_target_mappings()
    
    def _compute_normalization_params(self) -> Dict:
        """
        Computa média e desvio padrão de todos os dados para normalização.
        
        Returns:
            Dict com 'mean' e 'std'
        """
        all_values = []
        num_processed = 0
        num_errors = 0
        
        console.print(f"[cyan]Iniciando computação de normalização para {len(self.base_dataset)} amostras...[/cyan]")
        console.print(f"[cyan]Método de normalização: {self.normalization_method}[/cyan]")
        console.print(f"[cyan]Outputs a processar: {', '.join(self.alphagenome_outputs)}[/cyan]")
        console.print(f"[cyan]Modo de haplótipo: {self.haplotype_mode}[/cyan]")
        console.print(f"[cyan]Tamanho do trecho central: {self.window_center_size} bases[/cyan]")
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task(
                "Computando estatísticas para normalização...",
                total=len(self.base_dataset)
            )
            
            for idx in range(len(self.base_dataset)):
                try:
                    input_data, output_data = self.base_dataset[idx]
                    sample_id = output_data.get('sample_id', f'sample_{idx}')
                    
                    # Informar qual amostra está sendo processada
                    if idx % 10 == 0:  # A cada 10 amostras
                        console.print(f"[dim]  Processando amostra {idx + 1}/{len(self.base_dataset)}: {sample_id}[/dim]")
                    
                    processed_array = self._process_windows(input_data['windows'])
                    
                    if len(processed_array) > 0:
                        all_values.append(processed_array)
                        num_processed += 1
                    else:
                        console.print(f"[yellow]  ⚠ Amostra {sample_id} (idx={idx}): Nenhum dado válido encontrado[/yellow]")
                        num_errors += 1
                        
                except Exception as e:
                    num_errors += 1
                    console.print(f"[yellow]  ⚠ Erro ao processar amostra {idx}: {e}[/yellow]")
                progress.update(task, advance=1)
        
        if len(all_values) == 0:
            console.print(f"[red]ERRO: Nenhuma amostra válida para normalização![/red]")
            console.print(f"[red]  • Amostras processadas: {num_processed}[/red]")
            console.print(f"[red]  • Amostras com erro: {num_errors}[/red]")
            console.print(f"[red]  • Total esperado: {len(self.base_dataset)}[/red]")
            return {'method': self.normalization_method, 'mean': 0.0, 'std': 1.0}
        
        # Concatenar todos os arrays
        all_values = np.concatenate(all_values)
        
        # Calcular parâmetros dependendo do método escolhido
        params = {'method': self.normalization_method}
        
        if self.normalization_method == 'zscore':
            # Z-score: mean e std
            mean = float(np.mean(all_values))
            std = float(np.std(all_values))
            if std < 1e-8:
                std = 1.0
            params['mean'] = mean
            params['std'] = std
            
            console.print(f"\n[bold green]✓ Normalização Z-score Concluída:[/bold green]")
            console.print(f"  • Amostras processadas com sucesso: {num_processed}/{len(self.base_dataset)}")
            console.print(f"  • Amostras com erro: {num_errors}")
            console.print(f"  • Total de valores coletados: {len(all_values):,}")
            console.print(f"  • Média (mean): {mean:.6f}")
            console.print(f"  • Desvio padrão (std): {std:.6f}")
            
        elif self.normalization_method == 'minmax_keep_zero':
            # MinMax mantendo zeros: apenas o máximo dos valores não-zero
            nonzero_values = all_values[all_values > 0]
            if len(nonzero_values) > 0:
                xmax = float(nonzero_values.max())
            else:
                xmax = 1.0
            params['max'] = xmax
            
            console.print(f"\n[bold green]✓ Normalização MinMax (mantendo zeros) Concluída:[/bold green]")
            console.print(f"  • Amostras processadas com sucesso: {num_processed}/{len(self.base_dataset)}")
            console.print(f"  • Amostras com erro: {num_errors}")
            console.print(f"  • Total de valores coletados: {len(all_values):,}")
            console.print(f"  • Valores zeros: {(all_values == 0).sum():,}")
            console.print(f"  • Máximo não-zero: {xmax:.6f}")
            
        elif self.normalization_method == 'log':
            # Log: log1p do máximo
            nonzero_values = all_values[all_values > 0]
            if len(nonzero_values) > 0:
                xmax = float(nonzero_values.max())
                log_max = float(np.log1p(xmax))
            else:
                log_max = 1.0
            params['log_max'] = log_max
            
            console.print(f"\n[bold green]✓ Normalização Logarítmica Concluída:[/bold green]")
            console.print(f"  • Amostras processadas com sucesso: {num_processed}/{len(self.base_dataset)}")
            console.print(f"  • Amostras com erro: {num_errors}")
            console.print(f"  • Total de valores coletados: {len(all_values):,}")
            console.print(f"  • Valores zeros: {(all_values == 0).sum():,}")
            console.print(f"  • log1p(max): {log_max:.6f}")
            
        else:
            console.print(f"[red]Método de normalização desconhecido: {self.normalization_method}[/red]")
            console.print(f"[yellow]Usando zscore como fallback[/yellow]")
            mean = float(np.mean(all_values))
            std = float(np.std(all_values))
            if std < 1e-8:
                std = 1.0
            params['method'] = 'zscore'
            params['mean'] = mean
            params['std'] = std
        
        return params
    
    def _create_target_mappings(self):
        """Cria mapeamentos entre targets e índices."""
        console.print(f"\n[cyan]Criando mapeamento de classes...[/cyan]")
        unique_targets = set()
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task(
                "Escaneando targets...",
                total=len(self.base_dataset)
            )
            
            for idx in range(len(self.base_dataset)):
                try:
                    _, output_data = self.base_dataset[idx]
                    target = self._get_target_value(output_data)
                    if target is not None:
                        unique_targets.add(target)
                except Exception:
                    pass
                progress.update(task, advance=1)
        
        # Ordenar para consistência
        sorted_targets = sorted(list(unique_targets))
        
        self.target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
        self.idx_to_target = {idx: target for target, idx in self.target_to_idx.items()}
        
        console.print(f"[green]✓ Mapeamento de targets criado: {len(self.target_to_idx)} classes[/green]")
        console.print(f"[cyan]Classes: {sorted_targets}[/cyan]")
    
    def _get_target_value(self, output_data: Dict) -> Optional[str]:
        """
        Extrai valor do target dos dados de saída.
        
        Args:
            output_data: Dicionário de saída do dataset
            
        Returns:
            Target como string ou None se FROG likelihood
        """
        if self.prediction_target == 'superpopulation':
            return output_data.get('superpopulation')
        elif self.prediction_target == 'population':
            return output_data.get('population')
        elif self.prediction_target == 'frog_likelihood':
            return None  # Regressão, não há classes
        else:
            raise ValueError(f"prediction_target inválido: {self.prediction_target}")
    
    def _process_windows(self, windows: Dict) -> np.ndarray:
        """
        Processa todas as janelas e concatena em um único array.
        
        Args:
            windows: Dicionário de janelas do input_data
            
        Returns:
            Array numpy concatenado com todos os features
        """
        processed_arrays = []
        
        for window_name, window_data in windows.items():
            # Processar cada haplótipo
            if self.haplotype_mode in ['H1', 'H1+H2']:
                h1_array = self._process_haplotype(window_data.get('predictions_h1', {}))
                if h1_array is not None:
                    processed_arrays.append(h1_array)
            
            if self.haplotype_mode in ['H2', 'H1+H2']:
                h2_array = self._process_haplotype(window_data.get('predictions_h2', {}))
                if h2_array is not None:
                    processed_arrays.append(h2_array)
        
        if len(processed_arrays) == 0:
            # Retornar array vazio com dimensão correta
            return np.array([])
        
        # Concatenar todos os arrays
        result = np.concatenate(processed_arrays)
        
        # Garantir que resultado é 1D
        if result.ndim > 1:
            result = result.flatten()
        
        return result
    
    def _process_haplotype(self, predictions: Dict) -> Optional[np.ndarray]:
        """
        Processa predições de um haplótipo.
        
        Args:
            predictions: Dict com {output_type: array}
            
        Returns:
            Array processado ou None
        """
        arrays_to_concat = []
        
        for output_type in self.alphagenome_outputs:
            if output_type in predictions:
                array = predictions[output_type]
                
                # Garantir que é 1D
                if array.ndim > 1:
                    array = array.flatten()
                
                # Extrair trecho central
                center_array = self._extract_center(array)
                
                # Aplicar downsampling
                downsampled = self._downsample(center_array)
                
                arrays_to_concat.append(downsampled)
        
        if len(arrays_to_concat) == 0:
            return None
        
        return np.concatenate(arrays_to_concat)
    
    def _extract_center(self, array: np.ndarray) -> np.ndarray:
        """
        Extrai trecho central do array.
        
        Args:
            array: Array de predições (tamanho original ~1M)
            
        Returns:
            Trecho central
        """
        if self.window_center_size <= 0 or self.window_center_size >= len(array):
            return array
        
        center_idx = len(array) // 2
        half_size = self.window_center_size // 2
        
        start = max(0, center_idx - half_size)
        end = min(len(array), center_idx + half_size)
        
        return array[start:end]
    
    def _downsample(self, array: np.ndarray) -> np.ndarray:
        """
        Aplica downsampling no array.
        
        Args:
            array: Array a ser downsampled
            
        Returns:
            Array downsampled
        """
        if self.downsample_factor <= 1:
            return array
        
        return array[::self.downsample_factor]
    
    def __len__(self) -> int:
        return len(self.base_dataset)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Retorna item processado.
        
        Returns:
            Tupla (features_tensor, target_tensor)
        """
        input_data, output_data = self.base_dataset[idx]
        
        # Processar janelas
        features = self._process_windows(input_data['windows'])
        
        # Garantir que features é 1D (flatten se necessário)
        if features.ndim > 1:
            features = features.flatten()
        
        # Converter para tensor
        features_tensor = torch.FloatTensor(features)
        
        # Aplicar normalização conforme método escolhido
        method = self.normalization_params.get('method', 'zscore')
        
        if method == 'zscore':
            features_tensor = zscore_normalize(
                features_tensor, 
                self.normalization_params['mean'], 
                self.normalization_params['std']
            )
        elif method == 'minmax_keep_zero':
            features_tensor = minmax_keep_zero(
                features_tensor,
                self.normalization_params['max']
            )
        elif method == 'log':
            features_tensor = log_normalize(
                features_tensor,
                self.normalization_params['log_max']
            )
        else:
            # Fallback para zscore
            features_tensor = zscore_normalize(
                features_tensor,
                self.normalization_params.get('mean', 0.0),
                self.normalization_params.get('std', 1.0)
            )
        
        # Processar target
        if self.prediction_target == 'frog_likelihood':
            # Regressão: usar likelihood diretamente
            target = output_data.get('frog_likelihood', np.zeros(150))
            target_tensor = torch.FloatTensor(target)
        else:
            # Classificação: converter para índice (tensor escalar)
            target_value = self._get_target_value(output_data)
            if target_value in self.target_to_idx:
                target_idx = self.target_to_idx[target_value]
                target_tensor = torch.tensor(target_idx, dtype=torch.long)
                
                # Tainting em runtime (se habilitado) - apenas para classificação
                if self.config.get('debug', {}).get('taint_at_runtime', False):
                    num_classes = len(self.target_to_idx)
                    features_tensor = taint_sample(features_tensor, target_idx, num_classes)
            else:
                # Target desconhecido, usar -1
                target_tensor = torch.tensor(-1, dtype=torch.long)
        
        return features_tensor, target_tensor
    
    def get_num_classes(self) -> int:
        """Retorna número de classes."""
        return len(self.target_to_idx)
    
    def get_input_size(self) -> int:
        """Calcula tamanho da entrada da rede."""
        # Pegar primeira amostra válida para calcular
        for idx in range(len(self.base_dataset)):
            try:
                features, _ = self[idx]
                return features.shape[0]
            except Exception:
                continue
        
        # Fallback: calcular teoricamente
        num_windows = len(self.base_dataset[0][0]['windows'])
        num_outputs = len(self.alphagenome_outputs)
        num_haplotypes = 2 if self.haplotype_mode == 'H1+H2' else 1
        
        effective_size = self.window_center_size // self.downsample_factor
        
        return num_windows * num_outputs * num_haplotypes * effective_size


# ==============================================================================
# MODEL ARCHITECTURE
# ==============================================================================

class AncestryPredictor(nn.Module):
    """
    Rede neural para predição de ancestralidade.
    
    Arquitetura:
    - Camada de entrada (tamanho variável)
    - Camadas ocultas (configurável)
    - Camada de saída (softmax para classificação ou linear para regressão)
    """
    
    def __init__(self, config: Dict, input_size: int, num_classes: int):
        """
        Inicializa o modelo.
        
        Args:
            config: Configuração do YAML
            input_size: Tamanho da entrada
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(AncestryPredictor, self).__init__()
        
        self.config = config
        self.input_size = input_size
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Parâmetros de arquitetura
        hidden_layers = config['model']['hidden_layers']
        activation_type = config['model']['activation']
        dropout_rate = config['model']['dropout_rate']
        
        # Escolher função de ativação (para camadas intermediárias)
        if activation_type == 'relu':
            self.activation = nn.ReLU()
        elif activation_type == 'tanh':
            self.activation = nn.Tanh()
        elif activation_type == 'sigmoid':
            self.activation = nn.Sigmoid()
        else:
            raise ValueError(f"Ativação não suportada: {activation_type}")
        
        # Construir camadas
        # Arquitetura: camadas intermediárias (com ativação) + última hidden (LINEAR) + saída
        layers = []
        prev_size = input_size
        
        # Camadas intermediárias (todas exceto última)
        for i, hidden_size in enumerate(hidden_layers[:-1]):
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(self.activation)
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Última camada hidden (pré-softmax): sempre LINEAR (sem ativação)
        if len(hidden_layers) > 0:
            pre_softmax_size = hidden_layers[-1]
            layers.append(nn.Linear(prev_size, pre_softmax_size))
            # Sem ativação aqui - linear para softmax
            prev_size = pre_softmax_size
        
        # Camada de saída
        layers.append(nn.Linear(prev_size, num_classes))
        
        # Softmax para classificação (aplicado no forward)
        if self.is_classification:
            self.softmax = nn.Softmax(dim=1)
        
        self.network = nn.Sequential(*layers)
        
        # Inicializar pesos apropriadamente
        self._initialize_weights(activation_type)
        
        console.print(f"[green]Modelo criado:[/green]")
        console.print(f"  • Input size: {input_size}")
        console.print(f"  • Hidden layers: {hidden_layers}")
        console.print(f"  • Arquitetura: camadas intermediárias ({activation_type}) + pré-softmax (linear) + saída")
        console.print(f"  • Output size: {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward pass."""
        logits = self.network(x)
        
        if self.is_classification and not self.training:
            # Aplicar softmax apenas durante inferência
            return self.softmax(logits)
        
        return logits
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self, activation_type: str):
        """
        Inicializa pesos apropriadamente de acordo com tipo de ativação.
        
        Estratégia:
        - Camadas intermediárias (com ativação):
            * ReLU: He/Kaiming initialization (fan_in)
            * Tanh/Sigmoid: Xavier/Glorot initialization
        - Últimas 2 camadas (pré-softmax linear + saída): Xavier initialization
        - Bias: zeros (padrão recomendado)
        
        Args:
            activation_type: Tipo de função de ativação das camadas intermediárias
        """
        layer_count = 0
        total_layers = sum(1 for m in self.modules() if isinstance(m, nn.Linear))
        
        for m in self.modules():
            if isinstance(m, nn.Linear):
                layer_count += 1
                
                # Últimas 2 camadas (pré-softmax + saída): sempre Xavier
                # (independente da ativação das camadas anteriores)
                if layer_count >= total_layers - 1:
                    nn.init.xavier_normal_(m.weight)
                # Camadas intermediárias: depende da ativação configurada
                elif activation_type == 'relu':
                    # He initialization (Kaiming) para ReLU
                    # Usa fan_in para manter variância durante forward pass
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                elif activation_type in ['tanh', 'sigmoid']:
                    # Xavier/Glorot initialization para tanh/sigmoid
                    nn.init.xavier_normal_(m.weight)
                else:
                    # Fallback: Xavier para qualquer outra ativação
                    nn.init.xavier_normal_(m.weight)
                
                # Bias: inicializar com zeros (padrão recomendado)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        
        console.print(f"[green]✓ Pesos inicializados:[/green]")
        console.print(f"  • Camadas intermediárias: {activation_type} initialization")
        console.print(f"  • Pré-softmax + saída: Xavier initialization")
        console.print(f"  • Bias: zeros")


# ==============================================================================
# TRAINING AND EVALUATION
# ==============================================================================

class Trainer:
    """Classe para gerenciar treinamento e avaliação."""
    
    def __init__(
        self,
        model: nn.Module,
        train_loader: DataLoader,
        val_loader: DataLoader,
        config: Dict,
        device: torch.device,
        wandb_run: Optional[Any] = None
    ):
        """
        Inicializa trainer.
        
        Args:
            model: Modelo a treinar
            train_loader: DataLoader de treino
            val_loader: DataLoader de validação
            config: Configuração
            device: Device (CPU ou GPU)
            wandb_run: Run do W&B (opcional)
        """
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.config = config
        self.device = device
        self.wandb_run = wandb_run
        
        # Configurações de visualização/debug
        self.enable_visualization = config.get('debug', {}).get('enable_visualization', False)
        self.max_samples_per_epoch = config.get('debug', {}).get('max_samples_per_epoch', None)
        
        if self.enable_visualization:
            console.print(f"[yellow]⚠ Modo de visualização ativado! Batch size será forçado para 1.[/yellow]")
            console.print(f"[yellow]Pressione qualquer tecla na janela do gráfico para continuar.[/yellow]")
            plt.ion()  # Modo interativo
            self._key_pressed = False
        
        # Configurar otimizador
        optimizer_type = config['training']['optimizer'].lower()
        lr = config['training']['learning_rate']
        
        if optimizer_type == 'adam':
            self.optimizer = optim.Adam(model.parameters(), lr=lr)
        elif optimizer_type == 'adamw':
            self.optimizer = optim.AdamW(model.parameters(), lr=lr)
        elif optimizer_type == 'sgd':
            self.optimizer = optim.SGD(model.parameters(), lr=lr, momentum=0.9)
        else:
            raise ValueError(f"Otimizador não suportado: {optimizer_type}")
        
        # Configurar loss function
        loss_type = config['training']['loss_function']
        if loss_type == 'cross_entropy':
            self.criterion = nn.CrossEntropyLoss()
        elif loss_type == 'mse':
            self.criterion = nn.MSELoss()
        else:
            raise ValueError(f"Loss function não suportada: {loss_type}")
        
        # Histórico
        self.history = {
            'train_loss': [],
            'val_loss': [],
            'val_accuracy': [],
            'epoch': []
        }
        
        self.best_val_loss = float('inf')
        self.best_val_accuracy = 0.0
    
    def _on_key_press(self, event):
        """Callback para capturar tecla pressionada na janela do gráfico."""
        self._key_pressed = True
        plt.close()
    
    def _visualize_sample(self, features: torch.Tensor, targets: torch.Tensor, 
                         outputs: torch.Tensor, batch_idx: int, epoch: int):
        """
        Visualiza uma amostra de entrada e suas predições.
        
        Args:
            features: Tensor de entrada (1, num_features)
            targets: Target verdadeiro (1,)
            outputs: Saída da rede (1, num_classes)
            batch_idx: Índice do batch
            epoch: Número da época
        """
        # Converter para CPU e numpy
        features_np = features.cpu().detach().numpy().flatten()
        target_idx = targets.cpu().item()
        output_probs = torch.softmax(outputs, dim=1).cpu().detach().numpy()[0]
        predicted_idx = output_probs.argmax()
        
        # Obter nomes das classes (se disponível)
        class_names = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']  # Para superpopulation
        if target_idx < len(class_names):
            target_name = class_names[target_idx]
            predicted_name = class_names[predicted_idx]
        else:
            target_name = f"Class {target_idx}"
            predicted_name = f"Class {predicted_idx}"
        
        # Criar figura
        plt.clf()
        fig = plt.gcf()
        fig.set_size_inches(16, 8)
        
        # Plot das features de entrada
        plt.subplot(2, 1, 1)
        plt.plot(features_np, linewidth=0.5, alpha=0.7)
        plt.xlabel('Feature Index', fontsize=12)
        plt.ylabel('Feature Value', fontsize=12)
        plt.title(f'Época {epoch + 1} | Amostra {batch_idx + 1} | Input Features (n={len(features_np)})', 
                 fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        
        # Estatísticas das features
        stats_text = (f'Min: {features_np.min():.3f} | Max: {features_np.max():.3f} | '
                     f'Mean: {features_np.mean():.3f} | Std: {features_np.std():.3f}')
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Plot das probabilidades de saída
        plt.subplot(2, 1, 2)
        bars = plt.bar(range(len(output_probs)), output_probs, color='steelblue', alpha=0.7)
        bars[target_idx].set_color('green')
        bars[predicted_idx].set_edgecolor('red')
        bars[predicted_idx].set_linewidth(3)
        
        plt.xlabel('Class', fontsize=12)
        plt.ylabel('Probability', fontsize=12)
        plt.title('Network Output Probabilities', fontsize=14, fontweight='bold')
        plt.xticks(range(len(output_probs)), 
                  class_names[:len(output_probs)] if len(output_probs) <= len(class_names) 
                  else [str(i) for i in range(len(output_probs))])
        plt.grid(True, alpha=0.3, axis='y')
        plt.ylim([0, 1])
        
        # Texto com predição e target
        correct = "✓ CORRETO" if predicted_idx == target_idx else "✗ ERRADO"
        color = 'green' if predicted_idx == target_idx else 'red'
        result_text = (f'{correct}\n'
                      f'Target: {target_name} (classe {target_idx})\n'
                      f'Predito: {predicted_name} (classe {predicted_idx}, prob={output_probs[predicted_idx]:.3f})')
        
        plt.text(0.98, 0.98, result_text, transform=plt.gca().transAxes,
                fontsize=11, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
        
        # Adicionar instruções no canto superior direito da figura
        fig.text(0.98, 0.98, 'Pressione qualquer tecla para continuar...', 
                ha='right', va='top', fontsize=11, style='italic',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.6))
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)  # Dar espaço no topo
        
        # Conectar callback de tecla
        self._key_pressed = False
        cid = fig.canvas.mpl_connect('key_press_event', self._on_key_press)
        
        # Mostrar e aguardar tecla
        plt.show(block=True)
        
        # Desconectar callback
        fig.canvas.mpl_disconnect(cid)
    
    def train_epoch(self, epoch: int) -> float:
        """
        Treina por uma época.
        
        Returns:
            Loss média da época
        """
        self.model.train()
        total_loss = 0.0
        num_batches = 0
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task(
                f"[cyan]Época {epoch + 1} - Treino",
                total=len(self.train_loader)
            )
            
            for batch_idx, (features, targets) in enumerate(self.train_loader):
                # Verificar limite de amostras para visualização
                if self.enable_visualization and self.max_samples_per_epoch is not None:
                    if batch_idx >= self.max_samples_per_epoch:
                        break
                
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                
                # Forward pass
                self.optimizer.zero_grad()
                outputs = self.model(features)
                loss = self.criterion(outputs, targets)
                
                # Backward pass
                loss.backward()
                self.optimizer.step()
                
                total_loss += loss.item()
                num_batches += 1
                
                # Visualização interativa (se habilitada)
                if self.enable_visualization:
                    self._visualize_sample(features, targets, outputs, batch_idx, epoch)
                
                # Log no W&B
                if self.wandb_run and batch_idx % self.config['wandb']['log_frequency'] == 0:
                    self.wandb_run.log({
                        'batch_loss': loss.item(),
                        'epoch': epoch,
                        'batch': batch_idx
                    })
                
                progress.update(task, advance=1)
        
        avg_loss = total_loss / num_batches if num_batches > 0 else 0.0
        return avg_loss
    
    def validate(self, epoch: int) -> Tuple[float, float]:
        """
        Valida o modelo.
        
        Returns:
            Tupla (loss, accuracy)
        """
        self.model.eval()
        total_loss = 0.0
        all_predictions = []
        all_targets = []
        
        with torch.no_grad():
            for features, targets in self.val_loader:
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                
                outputs = self.model(features)
                loss = self.criterion(outputs, targets)
                
                total_loss += loss.item()
                
                # Para classificação, calcular accuracy
                if self.config['output']['prediction_target'] != 'frog_likelihood':
                    predictions = torch.argmax(outputs, dim=1)
                    all_predictions.extend(predictions.cpu().numpy())
                    all_targets.extend(targets.cpu().numpy())
        
        avg_loss = total_loss / len(self.val_loader)
        
        # Calcular accuracy
        if len(all_predictions) > 0:
            accuracy = accuracy_score(all_targets, all_predictions)
        else:
            accuracy = 0.0
        
        console.print(f"[green]Validação - Época {epoch + 1}:[/green] Loss={avg_loss:.4f}, Accuracy={accuracy:.4f}")
        
        return avg_loss, accuracy
    
    def train(self) -> Dict:
        """
        Executa loop de treinamento completo.
        
        Returns:
            Histórico de treinamento
        """
        num_epochs = self.config['training']['num_epochs']
        val_frequency = self.config['training']['validation_frequency']
        save_frequency = self.config['checkpointing']['save_frequency']
        
        console.print(Panel.fit(
            f"[bold cyan]Iniciando Treinamento[/bold cyan]\n"
            f"Épocas: {num_epochs}\n"
            f"Batch size: {self.config['training']['batch_size']}\n"
            f"Learning rate: {self.config['training']['learning_rate']}"
        ))
        
        for epoch in range(num_epochs):
            # Treinar
            train_loss = self.train_epoch(epoch)
            
            # Validar
            if (epoch + 1) % val_frequency == 0:
                val_loss, val_accuracy = self.validate(epoch)
                
                # Salvar histórico
                self.history['train_loss'].append(train_loss)
                self.history['val_loss'].append(val_loss)
                self.history['val_accuracy'].append(val_accuracy)
                self.history['epoch'].append(epoch + 1)
                
                # Log no W&B
                if self.wandb_run:
                    self.wandb_run.log({
                        'epoch': epoch + 1,
                        'train_loss': train_loss,
                        'val_loss': val_loss,
                        'val_accuracy': val_accuracy
                    })
                
                # Salvar melhor modelo
                if val_loss < self.best_val_loss:
                    self.best_val_loss = val_loss
                    self.save_checkpoint(epoch, 'best_loss')
                
                if val_accuracy > self.best_val_accuracy:
                    self.best_val_accuracy = val_accuracy
                    self.save_checkpoint(epoch, 'best_accuracy')
            
            # Salvar checkpoint periódico
            if (epoch + 1) % save_frequency == 0:
                self.save_checkpoint(epoch, f'epoch_{epoch + 1}')
        
        console.print("[bold green]✓ Treinamento concluído![/bold green]")
        return self.history
    
    def save_checkpoint(self, epoch: int, name: str):
        """Salva checkpoint do modelo."""
        checkpoint_dir = Path(self.config['checkpointing']['checkpoint_dir'])
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        
        checkpoint_path = checkpoint_dir / f"{name}.pt"
        
        torch.save({
            'epoch': epoch,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'best_val_loss': self.best_val_loss,
            'best_val_accuracy': self.best_val_accuracy,
            'history': self.history,
            'config': self.config
        }, checkpoint_path)
        
        console.print(f"[green]✓ Checkpoint salvo: {checkpoint_path}[/green]")


class Tester:
    """Classe para teste e avaliação."""
    
    def __init__(
        self,
        model: nn.Module,
        test_loader: DataLoader,
        dataset: ProcessedGenomicDataset,
        config: Dict,
        device: torch.device,
        wandb_run: Optional[Any] = None
    ):
        """Inicializa tester."""
        self.model = model
        self.test_loader = test_loader
        self.dataset = dataset
        self.config = config
        self.device = device
        self.wandb_run = wandb_run
    
    def test(self) -> Dict:
        """
        Executa teste e gera métricas.
        
        Returns:
            Dict com resultados
        """
        console.print(Panel.fit("[bold cyan]Executando Teste[/bold cyan]"))
        
        self.model.eval()
        all_predictions = []
        all_targets = []
        all_probs = []
        
        with torch.no_grad():
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                console=console
            ) as progress:
                task = progress.add_task(
                    "[cyan]Testando modelo...",
                    total=len(self.test_loader)
                )
                
                for features, targets in self.test_loader:
                    features = features.to(self.device, non_blocking=True)
                    targets = targets.to(self.device, non_blocking=True)
                    
                    outputs = self.model(features)
                    
                    if self.config['output']['prediction_target'] != 'frog_likelihood':
                        predictions = torch.argmax(outputs, dim=1)
                        all_predictions.extend(predictions.cpu().numpy())
                        all_targets.extend(targets.cpu().numpy())
                        all_probs.append(outputs.cpu().numpy())
                    
                    progress.update(task, advance=1)
        
        # Calcular métricas
        results = {}
        
        if len(all_predictions) > 0:
            # Preparar labels e nomes para todas as classes
            # (importante para incluir classes que podem não aparecer no conjunto de teste)
            labels = list(range(self.dataset.get_num_classes()))
            target_names = [self.dataset.idx_to_target[i] for i in labels]
            
            # Calcular métricas
            results['accuracy'] = accuracy_score(all_targets, all_predictions)
            results['precision'], results['recall'], results['f1'], _ = \
                precision_recall_fscore_support(all_targets, all_predictions, average='weighted', zero_division=0)
            results['confusion_matrix'] = confusion_matrix(all_targets, all_predictions, labels=labels)
            results['classification_report'] = classification_report(
                all_targets, all_predictions, 
                labels=labels,
                target_names=target_names, 
                zero_division=0
            )
            
            # Imprimir resultados
            self._print_results(results, target_names)
            
            # Log no W&B
            if self.wandb_run:
                self._log_wandb(results, all_targets, all_predictions, target_names)
        
        return results
    
    def _print_results(self, results: Dict, target_names: List[str]):
        """Imprime resultados formatados."""
        console.print("\n[bold cyan]═══ RESULTADOS DO TESTE ═══[/bold cyan]\n")
        
        # Métricas gerais
        table = Table(title="Métricas de Performance")
        table.add_column("Métrica", style="cyan")
        table.add_column("Valor", style="green")
        
        table.add_row("Accuracy", f"{results['accuracy']:.4f}")
        table.add_row("Precision (weighted)", f"{results['precision']:.4f}")
        table.add_row("Recall (weighted)", f"{results['recall']:.4f}")
        table.add_row("F1-Score (weighted)", f"{results['f1']:.4f}")
        
        console.print(table)
        
        # Classification report
        console.print("\n[bold]Classification Report:[/bold]")
        console.print(results['classification_report'])
        
        # Confusion Matrix
        console.print("\n[bold]Confusion Matrix:[/bold]")
        cm = results['confusion_matrix']
        
        cm_table = Table(title="Confusion Matrix")
        cm_table.add_column("True \\ Pred", style="cyan")
        for name in target_names:
            cm_table.add_column(name[:10], justify="right")
        
        for i, name in enumerate(target_names):
            row = [name[:10]] + [str(cm[i, j]) for j in range(len(target_names))]
            cm_table.add_row(*row)
        
        console.print(cm_table)
    
    def _log_wandb(self, results: Dict, all_targets: List, all_predictions: List, target_names: List[str]):
        """Loga resultados no W&B."""
        import wandb
        
        self.wandb_run.log({
            'test_accuracy': results['accuracy'],
            'test_precision': results['precision'],
            'test_recall': results['recall'],
            'test_f1': results['f1']
        })
        
        # Confusion matrix
        self.wandb_run.log({
            'confusion_matrix': wandb.plot.confusion_matrix(
                probs=None,
                y_true=all_targets,
                preds=all_predictions,
                class_names=target_names
            )
        })


# ==============================================================================
# DATASET CACHE FUNCTIONS
# ==============================================================================

def validate_cache(cache_dir: Path, config: Dict) -> bool:
    """
    Valida se cache existe e é compatível com configuração atual.
    
    Args:
        cache_dir: Diretório do cache
        config: Configuração atual
        
    Returns:
        True se cache é válido, False caso contrário
    """
    cache_dir = Path(cache_dir)
    
    # Verificar se diretório existe
    if not cache_dir.exists():
        return False
    
    # Verificar arquivo de flag que indica cache completo
    complete_flag = cache_dir / '.cache_complete'
    if not complete_flag.exists():
        console.print(f"[yellow]Cache incompleto: processo foi interrompido antes de concluir[/yellow]")
        return False
    
    # Verificar arquivos necessários
    required_files = [
        'metadata.json',
        'normalization_params.json',
        'splits.json',
        'train_data.pt',
        'val_data.pt',
        'test_data.pt'
    ]
    
    for filename in required_files:
        if not (cache_dir / filename).exists():
            console.print(f"[yellow]Cache incompleto: falta {filename}[/yellow]")
            return False
    
    # Carregar e verificar metadados
    try:
        with open(cache_dir / 'metadata.json', 'r') as f:
            metadata = json.load(f)
        
        # Verificar compatibilidade dos parâmetros críticos
        processing_params = metadata.get('processing_params', {})
        current_params = {
            'alphagenome_outputs': config['dataset_input']['alphagenome_outputs'],
            'haplotype_mode': config['dataset_input']['haplotype_mode'],
            'window_center_size': config['dataset_input']['window_center_size'],
            'downsample_factor': config['dataset_input']['downsample_factor'],
            'normalization_method': config['dataset_input'].get('normalization_method', 'zscore'),
            'dataset_dir': config['dataset_input']['dataset_dir'],
            'taint_at_cache_save': config.get('debug', {}).get('taint_at_cache_save', False),
        }
        
        for key, current_value in current_params.items():
            if key == 'dataset_dir':
                # Dataset dir deve existir
                if not Path(current_value).exists():
                    console.print(f"[yellow]Dataset dir não existe: {current_value}[/yellow]")
                    return False
                # Comparar path absoluto
                cached_dir = Path(metadata.get('dataset_dir', ''))
                if cached_dir.resolve() != Path(current_value).resolve():
                    console.print(f"[yellow]Dataset dir diferente: cache={metadata.get('dataset_dir')}, atual={current_value}[/yellow]")
                    return False
            else:
                cached_value = processing_params.get(key)
                if cached_value != current_value:
                    console.print(f"[yellow]Parâmetro {key} diferente: cache={cached_value}, atual={current_value}[/yellow]")
                    return False
        
        # Verificar splits
        splits_metadata = metadata.get('splits', {})
        current_splits = {
            'train_split': config['data_split']['train_split'],
            'val_split': config['data_split']['val_split'],
            'test_split': config['data_split']['test_split'],
            'random_seed': config['data_split']['random_seed']
        }
        
        if splits_metadata.get('random_seed') != current_splits['random_seed']:
            console.print(f"[yellow]Random seed diferente[/yellow]")
            return False
        
        console.print("[green]✓ Cache válido e compatível[/green]")
        return True
        
    except Exception as e:
        console.print(f"[yellow]Erro ao validar cache: {e}[/yellow]")
        return False


def save_processed_dataset(
    cache_dir: Path,
    processed_dataset: ProcessedGenomicDataset,
    train_indices: List[int],
    val_indices: List[int],
    test_indices: List[int],
    config: Dict
):
    """
    Salva dataset processado em cache de forma idempotente.
    
    Args:
        cache_dir: Diretório onde salvar cache
        processed_dataset: Dataset processado
        train_indices: Índices de treino
        val_indices: Índices de validação
        test_indices: Índices de teste
        config: Configuração usada
    """
    cache_dir = Path(cache_dir)
    
    # Criar diretório temporário para escrita atômica
    temp_cache_dir = cache_dir.parent / f"{cache_dir.name}_tmp_{os.getpid()}"
    
    # Limpar temp dir se existir (de alguma execução anterior interrompida)
    if temp_cache_dir.exists():
        console.print(f"[yellow]Limpando diretório temporário de execução anterior...[/yellow]")
        shutil.rmtree(temp_cache_dir)
    
    temp_cache_dir.mkdir(parents=True, exist_ok=True)
    
    console.print(f"\n[bold cyan]💾 Salvando Dataset Processado em Cache[/bold cyan]")
    console.print(f"  📁 Diretório: {cache_dir}")
    console.print(f"  ⚙️  Escrevendo em diretório temporário primeiro...")
    console.print(f"  📊 Amostras de treino: {len(train_indices)}")
    console.print(f"  📊 Amostras de validação: {len(val_indices)}")
    console.print(f"  📊 Amostras de teste: {len(test_indices)}")
    
    # Preparar dados de cada split
    train_data = []
    val_data = []
    test_data = []
    
    # Verificar se tainting está habilitado ao salvar cache
    taint_at_save = config.get('debug', {}).get('taint_at_cache_save', False)
    num_classes = processed_dataset.get_num_classes() if taint_at_save else 0
    
    # IMPORTANTE: Desabilitar temporariamente taint_at_runtime enquanto salvamos o cache
    # para evitar que dados sejam salvos com tainting não-intencional
    original_taint_runtime = processed_dataset.config.get('debug', {}).get('taint_at_runtime', False)
    if 'debug' not in processed_dataset.config:
        processed_dataset.config['debug'] = {}
    processed_dataset.config['debug']['taint_at_runtime'] = False
    
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Processando e salvando dados...", total=len(processed_dataset))
            
            for idx in range(len(processed_dataset)):
                features, target = processed_dataset[idx]
                
                # Aplicar tainting se habilitado (apenas para classificação)
                if taint_at_save and config['output']['prediction_target'] != 'frog_likelihood':
                    # Extrair classe do target tensor
                    if target.ndim == 0:  # Escalar
                        target_class = target.item()
                    else:
                        target_class = target[0].item()
                    
                    # Aplicar tainting
                    features = taint_sample(features, target_class, num_classes)
                
                if idx in train_indices:
                    train_data.append((features, target))
                elif idx in val_indices:
                    val_data.append((features, target))
                elif idx in test_indices:
                    test_data.append((features, target))
                
                progress.update(task, advance=1)
        
        # Salvar dados no diretório temporário
        console.print(f"  💾 Salvando train_data.pt ({len(train_data)} amostras)...")
        torch.save(train_data, temp_cache_dir / 'train_data.pt')
        
        console.print(f"  💾 Salvando val_data.pt ({len(val_data)} amostras)...")
        torch.save(val_data, temp_cache_dir / 'val_data.pt')
        
        console.print(f"  💾 Salvando test_data.pt ({len(test_data)} amostras)...")
        torch.save(test_data, temp_cache_dir / 'test_data.pt')
        
        # Salvar splits
        console.print(f"  💾 Salvando splits.json...")
        splits = {
            'train_indices': train_indices,
            'val_indices': val_indices,
            'test_indices': test_indices
        }
        with open(temp_cache_dir / 'splits.json', 'w') as f:
            json.dump(splits, f, indent=2)
        
        # Salvar normalization params
        console.print(f"  💾 Salvando normalization_params.json...")
        with open(temp_cache_dir / 'normalization_params.json', 'w') as f:
            json.dump(processed_dataset.normalization_params, f, indent=2)
        
        # Salvar metadados
        metadata = {
            'creation_date': datetime.now().isoformat(),
            'dataset_dir': config['dataset_input']['dataset_dir'],
            'processing_params': {
                'alphagenome_outputs': config['dataset_input']['alphagenome_outputs'],
                'haplotype_mode': config['dataset_input']['haplotype_mode'],
                'window_center_size': config['dataset_input']['window_center_size'],
                'downsample_factor': config['dataset_input']['downsample_factor'],
                'normalization_method': config['dataset_input'].get('normalization_method', 'zscore'),
                'taint_at_cache_save': config.get('debug', {}).get('taint_at_cache_save', False),
            },
            'splits': {
                'train_size': len(train_indices),
                'val_size': len(val_indices),
                'test_size': len(test_indices),
                'train_split': config['data_split']['train_split'],
                'val_split': config['data_split']['val_split'],
                'test_split': config['data_split']['test_split'],
                'random_seed': config['data_split']['random_seed']
            },
            'total_samples': len(processed_dataset),
            'num_classes': processed_dataset.get_num_classes(),
            'input_size': processed_dataset.get_input_size(),
            'prediction_target': config['output']['prediction_target']
        }
        console.print(f"  💾 Salvando metadata.json...")
        with open(temp_cache_dir / 'metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Criar arquivo de flag indicando cache completo
        console.print(f"  ✓ Marcando cache como completo...")
        (temp_cache_dir / '.cache_complete').touch()
        
        # Mover atomicamente do temp para o diretório final
        console.print(f"  🔄 Movendo cache para localização final...")
        
        # Se cache_dir já existir (cache antigo), remover
        if cache_dir.exists():
            console.print(f"  🗑️  Removendo cache antigo...")
            shutil.rmtree(cache_dir)
        
        # Renomear temp_cache_dir para cache_dir (operação atômica no mesmo filesystem)
        temp_cache_dir.rename(cache_dir)
        
        console.print(f"\n[bold green]✓ Cache Salvo com Sucesso![/bold green]")
        console.print(f"  📁 Localização: {cache_dir}")
        console.print(f"  📊 Train: {len(train_data)} amostras")
        console.print(f"  📊 Val: {len(val_data)} amostras")
        console.print(f"  📊 Test: {len(test_data)} amostras")
        console.print(f"  💡 Próximas execuções usarão este cache automaticamente!")
    
    finally:
        # Restaurar o valor original de taint_at_runtime
        processed_dataset.config['debug']['taint_at_runtime'] = original_taint_runtime


class CachedProcessedDataset(Dataset):
    """
    Dataset wrapper para dados carregados do cache.
    """
    
    def __init__(self, data_file: Path, target_to_idx: Dict, idx_to_target: Dict, config: Dict):
        """
        Inicializa dataset do cache.
        
        Args:
            data_file: Arquivo .pt com dados processados
            target_to_idx: Mapeamento target->índice
            idx_to_target: Mapeamento índice->target
            config: Configuração (necessário para taint_at_runtime)
        """
        console.print(f"[cyan]Carregando {data_file.name}...[/cyan]")
        self.data = torch.load(data_file)
        self.target_to_idx = target_to_idx
        self.idx_to_target = idx_to_target
        self.config = config
        console.print(f"[green]✓ {len(self.data)} samples carregados[/green]")
    
    def __len__(self) -> int:
        return len(self.data)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        features, target = self.data[idx]
        
        # Aplicar tainting em runtime (se habilitado)
        # Apenas para classificação (não para regressão/frog_likelihood)
        if (self.config.get('debug', {}).get('taint_at_runtime', False) and
            self.config['output']['prediction_target'] != 'frog_likelihood'):
            num_classes = len(self.target_to_idx)
            # Extrair classe do target tensor
            if target.ndim == 0:  # Escalar
                target_class = target.item()
            else:
                target_class = target[0].item()
            # Aplicar tainting
            features = taint_sample(features, target_class, num_classes)
        
        return features, target
    
    def get_num_classes(self) -> int:
        return len(self.target_to_idx)
    
    def get_input_size(self) -> int:
        if len(self.data) > 0:
            return self.data[0][0].shape[0]
        return 0


def load_processed_dataset(
    cache_dir: Path,
    config: Dict
) -> Tuple[CachedProcessedDataset, DataLoader, DataLoader, DataLoader]:
    """
    Carrega dataset processado do cache.
    
    Args:
        cache_dir: Diretório do cache
        config: Configuração atual
        
    Returns:
        Tupla (full_dataset, train_loader, val_loader, test_loader)
    """
    cache_dir = Path(cache_dir)
    
    console.print(Panel.fit(f"[bold cyan]Carregando Dataset do Cache[/bold cyan]\n{cache_dir}"))
    
    # Carregar metadados
    with open(cache_dir / 'metadata.json', 'r') as f:
        metadata = json.load(f)
    
    console.print(f"[green]Cache criado em: {metadata['creation_date']}[/green]")
    console.print(f"[green]Total de amostras: {metadata['total_samples']}[/green]")
    
    # Carregar normalization params (para referência)
    with open(cache_dir / 'normalization_params.json', 'r') as f:
        norm_params = json.load(f)
    
    # Mostrar parâmetros de normalização de acordo com o método
    method = norm_params.get('method', 'zscore')
    if method == 'zscore':
        console.print(f"[green]Normalização: zscore (mean={norm_params.get('mean', 0):.6f}, std={norm_params.get('std', 1):.6f})[/green]")
    elif method == 'minmax_keep_zero':
        console.print(f"[green]Normalização: minmax_keep_zero (max={norm_params.get('max', 1):.6f})[/green]")
    elif method == 'log':
        console.print(f"[green]Normalização: log (log_max={norm_params.get('log_max', 1):.6f})[/green]")
    else:
        console.print(f"[green]Normalização: {method}[/green]")
    
    # Criar mapeamentos de target (do metadata)
    # Para cached dataset, vamos reconstruir os mapeamentos baseados no prediction_target
    # Isso é simplificado - assumimos que os targets já estão convertidos para índices
    target_to_idx = {}
    idx_to_target = {}
    
    # Para classificação, criar mapeamentos dummy (os dados já estão como índices)
    num_classes = metadata.get('num_classes', 0)
    for i in range(num_classes):
        target_to_idx[str(i)] = i
        idx_to_target[i] = str(i)
    
    # Carregar datasets
    train_dataset = CachedProcessedDataset(cache_dir / 'train_data.pt', target_to_idx, idx_to_target, config)
    val_dataset = CachedProcessedDataset(cache_dir / 'val_data.pt', target_to_idx, idx_to_target, config)
    test_dataset = CachedProcessedDataset(cache_dir / 'test_data.pt', target_to_idx, idx_to_target, config)
    
    # Preparar para criar DataLoaders
    console.print("\n[cyan]⚙️  Criando DataLoaders...[/cyan]")
    batch_size = config['training']['batch_size']
    
    # Forçar batch_size=1 se visualização estiver habilitada
    if config.get('debug', {}).get('enable_visualization', False):
        batch_size = 1
        console.print(f"  • [yellow]Batch size forçado para 1 (visualização habilitada)[/yellow]")
    
    console.print(f"  • Inicializando workers paralelos (train: 4 workers, val/test: 2 workers)")
    console.print(f"  • Batch size: {batch_size}")
    console.print(f"  • Isso pode levar alguns segundos na primeira vez...")
    
    # Função collate para empilhar batches corretamente
    def collate_fn(batch):
        """Empilha batch de tuplas (features, target) em tensors."""
        # Como os dados vêm do cache já processados, apenas empilhar
        features_list, targets_list = zip(*batch)
        
        # Empilhar features: (batch_size, num_features)
        features_batch = torch.stack(features_list, dim=0)
        
        # Empilhar targets: (batch_size,)
        targets_batch = torch.stack(targets_list, dim=0)
        
        return features_batch, targets_batch
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=4,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=True
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=2,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=True
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=2,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=True
    )
    
    console.print(f"[green]✓ DataLoaders criados com sucesso![/green]\n")
    console.print(f"[green]Dataset splits carregados:[/green]")
    console.print(f"  • Treino: {len(train_dataset)} amostras")
    console.print(f"  • Validação: {len(val_dataset)} amostras")
    console.print(f"  • Teste: {len(test_dataset)} amostras")
    
    # Usar train_dataset como full_dataset (para compatibilidade com código existente)
    return train_dataset, train_loader, val_loader, test_loader


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

def load_config(config_path: Path) -> Dict:
    """Carrega configuração do YAML."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def save_config(config: Dict, output_path: Path):
    """Salva configuração em JSON."""
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)


def prepare_data(config: Dict) -> Tuple[Any, DataLoader, DataLoader, DataLoader]:
    """
    Prepara datasets e dataloaders.
    Tenta carregar do cache se disponível, senão processa e salva.
    
    Returns:
        Tupla (full_dataset, train_loader, val_loader, test_loader)
    """
    # Verificar se cache está configurado
    cache_dir = config['dataset_input'].get('processed_cache_dir')
    
    # Tentar carregar do cache
    if cache_dir is not None:
        cache_path = Path(cache_dir)
        if cache_path.exists() and validate_cache(cache_path, config):
            console.print(Panel.fit(
                "[bold green]Cache Encontrado![/bold green]\n"
                "Carregando dataset processado do cache para economizar tempo..."
            ))
            return load_processed_dataset(cache_path, config)
        elif cache_path.exists():
            console.print(Panel.fit(
                "[bold yellow]Cache Inválido ou Incompatível[/bold yellow]\n"
                "Parâmetros mudaram. Reprocessando dataset..."
            ))
        else:
            console.print(Panel.fit(
                "[bold cyan]Cache Não Encontrado[/bold cyan]\n"
                "Primeira execução. Processando e salvando dataset..."
            ))
    else:
        console.print(Panel.fit("[bold cyan]Cache Desabilitado[/bold cyan]\nProcessando dataset..."))
    
    # Processar dataset do zero
    console.print("[bold cyan]Carregando Dataset Base[/bold cyan]")
    
    # Carregar dataset base
    base_dataset = GenomicLongevityDataset(
        dataset_dir=Path(config['dataset_input']['dataset_dir']),
        load_predictions=True,
        load_sequences=False,
        cache_sequences=False
    )
    
    # Tentar carregar parâmetros de normalização do cache (se existirem e forem compatíveis)
    normalization_params = None
    if cache_dir is not None:
        cache_path = Path(cache_dir)
        norm_file = cache_path / 'normalization_params.json'
        metadata_file = cache_path / 'metadata.json'
        
        if norm_file.exists() and metadata_file.exists():
            try:
                # Verificar se parâmetros de processamento são compatíveis
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                processing_params = metadata.get('processing_params', {})
                current_params = {
                    'alphagenome_outputs': config['dataset_input']['alphagenome_outputs'],
                    'haplotype_mode': config['dataset_input']['haplotype_mode'],
                    'window_center_size': config['dataset_input']['window_center_size'],
                    'downsample_factor': config['dataset_input']['downsample_factor'],
                    'normalization_method': config['dataset_input'].get('normalization_method', 'zscore'),
                }
                
                params_match = all(
                    processing_params.get(k) == v 
                    for k, v in current_params.items()
                )
                
                if params_match:
                    with open(norm_file, 'r') as f:
                        normalization_params = json.load(f)
                    console.print(f"[green]✓ Parâmetros de normalização carregados do cache[/green]")
                    
                    # Mostrar parâmetros de acordo com o método
                    method = normalization_params.get('method', 'zscore')
                    if method == 'zscore':
                        console.print(f"  • Método: Z-score")
                        console.print(f"  • Média: {normalization_params.get('mean', 0):.6f}")
                        console.print(f"  • Desvio padrão: {normalization_params.get('std', 1):.6f}")
                    elif method == 'minmax_keep_zero':
                        console.print(f"  • Método: MinMax (mantém zeros)")
                        console.print(f"  • Máximo não-zero: {normalization_params.get('max', 1):.6f}")
                    elif method == 'log':
                        console.print(f"  • Método: Logarítmico")
                        console.print(f"  • log1p(max): {normalization_params.get('log_max', 1):.6f}")
                    else:
                        console.print(f"  • Método: {method}")
            except Exception as e:
                console.print(f"[yellow]⚠ Não foi possível carregar parâmetros de normalização: {e}[/yellow]")
                normalization_params = None
    
    # Criar dataset processado
    processed_dataset = ProcessedGenomicDataset(
        base_dataset=base_dataset,
        config=config,
        normalization_params=normalization_params,
        compute_normalization=(normalization_params is None)
    )
    
    # Salvar parâmetros de normalização no diretório de checkpoints (para referência)
    norm_path = Path(config['checkpointing']['checkpoint_dir']) / 'normalization_params.json'
    norm_path.parent.mkdir(parents=True, exist_ok=True)
    with open(norm_path, 'w') as f:
        json.dump(processed_dataset.normalization_params, f, indent=2)
    console.print(f"[green]✓ Parâmetros de normalização salvos em {norm_path}[/green]")
    
    # Também salvar no cache_dir para reutilização (mesmo que cache completo não exista ainda)
    if cache_dir is not None and normalization_params is None:  # Só se computamos agora
        cache_path = Path(cache_dir)
        cache_path.mkdir(parents=True, exist_ok=True)
        
        # Salvar metadados parciais para validação futura
        partial_metadata = {
            'creation_date': datetime.now().isoformat(),
            'processing_params': {
                'alphagenome_outputs': config['dataset_input']['alphagenome_outputs'],
                'haplotype_mode': config['dataset_input']['haplotype_mode'],
                'window_center_size': config['dataset_input']['window_center_size'],
                'downsample_factor': config['dataset_input']['downsample_factor'],
            }
        }
        with open(cache_path / 'metadata.json', 'w') as f:
            json.dump(partial_metadata, f, indent=2)
        
        with open(cache_path / 'normalization_params.json', 'w') as f:
            json.dump(processed_dataset.normalization_params, f, indent=2)
        
        console.print(f"[green]✓ Parâmetros de normalização também salvos no cache para reutilização[/green]")
    
    # Split dataset
    total_size = len(processed_dataset)
    train_size = int(config['data_split']['train_split'] * total_size)
    val_size = int(config['data_split']['val_split'] * total_size)
    test_size = total_size - train_size - val_size
    
    # Criar índices
    indices = list(range(total_size))
    if config['data_split']['random_seed'] is not None:
        np.random.seed(config['data_split']['random_seed'])
        np.random.shuffle(indices)
    
    train_indices = indices[:train_size]
    val_indices = indices[train_size:train_size + val_size]
    test_indices = indices[train_size + val_size:]
    
    console.print(f"[green]Dataset split:[/green]")
    console.print(f"  • Treino: {len(train_indices)} amostras")
    console.print(f"  • Validação: {len(val_indices)} amostras")
    console.print(f"  • Teste: {len(test_indices)} amostras")
    
    # Salvar cache se configurado
    if cache_dir is not None:
        cache_path = Path(cache_dir)
        save_processed_dataset(
            cache_path,
            processed_dataset,
            train_indices,
            val_indices,
            test_indices,
            config
        )
        
        # IMPORTANTE: Agora que o cache foi salvo, carregar dele para usar
        # o CachedProcessedDataset (rápido) em vez do ProcessedGenomicDataset (lento)
        console.print("\n[cyan]✓ Cache salvo! Recarregando do cache para treino rápido...[/cyan]")
        return load_processed_dataset(cache_path, config)
    
    # Se cache não está configurado, criar subsets do ProcessedGenomicDataset
    # (será lento porque processa on-the-fly)
    console.print("\n[yellow]⚠ Cache desabilitado: treino será mais lento (processamento on-the-fly)[/yellow]")
    train_dataset = Subset(processed_dataset, train_indices)
    val_dataset = Subset(processed_dataset, val_indices)
    test_dataset = Subset(processed_dataset, test_indices)
    
    # Preparar para criar DataLoaders
    console.print("\n[cyan]⚙️  Criando DataLoaders...[/cyan]")
    batch_size = config['training']['batch_size']
    
    # Forçar batch_size=1 se visualização estiver habilitada
    if config.get('debug', {}).get('enable_visualization', False):
        batch_size = 1
        console.print(f"  • [yellow]Batch size forçado para 1 (visualização habilitada)[/yellow]")
    
    console.print(f"  • Inicializando workers paralelos (train: 4 workers, val/test: 2 workers)")
    console.print(f"  • Batch size: {batch_size}")
    console.print(f"  • Isso pode levar alguns segundos na primeira vez...")
    
    # Função collate para empilhar batches corretamente
    def collate_fn(batch):
        """Empilha batch de tuplas (features, target) em tensors."""
        # Como os dados vêm do cache já processados, apenas empilhar
        features_list, targets_list = zip(*batch)
        
        # Empilhar features: (batch_size, num_features)
        features_batch = torch.stack(features_list, dim=0)
        
        # Empilhar targets: (batch_size,)
        targets_batch = torch.stack(targets_list, dim=0)
        
        return features_batch, targets_batch
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=4,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=True
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=2,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=True
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=2,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=True
    )
    
    return processed_dataset, train_loader, val_loader, test_loader


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description="Neural Ancestry Predictor - Predição de ancestralidade a partir de dados AlphaGenome"
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Caminho para arquivo de configuração YAML'
    )
    parser.add_argument(
        '--mode',
        type=str,
        choices=['train', 'test'],
        help='Modo de operação (sobrescreve config)'
    )
    
    args = parser.parse_args()
    
    # Carregar configuração
    config = load_config(Path(args.config))
    
    # Sobrescrever modo se fornecido
    if args.mode:
        config['mode'] = args.mode
    
    # Banner
    console.print(Panel.fit(
        "[bold cyan]Neural Ancestry Predictor[/bold cyan]\n"
        f"Modo: {config['mode']}\n"
        f"Target: {config['output']['prediction_target']}\n"
        f"Config: {args.config}",
        title="🧬 Genomics"
    ))
    
    # Configurar device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    console.print(f"[green]Device: {device}[/green]")
    
    # Preparar dados
    full_dataset, train_loader, val_loader, test_loader = prepare_data(config)
    
    # Calcular tamanhos
    input_size = full_dataset.get_input_size()
    num_classes = full_dataset.get_num_classes() if config['output']['prediction_target'] != 'frog_likelihood' else 150
    
    console.print(f"[green]Input size: {input_size}[/green]")
    console.print(f"[green]Number of classes: {num_classes}[/green]")
    
    # Criar modelo
    model = AncestryPredictor(config, input_size, num_classes).to(device)
    
    # Inicializar W&B
    wandb_run = None
    if config['wandb']['use_wandb']:
        try:
            import wandb
            wandb_run = wandb.init(
                project=config['wandb']['project_name'],
                name=config['wandb'].get('run_name'),
                config=config
            )
            console.print("[green]✓ Weights & Biases inicializado[/green]")
        except ImportError:
            console.print("[yellow]⚠ Weights & Biases não disponível. Instale com: pip install wandb[/yellow]")
        except Exception as e:
            console.print(f"[yellow]⚠ Erro ao inicializar W&B: {e}[/yellow]")
    
    # Modo de operação
    if config['mode'] == 'train':
        # Carregar checkpoint se fornecido
        if config['checkpointing'].get('load_checkpoint'):
            checkpoint_path = Path(config['checkpointing']['load_checkpoint'])
            if checkpoint_path.exists():
                console.print(f"[yellow]Carregando checkpoint: {checkpoint_path}[/yellow]")
                checkpoint = torch.load(checkpoint_path, map_location=device)
                model.load_state_dict(checkpoint['model_state_dict'])
        
        # Treinar
        trainer = Trainer(model, train_loader, val_loader, config, device, wandb_run)
        history = trainer.train()
        
        # Salvar histórico
        history_path = Path(config['checkpointing']['checkpoint_dir']) / 'training_history.json'
        with open(history_path, 'w') as f:
            json.dump(history, f, indent=2)
        console.print(f"[green]✓ Histórico salvo em {history_path}[/green]")
        
    elif config['mode'] == 'test':
        # Carregar melhor checkpoint
        checkpoint_path = Path(config['checkpointing']['checkpoint_dir']) / 'best_accuracy.pt'
        
        if not checkpoint_path.exists():
            console.print(f"[red]ERRO: Checkpoint não encontrado: {checkpoint_path}[/red]")
            sys.exit(1)
        
        console.print(f"[yellow]Carregando checkpoint: {checkpoint_path}[/yellow]")
        checkpoint = torch.load(checkpoint_path, map_location=device)
        model.load_state_dict(checkpoint['model_state_dict'])
        
        # Testar
        tester = Tester(model, test_loader, full_dataset, config, device, wandb_run)
        results = tester.test()
    
    # Finalizar W&B
    if wandb_run:
        wandb_run.finish()
    
    console.print("\n[bold green]✓ Execução concluída![/bold green]")


if __name__ == '__main__':
    main()

