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
import signal
import io
import contextlib
import time
import psutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings
from datetime import datetime

import numpy as np
from scipy import ndimage
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

# Global flag para capturar CTRL+C
class InterruptState:
    """Estado global para interrupção de treinamento."""
    interrupted = False

interrupt_state = InterruptState()

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


def block_reduce_2d(arr: np.ndarray, target_shape: Tuple[int, int], 
                    func: str = 'max') -> np.ndarray:
    """
    Redimensiona array 2D para target_shape usando agregação por bloco para downscale
    e interpolação para upscale, tratando altura e largura independentemente.
    
    Isso permite preservar picos (ex: taint) quando há downscale em uma dimensão
    e upscale em outra.
    
    Args:
        arr: Array 2D numpy de entrada
        target_shape: Tupla (altura, largura) do tamanho desejado
        func: Função de agregação para downscale - 'max', 'min' ou 'mean'
        
    Returns:
        Array 2D redimensionado para target_shape
    """
    h, w = arr.shape
    th, tw = target_shape
    
    # Selecionar função de agregação
    if func == 'max':
        agg_func = np.max
    elif func == 'min':
        agg_func = np.min
    elif func == 'mean':
        agg_func = np.mean
    else:
        raise ValueError(f"Função de agregação desconhecida: {func}. Use 'max', 'min' ou 'mean'.")
    
    # Passo 1: Processar LARGURA (colunas)
    if tw >= w:
        # Upscale na largura - usar interpolação
        temp = ndimage.zoom(arr, (1, tw / w), order=1)
    else:
        # Downscale na largura - usar block_reduce
        bw = w / tw
        temp = np.zeros((h, tw), dtype=arr.dtype)
        for j in range(tw):
            x_start = int(j * bw)
            x_end = int((j + 1) * bw)
            x_end = min(x_end, w)
            if x_end > x_start:
                temp[:, j] = agg_func(arr[:, x_start:x_end], axis=1)
    
    # Passo 2: Processar ALTURA (linhas)
    h_temp = temp.shape[0]
    if th >= h_temp:
        # Upscale na altura - usar interpolação
        result = ndimage.zoom(temp, (th / h_temp, 1), order=1)
    else:
        # Downscale na altura - usar block_reduce
        bh = h_temp / th
        result = np.zeros((th, tw), dtype=arr.dtype)
        for i in range(th):
            y_start = int(i * bh)
            y_end = int((i + 1) * bh)
            y_end = min(y_end, h_temp)
            if y_end > y_start:
                result[i, :] = agg_func(temp[y_start:y_end, :], axis=0)
    
    return result


def taint_sample(
    features: torch.Tensor, 
    target_class: int, 
    num_classes: int,
    taint_type: str = 'additive',
    taint_value: float = 1.0,
    taint_horizontal_size: int = 6,
    taint_vertical_size: int = 6,
    taint_horizontal_step: int = 100,
    taint_vertical_step: int = 6
) -> torch.Tensor:
    """
    Marca amostra com valor sentinela para debug - versão configurável para CNN2.
    
    MODOS DE TAINTING:
    
    1. additive (padrão): Adiciona bloco de taint sobre os dados originais
       - Mantém dados originais intactos
       - Coloca bloco de taint_value em posição específica da classe
       - Útil para testar se a rede consegue "ver" o sinal
    
    2. override: Zera TODA a entrada e coloca apenas o taint
       - Entrada fica toda zerada exceto o bloco de taint
       - Teste mais radical: a rede SÓ vê o taint
       - Se não aprender com isso, há problema sério na arquitetura
    
    ESTRATÉGIA DE POSICIONAMENTO:
    - A posição do bloco é calculada com base na classe:
      - start_col = target_class * taint_horizontal_step
      - start_row = target_class * taint_vertical_step
    - O tamanho do bloco é definido por:
      - largura = taint_horizontal_size
      - altura = taint_vertical_size
    
    Exemplo com horizontal_step=100, vertical_step=6, horizontal_size=6, vertical_size=6:
        - Classe 0: bloco em (row=0, col=0) até (row=6, col=6)
        - Classe 1: bloco em (row=6, col=100) até (row=12, col=106)
        - Classe 2: bloco em (row=12, col=200) até (row=18, col=206)
        - etc.
    
    Args:
        features: Tensor de entrada (1D ou 2D)
        target_class: Classe de saída (0 a num_classes-1)
        num_classes: Número total de classes
        taint_type: 'additive' (mantém dados, adiciona taint) ou 
                    'override' (zera tudo, só taint)
        taint_value: Valor do taint (default: 1.0)
        taint_horizontal_size: Largura do bloco de taint em colunas (default: 6)
        taint_vertical_size: Altura do bloco de taint em linhas (default: 6)
        taint_horizontal_step: Passo horizontal por classe (default: 100)
        taint_vertical_step: Passo vertical por classe (default: 6)
        
    Returns:
        Tensor com tainting aplicado
    """
    if taint_type == 'override':
        # Override: zerar toda a entrada primeiro
        features = torch.zeros_like(features)
    else:
        # Additive: manter dados originais
        features = features.clone()
    
    if features.ndim == 1:
        # Dados 1D: comportamento simples baseado em horizontal_step
        input_size = features.shape[0]
        position = int(target_class * taint_horizontal_step)
        end_position = min(position + taint_horizontal_size, input_size)
        if position < input_size:
            features[position:end_position] = taint_value
    
    elif features.ndim == 2:
        # Dados 2D: criar bloco visível para CNN
        num_rows, num_cols = features.shape
        
        # Posição do bloco baseada na classe e nos passos
        start_col = int(target_class * taint_horizontal_step)
        start_row = int(target_class * taint_vertical_step)
        
        # Garantir que o bloco está dentro dos limites
        end_row = min(start_row + taint_vertical_size, num_rows)
        end_col = min(start_col + taint_horizontal_size, num_cols)
        
        # Garantir que as posições iniciais estão dentro dos limites
        if start_row < num_rows and start_col < num_cols:
            features[start_row:end_row, start_col:end_col] = taint_value
    
    else:
        # Mais dimensões: flatten e aplicar
        flat = features.view(-1)
        input_size = flat.numel()
        position = int(target_class * taint_horizontal_step)
        end_position = min(position + taint_horizontal_size, input_size)
        if position < input_size:
            flat[position:end_position] = taint_value
    
    return features


def generate_experiment_name(config: Dict) -> str:
    """
    Gera nome do experimento baseado nos parâmetros de configuração.
    
    Formato NN: nn_<alphagenome_outputs>_<haplotype_mode>_<window_center_size>_
                <normalization_method>_<hidden_layers>_<activation>_<dropout_rate>_<optimizer>
    
    Formato CNN: cnn_<alphagenome_outputs>_<haplotype_mode>_<window_center_size>_
                 <normalization_method>_k<kernel>_f<filters>_s<stride>_p<padding>_[pool<size>_]
                 <hidden_layers>_<activation>_<dropout_rate>_<optimizer>
    
    Formato CNN2: cnn2_<alphagenome_outputs>_<haplotype_mode>_<window_center_size>_
                  <normalization_method>_s1k<kernel>f<filters>_s2f<filters>_s3f<filters>_
                  fc<hidden_size>_<hidden_layers>_<activation>_<dropout_rate>_<optimizer>
    
    Exemplo NN:   nn_atac_H1_1002_log_L100-40_relu_0.0_adam
    Exemplo CNN:  cnn_atac_H1_1002_log_k5x5_f20_s5_p0_L100-40_relu_0.0_adam
    Exemplo CNN2: cnn2_rna_seq_H1_32768_log_s1k6x32f16_s2f32_s3f64_fc128_L100-40_relu_0.4_adam
    
    Args:
        config: Dicionário de configuração
        
    Returns:
        Nome do experimento
    """
    # Tipo de rede
    model_type = config['model'].get('type', 'NN').lower()
    
    # Parâmetros de entrada de dados
    alphagenome_outputs = '_'.join(config['dataset_input']['alphagenome_outputs'])
    haplotype_mode = config['dataset_input']['haplotype_mode']
    window_center_size = config['dataset_input']['window_center_size']
    normalization_method = config['dataset_input'].get('normalization_method', 'zscore')
    
    # Parâmetros do modelo
    hidden_layers = config['model']['hidden_layers']
    hidden_layers_str = 'L' + '-'.join(map(str, hidden_layers))
    
    activation = config['model']['activation']
    dropout_rate = config['model']['dropout_rate']
    
    # Parâmetros de treinamento
    optimizer = config['training']['optimizer']
    
    # Construir parte específica do modelo
    if model_type == 'cnn':
        # Parâmetros CNN
        cnn_config = config['model']['cnn']
        kernel_size = cnn_config['kernel_size']
        kernel_str = f"k{kernel_size[0]}x{kernel_size[1]}"
        
        num_filters = cnn_config['num_filters']
        filters_str = f"f{num_filters}"
        
        stride = cnn_config['stride']
        if isinstance(stride, list):
            stride_str = f"s{stride[0]}x{stride[1]}"
        else:
            stride_str = f"s{stride}"
        
        padding = cnn_config['padding']
        if isinstance(padding, list):
            padding_str = f"p{padding[0]}x{padding[1]}"
        else:
            padding_str = f"p{padding}"
        
        # Pool (opcional)
        pool_size = cnn_config.get('pool_size')
        if pool_size is not None:
            pool_str = f"pool{pool_size[0]}x{pool_size[1]}_"
        else:
            pool_str = ""
        
        # Montar nome CNN
        experiment_name = (
            f"cnn_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
            f"{normalization_method}_{kernel_str}_{filters_str}_{stride_str}_{padding_str}_{pool_str}"
            f"{hidden_layers_str}_{activation}_{dropout_rate}_{optimizer}"
        )
    
    elif model_type == 'cnn2':
        # Parâmetros CNN2
        cnn2_config = config['model'].get('cnn2', {})
        
        # Stage 1
        num_filters_s1 = cnn2_config.get('num_filters_stage1', 16)
        kernel_s1 = cnn2_config.get('kernel_stage1', [6, 32])
        stage1_str = f"s1k{kernel_s1[0]}x{kernel_s1[1]}f{num_filters_s1}"
        
        # Stage 2
        num_filters_s2 = cnn2_config.get('num_filters_stage2', 32)
        stage2_str = f"s2f{num_filters_s2}"
        
        # Stage 3
        num_filters_s3 = cnn2_config.get('num_filters_stage3', 64)
        stage3_str = f"s3f{num_filters_s3}"
        
        # Global pooling type
        pool_type = cnn2_config.get('global_pool_type', 'max')
        pool_str = f"gp{pool_type}"
        
        # FC hidden size
        fc_hidden_size = cnn2_config.get('fc_hidden_size', 128)
        fc_str = f"fc{fc_hidden_size}"
        
        # Montar nome CNN2
        experiment_name = (
            f"cnn2_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
            f"{normalization_method}_{stage1_str}_{stage2_str}_{stage3_str}_{pool_str}_{fc_str}_"
            f"{hidden_layers_str}_{activation}_{dropout_rate}_{optimizer}"
        )
    
    else:
        # Montar nome NN (default)
        experiment_name = (
            f"nn_{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
            f"{normalization_method}_{hidden_layers_str}_{activation}_{dropout_rate}_{optimizer}"
        )
    
    return experiment_name


def generate_dataset_name(config: Dict) -> str:
    """
    Gera nome único para o dataset baseado apenas nos parâmetros de processamento.
    
    Experimentos com os mesmos parâmetros de dataset compartilharão o mesmo cache.
    
    Formato: <alphagenome_outputs>_<haplotype_mode>_<window_center_size>_
             ds<downsample_factor>_<normalization_method>_<balancing>_
             split<train>-<val>-<test>_seed<random_seed>
    
    Exemplo: atac_H1_1002_ds1_log_strat_split0.7-0.15-0.15_seed42
    
    Args:
        config: Dicionário de configuração
        
    Returns:
        Nome do dataset
    """
    # Parâmetros de entrada de dados
    alphagenome_outputs = '_'.join(config['dataset_input']['alphagenome_outputs'])
    haplotype_mode = config['dataset_input']['haplotype_mode']
    window_center_size = config['dataset_input']['window_center_size']
    downsample_factor = config['dataset_input']['downsample_factor']
    normalization_method = config['dataset_input'].get('normalization_method', 'zscore')
    
    # Parâmetros de split
    train_split = config['data_split']['train_split']
    val_split = config['data_split']['val_split']
    test_split = config['data_split']['test_split']
    random_seed = config['data_split']['random_seed']
    
    # Estratégia de balanceamento
    balancing_strategy = config['data_split'].get('balancing_strategy', 'stratified')
    balance_str = 'strat' if balancing_strategy == 'stratified' else 'shuf'
    
    # Montar nome do dataset
    dataset_name = (
        f"{alphagenome_outputs}_{haplotype_mode}_{window_center_size}_"
        f"ds{downsample_factor}_{normalization_method}_{balance_str}_"
        f"split{train_split}-{val_split}-{test_split}_seed{random_seed}"
    )
    
    return dataset_name


def get_dataset_cache_dir(config: Dict) -> Path:
    """
    Retorna o diretório onde o dataset deve ser cacheado.
    
    O cache do dataset é compartilhado entre experimentos com mesmos parâmetros de dados.
    
    Args:
        config: Dicionário de configuração
        
    Returns:
        Path do diretório de cache do dataset
    """
    base_cache_dir = Path(config['dataset_input']['processed_cache_dir'])
    dataset_name = generate_dataset_name(config)
    dataset_cache_dir = base_cache_dir / 'datasets' / dataset_name
    
    return dataset_cache_dir


def setup_experiment_dir(config: Dict, config_path: str) -> Path:
    """
    Cria e configura diretório do experimento (sem cache de dataset).
    
    Args:
        config: Dicionário de configuração
        config_path: Caminho do arquivo de configuração original
        
    Returns:
        Path do diretório do experimento
    """
    # Gerar nome do experimento
    experiment_name = generate_experiment_name(config)
    
    # Criar path do experimento
    base_cache_dir = Path(config['dataset_input']['processed_cache_dir'])
    experiment_dir = base_cache_dir / experiment_name
    
    # Criar diretórios
    experiment_dir.mkdir(parents=True, exist_ok=True)
    (experiment_dir / 'models').mkdir(exist_ok=True)
    
    # Copiar arquivo de configuração
    config_copy_path = experiment_dir / 'config.yaml'
    shutil.copy(config_path, config_copy_path)
    
    console.print(f"[green]📁 Diretório do experimento:[/green] {experiment_dir}")
    console.print(f"[green]   Nome:[/green] {experiment_name}")
    
    return experiment_dir


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
        Computa parâmetros de normalização por track (por gene/ontologia).
        
        Cada track é normalizada independentemente para evitar que tracks com
        valores altos dominem o treinamento.
        
        Returns:
            Dict com parâmetros de normalização por track
        """
        # Verificar se há valor pré-definido no config (não suportado para per-track)
        predefined_value = self.config['dataset_input'].get('normalization_value', 0)
        
        if predefined_value != 0 and predefined_value is not None:
            console.print(f"[yellow]⚠ AVISO: normalization_value não é suportado para normalização per-track[/yellow]")
            console.print(f"[yellow]  Ignorando e computando parâmetros por track...[/yellow]")
        
        # Coletar valores por track
        track_values = None  # Será inicializado no primeiro sample
        num_processed = 0
        num_errors = 0
        
        total_samples = len(self.base_dataset)
        console.print(f"[cyan]Iniciando computação de normalização POR TRACK...[/cyan]")
        console.print(f"[cyan]  • Amostras: {total_samples}[/cyan]")
        console.print(f"[cyan]  • Método: {self.normalization_method}[/cyan]")
        console.print(f"[cyan]  • Outputs: {', '.join(self.alphagenome_outputs)}[/cyan]")
        console.print(f"[cyan]  • Haplótipo: {self.haplotype_mode}[/cyan]")
        console.print(f"[cyan]  • Window center: {self.window_center_size} bases[/cyan]")
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task(
                "Coletando valores por track...",
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
                        # Inicializar track_values no primeiro sample válido
                        if track_values is None:
                            num_tracks = processed_array.shape[0]
                            track_values = [[] for _ in range(num_tracks)]
                            console.print(f"[cyan]  • Total de tracks detectadas: {num_tracks}[/cyan]")
                        
                        # Adicionar valores de cada track
                        for track_idx in range(processed_array.shape[0]):
                            track_values[track_idx].append(processed_array[track_idx])
                        
                        num_processed += 1
                    else:
                        console.print(f"[yellow]  ⚠ Amostra {sample_id} (idx={idx}): Nenhum dado válido encontrado[/yellow]")
                        num_errors += 1
                        
                except Exception as e:
                    num_errors += 1
                    console.print(f"[yellow]  ⚠ Erro ao processar amostra {idx}: {e}[/yellow]")
                progress.update(task, advance=1)
        
        if track_values is None or len(track_values) == 0:
            console.print(f"[red]ERRO: Nenhuma amostra válida para normalização![/red]")
            console.print(f"[red]  • Amostras processadas: {num_processed}[/red]")
            console.print(f"[red]  • Amostras com erro: {num_errors}[/red]")
            console.print(f"[red]  • Total esperado: {len(self.base_dataset)}[/red]")
            return {
                'method': self.normalization_method,
                'per_track': False,
                'mean': 0.0,
                'std': 1.0
            }
        
        num_tracks = len(track_values)
        console.print(f"\n[cyan]Computando parâmetros para {num_tracks} tracks...[/cyan]")
        
        # Concatenar valores de cada track e computar parâmetros
        track_params = []
        
        for track_idx in range(num_tracks):
            # Concatenar valores desta track de todas as amostras
            track_array = np.concatenate(track_values[track_idx])
            
            if self.normalization_method == 'zscore':
                mean = float(np.mean(track_array))
                std = float(np.std(track_array))
                if std < 1e-8:
                    std = 1.0
                track_params.append({'mean': mean, 'std': std})
                
            elif self.normalization_method == 'minmax_keep_zero':
                nonzero_values = track_array[track_array > 0]
                if len(nonzero_values) > 0:
                    xmax = float(nonzero_values.max())
                else:
                    xmax = 1.0
                track_params.append({'max': xmax})
                
            elif self.normalization_method == 'log':
                nonzero_values = track_array[track_array > 0]
                if len(nonzero_values) > 0:
                    xmax = float(nonzero_values.max())
                    log_max = float(np.log1p(xmax))
                else:
                    log_max = 1.0
                track_params.append({'log_max': log_max})
                
            else:
                # Fallback para zscore
                mean = float(np.mean(track_array))
                std = float(np.std(track_array))
                if std < 1e-8:
                    std = 1.0
                track_params.append({'mean': mean, 'std': std})
        
        # Construir resultado
        params = {
            'method': self.normalization_method,
            'per_track': True,
            'num_tracks': num_tracks,
            'track_params': track_params
        }
        
        # Exibir resumo
        if self.normalization_method == 'zscore':
            console.print(f"\n[bold green]✓ Normalização Z-score Por-Track Concluída:[/bold green]")
            console.print(f"  • Amostras processadas: {num_processed}/{len(self.base_dataset)}")
            console.print(f"  • Amostras com erro: {num_errors}")
            console.print(f"  • Tracks normalizadas: {num_tracks}")
            # Mostrar algumas tracks como exemplo
            for i in [0, 1, num_tracks-1]:
                if i < len(track_params):
                    console.print(f"  • Track {i}: mean={track_params[i]['mean']:.6f}, std={track_params[i]['std']:.6f}")
            
        elif self.normalization_method == 'minmax_keep_zero':
            console.print(f"\n[bold green]✓ Normalização MinMax Por-Track Concluída:[/bold green]")
            console.print(f"  • Amostras processadas: {num_processed}/{len(self.base_dataset)}")
            console.print(f"  • Amostras com erro: {num_errors}")
            console.print(f"  • Tracks normalizadas: {num_tracks}")
            for i in [0, 1, num_tracks-1]:
                if i < len(track_params):
                    console.print(f"  • Track {i}: max={track_params[i]['max']:.6f}")
            
        elif self.normalization_method == 'log':
            console.print(f"\n[bold green]✓ Normalização Logarítmica Por-Track Concluída:[/bold green]")
            console.print(f"  • Amostras processadas: {num_processed}/{len(self.base_dataset)}")
            console.print(f"  • Amostras com erro: {num_errors}")
            console.print(f"  • Tracks normalizadas: {num_tracks}")
            for i in [0, 1, num_tracks-1]:
                if i < len(track_params):
                    console.print(f"  • Track {i}: log1p(max)={track_params[i]['log_max']:.6f}")
        
        return params
    
    def _create_target_mappings(self):
        """Cria mapeamentos entre targets e índices."""
        
        # Prioridade 1: Tentar carregar dos metadados do dataset (rápido - sem I/O de dados)
        dataset_metadata = getattr(self.base_dataset, 'dataset_metadata', None)
        
        # Se não está no base_dataset, tentar carregar do arquivo
        if dataset_metadata is None:
            dataset_dir = Path(self.config['dataset_input']['dataset_dir'])
            metadata_file = dataset_dir / 'dataset_metadata.json'
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    dataset_metadata = json.load(f)
        
        classes_from_metadata = None
        
        if dataset_metadata:
            if self.prediction_target == 'superpopulation':
                # Tentar superpopulation_distribution primeiro
                superpop_dist = dataset_metadata.get('superpopulation_distribution', {})
                if superpop_dist:
                    classes_from_metadata = list(superpop_dist.keys())
                else:
                    # Fallback: extrair de individuals_pedigree
                    pedigree = dataset_metadata.get('individuals_pedigree', {})
                    if pedigree:
                        classes_from_metadata = list(set(
                            p.get('superpopulation') for p in pedigree.values() 
                            if p.get('superpopulation')
                        ))
            elif self.prediction_target == 'population':
                pop_dist = dataset_metadata.get('population_distribution', {})
                if pop_dist:
                    classes_from_metadata = list(pop_dist.keys())
                else:
                    # Fallback: extrair de individuals_pedigree
                    pedigree = dataset_metadata.get('individuals_pedigree', {})
                    if pedigree:
                        classes_from_metadata = list(set(
                            p.get('population') for p in pedigree.values() 
                            if p.get('population')
                        ))
        
        if classes_from_metadata:
            console.print(f"\n[cyan]Carregando classes dos metadados do dataset...[/cyan]")
            sorted_targets = sorted(classes_from_metadata)
            
            self.target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
            self.idx_to_target = {idx: target for target, idx in self.target_to_idx.items()}
            
            console.print(f"[green]✓ Mapeamento de targets criado: {len(self.target_to_idx)} classes[/green]")
            console.print(f"[cyan]Classes: {sorted_targets}[/cyan]")
            return
        
        # Prioridade 2: Verificar se há classes conhecidas no config
        known_classes = self.config.get('output', {}).get('known_classes')
        
        if known_classes is not None and len(known_classes) > 0:
            # Usar classes pré-definidas (fallback)
            console.print(f"\n[cyan]Usando classes conhecidas do config (fallback)...[/cyan]")
            console.print(f"[yellow]⚠ Recomendado: Remova 'known_classes' do YAML para usar metadados do dataset[/yellow]")
            sorted_targets = sorted(known_classes)
            
            self.target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
            self.idx_to_target = {idx: target for target, idx in self.target_to_idx.items()}
            
            console.print(f"[green]✓ Mapeamento de targets criado: {len(self.target_to_idx)} classes[/green]")
            console.print(f"[cyan]Classes: {sorted_targets}[/cyan]")
            return
        
        # Prioridade 3: Escanear o dataset (mais lento)
        console.print(f"\n[yellow]Classes não encontradas nos metadados ou config.[/yellow]")
        console.print(f"[cyan]Escaneando dataset para descobrir classes...[/cyan]")
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
        Processa todas as janelas e retorna matriz 2D.
        
        Cada linha da matriz representa um haplótipo de um único tipo de saída.
        
        Args:
            windows: Dicionário de janelas do input_data
            
        Returns:
            Array numpy 2D com shape [num_rows, effective_size]
            onde num_rows = num_windows * num_outputs * num_haplotypes
        """
        processed_rows = []
        
        for window_name, window_data in windows.items():
            # Processar cada haplótipo
            if self.haplotype_mode in ['H1', 'H1+H2']:
                h1_rows = self._process_haplotype(window_data.get('predictions_h1', {}))
                if h1_rows is not None:
                    # h1_rows pode ser 1D (um output) ou 2D (múltiplos outputs)
                    if h1_rows.ndim == 1:
                        processed_rows.append(h1_rows.reshape(1, -1))
                    else:
                        processed_rows.append(h1_rows)
            
            if self.haplotype_mode in ['H2', 'H1+H2']:
                h2_rows = self._process_haplotype(window_data.get('predictions_h2', {}))
                if h2_rows is not None:
                    if h2_rows.ndim == 1:
                        processed_rows.append(h2_rows.reshape(1, -1))
                    else:
                        processed_rows.append(h2_rows)
        
        if len(processed_rows) == 0:
            # Retornar array vazio 2D
            return np.array([[]]).reshape(0, 0)
        
        # Empilhar todas as linhas verticalmente para criar matriz 2D
        result = np.vstack(processed_rows)
        
        return result
    
    def _process_haplotype(self, predictions: Dict) -> Optional[np.ndarray]:
        """
        Processa predições de um haplótipo.
        
        Cada tipo de saída (ex: atac, rna_seq) pode ter múltiplas tracks.
        - Se array é 1D: uma única track
        - Se array é 2D com shape (bases, num_tracks): cada coluna é uma track
        
        Cada track gera uma linha separada no output.
        
        Args:
            predictions: Dict com {output_type: np.ndarray}
                         array pode ser 1D (1 track) ou 2D (múltiplas tracks)
            
        Returns:
            Array 2D com shape [num_rows, effective_size] ou None
            num_rows = total de tracks de todos os output_types
            Se apenas uma track, retorna array 1D com shape [effective_size]
        """
        rows = []
        
        for output_type in self.alphagenome_outputs:
            if output_type in predictions:
                array = predictions[output_type]
                
                # Verificar se é 2D com múltiplas colunas (tracks)
                if array.ndim == 2 and array.shape[1] > 1:
                    # Array 2D: cada coluna é uma track separada
                    num_tracks = array.shape[1]
                    
                    for track_idx in range(num_tracks):
                        # Extrair coluna (track específica)
                        track_array = array[:, track_idx]
                        
                        # Extrair trecho central
                        center_array = self._extract_center(track_array)
                        
                        # Aplicar downsampling
                        downsampled = self._downsample(center_array)
                        
                        rows.append(downsampled)
                else:
                    # Array 1D ou 2D com apenas 1 coluna: track única
                    if array.ndim > 1:
                        array = array.flatten()
                    
                    # Extrair trecho central
                    center_array = self._extract_center(array)
                    
                    # Aplicar downsampling
                    downsampled = self._downsample(center_array)
                    
                    rows.append(downsampled)
        
        if len(rows) == 0:
            return None
        
        # Se múltiplas tracks, empilhar como matriz 2D
        # Se apenas uma track, retornar como array 1D (será convertido para 2D em _process_windows)
        if len(rows) == 1:
            return rows[0]
        else:
            return np.vstack(rows)
    
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
            features_tensor tem shape [num_rows, effective_size] (2D)
        """
        input_data, output_data = self.base_dataset[idx]
        
        # Processar janelas (retorna matriz 2D)
        features = self._process_windows(input_data['windows'])
        
        # Converter para tensor 2D
        features_tensor = torch.FloatTensor(features)
        
        # Aplicar normalização conforme método escolhido
        method = self.normalization_params.get('method', 'zscore')
        per_track = self.normalization_params.get('per_track', False)
        
        if per_track:
            # Normalização por track (cada linha/track normalizada independentemente)
            track_params = self.normalization_params['track_params']
            normalized_rows = []
            
            for track_idx in range(features_tensor.shape[0]):
                track_tensor = features_tensor[track_idx:track_idx+1, :]  # Manter 2D
                params = track_params[track_idx]
                
                if method == 'zscore':
                    normalized = zscore_normalize(track_tensor, params['mean'], params['std'])
                elif method == 'minmax_keep_zero':
                    normalized = minmax_keep_zero(track_tensor, params['max'])
                elif method == 'log':
                    normalized = log_normalize(track_tensor, params['log_max'])
                else:
                    # Fallback para zscore
                    normalized = zscore_normalize(track_tensor, params.get('mean', 0.0), params.get('std', 1.0))
                
                normalized_rows.append(normalized)
            
            features_tensor = torch.cat(normalized_rows, dim=0)
        else:
            # Normalização global (backwards compatibility)
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
                    debug_config = self.config.get('debug', {})
                    features_tensor = taint_sample(
                        features_tensor, target_idx, num_classes,
                        taint_type=debug_config.get('taint_type', 'additive'),
                        taint_value=debug_config.get('taint_value', 1.0),
                        taint_horizontal_size=debug_config.get('taint_horizontal_size', 6),
                        taint_vertical_size=debug_config.get('taint_vertical_size', 6),
                        taint_horizontal_step=debug_config.get('taint_horizontal_step', 100),
                        taint_vertical_step=debug_config.get('taint_vertical_step', 6)
                    )
            else:
                # Target desconhecido, usar -1
                target_tensor = torch.tensor(-1, dtype=torch.long)
        
        return features_tensor, target_tensor
    
    def get_num_classes(self) -> int:
        """Retorna número de classes."""
        return len(self.target_to_idx)
    
    def get_input_shape(self) -> Tuple[int, int]:
        """
        Calcula shape da entrada da rede.
        
        Returns:
            Tupla (num_rows, effective_size)
            onde num_rows = num_windows * num_outputs * num_haplotypes
        """
        # Pegar primeira amostra válida para calcular
        for idx in range(len(self.base_dataset)):
            try:
                features, _ = self[idx]
                return tuple(features.shape)
            except Exception:
                continue
        
        # Fallback: calcular teoricamente
        num_windows = len(self.base_dataset[0][0]['windows'])
        num_outputs = len(self.alphagenome_outputs)
        num_haplotypes = 2 if self.haplotype_mode == 'H1+H2' else 1
        
        effective_size = self.window_center_size // self.downsample_factor
        num_rows = num_windows * num_outputs * num_haplotypes
        
        return (num_rows, effective_size)
    
    def get_input_size(self) -> int:
        """
        DEPRECATED: Retorna tamanho total da entrada (para retrocompatibilidade).
        Use get_input_shape() para obter shape 2D.
        """
        num_rows, effective_size = self.get_input_shape()
        return num_rows * effective_size


# ==============================================================================
# MODEL ARCHITECTURE
# ==============================================================================

class NNAncestryPredictor(nn.Module):
    """
    Rede neural totalmente conectada (NN) para predição de ancestralidade.
    
    Arquitetura:
    - Flatten (converte entrada 2D para 1D)
    - Camada de entrada (tamanho variável)
    - Camadas ocultas (configurável)
    - Camada de saída (softmax para classificação ou linear para regressão)
    """
    
    def __init__(self, config: Dict, input_shape: Tuple[int, int], num_classes: int):
        """
        Inicializa o modelo.
        
        Args:
            config: Configuração do YAML
            input_shape: Tupla (num_rows, effective_size) do shape de entrada 2D
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(NNAncestryPredictor, self).__init__()
        
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Determinar quais linhas usar (genes específicos)
        genes_to_use = config['dataset_input'].get('genes_to_use', None)
        tracks_per_gene = 6  # rna_seq tem 6 tracks
        
        # Ordem dos genes no dataset (carregada do cache metadata)
        # Fallback para ordem padrão se não disponível
        GENE_ORDER = config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        
        if genes_to_use:
            # Validar genes
            for gene in genes_to_use:
                if gene not in GENE_ORDER:
                    raise ValueError(f"Gene inválido: {gene}. Opções: {GENE_ORDER}")
            
            # Criar lista de índices de genes (mantendo ordem do dataset)
            gene_indices = [i for i, gene in enumerate(GENE_ORDER) if gene in genes_to_use]
            self.genes_selected = [GENE_ORDER[i] for i in gene_indices]
            
            # Calcular linhas a extrair (cada gene = 6 linhas consecutivas)
            self.rows_to_use = []
            for gene_idx in gene_indices:
                start_row = gene_idx * tracks_per_gene
                self.rows_to_use.extend(range(start_row, start_row + tracks_per_gene))
        else:
            # Default: usar todos os genes
            self.genes_selected = GENE_ORDER
            self.rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))
        
        # Ajustar input_shape para os genes selecionados
        num_rows_original = input_shape[0]
        effective_size = input_shape[1]
        num_rows = len(self.rows_to_use)
        self.input_shape = (num_rows, effective_size)
        self.input_size = num_rows * effective_size  # Tamanho após flatten
        
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
        # Arquitetura: camadas hidden (TODAS com ativação) + camada de saída (LINEAR antes do softmax)
        layers = []
        prev_size = self.input_size
        
        # TODAS as camadas hidden (com ativação configurada)
        for hidden_size in hidden_layers:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(self.activation)
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Camada de saída (sempre LINEAR, sem ativação - logits)
        layers.append(nn.Linear(prev_size, num_classes))
        
        self.network = nn.Sequential(*layers)
        
        # Inicializar pesos apropriadamente
        self._initialize_weights(activation_type)
        
        console.print(f"[green]Modelo NN criado:[/green]")
        console.print(f"  • Input shape original: {input_shape[0]} x {input_shape[1]}")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {len(self.rows_to_use)} linhas")
        console.print(f"  • Genes: {', '.join(self.genes_selected)}")
        console.print(f"  • Input shape efetivo: {self.input_shape[0]} x {self.input_shape[1]}")
        console.print(f"  • Input size (após flatten): {self.input_size}")
        console.print(f"  • Hidden layers: {hidden_layers}")
        console.print(f"  • Ativação: {activation_type} (em todas as camadas hidden)")
        console.print(f"  • Arquitetura: flatten → camadas hidden → saída (logits)")
        console.print(f"  • Output size: {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo logits (não probabilidades)
        """
        # Selecionar apenas as linhas dos genes escolhidos
        # [batch, num_rows_original, effective_size] -> [batch, num_rows_selected, effective_size]
        x = x[:, self.rows_to_use, :]
        
        # Flatten: [batch, num_rows_selected, effective_size] -> [batch, num_rows_selected * effective_size]
        x = x.view(x.size(0), -1)
        
        logits = self.network(x)
        return logits
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Calcula probabilidades de classe usando softmax sobre os logits.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo probabilidades
        """
        logits = self.forward(x)
        return torch.softmax(logits, dim=1)
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self, activation_type: str):
        """
        Inicializa pesos apropriadamente de acordo com tipo de ativação.
        
        Estratégia:
        - Camadas hidden (com ativação):
            * ReLU: He/Kaiming initialization (fan_in)
            * Tanh/Sigmoid: Xavier/Glorot initialization
        - Última camada (saída, linear antes do softmax): Xavier initialization
        - Bias: zeros (padrão recomendado)
        
        Args:
            activation_type: Tipo de função de ativação das camadas hidden
        """
        layer_count = 0
        total_layers = sum(1 for m in self.modules() if isinstance(m, nn.Linear))
        
        for m in self.modules():
            if isinstance(m, nn.Linear):
                layer_count += 1
                
                # Última camada (saída): sempre Xavier
                # (linear antes do softmax, sem ativação adicional)
                if layer_count == total_layers:
                    nn.init.xavier_normal_(m.weight)
                # Camadas hidden: depende da ativação configurada
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
        console.print(f"  • Camadas hidden: {activation_type} initialization")
        console.print(f"  • Camada de saída: Xavier initialization")
        console.print(f"  • Bias: zeros")


class CNNAncestryPredictor(nn.Module):
    """
    Rede neural convolucional (CNN) para predição de ancestralidade.
    
    Arquitetura:
    - Input: [batch, num_rows, effective_size] como imagem com 1 canal
    - Conv2D: camada convolucional com kernel configurável
    - Ativação (ReLU/Tanh/Sigmoid)
    - MaxPool2D (opcional)
    - Flatten
    - Camadas fully connected (hidden_layers)
    - Output linear → Softmax (classificação)
    """
    
    def __init__(self, config: Dict, input_shape: Tuple[int, int], num_classes: int):
        """
        Inicializa o modelo CNN.
        
        Args:
            config: Configuração do YAML
            input_shape: Tupla (num_rows, effective_size) do shape de entrada 2D
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(CNNAncestryPredictor, self).__init__()
        
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Determinar quais linhas usar (genes específicos)
        genes_to_use = config['dataset_input'].get('genes_to_use', None)
        tracks_per_gene = 6  # rna_seq tem 6 tracks
        
        # Ordem dos genes no dataset (carregada do cache metadata)
        # Fallback para ordem padrão se não disponível
        GENE_ORDER = config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        
        if genes_to_use:
            # Validar genes
            for gene in genes_to_use:
                if gene not in GENE_ORDER:
                    raise ValueError(f"Gene inválido: {gene}. Opções: {GENE_ORDER}")
            
            # Criar lista de índices de genes (mantendo ordem do dataset)
            gene_indices = [i for i, gene in enumerate(GENE_ORDER) if gene in genes_to_use]
            self.genes_selected = [GENE_ORDER[i] for i in gene_indices]
            
            # Calcular linhas a extrair (cada gene = 6 linhas consecutivas)
            self.rows_to_use = []
            for gene_idx in gene_indices:
                start_row = gene_idx * tracks_per_gene
                self.rows_to_use.extend(range(start_row, start_row + tracks_per_gene))
        else:
            # Default: usar todos os genes
            self.genes_selected = GENE_ORDER
            self.rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))
        
        # Ajustar input_shape para os genes selecionados
        num_rows_original = input_shape[0]
        effective_size = input_shape[1]
        num_rows = len(self.rows_to_use)
        self.input_shape = (num_rows, effective_size)
        
        # Parâmetros CNN
        cnn_config = config['model']['cnn']
        kernel_size = tuple(cnn_config['kernel_size'])  # [height, width]
        num_filters = cnn_config['num_filters']
        
        # Stride pode ser escalar ou lista [vertical, horizontal]
        stride_config = cnn_config['stride']
        if isinstance(stride_config, list):
            stride = tuple(stride_config)
        else:
            stride = (stride_config, stride_config)  # Usar mesmo valor para ambas dimensões
        
        # Padding pode ser escalar ou lista [vertical, horizontal]
        padding_config = cnn_config['padding']
        if isinstance(padding_config, list):
            padding = tuple(padding_config)
        else:
            padding = padding_config  # PyTorch aceita int ou tuple
        
        pool_size = cnn_config.get('pool_size')  # Pode ser None
        
        # Parâmetros gerais
        activation_type = config['model']['activation']
        dropout_rate = config['model']['dropout_rate']
        hidden_layers = config['model']['hidden_layers']
        
        # Escolher função de ativação
        if activation_type == 'relu':
            self.activation = nn.ReLU()
        elif activation_type == 'tanh':
            self.activation = nn.Tanh()
        elif activation_type == 'sigmoid':
            self.activation = nn.Sigmoid()
        else:
            raise ValueError(f"Ativação não suportada: {activation_type}")
        
        # Camada convolucional
        # Input: [batch, 1, num_rows, effective_size]
        # Output: [batch, num_filters, out_height, out_width]
        self.conv1 = nn.Conv2d(
            in_channels=1,
            out_channels=num_filters,
            kernel_size=kernel_size,
            stride=stride,
            padding=padding
        )
        
        # Calcular dimensões após convolução (usando input_shape ajustado)
        num_rows, effective_size = self.input_shape
        
        # Extrair padding como tupla
        if isinstance(padding, int):
            pad_h, pad_w = padding, padding
        else:
            pad_h, pad_w = padding
        
        # Extrair stride como tupla
        if isinstance(stride, int):
            stride_h, stride_w = stride, stride
        else:
            stride_h, stride_w = stride
        
        # Fórmula: out_size = (in_size + 2*padding - kernel_size) / stride + 1
        conv_out_h = (num_rows + 2*pad_h - kernel_size[0]) // stride_h + 1
        conv_out_w = (effective_size + 2*pad_w - kernel_size[1]) // stride_w + 1
        
        # Pooling (opcional)
        if pool_size is not None:
            pool_size = tuple(pool_size)
            self.pool = nn.MaxPool2d(kernel_size=pool_size)
            pool_out_h = conv_out_h // pool_size[0]
            pool_out_w = conv_out_w // pool_size[1]
        else:
            self.pool = None
            pool_out_h = conv_out_h
            pool_out_w = conv_out_w
        
        # Tamanho após flatten
        flattened_size = num_filters * pool_out_h * pool_out_w
        
        # Camadas fully connected
        fc_layers = []
        prev_size = flattened_size
        
        for hidden_size in hidden_layers:
            fc_layers.append(nn.Linear(prev_size, hidden_size))
            fc_layers.append(self.activation)
            if dropout_rate > 0:
                fc_layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Camada de saída (linear, sem ativação - logits)
        fc_layers.append(nn.Linear(prev_size, num_classes))
        
        self.fc_network = nn.Sequential(*fc_layers)
        
        # Inicializar pesos
        self._initialize_weights(activation_type)
        
        console.print(f"[green]Modelo CNN criado:[/green]")
        console.print(f"  • Input shape original: {input_shape[0]} x {input_shape[1]}")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {len(self.rows_to_use)} linhas")
        console.print(f"  • Genes: {', '.join(self.genes_selected)}")
        console.print(f"  • Input shape efetivo: {self.input_shape[0]} x {self.input_shape[1]} (1 canal)")
        console.print(f"  • Conv2D: {num_filters} filters, kernel={kernel_size}, stride={stride}, padding={padding}")
        console.print(f"  • Após Conv2D: {num_filters} x {conv_out_h} x {conv_out_w}")
        if pool_size:
            console.print(f"  • MaxPool2D: kernel={pool_size}")
            console.print(f"  • Após Pool: {num_filters} x {pool_out_h} x {pool_out_w}")
        console.print(f"  • Flatten size: {flattened_size}")
        console.print(f"  • FC hidden layers: {hidden_layers}")
        console.print(f"  • Ativação: {activation_type}")
        console.print(f"  • Output size: {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo logits (não probabilidades)
        """
        # Selecionar apenas as linhas dos genes escolhidos
        # [batch, num_rows_original, effective_size] -> [batch, num_rows_selected, effective_size]
        x = x[:, self.rows_to_use, :]
        
        # Adicionar dimensão de canal: [batch, num_rows_selected, effective_size] -> [batch, 1, num_rows_selected, effective_size]
        x = x.unsqueeze(1)
        
        # Convolução + Ativação
        x = self.conv1(x)
        x = self.activation(x)
        
        # Pooling (se habilitado)
        if self.pool is not None:
            x = self.pool(x)
        
        # Flatten
        x = x.view(x.size(0), -1)
        
        # Camadas fully connected
        logits = self.fc_network(x)
        return logits
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Calcula probabilidades de classe usando softmax sobre os logits.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo probabilidades
        """
        logits = self.forward(x)
        return torch.softmax(logits, dim=1)
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self, activation_type: str):
        """
        Inicializa pesos apropriadamente de acordo com tipo de ativação.
        
        Estratégia:
        - Conv2D: He/Kaiming para ReLU, Xavier para Tanh/Sigmoid
        - FC hidden layers: He/Kaiming para ReLU, Xavier para Tanh/Sigmoid
        - Última camada FC (saída): Xavier initialization
        - Bias: zeros
        
        Args:
            activation_type: Tipo de função de ativação
        """
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                # Inicialização convolucional
                if activation_type == 'relu':
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                elif activation_type in ['tanh', 'sigmoid']:
                    nn.init.xavier_normal_(m.weight)
                else:
                    nn.init.xavier_normal_(m.weight)
                
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            
            elif isinstance(m, nn.Linear):
                # Contar quantas camadas lineares existem
                linear_layers = [module for module in self.fc_network.modules() if isinstance(module, nn.Linear)]
                total_layers = len(linear_layers)
                layer_count = sum(1 for module in self.fc_network.modules() 
                                 if isinstance(module, nn.Linear) and id(module) <= id(m))
                
                # Última camada (saída): Xavier
                if layer_count == total_layers:
                    nn.init.xavier_normal_(m.weight)
                # Camadas hidden
                elif activation_type == 'relu':
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                elif activation_type in ['tanh', 'sigmoid']:
                    nn.init.xavier_normal_(m.weight)
                else:
                    nn.init.xavier_normal_(m.weight)
                
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        
        console.print(f"[green]✓ Pesos CNN inicializados:[/green]")
        console.print(f"  • Conv2D: {activation_type} initialization")
        console.print(f"  • FC hidden layers: {activation_type} initialization")
        console.print(f"  • Camada de saída: Xavier initialization")
        console.print(f"  • Bias: zeros")


class CNN2AncestryPredictor(nn.Module):
    """
    Rede neural convolucional avançada (CNN2) para predição de ancestralidade.
    
    Arquitetura multi-layer com global pooling, mais adequada para dados de
    genes com múltiplas tracks (ex: 11 genes × 6 ontologias = 66 tracks).
    
    Arquitetura:
    - Input: [batch, 1, num_rows, effective_size]
    - Conv2D Stage 1: Agrupa tracks de genes (ex: 6→1 por gene)
    - Conv2D Stage 2: Processa ao longo das bases
    - Conv2D Stage 3: Features mais abstratas
    - AdaptiveAvgPool2d: Global pooling → vetor fixo
    - FC: Classificação final
    
    Vantagens sobre CNN simples:
    - Menos parâmetros (~100-200k vs 5M+)
    - Hierarquia de features (padrões locais → globais)
    - Global pooling = invariante ao tamanho da entrada
    """
    
    def __init__(self, config: Dict, input_shape: Tuple[int, int], num_classes: int):
        """
        Inicializa o modelo CNN2.
        
        Args:
            config: Configuração do YAML
            input_shape: Tupla (num_rows, effective_size) do shape de entrada 2D
            num_classes: Número de classes (ou tamanho da saída para regressão)
        """
        super(CNN2AncestryPredictor, self).__init__()
        
        self.config = config
        self.num_classes = num_classes
        self.is_classification = config['output']['prediction_target'] != 'frog_likelihood'
        
        # Determinar quais linhas usar (genes específicos)
        genes_to_use = config['dataset_input'].get('genes_to_use', None)
        tracks_per_gene = 6  # rna_seq tem 6 tracks
        
        # Ordem dos genes no dataset (carregada do cache metadata)
        # Fallback para ordem padrão se não disponível
        GENE_ORDER = config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        
        if genes_to_use:
            # Validar genes
            for gene in genes_to_use:
                if gene not in GENE_ORDER:
                    raise ValueError(f"Gene inválido: {gene}. Opções: {GENE_ORDER}")
            
            # Criar lista de índices de genes (mantendo ordem do dataset)
            gene_indices = [i for i, gene in enumerate(GENE_ORDER) if gene in genes_to_use]
            self.genes_selected = [GENE_ORDER[i] for i in gene_indices]
            
            # Calcular linhas a extrair (cada gene = 6 linhas consecutivas)
            self.rows_to_use = []
            for gene_idx in gene_indices:
                start_row = gene_idx * tracks_per_gene
                self.rows_to_use.extend(range(start_row, start_row + tracks_per_gene))
        else:
            # Default: usar todos os genes
            self.genes_selected = GENE_ORDER
            self.rows_to_use = list(range(len(GENE_ORDER) * tracks_per_gene))
        
        # Ajustar input_shape para os genes selecionados
        num_rows_original = input_shape[0]
        effective_size = input_shape[1]
        num_rows = len(self.rows_to_use)
        self.input_shape = (num_rows, effective_size)
        
        # Ler parâmetros da configuração CNN2
        cnn2_config = config['model'].get('cnn2', {})
        
        # Stage 1: Primeira convolução (agrupa tracks)
        num_filters_s1 = cnn2_config.get('num_filters_stage1', 16)
        kernel_s1 = tuple(cnn2_config.get('kernel_stage1', [6, 32]))
        stride_s1 = tuple(cnn2_config.get('stride_stage1', [6, 32]))
        
        # Stage 2 e 3: Convoluções subsequentes
        num_filters_s2 = cnn2_config.get('num_filters_stage2', 32)
        num_filters_s3 = cnn2_config.get('num_filters_stage3', 64)
        kernel_s23 = tuple(cnn2_config.get('kernel_stages23', [1, 5]))
        stride_s23 = tuple(cnn2_config.get('stride_stages23', [1, 2]))
        padding_s23 = tuple(cnn2_config.get('padding_stages23', [0, 2]))
        
        # FC e dropout
        fc_hidden_size = cnn2_config.get('fc_hidden_size', 128)
        dropout_rate = config['model'].get('dropout_rate', 0.3)
        
        # Calcular dimensões após cada convolução
        # Stage 1
        conv1_h = (num_rows - kernel_s1[0]) // stride_s1[0] + 1
        conv1_w = (effective_size - kernel_s1[1]) // stride_s1[1] + 1
        
        # Stage 2
        conv2_h = (conv1_h + 2*padding_s23[0] - kernel_s23[0]) // stride_s23[0] + 1
        conv2_w = (conv1_w + 2*padding_s23[1] - kernel_s23[1]) // stride_s23[1] + 1
        
        # Stage 3
        conv3_h = (conv2_h + 2*padding_s23[0] - kernel_s23[0]) // stride_s23[0] + 1
        conv3_w = (conv2_w + 2*padding_s23[1] - kernel_s23[1]) // stride_s23[1] + 1
        
        # Bloco de features convolucionais
        self.features = nn.Sequential(
            # Stage 1: Agrupa tracks de genes
            nn.Conv2d(1, num_filters_s1, kernel_size=kernel_s1, stride=stride_s1),
            nn.ReLU(),
            
            # Stage 2: Processa ao longo das bases
            nn.Conv2d(num_filters_s1, num_filters_s2, 
                     kernel_size=kernel_s23, stride=stride_s23, padding=padding_s23),
            nn.ReLU(),
            
            # Stage 3: Features mais abstratas
            nn.Conv2d(num_filters_s2, num_filters_s3,
                     kernel_size=kernel_s23, stride=stride_s23, padding=padding_s23),
            nn.ReLU()
        )
        
        # Global pooling determinístico: usa MaxPool2d ou AvgPool2d com kernel fixo
        # Faz pooling apenas na dimensão da largura, mantendo a altura (genes)
        # NOTA: AdaptiveAvgPool2d não é determinístico em CUDA no backward pass
        pool_type = cnn2_config.get('global_pool_type', 'max')
        if pool_type == 'max':
            self.global_pool = nn.MaxPool2d(kernel_size=(1, conv3_w))
        else:
            self.global_pool = nn.AvgPool2d(kernel_size=(1, conv3_w))
        self.pool_type = pool_type
        
        # Dimensão após flatten
        flattened_size = num_filters_s3 * conv1_h * 1
        
        # Classificador
        self.classifier = nn.Sequential(
            nn.Linear(flattened_size, fc_hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(fc_hidden_size, num_classes)
        )
        
        # Inicializar pesos
        self._initialize_weights()
        
        # Imprimir informações do modelo
        console.print(f"[green]Modelo CNN2 criado:[/green]")
        console.print(f"  • Input shape original: {input_shape[0]} x {input_shape[1]}")
        console.print(f"  • Genes selecionados: {len(self.genes_selected)} genes → {len(self.rows_to_use)} linhas")
        console.print(f"  • Genes: {', '.join(self.genes_selected)}")
        console.print(f"  • Input shape efetivo: {num_rows} x {effective_size} (1 canal)")
        console.print(f"  • Stage 1: {num_filters_s1} filters, kernel={kernel_s1}, stride={stride_s1}")
        console.print(f"    → Output: {num_filters_s1} x {conv1_h} x {conv1_w}")
        console.print(f"  • Stage 2: {num_filters_s2} filters, kernel={kernel_s23}, stride={stride_s23}, padding={padding_s23}")
        console.print(f"    → Output: {num_filters_s2} x {conv2_h} x {conv2_w}")
        console.print(f"  • Stage 3: {num_filters_s3} filters, kernel={kernel_s23}, stride={stride_s23}, padding={padding_s23}")
        console.print(f"    → Output: {num_filters_s3} x {conv3_h} x {conv3_w}")
        console.print(f"  • Global Pool ({pool_type}): kernel=(1, {conv3_w}) → {num_filters_s3} x {conv1_h} x 1")
        console.print(f"  • Flatten size: {flattened_size}")
        console.print(f"  • FC: {flattened_size} → {fc_hidden_size} → {num_classes}")
        console.print(f"  • Dropout: {dropout_rate}")
        console.print(f"  • Total parameters: {self.count_parameters():,}")
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo logits (não probabilidades)
        """
        # Selecionar apenas as linhas dos genes escolhidos
        # [batch, num_rows_original, effective_size] -> [batch, num_rows_selected, effective_size]
        x = x[:, self.rows_to_use, :]
        
        # Adicionar dimensão de canal: [batch, num_rows_selected, effective_size] -> [batch, 1, num_rows_selected, effective_size]
        x = x.unsqueeze(1)
        
        # Features convolucionais
        x = self.features(x)
        
        # Global pooling
        x = self.global_pool(x)
        
        # Flatten
        x = x.view(x.size(0), -1)
        
        # Classificador
        logits = self.classifier(x)
        return logits
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Calcula probabilidades de classe usando softmax sobre os logits.
        
        Args:
            x: Input tensor com shape [batch, num_rows, effective_size] (2D)
               
        Returns:
            Output tensor com shape [batch, num_classes] contendo probabilidades
        """
        logits = self.forward(x)
        return torch.softmax(logits, dim=1)
    
    def count_parameters(self) -> int:
        """Conta número total de parâmetros treináveis."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
    
    def _initialize_weights(self):
        """
        Inicializa pesos usando He/Kaiming para ReLU.
        
        Estratégia:
        - Conv2D: Kaiming (He) initialization para ReLU
        - FC layers: Kaiming para camadas intermediárias, Xavier para saída
        - Bias: zeros
        """
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                # Última camada (saída) usa Xavier, outras usam Kaiming
                if m.out_features == self.num_classes:
                    nn.init.xavier_normal_(m.weight)
                else:
                    nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
        
        console.print(f"[green]✓ Pesos CNN2 inicializados:[/green]")
        console.print(f"  • Conv2D: Kaiming (He) initialization")
        console.print(f"  • FC hidden: Kaiming initialization")
        console.print(f"  • FC output: Xavier initialization")
        console.print(f"  • Bias: zeros")


# ==============================================================================
# INTERPRETABILITY: GRAD-CAM AND DEEPLIFT
# ==============================================================================

class GradCAM:
    """
    Grad-CAM (Gradient-weighted Class Activation Mapping) para CNNs.
    
    Calcula mapas de ativação que destacam regiões importantes para a predição
    de uma classe específica, usando gradientes da saída em relação às ativações
    da última camada convolucional.
    
    Referência: Selvaraju et al., 2017 - "Grad-CAM: Visual Explanations from 
    Deep Networks via Gradient-based Localization"
    
    Uso:
        gradcam = GradCAM(model, target_layer='auto')
        cam = gradcam.generate(input_tensor, target_class=None)  # None = classe predita
    """
    
    def __init__(self, model: nn.Module, target_layer: str = 'auto'):
        """
        Inicializa Grad-CAM.
        
        Args:
            model: Modelo CNN (CNNAncestryPredictor ou CNN2AncestryPredictor)
            target_layer: Camada alvo para extrair ativações
                         'auto' = detecta automaticamente a última conv
        """
        self.model = model
        self.target_layer = target_layer
        self.activations = None
        self.gradients = None
        self._hooks = []
        
        # Detectar tipo de modelo e camada alvo
        self._setup_target_layer()
    
    def _setup_target_layer(self):
        """Configura a camada alvo e registra hooks."""
        model_type = type(self.model).__name__
        
        if model_type == 'CNNAncestryPredictor':
            # CNN simples: usar conv1
            target = self.model.conv1
            self._layer_name = 'conv1'
        elif model_type == 'CNN2AncestryPredictor':
            # CNN2: usar última conv em features (índice 4 = Stage 3 Conv)
            # features = [Conv2d, ReLU, Conv2d, ReLU, Conv2d, ReLU]
            target = self.model.features[4]  # Stage 3 Conv2d
            self._layer_name = 'features[4] (Stage 3)'
        else:
            raise ValueError(f"Grad-CAM não suportado para modelo tipo: {model_type}")
        
        # Registrar hooks
        self._hooks.append(
            target.register_forward_hook(self._forward_hook)
        )
        self._hooks.append(
            target.register_full_backward_hook(self._backward_hook)
        )
        
        console.print(f"[cyan]Grad-CAM configurado para camada: {self._layer_name}[/cyan]")
    
    def _forward_hook(self, module, input, output):
        """Hook para capturar ativações no forward pass."""
        self.activations = output.detach()
    
    def _backward_hook(self, module, grad_input, grad_output):
        """Hook para capturar gradientes no backward pass."""
        self.gradients = grad_output[0].detach()
    
    def generate(self, input_tensor: torch.Tensor, target_class: Optional[int] = None) -> torch.Tensor:
        """
        Gera mapa Grad-CAM para uma entrada.
        
        Args:
            input_tensor: Tensor de entrada [batch, num_rows, effective_size]
            target_class: Classe alvo para gerar CAM. Se None, usa classe predita.
            
        Returns:
            Tensor com mapa de ativação [num_rows, effective_size] normalizado [0, 1]
        """
        # Garantir que modelo está em modo eval mas com gradientes
        was_training = self.model.training
        self.model.eval()
        
        # Forward pass
        input_tensor = input_tensor.requires_grad_(True)
        output = self.model(input_tensor)
        
        # Determinar classe alvo
        if target_class is None:
            target_class = output.argmax(dim=1).item()
        
        # Backward pass para a classe alvo
        # Estratégia inversa: propagar gradiente das outras classes para evitar
        # saturação quando a classe alvo tem probabilidade muito alta.
        # Gradiente = 0 para classe alvo, 1/(N-1) para as demais.
        # Isso captura "o que reduz as outras classes" ≈ "o que ativa a classe alvo"
        self.model.zero_grad()
        num_classes = output.shape[1]
        gradient = torch.ones_like(output) / (num_classes - 1)
        gradient[0, target_class] = 0
        output.backward(gradient=gradient, retain_graph=True)
        
        # Calcular pesos (média global dos gradientes por canal)
        # gradients shape: [batch, channels, height, width]
        weights = self.gradients.mean(dim=(2, 3), keepdim=True)  # [batch, channels, 1, 1]
        
        # Combinação ponderada das ativações
        # activations shape: [batch, channels, height, width]
        cam = (weights * self.activations).sum(dim=1, keepdim=True)  # [batch, 1, height, width]
        
        # ReLU para manter apenas contribuições positivas
        cam = torch.relu(cam)
        
        # Sem normalização - manter valores originais do CAM
        cam = cam.squeeze()  # [height, width]
        
        # Interpolar para tamanho da entrada original
        # Obter dimensão original da entrada (após seleção de genes)
        if hasattr(self.model, 'rows_to_use'):
            target_h = len(self.model.rows_to_use)
        else:
            target_h = input_tensor.shape[1]
        target_w = input_tensor.shape[2]
        
        # Redimensionar CAM para tamanho da entrada
        cam_resized = torch.nn.functional.interpolate(
            cam.unsqueeze(0).unsqueeze(0),
            size=(target_h, target_w),
            mode='bilinear',
            align_corners=False
        ).squeeze()
        
        # Restaurar estado do modelo
        if was_training:
            self.model.train()
        
        return cam_resized.cpu(), target_class
    
    def remove_hooks(self):
        """Remove todos os hooks registrados."""
        for hook in self._hooks:
            hook.remove()
        self._hooks = []


class DeepLIFT:
    """
    DeepLIFT (Deep Learning Important Features) para CNNs.
    
    Calcula atribuições de importância para cada feature de entrada,
    propagando diferenças de ativação (em relação a um baseline) através
    da rede usando regras de contribuição específicas.
    
    Implementação simplificada que usa gradientes × (input - baseline),
    uma aproximação válida para redes com ReLU (Rescale Rule).
    
    Referência: Shrikumar et al., 2017 - "Learning Important Features 
    Through Propagating Activation Differences"
    
    Uso:
        deeplift = DeepLIFT(model)
        attributions = deeplift.generate(input_tensor, baseline='zeros')
    """
    
    def __init__(self, model: nn.Module):
        """
        Inicializa DeepLIFT.
        
        Args:
            model: Modelo (NN, CNN ou CNN2)
        """
        self.model = model
        self._baseline_cache = None
        self._class_mean_cache: Dict[int, torch.Tensor] = {}  # Cache para média de atribuições por classe
        self._class_input_mean_cache: Dict[int, Tuple[torch.Tensor, int]] = {}  # Cache para (média das entradas, num_samples) por classe
    
    def _get_baseline(self, input_tensor: torch.Tensor, baseline_type: str,
                      dataset: Optional[Any] = None) -> torch.Tensor:
        """
        Obtém tensor baseline.
        
        Args:
            input_tensor: Tensor de entrada para obter shape e device
            baseline_type: 'zeros' ou 'mean'
            dataset: Dataset para calcular média (necessário se baseline_type='mean')
            
        Returns:
            Tensor baseline com mesmo shape que input_tensor
        """
        if baseline_type == 'zeros':
            return torch.zeros_like(input_tensor)
        
        elif baseline_type == 'mean':
            if self._baseline_cache is not None:
                return self._baseline_cache.to(input_tensor.device)
            
            if dataset is None:
                console.print("[yellow]⚠ Dataset não fornecido para baseline='mean', usando zeros[/yellow]")
                return torch.zeros_like(input_tensor)
            
            # Calcular média de todas as amostras
            console.print("\n\n[red]Calculando baseline (média do dataset)...[/red]\n\n")
            all_samples = []
            for i in range(min(len(dataset), 10000)):  # Limitar a 10000 amostras
                sample, _, _ = dataset[i]
                all_samples.append(sample)
            
            mean_baseline = torch.stack(all_samples).mean(dim=0)
            self._baseline_cache = mean_baseline
            return mean_baseline.unsqueeze(0).to(input_tensor.device)
        
        else:
            raise ValueError(f"Baseline type desconhecido: {baseline_type}")
    
    def generate(self, input_tensor: torch.Tensor, target_class: Optional[int] = None,
                 baseline_type: str = 'zeros', dataset: Optional[Any] = None) -> torch.Tensor:
        """
        Gera atribuições DeepLIFT para uma entrada.
        
        Usa a aproximação gradient × (input - baseline), que é equivalente
        ao DeepLIFT Rescale Rule para redes com ReLU.
        
        Args:
            input_tensor: Tensor de entrada [batch, num_rows, effective_size]
            target_class: Classe alvo. Se None, usa classe predita.
            baseline_type: 'zeros' ou 'mean'
            dataset: Dataset para calcular média (necessário se baseline_type='mean')
            
        Returns:
            Tensor de atribuições [num_rows, effective_size]
        """
        # Garantir modo eval
        was_training = self.model.training
        self.model.eval()
        
        # Obter baseline
        baseline = self._get_baseline(input_tensor, baseline_type, dataset)
        
        # Calcular diferença input - baseline
        delta = input_tensor - baseline
        
        # Habilitar gradientes
        input_tensor = input_tensor.clone().requires_grad_(True)
        
        # Forward pass
        output = self.model(input_tensor)
        
        # Determinar classe alvo
        if target_class is None:
            target_class = output.argmax(dim=1).item()
        
        # Backward pass
        self.model.zero_grad()
        one_hot = torch.zeros_like(output)
        one_hot[0, target_class] = 1
        output.backward(gradient=one_hot)
        
        # Atribuições = gradiente × delta (aproximação do DeepLIFT Rescale Rule)
        gradients = input_tensor.grad.detach()
        attributions = gradients * delta
        
        # Remover dimensão de batch
        attributions = attributions.squeeze(0)
        
        # Se modelo tem seleção de genes, aplicar máscara
        if hasattr(self.model, 'rows_to_use'):
            # Manter apenas as linhas usadas pelo modelo
            rows_to_use = self.model.rows_to_use
            # Criar tensor com zeros e preencher apenas as linhas usadas
            full_attr = torch.zeros_like(attributions)
            full_attr[rows_to_use, :] = attributions[rows_to_use, :]
            attributions = full_attr
        
        # Restaurar estado
        if was_training:
            self.model.train()
        
        return attributions.cpu(), target_class
    
    def generate_class_mean(self, target_class_idx: int, dataset: Any,
                            baseline_type: str = 'zeros') -> Tuple[torch.Tensor, torch.Tensor, int]:
        """
        Gera a média das atribuições DeepLIFT e das entradas para todas as amostras de uma classe.
        
        Args:
            target_class_idx: Índice da classe alvo
            dataset: Dataset para buscar amostras
            baseline_type: 'zeros' ou 'mean'
            
        Returns:
            Tupla (mean_attributions, mean_input, num_samples):
            - mean_attributions: Tensor com a média das atribuições [num_rows, effective_size]
            - mean_input: Tensor com a média das entradas [num_rows, effective_size]
            - num_samples: Número de amostras usadas
        """
        # Verificar cache
        if target_class_idx in self._class_mean_cache and target_class_idx in self._class_input_mean_cache:
            mean_input, num_samples = self._class_input_mean_cache[target_class_idx]
            return self._class_mean_cache[target_class_idx], mean_input, num_samples
        
        console.print(f"[cyan]Calculando média DeepLIFT para classe {target_class_idx}...[/cyan]")
        
        # Coletar todas as amostras da classe alvo
        class_samples = []
        for i in range(len(dataset)):
            sample, target, _ = dataset[i]
            if target == target_class_idx:
                class_samples.append(sample)
        
        if len(class_samples) == 0:
            console.print(f"[yellow]⚠ Nenhuma amostra encontrada para classe {target_class_idx}[/yellow]")
            return None, None, 0
        
        num_samples = len(class_samples)
        console.print(f"[cyan]  → {num_samples} amostras encontradas[/cyan]")
        
        # Calcular média das entradas
        mean_input = torch.stack(class_samples).mean(dim=0)
        
        # Calcular DeepLIFT para cada amostra
        all_attributions = []
        device = next(self.model.parameters()).device
        
        for i, sample in enumerate(class_samples):
            # Preparar entrada
            input_tensor = sample.unsqueeze(0).to(device)
            
            # Calcular atribuições para esta amostra
            attr, _ = self.generate(input_tensor, target_class=target_class_idx,
                                   baseline_type=baseline_type, dataset=dataset)
            all_attributions.append(attr)
            
            # Progress (a cada 10 amostras)
            if (i + 1) % 10 == 0:
                console.print(f"[dim]  → Processadas {i + 1}/{num_samples} amostras[/dim]")
        
        # Calcular média das atribuições
        mean_attributions = torch.stack(all_attributions).mean(dim=0)
        
        # Armazenar nos caches
        self._class_mean_cache[target_class_idx] = mean_attributions
        self._class_input_mean_cache[target_class_idx] = (mean_input, num_samples)
        
        console.print(f"[green]✓ Média DeepLIFT calculada para classe {target_class_idx}[/green]")
        
        return mean_attributions, mean_input, num_samples


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
        experiment_dir: Path,
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
            experiment_dir: Diretório do experimento
            wandb_run: Run do W&B (opcional)
        """
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.config = config
        self.device = device
        self.experiment_dir = experiment_dir
        self.wandb_run = wandb_run
        
        # Configurações de visualização/debug
        self.enable_visualization = config.get('debug', {}).get('enable_visualization', False)
        self.max_samples_per_epoch = config.get('debug', {}).get('max_samples_per_epoch', None)
        
        if self.enable_visualization:
            console.print(f"[yellow]⚠ Modo de visualização ativado! Batch size será forçado para 1.[/yellow]")
            console.print(f"[yellow]Pressione qualquer tecla na janela do gráfico para continuar. Pressione 'q' para sair.[/yellow]")
            plt.ion()  # Modo interativo
            self._key_pressed = False
            self._quit_requested = False
        
        # Configurar otimizador
        optimizer_type = config['training']['optimizer'].lower()
        lr = config['training']['learning_rate']
        weight_decay = config['training'].get('weight_decay', 0.0)
        
        if optimizer_type == 'adam':
            self.optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
        elif optimizer_type == 'adamw':
            self.optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
        elif optimizer_type == 'sgd':
            self.optimizer = optim.SGD(model.parameters(), lr=lr, momentum=0.9, weight_decay=weight_decay)
        else:
            raise ValueError(f"Otimizador não suportado: {optimizer_type}")
        
        if weight_decay > 0:
            console.print(f"[green]✓ Weight decay (L2): {weight_decay}[/green]")
        
        # Configurar learning rate scheduler
        self.scheduler = None
        if config['training'].get('lr_scheduler', {}).get('enabled', False):
            scheduler_config = config['training']['lr_scheduler']
            scheduler_type = scheduler_config.get('type', 'plateau').lower()
            
            if scheduler_type == 'plateau':
                self.scheduler = optim.lr_scheduler.ReduceLROnPlateau(
                    self.optimizer,
                    mode=scheduler_config.get('mode', 'min'),
                    factor=scheduler_config.get('factor', 0.5),
                    patience=scheduler_config.get('patience', 10),
                    min_lr=scheduler_config.get('min_lr', 1e-6)
                )
                console.print(f"[green]✓ LR Scheduler: ReduceLROnPlateau (mode={scheduler_config.get('mode', 'min')}, patience={scheduler_config.get('patience', 10)}, factor={scheduler_config.get('factor', 0.5)})[/green]")
            
            elif scheduler_type == 'step':
                self.scheduler = optim.lr_scheduler.StepLR(
                    self.optimizer,
                    step_size=scheduler_config.get('step_size', 30),
                    gamma=scheduler_config.get('gamma', 0.1)
                )
                console.print(f"[green]✓ LR Scheduler: StepLR (step_size={scheduler_config.get('step_size', 30)})[/green]")
            
            elif scheduler_type == 'cosine':
                T_max = scheduler_config.get('T_max', config['training']['num_epochs'])
                self.scheduler = optim.lr_scheduler.CosineAnnealingLR(
                    self.optimizer,
                    T_max=T_max,
                    eta_min=scheduler_config.get('eta_min', 1e-6)
                )
                console.print(f"[green]✓ LR Scheduler: CosineAnnealingLR (T_max={T_max})[/green]")
            
            elif scheduler_type == 'exponential':
                self.scheduler = optim.lr_scheduler.ExponentialLR(
                    self.optimizer,
                    gamma=scheduler_config.get('gamma', 0.95)
                )
                console.print(f"[green]✓ LR Scheduler: ExponentialLR (gamma={scheduler_config.get('gamma', 0.95)})[/green]")
            
            elif scheduler_type == 'multistep':
                self.scheduler = optim.lr_scheduler.MultiStepLR(
                    self.optimizer,
                    milestones=scheduler_config.get('milestones', [30, 60, 90]),
                    gamma=scheduler_config.get('gamma', 0.1)
                )
                console.print(f"[green]✓ LR Scheduler: MultiStepLR (milestones={scheduler_config.get('milestones', [30, 60, 90])})[/green]")
            
            elif scheduler_type == 'cosine_warm_restarts':
                T_0 = scheduler_config.get('T_0', 50)
                T_mult = scheduler_config.get('T_mult', 1)
                self.scheduler = optim.lr_scheduler.CosineAnnealingWarmRestarts(
                    self.optimizer,
                    T_0=T_0,
                    T_mult=T_mult,
                    eta_min=scheduler_config.get('eta_min', 1e-6)
                )
                console.print(f"[green]✓ LR Scheduler: CosineAnnealingWarmRestarts (T_0={T_0}, T_mult={T_mult})[/green]")
            
            else:
                console.print(f"[yellow]⚠ Scheduler type '{scheduler_type}' não reconhecido, continuando sem scheduler[/yellow]")
                self.scheduler = None
        
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
            'train_accuracy': [],
            'val_loss': [],
            'val_accuracy': [],
            'epoch': []
        }
        
        self.best_val_loss = float('inf')
        self.best_val_accuracy = 0.0
        self.best_val_loss_at_best_accuracy = float('inf')  # Para atualizar best_accuracy quando loss melhora
    
    def _on_key_press(self, event):
        """Callback para capturar tecla pressionada na janela do gráfico."""
        self._key_pressed = True
        if event.key == 'q':
            self._quit_requested = True
        plt.close()
    
    def _visualize_sample(self, features: torch.Tensor, targets: torch.Tensor, 
                         outputs: torch.Tensor, sample_idx: int, sample_id: str, epoch: int):
        """
        Visualiza uma amostra de entrada e suas predições.
        
        Args:
            features: Tensor de entrada with shape (1, num_rows, effective_size) for 2D
            targets: True target (1,)
            outputs: Network output (1, num_classes) - logits
            sample_idx: Sample index in the split
            sample_id: Sample ID (e.g., HG00138)
            epoch: Epoch number
        """
        # Convert to CPU and numpy
        features_cpu = features.cpu().detach()
        target_idx = targets.cpu().item()
        # Apply softmax over logits to get probabilities
        output_probs = torch.softmax(outputs, dim=1).cpu().detach().numpy()[0]
        predicted_idx = output_probs.argmax()
        
        # Get class names (if available)
        class_names = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']  # For superpopulation
        if target_idx < len(class_names):
            target_name = class_names[target_idx]
            predicted_name = class_names[predicted_idx]
        else:
            target_name = f"Class {target_idx}"
            predicted_name = f"Class {predicted_idx}"
        
        # Get gene names from config
        genes_to_use = self.config['dataset_input'].get('genes_to_use', None)
        GENE_ORDER = self.config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1",
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        if genes_to_use:
            # Maintain gene order from dataset
            gene_names = [gene for gene in GENE_ORDER if gene in genes_to_use]
        else:
            gene_names = GENE_ORDER
        tracks_per_gene = 6
        
        # Create figure
        plt.clf()
        fig = plt.gcf()
        fig.set_size_inches(16, 8)
        
        # Plot input features
        ax1 = plt.subplot(2, 1, 1)
        
        # Detect if input is 2D or 1D
        if features_cpu.ndim == 3 and features_cpu.shape[0] == 1:
            # 2D input: [1, num_rows, effective_size]
            img_data = features_cpu[0].numpy()  # [num_rows, effective_size]
            
            # Rescale for visualization
            viz_height = self.config.get('debug', {}).get('visualization', {}).get('height', 300)
            viz_width = self.config.get('debug', {}).get('visualization', {}).get('width', 600)
            
            # Calculate zoom factors
            zoom_factors = (viz_height / img_data.shape[0], viz_width / img_data.shape[1])
            img_resized = ndimage.zoom(img_data, zoom_factors, order=1)  # order=1 = bilinear
            
            # Normalize for visualization (0=black, 1=white)
            img_min, img_max = img_resized.min(), img_resized.max()
            if img_max > img_min:
                img_normalized = (img_resized - img_min) / (img_max - img_min)
            else:
                img_normalized = np.zeros_like(img_resized)
            
            # Plot as image
            plt.imshow(img_normalized, cmap='gray', aspect='auto', interpolation='nearest')
            plt.xlabel('Genomic Position (rescaled)', fontsize=12)
            plt.title(f'Epoch {epoch + 1} | Sample {sample_id} ({target_name}) | Input 2D ({img_data.shape[0]}x{img_data.shape[1]} → {viz_height}x{viz_width})', 
                     fontsize=14, fontweight='bold')
            plt.colorbar(label='Normalized Value')
            
            # Configure Y axis with gene names
            num_genes = len(gene_names)
            if img_data.shape[0] == num_genes * tracks_per_gene:
                pixels_per_row = viz_height / img_data.shape[0]
                
                # Major ticks at gene boundaries (before and after each gene)
                y_major_ticks = [i * tracks_per_gene * pixels_per_row for i in range(num_genes + 1)]
                ax1.set_yticks(y_major_ticks)
                ax1.set_yticklabels([''] * len(y_major_ticks))  # No labels on boundary ticks
                ax1.tick_params(axis='y', which='major', length=8, width=0.8)
                
                # Minor ticks at center of each gene for labels
                y_minor_ticks = [(i * tracks_per_gene + tracks_per_gene / 2) * pixels_per_row 
                                 for i in range(num_genes)]
                ax1.set_yticks(y_minor_ticks, minor=True)
                ax1.set_yticklabels(gene_names, minor=True, fontsize=10)
                ax1.tick_params(axis='y', which='minor', length=0)  # Hide minor tick marks
                
                ax1.set_ylabel('Genes', fontsize=12)
            else:
                ax1.set_ylabel('Tracks (rescaled)', fontsize=12)
        else:
            # 1D input (fallback for backwards compatibility)
            features_np = features_cpu.numpy().flatten()
            plt.plot(features_np, linewidth=0.5, alpha=0.7)
            plt.xlabel('Feature Index', fontsize=12)
            plt.ylabel('Feature Value', fontsize=12)
            plt.title(f'Epoch {epoch + 1} | Sample {sample_id} | Input Features (n={len(features_np)})', 
                     fontsize=14, fontweight='bold')
            plt.grid(True, alpha=0.3)
            
        # Plot output probabilities
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
        
        # Text with prediction and target
        correct = "✓ CORRECT" if predicted_idx == target_idx else "✗ WRONG"
        color = 'green' if predicted_idx == target_idx else 'red'
        result_text = (f'{correct}\n'
                      f'Target: {target_name} (class {target_idx})\n'
                      f'Predicted: {predicted_name} (class {predicted_idx}, prob={output_probs[predicted_idx]:.3f})')
        
        plt.text(0.98, 0.98, result_text, transform=plt.gca().transAxes,
                fontsize=11, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
        
        plt.tight_layout()
        
        # Connect key callback
        self._key_pressed = False
        cid = fig.canvas.mpl_connect('key_press_event', self._on_key_press)
        
        # Show and wait for key
        plt.show(block=True)
        
        # Disconnect callback
        fig.canvas.mpl_disconnect(cid)
    
    def train_epoch(self, epoch: int) -> float:
        """
        Treina por uma época.
        
        Returns:
            Loss média da época
        """
        self.model.train()
        total_loss = 0.0
        total_samples = 0
        
        # Debug: salvar detalhes na época 50
        debug_epoch_50 = (epoch + 1 == 50)
        debug_data = []
        
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
            
            for batch_idx, (features, targets, indices) in enumerate(self.train_loader):
                # Verificar limite de amostras para visualização
                if self.enable_visualization and self.max_samples_per_epoch is not None:
                    if batch_idx >= self.max_samples_per_epoch:
                        break
                
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                
                batch_size = targets.size(0)
                
                # Forward pass
                self.optimizer.zero_grad()
                outputs = self.model(features)
                loss = self.criterion(outputs, targets)
                
                # Debug: coletar dados individuais (antes do backward)
                if debug_epoch_50:
                    with torch.no_grad():
                        for i in range(batch_size):
                            target_val = targets[i].item()
                            output_logits = outputs[i].cpu().numpy()
                            # Calcular loss individual
                            individual_loss = torch.nn.functional.cross_entropy(
                                outputs[i:i+1], targets[i:i+1], reduction='none'
                            ).item()
                            debug_data.append({
                                'target': target_val,
                                'logits': output_logits,
                                'loss': individual_loss
                            })
                
                # Backward pass
                loss.backward()
                self.optimizer.step()
                
                # Acumular loss ponderada pelo tamanho do batch
                total_loss += loss.item() * batch_size
                total_samples += batch_size
                
                # Visualização interativa (se habilitada)
                if self.enable_visualization:
                    sample_idx = indices[0].item()  # batch_size=1 na visualização
                    sample_id = self.train_loader.dataset.get_sample_id(sample_idx)
                    self._visualize_sample(features, targets, outputs, sample_idx, sample_id, epoch)
                    # Check if user requested to quit
                    if self._quit_requested:
                        console.print("[yellow]⚠ Quit requested by user (pressed 'q')[/yellow]")
                        import sys
                        sys.exit(0)
                
                # Log no W&B
                if self.wandb_run and batch_idx % self.config['wandb']['log_frequency'] == 0:
                    self.wandb_run.log({
                        'batch_loss': loss.item(),
                        'epoch': epoch,
                        'batch': batch_idx
                    })
                
                progress.update(task, advance=1)
        
        avg_loss = total_loss / total_samples if total_samples > 0 else 0.0
        
        # Salvar debug na época 50
        if debug_epoch_50:
            with open('train_loss.txt', 'w') as f:
                f.write("Training Loss Debug - Época 50\n")
                f.write("=" * 80 + "\n")
                f.write(f"Total samples: {total_samples}\n")
                f.write(f"Average loss: {avg_loss:.6f}\n")
                f.write("=" * 80 + "\n\n")
                for idx, data in enumerate(debug_data):
                    logits = data['logits']
                    # Calcular softmax dos logits
                    logits_tensor = torch.tensor(logits) if not isinstance(logits, torch.Tensor) else logits
                    softmax_probs = torch.softmax(logits_tensor, dim=0).numpy()
                    
                    # Formatar valores com 3 casas decimais, alinhados (7 chars: sinal + 1 digito + ponto + 3 decimais)
                    logits_str = "[" + ", ".join(f"{v:7.3f}" for v in logits) + "]"
                    softmax_str = "[" + ", ".join(f"{v:7.3f}" for v in softmax_probs) + "]"
                    
                    f.write(f"Sample {idx+1}:\n")
                    f.write(f"  Target:  {data['target']}\n")
                    f.write(f"  Logits:  {logits_str}\n")
                    f.write(f"  Softmax: {softmax_str}\n")
                    f.write(f"  Loss:    {data['loss']:.6f}\n\n")
            console.print(f"[yellow]Debug: Training loss detalhada salva em train_loss.txt[/yellow]")
        
        return avg_loss
    
    def evaluate_train_accuracy(self) -> float:
        """
        Avalia acurácia no conjunto de treino sem atualizar pesos.
        
        O modelo é colocado em modo eval() (sem dropout, BatchNorm em inference),
        garantindo avaliação justa e consistente com a validação.
        
        Chamado apenas quando for logar métricas (não a cada época).
        
        Returns:
            Acurácia em porcentagem
        """
        self.model.eval()
        correct = 0
        total = 0
        
        with torch.no_grad():
            for features, targets, indices in self.train_loader:
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                
                outputs = self.model(features)
                _, predicted = torch.max(outputs.data, 1)
                total += targets.size(0)
                correct += (predicted == targets).sum().item()
        
        accuracy = correct / total if total > 0 else 0.0
        return accuracy
    
    def validate(self, epoch: int) -> Tuple[float, float]:
        """
        Valida o modelo.
        
        Returns:
            Tupla (loss, accuracy)
        """
        self.model.eval()
        total_loss = 0.0
        total_samples = 0
        all_predictions = []
        all_targets = []
        
        # Debug: salvar detalhes na época 50
        debug_epoch_50 = (epoch + 1 == 50)
        debug_data = []
        
        with torch.no_grad():
            for features, targets, indices in self.val_loader:
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                
                batch_size = targets.size(0)
                
                outputs = self.model(features)
                loss = self.criterion(outputs, targets)
                
                # Debug: coletar dados individuais
                if debug_epoch_50:
                    for i in range(batch_size):
                        target_val = targets[i].item()
                        output_logits = outputs[i].cpu().numpy()
                        # Calcular loss individual
                        individual_loss = torch.nn.functional.cross_entropy(
                            outputs[i:i+1], targets[i:i+1], reduction='none'
                        ).item()
                        debug_data.append({
                            'target': target_val,
                            'logits': output_logits,
                            'loss': individual_loss
                        })
                
                # Acumular loss ponderada pelo tamanho do batch
                total_loss += loss.item() * batch_size
                total_samples += batch_size
                
                # Para classificação, calcular accuracy
                if self.config['output']['prediction_target'] != 'frog_likelihood':
                    predictions = torch.argmax(outputs, dim=1)
                    all_predictions.extend(predictions.cpu().numpy())
                    all_targets.extend(targets.cpu().numpy())
        
        avg_loss = total_loss / total_samples if total_samples > 0 else 0.0
        
        # Salvar debug na época 50
        if debug_epoch_50:
            with open('val_loss.txt', 'w') as f:
                f.write("Validation Loss Debug - Época 50\n")
                f.write("=" * 80 + "\n")
                f.write(f"Total samples: {total_samples}\n")
                f.write(f"Average loss: {avg_loss:.6f}\n")
                f.write("=" * 80 + "\n\n")
                for idx, data in enumerate(debug_data):
                    f.write(f"Sample {idx+1}:\n")
                    f.write(f"  Target: {data['target']}\n")
                    f.write(f"  Logits: {data['logits']}\n")
                    f.write(f"  Loss: {data['loss']:.6f}\n\n")
            console.print(f"[yellow]Debug: Validation loss detalhada salva em val_loss.txt[/yellow]")
        
        # Calcular accuracy
        if len(all_predictions) > 0:
            accuracy = accuracy_score(all_targets, all_predictions)
        else:
            accuracy = 0.0
        
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
        
        # Imprimir arquitetura da rede
        console.print("\n[bold cyan]═══ ARQUITETURA DA REDE ═══[/bold cyan]")
        console.print(self.model)
        console.print(f"\n[bold green]Total de parâmetros treináveis: {self.model.count_parameters():,}[/bold green]")
        console.print()
        
        try:
            for epoch in range(num_epochs):
                epoch_start = time.time()
                # Verificar se houve interrupção (CTRL+C)
                if interrupt_state.interrupted:
                    console.print("\n[yellow]⚠ Treinamento interrompido pelo usuário (CTRL+C)[/yellow]")
                    break
                
                # Treinar
                process = psutil.Process(os.getpid())
                train_start = time.time()
                train_loss = self.train_epoch(epoch)
                train_duration = time.time() - train_start
                rss_mb = process.memory_info().rss / 1024 / 1024
                console.print(f"[dim cyan]Tempo treino: {train_duration:.2f}s | RSS: {rss_mb:.1f} MB[/dim cyan]")
                
                # Validar
                if (epoch + 1) % val_frequency == 0:
                    # Avaliar acurácia de treino (apenas quando validar)
                    eval_start = time.time()
                    train_accuracy = self.evaluate_train_accuracy()
                    eval_duration = time.time() - eval_start
                    rss_mb = process.memory_info().rss / 1024 / 1024
                    console.print(f"[dim cyan]Tempo avaliação treino: {eval_duration:.2f}s | RSS: {rss_mb:.1f} MB[/dim cyan]")

                    # Liberar memória dos datasets
                    # if hasattr(self.train_loader.dataset, 'unload_data'):
                    #     self.train_loader.dataset.unload_data()
                    
                    val_start = time.time()
                    val_loss, val_accuracy = self.validate(epoch)
                    val_duration = time.time() - val_start
                    rss_mb = process.memory_info().rss / 1024 / 1024
                    console.print(f"[dim cyan]Tempo validação: {val_duration:.2f}s | RSS: {rss_mb:.1f} MB[/dim cyan]")
                    
                    # Liberar memória dos datasets
                    # if hasattr(self.val_loader.dataset, 'unload_data'):
                    #     self.val_loader.dataset.unload_data()
                        
                    # Obter learning rate atual
                    current_lr = self.optimizer.param_groups[0]['lr']
                    
                    # Imprimir tudo em uma linha
                    console.print(
                        f"[cyan]Época {epoch + 1}:[/cyan] "
                        f"Train Acc={train_accuracy * 100:.2f}% | "
                        f"Val Acc={val_accuracy * 100:.2f}% | "
                        f"Train Loss={train_loss:.4f} | "
                        f"Val Loss={val_loss:.4f} | "
                        f"LR={current_lr:.2e}"
                    )
                    
                    # Salvar histórico
                    self.history['train_loss'].append(train_loss)
                    self.history['train_accuracy'].append(train_accuracy)
                    self.history['val_loss'].append(val_loss)
                    self.history['val_accuracy'].append(val_accuracy)
                    self.history['epoch'].append(epoch + 1)
                    
                    # Log no W&B
                    if self.wandb_run:
                        log_dict = {
                            'epoch': epoch + 1,
                            'train_loss': train_loss,
                            'train_accuracy': train_accuracy,
                            'val_loss': val_loss,
                            'val_accuracy': val_accuracy,
                            'learning_rate': current_lr
                        }
                        self.wandb_run.log(log_dict)
                    
                    # Atualizar e salvar melhores métricas (apenas se habilitado)
                    save_during_training = self.config['checkpointing'].get('save_during_training', True)
                    
                    if val_loss < self.best_val_loss:
                        self.best_val_loss = val_loss
                        if save_during_training:
                            self.save_checkpoint(epoch, 'best_loss')
                    
                    # Salvar best_accuracy se:
                    # 1. Accuracy melhorou, OU
                    # 2. Accuracy igual à melhor mas loss menor (logits mais confiantes)
                    if val_accuracy > self.best_val_accuracy:
                        self.best_val_accuracy = val_accuracy
                        self.best_val_loss_at_best_accuracy = val_loss
                        if save_during_training:
                            self.save_checkpoint(epoch, 'best_accuracy')
                    elif val_accuracy == self.best_val_accuracy and val_loss < self.best_val_loss_at_best_accuracy:
                        self.best_val_loss_at_best_accuracy = val_loss
                        if save_during_training:
                            self.save_checkpoint(epoch, 'best_accuracy')
                            console.print(f"[green]  ✓ best_accuracy atualizado (mesma acc, loss menor: {val_loss:.6f})[/green]")
                    
                    # Atualizar ReduceLROnPlateau (precisa da métrica, só quando valida)
                    if self.scheduler is not None and isinstance(self.scheduler, optim.lr_scheduler.ReduceLROnPlateau):
                        scheduler_config = self.config['training']['lr_scheduler']
                        if scheduler_config.get('mode', 'min') == 'min':
                            self.scheduler.step(val_loss)
                        else:  # mode == 'max'
                            self.scheduler.step(val_accuracy)
                
                # Atualizar schedulers baseados em época (TODA época, fora do bloco de validação)
                if self.scheduler is not None and not isinstance(self.scheduler, optim.lr_scheduler.ReduceLROnPlateau):
                    self.scheduler.step()
                            
                # Salvar checkpoint periódico (apenas se habilitado)
                save_during_training = self.config['checkpointing'].get('save_during_training', True)
                if save_during_training and (epoch + 1) % save_frequency == 0:
                    self.save_checkpoint(epoch, f'epoch_{epoch + 1}')
                
                epoch_duration = time.time() - epoch_start
                rss_mb = process.memory_info().rss / 1024 / 1024
                console.print(f"[dim cyan]Tempo total da época {epoch + 1}: {epoch_duration:.2f}s | RSS: {rss_mb:.1f} MB[/dim cyan]")
        
        except KeyboardInterrupt:
            console.print("\n[yellow]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/yellow]")
            console.print("[yellow]⚠ CTRL+C detectado - Finalizando treino graciosamente...[/yellow]")
            console.print("[yellow]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/yellow]")
            # Salvar checkpoint de interrupção
            if 'epoch' in locals():
                console.print(f"[yellow]Salvando checkpoint da época {epoch + 1}...[/yellow]")
                self.save_checkpoint(epoch, 'interrupted')
            console.print("[yellow]Treinamento interrompido. Checkpoint salvo.[/yellow]")
            # Marcar como interrompido (mas continuar para executar testes)
            interrupt_state.interrupted = True
        
        # Salvar checkpoint final (sempre)
        console.print("[yellow]Salvando checkpoint final...[/yellow]")
        self.save_checkpoint(num_epochs - 1, 'final')
        
        console.print("[bold green]✓ Treinamento concluído![/bold green]")
        return self.history
    
    def save_checkpoint(self, epoch: int, name: str):
        """Salva checkpoint do modelo."""
        checkpoint_dir = self.experiment_dir / 'models'
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
        wandb_run: Optional[Any] = None,
        dataset_name: str = "Test"
    ):
        """Inicializa tester.
        
        Args:
            model: Modelo a testar
            test_loader: DataLoader com dados
            dataset: Dataset completo
            config: Configuração
            device: Device (CPU ou GPU)
            wandb_run: Run do W&B (opcional)
            dataset_name: Nome do conjunto de dados sendo testado (ex: "Teste", "Treino", "Validação")
        """
        self.model = model
        self.test_loader = test_loader
        self.dataset = dataset
        self.config = config
        self.device = device
        self.wandb_run = wandb_run
        self.dataset_name = dataset_name
        
        # Debug/visualização
        self.enable_visualization = config.get('debug', {}).get('enable_visualization', False)
        self.max_samples = config.get('debug', {}).get('max_samples_per_epoch', None)
        
        # Configuração de interpretabilidade
        interp_config = config.get('debug', {}).get('interpretability', {})
        self.interpretability_enabled = (
            self.enable_visualization and 
            interp_config.get('enabled', False)
        )
        self.interp_method = interp_config.get('method', 'gradcam')
        self.interp_save_images = interp_config.get('save_images', True)
        self.interp_output_dir = interp_config.get('output_dir', 'interpretability_results')
        self.deeplift_baseline = interp_config.get('deeplift', {}).get('baseline', 'zeros')
        self.deeplift_target_class = interp_config.get('deeplift', {}).get('target_class', 'predicted')
        self.gradcam_target_class = interp_config.get('gradcam', {}).get('target_class', 'predicted')
        
        # Inicializar objetos de interpretabilidade
        self.gradcam = None
        self.deeplift = None
        
        if self.interpretability_enabled:
            model_type = type(model).__name__
            
            if model_type in ['CNNAncestryPredictor', 'CNN2AncestryPredictor']:
                if self.interp_method in ['gradcam', 'both']:
                    try:
                        self.gradcam = GradCAM(model)
                        console.print("[green]✓ Grad-CAM inicializado[/green]")
                    except Exception as e:
                        console.print(f"[yellow]⚠ Erro ao inicializar Grad-CAM: {e}[/yellow]")
                
                if self.interp_method in ['deeplift', 'both']:
                    self.deeplift = DeepLIFT(model)
                    console.print("[green]✓ DeepLIFT inicializado[/green]")
            else:
                console.print(f"[yellow]⚠ Interpretabilidade não suportada para {model_type}[/yellow]")
                self.interpretability_enabled = False
        
        # Forçar modo interativo do matplotlib se visualização habilitada
        if self.enable_visualization:
            plt.ion()
            self._key_pressed = False
            self._quit_requested = False
    
    def _on_key_press(self, event):
        """Callback para detectar tecla pressionada."""
        self._key_pressed = True
        if event.key == 'q':
            self._quit_requested = True
        plt.close()
    
    def _print_sample_info(self, sample_idx: int, sample_id: str, 
                          targets: torch.Tensor, outputs: torch.Tensor):
        """
        Imprime informações da amostra no terminal.
        
        Args:
            sample_idx: Índice da amostra no split
            sample_id: ID da amostra (e.g., HG00138)
            targets: Target verdadeiro (1,)
            outputs: Saída da rede (1, num_classes) - logits
        """
        target_idx = targets.cpu().item()
        logits = outputs.cpu().detach().numpy()[0]
        softmax_probs = torch.softmax(outputs, dim=1).cpu().detach().numpy()[0]
        
        # Formatar valores com 3 casas decimais, alinhados
        logits_str = "[" + ", ".join(f"{v:7.3f}" for v in logits) + "]"
        softmax_str = "[" + ", ".join(f"{v:7.3f}" for v in softmax_probs) + "]"
        
        # Obter nomes das classes
        class_names = self.config['output'].get('known_classes') or ["AFR", "AMR", "EAS", "EUR", "SAS"]
        target_name = class_names[target_idx] if target_idx < len(class_names) else str(target_idx)
        predicted_idx = softmax_probs.argmax()
        predicted_name = class_names[predicted_idx] if predicted_idx < len(class_names) else str(predicted_idx)
        correct = "✓" if target_idx == predicted_idx else "✗"
        
        console.print(f"\n[bold cyan]{'='*80}[/bold cyan]")
        console.print(f"[bold]Sample {sample_idx + 1}: {sample_id}[/bold]")
        console.print(f"  Target:    {target_idx} ({target_name})")
        console.print(f"  Predicted: {predicted_idx} ({predicted_name}) {correct}")
        console.print(f"  Logits:    {logits_str}")
        console.print(f"  Softmax:   {softmax_str}")
        console.print(f"[bold cyan]{'='*80}[/bold cyan]")
    
    def _visualize_sample(self, features: torch.Tensor, targets: torch.Tensor, 
                         outputs: torch.Tensor, sample_idx: int, sample_id: str):
        """
        Visualizes a test sample and its predictions.
        
        Args:
            features: Input tensor with shape (1, num_rows, effective_size) for 2D
            targets: True target (1,)
            outputs: Network output (1, num_classes) - logits
            sample_idx: Sample index in the split
            sample_id: Sample ID (e.g., HG00138)
        """
        # Convert to CPU and numpy
        features_cpu = features.cpu().detach()
        target_idx = targets.cpu().item()
        # Apply softmax over logits to get probabilities
        output_probs = torch.softmax(outputs, dim=1).cpu().detach().numpy()[0]
        predicted_idx = output_probs.argmax()
        
        # Get class names (if available)
        class_names = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']  # For superpopulation
        if target_idx < len(class_names):
            target_name = class_names[target_idx]
            predicted_name = class_names[predicted_idx]
        else:
            target_name = f"Class {target_idx}"
            predicted_name = f"Class {predicted_idx}"
        
        # Get gene names from config
        genes_to_use = self.config['dataset_input'].get('genes_to_use', None)
        GENE_ORDER = self.config['dataset_input'].get('gene_order', 
            ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1",
             "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"])
        if genes_to_use:
            # Maintain gene order from dataset
            gene_names = [gene for gene in GENE_ORDER if gene in genes_to_use]
        else:
            gene_names = GENE_ORDER
        tracks_per_gene = 6
        
        # Get gene window metadata for genomic coordinates
        gene_window_metadata = self.config['dataset_input'].get('gene_window_metadata', {})
        
        # Compute interpretability maps if enabled
        gradcam_map = None
        deeplift_map = None
        gradcam_target_class_idx = None
        gradcam_target_class_name = None
        deeplift_target_class_name = None
        
        if self.interpretability_enabled:
            if self.gradcam is not None:
                try:
                    # Determine target class for Grad-CAM
                    if self.gradcam_target_class == 'predicted':
                        gradcam_target_class_idx = predicted_idx
                        gradcam_target_class_name = predicted_name
                    else:
                        # Map class name to index
                        class_name_to_idx = {name: i for i, name in enumerate(class_names)}
                        if self.gradcam_target_class in class_name_to_idx:
                            gradcam_target_class_idx = class_name_to_idx[self.gradcam_target_class]
                            gradcam_target_class_name = self.gradcam_target_class
                        else:
                            console.print(f"[yellow]⚠ Classe '{self.gradcam_target_class}' não encontrada, usando predita[/yellow]")
                            gradcam_target_class_idx = predicted_idx
                            gradcam_target_class_name = predicted_name
                    
                    gradcam_map, _ = self.gradcam.generate(features.clone(), target_class=gradcam_target_class_idx)
                except Exception as e:
                    console.print(f"[yellow]⚠ Erro ao gerar Grad-CAM: {e}[/yellow]")
            
            # Variáveis para modo de média de classe
            deeplift_class_mean_mode = False
            deeplift_class_mean_input = None
            deeplift_class_mean_num_samples = 0
            
            if self.deeplift is not None:
                try:
                    # Determine target class for DeepLIFT
                    if self.deeplift_target_class == 'predicted':
                        deeplift_target_class_idx = predicted_idx
                        deeplift_target_class_name = predicted_name
                        # Calcular individualmente para cada amostra
                        deeplift_map, _ = self.deeplift.generate(
                            features.clone(), 
                            target_class=deeplift_target_class_idx,
                            baseline_type=self.deeplift_baseline,
                            dataset=self.dataset
                        )
                    else:
                        # Map class name to index
                        class_name_to_idx = {name: i for i, name in enumerate(class_names)}
                        if self.deeplift_target_class in class_name_to_idx:
                            deeplift_target_class_idx = class_name_to_idx[self.deeplift_target_class]
                            deeplift_target_class_name = self.deeplift_target_class
                        else:
                            console.print(f"[yellow]⚠ Classe '{self.deeplift_target_class}' não encontrada, usando predita[/yellow]")
                            deeplift_target_class_idx = predicted_idx
                            deeplift_target_class_name = predicted_name
                        
                        # Usar média de todas as amostras da classe alvo
                        deeplift_map, deeplift_class_mean_input, deeplift_class_mean_num_samples = self.deeplift.generate_class_mean(
                            target_class_idx=deeplift_target_class_idx,
                            dataset=self.dataset,
                            baseline_type=self.deeplift_baseline
                        )
                        deeplift_class_mean_mode = True
                except Exception as e:
                    console.print(f"[yellow]⚠ Erro ao gerar DeepLIFT: {e}[/yellow]")
        
        # Determine layout based on interpretability
        has_gradcam = gradcam_map is not None
        has_deeplift = deeplift_map is not None
        num_interp_panels = (1 if has_gradcam else 0) + (1 if has_deeplift else 0)
        
        # Calculate number of rows for subplot
        # Row 1: Input features
        # Row 2: Interpretability maps (if any)
        # Row 3 (or 2): Output probabilities
        if num_interp_panels > 0:
            num_rows_subplot = 3
            fig_height = 14  # Increased for better spacing between plots
        else:
            num_rows_subplot = 2
            fig_height = 8
        
        # Create figure
        plt.clf()
        fig = plt.gcf()
        fig.set_size_inches(16, fig_height)
        
        # Rescale parameters
        viz_height = self.config.get('debug', {}).get('visualization', {}).get('height', 300)
        viz_width = self.config.get('debug', {}).get('visualization', {}).get('width', 600)
        downsample_agg = self.config.get('debug', {}).get('visualization', {}).get('downsample_aggregation', 'max')
        
        # ─────────────────────────────────────────────────────────────────
        # Row 1: Input features (with CAM overlay if available)
        # ─────────────────────────────────────────────────────────────────
        ax1 = plt.subplot(num_rows_subplot, 1, 1)
        
        # Detect if input is 2D or 1D
        if features_cpu.ndim == 3 and features_cpu.shape[0] == 1:
            # 2D input: [1, num_rows, effective_size]
            # Se estamos no modo de média de classe, usar a média das entradas
            if deeplift_class_mean_mode and deeplift_class_mean_input is not None:
                img_data = deeplift_class_mean_input.numpy()  # [num_rows, effective_size]
                input_title = f'{self.dataset_name.upper()} | Class {deeplift_target_class_name} ({deeplift_class_mean_num_samples} samples) | Input 2D Mean ({img_data.shape[0]}x{img_data.shape[1]})'
            else:
                img_data = features_cpu[0].numpy()  # [num_rows, effective_size]
                input_title = f'{self.dataset_name.upper()} | Sample {sample_id} ({target_name}) | Input 2D ({img_data.shape[0]}x{img_data.shape[1]})'
            
            # Calculate zoom factors for later use
            zoom_factors = (viz_height / img_data.shape[0], viz_width / img_data.shape[1])
            
            # Usar block_reduce para preservar picos (ex: taint) ao invés de interpolação
            img_resized = block_reduce_2d(img_data, (viz_height, viz_width), func=downsample_agg)
            
            # Normalize for visualization (0=black, 1=white)
            img_min, img_max = img_resized.min(), img_resized.max()
            if img_max > img_min:
                img_normalized = (img_resized - img_min) / (img_max - img_min)
            else:
                img_normalized = np.zeros_like(img_resized)
            
            # Plot as image (without overlay - interpretability maps shown separately)
            plt.imshow(img_normalized, cmap='gray', aspect='auto', interpolation='nearest')
            
            plt.xlabel('Gene Position (rescaled)', fontsize=12)
            plt.title(input_title, fontsize=14, fontweight='bold')
            plt.colorbar(label='Normalized Value')
            
            # Configure Y axis with gene names
            num_genes = len(gene_names)
            if img_data.shape[0] == num_genes * tracks_per_gene:
                pixels_per_row = viz_height / img_data.shape[0]
                
                # Major ticks at gene boundaries (before and after each gene)
                y_major_ticks = [i * tracks_per_gene * pixels_per_row for i in range(num_genes + 1)]
                ax1.set_yticks(y_major_ticks)
                ax1.set_yticklabels([''] * len(y_major_ticks))  # No labels on boundary ticks
                ax1.tick_params(axis='y', which='major', length=8, width=0.8)
                
                # Minor ticks at center of each gene for labels
                y_minor_ticks = [(i * tracks_per_gene + tracks_per_gene / 2) * pixels_per_row 
                                 for i in range(num_genes)]
                ax1.set_yticks(y_minor_ticks, minor=True)
                ax1.set_yticklabels(gene_names, minor=True, fontsize=10)
                ax1.tick_params(axis='y', which='minor', length=0)  # Hide minor tick marks
                
                ax1.set_ylabel('Genes', fontsize=12)
            else:
                ax1.set_ylabel('Tracks (rescaled)', fontsize=12)
        else:
            # 1D input (fallback for backwards compatibility)
            features_np = features_cpu.numpy().flatten()
            plt.plot(features_np, linewidth=0.5, alpha=0.7)
            plt.xlabel('Feature Index', fontsize=12)
            plt.ylabel('Feature Value', fontsize=12)
            plt.title(f'{self.dataset_name.upper()} | Sample {sample_id} | Input Features (n={len(features_np)})', 
                     fontsize=14, fontweight='bold')
            plt.grid(True, alpha=0.3)
        
        # ─────────────────────────────────────────────────────────────────
        # Row 2: Interpretability maps (if enabled)
        # ─────────────────────────────────────────────────────────────────
        if num_interp_panels > 0:
            # Determine subplot positions
            if has_gradcam and has_deeplift:
                # Two panels side by side
                ax_gc = plt.subplot(num_rows_subplot, 2, 3)
                ax_dl = plt.subplot(num_rows_subplot, 2, 4)
            elif has_gradcam:
                ax_gc = plt.subplot(num_rows_subplot, 1, 2)
            elif has_deeplift:
                ax_dl = plt.subplot(num_rows_subplot, 1, 2)
            
            # Plot Grad-CAM
            if has_gradcam:
                cam_np = gradcam_map.numpy()
                # Usar block_reduce para preservar picos ao invés de interpolação
                cam_resized = block_reduce_2d(cam_np, (viz_height, viz_width), func=downsample_agg)
                
                plt.sca(ax_gc)
                im_gc = plt.imshow(cam_resized, cmap='hot', aspect='auto', interpolation='nearest',
                                   vmin=0.0, vmax=0.0025)
                plt.colorbar(im_gc, label='Activation')
                # Show which class is being visualized and individual's superpopulation
                gc_class_label = gradcam_target_class_name if gradcam_target_class_name else predicted_name
                plt.title(f'Grad-CAM: Important Regions for Class {gc_class_label} (Individual: {target_name})', fontsize=14, fontweight='bold')
                
                # Compute gene importance first (needed for xlabel)
                gene_importance = []
                for i in range(len(gene_names)):
                    start_row = i * tracks_per_gene
                    end_row = (i + 1) * tracks_per_gene
                    if end_row <= cam_np.shape[0]:
                        importance = cam_np[start_row:end_row, :].mean()
                        gene_importance.append((gene_names[i], importance))
                
                # Build xlabel with top genes included (same style as DeepLIFT)
                if gene_importance:
                    gene_importance.sort(key=lambda x: x[1], reverse=True)
                    top_genes_str = ', '.join([f"{g[0]}({g[1]:.5f})" for g in gene_importance[:3]])
                    plt.xlabel(f'Gene Position (rescaled) — Top genes: {top_genes_str}', fontsize=11)
                else:
                    plt.xlabel('Gene Position (rescaled)', fontsize=12)
                
                # Configure Y axis (same style as first plot)
                num_genes = len(gene_names)
                if features_cpu.shape[1] == num_genes * tracks_per_gene:
                    pixels_per_row = viz_height / features_cpu.shape[1]
                    
                    # Major ticks at gene boundaries
                    y_major_ticks = [i * tracks_per_gene * pixels_per_row for i in range(num_genes + 1)]
                    ax_gc.set_yticks(y_major_ticks)
                    ax_gc.set_yticklabels([''] * len(y_major_ticks))
                    ax_gc.tick_params(axis='y', which='major', length=8, width=0.8)
                    
                    # Minor ticks at center of each gene for labels
                    y_minor_ticks = [(i * tracks_per_gene + tracks_per_gene / 2) * pixels_per_row 
                                     for i in range(num_genes)]
                    ax_gc.set_yticks(y_minor_ticks, minor=True)
                    ax_gc.set_yticklabels(gene_names, minor=True, fontsize=10)
                    ax_gc.tick_params(axis='y', which='minor', length=0)
                    
                    ax_gc.set_ylabel('Genes', fontsize=12)
            
            # Plot DeepLIFT
            if has_deeplift:
                dl_np = deeplift_map.numpy()
                
                # Usar block_reduce para preservar picos ao invés de interpolação
                # Para DeepLIFT (valores positivos e negativos), usar max do valor absoluto
                # mas preservando o sinal
                if downsample_agg == 'max':
                    # Para preservar extremos (positivos e negativos), comparar por valor absoluto
                    dl_resized_pos = block_reduce_2d(np.maximum(dl_np, 0), (viz_height, viz_width), func='max')
                    dl_resized_neg = block_reduce_2d(np.minimum(dl_np, 0), (viz_height, viz_width), func='min')
                    # Combinar: usar o que tiver maior magnitude
                    dl_resized = np.where(np.abs(dl_resized_pos) >= np.abs(dl_resized_neg), 
                                          dl_resized_pos, dl_resized_neg)
                else:
                    dl_resized = block_reduce_2d(dl_np, (viz_height, viz_width), func=downsample_agg)
                
                plt.sca(ax_dl)
                
                # Colormap personalizado: preto no centro, clareando para os extremos
                # Azul claro (negativo) <- Preto (zero) -> Vermelho claro (positivo)
                from matplotlib.colors import LinearSegmentedColormap
                colors_diverging = [
                    (0.6, 0.8, 1.0),   # Azul claro (negativo máximo)
                    (0.0, 0.0, 0.0),   # Preto (zero)
                    (1.0, 0.7, 0.6)    # Vermelho claro (positivo máximo)
                ]
                cmap_black_center = LinearSegmentedColormap.from_list('black_center', colors_diverging)
                
                # Escala: dinâmica para modo de média de classe, fixa para outros casos
                if deeplift_class_mean_mode:
                    # Escala dinâmica baseada nos dados
                    dl_abs_max = max(abs(dl_resized.min()), abs(dl_resized.max()))
                    dl_vmin, dl_vmax = -dl_abs_max, dl_abs_max
                else:
                    # Escala fixa para comparação entre amostras
                    dl_vmin, dl_vmax = -0.10, 0.10
                
                im_dl = plt.imshow(dl_resized, cmap=cmap_black_center, aspect='auto', 
                                  interpolation='nearest', vmin=dl_vmin, vmax=dl_vmax)
                plt.colorbar(im_dl, label='Attribution (+ → class, - → not class)')
                
                # Título: diferente para modo de média de classe
                dl_class_label = deeplift_target_class_name if deeplift_target_class_name else predicted_name
                if deeplift_class_mean_mode:
                    plt.title(f'DeepLIFT: Mean Attribution for Class {dl_class_label} ({deeplift_class_mean_num_samples} samples)', fontsize=14, fontweight='bold')
                else:
                    plt.title(f'DeepLIFT: Feature Attribution for Class {dl_class_label} (Individual: {target_name})', fontsize=14, fontweight='bold')
                
                # Compute top 5 most active regions (highest positive values, not absolute)
                top_regions = []
                total_cols = dl_np.shape[1]
                
                for i in range(len(gene_names)):
                    gene_name = gene_names[i]
                    start_row = i * tracks_per_gene
                    end_row = (i + 1) * tracks_per_gene
                    if end_row <= dl_np.shape[0]:
                        gene_data = dl_np[start_row:end_row, :]
                        # Find maximum positive value and its position
                        max_val = gene_data.max()
                        if max_val > 0:
                            # Find position of maximum
                            max_idx = np.unravel_index(gene_data.argmax(), gene_data.shape)
                            col_idx = max_idx[1]
                            
                            # Calculate genomic coordinate
                            if gene_name in gene_window_metadata:
                                meta = gene_window_metadata[gene_name]
                                chrom = meta.get('chromosome', 'N/A')
                                start_pos = meta.get('start', 0)
                                window_size = meta.get('window_size', 0)
                                # Calculate position within the window
                                genomic_pos = start_pos + int((col_idx / total_cols) * window_size) if window_size > 0 else start_pos
                            else:
                                chrom = 'N/A'
                                genomic_pos = 0
                            
                            top_regions.append((gene_name, max_val, chrom, genomic_pos, col_idx))
                
                # Sort by value (highest positive first) and take top 5
                top_regions.sort(key=lambda x: x[1], reverse=True)
                top_5_regions = top_regions[:5]
                
                # Print details to terminal
                if top_5_regions:
                    console.print(f"\n[bold cyan]Top 5 Regiões Mais Ativas (DeepLIFT):[/bold cyan]")
                    for rank, (gene, val, chrom, pos, col) in enumerate(top_5_regions, 1):
                        console.print(f"  {rank}. [green]{gene}[/green]: valor = {val:.6f}, {chrom}: {pos:,}")
                
                # Build xlabel with top 5 regions
                if top_5_regions:
                    top_regions_str = ', '.join([f"{r[0]}({r[1]:.5f})" for r in top_5_regions])
                    plt.xlabel(f'Gene Position (rescaled) — Top regions: {top_regions_str}', fontsize=11)
                else:
                    plt.xlabel('Gene Position (rescaled)', fontsize=12)
                
                # Configure Y axis (same style as first plot)
                num_genes = len(gene_names)
                if features_cpu.shape[1] == num_genes * tracks_per_gene:
                    pixels_per_row = viz_height / features_cpu.shape[1]
                    
                    # Major ticks at gene boundaries
                    y_major_ticks = [i * tracks_per_gene * pixels_per_row for i in range(num_genes + 1)]
                    ax_dl.set_yticks(y_major_ticks)
                    ax_dl.set_yticklabels([''] * len(y_major_ticks))
                    ax_dl.tick_params(axis='y', which='major', length=8, width=0.8)
                    
                    # Minor ticks at center of each gene for labels
                    y_minor_ticks = [(i * tracks_per_gene + tracks_per_gene / 2) * pixels_per_row 
                                     for i in range(num_genes)]
                    ax_dl.set_yticks(y_minor_ticks, minor=True)
                    ax_dl.set_yticklabels(gene_names, minor=True, fontsize=10)
                    ax_dl.tick_params(axis='y', which='minor', length=0)
                    
                    ax_dl.set_ylabel('Genes', fontsize=12)
        
        # ─────────────────────────────────────────────────────────────────
        # Last Row: Output probabilities
        # ─────────────────────────────────────────────────────────────────
        ax_prob = plt.subplot(num_rows_subplot, 1, num_rows_subplot)
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
        
        # Text with prediction and target
        correct = "✓ CORRECT" if predicted_idx == target_idx else "✗ WRONG"
        color = 'green' if predicted_idx == target_idx else 'red'
        result_text = (f'{correct}\n'
                      f'Target: {target_name} (class {target_idx})\n'
                      f'Predicted: {predicted_name} (class {predicted_idx}, prob={output_probs[predicted_idx]:.3f})')
        
        plt.text(0.98, 0.98, result_text, transform=plt.gca().transAxes,
                fontsize=11, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
        
        # Adjust layout with more vertical spacing between subplots
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.35)  # Increase vertical spacing between plots
        
        # ─────────────────────────────────────────────────────────────────
        # Save image if enabled
        # ─────────────────────────────────────────────────────────────────
        if self.interpretability_enabled and self.interp_save_images:
            # Create output directory
            output_dir = Path(self.interp_output_dir)
            if not output_dir.is_absolute():
                # Make relative to config's processed_cache_dir or current directory
                cache_dir = self.config.get('dataset_input', {}).get('processed_cache_dir', '.')
                output_dir = Path(cache_dir) / output_dir
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate filename
            method_suffix = self.interp_method
            
            # Nome diferente para modo de média de classe
            if deeplift_class_mean_mode:
                filename = f"class_mean_{deeplift_target_class_name}_{deeplift_class_mean_num_samples}samples_{method_suffix}.png"
            else:
                correct_str = "correct" if predicted_idx == target_idx else "wrong"
                filename = f"{sample_id}_{predicted_name}_{correct_str}_{method_suffix}.png"
            
            filepath = output_dir / filename
            
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            console.print(f"[dim]Saved: {filepath}[/dim]")
        
        # Connect key callback
        self._key_pressed = False
        cid = fig.canvas.mpl_connect('key_press_event', self._on_key_press)
        
        # Show and wait for key
        plt.show(block=True)
        
        # Disconnect callback
        fig.canvas.mpl_disconnect(cid)
    
    def test(self) -> Dict:
        """
        Executa teste e gera métricas.
        
        Returns:
            Dict com resultados
        """
        console.print(Panel.fit(
            f"[bold cyan]Executando Teste[/bold cyan]\n"
            f"Conjunto: {self.dataset_name}"
        ))
        
        # Imprimir arquitetura da rede
        console.print("\n[bold cyan]═══ ARQUITETURA DA REDE ═══[/bold cyan]")
        console.print(self.model)
        console.print(f"\n[bold green]Total de parâmetros treináveis: {self.model.count_parameters():,}[/bold green]")
        console.print()
        
        if self.enable_visualization:
            console.print("[yellow]📊 Modo de visualização habilitado - mostrando gráficos interativos[/yellow]")
            if self.max_samples:
                console.print(f"[yellow]   Limitado a {self.max_samples} amostras[/yellow]")
            if self.interpretability_enabled:
                console.print(f"[cyan]🔍 Interpretabilidade habilitada - método: {self.interp_method}[/cyan]")
                if self.interp_save_images:
                    console.print(f"[cyan]   Salvando imagens em: {self.interp_output_dir}[/cyan]")
        
        self.model.eval()
        all_predictions = []
        all_targets = []
        all_probs = []
        
        # Se interpretabilidade está habilitada, precisamos de gradientes
        # então não usamos torch.no_grad() no loop de visualização
        if self.enable_visualization and self.interpretability_enabled:
            # Loop COM gradientes para interpretabilidade
            sample_count = 0
            for batch_idx, (features, targets, indices) in enumerate(self.test_loader):
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                
                # Forward pass normal (sem gradientes) para obter outputs
                with torch.no_grad():
                    outputs = self.model(features)
                
                # Visualizar amostra (a interpretabilidade faz seu próprio forward com gradientes)
                sample_idx = indices[0].item()  # batch_size=1 na visualização
                sample_id = self.dataset.get_sample_id(sample_idx)
                
                # Imprimir informações no terminal
                self._print_sample_info(sample_idx, sample_id, targets, outputs)
                
                self._visualize_sample(features.clone(), targets, outputs, sample_idx, sample_id)
                
                # Check if user requested to quit
                if self._quit_requested:
                    console.print("[yellow]⚠ Quit requested by user (pressed 'q')[/yellow]")
                    import sys
                    sys.exit(0)
                
                if self.config['output']['prediction_target'] != 'frog_likelihood':
                    predictions = torch.argmax(outputs, dim=1)
                    all_predictions.extend(predictions.cpu().numpy())
                    all_targets.extend(targets.cpu().numpy())
                    all_probs.append(outputs.cpu().numpy())
                
                sample_count += 1
                
                # Limitar número de amostras se especificado
                if self.max_samples is not None and sample_count >= self.max_samples:
                    console.print(f"[yellow]⚠ Limitado a {self.max_samples} amostras (debug)[/yellow]")
                    break
        
        elif self.enable_visualization:
            # Loop SEM interpretabilidade - visualização simples com torch.no_grad()
            with torch.no_grad():
                sample_count = 0
                for batch_idx, (features, targets, indices) in enumerate(self.test_loader):
                    features = features.to(self.device, non_blocking=True)
                    targets = targets.to(self.device, non_blocking=True)
                    
                    outputs = self.model(features)
                    
                    # Visualizar amostra
                    sample_idx = indices[0].item()  # batch_size=1 na visualização
                    sample_id = self.dataset.get_sample_id(sample_idx)
                    
                    # Imprimir informações no terminal
                    self._print_sample_info(sample_idx, sample_id, targets, outputs)
                    
                    self._visualize_sample(features, targets, outputs, sample_idx, sample_id)
                    
                    # Check if user requested to quit
                    if self._quit_requested:
                        console.print("[yellow]⚠ Quit requested by user (pressed 'q')[/yellow]")
                        import sys
                        sys.exit(0)
                    
                    if self.config['output']['prediction_target'] != 'frog_likelihood':
                        predictions = torch.argmax(outputs, dim=1)
                        all_predictions.extend(predictions.cpu().numpy())
                        all_targets.extend(targets.cpu().numpy())
                        all_probs.append(outputs.cpu().numpy())
                    
                    sample_count += 1
                    
                    # Limitar número de amostras se especificado
                    if self.max_samples is not None and sample_count >= self.max_samples:
                        console.print(f"[yellow]⚠ Limitado a {self.max_samples} amostras (debug)[/yellow]")
                        break
        
        else:
            # Modo normal sem visualização - com progress bar e torch.no_grad()
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
                    
                    for features, targets, indices in self.test_loader:
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
        console.print(f"\n[bold cyan]═══ RESULTADOS DO TESTE ({self.dataset_name.upper()}) ═══[/bold cyan]\n")
        
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
        # Listar nomes das classes
        console.print(f"[dim]Classes: {', '.join([f'{i}={name}' for i, name in enumerate(target_names)])}[/dim]")
        console.print(results['classification_report'])
        
        # Confusion Matrix
        console.print("\n[bold]Confusion Matrix:[/bold]")
        # Listar nomes das classes
        console.print(f"[dim]Classes: {', '.join([f'{i}={name}' for i, name in enumerate(target_names)])}[/dim]")
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
    
    # Se taint_at_cache_save=True e mode=train, forçar recriação do cache
    # Isso permite que o usuário controle os parâmetros de tainting sem verificar versão do cache
    taint_at_cache_save = config.get('debug', {}).get('taint_at_cache_save', False)
    mode = config.get('mode', 'train')
    if taint_at_cache_save and mode == 'train':
        console.print(f"[yellow]Cache será recriado: taint_at_cache_save=True em modo train[/yellow]")
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
        'splits_metadata.json',
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
        
        # Verificar versão do formato - requer formato v2 com metadados completos em splits_metadata.json
        format_version = metadata.get('format_version', 1)
        if format_version < 2:
            console.print(f"[yellow]Cache invalidado: formato legado (v{format_version}). Requer v2 com metadados completos.[/yellow]")
            console.print(f"[yellow]  O cache será recriado automaticamente.[/yellow]")
            return False
        
        # Verificar compatibilidade dos parâmetros críticos
        processing_params = metadata.get('processing_params', {})
        
        # Verificar formato de entrada (1D vs 2D)
        cached_input_shape = processing_params.get('input_shape', '1D')
        if cached_input_shape != '2D':
            console.print(f"[yellow]Cache invalidado: formato de entrada mudou de 1D para 2D[/yellow]")
            return False
        # Parâmetros básicos que sempre invalidam o cache
        # NOTA: taint_type, taint_value e outros parâmetros de tainting NÃO são verificados aqui
        # Se taint_at_cache_save=True em mode=train, o cache será recriado de qualquer forma
        current_params = {
            'alphagenome_outputs': config['dataset_input']['alphagenome_outputs'],
            'haplotype_mode': config['dataset_input']['haplotype_mode'],
            'window_center_size': config['dataset_input']['window_center_size'],
            'downsample_factor': config['dataset_input']['downsample_factor'],
            'normalization_method': config['dataset_input'].get('normalization_method', 'zscore'),
            'balancing_strategy': config['data_split'].get('balancing_strategy', 'stratified'),
            'dataset_dir': config['dataset_input']['dataset_dir'],
            'input_shape': '2D',
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


def get_gene_window_metadata(
    dataset_dir: Path,
    gene_order: List[str],
    window_center_size: int
) -> Optional[Dict]:
    """
    Obtém metadados das janelas genômicas com coordenadas ajustadas para window_center_size.
    
    Args:
        dataset_dir: Diretório do dataset base (contém individuals/)
        gene_order: Lista ordenada de genes
        window_center_size: Tamanho da janela central extraída
        
    Returns:
        Dicionário com metadados de cada gene, ou None se não for possível obter
    """
    if not gene_order:
        return None
    
    dataset_dir = Path(dataset_dir)
    individuals_dir = dataset_dir / 'individuals'
    
    if not individuals_dir.exists():
        console.print(f"[yellow]⚠ Diretório individuals não encontrado: {individuals_dir}[/yellow]")
        return None
    
    # Encontrar primeiro indivíduo válido
    individual_dirs = [d for d in individuals_dir.iterdir() if d.is_dir()]
    if not individual_dirs:
        console.print(f"[yellow]⚠ Nenhum indivíduo encontrado em {individuals_dir}[/yellow]")
        return None
    
    # Usar primeiro indivíduo para obter window_metadata
    first_individual = individual_dirs[0]
    metadata_file = first_individual / 'individual_metadata.json'
    
    if not metadata_file.exists():
        console.print(f"[yellow]⚠ individual_metadata.json não encontrado: {metadata_file}[/yellow]")
        return None
    
    try:
        with open(metadata_file, 'r') as f:
            individual_metadata = json.load(f)
        
        window_metadata = individual_metadata.get('window_metadata', {})
        if not window_metadata:
            console.print(f"[yellow]⚠ window_metadata não encontrado em {metadata_file}[/yellow]")
            return None
        
        gene_window_metadata = {}
        
        for gene_name in gene_order:
            if gene_name not in window_metadata:
                console.print(f"[yellow]⚠ Gene {gene_name} não encontrado em window_metadata[/yellow]")
                continue
            
            gene_meta = window_metadata[gene_name]
            original_start = gene_meta.get('start', 0)
            original_window_size = gene_meta.get('window_size', 0)
            chromosome = gene_meta.get('chromosome', '')
            
            # Calcular centro da janela original
            # original_window_size é o tamanho -1 (end - start), então temos original_window_size + 1 elementos
            center_idx = (original_window_size + 1) // 2
            genomic_center = original_start + center_idx
            
            # Calcular nova janela centrada
            half_size = window_center_size // 2
            new_start = genomic_center - half_size
            new_end = genomic_center + half_size - 1
            new_window_size = new_end - new_start  # = window_center_size - 1
            
            gene_window_metadata[gene_name] = {
                'chromosome': chromosome,
                'start': new_start,
                'end': new_end,
                'window_size': new_window_size
            }
        
        return gene_window_metadata
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro ao obter gene_window_metadata: {e}[/yellow]")
        return None


def save_processed_dataset(
    cache_dir: Path,
    processed_dataset: ProcessedGenomicDataset,
    train_indices: List[int],
    val_indices: List[int],
    test_indices: List[int],
    config: Dict
):
    """
    Salva dataset processado em cache com balanceamento estratificado sequencial.
    
    Os dados são organizados em grupos onde cada grupo contém exatamente uma amostra
    de cada classe (superpopulation), garantindo balanceamento. Amostras excedentes
    de classes majoritárias são descartadas.
    
    A ordem é: [classe1_sample1, classe2_sample1, classe3_sample1, 
                classe1_sample2, classe2_sample2, classe3_sample2, ...]
    
    Args:
        cache_dir: Diretório onde salvar cache
        processed_dataset: Dataset processado
        train_indices: Índices de treino (serão reorganizados)
        val_indices: Índices de validação (serão reorganizados)
        test_indices: Índices de teste (serão reorganizados)
        config: Configuração usada
    """
    import random
    
    cache_dir = Path(cache_dir)
    random_seed = config['data_split']['random_seed']
    
    # Criar diretório temporário para escrita atômica
    temp_cache_dir = cache_dir.parent / f"{cache_dir.name}_tmp_{os.getpid()}"
    
    # Limpar temp dir se existir (de alguma execução anterior interrompida)
    if temp_cache_dir.exists():
        console.print(f"[yellow]Limpando diretório temporário de execução anterior...[/yellow]")
        shutil.rmtree(temp_cache_dir)
    
    temp_cache_dir.mkdir(parents=True, exist_ok=True)
    
    console.print(f"\n[bold cyan]💾 Salvando Dataset Processado em Cache (Estratificado)[/bold cyan]")
    console.print(f"  📁 Diretório: {cache_dir}")
    console.print(f"  🎲 Random seed: {random_seed}")
    
    # Carregar metadados do dataset fonte para obter informações dos indivíduos
    dataset_dir = Path(config['dataset_input']['dataset_dir'])
    dataset_metadata_file = dataset_dir / 'dataset_metadata.json'
    if not dataset_metadata_file.exists():
        raise FileNotFoundError(f"dataset_metadata.json não encontrado em {dataset_dir}")
    
    with open(dataset_metadata_file, 'r') as f:
        source_metadata = json.load(f)
    
    individuals = source_metadata.get('individuals', [])
    individuals_pedigree = source_metadata.get('individuals_pedigree', {})
    
    if not individuals:
        raise ValueError("Lista de individuals vazia em dataset_metadata.json")
    
    console.print(f"  📋 Metadados carregados: {len(individuals)} indivíduos")
    
    def get_individual_metadata(idx: int) -> Dict:
        """Obtém metadados completos de um indivíduo pelo índice."""
        if idx >= len(individuals):
            return {"sample_id": f"sample_{idx}", "superpopulation": "UNK", "population": "UNK", "sex": 0}
        
        sample_id = individuals[idx]
        pedigree = individuals_pedigree.get(sample_id, {})
        return {
            "sample_id": sample_id,
            "superpopulation": pedigree.get("superpopulation", "UNK"),
            "population": pedigree.get("population", "UNK"),
            "sex": pedigree.get("sex", 0)
        }
    
    def stratified_interleave(indices: List[int], seed: int) -> Tuple[List[int], List[int], Dict[str, int]]:
        """
        Organiza índices em ordem estratificada sequencial.
        
        Retorna:
            - Lista de índices intercalados por classe
            - Lista de índices descartados (excesso de classes majoritárias)
            - Dicionário com contagem por classe
        """
        # Agrupar índices por classe (superpopulation)
        indices_by_class = {}
        for idx in indices:
            meta = get_individual_metadata(idx)
            superpop = meta['superpopulation']
            if superpop not in indices_by_class:
                indices_by_class[superpop] = []
            indices_by_class[superpop].append(idx)
        
        # Embaralhar cada classe com seed fixa para reprodutibilidade
        rng = random.Random(seed)
        for class_name in indices_by_class:
            rng.shuffle(indices_by_class[class_name])
        
        # Determinar tamanho mínimo entre todas as classes
        class_counts = {c: len(idxs) for c, idxs in indices_by_class.items()}
        min_count = min(class_counts.values())
        
        # Coletar índices descartados (excesso)
        discarded = []
        for class_name, idxs in indices_by_class.items():
            if len(idxs) > min_count:
                discarded.extend(idxs[min_count:])
        
        # Truncar cada classe ao tamanho mínimo
        for class_name in indices_by_class:
            indices_by_class[class_name] = indices_by_class[class_name][:min_count]
        
        # Intercalar: classe1[0], classe2[0], ..., classeN[0], classe1[1], ...
        # Ordenar classes alfabeticamente para determinismo
        sorted_classes = sorted(indices_by_class.keys())
        interleaved = []
        for i in range(min_count):
            for class_name in sorted_classes:
                interleaved.append(indices_by_class[class_name][i])
        
        return interleaved, discarded, class_counts
    
    def simple_shuffle(indices: List[int], seed: int) -> Tuple[List[int], Dict[str, int]]:
        """
        Embaralha índices de forma simples sem descartar dados.
        
        Retorna:
            - Lista de índices embaralhados
            - Dicionário com contagem por classe (para info)
        """
        # Contar classes para info
        class_counts = {}
        for idx in indices:
            meta = get_individual_metadata(idx)
            superpop = meta['superpopulation']
            class_counts[superpop] = class_counts.get(superpop, 0) + 1
        
        # Embaralhar com seed fixa para reprodutibilidade
        rng = random.Random(seed)
        shuffled = list(indices)
        rng.shuffle(shuffled)
        
        return shuffled, class_counts
    
    # Combinar todos os índices
    all_indices = list(train_indices) + list(val_indices) + list(test_indices)
    console.print(f"  📊 Total de amostras: {len(all_indices)}")
    
    # Ler estratégia de balanceamento do config
    balancing_strategy = config['data_split'].get('balancing_strategy', 'stratified')
    
    train_split = config['data_split']['train_split']
    val_split = config['data_split']['val_split']
    test_split = config['data_split']['test_split']
    
    if balancing_strategy == 'stratified':
        # Aplicar balanceamento estratificado (descarta excedentes)
        balanced_indices, discarded_indices, original_class_counts = stratified_interleave(all_indices, random_seed)
        
        num_classes = len(original_class_counts)
        samples_per_class = len(balanced_indices) // num_classes
        
        console.print(f"\n  [bold]Balanceamento Estratificado:[/bold]")
        console.print(f"  • Classes encontradas: {num_classes}")
        for class_name, count in sorted(original_class_counts.items()):
            discarded_count = count - samples_per_class
            status = f"[yellow](-{discarded_count})[/yellow]" if discarded_count > 0 else "[green](OK)[/green]"
            console.print(f"    - {class_name}: {count} → {samples_per_class} {status}")
        console.print(f"  • Amostras balanceadas: {len(balanced_indices)}")
        console.print(f"  • Amostras descartadas: {len(discarded_indices)}")
        
        # Calcular splits sobre os dados balanceados
        # Os grupos devem permanecer intactos (não quebrar grupos de N classes)
        total_groups = len(balanced_indices) // num_classes
        
        train_groups = int(total_groups * train_split)
        val_groups = int(total_groups * val_split)
        test_groups = total_groups - train_groups - val_groups  # Resto vai para teste
        
        # Converter grupos para índices
        train_end = train_groups * num_classes
        val_end = train_end + val_groups * num_classes
        
        new_train_indices = balanced_indices[:train_end]
        new_val_indices = balanced_indices[train_end:val_end]
        new_test_indices = balanced_indices[val_end:]
        
        console.print(f"\n  [bold]Splits após balanceamento:[/bold]")
        console.print(f"  • Train: {len(new_train_indices)} amostras ({train_groups} grupos)")
        console.print(f"  • Val: {len(new_val_indices)} amostras ({val_groups} grupos)")
        console.print(f"  • Test: {len(new_test_indices)} amostras ({test_groups} grupos)")
        
    else:  # shuffle
        # Randomização simples (preserva todos os dados)
        shuffled_indices, class_counts = simple_shuffle(all_indices, random_seed)
        discarded_indices = []  # Nenhum dado descartado
        
        console.print(f"\n  [bold]Randomização Simples (shuffle):[/bold]")
        console.print(f"  • Classes encontradas: {len(class_counts)}")
        for class_name, count in sorted(class_counts.items()):
            console.print(f"    - {class_name}: {count}")
        console.print(f"  • Total de amostras: {len(shuffled_indices)} [green](nenhuma descartada)[/green]")
        
        # Calcular splits simples
        total_samples = len(shuffled_indices)
        train_end = int(total_samples * train_split)
        val_end = train_end + int(total_samples * val_split)
        
        new_train_indices = shuffled_indices[:train_end]
        new_val_indices = shuffled_indices[train_end:val_end]
        new_test_indices = shuffled_indices[val_end:]
        
        console.print(f"\n  [bold]Splits:[/bold]")
        console.print(f"  • Train: {len(new_train_indices)} amostras")
        console.print(f"  • Val: {len(new_val_indices)} amostras")
        console.print(f"  • Test: {len(new_test_indices)} amostras")
    
    # Preparar dados de cada split
    train_data = []
    val_data = []
    test_data = []
    
    # Preparar metadados de cada split (na mesma ordem dos dados)
    train_metadata = []
    val_metadata = []
    test_metadata = []
    
    # Verificar se tainting está habilitado ao salvar cache
    taint_at_save = config.get('debug', {}).get('taint_at_cache_save', False)
    taint_num_classes = processed_dataset.get_num_classes() if taint_at_save else 0
    debug_config = config.get('debug', {})
    taint_type = debug_config.get('taint_type', 'additive')
    taint_value = debug_config.get('taint_value', 1.0)
    taint_horizontal_size = debug_config.get('taint_horizontal_size', 6)
    taint_vertical_size = debug_config.get('taint_vertical_size', 6)
    taint_horizontal_step = debug_config.get('taint_horizontal_step', 100)
    taint_vertical_step = debug_config.get('taint_vertical_step', 6)
    
    # IMPORTANTE: Desabilitar temporariamente taint_at_runtime enquanto salvamos o cache
    original_taint_runtime = processed_dataset.config.get('debug', {}).get('taint_at_runtime', False)
    if 'debug' not in processed_dataset.config:
        processed_dataset.config['debug'] = {}
    processed_dataset.config['debug']['taint_at_runtime'] = False
    
    try:
        # Processar dados de treino
        console.print(f"\n  📦 Processando dados de treino...")
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Train", total=len(new_train_indices))
            for idx in new_train_indices:
                features, target = processed_dataset[idx]
                
                # Aplicar tainting se habilitado
                if taint_at_save and config['output']['prediction_target'] != 'frog_likelihood':
                    if target.ndim == 0:
                        target_class = target.item()
                    else:
                        target_class = target[0].item()
                    features = taint_sample(
                        features, target_class, taint_num_classes,
                        taint_type=taint_type, taint_value=taint_value,
                        taint_horizontal_size=taint_horizontal_size,
                        taint_vertical_size=taint_vertical_size,
                        taint_horizontal_step=taint_horizontal_step,
                        taint_vertical_step=taint_vertical_step
                    )
                
                train_data.append((features, target))
                train_metadata.append(get_individual_metadata(idx))
                progress.update(task, advance=1)
        
        # Processar dados de validação
        console.print(f"  📦 Processando dados de validação...")
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Val", total=len(new_val_indices))
            for idx in new_val_indices:
                features, target = processed_dataset[idx]
                
                if taint_at_save and config['output']['prediction_target'] != 'frog_likelihood':
                    if target.ndim == 0:
                        target_class = target.item()
                    else:
                        target_class = target[0].item()
                    features = taint_sample(
                        features, target_class, taint_num_classes,
                        taint_type=taint_type, taint_value=taint_value,
                        taint_horizontal_size=taint_horizontal_size,
                        taint_vertical_size=taint_vertical_size,
                        taint_horizontal_step=taint_horizontal_step,
                        taint_vertical_step=taint_vertical_step
                    )
                
                val_data.append((features, target))
                val_metadata.append(get_individual_metadata(idx))
                progress.update(task, advance=1)
        
        # Processar dados de teste
        console.print(f"  📦 Processando dados de teste...")
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Test", total=len(new_test_indices))
            for idx in new_test_indices:
                features, target = processed_dataset[idx]
                
                if taint_at_save and config['output']['prediction_target'] != 'frog_likelihood':
                    if target.ndim == 0:
                        target_class = target.item()
                    else:
                        target_class = target[0].item()
                    features = taint_sample(
                        features, target_class, taint_num_classes,
                        taint_type=taint_type, taint_value=taint_value,
                        taint_horizontal_size=taint_horizontal_size,
                        taint_vertical_size=taint_vertical_size,
                        taint_horizontal_step=taint_horizontal_step,
                        taint_vertical_step=taint_vertical_step
                    )
                
                test_data.append((features, target))
                test_metadata.append(get_individual_metadata(idx))
                progress.update(task, advance=1)
        
        # Salvar dados no diretório temporário
        console.print(f"\n  💾 Salvando train_data.pt ({len(train_data)} amostras)...")
        torch.save(train_data, temp_cache_dir / 'train_data.pt')
        
        console.print(f"  💾 Salvando val_data.pt ({len(val_data)} amostras)...")
        torch.save(val_data, temp_cache_dir / 'val_data.pt')
        
        console.print(f"  💾 Salvando test_data.pt ({len(test_data)} amostras)...")
        torch.save(test_data, temp_cache_dir / 'test_data.pt')
        
        # Preparar metadados dos samples descartados
        discarded_metadata = [get_individual_metadata(idx) for idx in discarded_indices]
        
        # Salvar splits com metadados completos
        console.print(f"  💾 Salvando splits_metadata.json (com metadados completos)...")
        splits = {
            'format_version': 3,  # Nova versão com balanceamento estratificado
            'train': train_metadata,
            'val': val_metadata,
            'test': test_metadata,
            'discarded': discarded_metadata  # Samples descartados por balanceamento
        }
        with open(temp_cache_dir / 'splits_metadata.json', 'w') as f:
            json.dump(splits, f, indent=2)
        
        # Salvar normalization params
        console.print(f"  💾 Salvando normalization_params.json...")
        with open(temp_cache_dir / 'normalization_params.json', 'w') as f:
            json.dump(processed_dataset.normalization_params, f, indent=2)
        
        # Obter ordem dos genes do dataset base
        try:
            first_sample_data, _ = processed_dataset.base_dataset[0]
            gene_order = list(first_sample_data['windows'])
            console.print(f"  📊 Ordem dos genes detectada: {', '.join(gene_order)}")
        except Exception as e:
            console.print(f"[yellow]⚠ Não foi possível obter ordem dos genes: {e}[/yellow]")
            gene_order = None
        
        # Obter metadados das janelas genômicas com coordenadas ajustadas
        window_center_size = config['dataset_input']['window_center_size']
        dataset_dir = Path(config['dataset_input']['dataset_dir'])
        gene_window_metadata = get_gene_window_metadata(dataset_dir, gene_order, window_center_size)
        if gene_window_metadata:
            console.print(f"  📊 Metadados de janelas genômicas obtidos para {len(gene_window_metadata)} genes")
        
        # Construir balancing_info conforme estratégia
        if balancing_strategy == 'stratified':
            balanced_class_counts = {c: samples_per_class for c in original_class_counts.keys()}
            balancing_info = {
                'strategy': 'stratified',
                'method': 'stratified_sequential',
                'original_total': len(all_indices),
                'balanced_total': len(new_train_indices) + len(new_val_indices) + len(new_test_indices),
                'discarded_total': len(discarded_indices),
                'samples_per_class': samples_per_class,
                'original_class_counts': original_class_counts,
                'balanced_class_counts': balanced_class_counts,
            }
            splits_info = {
                'train_size': len(new_train_indices),
                'val_size': len(new_val_indices),
                'test_size': len(new_test_indices),
                'train_groups': train_groups,
                'val_groups': val_groups,
                'test_groups': test_groups,
                'train_split': config['data_split']['train_split'],
                'val_split': config['data_split']['val_split'],
                'test_split': config['data_split']['test_split'],
                'random_seed': random_seed
            }
            total_samples = len(new_train_indices) + len(new_val_indices) + len(new_test_indices)
        else:  # shuffle
            balancing_info = {
                'strategy': 'shuffle',
                'method': 'simple_random',
                'original_total': len(all_indices),
                'balanced_total': len(all_indices),  # Todos preservados
                'discarded_total': 0,
                'class_counts': class_counts,
            }
            splits_info = {
                'train_size': len(new_train_indices),
                'val_size': len(new_val_indices),
                'test_size': len(new_test_indices),
                'train_split': config['data_split']['train_split'],
                'val_split': config['data_split']['val_split'],
                'test_split': config['data_split']['test_split'],
                'random_seed': random_seed
            }
            total_samples = len(all_indices)
        
        # Salvar metadados
        metadata = {
            'creation_date': datetime.now().isoformat(),
            'format_version': 3,
            'dataset_dir': config['dataset_input']['dataset_dir'],
            'processing_params': {
                'alphagenome_outputs': config['dataset_input']['alphagenome_outputs'],
                'haplotype_mode': config['dataset_input']['haplotype_mode'],
                'window_center_size': config['dataset_input']['window_center_size'],
                'downsample_factor': config['dataset_input']['downsample_factor'],
                'normalization_method': config['dataset_input'].get('normalization_method', 'zscore'),
                'balancing_strategy': balancing_strategy,
                'input_shape': '2D',
            },
            'balancing_info': balancing_info,
            'splits': splits_info,
            'total_samples': total_samples,
            'num_classes': processed_dataset.get_num_classes(),
            'input_size': processed_dataset.get_input_size(),
            'prediction_target': config['output']['prediction_target'],
            'class_names': processed_dataset.idx_to_target,
            'gene_order': gene_order,
            'gene_window_metadata': gene_window_metadata,
            'tracks_per_gene': 6
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
    
    Suporta dois modos de carregamento:
    - preload: Carrega tudo na RAM (rápido, usa muita memória)
    - lazy: Carrega sob demanda (economiza memória, pequeno overhead I/O)
    """
    
    def __init__(
        self,
        data_file: Path,
        target_to_idx: Dict,
        idx_to_target: Dict,
        config: Dict,
        split_name: Optional[str] = None,
        length_hint: Optional[int] = None
    ):
        """
        Inicializa dataset do cache.
        
        Args:
            data_file: Arquivo .pt com dados processados
            target_to_idx: Mapeamento target->índice
            idx_to_target: Mapeamento índice->target
            config: Configuração (necessário para taint_at_runtime e loading_strategy)
            split_name: Nome do split (train/val/test) para resolver tamanho via metadata
            length_hint: Comprimento esperado (evita carregar o arquivo apenas para len)
        """
        self.data_file = data_file
        self.target_to_idx = target_to_idx
        self.idx_to_target = idx_to_target
        self.config = config
        self.split_name = (split_name or self._infer_split_name()).lower()
        self._length_hint = length_hint
        
        # Carregar metadados dos samples (sample_id, superpopulation, population, sex)
        self._sample_metadata = self._load_sample_metadata()
        
        # Ler configuração de carregamento
        data_loading_config = config.get('data_loading', {})
        self.loading_strategy = data_loading_config.get('loading_strategy', 'preload')
        self.cache_size = data_loading_config.get('cache_size', 100)
        
        console.print(f"[cyan]Preparando {data_file.name}...[/cyan]")
        
        if self.loading_strategy == 'preload':
            # Modo tradicional: carrega tudo na RAM
            self.data = torch.load(data_file)
            self._length = len(self.data)
            self._data_loaded = True
            console.print(f"[green]✓ {self._length} samples carregados (preload)[/green]")
        
        elif self.loading_strategy == 'lazy':
            # Modo lazy: carrega apenas metadados agora
            self._length = self._determine_length_without_loading()
            
            self._data_loaded = False
            self.data = None
            
            # Cache LRU para samples recentemente acessados
            self._cache = {}  # {idx: (features, target)}
            self._cache_order = []  # Lista para LRU
            
            console.print(f"[green]✓ {self._length} samples preparados (lazy loading)[/green]")
            console.print(f"  • Modo: Lazy loading (carrega sob demanda)")
            console.print(f"  • Cache LRU: {self.cache_size} samples")
            console.print(f"  • Use unload_data() para liberar memória entre épocas")
        
        else:
            raise ValueError(f"loading_strategy inválida: {self.loading_strategy}. Use 'preload' ou 'lazy'.")
    
    def __len__(self) -> int:
        return self._length
    
    def _load_sample_metadata(self) -> Optional[List[Dict]]:
        """
        Load sample metadata from splits_metadata.json.
        
        Supports two formats:
        - Format v2 (new): splits_metadata.json contains complete metadata per sample
        - Format v1 (legacy): splits_metadata.json contains only indices, needs dataset_metadata.json
        
        Returns:
            List of dicts with sample metadata, or None if loading fails
        """
        try:
            cache_dir = self.data_file.parent
            
            # Load splits_metadata.json
            splits_file = cache_dir / 'splits_metadata.json'
            if not splits_file.exists():
                console.print(f"[yellow]⚠ Sample metadata: splits_metadata.json não encontrado em {cache_dir}[/yellow]")
                return None
            with open(splits_file, 'r') as f:
                splits = json.load(f)
            
            # Check format version
            format_version = splits.get('format_version', 1)
            
            if format_version >= 2:
                # New format: splits_metadata.json contains complete metadata
                split_data = splits.get(self.split_name, [])
                if not split_data:
                    console.print(f"[yellow]⚠ Sample metadata: '{self.split_name}' não encontrado em splits_metadata.json[/yellow]")
                    return None
                return split_data
            else:
                # Legacy format: need to load from dataset_metadata.json
                indices_key = f"{self.split_name}_indices"
                split_indices = splits.get(indices_key, [])
                if not split_indices:
                    console.print(f"[yellow]⚠ Sample metadata: {indices_key} não encontrado em splits_metadata.json (formato legado)[/yellow]")
                    return None
                
                # Load dataset metadata to get individual IDs
                metadata_file = cache_dir / 'metadata.json'
                if not metadata_file.exists():
                    return None
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                dataset_dir_str = metadata.get('dataset_dir', '') or metadata.get('processing_params', {}).get('dataset_dir', '')
                if not dataset_dir_str:
                    dataset_dir_str = self.config.get('dataset_input', {}).get('dataset_dir', '')
                
                if not dataset_dir_str:
                    return None
                
                dataset_dir = Path(dataset_dir_str)
                dataset_metadata_file = dataset_dir / 'dataset_metadata.json'
                if not dataset_metadata_file.exists():
                    return None
                with open(dataset_metadata_file, 'r') as f:
                    dataset_metadata = json.load(f)
                
                individuals = dataset_metadata.get('individuals', [])
                individuals_pedigree = dataset_metadata.get('individuals_pedigree', {})
                
                # Build metadata list from indices
                # NOTE: Legacy format has ordering issues - indices may not match data order
                sample_metadata = []
                for idx in split_indices:
                    if idx < len(individuals):
                        sample_id = individuals[idx]
                        pedigree = individuals_pedigree.get(sample_id, {})
                        sample_metadata.append({
                            "sample_id": sample_id,
                            "superpopulation": pedigree.get("superpopulation", "UNK"),
                            "population": pedigree.get("population", "UNK"),
                            "sex": pedigree.get("sex", 0)
                        })
                    else:
                        sample_metadata.append({
                            "sample_id": f"sample_{idx}",
                            "superpopulation": "UNK",
                            "population": "UNK",
                            "sex": 0
                        })
                
                # Warning about legacy format ordering issues
                console.print(f"[yellow]⚠ Cache formato legado detectado - sample_ids podem não corresponder aos dados![/yellow]")
                console.print(f"[yellow]  Recomendado: Apague o cache e execute novamente para usar formato v2.[/yellow]")
                
                return sample_metadata
                
        except Exception as e:
            console.print(f"[yellow]⚠ Sample metadata: erro ao carregar - {e}[/yellow]")
            return None
    
    def get_sample_id(self, idx: int) -> str:
        """Get sample ID for given index in this split."""
        if self._sample_metadata is not None and idx < len(self._sample_metadata):
            return self._sample_metadata[idx].get('sample_id', f"#{idx + 1}")
        return f"#{idx + 1}"
    
    def get_sample_metadata(self, idx: int) -> Dict:
        """
        Get complete metadata for a sample.
        
        Args:
            idx: Index in this split
            
        Returns:
            Dict with sample_id, superpopulation, population, sex
        """
        if self._sample_metadata is not None and idx < len(self._sample_metadata):
            return self._sample_metadata[idx]
        return {
            "sample_id": f"#{idx + 1}",
            "superpopulation": "UNK",
            "population": "UNK",
            "sex": 0
        }

    def _infer_split_name(self) -> str:
        """
        Tenta inferir automaticamente o split a partir do nome do arquivo.
        Retorna 'train' como fallback seguro.
        """
        stem = self.data_file.stem.lower()
        mapping = {
            'train': 'train',
            'train_data': 'train',
            'training': 'train',
            'val': 'val',
            'validation': 'val',
            'val_data': 'val',
            'test': 'test',
            'test_data': 'test'
        }
        for key, value in mapping.items():
            if key in stem:
                return value
        return 'train'

    def _determine_length_without_loading(self) -> int:
        """
        Resolve o comprimento do dataset sem carregar o arquivo .pt completo.
        Utiliza hint direto, metadata.json ou splits_metadata.json como fallback.
        """
        if self._length_hint is not None:
            return self._length_hint
        
        metadata_path = self.data_file.parent / 'metadata.json'
        if metadata_path.exists():
            try:
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                split_key = f"{self.split_name}_size"
                split_sizes = metadata.get('splits', {})
                if split_key in split_sizes:
                    return split_sizes[split_key]
            except json.JSONDecodeError:
                console.print(f"[yellow]Aviso: metadata.json inválido em {metadata_path}[/yellow]")
        
        splits_path = self.data_file.parent / 'splits_metadata.json'
        if splits_path.exists():
            try:
                with open(splits_path, 'r') as f:
                    splits = json.load(f)
                indices_key = f"{self.split_name}_indices"
                if indices_key in splits:
                    return len(splits[indices_key])
            except json.JSONDecodeError:
                console.print(f"[yellow]Aviso: splits_metadata.json inválido em {splits_path}[/yellow]")
        
        raise RuntimeError(
            f"Não foi possível determinar o tamanho de {self.data_file.name} sem carregar o arquivo."
        )
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, int]:
        """
        Retorna sample do dataset.
        
        Para preload: acesso direto da RAM
        Para lazy: verifica cache, senão carrega do arquivo
        
        Returns:
            Tuple (features, target, idx) onde idx é o índice original no split
        """
        if self.loading_strategy == 'preload':
            # Modo preload: acesso direto
            features, target = self.data[idx]
        
        elif self.loading_strategy == 'lazy':
            # Modo lazy: verificar cache primeiro
            if idx in self._cache:
                features, target = self._cache[idx]
            else:
                # Cache miss: carregar arquivo se necessário
                if not self._data_loaded or self.data is None:
                    self.data = torch.load(self.data_file)
                    self._data_loaded = True
                
                features, target = self.data[idx]
                
                # Adicionar ao cache LRU (guardar referências, não clones)
                if self.cache_size > 0:
                    if len(self._cache) >= self.cache_size:
                        # Remover item mais antigo (LRU)
                        oldest_idx = self._cache_order.pop(0)
                        del self._cache[oldest_idx]
                    
                    # Não clonar - apenas guardar referência ao tensor em self.data
                    # Isso não duplica memória, apenas cria ponteiros
                    self._cache[idx] = (features, target)
                    self._cache_order.append(idx)
        
        # Aplicar tainting em runtime (se habilitado)
        # Apenas para classificação (não para regressão/frog_likelihood)
        if (self.config.get('debug', {}).get('taint_at_runtime', False) and
            self.config['output']['prediction_target'] != 'frog_likelihood'):
            num_classes = len(self.target_to_idx)
            debug_config = self.config.get('debug', {})
            # Extrair classe do target tensor
            if target.ndim == 0:  # Escalar
                target_class = target.item()
            else:
                target_class = target[0].item()
            # Aplicar tainting
            features = taint_sample(
                features, target_class, num_classes,
                taint_type=debug_config.get('taint_type', 'additive'),
                taint_value=debug_config.get('taint_value', 1.0),
                taint_horizontal_size=debug_config.get('taint_horizontal_size', 6),
                taint_vertical_size=debug_config.get('taint_vertical_size', 6),
                taint_horizontal_step=debug_config.get('taint_horizontal_step', 100),
                taint_vertical_step=debug_config.get('taint_vertical_step', 6)
            )
        
        return features, target, idx
    
    def unload_data(self):
        """
        Libera arquivo da memória (apenas para lazy loading).
        Chamado automaticamente após cada época para economizar RAM.
        O cache LRU é mantido para acesso rápido.
        """
        if self.loading_strategy == 'lazy' and self.data is not None:
            # Calcular tamanho aproximado liberado
            if isinstance(self.data, list) and len(self.data) > 0:
                # Estimativa: cada sample ~8.6 MB (66 x 32768 x 4 bytes)
                size_mb = len(self.data) * 8.6
                
                del self.data
                self.data = None
                self._data_loaded = False
                
                import gc
                gc.collect()
                
                # Log apenas em modo verbose - mostra quanto foi liberado
                console.print(f"[dim cyan]  → Liberados ~{size_mb:.0f} MB ({self.data_file.name})[/dim cyan]")
    
    def get_num_classes(self) -> int:
        return len(self.target_to_idx)
    
    def get_input_shape(self) -> Tuple[int, int]:
        """
        Retorna shape da entrada.
        
        Returns:
            Tupla (num_rows, effective_size) para dados 2D
        """
        # Garantir que dados estejam carregados
        if self.loading_strategy == 'lazy' and not self._data_loaded:
            self.data = torch.load(self.data_file)
            self._data_loaded = True
        
        if len(self.data) > 0:
            features_shape = self.data[0][0].shape
            if len(features_shape) == 2:
                return tuple(features_shape)
            else:
                # Fallback para dados 1D antigos (retrocompatibilidade)
                return (1, features_shape[0])
        return (0, 0)
    
    def get_input_size(self) -> int:
        """
        DEPRECATED: Retorna tamanho total da entrada (para retrocompatibilidade).
        """
        num_rows, effective_size = self.get_input_shape()
        return num_rows * effective_size


def load_processed_dataset(
    cache_dir: Path,
    config: Dict
) -> Tuple[CachedProcessedDataset, DataLoader, DataLoader, DataLoader, Dict]:
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
    
    # Para classificação, criar mapeamentos a partir dos nomes salvos (ou reconstruir se necessário)
    metadata_updated = False
    temp_base_dataset = None
    
    if 'class_names' in metadata:
        # Carregar mapeamentos reais do metadata (novo formato)
        class_names = metadata['class_names']
        # Converter chaves de string para int se necessário
        idx_to_target = {int(k): v for k, v in class_names.items()}
        target_to_idx = {v: int(k) for k, v in idx_to_target.items()}
    else:
        # Cache antigo sem class_names: reconstruir a partir do dataset base
        console.print("[yellow]⚠ Cache não contém nomes de classes. Reconstruindo...[/yellow]")
        
        # Carregar dataset base temporariamente apenas para obter as classes
        temp_base_dataset = GenomicLongevityDataset(
            dataset_dir=Path(config['dataset_input']['dataset_dir']),
            load_predictions=False,  # Não precisa carregar predições
            load_sequences=False,
            cache_sequences=False
        )
        
        # Obter classes do metadata do dataset (mesmo procedimento de _create_target_mappings)
        dataset_metadata = temp_base_dataset.dataset_metadata
        prediction_target = config['output']['prediction_target']
        classes_from_metadata = None
        
        if prediction_target == 'superpopulation':
            superpop_dist = dataset_metadata.get('superpopulation_distribution', {})
            if superpop_dist:
                classes_from_metadata = list(superpop_dist.keys())
        elif prediction_target == 'population':
            pop_dist = dataset_metadata.get('population_distribution', {})
            if pop_dist:
                classes_from_metadata = list(pop_dist.keys())
        
        if classes_from_metadata:
            # Criar mapeamentos (ordenados alfabeticamente para consistência)
            sorted_targets = sorted(classes_from_metadata)
            target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
            idx_to_target = {idx: target for target, idx in target_to_idx.items()}
            
            # Atualizar metadata do cache com os nomes das classes
            metadata['class_names'] = idx_to_target
            metadata_updated = True
            
            console.print(f"[green]✓ Nomes de classes reconstruídos[/green]")
            console.print(f"[cyan]Classes: {sorted_targets}[/cyan]")
        else:
            # Último recurso: criar mapeamentos dummy
            console.print("[yellow]⚠ Não foi possível obter nomes de classes. Usando índices.[/yellow]")
            num_classes = metadata.get('num_classes', 0)
            for i in range(num_classes):
                target_to_idx[str(i)] = i
                idx_to_target[i] = str(i)
    
    # Carregar ou reconstruir gene_order
    if 'gene_order' in metadata:
        gene_order = metadata['gene_order']
        console.print(f"[green]Ordem dos genes: {', '.join(gene_order) if gene_order else 'N/A'}[/green]")
    else:
        # Cache antigo sem gene_order: reconstruir a partir do dataset base
        console.print("[yellow]⚠ Cache não contém ordem dos genes. Reconstruindo...[/yellow]")
        
        # Carregar dataset base se ainda não foi carregado
        if temp_base_dataset is None:
            temp_base_dataset = GenomicLongevityDataset(
                dataset_dir=Path(config['dataset_input']['dataset_dir']),
                load_predictions=False,
                load_sequences=False,
                cache_sequences=False
            )
        
        try:
            # Obter ordem dos genes do primeiro sample
            first_sample_data, _ = temp_base_dataset[0]
            gene_order = list(first_sample_data['windows'])
            
            # Atualizar metadata do cache com a ordem dos genes
            metadata['gene_order'] = gene_order
            metadata['tracks_per_gene'] = 6  # rna_seq tem 6 tracks
            metadata_updated = True
            
            console.print(f"[green]✓ Ordem dos genes reconstruída[/green]")
            console.print(f"[cyan]Genes: {', '.join(gene_order)}[/cyan]")
        except Exception as e:
            console.print(f"[yellow]⚠ Não foi possível obter ordem dos genes: {e}[/yellow]")
            gene_order = None
    
    # Carregar ou reconstruir gene_window_metadata
    if 'gene_window_metadata' in metadata and metadata['gene_window_metadata']:
        gene_window_metadata = metadata['gene_window_metadata']
        console.print(f"[green]Metadados de janelas genômicas: {len(gene_window_metadata)} genes[/green]")
    else:
        # Cache antigo sem gene_window_metadata: reconstruir a partir do dataset base
        console.print("[yellow]⚠ Cache não contém metadados de janelas genômicas. Reconstruindo...[/yellow]")
        
        # Obter gene_order se ainda não foi obtido
        if gene_order is None and 'gene_order' in metadata:
            gene_order = metadata['gene_order']
        
        if gene_order:
            window_center_size = config['dataset_input']['window_center_size']
            dataset_dir = Path(config['dataset_input']['dataset_dir'])
            gene_window_metadata = get_gene_window_metadata(dataset_dir, gene_order, window_center_size)
            
            if gene_window_metadata:
                # Atualizar metadata do cache com os metadados de janelas
                metadata['gene_window_metadata'] = gene_window_metadata
                metadata_updated = True
                
                console.print(f"[green]✓ Metadados de janelas genômicas reconstruídos[/green]")
                console.print(f"[cyan]Genes com metadados: {', '.join(gene_window_metadata.keys())}[/cyan]")
            else:
                console.print(f"[yellow]⚠ Não foi possível obter metadados de janelas genômicas[/yellow]")
        else:
            console.print(f"[yellow]⚠ Não foi possível obter metadados de janelas: gene_order não disponível[/yellow]")
    
    # Salvar metadata atualizado se houve mudanças
    if metadata_updated:
        metadata_file = cache_dir / 'metadata.json'
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        console.print(f"[green]✓ Metadata do cache atualizado e salvo[/green]")
    
    split_sizes = metadata.get('splits', {})
    
    # Carregar datasets
    train_dataset = CachedProcessedDataset(
        cache_dir / 'train_data.pt',
        target_to_idx,
        idx_to_target,
        config,
        split_name='train',
        length_hint=split_sizes.get('train_size')
    )
    val_dataset = CachedProcessedDataset(
        cache_dir / 'val_data.pt',
        target_to_idx,
        idx_to_target,
        config,
        split_name='val',
        length_hint=split_sizes.get('val_size')
    )
    test_dataset = CachedProcessedDataset(
        cache_dir / 'test_data.pt',
        target_to_idx,
        idx_to_target,
        config,
        split_name='test',
        length_hint=split_sizes.get('test_size')
    )
    
    # Preparar para criar DataLoaders
    console.print("\n[cyan]⚙️  Criando DataLoaders...[/cyan]")
    batch_size = config['training']['batch_size']
    
    # Forçar batch_size=1 se visualização estiver habilitada
    if config.get('debug', {}).get('enable_visualization', False):
        batch_size = 1
        console.print(f"  • [yellow]Batch size forçado para 1 (visualização habilitada)[/yellow]")
    
    loading_strategy = config.get('data_loading', {}).get('loading_strategy', 'preload').lower()
    if loading_strategy == 'lazy':
        train_workers = 0
        val_test_workers = 0
        persistent_workers = False
        console.print(f"  • [yellow]Lazy loading detectado: num_workers ajustado automaticamente para 0 (evita múltiplas cópias em RAM)[/yellow]")
        console.print(f"  • Batch size: {batch_size}")
    else:
        train_workers = 4
        val_test_workers = 2
        persistent_workers = True
        console.print(f"  • Inicializando workers paralelos (train: {train_workers} workers, val/test: {val_test_workers} workers)")
        console.print(f"  • Batch size: {batch_size}")
        console.print(f"  • Isso pode levar alguns segundos na primeira vez...")
    
    # Função collate para empilhar batches corretamente
    def collate_fn(batch):
        """Empilha batch de tuplas (features, target, idx) em tensors."""
        # Como os dados vêm do cache já processados, apenas empilhar
        features_list, targets_list, indices_list = zip(*batch)
        
        # Empilhar features: (batch_size, num_features)
        features_batch = torch.stack(features_list, dim=0)
        
        # Empilhar targets: (batch_size,)
        targets_batch = torch.stack(targets_list, dim=0)
        
        # Empilhar indices: (batch_size,)
        indices_batch = torch.tensor(indices_list, dtype=torch.long)
        
        return features_batch, targets_batch, indices_batch
    
    # Criar generator para shuffle determinístico (se seed configurada)
    generator = None
    if config['data_split']['random_seed'] is not None:
        generator = torch.Generator()
        generator.manual_seed(config['data_split']['random_seed'])
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=train_workers,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if train_workers > 0 else False,
        generator=generator,
        worker_init_fn=worker_init_fn if train_workers > 0 else None
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=val_test_workers,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if val_test_workers > 0 else False,
        worker_init_fn=worker_init_fn if val_test_workers > 0 else None
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=val_test_workers,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if val_test_workers > 0 else False,
        worker_init_fn=worker_init_fn if val_test_workers > 0 else None
    )
    
    console.print(f"[green]✓ DataLoaders criados com sucesso![/green]\n")
    console.print(f"[green]Dataset splits carregados:[/green]")
    console.print(f"  • Treino: {len(train_dataset)} amostras")
    console.print(f"  • Validação: {len(val_dataset)} amostras")
    console.print(f"  • Teste: {len(test_dataset)} amostras")
    
    # Usar train_dataset como full_dataset (para compatibilidade com código existente)
    # Retornar também o metadata para que prepare_data possa extrair gene_order
    return train_dataset, train_loader, val_loader, test_loader, metadata


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


def prepare_data(config: Dict, experiment_dir: Path) -> Tuple[Any, DataLoader, DataLoader, DataLoader]:
    """
    Prepara datasets e dataloaders.
    Tenta carregar do cache compartilhado se disponível, senão processa e salva.
    
    IMPORTANTE: Normalização é computada APENAS com dados de treino+validação,
    nunca incluindo dados de teste (previne data leakage).
    
    Args:
        config: Configuração do experimento
        experiment_dir: Diretório do experimento (onde salvar resultados)
    
    Returns:
        Tupla (full_dataset, train_loader, val_loader, test_loader)
    """
    # Cache do dataset é compartilhado em datasets/
    dataset_cache_dir = get_dataset_cache_dir(config)
    cache_dir = dataset_cache_dir
    
    # Tentar carregar do cache
    if cache_dir is not None:
        cache_path = Path(cache_dir)
        if cache_path.exists() and validate_cache(cache_path, config):
            console.print(Panel.fit(
                "[bold green]Cache Encontrado![/bold green]\n"
                f"Carregando dataset do cache compartilhado: {cache_path.name}"
            ))
            full_dataset, train_loader, val_loader, test_loader, cache_metadata = load_processed_dataset(cache_path, config)
            
            # Passar gene_order do cache para o config
            if 'gene_order' in cache_metadata and cache_metadata['gene_order']:
                config['dataset_input']['gene_order'] = cache_metadata['gene_order']
                console.print(f"[green]Gene order carregado do cache: {', '.join(cache_metadata['gene_order'])}[/green]")
            
            # Passar gene_window_metadata do cache para o config
            if 'gene_window_metadata' in cache_metadata and cache_metadata['gene_window_metadata']:
                config['dataset_input']['gene_window_metadata'] = cache_metadata['gene_window_metadata']
            
            result = (full_dataset, train_loader, val_loader, test_loader)
            
            # Copiar normalization_params para o diretório do experimento (referência)
            (experiment_dir / 'models').mkdir(exist_ok=True)
            norm_source = cache_path / 'normalization_params.json'
            norm_dest = experiment_dir / 'models' / 'normalization_params.json'
            if norm_source.exists() and not norm_dest.exists():
                shutil.copy(norm_source, norm_dest)
                console.print(f"[green]✓ Parâmetros de normalização copiados para o experimento[/green]")
            
            return result
        elif cache_path.exists():
            console.print(Panel.fit(
                "[bold yellow]Cache Inválido ou Incompatível[/bold yellow]\n"
                "Parâmetros mudaram. Reprocessando dataset..."
            ))
        else:
            dataset_name = generate_dataset_name(config)
            console.print(Panel.fit(
                "[bold cyan]Cache Não Encontrado[/bold cyan]\n"
                f"Primeira execução. Processando e salvando em cache compartilhado:\n"
                f"datasets/{dataset_name}"
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
    
    # ==================================================================================
    # PASSO 1: FAZER SPLIT ANTES DA NORMALIZAÇÃO (prevenir data leakage)
    # ==================================================================================
    console.print("\n[bold cyan]📊 Dividindo Dataset (Train/Val/Test)[/bold cyan]")
    
    total_size = len(base_dataset)
    train_size = int(config['data_split']['train_split'] * total_size)
    val_size = int(config['data_split']['val_split'] * total_size)
    test_size = total_size - train_size - val_size
    
    # Criar índices
    indices = list(range(total_size))
    random_seed = config['data_split']['random_seed']
    
    # random_seed == -1 significa NÃO embaralhar (modo debug)
    # random_seed == None significa embaralhar aleatoriamente (não reprodutível)
    # random_seed >= 0 significa embaralhar com seed (reprodutível)
    if random_seed is not None and random_seed != -1:
        np.random.seed(random_seed)
        np.random.shuffle(indices)
        console.print(f"[cyan]  • Dados embaralhados com seed {random_seed}[/cyan]")
    elif random_seed == -1:
        console.print(f"[yellow]  • MODO DEBUG: Dados NÃO embaralhados (random_seed=-1)[/yellow]")
    else:
        np.random.shuffle(indices)
        console.print(f"[yellow]  • Dados embaralhados aleatoriamente (sem seed)[/yellow]")
    
    train_indices = indices[:train_size]
    val_indices = indices[train_size:train_size + val_size]
    test_indices = indices[train_size + val_size:]
    
    console.print(f"[green]  • Treino: {len(train_indices)} amostras ({len(train_indices)/total_size*100:.1f}%)[/green]")
    console.print(f"[green]  • Validação: {len(val_indices)} amostras ({len(val_indices)/total_size*100:.1f}%)[/green]")
    console.print(f"[green]  • Teste: {len(test_indices)} amostras ({len(test_indices)/total_size*100:.1f}%)[/green]")
    
    # ==================================================================================
    # PASSO 2: COMPUTAR NORMALIZAÇÃO APENAS COM TRAIN+VAL (nunca incluir test!)
    # ==================================================================================
    
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
                    console.print(f"\n[green]✓ Parâmetros de normalização carregados do cache[/green]")
                    
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
                        if normalization_params.get('per_track') or 'track_params' in normalization_params:
                            num_tracks = normalization_params.get('num_tracks', len(normalization_params.get('track_params', [])))
                            console.print(f"  • Per-track log_max ({num_tracks} tracks)")
                        else:
                            console.print(f"  • log1p(max): {normalization_params.get('log_max', 1):.6f}")
                    else:
                        console.print(f"  • Método: {method}")
            except Exception as e:
                console.print(f"[yellow]⚠ Não foi possível carregar parâmetros de normalização: {e}[/yellow]")
                normalization_params = None
    
    # Se não temos parâmetros cached, computar APENAS com train+val
    if normalization_params is None:
        console.print(f"\n[bold yellow]🔬 Computando Normalização (APENAS Train+Val)[/bold yellow]")
        console.print(f"[yellow]   ⚠ Dados de teste NÃO serão usados (prevenir data leakage)[/yellow]")
        
        # Criar subset train+val
        train_val_indices = train_indices + val_indices
        train_val_subset = Subset(base_dataset, train_val_indices)
        
        console.print(f"[cyan]   • Usando {len(train_val_subset)} amostras para normalização[/cyan]")
        console.print(f"[cyan]   • Excluindo {len(test_indices)} amostras de teste[/cyan]")
        
        # Criar dataset temporário apenas para computar normalização
        temp_dataset = ProcessedGenomicDataset(
            base_dataset=train_val_subset,
            config=config,
            normalization_params=None,
            compute_normalization=True
        )
        
        # Extrair parâmetros computados
        normalization_params = temp_dataset.normalization_params
        console.print(f"[green]   ✓ Normalização computada com sucesso[/green]")
    
    # ==================================================================================
    # PASSO 3: APLICAR NORMALIZAÇÃO AO DATASET COMPLETO (com parâmetros de train+val)
    # ==================================================================================
    console.print(f"\n[bold cyan]⚙️  Criando Dataset Processado[/bold cyan]")
    
    processed_dataset = ProcessedGenomicDataset(
        base_dataset=base_dataset,
        config=config,
        normalization_params=normalization_params,
        compute_normalization=False  # Nunca recomputar aqui!
    )
    
    # Salvar parâmetros de normalização no diretório models do experimento (para referência)
    norm_path = experiment_dir / 'models' / 'normalization_params.json'
    norm_path.parent.mkdir(parents=True, exist_ok=True)
    with open(norm_path, 'w') as f:
        json.dump(normalization_params, f, indent=2)
    console.print(f"[green]✓ Parâmetros de normalização salvos em {norm_path.name}[/green]")
    
    # Também salvar no cache_dir para reutilização
    if cache_dir is not None:
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
                'normalization_method': config['dataset_input'].get('normalization_method', 'zscore'),
            },
            'normalization_note': 'Computed using ONLY train+validation samples (test excluded to prevent data leakage)'
        }
        with open(cache_path / 'metadata.json', 'w') as f:
            json.dump(partial_metadata, f, indent=2)
        
        with open(cache_path / 'normalization_params.json', 'w') as f:
            json.dump(normalization_params, f, indent=2)
        
        console.print(f"[green]✓ Parâmetros salvos no cache para reutilização[/green]")
    
    # ==================================================================================
    # PASSO 4: SALVAR CACHE E CRIAR DATALOADERS
    # ==================================================================================
    
    console.print(f"\n[green]Dataset split (já definido):[/green]")
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
        full_dataset, train_loader, val_loader, test_loader, cache_metadata = load_processed_dataset(cache_path, config)
        
        # Passar gene_order do cache para o config
        if 'gene_order' in cache_metadata and cache_metadata['gene_order']:
            config['dataset_input']['gene_order'] = cache_metadata['gene_order']
            console.print(f"[green]Gene order carregado do cache: {', '.join(cache_metadata['gene_order'])}[/green]")
        
        # Passar gene_window_metadata do cache para o config
        if 'gene_window_metadata' in cache_metadata and cache_metadata['gene_window_metadata']:
            config['dataset_input']['gene_window_metadata'] = cache_metadata['gene_window_metadata']
        
        result = (full_dataset, train_loader, val_loader, test_loader)
        
        # Copiar normalization_params para o diretório do experimento (referência)
        (experiment_dir / 'models').mkdir(exist_ok=True)
        norm_source = cache_path / 'normalization_params.json'
        norm_dest = experiment_dir / 'models' / 'normalization_params.json'
        if norm_source.exists():
            shutil.copy(norm_source, norm_dest)
            console.print(f"[green]✓ Parâmetros de normalização copiados para {norm_dest}[/green]")
        
        return result
    
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
    
    loading_strategy = config.get('data_loading', {}).get('loading_strategy', 'preload').lower()
    if loading_strategy == 'lazy':
        train_workers = 0
        val_test_workers = 0
        persistent_workers = False
        console.print(f"  • [yellow]Lazy loading detectado: num_workers ajustado automaticamente para 0 (evita múltiplas cópias em RAM)[/yellow]")
        console.print(f"  • Batch size: {batch_size}")
    else:
        train_workers = 4
        val_test_workers = 2
        persistent_workers = True
        console.print(f"  • Inicializando workers paralelos (train: {train_workers} workers, val/test: {val_test_workers} workers)")
        console.print(f"  • Batch size: {batch_size}")
        console.print(f"  • Isso pode levar alguns segundos na primeira vez...")
    
    # Função collate para empilhar batches corretamente
    def collate_fn(batch):
        """Empilha batch de tuplas (features, target, idx) em tensors."""
        # Como os dados vêm do cache já processados, apenas empilhar
        features_list, targets_list, indices_list = zip(*batch)
        
        # Empilhar features: (batch_size, num_features)
        features_batch = torch.stack(features_list, dim=0)
        
        # Empilhar targets: (batch_size,)
        targets_batch = torch.stack(targets_list, dim=0)
        
        # Empilhar indices: (batch_size,)
        indices_batch = torch.tensor(indices_list, dtype=torch.long)
        
        return features_batch, targets_batch, indices_batch
    
    # Criar generator para shuffle determinístico (se seed configurada)
    generator = None
    if config['data_split']['random_seed'] is not None:
        generator = torch.Generator()
        generator.manual_seed(config['data_split']['random_seed'])
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=train_workers,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if train_workers > 0 else False,
        generator=generator,
        worker_init_fn=worker_init_fn if train_workers > 0 else None
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=val_test_workers,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if val_test_workers > 0 else False,
        worker_init_fn=worker_init_fn if val_test_workers > 0 else None
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=val_test_workers,
        pin_memory=True,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if val_test_workers > 0 else False,
        worker_init_fn=worker_init_fn if val_test_workers > 0 else None
    )
    
    return processed_dataset, train_loader, val_loader, test_loader


def run_test_and_save(
    model: nn.Module,
    loader: DataLoader,
    full_dataset: Any,
    config: Dict,
    device: torch.device,
    dataset_name: str,
    experiment_dir: Path
) -> Dict:
    """
    Executa teste e salva resultados em arquivos.
    
    Args:
        model: Modelo treinado
        loader: DataLoader com dados de teste
        full_dataset: Dataset completo
        config: Configuração
        device: Device (CPU ou GPU)
        dataset_name: Nome do conjunto ('train', 'val', 'test')
        experiment_dir: Diretório do experimento
        
    Returns:
        Dict com métricas do teste
    """
    # Criar tester
    tester = Tester(model, loader, full_dataset, config, device, None, dataset_name.capitalize())
    
    # Executar teste (a saída vai para o console normalmente)
    results = tester.test()
    
    # Converter arrays numpy para listas para serialização JSON
    results_serializable = {}
    for key, value in results.items():
        if isinstance(value, np.ndarray):
            results_serializable[key] = value.tolist()
        elif isinstance(value, (list, tuple)):
            # Converter elementos que possam ser numpy arrays
            results_serializable[key] = [
                item.tolist() if isinstance(item, np.ndarray) else item
                for item in value
            ]
        else:
            results_serializable[key] = value
    
    # Salvar resultados em JSON
    json_file = experiment_dir / f'{dataset_name}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results_serializable, f, indent=2)
    
    console.print(f"[green]✓ Resultados de {dataset_name} salvos em: {json_file}[/green]\n")
    
    return results


def summarize_experiments(config: Dict, sort_by: str = 'test_acc'):
    """
    Sumariza resultados de todos os experimentos e cria gráfico comparativo.
    
    Args:
        config: Configuração (para obter processed_cache_dir)
        sort_by: Métrica(s) para ordenação, separada por vírgula para ordenação composta
                 Ex: 'test_acc' ou 'val_acc,test_acc'
    """
    base_cache_dir = Path(config['dataset_input']['processed_cache_dir'])
    
    if not base_cache_dir.exists():
        console.print(f"[red]Erro: Diretório base não encontrado: {base_cache_dir}[/red]")
        return
    
    # Listar todos os subdiretórios (experimentos)
    experiments = []
    for exp_dir in base_cache_dir.iterdir():
        if exp_dir.is_dir():
            # Verificar se tem resultados
            train_json = exp_dir / 'train_results.json'
            val_json = exp_dir / 'val_results.json'
            test_json = exp_dir / 'test_results.json'
            
            if train_json.exists() and val_json.exists() and test_json.exists():
                try:
                    with open(train_json) as f:
                        train_results = json.load(f)
                    with open(val_json) as f:
                        val_results = json.load(f)
                    with open(test_json) as f:
                        test_results = json.load(f)
                    
                    experiments.append({
                        'name': exp_dir.name,
                        'train_accuracy': train_results.get('accuracy', 0),
                        'val_accuracy': val_results.get('accuracy', 0),
                        'test_accuracy': test_results.get('accuracy', 0)
                    })
                except Exception as e:
                    console.print(f"[yellow]⚠ Erro ao ler resultados de {exp_dir.name}: {e}[/yellow]")
    
    if not experiments:
        console.print("[yellow]Nenhum experimento completo encontrado![/yellow]")
        return
    
    console.print(f"[green]Encontrados {len(experiments)} experimentos completos[/green]")
    
    # Parsear métricas de ordenação (pode ter múltiplas separadas por vírgula)
    sort_metrics = [s.strip() for s in sort_by.split(',')]
    
    # Mapear nomes curtos para nomes de campos
    sort_key_map = {
        'train_acc': 'train_accuracy',
        'val_acc': 'val_accuracy',
        'test_acc': 'test_accuracy'
    }
    
    # Validar todas as métricas
    sort_keys = []
    for metric in sort_metrics:
        if metric not in sort_key_map:
            console.print(f"[red]Erro: Métrica inválida '{metric}'. Use: train_acc, val_acc ou test_acc[/red]")
            return
        sort_keys.append(sort_key_map[metric])
    
    # Ordenar experimentos (maior para menor para todas as métricas)
    # Para ordenação composta, usamos tupla como chave
    def sort_key_func(exp):
        # Retorna tupla negativa para ordenar decrescente (maior primeiro)
        return tuple(-exp[key] for key in sort_keys)
    
    experiments = sorted(experiments, key=sort_key_func)
    
    # Mensagem de ordenação
    if len(sort_metrics) == 1:
        console.print(f"[cyan]Ordenando por: {sort_metrics[0]} (maior para menor)[/cyan]")
    else:
        sort_desc = " → ".join(sort_metrics)
        console.print(f"[cyan]Ordenando por: {sort_desc} (prioridade da esquerda para direita, maior para menor)[/cyan]")
    
    # Criar gráfico
    import matplotlib.pyplot as plt
    import numpy as np
    
    exp_names = [exp['name'] for exp in experiments]
    train_accs = [exp['train_accuracy'] for exp in experiments]
    val_accs = [exp['val_accuracy'] for exp in experiments]
    test_accs = [exp['test_accuracy'] for exp in experiments]
    
    # Configurar gráfico
    x = np.arange(len(exp_names))
    width = 0.25
    
    fig, ax = plt.subplots(figsize=(max(12, len(exp_names) * 2), 8))
    
    bars1 = ax.bar(x - width, train_accs, width, label='Train', color='#1f77b4')
    bars2 = ax.bar(x, val_accs, width, label='Validation', color='#ff7f0e')
    bars3 = ax.bar(x + width, test_accs, width, label='Test', color='#2ca02c')
    
    # Configurar eixos
    ax.set_xlabel('Experimento', fontsize=12, fontweight='bold')
    ax.set_ylabel('Accuracy', fontsize=12, fontweight='bold')
    ax.set_title('Comparação de Experimentos - Accuracy', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(exp_names, rotation=45, ha='right', fontsize=8)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.0)
    
    # Adicionar valores nas barras
    def add_value_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=7)
    
    add_value_labels(bars1)
    add_value_labels(bars2)
    add_value_labels(bars3)
    
    plt.tight_layout()
    
    # Salvar gráfico
    output_path = base_cache_dir / 'experiments_summary.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    console.print(f"[green]✓ Gráfico salvo em: {output_path}[/green]")
    
    # Mostrar gráfico
    plt.show()
    
    # Imprimir tabela resumo
    table = Table(title="Resumo dos Experimentos")
    table.add_column("Experimento", style="cyan")
    table.add_column("Train Acc", justify="right", style="blue")
    table.add_column("Val Acc", justify="right", style="yellow")
    table.add_column("Test Acc", justify="right", style="green")
    
    for exp in experiments:
        table.add_row(
            exp['name'],
            f"{exp['train_accuracy']:.4f}",
            f"{exp['val_accuracy']:.4f}",
            f"{exp['test_accuracy']:.4f}"
        )
    
    console.print(table)


def worker_init_fn(worker_id: int):
    """
    Inicializa a seed de cada worker do DataLoader de forma determinística.
    
    Essencial para reprodutibilidade com num_workers > 0 e persistent_workers=True.
    Sem isso, cada worker pode ter seeds aleatórias diferentes em cada execução.
    
    Args:
        worker_id: ID do worker (0, 1, 2, ...)
    """
    import random
    
    # Obter a seed base do PyTorch (configurada por set_random_seeds)
    worker_seed = torch.initial_seed() % 2**32
    
    # Configurar seeds específicas para este worker
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def set_random_seeds(seed: int, strict_determinism: bool = True):
    """
    Configura todas as sementes randômicas para reprodutibilidade.
    
    Args:
        seed: Valor da semente randômica
        strict_determinism: Se True, garante determinismo total (mais lento).
                          Se False, determinismo parcial (mais rápido, ~99% reprodutível).
    """
    import random
    import os
    
    # Python hash seed (para dicts e sets)
    os.environ['PYTHONHASHSEED'] = str(seed)
    
    # Python random
    random.seed(seed)
    
    # NumPy
    np.random.seed(seed)
    
    # PyTorch CPU
    torch.manual_seed(seed)
    
    # PyTorch GPU
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # Para multi-GPU
        
        # Limpar cache da GPU para eliminar estado residual
        torch.cuda.empty_cache()
    
    # Configurar determinismo estrito do PyTorch
    if strict_determinism:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        
        # Forçar algoritmos determinísticos (PyTorch 1.8+)
        try:
            torch.use_deterministic_algorithms(True)
        except AttributeError:
            # PyTorch < 1.8
            torch.set_deterministic(True)
        
        # Configurar cuBLAS para operações determinísticas
        os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'
        
        console.print(f"[green]🎲 Semente randômica configurada: {seed} (determinismo ESTRITO - 100% reprodutível)[/green]")
        console.print(f"[yellow]   ⚠ Treinamento pode ser 10-30% mais lento devido ao determinismo estrito[/yellow]")
    else:
        # Permite operações não-determinísticas para melhor performance
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True
        
        try:
            torch.use_deterministic_algorithms(False)
        except AttributeError:
            try:
                torch.set_deterministic(False)
            except AttributeError:
                pass
        
        console.print(f"[green]🎲 Semente randômica configurada: {seed} (determinismo PARCIAL - ~99% reprodutível)[/green]")
        console.print(f"[green]   ✓ Performance otimizada, pequenas variações podem ocorrer[/green]")


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
    parser.add_argument(
        '--summarize_results',
        action='store_true',
        help='Sumariza resultados de todos os experimentos e gera gráfico comparativo'
    )
    parser.add_argument(
        '--sort_by',
        type=str,
        default='test_acc',
        help='Métrica(s) para ordenar experimentos (padrão: test_acc). '
             'Aceita ordenação simples (ex: val_acc) ou composta separada por vírgula '
             '(ex: val_acc,test_acc ordena por val_acc, depois test_acc como desempate)'
    )
    
    args = parser.parse_args()
    
    # Carregar configuração
    config = load_config(Path(args.config))
    
    # Limpeza agressiva de estado CUDA ANTES de qualquer coisa
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
        try:
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.reset_accumulated_memory_stats()
        except:
            pass  # Versões antigas do PyTorch podem não ter essas funções
        console.print(f"[green]✓ Estado CUDA limpo (pré-inicialização)[/green]")
    
    # Configurar semente randômica para reprodutibilidade (ANTES de qualquer operação)
    random_seed = config['data_split']['random_seed']
    if random_seed is not None and random_seed != -1:
        strict_determinism = config['data_split'].get('strict_determinism', True)
        set_random_seeds(random_seed, strict_determinism)
    elif random_seed == -1:
        console.print("[yellow]⚠ MODO DEBUG: random_seed=-1 (dados não serão embaralhados)[/yellow]")
        # Ainda configurar determinismo para outras operações
        strict_determinism = config['data_split'].get('strict_determinism', True)
        set_random_seeds(0, strict_determinism)  # Usar seed 0 para outras operações
    else:
        console.print("[yellow]⚠ Semente randômica não configurada - resultados não serão reprodutíveis[/yellow]")
    
    # Limpeza adicional APÓS configurar seeds
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
        console.print(f"[green]✓ Estado CUDA re-sincronizado (pós-seeds)[/green]")
    
    # Se flag --summarize_results, executar e sair
    if args.summarize_results:
        summarize_experiments(config, sort_by=args.sort_by)
        return
    
    # Sobrescrever modo se fornecido
    if args.mode:
        config['mode'] = args.mode
    
    # Configurar signal handler para CTRL+C (diferente para train vs test)
    if config['mode'] == 'train':
        def signal_handler(sig, frame):
            """Handler para capturar CTRL+C durante treinamento."""
            console.print("\n[yellow]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/yellow]")
            console.print("[yellow]⚠ CTRL+C detectado - Finalizando treino graciosamente...[/yellow]")
            console.print("[yellow]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/yellow]")
            interrupt_state.interrupted = True
    else:
        def signal_handler(sig, frame):
            """Handler para capturar CTRL+C durante teste - interrompe imediatamente."""
            console.print("\n[yellow]⚠ CTRL+C detectado - Interrompendo...[/yellow]")
            sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    
    # Banner
    banner_text = (
        "[bold cyan]Neural Ancestry Predictor[/bold cyan]\n"
        f"Modo: {config['mode']}\n"
    )
    
    # Adicionar info do conjunto de teste se modo for test
    if config['mode'] == 'test':
        test_dataset = config.get('test_dataset', 'test')
        banner_text += f"Conjunto: {test_dataset}\n"
    
    banner_text += (
        f"Target: {config['output']['prediction_target']}\n"
        f"Config: {args.config}"
    )
    
    console.print(Panel.fit(banner_text, title="🧬 Genomics"))
    
    # Configurar device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    console.print(f"[green]Device: {device}[/green]")
    
    # Setup do diretório do experimento (apenas para treino)
    if config['mode'] == 'train':
        experiment_dir = setup_experiment_dir(config, args.config)
    else:
        # Para teste, precisa reconstruir o experiment_dir a partir dos parâmetros
        experiment_name = generate_experiment_name(config)
        base_cache_dir = Path(config['dataset_input']['processed_cache_dir'])
        experiment_dir = base_cache_dir / experiment_name
        
        if not experiment_dir.exists():
            console.print(f"[red]Erro: Experimento não encontrado: {experiment_dir}[/red]")
            console.print("[yellow]Execute o treinamento primeiro![/yellow]")
            sys.exit(1)
    
    # Limpeza CUDA antes de preparar dados
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
    
    # Preparar dados
    full_dataset, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    
    # Calcular dimensões de entrada
    input_shape = full_dataset.get_input_shape()
    num_classes = full_dataset.get_num_classes() if config['output']['prediction_target'] != 'frog_likelihood' else 150
    
    console.print(f"[green]Input shape: {input_shape[0]} x {input_shape[1]} (2D)[/green]")
    console.print(f"[green]Number of classes: {num_classes}[/green]")
    
    # Limpeza CUDA antes de criar modelo
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
        console.print(f"[green]✓ Estado CUDA limpo (pré-modelo)[/green]")
    
    # Criar modelo baseado no tipo configurado
    model_type = config['model'].get('type', 'NN').upper()
    
    if model_type == 'NN':
        console.print(f"[cyan]Criando modelo: Neural Network (NN) totalmente conectada[/cyan]")
        model = NNAncestryPredictor(config, input_shape, num_classes).to(device)
    elif model_type == 'CNN':
        console.print(f"[cyan]Criando modelo: Convolutional Neural Network (CNN)[/cyan]")
        model = CNNAncestryPredictor(config, input_shape, num_classes).to(device)
    elif model_type == 'CNN2':
        console.print(f"[cyan]Criando modelo: CNN2 (Multi-layer with Global Pooling)[/cyan]")
        model = CNN2AncestryPredictor(config, input_shape, num_classes).to(device)
    else:
        raise ValueError(f"Tipo de modelo não suportado: {model_type}. Use 'NN', 'CNN' ou 'CNN2'.")
    
    # Inicializar W&B
    wandb_run = None
    if config['wandb']['use_wandb']:
        try:
            import wandb
            
            # Se run_name não for especificado, usar o nome do experimento
            run_name = config['wandb'].get('run_name')
            if run_name is None:
                run_name = experiment_dir.name  # Nome do diretório do experimento
            
            wandb_run = wandb.init(
                project=config['wandb']['project_name'],
                name=run_name,
                config=config
            )
            console.print("[green]✓ Weights & Biases inicializado[/green]")
            console.print(f"  • Run name: {run_name}")
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
                
                # Limpar memória GPU e carregar na CPU primeiro
                if device.type == 'cuda':
                    torch.cuda.empty_cache()
                
                checkpoint = torch.load(checkpoint_path, map_location='cpu')
                model.load_state_dict(checkpoint['model_state_dict'])
        
        # Limpeza CUDA final antes de treinar
        if device.type == 'cuda':
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
            console.print(f"[green]✓ Estado CUDA limpo (pré-treino)[/green]")
        
        # Treinar
        trainer = Trainer(model, train_loader, val_loader, config, device, experiment_dir, wandb_run)
        history = trainer.train()
        
        # Salvar histórico
        history_path = experiment_dir / 'models' / 'training_history.json'
        with open(history_path, 'w') as f:
            json.dump(history, f, indent=2)
        console.print(f"[green]✓ Histórico salvo em {history_path}[/green]")
        
        # Executar testes automáticos após o treino (executar MESMO se CTRL+C)
        console.print("\n[bold cyan]═══════════════════════════════════════════════[/bold cyan]")
        if interrupt_state.interrupted:
            console.print("[bold cyan]Executando Testes Após CTRL+C[/bold cyan]")
        else:
            console.print("[bold cyan]Executando Testes Automáticos Após Treinamento[/bold cyan]")
        console.print("[bold cyan]═══════════════════════════════════════════════[/bold cyan]\n")
        
        # Determinar qual checkpoint usar
        models_dir = experiment_dir / 'models'
        if (models_dir / 'best_accuracy.pt').exists():
            checkpoint_path = models_dir / 'best_accuracy.pt'
            console.print(f"[green]Carregando melhor modelo (accuracy): {checkpoint_path}[/green]")
        elif (models_dir / 'final.pt').exists():
            checkpoint_path = models_dir / 'final.pt'
            console.print(f"[green]Carregando modelo final: {checkpoint_path}[/green]")
        else:
            console.print("[yellow]⚠ Nenhum checkpoint encontrado para teste automático[/yellow]")
            checkpoint_path = None
        
        if checkpoint_path:
            # Limpar memória GPU antes de carregar checkpoint
            if device.type == 'cuda':
                torch.cuda.empty_cache()
            
            # Carregar checkpoint na CPU, depois aplicar ao modelo (já na GPU)
            checkpoint = torch.load(checkpoint_path, map_location='cpu')
            model.load_state_dict(checkpoint['model_state_dict'])
            
            # Teste no conjunto de treino
            console.print("\n[cyan]━━━ Testando no conjunto de TREINO ━━━[/cyan]")
            train_results = run_test_and_save(model, train_loader, full_dataset, config, device, 'train', experiment_dir)
            
            # Teste no conjunto de validação
            console.print("\n[cyan]━━━ Testando no conjunto de VALIDAÇÃO ━━━[/cyan]")
            val_results = run_test_and_save(model, val_loader, full_dataset, config, device, 'val', experiment_dir)
            
            # Teste no conjunto de teste
            console.print("\n[cyan]━━━ Testando no conjunto de TESTE ━━━[/cyan]")
            test_results = run_test_and_save(model, test_loader, full_dataset, config, device, 'test', experiment_dir)
            
            console.print("\n[bold green]✓ Testes automáticos concluídos![/bold green]")
        
    elif config['mode'] == 'test':
        # Carregar melhor checkpoint do experiment_dir
        models_dir = experiment_dir / 'models'
        checkpoint_path = models_dir / 'best_accuracy.pt'
        
        if not checkpoint_path.exists():
            # Tentar final.pt como fallback
            checkpoint_path = models_dir / 'final.pt'
            if not checkpoint_path.exists():
                console.print(f"[red]ERRO: Nenhum checkpoint encontrado em: {models_dir}[/red]")
                sys.exit(1)
        
        console.print(f"[yellow]Carregando checkpoint: {checkpoint_path}[/yellow]")
        
        # Limpar memória GPU antes de carregar checkpoint
        if device.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Carregar checkpoint na CPU, depois aplicar ao modelo (já na GPU)
        checkpoint = torch.load(checkpoint_path, map_location='cpu')
        model.load_state_dict(checkpoint['model_state_dict'])
        
        # Selecionar conjunto de dados para teste
        test_dataset_choice = config.get('test_dataset', 'test').lower()
        
        if test_dataset_choice == 'train':
            selected_loader = train_loader
            dataset_name = "Train"
        elif test_dataset_choice == 'val':
            selected_loader = val_loader
            dataset_name = "Validation"
        else:  # 'test' is the default
            selected_loader = test_loader
            dataset_name = "Test"
        
        console.print(f"[cyan]Testando no conjunto de: {dataset_name}[/cyan]")
        
        # Limpeza CUDA final antes de testar
        if device.type == 'cuda':
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
            console.print(f"[green]✓ Estado CUDA limpo (pré-teste)[/green]")
        
        # Testar
        tester = Tester(model, selected_loader, full_dataset, config, device, wandb_run, dataset_name)
        results = tester.test()
    
    # Finalizar W&B
    if wandb_run:
        wandb_run.finish()
    
    console.print("\n[bold green]✓ Execução concluída![/bold green]")


if __name__ == '__main__':
    main()

