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

# Importar dataset genômico
sys.path.insert(0, str(Path(__file__).parent.parent / "build_non_longevous_dataset"))
from genomic_dataset import GenomicLongevityDataset

console = Console()

# ==============================================================================
# DATA PROCESSING FUNCTIONS
# ==============================================================================

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
                    input_data, _ = self.base_dataset[idx]
                    processed_array = self._process_windows(input_data['windows'])
                    all_values.append(processed_array)
                except Exception as e:
                    console.print(f"[yellow]Aviso: Erro ao processar amostra {idx}: {e}[/yellow]")
                progress.update(task, advance=1)
        
        if len(all_values) == 0:
            console.print("[red]ERRO: Nenhuma amostra válida para normalização![/red]")
            return {'mean': 0.0, 'std': 1.0}
        
        # Concatenar todos os arrays
        all_values = np.concatenate(all_values)
        
        mean = float(np.mean(all_values))
        std = float(np.std(all_values))
        
        # Evitar divisão por zero
        if std < 1e-8:
            std = 1.0
        
        console.print(f"[green]Normalização: mean={mean:.6f}, std={std:.6f}[/green]")
        
        return {'mean': mean, 'std': std}
    
    def _create_target_mappings(self):
        """Cria mapeamentos entre targets e índices."""
        unique_targets = set()
        
        for idx in range(len(self.base_dataset)):
            try:
                _, output_data = self.base_dataset[idx]
                target = self._get_target_value(output_data)
                if target is not None:
                    unique_targets.add(target)
            except Exception:
                continue
        
        # Ordenar para consistência
        sorted_targets = sorted(list(unique_targets))
        
        self.target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
        self.idx_to_target = {idx: target for target, idx in self.target_to_idx.items()}
        
        console.print(f"[green]Mapeamento de targets criado: {len(self.target_to_idx)} classes[/green]")
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
        return np.concatenate(processed_arrays)
    
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
        
        # Normalizar
        features = (features - self.normalization_params['mean']) / self.normalization_params['std']
        
        # Converter para tensor
        features_tensor = torch.FloatTensor(features)
        
        # Processar target
        if self.prediction_target == 'frog_likelihood':
            # Regressão: usar likelihood diretamente
            target = output_data.get('frog_likelihood', np.zeros(150))
            target_tensor = torch.FloatTensor(target)
        else:
            # Classificação: converter para índice
            target_value = self._get_target_value(output_data)
            if target_value in self.target_to_idx:
                target_idx = self.target_to_idx[target_value]
                target_tensor = torch.LongTensor([target_idx])
            else:
                # Target desconhecido, usar -1
                target_tensor = torch.LongTensor([-1])
        
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
        
        # Escolher função de ativação
        if activation_type == 'relu':
            self.activation = nn.ReLU()
        elif activation_type == 'tanh':
            self.activation = nn.Tanh()
        elif activation_type == 'sigmoid':
            self.activation = nn.Sigmoid()
        else:
            raise ValueError(f"Ativação não suportada: {activation_type}")
        
        # Construir camadas
        layers = []
        prev_size = input_size
        
        for hidden_size in hidden_layers:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(self.activation)
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            prev_size = hidden_size
        
        # Camada de saída
        layers.append(nn.Linear(prev_size, num_classes))
        
        # Softmax para classificação (aplicado no forward)
        if self.is_classification:
            self.softmax = nn.Softmax(dim=1)
        
        self.network = nn.Sequential(*layers)
        
        console.print(f"[green]Modelo criado:[/green]")
        console.print(f"  • Input size: {input_size}")
        console.print(f"  • Hidden layers: {hidden_layers}")
        console.print(f"  • Output size: {num_classes}")
        console.print(f"  • Activation: {activation_type}")
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
                features = features.to(self.device)
                targets = targets.to(self.device)
                
                # Squeeze targets se necessário
                if targets.dim() > 1 and targets.shape[1] == 1:
                    targets = targets.squeeze(1)
                
                # Forward pass
                self.optimizer.zero_grad()
                outputs = self.model(features)
                loss = self.criterion(outputs, targets)
                
                # Backward pass
                loss.backward()
                self.optimizer.step()
                
                total_loss += loss.item()
                num_batches += 1
                
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
                features = features.to(self.device)
                targets = targets.to(self.device)
                
                # Squeeze targets se necessário
                if targets.dim() > 1 and targets.shape[1] == 1:
                    targets = targets.squeeze(1)
                
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
                    features = features.to(self.device)
                    targets = targets.to(self.device)
                    
                    if targets.dim() > 1 and targets.shape[1] == 1:
                        targets = targets.squeeze(1)
                    
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
            results['accuracy'] = accuracy_score(all_targets, all_predictions)
            results['precision'], results['recall'], results['f1'], _ = \
                precision_recall_fscore_support(all_targets, all_predictions, average='weighted', zero_division=0)
            results['confusion_matrix'] = confusion_matrix(all_targets, all_predictions)
            
            # Classification report
            target_names = [self.dataset.idx_to_target[i] for i in range(self.dataset.get_num_classes())]
            results['classification_report'] = classification_report(
                all_targets, all_predictions, target_names=target_names, zero_division=0
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
            'dataset_dir': config['dataset_input']['dataset_dir'],
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
    Salva dataset processado em cache.
    
    Args:
        cache_dir: Diretório onde salvar cache
        processed_dataset: Dataset processado
        train_indices: Índices de treino
        val_indices: Índices de validação
        test_indices: Índices de teste
        config: Configuração usada
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    console.print(f"[cyan]Salvando cache em {cache_dir}...[/cyan]")
    
    # Preparar dados de cada split
    train_data = []
    val_data = []
    test_data = []
    
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
            
            if idx in train_indices:
                train_data.append((features, target))
            elif idx in val_indices:
                val_data.append((features, target))
            elif idx in test_indices:
                test_data.append((features, target))
            
            progress.update(task, advance=1)
    
    # Salvar dados
    torch.save(train_data, cache_dir / 'train_data.pt')
    torch.save(val_data, cache_dir / 'val_data.pt')
    torch.save(test_data, cache_dir / 'test_data.pt')
    
    # Salvar splits
    splits = {
        'train_indices': train_indices,
        'val_indices': val_indices,
        'test_indices': test_indices
    }
    with open(cache_dir / 'splits.json', 'w') as f:
        json.dump(splits, f, indent=2)
    
    # Salvar normalization params
    with open(cache_dir / 'normalization_params.json', 'w') as f:
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
    with open(cache_dir / 'metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    console.print(f"[green]✓ Cache salvo com sucesso em {cache_dir}[/green]")
    console.print(f"  • Train: {len(train_data)} samples")
    console.print(f"  • Val: {len(val_data)} samples")
    console.print(f"  • Test: {len(test_data)} samples")


class CachedProcessedDataset(Dataset):
    """
    Dataset wrapper para dados carregados do cache.
    """
    
    def __init__(self, data_file: Path, target_to_idx: Dict, idx_to_target: Dict):
        """
        Inicializa dataset do cache.
        
        Args:
            data_file: Arquivo .pt com dados processados
            target_to_idx: Mapeamento target->índice
            idx_to_target: Mapeamento índice->target
        """
        console.print(f"[cyan]Carregando {data_file.name}...[/cyan]")
        self.data = torch.load(data_file)
        self.target_to_idx = target_to_idx
        self.idx_to_target = idx_to_target
        console.print(f"[green]✓ {len(self.data)} samples carregados[/green]")
    
    def __len__(self) -> int:
        return len(self.data)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.data[idx]
    
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
    
    console.print(f"[green]Normalização: mean={norm_params['mean']:.6f}, std={norm_params['std']:.6f}[/green]")
    
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
    train_dataset = CachedProcessedDataset(cache_dir / 'train_data.pt', target_to_idx, idx_to_target)
    val_dataset = CachedProcessedDataset(cache_dir / 'val_data.pt', target_to_idx, idx_to_target)
    test_dataset = CachedProcessedDataset(cache_dir / 'test_data.pt', target_to_idx, idx_to_target)
    
    # Criar dataloaders
    batch_size = config['training']['batch_size']
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=0
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0
    )
    
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
    
    # Criar dataset processado
    processed_dataset = ProcessedGenomicDataset(
        base_dataset=base_dataset,
        config=config,
        normalization_params=None,
        compute_normalization=True
    )
    
    # Salvar parâmetros de normalização no diretório de checkpoints (para referência)
    norm_path = Path(config['checkpointing']['checkpoint_dir']) / 'normalization_params.json'
    norm_path.parent.mkdir(parents=True, exist_ok=True)
    with open(norm_path, 'w') as f:
        json.dump(processed_dataset.normalization_params, f, indent=2)
    console.print(f"[green]✓ Parâmetros de normalização salvos em {norm_path}[/green]")
    
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
    
    # Criar subsets
    train_dataset = Subset(processed_dataset, train_indices)
    val_dataset = Subset(processed_dataset, val_indices)
    test_dataset = Subset(processed_dataset, test_indices)
    
    # Criar dataloaders
    batch_size = config['training']['batch_size']
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=0
    )
    
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0
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

