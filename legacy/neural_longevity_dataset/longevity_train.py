#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script de Treinamento - Marcadores de Longevidade

Treina uma rede neural profunda para classificar longevidade
usando features genômicas e predições do AlphaGenome.
"""

import argparse
import sys
from pathlib import Path
import json
import yaml

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, roc_auc_score, confusion_matrix

import numpy as np
from tqdm import tqdm
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from neural_longevity_dataset import LongevityDataset

console = Console()


# ═══════════════════════════════════════════════════════════════════
# Modelos
# ═══════════════════════════════════════════════════════════════════

class LongevityCNN(nn.Module):
    """
    CNN 1D para classificação de longevidade.
    
    Processa sequência DNA com convoluções e funde com features AlphaGenome.
    """
    
    def __init__(self, sequence_length: int = 2048, alphagenome_features: int = 35):
        super().__init__()
        
        # DNA sequence processing (Conv1D)
        self.conv1 = nn.Conv1d(4, 64, kernel_size=7, padding=3)
        self.bn1 = nn.BatchNorm1d(64)
        self.pool1 = nn.MaxPool1d(2)
        
        self.conv2 = nn.Conv1d(64, 128, kernel_size=5, padding=2)
        self.bn2 = nn.BatchNorm1d(128)
        self.pool2 = nn.MaxPool1d(2)
        
        self.conv3 = nn.Conv1d(128, 256, kernel_size=3, padding=1)
        self.bn3 = nn.BatchNorm1d(256)
        
        # Global pooling
        self.global_pool = nn.AdaptiveMaxPool1d(1)
        
        # AlphaGenome features processing
        self.fc_alpha1 = nn.Linear(alphagenome_features, 128)
        self.bn_alpha1 = nn.BatchNorm1d(128)
        self.fc_alpha2 = nn.Linear(128, 128)
        self.bn_alpha2 = nn.BatchNorm1d(128)
        
        # Fusion and classification
        self.fc1 = nn.Linear(256 + 128 + 1, 512)  # seq + alpha + position
        self.bn_fc1 = nn.BatchNorm1d(512)
        self.dropout1 = nn.Dropout(0.5)
        
        self.fc2 = nn.Linear(512, 256)
        self.bn_fc2 = nn.BatchNorm1d(256)
        self.dropout2 = nn.Dropout(0.3)
        
        self.fc3 = nn.Linear(256, 1)
    
    def forward(self, sequence, position, alphagenome):
        """
        Args:
            sequence: (batch, 4, seq_len) - one-hot encoded DNA
            position: (batch, 1) - normalized position
            alphagenome: (batch, n_features) - AlphaGenome features
            
        Returns:
            (batch, 1) - probability of longevity
        """
        # Process DNA sequence
        x_seq = F.relu(self.bn1(self.conv1(sequence)))
        x_seq = self.pool1(x_seq)
        
        x_seq = F.relu(self.bn2(self.conv2(x_seq)))
        x_seq = self.pool2(x_seq)
        
        x_seq = F.relu(self.bn3(self.conv3(x_seq)))
        x_seq = self.global_pool(x_seq).squeeze(-1)  # (batch, 256)
        
        # Process AlphaGenome features
        x_alpha = F.relu(self.bn_alpha1(self.fc_alpha1(alphagenome)))
        x_alpha = F.relu(self.bn_alpha2(self.fc_alpha2(x_alpha)))  # (batch, 128)
        
        # Concatenate all features
        x = torch.cat([x_seq, x_alpha, position], dim=1)
        
        # Classification layers
        x = F.relu(self.bn_fc1(self.fc1(x)))
        x = self.dropout1(x)
        
        x = F.relu(self.bn_fc2(self.fc2(x)))
        x = self.dropout2(x)
        
        x = torch.sigmoid(self.fc3(x))
        
        return x


# ═══════════════════════════════════════════════════════════════════
# Treinamento
# ═══════════════════════════════════════════════════════════════════

class LongevityTrainer:
    """Orquestrador de treinamento."""
    
    def __init__(self, config: dict):
        self.config = config
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        console.print(f"[cyan]Dispositivo: {self.device}[/cyan]")
        
        # Carregar datasets
        self.train_dataset = LongevityDataset(Path(config['data']['train_file']))
        self.val_dataset = LongevityDataset(Path(config['data']['val_file']))
        self.test_dataset = LongevityDataset(Path(config['data']['test_file']))
        
        # DataLoaders
        self.train_loader = DataLoader(
            self.train_dataset,
            batch_size=config['training']['batch_size'],
            shuffle=True,
            num_workers=config['training']['num_workers']
        )
        
        self.val_loader = DataLoader(
            self.val_dataset,
            batch_size=config['training']['batch_size'],
            shuffle=False,
            num_workers=config['training']['num_workers']
        )
        
        self.test_loader = DataLoader(
            self.test_dataset,
            batch_size=config['training']['batch_size'],
            shuffle=False,
            num_workers=config['training']['num_workers']
        )
        
        # Modelo
        self.model = LongevityCNN(
            sequence_length=config['model']['sequence_length'],
            alphagenome_features=config['model']['alphagenome_features']
        ).to(self.device)
        
        console.print(f"[green]✓ Modelo criado: {sum(p.numel() for p in self.model.parameters()):,} parâmetros[/green]")
        
        # Loss e otimizador
        self.criterion = nn.BCELoss()
        self.optimizer = optim.Adam(
            self.model.parameters(),
            lr=config['training']['learning_rate'],
            weight_decay=config['training']['weight_decay']
        )
        
        # Learning rate scheduler
        self.scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer,
            mode='min',
            factor=0.5,
            patience=5,
            verbose=True
        )
        
        # Histórico
        self.history = {
            'train_loss': [],
            'val_loss': [],
            'val_accuracy': [],
            'val_f1': []
        }
        
        self.best_val_loss = float('inf')
        self.best_model_path = Path(config['output']['model_dir']) / 'best_model.pt'
        self.best_model_path.parent.mkdir(parents=True, exist_ok=True)
    
    def train_epoch(self) -> float:
        """Treina uma época."""
        self.model.train()
        total_loss = 0.0
        
        for features, labels in tqdm(self.train_loader, desc="Treinando"):
            # Move para device
            sequence = features['sequence'].to(self.device)
            position = features['position'].to(self.device)
            alphagenome = features['alphagenome'].to(self.device)
            labels = labels.float().unsqueeze(1).to(self.device)
            
            # Forward
            self.optimizer.zero_grad()
            outputs = self.model(sequence, position, alphagenome)
            loss = self.criterion(outputs, labels)
            
            # Backward
            loss.backward()
            self.optimizer.step()
            
            total_loss += loss.item()
        
        return total_loss / len(self.train_loader)
    
    def validate(self) -> tuple:
        """Valida o modelo."""
        self.model.eval()
        total_loss = 0.0
        all_preds = []
        all_labels = []
        all_probs = []
        
        with torch.no_grad():
            for features, labels in tqdm(self.val_loader, desc="Validando"):
                sequence = features['sequence'].to(self.device)
                position = features['position'].to(self.device)
                alphagenome = features['alphagenome'].to(self.device)
                labels_t = labels.float().unsqueeze(1).to(self.device)
                
                outputs = self.model(sequence, position, alphagenome)
                loss = self.criterion(outputs, labels_t)
                
                total_loss += loss.item()
                
                preds = (outputs > 0.5).cpu().numpy().flatten()
                probs = outputs.cpu().numpy().flatten()
                
                all_preds.extend(preds)
                all_labels.extend(labels.numpy())
                all_probs.extend(probs)
        
        # Métricas
        avg_loss = total_loss / len(self.val_loader)
        accuracy = accuracy_score(all_labels, all_preds)
        precision, recall, f1, _ = precision_recall_fscore_support(
            all_labels, all_preds, average='binary'
        )
        auc_roc = roc_auc_score(all_labels, all_probs)
        
        return avg_loss, accuracy, precision, recall, f1, auc_roc
    
    def train(self):
        """Treinamento completo."""
        console.print(Panel.fit(
            "[bold cyan]Iniciando Treinamento[/bold cyan]\n"
            f"Épocas: {self.config['training']['epochs']}\n"
            f"Batch size: {self.config['training']['batch_size']}\n"
            f"Learning rate: {self.config['training']['learning_rate']}",
            border_style="cyan"
        ))
        
        for epoch in range(1, self.config['training']['epochs'] + 1):
            console.print(f"\n[bold cyan]Época {epoch}/{self.config['training']['epochs']}[/bold cyan]")
            
            # Treinar
            train_loss = self.train_epoch()
            
            # Validar
            val_loss, accuracy, precision, recall, f1, auc_roc = self.validate()
            
            # Atualizar scheduler
            self.scheduler.step(val_loss)
            
            # Histórico
            self.history['train_loss'].append(train_loss)
            self.history['val_loss'].append(val_loss)
            self.history['val_accuracy'].append(accuracy)
            self.history['val_f1'].append(f1)
            
            # Exibir métricas
            table = Table(title=f"Época {epoch}")
            table.add_column("Métrica", style="cyan")
            table.add_column("Valor", style="green")
            
            table.add_row("Train Loss", f"{train_loss:.4f}")
            table.add_row("Val Loss", f"{val_loss:.4f}")
            table.add_row("Accuracy", f"{accuracy:.4f}")
            table.add_row("Precision", f"{precision:.4f}")
            table.add_row("Recall", f"{recall:.4f}")
            table.add_row("F1-Score", f"{f1:.4f}")
            table.add_row("AUC-ROC", f"{auc_roc:.4f}")
            
            console.print(table)
            
            # Salvar melhor modelo
            if val_loss < self.best_val_loss:
                self.best_val_loss = val_loss
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': self.model.state_dict(),
                    'optimizer_state_dict': self.optimizer.state_dict(),
                    'val_loss': val_loss,
                    'metrics': {
                        'accuracy': accuracy,
                        'precision': precision,
                        'recall': recall,
                        'f1': f1,
                        'auc_roc': auc_roc
                    }
                }, self.best_model_path)
                console.print(f"[green]✓ Melhor modelo salvo![/green]")
        
        console.print("\n[bold green]✓ Treinamento concluído![/bold green]")
        
        # Salvar histórico
        history_file = Path(self.config['output']['model_dir']) / 'history.json'
        with open(history_file, 'w') as f:
            json.dump(self.history, f, indent=2)
        
        console.print(f"[green]✓ Histórico salvo: {history_file}[/green]")
    
    def test(self):
        """Avalia no conjunto de teste."""
        console.print("\n[bold cyan]Avaliação no Conjunto de Teste[/bold cyan]")
        
        # Carregar melhor modelo
        checkpoint = torch.load(self.best_model_path)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        
        self.model.eval()
        all_preds = []
        all_labels = []
        all_probs = []
        
        with torch.no_grad():
            for features, labels in tqdm(self.test_loader, desc="Testando"):
                sequence = features['sequence'].to(self.device)
                position = features['position'].to(self.device)
                alphagenome = features['alphagenome'].to(self.device)
                
                outputs = self.model(sequence, position, alphagenome)
                
                preds = (outputs > 0.5).cpu().numpy().flatten()
                probs = outputs.cpu().numpy().flatten()
                
                all_preds.extend(preds)
                all_labels.extend(labels.numpy())
                all_probs.extend(probs)
        
        # Métricas
        accuracy = accuracy_score(all_labels, all_preds)
        precision, recall, f1, _ = precision_recall_fscore_support(
            all_labels, all_preds, average='binary'
        )
        auc_roc = roc_auc_score(all_labels, all_probs)
        cm = confusion_matrix(all_labels, all_preds)
        
        # Exibir
        table = Table(title="Resultados no Teste")
        table.add_column("Métrica", style="cyan")
        table.add_column("Valor", style="green")
        
        table.add_row("Accuracy", f"{accuracy:.4f}")
        table.add_row("Precision", f"{precision:.4f}")
        table.add_row("Recall", f"{recall:.4f}")
        table.add_row("F1-Score", f"{f1:.4f}")
        table.add_row("AUC-ROC", f"{auc_roc:.4f}")
        
        console.print(table)
        
        console.print("\n[cyan]Matriz de Confusão:[/cyan]")
        console.print(f"  TN={cm[0,0]:4d}  FP={cm[0,1]:4d}")
        console.print(f"  FN={cm[1,0]:4d}  TP={cm[1,1]:4d}")
        
        # Salvar resultados
        results = {
            'accuracy': float(accuracy),
            'precision': float(precision),
            'recall': float(recall),
            'f1': float(f1),
            'auc_roc': float(auc_roc),
            'confusion_matrix': cm.tolist()
        }
        
        results_file = Path(self.config['output']['model_dir']) / 'test_results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        console.print(f"\n[green]✓ Resultados salvos: {results_file}[/green]")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Treinamento - Marcadores de Longevidade')
    parser.add_argument('--config', type=Path, required=True,
                       help='Arquivo de configuração YAML')
    parser.add_argument('--test-only', action='store_true',
                       help='Apenas testar modelo existente')
    
    args = parser.parse_args()
    
    # Carregar configuração
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Criar trainer
    trainer = LongevityTrainer(config)
    
    if args.test_only:
        trainer.test()
    else:
        trainer.train()
        trainer.test()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n[yellow]⚠ Treinamento interrompido[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[red]✗ Erro: {e}[/red]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")
        sys.exit(1)

