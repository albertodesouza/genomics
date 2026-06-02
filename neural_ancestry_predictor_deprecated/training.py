from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import psutil
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
import matplotlib.pyplot as plt

from neural_ancestry_predictor_deprecated.interrupts import interrupt_state


console = Console()

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
            'epoch': [],
            'train_data_time': [],
            'train_compute_time': [],
            'train_data_fraction': [],
            'val_data_time': [],
            'val_compute_time': [],
            'val_data_fraction': [],
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
        class_names = self.dataset.get_class_names()
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
            gene_names_original = [gene for gene in GENE_ORDER if gene in genes_to_use]
        else:
            gene_names_original = GENE_ORDER
        tracks_per_gene = 6
        
        # Check if alphabetical ordering is requested for visualization
        viz_gene_order = self.config.get('debug', {}).get('visualization', {}).get('gene_order', 'dataset')
        if viz_gene_order == 'alphabetical':
            gene_names = sorted(gene_names_original)
            # Create reordering indices: map from new order to original order
            gene_reorder_indices = [gene_names_original.index(g) for g in gene_names]
        else:
            gene_names = gene_names_original
            gene_reorder_indices = None
        
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
            
            # Reorder genes if alphabetical ordering is requested
            if gene_reorder_indices is not None:
                reordered_rows = []
                for new_idx in range(len(gene_reorder_indices)):
                    orig_idx = gene_reorder_indices[new_idx]
                    start_row = orig_idx * tracks_per_gene
                    end_row = (orig_idx + 1) * tracks_per_gene
                    reordered_rows.append(img_data[start_row:end_row, :])
                img_data = np.vstack(reordered_rows)
            
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
            plt.xlabel('Genomic Position (rescaled)', fontsize=20)
            plt.title(f'Epoch {epoch + 1} | Sample {sample_id} ({target_name}) | Input 2D ({img_data.shape[0]}x{img_data.shape[1]} → {viz_height}x{viz_width})', 
                     fontsize=18, fontweight='bold')
            cbar = plt.colorbar()
            cbar.set_label('Normalized Value', fontsize=19)
            cbar.ax.tick_params(labelsize=17)
            
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
                ax1.set_yticklabels(gene_names, minor=True, fontsize=17)
                ax1.tick_params(axis='y', which='minor', length=0)  # Hide minor tick marks
                
                ax1.set_ylabel('Genes', fontsize=20)
            else:
                ax1.set_ylabel('Tracks (rescaled)', fontsize=20)
        else:
            # 1D input (fallback for backwards compatibility)
            features_np = features_cpu.numpy().flatten()
            plt.plot(features_np, linewidth=0.5, alpha=0.7)
            plt.xlabel('Feature Index', fontsize=20)
            plt.ylabel('Feature Value', fontsize=20)
            plt.title(f'Epoch {epoch + 1} | Sample {sample_id} | Input Features (n={len(features_np)})', 
                     fontsize=18, fontweight='bold')
            plt.grid(True, alpha=0.3)
            
        # Plot output probabilities
        plt.subplot(2, 1, 2)
        bars = plt.bar(range(len(output_probs)), output_probs, color='steelblue', alpha=0.7)
        bars[target_idx].set_color('green')
        bars[predicted_idx].set_edgecolor('red')
        bars[predicted_idx].set_linewidth(3)
        
        plt.xlabel('Class', fontsize=20)
        plt.ylabel('Probability', fontsize=20)
        plt.title('Network Output Probabilities', fontsize=18, fontweight='bold')
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
                fontsize=19, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
        
        plt.tight_layout()
        
        # Connect key callback
        self._key_pressed = False
        cid = fig.canvas.mpl_connect('key_press_event', self._on_key_press)
        
        # Show and wait for key
        plt.show(block=True)
        
        # Disconnect callback
        fig.canvas.mpl_disconnect(cid)
    
    def train_epoch(self, epoch: int) -> Tuple[float, Dict[str, float]]:
        """
        Treina por uma época.
        
        Returns:
            Loss média da época
        """
        self.model.train()
        total_loss = 0.0
        total_samples = 0
        total_data_time = 0.0
        total_compute_time = 0.0
        
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

            data_wait_start = time.perf_counter()
            for batch_idx, (features, targets, indices) in enumerate(self.train_loader):
                data_ready_time = time.perf_counter()
                total_data_time += data_ready_time - data_wait_start

                # Verificar limite de amostras para visualização
                if self.enable_visualization and self.max_samples_per_epoch is not None:
                    if batch_idx >= self.max_samples_per_epoch:
                        break

                compute_start = time.perf_counter()
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

                total_compute_time += time.perf_counter() - compute_start
                data_wait_start = time.perf_counter()
                progress.update(task, advance=1)
        
        avg_loss = total_loss / total_samples if total_samples > 0 else 0.0
        total_train_time = total_data_time + total_compute_time
        timing_metrics = {
            'data_time': total_data_time,
            'compute_time': total_compute_time,
            'data_fraction': (total_data_time / total_train_time) if total_train_time > 0 else 0.0,
        }
        
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
        
        return avg_loss, timing_metrics
    
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
    
    def validate(self, epoch: int) -> Tuple[float, float, Dict[str, float]]:
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
        total_data_time = 0.0
        total_compute_time = 0.0
        
        # Debug: salvar detalhes na época 50
        debug_epoch_50 = (epoch + 1 == 50)
        debug_data = []
        
        with torch.no_grad():
            data_wait_start = time.perf_counter()
            for features, targets, indices in self.val_loader:
                data_ready_time = time.perf_counter()
                total_data_time += data_ready_time - data_wait_start
                compute_start = time.perf_counter()
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
                total_compute_time += time.perf_counter() - compute_start
                data_wait_start = time.perf_counter()
        
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
            accuracy = float((np.asarray(all_targets) == np.asarray(all_predictions)).mean())
        else:
            accuracy = 0.0

        total_val_time = total_data_time + total_compute_time
        timing_metrics = {
            'data_time': total_data_time,
            'compute_time': total_compute_time,
            'data_fraction': (total_data_time / total_val_time) if total_val_time > 0 else 0.0,
        }

        return avg_loss, accuracy, timing_metrics
    
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
                train_loss, train_timing = self.train_epoch(epoch)
                train_duration = time.time() - train_start
                rss_mb = process.memory_info().rss / 1024 / 1024
                console.print(f"[dim cyan]Tempo treino: {train_duration:.2f}s | RSS: {rss_mb:.1f} MB[/dim cyan]")
                console.print(
                    f"[dim cyan]  Train timing: data={train_timing['data_time']:.2f}s | "
                    f"compute={train_timing['compute_time']:.2f}s | "
                    f"data%={train_timing['data_fraction'] * 100:.1f}%[/dim cyan]"
                )
                
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
                    val_loss, val_accuracy, val_timing = self.validate(epoch)
                    val_duration = time.time() - val_start
                    rss_mb = process.memory_info().rss / 1024 / 1024
                    console.print(f"[dim cyan]Tempo validação: {val_duration:.2f}s | RSS: {rss_mb:.1f} MB[/dim cyan]")
                    console.print(
                        f"[dim cyan]  Val timing: data={val_timing['data_time']:.2f}s | "
                        f"compute={val_timing['compute_time']:.2f}s | "
                        f"data%={val_timing['data_fraction'] * 100:.1f}%[/dim cyan]"
                    )
                    
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
                    self.history['train_data_time'].append(train_timing['data_time'])
                    self.history['train_compute_time'].append(train_timing['compute_time'])
                    self.history['train_data_fraction'].append(train_timing['data_fraction'])
                    self.history['val_data_time'].append(val_timing['data_time'])
                    self.history['val_compute_time'].append(val_timing['compute_time'])
                    self.history['val_data_fraction'].append(val_timing['data_fraction'])
                    self.history['epoch'].append(epoch + 1)
                    
                    # Log no W&B
                    if self.wandb_run:
                        log_dict = {
                            'epoch': epoch + 1,
                            'train_loss': train_loss,
                            'train_accuracy': train_accuracy,
                            'val_loss': val_loss,
                            'val_accuracy': val_accuracy,
                            'learning_rate': current_lr,
                            'train/data_time': train_timing['data_time'],
                            'train/compute_time': train_timing['compute_time'],
                            'train/data_fraction': train_timing['data_fraction'],
                            'val/data_time': val_timing['data_time'],
                            'val/compute_time': val_timing['compute_time'],
                            'val/data_fraction': val_timing['data_fraction'],
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

