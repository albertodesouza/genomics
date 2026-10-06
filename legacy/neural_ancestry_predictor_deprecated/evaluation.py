from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_fscore_support
from torch.utils.data import DataLoader
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
from rich.table import Table
import matplotlib.pyplot as plt

from neural_ancestry_predictor_deprecated.data_pipeline import ProcessedGenomicDataset, block_reduce_2d
from neural_ancestry_predictor_deprecated.interpretability import DeepLIFT, GradCAM


console = Console()

def extract_dna_sequence(
    dataset_dir: Path,
    sample_id: str,
    gene_name: str,
    center_position: int,
    window_center_size: int,
    sequence_length: int = 1000,
    haplotype: str = 'H1'
) -> Optional[Tuple[str, str]]:
    """
    Extrai uma sequência de DNA centrada em uma posição específica.
    
    Args:
        dataset_dir: Diretório base do dataset
        sample_id: ID do indivíduo (ex: HG02635)
        gene_name: Nome do gene (ex: DDB1)
        center_position: Posição genômica central (para o header)
        window_center_size: Tamanho da janela central usada no processamento
        sequence_length: Comprimento da sequência a extrair (default 1000bp)
        haplotype: Haplótipo a usar ('H1' ou 'H2')
        
    Returns:
        Tupla (header, sequence) ou None se não encontrar
    """
    dataset_dir = Path(dataset_dir)
    
    # Caminho para o arquivo FASTA
    fasta_path = dataset_dir / 'individuals' / sample_id / 'windows' / gene_name / f'{sample_id}.{haplotype}.window.fixed.fa'
    
    if not fasta_path.exists():
        console.print(f"[yellow]⚠ Arquivo FASTA não encontrado: {fasta_path}[/yellow]")
        return None
    
    try:
        # Ler arquivo FASTA
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
        
        # Primeira linha é o header, resto é a sequência
        if len(lines) < 2:
            return None
        
        # Concatenar todas as linhas de sequência (removendo newlines)
        sequence = ''.join(line.strip() for line in lines[1:])
        
        # A sequência no FASTA é a janela original completa
        # Precisamos encontrar o centro correspondente ao window_center_size
        original_length = len(sequence)
        original_center = original_length // 2
        
        # O window_center_size é o tamanho da janela central extraída
        # O centro da janela central corresponde ao centro da sequência original
        half_window = window_center_size // 2
        
        # Extrair a sequência central de sequence_length bases
        half_seq = sequence_length // 2
        
        # A posição no centro da janela central
        # center_position é a posição genômica, mas precisamos trabalhar com índices na sequência
        # O índice na sequência é relativo ao centro
        start_idx = original_center - half_seq
        end_idx = original_center + half_seq
        
        # Garantir que estamos dentro dos limites
        start_idx = max(0, start_idx)
        end_idx = min(original_length, end_idx)
        
        extracted_seq = sequence[start_idx:end_idx]
        
        # Criar header FASTA
        header = f">{sample_id}_{haplotype}_{gene_name}_center_{center_position}"
        
        return header, extracted_seq
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro ao ler FASTA: {e}[/yellow]")
        return None




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
        self.deeplift_fasta_length = interp_config.get('deeplift', {}).get('fasta_length', 1000)
        self.top_regions_mode = interp_config.get('deeplift', {}).get('top_regions_mode', 'per_gene')
        self.top_regions_count = interp_config.get('deeplift', {}).get('top_regions_count', 5)
        self.min_distance_bp = interp_config.get('deeplift', {}).get('min_distance_bp', 50)
        self.show_top_regions_xlabel = interp_config.get('deeplift', {}).get('show_top_regions_xlabel', True)
        self.gradcam_target_class = interp_config.get('gradcam', {}).get('target_class', 'predicted')
        self.highlight_positions_bed = interp_config.get('highlight_positions_bed', None)
        self._highlight_positions_cache = None
        
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
        class_names = self.dataset.get_class_names()
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
    
    def _load_highlight_positions(
        self,
        gene_names: List[str],
        total_cols: int,
        gene_window_metadata: Dict
    ) -> List[Tuple[str, int]]:
        """
        Carrega posições de destaque de um arquivo BED e converte para
        coordenadas de plotagem (gene_name, col_idx).
        
        Retorna lista de (gene_name, col_idx) para cada posição dentro
        de uma janela gênica válida.
        """
        if self._highlight_positions_cache is not None:
            return self._highlight_positions_cache
        
        if not self.highlight_positions_bed:
            self._highlight_positions_cache = []
            return self._highlight_positions_cache
        
        bed_path = Path(self.highlight_positions_bed)
        if not bed_path.exists():
            console.print(f"[yellow]⚠ Arquivo BED de highlight não encontrado: {bed_path}[/yellow]")
            self._highlight_positions_cache = []
            return self._highlight_positions_cache
        
        positions = []
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                chrom = fields[0]
                genomic_pos = int(fields[1])
                gene_name_bed = fields[4]
                
                if gene_name_bed not in gene_names:
                    continue
                if gene_name_bed not in gene_window_metadata:
                    continue
                
                meta = gene_window_metadata[gene_name_bed]
                start_pos = meta.get('start', 0)
                window_size = meta.get('window_size', 0)
                if window_size <= 0:
                    continue
                
                col_idx = int(((genomic_pos - start_pos) / window_size) * total_cols)
                if 0 <= col_idx < total_cols:
                    positions.append((gene_name_bed, col_idx))
        
        self._highlight_positions_cache = positions
        if positions:
            console.print(f"[green]✓ {len(positions)} posições de highlight carregadas de {bed_path.name}[/green]")
        return positions
    
    def _get_genomic_coord(
        self,
        gene_name: str,
        col_idx: int,
        total_cols: int,
        gene_window_metadata: Dict
    ) -> Tuple[int, str]:
        """
        Calcula a coordenada genômica para uma posição (col_idx) dentro de um gene.
        
        Returns:
            Tuple[genomic_pos, chrom]
        """
        if gene_name in gene_window_metadata:
            meta = gene_window_metadata[gene_name]
            chrom = meta.get('chromosome', 'N/A')
            start_pos = meta.get('start', 0)
            window_size = meta.get('window_size', 0)
            genomic_pos = start_pos + int((col_idx / total_cols) * window_size) if window_size > 0 else start_pos
        else:
            chrom = 'N/A'
            genomic_pos = 0
        return genomic_pos, chrom
    
    def _compute_top_regions_per_gene(
        self,
        dl_np: np.ndarray,
        gene_names: List[str],
        tracks_per_gene: int,
        gene_window_metadata: Dict,
        top_n: int
    ) -> List[Tuple]:
        """
        Modo per_gene: encontra 1 ponto máximo por gene e ordena os genes.
        
        Cada gene pode aparecer no máximo 1 vez nos resultados.
        Este é o modo original, compatível com publicações anteriores.
        
        Returns:
            Lista de tuplas (gene_name, max_val, chrom, genomic_pos, col_idx)
        """
        top_regions = []
        total_cols = dl_np.shape[1]
        
        for i in range(len(gene_names)):
            gene_name = gene_names[i]
            start_row = i * tracks_per_gene
            end_row = (i + 1) * tracks_per_gene
            if end_row <= dl_np.shape[0]:
                gene_data = dl_np[start_row:end_row, :]
                # Encontra valor máximo positivo e sua posição
                max_val = gene_data.max()
                if max_val > 0:
                    max_idx = np.unravel_index(gene_data.argmax(), gene_data.shape)
                    col_idx = max_idx[1]
                    genomic_pos, chrom = self._get_genomic_coord(
                        gene_name, col_idx, total_cols, gene_window_metadata)
                    top_regions.append((gene_name, max_val, chrom, genomic_pos, col_idx))
        
        # Ordena por valor (maior primeiro) e retorna top N
        top_regions.sort(key=lambda x: x[1], reverse=True)
        return top_regions[:top_n]
    
    def _compute_top_regions_global(
        self,
        dl_np: np.ndarray,
        gene_names: List[str],
        tracks_per_gene: int,
        gene_window_metadata: Dict,
        top_n: int,
        min_distance_bp: int = 50
    ) -> List[Tuple]:
        """
        Modo global: considera cada posição (row, col) da matriz como candidato.
        
        Permite múltiplos pontos por gene. Útil para análise mais detalhada
        quando um gene tem várias regiões de interesse.
        
        Args:
            dl_np: DeepLIFT attribution map
            gene_names: List of gene names
            tracks_per_gene: Number of tracks per gene
            gene_window_metadata: Metadata for each gene window
            top_n: Number of top regions to return
            min_distance_bp: Minimum distance in base pairs between selected regions
        
        Returns:
            Lista de tuplas (gene_name, val, chrom, genomic_pos, col_idx, row_idx, track_idx)
        """
        candidates = []
        total_cols = dl_np.shape[1]
        
        for row_idx in range(dl_np.shape[0]):
            gene_idx = row_idx // tracks_per_gene
            if gene_idx >= len(gene_names):
                continue
            gene_name = gene_names[gene_idx]
            track_idx = row_idx % tracks_per_gene
            
            for col_idx in range(total_cols):
                val = dl_np[row_idx, col_idx]
                if val > 0:
                    genomic_pos, chrom = self._get_genomic_coord(
                        gene_name, col_idx, total_cols, gene_window_metadata)
                    # Inclui track_idx para identificação única
                    candidates.append((gene_name, val, chrom, genomic_pos, col_idx, row_idx, track_idx))
        
        # Ordena por valor (maior primeiro)
        candidates.sort(key=lambda x: x[1], reverse=True)
        
        # Filtrar candidatos muito próximos (manter apenas os de maior valor)
        selected = []
        for candidate in candidates:
            gene_name, val, chrom, genomic_pos, col_idx, row_idx, track_idx = candidate
            
            # Verificar se está muito próximo de algum ponto já selecionado
            too_close = False
            for selected_region in selected:
                sel_gene, sel_val, sel_chrom, sel_genomic_pos, sel_col_idx, sel_row_idx, sel_track_idx = selected_region
                
                # Só verificar distância se for o mesmo cromossomo
                if chrom == sel_chrom:
                    distance = abs(genomic_pos - sel_genomic_pos)
                    if distance < min_distance_bp:
                        too_close = True
                        break
            
            if not too_close:
                selected.append(candidate)
                if len(selected) >= top_n:
                    break
        
        return selected
    
    def _plot_deeplift_track_profile(
        self,
        dl_np: np.ndarray,
        gene_names: List[str],
        tracks_per_gene: int,
        window_center_size: int,
        title_suffix: str = ""
    ):
        """
        Plots an interactive line graph showing DeepLIFT values along the track
        with maximum activation for the top 5 most active genes.
        
        Args:
            dl_np: DeepLIFT attribution map (num_rows x num_cols)
            gene_names: List of gene names
            tracks_per_gene: Number of tracks per gene (typically 6)
            window_center_size: Size of the genomic window for x-axis
            title_suffix: Optional suffix for the plot title
        """
        # X-axis: positions from 0 to window_center_size
        num_cols = dl_np.shape[1]
        x_positions = np.linspace(0, window_center_size, num_cols)
        
        # First pass: calculate max value for each gene to find top 5
        gene_max_values = []
        for gene_idx, gene_name in enumerate(gene_names):
            start_row = gene_idx * tracks_per_gene
            end_row = (gene_idx + 1) * tracks_per_gene
            
            if end_row > dl_np.shape[0]:
                continue
            
            gene_data = dl_np[start_row:end_row, :]
            max_val = gene_data.max()
            max_track_idx = gene_data.max(axis=1).argmax()
            
            gene_max_values.append({
                'gene_idx': gene_idx,
                'gene_name': gene_name,
                'max_val': max_val,
                'max_track_idx': max_track_idx,
                'start_row': start_row,
                'end_row': end_row
            })
        
        # Sort by max value and take top 5
        gene_max_values.sort(key=lambda x: x['max_val'], reverse=True)
        top_5_genes = gene_max_values[:5]
        
        # Create a new figure with interactive toolbar
        fig, ax = plt.subplots(figsize=(14, 8))
        fig.canvas.manager.set_window_title('DeepLIFT Track Profile - Top 5 Genes')
        
        # Use a colormap with distinct colors for 5 genes
        cmap = plt.cm.get_cmap('tab10')
        colors = [cmap(i) for i in range(5)]
        
        # Plot one line per top gene (using the track with maximum value)
        for plot_idx, gene_info in enumerate(top_5_genes):
            gene_name = gene_info['gene_name']
            max_track_idx = gene_info['max_track_idx']
            start_row = gene_info['start_row']
            end_row = gene_info['end_row']
            
            # Get the values along the track with maximum
            gene_data = dl_np[start_row:end_row, :]
            track_values = gene_data[max_track_idx, :]
            
            # Legend shows gene name and track index (0-5)
            label = f"{gene_name} (track {max_track_idx})"
            ax.plot(x_positions, track_values, color=colors[plot_idx], 
                   linewidth=1.5, label=label, alpha=0.8)
        
        # Configure axes
        ax.set_xlabel('Genomic Position (bp)', fontsize=20)
        ax.set_ylabel('DeepLIFT Attribution', fontsize=20)
        ax.set_title(f'DeepLIFT Track Profile: Top 5 Most Active Genes{title_suffix}', 
                    fontsize=18, fontweight='bold')
        
        # Grid and legend
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
        
        # Legend outside plot
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=17, 
                 title='Gene (Track Index)', title_fontsize=19)
        
        # Adjust layout to accommodate legend
        plt.tight_layout()
        plt.subplots_adjust(right=0.82)
        
        # Enable interactive mode with toolbar (zoom, pan, home)
        # The toolbar is automatically shown by matplotlib
        
        # Show the figure (non-blocking)
        plt.show(block=False)
    
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
        
        # Get class names from dataset mapping (supports superpopulation, population, derived_targets, etc.)
        if hasattr(self.dataset, 'idx_to_target') and getattr(self.dataset, 'idx_to_target', None):
            _nc = self.dataset.get_num_classes()
            class_names = [self.dataset.idx_to_target[i] for i in range(_nc)]
        else:
            class_names = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']  # fallback legado
        
        # ─────────────────────────────────────────────────────────────────
        # Idempotency check: skip if output already exists
        # ─────────────────────────────────────────────────────────────────
        if self.interpretability_enabled and self.interp_save_images:
            output_dir = Path(self.interp_output_dir)
            if not output_dir.is_absolute():
                cache_dir = self.config.get('dataset_input', {}).get('processed_cache_dir', '.')
                output_dir = Path(cache_dir) / output_dir
            
            # Determine expected filename with interpretability parameters
            method_suffix = self.interp_method
            is_deeplift = self.interp_method == 'deeplift'
            
            # Build interpretability suffix for filename matching
            if is_deeplift:
                interp_suffix = f"{self.top_regions_mode}_top{self.top_regions_count}_{self.deeplift_fasta_length}bp_dist{self.min_distance_bp}bp_base_{self.deeplift_baseline}"
            else:
                interp_suffix = ""
            
            if is_deeplift and self.deeplift_target_class != 'predicted':
                # Class mean mode - check if class mean file exists
                target_class_name = self.deeplift_target_class
                target_class_file = str(target_class_name).replace(' ', '_')
                # We need to know the number of samples, but we can check for any file matching the pattern
                import glob
                pattern = str(output_dir / f"class_mean_{target_class_file}_*samples_{interp_suffix}_{method_suffix}.png")
                existing_files = glob.glob(pattern)
                if existing_files:
                    console.print(f"[dim]⏭ Skipping visualization (already exists): {existing_files[0]}[/dim]")
                    return
            else:
                # Individual mode
                if target_idx < len(class_names):
                    predicted_name_check = class_names[predicted_idx]
                else:
                    predicted_name_check = f"Class {predicted_idx}"
                correct_str = "correct" if predicted_idx == target_idx else "wrong"
                if interp_suffix:
                    filename = f"{sample_id}_{predicted_name_check}_{correct_str}_{interp_suffix}_{method_suffix}.png"
                else:
                    filename = f"{sample_id}_{predicted_name_check}_{correct_str}_{method_suffix}.png"
                filepath = output_dir / filename
                if filepath.exists():
                    console.print(f"[dim]⏭ Skipping visualization (already exists): {filepath}[/dim]")
                    return
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
        
        # Get gene window metadata for genomic coordinates
        gene_window_metadata = self.config['dataset_input'].get('gene_window_metadata', {})
        
        # Compute interpretability maps if enabled
        gradcam_map = None
        deeplift_map = None
        gradcam_target_class_idx = None
        gradcam_target_class_name = None
        deeplift_target_class_name = None
        deeplift_class_mean_mode = False
        deeplift_class_mean_input = None
        deeplift_class_mean_num_samples = 0
        
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
        
        # Nome de classe seguro para nomes de arquivo (evita espaços, ex.: "strong pigmentation")
        deeplift_target_class_file_slug = (
            deeplift_target_class_name.replace(' ', '_') if deeplift_target_class_name else ''
        )
        
        # Determine layout based on interpretability
        has_gradcam = gradcam_map is not None
        has_deeplift = deeplift_map is not None
        num_interp_panels = (1 if has_gradcam else 0) + (1 if has_deeplift else 0)
        
        # Calculate number of rows for subplot
        # Row 1: Input features
        # Row 2: Interpretability maps (if any)
        # Row 3 (or 2): Output probabilities (NOT shown for DeepLIFT)
        # For DeepLIFT: only 2 rows (Input + DeepLIFT), probabilities removed
        is_deeplift_mode = self.interpretability_enabled and self.interp_method == 'deeplift' and has_deeplift
        
        if is_deeplift_mode:
            num_rows_subplot = 2
            fig_height = 10  # Just Input + DeepLIFT
        elif num_interp_panels > 0:
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
        
        # Limpar dataset_name: remover conteúdo entre parênteses e adicionar "SET"
        import re
        clean_dataset_name = re.sub(r'\s*\([^)]*\)', '', self.dataset_name).strip().upper() + ' SET'
        
        # Detect if input is 2D or 1D
        if features_cpu.ndim == 3 and features_cpu.shape[0] == 1:
            # 2D input: [1, num_rows, effective_size]
            # Se estamos no modo de média de classe, usar a média das entradas
            if deeplift_class_mean_mode and deeplift_class_mean_input is not None:
                img_data = deeplift_class_mean_input.numpy()  # [num_rows, effective_size]
                input_title = f'{clean_dataset_name} | Class {deeplift_target_class_name} ({deeplift_class_mean_num_samples} samples) | Input 2D Mean ({img_data.shape[0]}x{img_data.shape[1]})'
            else:
                img_data = features_cpu[0].numpy()  # [num_rows, effective_size]
                input_title = f'{clean_dataset_name} | Sample {sample_id} ({target_name}) | Input 2D ({img_data.shape[0]}x{img_data.shape[1]})'
            
            # Reorder genes if alphabetical ordering is requested
            if gene_reorder_indices is not None:
                # Reorder rows: each gene has tracks_per_gene rows
                reordered_rows = []
                for new_idx in range(len(gene_reorder_indices)):
                    orig_idx = gene_reorder_indices[new_idx]
                    start_row = orig_idx * tracks_per_gene
                    end_row = (orig_idx + 1) * tracks_per_gene
                    reordered_rows.append(img_data[start_row:end_row, :])
                img_data = np.vstack(reordered_rows)
            
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
            
            plt.xlabel('Gene Position', fontsize=20)
            plt.title(input_title, fontsize=18, fontweight='bold')
            cbar = plt.colorbar()
            cbar.set_label('Normalized Value', fontsize=19)
            cbar.ax.tick_params(labelsize=17)
            
            # Configure X axis to show original scale (0 to window_center_size)
            window_center_size = self.config['dataset_input']['window_center_size']
            num_xticks = 5
            xtick_positions = np.linspace(0, viz_width - 1, num_xticks)
            xtick_labels = [f'{int(x)}' for x in np.linspace(0, window_center_size, num_xticks)]
            ax1.set_xticks(xtick_positions)
            ax1.set_xticklabels(xtick_labels, fontsize=17)
            
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
            plt.title(f'{clean_dataset_name} | Sample {sample_id} | Input Features (n={len(features_np)})', 
                     fontsize=18, fontweight='bold')
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
                
                # Reorder genes if alphabetical ordering is requested
                if gene_reorder_indices is not None:
                    reordered_rows = []
                    for new_idx in range(len(gene_reorder_indices)):
                        orig_idx = gene_reorder_indices[new_idx]
                        start_row = orig_idx * tracks_per_gene
                        end_row = (orig_idx + 1) * tracks_per_gene
                        reordered_rows.append(cam_np[start_row:end_row, :])
                    cam_np = np.vstack(reordered_rows)
                
                # Usar block_reduce para preservar picos ao invés de interpolação
                cam_resized = block_reduce_2d(cam_np, (viz_height, viz_width), func=downsample_agg)
                
                plt.sca(ax_gc)
                im_gc = plt.imshow(cam_resized, cmap='hot', aspect='auto', interpolation='nearest',
                                   vmin=0.0, vmax=0.0025)
                cbar_gc = plt.colorbar(im_gc)
                cbar_gc.set_label('Activation', fontsize=19)
                cbar_gc.ax.tick_params(labelsize=17)
                # Show which class is being visualized and individual's superpopulation
                gc_class_label = gradcam_target_class_name if gradcam_target_class_name else predicted_name
                plt.title(f'Grad-CAM: Important Regions for Class {gc_class_label} (Individual: {target_name})', fontsize=18, fontweight='bold')
                
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
                    plt.xlabel(f'Gene Position — Top genes: {top_genes_str}', fontsize=19)
                else:
                    plt.xlabel('Gene Position', fontsize=20)
                
                # Configure X axis to show original scale (0 to window_center_size)
                window_center_size = self.config['dataset_input']['window_center_size']
                num_xticks = 5
                xtick_positions = np.linspace(0, viz_width - 1, num_xticks)
                xtick_labels = [f'{int(x)}' for x in np.linspace(0, window_center_size, num_xticks)]
                ax_gc.set_xticks(xtick_positions)
                ax_gc.set_xticklabels(xtick_labels, fontsize=17)
                
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
                    ax_gc.set_yticklabels(gene_names, minor=True, fontsize=17)
                    ax_gc.tick_params(axis='y', which='minor', length=0)
                    
                    ax_gc.set_ylabel('Genes', fontsize=20)
            
            # Plot DeepLIFT
            if has_deeplift:
                dl_np = deeplift_map.numpy()
                
                # Reorder genes if alphabetical ordering is requested
                if gene_reorder_indices is not None:
                    reordered_rows = []
                    for new_idx in range(len(gene_reorder_indices)):
                        orig_idx = gene_reorder_indices[new_idx]
                        start_row = orig_idx * tracks_per_gene
                        end_row = (orig_idx + 1) * tracks_per_gene
                        reordered_rows.append(dl_np[start_row:end_row, :])
                    dl_np = np.vstack(reordered_rows)
                
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
                
                # Colormap personalizado: cor configurável no centro, cores vivas nos extremos
                # Azul brilhante (negativo) <- Branco/Preto (zero) -> Vermelho brilhante (positivo)
                from matplotlib.colors import LinearSegmentedColormap
                colormap_center = self.config.get('debug', {}).get('interpretability', {}).get('deeplift', {}).get('colormap_center', 'white')
                if colormap_center == 'black':
                    center_color = (0.0, 0.0, 0.0)  # Preto
                else:
                    center_color = (1.0, 1.0, 1.0)  # Branco (padrão)
                colors_diverging = [
                    (0.0, 0.5, 1.0),   # Azul brilhante (negativo máximo)
                    center_color,      # Cor central (zero)
                    (1.0, 0.2, 0.0)    # Vermelho brilhante (positivo máximo)
                ]
                cmap_diverging = LinearSegmentedColormap.from_list('diverging_center', colors_diverging)
                
                # Escala: dinâmica para modo de média de classe, fixa para outros casos
                if deeplift_class_mean_mode:
                    # Escala dinâmica baseada nos dados
                    dl_abs_max = max(abs(dl_resized.min()), abs(dl_resized.max()))
                    dl_vmin, dl_vmax = -dl_abs_max, dl_abs_max
                else:
                    # Escala fixa para comparação entre amostras
                    dl_vmin, dl_vmax = -0.10, 0.10
                
                # Aplicar transformação gamma para comprimir valores pequenos
                # gamma > 1.0: valores pequenos ficam mais próximos do centro (branco/preto)
                # gamma = 1.0: comportamento linear (sem transformação)
                colormap_gamma = self.config.get('debug', {}).get('interpretability', {}).get('deeplift', {}).get('colormap_gamma', 2.0)
                if colormap_gamma != 1.0:
                    # Normalizar para [-1, 1], aplicar gamma preservando sinal, depois desnormalizar
                    dl_normalized = dl_resized / max(abs(dl_vmin), abs(dl_vmax))
                    dl_display = np.sign(dl_normalized) * (np.abs(dl_normalized) ** (1.0 / colormap_gamma))
                    dl_display = dl_display * max(abs(dl_vmin), abs(dl_vmax))
                else:
                    dl_display = dl_resized
                
                im_dl = plt.imshow(dl_display, cmap=cmap_diverging, aspect='auto', 
                                  interpolation='nearest', vmin=dl_vmin, vmax=dl_vmax)
                cbar_dl = plt.colorbar(im_dl)
                cbar_dl.set_label('Attribution (+ → class, - → not class)', fontsize=19)
                cbar_dl.ax.tick_params(labelsize=17)
                
                # Título: diferente para modo de média de classe
                dl_class_label = deeplift_target_class_name if deeplift_target_class_name else predicted_name
                if deeplift_class_mean_mode:
                    plt.title(f'DeepLIFT: Mean Attribution for Class {dl_class_label} ({deeplift_class_mean_num_samples} samples)', fontsize=18, fontweight='bold')
                else:
                    plt.title(f'DeepLIFT: Feature Attribution for Class {dl_class_label} (Individual: {target_name})', fontsize=18, fontweight='bold')
                
                # Compute top N most active regions using configured mode
                is_global_mode = self.top_regions_mode == 'global'
                top_n = self.top_regions_count
                
                if is_global_mode:
                    # Modo global: considera cada posição da matriz como candidato
                    # Permite múltiplos pontos por gene
                    top_regions_raw = self._compute_top_regions_global(
                        dl_np, gene_names, tracks_per_gene, gene_window_metadata, top_n,
                        min_distance_bp=self.min_distance_bp)
                else:
                    # Modo per_gene: 1 ponto máximo por gene (modo original)
                    top_regions_raw = self._compute_top_regions_per_gene(
                        dl_np, gene_names, tracks_per_gene, gene_window_metadata, top_n)
                
                # Normalizar formato para (gene_name, val, chrom, genomic_pos, col_idx)
                # O modo global retorna tuplas com campos extras (row_idx, track_idx)
                if is_global_mode and top_regions_raw and len(top_regions_raw[0]) > 5:
                    top_5_regions = [(r[0], r[1], r[2], r[3], r[4]) for r in top_regions_raw]
                    top_regions_with_track = top_regions_raw  # Manter versão completa para output
                else:
                    top_5_regions = top_regions_raw
                    top_regions_with_track = None
                
                # Print details to terminal
                if top_5_regions:
                    mode_label = "global" if is_global_mode else "per_gene"
                    console.print(f"\n[bold cyan]Top {top_n} Regiões Mais Ativas (DeepLIFT, modo={mode_label}):[/bold cyan]")
                    for rank, region_tuple in enumerate(top_5_regions, 1):
                        gene, val, chrom, pos, col = region_tuple[:5]
                        # No modo global, mostrar track se disponível
                        if top_regions_with_track and rank <= len(top_regions_with_track):
                            track_idx = top_regions_with_track[rank-1][6] if len(top_regions_with_track[rank-1]) > 6 else None
                            if track_idx is not None:
                                console.print(f"  {rank}. [green]{gene}[/green] (track {track_idx}): valor = {val:.6f}, {chrom}: {pos:,}")
                            else:
                                console.print(f"  {rank}. [green]{gene}[/green]: valor = {val:.6f}, {chrom}: {pos:,}")
                        else:
                            console.print(f"  {rank}. [green]{gene}[/green]: valor = {val:.6f}, {chrom}: {pos:,}")
                    
                    # Save to file if save_images is enabled
                    if self.interp_save_images:
                        # Create output directory
                        top_regions_output_dir = Path(self.interp_output_dir)
                        if not top_regions_output_dir.is_absolute():
                            cache_dir = self.config.get('dataset_input', {}).get('processed_cache_dir', '.')
                            top_regions_output_dir = Path(cache_dir) / top_regions_output_dir
                        top_regions_output_dir.mkdir(parents=True, exist_ok=True)
                        
                        # Generate filename with all interpretability parameters
                        fasta_length = self.deeplift_fasta_length
                        interp_suffix = f"{self.top_regions_mode}_top{top_n}_{fasta_length}bp_dist{self.min_distance_bp}bp_base_{self.deeplift_baseline}"
                        if deeplift_class_mean_mode:
                            txt_filename = f"top_regions_class_mean_{deeplift_target_class_file_slug}_{deeplift_class_mean_num_samples}samples_{interp_suffix}_deeplift.txt"
                        else:
                            correct_str = "correct" if predicted_idx == target_idx else "wrong"
                            txt_filename = f"top_regions_{sample_id}_{predicted_name}_{correct_str}_{interp_suffix}_deeplift.txt"
                        
                        txt_filepath = top_regions_output_dir / txt_filename
                        
                        # Find individuals with max DeepLIFT in each region (only in class_mean mode)
                        max_individuals = {}
                        if deeplift_class_mean_mode and self.deeplift is not None:
                            max_individuals = self.deeplift.find_max_individuals_for_regions(
                                top_regions=top_5_regions,
                                target_class_idx=deeplift_target_class_idx,
                                dataset=self.dataset,
                                baseline_type=self.deeplift_baseline,
                                tracks_per_gene=tracks_per_gene,
                                gene_names=gene_names
                            )
                        
                        # Get dataset_dir and window_center_size for DNA extraction
                        dataset_dir = Path(self.config['dataset_input']['dataset_dir'])
                        window_center_size = self.config['dataset_input']['window_center_size']
                        
                        with open(txt_filepath, 'w') as f:
                            mode_desc = "modo global" if is_global_mode else "modo per_gene"
                            f.write(f"Top {top_n} Regiões Mais Ativas (DeepLIFT, {mode_desc})\n")
                            f.write(f"{'=' * 80}\n")
                            if deeplift_class_mean_mode:
                                f.write(f"Classe: {deeplift_target_class_name} ({deeplift_class_mean_num_samples} amostras)\n")
                            else:
                                f.write(f"Sample: {sample_id}\n")
                                f.write(f"Target: {target_name}\n")
                            f.write(f"{'=' * 80}\n\n")
                            
                            for rank, region_tuple in enumerate(top_5_regions, 1):
                                gene, val, chrom, pos, col = region_tuple[:5]
                                
                                # No modo global, mostrar track se disponível
                                track_info = ""
                                region_key = gene  # Chave para max_individuals
                                if top_regions_with_track and rank <= len(top_regions_with_track):
                                    raw_region = top_regions_with_track[rank-1]
                                    if len(raw_region) > 6:
                                        track_idx = raw_region[6]
                                        track_info = f" (track {track_idx})"
                                        # Chave única para modo global: gene_col_track
                                        region_key = f"{gene}_{col}_{track_idx}"
                                
                                f.write(f"{rank}. {gene}{track_info}: valor_medio = {val:.6f}, {chrom}: {pos:,}\n")
                                
                                # Add individual details if available
                                if region_key in max_individuals or gene in max_individuals:
                                    ind_info = max_individuals.get(region_key) or max_individuals.get(gene)
                                    ind_info = max_individuals[gene]
                                    ind_sample_id = ind_info['sample_id']
                                    ind_superpop = ind_info.get('superpopulation', 'UNK')
                                    ind_pop = ind_info.get('population', 'UNK')
                                    ind_max_val = ind_info['max_value']
                                    all_vals = ind_info['all_region_values']
                                    
                                    f.write(f"   Indivíduo com maior valor: {ind_sample_id} ({ind_superpop}/{ind_pop}, valor = {ind_max_val:.6f})\n")
                                    
                                    # Valores em todas as N regiões
                                    vals_str = ', '.join([f"{g}={v:.5f}" for g, v in all_vals.items()])
                                    f.write(f"   Valores nas {top_n} regiões: {vals_str}\n")
                                    
                                    # Extract DNA sequences for both haplotypes
                                    for haplotype in ['H1', 'H2']:
                                        dna_result = extract_dna_sequence(
                                            dataset_dir=dataset_dir,
                                            sample_id=ind_sample_id,
                                            gene_name=gene,
                                            center_position=pos,
                                            window_center_size=window_center_size,
                                            sequence_length=fasta_length,
                                            haplotype=haplotype
                                        )
                                        
                                        if dna_result:
                                            header, sequence = dna_result
                                            f.write(f"   DNA {haplotype} ({fasta_length}bp centradas em {chrom}:{pos:,}):\n")
                                            f.write(f"   {header}\n")
                                            # Write sequence in lines of 60 chars
                                            for i in range(0, len(sequence), 60):
                                                f.write(f"   {sequence[i:i+60]}\n")
                                        else:
                                            f.write(f"   DNA {haplotype}: não disponível\n")
                                
                                f.write("\n")
                        
                        console.print(f"[dim]Saved top regions: {txt_filepath}[/dim]")
                    
                    # Draw green circles centered on each top region
                    # Use Ellipse to compensate for aspect ratio distortion
                    from matplotlib.patches import Ellipse
                    total_cols = dl_np.shape[1]  # Needed for position calculation
                    num_genes = len(gene_names)
                    height_per_gene = viz_height / num_genes
                    
                    # Calculate aspect ratio to make circles appear circular
                    # With aspect='auto', the image is stretched to fit the axes
                    # If viz_width > viz_height, circles appear elongated horizontally
                    # To compensate: make ellipse width larger in data units
                    aspect_ratio = (viz_width / viz_height) / 3.5
                    circle_height = height_per_gene  # Diameter in Y direction = height per gene
                    circle_width = circle_height * aspect_ratio  # Multiply to compensate for horizontal stretch
                    
                    for gene_name_region, max_val, chrom, genomic_pos, col_idx in top_5_regions:
                        # Find gene index
                        if gene_name_region in gene_names:
                            gene_idx = gene_names.index(gene_name_region)
                            # X position: col_idx scaled to viz_width
                            x_pos = (col_idx / total_cols) * viz_width
                            # Y position: center of the gene
                            y_pos = (gene_idx + 0.5) * height_per_gene
                            
                            # Create hollow green ellipse that appears as a circle
                            ellipse_dl = Ellipse((x_pos, y_pos), circle_width, circle_height, 
                                            fill=False, edgecolor='lime', linewidth=2.5)
                            ax_dl.add_patch(ellipse_dl)

                            # Add matching circle on the input plot
                            ellipse_in = Ellipse((x_pos, y_pos), circle_width, circle_height, 
                                            fill=False, edgecolor='lime', linewidth=2.5)
                            ax1.add_patch(ellipse_in)
                
                # Draw orange circles at BED highlight positions
                if self.highlight_positions_bed:
                    highlight_positions = self._load_highlight_positions(
                        gene_names, total_cols, gene_window_metadata)
                    
                    if highlight_positions:
                        from matplotlib.patches import Ellipse as EllipseHL
                        hl_num_genes = len(gene_names)
                        hl_height_per_gene = viz_height / hl_num_genes
                        hl_aspect_ratio = (viz_width / viz_height) / 3.5
                        hl_circle_height = hl_height_per_gene
                        hl_circle_width = hl_circle_height * hl_aspect_ratio
                        
                        for hl_gene_name, hl_col_idx in highlight_positions:
                            if hl_gene_name in gene_names:
                                hl_gene_idx = gene_names.index(hl_gene_name)
                                hl_x_pos = (hl_col_idx / total_cols) * viz_width
                                hl_y_pos = (hl_gene_idx + 0.5) * hl_height_per_gene
                                
                                hl_ellipse_dl = EllipseHL(
                                    (hl_x_pos, hl_y_pos), hl_circle_width, hl_circle_height,
                                    fill=False, edgecolor='orange', linewidth=2.5)
                                ax_dl.add_patch(hl_ellipse_dl)
                                
                                hl_ellipse_in = EllipseHL(
                                    (hl_x_pos, hl_y_pos), hl_circle_width, hl_circle_height,
                                    fill=False, edgecolor='orange', linewidth=2.5)
                                ax1.add_patch(hl_ellipse_in)
                
                # Optional: include top regions in X-axis label
                if self.show_top_regions_xlabel and top_5_regions:
                    top_regions_labels = []
                    for idx, region_tuple in enumerate(top_5_regions):
                        gene = region_tuple[0]
                        val = region_tuple[1]
                        label = f"{gene}({val:.5f})"
                        if top_regions_with_track and idx < len(top_regions_with_track):
                            raw_region = top_regions_with_track[idx]
                            if len(raw_region) > 6:
                                track_idx = raw_region[6]
                                label = f"{gene}(t{track_idx},{val:.5f})"
                        top_regions_labels.append(label)
                    top_regions_str = ', '.join(top_regions_labels)
                    plt.xlabel(f"Gene Position\nTop regions: {top_regions_str}", fontsize=19)
                else:
                    plt.xlabel('Gene Position', fontsize=20)
                
                # Configure X axis to show original scale (0 to window_center_size)
                window_center_size = self.config['dataset_input']['window_center_size']
                num_xticks = 5
                xtick_positions = np.linspace(0, viz_width - 1, num_xticks)
                xtick_labels = [f'{int(x)}' for x in np.linspace(0, window_center_size, num_xticks)]
                ax_dl.set_xticks(xtick_positions)
                ax_dl.set_xticklabels(xtick_labels, fontsize=17)
                
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
                    ax_dl.set_yticklabels(gene_names, minor=True, fontsize=17)
                    ax_dl.tick_params(axis='y', which='minor', length=0)
                    
                    ax_dl.set_ylabel('Genes', fontsize=20)
        
        # ─────────────────────────────────────────────────────────────────
        # Last Row: Output probabilities (NOT shown for DeepLIFT mode)
        # ─────────────────────────────────────────────────────────────────
        if not is_deeplift_mode:
            ax_prob = plt.subplot(num_rows_subplot, 1, num_rows_subplot)
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
            
            # Generate filename with interpretability parameters
            method_suffix = self.interp_method
            
            # Build interpretability suffix with key parameters
            interp_params = []
            if is_deeplift_mode:
                interp_params.append(f"{self.top_regions_mode}")
                interp_params.append(f"top{self.top_regions_count}")
                interp_params.append(f"{self.deeplift_fasta_length}bp")
                interp_params.append(f"dist{self.min_distance_bp}bp")
                interp_params.append(f"base_{self.deeplift_baseline}")
            interp_suffix = "_".join(interp_params) if interp_params else ""
            
            # Nome diferente para modo de média de classe
            if deeplift_class_mean_mode:
                if interp_suffix:
                    filename = f"class_mean_{deeplift_target_class_file_slug}_{deeplift_class_mean_num_samples}samples_{interp_suffix}_{method_suffix}.png"
                else:
                    filename = f"class_mean_{deeplift_target_class_file_slug}_{deeplift_class_mean_num_samples}samples_{method_suffix}.png"
            else:
                correct_str = "correct" if predicted_idx == target_idx else "wrong"
                if interp_suffix:
                    filename = f"{sample_id}_{predicted_name}_{correct_str}_{interp_suffix}_{method_suffix}.png"
                else:
                    filename = f"{sample_id}_{predicted_name}_{correct_str}_{method_suffix}.png"
            
            filepath = output_dir / filename
            
            plt.savefig(filepath, dpi=150, bbox_inches='tight')
            console.print(f"[dim]Saved: {filepath}[/dim]")
        
        # ─────────────────────────────────────────────────────────────────
        # DeepLIFT Track Profile: separate interactive line graph
        # ─────────────────────────────────────────────────────────────────
        if is_deeplift_mode and deeplift_map is not None:
            # Get window_center_size for x-axis
            window_center_size = self.config['dataset_input']['window_center_size']
            
            # Determine title suffix based on mode
            if deeplift_class_mean_mode:
                title_suffix = f" | Class {deeplift_target_class_name} ({deeplift_class_mean_num_samples} samples)"
            else:
                title_suffix = f" | {sample_id} ({target_name})"
            
            # Plot the track profile (non-blocking, separate window)
            self._plot_deeplift_track_profile(
                dl_np=deeplift_map if isinstance(deeplift_map, np.ndarray) else deeplift_map.cpu().numpy(),
                gene_names=gene_names,
                tracks_per_gene=tracks_per_gene,
                window_center_size=window_center_size,
                title_suffix=title_suffix
            )
        
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


