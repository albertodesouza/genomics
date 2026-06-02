from __future__ import annotations

import json
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import ndimage
import torch
from torch.utils.data import Dataset, DataLoader, Subset
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn

from genomics_pipeline.splitting import build_family_groups, extract_family_links

sys.path.insert(0, str(Path(__file__).parent.parent / "build_non_longevous_dataset"))
from genomic_dataset import GenomicLongevityDataset
from neural_ancestry_predictor_deprecated.config import generate_dataset_name, get_dataset_cache_dir

console = Console()

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
        self.ontology_terms = config['dataset_input'].get('ontology_terms', None)
        self.haplotype_mode = config['dataset_input']['haplotype_mode']
        self.window_center_size = config['dataset_input']['window_center_size']
        self.downsample_factor = config['dataset_input']['downsample_factor']
        self.prediction_target = config['output']['prediction_target']
        self.normalization_method = config['dataset_input'].get('normalization_method', 'zscore')
        self.derived_targets = config.get('output', {}).get('derived_targets', {})
        self._derived_target_config = self.derived_targets.get(self.prediction_target, {})
        self.valid_sample_indices = list(range(len(self.base_dataset)))
        self._valid_sample_index_set = set(self.valid_sample_indices)
        
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

    def _filter_track_indices_by_ontology(self, output_type: str, track_metadata: List[Dict]) -> List[int]:
        """Return column indices matching the requested ontology CURIEs."""
        if not self.ontology_terms:
            return list(range(len(track_metadata)))

        requested = set(self.ontology_terms)
        indices = [
            idx for idx, metadata in enumerate(track_metadata)
            if metadata.get('ontology_curie') in requested
        ]

        if not indices:
            available = sorted({metadata.get('ontology_curie', '') for metadata in track_metadata if metadata.get('ontology_curie')})
            raise ValueError(
                f"Nenhuma track encontrada para ontology_terms={sorted(requested)} em {output_type}. "
                f"Disponiveis: {available}"
            )

        return indices
    
    def _compute_normalization_params(self) -> Dict:
        """
        Computa parâmetros de normalização por track (por gene/ontologia).
        
        Cada track é normalizada independentemente para evitar que tracks com
        valores altos dominem o treinamento.
        
        Usa estatísticas online (streaming) para evitar acumular todos os dados
        em memória.  Para log/minmax apenas o máximo por track é mantido; para
        zscore usa-se o algoritmo de Welford (count, mean, M2 por track).
        Custo de memória: O(num_tracks) em vez de O(num_tracks × N × W).
        
        Returns:
            Dict com parâmetros de normalização por track
        """
        # Verificar se há valor pré-definido no config (não suportado para per-track)
        predefined_value = self.config['dataset_input'].get('normalization_value', 0)
        
        if predefined_value != 0 and predefined_value is not None:
            console.print(f"[yellow]⚠ AVISO: normalization_value não é suportado para normalização per-track[/yellow]")
            console.print(f"[yellow]  Ignorando e computando parâmetros por track...[/yellow]")
        
        num_processed = 0
        num_errors = 0
        
        total_samples = len(self.base_dataset)
        console.print(f"[cyan]Iniciando computação de normalização POR TRACK (streaming)...[/cyan]")
        console.print(f"[cyan]  • Amostras: {total_samples}[/cyan]")
        console.print(f"[cyan]  • Método: {self.normalization_method}[/cyan]")
        console.print(f"[cyan]  • Outputs: {', '.join(self.alphagenome_outputs)}[/cyan]")
        console.print(f"[cyan]  • Haplótipo: {self.haplotype_mode}[/cyan]")
        console.print(f"[cyan]  • Window center: {self.window_center_size} bases[/cyan]")
        
        # ─── Online / streaming accumulators — O(num_tracks) memory ───
        # log / minmax_keep_zero: apenas o máximo não-zero por track
        track_max   = None   # np.ndarray [num_tracks]
        # zscore: algoritmo de Welford (count, mean, M2 por track)
        track_count = None   # np.ndarray [num_tracks]
        track_mean  = None   # np.ndarray [num_tracks]
        track_M2    = None   # np.ndarray [num_tracks]
        num_tracks  = None

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task(
                "Calculando estatísticas por track (streaming)...",
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
                        # Inicializar acumuladores no primeiro sample válido
                        if num_tracks is None:
                            num_tracks = processed_array.shape[0]
                            console.print(f"[cyan]  • Total de tracks detectadas: {num_tracks}[/cyan]")
                            if self.normalization_method == 'zscore':
                                track_count = np.zeros(num_tracks, dtype=np.float64)
                                track_mean  = np.zeros(num_tracks, dtype=np.float64)
                                track_M2    = np.zeros(num_tracks, dtype=np.float64)
                            else:
                                track_max = np.zeros(num_tracks, dtype=np.float64)
                        
                        # Actualizar acumuladores por track (sem guardar dados brutos)
                        for track_idx in range(num_tracks):
                            row = processed_array[track_idx].astype(np.float64)
                            
                            if self.normalization_method == 'zscore':
                                # Welford online (parallel batch update)
                                n          = row.size
                                batch_mean = row.mean()
                                batch_var  = row.var()
                                old_count  = track_count[track_idx]
                                new_count  = old_count + n
                                delta      = batch_mean - track_mean[track_idx]
                                track_mean[track_idx]  += delta * n / new_count
                                track_M2[track_idx]    += batch_var * n + delta * delta * old_count * n / new_count
                                track_count[track_idx]  = new_count
                            else:
                                # log / minmax_keep_zero: apenas o máximo não-zero
                                nonzero = row[row > 0]
                                if len(nonzero) > 0:
                                    track_max[track_idx] = max(track_max[track_idx], float(nonzero.max()))
                        
                        num_processed += 1
                    else:
                        console.print(f"[yellow]  ⚠ Amostra {sample_id} (idx={idx}): Nenhum dado válido encontrado[/yellow]")
                        num_errors += 1
                        
                except Exception as e:
                    num_errors += 1
                    console.print(f"[yellow]  ⚠ Erro ao processar amostra {idx}: {e}[/yellow]")
                progress.update(task, advance=1)
        
        if num_tracks is None:
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
        
        console.print(f"\n[cyan]Computando parâmetros para {num_tracks} tracks...[/cyan]")
        
        # Converter acumuladores em parâmetros por track
        track_params = []
        
        for track_idx in range(num_tracks):
            if self.normalization_method == 'zscore':
                mean = float(track_mean[track_idx])
                variance = float(track_M2[track_idx] / track_count[track_idx]) if track_count[track_idx] > 0 else 0.0
                std = float(np.sqrt(variance))
                if std < 1e-8:
                    std = 1.0
                track_params.append({'mean': mean, 'std': std})
                
            elif self.normalization_method == 'minmax_keep_zero':
                xmax = float(track_max[track_idx]) if track_max[track_idx] > 0 else 1.0
                track_params.append({'max': xmax})
                
            elif self.normalization_method == 'log':
                xmax = float(track_max[track_idx]) if track_max[track_idx] > 0 else 0.0
                log_max = float(np.log1p(xmax)) if xmax > 0 else 1.0
                track_params.append({'log_max': log_max})
                
            else:
                # Fallback para zscore
                mean = float(track_mean[track_idx]) if track_mean is not None else 0.0
                variance = (float(track_M2[track_idx] / track_count[track_idx])
                            if track_count is not None and track_count[track_idx] > 0 else 0.0)
                std = float(np.sqrt(variance))
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
        
        known_classes = self.config.get('output', {}).get('known_classes')
        if known_classes is not None and len(known_classes) > 0:
            console.print(f"\n[cyan]Usando classes conhecidas do config...[/cyan]")
            sorted_targets = list(known_classes)
            self.target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
            self.idx_to_target = {idx: target for target, idx in self.target_to_idx.items()}
            self._discover_valid_sample_indices()
            console.print(f"[green]✓ Mapeamento de targets criado: {len(self.target_to_idx)} classes[/green]")
            console.print(f"[cyan]Classes: {sorted_targets}[/cyan]")
            return

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
            elif self.prediction_target in self.derived_targets:
                pedigree = dataset_metadata.get('individuals_pedigree', {})
                if pedigree:
                    classes_from_metadata = sorted({
                        target for target in
                        (self._get_target_value(pedigree_data) for pedigree_data in pedigree.values())
                        if target is not None
                    })
        
        if classes_from_metadata:
            console.print(f"\n[cyan]Carregando classes dos metadados do dataset...[/cyan]")
            sorted_targets = sorted(classes_from_metadata)
            
            self.target_to_idx = {target: idx for idx, target in enumerate(sorted_targets)}
            self.idx_to_target = {idx: target for target, idx in self.target_to_idx.items()}
            self._discover_valid_sample_indices()
            
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
        self._discover_valid_sample_indices()
        
        console.print(f"[green]✓ Mapeamento de targets criado: {len(self.target_to_idx)} classes[/green]")
        console.print(f"[cyan]Classes: {sorted_targets}[/cyan]")

    def _discover_valid_sample_indices(self):
        """Descobre quais amostras possuem target válido."""
        if self.prediction_target == 'frog_likelihood':
            self.valid_sample_indices = list(range(len(self.base_dataset)))
            self._valid_sample_index_set = set(self.valid_sample_indices)
            return

        dataset_metadata = getattr(self.base_dataset, 'dataset_metadata', None)
        if dataset_metadata is None:
            dataset_dir = Path(self.config['dataset_input']['dataset_dir'])
            metadata_file = dataset_dir / 'dataset_metadata.json'
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    dataset_metadata = json.load(f)

        individuals = []
        pedigree = {}
        if dataset_metadata:
            individuals = dataset_metadata.get('individuals', [])
            pedigree = dataset_metadata.get('individuals_pedigree', {})

        if individuals and pedigree:
            valid_indices = []
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console
            ) as progress:
                task = progress.add_task(
                    f"Filtrando amostras válidas para {self.prediction_target} (metadata)...",
                    total=len(individuals)
                )

                for idx, sample_id in enumerate(individuals):
                    sample_metadata = pedigree.get(sample_id, {})
                    target = self._get_target_value(sample_metadata)
                    if target in self.target_to_idx:
                        valid_indices.append(idx)
                    progress.update(task, advance=1)

            self.valid_sample_indices = valid_indices
            self._valid_sample_index_set = set(valid_indices)
            return

        valid_indices = []
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task(
                f"Filtrando amostras válidas para {self.prediction_target}...",
                total=len(self.base_dataset)
            )

            for idx in range(len(self.base_dataset)):
                try:
                    _, output_data = self.base_dataset[idx]
                    target = self._get_target_value(output_data)
                    if target in self.target_to_idx:
                        valid_indices.append(idx)
                except Exception:
                    pass
                progress.update(task, advance=1)

        self.valid_sample_indices = valid_indices
        self._valid_sample_index_set = set(valid_indices)
    
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
        elif self.prediction_target in self.derived_targets:
            return self._get_derived_target_value(output_data)
        elif self.prediction_target == 'frog_likelihood':
            return None  # Regressão, não há classes
        else:
            raise ValueError(f"prediction_target inválido: {self.prediction_target}")

    def _get_derived_target_value(self, output_data: Dict) -> Optional[str]:
        """Extrai target derivado a partir de outro campo categórico."""
        source_field = self._derived_target_config.get('source_field')
        class_map = self._derived_target_config.get('class_map', {})
        exclude_unmapped = self._derived_target_config.get('exclude_unmapped', False)

        if not source_field:
            raise ValueError(f"derived_targets.{self.prediction_target}.source_field não configurado")

        source_value = output_data.get(source_field)
        if source_value is None:
            return None

        for class_name, source_values in class_map.items():
            if source_value in source_values:
                return class_name

        if exclude_unmapped:
            return None

        return source_value
    
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
                h1_rows = self._process_haplotype(
                    window_data.get('predictions_h1', {}),
                    window_data.get('prediction_metadata_h1', {}),
                )
                if h1_rows is not None:
                    # h1_rows pode ser 1D (um output) ou 2D (múltiplos outputs)
                    if h1_rows.ndim == 1:
                        processed_rows.append(h1_rows.reshape(1, -1))
                    else:
                        processed_rows.append(h1_rows)
            
            if self.haplotype_mode in ['H2', 'H1+H2']:
                h2_rows = self._process_haplotype(
                    window_data.get('predictions_h2', {}),
                    window_data.get('prediction_metadata_h2', {}),
                )
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
    
    def _process_haplotype(self, predictions: Dict, prediction_metadata: Optional[Dict[str, List[Dict]]] = None) -> Optional[np.ndarray]:
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
                track_metadata = None
                if prediction_metadata is not None:
                    track_metadata = prediction_metadata.get(output_type)
                
                # Verificar se é 2D com múltiplas colunas (tracks)
                if array.ndim == 2 and array.shape[1] > 1:
                    # Array 2D: cada coluna é uma track separada
                    num_tracks = array.shape[1]
                    selected_track_indices = list(range(num_tracks))
                    if track_metadata is not None:
                        selected_track_indices = self._filter_track_indices_by_ontology(output_type, track_metadata)
                    
                    for track_idx in selected_track_indices:
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
        return len(self.valid_sample_indices)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Retorna item processado.
        
        Returns:
            Tupla (features_tensor, target_tensor)
            features_tensor tem shape [num_rows, effective_size] (2D)
        """
        if idx < 0 or idx >= len(self.valid_sample_indices):
            raise IndexError(f"Índice fora do range para dataset processado: {idx}")

        base_idx = self.valid_sample_indices[idx]
        input_data, output_data = self.base_dataset[base_idx]
        
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
    
    def get_class_names(self) -> List[str]:
        """Lista de nomes de classe na ordem do índice (0 .. num_classes-1)."""
        n = self.get_num_classes()
        return [self.idx_to_target[i] for i in range(n)]
    
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
            'family_split_mode': config['data_split'].get('family_split_mode', 'family_aware'),
            'balancing_strategy': config['data_split'].get('balancing_strategy', 'stratified'),
            'dataset_dir': config['dataset_input']['dataset_dir'],
            'prediction_target': config['output']['prediction_target'],
            'derived_target_config': config.get('output', {}).get('derived_targets', {}).get(config['output']['prediction_target']),
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
        if idx < 0 or idx >= len(processed_dataset.valid_sample_indices):
            return {"sample_id": f"sample_{idx}", "superpopulation": "UNK", "population": "UNK", "sex": 0}

        base_idx = processed_dataset.valid_sample_indices[idx]
        if base_idx >= len(individuals):
            return {"sample_id": f"sample_{base_idx}", "superpopulation": "UNK", "population": "UNK", "sex": 0}

        sample_id = individuals[base_idx]
        pedigree = individuals_pedigree.get(sample_id, {})
        _, output_data = processed_dataset.base_dataset[base_idx]
        return {
            "sample_id": sample_id,
            "superpopulation": pedigree.get("superpopulation", "UNK"),
            "population": pedigree.get("population", "UNK"),
            "sex": pedigree.get("sex", 0),
            "target": processed_dataset._get_target_value(output_data)
        }

    def get_class_label(idx: int) -> str:
        meta = get_individual_metadata(idx)
        return meta.get('target') or 'UNK'
    
    def stratified_interleave(indices: List[int], seed: int) -> Tuple[List[int], List[int], Dict[str, int]]:
        """
        Organiza índices em ordem estratificada sequencial.
        
        Retorna:
            - Lista de índices intercalados por classe
            - Lista de índices descartados (excesso de classes majoritárias)
            - Dicionário com contagem por classe
        """
        # Agrupar índices por classe do target atual
        indices_by_class = {}
        for idx in indices:
            class_name = get_class_label(idx)
            if class_name not in indices_by_class:
                indices_by_class[class_name] = []
            indices_by_class[class_name].append(idx)
        
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
            class_name = get_class_label(idx)
            class_counts[class_name] = class_counts.get(class_name, 0) + 1
        
        # Embaralhar com seed fixa para reprodutibilidade
        rng = random.Random(seed)
        shuffled = list(indices)
        rng.shuffle(shuffled)
        
        return shuffled, class_counts
    
    # Combinar todos os índices apenas para metadados/contagens.
    # Os splits já foram calculados em prepare_data() e precisam ser preservados aqui,
    # especialmente quando o modo family_aware mantém famílias no mesmo split.
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
        # Para shuffle, os splits já foram definidos anteriormente e devem ser mantidos.
        # Aqui registramos apenas a distribuição de classes para metadados do cache.
        class_counts = {}
        for idx in all_indices:
            class_name = get_class_label(idx)
            class_counts[class_name] = class_counts.get(class_name, 0) + 1
        discarded_indices = []  # Nenhum dado descartado
        new_train_indices = list(train_indices)
        new_val_indices = list(val_indices)
        new_test_indices = list(test_indices)
        
        console.print(f"\n  [bold]Randomização Simples (shuffle):[/bold]")
        console.print(f"  • Classes encontradas: {len(class_counts)}")
        for class_name, count in sorted(class_counts.items()):
            console.print(f"    - {class_name}: {count}")
        console.print(f"  • Total de amostras: {len(all_indices)} [green](nenhuma descartada)[/green]")
        console.print(f"  • Splits preservados da etapa anterior (family-aware quando aplicável)")
        
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
        train_dir = temp_cache_dir / 'train'
        train_dir.mkdir(exist_ok=True)
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Train", total=len(new_train_indices))
            for i, idx in enumerate(new_train_indices):
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
                
                # Save each sample to disk individually to save memory
                torch.save((features, target), train_dir / f"{i}.pt")
                train_metadata.append(get_individual_metadata(idx))
                progress.update(task, advance=1)
        
        # Processar dados de validação
        console.print(f"  📦 Processando dados de validação...")
        val_dir = temp_cache_dir / 'val'
        val_dir.mkdir(exist_ok=True)
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Val", total=len(new_val_indices))
            for i, idx in enumerate(new_val_indices):
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
                
                torch.save((features, target), val_dir / f"{i}.pt")
                val_metadata.append(get_individual_metadata(idx))
                progress.update(task, advance=1)
        
        # Processar dados de teste
        console.print(f"  📦 Processando dados de teste...")
        test_dir = temp_cache_dir / 'test'
        test_dir.mkdir(exist_ok=True)
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Test", total=len(new_test_indices))
            for i, idx in enumerate(new_test_indices):
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
                
                torch.save((features, target), test_dir / f"{i}.pt")
                test_metadata.append(get_individual_metadata(idx))
                progress.update(task, advance=1)
        
        # Dados foram salvos como chunks individuais para poupar memória.
        # Criamos arquivos placeholder apenas para manter compatibilidade com verificadores antigos
        console.print(f"\n  💾 Criando referências de compatibilidade (.pt)...")
        with open(temp_cache_dir / 'train_data.pt', 'w') as f: f.write("chunked_format")
        with open(temp_cache_dir / 'val_data.pt', 'w') as f: f.write("chunked_format")
        with open(temp_cache_dir / 'test_data.pt', 'w') as f: f.write("chunked_format")
        
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
                'family_split_mode': config['data_split'].get('family_split_mode', 'family_aware'),
                'balancing_strategy': balancing_strategy,
                'prediction_target': config['output']['prediction_target'],
                'derived_target_config': config.get('output', {}).get('derived_targets', {}).get(config['output']['prediction_target']),
                'input_shape': '2D',
            },
            'balancing_info': balancing_info,
            'family_split_info': config['data_split'].get('_family_split_info', {}),
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
            self._length = self._determine_length_without_loading()
            split_dir = self.data_file.parent / self.split_name
            
            if split_dir.exists() and split_dir.is_dir():
                self.data = []
                for i in range(self._length):
                    self.data.append(torch.load(split_dir / f"{i}.pt", weights_only=False))
            else:
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
                # Cache miss: carregar arquivo
                split_dir = self.data_file.parent / self.split_name
                if split_dir.exists() and split_dir.is_dir():
                    # Formato chunked (novo): carrega sample individual
                    features, target = torch.load(split_dir / f"{idx}.pt", weights_only=False)
                else:
                    # Formato legado (antigo): carrega arquivo inteiro se necessário
                    if not self._data_loaded or self.data is None:
                        self.data = torch.load(self.data_file, weights_only=False)
                        self._data_loaded = True
                    features, target = self.data[idx]
                
                # Adicionar ao cache LRU
                if self.cache_size > 0:
                    if len(self._cache) >= self.cache_size:
                        # Remover item mais antigo (LRU)
                        oldest_idx = self._cache_order.pop(0)
                        del self._cache[oldest_idx]
                    
                    # Guardar no cache LRU
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
    
    def get_class_names(self) -> List[str]:
        """Lista de nomes de classe na ordem do índice (0 .. num_classes-1)."""
        n = self.get_num_classes()
        return [self.idx_to_target[i] for i in range(n)]
    
    def get_input_shape(self) -> Tuple[int, int]:
        """
        Retorna shape da entrada.
        
        Returns:
            Tupla (num_rows, effective_size) para dados 2D
        """
        features_shape = None
        
        split_dir = self.data_file.parent / self.split_name
        if split_dir.exists() and split_dir.is_dir():
            if self._determine_length_without_loading() > 0:
                features, _ = torch.load(split_dir / "0.pt", weights_only=False)
                features_shape = features.shape
        else:
            # Format legado
            if self.loading_strategy == 'lazy' and not self._data_loaded:
                self.data = torch.load(self.data_file, weights_only=False)
                self._data_loaded = True
            if len(self.data) > 0:
                features_shape = self.data[0][0].shape

        if features_shape is not None:
            if len(features_shape) == 2:
                return tuple(features_shape)
            else:
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
        elif prediction_target in config.get('output', {}).get('derived_targets', {}):
            derived_cfg = config['output']['derived_targets'][prediction_target]
            source_field = derived_cfg.get('source_field')
            class_map = derived_cfg.get('class_map', {})
            exclude_unmapped = derived_cfg.get('exclude_unmapped', False)
            pedigree = dataset_metadata.get('individuals_pedigree', {})
            classes = set()
            for p in pedigree.values():
                source_value = p.get(source_field)
                if source_value is None:
                    continue
                matched = None
                for class_name, source_values in class_map.items():
                    if source_value in source_values:
                        matched = class_name
                        break
                if matched is not None:
                    classes.add(matched)
                elif not exclude_unmapped:
                    classes.add(source_value)
            if classes:
                classes_from_metadata = list(classes)
        
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

    model_type = config['model'].get('type', 'NN').upper()
    use_sklearn_baseline = model_type in SKLEARN_BASELINE_TYPES
    pin_memory = not use_sklearn_baseline
    if use_sklearn_baseline:
        train_workers = 0
        val_test_workers = 0
        persistent_workers = False
        console.print("  • [yellow]Sklearn baseline detectado: desabilitando pin_memory e workers paralelos[/yellow]")
    
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=train_workers,
        pin_memory=pin_memory,
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
        pin_memory=pin_memory,
        collate_fn=collate_fn,
        persistent_workers=persistent_workers if val_test_workers > 0 else False,
        worker_init_fn=worker_init_fn if val_test_workers > 0 else None
    )
    
    test_loader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=val_test_workers,
        pin_memory=pin_memory,
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

def _extract_family_links(pedigree: Dict) -> List[str]:
    """Extrai possíveis ligações familiares a partir do pedigree."""
    return extract_family_links(pedigree)


def build_family_aware_sample_groups(base_dataset: GenomicLongevityDataset, config: Dict) -> Tuple[List[List[int]], Dict[str, Any]]:
    """Constrói grupos de amostras que devem permanecer no mesmo split."""
    split_cfg = config.get('data_split', {})
    family_mode = split_cfg.get('family_split_mode', 'family_aware')

    dataset_metadata = getattr(base_dataset, 'dataset_metadata', {}) or {}
    individuals = dataset_metadata.get('individuals', [])
    pedigree_map = dataset_metadata.get('individuals_pedigree', {})
    if not individuals:
        individuals = [str(idx) for idx in range(len(base_dataset))]
    return build_family_groups(
        sample_ids=individuals,
        pedigree_map=pedigree_map,
        individuals_dir=Path(base_dataset.dataset_dir) / 'individuals',
        family_split_mode=family_mode,
    )


def build_valid_sample_index_map(base_dataset: GenomicLongevityDataset, config: Dict) -> List[int]:
    """Retorna os índices do dataset base válidos para o target configurado."""
    probe_dataset = ProcessedGenomicDataset(
        base_dataset=base_dataset,
        config=config,
        normalization_params={'mean': 0.0, 'std': 1.0},
        compute_normalization=False,
    )
    return list(probe_dataset.valid_sample_indices)


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
                shutil.copyfile(norm_source, norm_dest)
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

    valid_base_indices = build_valid_sample_index_map(base_dataset, config)
    valid_base_index_set = set(valid_base_indices)
    processed_idx_by_base_idx = {base_idx: processed_idx for processed_idx, base_idx in enumerate(valid_base_indices)}
    console.print(f"[cyan]Amostras válidas para {config['output']['prediction_target']}: {len(valid_base_indices)} / {len(base_dataset)}[/cyan]")
    
    # ==================================================================================
    # PASSO 1: FAZER SPLIT ANTES DA NORMALIZAÇÃO (prevenir data leakage)
    # ==================================================================================
    console.print("\n[bold cyan]📊 Dividindo Dataset (Train/Val/Test)[/bold cyan]")
    
    all_sample_groups, family_split_info = build_family_aware_sample_groups(base_dataset, config)
    sample_groups = []
    for group in all_sample_groups:
        filtered_group = [base_idx for base_idx in group if base_idx in valid_base_index_set]
        if filtered_group:
            sample_groups.append(filtered_group)

    total_size = sum(len(group) for group in sample_groups)
    total_groups = len(sample_groups)

    if total_groups == 0:
        raise ValueError(f"Nenhuma amostra válida encontrada para target {config['output']['prediction_target']}")

    train_group_count = int(config['data_split']['train_split'] * total_groups)
    val_group_count = int(config['data_split']['val_split'] * total_groups)
    test_group_count = total_groups - train_group_count - val_group_count

    # Criar índices de grupos
    indices = list(range(total_groups))
    random_seed = config['data_split']['random_seed']
    
    # random_seed == -1 significa NÃO embaralhar (modo debug)
    # random_seed == None significa embaralhar aleatoriamente (não reprodutível)
    # random_seed >= 0 significa embaralhar com seed (reprodutível)
    if random_seed is not None and random_seed != -1:
        np.random.seed(random_seed)
        np.random.shuffle(indices)
        console.print(f"[cyan]  • Grupos familiares embaralhados com seed {random_seed}[/cyan]")
    elif random_seed == -1:
        console.print(f"[yellow]  • MODO DEBUG: Grupos familiares NÃO embaralhados (random_seed=-1)[/yellow]")
    else:
        np.random.shuffle(indices)
        console.print(f"[yellow]  • Grupos familiares embaralhados aleatoriamente (sem seed)[/yellow]")

    train_group_indices = indices[:train_group_count]
    val_group_indices = indices[train_group_count:train_group_count + val_group_count]
    test_group_indices = indices[train_group_count + val_group_count:]

    train_base_indices = [idx for group_idx in train_group_indices for idx in sample_groups[group_idx]]
    val_base_indices = [idx for group_idx in val_group_indices for idx in sample_groups[group_idx]]
    test_base_indices = [idx for group_idx in test_group_indices for idx in sample_groups[group_idx]]

    train_indices = [processed_idx_by_base_idx[idx] for idx in train_base_indices]
    val_indices = [processed_idx_by_base_idx[idx] for idx in val_base_indices]
    test_indices = [processed_idx_by_base_idx[idx] for idx in test_base_indices]

    config['data_split']['_family_split_info'] = family_split_info
    console.print(f"[cyan]  • Modo de split familiar: {family_split_info['family_split_mode']}[/cyan]")
    console.print(f"[cyan]  • Fonte dos grupos familiares: {family_split_info.get('grouping_source', 'unknown')}[/cyan]")
    console.print(f"[cyan]  • Grupos válidos: {len(sample_groups)} (famílias={sum(1 for group in sample_groups if len(group) > 1)}, singletons={sum(1 for group in sample_groups if len(group) == 1)})[/cyan]")
    if test_group_count < 0:
        raise ValueError("Configuração inválida de train/val/test para número de grupos familiares")
    
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
        train_val_indices = train_base_indices + val_base_indices
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
            shutil.copyfile(norm_source, norm_dest)
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


# ==============================================================================
# SKLEARN BASELINES (PCA + Linear SVM / Random Forest / XGBoost)


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
