#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_processed_dataset.py

Programa para verificar e comparar dados processados em cache com dados originais do AlphaGenome.
Visualiza tracks (genes × ontologias RNA-seq) para detectar bugs no pipeline.

Uso:
    python3 verify_processed_dataset.py --config configs/verify_processed_dataset.yaml

Author: ChatGPT (for Alberto)
Created: 2025-11-23
Updated: 2025-11-23 - Adicionado filtro por gene, navegação interativa, config YAML
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import warnings

import numpy as np
import torch
import yaml
import matplotlib.pyplot as plt
from matplotlib.backend_bases import KeyEvent
from scipy.stats import pearsonr
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

# Importar funções de normalização do neural_ancestry_predictor
sys.path.insert(0, str(Path(__file__).parent))
from neural_ancestry_predictor import log_normalize, minmax_keep_zero

# Try to import AlphaGenome API client
try:
    from alphagenome.models import dna_client
    ALPHAGENOME_AVAILABLE = True
except ImportError:
    ALPHAGENOME_AVAILABLE = False
    dna_client = None

console = Console()

NUM_ONTOLOGIES = 6  # Número de ontologias RNA-seq por gene


def get_genes_in_dataset_order(dataset_dir: Path, sample_id: str) -> List[str]:
    """
    Obtém a lista de genes na mesma ordem que o GenomicLongevityDataset processa.
    
    IMPORTANTE: A ordem dos genes no cache depende da ordem em que o GenomicLongevityDataset
    itera sobre o dicionário 'windows', que vem do individual_metadata.json.
    Esta ordem NÃO é alfabética!
    
    Args:
        dataset_dir: Diretório do dataset
        sample_id: ID do indivíduo
        
    Returns:
        Lista de genes na ordem correta (mesma do cache)
    """
    # Ler ordem do individual_metadata.json (mesma que GenomicLongevityDataset usa)
    metadata_file = dataset_dir / 'individuals' / sample_id / 'individual_metadata.json'
    
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata não encontrado: {metadata_file}")
    
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    
    genes = metadata.get('windows', [])
    
    if not genes:
        raise ValueError(f"Nenhuma janela encontrada no metadata de {sample_id}")
    
    return genes


def load_haplotype_sequence(
    dataset_dir: Path,
    sample_id: str,
    gene_name: str,
    haplotype: str = "H1"
) -> Optional[str]:
    """
    Load personalized haplotype sequence from dataset_dir.
    
    Args:
        dataset_dir: Dataset directory
        sample_id: Sample ID
        gene_name: Gene name
        haplotype: "H1" or "H2"
    
    Returns:
        DNA sequence string or None
    """
    seq_file = dataset_dir / "individuals" / sample_id / "windows" / gene_name / f"{sample_id}.{haplotype}.window.fixed.fa"
    
    if not seq_file.exists():
        return None
    
    # Read FASTA file
    sequence = ""
    with open(seq_file) as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    
    return sequence


def _extract_reference_sequence_from_interval(
    interval,
    gene_name: str,
) -> Optional[str]:
    """
    Extrai a sequência de referência do intervalo usando samtools faidx.
    
    Args:
        interval: AlphaGenome Interval object
        gene_name: Nome do gene
        fasta_path: Caminho para o arquivo FASTA de referência
    
    Returns:
        Sequência de referência como string, ou None se falhar
    """
    fasta_path = Path("/dados/GENOMICS_DATA/top3/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    try:
        import subprocess
        
        if not fasta_path.exists():
            console.print(f"[yellow]  Reference FASTA not found at {fasta_path}, skipping FASTA export[/yellow]")
            return None
        
        # Alguns objetos Interval usam .chrom, outros .chromosome; tentamos ambos
        chrom = getattr(interval, "chrom", None)
        if chrom is None:
            chrom = getattr(interval, "chromosome", None)
        if chrom is None:
            raise AttributeError("Interval object has no 'chrom' or 'chromosome' attribute")
        
        # AlphaGenome usa coordenadas 0-based [start, end),
        # samtools faidx usa 1-based inclusivo: start+1 .. end
        start_1based = interval.start + 1
        end_1based = interval.end
        region = f"{chrom}:{start_1based}-{end_1based}"
        
        console.print(f"[cyan]  Extracting reference interval with samtools faidx: {region}[/cyan]")
        
        interval_fasta_result = subprocess.run(
            ["samtools", "faidx", str(fasta_path), region],
            capture_output=True,
            text=True
        )
        
        if interval_fasta_result.returncode != 0:
            console.print(f"[red]  samtools faidx failed: {interval_fasta_result.stderr}[/red]")
            return None
        
        # Salvar FASTA em arquivo
        interval_fasta_output_path = Path("reference_interval.fasta")
        with open(interval_fasta_output_path, "w") as f:
            f.write(interval_fasta_result.stdout)
        console.print(f"[green]  ✓ Saved interval FASTA to: {interval_fasta_output_path}[/green]")
        
        # Converter FASTA para sequência contínua
        fasta_text = interval_fasta_result.stdout
        seq_lines = []
        for line in fasta_text.split("\n"):
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line.strip())
        
        seq = "".join(seq_lines)
               
        if seq is None:
            console.print("[red]  Cannot call predict_sequence: reference extraction failed[/red]")
            sys.exit(1)

        return seq
        
    except Exception as e:
        console.print(f"[red]  Error extracting reference sequence: {e}[/red]")
        import traceback
        traceback.print_exc()
        sys.exit(1)


# def _load_individual_haplotype_sequence(
#     sample_id: str,
#     gene_name: str,
#     target_length: int,
#     config: Dict
# ) -> Optional[str]:
#     """
#     Carrega a sequência do haplótipo H1 do indivíduo de dataset_dir.
    
#     Args:
#         sample_id: ID do indivíduo
#         gene_name: Nome do gene
#         dataset_dir: Diretório raiz do dataset
#         target_length: Tamanho da janela a extrair do centro
    
#     Returns:
#         Sequência do haplótipo como string, ou None se falhar
#     """
#     dataset_dir = Path(config.get('dataset_dir', '/dados/GENOMICS_DATA/top3/non_longevous_results_genes'))
#     try:
#         console.print(f"[cyan]  Loading individual's haplotype sequence from dataset_dir...[/cyan]")
        
#         # Caminho para a sequência do haplótipo H1
#         # Padrão: dataset_dir/individuals/SAMPLE_ID/windows/GENE_NAME/SAMPLE_ID.H1.window.fixed.fa
#         haplotype_file = dataset_dir / "individuals" / sample_id / "windows" / gene_name / f"{sample_id}.H1.window.fixed.fa"
        
#         if not haplotype_file.exists():
#             console.print(f"[yellow]  Haplotype file not found: {haplotype_file}[/yellow]")
#             return None
        
#         # Ler o arquivo FASTA do haplótipo
#         console.print(f"[cyan]  Reading: {haplotype_file}[/cyan]")
        
#         with open(haplotype_file, 'r') as f:
#             haplotype_lines = []
#             for line in f:
#                 if not line or line.startswith(">"):
#                     continue
#                 haplotype_lines.append(line.strip())
        
#         haplotype_seq = "".join(haplotype_lines)
        
#         # Extrair a janela central correspondente ao target_length
#         seq_length = len(haplotype_seq)
        
#         console.print(f"[dim]  Haplotype full length: {seq_length} bp[/dim]")
#         console.print(f"[dim]  Target length: {target_length} bp[/dim]")
        
#         # Extrair centro
#         center_idx = seq_length // 2
#         half_size = target_length // 2
#         start_idx = max(0, center_idx - half_size)
#         end_idx = min(seq_length, start_idx + target_length)
        
#         individual_seq = haplotype_seq[start_idx:end_idx]
        
#         console.print(f"[green]  ✓ Loaded individual sequence: {len(individual_seq)} bp (from position {start_idx} to {end_idx})[/green]")
#         console.print(f"[green]  ✓ Using {sample_id}'s H1 haplotype sequence[/green]")
        
#         return individual_seq
        
#     except Exception as e:
#         console.print(f"[yellow]  Warning: Could not load individual sequence: {e}[/yellow]")
#         import traceback
#         traceback.print_exc()
#         return None


def predict_with_alphagenome(
    gene_name: str,
    config: Dict,
    ontology_terms: List[str],
    window_size_key: str = "SEQUENCE_LENGTH_16KB",
    return_full_output: bool = False,
    sample_id: Optional[str] = None
) -> Optional[np.ndarray]:
    """
    Get predictions from AlphaGenome API usando predict_interval com genoma de referência (como no Colab).
    
    Args:
        gene_name: Gene symbol (e.g., 'MC1R')
        config: Configuration dict
        ontology_terms: List of ontology term CURIEs
        window_size_key: Tamanho da janela (e.g., 'SEQUENCE_LENGTH_16KB')
        return_full_output: If True, return (output_object, gtf) tuple instead of just values array
        sample_id: Optional sample ID to load individual's haplotype sequence instead of reference
    
    Returns:
        If return_full_output=False: Array shape [sequence_length, num_ontologies] with RNA-seq predictions or None
        If return_full_output=True: Tuple (output_object, gtf) for Colab-style plotting
    """
    if not ALPHAGENOME_AVAILABLE:
        raise ImportError("alphagenome package not available. Install with: pip install alphagenome")
    
    try:
        import pandas as pd
        from alphagenome.data import gene_annotation, genome
        import subprocess
        from pathlib import Path
    except ImportError as e:
        console.print(f"[red]Error importing required modules: {e}[/red]")
        return None
    
    api_key = config['alphagenome_api'].get('api_key') or os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError("AlphaGenome API key not provided. Set alphagenome_api.api_key in config or ALPHAGENOME_API_KEY environment variable")
    
    # Initialize client
    client = dna_client.create(api_key)
    
    console.print(f"[cyan]  Loading GTF for gene {gene_name}...[/cyan]")
    
    try:
        # Load GTF from local cache
        gtf_cache_path = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes/gtf_cache.feather")
        if not gtf_cache_path.exists():
            console.print(f"[yellow]  GTF cache not found, downloading from internet...[/yellow]")
            gtf = pd.read_feather(
                'https://storage.googleapis.com/alphagenome/reference/gencode/'
                'hg38/gencode.v46.annotation.gtf.gz.feather'
            )
        else:
            console.print(f"[green]  Loading GTF from cache: {gtf_cache_path}[/green]")
            gtf = pd.read_feather(gtf_cache_path)
        
        # Get gene interval
        interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_name)
        console.print(f"[green]  ✓ Found gene interval: {interval}[/green]")
               
        # Resize to specified window size (como no Colab: interval.resize(dna_client.SEQUENCE_LENGTH_16KB))
        window_size_enum = getattr(dna_client, window_size_key, None)
        if window_size_enum is None:
            console.print(f"[red]  Invalid window_size_key: {window_size_key}[/red]")
            return None
        
        interval = interval.resize(window_size_enum)
        console.print(f"[green]  ✓ Resized to {window_size_key}: {interval}[/green]")
 
        requested_outputs = [dna_client.OutputType.RNA_SEQ]
        console.print(f"[cyan]  Calling AlphaGenome API...[/cyan]")
        console.print(f"[dim]  Ontology terms: {ontology_terms}[/dim]")
        
        # 1. Extrair sequência de referência do intervalo
        seq = _extract_reference_sequence_from_interval(interval, gene_name)
         
        # 2. Se sample_id fornecido, sobrescrever com sequência do indivíduo
        # seq = _load_individual_haplotype_sequence(sample_id, gene_name, interval.end - interval.start, config)
     
        console.print(f"[cyan]  Calling AlphaGenome API (predict_sequence) on FASTA interval...[/cyan]")
        start_time = time.time()
        output = client.predict_sequence(
            seq,
            interval=interval,
            requested_outputs=requested_outputs,
            ontology_terms=ontology_terms
        )
        elapsed = time.time() - start_time

        # fim do NOVO TRECHO: usa o FASTA extraído (interval_fasta_result)
        # --------------------------------------------------------------
        
        # --------------------------------------------------------------
        # TRECHO ANTIGO: usa o interval (interval)
        # --------------------------------------------------------------
        # console.print(f"[cyan]  Calling AlphaGenome API (predict_interval) on interval...[/cyan]")
        # start_time = time.time()
        # output = client.predict_interval(
        #     interval=interval,
        #     requested_outputs=requested_outputs,
        #     ontology_terms=ontology_terms
        # )
        # elapsed = time.time() - start_time
        
        # Fim do trecho antigo: usa o interval (interval)
        # --------------------------------------------------------------

        console.print(f"[cyan]  API call completed in {elapsed:.1f}s[/cyan]")
        
        # Extract RNA-seq data
        rna_data = output.rna_seq
        if rna_data is None:
            console.print(f"[red]  No RNA-seq data returned from API[/red]")
            return None
        
        # Convert to numpy array [sequence_length, num_ontologies]
        values = np.array(rna_data.values)
        
        console.print(f"[green]  ✓ Received predictions shape: {values.shape}[/green]")
        
        # Apply rate limiting
        delay = config['alphagenome_api'].get('rate_limit_delay', 0.5)
        if delay > 0:
            time.sleep(delay)
        
        # Return full output or just values
        if return_full_output:
            return (output, gtf)
        else:
            return values
        
    except Exception as e:
        console.print(f"[red]  API error: {e}[/red]")
        import traceback
        traceback.print_exc()
        return None


# ==============================================================================
# CLASSE DE NAVEGAÇÃO INTERATIVA
# ==============================================================================

class InteractiveViewer:
    """Gerencia navegação interativa entre amostras."""
    
    def __init__(self, config: Dict, splits_data: Dict, dataset_metadata: Dict):
        """
        Inicializa o viewer interativo.
        
        Args:
            config: Configuração do programa
            splits_data: Dados dos splits (train/val/test indices)
            dataset_metadata: Metadata do dataset (lista de indivíduos)
        """
        self.config = config
        self.splits_data = splits_data
        self.dataset_metadata = dataset_metadata
        self.current_index = config.get('index', 0)
        self.current_split = config.get('split', 'test')
        self.should_exit = False
        self.fig = None
        
        # Obter lista de índices do split atual
        split_key = f"{self.current_split}_indices"
        self.split_indices = splits_data[split_key]
        self.max_index = len(self.split_indices) - 1
        
        console.print(f"[cyan]Modo interativo ativado:[/cyan]")
        console.print(f"  • Split: {self.current_split}")
        console.print(f"  • Total de amostras: {len(self.split_indices)}")
        console.print(f"  • Índice inicial: {self.current_index}")
        console.print(f"[yellow]  • Use ← → para navegar, 'q' para sair[/yellow]\n")
    
    def on_key_press(self, event: KeyEvent):
        """Handler para eventos de teclado."""
        if event.key == 'right':
            # Próxima amostra
            if self.current_index < self.max_index:
                self.current_index += 1
                console.print(f"[cyan]→ Avançando para índice {self.current_index}[/cyan]")
                plt.close(self.fig)
            else:
                console.print(f"[yellow]⚠ Já está na última amostra do split ({self.max_index})[/yellow]")
        
        elif event.key == 'left':
            # Amostra anterior
            if self.current_index > 0:
                self.current_index -= 1
                console.print(f"[cyan]← Retrocedendo para índice {self.current_index}[/cyan]")
                plt.close(self.fig)
            else:
                console.print(f"[yellow]⚠ Já está na primeira amostra do split (0)[/yellow]")
        
        elif event.key == 'q':
            # Sair
            console.print(f"[yellow]Saindo...[/yellow]")
            self.should_exit = True
            plt.close(self.fig)
    
    def get_current_sample_id(self) -> str:
        """Retorna o sample_id da amostra atual."""
        global_index = self.split_indices[self.current_index]
        individuals = self.dataset_metadata['individuals']
        return individuals[global_index]


# ==============================================================================
# FUNÇÕES DE CARREGAMENTO DE DADOS
# ==============================================================================

def load_cache_data(
    cache_dir: Path,
    split: str,
    index: int,
    gene_filter: Optional[Union[str, List[str]]] = None
) -> Tuple[np.ndarray, Dict, int]:
    """
    Carrega dados do cache processado.
    
    Args:
        cache_dir: Diretório do cache
        split: Conjunto de dados (train/val/test)
        index: Índice do indivíduo no split
        gene_filter: Gene(s) a filtrar. None = todos, string = um gene, lista = múltiplos genes
        
    Returns:
        Tupla (features_array, metadata, global_index)
        - features_array: Array numpy shape [num_genes*6, window_center_size]
        - metadata: Dicionário com metadados do cache
        - global_index: Índice global do indivíduo (0-77)
    """
    cache_dir = Path(cache_dir)
    
    # Carregar metadados
    metadata_file = cache_dir / 'metadata.json'
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata não encontrado: {metadata_file}")
    
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    
    # Carregar splits
    splits_file = cache_dir / 'splits.json'
    if not splits_file.exists():
        raise FileNotFoundError(f"Splits não encontrado: {splits_file}")
    
    with open(splits_file, 'r') as f:
        splits = json.load(f)
    
    # Obter índices do split selecionado
    split_key = f"{split}_indices"
    if split_key not in splits:
        raise ValueError(f"Split '{split}' não encontrado. Opções: train, val, test")
    
    split_indices = splits[split_key]
    
    if index >= len(split_indices):
        raise ValueError(f"Índice {index} fora do range para split '{split}' (max: {len(split_indices)-1})")
    
    global_index = split_indices[index]
    
    # Carregar dados do split
    data_file = cache_dir / f'{split}_data.pt'
    if not data_file.exists():
        raise FileNotFoundError(f"Arquivo de dados não encontrado: {data_file}")
    
    data = torch.load(data_file, map_location='cpu')
    
    # Extrair features e target do indivíduo
    features, target = data[index]
    features_array = features.numpy()  # Shape: [66, window_center_size]
    
    # NÃO filtramos aqui - o filtro será aplicado depois que soubermos o sample_id
    # (necessário para obter a ordem correta dos genes)
    
    return features_array, metadata, global_index


def filter_genes_from_features(
    features: np.ndarray,
    gene_filter: Union[str, List[str]],
    genes_in_order: List[str]
) -> Tuple[np.ndarray, List[str]]:
    """
    Filtra genes específicos do array de features.
    
    Args:
        features: Array shape [66, window_center_size]
        gene_filter: Gene(s) a filtrar
        genes_in_order: Lista de genes na ordem que aparecem nas features
        
    Returns:
        Tupla (filtered_features, genes_kept):
        - filtered_features: Array filtrado shape [num_genes*6, window_center_size]
        - genes_kept: Lista de genes mantidos
    """
    # Determinar quais genes filtrar
    if isinstance(gene_filter, str):
        genes_to_keep = [gene_filter]
    else:
        genes_to_keep = gene_filter
    
    # Verificar genes válidos
    for gene in genes_to_keep:
        if gene not in genes_in_order:
            raise ValueError(f"Gene '{gene}' não encontrado. Opções: {genes_in_order}")
    
    # Extrair índices das tracks correspondentes
    track_indices = []
    for gene in genes_to_keep:
        gene_idx = genes_in_order.index(gene)
        start_track = gene_idx * NUM_ONTOLOGIES
        end_track = start_track + NUM_ONTOLOGIES
        track_indices.extend(range(start_track, end_track))
    
    # Filtrar features
    filtered_features = features[track_indices, :]
    
    return filtered_features, genes_to_keep


def get_sample_id_from_index(dataset_dir: Path, global_index: int) -> str:
    """
    Mapeia índice global para sample_id.
    
    Args:
        dataset_dir: Diretório do dataset original
        global_index: Índice global do indivíduo (0-77)
        
    Returns:
        sample_id (ex: "HG00120")
    """
    dataset_metadata_file = dataset_dir / 'dataset_metadata.json'
    if not dataset_metadata_file.exists():
        raise FileNotFoundError(f"Dataset metadata não encontrado: {dataset_metadata_file}")
    
    with open(dataset_metadata_file, 'r') as f:
        dataset_metadata = json.load(f)
    
    individuals = dataset_metadata['individuals']
    
    if global_index >= len(individuals):
        raise ValueError(f"Índice global {global_index} fora do range (max: {len(individuals)-1})")
    
    sample_id = individuals[global_index]
    
    return sample_id


def load_alphagenome_data(
    dataset_dir: Path,
    sample_id: str,
    window_center_size: int,
    gene_filter: Optional[Union[str, List[str]]] = None,
    config: Optional[Dict] = None
) -> Tuple[np.ndarray, List[str]]:
    """
    Load AlphaGenome data from files OR API based on config.
    
    Args:
        dataset_dir: Diretório do dataset original
        sample_id: ID do indivíduo (ex: "HG00120")
        window_center_size: Tamanho do trecho central a extrair
        gene_filter: Gene(s) a carregar. None = todos, string = um gene, lista = múltiplos genes
        config: Configuration dict (optional, for API mode)
        
    Returns:
        Tupla (array, gene_list):
        - array: Array numpy shape [num_genes*6, window_center_size] com dados
        - gene_list: Lista de genes carregados
    """
    # Check if API mode is enabled
    use_api = config and config.get('alphagenome_api', {}).get('enabled', False)
    
    if use_api:
        return _load_from_alphagenome_api(
            dataset_dir, sample_id, window_center_size, gene_filter, config
        )
    else:
        return _load_from_files(
            dataset_dir, sample_id, window_center_size, gene_filter
        )


def _load_from_alphagenome_api(
    dataset_dir: Path,
    sample_id: str,
    window_center_size: int,
    gene_filter: Optional[Union[str, List[str]]],
    config: Dict
) -> Tuple[np.ndarray, List[str]]:
    """Load data from AlphaGenome API."""
    console.print(f"[bold cyan]🌐 AlphaGenome API Mode Enabled[/bold cyan]")
    
    # Get gene order and filter
    genes_in_order = get_genes_in_dataset_order(dataset_dir, sample_id)
    
    if gene_filter is None:
        genes_to_load = genes_in_order
    elif isinstance(gene_filter, str):
        if gene_filter not in genes_in_order:
            raise ValueError(f"Gene '{gene_filter}' não encontrado. Opções: {genes_in_order}")
        genes_to_load = [gene_filter]
    elif isinstance(gene_filter, list):
        for gene in gene_filter:
            if gene not in genes_in_order:
                raise ValueError(f"Gene '{gene}' não encontrado. Opções: {genes_in_order}")
        genes_to_load = [g for g in genes_in_order if g in gene_filter]
    else:
        raise ValueError(f"gene_filter inválido: {gene_filter}")
    
    ontology_terms = config['alphagenome_api']['ontology_terms']
    
    # Determine window size from window_center_size
    # Map common sizes to AlphaGenome constants
    window_size_map = {
        10000: "SEQUENCE_LENGTH_16KB",      # 10k fits in 16k
        100000: "SEQUENCE_LENGTH_128KB",    # 100k fits in 128k  
        524288: "SEQUENCE_LENGTH_512KB",    # 512k exact
        1048576: "SEQUENCE_LENGTH_1MB",     # 1MB exact
    }
    
    # Use window_center_size_key from config if available, otherwise infer from window_center_size
    if 'window_center_size_key' in config:
        window_size_key = config['window_center_size_key']
    else:
        window_size_key = window_size_map.get(window_center_size, "SEQUENCE_LENGTH_512KB")
    
    all_gene_data = []
    
    console.print(f"[cyan]Carregando {len(genes_to_load)} gene(s) via API (reference genome)...[/cyan]")
    console.print(f"[dim]  Window size: {window_size_key} (center: {window_center_size} bp)[/dim]")
    
    for gene_name in genes_to_load:
        # Call AlphaGenome API with individual's sequence from dataset_dir
        values = predict_with_alphagenome(
            gene_name=gene_name,
            config=config,
            ontology_terms=ontology_terms,
            window_size_key=window_size_key,
            sample_id=sample_id  # Pass sample_id to load individual's sequence
        )
        if values is None:
            raise RuntimeError(f"API prediction failed for {gene_name}")
        
        # Process each ontology (same as file-based)
        gene_tracks = []
        for ont_idx in range(values.shape[1]):
            track_array = values[:, ont_idx]
            
            # Extract center if needed
            if len(track_array) > window_center_size:
                sequence_length = len(track_array)
                center_idx = sequence_length // 2
                half_size = window_center_size // 2
                start_idx = max(0, center_idx - half_size)
                end_idx = min(sequence_length, center_idx + half_size)
                central_track = track_array[start_idx:end_idx]
            else:
                central_track = track_array
            
            gene_tracks.append(central_track)
        
        gene_array = np.array(gene_tracks)
        all_gene_data.append(gene_array)
    
    result = np.concatenate(all_gene_data, axis=0)
    console.print(f"[bold green]✓ API data loaded: shape {result.shape}[/bold green]")
    return result, genes_to_load


def load_raw_alphagenome_data(
    dataset_dir: Path,
    sample_id: str,
    gene_name: str,
    center_bp: int,
    source: str,
    config: Dict
) -> tuple:
    """
    Load raw AlphaGenome data without normalization.
    
    Args:
        dataset_dir: Dataset directory
        sample_id: Sample ID
        gene_name: Gene name
        center_bp: Number of base pairs to extract from center
        source: "files" or "api"
        config: Full config dict
    
    Returns:
        Array [6, center_bp] with raw RNA-seq values
    """
    if source == "api":
        # Call API with reference genome (como no Colab)
        ontology_terms = config['alphagenome_api']['ontology_terms']
        
        # Determine window size key from center_bp
        # Map common sizes to AlphaGenome constants
        window_size_map = {
            2048: "SEQUENCE_LENGTH_2KB",
            4096: "SEQUENCE_LENGTH_4KB",
            8192: "SEQUENCE_LENGTH_8KB",
            16384: "SEQUENCE_LENGTH_16KB",
            32768: "SEQUENCE_LENGTH_32KB",
            65536: "SEQUENCE_LENGTH_64KB",
            131072: "SEQUENCE_LENGTH_128KB",
            262144: "SEQUENCE_LENGTH_256KB",
            524288: "SEQUENCE_LENGTH_512KB",
            1048576: "SEQUENCE_LENGTH_1MB",
        }
        
        window_size_key = window_size_map.get(center_bp, "SEQUENCE_LENGTH_16KB")
        console.print(f"[cyan]  Using window size: {window_size_key} ({center_bp} bp)[/cyan]")
        
        result = predict_with_alphagenome(
            gene_name=gene_name,
            config=config,
            ontology_terms=ontology_terms,
            window_size_key=window_size_key,
            return_full_output=True,  # Get full output for Colab-style plotting
            sample_id=sample_id  # Pass sample_id to load individual's sequence
        )
        if result is None:
            raise RuntimeError(f"API prediction failed for {gene_name}")
        
        output, gtf = result
        values = np.array(output.rna_seq.values)
        
        # Convert to expected format [6, center_bp]
        gene_tracks = []
        for ont_idx in range(values.shape[1]):
            gene_tracks.append(values[:, ont_idx])
        
        return (np.array(gene_tracks), output, gtf)  # Shape: [6, center_bp], plus AlphaGenome objects
        
    else:
        # Read from .npz files
        gene_dir = dataset_dir / 'individuals' / sample_id / 'windows' / gene_name / 'predictions_H1'
        rna_seq_file = gene_dir / 'rna_seq.npz'
        if not rna_seq_file.exists():
            raise FileNotFoundError(f"File not found: {rna_seq_file}")
        
        data = np.load(rna_seq_file)
        values = data['values']  # Shape: (sequence_length, 6)
        
        # Extract center window for each ontology
        gene_tracks = []
        for ont_idx in range(values.shape[1]):
            track_array = values[:, ont_idx]
            
            # Extract center
            sequence_length = len(track_array)
            center_idx = sequence_length // 2
            half_size = center_bp // 2
            start_idx = max(0, center_idx - half_size)
            end_idx = min(sequence_length, center_idx + half_size)
            
            central_track = track_array[start_idx:end_idx]
            gene_tracks.append(central_track)
        
        return (np.array(gene_tracks), None, None)  # Shape: [6, center_bp], no AlphaGenome objects in file mode


def _load_from_files(
    dataset_dir: Path,
    sample_id: str,
    window_center_size: int,
    gene_filter: Optional[Union[str, List[str]]] = None
) -> Tuple[np.ndarray, List[str]]:
    """
    Load data from .npz files in dataset_dir (original implementation).
    
    Args:
        dataset_dir: Diretório do dataset original
        sample_id: ID do indivíduo (ex: "HG00120")
        window_center_size: Tamanho do trecho central a extrair
        gene_filter: Gene(s) a carregar. None = todos, string = um gene, lista = múltiplos genes
        
    Returns:
        Tupla (array, gene_list):
        - array: Array numpy shape [num_genes*6, window_center_size] com dados brutos
        - gene_list: Lista de genes carregados
    """
    dataset_dir = Path(dataset_dir)
    individual_dir = dataset_dir / 'individuals' / sample_id
    
    if not individual_dir.exists():
        raise FileNotFoundError(f"Diretório do indivíduo não encontrado: {individual_dir}")
    
    # Obter genes na mesma ordem que o dataset processa
    genes_in_order = get_genes_in_dataset_order(dataset_dir, sample_id)
    
    # Determinar quais genes carregar (SEMPRE NA ORDEM DO DATASET)
    if gene_filter is None:
        genes_to_load = genes_in_order
    elif isinstance(gene_filter, str):
        if gene_filter not in genes_in_order:
            raise ValueError(f"Gene '{gene_filter}' não encontrado. Opções: {genes_in_order}")
        genes_to_load = [gene_filter]
    elif isinstance(gene_filter, list):
        # Validar genes
        for gene in gene_filter:
            if gene not in genes_in_order:
                raise ValueError(f"Gene '{gene}' não encontrado. Opções: {genes_in_order}")
        # Manter a ordem do dataset, não a ordem do gene_filter!
        genes_to_load = [g for g in genes_in_order if g in gene_filter]
    else:
        raise ValueError(f"gene_filter inválido: {gene_filter}")
    
    all_gene_data = []
    
    for gene_name in genes_to_load:
        gene_dir = individual_dir / 'windows' / gene_name / 'predictions_H1'
        rna_seq_file = gene_dir / 'rna_seq.npz'
        
        if not rna_seq_file.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {rna_seq_file}")
        
        # Carregar .npz
        data = np.load(rna_seq_file)
        
        if 'values' not in data:
            raise ValueError(f"Key 'values' não encontrada em {rna_seq_file}")
        
        values = data['values']  # Shape: (1048576, 6) ou similar
        
        # Processar CADA COLUNA (ontologia) separadamente
        # (mesmo método usado em neural_ancestry_predictor.py)
        num_ontologies = values.shape[1]
        gene_tracks = []
        
        for ont_idx in range(num_ontologies):
            # Extrair coluna (uma ontologia por vez)
            track_array = values[:, ont_idx]
            
            # Extrair trecho central desta coluna
            sequence_length = len(track_array)
            center_idx = sequence_length // 2
            half_size = window_center_size // 2
            start_idx = max(0, center_idx - half_size)
            end_idx = min(sequence_length, center_idx + half_size)
            
            central_track = track_array[start_idx:end_idx]
            gene_tracks.append(central_track)
        
        # Empilhar as tracks como linhas (shape: [6, window_center_size])
        gene_array = np.array(gene_tracks)
        all_gene_data.append(gene_array)
    
    # Concatenar todos os genes: (num_genes × 6 ontologias, window_center_size)
    alphagenome_array = np.concatenate(all_gene_data, axis=0)
    
    return alphagenome_array, genes_to_load


def apply_normalization(
    data: np.ndarray,
    method: str,
    norm_params: Dict
) -> np.ndarray:
    """
    Aplica normalização aos dados do AlphaGenome.
    
    Args:
        data: Array numpy com dados brutos
        method: Método de normalização ('log', 'minmax_keep_zero', 'zscore')
        norm_params: Parâmetros de normalização (do metadata.json)
        
    Returns:
        Array numpy normalizado
    """
    # Converter para tensor PyTorch para usar as funções de normalização
    data_tensor = torch.from_numpy(data).float()
    
    if method == 'log':
        log_max = norm_params.get('log_max')
        if log_max is None:
            raise ValueError("log_max não encontrado nos parâmetros de normalização")
        normalized = log_normalize(data_tensor, log_max)
    
    elif method == 'minmax_keep_zero':
        max_val = norm_params.get('max')
        if max_val is None:
            raise ValueError("max não encontrado nos parâmetros de normalização")
        normalized = minmax_keep_zero(data_tensor, max_val)
    
    elif method == 'zscore':
        mean = norm_params.get('mean')
        std = norm_params.get('std')
        if mean is None or std is None:
            raise ValueError("mean ou std não encontrado nos parâmetros de normalização")
        normalized = (data_tensor - mean) / std
    
    else:
        raise ValueError(f"Método de normalização desconhecido: {method}")
    
    normalized_array = normalized.numpy()
    
    return normalized_array


# ==============================================================================
# FUNÇÕES DE ANÁLISE E VISUALIZAÇÃO
# ==============================================================================

def compute_metrics(cache_data: np.ndarray, alphagenome_data: np.ndarray, verbose: bool = True) -> Dict:
    """
    Calcula métricas de comparação entre cache e AlphaGenome.
    
    Args:
        cache_data: Array do cache [num_tracks, window_center_size]
        alphagenome_data: Array do AlphaGenome [num_tracks, window_center_size]
        verbose: Exibir métricas detalhadas
        
    Returns:
        Dicionário com métricas
    """
    num_tracks = cache_data.shape[0]
    
    # Métricas por track
    mae_per_track = []
    corr_per_track = []
    
    for track_idx in range(num_tracks):
        cache_track = cache_data[track_idx, :]
        alpha_track = alphagenome_data[track_idx, :]
        
        # MAE (Mean Absolute Error)
        mae = np.mean(np.abs(cache_track - alpha_track))
        mae_per_track.append(mae)
        
        # Correlação de Pearson
        try:
            corr, _ = pearsonr(cache_track, alpha_track)
            corr_per_track.append(corr)
        except:
            corr_per_track.append(np.nan)
    
    mae_per_track = np.array(mae_per_track)
    corr_per_track = np.array(corr_per_track)
    
    # Estatísticas globais
    metrics = {
        'mae_per_track': mae_per_track,
        'corr_per_track': corr_per_track,
        'mae_mean': np.mean(mae_per_track),
        'mae_max': np.max(mae_per_track),
        'mae_min': np.min(mae_per_track),
        'mae_max_track': np.argmax(mae_per_track),
        'corr_mean': np.nanmean(corr_per_track),
        'corr_min': np.nanmin(corr_per_track),
        'corr_max': np.nanmax(corr_per_track),
        'corr_min_track': np.nanargmin(corr_per_track),
    }
    
    # Exibir estatísticas se verbose
    if verbose:
        console.print("\n[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]")
        console.print("[bold cyan]           MÉTRICAS DE COMPARAÇÃO                      [/bold cyan]")
        console.print("[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]")
        
        table = Table(show_header=True, header_style="bold magenta")
        table.add_column("Métrica", style="cyan")
        table.add_column("Valor", justify="right", style="green")
        
        table.add_row("MAE Média (global)", f"{metrics['mae_mean']:.6f}")
        table.add_row("MAE Máximo", f"{metrics['mae_max']:.6f}")
        table.add_row("MAE Mínimo", f"{metrics['mae_min']:.6f}")
        table.add_row("Track com maior MAE", f"{metrics['mae_max_track']} (gene {metrics['mae_max_track']//6}, ont {metrics['mae_max_track']%6})")
        table.add_row("", "")
        table.add_row("Correlação Média", f"{metrics['corr_mean']:.6f}")
        table.add_row("Correlação Mínima", f"{metrics['corr_min']:.6f}")
        table.add_row("Correlação Máxima", f"{metrics['corr_max']:.6f}")
        table.add_row("Track com menor corr.", f"{metrics['corr_min_track']} (gene {metrics['corr_min_track']//6}, ont {metrics['corr_min_track']%6})")
        
        console.print(table)
        console.print("[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]\n")
    
    return metrics


def plot_comparison(
    cache_data: np.ndarray,
    alphagenome_data: np.ndarray,
    sample_id: str,
    metrics: Dict,
    genes_displayed: List[str],
    config: Dict,
    viewer: Optional[InteractiveViewer] = None
) -> plt.Figure:
    """
    Plota gráfico de comparação com subplots independentes para cada track.
    
    Args:
        cache_data: Array do cache [num_tracks, window_center_size]
        alphagenome_data: Array do AlphaGenome [num_tracks, window_center_size]
        sample_id: ID do indivíduo
        metrics: Dicionário com métricas de comparação
        genes_displayed: Lista de genes exibidos
        config: Configuração do programa
        viewer: Viewer interativo (opcional)
        
    Returns:
        Figura matplotlib
    """
    num_tracks, window_size = cache_data.shape
    
    # Subsampling para visualização (plotar 1 a cada N pontos)
    subsample = max(1, window_size // 2000)  # Máximo de 2000 pontos no gráfico
    x_positions = np.arange(0, window_size, subsample)
    
    # Criar figura com subplots (um para cada track)
    fig, axes = plt.subplots(num_tracks, 1, figsize=(16, num_tracks * 1.5), sharex=True)
    
    # Se apenas uma track, axes não é array
    if num_tracks == 1:
        axes = [axes]
    
    # Plotar cada track em seu próprio subplot
    for track_idx in range(num_tracks):
        ax = axes[track_idx]
        
        # Obter dados de ambas as fontes
        cache_track = cache_data[track_idx, ::subsample]
        alpha_track = alphagenome_data[track_idx, ::subsample]
        
        # Plotar
        ax.plot(x_positions, cache_track, 
                color='blue', linewidth=1.0, alpha=0.7, label='Cache')
        ax.plot(x_positions, alpha_track,
                color='red', linewidth=1.0, linestyle='--', alpha=0.7, label='AlphaGenome')
        
        # Criar label com gene e ontologia
        gene_idx = track_idx // NUM_ONTOLOGIES
        ont_idx = track_idx % NUM_ONTOLOGIES
        if gene_idx < len(genes_displayed):
            gene_name = genes_displayed[gene_idx]
            ylabel = f"{gene_name}\nOnt{ont_idx}"
        else:
            ylabel = f"Track {track_idx}"
        
        # Configurar subplot
        ax.set_ylabel(ylabel, fontsize=9, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='upper right', fontsize=8)
        
        # Mostrar MAE e correlação para esta track
        track_mae = metrics['mae_per_track'][track_idx]
        track_corr = metrics['corr_per_track'][track_idx]
        ax.text(0.02, 0.98, f'MAE: {track_mae:.4f} | Corr: {track_corr:.3f}',
                transform=ax.transAxes,
                fontsize=7,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Adicionar linha horizontal entre genes
        if ont_idx == NUM_ONTOLOGIES - 1 and gene_idx < len(genes_displayed) - 1:
            ax.axhline(y=ax.get_ylim()[0], color='gray', linewidth=2, alpha=0.8)
    
    # Título principal
    gene_info = f"Genes: {', '.join(genes_displayed)}" if len(genes_displayed) <= 3 else f"{len(genes_displayed)} genes"
    title = f'Sample {sample_id} - Cache vs AlphaGenome ({num_tracks} tracks)\n'
    title += f'{gene_info} | MAE médio: {metrics["mae_mean"]:.6f} | Corr média: {metrics["corr_mean"]:.6f}'
    
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    # Label do eixo X apenas no último subplot
    axes[-1].set_xlabel('Posição na janela (bp)', fontsize=12, fontweight='bold')
    
    # Texto com instruções de navegação (se modo interativo)
    if viewer and config.get('show_navigation_help', True):
        fig.text(0.5, 0.01, '← Anterior | → Próxima | Q Sair',
                ha='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    # Nota explicativa
    note_text = 'Nota: Dados do cache já estão normalizados. Dados do AlphaGenome foram normalizados com mesmos parâmetros para comparação.'
    fig.text(0.5, 0.98, note_text,
            ha='center', fontsize=8, style='italic',
            transform=fig.transFigure)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    
    return fig


def plot_raw_data_colab_style(
    output,
    gene_name: str,
    sample_id: str,
    gtf,  # pandas DataFrame
    config: Dict,
    viewer: Optional[InteractiveViewer] = None
) -> plt.Figure:
    """
    Plot using AlphaGenome's plot_components (Colab style).
    
    Args:
        output: AlphaGenome output object with rna_seq attribute
        gene_name: Gene name
        sample_id: Sample ID
        gtf: GTF dataframe
        config: Configuration dict
        viewer: Optional interactive viewer
    
    Returns:
        Figure object
    """
    try:
        from alphagenome.data import transcript as transcript_utils
        from alphagenome.data import gene_annotation
        from alphagenome.visualization import plot_components
    except ImportError as e:
        console.print(f"[red]Error importing visualization modules: {e}[/red]")
        console.print(f"[yellow]Falling back to simple plot...[/yellow]")
        # Fallback to simple plot
        values = np.array(output.rna_seq.values)
        gene_tracks = []
        for ont_idx in range(values.shape[1]):
            gene_tracks.append(values[:, ont_idx])
        raw_data = np.array(gene_tracks)
        return plot_raw_data(raw_data, sample_id, gene_name, config, viewer)
    
    console.print(f"[cyan]Creating Colab-style visualization...[/cyan]")
    
    # Extract transcripts
    gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
    try:
        # Try to use MANE select if available
        gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
    except (ImportError, AttributeError):
        # If not available, use longest transcript
        pass
    
    transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)
    transcripts = transcript_extractor.extract(output.rna_seq.interval)
    
    console.print(f"[green]  ✓ Extracted {len(transcripts)} transcript(s)[/green]")
    
    # Create plot using AlphaGenome's plot_components (como no Colab)
    fig = plot_components.plot(
        components=[
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Tracks(output.rna_seq),
        ],
        interval=output.rna_seq.interval,
    )
    
    # Add title
    plt.suptitle(f"AlphaGenome RNA-seq (Reference Genome): {gene_name} ({sample_id})", 
                 fontsize=14, y=0.98, fontweight='bold')
    
    return fig


def plot_raw_data(
    raw_data: np.ndarray,
    sample_id: str,
    gene_name: str,
    config: Dict,
    viewer: Optional[InteractiveViewer] = None
) -> plt.Figure:
    """
    Plot raw AlphaGenome data without comparison (simple subplot style).
    
    Args:
        raw_data: Array [6, length] with raw values
        sample_id: Sample ID
        gene_name: Gene name
        config: Config dict
        viewer: Interactive viewer (optional)
    
    Returns:
        Matplotlib figure
    """
    num_tracks = raw_data.shape[0]
    fig, axes = plt.subplots(num_tracks, 1, figsize=(14, 2*num_tracks), sharex=True)
    
    if num_tracks == 1:
        axes = [axes]
    
    fig.suptitle(f'AlphaGenome Raw Data - {sample_id} - {gene_name}', 
                 fontsize=14, fontweight='bold')
    
    for idx in range(num_tracks):
        ax = axes[idx]
        track_data = raw_data[idx, :]
        
        # Plot raw data
        ax.plot(track_data, linewidth=0.8, color='blue', alpha=0.8)
        
        # Labels
        ax.set_ylabel(f'Ontology {idx}', fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # Show value range
        vmin, vmax = track_data.min(), track_data.max()
        ax.text(0.99, 0.95, f'[{vmin:.2e}, {vmax:.2e}]',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=8, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    axes[-1].set_xlabel('Position (bp)', fontsize=10)
    
    plt.tight_layout()
    return fig


# ==============================================================================
# FUNÇÃO PRINCIPAL
# ==============================================================================

def process_sample_raw_mode(
    config: Dict,
    index: int,
    dataset_metadata: Dict,
    viewer: Optional[InteractiveViewer] = None
) -> bool:
    """Process sample in raw visualization mode (no cache comparison)."""
    
    console.print(f"\n[bold cyan]RAW MODE: AlphaGenome Only (No Normalization)[/bold cyan]")
    
    # Get sample_id from dataset_dir
    dataset_dir = Path(config['dataset_dir'])
    
    # Determine sample_id
    if config.get('sample_id'):
        sample_id = config['sample_id']
    else:
        # Get from dataset metadata
        sample_id = dataset_metadata['individuals'][index]
    
    console.print(f"[green]✓ Sample: {sample_id}[/green]")
    
    # Get gene filter
    gene_filter = config.get('gene_filter')
    if isinstance(gene_filter, str):
        genes_to_plot = [gene_filter]
    elif isinstance(gene_filter, list):
        genes_to_plot = gene_filter
    else:
        # Load all genes
        genes_to_plot = get_genes_in_dataset_order(dataset_dir, sample_id)
    
    # Load raw data for each gene
    raw_mode_config = config['raw_mode']
    source = raw_mode_config['source']
    window_size_key = raw_mode_config['window_size_key']
    
    # Convert window_size_key to actual bp
    window_size_map = {
        'SEQUENCE_LENGTH_16KB': 16384,
        'SEQUENCE_LENGTH_100KB': 102400,
        'SEQUENCE_LENGTH_500KB': 512000,
        'SEQUENCE_LENGTH_1MB': 1048576
    }
    center_bp = window_size_map.get(window_size_key, 102400)
    
    for gene_name in genes_to_plot:
        console.print(f"\n[cyan]Loading raw data for {gene_name}...[/cyan]")
        
        raw_data, output, gtf = load_raw_alphagenome_data(
            dataset_dir, sample_id, gene_name, center_bp, source, config
        )
        
        # Plot with Colab-style if we have output object, otherwise simple plot
        if output is not None and gtf is not None:
            fig = plot_raw_data_colab_style(output, gene_name, sample_id, gtf, config, viewer)
        else:
            fig = plot_raw_data(raw_data, sample_id, gene_name, config, viewer)
        
        # Save if requested
        if config.get('save_plots', False):
            output_dir = config.get('output_dir')
            if output_dir:
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)
                filename = f"{config.get('output_prefix', 'raw')}_{sample_id}_{gene_name}.png"
                fig.savefig(output_path / filename, dpi=150, bbox_inches='tight')
                console.print(f"[green]✓ Saved: {output_path / filename}[/green]")
        
        plt.show()
    
    return True


def process_sample(
    config: Dict,
    index: int,
    splits_data: Dict,
    dataset_metadata: Dict,
    norm_params: Dict,
    viewer: Optional[InteractiveViewer] = None
) -> Optional[plt.Figure]:
    """
    Processa e visualiza uma amostra.
    
    Args:
        config: Configuração do programa
        index: Índice da amostra no split
        splits_data: Dados dos splits
        dataset_metadata: Metadata do dataset
        norm_params: Parâmetros de normalização
        viewer: Viewer interativo (opcional)
        
    Returns:
        Figura matplotlib ou None se erro
    """
    # Check if raw mode is enabled
    if config.get('raw_mode', {}).get('enabled', False):
        return process_sample_raw_mode(config, index, dataset_metadata, viewer)
    
    try:
        # Carregar dados do cache (todas as tracks)
        cache_features, cache_metadata, global_index = load_cache_data(
            Path(config['cache_dir']),
            config['split'],
            index,
            None  # Sem filtro ainda
        )
        
        # Obter sample_id
        sample_id = get_sample_id_from_index(Path(config['dataset_dir']), global_index)
        console.print(f"[green]✓ Sample: {sample_id} (índice {index}, global {global_index})[/green]")
        
        # Obter ordem dos genes no dataset
        genes_in_order = get_genes_in_dataset_order(Path(config['dataset_dir']), sample_id)
        
        # Aplicar filtro de genes ao cache se necessário
        gene_filter = config.get('gene_filter')
        genes_loaded = genes_in_order  # Por padrão, todos os genes
        if gene_filter is not None:
            cache_features, genes_loaded = filter_genes_from_features(
                cache_features, gene_filter, genes_in_order
            )
        
        # Carregar dados originais do AlphaGenome
        window_center_size = cache_metadata['processing_params']['window_center_size']
        alphagenome_features, _ = load_alphagenome_data(
            Path(config['dataset_dir']),
            sample_id,
            window_center_size,
            gene_filter,
            config  # Pass config for API mode
        )
        
        # Aplicar mesma normalização aos dados do AlphaGenome
        normalization_method = cache_metadata['processing_params']['normalization_method']
        alphagenome_normalized = apply_normalization(
            alphagenome_features, normalization_method, norm_params
        )
        
        # Calcular métricas de comparação
        verbose_metrics = config.get('verbose_metrics', True)
        metrics = compute_metrics(cache_features, alphagenome_normalized, verbose=verbose_metrics)
        
        # Visualizar
        fig = plot_comparison(
            cache_features, alphagenome_normalized, sample_id, metrics,
            genes_loaded, config, viewer
        )
        
        # Salvar gráfico se configurado
        if config.get('save_plots', False) and config.get('output_dir'):
            output_dir = Path(config['output_dir'])
            output_dir.mkdir(parents=True, exist_ok=True)
            output_prefix = config.get('output_prefix', 'verify')
            output_file = output_dir / f"{output_prefix}_{sample_id}.png"
            fig.savefig(output_file, dpi=150, bbox_inches='tight')
            console.print(f"[green]✓ Gráfico salvo em: {output_file}[/green]")
        
        return fig
    
    except Exception as e:
        console.print(f"[red]ERRO ao processar amostra {index}: {e}[/red]")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description='Verifica e compara dados processados em cache com dados originais do AlphaGenome',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--config',
        type=Path,
        required=True,
        help='Arquivo de configuração YAML'
    )
    
    args = parser.parse_args()
    
    # Carregar configuração
    if not args.config.exists():
        console.print(f"[red]ERRO: Arquivo de configuração não existe: {args.config}[/red]")
        sys.exit(1)
    
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    # Validate raw_mode configuration if enabled
    if config.get('raw_mode', {}).get('enabled', False):
        raw_mode = config['raw_mode']
        
        # Validate window_size_key
        valid_window_sizes = ['SEQUENCE_LENGTH_16KB', 'SEQUENCE_LENGTH_100KB', 
                             'SEQUENCE_LENGTH_500KB', 'SEQUENCE_LENGTH_1MB']
        window_size_key = raw_mode.get('window_size_key', 'SEQUENCE_LENGTH_100KB')
        if window_size_key not in valid_window_sizes:
            console.print(f"[red]ERRO: window_size_key inválido: {window_size_key}[/red]")
            console.print(f"[red]Opções válidas: {', '.join(valid_window_sizes)}[/red]")
            sys.exit(1)
        
        # Validate source
        valid_sources = ['files', 'api']
        source = raw_mode.get('source', 'files')
        if source not in valid_sources:
            console.print(f"[red]ERRO: source inválido: {source}[/red]")
            console.print(f"[red]Opções válidas: {', '.join(valid_sources)}[/red]")
            sys.exit(1)
        
        # If API source, check API configuration
        if source == 'api':
            if not config.get('alphagenome_api', {}).get('enabled', False):
                console.print(f"[red]ERRO: raw_mode.source='api' requer alphagenome_api.enabled=true[/red]")
                sys.exit(1)
            
            api_key = config['alphagenome_api'].get('api_key') or os.environ.get('ALPHAGENOME_API_KEY')
            if not api_key:
                console.print(f"[red]ERRO: API key não configurada[/red]")
                console.print(f"[yellow]Configure alphagenome_api.api_key ou ALPHAGENOME_API_KEY env var[/yellow]")
                sys.exit(1)
    
    console.print("\n[bold green]════════════════════════════════════════════════════════════[/bold green]")
    console.print("[bold green]       VERIFICAÇÃO DE DATASET PROCESSADO                   [/bold green]")
    console.print("[bold green]════════════════════════════════════════════════════════════[/bold green]\n")
    
    # Validar diretórios
    dataset_dir = Path(config['dataset_dir'])
    
    # Raw mode doesn't require cache_dir
    if config.get('raw_mode', {}).get('enabled', False):
        if not dataset_dir.exists():
            console.print(f"[red]ERRO: Dataset dir não existe: {dataset_dir}[/red]")
            sys.exit(1)
        
        # Skip cache loading in raw mode
        cache_dir = None
    else:
        cache_dir = Path(config['cache_dir'])
        
        if not cache_dir.exists():
            console.print(f"[red]ERRO: Cache dir não existe: {cache_dir}[/red]")
            sys.exit(1)
        
        if not dataset_dir.exists():
            console.print(f"[red]ERRO: Dataset dir não existe: {dataset_dir}[/red]")
            sys.exit(1)
    
    # Carregar dados necessários
    try:
        # Carregar dataset metadata (always needed)
        with open(dataset_dir / 'dataset_metadata.json', 'r') as f:
            dataset_metadata = json.load(f)
        
        # Raw mode: different processing path
        if config.get('raw_mode', {}).get('enabled', False):
            # Raw mode doesn't need splits, norm_params
            index = config.get('index', 0)
            success = process_sample_raw_mode(config, index, dataset_metadata)
            
            if not success:
                console.print("[red]Erro ao processar amostra em raw mode[/red]")
                sys.exit(1)
            
            console.print("\n[bold green]✓ Raw mode visualization concluída![/bold green]\n")
        else:
            # Normal comparison mode
            # Carregar splits
            with open(cache_dir / 'splits.json', 'r') as f:
                splits_data = json.load(f)
            
            # Carregar normalization params
            with open(cache_dir / 'normalization_params.json', 'r') as f:
                norm_params = json.load(f)
            
            # Modo interativo ou não
            if config.get('interactive_mode', False):
                viewer = InteractiveViewer(config, splits_data, dataset_metadata)
                
                # Loop interativo
                while not viewer.should_exit:
                    fig = process_sample(config, viewer.current_index, splits_data, 
                                       dataset_metadata, norm_params, viewer)
                    
                    if fig is not None:
                        viewer.fig = fig
                        # Conectar handler de teclado
                        fig.canvas.mpl_connect('key_press_event', viewer.on_key_press)
                        plt.show()
                    else:
                        # Erro ao processar, sair
                        break
                
                console.print("\n[bold green]✓ Sessão interativa encerrada[/bold green]\n")
            
            else:
                # Modo não-interativo: processa apenas uma amostra
                index = config.get('index', 0)
                fig = process_sample(config, index, splits_data, dataset_metadata, norm_params)
                
                if fig is not None:
                    plt.show()
                    console.print("\n[bold green]✓ Verificação concluída com sucesso![/bold green]\n")
                else:
                    sys.exit(1)
    
    except Exception as e:
        console.print(f"\n[bold red]ERRO: {e}[/bold red]")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
