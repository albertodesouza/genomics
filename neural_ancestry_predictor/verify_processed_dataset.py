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


# ═══════════════════════════════════════════════════════════════════
# Helper Function: Center Window Extraction
# ═══════════════════════════════════════════════════════════════════

def extract_center_window(
    data: np.ndarray,
    target_size: int,
    axis: int = 0
) -> np.ndarray:
    """
    Extrai janela central de um array numpy de forma consistente.
    
    IMPORTANTE: Usa sempre a fórmula end_idx = start_idx + target_size
    para garantir que o tamanho extraído seja EXATAMENTE target_size,
    mesmo com números ímpares.
    
    Args:
        data: Array numpy de entrada
        target_size: Tamanho desejado da janela central
        axis: Eixo ao longo do qual extrair (default: 0)
    
    Returns:
        Array com janela central extraída
        
    Examples:
        # Para array 1D ou extração ao longo do primeiro eixo
        >>> data = np.arange(1000000)
        >>> window = extract_center_window(data, 500000)
        >>> window.shape
        (500000,)
        
        # Para array 2D, extraindo ao longo de colunas
        >>> data = np.random.rand(6, 1048576)
        >>> window = extract_center_window(data, 524288, axis=1)
        >>> window.shape
        (6, 524288)
    """
    sequence_length = data.shape[axis]
    
    # Se já tem o tamanho correto, retorna direto
    if sequence_length == target_size:
        return data
    
    # Se é menor, retorna o que tem (com warning)
    if sequence_length < target_size:
        console.print(f"[yellow]  Warning: sequence_length ({sequence_length}) < target_size ({target_size})[/yellow]")
        console.print(f"[yellow]  Using full sequence[/yellow]")
        return data
    
    # Calcular índices da janela central
    # SEMPRE usar: end_idx = start_idx + target_size (não center_idx + half_size!)
    center_idx = sequence_length // 2
    half_size = target_size // 2
    start_idx = center_idx - half_size
    end_idx = start_idx + target_size
    
    # Garantir que não ultrapassa os limites
    if start_idx < 0:
        start_idx = 0
        end_idx = target_size
    elif end_idx > sequence_length:
        end_idx = sequence_length
        start_idx = end_idx - target_size
    
    # Extrair ao longo do eixo especificado
    if axis == 0:
        return data[start_idx:end_idx]
    elif axis == 1:
        return data[:, start_idx:end_idx]
    elif axis == 2:
        return data[:, :, start_idx:end_idx]
    else:
        raise ValueError(f"Axis {axis} not supported")


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
    fai_path = Path(str(fasta_path) + ".fai")
    
    try:
        import subprocess
        
        if not fasta_path.exists():
            console.print(f"[yellow]  Reference FASTA not found at {fasta_path}, skipping FASTA export[/yellow]")
            return None
        
        # Ler tamanhos dos cromossomos do .fai
        chrom_sizes = {}
        with open(fai_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
        
        # Alguns objetos Interval usam .chrom, outros .chromosome; tentamos ambos
        chrom = getattr(interval, "chrom", None)
        if chrom is None:
            chrom = getattr(interval, "chromosome", None)
        if chrom is None:
            raise AttributeError("Interval object has no 'chrom' or 'chromosome' attribute")
        
        chrom_size = chrom_sizes.get(chrom)
        
        # Coordenadas solicitadas (0-based)
        requested_start = interval.start
        requested_end = interval.end
        
        # Clipar coordenadas para valores válidos
        actual_start = max(0, requested_start)
        actual_end = requested_end
        if chrom_size is not None:
            actual_end = min(requested_end, chrom_size)
        
        # AlphaGenome usa coordenadas 0-based [start, end),
        # samtools faidx usa 1-based inclusivo: start+1 .. end
        start_1based = actual_start + 1
        end_1based = actual_end
        region = f"{chrom}:{start_1based}-{end_1based}"
        
        console.print(f"[cyan]  Extracting reference interval with samtools faidx: {region}[/cyan]")
        
        if requested_start < 0:
            console.print(f"[yellow]  Warning: Requested start was negative ({requested_start}), clipped to {actual_start}[/yellow]")
        if chrom_size and requested_end > chrom_size:
            console.print(f"[yellow]  Warning: Requested end exceeded chromosome size ({requested_end} > {chrom_size}), clipped to {actual_end}[/yellow]")
        
        interval_fasta_result = subprocess.run(
            ["samtools", "faidx", str(fasta_path), region],
            capture_output=True,
            text=True
        )
        
        if interval_fasta_result.returncode != 0:
            console.print(f"[red]  samtools faidx failed: {interval_fasta_result.stderr}[/red]")
            return None
        
        # Converter FASTA para sequência contínua
        fasta_text = interval_fasta_result.stdout
        seq_lines = []
        for line in fasta_text.split("\n"):
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line.strip())
        
        seq_raw = "".join(seq_lines)
        
        # Aplicar padding se necessário
        expected_length = requested_end - requested_start
        actual_length = actual_end - actual_start
        
        # Calcular N's no início (quando próximo ao início do cromossomo)
        n_prefix = 0
        if requested_start < 0:
            n_prefix = actual_start - requested_start  # Ex: 0 - (-100) = 100
        
        # Calcular N's no final (quando próximo ao final do cromossomo)
        n_suffix = 0
        if actual_length < expected_length:
            n_suffix = expected_length - actual_length - n_prefix
        
        # Construir sequência com padding
        seq = 'N' * n_prefix + seq_raw + 'N' * n_suffix
        
        if n_prefix > 0 or n_suffix > 0:
            console.print(f"[yellow]  Applied padding: {n_prefix} N's at start, {n_suffix} N's at end[/yellow]")
            console.print(f"[yellow]  Final sequence length: {len(seq)} bp (expected: {expected_length})[/yellow]")
        
        # Salvar FASTA em arquivo com sequência padded
        interval_fasta_output_path = Path("reference_interval.fasta")
        with open(interval_fasta_output_path, "w") as f:
            f.write(f">{chrom}:{requested_start}-{requested_end}\n")
            # Wrap at 60 characters per line
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
        console.print(f"[green]  ✓ Saved interval FASTA to: {interval_fasta_output_path}[/green]")
               
        if seq is None or len(seq) == 0:
            console.print("[red]  Cannot call predict_sequence: reference extraction failed[/red]")
            sys.exit(1)

        return seq
        
    except Exception as e:
        console.print(f"[red]  Error extracting reference sequence: {e}[/red]")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def _load_individual_haplotype_sequence(
    sample_id: str,
    gene_name: str,
    target_length: int,
    config: Dict
) -> Optional[str]:
    """
    Carrega a sequência do haplótipo H1 do indivíduo de dataset_dir.
    
    Args:
        sample_id: ID do indivíduo
        gene_name: Nome do gene
        target_length: Tamanho da janela a extrair do centro
        config: Configuration dict
    
    Returns:
        Sequência do haplótipo como string, ou None se falhar
    """
    dataset_dir = Path(config.get('dataset_dir', '/dados/GENOMICS_DATA/top3/non_longevous_results_genes'))
    try:
        console.print(f"[cyan]  Loading individual's haplotype sequence from dataset_dir...[/cyan]")
        
        # Caminho para a sequência do haplótipo H1
        # Padrão: dataset_dir/individuals/SAMPLE_ID/windows/GENE_NAME/SAMPLE_ID.H1.window.fixed.fa
        haplotype_file = dataset_dir / "individuals" / sample_id / "windows" / gene_name / f"{sample_id}.H1.window.fixed.fa"
        
        if not haplotype_file.exists():
            console.print(f"[yellow]  Haplotype file not found: {haplotype_file}[/yellow]")
            return None
        
        # Ler o arquivo FASTA do haplótipo
        console.print(f"[cyan]  Reading: {haplotype_file}[/cyan]")
        
        with open(haplotype_file, 'r') as f:
            haplotype_lines = []
            for line in f:
                if not line or line.startswith(">"):
                    continue
                haplotype_lines.append(line.strip())
        
        haplotype_seq = "".join(haplotype_lines)
        
        # Extrair a janela central correspondente ao target_length
        seq_length = len(haplotype_seq)
        
        console.print(f"[dim]  Haplotype full length: {seq_length} bp[/dim]")
        console.print(f"[dim]  Target length: {target_length} bp[/dim]")
        
        # Extrair centro usando função centralizada
        # Converter string para array, extrair, converter de volta
        seq_array = np.array(list(haplotype_seq))
        seq_array_center = extract_center_window(seq_array, target_length, axis=0)
        individual_seq = ''.join(seq_array_center)
        
        console.print(f"[green]  ✓ Loaded individual sequence: {len(individual_seq)} bp[/green]")
        console.print(f"[green]  ✓ Using {sample_id}'s H1 haplotype sequence[/green]")
        
        return individual_seq
        
    except Exception as e:
        console.print(f"[yellow]  Warning: Could not load individual sequence: {e}[/yellow]")
        import traceback
        traceback.print_exc()
        return None


def _detect_window_size_key_from_length(length: int) -> str:
    """
    Detecta a window_size_key correspondente ao tamanho em bp.
    
    Args:
        length: Tamanho em base pairs
    
    Returns:
        window_size_key correspondente (ex: "SEQUENCE_LENGTH_1MB")
    """
    length_to_key = {
        2048: "SEQUENCE_LENGTH_2KB",
        4096: "SEQUENCE_LENGTH_4KB",
        8192: "SEQUENCE_LENGTH_8KB",
        16384: "SEQUENCE_LENGTH_16KB",
        32768: "SEQUENCE_LENGTH_32KB",
        65536: "SEQUENCE_LENGTH_64KB",
        102400: "SEQUENCE_LENGTH_100KB",
        131072: "SEQUENCE_LENGTH_100KB",  # AlphaGenome usa 100KB para 128KB
        262144: "SEQUENCE_LENGTH_256KB",
        512000: "SEQUENCE_LENGTH_500KB",
        524288: "SEQUENCE_LENGTH_500KB",  # AlphaGenome usa 500KB para 512KB (524288 bp)
        1048576: "SEQUENCE_LENGTH_1MB",
    }
    
    key = length_to_key.get(length)
    if key is None:
        raise ValueError(f"Unsupported length: {length}. Supported: {list(length_to_key.keys())}")
    
    return key


def _load_dataset_dir_predictions(
    dataset_dir: Path,
    sample_id: str,
    gene_name: str,
    window_size_key: Optional[str] = None
) -> Tuple[np.ndarray, int]:
    """
    Carrega predições do dataset_dir.
    Se window_size_key for None, retorna dados completos.
    Se window_size_key for especificado, extrai janela central.
    
    Args:
        dataset_dir: Diretório raiz do dataset
        sample_id: ID do indivíduo
        gene_name: Nome do gene
        window_size_key: Chave do tamanho da janela (ex: "SEQUENCE_LENGTH_16KB"), 
                        ou None para retornar dados completos
    
    Returns:
        Tuple (values, full_length):
        - values: Array [window_size, 6] (se window_size_key) ou [full_length, 6] (se None)
        - full_length: Tamanho completo original dos dados
    """
    # Caminho para o arquivo de predições
    predictions_file = (dataset_dir / "individuals" / sample_id / "windows" / 
                       gene_name / "predictions_H1" / "rna_seq.npz")
    
    if not predictions_file.exists():
        raise FileNotFoundError(f"Predictions file not found: {predictions_file}")
    
    console.print(f"[cyan]  Loading predictions from: {predictions_file}[/cyan]")
    
    # Carregar arquivo .npz
    data = np.load(predictions_file)
    values = data['values']  # Shape esperado: [sequence_length, 6]
    
    full_length = values.shape[0]
    console.print(f"[dim]  Loaded shape: {values.shape} (full length: {full_length} bp)[/dim]")
    
    # Se window_size_key não especificado, retornar dados completos
    if window_size_key is None:
        console.print(f"[green]  ✓ Returning full data (no windowing)[/green]")
        return values, full_length
    
    # Caso contrário, extrair janela central
    # Mapa de window_size_key para número de bases
    window_size_map = {
        "SEQUENCE_LENGTH_2KB": 2048,
        "SEQUENCE_LENGTH_4KB": 4096,
        "SEQUENCE_LENGTH_8KB": 8192,
        "SEQUENCE_LENGTH_16KB": 16384,
        "SEQUENCE_LENGTH_32KB": 32768,
        "SEQUENCE_LENGTH_64KB": 65536,
        "SEQUENCE_LENGTH_100KB": 131072,  # AlphaGenome: 100KB = 131072 bp (128 KB)
        "SEQUENCE_LENGTH_128KB": 131072,
        "SEQUENCE_LENGTH_256KB": 262144,
        "SEQUENCE_LENGTH_500KB": 524288,  # AlphaGenome: 500KB = 524288 bp (512 KB)
        "SEQUENCE_LENGTH_512KB": 524288,  # Alias (use 500KB)
        "SEQUENCE_LENGTH_1MB": 1048576,
    }
    
    target_length = window_size_map.get(window_size_key)
    if target_length is None:
        raise ValueError(f"Invalid window_size_key: {window_size_key}")
    
    # Extrair janela central usando função centralizada
    window_values = extract_center_window(values, target_length, axis=0)
    
    console.print(f"[green]  ✓ Extracted center window: {window_values.shape}[/green]")
    
    return window_values, full_length


def predict_with_alphagenome_interval(
    gene_name: str,
    config: Dict,
    ontology_terms: List[str],
    window_size_key: str = "SEQUENCE_LENGTH_16KB",
    return_full_output: bool = False
) -> Optional[np.ndarray]:
    """
    Get predictions from AlphaGenome API usando predict_interval (sem extrair FASTA).
    
    Esta função usa predict_interval diretamente, deixando o AlphaGenome buscar
    a sequência do genoma de referência internamente. 
    
    Comparado a predict_with_alphagenome, esta função:
    - NÃO extrai o FASTA do genoma de referência
    - Chama client.predict_interval() diretamente
    - Deixa o AlphaGenome buscar a sequência internamente
    
    Args:
        gene_name: Nome do gene (e.g., 'MC1R')
        config: Dicionário de configuração
        ontology_terms: Lista de termos de ontologia (tissue/cell types)
        window_size_key: Tamanho da janela (e.g., 'SEQUENCE_LENGTH_16KB')
        return_full_output: Se True, retorna (values, output, gtf) ao invés de apenas values
    
    Returns:
        Array numpy com predições [window_size, 6] ou None se falhar
        Se return_full_output=True, retorna (values, output, gtf)
    """
    if not ALPHAGENOME_AVAILABLE:
        raise ImportError("alphagenome package not available")
    
    try:
        import pandas as pd
        from alphagenome.data import gene_annotation
        import time
    except ImportError as e:
        console.print(f"[red]Error importing required modules: {e}[/red]")
        return None
    
    api_key = config['alphagenome_api'].get('api_key') or os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError("AlphaGenome API key not provided")
    
    # Initialize client
    client = dna_client.create(api_key)
    
    console.print(f"[cyan]  Loading GTF for gene {gene_name}...[/cyan]")
    
    try:
        # Load GTF from local cache
        gtf_cache_path = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes/gtf_cache.feather")
        if not gtf_cache_path.exists():
            console.print(f"[red]GTF cache not found: {gtf_cache_path}[/red]")
            return None
        
        console.print(f"[dim]  Loading GTF from cache: {gtf_cache_path}[/dim]")
        gtf = pd.read_feather(gtf_cache_path)
        
        # Get gene interval
        interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_name)
        console.print(f"[dim]  ✓ Found gene interval: {interval}[/dim]")
        
        # Resize to desired window
        window_size_enum = getattr(dna_client, window_size_key)
        interval = interval.resize(window_size_enum)
        console.print(f"[dim]  ✓ Resized to {window_size_key}: {interval}[/dim]")
        
        # Get requested outputs
        requested_outputs = [dna_client.OutputType.RNA_SEQ]
        
        console.print(f"[cyan]  Calling AlphaGenome API (predict_interval)...[/cyan]")
        console.print(f"[dim]  Ontology terms: {ontology_terms}[/dim]")
        
        # Call AlphaGenome API using predict_interval
        start_time = time.time()
        output = client.predict_interval(
            interval=interval,
            requested_outputs=requested_outputs,
            ontology_terms=ontology_terms
        )
        elapsed = time.time() - start_time
        
        console.print(f"[cyan]  API call completed in {elapsed:.1f}s[/cyan]")
        
        # Extract values from output
        if hasattr(output, 'rna_seq') and hasattr(output.rna_seq, 'values'):
            values = output.rna_seq.values
            console.print(f"[dim]  ✓ Received predictions shape: {values.shape}[/dim]")
            
            # Apply rate limiting
            delay = config['alphagenome_api'].get('rate_limit_delay', 0.5)
            if delay > 0:
                time.sleep(delay)
            
            if return_full_output:
                return (values, output, gtf)
            return values
        else:
            console.print(f"[red]No RNA-seq data in output[/red]")
            return None
            
    except Exception as e:
        console.print(f"[red]Error in predict_with_alphagenome_interval: {e}[/red]")
        import traceback
        traceback.print_exc()
        return None


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
        
        # 1. Extrair sequência de referência do intervalo - usado nos comparison_mode: com alphagenome_ref   
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
        100000: "SEQUENCE_LENGTH_100KB",    # 100k = 131072 bp (AlphaGenome)
        524288: "SEQUENCE_LENGTH_500KB",    # 512k = 524288 bp (AlphaGenome)
        1048576: "SEQUENCE_LENGTH_1MB",     # 1MB exact
    }
    
    # Use window_center_size_key from config if available, otherwise infer from window_center_size
    if 'window_center_size_key' in config:
        window_size_key = config['window_center_size_key']
    else:
        window_size_key = window_size_map.get(window_center_size, "SEQUENCE_LENGTH_500KB")
    
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
            
            # Extract center if needed using centralized function
            central_track = extract_center_window(track_array, window_center_size, axis=0)
            
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
            131072: "SEQUENCE_LENGTH_100KB",  # AlphaGenome: 100KB = 131072 bp
            262144: "SEQUENCE_LENGTH_256KB",
            524288: "SEQUENCE_LENGTH_500KB",  # AlphaGenome: 500KB = 524288 bp
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
            
            # Extract center using centralized function
            central_track = extract_center_window(track_array, center_bp, axis=0)
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
            
            # Extrair trecho central usando função centralizada
            central_track = extract_center_window(track_array, window_center_size, axis=0)
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
    viewer: Optional[InteractiveViewer] = None,
    label1: str = "Cache",
    label2: str = "AlphaGenome"
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
        label1: Label for first dataset
        label2: Label for second dataset
        
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
                color='blue', linewidth=1.0, alpha=0.7, label=label1)
        ax.plot(x_positions, alpha_track,
                color='red', linewidth=1.0, linestyle='--', alpha=0.7, label=label2)
        
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


def process_sample_comparison_mode(
    config: Dict,
    index: int,
    dataset_metadata: Dict,
    splits_data: Optional[Dict] = None,
    norm_params: Optional[Dict] = None,
    viewer: Optional[InteractiveViewer] = None
) -> bool:
    """
    Process sample in comparison mode (3 modes available).
    
    Modes:
    - alphagenome_ref_x_dataset_dir: Compare AlphaGenome (reference) vs dataset_dir
    - alphagenome_ind_x_dataset_dir: Compare AlphaGenome (individual) vs dataset_dir
    - dataset_dir_x_cache_dir: Compare dataset_dir (normalized) vs cache_dir
    """
    comparison_mode = config.get('comparison_mode', 'dataset_dir_x_cache_dir')
    dataset_dir = Path(config['dataset_dir'])
    
    # Determine sample_id
    if config.get('sample_id'):
        sample_id = config['sample_id']
    else:
        # If comparison_mode is dataset_dir_x_cache_dir, use index from splits
        if comparison_mode == 'dataset_dir_x_cache_dir' and splits_data:
            split = config.get('split', 'train')
            # splits_data has keys like 'train_indices', 'val_indices', 'test_indices'
            split_key = f"{split}_indices"
            split_indices = splits_data.get(split_key, [])
            if index >= len(split_indices):
                console.print(f"[red]Index {index} out of range for split {split} (size: {len(split_indices)})[/red]")
                return False
            # Get the actual index in dataset_metadata['individuals']
            dataset_index = split_indices[index]
            sample_id = dataset_metadata['individuals'][dataset_index]
        else:
            sample_id = dataset_metadata['individuals'][index]
    
    console.print(f"\n[bold cyan]COMPARISON MODE: {comparison_mode}[/bold cyan]")
    console.print(f"[green]✓ Sample: {sample_id}[/green]")
    
    # Get gene filter
    gene_filter = config.get('gene_filter')
    genes_in_order = get_genes_in_dataset_order(dataset_dir, sample_id)
    
    if gene_filter is None:
        genes_to_process = genes_in_order
    elif isinstance(gene_filter, str):
        genes_to_process = [gene_filter]
    elif isinstance(gene_filter, list):
        genes_to_process = gene_filter
    else:
        raise ValueError(f"Invalid gene_filter: {gene_filter}")
    
    # Get window size
    window_size_key = config.get('raw_mode', {}).get('window_size_key', 'SEQUENCE_LENGTH_16KB')
    
    for gene_name in genes_to_process:
        console.print(f"\n[cyan]Processing {gene_name}...[/cyan]")
        
        try:
            if comparison_mode == "alphagenome_x_alphagenome_ref":
                # Mode: AlphaGenome (predict_interval) vs AlphaGenome (predict_sequence com FASTA)
                # Ambos devem usar o mesmo tamanho do dataset_dir, mas visualizar apenas window_size_key
                
                console.print(f"[dim]  Detecting dataset_dir window size...[/dim]")
                ontology_terms = config['alphagenome_api']['ontology_terms']
                
                # Primeiro, carregar dataset_dir SEM extrair janela, só para detectar tamanho
                data1_full, full_length = _load_dataset_dir_predictions(
                    dataset_dir, sample_id, gene_name, window_size_key=None
                )
                
                # Detectar o window_size_key usado no dataset_dir
                window_size_key_dataset = _detect_window_size_key_from_length(full_length)
                console.print(f"[cyan]  Dataset was generated with: {window_size_key_dataset} ({full_length} bp)[/cyan]")
                console.print(f"[cyan]  Will visualize center: {window_size_key}[/cyan]")
                
                # AlphaGenome usando predict_interval (sem FASTA) - COM O TAMANHO DO DATASET
                console.print(f"[dim]  Loading AlphaGenome via predict_interval ({window_size_key_dataset})...[/dim]")
                result_interval = predict_with_alphagenome_interval(
                    gene_name=gene_name,
                    config=config,
                    ontology_terms=ontology_terms,
                    window_size_key=window_size_key_dataset,  # USAR TAMANHO DO DATASET
                    return_full_output=False
                )
                
                if result_interval is None:
                    console.print(f"[red]  Failed to load AlphaGenome data (predict_interval)[/red]")
                    continue
                
                # AlphaGenome usando predict_sequence (com FASTA extraído) - COM O TAMANHO DO DATASET
                console.print(f"[dim]  Loading AlphaGenome via predict_sequence ({window_size_key_dataset})...[/dim]")
                result_fasta = predict_with_alphagenome(
                    gene_name=gene_name,
                    config=config,
                    ontology_terms=ontology_terms,
                    window_size_key=window_size_key_dataset,  # USAR TAMANHO DO DATASET
                    return_full_output=False,
                    sample_id=None  # Reference genome
                )
                
                if result_fasta is None:
                    console.print(f"[red]  Failed to load AlphaGenome data (predict_sequence)[/red]")
                    continue
                
                # result_interval: [full_length, 6]
                # result_fasta: [full_length, 6]
                
                # Transpor ambos para [6, full_length]
                data_interval_T = np.array(result_interval).T  # [6, full_length]
                data_fasta_T = np.array(result_fasta).T  # [6, full_length]
                
                # EXTRAIR JANELA CENTRAL APENAS PARA VISUALIZAÇÃO
                window_size_map = {
                    "SEQUENCE_LENGTH_2KB": 2048,
                    "SEQUENCE_LENGTH_16KB": 16384,
                    "SEQUENCE_LENGTH_100KB": 102400,
                    "SEQUENCE_LENGTH_500KB": 512000,
                    "SEQUENCE_LENGTH_1MB": 1048576
                }
                viz_length = window_size_map.get(window_size_key, 16384)
                
                console.print(f"[dim]  Extracting center window for visualization: {viz_length} bp[/dim]")
                data1 = extract_center_window(data_interval_T, viz_length, axis=1)
                data2 = extract_center_window(data_fasta_T, viz_length, axis=1)
                
                console.print(f"[dim]  data1 shape: {data1.shape}, data2 shape: {data2.shape}[/dim]")
                
                label1 = "AlphaGenome (predict_interval)"
                label2 = "AlphaGenome (predict_sequence/FASTA)"
                
            elif comparison_mode == "alphagenome_ref_x_dataset_dir":
                # Mode 1: AlphaGenome (reference) vs dataset_dir
                # USAR TAMANHO DO DATASET, VISUALIZAR window_size_key
                
                console.print(f"[dim]  Detecting dataset_dir window size...[/dim]")
                ontology_terms = config['alphagenome_api']['ontology_terms']
                
                # Carregar dataset_dir SEM extrair janela
                data_dataset_full, full_length = _load_dataset_dir_predictions(
                    dataset_dir, sample_id, gene_name, window_size_key=None
                )
                
                # Detectar o window_size_key usado no dataset_dir
                window_size_key_dataset = _detect_window_size_key_from_length(full_length)
                console.print(f"[cyan]  Dataset was generated with: {window_size_key_dataset} ({full_length} bp)[/cyan]")
                console.print(f"[cyan]  Will visualize center: {window_size_key}[/cyan]")
                
                # AlphaGenome (reference) - COM O TAMANHO DO DATASET
                console.print(f"[dim]  Loading AlphaGenome reference ({window_size_key_dataset})...[/dim]")
                result = predict_with_alphagenome(
                    gene_name=gene_name,
                    config=config,
                    ontology_terms=ontology_terms,
                    window_size_key=window_size_key_dataset,  # USAR TAMANHO DO DATASET
                    return_full_output=False,
                    sample_id=None  # Use reference genome
                )
                
                if result is None:
                    console.print(f"[red]  Failed to load AlphaGenome data[/red]")
                    continue
                
                # result: [full_length, 6]
                # data_dataset_full: [full_length, 6]
                alphagenome_data = np.array(result)
                
                # Transpor ambos para [6, full_length]
                alphagenome_data_T = alphagenome_data.T
                data_dataset_full_T = data_dataset_full.T
                
                # EXTRAIR JANELA CENTRAL APENAS PARA VISUALIZAÇÃO
                window_size_map = {
                    "SEQUENCE_LENGTH_2KB": 2048,
                    "SEQUENCE_LENGTH_16KB": 16384,
                    "SEQUENCE_LENGTH_100KB": 102400,
                    "SEQUENCE_LENGTH_500KB": 512000,
                    "SEQUENCE_LENGTH_1MB": 1048576
                }
                viz_length = window_size_map.get(window_size_key, 16384)
                
                console.print(f"[dim]  Extracting center window for visualization: {viz_length} bp[/dim]")
                data1 = extract_center_window(alphagenome_data_T, viz_length, axis=1)
                data2 = extract_center_window(data_dataset_full_T, viz_length, axis=1)
                
                console.print(f"[dim]  AlphaGenome shape: {data1.shape}, Dataset shape: {data2.shape}[/dim]")
                
                label1 = "AlphaGenome (Ref)"
                label2 = "Dataset Dir"
                
            elif comparison_mode == "alphagenome_ind_x_dataset_dir":
                # Mode 2: AlphaGenome (individual) vs dataset_dir
                # USAR TAMANHO DO DATASET, VISUALIZAR window_size_key
                
                console.print(f"[dim]  Detecting dataset_dir window size...[/dim]")
                ontology_terms = config['alphagenome_api']['ontology_terms']
                
                # Carregar dataset_dir SEM extrair janela
                data_dataset_full, full_length = _load_dataset_dir_predictions(
                    dataset_dir, sample_id, gene_name, window_size_key=None
                )
                
                # Detectar o window_size_key usado no dataset_dir
                window_size_key_dataset = _detect_window_size_key_from_length(full_length)
                console.print(f"[cyan]  Dataset was generated with: {window_size_key_dataset} ({full_length} bp)[/cyan]")
                console.print(f"[cyan]  Will visualize center: {window_size_key}[/cyan]")
                
                # AlphaGenome (individual) - COM O TAMANHO DO DATASET
                console.print(f"[dim]  Loading AlphaGenome individual ({window_size_key_dataset})...[/dim]")
                result = predict_with_alphagenome(
                    gene_name=gene_name,
                    config=config,
                    ontology_terms=ontology_terms,
                    window_size_key=window_size_key_dataset,  # USAR TAMANHO DO DATASET
                    return_full_output=False,
                    sample_id=sample_id  # Use individual's genome
                )
                
                if result is None:
                    console.print(f"[red]  Failed to load AlphaGenome data[/red]")
                    continue
                
                # result: [full_length, 6]
                # data_dataset_full: [full_length, 6]
                alphagenome_data = np.array(result)
                
                # Transpor ambos para [6, full_length]
                alphagenome_data_T = alphagenome_data.T
                data_dataset_full_T = data_dataset_full.T
                
                # EXTRAIR JANELA CENTRAL APENAS PARA VISUALIZAÇÃO
                window_size_map = {
                    "SEQUENCE_LENGTH_2KB": 2048,
                    "SEQUENCE_LENGTH_16KB": 16384,
                    "SEQUENCE_LENGTH_100KB": 102400,
                    "SEQUENCE_LENGTH_500KB": 512000,
                    "SEQUENCE_LENGTH_1MB": 1048576
                }
                viz_length = window_size_map.get(window_size_key, 16384)
                
                console.print(f"[dim]  Extracting center window for visualization: {viz_length} bp[/dim]")
                data1 = extract_center_window(alphagenome_data_T, viz_length, axis=1)
                data2 = extract_center_window(data_dataset_full_T, viz_length, axis=1)
                
                console.print(f"[dim]  AlphaGenome shape: {data1.shape}, Dataset shape: {data2.shape}[/dim]")
                
                label1 = f"AlphaGenome ({sample_id})"
                label2 = "Dataset Dir"
                
            elif comparison_mode == "dataset_dir_x_cache_dir":
                # Mode 3: dataset_dir (normalized) vs cache_dir
                console.print(f"[dim]  Loading dataset_dir data...[/dim]")
                
                # Load dataset_dir data
                dataset_data, _ = _load_dataset_dir_predictions(
                    dataset_dir, sample_id, gene_name, window_size_key
                )
                
                # Apply normalization
                console.print(f"[dim]  Applying normalization...[/dim]")
                norm_result = apply_normalization(
                    dataset_data,  # Already numpy array
                    norm_params['method'],
                    norm_params
                )
                # Convert to numpy if it's a tensor
                if torch.is_tensor(norm_result):
                    dataset_data_normalized = norm_result.numpy()
                else:
                    dataset_data_normalized = norm_result
                
                # Load cache data
                console.print(f"[dim]  Loading cache data...[/dim]")
                cache_dir = Path(config['cache_dir'])
                split = config['split']
                
                # Find index in split
                if splits_data is None:
                    splits_file = cache_dir / 'splits.json'
                    with open(splits_file) as f:
                        splits_data = json.load(f)
                
                # splits_data has keys like 'train_indices', not 'train'
                # We need to find which split contains the sample and its position
                # The sample_id was determined from dataset_metadata['individuals'][dataset_index]
                # So we need to find dataset_index in one of the split_indices arrays
                split_found = None
                cache_index = None
                
                # Get dataset_index from sample_id
                dataset_index = dataset_metadata['individuals'].index(sample_id)
                
                for split_name in ['train', 'val', 'test']:
                    split_key = f"{split_name}_indices"
                    if split_key in splits_data:
                        split_indices = splits_data[split_key]
                        try:
                            cache_index = split_indices.index(dataset_index)
                            split_found = split_name
                            break
                        except ValueError:
                            continue
                
                if split_found is None:
                    console.print(f"[red]  Sample {sample_id} (index {dataset_index}) not found in any split[/red]")
                    continue
                
                console.print(f"[dim]  Sample found in {split_found} split at cache index {cache_index}[/dim]")
                split = split_found
                
                # Load cache data (returns all genes, need to filter)
                cache_features_all, cache_metadata, _ = load_cache_data(
                    cache_dir, split, cache_index, gene_filter=None  # Load all genes first
                )
                
                console.print(f"[dim]  cache_features_all shape: {cache_features_all.shape}[/dim]")
                
                # Find the gene index in the cache
                # If genes_in_order is not in metadata, get it from dataset_dir
                if 'genes_in_order' in cache_metadata:
                    genes_in_cache = cache_metadata['genes_in_order']
                else:
                    genes_in_cache = get_genes_in_dataset_order(dataset_dir, sample_id)
                
                console.print(f"[dim]  Genes in cache: {genes_in_cache}[/dim]")
                
                try:
                    gene_idx_in_cache = genes_in_cache.index(gene_name)
                except ValueError:
                    console.print(f"[red]  Gene {gene_name} not found in cache[/red]")
                    continue
                
                # Extract the 6 tracks for this gene
                start_track = gene_idx_in_cache * 6
                end_track = start_track + 6
                cache_features = cache_features_all[start_track:end_track, :]
                
                console.print(f"[dim]  Extracted cache tracks [{start_track}:{end_track}], shape: {cache_features.shape}[/dim]")
                
                # Transpose both to [6, window_size] for plotting
                dataset_data_normalized_T = dataset_data_normalized.T  # [6, window_size]
                
                console.print(f"[dim]  dataset_data_normalized_T type: {type(dataset_data_normalized_T)}, cache_features type: {type(cache_features)}[/dim]")
                
                # Make sure cache_features is a numpy array
                if torch.is_tensor(cache_features):
                    cache_features_np = cache_features.numpy()
                else:
                    cache_features_np = cache_features
                
                # Extract central region from dataset_data to match cache window_size
                dataset_window_size = dataset_data_normalized_T.shape[1]  # e.g., 102400
                cache_window_size = cache_features_np.shape[1]  # e.g., 10000
                
                console.print(f"[dim]  Dataset window: {dataset_window_size}, Cache window: {cache_window_size}[/dim]")
                
                if dataset_window_size > cache_window_size:
                    # Extract center from dataset_data using centralized function
                    dataset_data_trimmed = extract_center_window(dataset_data_normalized_T, cache_window_size, axis=1)
                    console.print(f"[dim]  Trimmed dataset_data to cache window size: {dataset_data_trimmed.shape}[/dim]")
                else:
                    dataset_data_trimmed = dataset_data_normalized_T
                
                data1 = dataset_data_trimmed  # [6, cache_window_size]
                data2 = cache_features_np  # [6, cache_window_size]
                
                label1 = "Dataset Dir (normalized)"
                label2 = "Cache Dir"
                
            else:
                console.print(f"[red]  Invalid comparison_mode: {comparison_mode}[/red]")
                continue
            
            # Compute metrics first
            console.print(f"[dim]  data1 shape: {data1.shape}, data2 shape: {data2.shape}[/dim]")
            metrics = compute_metrics(data1, data2, config.get('verbose_metrics', True))
            
            # Plot comparison
            fig = plot_comparison(
                data1, data2, sample_id, metrics, [gene_name], config, viewer,
                label1=label1, label2=label2
            )
            
            # Save if requested
            if config.get('save_plots', False):
                output_dir = config.get('output_dir')
                if output_dir:
                    output_path = Path(output_dir)
                    output_path.mkdir(parents=True, exist_ok=True)
                    filename = f"{config.get('output_prefix', 'compare')}_{sample_id}_{gene_name}.png"
                    fig.savefig(output_path / filename, dpi=150, bbox_inches='tight')
                    console.print(f"[green]✓ Saved: {output_path / filename}[/green]")
            
            plt.show()
            
        except Exception as e:
            console.print(f"[red]  Error processing {gene_name}: {e}[/red]")
            import traceback
            traceback.print_exc()
            continue
    
    console.print(f"\n[green]✓ Comparison mode processing completed![/green]")
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
    
    # Check for comparison mode first
    comparison_mode = config.get('comparison_mode')
    if comparison_mode:
        # Validate comparison_mode
        valid_modes = [
            'alphagenome_x_alphagenome_ref',
            'alphagenome_ref_x_dataset_dir', 
            'alphagenome_ind_x_dataset_dir', 
            'dataset_dir_x_cache_dir'
        ]
        if comparison_mode not in valid_modes:
            console.print(f"[red]Invalid comparison_mode: {comparison_mode}[/red]")
            console.print(f"[yellow]Valid modes: {', '.join(valid_modes)}[/yellow]")
            sys.exit(1)
        
        # Validate required configs for each mode
        if 'alphagenome' in comparison_mode:
            if not config.get('alphagenome_api', {}).get('enabled', False):
                console.print(f"[red]AlphaGenome API must be enabled for mode: {comparison_mode}[/red]")
                sys.exit(1)
        
        if comparison_mode == 'dataset_dir_x_cache_dir':
            cache_dir = Path(config.get('cache_dir', ''))
            if not cache_dir or not cache_dir.exists():
                console.print(f"[red]cache_dir must exist for mode: {comparison_mode}[/red]")
                console.print(f"[yellow]cache_dir: {cache_dir}[/yellow]")
                sys.exit(1)
        
        # Load dataset metadata
        dataset_dir = Path(config['dataset_dir'])
        dataset_metadata_file = dataset_dir / 'dataset_metadata.json'
        with open(dataset_metadata_file) as f:
            dataset_metadata = json.load(f)
        
        # Load splits and norm_params if needed
        splits_data = None
        norm_params = None
        if comparison_mode == 'dataset_dir_x_cache_dir':
            cache_dir = Path(config['cache_dir'])
            splits_file = cache_dir / 'splits.json'
            with open(splits_file) as f:
                splits_data = json.load(f)
            
            norm_params_file = cache_dir / 'normalization_params.json'
            with open(norm_params_file) as f:
                norm_params = json.load(f)
        
        # Process sample in comparison mode
        index = config.get('index', 0)
        success = process_sample_comparison_mode(
            config, index, dataset_metadata, 
            splits_data=splits_data,
            norm_params=norm_params
        )
        
        if success:
            console.print("\n[bold green]✓ Comparison completed successfully![/bold green]\n")
        else:
            sys.exit(1)
        
        return
    
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
