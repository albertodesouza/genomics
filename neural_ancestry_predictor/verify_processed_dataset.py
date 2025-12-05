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
import shutil
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

# ═══════════════════════════════════════════════════════════════════════════
# Constantes Globais - Tamanhos de Janela AlphaGenome
# ═══════════════════════════════════════════════════════════════════════════

# Tamanhos de janela suportados pelo AlphaGenome (em pares de base)
WINDOW_SIZE_4KB = 4096
WINDOW_SIZE_8KB = 8192
WINDOW_SIZE_16KB = 16384
WINDOW_SIZE_32KB = 32768
WINDOW_SIZE_64KB = 65536
WINDOW_SIZE_128KB = 131072
WINDOW_SIZE_256KB = 262144
WINDOW_SIZE_512KB = 524288   # AlphaGenome: "500KB" = 524288 bp (512 KB real)
WINDOW_SIZE_1MB = 1048576

# Mapeamento de nomes de sequência para tamanhos (usado para validação)
SEQUENCE_LENGTH_MAP = {
    "SEQUENCE_LENGTH_4KB": WINDOW_SIZE_4KB,
    "SEQUENCE_LENGTH_8KB": WINDOW_SIZE_8KB,
    "SEQUENCE_LENGTH_16KB": WINDOW_SIZE_16KB,
    "SEQUENCE_LENGTH_32KB": WINDOW_SIZE_32KB,
    "SEQUENCE_LENGTH_64KB": WINDOW_SIZE_64KB,
    "SEQUENCE_LENGTH_128KB": WINDOW_SIZE_128KB,
    "SEQUENCE_LENGTH_256KB": WINDOW_SIZE_256KB,
    "SEQUENCE_LENGTH_500KB": WINDOW_SIZE_512KB,
    "SEQUENCE_LENGTH_512KB": WINDOW_SIZE_512KB,  # Alias
    "SEQUENCE_LENGTH_1MB": WINDOW_SIZE_1MB,
}

# Mapeamento reverso: tamanho -> nome de sequência
SIZE_TO_SEQUENCE_NAME = {
    WINDOW_SIZE_4KB: "SEQUENCE_LENGTH_4KB",
    WINDOW_SIZE_8KB: "SEQUENCE_LENGTH_8KB",
    WINDOW_SIZE_16KB: "SEQUENCE_LENGTH_16KB",
    WINDOW_SIZE_32KB: "SEQUENCE_LENGTH_32KB",
    WINDOW_SIZE_64KB: "SEQUENCE_LENGTH_64KB",
    WINDOW_SIZE_128KB: "SEQUENCE_LENGTH_128KB",
    WINDOW_SIZE_256KB: "SEQUENCE_LENGTH_256KB",
    WINDOW_SIZE_512KB: "SEQUENCE_LENGTH_500KB",
    WINDOW_SIZE_1MB: "SEQUENCE_LENGTH_1MB",
}

NUM_ONTOLOGIES = 6  # Número de ontologias RNA-seq por gene (valor padrão)

# ═══════════════════════════════════════════════════════════════════════════
# Helper Function: Generate Ontology Labels from Dataset Metadata
# ═══════════════════════════════════════════════════════════════════════════

def get_ontology_terms_from_metadata(dataset_metadata: Dict) -> List[str]:
    """
    Extrai lista de ontology terms (CURIEs) dos metadados do dataset.
    
    Args:
        dataset_metadata: Dicionário com metadados do dataset
        
    Returns:
        Lista de ontology CURIEs ordenados (ex: ["CL:0000346", "CL:1000458", ...])
    """
    ontology_details = dataset_metadata.get('ontology_details', {})
    
    if ontology_details:
        # Retornar CURIEs ordenados das ontologies
        return sorted(ontology_details.keys())
    
    # Fallback: tentar ler do campo 'ontologies' (pode ter espaços indevidos, então limpar)
    ontologies = dataset_metadata.get('ontologies', [])
    if ontologies:
        return sorted([ont.strip() for ont in ontologies if ont.strip()])
    
    # Fallback final: retornar lista padrão se nada encontrado
    console.print("[yellow]⚠ Ontologias não encontradas nos metadados, usando lista padrão[/yellow]")
    return ["CL:1000458", "CL:0000346", "CL:2000092"]


def generate_ontology_labels(dataset_metadata: Dict) -> List[str]:
    """
    Gera labels de ontologias dinamicamente a partir dos metadados do dataset.
    
    Args:
        dataset_metadata: Dicionário com metadados do dataset (deve conter 'ontology_details')
        
    Returns:
        Lista de labels formatados para visualização (com strand + e -)
        
    Example:
        >>> labels = generate_ontology_labels(metadata)
        >>> labels[0]
        'CL:1000458 (+)\\nmelanocyte of skin'
    """
    ontology_details = dataset_metadata.get('ontology_details', {})
    
    if not ontology_details:
        # Fallback para labels padrão se metadados não tiverem ontology_details
        console.print("[yellow]⚠ Campo 'ontology_details' não encontrado nos metadados, usando labels padrão[/yellow]")
        return [
            "CL:1000458 (+)\nMelanocyte",
            "CL:0000346 (+)\nDermal Papilla",
            "CL:2000092 (+)\nKeratinocyte",
            "CL:1000458 (-)\nMelanocyte",
            "CL:0000346 (-)\nDermal Papilla",
            "CL:2000092 (-)\nKeratinocyte"
        ]
    
    # Ordenar ontologias para ordem consistente
    sorted_ontologies = sorted(ontology_details.keys())
    
    labels = []
    # Primeiro todas as ontologias com strand +
    for ontology_curie in sorted_ontologies:
        details = ontology_details[ontology_curie]
        biosample_name = details.get('biosample_name', ontology_curie)
        # Simplificar nome se muito longo
        if len(biosample_name) > 30:
            # Pegar apenas a primeira parte antes de vírgula ou "of"
            biosample_name = biosample_name.split(',')[0].split(' of ')[0]
        label = f"{ontology_curie} (+)\n{biosample_name}"
        labels.append(label)
    
    # Depois todas as ontologias com strand -
    for ontology_curie in sorted_ontologies:
        details = ontology_details[ontology_curie]
        biosample_name = details.get('biosample_name', ontology_curie)
        if len(biosample_name) > 30:
            biosample_name = biosample_name.split(',')[0].split(' of ')[0]
        label = f"{ontology_curie} (-)\n{biosample_name}"
        labels.append(label)
    
    return labels


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


def _get_alphagenome_predictions_for_individual(
    sample_id: str,
    gene_name: str,
    window_size_key: str,
    config: Dict,
    ontology_terms: List[str]
) -> Optional[np.ndarray]:
    """
    Usa build_window_and_predict.py como biblioteca para gerar predições
    AlphaGenome corretas para o indivíduo.
    
    Args:
        sample_id: ID do indivíduo (e.g., HG02445)
        gene_name: Nome do gene (e.g., MC1R)
        window_size_key: Tamanho da janela (e.g., SEQUENCE_LENGTH_500KB)
        config: Configuration dict
        ontology_terms: Lista de ontology terms
    
    Returns:
        Array [window_size, 6] com predições ou None se falhar
    """
    try:
        import argparse
        from alphagenome.data import gene_annotation
        
        # Import build_window_and_predict dentro da função para evitar erros de import
        sys.path.insert(0, str(Path(__file__).parent.parent / "build_non_longevous_dataset"))
        import build_window_and_predict
        
        console.print(f"[cyan]  Using build_window_and_predict.py to generate AlphaGenome data...[/cyan]")
        
        # Criar diretório temporário
        tmp_outdir = Path("/tmp/GENOMICS_DATA/top3")
        tmp_outdir.mkdir(parents=True, exist_ok=True)
        
        # Inferir paths a partir do dataset_dir
        # Assumindo estrutura padrão do GENOMICS_DATA
        dataset_dir = Path(config.get('dataset_dir', '/dados/GENOMICS_DATA/top3/non_longevous_results_genes'))
        genomics_data_root = dataset_dir.parent  # /dados/GENOMICS_DATA/top3
        
        ref_fasta = genomics_data_root / "refs" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
        vcf_pattern = str(genomics_data_root / "longevity_dataset" / "vcf_chromosomes" / "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
        
        if not ref_fasta.exists():
            console.print(f"[red]  Reference FASTA not found: {ref_fasta}[/red]")
            console.print(f"[yellow]  Inferred from dataset_dir: {dataset_dir}[/yellow]")
            console.print(f"[yellow]  Tried: {ref_fasta}[/yellow]")
            return None
        
        # Carregar GTF para obter intervalo do gene
        gtf_cache_path = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes/gtf_cache.feather")
        if not gtf_cache_path.exists():
            console.print(f"[red]  GTF cache not found: {gtf_cache_path}[/red]")
            return None
        
        import pandas as pd
        gtf = pd.read_feather(gtf_cache_path)
        interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_name)
        
        # Mapear window_size_key para número de bases usando constantes globais
        window_size_map = SEQUENCE_LENGTH_MAP
        window_size = window_size_map.get(window_size_key)
        if window_size is None:
            console.print(f"[red]  Invalid window_size_key: {window_size_key}[/red]")
            return None
        
        # Resize interval
        interval = interval.resize(window_size)
        
        # Detectar prefixo do cromossomo
        ref_fai = Path(str(ref_fasta) + ".fai")
        prefix = build_window_and_predict.detect_chr_prefix(ref_fai)
        chrom = build_window_and_predict.coerce_chromosome_name(interval.chromosome, prefix)
        
        # Resolver caminho do VCF para este cromossomo
        vcf_path = Path(vcf_pattern.replace('{chrom}', chrom))
        if not vcf_path.exists():
            console.print(f"[red]  VCF not found: {vcf_path}[/red]")
            return None
        
        # Criar objeto args mock para process_window
        args = argparse.Namespace(
            predict=True,
            api_key=config['alphagenome_api'].get('api_key') or os.environ.get('ALPHAGENOME_API_KEY'),
            outputs="RNA_SEQ",
            ontology=",".join(ontology_terms),
            all_tissues=False,
            skip_h2=True,
            also_iupac=False,
            reference_only=False,
            window_size=window_size,
            api_rate_limit_delay=config['alphagenome_api'].get('rate_limit_delay', 0.5)
        )
        
        console.print(f"[dim]  Calling process_window for {sample_id} / {gene_name}...[/dim]")
        console.print(f"[dim]  Interval: {chrom}:{interval.start}-{interval.end} ({window_size} bp)[/dim]")
        console.print(f"[dim]  Temporary output: {tmp_outdir}[/dim]")
        
        # Chamar process_window
        case_dir = build_window_and_predict.process_window(
            sample=sample_id,
            target_name=gene_name,
            chrom=chrom,
            start=interval.start,
            end=interval.end,
            ref_fa=ref_fasta,
            vcf_path=vcf_path,
            outdir=tmp_outdir,
            window_size=window_size,
            args=args
        )
        
        # Carregar resultado do .npz
        predictions_h1_dir = case_dir / "predictions_H1"
        rna_seq_file = predictions_h1_dir / "rna_seq.npz"
        
        if not rna_seq_file.exists():
            console.print(f"[red]  Predictions file not found: {rna_seq_file}[/red]")
            return None
        
        console.print(f"[green]  ✓ Loading predictions from: {rna_seq_file}[/green]")
        data = np.load(rna_seq_file)
        values = data['values']  # Shape: [window_size, 6]
        
        console.print(f"[green]  ✓ Loaded predictions shape: {values.shape}[/green]")
        
        return values
        
    except Exception as e:
        console.print(f"[red]  Error in _get_alphagenome_predictions_for_individual: {e}[/red]")
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
    length_to_key = SIZE_TO_SEQUENCE_NAME
    
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
    # Usar mapeamento global de constantes
    target_length = SEQUENCE_LENGTH_MAP.get(window_size_key)
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
            splits_data: Dados dos splits (train/val/test indices ou metadados completos)
            dataset_metadata: Metadata do dataset (lista de indivíduos) - pode ser None para formato v2
        """
        self.config = config
        self.splits_data = splits_data
        self.dataset_metadata = dataset_metadata
        self.current_index = config.get('index', 0)
        self.current_split = config.get('split', 'test')
        self.should_exit = False
        self.fig = None
        
        # Detectar formato do splits_data
        self.format_version = splits_data.get('format_version', 1)
        
        # Obter tamanho do split atual
        self.split_size = get_split_sample_count(splits_data, self.current_split, self.format_version)
        self.max_index = self.split_size - 1
        
        console.print(f"[cyan]Modo interativo ativado:[/cyan]")
        console.print(f"  • Split: {self.current_split}")
        console.print(f"  • Total de amostras: {self.split_size}")
        console.print(f"  • Índice inicial: {self.current_index}")
        console.print(f"  • Formato splits_metadata.json: v{self.format_version}")
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
        sample_meta = get_sample_info(
            self.splits_data, self.current_split, self.current_index,
            self.format_version, self.dataset_metadata
        )
        return sample_meta.get('sample_id', f"#{self.current_index + 1}")
    
    def get_current_sample_metadata(self) -> Dict:
        """Retorna os metadados completos da amostra atual."""
        return get_sample_info(
            self.splits_data, self.current_split, self.current_index,
            self.format_version, self.dataset_metadata
        )


class InteractiveComparisonViewer:
    """Gerencia navegação interativa de comparação entre dois indivíduos."""
    
    def __init__(self, config: Dict, splits_data: Dict, dataset_metadata: Dict, 
                 metadata_csv_path: str):
        self.config = config
        self.splits_data = splits_data
        self.dataset_metadata = dataset_metadata
        self.current_index_1 = config.get('index', 0)  # Primeiro indivíduo
        self.current_index_2 = config.get('index', 0)  # Segundo indivíduo (começa igual)
        self.current_gene_index = 0  # Índice do gene atual
        self.current_split = config.get('split', 'test')
        self.should_exit = False
        self.fig = None
        
        # Detectar formato do splits_data
        self.format_version = splits_data.get('format_version', 1)
        
        # Ler lista de genes disponíveis dos metadados do dataset
        self.available_genes = dataset_metadata.get('genes', []) if dataset_metadata else []
        if not self.available_genes:
            # Erro: campo 'genes' não encontrado nos metadados
            console.print("[red]❌ ERRO: Campo 'genes' não encontrado nos metadados do dataset![/red]")
            console.print("[yellow]Os metadados do dataset precisam ser regenerados com a versão atualizada do dataset_builder.py[/yellow]")
            console.print("[cyan]Execute: python3 build_non_longevous_dataset.py --config <seu_config.yaml>[/cyan]")
            console.print("[cyan]Com generate_dataset_metadata: true no config[/cyan]")
            raise ValueError(
                "Campo 'genes' não encontrado nos metadados do dataset. "
                "Regenere os metadados do dataset com a versão atualizada do dataset builder"
            )
        
        # Carregar metadata de superpopulação
        self.superpopulation_map = self._load_superpopulation_map(metadata_csv_path)
        
        # Obter tamanho do split atual
        self.split_size = get_split_sample_count(splits_data, self.current_split, self.format_version)
        self.max_index = self.split_size - 1
        self.max_gene_index = len(self.available_genes) - 1
        
    def _load_superpopulation_map(self, csv_path: str) -> Dict[str, Tuple[str, str]]:
        """Carrega mapeamento sample_id -> (population, superpopulation)."""
        import pandas as pd
        df = pd.read_csv(csv_path)
        return {
            row['SampleID']: (row['Population'], row['Superpopulation'])
            for _, row in df.iterrows()
        }
    
    def on_key_press(self, event: KeyEvent):
        """Handler para eventos de teclado."""
        if event.key == 'right':
            # Avançar ambos os indivíduos
            if self.current_index_1 < self.max_index:
                self.current_index_1 += 1
                self.current_index_2 += 1
                if self.current_index_2 > self.max_index:
                    self.current_index_2 = self.max_index
                console.print(f"[cyan]→ Avançando ambos: Ind1={self.current_index_1}, Ind2={self.current_index_2}[/cyan]")
                plt.close(self.fig)
        
        elif event.key == 'left':
            # Retroceder ambos os indivíduos
            if self.current_index_1 > 0:
                self.current_index_1 -= 1
                self.current_index_2 = max(0, self.current_index_2 - 1)
                console.print(f"[cyan]← Retrocedendo ambos: Ind1={self.current_index_1}, Ind2={self.current_index_2}[/cyan]")
                plt.close(self.fig)
        
        elif event.key == 'd':
            # Avançar apenas o segundo indivíduo
            if self.current_index_2 < self.max_index:
                self.current_index_2 += 1
                console.print(f"[cyan]→ Avançando Ind2: {self.current_index_2}[/cyan]")
                plt.close(self.fig)
        
        elif event.key == 'a':
            # Retroceder apenas o segundo indivíduo
            if self.current_index_2 > 0:
                self.current_index_2 -= 1
                console.print(f"[cyan]← Retrocedendo Ind2: {self.current_index_2}[/cyan]")
                plt.close(self.fig)
        
        elif event.key == 'w':
            # Avançar gene
            if self.current_gene_index < self.max_gene_index:
                self.current_gene_index += 1
                console.print(f"[cyan]↑ Gene: {self.available_genes[self.current_gene_index]}[/cyan]")
                plt.close(self.fig)
        
        elif event.key == 'z':
            # Retroceder gene
            if self.current_gene_index > 0:
                self.current_gene_index -= 1
                console.print(f"[cyan]↓ Gene: {self.available_genes[self.current_gene_index]}[/cyan]")
                plt.close(self.fig)
        
        elif event.key == 'q':
            # Sair
            console.print(f"[yellow]Saindo...[/yellow]")
            self.should_exit = True
            plt.close(self.fig)
    
    def get_current_gene(self) -> str:
        """Retorna o gene atual."""
        return self.available_genes[self.current_gene_index]
    
    def get_sample_info(self, index: int) -> Tuple[str, str, str]:
        """Retorna (sample_id, population, superpopulation) para um índice."""
        sample_meta = get_sample_info(
            self.splits_data, self.current_split, index,
            self.format_version, self.dataset_metadata
        )
        sample_id = sample_meta.get('sample_id', f"#{index + 1}")
        
        # Tentar obter do mapa de superpopulação primeiro (mais confiável para formato legado)
        pop, superpop = self.superpopulation_map.get(sample_id, (None, None))
        
        # Se não encontrou no mapa, usar os metadados do splits_metadata.json (formato v2)
        if pop is None:
            pop = sample_meta.get('population', 'Unknown')
        if superpop is None:
            superpop = sample_meta.get('superpopulation', 'Unknown')
        
        return sample_id, pop, superpop


# ==============================================================================
# FUNÇÕES DE CARREGAMENTO DE DADOS
# ==============================================================================

def load_splits_metadata(cache_dir: Path) -> Tuple[Dict, int]:
    """
    Carrega dados de splits_metadata.json detectando automaticamente o formato.
    
    Args:
        cache_dir: Diretório do cache
        
    Returns:
        Tupla (splits_data, format_version)
    """
    splits_file = cache_dir / 'splits_metadata.json'
    if not splits_file.exists():
        raise FileNotFoundError(f"Splits não encontrado: {splits_file}")
    
    with open(splits_file, 'r') as f:
        splits = json.load(f)
    
    format_version = splits.get('format_version', 1)
    return splits, format_version


def get_sample_info(
    splits: Dict,
    split: str,
    index: int,
    format_version: int,
    dataset_metadata: Optional[Dict] = None
) -> Dict:
    """
    Obtém metadados de um sample a partir de splits_metadata.json.
    
    Args:
        splits: Dados carregados de splits_metadata.json
        split: Nome do split (train/val/test)
        index: Índice no split
        format_version: Versão do formato (1=legado, 2=novo)
        dataset_metadata: Metadados do dataset fonte (necessário para formato v1)
        
    Returns:
        Dict com sample_id, superpopulation, population, sex
    """
    if format_version >= 2:
        # Novo formato: metadados completos diretamente em splits_metadata.json
        split_data = splits.get(split, [])
        if index >= len(split_data):
            return {"sample_id": f"#{index + 1}", "superpopulation": "UNK", "population": "UNK", "sex": 0}
        return split_data[index]
    else:
        # Formato legado: precisa de dataset_metadata
        split_key = f"{split}_indices"
        split_indices = splits.get(split_key, [])
        if index >= len(split_indices):
            return {"sample_id": f"#{index + 1}", "superpopulation": "UNK", "population": "UNK", "sex": 0}
        
        global_index = split_indices[index]
        
        if dataset_metadata is None:
            return {"sample_id": f"sample_{global_index}", "superpopulation": "UNK", "population": "UNK", "sex": 0}
        
        individuals = dataset_metadata.get('individuals', [])
        individuals_pedigree = dataset_metadata.get('individuals_pedigree', {})
        
        if global_index >= len(individuals):
            return {"sample_id": f"sample_{global_index}", "superpopulation": "UNK", "population": "UNK", "sex": 0}
        
        sample_id = individuals[global_index]
        pedigree = individuals_pedigree.get(sample_id, {})
        
        return {
            "sample_id": sample_id,
            "superpopulation": pedigree.get("superpopulation", "UNK"),
            "population": pedigree.get("population", "UNK"),
            "sex": pedigree.get("sex", 0)
        }


def get_split_sample_count(splits: Dict, split: str, format_version: int) -> int:
    """
    Obtém o tamanho de um split.
    
    Args:
        splits: Dados carregados de splits_metadata.json
        split: Nome do split (train/val/test)
        format_version: Versão do formato
        
    Returns:
        Número de samples no split
    """
    if format_version >= 2:
        return len(splits.get(split, []))
    else:
        split_key = f"{split}_indices"
        return len(splits.get(split_key, []))


def load_cache_data(
    cache_dir: Path,
    split: str,
    index: int,
    gene_filter: Optional[Union[str, List[str]]] = None
) -> Tuple[np.ndarray, Dict, Optional[Dict]]:
    """
    Carrega dados do cache processado.
    
    Args:
        cache_dir: Diretório do cache
        split: Conjunto de dados (train/val/test)
        index: Índice do indivíduo no split
        gene_filter: Gene(s) a filtrar. None = todos, string = um gene, lista = múltiplos genes
        
    Returns:
        Tupla (features_array, metadata, sample_metadata)
        - features_array: Array numpy shape [num_genes*6, window_center_size]
        - metadata: Dicionário com metadados do cache
        - sample_metadata: Dicionário com metadados do sample (sample_id, superpopulation, etc)
    """
    cache_dir = Path(cache_dir)
    
    # Carregar metadados do cache
    metadata_file = cache_dir / 'metadata.json'
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata não encontrado: {metadata_file}")
    
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    
    # Carregar splits
    splits, format_version = load_splits_metadata(cache_dir)
    
    # Verificar tamanho do split
    split_size = get_split_sample_count(splits, split, format_version)
    if split_size == 0:
        raise ValueError(f"Split '{split}' não encontrado ou vazio. Opções: train, val, test")
    
    if index >= split_size:
        raise ValueError(f"Índice {index} fora do range para split '{split}' (max: {split_size - 1})")
    
    # Carregar dados do split
    data_file = cache_dir / f'{split}_data.pt'
    if not data_file.exists():
        raise FileNotFoundError(f"Arquivo de dados não encontrado: {data_file}")
    
    data = torch.load(data_file, map_location='cpu')
    
    # Extrair features e target do indivíduo
    features, target = data[index]
    features_array = features.numpy()  # Shape: [66, window_center_size]
    
    # Obter metadados do sample
    # Para formato legado, precisamos carregar dataset_metadata
    dataset_metadata = None
    if format_version < 2:
        dataset_dir_str = metadata.get('dataset_dir', '')
        if dataset_dir_str:
            dataset_metadata_file = Path(dataset_dir_str) / 'dataset_metadata.json'
            if dataset_metadata_file.exists():
                with open(dataset_metadata_file, 'r') as f:
                    dataset_metadata = json.load(f)
    
    sample_metadata = get_sample_info(splits, split, index, format_version, dataset_metadata)
    
    return features_array, metadata, sample_metadata


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
        # Carregar dataset_metadata para passar à função
        metadata_file = dataset_dir / 'dataset_metadata.json'
        dataset_metadata = None
        if metadata_file.exists():
            import json
            with open(metadata_file, 'r') as f:
                dataset_metadata = json.load(f)
        
        return _load_from_alphagenome_api(
            dataset_dir, sample_id, window_center_size, gene_filter, config, dataset_metadata
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
    config: Dict,
    dataset_metadata: Optional[Dict] = None
) -> Tuple[np.ndarray, List[str]]:
    """Load data from AlphaGenome API."""
    console.print(f"[bold cyan]🌐 AlphaGenome API Mode Enabled[/bold cyan]")
    
    # Carregar dataset_metadata se não fornecido
    if dataset_metadata is None:
        metadata_file = dataset_dir / 'dataset_metadata.json'
        if metadata_file.exists():
            import json
            with open(metadata_file, 'r') as f:
                dataset_metadata = json.load(f)
        else:
            dataset_metadata = {}
    
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
    
    # Obter ontology terms dos metadados do dataset
    ontology_terms = get_ontology_terms_from_metadata(dataset_metadata)
    
    # Determine window size from window_center_size
    # Map common sizes to AlphaGenome constants (usando constantes globais)
    window_size_map = SIZE_TO_SEQUENCE_NAME
    
    # Use window_center_size_key from config if available, otherwise infer from window_center_size
    if 'window_center_size_key' in config:
        window_size_key = config['window_center_size_key']
    else:
        window_size_key = window_size_map.get(window_center_size, "SEQUENCE_LENGTH_16KB")
    
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
    config: Dict,
    dataset_metadata: Optional[Dict] = None
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
    # Carregar dataset_metadata se não fornecido
    if dataset_metadata is None:
        metadata_file = dataset_dir / 'dataset_metadata.json'
        if metadata_file.exists():
            import json
            with open(metadata_file, 'r') as f:
                dataset_metadata = json.load(f)
        else:
            dataset_metadata = {}
    
    if source == "api":
        # Call API with reference genome (como no Colab)
        ontology_terms = get_ontology_terms_from_metadata(dataset_metadata)
        
        # Determinar window size key usando mapeamento global
        # Criar dicionário completo com tamanhos adicionais que não estão no SIZE_TO_SEQUENCE_NAME
        extended_size_map = SIZE_TO_SEQUENCE_NAME
        
        window_size_key = extended_size_map.get(center_bp, "SEQUENCE_LENGTH_16KB")
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
    Suporta normalização global e por-track.
    
    Args:
        data: Array numpy com dados brutos
        method: Método de normalização ('log', 'minmax_keep_zero', 'zscore')
        norm_params: Parâmetros de normalização (do metadata.json)
        
    Returns:
        Array numpy normalizado
    """
    # Verificar se é normalização per-track
    per_track = norm_params.get('per_track', False)
    
    if per_track:
        # Normalização por track (cada linha normalizada independentemente)
        track_params = norm_params['track_params']
        normalized_rows = []
        
        for track_idx in range(data.shape[0]):
            track_data = data[track_idx:track_idx+1, :]  # Manter 2D [1, width]
            params = track_params[track_idx]
            
            if method == 'log':
                log_max = params['log_max']
                track_tensor = torch.from_numpy(track_data).float()
                track_norm = log_normalize(track_tensor, log_max)
                normalized_rows.append(track_norm.numpy())
                
            elif method == 'minmax_keep_zero':
                xmax = params['max']
                track_tensor = torch.from_numpy(track_data).float()
                track_norm = minmax_keep_zero(track_tensor, xmax)
                normalized_rows.append(track_norm.numpy())
                
            elif method == 'zscore':
                mean = params['mean']
                std = params['std']
                track_tensor = torch.from_numpy(track_data).float()
                track_norm = (track_tensor - mean) / std if std > 1e-8 else track_tensor
                normalized_rows.append(track_norm.numpy())
                
            else:
                raise ValueError(f"Método de normalização desconhecido: {method}")
        
        return np.vstack(normalized_rows)
        
    else:
        # Normalização global (backwards compatibility)
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
    label2: str = "AlphaGenome",
    dataset_metadata: Optional[Dict] = None
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
        
        # Criar label apenas com ontologia (gene já está no título)
        # Gerar labels dinamicamente dos metadados
        if dataset_metadata:
            ontology_labels = generate_ontology_labels(dataset_metadata)
        else:
            ontology_labels = generate_ontology_labels({})  # Usa fallback
        
        gene_idx = track_idx // len(ontology_labels)
        ont_idx = track_idx % len(ontology_labels)
        if ont_idx < len(ontology_labels):
            ylabel = ontology_labels[ont_idx]
        else:
            ylabel = f"Track {track_idx}"
        
        # Configurar subplot
        ax.set_ylabel(ylabel, fontsize=5, fontweight='bold')
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


def plot_individual_comparison(
    features_1: np.ndarray,
    features_2: np.ndarray,
    sample_id_1: str,
    pop_1: str,
    superpop_1: str,
    sample_id_2: str,
    pop_2: str,
    superpop_2: str,
    gene_name: str,
    config: Dict,
    viewer: Optional['InteractiveComparisonViewer'] = None,
    dataset_metadata: Optional[Dict] = None
) -> plt.Figure:
    """
    Plota comparação entre dois indivíduos.
    
    Args:
        features_1: Array [6, window_size] do primeiro indivíduo
        features_2: Array [6, window_size] do segundo indivíduo
        sample_id_1, pop_1, superpop_1: Informações do primeiro indivíduo
        sample_id_2, pop_2, superpop_2: Informações do segundo indivíduo
        gene_name: Nome do gene sendo visualizado
        config: Configuração
        viewer: Viewer interativo
    
    Returns:
        Figure do matplotlib
    """
    # Usar labels globais das ontologias
    
    num_tracks = features_1.shape[0]  # Deve ser 6
    
    # Criar subplots (6 linhas x 1 coluna para 6 tracks empilhadas verticalmente)
    fig, axes = plt.subplots(6, 1, figsize=(14, 12))
    
    # Título principal com gene
    label_1 = f"{sample_id_1} ({pop_1}/{superpop_1})"
    label_2 = f"{sample_id_2} ({pop_2}/{superpop_2})"
    
    fig.suptitle(
        f"Gene: {gene_name} - Comparação: {label_1} vs {label_2}",
        fontsize=16, fontweight='bold'
    )
    
    # Plotar cada track
    for i in range(num_tracks):
        ax = axes[i]
        
        # Plotar primeiro indivíduo (azul sólido)
        ax.plot(features_1[i], color='blue', linewidth=1.5, label=label_1, alpha=0.8)
        
        # Plotar segundo indivíduo (vermelho tracejado)
        ax.plot(features_2[i], color='red', linewidth=1.5, linestyle='--', 
                label=label_2, alpha=0.8)
        
        # Usar nome da ontologia como ylabel (gerar dinamicamente dos metadados)
        if dataset_metadata:
            ontology_labels = generate_ontology_labels(dataset_metadata)
        else:
            ontology_labels = generate_ontology_labels({})  # Usa fallback
        
        ylabel = ontology_labels[i] if i < len(ontology_labels) else f"Track {i}"
        ax.set_ylabel(ylabel, fontsize=5, fontweight='bold')
        ax.set_xlabel('Position', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='best')
    
    # Ajustar layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    # Adicionar instruções de navegação
    if config.get('show_navigation_help', True) and viewer is not None:
        nav_text = (
            "Navegação: ← → (ambos indivíduos) | "
            "A D (segundo indivíduo) | "
            "W Z (genes) | Q (sair)"
        )
        fig.text(0.5, 0.01, nav_text, ha='center', fontsize=9, 
                style='italic', color='gray')
    
    return fig


def plot_raw_data_colab_style(
    output,
    gene_name: str,
    sample_id: str,
    gtf,  # pandas DataFrame
    config: Dict,
    viewer: Optional[InteractiveViewer] = None,
    dataset_metadata: Optional[Dict] = None
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
        # dataset_metadata já está disponível como parâmetro da função
        return plot_raw_data(raw_data, sample_id, gene_name, config, viewer, dataset_metadata)
    
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
    viewer: Optional[InteractiveViewer] = None,
    dataset_metadata: Optional[Dict] = None
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
        
        # Labels - gerar dinamicamente dos metadados do dataset
        if dataset_metadata:
            ontology_labels = generate_ontology_labels(dataset_metadata)
        else:
            ontology_labels = generate_ontology_labels({})  # Usa fallback
        
        ylabel = ontology_labels[idx] if idx < len(ontology_labels) else f'Track {idx}'
        ax.set_ylabel(ylabel, fontsize=5, fontweight='bold')
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
    
    # Convert window_size_key to actual bp (usando constantes globais)
    window_size_map = SEQUENCE_LENGTH_MAP
    center_bp = window_size_map.get(window_size_key, 16384)
    
    for gene_name in genes_to_plot:
        console.print(f"\n[cyan]Loading raw data for {gene_name}...[/cyan]")
        
        raw_data, output, gtf = load_raw_alphagenome_data(
            dataset_dir, sample_id, gene_name, center_bp, source, config, dataset_metadata
        )
        
        # Plot with Colab-style if we have output object, otherwise simple plot
        if output is not None and gtf is not None:
            fig = plot_raw_data_colab_style(output, gene_name, sample_id, gtf, config, viewer, dataset_metadata)
        else:
            # dataset_metadata já está disponível como parâmetro da função
            fig = plot_raw_data(raw_data, sample_id, gene_name, config, viewer, dataset_metadata)
        
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
    # If viewer is navigating, use index to get sample_id (ignore fixed sample_id in config)
    if viewer is not None and splits_data:
        split = config.get('split', 'train')
        format_version = splits_data.get('format_version', 1)
        sample_meta = get_sample_info(splits_data, split, index, format_version, dataset_metadata)
        sample_id = sample_meta.get('sample_id', f"#{index + 1}")
    elif config.get('sample_id'):
        sample_id = config['sample_id']
    else:
        # If comparison_mode is dataset_dir_x_cache_dir, use index from splits
        if comparison_mode == 'dataset_dir_x_cache_dir' and splits_data:
            split = config.get('split', 'train')
            format_version = splits_data.get('format_version', 1)
            split_size = get_split_sample_count(splits_data, split, format_version)
            
            if index >= split_size:
                console.print(f"[red]Index {index} out of range for split {split} (size: {split_size})[/red]")
                return False
            
            # Get sample_id using the helper function
            sample_meta = get_sample_info(splits_data, split, index, format_version, dataset_metadata)
            sample_id = sample_meta.get('sample_id', f"#{index + 1}")
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
                ontology_terms = get_ontology_terms_from_metadata(dataset_metadata)
                
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
                # Usar constante global
                window_size_map = SEQUENCE_LENGTH_MAP
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
                ontology_terms = get_ontology_terms_from_metadata(dataset_metadata)
                
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
                # Usar constante global
                window_size_map = SEQUENCE_LENGTH_MAP
                viz_length = window_size_map.get(window_size_key, 16384)
                
                console.print(f"[dim]  Extracting center window for visualization: {viz_length} bp[/dim]")
                data1 = extract_center_window(alphagenome_data_T, viz_length, axis=1)
                data2 = extract_center_window(data_dataset_full_T, viz_length, axis=1)
                
                console.print(f"[dim]  AlphaGenome shape: {data1.shape}, Dataset shape: {data2.shape}[/dim]")
                
                label1 = "AlphaGenome (Ref)"
                label2 = "Dataset Dir"
                
            elif comparison_mode == "alphagenome_ind_x_dataset_dir":
                # Mode 2: AlphaGenome (individual) vs dataset_dir
                # USA build_window_and_predict.py como biblioteca
                
                console.print(f"[dim]  Detecting dataset_dir window size...[/dim]")
                ontology_terms = get_ontology_terms_from_metadata(dataset_metadata)
                
                # Carregar dataset_dir SEM extrair janela
                data_dataset_full, full_length = _load_dataset_dir_predictions(
                    dataset_dir, sample_id, gene_name, window_size_key=None
                )
                
                # Detectar o window_size_key usado no dataset_dir
                window_size_key_dataset = _detect_window_size_key_from_length(full_length)
                console.print(f"[cyan]  Dataset was generated with: {window_size_key_dataset} ({full_length} bp)[/cyan]")
                console.print(f"[cyan]  Will visualize center: {window_size_key}[/cyan]")
                
                # AlphaGenome (individual) - usando build_window_and_predict.py
                console.print(f"[dim]  Generating AlphaGenome predictions for {sample_id}...[/dim]")
                alphagenome_values = _get_alphagenome_predictions_for_individual(
                    sample_id=sample_id,
                    gene_name=gene_name,
                    window_size_key=window_size_key_dataset,
                    config=config,
                    ontology_terms=ontology_terms
                )
                
                if alphagenome_values is None:
                    console.print(f"[red]  Failed to generate AlphaGenome data[/red]")
                    continue
                
                # alphagenome_values: [full_length, 6]
                # data_dataset_full: [full_length, 6]
                
                # Transpor ambos para [6, full_length]
                alphagenome_data_T = alphagenome_values.T
                data_dataset_full_T = data_dataset_full.T
                
                # EXTRAIR JANELA CENTRAL APENAS PARA VISUALIZAÇÃO
                # Usar constante global
                window_size_map = SEQUENCE_LENGTH_MAP
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
                
                # Load dataset_dir data - returns [window_size, 6]
                dataset_data, _ = _load_dataset_dir_predictions(
                    dataset_dir, sample_id, gene_name, window_size_key
                )
                
                # Transpose to [6, window_size] for normalization
                dataset_data_T = dataset_data.T  # [6, window_size]
                
                # Get gene index to select correct track_params
                gene_idx = genes_in_order.index(gene_name)
                track_offset = gene_idx * 6  # 6 tracks per gene
                
                # Create subset of norm_params for this gene's tracks
                gene_norm_params = {
                    'method': norm_params['method'],
                    'per_track': norm_params.get('per_track', False)
                }
                if gene_norm_params['per_track'] and 'track_params' in norm_params:
                    # Extract only the track_params for this gene
                    all_track_params = norm_params['track_params']
                    gene_norm_params['track_params'] = all_track_params[track_offset:track_offset + 6]
                    console.print(f"[dim]  Using track_params[{track_offset}:{track_offset + 6}] for gene {gene_name}[/dim]")
                else:
                    # Copy global params
                    for key in ['log_max', 'max', 'mean', 'std']:
                        if key in norm_params:
                            gene_norm_params[key] = norm_params[key]
                
                # Apply normalization
                console.print(f"[dim]  Applying normalization...[/dim]")
                norm_result = apply_normalization(
                    dataset_data_T,  # [6, window_size]
                    gene_norm_params['method'],
                    gene_norm_params
                )
                # Convert to numpy if it's a tensor
                if torch.is_tensor(norm_result):
                    dataset_data_normalized = norm_result.numpy()  # [6, window_size]
                else:
                    dataset_data_normalized = norm_result  # [6, window_size]
                
                # Load cache data
                console.print(f"[dim]  Loading cache data...[/dim]")
                cache_dir = Path(config['cache_dir'])
                split = config['split']
                
                # Find index in split
                if splits_data is None:
                    splits_file = cache_dir / 'splits_metadata.json'
                    with open(splits_file) as f:
                        splits_data = json.load(f)
                
                format_version = splits_data.get('format_version', 1)
                split_found = None
                cache_index = None
                
                if format_version >= 2:
                    # New format: search for sample_id in split metadata
                    for split_name in ['train', 'val', 'test']:
                        split_data = splits_data.get(split_name, [])
                        for idx, sample_meta in enumerate(split_data):
                            if sample_meta.get('sample_id') == sample_id:
                                cache_index = idx
                                split_found = split_name
                                break
                        if split_found:
                            break
                else:
                    # Legacy format: search by dataset index
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
                    console.print(f"[red]  Sample {sample_id} not found in any split[/red]")
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
                
                # dataset_data_normalized is already [6, window_size] after normalization fix
                console.print(f"[dim]  dataset_data_normalized type: {type(dataset_data_normalized)}, shape: {dataset_data_normalized.shape}[/dim]")
                console.print(f"[dim]  cache_features type: {type(cache_features)}[/dim]")
                
                # Make sure cache_features is a numpy array
                if torch.is_tensor(cache_features):
                    cache_features_np = cache_features.numpy()
                else:
                    cache_features_np = cache_features
                
                # Match window sizes by extracting central region from the larger one
                dataset_window_size = dataset_data_normalized.shape[1]  # e.g., 16384
                cache_window_size = cache_features_np.shape[1]  # e.g., 32768
                
                console.print(f"[dim]  Dataset window: {dataset_window_size}, Cache window: {cache_window_size}[/dim]")
                
                if dataset_window_size > cache_window_size:
                    # Extract center from dataset_data to match cache
                    dataset_data_trimmed = extract_center_window(dataset_data_normalized, cache_window_size, axis=1)
                    cache_data_trimmed = cache_features_np
                    console.print(f"[dim]  Trimmed dataset_data to cache window size: {dataset_data_trimmed.shape}[/dim]")
                elif cache_window_size > dataset_window_size:
                    # Extract center from cache to match dataset_data
                    cache_data_trimmed = extract_center_window(cache_features_np, dataset_window_size, axis=1)
                    dataset_data_trimmed = dataset_data_normalized
                    console.print(f"[dim]  Trimmed cache to dataset window size: {cache_data_trimmed.shape}[/dim]")
                else:
                    dataset_data_trimmed = dataset_data_normalized
                    cache_data_trimmed = cache_features_np
                
                data1 = dataset_data_trimmed  # [6, window_size]
                data2 = cache_data_trimmed  # [6, window_size]
                
                label1 = "Dataset Dir (normalized)"
                label2 = "Cache Dir"
                
            else:
                console.print(f"[red]  Invalid comparison_mode: {comparison_mode}[/red]")
                continue
            
            # Compute metrics first
            console.print(f"[dim]  data1 shape: {data1.shape}, data2 shape: {data2.shape}[/dim]")
            metrics = compute_metrics(data1, data2, config.get('verbose_metrics', True))
            
            # Plot comparison
            # dataset_metadata já está disponível como parâmetro da função
            fig = plot_comparison(
                data1, data2, sample_id, metrics, [gene_name], config, viewer,
                label1=label1, label2=label2, dataset_metadata=dataset_metadata
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
            
            # Interactive mode: store figure in viewer, don't show yet
            if viewer is not None:
                viewer.fig = fig
                # Return after first gene in interactive mode
                return True
            else:
            plt.show()
            
        except Exception as e:
            console.print(f"[red]  Error processing {gene_name}: {e}[/red]")
            import traceback
            traceback.print_exc()
            continue
    
    # Cleanup: remover arquivos temporários
    tmp_dir = Path("/tmp/GENOMICS_DATA/top3")
    if tmp_dir.exists():
        try:
            console.print(f"[dim]  Cleaning up temporary files: {tmp_dir}[/dim]")
            shutil.rmtree(tmp_dir)
            console.print(f"[green]  ✓ Cleanup completed[/green]")
        except Exception as e:
            console.print(f"[yellow]  Warning: Could not clean up {tmp_dir}: {e}[/yellow]")
    
    console.print(f"\n[green]✓ Comparison mode processing completed![/green]")
    return True


def process_comparison_sample(
    config: Dict,
    index_1: int,
    index_2: int,
    gene_name: str,
    splits_data: Dict,
    dataset_metadata: Dict,
    viewer: 'InteractiveComparisonViewer'
) -> Optional[plt.Figure]:
    """
    Processa e visualiza comparação entre dois indivíduos.
    
    Args:
        config: Configuração
        index_1: Índice do primeiro indivíduo no split
        index_2: Índice do segundo indivíduo no split
        gene_name: Nome do gene a visualizar
        splits_data: Dados dos splits
        dataset_metadata: Metadata do dataset
        viewer: Viewer de comparação
    
    Returns:
        Figure do matplotlib ou None se erro
    """
    cache_dir = Path(config['cache_dir'])
    split = config.get('split', 'test')
    
    try:
        # Carregar dados do primeiro indivíduo
        features_1_raw, metadata_1, global_index_1 = load_cache_data(
            cache_dir, split, index_1, gene_filter=None
        )
        sample_id_1, pop_1, superpop_1 = viewer.get_sample_info(index_1)
        
        # Carregar dados do segundo indivíduo
        features_2_raw, metadata_2, global_index_2 = load_cache_data(
            cache_dir, split, index_2, gene_filter=None
        )
        sample_id_2, pop_2, superpop_2 = viewer.get_sample_info(index_2)
        
        # Obter lista de genes ordenada (usar a mesma da classe viewer)
        genes_in_order = viewer.available_genes
        
        # Filtrar por gene específico
        features_1, _ = filter_genes_from_features(
            features_1_raw, gene_name, genes_in_order
        )
        features_2, _ = filter_genes_from_features(
            features_2_raw, gene_name, genes_in_order
        )
        
        # Plotar comparação
        # dataset_metadata já está disponível como parâmetro da função
        fig = plot_individual_comparison(
            features_1, features_2,
            sample_id_1, pop_1, superpop_1,
            sample_id_2, pop_2, superpop_2,
            gene_name, config, viewer, dataset_metadata
        )
        
        return fig
        
    except Exception as e:
        console.print(f"[red]Erro ao processar comparação: {e}[/red]")
        import traceback
        traceback.print_exc()
        return None


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
        # dataset_metadata já está disponível como parâmetro da função
        fig = plot_comparison(
            cache_features, alphagenome_normalized, sample_id, metrics,
            genes_loaded, config, viewer, dataset_metadata=dataset_metadata
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
        
        # Validate window_size_key AlphaGenome
        valid_window_sizes = ['SEQUENCE_LENGTH_16KB', 'SEQUENCE_LENGTH_100KB', 
                             'SEQUENCE_LENGTH_500KB', 'SEQUENCE_LENGTH_1MB']
        window_size_key = raw_mode.get('window_size_key', 'SEQUENCE_LENGTH_16KB')
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
            splits_file = cache_dir / 'splits_metadata.json'
            with open(splits_file) as f:
                splits_data = json.load(f)
            
            norm_params_file = cache_dir / 'normalization_params.json'
            with open(norm_params_file) as f:
                norm_params = json.load(f)
        
        # Process sample in comparison mode
        if config.get('interactive_mode', False):
            # Create splits_data if not available (for modes without cache)
            if splits_data is None:
                # Create fake splits_data based on dataset_metadata individuals
                individuals = dataset_metadata.get('individuals', [])
                individuals_pedigree = dataset_metadata.get('individuals_pedigree', {})
                
                # Build metadata for all individuals
                all_samples = []
                for ind in individuals:
                    pedigree = individuals_pedigree.get(ind, {})
                    all_samples.append({
                        'sample_id': ind,
                        'superpopulation': pedigree.get('superpopulation', 'UNK'),
                        'population': pedigree.get('population', 'UNK'),
                        'sex': pedigree.get('sex', 0)
                    })
                
                splits_data = {
                    'format_version': 2,
                    'train': all_samples,
                    'val': [],
                    'test': []
                }
                # Force split to 'train' for navigation
                config['split'] = 'train'
            
            # Check interactive_comparison_mode
            interactive_mode_type = config.get('interactive_comparison_mode', 'single')
            
            if interactive_mode_type == 'comparison':
                # Modo de comparação entre dois indivíduos
                metadata_csv = config.get('data_sources', {}).get('metadata_csv')
                if metadata_csv is None:
                    metadata_csv = "../docs/1000_genomes_metadata.csv"
                
                viewer = InteractiveComparisonViewer(
                    config, splits_data, dataset_metadata, metadata_csv
                )
                
                console.print(f"[cyan]Modo de comparação interativa:[/cyan]")
                console.print(f"  • Split: {viewer.current_split}")
                console.print(f"  • Total de amostras: {viewer.split_size}")
                console.print(f"  • Genes: {viewer.available_genes}")
                console.print(f"[yellow]  • ← → (ambos), A D (ind2), W Z (genes), Q (sair)[/yellow]\n")
                
                # Interactive loop for comparison mode
                while not viewer.should_exit:
                    current_gene = viewer.get_current_gene()
                    
                    # Use process_sample_comparison_mode for each individual separately
                    # and combine in a comparison view
                    fig = process_comparison_sample(
                        config, viewer.current_index_1, viewer.current_index_2,
                        current_gene, splits_data, dataset_metadata, viewer
                    )
                    
                    if fig is not None:
                        viewer.fig = fig
                        fig.canvas.mpl_connect('key_press_event', viewer.on_key_press)
                        plt.show()
                    else:
                        break
                
                console.print("\n[bold green]✓ Interactive comparison session ended[/bold green]\n")
            else:
                # Modo single - um indivíduo por vez
                viewer = InteractiveViewer(config, splits_data, dataset_metadata)
                
                console.print(f"[yellow]Modo interativo: Use ← → para navegar, 'q' para sair[/yellow]\n")
                
                # Interactive loop
                while not viewer.should_exit:
                    success = process_sample_comparison_mode(
                        config, viewer.current_index, dataset_metadata,
                        splits_data=splits_data,
                        norm_params=norm_params,
                        viewer=viewer
                    )
                    
                    if not success:
                        break
                    
                    # Connect keyboard handler if figure exists
                    if viewer.fig is not None:
                        viewer.fig.canvas.mpl_connect('key_press_event', viewer.on_key_press)
                        plt.show()
                
                console.print("\n[bold green]✓ Interactive session ended[/bold green]\n")
            
            return
        else:
            # Non-interactive: process single sample
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
            with open(cache_dir / 'splits_metadata.json', 'r') as f:
                splits_data = json.load(f)
            
            # Carregar normalization params
            with open(cache_dir / 'normalization_params.json', 'r') as f:
                norm_params = json.load(f)
            
            # Modo interativo ou não
            if config.get('interactive_mode', False):
                # Verificar qual tipo de modo interativo
                comparison_mode = config.get('interactive_comparison_mode', 'single')
                
                if comparison_mode == 'comparison':
                    # NOVO: Modo de comparação entre indivíduos
                    metadata_csv = config.get('data_sources', {}).get('metadata_csv')
                    if metadata_csv is None:
                        metadata_csv = "../docs/1000_genomes_metadata.csv"
                    
                    viewer = InteractiveComparisonViewer(
                        config, splits_data, dataset_metadata, metadata_csv
                    )
                    
                    console.print(f"[cyan]Modo de comparação interativa ativado:[/cyan]")
                    console.print(f"  • Split: {viewer.current_split}")
                    console.print(f"  • Total de amostras: {len(viewer.split_indices)}")
                    console.print(f"  • Total de genes: {len(viewer.available_genes)}")
                    console.print(f"[yellow]  • Use ← → (ambos), A D (ind2), W Z (genes), Q (sair)[/yellow]\n")
                    
                    # Loop interativo
                    while not viewer.should_exit:
                        current_gene = viewer.get_current_gene()
                        fig = process_comparison_sample(
                            config, viewer.current_index_1, viewer.current_index_2,
                            current_gene, splits_data, dataset_metadata, viewer
                        )
                        
                        if fig is not None:
                            viewer.fig = fig
                            fig.canvas.mpl_connect('key_press_event', viewer.on_key_press)
                            plt.show()
                        else:
                            break
                    
                    console.print("\n[bold green]✓ Sessão de comparação interativa encerrada[/bold green]\n")
                
                else:
                    # Modo interativo existente (single)
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
