#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Neural Longevity Dataset Builder

Módulo para construção de dataset de marcadores genéticos e epigenéticos
de longevidade usando AlphaGenome e PyTorch.

Pipeline:
1. Download de genomas (1000 Genomes ou outras fontes)
2. Identificação de variantes vs referência
3. Seleção de pontos centrais (variantes)
4. Extração de sequências FASTA
5. Processamento com AlphaGenome
6. Construção de dataset PyTorch

Autor: IA Neuro-Simbólica para Longevidade
Data: Outubro 2025
"""

import argparse
import sys
import yaml
import json
import pickle
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
import subprocess as sp
import hashlib
from collections import defaultdict
import gzip

import numpy as np
import pandas as pd
from tqdm import tqdm

# PyTorch
import torch
from torch.utils.data import Dataset, DataLoader, random_split

# Rich para terminal
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint

console = Console()


# ═══════════════════════════════════════════════════════════════════
# Data Classes
# ═══════════════════════════════════════════════════════════════════

@dataclass
class GenomicVariant:
    """Representa uma variante genômica."""
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    quality: float
    depth: int
    allele_frequency: float = 0.0
    filter_status: str = "PASS"
    variant_type: str = "SNV"  # SNV, INSERTION, DELETION
    
    def to_dict(self) -> Dict:
        return {
            'chromosome': self.chromosome,
            'position': self.position,
            'ref_allele': self.ref_allele,
            'alt_allele': self.alt_allele,
            'quality': self.quality,
            'depth': self.depth,
            'allele_frequency': self.allele_frequency,
            'filter_status': self.filter_status,
            'variant_type': self.variant_type
        }


@dataclass
class CentralPoint:
    """Ponto central para extração de sequência."""
    variant: GenomicVariant
    importance_score: float = 0.0
    selected_for_dataset: bool = False


@dataclass
class SequenceRecord:
    """Registro de sequência extraída."""
    sample_id: str
    central_point: CentralPoint
    sequence: str
    fasta_file: Path
    label: int  # 1=longevo, 0=não-longevo
    
    def to_dict(self) -> Dict:
        return {
            'sample_id': self.sample_id,
            'variant': self.central_point.variant.to_dict(),
            'sequence': self.sequence,
            'sequence_length': len(self.sequence),
            'label': self.label
        }


@dataclass
class AlphaGenomeResult:
    """Resultado de predição do AlphaGenome."""
    sequence_record: SequenceRecord
    predictions: Dict[str, np.ndarray]  # output_type -> array de predições
    metadata: Dict[str, pd.DataFrame]  # output_type -> DataFrame de metadados
    
    def aggregate_statistics(self) -> Dict[str, Dict[str, float]]:
        """Agrega estatísticas das predições."""
        stats = {}
        for output_type, pred_array in self.predictions.items():
            stats[output_type] = {
                'mean': float(np.mean(pred_array)),
                'std': float(np.std(pred_array)),
                'min': float(np.min(pred_array)),
                'max': float(np.max(pred_array)),
                'median': float(np.median(pred_array))
            }
        return stats


# ═══════════════════════════════════════════════════════════════════
# PyTorch Dataset
# ═══════════════════════════════════════════════════════════════════

class LongevityDataset(Dataset):
    """
    Dataset PyTorch para marcadores de longevidade.
    
    Cada amostra contém:
    - Sequência DNA (one-hot encoded)
    - Posição genômica
    - Predições AlphaGenome
    - Label (longevo=1, não-longevo=0)
    - Metadados
    """
    
    def __init__(self, data_file: Path, transform=None):
        """
        Args:
            data_file: Arquivo pickle com dados processados
            transform: Transformações opcionais
        """
        self.data_file = data_file
        self.transform = transform
        
        # Carregar dados
        with open(data_file, 'rb') as f:
            data = pickle.load(f)
        
        self.sequences = data['sequences']
        self.positions = data['positions']
        self.alphagenome_features = data['alphagenome_features']
        self.labels = data['labels']
        self.metadata = data['metadata']
        
        console.print(f"[green]✓ Dataset carregado: {len(self)} amostras[/green]")
    
    def __len__(self) -> int:
        return len(self.labels)
    
    def __getitem__(self, idx: int) -> Tuple[Dict[str, torch.Tensor], int]:
        """
        Retorna uma amostra do dataset.
        
        Returns:
            features: Dict com 'sequence', 'position', 'alphagenome'
            label: 0 ou 1
        """
        # Sequência (one-hot encoded)
        sequence = torch.FloatTensor(self.sequences[idx])
        
        # Posição normalizada
        position = torch.FloatTensor([self.positions[idx]])
        
        # Features AlphaGenome
        alphagenome = torch.FloatTensor(self.alphagenome_features[idx])
        
        # Label
        label = int(self.labels[idx])
        
        features = {
            'sequence': sequence,
            'position': position,
            'alphagenome': alphagenome,
            'metadata': self.metadata[idx]
        }
        
        if self.transform:
            features = self.transform(features)
        
        return features, label
    
    def get_metadata(self, idx: int) -> Dict:
        """Retorna metadados de uma amostra."""
        return self.metadata[idx]
    
    def get_class_distribution(self) -> Dict[int, int]:
        """Retorna distribuição de classes."""
        unique, counts = np.unique(self.labels, return_counts=True)
        return dict(zip(unique, counts))


# ═══════════════════════════════════════════════════════════════════
# Processador Principal
# ═══════════════════════════════════════════════════════════════════

class LongevityDatasetBuilder:
    """
    Construtor do dataset de longevidade.
    
    Orquestra todo o pipeline de construção do dataset.
    """
    
    def __init__(self, config_path: Path):
        """
        Args:
            config_path: Caminho para arquivo YAML de configuração
        """
        self.config_path = config_path
        
        # Carregar configuração
        with open(config_path) as f:
            self.config = yaml.safe_load(f)
        
        # Diretórios (usar diretório de execução atual)
        # output_dir é relativo ao diretório onde o comando é executado
        self.work_dir = Path.cwd()  # Diretório de execução
        self.output_dir = self.work_dir / self.config['project']['output_dir']
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        console.print(f"[cyan]Diretório de trabalho: {self.work_dir}[/cyan]")
        console.print(f"[cyan]Saída em: {self.output_dir}[/cyan]")
        
        # Cache dir também relativo ao work_dir se for caminho relativo
        cache_path = Path(self.config['alphagenome']['cache_dir'])
        if cache_path.is_absolute():
            self.cache_dir = cache_path
        else:
            self.cache_dir = self.work_dir / cache_path
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Estado do processamento
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.state = self._load_checkpoint()
        
        console.print(Panel.fit(
            f"[bold cyan]{self.config['project']['name']}[/bold cyan]\n"
            f"[dim]{self.config['project']['description']}[/dim]",
            border_style="cyan"
        ))
    
    def _load_checkpoint(self) -> Dict:
        """Carrega checkpoint se existir."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file) as f:
                return json.load(f)
        return {
            'samples_downloaded': [],
            'variants_extracted': [],
            'central_points_selected': False,
            'sequences_extracted': [],
            'alphagenome_processed': [],
            'dataset_built': False
        }
    
    def _save_checkpoint(self):
        """Salva checkpoint."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.state, f, indent=2)
    
    # ───────────────────────────────────────────────────────────────
    # Passo 1: Download de Amostras
    # ───────────────────────────────────────────────────────────────
    
    def download_samples(self):
        """
        Baixa VCFs de amostras do 1000 Genomes.
        
        Por enquanto, simula obtendo lista de amostras disponíveis.
        """
        console.print("\n[bold cyan]Passo 1: Download de Amostras[/bold cyan]")
        
        is_dry_run = self.config['debug']['dry_run']
        
        if is_dry_run:
            console.print("[yellow]⚠ Modo dry-run: simulando downloads[/yellow]")
        
        # Obter lista de amostras
        samples_file = self.output_dir / "samples_list.txt"
        
        if not samples_file.exists():
            if not is_dry_run:
                console.print("[cyan]Baixando lista de amostras do 1000 Genomes...[/cyan]")
            self._download_samples_list(samples_file)
        
        # Ler lista (filtrar comentários e linhas vazias)
        with open(samples_file) as f:
            all_samples = [
                line.strip() for line in f 
                if line.strip() and not line.startswith('#')
            ]
        
        console.print(f"[green]✓ {len(all_samples)} amostras disponíveis[/green]")
        
        # Selecionar amostras longevas e não-longevas
        longevous_range = self.config['data_sources']['longevous']['sample_range']
        non_longevous_range = self.config['data_sources']['non_longevous']['sample_range']
        
        longevous_samples = all_samples[longevous_range[0]:longevous_range[1]]
        non_longevous_samples = all_samples[non_longevous_range[0]:non_longevous_range[1]]
        
        console.print(f"[cyan]Amostras longevas: {len(longevous_samples)}[/cyan]")
        console.print(f"[cyan]Amostras não-longevas: {len(non_longevous_samples)}[/cyan]")
        
        # Salvar listas (mesmo em dry-run, para próximos passos)
        (self.output_dir / "longevous_samples.txt").write_text('\n'.join(longevous_samples))
        (self.output_dir / "non_longevous_samples.txt").write_text('\n'.join(non_longevous_samples))
        
        # Baixar VCFs (placeholder - implementação depende de acesso real)
        if not is_dry_run:
            self._download_vcfs(longevous_samples, label="longevous")
            self._download_vcfs(non_longevous_samples, label="non_longevous")
        else:
            console.print("[yellow]⚠ Dry-run: pulando download de VCFs[/yellow]")
        
        console.print("[green]✓ Download de amostras concluído[/green]")
    
    def _download_samples_list(self, output_file: Path):
        """Baixa lista de amostras do 1000 Genomes."""
        # URL da lista de amostras
        samples_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
        
        try:
            # Baixar arquivo temporário
            temp_file = output_file.parent / "samples_raw.tsv"
            cmd = ['wget', '-O', str(temp_file), samples_url]
            sp.run(cmd, check=True, capture_output=True)
            
            # Parsear TSV e extrair IDs de amostra (coluna 9: SAMPLE_NAME)
            sample_ids = set()  # Usar set para evitar duplicatas
            with open(temp_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) > 9:  # Garantir que tem colunas suficientes
                        sample_id = fields[9].strip()  # Coluna SAMPLE_NAME (índice 9)
                        if sample_id and sample_id != 'SAMPLE_NAME':  # Pular header
                            sample_ids.add(sample_id)
            
            # Salvar lista limpa de IDs únicos
            sorted_ids = sorted(sample_ids)
            output_file.write_text('\n'.join(sorted_ids))
            
            # Remover arquivo temporário
            temp_file.unlink()
            
            console.print(f"[green]✓ Lista de amostras baixada e processada[/green]")
            console.print(f"[green]  {len(sorted_ids)} IDs únicos extraídos[/green]")
            
        except Exception as e:
            console.print(f"[yellow]⚠ Erro ao baixar lista: {e}[/yellow]")
            # Criar lista simulada para desenvolvimento
            self._create_simulated_samples_list(output_file)
    
    def _create_simulated_samples_list(self, output_file: Path):
        """Cria lista simulada de amostras para desenvolvimento."""
        console.print("[yellow]Criando lista simulada de amostras...[/yellow]")
        # Criar IDs de amostras simuladas (sem comentários)
        samples = [f"SAMPLE_{i:04d}" for i in range(1000)]
        output_file.write_text('\n'.join(samples))
        console.print(f"[green]✓ {len(samples)} amostras simuladas criadas[/green]")
    
    def _download_vcfs(self, samples: List[str], label: str):
        """
        Baixa VCFs das amostras do 1000 Genomes (idempotente).
        
        Se o VCF já existir, não re-baixa.
        """
        vcf_dir = self.output_dir / f"vcfs_{label}"
        vcf_dir.mkdir(exist_ok=True)
        
        console.print(f"[cyan]Preparando VCFs para {len(samples)} amostras ({label})...[/cyan]")
        
        # URL base para VCFs do 1000 Genomes (30x coverage)
        vcf_base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"
        
        downloaded = 0
        skipped = 0
        failed = 0
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task(f"Baixando VCFs ({label})...", total=len(samples))
            
            for sample_id in samples:
                vcf_file = vcf_dir / f"{sample_id}.vcf.gz"
                vcf_index = vcf_dir / f"{sample_id}.vcf.gz.tbi"
                
                # Idempotência: pular se já existe
                if vcf_file.exists() and vcf_index.exists():
                    skipped += 1
                    progress.advance(task)
                    continue
                
                # Tentar baixar VCF chr por chr e concatenar, ou pegar merged
                # Por simplicidade, vamos baixar o VCF merged se disponível
                # Formato: sample_id subdividido por cromossomo ou merged
                
                # URL exemplo: .../CHR{chr}/{sample_id}.{chr}.vcf.gz
                # Para simplificar, vamos usar a estratégia de pular o download
                # e criar um placeholder que indica que deveria ser baixado
                
                # Como os VCFs são muito grandes e específicos por sample,
                # vamos criar apenas a estrutura e avisar ao usuário
                console.print(f"[yellow]  ⚠ VCF para {sample_id} deve ser obtido manualmente[/yellow]")
                console.print(f"[yellow]    URL base: {vcf_base_url}[/yellow]")
                
                failed += 1
                progress.advance(task)
        
        console.print(f"[green]✓ Processamento de VCFs ({label}):[/green]")
        console.print(f"[green]  Baixados: {downloaded}, Já existiam: {skipped}, Pulados: {failed}[/green]")
        
        if failed > 0:
            console.print(f"[yellow]⚠ {failed} VCFs não foram baixados automaticamente[/yellow]")
            console.print(f"[yellow]  Para obter VCFs reais:[/yellow]")
            console.print(f"[yellow]  1. Use genomes_analyzer.py para processar FASTQs, OU[/yellow]")
            console.print(f"[yellow]  2. Use neural_integration.py para extrair de VCFs existentes[/yellow]")
            console.print(f"[yellow]  3. O pipeline continuará com pontos centrais simulados[/yellow]")
    
    # ───────────────────────────────────────────────────────────────
    # Passo 2: Extração de Variantes
    # ───────────────────────────────────────────────────────────────
    
    def extract_variants(self, sample_id: str, vcf_path: Path) -> List[GenomicVariant]:
        """
        Extrai variantes de um VCF.
        
        Args:
            sample_id: ID da amostra
            vcf_path: Caminho para VCF
            
        Returns:
            Lista de variantes filtradas
        """
        variants = []
        filters = self.config['variant_selection']['filters']
        
        try:
            # Usar bcftools para extrair variantes
            cmd = [
                'bcftools', 'view',
                '-f', 'PASS' if filters['filter_pass_only'] else '.',
                '-i', f'QUAL>={filters["min_quality"]} && FORMAT/DP>={filters["min_depth"]}',
                str(vcf_path)
            ]
            
            result = sp.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse VCF
            for line in result.stdout.split('\n'):
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.split('\t')
                if len(fields) < 10:
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4].split(',')[0]  # Primeira alternativa
                qual = float(fields[5]) if fields[5] != '.' else 0.0
                filter_status = fields[6]
                
                # Parse INFO para AF
                info = dict(item.split('=', 1) for item in fields[7].split(';') if '=' in item)
                af = float(info.get('AF', '0.0').split(',')[0])
                
                # Parse FORMAT para DP
                format_fields = fields[8].split(':')
                format_values = fields[9].split(':')
                format_dict = dict(zip(format_fields, format_values))
                dp = int(format_dict.get('DP', '0'))
                
                # Filtros
                if filters['exclude_common'] and af > filters['max_allele_frequency']:
                    continue
                
                # Determinar tipo de variante
                if len(ref) == 1 and len(alt) == 1:
                    var_type = "SNV"
                elif len(alt) > len(ref):
                    var_type = "INSERTION"
                else:
                    var_type = "DELETION"
                
                variant = GenomicVariant(
                    chromosome=chrom,
                    position=pos,
                    ref_allele=ref,
                    alt_allele=alt,
                    quality=qual,
                    depth=dp,
                    allele_frequency=af,
                    filter_status=filter_status,
                    variant_type=var_type
                )
                
                variants.append(variant)
        
        except Exception as e:
            console.print(f"[red]✗ Erro ao extrair variantes de {vcf_path}: {e}[/red]")
        
        return variants
    
    # ───────────────────────────────────────────────────────────────
    # Passo 3: Seleção de Pontos Centrais
    # ───────────────────────────────────────────────────────────────
    
    def select_central_points(self) -> List[CentralPoint]:
        """
        Seleciona pontos centrais para construção do dataset.
        
        Estratégia inicial: selecionar N variantes da primeira pessoa longeva.
        
        Returns:
            Lista de pontos centrais selecionados
        """
        console.print("\n[bold cyan]Passo 3: Seleção de Pontos Centrais[/bold cyan]")
        
        strategy = self.config['variant_selection']['initial_strategy']
        n_points = self.config['variant_selection']['n_central_points']
        
        if strategy == "first_longevous_sample":
            # Obter primeira amostra longeva
            longevous_samples_file = self.output_dir / "longevous_samples.txt"
            
            # Verificar se arquivo existe (pode não existir em dry-run)
            if not longevous_samples_file.exists():
                console.print(f"[yellow]⚠ Arquivo não encontrado: {longevous_samples_file}[/yellow]")
                console.print(f"[yellow]  Criando pontos centrais simulados...[/yellow]")
                return self._create_simulated_central_points(n_points)
            
            with open(longevous_samples_file) as f:
                first_sample = f.readline().strip()
            
            console.print(f"[cyan]Selecionando variantes de: {first_sample}[/cyan]")
            
            # Obter VCF (assumindo que existe ou foi gerado)
            vcf_path = self.output_dir / f"vcfs_longevous/{first_sample}.vcf.gz"
            
            if not vcf_path.exists():
                console.print(f"[yellow]⚠ VCF não encontrado: {vcf_path}[/yellow]")
                console.print(f"[yellow]  Criando pontos centrais simulados...[/yellow]")
                return self._create_simulated_central_points(n_points)
            
            # Extrair variantes
            variants = self.extract_variants(first_sample, vcf_path)
            
            console.print(f"[green]✓ {len(variants)} variantes encontradas[/green]")
            
            # Selecionar top N variantes (por qualidade)
            variants_sorted = sorted(variants, key=lambda v: v.quality, reverse=True)
            selected_variants = variants_sorted[:n_points]
            
            central_points = [
                CentralPoint(variant=v, importance_score=v.quality, selected_for_dataset=True)
                for v in selected_variants
            ]
            
            # Salvar pontos centrais
            self._save_central_points(central_points)
            
            console.print(f"[green]✓ {len(central_points)} pontos centrais selecionados[/green]")
            
            return central_points
        
        else:
            raise NotImplementedError(f"Estratégia '{strategy}' não implementada")
    
    def _create_simulated_central_points(self, n_points: int) -> List[CentralPoint]:
        """Cria pontos centrais simulados para desenvolvimento."""
        central_points = []
        for i in range(n_points):
            variant = GenomicVariant(
                chromosome=f"chr{(i % 22) + 1}",
                position=1000000 + i * 100000,
                ref_allele="A",
                alt_allele="G",
                quality=30.0,
                depth=20,
                allele_frequency=0.1
            )
            central_points.append(
                CentralPoint(variant=variant, importance_score=1.0, selected_for_dataset=True)
            )
        return central_points
    
    def _save_central_points(self, points: List[CentralPoint]):
        """Salva pontos centrais em arquivo."""
        output_file = self.output_dir / "central_points.json"
        data = [
            {
                'variant': p.variant.to_dict(),
                'importance_score': p.importance_score,
                'selected': p.selected_for_dataset
            }
            for p in points
        ]
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        console.print(f"[green]✓ Pontos centrais salvos em: {output_file}[/green]")
    
    def _load_central_points(self) -> List[CentralPoint]:
        """Carrega pontos centrais de arquivo."""
        input_file = self.output_dir / "central_points.json"
        with open(input_file) as f:
            data = json.load(f)
        
        points = []
        for item in data:
            var_dict = item['variant']
            variant = GenomicVariant(**var_dict)
            point = CentralPoint(
                variant=variant,
                importance_score=item['importance_score'],
                selected_for_dataset=item['selected']
            )
            points.append(point)
        
        return points
    
    # ───────────────────────────────────────────────────────────────
    # Passo 4: Extração de Sequências FASTA
    # ───────────────────────────────────────────────────────────────
    
    def extract_sequences(self, central_points: List[CentralPoint]) -> List[SequenceRecord]:
        """
        Extrai sequências FASTA centradas nos pontos selecionados.
        
        Args:
            central_points: Pontos centrais para extração
            
        Returns:
            Lista de registros de sequência
        """
        console.print("\n[bold cyan]Passo 4: Extração de Sequências FASTA[/bold cyan]")
        
        if self.config['debug']['dry_run']:
            console.print("[yellow]⚠ Modo dry-run: simulando extração de sequências[/yellow]")
            console.print(f"[cyan]Seria extraído: {len(central_points)} pontos centrais × N amostras[/cyan]")
            return []
        
        window_size = self.config['sequence_extraction']['window_size']
        
        # Resolver caminho da referência (relativo ao work_dir ou absoluto)
        ref_path = Path(self.config['data_sources']['reference']['fasta'])
        if ref_path.is_absolute():
            ref_fasta = ref_path
        else:
            ref_fasta = self.work_dir / ref_path
        
        if not ref_fasta.exists():
            console.print(f"[red]✗ Genoma de referência não encontrado: {ref_fasta}[/red]")
            return []
        
        # Obter todas as amostras
        longevous_samples = self._load_sample_list("longevous_samples.txt", label=1)
        non_longevous_samples = self._load_sample_list("non_longevous_samples.txt", label=0)
        all_samples = longevous_samples + non_longevous_samples
        
        console.print(f"[cyan]Extraindo sequências para {len(all_samples)} amostras[/cyan]")
        console.print(f"[cyan]  {len(central_points)} pontos centrais × {len(all_samples)} amostras = {len(central_points) * len(all_samples)} sequências[/cyan]")
        
        records = []
        sequences_dir = self.output_dir / "sequences"
        sequences_dir.mkdir(exist_ok=True)
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task("Extraindo sequências...", total=len(all_samples) * len(central_points))
            
            for sample_id, label in all_samples:
                for point in central_points:
                    # Calcular região
                    chrom = point.variant.chromosome
                    center_pos = point.variant.position
                    start = max(1, center_pos - window_size // 2)
                    end = start + window_size
                    
                    region = f"{chrom}:{start}-{end}"
                    
                    # Extrair com samtools
                    fasta_file = sequences_dir / f"{sample_id}_{chrom}_{center_pos}.fasta"
                    
                    try:
                        cmd = ['samtools', 'faidx', str(ref_fasta), region]
                        result = sp.run(cmd, capture_output=True, text=True, check=True)
                        
                        # Parse sequência
                        lines = result.stdout.split('\n')
                        header = f">{sample_id}_{chrom}:{start}-{end}"
                        sequence = ''.join(lines[1:]).replace('\n', '').upper()
                        
                        # Salvar FASTA
                        with open(fasta_file, 'w') as f:
                            f.write(f"{header}\n{sequence}\n")
                        
                        # Criar registro
                        record = SequenceRecord(
                            sample_id=sample_id,
                            central_point=point,
                            sequence=sequence,
                            fasta_file=fasta_file,
                            label=label
                        )
                        records.append(record)
                    
                    except Exception as e:
                        console.print(f"[yellow]⚠ Erro ao extrair {region}: {e}[/yellow]")
                    
                    progress.advance(task)
        
        console.print(f"[green]✓ {len(records)} sequências extraídas[/green]")
        
        # Salvar índice de sequências
        self._save_sequence_index(records)
        
        return records
    
    def _load_sample_list(self, filename: str, label: int) -> List[Tuple[str, int]]:
        """Carrega lista de amostras com labels."""
        samples_file = self.output_dir / filename
        if not samples_file.exists():
            return []
        
        with open(samples_file) as f:
            samples = [(line.strip(), label) for line in f if line.strip()]
        
        return samples
    
    def _save_sequence_index(self, records: List[SequenceRecord]):
        """Salva índice de sequências extraídas."""
        index_file = self.output_dir / "sequences_index.json"
        data = [rec.to_dict() | {'fasta_file': str(rec.fasta_file)} for rec in records]
        with open(index_file, 'w') as f:
            json.dump(data, f, indent=2)
        console.print(f"[green]✓ Índice de sequências salvo: {index_file}[/green]")
    
    # ───────────────────────────────────────────────────────────────
    # Passo 5: Processamento com AlphaGenome
    # ───────────────────────────────────────────────────────────────
    
    def run_alphagenome(self, records: List[SequenceRecord]) -> List[AlphaGenomeResult]:
        """
        Processa sequências com AlphaGenome.
        
        Args:
            records: Registros de sequências
            
        Returns:
            Lista de resultados AlphaGenome
        """
        console.print("\n[bold cyan]Passo 5: Processamento com AlphaGenome[/bold cyan]")
        
        if self.config['debug']['dry_run']:
            console.print("[yellow]⚠ Modo dry-run: simulando processamento com AlphaGenome[/yellow]")
            console.print(f"[cyan]Seria processado: {len(records)} sequências[/cyan]")
            return []
        
        api_key = self.config['alphagenome']['api_key']
        if not api_key:
            console.print("[red]✗ API key do AlphaGenome não configurada![/red]")
            console.print("[yellow]  Configure alphagenome.api_key em longevity_config.yaml[/yellow]")
            return []
        
        console.print(f"[cyan]Processando {len(records)} sequências com AlphaGenome...[/cyan]")
        
        # Usar neural_module.py para processar
        results = []
        
        for record in tqdm(records, desc="AlphaGenome"):
            result = self._process_with_alphagenome(record, api_key)
            if result:
                results.append(result)
        
        console.print(f"[green]✓ {len(results)} sequências processadas[/green]")
        
        return results
    
    def _process_with_alphagenome(self, record: SequenceRecord, api_key: str) -> Optional[AlphaGenomeResult]:
        """Processa uma sequência com AlphaGenome."""
        # Gerar hash para cache
        seq_hash = hashlib.md5(record.sequence.encode()).hexdigest()
        cache_file = self.cache_dir / f"{seq_hash}.pkl"
        
        # Verificar cache
        if cache_file.exists() and self.config['alphagenome']['cache_results']:
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
        
        # Processar com neural_module.py
        try:
            output_dir = self.cache_dir / f"temp_{seq_hash}"
            output_dir.mkdir(exist_ok=True)
            
            cmd = [
                'python', 'neural_module.py',
                '-i', str(record.fasta_file),
                '-k', api_key,
                '-o', str(output_dir),
                '--outputs'] + self.config['alphagenome']['outputs']
            
            sp.run(cmd, check=True, capture_output=True)
            
            # Carregar resultados
            report_file = output_dir / "analysis_report.json"
            with open(report_file) as f:
                report = json.load(f)
            
            # TODO: Parse predições e metadados dos arquivos gerados
            # Por enquanto, placeholder
            predictions = {}
            metadata = {}
            
            result = AlphaGenomeResult(
                sequence_record=record,
                predictions=predictions,
                metadata=metadata
            )
            
            # Salvar em cache
            with open(cache_file, 'wb') as f:
                pickle.dump(result, f)
            
            return result
        
        except Exception as e:
            console.print(f"[yellow]⚠ Erro ao processar com AlphaGenome: {e}[/yellow]")
            return None
    
    # ───────────────────────────────────────────────────────────────
    # Passo 6: Construção do Dataset PyTorch
    # ───────────────────────────────────────────────────────────────
    
    def build_dataset(self, results: List[AlphaGenomeResult]):
        """
        Constrói dataset PyTorch a partir dos resultados.
        
        Args:
            results: Resultados do AlphaGenome
        """
        console.print("\n[bold cyan]Passo 6: Construção do Dataset PyTorch[/bold cyan]")
        
        if self.config['debug']['dry_run']:
            console.print("[yellow]⚠ Modo dry-run: simulando construção do dataset[/yellow]")
            console.print(f"[cyan]Seria criado: train.pkl, val.pkl, test.pkl[/cyan]")
            return
        
        # TODO: Implementar construção completa
        console.print("[yellow]⚠ Construção de dataset em desenvolvimento...[/yellow]")
    
    # ───────────────────────────────────────────────────────────────
    # Pipeline Completo
    # ───────────────────────────────────────────────────────────────
    
    def run_pipeline(self):
        """Executa pipeline completo."""
        console.print("\n[bold green]Iniciando Pipeline de Construção do Dataset[/bold green]\n")
        
        steps = self.config['pipeline']['steps']
        
        try:
            # Passo 1: Download
            if steps['download_samples']:
                self.download_samples()
            
            # Passo 2: Seleção de pontos centrais
            if steps['select_central_points']:
                central_points = self.select_central_points()
            else:
                central_points = self._load_central_points()
            
            # Passo 3: Extração de sequências
            if steps['extract_sequences']:
                records = self.extract_sequences(central_points)
            else:
                # Carregar de índice
                records = []
            
            # Passo 4: AlphaGenome
            if steps['run_alphagenome']:
                results = self.run_alphagenome(records)
            else:
                results = []
            
            # Passo 5: Build dataset
            if steps['build_dataset']:
                self.build_dataset(results)
            
            console.print("\n[bold green]✓ Pipeline concluído com sucesso![/bold green]")
        
        except KeyboardInterrupt:
            console.print("\n[yellow]⚠ Pipeline interrompido pelo usuário[/yellow]")
            self._save_checkpoint()
        except Exception as e:
            console.print(f"\n[red]✗ Erro no pipeline: {e}[/red]")
            raise


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Construtor de Dataset para Marcadores de Longevidade',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:

1. Construir dataset completo:
   %(prog)s --config longevity_config.yaml

2. Apenas download de amostras:
   %(prog)s --config longevity_config.yaml --steps download_samples

3. Modo dry-run (simulação):
   %(prog)s --config longevity_config.yaml --dry-run
        """
    )
    
    parser.add_argument('--config', type=Path, required=True,
                       help='Arquivo de configuração YAML')
    parser.add_argument('--steps', nargs='+',
                       help='Executar apenas etapas específicas')
    parser.add_argument('--dry-run', action='store_true',
                       help='Modo simulação (não executa, apenas mostra o que faria)')
    
    args = parser.parse_args()
    
    # Criar builder
    builder = LongevityDatasetBuilder(args.config)
    
    # Dry-run
    if args.dry_run:
        builder.config['debug']['dry_run'] = True
    
    # Executar pipeline
    builder.run_pipeline()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n[yellow]⚠ Interrompido pelo usuário[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[red]✗ Erro fatal: {e}[/red]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")
        sys.exit(1)

