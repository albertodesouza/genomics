#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Neural Longevity Dataset Builder
===================================

Ferramenta de linha de comando para montar datasets de marcadores de
longevidade a partir de genomas individuais, executando o pipeline:

1. Download de genomas (ex.: 1000 Genomes High Coverage) a partir do ENA
2. Chamada de variantes contra uma referência (GRCh38) com *bcftools*
3. Seleção dos pontos centrais (variantes) que serão usados no dataset
4. Extração de janelas FASTA centradas nessas variantes, aplicando o alelo ALT
5. Processamento das sequências com o AlphaGenome para gerar *features*
6. Consolidação de um dataset PyTorch balanceado em *train/val/test*

Uso rápido
----------
1. Prepare um arquivo YAML de configuração (ex.: ``longevity_config.yaml``)
   contendo caminhos de saída, opções de AlphaGenome e parâmetros do dataset.
2. Garanta que as dependências de sistema estejam instaladas (``bcftools``,
   ``samtools``, AlphaGenome CLI, PyTorch, etc.).
3. Execute:

   ``python3 neural_longevity_dataset.py --config longevity_config.yaml``

   O programa respeita as etapas definidas em ``pipeline.steps`` dentro do YAML.

Selecionando etapas específicas
-------------------------------
Você pode rodar etapas isoladas informando ``--steps`` (ex.: apenas download e
extração de sequências):

``python3 neural_longevity_dataset.py --config longevity_config.yaml \
    --steps download_samples extract_sequences``

O argumento ``--dry-run`` mostra o que seria executado sem alterar arquivos.

Estrutura de saída
------------------
- ``cram/``: CRAM/CRAI baixados do ENA
- ``vcf/``: VCFs gerados via *bcftools*
- ``sequences/``: FASTAs centrados nas variantes com alelo aplicado por amostra
- ``alphagenome/``: caches das predições e estatísticas agregadas
- ``sequences_index.json``: índice consolidado das sequências extraídas e seus metadados
- ``central_points.json``: pontos centrais selecionados para o dataset
- ``torch_dataset/``: arquivos ``train.pkl``, ``validation.pkl`` e ``test.pkl``
  acompanhados de ``samples.csv`` com um resumo tabular das amostras

Metadados dos pontos centrais
-----------------------------
O arquivo ``central_points.json`` é gravado diretamente em ``project.output_dir``
(por exemplo, ``/dados/GENOMICS_DATA/top3/<nome_do_projeto>/central_points.json``)
logo após a etapa ``select_central_points``. Cada elemento da lista possui a
seguinte estrutura:

``variant``
    Objeto com as informações da variante genômica que origina o ponto central,
    contendo ``chromosome`` (cromossomo), ``position`` (posição 1-indexada),
    ``ref_allele`` e ``alt_allele`` (alelos de referência e alternativo), além de
    métricas opcionais como ``quality`` (QUAL do VCF), ``depth`` (profundidade de
    leitura), ``allele_frequency`` (frequência estimada), ``filter_status``
    (status da coluna FILTER) e ``variant_type`` (SNV, INSERTION ou DELETION).
``selected``
    Flag booleana indicando se o ponto central foi efetivamente escolhido para
    compor o dataset balanceado.
``source_sample_id``
    Identificador da amostra de onde a variante foi obtida. Para pontos
    simulados esse campo fica como ``null``.

**Importante**: Este script **deve** ser executado a partir do diretório
``/dados/GENOMICS_DATA/top3``. Todos os dados baixados ou gerados são
armazenados em subdiretórios dessa pasta. Use caminhos absolutos para
referenciar arquivos de configuração localizados no repositório em
``~/genomics``.

Autor: Alberto F. De Souza
Data: Outubro 2025
"""

import argparse
import sys
import yaml
import json
import pickle
import random
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any, Sequence, Set
from dataclasses import dataclass, field
import subprocess as sp
import hashlib
from collections import defaultdict
import gzip
import tempfile
import shutil
from urllib.parse import quote, unquote
import shlex

import numpy as np
import pandas as pd
from tqdm import tqdm

# PyTorch
import torch
from torch.utils.data import Dataset, DataLoader, random_split

# Rich para terminal
from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    DownloadColumn,
    TransferSpeedColumn,
)
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint

console = Console()


# Diretório raiz obrigatório para dados persistentes do projeto
DATA_ROOT = Path("/dados/GENOMICS_DATA/top3").resolve()


def _format_command(command: Sequence[Any]) -> str:
    """Formata uma lista de argumentos para exibição em shell."""

    return " ".join(shlex.quote(str(arg)) for arg in command)


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
    source_sample_id: Optional[str] = None

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
            'variant_type': self.variant_type,
            'source_sample_id': self.source_sample_id,
        }


@dataclass
class CentralPoint:
    """Ponto central para extração de sequência."""
    variant: GenomicVariant
    selected_for_dataset: bool = False
    source_sample_id: Optional[str] = None
    importance_score: Optional[float] = None

    def __post_init__(self):
        if self.source_sample_id and not self.variant.source_sample_id:
            self.variant.source_sample_id = self.source_sample_id
        elif not self.source_sample_id and self.variant.source_sample_id:
            self.source_sample_id = self.variant.source_sample_id

    def to_dict(self) -> Dict[str, Any]:
        variant_payload = self.variant.to_dict()
        if self.source_sample_id and not variant_payload.get('source_sample_id'):
            variant_payload['source_sample_id'] = self.source_sample_id

        payload: Dict[str, Any] = {
            'variant': variant_payload,
            'selected': self.selected_for_dataset,
        }

        if self.importance_score is not None:
            payload['importance_score'] = self.importance_score

        return payload

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CentralPoint":
        variant_data = data.get('variant', {})
        variant = GenomicVariant(
            chromosome=variant_data.get('chromosome'),
            position=variant_data.get('position'),
            ref_allele=variant_data.get('ref_allele'),
            alt_allele=variant_data.get('alt_allele'),
            quality=variant_data.get('quality', 0.0),
            depth=variant_data.get('depth', 0),
            allele_frequency=variant_data.get('allele_frequency', 0.0),
            filter_status=variant_data.get('filter_status', 'PASS'),
            variant_type=variant_data.get('variant_type', 'SNV'),
            source_sample_id=variant_data.get('source_sample_id'),
        )

        return cls(
            variant=variant,
            selected_for_dataset=data.get('selected', False),
            source_sample_id=data.get('source_sample_id', variant.source_sample_id),
            importance_score=data.get('importance_score'),
        )


@dataclass
class SequenceRecord:
    """Registro de sequência extraída."""
    sample_id: str
    central_point: CentralPoint
    sequence: str
    GRCh38_s: str  # Sequência original de referência (mesmo tamanho que sequence)
    fasta_file: Path
    label: int  # 1=longevo, 0=não-longevo
    genotype: str = "0/0"
    variant_present: bool = False

    def to_dict(self) -> Dict:
        return {
            'sample_id': self.sample_id,
            'variant': self.central_point.variant.to_dict(),
            'sequence': self.sequence,
            'GRCh38_s': self.GRCh38_s,
            'sequence_length': len(self.sequence),
            'label': self.label,
            'genotype': self.genotype,
            'variant_present': self.variant_present,
            'fasta_file': str(self.fasta_file),
        }

    @staticmethod
    def from_dict(data: Dict[str, Any]) -> "SequenceRecord":
        central_point_data = {
            'variant': data.get('variant', {}),
            'selected': data.get('selected', True),
            'source_sample_id': data.get(
                'central_point_source_sample_id',
                data.get('variant', {}).get('source_sample_id'),
            ),
        }

        central_point = CentralPoint.from_dict(central_point_data)

        fasta_path = data.get('fasta_file')
        if not fasta_path:
            raise ValueError('sequence record missing fasta_file path')
        if not data.get('sample_id'):
            raise ValueError('sequence record missing sample_id')

        return SequenceRecord(
            sample_id=data.get('sample_id'),
            central_point=central_point,
            sequence=data.get('sequence', ''),
            GRCh38_s=data.get('GRCh38_s', ''),
            fasta_file=Path(fasta_path),
            label=data.get('label', 0),
            genotype=data.get('genotype', '0/0'),
            variant_present=data.get('variant_present', False)
        )


@dataclass
class AlphaGenomeResult:
    """Resultado de predição do AlphaGenome."""
    sequence_record: SequenceRecord
    predictions: Dict[str, np.ndarray]  # output_type -> array de predições
    metadata: Dict[str, pd.DataFrame]  # output_type -> DataFrame de metadados
    statistics: Dict[str, Dict[str, float]] = field(default_factory=dict)
    feature_vector: Optional[np.ndarray] = None

    def aggregate_statistics(self) -> Dict[str, Dict[str, float]]:
        """Agrega estatísticas das predições."""
        if self.statistics:
            return self.statistics

        stats = {}
        for output_type, pred_array in self.predictions.items():
            if pred_array is None:
                continue
            try:
                arr = np.asarray(pred_array)
            except Exception:
                continue

            if arr.size == 0:
                continue

            stats[output_type] = {
                'mean': float(np.mean(arr)),
                'std': float(np.std(arr)),
                'min': float(np.min(arr)),
                'max': float(np.max(arr)),
                'median': float(np.median(arr))
            }

        self.statistics = stats
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

        self.config_dir = self.config_path.parent.resolve()
        
        # Diretórios: execução deve ocorrer dentro de DATA_ROOT
        self.work_dir = Path.cwd().resolve()
        self.data_root = DATA_ROOT

        if not self.data_root.exists():
            self.data_root.mkdir(parents=True, exist_ok=True)

        if self.work_dir != self.data_root:
            raise RuntimeError(
                "neural_longevity_dataset.py must be executed from /dados/GENOMICS_DATA/top3. "
                f"Current working directory: {self.work_dir}"
            )

        output_cfg = Path(self.config['project']['output_dir'])
        if output_cfg.is_absolute():
            output_resolved = output_cfg.resolve()
            if not str(output_resolved).startswith(str(self.data_root)):
                raise ValueError(
                    "project.output_dir must be inside /dados/GENOMICS_DATA/top3"
                )
            self.output_dir = output_resolved
        else:
            self.output_dir = (self.data_root / output_cfg).resolve()

        self.output_dir.mkdir(parents=True, exist_ok=True)

        console.print(f"[cyan]Diretório de trabalho: {self.work_dir}[/cyan]")
        console.print(f"[cyan]Saída em: {self.output_dir}[/cyan]")

        # Cache dir também relativo ao work_dir se for caminho relativo
        cache_path = Path(self.config['alphagenome']['cache_dir'])
        if cache_path.is_absolute():
            self.cache_dir = cache_path
        else:
            self.cache_dir = self.data_root / cache_path
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Lista ordenada das etapas do pipeline
        self._step_order = [
            'download_samples',
            'select_central_points',
            'extract_sequences',
            'run_alphagenome',
            'build_dataset'
        ]

        # Estado do processamento
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.state = self._load_checkpoint()
        self._variant_cache_runtime: Dict[str, Dict[Tuple[str, int], Dict[str, Any]]] = {}
        self._alphagenome_client = None
        self._chrom_lengths: Optional[Dict[str, int]] = None
        self._remote_size_cache: Dict[str, Optional[int]] = {}
        
        console.print(Panel.fit(
            f"[bold cyan]{self.config['project']['name']}[/bold cyan]\n"
            f"[dim]{self.config['project']['description']}[/dim]",
            border_style="cyan"
        ))
    
    def _load_checkpoint(self) -> Dict:
        """Carrega checkpoint se existir."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file) as f:
                state = json.load(f)
                return self._ensure_state_keys(state)
        return self._ensure_state_keys({
            'samples_downloaded': [],
            'variants_extracted': [],
            'central_points_selected': False,
            'sequences_extracted': [],
            'alphagenome_processed': [],
            'dataset_built': False
        })

    def _save_checkpoint(self):
        """Salva checkpoint."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.state, f, indent=2)

    # ───────────────────────────────────────────────────────────────
    # Configuração de etapas
    # ───────────────────────────────────────────────────────────────

    def override_steps(self, requested_steps: Sequence[str]):
        """Sobrescreve as etapas do pipeline conforme argumento CLI."""

        if not requested_steps:
            return

        valid_steps = set(self._step_order)
        normalized: List[str] = []
        for step in requested_steps:
            step_name = step.strip()
            if step_name not in valid_steps:
                valid_list = ', '.join(sorted(valid_steps))
                raise ValueError(
                    f"Invalid step '{step_name}'. Valid options: {valid_list}"
                )
            normalized.append(step_name)

        steps_cfg = self.config.setdefault('pipeline', {}).setdefault('steps', {})

        for step_name in self._step_order:
            steps_cfg[step_name] = False
        for step_name in normalized:
            steps_cfg[step_name] = True

        console.print(
            "[cyan]Executing subset of steps:[/cyan] " + ", ".join(normalized)
        )

    def _ensure_state_keys(self, state: Dict[str, Any]) -> Dict[str, Any]:
        """Garante que o estado contenha todas as chaves necessárias."""
        defaults = {
            'samples_downloaded': [],
            'variants_extracted': [],
            'central_points_selected': False,
            'sequences_extracted': [],
            'alphagenome_processed': [],
            'dataset_built': False,
            'sample_metadata': {},
            'vcf_paths': {},
            'variant_cache': {},
            'downloaded_chromosomes': [],  # NOVO: cromossomos baixados no modo multisample
            'vcf_chromosome_dir': None  # NOVO: diretório dos VCFs por cromossomo
        }

        for key, value in defaults.items():
            if key not in state:
                state[key] = value

        return state

    # ───────────────────────────────────────────────────────────────
    # Utilitários internos
    # ───────────────────────────────────────────────────────────────

    def _normalize_ena_url(self, url: str) -> str:
        """Normaliza URLs do ENA para um esquema HTTPS utilizável."""
        if not url:
            return url

        cleaned = url.strip()

        if cleaned.startswith('ftp://'):
            return cleaned.replace('ftp://', 'https://', 1)
        if cleaned.startswith('http://'):
            return cleaned.replace('http://', 'https://', 1)
        if cleaned.startswith('https://'):
            return cleaned

        # Muitos registros "submitted_ftp" retornam sem esquema (ex.: ftp.sra...).
        # Nesses casos, forçamos HTTPS explícito.
        return f"https://{cleaned}"

    def _group_cram_crai_pairs(self, urls: Sequence[str]) -> List[Dict[str, str]]:
        """Agrupa URLs por par CRAM/CRAI compartilhando o mesmo basename."""

        groups: Dict[str, Dict[str, str]] = {}
        order: List[str] = []
        
        for url in urls:
            lowered = url.lower()
            
            # Verificar primeiro se termina com .cram.crai (arquivo de índice)
            if lowered.endswith('.cram.crai'):
                # Remover .crai (5 caracteres), deixando .cram no final
                base_key = url[:-5]
                if base_key not in groups:
                    groups[base_key] = {'basename': base_key.rsplit('/', 1)[-1]}
                    order.append(base_key)
                groups[base_key]['crai_url'] = url
            elif lowered.endswith('.cram'):
                # Arquivo CRAM
                base_key = url  # Manter .cram no base_key
                if base_key not in groups:
                    groups[base_key] = {'basename': base_key.rsplit('/', 1)[-1]}
                    order.append(base_key)
                groups[base_key]['cram_url'] = url
            elif lowered.endswith('.crai'):
                # Arquivo .crai sem .cram antes (caso raro)
                base_key = url[:-1] + 'm'  # Trocar .crai por .cram
                if base_key not in groups:
                    groups[base_key] = {'basename': base_key.rsplit('/', 1)[-1]}
                    order.append(base_key)
                groups[base_key]['crai_url'] = url

        pairs: List[Dict[str, str]] = []
        for base_key in order:
            entry = groups.get(base_key, {})
            cram = entry.get('cram_url')
            crai = entry.get('crai_url')
            if cram and crai:
                pairs.append({
                    'basename': entry.get('basename', base_key.rsplit('/', 1)[-1]),
                    'cram_url': cram,
                    'crai_url': crai,
                })

        return pairs

    def _is_chromosomal_shard(self, basename: str) -> bool:
        """Detecta shards por cromossomo a partir do nome base do arquivo."""

        lower = basename.lower()
        chromosome_patterns = [
            r'(?:^|[_.-])chr(?:[0-9]{1,2}|[xy]|mt|m)(?:[_.-]|$)',
            r'(?:^|[_.-])chrom(?:osome)?[_.-]?(?:[0-9]{1,2}|[xy]|mt)(?:[_.-]|$)',
            r'(?:^|[_.-])chr(?:un|random)(?:[_.-]|$)'
        ]

        for pattern in chromosome_patterns:
            if re.search(pattern, lower):
                return True

        return False

    def _fetch_remote_size(self, url: str) -> Optional[int]:
        """Obtém tamanho remoto via HEAD, com cache em memória."""

        if url in self._remote_size_cache:
            return self._remote_size_cache[url]

        try:
            from urllib.request import Request, urlopen
        except Exception:  # pragma: no cover - fallback improvável
            self._remote_size_cache[url] = None
            return None

        try:
            try:
                request = Request(url, method='HEAD')
            except TypeError:  # Python < 3.8 compat
                request = Request(url)

                def _get_method():
                    return 'HEAD'

                request.get_method = _get_method  # type: ignore[attr-defined]

            with urlopen(request) as response:  # type: ignore[call-arg]
                length = response.headers.get('Content-Length')
        except Exception:
            self._remote_size_cache[url] = None
            return None

        if not length:
            self._remote_size_cache[url] = None
            return None

        try:
            size = int(length)
        except (TypeError, ValueError):
            self._remote_size_cache[url] = None
            return None

        self._remote_size_cache[url] = size
        return size

    def _select_cram_crai_pair(
        self,
        pairs: Sequence[Dict[str, str]],
    ) -> Tuple[Optional[Dict[str, str]], bool]:
        """Seleciona o melhor par CRAM/CRAI disponível.

        Retorna uma tupla com o par selecionado e um sinalizador indicando se
        apenas shards por cromossomo estavam disponíveis.
        """

        if not pairs:
            return None, False

        whole_genome: List[Dict[str, str]] = []
        chrom_shards: List[Dict[str, str]] = []

        for pair in pairs:
            basename = pair.get('basename', '')
            if self._is_chromosomal_shard(basename):
                chrom_shards.append(pair)
            else:
                whole_genome.append(pair)

        def choose_best(candidates: Sequence[Dict[str, str]]) -> Dict[str, str]:
            if len(candidates) == 1:
                return candidates[0]

            selected = candidates[0]
            best_size = self._fetch_remote_size(selected['cram_url'])

            for candidate in candidates[1:]:
                size = self._fetch_remote_size(candidate['cram_url'])
                if size is None:
                    continue
                if best_size is None or size > best_size:
                    selected = candidate
                    best_size = size

            return selected

        if whole_genome:
            return choose_best(whole_genome), False

        return choose_best(chrom_shards or pairs), True

    def _get_processed_samples(self) -> set:
        """Retorna conjunto com IDs já processados."""
        processed = set()
        for entry in self.state.get('samples_downloaded', []):
            if isinstance(entry, dict):
                processed.add(entry.get('sample_id'))
            elif isinstance(entry, str):
                processed.add(entry)
        return processed

    def _fetch_project_runs(self, project_id: str, max_records: int, batch_size: int = 500) -> pd.DataFrame:
        """Obtém tabela de runs do ENA Portal API."""
        fields = "run_accession,sample_accession,sample_alias,submitted_format,submitted_ftp,first_public"
        collected: List[pd.DataFrame] = []
        total = 0
        offset = 0
        
        first_batch = True

        while total < max_records:
            limit = min(batch_size, max_records - total)
            query = quote(f'study_accession="{project_id}" AND submitted_format="cram"')
            url = (
                "https://www.ebi.ac.uk/ena/portal/api/search?"
                f"result=read_run&format=tsv&fields={fields}&limit={limit}&offset={offset}&query={query}"
            )
            
            # Imprimir URL da API na primeira consulta
            if first_batch:
                # Decodificar a URL para exibição legível (pode copiar e colar no browser)
                url_legivel = unquote(url)
                console.print(f"\n[bold cyan]Consultando API do ENA:[/bold cyan]")
                console.print(f"[cyan]URL (copie e cole no browser):[/cyan]")
                console.print(f"[blue]{url_legivel}[/blue]")
                console.print(f"\n[cyan]Projeto: {project_id}[/cyan]")
                console.print(f"[cyan]Formato: CRAM[/cyan]")
                console.print(f"[cyan]Registros solicitados: {max_records}[/cyan]\n")
                first_batch = False

            try:
                df = pd.read_csv(url, sep='\t')
            except Exception as exc:
                console.print(f"[red]✗ Erro ao consultar ENA: {exc}[/red]")
                break

            if df.empty:
                break

            collected.append(df)
            total += len(df)
            offset += len(df)

            if len(df) < limit:
                break

        if not collected:
            return pd.DataFrame(columns=fields.split(','))

        combined = pd.concat(collected, ignore_index=True)
        
        # CORREÇÃO: Usar run_accession quando sample_alias está vazio
        if 'sample_alias' in combined.columns and 'run_accession' in combined.columns:
            # Preencher sample_alias vazios com run_accession
            combined['sample_alias'] = combined['sample_alias'].fillna(combined['run_accession'])
            # Remover registros que não têm nem sample_alias nem run_accession
            combined = combined.dropna(subset=['sample_alias'])
            combined = combined.drop_duplicates(subset=['sample_alias'], keep='first')
        elif 'run_accession' in combined.columns:
            # Se não houver sample_alias, usar run_accession
            combined['sample_alias'] = combined['run_accession']
            combined = combined.dropna(subset=['sample_alias'])
            combined = combined.drop_duplicates(subset=['sample_alias'], keep='first')

        return combined

    def _prepare_sample_records(self, df: pd.DataFrame, sample_range: List[int], label: int) -> List[Dict[str, Any]]:
        """Seleciona registros de acordo com o intervalo configurado."""
        if df.empty:
            return []

        start, end = sample_range
        subset = df.iloc[start:end]
        records: List[Dict[str, Any]] = []

        for idx, row in subset.iterrows():
            sample_alias = str(row.get('sample_alias') or '').strip()
            if not sample_alias:
                continue

            submitted = str(row.get('submitted_ftp') or '').strip()
            if not submitted:
                continue

            urls = [self._normalize_ena_url(part.strip()) for part in submitted.split(';') if part.strip()]
            pairs = self._group_cram_crai_pairs(urls)
            selected_pair, only_shards = self._select_cram_crai_pair(pairs)

            if not selected_pair:
                continue

            if only_shards:
                run_accession = str(row.get('run_accession', '')).strip()
                chosen_name = selected_pair.get('basename', 'desconhecido')
                console.print(
                    "[yellow]⚠ Nenhum CRAM de genoma completo para amostra"
                    f" {sample_alias} (run {run_accession}). Utilizando shard {chosen_name}.[/yellow]"
                )

            records.append({
                'sample_id': sample_alias,
                'run_accession': str(row.get('run_accession', '')).strip(),
                'label': label,
                'cram_url': selected_pair['cram_url'],
                'crai_url': selected_pair['crai_url']
            })

        return records

    def _download_file(self, url: str, destination: Path) -> bool:
        """Baixa arquivo exibindo progresso detalhado (idempotente)."""
        destination.parent.mkdir(parents=True, exist_ok=True)

        if destination.exists() and destination.stat().st_size > 0:
            console.print(
                f"[green]• Reutilizando download existente:[/green] "
                f"{destination.name} ({destination.stat().st_size / (1024**2):.1f} MB)"
            )
            return True

        console.print(
            f"[cyan]⇣ Baixando {destination.name}[/cyan]\n"
            f"    Origem: {url}\n"
            f"    Destino: {destination}"
        )

        try:
            from urllib.request import urlopen
        except Exception as exc:  # pragma: no cover - fallback improvável
            console.print(
                "[yellow]⚠ urllib indisponível, tentando wget clássico (sem progresso)[/yellow]"
            )
            try:
                sp.run(['wget', '-O', str(destination), url], check=True)
                return True
            except Exception as wget_exc:
                console.print(f"[yellow]⚠ Falha ao baixar {url}: {wget_exc}[/yellow]")
                if destination.exists():
                    destination.unlink(missing_ok=True)
                return False

        try:
            response = urlopen(url)
            total = int(response.headers.get('Content-Length', 0)) or None
        except Exception as exc:
            console.print(f"[yellow]⚠ Falha ao iniciar download de {url}: {exc}[/yellow]")
            return False

        chunk_size = 1024 * 1024  # 1 MiB

        try:
            with response, open(destination, 'wb') as fh, Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                DownloadColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
                console=console,
            ) as progress:
                task_description = f"Baixando {destination.name}"
                task = progress.add_task(task_description, total=total)

                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    fh.write(chunk)
                    progress.update(task, advance=len(chunk))

            # Rich Progress garante flush na conclusão; informar tamanho final
            size_mb = destination.stat().st_size / (1024 ** 2)
            console.print(f"[green]✓ Download concluído:[/green] {destination.name} ({size_mb:.1f} MB)")
            return True
        except Exception as exc:
            console.print(f"[yellow]⚠ Falha ao baixar {url}: {exc}[/yellow]")
            if destination.exists():
                destination.unlink(missing_ok=True)
            return False

    def _ensure_reference_indices(self, fasta_path: Path):
        """Garante que índices do FASTA existam."""
        fai_path = fasta_path.with_suffix(fasta_path.suffix + '.fai')
        if not fai_path.exists():
            console.print(f"[cyan]Indexando referência: {fasta_path}[/cyan]")
            sp.run(['samtools', 'faidx', str(fasta_path)], check=True)
    
    def _download_reference_genome(self, target_path: Path) -> bool:
        """
        Baixa o genoma de referência GRCh38 do 1000 Genomes se não existir.
        
        Args:
            target_path: Caminho onde salvar o arquivo FASTA
            
        Returns:
            True se sucesso, False caso contrário
        """
        console.print("\n" + "="*70)
        console.print("[bold yellow]GENOMA DE REFERÊNCIA NÃO ENCONTRADO[/bold yellow]")
        console.print("="*70 + "\n")
        
        console.print("[cyan]O arquivo de referência configurado não existe.[/cyan]")
        console.print("[cyan]Iniciando download automático do GRCh38 (1000 Genomes)...[/cyan]\n")
        
        # URL oficial do 1000 Genomes para GRCh38
        base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
        fasta_url = base_url + "GRCh38_full_analysis_set_plus_decoy_hla.fa"
        fai_url = fasta_url + ".fai"
        
        console.print("[bold cyan]Informações do Arquivo:[/bold cyan]")
        console.print(f"  • Nome: GRCh38_full_analysis_set_plus_decoy_hla.fa")
        console.print(f"  • Fonte: 1000 Genomes Project")
        console.print(f"  • URL: [blue]{fasta_url}[/blue]")
        console.print(f"  • Tamanho: ~3.2 GB (compactado)")
        console.print(f"  • Destino: {target_path}")
        console.print(f"  • Tempo estimado: 5-30 minutos (dependendo da conexão)\n")
        
        console.print("[yellow]⚠ Este é um arquivo grande. O download pode levar algum tempo.[/yellow]\n")
        
        # Criar diretório se não existir
        target_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Baixar FASTA
        console.print("[bold cyan]Etapa 1/2: Baixando arquivo FASTA...[/bold cyan]")
        fasta_ok = self._download_file(fasta_url, target_path)
        
        if not fasta_ok:
            console.print("[red]✗ Falha ao baixar o genoma de referência.[/red]")
            return False
        
        # Baixar índice .fai
        console.print("\n[bold cyan]Etapa 2/2: Baixando índice FASTA (.fai)...[/bold cyan]")
        fai_path = target_path.with_suffix(target_path.suffix + '.fai')
        fai_ok = self._download_file(fai_url, fai_path)
        
        if not fai_ok:
            console.print("[yellow]⚠ Índice não disponível para download, criando localmente...[/yellow]")
            try:
                self._ensure_reference_indices(target_path)
                console.print("[green]✓ Índice criado com sucesso[/green]")
            except Exception as e:
                console.print(f"[red]✗ Falha ao criar índice: {e}[/red]")
                return False
        
        console.print("\n" + "="*70)
        console.print("[green]✓ Genoma de referência baixado e indexado com sucesso![/green]")
        console.print("="*70 + "\n")
        
        return True

    def _get_reference_fasta(self) -> Optional[Path]:
        """Retorna caminho absoluto do FASTA de referência."""
        ref_path = Path(self.config['data_sources']['reference']['fasta'])

        candidates: List[Path] = []
        if ref_path.is_absolute():
            candidates.append(ref_path)
        else:
            # 1) Caminho relativo ao diretório obrigatório de dados
            candidates.append((DATA_ROOT / ref_path).resolve())
            # 2) Caminho relativo ao diretório do arquivo de configuração
            candidates.append((self.config_dir / ref_path).resolve())
            # 3) Caminho relativo ao diretório de trabalho atual (por segurança)
            candidates.append((Path.cwd() / ref_path).resolve())

        for candidate in candidates:
            if candidate.exists():
                self._ensure_reference_indices(candidate)
                return candidate

        # Nenhum candidato encontrado - tentar download automático
        tried_paths = "\n".join(str(path) for path in candidates) if candidates else str(ref_path)
        console.print("[yellow]⚠ Referência não encontrada. Caminhos verificados:[/yellow]\n" + tried_paths + "\n")
        
        # Usar o primeiro candidato como destino do download (geralmente DATA_ROOT)
        download_target = candidates[0] if candidates else (DATA_ROOT / ref_path).resolve()
        
        console.print(f"[cyan]Tentando baixar automaticamente para: {download_target}[/cyan]\n")
        
        if self._download_reference_genome(download_target):
            self._ensure_reference_indices(download_target)
            return download_target
        
        console.print("[red]✗ Não foi possível obter o genoma de referência.[/red]")
        return None

    def _generate_vcf_for_sample(self, sample_id: str, cram_path: Path, vcf_dir: Path, ref_fasta: Path) -> Optional[Path]:
        """Executa bcftools mpileup+call para gerar VCF."""
        vcf_dir.mkdir(parents=True, exist_ok=True)
        vcf_path = vcf_dir / f"{sample_id}.vcf.gz"

        if vcf_path.exists() and (vcf_path.with_suffix('.vcf.gz.tbi')).exists():
            return vcf_path

        console.print(f"[cyan]Chamando variantes para {sample_id}...[/cyan]")

        try:
            mpileup_cmd = ['bcftools', 'mpileup', '-Ou', '-f', str(ref_fasta), str(cram_path)]
            call_cmd = ['bcftools', 'call', '-mv', '-Oz', '-o', str(vcf_path)]

            mpileup_stderr = ''
            call_stderr = ''

            with sp.Popen(
                mpileup_cmd,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                text=True,
            ) as mpileup_proc:
                with sp.Popen(
                    call_cmd,
                    stdin=mpileup_proc.stdout,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    text=True,
                ) as call_proc:
                    mpileup_proc.stdout.close()
                    _, call_stderr = call_proc.communicate()
                    mpileup_stderr = mpileup_proc.stderr.read()
                    mpileup_return = mpileup_proc.wait()
                    call_return = call_proc.returncode

            if mpileup_stderr and mpileup_stderr.strip():
                console.print(
                    f"\n[yellow]⚠ bcftools mpileup warnings for {sample_id}:[/yellow]\n"
                    f"{mpileup_stderr.strip()}\n"
                )

            if call_stderr and call_stderr.strip():
                console.print(
                    f"\n[yellow]⚠ bcftools call warnings for {sample_id}:[/yellow]\n"
                    f"{call_stderr.strip()}\n"
                )

            if mpileup_return != 0 or call_return != 0:
                raise RuntimeError(
                    f"bcftools pipeline failed (mpileup exit {mpileup_return}, "
                    f"call exit {call_return})"
                )

            sp.run(['tabix', '-p', 'vcf', str(vcf_path)], check=True)
            return vcf_path
        except Exception as exc:
            console.print(f"[red]✗ Erro ao gerar VCF para {sample_id}: {exc}[/red]")
            if vcf_path.exists():
                vcf_path.unlink(missing_ok=True)
            return None

    def _update_sample_state(self, record: Dict[str, Any], cram_path: Path, crai_path: Path, vcf_path: Optional[Path]):
        """Atualiza metadados de uma amostra no checkpoint."""
        sample_id = record['sample_id']
        label = record['label']
        entry = {
            'sample_id': sample_id,
            'label': label,
            'run_accession': record.get('run_accession'),
            'cram_path': str(cram_path),
            'crai_path': str(crai_path),
            'vcf_path': str(vcf_path) if vcf_path else None
        }

        processed = self._get_processed_samples()
        if sample_id not in processed:
            self.state['samples_downloaded'].append(entry)

        self.state['sample_metadata'][sample_id] = entry
        if vcf_path:
            self.state['vcf_paths'][sample_id] = str(vcf_path)
            if sample_id not in self.state['variants_extracted']:
                self.state['variants_extracted'].append(sample_id)

    def _write_sample_list(self, samples: List[Dict[str, Any]], filename: str):
        """Salva lista de IDs de amostras."""
        sample_ids = [item['sample_id'] for item in samples]
        (self.output_dir / filename).write_text('\n'.join(sample_ids))

    def _create_placeholder_samples(self, filename: str, count: int = 50):
        """Cria lista simulada caso não seja possível baixar dados reais."""
        samples = [f"SIMULATED_{i:04d}" for i in range(count)]
        (self.output_dir / filename).write_text('\n'.join(samples))
        console.print(f"[yellow]⚠ Lista simulada gerada em {filename}[/yellow]")

    def _load_variants_for_sample(self, sample_id: str, vcf_path: Path) -> Dict[Tuple[str, int], Dict[str, Any]]:
        """Carrega variantes de um VCF (cacheado em memória)."""
        if sample_id in self._variant_cache_runtime:
            return self._variant_cache_runtime[sample_id]

        variants: Dict[Tuple[str, int], Dict[str, Any]] = {}
        if not vcf_path or not vcf_path.exists():
            self._variant_cache_runtime[sample_id] = variants
            return variants

        open_func = gzip.open if vcf_path.suffix == '.gz' else open

        try:
            with open_func(vcf_path, 'rt') as handle:
                for line in handle:
                    if not line or line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 10:
                        continue

                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alts = fields[4].split(',') if fields[4] else []
                    qual = float(fields[5]) if fields[5] not in {'.', ''} else None
                    filter_status = fields[6]
                    info_pairs = [item.split('=', 1) for item in fields[7].split(';') if '=' in item]
                    info_dict = {k: v for k, v in info_pairs}
                    format_keys = fields[8].split(':')
                    sample_values = fields[9].split(':')
                    format_dict = dict(zip(format_keys, sample_values))

                    af_raw = info_dict.get('AF', '0')
                    try:
                        allele_frequency = float(af_raw.split(',')[0])
                    except ValueError:
                        allele_frequency = 0.0

                    try:
                        depth = int(format_dict.get('DP', '0'))
                    except ValueError:
                        depth = 0

                    genotype = format_dict.get('GT', '0/0')

                    variants[(chrom, pos)] = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alts': alts,
                        'qual': qual,
                        'filter': filter_status,
                        'allele_frequency': allele_frequency,
                        'depth': depth,
                        'genotype': genotype
                    }
        except Exception as exc:
            console.print(f"[yellow]⚠ Falha ao ler VCF de {sample_id}: {exc}[/yellow]")

        self._variant_cache_runtime[sample_id] = variants
        return variants

    def _resolve_sample_allele(self, variant_data: Dict[str, Any], use_alt: bool) -> Tuple[str, str, bool]:
        """Determina alelo a ser usado com base no genótipo."""
        if not variant_data:
            return '', '0/0', False

        genotype = variant_data.get('genotype', '0/0')
        ref = variant_data.get('ref', '')
        alts = variant_data.get('alts', [])

        allele_indices = []
        for token in genotype.replace('|', '/').split('/'):
            try:
                allele_indices.append(int(token))
            except ValueError:
                continue

        has_alt = any(idx > 0 for idx in allele_indices)

        if use_alt and has_alt and alts:
            idx = next((i for i in allele_indices if i > 0), 1)
            alt_index = idx - 1 if idx - 1 < len(alts) else 0
            return alts[alt_index], genotype, True

        return ref, genotype, has_alt

    def _apply_variant_to_sequence(self, sequence: str, ref: str, alt: str, center_idx: int, window_size: int) -> str:
        """Aplica substituição de alelo mantendo comprimento da janela."""
        if not alt:
            alt = ref

        if center_idx < 0 or center_idx >= len(sequence):
            seq = sequence[:window_size]
            if len(seq) < window_size:
                seq += 'N' * (window_size - len(seq))
            return seq

        ref_len = len(ref) if ref else 1
        alt_len = len(alt) if alt else ref_len

        end_idx = center_idx + ref_len
        ref_segment = sequence[center_idx:end_idx]

        new_seq = sequence[:center_idx] + alt + sequence[end_idx:]

        if len(new_seq) > window_size:
            new_seq = new_seq[:window_size]
        elif len(new_seq) < window_size:
            deficit = window_size - len(new_seq)
            new_seq = new_seq + 'N' * deficit

        return new_seq

    def _get_alphagenome_analyzer(self, api_key: str):
        """Inicializa (se necessário) o cliente AlphaGenome."""
        if self._alphagenome_client is False:
            return None

        if self._alphagenome_client is None:
            try:
                from neural_module.neural_module import AlphaGenomeAnalyzer, DEFAULT_CONFIG
            except ImportError as exc:
                console.print(f"[yellow]⚠ AlphaGenome não disponível: {exc}[/yellow]")
                self._alphagenome_client = False
                return None

            config = DEFAULT_CONFIG.copy()
            config['default_outputs'] = self.config['alphagenome']['outputs']
            config['save_metadata'] = False
            config['use_advanced_viz'] = False

            analyzer = AlphaGenomeAnalyzer(api_key, config)
            if not analyzer.initialize():
                self._alphagenome_client = False
                return None

            self._alphagenome_client = analyzer

        return self._alphagenome_client

    def _convert_output_to_array(self, output_data: Any) -> Optional[np.ndarray]:
        """Extrai array NumPy de um output do AlphaGenome."""
        if output_data is None:
            return None

        candidates = ['values', 'scores', 'data']
        for attr in candidates:
            value = getattr(output_data, attr, None)
            if value is None:
                continue
            try:
                arr = np.asarray(value)
                if arr.size > 0:
                    return arr
            except Exception:
                continue

        if hasattr(output_data, 'to_numpy'):
            try:
                arr = np.asarray(output_data.to_numpy())
                if arr.size > 0:
                    return arr
            except Exception:
                pass

        try:
            arr = np.asarray(output_data)
            if arr.size > 0:
                return arr
        except Exception:
            return None

        return None

    def _build_alphagenome_feature_vector(self, stats: Dict[str, Dict[str, float]]) -> np.ndarray:
        """Constrói vetor de features a partir das estatísticas configuradas."""
        outputs = self.config['alphagenome']['outputs']
        stat_order = self.config['dataset']['features']['alphagenome_predictions']['statistics']
        feature_values: List[float] = []

        for output_name in outputs:
            output_stats = stats.get(output_name, {})
            for stat_name in stat_order:
                value = output_stats.get(stat_name)
                feature_values.append(float(value) if value is not None else float('nan'))

        return np.asarray(feature_values, dtype=float)

    def _load_chrom_lengths(self) -> Dict[str, int]:
        """Carrega comprimentos dos cromossomos a partir do arquivo FAI."""
        if self._chrom_lengths is not None:
            return self._chrom_lengths

        ref_fasta = self._get_reference_fasta()
        if not ref_fasta:
            self._chrom_lengths = {}
            return {}

        fai_path = Path(str(ref_fasta) + '.fai')
        if not fai_path.exists():
            self._ensure_reference_indices(ref_fasta)

        lengths: Dict[str, int] = {}
        try:
            with open(fai_path) as fh:
                for line in fh:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        lengths[parts[0]] = int(parts[1])
        except Exception as exc:
            console.print(f"[yellow]⚠ Não foi possível carregar comprimentos cromossômicos: {exc}[/yellow]")

        self._chrom_lengths = lengths
        return lengths

    def _normalize_position(self, chrom: str, position: int) -> float:
        """Normaliza posição genômica para [0,1] usando comprimento do cromossomo."""
        lengths = self._load_chrom_lengths()
        chrom_len = lengths.get(chrom)
        if chrom_len is None and chrom.startswith('chr'):
            chrom_len = lengths.get(chrom.replace('chr', ''))
        if not chrom_len:
            return float(position)
        return float(position) / float(chrom_len)

    def _encode_sequence_one_hot(self, sequence: str) -> np.ndarray:
        """Codifica sequência em one-hot (A,C,G,T)."""
        mapping = {
            'A': np.array([1, 0, 0, 0], dtype=np.float32),
            'C': np.array([0, 1, 0, 0], dtype=np.float32),
            'G': np.array([0, 0, 1, 0], dtype=np.float32),
            'T': np.array([0, 0, 0, 1], dtype=np.float32),
        }
        seq_array = np.zeros((len(sequence), 4), dtype=np.float32)
        for idx, base in enumerate(sequence.upper()):
            seq_array[idx] = mapping.get(base, np.zeros(4, dtype=np.float32))
        return seq_array
    
    # ───────────────────────────────────────────────────────────────
    # Passo 1: Download de Amostras
    # ───────────────────────────────────────────────────────────────
    
    def download_samples(self):
        """
        Baixa dados genômicos conforme modo configurado.
        
        Modos suportados:
        - "vcf_multisample": Baixa VCFs filtrados por cromossomo (RECOMENDADO - 3,202 amostras)
        - "cram": Baixa CRAM/CRAI do ENA e gera VCFs individuais
        """
        console.print("\n[bold cyan]Passo 1: Download de Amostras[/bold cyan]")

        is_dry_run = self.config['debug']['dry_run']
        download_mode = self.config['data_sources'].get('download_mode', 'vcf_multisample')

        if is_dry_run:
            console.print("[yellow]⚠ Modo dry-run: simulando downloads[/yellow]")

        console.print(f"[cyan]Modo de download: [bold]{download_mode.upper()}[/bold][/cyan]\n")

        if download_mode == 'vcf_multisample':
            return self._download_samples_vcf_multisample_mode()
        elif download_mode == 'cram':
            longevous_range = self.config['data_sources']['longevous']['sample_range']
            non_longevous_range = self.config['data_sources']['non_longevous']['sample_range']
            return self._download_samples_cram_mode(longevous_range, non_longevous_range, is_dry_run)
        else:
            raise ValueError(f"Modo de download '{download_mode}' não reconhecido. Use 'vcf_multisample' ou 'cram'.")
    
    def _download_samples_cram_mode(self, longevous_range, non_longevous_range, is_dry_run):
        """Modo CRAM: Baixa CRAM/CRAI do ENA e gera VCFs."""
        
        project_id = self.config['data_sources']['cram_source'].get('ena_project', 'PRJEB31736')

        max_records = max(longevous_range[1], non_longevous_range[1])
        runs_df = self._fetch_project_runs(project_id, max_records)
        
        if not runs_df.empty:
            console.print(f"[green]✓ {len(runs_df)} registros CRAM encontrados no ENA[/green]")
            console.print(f"[cyan]  Selecionando amostras longevas: índices {longevous_range[0]} a {longevous_range[1]-1}[/cyan]")
            console.print(f"[cyan]  Selecionando amostras não-longevas: índices {non_longevous_range[0]} a {non_longevous_range[1]-1}[/cyan]\n")

        if runs_df.empty:
            console.print("[red]✗ Não foi possível obter lista de amostras reais.[/red]")
            self._create_placeholder_samples("longevous_samples.txt", longevous_range[1] - longevous_range[0])
            self._create_placeholder_samples("non_longevous_samples.txt", non_longevous_range[1] - non_longevous_range[0])
            return

        longevous_records = self._prepare_sample_records(runs_df, longevous_range, label=1)
        non_longevous_records = self._prepare_sample_records(runs_df, non_longevous_range, label=0)

        console.print(f"[green]✓ {len(longevous_records)} registros longevos selecionados[/green]")
        console.print(f"[green]✓ {len(non_longevous_records)} registros não-longevos selecionados[/green]")
        
        # Seção informativa sobre os dados
        all_records_preview = longevous_records + non_longevous_records
        if all_records_preview:
            console.print("\n" + "="*70)
            console.print("[bold yellow]INFORMAÇÕES SOBRE OS DADOS A SEREM BAIXADOS[/bold yellow]")
            console.print("="*70 + "\n")
            
            console.print("[bold cyan]Origem dos Dados:[/bold cyan]")
            console.print("  • Projeto: [bold]1000 Genomes High Coverage[/bold]")
            console.print("  • Repositório: European Nucleotide Archive (ENA)")
            console.print(f"  • Código do projeto: [bold]{project_id}[/bold]")
            console.print(f"  • URL do projeto: [blue]https://www.ebi.ac.uk/ena/browser/view/{project_id}[/blue]\n")
            
            console.print("[bold cyan]Formato dos Arquivos:[/bold cyan]")
            console.print("  • Tipo: CRAM (Compressed Reference-oriented Alignment Map)")
            console.print("  • Conteúdo: Sequenciamento completo do genoma (WGS) de indivíduos")
            console.print("  • Cobertura: ~30x (alta qualidade)")
            console.print("  • Uso: Extração de variantes genéticas (SNPs, indels) de todos os cromossomos\n")
            
            console.print("[bold cyan]Estatísticas dos Downloads:[/bold cyan]")
            console.print(f"  • Total de indivíduos: [bold]{len(all_records_preview)}[/bold]")
            console.print(f"  • Longevos (controle positivo): {len(longevous_records)}")
            console.print(f"  • Não-longevos (controle negativo): {len(non_longevous_records)}")
            console.print("  • Tamanho médio por indivíduo: ~15-20 GB (CRAM + CRAI)")
            console.print(f"  • Espaço em disco necessário: ~{len(all_records_preview) * 18} GB\n")
            
            console.print("="*70)
            console.print("[bold yellow]DETALHES DOS INDIVÍDUOS SELECIONADOS[/bold yellow]")
            console.print("="*70 + "\n")
            
            for i, record in enumerate(all_records_preview, 1):
                label_text = "LONGEVO" if record['label'] == 1 else "NÃO-LONGEVO"
                run_acc = record.get('run_accession', 'N/A')
                
                console.print(f"[bold cyan]┌─ Indivíduo {i}/{len(all_records_preview)} ({label_text})[/bold cyan]")
                console.print(f"[cyan]│[/cyan]")
                console.print(f"[cyan]│[/cyan] [bold]Identificadores:[/bold]")
                console.print(f"[cyan]│[/cyan]   • Sample ID (ENA): {record['sample_id']}")
                console.print(f"[cyan]│[/cyan]   • Run Accession: {run_acc}")
                console.print(f"[cyan]│[/cyan]")
                console.print(f"[cyan]│[/cyan] [bold]Inspecionar no ENA:[/bold]")
                console.print(f"[cyan]│[/cyan]   • Página do run: [blue]https://www.ebi.ac.uk/ena/browser/view/{run_acc}[/blue]")
                console.print(f"[cyan]│[/cyan]   • Página da amostra: [blue]https://www.ebi.ac.uk/ena/browser/view/{record['sample_id']}[/blue]")
                console.print(f"[cyan]│[/cyan]")
                console.print(f"[cyan]│[/cyan] [bold]URLs de Download:[/bold]")
                console.print(f"[cyan]│[/cyan]   • CRAM: [dim]{record['cram_url']}[/dim]")
                console.print(f"[cyan]│[/cyan]   • CRAI: [dim]{record['crai_url']}[/dim]")
                console.print(f"[cyan]└{'─'*68}[/cyan]\n")
            
            console.print("="*70)
            console.print("[bold yellow]PRÓXIMOS PASSOS[/bold yellow]")
            console.print("="*70)
            console.print("1. Download dos arquivos CRAM e CRAI (~20 GB por indivíduo)")
            console.print("2. Chamada de variantes usando bcftools mpileup + call")
            console.print("3. Geração de arquivos VCF com todas as variantes genéticas")
            console.print("4. Seleção de pontos centrais (variantes de interesse)")
            console.print("5. Extração de sequências FASTA ao redor das variantes")
            console.print("6. Processamento com AlphaGenome para predições epigenéticas")
            console.print("7. Construção do dataset PyTorch para treinamento\n")
            console.print("="*70 + "\n")
            
            console.print("[yellow]⏳ Iniciando downloads... (isto pode levar várias horas)[/yellow]\n")

        self._write_sample_list(longevous_records, "longevous_samples.txt")
        self._write_sample_list(non_longevous_records, "non_longevous_samples.txt")

        if is_dry_run:
            console.print("[yellow]⚠ Dry-run: download real e chamada de variantes foram pulados[/yellow]")
            return

        ref_fasta = self._get_reference_fasta()
        if not ref_fasta:
            return

        all_records = longevous_records + non_longevous_records
        processed = self._get_processed_samples()

        cram_root = self.output_dir / "cram"
        vcf_root = self.output_dir / "vcf"
        cram_root.mkdir(exist_ok=True)
        vcf_root.mkdir(exist_ok=True)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Processando amostras...", total=len(all_records))

            for record in all_records:
                sample_id = record['sample_id']
                label = record['label']
                label_dir = 'longevous' if label == 1 else 'non_longevous'

                cram_dir = cram_root / label_dir
                vcf_dir = vcf_root / label_dir
                cram_dir.mkdir(exist_ok=True)
                vcf_dir.mkdir(exist_ok=True)

                cram_path = cram_dir / f"{sample_id}.cram"
                crai_path = cram_dir / f"{sample_id}.cram.crai"
                existing_vcf = vcf_dir / f"{sample_id}.vcf.gz"
                existing_vcf_index = existing_vcf.with_suffix('.vcf.gz.tbi')

                if (
                    sample_id in processed
                    and cram_path.exists()
                    and crai_path.exists()
                    and existing_vcf.exists()
                    and existing_vcf_index.exists()
                ):
                    self._update_sample_state(record, cram_path, crai_path, existing_vcf)
                    progress.advance(task)
                    continue

                if sample_id not in processed:
                    cram_ok = self._download_file(record['cram_url'], cram_path)
                    crai_ok = self._download_file(record['crai_url'], crai_path)
                else:
                    cram_ok = cram_path.exists()
                    crai_ok = crai_path.exists()

                vcf_path = None
                if cram_ok and crai_ok:
                    vcf_path = self._generate_vcf_for_sample(sample_id, cram_path, vcf_dir, ref_fasta)
                else:
                    console.print(f"[yellow]⚠ Pulando chamada de variantes para {sample_id} (download incompleto)[/yellow]")

                self._update_sample_state(record, cram_path, crai_path, vcf_path)
                processed.add(sample_id)
                progress.advance(task)

        self._save_checkpoint()
        console.print("[green]✓ Download e preparação das amostras concluídos[/green]")
    
    def _download_samples_vcf_mode(self, longevous_range, non_longevous_range, is_dry_run):
        """Modo VCF: Baixa VCFs já processados do IGSR (1000 Genomes)."""
        
        vcf_config = self.config['data_sources']['vcf_source']
        base_url = vcf_config['base_url']
        vcf_suffix = vcf_config['vcf_suffix']
        index_suffix = vcf_config['index_suffix']
        
        # Primeiro, precisamos obter a lista de sample IDs disponíveis
        # Para o modo VCF, vamos consultar o ENA para mapear índices para sample IDs
        project_id = self.config['data_sources']['cram_source'].get('ena_project', 'PRJEB31736')
        max_records = max(longevous_range[1], non_longevous_range[1])
        
        console.print(f"[cyan]Consultando ENA para mapear sample IDs...[/cyan]")
        runs_df = self._fetch_project_runs(project_id, max_records)
        
        if runs_df.empty:
            console.print("[red]✗ Não foi possível obter lista de amostras do ENA.[/red]")
            return
        
        console.print(f"[green]✓ {len(runs_df)} amostras encontradas no ENA[/green]")
        console.print(f"[cyan]  Selecionando amostras longevas: índices {longevous_range[0]} a {longevous_range[1]-1}[/cyan]")
        console.print(f"[cyan]  Selecionando amostras não-longevas: índices {non_longevous_range[0]} a {non_longevous_range[1]-1}[/cyan]\n")
        
        # Preparar registros
        longevous_records = self._prepare_vcf_records(runs_df, longevous_range, label=1, base_url=base_url, vcf_suffix=vcf_suffix, index_suffix=index_suffix)
        non_longevous_records = self._prepare_vcf_records(runs_df, non_longevous_range, label=0, base_url=base_url, vcf_suffix=vcf_suffix, index_suffix=index_suffix)
        
        console.print(f"[green]✓ {len(longevous_records)} registros longevos preparados[/green]")
        console.print(f"[green]✓ {len(non_longevous_records)} registros não-longevos preparados[/green]")
        
        # Seção informativa
        all_records = longevous_records + non_longevous_records
        if all_records:
            console.print("\n" + "="*70)
            console.print("[bold yellow]INFORMAÇÕES SOBRE OS VCFs A SEREM BAIXADOS[/bold yellow]")
            console.print("="*70 + "\n")
            
            console.print("[bold cyan]Origem dos Dados:[/bold cyan]")
            console.print("  • Projeto: [bold]1000 Genomes High Coverage[/bold]")
            console.print("  • Repositório: IGSR (International Genome Sample Resource)")
            console.print(f"  • URL base: [blue]{base_url}[/blue]\n")
            
            console.print("[bold cyan]Formato dos Arquivos:[/bold cyan]")
            console.print("  • Tipo: VCF 4.2 (Variant Call Format)")
            console.print("  • Processamento: GATK HaplotypeCaller + VQSR")
            console.print("  • Conteúdo: Variantes genéticas já chamadas e filtradas")
            console.print("  • Qualidade: Alta (pipeline GATK best-practices)\n")
            
            console.print("[bold cyan]Estatísticas dos Downloads:[/bold cyan]")
            console.print(f"  • Total de indivíduos: [bold]{len(all_records)}[/bold]")
            console.print(f"  • Longevos: {len(longevous_records)}")
            console.print(f"  • Não-longevos: {len(non_longevous_records)}")
            console.print("  • Tamanho médio por VCF: ~200-500 MB (comprimido)")
            console.print(f"  • Espaço em disco necessário: ~{len(all_records) * 0.35:.1f} GB")
            console.print("  • [bold green]Vantagem vs CRAM: 100x mais rápido, 100x menor![/bold green]\n")
            
            console.print("="*70)
            console.print("[bold yellow]DETALHES DOS VCFs SELECIONADOS[/bold yellow]")
            console.print("="*70 + "\n")
            
            for i, record in enumerate(all_records, 1):
                label_text = "LONGEVO" if record['label'] == 1 else "NÃO-LONGEVO"
                run_acc = record.get('run_accession', 'N/A')
                
                console.print(f"[bold cyan]┌─ Indivíduo {i}/{len(all_records)} ({label_text})[/bold cyan]")
                console.print(f"[cyan]│[/cyan]")
                console.print(f"[cyan]│[/cyan] [bold]Identificadores:[/bold]")
                console.print(f"[cyan]│[/cyan]   • Sample ID: {record['sample_id']}")
                console.print(f"[cyan]│[/cyan]   • Run Accession: {run_acc}")
                console.print(f"[cyan]│[/cyan]")
                console.print(f"[cyan]│[/cyan] [bold]Inspecionar no ENA:[/bold]")
                console.print(f"[cyan]│[/cyan]   • Página do run: [blue]https://www.ebi.ac.uk/ena/browser/view/{run_acc}[/blue]")
                console.print(f"[cyan]│[/cyan]   • Página da amostra: [blue]https://www.ebi.ac.uk/ena/browser/view/{record['sample_id']}[/blue]")
                console.print(f"[cyan]│[/cyan]")
                console.print(f"[cyan]│[/cyan] [bold]URLs de Download:[/bold]")
                console.print(f"[cyan]│[/cyan]   • VCF: [dim]{record['vcf_url']}[/dim]")
                console.print(f"[cyan]│[/cyan]   • Index: [dim]{record['index_url']}[/dim]")
                console.print(f"[cyan]└{'─'*68}[/cyan]\n")
            
            console.print("="*70)
            console.print("[bold yellow]VANTAGENS DO MODO VCF[/bold yellow]")
            console.print("="*70)
            console.print("✓ Download 100x mais rápido (~30 min vs ~20 horas)")
            console.print("✓ Arquivos 100x menores (~200 MB vs ~20 GB por indivíduo)")
            console.print("✓ Variantes já chamadas com GATK (qualidade superior)")
            console.print("✓ VQSR aplicado (filtragem de qualidade profissional)")
            console.print("✓ Economia de tempo e espaço em disco")
            console.print("✓ Pipeline validado pelo consórcio 1000 Genomes\n")
            console.print("="*70 + "\n")
            
            console.print("[yellow]⏳ Iniciando downloads de VCFs... (muito mais rápido que CRAMs!)[/yellow]\n")
        
        self._write_sample_list(longevous_records, "longevous_samples.txt")
        self._write_sample_list(non_longevous_records, "non_longevous_samples.txt")
        
        if is_dry_run:
            console.print("[yellow]⚠ Dry-run: download real foi pulado[/yellow]")
            return
        
        # Baixar VCFs
        vcf_root = self.output_dir / "vcf"
        vcf_root.mkdir(exist_ok=True)
        
        processed = self._get_processed_samples()
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task("Baixando VCFs...", total=len(all_records))
            
            for record in all_records:
                sample_id = record['sample_id']
                label = record['label']
                label_dir = 'longevous' if label == 1 else 'non_longevous'
                
                vcf_dir = vcf_root / label_dir
                vcf_dir.mkdir(exist_ok=True)
                
                vcf_path = vcf_dir / f"{sample_id}.vcf.gz"
                vcf_index_path = vcf_dir / f"{sample_id}.vcf.gz.tbi"
                
                # Verificar se já existe
                if sample_id in processed and vcf_path.exists() and vcf_index_path.exists():
                    self._update_sample_state_vcf(record, vcf_path)
                    progress.advance(task)
                    continue
                
                # Baixar VCF e índice
                vcf_ok = self._download_file(record['vcf_url'], vcf_path)
                index_ok = self._download_file(record['index_url'], vcf_index_path)
                
                if vcf_ok and index_ok:
                    self._update_sample_state_vcf(record, vcf_path)
                else:
                    console.print(f"[yellow]⚠ Falha no download de {sample_id}[/yellow]")
                
                processed.add(sample_id)
                progress.advance(task)
        
        self._save_checkpoint()
        console.print("[green]✓ Download de VCFs concluído![/green]")
    
    def _prepare_vcf_records(self, df: pd.DataFrame, sample_range: List[int], label: int, base_url: str, vcf_suffix: str, index_suffix: str) -> List[Dict[str, Any]]:
        """Prepara registros para download de VCFs do IGSR."""
        if df.empty:
            return []
        
        start, end = sample_range
        subset = df.iloc[start:end]
        records: List[Dict[str, Any]] = []
        
        for idx, row in subset.iterrows():
            sample_id = str(row.get('sample_alias') or '').strip()
            if not sample_id:
                continue
            
            run_acc = str(row.get('run_accession', '')).strip()
            
            # Construir URLs do IGSR
            vcf_url = base_url + vcf_suffix.format(sample_id=sample_id)
            index_url = base_url + index_suffix.format(sample_id=sample_id)
            
            records.append({
                'sample_id': sample_id,
                'run_accession': run_acc,
                'label': label,
                'vcf_url': vcf_url,
                'index_url': index_url
            })
        
        return records
    
    def _update_sample_state_vcf(self, record: Dict[str, Any], vcf_path: Path):
        """Atualiza estado de amostra no modo VCF."""
        sample_id = record['sample_id']
        label = record['label']
        
        entry = {
            'sample_id': sample_id,
            'label': label,
            'run_accession': record.get('run_accession'),
            'vcf_path': str(vcf_path),
            'download_mode': 'vcf'
        }
        
        processed = self._get_processed_samples()
        if sample_id not in processed:
            self.state['samples_downloaded'].append(entry)
        
        self.state['sample_metadata'][sample_id] = entry
        self.state['vcf_paths'][sample_id] = str(vcf_path)
        
        if sample_id not in self.state['variants_extracted']:
            self.state['variants_extracted'].append(sample_id)
    
    def _download_samples_vcf_multisample_mode(self):
        """Modo VCF Multi-sample: Baixa VCFs filtrados por cromossomo (todos os 3,202 indivíduos)."""
        
        vcf_config = self.config['data_sources']['vcf_multisample_source']
        base_url = vcf_config['base_url']
        filename_pattern = vcf_config['filename_pattern']
        chromosomes = vcf_config['chromosomes']
        
        console.print("\n" + "="*70)
        console.print("[bold yellow]MODO VCF MULTI-SAMPLE (CROMOSSOMOS COMPLETOS)[/bold yellow]")
        console.print("="*70 + "\n")
        
        console.print("[bold cyan]Origem dos Dados:[/bold cyan]")
        console.print("  • Projeto: [bold]1000 Genomes High Coverage (30x)[/bold]")
        console.print("  • Repositório: IGSR - Filtered & Phased Variants")
        console.print(f"  • URL base: [blue]{base_url}[/blue]")
        console.print(f"  • Total de amostras no VCF: [bold]3,202 indivíduos[/bold]")
        console.print(f"  • Cromossomos a baixar: {len(chromosomes)}\n")
        
        console.print("[bold cyan]Formato dos Arquivos:[/bold cyan]")
        console.print("  • Tipo: VCF 4.2 multi-sample (todas as amostras por cromossomo)")
        console.print("  • Processamento: GATK HaplotypeCaller + VQSR + Phased")
        console.print("  • Filtros aplicados: PASS, MAC≥2, HWE, Mendel errors")
        console.print("  • Conteúdo: ~73.5 milhões de variantes de alta qualidade")
        console.print("  • Qualidade: Alta (pipeline NYGC/Broad best-practices)\n")
        
        # Estimativa de tamanho
        avg_size_gb = {
            'chr1': 5.8, 'chr2': 5.2, 'chr3': 4.3, 'chr4': 4.1, 'chr5': 3.9,
            'chr6': 3.7, 'chr7': 3.5, 'chr8': 3.3, 'chr9': 2.8, 'chr10': 3.0,
            'chr11': 3.1, 'chr12': 2.9, 'chr13': 2.0, 'chr14': 2.0, 'chr15': 1.9,
            'chr16': 2.1, 'chr17': 1.9, 'chr18': 1.7, 'chr19': 1.5, 'chr20': 1.4,
            'chr21': 0.9, 'chr22': 0.9, 'chrX': 2.8
        }
        total_size = sum(avg_size_gb.get(c, 3.0) for c in chromosomes)
        
        console.print("[bold cyan]Estatísticas dos Downloads:[/bold cyan]")
        console.print(f"  • Cromossomos selecionados: [bold]{len(chromosomes)}[/bold]")
        console.print(f"  • Tamanho estimado total: ~{total_size:.1f} GB")
        console.print(f"  • Tempo estimado: ~{total_size * 2:.0f}-{total_size * 5:.0f} minutos")
        console.print(f"  • Variantes totais: ~73.5M (genoma), ~0.6M (chr21), ~0.7M (chr22)\n")
        
        console.print("="*70)
        console.print("[bold yellow]CROMOSSOMOS A BAIXAR[/bold yellow]")
        console.print("="*70 + "\n")
        
        vcf_dir = self.output_dir / "vcf_chromosomes"
        vcf_dir.mkdir(parents=True, exist_ok=True)
        
        downloaded_chroms = []
        
        for i, chrom in enumerate(chromosomes, 1):
            # Tratamento especial para chrX que tem .v2 no nome
            if chrom == 'chrX':
                filename = filename_pattern.replace('.vcf.gz', '.v2.vcf.gz').format(chrom=chrom)
            else:
                filename = filename_pattern.format(chrom=chrom)
            
            vcf_url = base_url + filename
            vcf_path = vcf_dir / filename
            index_path = vcf_path.with_suffix(vcf_path.suffix + '.tbi')
            
            size_est = avg_size_gb.get(chrom, 3.0)
            
            console.print(f"[bold cyan]┌─ Cromossomo {i}/{len(chromosomes)}: {chrom}[/bold cyan]")
            console.print(f"[cyan]│[/cyan]")
            console.print(f"[cyan]│[/cyan] [bold]Arquivo:[/bold] {filename}")
            console.print(f"[cyan]│[/cyan] [bold]URL:[/bold]")
            console.print(f"[cyan]│[/cyan]   [blue]{vcf_url}[/blue]")
            console.print(f"[cyan]│[/cyan] [bold]Tamanho estimado:[/bold] ~{size_est:.1f} GB")
            console.print(f"[cyan]│[/cyan] [bold]Destino:[/bold]")
            console.print(f"[cyan]│[/cyan]   {vcf_path}")
            console.print(f"[cyan]└{'─'*68}[/cyan]\n")
            
            # Verificar se já existe
            if vcf_path.exists() and index_path.exists():
                actual_size = vcf_path.stat().st_size / (1024**3)
                console.print(f"[green]• Reutilizando download existente: {chrom} ({actual_size:.2f} GB)[/green]\n")
                downloaded_chroms.append(chrom)
                continue
            
            # Baixar VCF
            console.print(f"[yellow]⏳ Baixando {chrom}...[/yellow]")
            vcf_ok = self._download_file(vcf_url, vcf_path)
            
            if vcf_ok:
                # Baixar índice
                index_url = vcf_url + '.tbi'
                index_ok = self._download_file(index_url, index_path)
                
                if index_ok:
                    console.print(f"[green]✓ Download completo: {chrom}[/green]\n")
                    downloaded_chroms.append(chrom)
                else:
                    console.print(f"[yellow]⚠ Índice não disponível para {chrom}, criando localmente...[/yellow]")
                    try:
                        sp.run(['tabix', '-p', 'vcf', str(vcf_path)], check=True)
                        console.print(f"[green]✓ Índice criado localmente para {chrom}[/green]\n")
                        downloaded_chroms.append(chrom)
                    except Exception as e:
                        console.print(f"[red]✗ Falha ao criar índice para {chrom}: {e}[/red]\n")
            else:
                console.print(f"[red]✗ Falha no download de {chrom}[/red]\n")
        
        # Salvar estado
        self.state['downloaded_chromosomes'] = downloaded_chroms
        self.state['vcf_chromosome_dir'] = str(vcf_dir)
        self._save_checkpoint()
        
        console.print("="*70)
        console.print(f"[green]✓ Download concluído: {len(downloaded_chroms)}/{len(chromosomes)} cromossomos[/green]")
        console.print("="*70 + "\n")
        
        if len(downloaded_chroms) < len(chromosomes):
            console.print(f"[yellow]⚠ Alguns cromossomos não foram baixados. Verifique os erros acima.[/yellow]\n")
        
        return downloaded_chroms
    
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
            # Montar filtros dinamicamente (nem todo VCF possui FORMAT/DP)
            filter_conditions: List[str] = []

            min_qual = filters.get('min_quality', 0)
            if min_qual:
                filter_conditions.append(f'QUAL>={min_qual}')

            min_depth = filters.get('min_depth', 0)
            if min_depth:
                if self._vcf_has_format_field(vcf_path, 'DP'):
                    filter_conditions.append(f'FMT/DP>={min_depth}')
                else:
                    console.print(
                        f"[yellow]⚠ FORMAT/DP ausente no header de {vcf_path}. "
                        "Ignorando filtro por profundidade.[/yellow]"
                    )

            # Usar bcftools para extrair variantes
            cmd = ['bcftools', 'view']

            if filters.get('filter_pass_only', False):
                # Muitos VCFs gerados via mpileup não anotam FILTER=PASS explicitamente,
                # mantendo o campo como '.' (não filtrado). Para contemplar esse caso,
                # aceitamos tanto PASS quanto '.' ao aplicar o filtro.
                cmd.extend(['-f', '.,PASS'])

            if filter_conditions:
                cmd.extend(['-i', ' && '.join(filter_conditions)])

            cmd.append(str(vcf_path))

            command_str = _format_command(cmd)
            console.print(f"[blue]$ {command_str}[/blue]")

            repair_attempted = False
            while True:
                try:
                    result = sp.run(cmd, capture_output=True, text=True, check=True)
                    break
                except sp.CalledProcessError as e:
                    stderr = e.stderr.strip() if e.stderr else str(e)
                    if (
                        not repair_attempted
                        and 'No BGZF EOF marker' in stderr
                    ):
                        repair_attempted = True
                        console.print(
                            "[yellow]⚠ bcftools detectou EOF ausente em BGZF. "
                            "Tentando recompactar o VCF automaticamente...[/yellow]"
                        )
                        fixed_path = self._repair_truncated_bgzf(vcf_path)
                        if fixed_path and fixed_path != vcf_path:
                            vcf_path = fixed_path
                        if fixed_path:
                            cmd[-1] = str(vcf_path)
                            continue
                    raise
            
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
                    variant_type=var_type,
                    source_sample_id=sample_id,
                )

                variants.append(variant)
        
        except sp.CalledProcessError as e:
            stderr = e.stderr.strip() if e.stderr else str(e)
            console.print(
                f"[red]✗ Erro ao extrair variantes de {vcf_path}: {stderr}[/red]"
            )
            console.print(f"[red]  Comando executado: {_format_command(cmd)}[/red]")
        except Exception as e:
            console.print(
                f"[red]✗ Erro inesperado ao extrair variantes de {vcf_path}: {e}[/red]"
            )

        return variants

    def _vcf_has_format_field(self, vcf_path: Path, field_id: str) -> bool:
        """Verifica se o VCF declara um campo FORMAT específico no header."""

        opener = gzip.open if vcf_path.suffix.endswith('gz') else open
        try:
            with opener(vcf_path, 'rt') as handle:
                for line in handle:
                    if not line.startswith('##'):
                        break
                    if line.startswith('##FORMAT=') and f'ID={field_id}' in line:
                        return True
        except OSError:
            # Caso não seja possível ler, deixe o bcftools reportar o problema
            return False

        return False

    def _repair_truncated_bgzf(self, vcf_path: Path) -> Optional[Path]:
        """Tenta recompactar um VCF que está sem o marcador EOF do BGZF."""

        if not vcf_path.exists():
            return None

        temp_plain: Optional[Path] = None
        try:
            with tempfile.NamedTemporaryFile(
                delete=False,
                suffix='.vcf',
                dir=str(vcf_path.parent)
            ) as tmp_handle:
                temp_plain = Path(tmp_handle.name)
                with gzip.open(vcf_path, 'rb') as src:
                    shutil.copyfileobj(src, tmp_handle)

            bgzip_cmd = ['bgzip', '-f', str(temp_plain)]
            sp.run(bgzip_cmd, check=True)

            recompressed_path = temp_plain.with_suffix(temp_plain.suffix + '.gz')

            if self.config.get('debug', {}).get('save_intermediate', False):
                backup_path = vcf_path.with_name(vcf_path.name + '.broken')
                try:
                    shutil.move(str(vcf_path), str(backup_path))
                except Exception:
                    vcf_path.unlink(missing_ok=True)
            else:
                vcf_path.unlink(missing_ok=True)

            shutil.move(str(recompressed_path), str(vcf_path))

            try:
                sp.run(['tabix', '-f', '-p', 'vcf', str(vcf_path)], check=False)
            except FileNotFoundError:
                pass

            console.print(f"[green]✓ VCF reparado: {vcf_path.name}[/green]")
            return vcf_path
        except Exception as exc:
            console.print(f"[red]✗ Falha ao recompactar {vcf_path}: {exc}[/red]")
            return None
        finally:
            if temp_plain and temp_plain.exists():
                temp_plain.unlink(missing_ok=True)
    
    # ───────────────────────────────────────────────────────────────
    # Passo 3: Seleção de Pontos Centrais
    # ───────────────────────────────────────────────────────────────
    
    def select_central_points(self) -> List[CentralPoint]:
        """
        Seleciona pontos centrais conforme modo de download.
        
        Returns:
            Lista de pontos centrais selecionados
        """
        download_mode = self.config['data_sources'].get('download_mode', 'vcf_multisample')
        
        if download_mode == 'vcf_multisample':
            return self._select_central_points_multisample()
        else:
            return self._select_central_points_individual()
    
    def _select_central_points_individual(self) -> List[CentralPoint]:
        """Método original de seleção de pontos centrais (modo CRAM/individual)."""
        console.print("\n[bold cyan]Passo 3: Seleção de Pontos Centrais (Individual)[/bold cyan]")
        
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
            
            # Obter VCF gerado automaticamente
            vcf_path_str = self.state['vcf_paths'].get(first_sample)
            if vcf_path_str:
                vcf_path = Path(vcf_path_str)
            else:
                vcf_path = self.output_dir / f"vcf/longevous/{first_sample}.vcf.gz"

            if not vcf_path.exists():
                console.print(f"[yellow]⚠ VCF não encontrado para {first_sample}: {vcf_path}[/yellow]")
                console.print(f"[yellow]  Criando pontos centrais simulados...[/yellow]")
                return self._create_simulated_central_points(n_points)
            
            # Extrair variantes
            variants = self.extract_variants(first_sample, vcf_path)
            
            console.print(f"[green]✓ {len(variants)} variantes encontradas[/green]")
            
            # Selecionar top N variantes (por qualidade)
            variants_sorted = sorted(variants, key=lambda v: v.quality, reverse=True)
            selected_variants = variants_sorted[:n_points]
            
            central_points = [
                CentralPoint(
                    variant=v,
                    selected_for_dataset=True,
                    source_sample_id=first_sample,
                )
                for v in selected_variants
            ]
            
            # Salvar pontos centrais
            self._save_central_points(central_points)
            
            console.print(f"[green]✓ {len(central_points)} pontos centrais selecionados[/green]")

            return central_points

        elif strategy == "random_rotation_longevous_samples":
            longevous_samples_file = self.output_dir / "longevous_samples.txt"

            if not longevous_samples_file.exists():
                console.print(f"[yellow]⚠ Arquivo não encontrado: {longevous_samples_file}[/yellow]")
                console.print(f"[yellow]  Criando pontos centrais simulados...[/yellow]")
                return self._create_simulated_central_points(n_points)

            with open(longevous_samples_file) as f:
                longevous_samples = [line.strip() for line in f if line.strip()]

            if not longevous_samples:
                console.print("[yellow]⚠ Nenhuma amostra longeva listada para rotação aleatória[/yellow]")
                console.print(f"[yellow]  Criando pontos centrais simulados...[/yellow]")
                return self._create_simulated_central_points(n_points)

            console.print(
                f"[cyan]Selecionando variantes aleatórias entre {len(longevous_samples)} longevos[/cyan]"
            )

            rng_seed = self.config['variant_selection'].get('random_seed')
            rng = random.Random(rng_seed)

            sample_variants_pool: Dict[str, List[GenomicVariant]] = {}
            exhausted_samples: Set[str] = set()
            central_points: List[CentralPoint] = []

            total_samples = len(longevous_samples)
            sample_index = 0
            central_points_counter = 0

            while central_points_counter < n_points:
                if len(exhausted_samples) == total_samples:
                    break

                sample_id = longevous_samples[sample_index]

                # Avança até encontrar um longevo ainda disponível
                cycle_advance = 0
                while sample_id in exhausted_samples and cycle_advance < total_samples:
                    sample_index = (sample_index + 1) % total_samples
                    sample_id = longevous_samples[sample_index]
                    cycle_advance += 1

                if cycle_advance >= total_samples and sample_id in exhausted_samples:
                    break

                if sample_id not in sample_variants_pool:
                    vcf_path_str = self.state['vcf_paths'].get(sample_id)
                    if vcf_path_str:
                        vcf_path = Path(vcf_path_str)
                    else:
                        vcf_path = self.output_dir / f"vcf/longevous/{sample_id}.vcf.gz"

                    if not vcf_path.exists():
                        console.print(
                            f"[yellow]⚠ VCF não encontrado para {sample_id}: {vcf_path}[/yellow]"
                        )
                        exhausted_samples.add(sample_id)
                        sample_index = (sample_index + 1) % total_samples
                        continue

                    variants = self.extract_variants(sample_id, vcf_path)
                    if not variants:
                        console.print(
                            f"[yellow]⚠ Nenhuma variante elegível para {sample_id}[/yellow]"
                        )
                        exhausted_samples.add(sample_id)
                        sample_index = (sample_index + 1) % total_samples
                        continue

                    sample_variants_pool[sample_id] = variants

                variants_pool = sample_variants_pool.get(sample_id, [])
                if not variants_pool:
                    exhausted_samples.add(sample_id)
                    sample_index = (sample_index + 1) % total_samples
                    continue

                variant_index = rng.randrange(len(variants_pool))
                variant = variants_pool.pop(variant_index)
                central_points.append(
                    CentralPoint(
                        variant=variant,
                        selected_for_dataset=True,
                        source_sample_id=sample_id,
                    )
                )

                central_points_counter += 1

                if not variants_pool:
                    exhausted_samples.add(sample_id)

                sample_index = (sample_index + 1) % total_samples

            if central_points_counter < n_points:
                remaining = n_points - central_points_counter
                console.print(
                    f"[yellow]⚠ Apenas {central_points_counter} pontos centrais reais selecionados."
                    f" Completando com {remaining} pontos simulados.[/yellow]"
                )
                central_points.extend(self._create_simulated_central_points(remaining))

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
                CentralPoint(
                    variant=variant,
                    selected_for_dataset=True,
                    source_sample_id=None,
                )
            )
        return central_points
    
    def _select_central_points_multisample(self) -> List[CentralPoint]:
        """
        Seleciona pontos centrais dos VCFs multi-sample por cromossomo.
        
        Estratégia: Percorre cada cromossomo e seleciona N variantes 
        sequencialmente ou aleatoriamente.
        """
        console.print("\n[bold cyan]Passo 2: Seleção de Pontos Centrais (Multi-sample)[/bold cyan]")
        
        variant_config = self.config['data_sources']['variant_sampling']
        n_per_chrom = variant_config['n_variants_per_chromosome']
        selection_method = variant_config.get('selection_method', 'random')
        random_seed = variant_config.get('random_seed', 42)
        min_quality = variant_config.get('min_quality', 30)
        only_pass = variant_config.get('only_pass', True)
        show_counts = variant_config.get('show_variant_counts', False)
        
        downloaded_chroms = self.state.get('downloaded_chromosomes', [])
        vcf_dir = Path(self.state.get('vcf_chromosome_dir', self.output_dir / 'vcf_chromosomes'))
        
        if not downloaded_chroms:
            console.print("[red]✗ Nenhum cromossomo baixado. Execute download_samples primeiro.[/red]")
            return []
        
        console.print(f"[cyan]Método de seleção: [bold]{selection_method}[/bold][/cyan]")
        console.print(f"[cyan]Variantes por cromossomo: [bold]{n_per_chrom}[/bold][/cyan]")
        console.print(f"[cyan]Nota: VCFs do IGSR já são filtrados (alta qualidade)[/cyan]\n")
        
        # Contar variantes disponíveis em cada cromossomo (opcional)
        if show_counts:
            console.print("[bold cyan]Contando variantes disponíveis em cada cromossomo...[/bold cyan]")
            total_variants_available = 0
            variant_counts = {}
            
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                console=console
            ) as progress:
                count_task = progress.add_task("Contando variantes...", total=len(downloaded_chroms))
                
                for chrom in downloaded_chroms:
                    vcf_files = list(vcf_dir.glob(f"*{chrom}.*.vcf.gz"))
                    if vcf_files:
                        vcf_path = vcf_files[0]
                        try:
                            # Contar linhas do VCF (excluindo cabeçalho)
                            count_cmd = f"bcftools view -H {vcf_path} | wc -l"
                            result = sp.run(count_cmd, capture_output=True, text=True, shell=True, check=True)
                            count = int(result.stdout.strip())
                            variant_counts[chrom] = count
                            total_variants_available += count
                        except Exception as e:
                            variant_counts[chrom] = 0
                            console.print(f"[yellow]⚠ Erro ao contar {chrom}: {e}[/yellow]")
                    else:
                        variant_counts[chrom] = 0
                    
                    progress.advance(count_task)
            
            # Exibir resumo das variantes
            console.print(f"\n[bold cyan]Variantes Disponíveis por Cromossomo:[/bold cyan]")
            console.print("─" * 70)
            
            for chrom in downloaded_chroms:
                count = variant_counts.get(chrom, 0)
                count_str = f"{count:,}".replace(',', '.')
                bar_width = min(40, int(count / max(variant_counts.values()) * 40)) if count > 0 else 0
                bar = "█" * bar_width
                console.print(f"  {chrom:>6}: {count_str:>12} variantes  {bar}")
            
            console.print("─" * 70)
            total_str = f"{total_variants_available:,}".replace(',', '.')
            console.print(f"  [bold]TOTAL: {total_str:>12} variantes[/bold]")
            
            # Calcular quantas serão selecionadas
            total_to_select = min(len(downloaded_chroms) * n_per_chrom, total_variants_available)
            select_str = f"{total_to_select:,}".replace(',', '.')
            percent = (total_to_select / total_variants_available * 100) if total_variants_available > 0 else 0
            console.print(f"  [green]Selecionando: {select_str} ({percent:.4f}%)[/green]\n")
        
        all_central_points: List[CentralPoint] = []
        rng = random.Random(random_seed) if selection_method == 'random' else None
        
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            task = progress.add_task(
                f"Selecionando variantes de {len(downloaded_chroms)} cromossomos...",
                total=len(downloaded_chroms)
            )
            
            for chrom in downloaded_chroms:
                # Busca específica com ponto após cromossomo para evitar chr1 pegar chr10-19
                vcf_files = list(vcf_dir.glob(f"*{chrom}.*.vcf.gz"))
                if not vcf_files:
                    console.print(f"[yellow]⚠ VCF não encontrado para {chrom}[/yellow]")
                    progress.advance(task)
                    continue
                
                vcf_path = vcf_files[0]
                
                # Extrair variantes do cromossomo de forma eficiente
                # Nota: VCFs filtrados do IGSR já estão filtrados (sem QUAL, FILTER=.)
                
                try:
                    # Usar pipeline bash eficiente para extrair N variantes
                    if selection_method == 'random':
                        # Usar shuf simples (sem random-source que pode não estar disponível)
                        cmd = f"bcftools view -H {vcf_path} | shuf -n {n_per_chrom}"
                        result = sp.run(cmd, capture_output=True, text=True, shell=True, check=True)
                        lines = [l for l in result.stdout.strip().split('\n') if l]
                    else:
                        # Sequential: pegar primeiras N - mais confiável
                        cmd = f"bcftools view -H {vcf_path} | head -n {n_per_chrom}"
                        result = sp.run(cmd, capture_output=True, text=True, shell=True, check=True)
                        lines = [l for l in result.stdout.strip().split('\n') if l]
                    
                    if not lines:
                        console.print(f"[yellow]⚠ Nenhuma variante selecionada em {chrom}[/yellow]")
                        progress.advance(task)
                        continue
                    
                    selected_lines = lines
                    
                    # Parse variantes
                    for line in selected_lines:
                        fields = line.split('\t')
                        if len(fields) < 8:  # Precisamos até INFO (coluna 8)
                            continue
                        
                        chr_name = fields[0]
                        pos = int(fields[1])
                        ref = fields[3]
                        alt = fields[4].split(',')[0]
                        # QUAL e FILTER são '.' nestes VCFs (já filtrados)
                        qual = 0.0
                        filter_status = '.'
                        
                        # Determinar tipo
                        if len(ref) == 1 and len(alt) == 1:
                            var_type = "SNV"
                        elif len(alt) > len(ref):
                            var_type = "INSERTION"
                        else:
                            var_type = "DELETION"
                        
                        variant = GenomicVariant(
                            chromosome=chr_name,
                            position=pos,
                            ref_allele=ref,
                            alt_allele=alt,
                            quality=qual,
                            depth=0,  # Não disponível no VCF multi-sample
                            allele_frequency=0.0,  # Pode ser extraído do INFO se necessário
                            filter_status=filter_status,
                            variant_type=var_type,
                            source_sample_id=f"multisample_{chrom}"
                        )
                        
                        central_point = CentralPoint(
                            variant=variant,
                            selected_for_dataset=True,
                            source_sample_id=f"multisample_{chrom}"
                        )
                        
                        all_central_points.append(central_point)
                    
                    console.print(f"[green]✓ {chrom}: {len(selected_lines)} variantes selecionadas[/green]")
                    
                except Exception as e:
                    console.print(f"[red]✗ Erro ao processar {chrom}: {e}[/red]")
                
                progress.advance(task)
        
        # Salvar pontos centrais
        self._save_central_points(all_central_points)
        
        console.print(f"\n[green]✓ Total: {len(all_central_points)} pontos centrais selecionados[/green]")
        console.print(f"[cyan]  Distribuição: ~{len(all_central_points) // len(downloaded_chroms) if downloaded_chroms else 0} por cromossomo[/cyan]\n")
        
        return all_central_points

    def _save_central_points(self, points: List[CentralPoint]):
        """Salva pontos centrais em arquivo."""
        output_file = self.output_dir / "central_points.json"
        data = [p.to_dict() for p in points]
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        console.print(f"[green]✓ Pontos centrais salvos em: {output_file}[/green]")
        self.state['central_points_selected'] = True
        self._save_checkpoint()

    def _load_central_points(self) -> List[CentralPoint]:
        """Carrega pontos centrais de arquivo."""
        input_file = self.output_dir / "central_points.json"
        with open(input_file) as f:
            data = json.load(f)

        points = []
        for item in data:
            points.append(CentralPoint.from_dict(item))

        return points
    
    # ───────────────────────────────────────────────────────────────
    # Passo 4: Extração de Sequências FASTA
    # ───────────────────────────────────────────────────────────────
    
    def extract_sequences(self, central_points: List[CentralPoint]) -> List[SequenceRecord]:
        """
        Extrai sequências FASTA centradas nos pontos selecionados.
        Despacha para método apropriado conforme modo de download.
        
        Args:
            central_points: Pontos centrais para extração
            
        Returns:
            Lista de registros de sequência
        """
        console.print("\n[bold cyan]Passo 4: Extração de Sequências FASTA[/bold cyan]")
        
        if self.config['debug']['dry_run']:
            console.print("[yellow]⚠ Modo dry-run: simulando extração de sequências[/yellow]")
            console.print(f"[cyan]Seria extraído: {len(central_points)} pontos centrais[/cyan]")
            return []
        
        # Detectar modo de download
        download_mode = self.config['data_sources'].get('download_mode', 'cram')
        
        if download_mode == 'vcf_multisample':
            return self._extract_sequences_multisample(central_points)
        else:
            return self._extract_sequences_individual(central_points)
    
    def _extract_sequences_individual(self, central_points: List[CentralPoint]) -> List[SequenceRecord]:
        """
        Extração de sequências no modo CRAM (amostras individuais).
        Cada amostra × variante gera uma sequência.
        """
        window_size = self.config['sequence_extraction']['window_size']
        if window_size % 2 != 0:
            raise ValueError("sequence_extraction.window_size deve ser um número par")

        ref_fasta = self._get_reference_fasta()
        if not ref_fasta:
            return []

        longevous_samples = self._load_sample_list("longevous_samples.txt", label=1)
        non_longevous_samples = self._load_sample_list("non_longevous_samples.txt", label=0)
        all_samples = longevous_samples + non_longevous_samples

        console.print(f"[cyan]Extraindo sequências para {len(all_samples)} amostras[/cyan]")
        console.print(f"[cyan]  {len(central_points)} pontos centrais × {len(all_samples)} amostras = {len(central_points) * len(all_samples)} sequências[/cyan]")

        use_alt = self.config['sequence_extraction'].get('use_alternate_allele', True)
        center_on_variant = self.config['sequence_extraction'].get('center_on_variant', True)
        sample_metadata = self.state.get('sample_metadata', {})

        records: List[SequenceRecord] = []
        sequences_dir = self.output_dir / "sequences"
        sequences_dir.mkdir(exist_ok=True)

        total_sequences = len(all_samples) * len(central_points)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task("Extraindo sequências...", total=total_sequences)

            for sample_id, label in all_samples:
                meta = sample_metadata.get(sample_id, {})
                vcf_path_str = meta.get('vcf_path')
                vcf_path = Path(vcf_path_str) if vcf_path_str else None
                variant_map = self._load_variants_for_sample(sample_id, vcf_path) if vcf_path else {}

                sample_dir = sequences_dir / sample_id
                sample_dir.mkdir(exist_ok=True)

                for point in central_points:
                    chrom = point.variant.chromosome
                    center_pos = point.variant.position
                    start = max(1, center_pos - window_size // 2) if center_on_variant else max(1, center_pos)
                    end = start + window_size - 1

                    variant_data = variant_map.get((chrom, center_pos))
                    allele_candidate, genotype, has_alt = self._resolve_sample_allele(variant_data, use_alt)

                    ref_allele = variant_data.get('ref', point.variant.ref_allele) if variant_data else point.variant.ref_allele

                    if has_alt and use_alt:
                        diff = len(allele_candidate) - len(ref_allele)
                        if diff > 0:
                            end += diff
                        elif diff < 0:
                            end += abs(diff)

                    region = f"{chrom}:{start}-{end}"
                    allele_used = allele_candidate if (use_alt and has_alt) else ref_allele

                    fasta_file = sample_dir / f"{sample_id}_{chrom}_{center_pos}_{ref_allele}>{allele_used}_w{window_size}.fasta"

                    try:
                        cmd = ['samtools', 'faidx', str(ref_fasta), region]
                        result = sp.run(cmd, capture_output=True, text=True, check=True)
                        sequence_ref = ''.join(result.stdout.split('\n')[1:]).upper()
                        center_idx = center_pos - start

                        # Capturar sequência de referência original (GRCh38)
                        grch38_seq = sequence_ref[:window_size]
                        if len(grch38_seq) < window_size:
                            grch38_seq += 'N' * (window_size - len(grch38_seq))
                        
                        if use_alt and has_alt:
                            sequence_final = self._apply_variant_to_sequence(sequence_ref, ref_allele, allele_candidate, center_idx, window_size)
                            allele_used = allele_candidate
                        else:
                            sequence_final = grch38_seq
                            allele_used = ref_allele

                        header = f">{sample_id}|{chrom}:{center_pos}|{ref_allele}>{allele_used}|win={window_size}"
                        with open(fasta_file, 'w') as fasta_handle:
                            fasta_handle.write(header + '\n')
                            for i in range(0, len(sequence_final), 80):
                                fasta_handle.write(sequence_final[i:i+80] + '\n')

                        record = SequenceRecord(
                            sample_id=sample_id,
                            central_point=point,
                            sequence=sequence_final,
                            GRCh38_s=grch38_seq,
                            fasta_file=fasta_file,
                            label=label,
                            genotype=genotype,
                            variant_present=has_alt
                        )
                        records.append(record)

                    except Exception as e:
                        console.print(f"[yellow]⚠ Erro ao extrair {region} para {sample_id}: {e}[/yellow]")

                    progress.advance(task)

        console.print(f"[green]✓ {len(records)} sequências extraídas[/green]")

        # Salvar índice de sequências
        self._save_sequence_index(records)
        self.state['sequences_extracted'] = [str(rec.fasta_file) for rec in records]
        self._save_checkpoint()

        return records
    
    def _extract_sequences_multisample(self, central_points: List[CentralPoint]) -> List[SequenceRecord]:
        """
        Extração de sequências no modo VCF multisample.
        Extrai uma sequência para cada variante dos pontos centrais.
        """
        console.print(f"[cyan]Modo: VCF Multisample[/cyan]")
        console.print(f"[cyan]Extraindo sequências para {len(central_points)} variantes[/cyan]")
        
        window_size = self.config['sequence_extraction']['window_size']
        if window_size % 2 != 0:
            raise ValueError("sequence_extraction.window_size deve ser um número par")
        
        ref_fasta = self._get_reference_fasta()
        if not ref_fasta:
            return []
        
        use_alt = self.config['sequence_extraction'].get('use_alternate_allele', True)
        center_on_variant = self.config['sequence_extraction'].get('center_on_variant', True)
        
        records: List[SequenceRecord] = []
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
            task = progress.add_task("Extraindo sequências...", total=len(central_points))
            
            for point in central_points:
                variant = point.variant
                chrom = variant.chromosome
                center_pos = variant.position
                ref_allele = variant.ref_allele
                alt_allele = variant.alt_allele
                
                # Calcular região
                start = max(1, center_pos - window_size // 2) if center_on_variant else max(1, center_pos)
                end = start + window_size - 1
                
                # Ajustar para indels
                if use_alt:
                    diff = len(alt_allele) - len(ref_allele)
                    if diff > 0:
                        end += diff
                    elif diff < 0:
                        end += abs(diff)
                
                region = f"{chrom}:{start}-{end}"
                allele_used = alt_allele if use_alt else ref_allele
                
                # Nome do arquivo (simplificado - apenas cromossomo e posição)
                seq_id = f"{chrom}_{center_pos}_w{window_size}"
                fasta_file = sequences_dir / f"{seq_id}.fasta"
                
                try:
                    # Extrair sequência de referência usando samtools
                    cmd = ['samtools', 'faidx', str(ref_fasta), region]
                    result = sp.run(cmd, capture_output=True, text=True, check=True)
                    sequence_ref = ''.join(result.stdout.split('\n')[1:]).upper()
                    
                    # Capturar sequência de referência original (GRCh38)
                    grch38_seq = sequence_ref[:window_size]
                    if len(grch38_seq) < window_size:
                        grch38_seq += 'N' * (window_size - len(grch38_seq))
                    
                    # Aplicar alelo alternativo se configurado
                    center_idx = center_pos - start
                    
                    if use_alt:
                        sequence_final = self._apply_variant_to_sequence(
                            sequence_ref, ref_allele, alt_allele, center_idx, window_size
                        )
                    else:
                        sequence_final = grch38_seq
                    
                    # Escrever FASTA
                    header = f">multisample|{chrom}:{center_pos}|win={window_size}"
                    with open(fasta_file, 'w') as fasta_handle:
                        fasta_handle.write(header + '\n')
                        for i in range(0, len(sequence_final), 80):
                            fasta_handle.write(sequence_final[i:i+80] + '\n')
                    
                    # Criar registro
                    record = SequenceRecord(
                        sample_id="multisample",
                        central_point=point,
                        sequence=sequence_final,
                        GRCh38_s=grch38_seq,
                        fasta_file=fasta_file,
                        label=1,  # Placeholder para modo multisample
                        genotype=None,
                        variant_present=True
                    )
                    records.append(record)
                    
                except Exception as e:
                    console.print(f"[yellow]⚠ Erro ao extrair {region}: {e}[/yellow]")
                
                progress.advance(task)
        
        console.print(f"[green]✓ {len(records)} sequências extraídas[/green]")
        
        # Salvar índice de sequências
        self._save_sequence_index(records)
        self.state['sequences_extracted'] = [str(rec.fasta_file) for rec in records]
        self._save_checkpoint()
        
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
        data = []
        for rec in records:
            base = rec.to_dict()
            base['selected'] = rec.central_point.selected_for_dataset
            data.append(base)
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
        self._save_checkpoint()

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

        analyzer = self._get_alphagenome_analyzer(api_key)
        if analyzer is None:
            return None

        outputs_requested = self.config['alphagenome']['outputs']
        ontology_terms = self.config['alphagenome'].get('ontology_terms')

        chrom = record.central_point.variant.chromosome
        center_pos = record.central_point.variant.position
        start_pos = max(0, center_pos - len(record.sequence) // 2)
        seq_id = f"{record.sample_id}_{chrom}_{center_pos}"

        try:
            ag_result = analyzer.predict_sequence(
                sequence=record.sequence,
                seq_id=seq_id,
                chromosome=chrom,
                start=int(start_pos),
                requested_outputs=outputs_requested,
                ontology_terms=ontology_terms
            )

            if not ag_result:
                return None

            outputs_obj = ag_result.get('outputs')
            requested = ag_result.get('requested_outputs', outputs_requested)

            predictions: Dict[str, Optional[np.ndarray]] = {}
            metadata: Dict[str, pd.DataFrame] = {}

            for output_name in requested:
                data_obj = getattr(outputs_obj, output_name.lower(), None) if outputs_obj else None
                array = self._convert_output_to_array(data_obj)
                predictions[output_name] = array

                meta_obj = getattr(data_obj, 'metadata', None) if data_obj is not None else None
                if meta_obj is not None:
                    if isinstance(meta_obj, pd.DataFrame):
                        metadata[output_name] = meta_obj
                    else:
                        try:
                            metadata[output_name] = pd.DataFrame(meta_obj)
                        except Exception:
                            continue

            result = AlphaGenomeResult(
                sequence_record=record,
                predictions=predictions,
                metadata=metadata
            )

            stats = result.aggregate_statistics()
            result.statistics = stats
            result.feature_vector = self._build_alphagenome_feature_vector(stats)

            if self.config['alphagenome']['cache_results']:
                with open(cache_file, 'wb') as f:
                    pickle.dump(result, f)

            processed_files = self.state.get('alphagenome_processed', [])
            fasta_path_str = str(record.fasta_file)
            if fasta_path_str not in processed_files:
                processed_files.append(fasta_path_str)
                self.state['alphagenome_processed'] = processed_files

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

        if not results:
            console.print("[yellow]⚠ Nenhum resultado AlphaGenome disponível para construção do dataset[/yellow]")
            return

        dataset_dir = self.output_dir / "torch_dataset"
        dataset_dir.mkdir(exist_ok=True)

        outputs = self.config['alphagenome']['outputs']
        stat_order = self.config['dataset']['features']['alphagenome_predictions']['statistics']
        expected_feature_len = len(outputs) * len(stat_order)
        normalize_position = self.config['dataset']['features']['position'].get('normalize', False)
        metadata_fields = self.config['dataset']['features']['metadata'].get('fields', [])

        sequences_encoded = []
        positions = []
        alpha_features = []
        labels = []
        metadata_entries = []

        for res in results:
            record = res.sequence_record
            seq_encoded = self._encode_sequence_one_hot(record.sequence)
            sequences_encoded.append(seq_encoded)

            pos_value = float(record.central_point.variant.position)
            if normalize_position:
                pos_value = self._normalize_position(record.central_point.variant.chromosome, record.central_point.variant.position)
            positions.append(pos_value)

            feature_vec = res.feature_vector
            if feature_vec is None or feature_vec.size == 0:
                feature_vec = self._build_alphagenome_feature_vector(res.aggregate_statistics())
            feature_vec = np.asarray(feature_vec, dtype=np.float32)
            if feature_vec.size != expected_feature_len:
                padded = np.full(expected_feature_len, np.nan, dtype=np.float32)
                length = min(len(padded), feature_vec.size)
                padded[:length] = feature_vec[:length]
                feature_vec = padded
            alpha_features.append(feature_vec)

            labels.append(int(record.label))

            meta_entry: Dict[str, Any] = {}
            for field_name in metadata_fields:
                if field_name == 'sample_id':
                    meta_entry[field_name] = record.sample_id
                elif field_name == 'chromosome':
                    meta_entry[field_name] = record.central_point.variant.chromosome
                elif field_name == 'position':
                    meta_entry[field_name] = record.central_point.variant.position
                elif field_name == 'ref_allele':
                    meta_entry[field_name] = record.central_point.variant.ref_allele
                elif field_name == 'alt_allele':
                    meta_entry[field_name] = record.central_point.variant.alt_allele
                elif field_name == 'variant_type':
                    meta_entry[field_name] = record.central_point.variant.variant_type
                elif field_name == 'genotype':
                    meta_entry[field_name] = record.genotype
                elif field_name == 'allele_used':
                    meta_entry[field_name] = record.allele_used
                else:
                    meta_entry[field_name] = getattr(record, field_name, None)
            metadata_entries.append(meta_entry)

        sequences_array = np.stack(sequences_encoded)
        positions_array = np.asarray(positions, dtype=np.float32)
        alpha_array = np.stack(alpha_features)
        labels_array = np.asarray(labels, dtype=np.int64)

        rng = random.Random(self.config['dataset']['random_seed'])
        splits_cfg = self.config['dataset']['splits']
        balance = self.config['dataset'].get('balance_classes', False)

        indices = list(range(len(results)))

        if balance:
            class_groups: Dict[int, List[int]] = defaultdict(list)
            for idx, label in enumerate(labels_array):
                class_groups[int(label)].append(idx)

            splits = {'train': [], 'validation': [], 'test': []}
            for label_value, group_indices in class_groups.items():
                rng.shuffle(group_indices)
                n = len(group_indices)
                train_count = int(n * splits_cfg['train'])
                val_count = int(n * splits_cfg['validation'])
                test_count = n - train_count - val_count

                splits['train'].extend(group_indices[:train_count])
                splits['validation'].extend(group_indices[train_count:train_count + val_count])
                splits['test'].extend(group_indices[train_count + val_count:train_count + val_count + test_count])
        else:
            rng.shuffle(indices)
            train_end = int(len(indices) * splits_cfg['train'])
            val_end = train_end + int(len(indices) * splits_cfg['validation'])
            splits = {
                'train': indices[:train_end],
                'validation': indices[train_end:val_end],
                'test': indices[val_end:]
            }

        def build_split(split_indices: List[int]) -> Dict[str, Any]:
            split_sequences = np.stack([sequences_array[i] for i in split_indices]) if split_indices else np.empty((0, sequences_array.shape[1], sequences_array.shape[2]), dtype=np.float32)
            split_positions = np.asarray([positions_array[i] for i in split_indices], dtype=np.float32)
            split_alpha = np.stack([alpha_array[i] for i in split_indices]) if split_indices else np.empty((0, alpha_array.shape[1]), dtype=np.float32)
            split_labels = np.asarray([labels_array[i] for i in split_indices], dtype=np.int64)
            split_metadata = [metadata_entries[i] for i in split_indices]

            return {
                'sequences': split_sequences,
                'positions': split_positions,
                'alphagenome_features': split_alpha,
                'labels': split_labels,
                'metadata': split_metadata
            }

        for split_name, split_indices in splits.items():
            split_data = build_split(split_indices)
            output_file = dataset_dir / f"{split_name}.pkl"
            with open(output_file, 'wb') as fh:
                pickle.dump(split_data, fh)
            console.print(f"[green]✓ Dataset salvo: {output_file} ({len(split_indices)} amostras)")

        # CSV resumo
        summary_rows = []
        for idx, res in enumerate(results):
            record = res.sequence_record
            summary_rows.append({
                'index': idx,
                'sample_id': record.sample_id,
                'label': record.label,
                'chromosome': record.central_point.variant.chromosome,
                'position': record.central_point.variant.position,
                'ref': record.central_point.variant.ref_allele,
                'alt': record.central_point.variant.alt_allele,
                'genotype': record.genotype,
                'variant_present': record.variant_present,
                'fasta_path': str(record.fasta_file)
            })

        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(dataset_dir / 'samples.csv', index=False)

        self.state['dataset_built'] = True
        self._save_checkpoint()
        console.print(f"[green]✓ Dataset PyTorch construído em {dataset_dir}[/green]")
    
    # ───────────────────────────────────────────────────────────────
    # Pipeline Completo
    # ───────────────────────────────────────────────────────────────
    
    def run_pipeline(self):
        """Executa pipeline completo."""
        console.print("\n[bold green]Iniciando Pipeline de Construção do Dataset[/bold green]\n")

        steps = self.config['pipeline']['steps']
        active_steps = [step for step in self._step_order if steps.get(step)]

        try:
            # Passo 1: Download
            if steps['download_samples']:
                self.download_samples()
            else:
                console.print("[yellow]⏭ Passo 1 ignorado (download_samples desativado)[/yellow]")

            # Passo 2: Seleção de pontos centrais
            if steps['select_central_points']:
                central_points = self.select_central_points()
            elif any(steps.get(step) for step in ('extract_sequences', 'run_alphagenome', 'build_dataset')):
                points_file = self.output_dir / "central_points.json"
                if points_file.exists():
                    central_points = self._load_central_points()
                    console.print("[yellow]⏭ Passo 2 ignorado (utilizando pontos centrais existentes)[/yellow]")
                else:
                    raise RuntimeError(
                        "Central points not available. Run the select_central_points step first."
                    )
            else:
                central_points = []
                console.print("[yellow]⏭ Passo 2 ignorado (não requerido pelas etapas selecionadas)[/yellow]")

            # Passo 3: Extração de sequências
            if steps['extract_sequences']:
                records = self.extract_sequences(central_points)
            elif any(steps.get(step) for step in ('run_alphagenome', 'build_dataset')):
                index_file = self.output_dir / "sequences_index.json"
                if index_file.exists():
                    with open(index_file) as f:
                        records_data = json.load(f)
                    records = [SequenceRecord.from_dict(item) for item in records_data]
                    console.print("[yellow]⏭ Passo 3 ignorado (reutilizando sequências existentes)[/yellow]")
                else:
                    raise RuntimeError(
                        "Sequences not available. Run the extract_sequences step first."
                    )
            else:
                records = []
                console.print("[yellow]⏭ Passo 3 ignorado (não requerido pelas etapas selecionadas)[/yellow]")

            # Passo 4: AlphaGenome
            if steps['run_alphagenome']:
                results = self.run_alphagenome(records)
            else:
                results = []
                console.print("[yellow]⏭ Passo 4 ignorado (run_alphagenome desativado)[/yellow]")

            # Passo 5: Build dataset
            if steps['build_dataset']:
                self.build_dataset(results)
            else:
                console.print("[yellow]⏭ Passo 5 ignorado (build_dataset desativado)[/yellow]")

            if active_steps:
                console.print("\n[bold green]✓ Pipeline concluído com sucesso![/bold green]")
            else:
                console.print("\n[yellow]⚠ Nenhuma etapa foi selecionada para execução.[/yellow]")

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

    # Sobrescrever etapas específicas, se fornecido
    builder.override_steps(args.steps)

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

