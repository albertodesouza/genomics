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

   ``python neural_longevity_dataset.py --config longevity_config.yaml``

   O programa respeita as etapas definidas em ``pipeline.steps`` dentro do YAML.

Selecionando etapas específicas
-------------------------------
Você pode rodar etapas isoladas informando ``--steps`` (ex.: apenas download e
extração de sequências):

``python neural_longevity_dataset.py --config longevity_config.yaml \
    --steps download_samples extract_sequences``

O argumento ``--dry-run`` mostra o que seria executado sem alterar arquivos.

Estrutura de saída
------------------
- ``cram/``: CRAM/CRAI baixados do ENA
- ``vcf/``: VCFs gerados via *bcftools*
- ``windows/``: FASTAs centrados nas variantes, com metadados por amostra
- ``alphagenome/``: caches das predições e estatísticas agregadas
- ``dataset/``: arquivos ``train.pkl``, ``validation.pkl`` e ``test.pkl``
  acompanhados de ``samples.csv`` para inspeção manual

Autor: IA Neuro-Simbólica para Longevidade
Data: Outubro 2025
"""

import argparse
import sys
import yaml
import json
import pickle
import random
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
import subprocess as sp
import hashlib
from collections import defaultdict
import gzip
from urllib.parse import quote

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


# Diretório raiz obrigatório para dados persistentes do projeto
DATA_ROOT = Path("/dados/GENOMICS_DATA/top3").resolve()


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
    genotype: str = "0/0"
    variant_present: bool = False
    allele_used: str = "REF"
    vcf_path: Optional[Path] = None
    quality: Optional[float] = None
    depth: Optional[int] = None
    allele_frequency: Optional[float] = None

    def to_dict(self) -> Dict:
        return {
            'sample_id': self.sample_id,
            'variant': self.central_point.variant.to_dict(),
            'sequence': self.sequence,
            'sequence_length': len(self.sequence),
            'label': self.label,
            'genotype': self.genotype,
            'variant_present': self.variant_present,
            'allele_used': self.allele_used,
            'quality': self.quality,
            'depth': self.depth,
            'allele_frequency': self.allele_frequency
        }


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
        
        # Estado do processamento
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.state = self._load_checkpoint()
        self._variant_cache_runtime: Dict[str, Dict[Tuple[str, int], Dict[str, Any]]] = {}
        self._alphagenome_client = None
        self._chrom_lengths: Optional[Dict[str, int]] = None
        
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
            'variant_cache': {}
        }

        for key, value in defaults.items():
            if key not in state:
                state[key] = value

        return state

    # ───────────────────────────────────────────────────────────────
    # Utilitários internos
    # ───────────────────────────────────────────────────────────────

    def _normalize_ena_url(self, url: str) -> str:
        """Converte URLs FTP do ENA para HTTPS."""
        if not url:
            return url
        if url.startswith('ftp://'):
            return url.replace('ftp://', 'https://', 1)
        if url.startswith('http://'):
            return url.replace('http://', 'https://', 1)
        return url

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

        while total < max_records:
            limit = min(batch_size, max_records - total)
            query = quote(f'study_accession="{project_id}" AND submitted_format="cram"')
            url = (
                "https://www.ebi.ac.uk/ena/portal/api/search?"
                f"result=read_run&format=tsv&fields={fields}&limit={limit}&offset={offset}&query={query}"
            )

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
        if 'sample_alias' in combined.columns:
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

        for _, row in subset.iterrows():
            sample_alias = str(row.get('sample_alias') or '').strip()
            if not sample_alias:
                continue

            submitted = str(row.get('submitted_ftp') or '').strip()
            if not submitted:
                continue

            urls = [self._normalize_ena_url(part.strip()) for part in submitted.split(';') if part.strip()]
            cram_url = next((u for u in urls if u.endswith('.cram')), None)
            crai_url = next((u for u in urls if u.endswith('.crai')), None)

            if not cram_url or not crai_url:
                continue

            records.append({
                'sample_id': sample_alias,
                'run_accession': str(row.get('run_accession', '')).strip(),
                'label': label,
                'cram_url': cram_url,
                'crai_url': crai_url
            })

        return records

    def _download_file(self, url: str, destination: Path) -> bool:
        """Baixa arquivo com wget (idempotente)."""
        destination.parent.mkdir(parents=True, exist_ok=True)

        if destination.exists() and destination.stat().st_size > 0:
            return True

        try:
            cmd = ['wget', '-q', '-O', str(destination), url]
            sp.run(cmd, check=True)
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

    def _get_reference_fasta(self) -> Optional[Path]:
        """Retorna caminho absoluto do FASTA de referência."""
        ref_path = Path(self.config['data_sources']['reference']['fasta'])
        if not ref_path.is_absolute():
            ref_fasta = self.work_dir / ref_path
        else:
            ref_fasta = ref_path

        if not ref_fasta.exists():
            console.print(f"[red]✗ Referência não encontrada: {ref_fasta}[/red]")
            return None

        self._ensure_reference_indices(ref_fasta)
        return ref_fasta

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

            with sp.Popen(mpileup_cmd, stdout=sp.PIPE) as mpileup_proc:
                with sp.Popen(call_cmd, stdin=mpileup_proc.stdout) as call_proc:
                    mpileup_proc.stdout.close()
                    call_proc.communicate()

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
                from neural_module import AlphaGenomeAnalyzer, DEFAULT_CONFIG
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
        Baixa CRAM/CRAI do ENA e gera VCFs correspondentes.
        """
        console.print("\n[bold cyan]Passo 1: Download de Amostras[/bold cyan]")

        is_dry_run = self.config['debug']['dry_run']

        if is_dry_run:
            console.print("[yellow]⚠ Modo dry-run: simulando downloads[/yellow]")

        project_id = self.config['data_sources']['longevous'].get('ena_project', 'PRJEB31736')
        longevous_range = self.config['data_sources']['longevous']['sample_range']
        non_longevous_range = self.config['data_sources']['non_longevous']['sample_range']

        max_records = max(longevous_range[1], non_longevous_range[1])
        console.print(f"[cyan]Consultando ENA ({project_id}) para {max_records} registros...[/cyan]")
        runs_df = self._fetch_project_runs(project_id, max_records)

        if runs_df.empty:
            console.print("[red]✗ Não foi possível obter lista de amostras reais.[/red]")
            self._create_placeholder_samples("longevous_samples.txt", longevous_range[1] - longevous_range[0])
            self._create_placeholder_samples("non_longevous_samples.txt", non_longevous_range[1] - non_longevous_range[0])
            return

        longevous_records = self._prepare_sample_records(runs_df, longevous_range, label=1)
        non_longevous_records = self._prepare_sample_records(runs_df, non_longevous_range, label=0)

        console.print(f"[green]✓ {len(longevous_records)} registros longevos selecionados[/green]")
        console.print(f"[green]✓ {len(non_longevous_records)} registros não-longevos selecionados[/green]")

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
        self.state['central_points_selected'] = True
        self._save_checkpoint()
    
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

                        if use_alt and has_alt:
                            sequence_final = self._apply_variant_to_sequence(sequence_ref, ref_allele, allele_candidate, center_idx, window_size)
                            allele_used = allele_candidate
                        else:
                            sequence_final = sequence_ref[:window_size]
                            if len(sequence_final) < window_size:
                                sequence_final += 'N' * (window_size - len(sequence_final))
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
                            fasta_file=fasta_file,
                            label=label,
                            genotype=genotype,
                            variant_present=has_alt,
                            allele_used=allele_used,
                            vcf_path=vcf_path,
                            quality=variant_data.get('qual') if variant_data else None,
                            depth=variant_data.get('depth') if variant_data else None,
                            allele_frequency=variant_data.get('allele_frequency') if variant_data else None
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

