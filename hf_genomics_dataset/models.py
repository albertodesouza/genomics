from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass(frozen=True)
class SourceDataset:
    name: str
    path: Path


@dataclass(frozen=True)
class ConversionOptions:
    chunk_size: int = 0
    include_sequences: bool = True


@dataclass
class SampleRecord:
    sample_id: str
    family_id: Optional[str]
    sex: Optional[int]
    population: Optional[str]
    superpopulation: Optional[str]
    longevity: Optional[bool]
    frog_likelihood: List[float]
    frog_population_names: List[str]
    available_genes: List[str]
    created_at: Optional[str]
    last_updated: Optional[str]
    source_dataset: str
    source_path: str
    source_sample_path: str


@dataclass
class GeneWindowRecord:
    sample_id: str
    gene: str
    source_dataset: str
    chromosome: Optional[str]
    start: Optional[int]
    end: Optional[int]
    window_size: Optional[int]
    window_type: Optional[str]
    ref_sequence: Optional[str] = None
    h1_sequence: Optional[str] = None
    h2_sequence: Optional[str] = None
    outputs: List[str] = field(default_factory=list)
    ontologies: List[str] = field(default_factory=list)
    predictions_h1: Dict[str, List[List[float]]] = field(default_factory=dict)
    predictions_h2: Dict[str, List[List[float]]] = field(default_factory=dict)
    source_path: str = ""
    source_sample_path: str = ""
    source_window_path: str = ""
    source_metadata_path: str = ""


@dataclass
class BuildSummary:
    source_datasets: List[str]
    total_samples: int
    total_gene_windows: int
    genes: List[str]
    duplicate_samples_skipped: int = 0
    duplicate_gene_windows_skipped: int = 0
    warnings: List[str] = field(default_factory=list)
