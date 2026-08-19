"""Shared state for the promoter-location / knockout-experiment sections (CAGE lookups,
sequence-scramble knockouts, and their plots) -- the fixed bundle of config/model/annotation
state established once and threaded through every function in `cage.py`, `knockout.py`, and the
knockout/CAGE plots in `plotting.py`.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


@dataclass
class KnockoutContext:
    ag_client: Any
    organism: Any
    config: Any                 # KO_CONFIG: the dita_no_masks predictor config
    genes: Sequence[str]
    ontology_terms: Sequence[str]
    model: Any                  # KO_MODEL
    device: Any
    dataset_dir: Path
    full_ds: Any                 # ko_full_ds: real ProcessedGenomicDataset (not the cache reader)
    normalization_params: Any
    strong_idx: int
    weak_idx: int
    class_names: list
    class_colors: dict
    target_idx: int
    other_idx: int
    cache_dir: Path              # AlphaGenome prediction cache
    fig_dir: Path
    tss_df_mane: Any
    tss_df_coding: Any
    gtf_gene_rows: Any
    gene_id_map: Any
    transcript_extractor_mane: Any
    transcript_extractor_coding: Any
