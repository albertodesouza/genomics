"""Shared experiment/client/split bootstrapping for the `0N_*.ipynb` notebooks
(the split of what used to be one `simplified_notebook.ipynb`).

Kept as explicit, composable calls -- `load_experiment()` (no API key needed) and
`create_client()` (needs `ALPHAGENOME_API_KEY`) are separate so notebooks that only read
from `notebooks/.cache/predictions/` (no fresh AlphaGenome calls) don't need to load an
API key at all.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from . import annotations as annotations_utils
from . import data as data_utils
from .family_split import family_aware_train_test_split


@dataclass
class ExperimentContext:
    notebook_dir: Path
    GENES: List[str]
    ONTOLOGY_TERMS: List[str]
    CLASS_MAP: Dict[str, List[str]]
    WINDOW_CENTER_SIZE: int
    gene_rows: Dict[str, dict]
    individuals: Dict[str, dict]
    samples_by_class: Dict[str, List[str]]
    windows_df: pd.DataFrame
    gtf: pd.DataFrame
    ANNOTATIONS_CACHE_DIR: Path
    REFERENCE_PREDICTIONS_CACHE_DIR: Path
    CONSENSUS_CACHE_DIR: Path
    PREDICTIONS_CACHE_DIR: Path


def _build_windows_df(gene_rows: Dict[str, dict], expected_window_size: int) -> pd.DataFrame:
    windows_df = pd.DataFrame([
        {
            "gene": gene,
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "window_size": row["end"] - row["start"],
        }
        for gene, row in gene_rows.items()
    ]).sort_values("gene").reset_index(drop=True)
    assert (windows_df["window_size"] == expected_window_size).all()
    return windows_df


def load_experiment(notebook_dir: Path) -> ExperimentContext:
    """Loads everything downstream notebook cells expect as top-level globals: gene/ontology/
    class config, individuals, `windows_df`, `gtf`, and the notebook's `.cache/` subdirectories.
    Makes no AlphaGenome API calls -- see `create_client()` for that."""
    from alphagenome.models import dna_client

    notebook_dir = Path(notebook_dir)
    annotations_cache_dir = notebook_dir / ".cache" / "annotations"
    reference_predictions_cache_dir = notebook_dir / ".cache" / "reference_predictions"
    consensus_cache_dir = notebook_dir / ".cache" / "consensus"
    predictions_cache_dir = notebook_dir / ".cache" / "predictions"

    experiment = data_utils.load_experiment()
    gene_rows = data_utils.load_genes()
    individuals = data_utils.load_individuals()
    samples_by_class = data_utils.samples_by_class(individuals, experiment["class_map"])

    windows_df = _build_windows_df(gene_rows, dna_client.SEQUENCE_LENGTH_500KB)
    gtf = annotations_utils.load_gtf(annotations_cache_dir)

    return ExperimentContext(
        notebook_dir=notebook_dir,
        GENES=experiment["genes"],
        ONTOLOGY_TERMS=experiment["ontology_terms"],
        CLASS_MAP=experiment["class_map"],
        WINDOW_CENTER_SIZE=experiment["window_center_size"],
        gene_rows=gene_rows,
        individuals=individuals,
        samples_by_class=samples_by_class,
        windows_df=windows_df,
        gtf=gtf,
        ANNOTATIONS_CACHE_DIR=annotations_cache_dir,
        REFERENCE_PREDICTIONS_CACHE_DIR=reference_predictions_cache_dir,
        CONSENSUS_CACHE_DIR=consensus_cache_dir,
        PREDICTIONS_CACHE_DIR=predictions_cache_dir,
    )


def create_client():
    """Picks up ALPHAGENOME_API_KEY from either the shell environment or ~/.env (this repo's
    convention for storing API keys outside of chat/context -- see the `credentials` skill) and
    returns (client, ORGANISM)."""
    from alphagenome.models import dna_client
    from dotenv import load_dotenv

    load_dotenv(Path.home() / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise RuntimeError(
            "ALPHAGENOME_API_KEY not found in the environment or ~/.env. "
            "Get a key at https://www.alphagenomedocs.com/ and add it with:\n"
            '  printf "Enter ALPHAGENOME_API_KEY (typing hidden): " && read -s val && echo '
            '&& echo "ALPHAGENOME_API_KEY=$val" >> ~/.env'
        )
    client = dna_client.create(api_key=api_key)
    return client, dna_client.Organism.HOMO_SAPIENS


def get_train_test_split(
    individuals: Dict[str, dict], class_map: Dict[str, List[str]]
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Fixed-config wrapper around `family_aware_train_test_split` (test_size=0.2,
    random_seed=13) so every notebook in this split gets the identical train/test split."""
    return family_aware_train_test_split(individuals, class_map, test_size=0.2, random_seed=13)
