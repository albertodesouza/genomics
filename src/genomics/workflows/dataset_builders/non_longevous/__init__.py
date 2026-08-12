"""
build_non_longevous_dataset

Pipeline to build datasets of non-longevous individuals from 1000 Genomes Project.
Now with full PyTorch Dataset support!

Modules:
    - frog_ancestry_parser: FROGAncestryCalc likelihood parser
    - dataset_builder: Per-individual and global metadata builder
    - genomic_dataset: PyTorch Dataset for genomic data
    - build_non_longevous_dataset: Main pipeline
    - build_window_and_predict: Window extraction and AlphaGenome predictions

Example:
    >>> from build_non_longevous_dataset.genomic_dataset import GenomicLongevityDataset
    >>> dataset = GenomicLongevityDataset('non_longevous_results')
    >>> print(f"Total: {len(dataset)} individuals")

Author: Alberto F. De Souza
Version: 2.0.0 (with PyTorch Dataset)
Last updated: 2025-11-11
"""

__version__ = "2.0.0"
__author__ = "Alberto F. De Souza"

# Importar classes principais para facilitar uso
try:
    from .genomic_dataset import (
        GenomicLongevityDataset,
        collate_fn_variable_windows
    )
    from .frog_ancestry_parser import FROGAncestryData
    from .dataset_builder import (
        IndividualDatasetBuilder,
        DatasetMetadataBuilder
    )
    
    __all__ = [
        'GenomicLongevityDataset',
        'collate_fn_variable_windows',
        'FROGAncestryData',
        'IndividualDatasetBuilder',
        'DatasetMetadataBuilder',
    ]
    
except ImportError as e:
    # Imports relativos podem falhar se executado como script
    # Neste caso, as classes ainda podem ser importadas diretamente
    __all__ = []

