"""Reports derived from VCF, coverage, and gene annotations."""

from .legacy import (
    combine_pairwise_stats,
    gene_list_from_gtf,
    gene_presence_reports,
    pairwise_comparisons,
    paternity_analysis,
    trio_denovo_report,
)

__all__ = [
    "combine_pairwise_stats",
    "gene_list_from_gtf",
    "gene_presence_reports",
    "pairwise_comparisons",
    "paternity_analysis",
    "trio_denovo_report",
]
