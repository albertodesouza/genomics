"""Interactive, dataset-agnostic genomics visualizer.

A single in-process web server (``genomics visualize``) with a static single-page frontend for
browsing canonical-layout datasets: cohort/sample metadata, AlphaGenome prediction tracks in
genomic, haplotype or training-alignment coordinates, haplotype sequences against the reference,
and experiment runs. Heavy dependencies (pandas for GTF annotations, bcftools for the training
alignment axis) are optional and loaded lazily.
"""
