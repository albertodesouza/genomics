# -*- coding: utf-8 -*-
"""
reference_realign_mapper.py
============================

`ReferenceRealignMapper`: alignment_mapping="reference_realign".

Unlike `BcftoolsChainMapper` (alignment_mapping="bcftools_chain"), which builds a cohort-wide
"expanded axis" wide enough to fit every INDEL seen across the whole selected sample set (one
population VCF scan + chain-file parse per gene, shared across samples), this mapper realigns
each individual independently, straight onto true reference-genome coordinates: no cohort-wide
axis, no chain-file parsing, no population VCF scan. An inserted stretch collapses via `max()`
onto its single reference anchor slot; a deleted stretch is zero-filled. Output is always exactly
`window_center_size` long, in true reference coordinates.

This is a thin per-sample orchestration layer over `analysis.indel_track_realignment`
(`load_haplotype_variants`/`realign_consensus_values`), which does the actual position-mapping
math and already has its own unit tests.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from genomics.predictors.genotype_based.analysis.indel_track_realignment import (
    load_haplotype_variants,
    realign_consensus_values,
)

REFERENCE_REALIGN_MAPPER_VERSION = "reference_realign_mapper_v1"


class ReferenceRealignMapper:
    """Per-individual, no-cohort-axis realignment onto true reference coordinates.

    Parameters
    ----------
    dataset_dir : Path
        Root of the canonical dataset (has `references/windows/<gene>/window_metadata.json`).
    consensus_dataset_dir : Path
        Root containing `individuals/<sample>/windows/<gene>/<sample>.window.consensus_ready.vcf.gz`
        (same file `BcftoolsChainMapper` already consumes -- typically equal to `dataset_dir`).
    window_center_size : int
        Fixed output length after centering on the gene's reference midpoint.
    """

    def __init__(self, *, dataset_dir: Path, consensus_dataset_dir: Path, window_center_size: int):
        self.dataset_dir = Path(dataset_dir)
        self.consensus_dataset_dir = Path(consensus_dataset_dir)
        self.window_center_size = int(window_center_size)
        self._window_meta_cache: Dict[str, Dict] = {}
        self._variants_path_cache: Dict[Tuple[str, str], Optional[Path]] = {}
        self._variants_cache: Dict[Tuple[str, str, str], List[Tuple[int, str, str]]] = {}
        self.profile_stats = {"tracks_built": 0, "missing_vcf": 0}

    def _window_meta(self, gene: str) -> Dict:
        cached = self._window_meta_cache.get(gene)
        if cached is None:
            meta_path = self.dataset_dir / "references" / "windows" / gene / "window_metadata.json"
            with open(meta_path) as f:
                cached = json.load(f)
            self._window_meta_cache[gene] = cached
        return cached

    def _center_crop_bounds(self, full_ref_length: int) -> Tuple[int, int]:
        """Same centering formula as `DynamicIndelAligner._alignment_ref_span`, duplicated here
        (deliberately, it's ~6 lines) to keep this mapper fully decoupled from the cohort-axis
        aligner -- there is no expanded axis to slice through in this method.
        """
        size = min(self.window_center_size, full_ref_length)
        center_offset = full_ref_length // 2
        start = max(0, center_offset - size // 2)
        end = min(full_ref_length, start + size)
        if end - start < size:
            start = max(0, end - size)
        return start, end

    def _consensus_vcf_path(self, gene: str, sample_id: str) -> Optional[Path]:
        key = (gene, sample_id)
        if key not in self._variants_path_cache:
            path = (
                self.consensus_dataset_dir / "individuals" / sample_id / "windows" / gene
                / f"{sample_id}.window.consensus_ready.vcf.gz"
            )
            self._variants_path_cache[key] = path if path.exists() else None
        return self._variants_path_cache[key]

    def _variants_for(self, gene: str, sample_id: str, haplotype: str) -> Optional[List[Tuple[int, str, str]]]:
        key = (gene, sample_id, haplotype)
        cached = self._variants_cache.get(key)
        if cached is not None:
            return cached
        vcf_path = self._consensus_vcf_path(gene, sample_id)
        if vcf_path is None:
            return None
        variants = load_haplotype_variants(vcf_path, haplotype)
        self._variants_cache[key] = variants
        return variants

    def get_haplotype_track(
        self, gene: str, sample_id: str, haplotype: str, row: np.ndarray,
    ) -> Optional[np.ndarray]:
        """Realign one raw AlphaGenome track onto true reference coordinates, then center-crop.

        `row` is the raw per-position prediction for this (sample, gene, haplotype, track) --
        already loaded by the caller from the base dataset's `predictions` dict (full un-cropped
        window length, i.e. the same length as `ref.window.fa`). No disk I/O happens here beyond
        the small per-sample, per-window consensus VCF read.

        Returns a `window_center_size`-length float32 array, or None if the sample's consensus
        VCF for this gene isn't present on disk.
        """
        variants = self._variants_for(gene, sample_id, haplotype)
        if variants is None:
            self.profile_stats["missing_vcf"] += 1
            return None

        window_meta = self._window_meta(gene)
        start_1based = int(window_meta["start"])
        end_1based = int(window_meta["end"])
        full_ref_length = end_1based - start_1based + 1
        if full_ref_length != len(row):
            raise ValueError(
                f"{gene}: full_ref_length={full_ref_length} não corresponde a len(row)={len(row)}"
            )

        realigned_full = realign_consensus_values(
            consensus_values=np.asarray(row, dtype=np.float32).reshape(-1, 1),
            window_start=start_1based,
            window_width=full_ref_length,
            variants=variants,
        )[:, 0]

        crop_start, crop_end = self._center_crop_bounds(full_ref_length)
        self.profile_stats["tracks_built"] += 1
        return realigned_full[crop_start:crop_end].astype(np.float32, copy=False)
