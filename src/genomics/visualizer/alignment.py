"""Training-axis (bcftools_chain expanded alignment) access for the visualizer.

Wraps ``DynamicIndelAligner``/``BcftoolsChainMapper`` with the same sample set and centre window
the CNN training and the legacy track viewer use, so their persistent caches are shared. Entries
are converted once to compact numpy index arrays instead of re-walking the JSON index lists on
every request.
"""
from __future__ import annotations

import shutil
import threading
from typing import Any, Dict, Optional, Tuple

import numpy as np

from genomics.visualizer.cache import LRUCache
from genomics.visualizer.datasets import Dataset


class AlignmentService:
    def __init__(self, model_window: int = 32768, cache_bytes: int = 512 << 20, consensus_dirs: Optional[Dict[str, str]] = None):
        self.model_window = int(model_window)
        self.entries = LRUCache(cache_bytes)
        self.consensus_dirs = consensus_dirs or {}
        self._contexts: Dict[str, Tuple[Any, Any]] = {}
        self._axes: Dict[Tuple[str, str], Dict[str, Any]] = {}
        self._built: set = set()
        self._lock = threading.Lock()
        self._gene_locks: Dict[Tuple[str, str], threading.Lock] = {}

    @staticmethod
    def available() -> Dict[str, Any]:
        return {"bcftools": shutil.which("bcftools") is not None}

    def _context(self, dataset: Dataset):
        with self._lock:
            ctx = self._contexts.get(dataset.id)
            if ctx is None:
                from pathlib import Path

                from genomics.predictors.genotype_based.alignment.bcftools_chain_mapper import BcftoolsChainMapper
                from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner

                sample_ids = [row["sample_id"] for row in dataset.samples]
                aligner = DynamicIndelAligner(dataset.path, selected_sample_ids=sample_ids, center_window_size=self.model_window)
                consensus = Path(self.consensus_dirs.get(dataset.id) or dataset.path)
                mapper = BcftoolsChainMapper(dataset_dir=dataset.path, consensus_dataset_dir=consensus, aligner=aligner)
                ctx = (aligner, mapper)
                self._contexts[dataset.id] = ctx
            return ctx

    def _gene_lock(self, dataset: Dataset, gene: str) -> threading.Lock:
        with self._lock:
            return self._gene_locks.setdefault((dataset.id, gene), threading.Lock())

    def ensure_axis(self, dataset: Dataset, gene: str) -> None:
        key = (dataset.id, gene)
        if key in self._built:
            return
        with self._gene_lock(dataset, gene):
            if key in self._built:
                return
            aligner, _mapper = self._context(dataset)
            aligner.build_alignment_axis_for_gene(gene, [row["sample_id"] for row in dataset.samples])
            self._built.add(key)

    def has_axis(self, dataset: Dataset, gene: str) -> bool:
        return (dataset.id, gene) in self._axes

    def axis(self, dataset: Dataset, gene: str) -> Dict[str, Any]:
        key = (dataset.id, gene)
        cached = self._axes.get(key)
        if cached is not None:
            return cached
        self.ensure_axis(dataset, gene)
        aligner, _mapper = self._context(dataset)
        axis = aligner.get_alignment_axis(gene)
        expanded_length = int(axis["expanded_length"])
        ref_of_expanded = np.full(expanded_length, -1, dtype=np.int64)
        for ref_idx, exp_idx in axis["expanded_index_map"].items():
            if 0 <= int(exp_idx) < expanded_length:
                ref_of_expanded[int(exp_idx)] = int(ref_idx)
        center = aligner.get_reference_centered_expanded_slice(gene, self.model_window)
        result = {
            "expanded_length": expanded_length,
            "ref_length": int(axis["ref_length"]),
            "ref_start_offset": int(axis["ref_start_offset"]),
            "insertion_slots": np.nonzero(ref_of_expanded < 0)[0].astype(np.int64),
            "ref_of_expanded": ref_of_expanded,
            "model_window": {"start": int(center["expanded_start"]), "end": int(center["expanded_end"])},
            "sample_set_key": axis.get("sample_set_key"),
        }
        self._axes[key] = result
        return result

    def axis_key(self, dataset: Dataset, gene: str) -> str:
        axis = self.axis(dataset, gene)
        return f"{axis['sample_set_key']}:{axis['expanded_length']}:{axis['ref_start_offset']}"

    def entry_arrays(self, dataset: Dataset, gene: str, sample: str, haplotype: str) -> Tuple[np.ndarray, np.ndarray]:
        """``(expanded_indices, source_indices)`` for one haplotype on the training axis."""
        if haplotype not in ("H1", "H2"):
            raise ValueError("Aligned coordinates need haplotype H1 or H2")

        def load() -> Tuple[np.ndarray, np.ndarray]:
            # Cached entries load without bcftools; only cache misses shell out to it.
            self.ensure_axis(dataset, gene)
            _aligner, mapper = self._context(dataset)
            entry = mapper.get_haplotype_entry(gene, sample, haplotype)
            if entry is None:
                raise FileNotFoundError(f"No alignment entry for {sample} {haplotype} ({gene})")
            expanded = np.asarray(entry.get("expanded_indices", []), dtype=np.int64)
            source = np.asarray(entry.get("copy_from_indices", []), dtype=np.int64)
            if str(entry.get("mapping_method")) != "bcftools_chain":
                source = source + int(entry.get("source_start_idx", 0))
            return expanded, source

        return self.entries.get_or_load(("entry", dataset.id, gene, sample, haplotype), load)

    def reference_offset_to_expanded(self, dataset: Dataset, gene: str) -> np.ndarray:
        """Map full-window reference offsets to expanded indices (-1 outside the axis)."""
        axis = self.axis(dataset, gene)
        key = ("ref2exp", dataset.id, gene)
        hit = self.entries.get(key)
        if hit is not None:
            return hit
        ref_of_expanded = axis["ref_of_expanded"]
        size = axis["ref_start_offset"] + axis["ref_length"]
        out = np.full(max(size, 1), -1, dtype=np.int64)
        valid = ref_of_expanded >= 0
        out[ref_of_expanded[valid] + axis["ref_start_offset"]] = np.nonzero(valid)[0]
        self.entries.put(key, out)
        return out
