"""Helpers for manifest-based INDEL-aware alignment."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional


class IndelManifest:
    """Read-only helper around a JSON manifest with expanded-axis mappings."""

    def __init__(self, manifest_path: str | Path):
        self.path = Path(manifest_path)
        if not self.path.exists():
            raise FileNotFoundError(f"Manifest de INDELs não encontrado: {self.path}")

        with open(self.path) as f:
            self.data = json.load(f)

        self.genes = self.data.get("genes", {})

    def get_gene_manifest(self, gene_name: str) -> Dict[str, Any]:
        gene_manifest = self.genes.get(gene_name)
        if gene_manifest is None:
            raise KeyError(f"Gene ausente no manifest de INDELs: {gene_name}")
        return gene_manifest

    def get_haplotype_entry(self, gene_name: str, sample_id: str, haplotype: str) -> Optional[Dict[str, Any]]:
        gene_manifest = self.get_gene_manifest(gene_name)
        sample_entry = gene_manifest.get("samples", {}).get(sample_id)
        if sample_entry is None:
            return None
        return sample_entry.get(haplotype)

    def get_expanded_length(self, gene_name: str) -> int:
        gene_manifest = self.get_gene_manifest(gene_name)
        expanded_length = gene_manifest.get("expanded_length")
        if expanded_length is None:
            raise KeyError(f"expanded_length ausente para gene {gene_name}")
        return int(expanded_length)
