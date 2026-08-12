from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    gene_id: str

    @property
    def bcftools_region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass(frozen=True)
class SampleMetadata:
    sample_id: str
    target: str
    family_id: Optional[str] = None
    population: Optional[str] = None
    superpopulation: Optional[str] = None
    extra: Optional[Dict] = None


@dataclass(frozen=True)
class VariantToken:
    chrom: str
    position: int
    position_relative: int
    gene_id: str
    haplotype: str
    reference_allele: str
    alternate_allele: str
    variant_type: str
    length: int


def chromosome_sort_key(chrom: str) -> tuple[int, str]:
    text = str(chrom)
    token = text[3:] if text.startswith("chr") else text
    if token.isdigit():
        return int(token), ""
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return special.get(token.upper(), 10_000), token


def sort_tokens(tokens: List[VariantToken]) -> List[VariantToken]:
    return sorted(tokens, key=lambda t: (chromosome_sort_key(t.chrom), int(t.position), t.gene_id, t.haplotype))
