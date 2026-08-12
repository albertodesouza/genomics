from __future__ import annotations

import torch

from .constants import BASE_PAD_ID, BASE_TO_ID


def encode_allele(allele: str, l_max: int) -> torch.Tensor:
    ids = [BASE_TO_ID.get(base.upper(), BASE_TO_ID["N"]) for base in str(allele or "N")[:l_max]]
    if len(ids) < l_max:
        ids.extend([BASE_PAD_ID] * (l_max - len(ids)))
    return torch.tensor(ids, dtype=torch.long)


def classify_variant(ref: str, alt: str) -> str:
    ref = str(ref).upper()
    alt = str(alt).upper()
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(alt) > len(ref):
        return "INS"
    if len(alt) < len(ref):
        return "DEL"
    return "SNP"


def normalize_length(length: int, max_indel_size: int) -> float:
    denom = max(float(max_indel_size), 1.0)
    value = float(length) / denom
    return max(min(value, 1.0), -1.0)
