from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import torch
from torch.utils.data import Dataset


class VariantTokenDataset(Dataset):
    def __init__(
        self,
        processed_dir: str | Path,
        split: str,
        max_sequence_length: Optional[int] = None,
        truncate_policy: str = "error",
        loading_strategy: str = "lazy",
    ):
        self.processed_dir = Path(processed_dir)
        self.split = split
        self.max_sequence_length = max_sequence_length
        self.truncate_policy = truncate_policy
        self.loading_strategy = loading_strategy
        with open(self.processed_dir / "metadata.json") as f:
            self.metadata = json.load(f)
        with open(self.processed_dir / "sample_index.json") as f:
            index_payload = json.load(f)
        with open(self.processed_dir / "splits.json") as f:
            splits = json.load(f)
        split_ids = set(splits.get(split, []))
        self.samples = [row for row in index_payload.get("samples", []) if row.get("sample_id") in split_ids]
        if not self.samples:
            self.samples = []
        self.gene_vocab = self._load_gene_vocab()
        self._data = None
        if self.loading_strategy == "preload":
            self._data = [self._load_item(row) for row in self.samples]

    def _load_gene_vocab(self) -> Dict[str, int]:
        with open(self.processed_dir / "gene_vocab.json") as f:
            return json.load(f)

    def __len__(self) -> int:
        return len(self.samples)

    def get_num_classes(self) -> int:
        return len(self.metadata.get("classes", []))

    def get_class_names(self) -> List[str]:
        return list(self.metadata.get("classes", []))

    @property
    def idx_to_target(self) -> Dict[int, str]:
        return {idx: name for idx, name in enumerate(self.get_class_names())}

    def _truncate(self, item: Dict) -> Dict:
        max_len = self.max_sequence_length
        if max_len is None:
            return item
        n = int(item.get("num_tokens", item["variant_type"].shape[0]))
        if n <= max_len:
            return item
        if self.truncate_policy == "error":
            raise ValueError(
                f"Amostra {item.get('sample_id')} tem {n} tokens > max_sequence_length={max_len}; "
                "use truncate_policy keep_first/keep_last ou reduza as regioes"
            )
        if self.truncate_policy == "keep_last":
            sl = slice(n - max_len, n)
        else:
            sl = slice(0, max_len)
        out = dict(item)
        for key in ("variant_type", "haplotype", "gene", "length_norm", "position_relative", "position", "ref_allele", "alt_allele"):
            out[key] = item[key][sl]
        out["num_tokens"] = max_len
        return out

    def _load_item(self, row: Dict) -> Dict:
        path = self.processed_dir / row["path"]
        item = torch.load(path, map_location="cpu", weights_only=True)
        return self._truncate(item)

    def __getitem__(self, idx: int) -> Dict:
        if self._data is not None:
            return self._data[idx]
        return self._load_item(self.samples[idx])


def _pad_1d(values: List[torch.Tensor], pad_value: int = 0) -> torch.Tensor:
    max_len = max((v.shape[0] for v in values), default=0)
    out = torch.full((len(values), max_len), pad_value, dtype=values[0].dtype if values else torch.long)
    for i, value in enumerate(values):
        out[i, : value.shape[0]] = value
    return out


def _pad_2d(values: List[torch.Tensor], pad_value: int | float = 0) -> torch.Tensor:
    max_len = max((v.shape[0] for v in values), default=0)
    width = values[0].shape[1] if values else 0
    dtype = values[0].dtype if values else torch.long
    out = torch.full((len(values), max_len, width), pad_value, dtype=dtype)
    for i, value in enumerate(values):
        out[i, : value.shape[0], :] = value
    return out


def collate_variant_tokens(batch: List[Dict]) -> Dict[str, torch.Tensor | List[str]]:
    if not batch:
        raise ValueError("Batch vazio")
    lengths = torch.tensor([int(item["variant_type"].shape[0]) for item in batch], dtype=torch.long)
    max_len = int(lengths.max().item()) if len(lengths) else 0
    attention_mask = torch.zeros((len(batch), max_len + 1), dtype=torch.bool)
    attention_mask[:, 0] = True
    for i, length in enumerate(lengths.tolist()):
        attention_mask[i, 1: 1 + length] = True
    l_max = 0
    for item in batch:
        if item["ref_allele"].ndim == 2 and item["ref_allele"].shape[1] > 0:
            l_max = int(item["ref_allele"].shape[1])
            break
    if l_max == 0:
        l_max = 16
    ref_values = [item["ref_allele"] if item["ref_allele"].shape[1] else torch.empty((0, l_max), dtype=torch.long) for item in batch]
    alt_values = [item["alt_allele"] if item["alt_allele"].shape[1] else torch.empty((0, l_max), dtype=torch.long) for item in batch]
    return {
        "variant_type": _pad_1d([item["variant_type"] for item in batch]),
        "haplotype": _pad_1d([item["haplotype"] for item in batch]),
        "gene": _pad_1d([item["gene"] for item in batch]),
        "length_norm": _pad_2d([item["length_norm"] for item in batch], pad_value=0.0),
        "position_relative": _pad_1d([item["position_relative"] for item in batch]),
        "position": _pad_1d([item["position"] for item in batch]),
        "ref_allele": _pad_2d(ref_values, pad_value=5),
        "alt_allele": _pad_2d(alt_values, pad_value=5),
        "attention_mask": attention_mask,
        "targets": torch.stack([item["target"] for item in batch]),
        "lengths": lengths,
        "sample_ids": [str(item.get("sample_id", "")) for item in batch],
    }
