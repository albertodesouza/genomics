from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .config_io import load_json


@dataclass(frozen=True)
class SampleRecord:
    sample_id: str
    target: Optional[str] = None
    family_id: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


@dataclass(frozen=True)
class SplitSpec:
    train_split: float = 0.7
    val_split: float = 0.15
    test_split: float = 0.15
    random_seed: Optional[int] = 13
    family_split_mode: str = "family_aware"
    balancing_strategy: str = "shuffle"


def extract_family_links(pedigree: Dict[str, Any]) -> List[str]:
    related: List[str] = []
    for key, value in pedigree.items():
        key_lower = str(key).lower()
        if not any(token in key_lower for token in ["father", "mother", "parent", "child", "sibling", "family"]):
            continue
        if value in [None, "", "0"]:
            continue
        if isinstance(value, (list, tuple, set)):
            related.extend(str(v) for v in value if v not in [None, "", "0"])
        else:
            related.append(str(value))
    return related


def _singleton_groups(sample_ids: Sequence[str]) -> List[List[int]]:
    return [[idx] for idx in range(len(sample_ids))]


def _group_info(groups: List[List[int]], family_mode: str, source: str) -> Dict[str, Any]:
    return {
        "family_split_mode": family_mode,
        "grouping_source": source,
        "num_groups": len(groups),
        "num_family_groups": sum(1 for group in groups if len(group) > 1),
        "num_singletons": sum(1 for group in groups if len(group) == 1),
    }


def build_family_groups(
    *,
    sample_ids: Sequence[str],
    pedigree_map: Optional[Dict[str, Dict[str, Any]]] = None,
    individuals_dir: Optional[Path] = None,
    family_split_mode: str = "family_aware",
) -> Tuple[List[List[int]], Dict[str, Any]]:
    """Build groups of sample indices that must stay in the same split."""
    sample_ids = [str(sample_id) for sample_id in sample_ids]
    if not sample_ids:
        return [], _group_info([], family_split_mode, "individual")

    if family_split_mode == "ignore":
        groups = _singleton_groups(sample_ids)
        return groups, _group_info(groups, family_split_mode, "individual")

    sample_to_idx = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    family_ids: Dict[str, str] = {}

    if individuals_dir is not None:
        individuals_dir = Path(individuals_dir)
        for sample_id in sample_ids:
            metadata_file = individuals_dir / sample_id / "individual_metadata.json"
            if not metadata_file.exists():
                continue
            try:
                family_id = load_json(metadata_file).get("family_id")
            except Exception:
                continue
            if family_id not in [None, "", "0"]:
                family_ids[sample_id] = str(family_id)

    if family_ids:
        groups_by_family_id: Dict[str, List[int]] = {}
        for sample_id, idx in sample_to_idx.items():
            family_id = family_ids.get(sample_id, sample_id)
            groups_by_family_id.setdefault(family_id, []).append(idx)
        groups = list(groups_by_family_id.values())
        return groups, _group_info(groups, family_split_mode, "family_id")

    pedigree_map = pedigree_map or {}
    if not pedigree_map:
        groups = _singleton_groups(sample_ids)
        return groups, _group_info(groups, family_split_mode, "individual")

    parent = {sample_id: sample_id for sample_id in sample_ids}

    def find(sample_id: str) -> str:
        while parent[sample_id] != sample_id:
            parent[sample_id] = parent[parent[sample_id]]
            sample_id = parent[sample_id]
        return sample_id

    def union(left: str, right: str) -> None:
        if left not in parent or right not in parent:
            return
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for sample_id, pedigree in pedigree_map.items():
        sample_id = str(sample_id)
        if sample_id not in parent:
            continue
        for related_id in extract_family_links(pedigree or {}):
            if related_id in parent:
                union(sample_id, related_id)

    groups_by_root: Dict[str, List[int]] = {}
    for sample_id, idx in sample_to_idx.items():
        groups_by_root.setdefault(find(sample_id), []).append(idx)
    groups = list(groups_by_root.values())
    return groups, _group_info(groups, family_split_mode, "pedigree")


def build_sample_record_groups(records: Sequence[SampleRecord], family_split_mode: str = "family_aware") -> Tuple[List[List[int]], Dict[str, Any]]:
    if family_split_mode == "ignore":
        groups = _singleton_groups([record.sample_id for record in records])
        return groups, _group_info(groups, family_split_mode, "individual")

    groups_by_key: Dict[str, List[int]] = {}
    for idx, record in enumerate(records):
        family_id = record.family_id
        key = family_id if family_id not in (None, "", "0") else record.sample_id
        groups_by_key.setdefault(str(key), []).append(idx)
    groups = list(groups_by_key.values())
    source = "family_id" if any(len(group) > 1 for group in groups) else "individual"
    return groups, _group_info(groups, family_split_mode, source)


def split_groups(groups: Sequence[Sequence[int]], spec: SplitSpec) -> Dict[str, List[int]]:
    shuffled = [list(group) for group in groups]
    rng = random.Random(spec.random_seed)
    rng.shuffle(shuffled)
    n_train = int(spec.train_split * len(shuffled))
    n_val = int(spec.val_split * len(shuffled))
    return {
        "train": [idx for group in shuffled[:n_train] for idx in group],
        "val": [idx for group in shuffled[n_train:n_train + n_val] for idx in group],
        "test": [idx for group in shuffled[n_train + n_val:] for idx in group],
    }


def split_sample_records(records: Sequence[SampleRecord], spec: SplitSpec) -> Tuple[Dict[str, List[str]], Dict[str, Any]]:
    groups, info = build_sample_record_groups(records, spec.family_split_mode)
    split_indices = split_groups(groups, spec)
    return {
        split: [records[idx].sample_id for idx in indices]
        for split, indices in split_indices.items()
    }, info


def records_from_objects(samples: Iterable[Any]) -> List[SampleRecord]:
    records: List[SampleRecord] = []
    for sample in samples:
        if isinstance(sample, dict):
            sample_id = sample.get("sample_id") or sample.get("id")
            target = sample.get("target")
            family_id = sample.get("family_id")
            metadata = dict(sample)
        else:
            sample_id = getattr(sample, "sample_id")
            target = getattr(sample, "target", None)
            family_id = getattr(sample, "family_id", None)
            metadata = getattr(sample, "metadata", None)
        records.append(SampleRecord(str(sample_id), None if target is None else str(target), None if family_id is None else str(family_id), metadata))
    return records
