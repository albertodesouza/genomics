from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional


@dataclass(frozen=True)
class DatasetRef:
    dataset_id: str
    version: str
    path: Path
    role: str
    metadata_path: Path
    layout_version: Optional[int] = None
    source_dataset_id: Optional[str] = None


def data_root() -> Path:
    return Path(os.environ.get("GENOMICS_DATA_ROOT", "/dados/GENOMICS_DATA")).resolve()


def _dataset_defs() -> Dict[str, Dict[str, object]]:
    root = data_root()
    return {
        "1kg_high_coverage": {
            "version": "v1",
            "path": root / "v1" / "1kG_high_coverage",
            "role": "canonical",
            "layout_version": 1,
        },
        "legacy_top3_1kg_high_coverage": {
            "version": "legacy",
            "path": root / "top3" / "non_longevous_results_genes_1000_all",
            "role": "legacy",
            "layout_version": None,
            "source_dataset_id": "1kg_high_coverage",
        },
        "variant_transformer_superpopulation": {
            "version": "current",
            "path": root / "variant_transformer" / "superpopulation",
            "role": "derived",
            "layout_version": None,
            "source_dataset_id": "1kg_high_coverage",
        },
        "variant_transformer_superpopulation_32k": {
            "version": "current",
            "path": root / "variant_transformer" / "superpopulation_32k",
            "role": "derived",
            "layout_version": None,
            "source_dataset_id": "1kg_high_coverage",
        },
        "variant_transformer_pigmentation_binary": {
            "version": "current",
            "path": root / "variant_transformer" / "pigmentation_binary",
            "role": "derived",
            "layout_version": None,
            "source_dataset_id": "1kg_high_coverage",
        },
    }


def registered_dataset_ids() -> tuple[str, ...]:
    return tuple(sorted(_dataset_defs()))


def resolve_dataset(dataset_id: str, version: str = "current") -> DatasetRef:
    defs = _dataset_defs()
    if dataset_id not in defs:
        known = ", ".join(sorted(defs))
        raise KeyError(f"Dataset id desconhecido: {dataset_id}. Conhecidos: {known}")
    item = defs[dataset_id]
    item_version = str(item["version"])
    if version not in ("current", item_version):
        raise KeyError(f"Versao '{version}' indisponivel para {dataset_id}; versao registrada: {item_version}")
    path = Path(item["path"]).resolve()
    return DatasetRef(
        dataset_id=dataset_id,
        version=item_version,
        path=path,
        role=str(item["role"]),
        metadata_path=path / "dataset_metadata.json",
        layout_version=item.get("layout_version"),
        source_dataset_id=item.get("source_dataset_id"),
    )


def resolve_derived_dataset(name: str) -> DatasetRef:
    if name.startswith("variant_transformer_"):
        return resolve_dataset(name)
    return resolve_dataset(f"variant_transformer_{name}")


def canonical_dataset_ref() -> DatasetRef:
    return resolve_dataset("1kg_high_coverage")


def legacy_top3_ref() -> DatasetRef:
    return resolve_dataset("legacy_top3_1kg_high_coverage")


def is_same_or_child(path: Path, root: Path) -> bool:
    path = Path(path).resolve()
    root = Path(root).resolve()
    return path == root or root in path.parents


def classify_dataset_path(path: Optional[str | Path]) -> str:
    if not path:
        return "missing-data-path"
    path = Path(path).resolve()
    defs = [resolve_dataset(dataset_id) for dataset_id in registered_dataset_ids()]
    for ref in defs:
        if is_same_or_child(path, ref.path):
            if ref.role == "canonical":
                return "canonical-ok"
            if ref.role == "legacy":
                return "legacy-top3"
            if ref.role == "derived":
                return "derived-ok"
            return ref.role
    if is_same_or_child(path, data_root() / "top3"):
        return "legacy-top3"
    if is_same_or_child(path, data_root()):
        return "dados-unregistered"
    return "external-data-path"


def is_canonical_dataset(path: Path) -> bool:
    return is_same_or_child(path, canonical_dataset_ref().path)


def is_legacy_top3_path(path: Path) -> bool:
    return is_same_or_child(path, data_root() / "top3")
