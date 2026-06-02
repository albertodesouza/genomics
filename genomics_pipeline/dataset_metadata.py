from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .config_io import load_json


def load_dataset_metadata(dataset_dir: Path) -> Dict[str, Any]:
    dataset_dir = Path(dataset_dir)
    metadata_path = dataset_dir / "dataset_metadata.json"
    if not metadata_path.exists():
        raise FileNotFoundError(f"dataset_metadata.json nao encontrado: {metadata_path}")
    return load_json(metadata_path)


def load_layout_metadata(dataset_dir: Path) -> Dict[str, Any]:
    layout_path = Path(dataset_dir) / "layout_metadata.json"
    if not layout_path.exists():
        return {}
    return load_json(layout_path)


def get_individuals(metadata: Dict[str, Any]) -> List[str]:
    return [str(sample_id) for sample_id in metadata.get("individuals", [])]


def get_pedigree(metadata: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    pedigree = metadata.get("individuals_pedigree", {}) or {}
    return {str(k): dict(v or {}) for k, v in pedigree.items()}


def get_gene_order(metadata: Dict[str, Any]) -> List[str]:
    return [str(gene) for gene in metadata.get("genes", [])]


def load_window_catalog(dataset_dir: Path, metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Dict[str, Any]]:
    dataset_dir = Path(dataset_dir)
    metadata = metadata if metadata is not None else load_dataset_metadata(dataset_dir)
    catalog = dict(metadata.get("window_catalog", {}) or {})
    if catalog:
        return {str(k): dict(v or {}) for k, v in catalog.items()}

    ref_windows = dataset_dir / "references" / "windows"
    for window_meta_path in sorted(ref_windows.glob("*/window_metadata.json")):
        catalog[window_meta_path.parent.name] = load_json(window_meta_path)
    return {str(k): dict(v or {}) for k, v in catalog.items()}


def iter_sample_metadata(dataset_dir: Path, metadata: Optional[Dict[str, Any]] = None) -> Iterable[Dict[str, Any]]:
    dataset_dir = Path(dataset_dir)
    metadata = metadata if metadata is not None else load_dataset_metadata(dataset_dir)
    pedigree = get_pedigree(metadata)
    for sample_id in get_individuals(metadata):
        row = dict(pedigree.get(sample_id, {}) or {})
        row.setdefault("sample_id", sample_id)
        individual_metadata_path = dataset_dir / "individuals" / sample_id / "individual_metadata.json"
        if individual_metadata_path.exists():
            loaded = load_json(individual_metadata_path)
            loaded.update(row)
            row = loaded
            row.setdefault("sample_id", sample_id)
        yield row


def get_sample_target(sample_metadata: Dict[str, Any], target: str) -> Optional[str]:
    value = sample_metadata.get(target)
    if value is None:
        value = sample_metadata.get("target")
    if value in (None, ""):
        return None
    return str(value)


def load_vcf_sources_from_dataset(dataset_dir: Optional[Path]) -> Dict[str, str]:
    if dataset_dir is None:
        return {}
    dataset_dir = Path(dataset_dir)
    metadata_path = dataset_dir / "dataset_metadata.json"
    if not metadata_path.exists():
        return {}
    catalog = load_window_catalog(dataset_dir, load_json(metadata_path))
    sources: Dict[str, str] = {}
    for gene_id, item in catalog.items():
        raw_source = item.get("raw_variant_source") or {}
        vcf_path = raw_source.get("vcf_path")
        if not vcf_path:
            continue
        sources[str(gene_id)] = str(vcf_path)
        chrom = raw_source.get("chromosome") or item.get("chromosome") or item.get("chrom")
        if chrom:
            chrom_text = str(chrom)
            stripped = chrom_text[3:] if chrom_text.startswith("chr") else chrom_text
            sources.setdefault(chrom_text, str(vcf_path))
            sources.setdefault(stripped, str(vcf_path))
    return sources
