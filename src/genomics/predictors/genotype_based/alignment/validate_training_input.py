from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

from genomics.predictors.genotype_based.alignment.bcftools_chain_mapper import BcftoolsChainMapper, read_fasta_sequence
from genomics.predictors.genotype_based.config import PipelineConfig, load_config
from genomics.predictors.genotype_based.data.processed_dataset import ProcessedGenomicDataset
from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner
from genomics.predictors.genotype_based.data.genomic_dataset import GenomicDataset


def _resolve_sample_ids(config: PipelineConfig, base_dataset: GenomicDataset, sample_limit: Optional[int], all_samples: bool) -> List[str]:
    individuals = list(getattr(base_dataset, "dataset_metadata", {}).get("individuals", []))
    selected = []
    explicit = set(config.dataset_input.sample_ids or [])
    if explicit:
        selected = [sid for sid in individuals if sid in explicit]
    else:
        selected = [str(sid) for sid in individuals]
    if all_samples:
        return selected
    limit = sample_limit if sample_limit is not None and sample_limit > 0 else 5
    return selected[:limit]


def _prediction_array(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Prediction NPZ ausente: {path}")
    with np.load(path) as data:
        if "values" in data:
            return np.asarray(data["values"])
        keys = [k for k in data.files if not k.endswith("_metadata")]
        if not keys:
            raise ValueError(f"NPZ sem array de predicao: {path}")
        return np.asarray(data[keys[0]])


def _validate_entry_bounds(entry: Dict[str, object], prediction_length: int, expanded_length: int) -> Dict[str, int]:
    copy_from = [int(v) for v in entry.get("copy_from_indices", [])]
    copy_to = [int(v) for v in entry.get("expanded_indices", [])]
    if len(copy_from) != len(copy_to):
        raise ValueError(f"copy_from_indices e expanded_indices tem tamanhos diferentes: {len(copy_from)} != {len(copy_to)}")
    source_oob = sum(1 for v in copy_from if v < 0 or v >= prediction_length)
    target_oob = sum(1 for v in copy_to if v < 0 or v >= expanded_length)
    duplicate_targets = len(copy_to) - len(set(copy_to))
    if source_oob or target_oob or duplicate_targets:
        raise ValueError(
            "Entry invalida: "
            f"source_oob={source_oob}, target_oob={target_oob}, duplicate_targets={duplicate_targets}"
        )
    return {
        "copy_positions": len(copy_from),
        "source_oob": source_oob,
        "target_oob": target_oob,
        "duplicate_targets": duplicate_targets,
        "insertions": len(set(int(v) for v in entry.get("insertion_indices", []))),
        "deletions": len(set(int(v) for v in entry.get("deletion_indices", []))),
    }


def _validate_haplotype(
    *,
    mapper: BcftoolsChainMapper,
    consensus_root: Path,
    dataset_dir: Path,
    gene: str,
    sample_id: str,
    haplotype: str,
    output_names: Sequence[str],
    expanded_length: int,
) -> Dict[str, object]:
    case_dir = consensus_root / "individuals" / sample_id / "windows" / gene
    ref_path = consensus_root / "references" / "windows" / gene / "ref.window.fa"
    fixed_path = case_dir / f"{sample_id}.{haplotype}.window.fixed.fa"
    raw_path = case_dir / f"{sample_id}.{haplotype}.window.raw.fa"
    vcf_path = case_dir / f"{sample_id}.window.consensus_ready.vcf.gz"
    required = [ref_path, fixed_path, raw_path, vcf_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Arquivos obrigatorios ausentes para {gene}/{sample_id}/{haplotype}: {missing}")

    fixed_len = len(read_fasta_sequence(fixed_path))
    entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
    if entry is None:
        raise RuntimeError(f"Entry ausente para {gene}/{sample_id}/{haplotype}")
    if entry.get("mapping_method") != "bcftools_chain":
        raise ValueError(f"Entry nao usa bcftools_chain: {entry.get('mapping_method')}")

    npz_summaries = []
    for output in output_names:
        npz_path = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}" / f"{output}.npz"
        array = _prediction_array(npz_path)
        if int(array.shape[0]) != fixed_len:
            raise ValueError(f"NPZ length != fixed FASTA length para {npz_path}: {array.shape[0]} != {fixed_len}")
        bounds = _validate_entry_bounds(entry, int(array.shape[0]), expanded_length)
        finite_count = int(np.isfinite(array).sum())
        if finite_count == 0:
            raise ValueError(f"NPZ sem valores finitos: {npz_path}")
        npz_summaries.append({
            "output": output,
            "shape": list(array.shape),
            "finite_count": finite_count,
            **bounds,
        })

    return {
        "gene": gene,
        "sample_id": sample_id,
        "haplotype": haplotype,
        "fixed_len": fixed_len,
        "copy_positions": len(entry.get("copy_from_indices", [])),
        "insertions": len(entry.get("insertion_indices", [])),
        "deletions": len(entry.get("deletion_indices", [])),
        "npz": npz_summaries,
    }


def _validate_tensors(processed: ProcessedGenomicDataset, max_items: int) -> Dict[str, object]:
    shape_counts: Counter = Counter()
    checked = 0
    for idx in range(min(max_items, len(processed))):
        features, target = processed[idx]
        if not isinstance(features, torch.Tensor):
            raise TypeError(f"features nao e torch.Tensor no item {idx}: {type(features)}")
        if torch.isnan(features).any():
            raise ValueError(f"Tensor contem NaN no item {idx}")
        if torch.isinf(features).any():
            raise ValueError(f"Tensor contem Inf no item {idx}")
        if torch.count_nonzero(features).item() == 0:
            raise ValueError(f"Tensor todo zero no item {idx}")
        shape_counts[tuple(int(v) for v in features.shape)] += 1
        checked += 1
    return {"items_checked": checked, "shape_counts": {str(k): v for k, v in shape_counts.items()}}


def validate(config_path: Path, sample_limit: Optional[int], all_samples: bool, max_tensor_items: int) -> Dict[str, object]:
    config = load_config(config_path)
    if config.dataset_input.alignment_mapping != "bcftools_chain":
        raise ValueError("Este validador espera dataset_input.alignment_mapping='bcftools_chain'")
    if not config.dataset_input.consensus_dataset_dir:
        raise ValueError("dataset_input.consensus_dataset_dir e obrigatorio")

    dataset_dir = Path(config.dataset_input.dataset_dir)
    consensus_root = Path(config.dataset_input.consensus_dataset_dir)
    base_dataset = GenomicDataset(dataset_dir=dataset_dir, load_predictions=True, load_sequences=False, cache_sequences=False)
    sample_ids = _resolve_sample_ids(config, base_dataset, sample_limit, all_samples)
    genes = list(config.dataset_input.genes_to_use or [])
    if not genes:
        genes = list(getattr(base_dataset, "dataset_metadata", {}).get("genes", []))
    output_names = [str(name).split(".")[-1].lower() for name in config.dataset_input.alphagenome_outputs]

    aligner = DynamicIndelAligner(dataset_dir, selected_sample_ids=sample_ids, center_window_size=config.dataset_input.window_center_size)
    mapper = BcftoolsChainMapper(dataset_dir=dataset_dir, consensus_dataset_dir=consensus_root, aligner=aligner)

    hap_summary = []
    axis_summary = {}
    failures = []
    for gene in genes:
        aligner.build_alignment_axis_for_gene(gene, sample_ids)
        axis = aligner.get_alignment_axis(gene)
        expanded_length = int(axis["expanded_length"])
        axis_summary[gene] = {
            "expanded_length": expanded_length,
            "ref_length": int(axis.get("ref_length", 0)),
            "ref_start_offset": int(axis.get("ref_start_offset", 0)),
        }
        for sample_id in sample_ids:
            for haplotype in ("H1", "H2"):
                try:
                    hap_summary.append(_validate_haplotype(
                        mapper=mapper,
                        consensus_root=consensus_root,
                        dataset_dir=dataset_dir,
                        gene=gene,
                        sample_id=sample_id,
                        haplotype=haplotype,
                        output_names=output_names,
                        expanded_length=expanded_length,
                    ))
                except Exception as exc:
                    failures.append({"gene": gene, "sample_id": sample_id, "haplotype": haplotype, "error": str(exc)})

    if failures:
        return {
            "ok": False,
            "config": str(config_path),
            "samples_checked": len(sample_ids),
            "genes_checked": len(genes),
            "failures": failures[:50],
            "failure_count": len(failures),
        }

    processed = ProcessedGenomicDataset(base_dataset, config, compute_normalization=False)
    processed.prepare_alignment_cache(sample_ids)
    tensor_summary = _validate_tensors(processed, max_tensor_items)

    copy_counts = [int(row["copy_positions"]) for row in hap_summary]
    return {
        "ok": True,
        "config": str(config_path),
        "samples_checked": len(sample_ids),
        "genes_checked": len(genes),
        "haplotypes_checked": len(hap_summary),
        "outputs": list(output_names),
        "axis_summary": axis_summary,
        "copy_positions": {
            "min": min(copy_counts) if copy_counts else None,
            "max": max(copy_counts) if copy_counts else None,
        },
        "tensor_summary": tensor_summary,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate bcftools_chain training inputs before CNN training")
    parser.add_argument("config", type=Path)
    parser.add_argument("--sample-limit", type=int, default=5)
    parser.add_argument("--all-samples", action="store_true")
    parser.add_argument("--max-tensor-items", type=int, default=3)
    args = parser.parse_args()
    result = validate(args.config.resolve(), args.sample_limit, args.all_samples, args.max_tensor_items)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if not result.get("ok"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
