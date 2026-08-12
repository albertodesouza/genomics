from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from genomics.predictors.genotype_based.alignment.bcftools_chain_mapper import BcftoolsChainMapper
from genomics.predictors.genotype_based.config import PipelineConfig, load_config
from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner
from genomics.predictors.genotype_based.data.genomic_dataset import GenomicDataset


console = Console()

_WORKER_DATASET_DIR: Optional[Path] = None
_WORKER_CONSENSUS_DIR: Optional[Path] = None
_WORKER_SAMPLE_IDS: Optional[List[str]] = None
_WORKER_CENTER_WINDOW_SIZE: int = 32768
_WORKER_ALIGNER: Optional[DynamicIndelAligner] = None
_WORKER_MAPPER: Optional[BcftoolsChainMapper] = None


def _resolve_sample_ids(config: PipelineConfig, base_dataset: GenomicDataset, start: int, limit: Optional[int], all_samples: bool) -> List[str]:
    individuals = [str(x) for x in getattr(base_dataset, "dataset_metadata", {}).get("individuals", [])]
    selected = set(config.dataset_input.sample_ids or [])
    if selected:
        individuals = [sid for sid in individuals if sid in selected]
    start = max(int(start), 0)
    if all_samples:
        return individuals[start:]
    if limit is None or limit <= 0:
        return individuals[start:]
    return individuals[start:start + int(limit)]


def _resolve_genes(config: PipelineConfig, base_dataset: GenomicDataset, genes_arg: Optional[str]) -> List[str]:
    if genes_arg:
        requested = [item.strip() for item in genes_arg.split(",") if item.strip()]
        return requested
    genes = list(config.dataset_input.genes_to_use or [])
    if genes:
        return genes
    return [str(x) for x in getattr(base_dataset, "dataset_metadata", {}).get("genes", [])]


def _chunks(items: Sequence[Tuple[str, str]], size: int) -> List[List[Tuple[str, str]]]:
    size = max(int(size), 1)
    return [list(items[i:i + size]) for i in range(0, len(items), size)]


def _init_worker(dataset_dir: str, consensus_dir: str, sample_ids: List[str], center_window_size: int) -> None:
    global _WORKER_DATASET_DIR, _WORKER_CONSENSUS_DIR, _WORKER_SAMPLE_IDS, _WORKER_CENTER_WINDOW_SIZE
    global _WORKER_ALIGNER, _WORKER_MAPPER
    _WORKER_DATASET_DIR = Path(dataset_dir)
    _WORKER_CONSENSUS_DIR = Path(consensus_dir)
    _WORKER_SAMPLE_IDS = list(sample_ids)
    _WORKER_CENTER_WINDOW_SIZE = int(center_window_size)
    _WORKER_ALIGNER = None
    _WORKER_MAPPER = None


def _get_worker_mapper() -> Tuple[DynamicIndelAligner, BcftoolsChainMapper]:
    global _WORKER_ALIGNER, _WORKER_MAPPER
    if _WORKER_DATASET_DIR is None or _WORKER_CONSENSUS_DIR is None or _WORKER_SAMPLE_IDS is None:
        raise RuntimeError("Worker nao inicializado")
    if _WORKER_ALIGNER is None:
        _WORKER_ALIGNER = DynamicIndelAligner(
            _WORKER_DATASET_DIR,
            selected_sample_ids=_WORKER_SAMPLE_IDS,
            center_window_size=_WORKER_CENTER_WINDOW_SIZE,
            sample_cache_limit=2,
            sample_prefetch_size=1,
        )
    if _WORKER_MAPPER is None:
        _WORKER_MAPPER = BcftoolsChainMapper(
            dataset_dir=_WORKER_DATASET_DIR,
            consensus_dataset_dir=_WORKER_CONSENSUS_DIR,
            aligner=_WORKER_ALIGNER,
        )
    return _WORKER_ALIGNER, _WORKER_MAPPER


def _process_chunk(chunk: List[Tuple[str, str]]) -> Dict[str, object]:
    aligner, mapper = _get_worker_mapper()
    done = 0
    failures = []
    built_axes = set()
    for gene, sample_id in chunk:
        try:
            if gene not in built_axes:
                aligner.build_alignment_axis_for_gene(gene, _WORKER_SAMPLE_IDS or [])
                built_axes.add(gene)
            for haplotype in ("H1", "H2"):
                entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
                if entry is None:
                    raise RuntimeError("entry ausente")
                done += 1
        except Exception as exc:
            failures.append({"gene": gene, "sample_id": sample_id, "error": str(exc)})
    mapper._entry_cache.clear()
    for gene, _sample_id in chunk:
        cache = aligner._sample_entry_cache.get(gene)
        if cache is not None:
            cache.clear()
    return {"done": done, "failures": failures, "pid": os.getpid()}


def precompute(
    config_path: Path,
    workers: int,
    chunk_size: int,
    start: int,
    limit: Optional[int],
    all_samples: bool,
    genes_arg: Optional[str],
) -> Dict[str, object]:
    config = load_config(config_path)
    if config.dataset_input.alignment_mapping != "bcftools_chain":
        raise ValueError("A config deve usar dataset_input.alignment_mapping='bcftools_chain'")
    if not config.dataset_input.consensus_dataset_dir:
        raise ValueError("dataset_input.consensus_dataset_dir e obrigatorio")

    dataset_dir = Path(config.dataset_input.dataset_dir)
    consensus_dir = Path(config.dataset_input.consensus_dataset_dir)
    base_dataset = GenomicDataset(dataset_dir=dataset_dir, load_predictions=False, load_sequences=False, cache_sequences=False)
    sample_ids = _resolve_sample_ids(config, base_dataset, start=start, limit=limit, all_samples=all_samples)
    genes = _resolve_genes(config, base_dataset, genes_arg)
    if not sample_ids:
        raise ValueError("Nenhum sample_id selecionado")
    if not genes:
        raise ValueError("Nenhum gene selecionado")

    pairs = [(gene, sample_id) for gene in genes for sample_id in sample_ids]
    chunks = _chunks(pairs, chunk_size)
    total_entries = len(pairs) * 2
    workers = max(int(workers), 1)

    console.print(f"[cyan]Precomputando bcftools_chain[/cyan]: {len(genes)} genes, {len(sample_ids)} samples, {total_entries} entries")
    console.print(f"[cyan]Workers:[/cyan] {workers} | [cyan]chunk_size:[/cyan] {chunk_size}")

    failures = []
    completed = 0
    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("bcftools_chain entries", total=total_entries)
        if workers == 1:
            _init_worker(str(dataset_dir), str(consensus_dir), sample_ids, config.dataset_input.window_center_size)
            for chunk in chunks:
                result = _process_chunk(chunk)
                completed += int(result["done"])
                failures.extend(result["failures"])
                progress.update(task, advance=int(result["done"]))
        else:
            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_init_worker,
                initargs=(str(dataset_dir), str(consensus_dir), sample_ids, config.dataset_input.window_center_size),
            ) as executor:
                chunk_iter = iter(chunks)
                pending = set()
                max_pending = max(workers * 2, 1)
                for _ in range(min(max_pending, len(chunks))):
                    pending.add(executor.submit(_process_chunk, next(chunk_iter)))
                while pending:
                    done_futures, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done_futures:
                        result = future.result()
                        completed += int(result["done"])
                        failures.extend(result["failures"])
                        progress.update(task, advance=int(result["done"]))
                        try:
                            pending.add(executor.submit(_process_chunk, next(chunk_iter)))
                        except StopIteration:
                            pass

    return {
        "ok": not failures,
        "config": str(config_path),
        "genes": genes,
        "sample_count": len(sample_ids),
        "expected_entries": total_entries,
        "completed_entries": completed,
        "failure_count": len(failures),
        "failures": failures[:50],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Precompute bcftools_chain mapper cache in parallel")
    parser.add_argument("config", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--chunk-size", type=int, default=25)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--all-samples", action="store_true")
    parser.add_argument("--genes", default=None, help="Comma-separated genes; default uses config genes")
    args = parser.parse_args()
    result = precompute(args.config.resolve(), args.workers, args.chunk_size, args.start, args.limit, args.all_samples, args.genes)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if not result.get("ok"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
