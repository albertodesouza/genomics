# -*- coding: utf-8 -*-
"""
data_pipeline.py — Orquestra cache, DataLoaders e prepare_data().
"""
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import torch
from rich.console import Console
from torch.utils.data import DataLoader, Subset

from genotype_based_predictor.config import PipelineConfig, get_dataset_cache_dir, generate_dataset_name
from genotype_based_predictor.dataset import CachedProcessedDataset, ProcessedGenomicDataset
from genotype_based_predictor.dataset_layout import is_dataset_dir
from genotype_based_predictor.data_splitting import (
    build_family_aware_sample_groups,
    build_valid_sample_index_map,
)
from genotype_based_predictor.utils import worker_init_fn

console = Console()


def _load_external_normalization_params(config: PipelineConfig) -> Optional[Dict]:
    norm_path = config.dataset_input.normalization_params_path
    if not norm_path:
        return None
    norm_file = Path(norm_path)
    if not norm_file.exists():
        raise FileNotFoundError(f"normalization_params.json não encontrado: {norm_file}")
    with open(norm_file) as f:
        params = json.load(f)
    console.print(f"[green]Usando normalização pré-computada:[/green] {norm_file}")
    return params


def _build_view_definition(config: PipelineConfig) -> Dict[str, Any]:
    di = config.dataset_input
    return {
        "dataset_dir": str(Path(di.dataset_dir).resolve()),
        "sample_ids": di.sample_ids,
        "sample_ids_path": di.sample_ids_path,
        "superpopulations_to_use": di.superpopulations_to_use,
        "populations_to_use": di.populations_to_use,
        "genes_to_use": di.genes_to_use,
        "alphagenome_outputs": di.alphagenome_outputs,
        "haplotype_mode": di.haplotype_mode,
        "window_center_size": di.window_center_size,
        "downsample_factor": di.downsample_factor,
        "normalization_method": di.normalization_method,
        "selected_track_index": di.selected_track_index,
        "indel_include_valid_mask": di.indel_include_valid_mask,
        "tensor_layout": di.tensor_layout,
    }


def _build_resolved_view_definition(config: PipelineConfig, processed_dataset: ProcessedGenomicDataset) -> Dict[str, Any]:
    view = _build_view_definition(config)
    sample_ids = []
    for base_idx in processed_dataset.valid_sample_indices:
        if base_idx < len(processed_dataset.individuals):
            sample_ids.append(processed_dataset.individuals[base_idx])
    view["resolved_sample_ids"] = sample_ids
    view["resolved_genes"] = sorted(processed_dataset.genes_to_use) if processed_dataset.genes_to_use else None
    return view


def _resolve_runtime_dataset_dir(config: PipelineConfig) -> Path:
    dataset_dir = Path(config.dataset_input.dataset_dir)
    if is_dataset_dir(dataset_dir):
        return dataset_dir

    raise FileNotFoundError(
        f"dataset_input.dataset_dir deve apontar para um dataset válido: {dataset_dir}"
    )


def _load_base_dataset(runtime_dataset_dir: Path):
    if not is_dataset_dir(runtime_dataset_dir):
        raise FileNotFoundError(f"Dataset inválido: {runtime_dataset_dir}")

    from genotype_based_predictor.genomic_dataset import GenomicDataset

    return GenomicDataset(
        dataset_dir=runtime_dataset_dir,
        load_predictions=True,
        load_sequences=False,
        cache_sequences=False,
    )


# ---------------------------------------------------------------------------
# Cache validation
# ---------------------------------------------------------------------------

def validate_cache(cache_dir: Path, config: PipelineConfig) -> bool:
    """Verifica se o cache é válido e compatível com a configuração atual."""
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return False

    if config.debug.taint_at_cache_save and config.mode == "train":
        console.print("[yellow]Cache recriado: taint_at_cache_save=True[/yellow]")
        return False

    if not (cache_dir / ".cache_complete").exists():
        console.print("[yellow]Cache incompleto (processo interrompido)[/yellow]")
        return False

    required_files = ["metadata.json", "normalization_params.json", "split_index.json", "view_definition.json"]
    if config.dataset_input.cache_processed_tensors:
        required_files.extend(["train_data.pt", "val_data.pt", "test_data.pt"])
    for f in required_files:
        if not (cache_dir / f).exists():
            console.print(f"[yellow]Cache incompleto: falta {f}[/yellow]")
            return False

    try:
        with open(cache_dir / "metadata.json") as f:
            meta = json.load(f)

        if meta.get("format_version", 1) < 2:
            console.print("[yellow]Cache legado (v1) — recriando...[/yellow]")
            return False

        with open(cache_dir / "view_definition.json") as f:
            cache_view = json.load(f)
        if cache_view.get("requested_view") != _build_view_definition(config):
            console.print("[yellow]Cache inválido: definição lógica da view mudou[/yellow]")
            return False

        pp = meta.get("processing_params", {})
        checks = {
            "alphagenome_outputs": config.dataset_input.alphagenome_outputs,
            "haplotype_mode": config.dataset_input.haplotype_mode,
            "window_center_size": config.dataset_input.window_center_size,
            "downsample_factor": config.dataset_input.downsample_factor,
            "normalization_method": config.dataset_input.normalization_method,
            "selected_track_index": config.dataset_input.selected_track_index,
            "indel_include_valid_mask": config.dataset_input.indel_include_valid_mask,
            "tensor_layout": config.dataset_input.tensor_layout,
            "prediction_target": config.output.prediction_target,
            "input_shape": "3D_haplotype_channels",
        }
        for k, v in checks.items():
            if pp.get(k) != v:
                console.print(f"[yellow]Cache inválido: {k} mudou ({pp.get(k)} → {v})[/yellow]")
                return False

        if meta.get("splits", {}).get("random_seed") != config.data_split.random_seed:
            console.print("[yellow]Cache inválido: random_seed diferente[/yellow]")
            return False

        console.print("[green]✓ Cache válido e compatível[/green]")
        return True
    except Exception as e:
        console.print(f"[yellow]Erro ao validar cache: {e}[/yellow]")
        return False


# ---------------------------------------------------------------------------
# Gene window metadata
# ---------------------------------------------------------------------------

def get_gene_window_metadata(dataset_dir: Path, gene_order: list, window_center_size: int) -> Optional[Dict]:
    """Obtém metadados das janelas genômicas ajustados para window_center_size."""
    dataset_dir = Path(dataset_dir)
    ref_dataset_meta = dataset_dir / "dataset_metadata.json"
    if ref_dataset_meta.exists():
        try:
            with open(ref_dataset_meta) as f:
                dataset_meta = json.load(f)
            if dataset_meta.get("window_catalog"):
                catalog = dataset_meta.get("window_catalog", {})
                result = {}
                for gene in gene_order:
                    gm = catalog.get(gene)
                    if not gm:
                        continue
                    orig_start = gm.get("start", 0)
                    orig_size = gm.get("window_size", 0)
                    chrom = gm.get("chromosome", "")
                    center = orig_start + (orig_size + 1) // 2
                    half = window_center_size // 2
                    new_start = center - half
                    new_end = center + half - 1
                    result[gene] = {"chromosome": chrom, "start": new_start, "end": new_end, "window_size": new_end - new_start}
                if result:
                    return result
        except Exception:
            pass

    individuals_dir = dataset_dir / "individuals"
    if not individuals_dir.exists():
        return None

    individual_dirs = [d for d in individuals_dir.iterdir() if d.is_dir()]
    if not individual_dirs:
        return None

    meta_file = individual_dirs[0] / "individual_metadata.json"
    if not meta_file.exists():
        return None

    try:
        with open(meta_file) as f:
            ind_meta = json.load(f)
        window_metadata = ind_meta.get("window_metadata", {})
        result = {}
        for gene in gene_order:
            if gene not in window_metadata:
                continue
            gm = window_metadata[gene]
            orig_start = gm.get("start", 0)
            orig_size = gm.get("window_size", 0)
            chrom = gm.get("chromosome", "")
            center = orig_start + (orig_size + 1) // 2
            half = window_center_size // 2
            new_start = center - half
            new_end = center + half - 1
            result[gene] = {"chromosome": chrom, "start": new_start, "end": new_end, "window_size": new_end - new_start}
        return result
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Save / load cache
# ---------------------------------------------------------------------------

def save_processed_dataset(cache_dir: Path, processed_dataset: ProcessedGenomicDataset,
                           train_indices, val_indices, test_indices, config: PipelineConfig):
    """Salva dataset processado em cache com organização estratificada."""
    import random as _random
    cache_dir = Path(cache_dir)
    random_seed = config.data_split.random_seed
    temp_dir = cache_dir.parent / f"{cache_dir.name}_tmp_{os.getpid()}"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)

    console.print(f"\n[bold cyan]💾 Salvando cache: {cache_dir.name}[/bold cyan]")

    dataset_dir = Path(config.dataset_input.dataset_dir)
    with open(dataset_dir / "dataset_metadata.json") as f:
        src_meta = json.load(f)

    individuals = src_meta.get("individuals", [])
    pedigree = src_meta.get("individuals_pedigree", {})

    # Desativar taint_at_runtime durante o salvamento
    original_taint = processed_dataset.config.debug.taint_at_runtime
    processed_dataset.config.debug.taint_at_runtime = False

    try:
        class_names = {str(i): t for i, t in processed_dataset.idx_to_target.items()}
        gene_order = None
        try:
            sample_data, _ = processed_dataset.base_dataset[0]
            gene_order = list(sample_data["windows"].keys())
        except Exception:
            pass

        gene_window_metadata = None
        if gene_order:
            gene_window_metadata = get_gene_window_metadata(
                dataset_dir, gene_order, config.dataset_input.window_center_size)

        cache_tensors = config.dataset_input.cache_processed_tensors
        splits_meta = {"format_version": 2}
        split_index = {"format_version": 3}
        total_samples = 0

        for split_name, indices, filename in [
            ("train", train_indices, "train_data.pt"),
            ("val", val_indices, "val_data.pt"),
            ("test", test_indices, "test_data.pt"),
        ]:
            data = [] if cache_tensors else None
            split_sample_meta = []

            for idx in indices:
                try:
                    base_idx = processed_dataset.valid_sample_indices[idx]
                    sid = individuals[base_idx] if base_idx < len(individuals) else f"sample_{base_idx}"
                    ped = pedigree.get(sid, {})
                    split_sample_meta.append({
                        "sample_id": sid,
                        "superpopulation": ped.get("superpopulation", "UNK"),
                        "population": ped.get("population", "UNK"),
                        "sex": ped.get("sex", 0),
                    })
                    if cache_tensors:
                        features, target = processed_dataset[idx]
                        data.append((features, target))
                except Exception as e:
                    console.print(f"[yellow]Erro ao processar idx={idx}: {e}[/yellow]")

            if cache_tensors:
                torch.save(data, temp_dir / filename)
            splits_meta[split_name] = split_sample_meta
            split_index[split_name] = [item["sample_id"] for item in split_sample_meta]
            total_samples += len(split_sample_meta)
            if cache_tensors:
                console.print(f"  ✓ {split_name}: {len(split_sample_meta)} amostras salvas")
            else:
                console.print(f"  ✓ {split_name}: {len(split_sample_meta)} amostras indexadas (sem .pt)")

        if config.dataset_input.keep_split_metadata:
            with open(temp_dir / "splits_metadata.json", "w") as f:
                json.dump(splits_meta, f, indent=2)
        with open(temp_dir / "split_index.json", "w") as f:
            json.dump(split_index, f, indent=2)

        # Calcular normalization_params
        norm_params = processed_dataset.normalization_params
        with open(temp_dir / "normalization_params.json", "w") as f:
            json.dump(norm_params, f, indent=2)
        with open(temp_dir / "view_definition.json", "w") as f:
            json.dump({
                "requested_view": _build_view_definition(config),
                "resolved_view": _build_resolved_view_definition(config, processed_dataset),
            }, f, indent=2)

        # Salvar metadata.json
        meta = {
            "creation_date": datetime.now().isoformat(),
            "format_version": 2,
            "total_samples": total_samples,
            "num_classes": processed_dataset.get_num_classes(),
            "class_names": class_names,
            "dataset_dir": str(dataset_dir.resolve()),
            "gene_order": gene_order,
            "tracks_per_gene": 6,
            "gene_window_metadata": gene_window_metadata,
            "processing_params": {
                "alphagenome_outputs": config.dataset_input.alphagenome_outputs,
                "haplotype_mode": config.dataset_input.haplotype_mode,
                "window_center_size": config.dataset_input.window_center_size,
                "downsample_factor": config.dataset_input.downsample_factor,
                "normalization_method": config.dataset_input.normalization_method,
                "selected_track_index": config.dataset_input.selected_track_index,
                "indel_include_valid_mask": config.dataset_input.indel_include_valid_mask,
                "indel_neutral_value": config.dataset_input.indel_neutral_value,
                "tensor_layout": config.dataset_input.tensor_layout,
                "cache_processed_tensors": config.dataset_input.cache_processed_tensors,
                "runtime_dataset_dir": str(dataset_dir.resolve()),
                "family_split_mode": config.data_split.family_split_mode,
                "balancing_strategy": config.data_split.balancing_strategy,
                "dataset_dir": str(dataset_dir.resolve()),
                "prediction_target": config.output.prediction_target,
                "input_shape": "3D_haplotype_channels",
            },
            "splits": {
                "train_size": len(splits_meta.get("train", [])),
                "val_size": len(splits_meta.get("val", [])),
                "test_size": len(splits_meta.get("test", [])),
                "random_seed": config.data_split.random_seed,
            },
        }
        with open(temp_dir / "metadata.json", "w") as f:
            json.dump(meta, f, indent=2)

        (temp_dir / ".cache_complete").touch()

        if cache_dir.exists():
            shutil.rmtree(cache_dir)
        temp_dir.rename(cache_dir)
        console.print(f"[bold green]✓ Cache salvo em {cache_dir}[/bold green]")

    finally:
        processed_dataset.config.debug.taint_at_runtime = original_taint


def _collate_fn(batch):
    if len(batch[0]) == 2:
        features_list, targets_list = zip(*batch)
        indices_list = tuple(range(len(batch)))
    else:
        features_list, targets_list, indices_list = zip(*batch)
    return (
        torch.stack(features_list),
        torch.stack(targets_list),
        torch.tensor(indices_list, dtype=torch.long),
    )


def _make_data_loaders(train_ds, val_ds, test_ds, config: PipelineConfig):
    """Cria DataLoaders para train/val/test."""
    batch_size = config.training.batch_size
    if config.debug.enable_visualization:
        batch_size = 1

    from genotype_based_predictor.models import SKLEARN_BASELINE_TYPES
    model_type = config.model.type.upper()
    use_sklearn = model_type in SKLEARN_BASELINE_TYPES
    lazy = config.data_loading.loading_strategy.lower() == "lazy"
    pin_mem = not use_sklearn

    if lazy or use_sklearn:
        train_w = val_w = 0
        persistent = False
    else:
        train_w, val_w = 4, 2
        persistent = True

    seed = config.data_split.random_seed
    gen = None
    if seed is not None:
        gen = torch.Generator()
        gen.manual_seed(seed)

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              num_workers=train_w, pin_memory=pin_mem, collate_fn=_collate_fn,
                              persistent_workers=persistent and train_w > 0,
                              generator=gen,
                              worker_init_fn=worker_init_fn if train_w > 0 else None)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False,
                            num_workers=val_w, pin_memory=pin_mem, collate_fn=_collate_fn,
                            persistent_workers=persistent and val_w > 0,
                            worker_init_fn=worker_init_fn if val_w > 0 else None)
    test_loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False,
                             num_workers=val_w, pin_memory=pin_mem, collate_fn=_collate_fn,
                             persistent_workers=persistent and val_w > 0,
                             worker_init_fn=worker_init_fn if val_w > 0 else None)
    return train_loader, val_loader, test_loader


def _load_normalization_params(cache_dir: Path) -> Dict:
    with open(cache_dir / "normalization_params.json") as f:
        return json.load(f)


def _load_split_index(cache_dir: Path) -> Dict[str, Any]:
    with open(cache_dir / "split_index.json") as f:
        return json.load(f)


def _make_runtime_processed_datasets(runtime_dataset_dir: Path, cache_dir: Path, config: PipelineConfig):
    base_dataset = _load_base_dataset(runtime_dataset_dir)
    norm_params = _load_normalization_params(cache_dir)
    processed_ds = ProcessedGenomicDataset(
        base_dataset,
        config,
        normalization_params=norm_params,
        compute_normalization=False,
    )

    split_index = _load_split_index(cache_dir)
    dataset_metadata = getattr(base_dataset, "dataset_metadata", {}) or {}
    individuals = dataset_metadata.get("individuals", [])
    sample_to_processed_idx = {}
    for proc_idx, base_idx in enumerate(processed_ds.valid_sample_indices):
        if base_idx < len(individuals):
            sample_to_processed_idx[individuals[base_idx]] = proc_idx

    train_indices = [sample_to_processed_idx[sid] for sid in split_index.get("train", []) if sid in sample_to_processed_idx]
    val_indices = [sample_to_processed_idx[sid] for sid in split_index.get("val", []) if sid in sample_to_processed_idx]
    test_indices = [sample_to_processed_idx[sid] for sid in split_index.get("test", []) if sid in sample_to_processed_idx]

    train_ds = Subset(processed_ds, train_indices)
    val_ds = Subset(processed_ds, val_indices)
    test_ds = Subset(processed_ds, test_indices)
    return processed_ds, train_ds, val_ds, test_ds


def load_processed_dataset(cache_dir: Path, config: PipelineConfig):
    """Carrega dataset processado do cache e cria DataLoaders."""
    cache_dir = Path(cache_dir)
    with open(cache_dir / "metadata.json") as f:
        meta = json.load(f)

    console.print(f"[green]Cache: {meta['creation_date']} | {meta['total_samples']} amostras[/green]")

    class_names = meta.get("class_names", {})
    idx_to_target = {int(k): v for k, v in class_names.items()}
    target_to_idx = {v: int(k) for k, v in idx_to_target.items()}

    if not meta.get("processing_params", {}).get("cache_processed_tensors", True):
        runtime_dataset_dir = Path(
            meta.get("processing_params", {}).get("runtime_dataset_dir")
            or config.dataset_input.dataset_dir
        )
        full_ds, train_ds, val_ds, test_ds = _make_runtime_processed_datasets(runtime_dataset_dir, cache_dir, config)
        train_l, val_l, test_l = _make_data_loaders(train_ds, val_ds, test_ds, config)
        console.print(f"[green]✓ Train:{len(train_ds)} Val:{len(val_ds)} Test:{len(test_ds)} (runtime)[/green]")
        return full_ds, train_l, val_l, test_l, meta

    split_sizes = meta.get("splits", {})
    kwargs = dict(target_to_idx=target_to_idx, idx_to_target=idx_to_target, config=config)
    train_ds = CachedProcessedDataset(cache_dir / "train_data.pt", split_name="train",
                                      length_hint=split_sizes.get("train_size"), **kwargs)
    val_ds = CachedProcessedDataset(cache_dir / "val_data.pt", split_name="val",
                                    length_hint=split_sizes.get("val_size"), **kwargs)
    test_ds = CachedProcessedDataset(cache_dir / "test_data.pt", split_name="test",
                                     length_hint=split_sizes.get("test_size"), **kwargs)

    train_l, val_l, test_l = _make_data_loaders(train_ds, val_ds, test_ds, config)
    console.print(f"[green]✓ Train:{len(train_ds)} Val:{len(val_ds)} Test:{len(test_ds)}[/green]")
    return train_ds, train_l, val_l, test_l, meta


# ---------------------------------------------------------------------------
# prepare_data
# ---------------------------------------------------------------------------

def prepare_data(config: PipelineConfig, experiment_dir: Path):
    """
    Ponto de entrada principal para preparação de dados.

    Tenta carregar do cache. Se não existir ou for inválido, processa
    do zero e salva o cache para reutilização futura.

    Normalização é computada APENAS com train+val para evitar data leakage.
    """
    runtime_dataset_dir = _resolve_runtime_dataset_dir(config)

    cache_dir = get_dataset_cache_dir(config)
    cache_path = Path(cache_dir)

    if cache_path.exists() and validate_cache(cache_path, config):
        console.print(f"[bold green]Cache encontrado: {cache_path.name}[/bold green]")
        full_ds, train_l, val_l, test_l, cache_meta = load_processed_dataset(cache_path, config)
        if cache_meta.get("gene_order"):
            config.dataset_input.gene_order = cache_meta["gene_order"]
        if cache_meta.get("gene_window_metadata"):
            config.dataset_input.gene_window_metadata = cache_meta["gene_window_metadata"]
        (experiment_dir / "models").mkdir(exist_ok=True)
        norm_src = cache_path / "normalization_params.json"
        norm_dst = experiment_dir / "models" / "normalization_params.json"
        if norm_src.exists() and not norm_dst.exists():
            shutil.copyfile(norm_src, norm_dst)
        return full_ds, train_l, val_l, test_l

    # Processar do zero
    console.print("[bold cyan]Processando dataset do zero...[/bold cyan]")
    base_dataset = _load_base_dataset(runtime_dataset_dir)

    valid_base_indices = build_valid_sample_index_map(base_dataset, config)
    valid_set = set(valid_base_indices)
    proc_idx_by_base = {bi: pi for pi, bi in enumerate(valid_base_indices)}

    all_groups, fam_info = build_family_aware_sample_groups(base_dataset, config)
    sample_groups = [[bi for bi in g if bi in valid_set] for g in all_groups]
    sample_groups = [g for g in sample_groups if g]

    n_groups = len(sample_groups)
    n_train = int(config.data_split.train_split * n_groups)
    n_val = int(config.data_split.val_split * n_groups)
    indices = list(range(n_groups))

    seed = config.data_split.random_seed
    if seed is not None and seed != -1:
        np.random.seed(seed)
        np.random.shuffle(indices)

    train_g = indices[:n_train]
    val_g = indices[n_train:n_train + n_val]
    test_g = indices[n_train + n_val:]

    def flatten(group_idxs):
        return [proc_idx_by_base[bi] for gi in group_idxs for bi in sample_groups[gi]]

    train_idx, val_idx, test_idx = flatten(train_g), flatten(val_g), flatten(test_g)

    # Normalização apenas com train+val, a menos que um arquivo compatível já exista.
    norm_params = _load_external_normalization_params(config)
    if norm_params is None:
        train_val_base = [bi for gi in (train_g + val_g) for bi in sample_groups[gi]]
        tv_subset = Subset(base_dataset, train_val_base)
        temp_ds = ProcessedGenomicDataset(tv_subset, config, normalization_params=None, compute_normalization=True)
        norm_params = temp_ds.normalization_params

    processed_ds = ProcessedGenomicDataset(base_dataset, config, normalization_params=norm_params, compute_normalization=False)

    norm_path = experiment_dir / "models" / "normalization_params.json"
    norm_path.parent.mkdir(parents=True, exist_ok=True)
    with open(norm_path, "w") as f:
        json.dump(norm_params, f, indent=2)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    save_processed_dataset(cache_path, processed_ds, train_idx, val_idx, test_idx, config)

    full_ds, train_l, val_l, test_l, cache_meta = load_processed_dataset(cache_path, config)
    if cache_meta.get("gene_order"):
        config.dataset_input.gene_order = cache_meta["gene_order"]
    if cache_meta.get("gene_window_metadata"):
        config.dataset_input.gene_window_metadata = cache_meta["gene_window_metadata"]

    shutil.copyfile(cache_path / "normalization_params.json", norm_path)
    return full_ds, train_l, val_l, test_l
