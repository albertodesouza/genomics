from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import torch
from matplotlib.colors import LinearSegmentedColormap
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from genotype_based_predictor.config import generate_experiment_name, load_config
from genotype_based_predictor.data.pipeline import prepare_data
from genotype_based_predictor.analysis.interpretability import (
    DeepLIFT,
    build_haplotype_track_layout,
    find_top_deeplift_windows,
    summarize_deeplift_by_track,
)
from genotype_based_predictor.alignment.dynamic_indel_alignment import (
    DynamicIndelAligner,
    _allele_sequence,
    _classify_length_behavior,
    _parse_gt,
    _resolve_ref_index,
    _stream_gene_variants,
    _stream_gene_variants_from_gzip,
)
from genotype_based_predictor.models import CNN2AncestryPredictor, CNNAncestryPredictor, NNAncestryPredictor
from genotype_based_predictor.utils import set_random_seeds

console = Console()


def _build_model(config: Any, dataset: Any) -> torch.nn.Module:
    input_shape = dataset.get_input_shape()
    num_classes = dataset.get_num_classes()
    model_type = config.model.type.upper()
    if model_type == "NN":
        return NNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN":
        return CNNAncestryPredictor(config, input_shape, num_classes)
    if model_type == "CNN2":
        return CNN2AncestryPredictor(config, input_shape, num_classes)
    raise ValueError(f"DeepLIFT suporta apenas modelos PyTorch NN/CNN/CNN2, recebeu {config.model.type}")


def _select_dataset(split: str, train_loader: Any, val_loader: Any, test_loader: Any) -> Any:
    split = split.lower()
    if split == "train":
        return train_loader.dataset
    if split == "val":
        return val_loader.dataset
    if split == "test":
        return test_loader.dataset
    raise ValueError("split deve ser train, val ou test")


def _unpack_item(item: Any):
    if isinstance(item, (tuple, list)) and len(item) >= 3:
        return item[0], item[1], item[2]
    if isinstance(item, (tuple, list)) and len(item) == 2:
        return item[0], item[1], None
    raise ValueError("Dataset item deve ser (features, target) ou (features, target, idx)")


def _target_to_int(target: Any) -> Optional[int]:
    if isinstance(target, torch.Tensor):
        if target.numel() != 1:
            return None
        return int(target.item())
    try:
        return int(target)
    except Exception:
        return None


def _sample_id_for(dataset: Any, idx: int, meta_idx: Any) -> str:
    candidate = meta_idx
    if isinstance(candidate, torch.Tensor) and candidate.numel() == 1:
        candidate = int(candidate.item())
    if isinstance(candidate, int) and hasattr(dataset, "get_sample_id"):
        return str(dataset.get_sample_id(candidate))
    if hasattr(dataset, "get_sample_id"):
        return str(dataset.get_sample_id(idx))
    return f"#{idx + 1}"


def _class_arg_to_idx(class_arg: str, class_name_to_idx: Dict[str, int]) -> int:
    if class_arg in class_name_to_idx:
        return class_name_to_idx[class_arg]
    return int(class_arg)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def _write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with open(path, "w") as f:
            f.write("")
        return
    fieldnames: List[str] = []
    seen = set()
    for row in rows:
        for key in row.keys():
            if key not in seen and key != "metadata":
                fieldnames.append(key)
                seen.add(key)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: v for k, v in row.items() if k in seen})


def _aggregate_metric(rows: Iterable[Dict[str, Any]], key: str) -> List[Dict[str, Any]]:
    acc: Dict[tuple, Dict[str, Any]] = {}
    for row in rows:
        group_key = tuple(row[k] for k in key.split("|"))
        if group_key not in acc:
            entry = {k: row[k] for k in key.split("|")}
            entry.update({"sum": 0.0, "count": 0})
            acc[group_key] = entry
        entry = acc[group_key]
        entry["sum"] += float(row["mean_abs_attr"])
        entry["count"] += 1
    out = []
    for entry in acc.values():
        count = max(int(entry.pop("count")), 1)
        total = float(entry.pop("sum"))
        entry["mean_abs_attr"] = total / count
        entry["n_rows"] = count
        out.append(entry)
    out.sort(key=lambda r: r["mean_abs_attr"], reverse=True)
    return out


def _plot_layout(config: Any) -> Dict[str, Any]:
    genes = list(config.dataset_input.genes_to_use or config.dataset_input.gene_order or [])
    signal_track_count = 2 * len(config.dataset_input.ontology_terms or []) if config.dataset_input.ontology_terms else 1
    mask_track_count = 2
    if config.dataset_input.indel_include_valid_mask:
        mask_track_count += 1
    if config.dataset_input.indel_include_snp_mask:
        mask_track_count += 1
    feature_mode = getattr(config.dataset_input, "feature_mode", "signals_and_masks")
    if feature_mode == "signals_only":
        mask_track_count = 0
    elif feature_mode == "masks_only":
        signal_track_count = 0
    tracks_per_gene = signal_track_count + mask_track_count
    return {
        "genes": genes,
        "tracks_per_gene": tracks_per_gene,
        "rows_per_hap": len(genes) * tracks_per_gene,
        "signal_track_count": signal_track_count,
        "mask_track_count": mask_track_count,
    }


def _flatten_haps_for_plot(x: torch.Tensor, config: Any) -> np.ndarray:
    layout = _plot_layout(config)
    genes = layout["genes"]
    tracks_per_gene = layout["tracks_per_gene"]
    rows_per_hap = layout["rows_per_hap"]
    arr = x.detach().cpu()
    if arr.ndim == 3:
        if genes and arr.shape[0] == 2 and arr.shape[1] == rows_per_hap:
            blocks = []
            for gene_idx in range(len(genes)):
                start = gene_idx * tracks_per_gene
                end = start + tracks_per_gene
                blocks.append(arr[0, start:end, :])
                blocks.append(arr[1, start:end, :])
            arr = torch.cat(blocks, dim=0)
        else:
            arr = arr.reshape(arr.shape[0] * arr.shape[1], arr.shape[2])
    elif arr.ndim != 2:
        raise ValueError(f"Esperado tensor 2D/3D para plot, recebeu {tuple(arr.shape)}")
    return arr.numpy()


def _split_signal_and_masks(matrix: np.ndarray, config: Any) -> Tuple[np.ndarray, np.ndarray]:
    layout = _plot_layout(config)
    tracks_per_gene = layout["tracks_per_gene"]
    signal_track_count = layout["signal_track_count"]
    mask_track_count = layout["mask_track_count"]
    signal = matrix.copy()
    masks = np.full_like(matrix, np.nan, dtype=np.float32)
    if mask_track_count > 0 and tracks_per_gene > 0:
        for row_idx in range(matrix.shape[0]):
            if row_idx % tracks_per_gene >= signal_track_count:
                signal[row_idx, :] = np.nan
                masks[row_idx, :] = matrix[row_idx, :]
    return signal, masks


def _save_raw_pixel_images(mean_input: torch.Tensor, mean_attr: torch.Tensor, config: Any, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    input_np = _flatten_haps_for_plot(mean_input, config)
    attr_np = _flatten_haps_for_plot(mean_attr, config)
    signal_input, mask_input = _split_signal_and_masks(input_np, config)
    signal_attr, mask_attr = _split_signal_and_masks(attr_np, config)
    has_masks = _plot_layout(config)["mask_track_count"] > 0

    attr_lim = float(np.nanmax(np.abs(signal_attr))) if np.isfinite(signal_attr).any() else 1.0
    mask_attr_lim = float(np.nanmax(np.abs(mask_attr))) if np.isfinite(mask_attr).any() else 1.0
    if attr_lim <= 0:
        attr_lim = 1.0
    if mask_attr_lim <= 0:
        mask_attr_lim = 1.0
    mask_cmap = LinearSegmentedColormap.from_list("mask_red_to_blue", ["#d73027", "#f7f7f7", "#2166ac"])

    plt.imsave(out_dir / "raw_input_signal_gray.png", signal_input, cmap="gray", vmin=0)
    plt.imsave(out_dir / "raw_deeplift_signal_bwr.png", signal_attr, cmap="bwr", vmin=-attr_lim, vmax=attr_lim)
    if has_masks:
        plt.imsave(out_dir / "raw_input_masks_0_1.png", np.clip(mask_input, 0.0, 1.0), cmap=mask_cmap, vmin=0, vmax=1)
        plt.imsave(out_dir / "raw_deeplift_masks_bwr.png", mask_attr, cmap="bwr", vmin=-mask_attr_lim, vmax=mask_attr_lim)

    np.save(out_dir / "raw_input_matrix.npy", input_np)
    np.save(out_dir / "raw_deeplift_matrix.npy", attr_np)


def _to_scalar_int(value: Any) -> Optional[int]:
    if isinstance(value, torch.Tensor):
        if value.numel() != 1:
            return None
        return int(value.item())
    if isinstance(value, (int, np.integer)):
        return int(value)
    return None


def _metadata_for_dataset_item(dataset: Any, idx: int, meta_idx: Any) -> Dict[str, Any]:
    item_idx = _to_scalar_int(meta_idx)
    if item_idx is None:
        item_idx = idx

    metadata: Dict[str, Any] = {}
    if hasattr(dataset, "get_sample_metadata"):
        try:
            metadata.update(dict(dataset.get_sample_metadata(item_idx) or {}))
        except Exception:
            pass

    if hasattr(dataset, "indices") and hasattr(dataset, "dataset"):
        try:
            inner_idx = int(dataset.indices[idx])
            nested = _metadata_for_dataset_item(dataset.dataset, inner_idx, inner_idx)
            nested.update(metadata)
            metadata = nested
        except Exception:
            pass

    if not metadata.get("superpopulation") and hasattr(dataset, "valid_sample_indices") and hasattr(dataset, "_sample_id_for_base_index"):
        try:
            base_idx = int(dataset.valid_sample_indices[item_idx])
            sample_id = dataset._sample_id_for_base_index(base_idx)
            if sample_id:
                metadata["sample_id"] = sample_id
                pedigree = getattr(dataset, "dataset_metadata", {}).get("individuals_pedigree", {}) or {}
                metadata.update(pedigree.get(sample_id, {}) or {})
        except Exception:
            pass

    if "sample_id" not in metadata:
        metadata["sample_id"] = _sample_id_for(dataset, idx, meta_idx)
    if not metadata.get("superpopulation"):
        sample_id = metadata.get("sample_id")
        try:
            split_metadata_cache = getattr(dataset, "_opencode_split_metadata_by_sample", None)
            if split_metadata_cache is None:
                split_metadata_cache = {}
                data_file = getattr(dataset, "data_file", None)
                split_name = getattr(dataset, "split_name", None)
                if data_file is not None and split_name is not None:
                    splits_file = Path(data_file).parent / "splits_metadata.json"
                    if splits_file.exists():
                        with open(splits_file) as f:
                            splits = json.load(f)
                        for item in splits.get(str(split_name), []) or []:
                            sid = item.get("sample_id")
                            if sid:
                                split_metadata_cache[str(sid)] = item
                setattr(dataset, "_opencode_split_metadata_by_sample", split_metadata_cache)
            if sample_id and str(sample_id) in split_metadata_cache:
                metadata.update(split_metadata_cache[str(sample_id)] or {})
        except Exception:
            pass
    if not metadata.get("superpopulation"):
        sample_id = metadata.get("sample_id")
        try:
            pedigree_cache = getattr(dataset, "_opencode_pedigree_by_sample", None)
            if pedigree_cache is None:
                pedigree_cache = {}
                config = getattr(dataset, "config", None)
                dataset_dir = Path(config.dataset_input.dataset_dir) if config is not None else None
                if dataset_dir is not None:
                    metadata_file = dataset_dir / "dataset_metadata.json"
                    if metadata_file.exists():
                        with open(metadata_file) as f:
                            payload = json.load(f)
                        pedigree_cache = payload.get("individuals_pedigree", {}) or {}
                setattr(dataset, "_opencode_pedigree_by_sample", pedigree_cache)
            if sample_id and str(sample_id) in pedigree_cache:
                metadata.update(pedigree_cache[str(sample_id)] or {})
                metadata.setdefault("sample_id", sample_id)
        except Exception:
            pass
    return metadata


def _ordered_superpopulations(superpops: Iterable[str]) -> List[str]:
    preferred = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    seen = {str(s) for s in superpops if str(s)}
    ordered = [s for s in preferred if s in seen]
    ordered.extend(sorted(seen - set(ordered)))
    return ordered


def _ordered_class_groups(groups: Iterable[str], class_names: List[str]) -> List[str]:
    seen = {str(group) for group in groups if str(group)}
    ordered = [name for name in class_names if name in seen]
    ordered.extend(sorted(seen - set(ordered)))
    return ordered


def _class_name_for_target(target: Any, class_names: List[str]) -> str:
    target_idx = _target_to_int(target)
    if target_idx is not None and 0 <= target_idx < len(class_names):
        return class_names[target_idx]
    return "UNK" if target_idx is None else str(target_idx)


def _collect_superpopulation_mean_inputs(dataset: Any, split: str) -> Tuple[Dict[str, torch.Tensor], Dict[str, int], List[Dict[str, Any]]]:
    sums: Dict[str, torch.Tensor] = {}
    counts: Dict[str, int] = {}
    samples: List[Dict[str, Any]] = []

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(f"Calculando medias RNA-Seq por superpopulacao ({split})...", total=len(dataset))
        for idx in range(len(dataset)):
            features, target, meta_idx = _unpack_item(dataset[idx])
            metadata = _metadata_for_dataset_item(dataset, idx, meta_idx)
            superpop = str(metadata.get("superpopulation") or metadata.get("super_population") or "UNK")
            tensor = features.detach().cpu().to(torch.float64)
            if superpop not in sums:
                sums[superpop] = tensor.clone()
                counts[superpop] = 1
            else:
                sums[superpop] += tensor
                counts[superpop] += 1
            samples.append({
                "dataset_index": idx,
                "sample_id": metadata.get("sample_id", _sample_id_for(dataset, idx, meta_idx)),
                group_label: superpop,
                "population": metadata.get("population", "UNK"),
                "target": _target_to_int(target),
            })
            progress.update(task, advance=1)

    means = {superpop: (sums[superpop] / max(counts[superpop], 1)).to(torch.float32) for superpop in sums}
    return means, counts, samples


def _collect_class_mean_inputs(dataset: Any, split: str, class_names: List[str]) -> Tuple[Dict[str, torch.Tensor], Dict[str, int], List[Dict[str, Any]]]:
    sums: Dict[str, torch.Tensor] = {}
    counts: Dict[str, int] = {}
    samples: List[Dict[str, Any]] = []

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(f"Calculando medias RNA-Seq por classe ({split})...", total=len(dataset))
        for idx in range(len(dataset)):
            features, target, meta_idx = _unpack_item(dataset[idx])
            metadata = _metadata_for_dataset_item(dataset, idx, meta_idx)
            class_name = _class_name_for_target(target, class_names)
            tensor = features.detach().cpu().to(torch.float64)
            if class_name not in sums:
                sums[class_name] = tensor.clone()
                counts[class_name] = 1
            else:
                sums[class_name] += tensor
                counts[class_name] += 1
            samples.append({
                "dataset_index": idx,
                "sample_id": metadata.get("sample_id", _sample_id_for(dataset, idx, meta_idx)),
                "class_name": class_name,
                "superpopulation": metadata.get("superpopulation", metadata.get("super_population", "UNK")),
                "population": metadata.get("population", "UNK"),
                "target": _target_to_int(target),
            })
            progress.update(task, advance=1)

    means = {class_name: (sums[class_name] / max(counts[class_name], 1)).to(torch.float32) for class_name in sums}
    return means, counts, samples


def _signal_descriptors_for_tensor(config: Any, tensor: torch.Tensor) -> List[Dict[str, Any]]:
    arr = tensor.detach().cpu()
    if arr.ndim == 2:
        arr = arr.unsqueeze(0)
    if arr.ndim != 3:
        raise ValueError(f"Esperado tensor [hap,row,L] ou [row,L], recebido {tuple(tensor.shape)}")

    genes = list(config.dataset_input.genes_to_use or config.dataset_input.gene_order or [])
    if not genes:
        raise ValueError("genes_to_use/gene_order vazio; nao e possivel nomear as linhas de RNA-Seq")
    tracks_per_gene = arr.shape[1] // len(genes)
    if tracks_per_gene <= 0:
        raise ValueError(f"Tensor com {arr.shape[1]} linhas nao comporta {len(genes)} genes")

    ontology_terms = list(config.dataset_input.ontology_terms or [])
    if ontology_terms:
        signal_tracks = []
        for strand in ("+", "-"):
            for ontology in ontology_terms:
                signal_tracks.append({"ontology_curie": ontology, "strand": strand})
    else:
        signal_tracks = [{"ontology_curie": "signal", "strand": None}]
    signal_tracks = signal_tracks[:tracks_per_gene]

    haplotypes = ["H1", "H2"] if arr.shape[0] == 2 else [f"H{i + 1}" for i in range(arr.shape[0])]
    descriptors: List[Dict[str, Any]] = []
    for gene_idx, gene_name in enumerate(genes):
        for hap_idx, haplotype in enumerate(haplotypes):
            for signal_idx, signal_meta in enumerate(signal_tracks):
                descriptors.append({
                    "gene_index": gene_idx,
                    "gene_name": gene_name,
                    "haplotype_index": hap_idx,
                    "haplotype": haplotype,
                    "signal_index": signal_idx,
                    "row_index": gene_idx * tracks_per_gene + signal_idx,
                    "ontology_curie": signal_meta["ontology_curie"],
                    "strand": signal_meta["strand"],
                })
    return descriptors


def _safe_filename(value: str, max_len: int = 180) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.+-]+", "_", value).strip("_")
    if not safe:
        safe = "track"
    return safe[:max_len]


def _save_superpopulation_rnaseq_visualizations(
    means_by_superpop: Dict[str, torch.Tensor],
    counts_by_superpop: Dict[str, int],
    config: Any,
    split: str,
    out_dir: Path,
    delta_reference: bool = False,
    top_abs_markers: int = 20,
    group_label: str = "superpopulation",
    group_order: Optional[List[str]] = None,
) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from PIL import Image, ImageDraw

    if not means_by_superpop:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    line_dir = out_dir / "line_plots"
    line_dir.mkdir(parents=True, exist_ok=True)

    superpops = group_order or _ordered_superpopulations(means_by_superpop.keys())
    reference_tensor = means_by_superpop[superpops[0]]
    descriptors = _signal_descriptors_for_tensor(config, reference_tensor)
    length = int(reference_tensor.shape[-1])

    rows = []
    row_metadata = []
    for descriptor_idx, descriptor in enumerate(descriptors):
        hap_idx = int(descriptor["haplotype_index"])
        row_idx = int(descriptor["row_index"])
        for superpop in superpops:
            tensor = means_by_superpop[superpop]
            arr = tensor.detach().cpu()
            if arr.ndim == 2:
                arr = arr.unsqueeze(0)
            if hap_idx >= arr.shape[0] or row_idx >= arr.shape[1]:
                values = np.full(length, np.nan, dtype=np.float32)
            else:
                values = arr[hap_idx, row_idx, :].numpy().astype(np.float32, copy=False)
            rows.append(values)
            row_metadata.append({
                "matrix_row": len(row_metadata),
                "descriptor_index": descriptor_idx,
                "superpopulation": superpop,
                "num_samples": counts_by_superpop.get(superpop, 0),
                **{k: v for k, v in descriptor.items() if k != "haplotype_index"},
            })

    matrix = np.stack(rows, axis=0) if rows else np.zeros((0, length), dtype=np.float32)
    finite = matrix[np.isfinite(matrix)]
    vmin = float(finite.min()) if finite.size else 0.0
    vmax = float(finite.max()) if finite.size else 1.0
    if vmin == vmax:
        vmax = vmin + 1.0

    plt.imsave(
        out_dir / "superpopulation_rnaseq_grayscale.png",
        matrix,
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
        pil_kwargs={"compress_level": 0},
    )
    np.save(out_dir / "superpopulation_rnaseq_grayscale.npy", matrix)
    _write_json(out_dir / "superpopulation_rnaseq_rows.json", row_metadata)
    _write_json(out_dir / "superpopulation_rnaseq_metadata.json", {
        "split": split,
        f"{group_label}_order": superpops,
        f"counts_by_{group_label}": counts_by_superpop,
        "matrix_shape": list(matrix.shape),
        "grayscale_vmin": vmin,
        "grayscale_vmax": vmax,
        "row_grouping": "For each ontology/strand/gene/haplotype descriptor, rows are consecutive superpopulations in superpopulation_order.",
        "image_policy": "No interpolation, no positional downsampling, PNG compress_level=0; raw values are also saved as .npy.",
    })

    dark_cyan = (3, 18, 18)
    separator_height = 2
    if delta_reference:
        abs_max = float(np.nanmax(np.abs(matrix))) if np.isfinite(matrix).any() else 1.0
        if abs_max <= 0:
            abs_max = 1.0
        centered = np.clip((matrix + abs_max) / (2.0 * abs_max), 0.0, 1.0)
        centered = np.nan_to_num(centered, nan=0.5, posinf=1.0, neginf=0.0)
        gray_u8 = (centered * 255.0 + 0.5).astype(np.uint8)

        rendered = []
        rendered_rows = []
        separators = []
        matrix_to_rendered_row: Dict[int, int] = {}
        height_so_far = 0
        for idx, row_meta in enumerate(row_metadata):
            data_rgb = np.repeat(gray_u8[idx:idx + 1, :, None], 3, axis=2)
            rendered.append(data_rgb)
            rendered_rows.append({**row_meta, "rendered_data_row": height_so_far})
            matrix_to_rendered_row[int(row_meta["matrix_row"])] = height_so_far
            height_so_far += 1

            is_last = idx == len(row_metadata) - 1
            next_meta = None if is_last else row_metadata[idx + 1]
            descriptor_changed = is_last or next_meta.get("descriptor_index") != row_meta.get("descriptor_index")
            if not descriptor_changed:
                continue

            rendered.append(np.full((separator_height, gray_u8.shape[1], 3), dark_cyan, dtype=np.uint8))
            separators.append({
                "after_matrix_row": int(row_meta["matrix_row"]),
                "rendered_separator_start_row": height_so_far,
                "height": separator_height,
                "color_rgb": dark_cyan,
                "reason": "ontology_strand_boundary",
                "gene_name": row_meta.get("gene_name"),
                "haplotype": row_meta.get("haplotype"),
                "ontology_curie": row_meta.get("ontology_curie"),
                "strand": row_meta.get("strand"),
            })
            height_so_far += separator_height

        image = np.concatenate(rendered, axis=0)
        pil = Image.fromarray(image, mode="RGB")
        draw = ImageDraw.Draw(pil)
        abs_matrix = np.nan_to_num(np.abs(matrix), nan=-1.0, posinf=-1.0, neginf=-1.0)
        flat = abs_matrix.ravel()
        top_k = min(max(int(top_abs_markers), 0), flat.size)
        if top_k > 0:
            top_indices = np.argpartition(flat, -top_k)[-top_k:]
            top_indices = top_indices[np.argsort(flat[top_indices])[::-1]]
        else:
            top_indices = np.array([], dtype=np.int64)

        markers = []
        green = (0, 255, 0)
        radius = 7
        for rank, flat_idx in enumerate(top_indices.tolist(), start=1):
            row_idx, x = np.unravel_index(flat_idx, matrix.shape)
            y = matrix_to_rendered_row.get(int(row_idx))
            if y is None:
                continue
            draw.ellipse(
                [(int(x) - radius, int(y) - radius), (int(x) + radius, int(y) + radius)],
                outline=green,
                width=2,
            )
            row_meta = row_metadata[int(row_idx)]
            markers.append({
                "rank": rank,
                "matrix_row": int(row_idx),
                "x": int(x),
                "rendered_y": int(y),
                "value": float(matrix[row_idx, x]),
                "abs_value": float(abs_matrix[row_idx, x]),
                "gene_name": row_meta.get("gene_name"),
                "haplotype": row_meta.get("haplotype"),
                "ontology_curie": row_meta.get("ontology_curie"),
                "strand": row_meta.get("strand"),
                group_label: row_meta.get(group_label),
                "num_samples": row_meta.get("num_samples"),
            })

        image_path = out_dir / "superpopulation_rnaseq_delta_reference_dark_cyan_top_abs_deviation.png"
        pil.save(image_path, compress_level=0)
        _write_json(out_dir / "superpopulation_rnaseq_delta_reference_top_abs_deviation_markers.json", markers)
        _write_json(out_dir / "superpopulation_rnaseq_delta_reference_dark_cyan_separators.json", separators)
        _write_json(out_dir / "superpopulation_rnaseq_delta_reference_rows_with_rendered_indices.json", rendered_rows)
        _write_json(out_dir / "superpopulation_rnaseq_delta_reference_dark_cyan_top_abs_deviation_legend.json", {
            "split": split,
            "source_matrix": str(out_dir / "superpopulation_rnaseq_grayscale.npy"),
            "source_rows": str(out_dir / "superpopulation_rnaseq_rows.json"),
            "output_image": str(image_path),
            "image_shape": list(np.asarray(pil).shape),
            "value_encoding": {
                "black": -abs_max,
                "middle_gray": 0.0,
                "white": abs_max,
                "abs_max": abs_max,
            },
            "separator_color_rgb": dark_cyan,
            "separator_height_px": separator_height,
            "marker_color_rgb": green,
            "marker_count": len(markers),
            "top_abs_markers_requested": top_abs_markers,
            "policy": "All data rows are grayscale and centered at zero. No data row is resampled or interpolated. Dark-cyan separator rows are inserted after ontology/strand blocks. Green circles mark the largest absolute deviations. PNG compress_level=0.",
        })

    previous_simplify = mpl.rcParams.get("path.simplify", False)
    mpl.rcParams["path.simplify"] = False
    dpi = 100
    try:
        for descriptor_idx, descriptor in enumerate(descriptors):
            hap_idx = int(descriptor["haplotype_index"])
            row_idx = int(descriptor["row_index"])
            series = []
            for superpop in superpops:
                tensor = means_by_superpop[superpop]
                arr = tensor.detach().cpu()
                if arr.ndim == 2:
                    arr = arr.unsqueeze(0)
                if hap_idx >= arr.shape[0] or row_idx >= arr.shape[1]:
                    values = np.full(length, np.nan, dtype=np.float32)
                else:
                    values = arr[hap_idx, row_idx, :].numpy().astype(np.float32, copy=False)
                series.append((superpop, values))

            finite_values = np.concatenate([values[np.isfinite(values)] for _, values in series if np.isfinite(values).any()])
            ymin = float(finite_values.min()) if finite_values.size else 0.0
            ymax = float(finite_values.max()) if finite_values.size else 1.0
            if ymin == ymax:
                ymin -= 0.5
                ymax += 0.5

            fig_width = max(length / dpi, 8.0)
            fig_height = max(1.8 * len(superpops), 3.0)
            fig, axes = plt.subplots(len(superpops), 1, figsize=(fig_width, fig_height), sharex=True, squeeze=False)
            x = np.arange(length, dtype=np.int32)
            for ax, (superpop, values) in zip(axes[:, 0], series):
                ax.plot(x, values, color="black", linewidth=0.7, antialiased=False)
                ax.set_ylim(ymin, ymax)
                ax.set_ylabel(f"{superpop}\nn={counts_by_superpop.get(superpop, 0)}", rotation=0, ha="right", va="center", fontsize=8)
                ax.grid(False)
            axes[-1, 0].set_xlabel("Aligned gene position")
            title = f"{descriptor['gene_name']} | {descriptor['haplotype']} | {descriptor['ontology_curie']} | strand={descriptor['strand']}"
            fig.suptitle(title, fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            filename_stem = _safe_filename(
                f"{descriptor_idx:05d}_{descriptor['gene_name']}_{descriptor['haplotype']}_{descriptor['ontology_curie']}_{descriptor['strand']}"
            )
            fig.savefig(line_dir / f"{filename_stem}.png", dpi=dpi, pil_kwargs={"compress_level": 0})
            plt.close(fig)
    finally:
        mpl.rcParams["path.simplify"] = previous_simplify


def _variant_row_index(gene_idx: int, haplotype: str, track_idx: int, tracks_per_gene: int) -> int:
    hap_offset = 0 if haplotype == "H1" else tracks_per_gene
    return gene_idx * 2 * tracks_per_gene + hap_offset + track_idx


def _plot_row_index_for_window(win: Dict[str, Any], genes: List[str], rows_per_hap: int, tracks_per_gene: int) -> int:
    row_index = int(win.get("row_index", 0))
    if genes and rows_per_hap > 0:
        gene_idx = int(win.get("gene_index", row_index // tracks_per_gene))
        track_idx = int(win.get("track_index_in_gene", row_index % tracks_per_gene))
        hap_offset = tracks_per_gene if str(win.get("haplotype")) == "H2" else 0
        return gene_idx * 2 * tracks_per_gene + hap_offset + track_idx
    if str(win.get("haplotype")) == "H2":
        return row_index + rows_per_hap
    return row_index


def _variant_marker_label(marker: Dict[str, Any]) -> str:
    variant_id = str(marker.get("id") or marker.get("variant_id") or "").strip()
    if variant_id and variant_id != ".":
        return variant_id
    chrom = str(marker.get("chrom") or "").strip()
    pos = marker.get("pos_1based")
    if chrom and pos is not None:
        return f"{chrom}:{pos}"
    gene = str(marker.get("gene_name") or "")
    x = marker.get("x")
    return f"{gene}:x{x}" if gene else f"x{x}"


def _select_variant_markers_near_windows(
    variant_markers: List[Dict[str, Any]],
    top_windows: List[Dict[str, Any]],
    genes: List[str],
    rows_per_hap: int,
    tracks_per_gene: int,
    max_windows: int = 10,
) -> List[Dict[str, Any]]:
    drawable = [m for m in variant_markers if "x" in m and "row_index" in m]
    selected: List[Dict[str, Any]] = []
    used_keys = set()

    for win in top_windows[:max_windows]:
        win_gene = win.get("gene_name")
        win_haplotype = str(win.get("haplotype", ""))
        win_center = int(win.get("window_center", 0))
        win_row = _plot_row_index_for_window(win, genes, rows_per_hap, tracks_per_gene)
        candidates = [
            marker for marker in drawable
            if marker.get("gene_name") == win_gene and str(marker.get("haplotype", "")) == win_haplotype
        ]
        if not candidates:
            continue

        marker = min(
            candidates,
            key=lambda item: (
                abs(int(item["x"]) - win_center),
                abs(int(item["row_index"]) - win_row),
                -float(item.get("frequency", 0.0) or 0.0),
            ),
        )
        key = (
            marker.get("gene_name"),
            marker.get("haplotype"),
            int(marker["x"]),
            int(marker["row_index"]),
            marker.get("variant_type"),
            marker.get("id"),
            marker.get("pos_1based"),
        )
        if key in used_keys:
            continue
        used_keys.add(key)
        enriched = dict(marker)
        enriched["nearest_top_window_center"] = win_center
        enriched["nearest_top_window_row_index"] = win_row
        enriched["nearest_top_window_distance"] = abs(int(marker["x"]) - win_center)
        selected.append(enriched)

    return selected


def _collect_indel_markers_from_mean_input(
    *,
    mean_input: torch.Tensor,
    config: Any,
    max_markers: int,
    min_frequency: float,
    genes_filter: Optional[set[str]] = None,
) -> List[Dict[str, Any]]:
    matrix = _flatten_haps_for_plot(mean_input, config)
    layout = _plot_layout(config)
    genes = layout["genes"]
    tracks_per_gene = layout["tracks_per_gene"]
    signal_track_count = layout["signal_track_count"]
    markers: List[Dict[str, Any]] = []
    if not genes:
        return markers

    for gene_idx, gene_name in enumerate(genes):
        if genes_filter and gene_name not in genes_filter:
            continue
        for hap_idx, haplotype in enumerate(("H1", "H2")):
            block_offset = gene_idx * 2 * tracks_per_gene + hap_idx * tracks_per_gene
            for track_idx, variant_type in ((signal_track_count, "INDEL_ins"), (signal_track_count + 1, "INDEL_del")):
                row_index = block_offset + track_idx
                if row_index >= matrix.shape[0]:
                    continue
                row = matrix[row_index]
                row_values = np.nan_to_num(row, nan=0.0)
                if min_frequency <= 0:
                    positions = np.where(row_values > 0)[0]
                else:
                    positions = np.where(row_values >= min_frequency)[0]
                for x in positions.tolist():
                    frequency = float(row_values[x])
                    markers.append({
                        "gene_name": gene_name,
                        "gene_index": gene_idx,
                        "haplotype": haplotype,
                        "x": int(x),
                        "variant_type": variant_type,
                        "count": None,
                        "frequency": frequency,
                        "track_idx": track_idx,
                        "row_index": int(row_index),
                        "source": "mean_input_indel_mask",
                    })
    markers.sort(key=lambda item: (item["frequency"], item["gene_name"], item["x"]), reverse=True)
    return markers if max_markers <= 0 else markers[:max_markers]


def _variant_cache_path(
    out_dir: Path,
    sample_ids: List[str],
    variant_types: str,
    min_frequency: float,
    max_markers: int,
    genes_filter: Optional[set[str]],
) -> Path:
    payload = json.dumps({
        "cache_version": 3,
        "sample_ids": sorted(sample_ids),
        "variant_types": variant_types,
        "min_frequency": min_frequency,
        "max_markers": max_markers,
        "genes_filter": sorted(genes_filter or []),
    }, sort_keys=True)
    digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]
    return out_dir / "variant_marker_cache" / f"variant_markers_{digest}.json"


def _stream_variants_with_gzip_fallback(vcf_path: str, sample_ids: List[str], region: str):
    try:
        yield from _stream_gene_variants(vcf_path, sample_ids, region)
    except subprocess.CalledProcessError as exc:
        stderr = str(getattr(exc, "stderr", "") or "")
        if "BGZF EOF" not in stderr and "may be truncated" not in stderr:
            raise
        console.print(
            "[yellow]bcftools recusou o VCF por EOF BGZF ausente; "
            f"usando fallback gzip para {Path(vcf_path).name} {region}[/yellow]"
        )
        yield from _stream_gene_variants_from_gzip(vcf_path, sample_ids, region)


def _collect_vcf_variant_markers(
    *,
    config: Any,
    sample_ids: List[str],
    max_markers: int,
    min_frequency: float,
    include_snps: bool,
    include_indels: bool,
    genes_filter: Optional[set[str]] = None,
) -> List[Dict[str, Any]]:
    if not sample_ids:
        return []
    genes = list(config.dataset_input.genes_to_use or config.dataset_input.gene_order or [])
    if not genes:
        return []
    ontology_count = len(config.dataset_input.ontology_terms or []) or 1
    signal_track_count = 2 * ontology_count
    tracks_per_gene = signal_track_count + 3
    aligner = DynamicIndelAligner(
        Path(config.dataset_input.dataset_dir),
        selected_sample_ids=set(sample_ids),
        center_window_size=config.dataset_input.window_center_size,
    )
    counts: Dict[Tuple[str, str, int, str], Dict[str, Any]] = {}
    summary: Dict[str, Dict[str, Any]] = {}
    sample_count = max(len(sample_ids), 1)

    for gene_idx, gene_name in enumerate(genes):
        if genes_filter and gene_name not in genes_filter:
            continue
        gene_summary = summary.setdefault(gene_name, {
            "rows_seen": 0,
            "rows_after_type_filter": 0,
            "rows_inside_center_slice": 0,
            "non_ref_alleles": 0,
            "markers_added": 0,
            "error": None,
        })
        try:
            axis = aligner.get_alignment_axis(gene_name)
            center_slice = aligner.get_reference_centered_expanded_slice(gene_name, config.dataset_input.window_center_size)
            gene_payload = aligner._load_gene_mapping(gene_name)
        except Exception as exc:
            gene_summary["error"] = str(exc)
            console.print(f"[yellow]Variantes: pulando {gene_name}: {exc}[/yellow]")
            continue

        expanded_index_map = {int(k): int(v) for k, v in axis.get("expanded_index_map", {}).items()}
        expanded_start = int(center_slice["expanded_start"])
        expanded_end = int(center_slice["expanded_end"])
        ref_start_offset = int(axis.get("ref_start_offset", 0))
        full_start_1based = int(gene_payload.get("full_start_1based", gene_payload.get("start_1based", 0)))
        ref_sequence = str(gene_payload.get("ref_sequence", ""))
        ref_length = int(axis.get("ref_length", 0))

        try:
            rows = _stream_variants_with_gzip_fallback(gene_payload["vcf_path"], sample_ids, gene_payload["region"])
            for row in rows:
                gene_summary["rows_seen"] += 1
                ref = str(row["ref"])
                alt = str(row["alt"])
                pos0_full = int(row["pos_1based"]) - full_start_1based
                pos0 = _resolve_ref_index(ref_sequence, pos0_full - ref_start_offset, ref)
                if pos0 < 0 or pos0 >= ref_length:
                    continue
                behavior = _classify_length_behavior(ref, alt)
                if behavior == "length_change" and not include_indels:
                    continue
                if behavior != "length_change" and not include_snps:
                    continue
                gene_summary["rows_after_type_filter"] += 1
                ref_expanded = expanded_index_map.get(pos0)
                if ref_expanded is None or not (expanded_start <= ref_expanded < expanded_end):
                    continue
                gene_summary["rows_inside_center_slice"] += 1

                for sample_idx, gt in enumerate(row.get("gts", [])):
                    if sample_idx >= len(sample_ids):
                        break
                    gt_left, gt_right = _parse_gt(str(gt))
                    for haplotype, token in (("H1", gt_left), ("H2", gt_right)):
                        if token in {"0", ".", "nan", "None"}:
                            continue
                        allele = _allele_sequence(ref, alt, token)
                        if allele == ref:
                            continue
                        gene_summary["non_ref_alleles"] += 1
                        if behavior == "length_change":
                            entry = aligner.get_haplotype_entry(gene_name, sample_ids[sample_idx], haplotype)
                            if entry is None:
                                continue
                            if len(allele) > len(ref):
                                positions = [int(p) for p in entry.get("insertion_indices", [])]
                                variant_type = "INDEL_ins"
                                track_idx = signal_track_count
                            elif len(allele) < len(ref):
                                positions = [int(p) for p in entry.get("deletion_indices", [])]
                                variant_type = "INDEL_del"
                                track_idx = signal_track_count + 1
                            else:
                                positions = [ref_expanded]
                                variant_type = "SNP"
                                track_idx = 0
                        else:
                            positions = [ref_expanded]
                            variant_type = "SNP"
                            track_idx = 0

                        for expanded_pos in positions:
                            if not (expanded_start <= expanded_pos < expanded_end):
                                continue
                            x = int(expanded_pos - expanded_start)
                            key = (gene_name, haplotype, x, variant_type)
                            marker = counts.setdefault(key, {
                                "gene_name": gene_name,
                                "gene_index": gene_idx,
                                "haplotype": haplotype,
                                "x": x,
                                "chrom": row.get("chrom"),
                                "pos_1based": int(row["pos_1based"]),
                                "id": row.get("id"),
                                "ref": ref,
                                "alt": alt,
                                "variant_type": variant_type,
                                "count": 0,
                                "track_idx": track_idx,
                            })
                            marker["count"] += 1
                            gene_summary["markers_added"] += 1
        except Exception as exc:
            gene_summary["error"] = str(exc)
            console.print(f"[yellow]Variantes: erro em {gene_name}: {exc}[/yellow]")

    markers = []
    for marker in counts.values():
        frequency = float(marker["count"]) / sample_count
        if frequency < min_frequency:
            continue
        row_index = _variant_row_index(
            int(marker["gene_index"]),
            str(marker["haplotype"]),
            int(marker["track_idx"]),
            tracks_per_gene,
        )
        marker["frequency"] = frequency
        marker["row_index"] = row_index
        markers.append(marker)
    markers.sort(key=lambda item: (item["frequency"], item["count"]), reverse=True)
    selected = markers if max_markers <= 0 else markers[:max_markers]
    for marker in selected:
        marker["_summary_by_gene"] = summary.get(marker.get("gene_name"), {})
    selected.append({"_variant_collection_summary": summary})
    return selected


def _plot_mean_deeplift(
    mean_input: torch.Tensor,
    mean_attr: torch.Tensor,
    top_windows: List[Dict[str, Any]],
    config: Any,
    class_label: str,
    split: str,
    num_samples: int,
    output_path: Path,
    show_track_separators: bool = False,
    variant_markers: Optional[List[Dict[str, Any]]] = None,
    deeplift_cmap: str = "bwr",
) -> None:
    import matplotlib.pyplot as plt

    layout = _plot_layout(config)
    genes = layout["genes"]
    tracks_per_gene = layout["tracks_per_gene"]
    rows_per_hap = layout["rows_per_hap"]
    signal_track_count = layout["signal_track_count"]
    mask_track_count = layout["mask_track_count"]
    ontology_terms = list(config.dataset_input.ontology_terms or [])
    signal_track_labels = []
    if ontology_terms:
        for strand in ("+", "-"):
            for ontology in ontology_terms:
                signal_track_labels.append(f"{len(signal_track_labels)}:{ontology}{strand}")
    else:
        signal_track_labels.append("0:signal")

    input_np = _flatten_haps_for_plot(mean_input, config)
    attr_np = _flatten_haps_for_plot(mean_attr, config)
    signal_input_np, mask_input_np = _split_signal_and_masks(input_np, config)
    signal_attr_np, mask_attr_np = _split_signal_and_masks(attr_np, config)
    mask_input_np = np.clip(mask_input_np, 0.0, 1.0)
    abs_max = float(np.nanmax(np.abs(attr_np))) if attr_np.size else 0.0
    attr_lim = abs_max if abs_max > 0 else 1.0

    yticks = []
    ylabels = []
    if genes and input_np.shape[0] >= rows_per_hap:
        track_labels = list(signal_track_labels)
        mask_labels = ["ins", "del"]
        if config.dataset_input.indel_include_valid_mask:
            mask_labels.append("valid")
        if config.dataset_input.indel_include_snp_mask:
            mask_labels.append("snp")
        track_labels.extend(f"{signal_track_count + idx}:{name}" for idx, name in enumerate(mask_labels[:mask_track_count]))
        for gene_idx, gene in enumerate(genes):
            gene_offset = gene_idx * 2 * tracks_per_gene
            for hap_idx, hap in enumerate(("H1", "H2")):
                hap_offset = gene_offset + hap_idx * tracks_per_gene
                for track_idx, label in enumerate(track_labels):
                    yticks.append(hap_offset + track_idx)
                    ylabels.append(f"{hap}:{gene}:{label}")

    has_masks = mask_track_count > 0
    fig = plt.figure(figsize=(20, 18), constrained_layout=True)
    gs = fig.add_gridspec(2, 3 if has_masks else 2, width_ratios=[40, 1.4, 1.4] if has_masks else [40, 1.4])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    cax_signal = fig.add_subplot(gs[0, 1])
    cax_attr = fig.add_subplot(gs[1, 1])
    cax_mask = fig.add_subplot(gs[0, 2]) if has_masks else None
    cax_mask_attr = fig.add_subplot(gs[1, 2]) if has_masks else None
    gray_cmap = plt.colormaps["gray"].copy()
    gray_cmap.set_bad(alpha=0.0)
    bwr_cmap = plt.colormaps["bwr"].copy()
    bwr_cmap.set_bad(alpha=0.0)
    mask_cmap = LinearSegmentedColormap.from_list("mask_red_to_blue", ["#d73027", "#f7f7f7", "#2166ac"]).copy()
    mask_cmap.set_bad(alpha=0.0)

    im0 = axes[0].imshow(signal_input_np, aspect="auto", cmap=gray_cmap, interpolation="nearest", vmin=0)
    if has_masks:
        im0_mask = axes[0].imshow(mask_input_np, aspect="auto", cmap=mask_cmap, interpolation="nearest", vmin=0, vmax=1)
    axes[0].set_title(f"{split.upper()} SET | Class {class_label} ({num_samples} samples) | Input 2D Mean ({input_np.shape[0]}x{input_np.shape[1]})", fontsize=16, fontweight="bold")
    axes[0].set_ylabel("Tracks")
    axes[0].set_xlabel("Gene Position")
    if yticks:
        axes[0].set_yticks(yticks)
        axes[0].set_yticklabels(ylabels, fontsize=4.5)
    fig.colorbar(im0, cax=cax_signal, label="Ontology/strand tracks (gray)")
    if has_masks and cax_mask is not None:
        fig.colorbar(im0_mask, cax=cax_mask, label="Mask tracks (0 red -> 1 blue)")

    if deeplift_cmap == "gray_abs":
        attr_to_plot = np.abs(signal_attr_np)
        im1 = axes[1].imshow(attr_to_plot, aspect="auto", cmap=gray_cmap, interpolation="nearest", vmin=0, vmax=attr_lim)
        attr_colorbar_label = "Ontology/strand attribution magnitude"
    else:
        im1 = axes[1].imshow(signal_attr_np, aspect="auto", cmap=bwr_cmap, interpolation="nearest", vmin=-attr_lim, vmax=attr_lim)
        attr_colorbar_label = "Ontology/strand signed attribution"
    mask_attr_abs = np.abs(mask_attr_np)
    mask_attr_max = float(np.nanmax(mask_attr_abs)) if np.isfinite(mask_attr_abs).any() else 1.0
    if mask_attr_max <= 0:
        mask_attr_max = 1.0
    if has_masks:
        im1_mask = axes[1].imshow(mask_attr_np, aspect="auto", cmap=bwr_cmap, interpolation="nearest", vmin=-mask_attr_max, vmax=mask_attr_max)
    axes[1].set_title(f"DeepLIFT: Mean Attribution for Class {class_label} ({num_samples} samples)", fontsize=16, fontweight="bold")
    axes[1].set_ylabel("Tracks")
    axes[1].set_xlabel("Gene Position")
    if yticks:
        axes[1].set_yticks(yticks)
        axes[1].set_yticklabels(ylabels, fontsize=4.5)
    fig.colorbar(im1, cax=cax_attr, label=attr_colorbar_label)
    if has_masks and cax_mask_attr is not None:
        fig.colorbar(im1_mask, cax=cax_mask_attr, label="Mask attribution")

    if show_track_separators and genes and tracks_per_gene > 0:
        for ax in axes:
            for gene_idx in range(len(genes)):
                for block_offset in (gene_idx * 2 * tracks_per_gene, gene_idx * 2 * tracks_per_gene + tracks_per_gene):
                    if has_masks and signal_track_count > 0:
                        ax.axhline(block_offset + signal_track_count - 0.5, color="#f6f6f6", linewidth=0.35, alpha=0.7)
                    ax.axhline(block_offset + tracks_per_gene - 0.5, color="#dddddd", linewidth=0.55, alpha=0.7)
    if genes and tracks_per_gene > 0:
        legend_text = "Track order per haplotype: signals = " + ", ".join(signal_track_labels)
        if has_masks:
            mask_labels = ["ins", "del"]
            if config.dataset_input.indel_include_valid_mask:
                mask_labels.append("valid")
            if config.dataset_input.indel_include_snp_mask:
                mask_labels.append("snp")
            legend_text += "; masks = " + ", ".join(
                f"{signal_track_count + idx}:{name}" for idx, name in enumerate(mask_labels[:mask_track_count])
            )
        fig.text(0.02, 0.035, legend_text, fontsize=10)

    for win in top_windows[:10]:
        row_index = _plot_row_index_for_window(win, genes, rows_per_hap, tracks_per_gene)
        x = int(win.get("window_center", 0))
        for ax in axes:
            ax.scatter([x], [row_index], facecolors="none", edgecolors="lime", s=160, linewidths=2)

    if variant_markers:
        drawable_markers = _select_variant_markers_near_windows(
            variant_markers,
            top_windows,
            genes,
            rows_per_hap,
            tracks_per_gene,
            max_windows=10,
        )
        for marker in drawable_markers:
            color = "#ff00ff" if marker.get("variant_type") == "SNP" else "#00e5ff"
            x = int(marker["x"])
            y = int(marker["row_index"])
            axes[1].scatter(
                [x],
                [y],
                facecolors="none",
                edgecolors=color,
                s=80,
                linewidths=1.2,
                alpha=0.85,
            )
            axes[1].annotate(
                _variant_marker_label(marker),
                xy=(x, y),
                xytext=(6, 4),
                textcoords="offset points",
                fontsize=7,
                color=color,
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": color, "alpha": 0.75, "lw": 0.6},
            )
        if drawable_markers:
            axes[1].scatter([], [], facecolors="none", edgecolors="#ff00ff", s=80, linewidths=1.2, label="SNP")
            axes[1].scatter([], [], facecolors="none", edgecolors="#00e5ff", s=80, linewidths=1.2, label="INDEL")
            axes[1].legend(loc="lower right", fontsize=8, frameon=True, title="Nearest variants")

    top_labels = []
    for win in top_windows[:5]:
        gene = win.get("gene_name")
        track = win.get("track_name")
        value = float(win.get("mean_abs_attr", 0.0))
        top_labels.append(f"{gene}/{track}({value:.5f})")
    if top_labels:
        fig.text(0.02, 0.01, "Top regions: " + ", ".join(top_labels), fontsize=12)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def _window_too_close(candidate: Dict[str, Any], selected: List[Dict[str, Any]], min_distance_bp: int) -> bool:
    if min_distance_bp <= 0:
        return False
    cand_gene = candidate.get("gene_name")
    cand_center = int(candidate.get("window_center", 0))
    for row in selected:
        if row.get("gene_name") != cand_gene:
            continue
        if abs(cand_center - int(row.get("window_center", 0))) < min_distance_bp:
            return True
    return False


def _select_top_regions(
    attributions: torch.Tensor,
    config: Any,
    window_size: int,
    top_k: int,
    mode: str,
    min_distance_bp: int,
) -> List[Dict[str, Any]]:
    if mode == "global":
        candidates = find_top_deeplift_windows(
            attributions,
            config,
            window_size=window_size,
            top_k=max(top_k * 50, top_k),
        )
        selected: List[Dict[str, Any]] = []
        for candidate in candidates:
            if _window_too_close(candidate, selected, min_distance_bp):
                continue
            selected.append(candidate)
            if len(selected) >= top_k:
                break
        return selected

    attrs = attributions.detach().cpu()
    if attrs.ndim == 2:
        attrs = attrs.unsqueeze(0)
    if attrs.ndim != 3:
        raise ValueError(f"Esperado attributions [hap,row,L] ou [row,L], recebido {tuple(attrs.shape)}")

    layout = build_haplotype_track_layout(config)
    hap_names = ["H1", "H2"] if attrs.shape[0] == 2 else [f"H{i + 1}" for i in range(attrs.shape[0])]
    length = attrs.shape[-1]
    actual_window = min(window_size, length)
    best_by_gene: Dict[str, Dict[str, Any]] = {}

    for hap_idx, hap_name in enumerate(hap_names):
        for row_meta in layout[: attrs.shape[1]]:
            row = attrs[hap_idx, row_meta["row_index"]].abs().view(1, 1, -1)
            scores = torch.nn.functional.avg_pool1d(row, kernel_size=actual_window, stride=1).view(-1)
            if scores.numel() == 0:
                continue
            value, start = torch.max(scores, dim=0)
            candidate = {
                "haplotype": hap_name,
                **row_meta,
                "window_start": int(start.item()),
                "window_end": int(start.item() + actual_window),
                "window_center": int(start.item() + actual_window // 2),
                "window_size": int(actual_window),
                "mean_abs_attr": float(value.item()),
            }
            gene = str(candidate.get("gene_name"))
            current = best_by_gene.get(gene)
            if current is None or candidate["mean_abs_attr"] > current["mean_abs_attr"]:
                best_by_gene[gene] = candidate

    selected = sorted(best_by_gene.values(), key=lambda r: r["mean_abs_attr"], reverse=True)
    return selected[:top_k]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DeepLIFT interpretability for genotype_based_predictor")
    parser.add_argument("config_path", type=str, help="Path to the YAML config used for training")
    parser.add_argument("--checkpoint", type=str, default="best_accuracy.pt", help="Checkpoint file or absolute path")
    parser.add_argument("--split", choices=["train", "val", "test"], default="test")
    parser.add_argument("--target-class", type=str, default="predicted", help="predicted, true, class index, or class name")
    parser.add_argument("--filter-class", type=str, default=None, help="Use only samples whose true class matches this class name/index")
    parser.add_argument("--baseline", choices=["zeros", "mean"], default="zeros")
    parser.add_argument("--max-samples", type=int, default=50)
    parser.add_argument("--window-size", type=int, default=512)
    parser.add_argument("--top-k-windows", type=int, default=100)
    parser.add_argument("--top-regions-mode", choices=["per_gene", "global"], default="global")
    parser.add_argument("--min-distance-bp", type=int, default=0)
    parser.add_argument("--out-dir", type=str, default=None)
    parser.add_argument("--plot", action="store_true", help="Save a PNG similar to the legacy DeepLIFT class-mean panel")
    parser.add_argument("--deeplift-cmap", choices=["bwr", "gray_abs"], default="bwr", help="Colormap for DeepLIFT panel: signed red/blue or absolute grayscale")
    parser.add_argument("--plot-superpopulation-rnaseq", action="store_true", help="Save grayscale and line-plot RNA-Seq means for each superpopulation in the selected split")
    parser.add_argument("--plot-class-rnaseq", action="store_true", help="Save grayscale and line-plot RNA-Seq means for each target class in the selected split")
    parser.add_argument("--superpopulation-rnaseq-delta-reference", action="store_true", help="Treat superpopulation RNA-Seq values as sample-reference deltas; save centered grayscale map with top absolute deviations")
    parser.add_argument("--reference-predictions-dataset-dir", type=str, default=None, help="Override dataset_input.reference_predictions_dataset_dir for delta-reference cache lookup")
    parser.add_argument("--reference-predictions-sample-id", type=str, default=None, help="Override dataset_input.reference_predictions_sample_id for delta-reference cache lookup")
    parser.add_argument("--superpopulation-rnaseq-normalization-method", choices=["zscore", "minmax_keep_zero", "none"], default=None, help="Override dataset_input.normalization_method for superpopulation RNA-Seq cache lookup")
    parser.add_argument("--superpopulation-rnaseq-top-abs-markers", type=int, default=20, help="Number of green circles for largest absolute delta-reference deviations")
    parser.add_argument("--plot-top-from", choices=["mean", "samples"], default="mean", help="Circle top windows from the mean map or from per-sample windows")
    parser.add_argument("--save-raw-pixels", action="store_true", help="Save raw matrix PNGs with no axes/text/resize plus .npy matrices")
    parser.add_argument("--show-track-separators", action="store_true", help="Draw horizontal separators between signal/mask/gene blocks")
    parser.add_argument("--show-variant-markers", action="store_true", help="Draw variant markers on the DeepLIFT panel using the aligned cache axis")
    parser.add_argument("--variant-types", choices=["indel", "snp", "all"], default="indel", help="Variant marker source: indel is fast from cached masks; snp/all may query VCF once and cache")
    parser.add_argument("--variant-marker-min-frequency", type=float, default=0.0, help="Minimum fraction of selected samples carrying a variant marker")
    parser.add_argument("--variant-marker-max", type=int, default=0, help="Maximum number of variant markers to draw; 0 means no limit")
    parser.add_argument("--variant-marker-genes", type=str, default=None, help="Comma-separated genes for variant markers only, e.g. SLC45A2")
    parser.add_argument("--variant-marker-force-refresh", action="store_true", help="Recompute variant marker cache only for the requested marker settings")
    args = parser.parse_args()

    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    if args.superpopulation_rnaseq_normalization_method:
        config.dataset_input.normalization_method = args.superpopulation_rnaseq_normalization_method
    if args.superpopulation_rnaseq_delta_reference:
        config.dataset_input.alphagenome_signal_transform = "delta_reference"
        if args.reference_predictions_dataset_dir:
            config.dataset_input.reference_predictions_dataset_dir = str(Path(args.reference_predictions_dataset_dir).resolve())
        if args.reference_predictions_sample_id is not None:
            config.dataset_input.reference_predictions_sample_id = args.reference_predictions_sample_id
        if not config.dataset_input.reference_predictions_dataset_dir:
            raise ValueError(
                "--superpopulation-rnaseq-delta-reference requer --reference-predictions-dataset-dir "
                "quando a config nao define reference_predictions_dataset_dir"
            )
    if config.data_split.random_seed is not None and config.data_split.random_seed != -1:
        set_random_seeds(config.data_split.random_seed, config.data_split.strict_determinism)

    experiment_dir = Path(config.dataset_input.processed_cache_dir) / generate_experiment_name(config)
    experiment_dir.mkdir(parents=True, exist_ok=True)
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    dataset = _select_dataset(args.split, train_loader, val_loader, test_loader)

    if args.plot_superpopulation_rnaseq:
        out_dir = Path(args.out_dir) if args.out_dir else experiment_dir / "interpretability" / f"superpopulation_rnaseq_{args.split}"
        out_dir.mkdir(parents=True, exist_ok=True)
        means_by_superpop, counts_by_superpop, superpop_samples = _collect_superpopulation_mean_inputs(dataset, args.split)
        _save_superpopulation_rnaseq_visualizations(
            means_by_superpop=means_by_superpop,
            counts_by_superpop=counts_by_superpop,
            config=config,
            split=args.split,
            out_dir=out_dir,
            delta_reference=args.superpopulation_rnaseq_delta_reference,
            top_abs_markers=args.superpopulation_rnaseq_top_abs_markers,
        )
        _write_json(out_dir / "samples.json", superpop_samples)
        _write_json(out_dir / "run_metadata.json", {
            "config_path": str(config_path),
            "experiment_dir": str(experiment_dir),
            "split": args.split,
            "mode": "superpopulation_rnaseq",
            "superpopulation_rnaseq_delta_reference": args.superpopulation_rnaseq_delta_reference,
            "reference_predictions_dataset_dir": config.dataset_input.reference_predictions_dataset_dir,
            "reference_predictions_sample_id": config.dataset_input.reference_predictions_sample_id,
            "superpopulation_rnaseq_top_abs_markers": args.superpopulation_rnaseq_top_abs_markers,
            "num_samples": len(superpop_samples),
                f"counts_by_{group_label}": counts_by_superpop,
            "deeplift_computed": False,
        })
        console.print(f"[green]Visualizacoes RNA-Seq por superpopulacao salvas em:[/green] {out_dir}")
        return

    if args.plot_class_rnaseq:
        out_dir = Path(args.out_dir) if args.out_dir else experiment_dir / "interpretability" / f"class_rnaseq_{args.split}"
        out_dir.mkdir(parents=True, exist_ok=True)
        class_names = full_ds.get_class_names() if hasattr(full_ds, "get_class_names") else []
        means_by_class, counts_by_class, class_samples = _collect_class_mean_inputs(dataset, args.split, class_names)
        _save_superpopulation_rnaseq_visualizations(
            means_by_superpop=means_by_class,
            counts_by_superpop=counts_by_class,
            config=config,
            split=args.split,
            out_dir=out_dir,
            delta_reference=args.superpopulation_rnaseq_delta_reference,
            top_abs_markers=args.superpopulation_rnaseq_top_abs_markers,
            group_label="class_name",
            group_order=_ordered_class_groups(means_by_class.keys(), class_names),
        )
        _write_json(out_dir / "samples.json", class_samples)
        _write_json(out_dir / "run_metadata.json", {
            "config_path": str(config_path),
            "experiment_dir": str(experiment_dir),
            "split": args.split,
            "mode": "class_rnaseq",
            "class_names": class_names,
            "num_samples": len(class_samples),
            "counts_by_class_name": counts_by_class,
            "deeplift_computed": False,
        })
        console.print(f"[green]Visualizacoes RNA-Seq por classe salvas em:[/green] {out_dir}")
        return

    checkpoint_path = Path(args.checkpoint)
    if not checkpoint_path.is_absolute():
        checkpoint_path = experiment_dir / "models" / args.checkpoint
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint nao encontrado: {checkpoint_path}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = _build_model(config, full_ds).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device)
    state = checkpoint.get("model_state_dict", checkpoint)
    model.load_state_dict(state)
    model.eval()

    checkpoint_slug = checkpoint_path.stem
    out_dir = Path(args.out_dir) if args.out_dir else experiment_dir / "interpretability" / f"deeplift_{checkpoint_slug}_{args.split}"
    out_dir.mkdir(parents=True, exist_ok=True)

    class_names = full_ds.get_class_names() if hasattr(full_ds, "get_class_names") else []
    class_name_to_idx = {name: idx for idx, name in enumerate(class_names)}
    deeplift = DeepLIFT(model)
    filter_class_idx = _class_arg_to_idx(args.filter_class, class_name_to_idx) if args.filter_class is not None else None
    variant_marker_genes = None
    if args.variant_marker_genes:
        variant_marker_genes = {g.strip() for g in args.variant_marker_genes.split(",") if g.strip()}

    per_sample_summaries: List[Dict[str, Any]] = []
    all_track_rows: List[Dict[str, Any]] = []
    all_window_rows: List[Dict[str, Any]] = []
    input_sum: Optional[torch.Tensor] = None
    attr_sum: Optional[torch.Tensor] = None
    plotted_class_label: Optional[str] = None
    selected_sample_ids: List[str] = []
    max_samples = args.max_samples if args.max_samples and args.max_samples > 0 else len(dataset)

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Calculando DeepLIFT...", total=max_samples)
        for idx in range(len(dataset)):
            if len(per_sample_summaries) >= max_samples:
                break
            features, target, meta_idx = _unpack_item(dataset[idx])
            true_idx = _target_to_int(target)
            if filter_class_idx is not None and true_idx != filter_class_idx:
                continue
            sample_id = _sample_id_for(dataset, idx, meta_idx)
            input_tensor = features.unsqueeze(0).to(device)
            with torch.no_grad():
                logits = model(input_tensor)
                predicted_idx = int(logits.argmax(dim=1).item())

            if args.target_class == "predicted":
                target_idx = predicted_idx
            elif args.target_class == "true":
                if true_idx is None:
                    raise ValueError("--target-class true requer target escalar")
                target_idx = true_idx
            elif args.target_class in class_name_to_idx:
                target_idx = class_name_to_idx[args.target_class]
            else:
                target_idx = int(args.target_class)

            attrs, _ = deeplift.generate(input_tensor, target_class=target_idx, baseline_type=args.baseline, dataset=dataset)
            if input_sum is None:
                input_sum = features.detach().cpu().clone()
                attr_sum = attrs.detach().cpu().clone()
            else:
                input_sum += features.detach().cpu()
                attr_sum += attrs.detach().cpu()
            track_rows = summarize_deeplift_by_track(attrs, config)
            window_rows = _select_top_regions(
                attrs,
                config,
                window_size=args.window_size,
                top_k=args.top_k_windows,
                mode=args.top_regions_mode,
                min_distance_bp=args.min_distance_bp,
            )

            target_name = class_names[target_idx] if 0 <= target_idx < len(class_names) else str(target_idx)
            if plotted_class_label is None:
                plotted_class_label = target_name if args.target_class not in {"predicted", "true"} else args.target_class
            pred_name = class_names[predicted_idx] if 0 <= predicted_idx < len(class_names) else str(predicted_idx)
            true_name = class_names[true_idx] if true_idx is not None and 0 <= true_idx < len(class_names) else str(true_idx)
            per_sample_summaries.append({
                "sample_index": idx,
                "sample_id": sample_id,
                "true_class": true_name,
                "predicted_class": pred_name,
                "deeplift_target_class": target_name,
            })
            selected_sample_ids.append(sample_id)

            for row in track_rows:
                row.update({"sample_index": idx, "sample_id": sample_id, "deeplift_target_class": target_name})
            for row in window_rows:
                row.update({"sample_index": idx, "sample_id": sample_id, "deeplift_target_class": target_name})
            all_track_rows.extend(track_rows)
            all_window_rows.extend(window_rows)
            progress.update(task, advance=1)

    n_samples = len(per_sample_summaries)
    if n_samples == 0:
        raise RuntimeError("Nenhuma amostra encontrada para os filtros solicitados")

    aggregate_by_track = _aggregate_metric(
        all_track_rows,
        "haplotype|gene_name|track_type|track_name|ontology_curie|strand",
    )
    aggregate_by_gene = _aggregate_metric(all_track_rows, "haplotype|gene_name")
    aggregate_by_track_type = _aggregate_metric(all_track_rows, "track_type|track_name")

    _write_json(out_dir / "run_metadata.json", {
        "config_path": str(config_path),
        "experiment_dir": str(experiment_dir),
        "checkpoint_path": str(checkpoint_path),
        "split": args.split,
        "target_class": args.target_class,
        "filter_class": args.filter_class,
        "baseline": args.baseline,
        "num_samples": n_samples,
        "window_size": args.window_size,
        "top_k_windows_per_sample": args.top_k_windows,
        "top_regions_mode": args.top_regions_mode,
        "min_distance_bp": args.min_distance_bp,
        "save_raw_pixels": args.save_raw_pixels,
        "show_track_separators": args.show_track_separators,
        "plot_superpopulation_rnaseq": args.plot_superpopulation_rnaseq,
        "show_variant_markers": args.show_variant_markers,
        "variant_types": args.variant_types,
        "variant_marker_min_frequency": args.variant_marker_min_frequency,
        "variant_marker_max": args.variant_marker_max,
        "variant_marker_genes": sorted(variant_marker_genes or []),
        "variant_marker_force_refresh": args.variant_marker_force_refresh,
    })
    _write_json(out_dir / "samples.json", per_sample_summaries)
    _write_json(out_dir / "aggregate_by_track.json", aggregate_by_track)
    _write_json(out_dir / "aggregate_by_gene.json", aggregate_by_gene)
    _write_json(out_dir / "aggregate_by_track_type.json", aggregate_by_track_type)
    _write_json(out_dir / "top_windows.json", all_window_rows)

    _write_csv(out_dir / "aggregate_by_track.csv", aggregate_by_track)
    _write_csv(out_dir / "aggregate_by_gene.csv", aggregate_by_gene)
    _write_csv(out_dir / "aggregate_by_track_type.csv", aggregate_by_track_type)
    _write_csv(out_dir / "top_windows.csv", all_window_rows)

    if args.plot and input_sum is not None and attr_sum is not None:
        mean_input = input_sum / n_samples
        mean_attr = attr_sum / n_samples
        variant_markers = None
        if args.show_variant_markers:
            if args.variant_types == "indel":
                if getattr(config.dataset_input, "feature_mode", "signals_and_masks") == "signals_only":
                    raise ValueError(
                        "--variant-types indel usa tracks de mascara no tensor, mas feature_mode='signals_only'. "
                        "Use --variant-types snp ou --variant-types all para coletar variantes via VCF."
                    )
                variant_markers = _collect_indel_markers_from_mean_input(
                    mean_input=mean_input,
                    config=config,
                    max_markers=args.variant_marker_max,
                    min_frequency=args.variant_marker_min_frequency,
                    genes_filter=variant_marker_genes,
                )
            else:
                cache_path = _variant_cache_path(
                    out_dir,
                    selected_sample_ids,
                    args.variant_types,
                    args.variant_marker_min_frequency,
                    args.variant_marker_max,
                    variant_marker_genes,
                )
                if cache_path.exists() and not args.variant_marker_force_refresh:
                    with open(cache_path) as f:
                        variant_markers = json.load(f)
                    console.print(f"[green]Marcadores de variantes carregados do cache:[/green] {cache_path}")
                else:
                    variant_markers = _collect_vcf_variant_markers(
                        config=config,
                        sample_ids=selected_sample_ids,
                        max_markers=args.variant_marker_max,
                        min_frequency=args.variant_marker_min_frequency,
                        include_snps=args.variant_types in {"snp", "all"},
                        include_indels=args.variant_types == "all",
                        genes_filter=variant_marker_genes,
                    )
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    _write_json(cache_path, variant_markers)
            _write_json(out_dir / "variant_markers.json", variant_markers)
            for item in variant_markers:
                if "_variant_collection_summary" in item:
                    _write_json(out_dir / "variant_collection_summary.json", item["_variant_collection_summary"])
                    break
        plot_windows = all_window_rows
        if args.plot_top_from == "mean":
            plot_windows = _select_top_regions(
                mean_attr,
                config,
                window_size=args.window_size,
                top_k=min(args.top_k_windows, 100),
                mode=args.top_regions_mode,
                min_distance_bp=args.min_distance_bp,
            )
        _plot_mean_deeplift(
            mean_input=mean_input,
            mean_attr=mean_attr,
            top_windows=plot_windows,
            config=config,
            class_label=plotted_class_label or args.target_class,
            split=args.split,
            num_samples=n_samples,
            output_path=out_dir / "deeplift_mean.png",
            show_track_separators=args.show_track_separators,
            variant_markers=variant_markers,
            deeplift_cmap=args.deeplift_cmap,
        )
        if args.save_raw_pixels:
            _save_raw_pixel_images(mean_input, mean_attr, config, out_dir / "raw_pixels")

    console.print(f"[green]DeepLIFT salvo em:[/green] {out_dir}")


if __name__ == "__main__":
    main()
