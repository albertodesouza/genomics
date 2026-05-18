from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import torch
from matplotlib.colors import LinearSegmentedColormap
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from genotype_based_predictor.config import generate_experiment_name, load_config
from genotype_based_predictor.data_pipeline import prepare_data
from genotype_based_predictor.interpretability import (
    DeepLIFT,
    find_top_deeplift_windows,
    summarize_deeplift_by_track,
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
    tracks_per_gene = 2 * len(config.dataset_input.ontology_terms or []) + 3
    return {
        "genes": genes,
        "tracks_per_gene": tracks_per_gene,
        "rows_per_hap": len(genes) * tracks_per_gene,
        "signal_track_count": tracks_per_gene - 3,
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
    signal = matrix.copy()
    masks = np.full_like(matrix, np.nan, dtype=np.float32)
    if tracks_per_gene >= 4:
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

    attr_lim = float(np.nanmax(np.abs(signal_attr))) if np.isfinite(signal_attr).any() else 1.0
    mask_attr_lim = float(np.nanmax(np.abs(mask_attr))) if np.isfinite(mask_attr).any() else 1.0
    if attr_lim <= 0:
        attr_lim = 1.0
    if mask_attr_lim <= 0:
        mask_attr_lim = 1.0
    mask_cmap = LinearSegmentedColormap.from_list("mask_red_to_blue", ["#d73027", "#f7f7f7", "#2166ac"])

    plt.imsave(out_dir / "raw_input_signal_gray.png", signal_input, cmap="gray", vmin=0)
    plt.imsave(out_dir / "raw_input_masks_0_1.png", np.clip(mask_input, 0.0, 1.0), cmap=mask_cmap, vmin=0, vmax=1)
    plt.imsave(out_dir / "raw_deeplift_signal_bwr.png", signal_attr, cmap="bwr", vmin=-attr_lim, vmax=attr_lim)
    plt.imsave(out_dir / "raw_deeplift_masks_bwr.png", mask_attr, cmap="bwr", vmin=-mask_attr_lim, vmax=mask_attr_lim)

    np.save(out_dir / "raw_input_matrix.npy", input_np)
    np.save(out_dir / "raw_deeplift_matrix.npy", attr_np)


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
) -> None:
    import matplotlib.pyplot as plt

    layout = _plot_layout(config)
    genes = layout["genes"]
    tracks_per_gene = layout["tracks_per_gene"]
    rows_per_hap = layout["rows_per_hap"]
    signal_track_count = layout["signal_track_count"]
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
        track_labels = list(signal_track_labels) + [
            f"{signal_track_count}:ins",
            f"{signal_track_count + 1}:del",
            f"{signal_track_count + 2}:valid",
        ]
        for gene_idx, gene in enumerate(genes):
            gene_offset = gene_idx * 2 * tracks_per_gene
            for hap_idx, hap in enumerate(("H1", "H2")):
                hap_offset = gene_offset + hap_idx * tracks_per_gene
                for track_idx, label in enumerate(track_labels):
                    yticks.append(hap_offset + track_idx)
                    ylabels.append(f"{hap}:{gene}:{label}")

    fig = plt.figure(figsize=(20, 18), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[40, 1.4, 1.4])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])]
    cax_signal = fig.add_subplot(gs[0, 1])
    cax_mask = fig.add_subplot(gs[0, 2])
    cax_attr = fig.add_subplot(gs[1, 1])
    cax_mask_attr = fig.add_subplot(gs[1, 2])
    gray_cmap = plt.colormaps["gray"].copy()
    gray_cmap.set_bad(alpha=0.0)
    bwr_cmap = plt.colormaps["bwr"].copy()
    bwr_cmap.set_bad(alpha=0.0)
    mask_cmap = LinearSegmentedColormap.from_list("mask_red_to_blue", ["#d73027", "#f7f7f7", "#2166ac"]).copy()
    mask_cmap.set_bad(alpha=0.0)

    im0 = axes[0].imshow(signal_input_np, aspect="auto", cmap=gray_cmap, interpolation="nearest", vmin=0)
    im0_mask = axes[0].imshow(mask_input_np, aspect="auto", cmap=mask_cmap, interpolation="nearest", vmin=0, vmax=1)
    axes[0].set_title(f"{split.upper()} SET | Class {class_label} ({num_samples} samples) | Input 2D Mean ({input_np.shape[0]}x{input_np.shape[1]})", fontsize=16, fontweight="bold")
    axes[0].set_ylabel("Tracks")
    axes[0].set_xlabel("Gene Position")
    if yticks:
        axes[0].set_yticks(yticks)
        axes[0].set_yticklabels(ylabels, fontsize=4.5)
    fig.colorbar(im0, cax=cax_signal, label="Ontology/strand tracks (gray)")
    fig.colorbar(im0_mask, cax=cax_mask, label="Mask tracks (0 red -> 1 blue)")

    im1 = axes[1].imshow(signal_attr_np, aspect="auto", cmap=bwr_cmap, interpolation="nearest", vmin=-attr_lim, vmax=attr_lim)
    mask_attr_abs = np.abs(mask_attr_np)
    mask_attr_max = float(np.nanmax(mask_attr_abs)) if np.isfinite(mask_attr_abs).any() else 1.0
    if mask_attr_max <= 0:
        mask_attr_max = 1.0
    im1_mask = axes[1].imshow(mask_attr_np, aspect="auto", cmap=bwr_cmap, interpolation="nearest", vmin=-mask_attr_max, vmax=mask_attr_max)
    axes[1].set_title(f"DeepLIFT: Mean Attribution for Class {class_label} ({num_samples} samples)", fontsize=16, fontweight="bold")
    axes[1].set_ylabel("Tracks")
    axes[1].set_xlabel("Gene Position")
    if yticks:
        axes[1].set_yticks(yticks)
        axes[1].set_yticklabels(ylabels, fontsize=4.5)
    fig.colorbar(im1, cax=cax_attr, label="Ontology/strand attribution")
    fig.colorbar(im1_mask, cax=cax_mask_attr, label="Mask attribution")

    if show_track_separators and genes and tracks_per_gene >= 4:
        for ax in axes:
            for gene_idx in range(len(genes)):
                for block_offset in (gene_idx * 2 * tracks_per_gene, gene_idx * 2 * tracks_per_gene + tracks_per_gene):
                    ax.axhline(block_offset + signal_track_count - 0.5, color="#f6f6f6", linewidth=0.35, alpha=0.7)
                    ax.axhline(block_offset + tracks_per_gene - 0.5, color="#dddddd", linewidth=0.55, alpha=0.7)
    if genes and tracks_per_gene >= 4:
        legend_text = (
            "Track order per haplotype: "
            f"signals = {', '.join(signal_track_labels)}; masks = "
            f"{signal_track_count}:ins, {signal_track_count + 1}:del, {signal_track_count + 2}:valid"
        )
        fig.text(0.02, 0.035, legend_text, fontsize=10)

    for win in top_windows[:10]:
        row_index = int(win.get("row_index", 0))
        if genes and rows_per_hap > 0:
            gene_idx = int(win.get("gene_index", row_index // tracks_per_gene))
            track_idx = int(win.get("track_index_in_gene", row_index % tracks_per_gene))
            hap_offset = tracks_per_gene if str(win.get("haplotype")) == "H2" else 0
            row_index = gene_idx * 2 * tracks_per_gene + hap_offset + track_idx
        elif str(win.get("haplotype")) == "H2":
            row_index += rows_per_hap
        x = int(win.get("window_center", 0))
        for ax in axes:
            ax.scatter([x], [row_index], facecolors="none", edgecolors="lime", s=160, linewidths=2)

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

    candidates = find_top_deeplift_windows(
        attributions,
        config,
        window_size=window_size,
        top_k=max(top_k * 50, top_k),
    )
    best_by_gene: Dict[str, Dict[str, Any]] = {}
    for candidate in candidates:
        gene = str(candidate.get("gene_name"))
        if gene not in best_by_gene:
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
    parser.add_argument("--plot-top-from", choices=["mean", "samples"], default="mean", help="Circle top windows from the mean map or from per-sample windows")
    parser.add_argument("--save-raw-pixels", action="store_true", help="Save raw matrix PNGs with no axes/text/resize plus .npy matrices")
    parser.add_argument("--show-track-separators", action="store_true", help="Draw horizontal separators between signal/mask/gene blocks")
    args = parser.parse_args()

    config_path = Path(args.config_path).resolve()
    config = load_config(config_path)
    if config.data_split.random_seed is not None and config.data_split.random_seed != -1:
        set_random_seeds(config.data_split.random_seed, config.data_split.strict_determinism)

    experiment_dir = Path(config.dataset_input.processed_cache_dir) / generate_experiment_name(config)
    experiment_dir.mkdir(parents=True, exist_ok=True)
    full_ds, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    dataset = _select_dataset(args.split, train_loader, val_loader, test_loader)

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

    per_sample_summaries: List[Dict[str, Any]] = []
    all_track_rows: List[Dict[str, Any]] = []
    all_window_rows: List[Dict[str, Any]] = []
    input_sum: Optional[torch.Tensor] = None
    attr_sum: Optional[torch.Tensor] = None
    plotted_class_label: Optional[str] = None
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
        )
        if args.save_raw_pixels:
            _save_raw_pixel_images(mean_input, mean_attr, config, out_dir / "raw_pixels")

    console.print(f"[green]DeepLIFT salvo em:[/green] {out_dir}")


if __name__ == "__main__":
    main()
