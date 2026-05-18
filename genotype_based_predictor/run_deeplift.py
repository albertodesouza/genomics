from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import torch
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


def _plot_mean_deeplift(
    mean_input: torch.Tensor,
    mean_attr: torch.Tensor,
    top_windows: List[Dict[str, Any]],
    config: Any,
    class_label: str,
    split: str,
    num_samples: int,
    output_path: Path,
) -> None:
    import matplotlib.pyplot as plt

    def flatten_haps(x: torch.Tensor) -> np.ndarray:
        arr = x.detach().cpu()
        if arr.ndim == 3:
            arr = arr.reshape(arr.shape[0] * arr.shape[1], arr.shape[2])
        elif arr.ndim != 2:
            raise ValueError(f"Esperado tensor 2D/3D para plot, recebeu {tuple(arr.shape)}")
        return arr.numpy()

    input_np = flatten_haps(mean_input)
    attr_np = flatten_haps(mean_attr)
    abs_max = float(np.nanmax(np.abs(attr_np))) if attr_np.size else 0.0
    attr_lim = abs_max if abs_max > 0 else 1.0

    genes = list(config.dataset_input.genes_to_use or config.dataset_input.gene_order or [])
    tracks_per_gene = 2 * len(config.dataset_input.ontology_terms or []) + 3
    rows_per_hap = len(genes) * tracks_per_gene
    yticks = []
    ylabels = []
    if genes and input_np.shape[0] >= rows_per_hap:
        for hap_idx in range(max(input_np.shape[0] // rows_per_hap, 1)):
            prefix = f"H{hap_idx + 1}:" if input_np.shape[0] > rows_per_hap else ""
            offset = hap_idx * rows_per_hap
            for gene_idx, gene in enumerate(genes):
                yticks.append(offset + gene_idx * tracks_per_gene + (tracks_per_gene - 1) / 2)
                ylabels.append(f"{prefix}{gene}")

    fig, axes = plt.subplots(2, 1, figsize=(16, 10), constrained_layout=True)
    im0 = axes[0].imshow(input_np, aspect="auto", cmap="gray", interpolation="nearest", vmin=0)
    axes[0].set_title(f"{split.upper()} SET | Class {class_label} ({num_samples} samples) | Input 2D Mean ({input_np.shape[0]}x{input_np.shape[1]})", fontsize=16, fontweight="bold")
    axes[0].set_ylabel("Genes / tracks")
    axes[0].set_xlabel("Gene Position")
    if yticks:
        axes[0].set_yticks(yticks)
        axes[0].set_yticklabels(ylabels, fontsize=8)
    fig.colorbar(im0, ax=axes[0], label="Normalized Value")

    im1 = axes[1].imshow(attr_np, aspect="auto", cmap="bwr", interpolation="nearest", vmin=-attr_lim, vmax=attr_lim)
    axes[1].set_title(f"DeepLIFT: Mean Attribution for Class {class_label} ({num_samples} samples)", fontsize=16, fontweight="bold")
    axes[1].set_ylabel("Genes / tracks")
    axes[1].set_xlabel("Gene Position")
    if yticks:
        axes[1].set_yticks(yticks)
        axes[1].set_yticklabels(ylabels, fontsize=8)
    fig.colorbar(im1, ax=axes[1], label="Attribution (+ class, - not class)")

    for win in top_windows[:10]:
        row_index = int(win.get("row_index", 0))
        if str(win.get("haplotype")) == "H2":
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DeepLIFT interpretability for genotype_based_predictor")
    parser.add_argument("config_path", type=str, help="Path to the YAML config used for training")
    parser.add_argument("--checkpoint", type=str, default="best_accuracy.pt", help="Checkpoint file or absolute path")
    parser.add_argument("--split", choices=["train", "val", "test"], default="test")
    parser.add_argument("--target-class", type=str, default="predicted", help="predicted, true, class index, or class name")
    parser.add_argument("--baseline", choices=["zeros", "mean"], default="zeros")
    parser.add_argument("--max-samples", type=int, default=50)
    parser.add_argument("--window-size", type=int, default=512)
    parser.add_argument("--top-k-windows", type=int, default=100)
    parser.add_argument("--out-dir", type=str, default=None)
    parser.add_argument("--plot", action="store_true", help="Save a PNG similar to the legacy DeepLIFT class-mean panel")
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

    per_sample_summaries: List[Dict[str, Any]] = []
    all_track_rows: List[Dict[str, Any]] = []
    all_window_rows: List[Dict[str, Any]] = []
    input_sum: Optional[torch.Tensor] = None
    attr_sum: Optional[torch.Tensor] = None
    plotted_class_label: Optional[str] = None
    n_samples = min(args.max_samples, len(dataset)) if args.max_samples and args.max_samples > 0 else len(dataset)

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Calculando DeepLIFT...", total=n_samples)
        for idx in range(n_samples):
            features, target, meta_idx = _unpack_item(dataset[idx])
            sample_id = _sample_id_for(dataset, idx, meta_idx)
            input_tensor = features.unsqueeze(0).to(device)
            with torch.no_grad():
                logits = model(input_tensor)
                predicted_idx = int(logits.argmax(dim=1).item())
            true_idx = _target_to_int(target)

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
            window_rows = find_top_deeplift_windows(attrs, config, window_size=args.window_size, top_k=args.top_k_windows)

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
        "baseline": args.baseline,
        "num_samples": n_samples,
        "window_size": args.window_size,
        "top_k_windows_per_sample": args.top_k_windows,
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
        mean_input = input_sum / max(n_samples, 1)
        mean_attr = attr_sum / max(n_samples, 1)
        _plot_mean_deeplift(
            mean_input=mean_input,
            mean_attr=mean_attr,
            top_windows=all_window_rows,
            config=config,
            class_label=plotted_class_label or args.target_class,
            split=args.split,
            num_samples=n_samples,
            output_path=out_dir / "deeplift_mean.png",
        )

    console.print(f"[green]DeepLIFT salvo em:[/green] {out_dir}")


if __name__ == "__main__":
    main()
