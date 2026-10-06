from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np


def _load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _load_tsv(path: Path) -> List[dict]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _float_value(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def discover_model_dirs(ml_dir: Path) -> List[Path]:
    if not ml_dir.exists():
        raise FileNotFoundError(f"ML directory not found: {ml_dir}")
    return sorted(path for path in ml_dir.iterdir() if path.is_dir() and (path / "metrics.json").exists())


def metric_labels(metrics: dict) -> List[str]:
    report = metrics.get("classification_report", "")
    labels = []
    for line in str(report).splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[0] not in {"accuracy", "macro", "weighted"}:
            labels.append(parts[0])
    matrix = metrics.get("confusion_matrix") or []
    if labels and len(labels) == len(matrix):
        return labels
    return [str(i) for i in range(len(matrix))]


def plot_confusion_matrix(metrics: dict, output: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    matrix = np.asarray(metrics.get("confusion_matrix") or [], dtype=float)
    if matrix.size == 0:
        return
    labels = metric_labels(metrics)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.1), max(5, len(labels) * 0.95)))
    im = ax.imshow(matrix, cmap="Blues")
    ax.set_title(title)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    threshold = matrix.max() / 2.0 if matrix.size else 0.0
    for row_idx in range(matrix.shape[0]):
        for col_idx in range(matrix.shape[1]):
            value = matrix[row_idx, col_idx]
            ax.text(
                col_idx,
                row_idx,
                str(int(value)),
                ha="center",
                va="center",
                color="white" if value > threshold else "black",
            )
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_feature_importance(rows: Sequence[dict], output: Path, title: str, top: int) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    ranked = sorted(rows, key=lambda row: _float_value(row.get("importance")), reverse=True)[:top]
    if not ranked:
        return
    labels = [str(row.get("rsid", "")) for row in ranked][::-1]
    values = [_float_value(row.get("importance")) for row in ranked][::-1]
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, max(4.5, len(labels) * 0.28)))
    ax.barh(np.arange(len(labels)), values, color="#4059ad")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Importance")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_ablation_curves(rows: Sequence[dict], output: Path, metrics: Sequence[str]) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    if not rows:
        return
    models = sorted({str(row.get("model", "")) for row in rows if row.get("model")})
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(metrics), 1, figsize=(8.5, max(4.2, 3.2 * len(metrics))), squeeze=False)
    for ax, metric in zip(axes[:, 0], metrics):
        for model in models:
            model_rows = [row for row in rows if row.get("model") == model]
            model_rows.sort(key=lambda row: _float_value(row.get("remove_top")))
            x = [_float_value(row.get("remove_top")) for row in model_rows]
            y = [_float_value(row.get(metric)) for row in model_rows]
            ax.plot(x, y, marker="o", linewidth=2, label=model)
        ax.set_ylim(0.0, 1.02)
        ax.set_xlabel("Top AIMs removed")
        ax.set_ylabel(metric.replace("_", " "))
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle("AIM ablation robustness", y=0.995)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_ml_dir(ml_dir: Path, output_dir: Path, splits: Iterable[str], top_features: int) -> List[Path]:
    created = []
    for model_dir in discover_model_dirs(ml_dir):
        model_name = model_dir.name
        metrics_path = model_dir / "metrics.json"
        payload = _load_json(metrics_path)
        split_metrics: Dict[str, dict] = payload.get("metrics", {})
        for split_name in splits:
            metrics = split_metrics.get(split_name)
            if not metrics:
                continue
            output = output_dir / model_name / f"confusion_matrix_{split_name}.png"
            plot_confusion_matrix(metrics, output, f"{model_name} {split_name} confusion matrix")
            created.append(output)
        importance_path = model_dir / "feature_importance.tsv"
        if importance_path.exists():
            output = output_dir / model_name / f"feature_importance_top{top_features}.png"
            plot_feature_importance(
                _load_tsv(importance_path),
                output,
                f"{model_name} top {top_features} AIM importances",
                top_features,
            )
            created.append(output)
    return created


def plot_ablation_dir(ablation_dir: Path, output_dir: Path, metrics: Sequence[str]) -> List[Path]:
    path = ablation_dir / "ablation.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Ablation TSV not found: {path}")
    output = output_dir / "ablation_curves.png"
    plot_ablation_curves(_load_tsv(path), output, metrics)
    return [output]


def run_plots(args: argparse.Namespace) -> int:
    try:
        import matplotlib  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "Missing matplotlib dependency. Install with `python3 -m pip install -e \".[snp-ancestry]\"` or `python3 -m pip install matplotlib`."
        ) from exc

    if args.ml_dir is None and args.ablation_dir is None:
        raise SystemExit("At least one of --ml-dir or --ablation-dir is required")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    created: List[Path] = []
    if args.ml_dir is not None:
        created.extend(plot_ml_dir(args.ml_dir, args.output_dir / "ml", args.splits, args.top_features))
    if args.ablation_dir is not None:
        created.extend(plot_ablation_dir(args.ablation_dir, args.output_dir / "ablation", args.metrics))

    for path in created:
        print(f"Plot: {path}")
    print(f"Output: {args.output_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot SNP ancestry ML and AIM-ablation artifacts")
    parser.add_argument("--ml-dir", type=Path, default=None, help="Directory produced by `genomics snp-ancestry train-ml`")
    parser.add_argument("--ablation-dir", type=Path, default=None, help="Directory produced by `genomics snp-ancestry ablate`")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--splits", nargs="+", choices=["train", "val", "test"], default=["test", "val", "train"])
    parser.add_argument("--top-features", type=int, default=50)
    parser.add_argument("--metrics", nargs="+", choices=["accuracy", "balanced_accuracy", "macro_f1", "weighted_f1"], default=["accuracy", "balanced_accuracy", "macro_f1"])
    parser.set_defaults(func=run_plots)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
