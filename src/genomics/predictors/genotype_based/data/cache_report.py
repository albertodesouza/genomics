from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

import numpy as np


def _read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    return payload if isinstance(payload, dict) else {}


def _sample_records(dataset_metadata: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    records: Dict[str, Dict[str, Any]] = {}
    for key in ("individual_metadata", "sample_metadata", "samples_metadata", "pedigree"):
        value = dataset_metadata.get(key)
        if isinstance(value, dict):
            for sample_id, record in value.items():
                if isinstance(record, dict):
                    records[str(sample_id)] = record
    individuals = dataset_metadata.get("individuals")
    if isinstance(individuals, list):
        for item in individuals:
            if isinstance(item, dict):
                sample_id = item.get("sample_id") or item.get("SampleID") or item.get("id")
                if sample_id:
                    records.setdefault(str(sample_id), item)
    return records


def _load_sample_records(dataset_dir: Path, dataset_metadata: Mapping[str, Any], sample_ids: Iterable[str]) -> Dict[str, Dict[str, Any]]:
    records = _sample_records(dataset_metadata)
    individuals_dir = Path(dataset_dir) / "individuals"
    for sample_id in sample_ids:
        sample_id = str(sample_id)
        if sample_id in records:
            continue
        metadata_path = individuals_dir / sample_id / "individual_metadata.json"
        if metadata_path.exists():
            records[sample_id] = _read_json(metadata_path)
    return records


def _metadata_value(record: Mapping[str, Any], *keys: str, default: str = "UNK") -> str:
    for key in keys:
        value = record.get(key)
        if value not in (None, ""):
            return str(value)
    return default


def _normalize_metadata_label(field: str, value: str) -> str:
    if field == "sex":
        mapping = {
            "1": "male",
            "1.0": "male",
            "M": "male",
            "m": "male",
            "male": "male",
            "Male": "male",
            "2": "female",
            "2.0": "female",
            "F": "female",
            "f": "female",
            "female": "female",
            "Female": "female",
            "0": "unknown",
            "0.0": "unknown",
            "UNK": "unknown",
        }
        return mapping.get(str(value), str(value))
    return str(value)


def _counter_for(samples: Iterable[str], records: Mapping[str, Mapping[str, Any]], field: str, *keys: str) -> Dict[str, int]:
    counter: Counter[str] = Counter()
    for sample_id in samples:
        value = _metadata_value(records.get(sample_id, {}), *keys)
        counter[_normalize_metadata_label(field, value)] += 1
    return dict(sorted(counter.items()))


def build_cache_report(cache_dir: Path) -> Dict[str, Any]:
    cache_dir = Path(cache_dir)
    metadata = _read_json(cache_dir / "metadata.json")
    split_index = _read_json(cache_dir / "split_index.json")
    view_definition = _read_json(cache_dir / "view_definition.json") if (cache_dir / "view_definition.json").exists() else {}
    dataset_dir = Path(metadata.get("dataset_dir") or metadata.get("processing_params", {}).get("dataset_dir") or "")
    dataset_metadata = _read_json(dataset_dir / "dataset_metadata.json") if dataset_dir and (dataset_dir / "dataset_metadata.json").exists() else {}
    splits = {name: [str(item) for item in split_index.get(name, [])] for name in ("train", "val", "test")}
    all_samples = [sample_id for samples in splits.values() for sample_id in samples]
    records = _load_sample_records(dataset_dir, dataset_metadata, all_samples)

    split_stats = {}
    for split_name, samples in splits.items():
        split_stats[split_name] = {
            "count": len(samples),
            "superpopulation": _counter_for(samples, records, "superpopulation", "superpopulation", "Superpopulation"),
            "population": _counter_for(samples, records, "population", "population", "Population"),
            "sex": _counter_for(samples, records, "sex", "sex", "Sex"),
        }

    class_names = metadata.get("class_names") or {}
    class_labels = [class_names[key] for key in sorted(class_names, key=lambda item: int(item))] if isinstance(class_names, dict) else []
    report = {
        "cache_dir": str(cache_dir),
        "dataset_dir": str(dataset_dir) if dataset_dir else None,
        "created_at": metadata.get("creation_date"),
        "total_samples": metadata.get("total_samples", len(all_samples)),
        "model_input_shape": metadata.get("model_input_shape"),
        "num_classes": metadata.get("num_classes"),
        "class_names": class_labels,
        "gene_count": len(metadata.get("gene_order") or []),
        "genes": metadata.get("gene_order") or [],
        "tracks_per_gene": metadata.get("tracks_per_gene"),
        "processing_params": metadata.get("processing_params", {}),
        "requested_view": view_definition.get("requested_view", {}),
        "target_definition": view_definition.get("target_definition", {}),
        "splits": split_stats,
    }
    return report


def save_cache_report(cache_dir: Path, console: Optional[Any] = None) -> Path:
    cache_dir = Path(cache_dir)
    report = build_cache_report(cache_dir)
    report_dir = cache_dir / "report"
    plots_dir = report_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / "dataset_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    _write_markdown_report(report, report_dir / "dataset_report.md")
    _save_plots(report, plots_dir, console)
    if console is not None:
        console.print(f"[green]Dataset report salvo:[/green] {report_dir}")
    return report_path


def _write_markdown_report(report: Mapping[str, Any], path: Path) -> None:
    lines = [
        "# Genotype Dataset Cache Report",
        "",
        f"Cache: `{report.get('cache_dir')}`",
        f"Dataset: `{report.get('dataset_dir')}`",
        f"Created at: `{report.get('created_at')}`",
        "",
        "## Summary",
        "",
        f"- Total samples: {report.get('total_samples')}",
        f"- Input shape: `{report.get('model_input_shape')}`",
        f"- Classes: {', '.join(report.get('class_names') or [])}",
        f"- Genes: {report.get('gene_count')}",
        f"- Tracks per gene: {report.get('tracks_per_gene')}",
        "",
        "## Splits",
        "",
        "| Split | Count |",
        "|---|---:|",
    ]
    for split_name, stats in (report.get("splits") or {}).items():
        lines.append(f"| {split_name} | {stats.get('count', 0)} |")
    lines.extend(["", "## Plots", "", "See `plots/` for PNG summaries.", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def _save_plots(report: Mapping[str, Any], plots_dir: Path, console: Optional[Any]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError:
        if console is not None:
            console.print("[yellow]matplotlib indisponível; plots de dataset não foram gerados.[/yellow]")
        return

    splits = report.get("splits") or {}
    split_names = list(splits)
    counts = [int(splits[name].get("count", 0)) for name in split_names]
    _bar_plot(plt, split_names, counts, "Examples per split", "Examples", plots_dir / "split_counts.png")

    for field in ("superpopulation", "population", "sex"):
        labels = sorted({label for stats in splits.values() for label in (stats.get(field) or {})})
        if not labels:
            continue
        fig, ax = plt.subplots(figsize=(max(9, len(labels) * 0.9), 5.5))
        x = np.arange(len(labels))
        width = min(0.8 / max(len(split_names), 1), 0.25)
        for split_name in split_names:
            offset = (split_names.index(split_name) - (len(split_names) - 1) / 2) * width
            values = [int((splits[split_name].get(field) or {}).get(label, 0)) for label in labels]
            bars = ax.bar(x + offset, values, width=width, label=split_name)
            _annotate_bars(ax, bars)
        ax.set_title(f"{field.title()} distribution by split")
        ax.set_ylabel("Examples")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / f"{field}_by_split.png", dpi=160)
        plt.close(fig)

    genes = report.get("genes") or []
    if genes:
        _bar_plot(plt, genes, [1] * len(genes), "Genes included in view", "Included", plots_dir / "genes_included.png")


def _bar_plot(plt: Any, labels: List[str], values: List[int], title: str, ylabel: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 0.7), 4.5))
    bars = ax.bar(labels, values)
    _annotate_bars(ax, bars)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def _annotate_bars(ax: Any, bars: Any) -> None:
    for bar in bars:
        height = bar.get_height()
        if height == 0:
            continue
        ax.annotate(
            f"{int(height)}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
        )
