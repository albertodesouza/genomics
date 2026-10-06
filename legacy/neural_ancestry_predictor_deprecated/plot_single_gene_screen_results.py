#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt


CURATED_GENE_ORDER = [
    "MC1R",
    "TYRP1",
    "TYR",
    "SLC45A2",
    "DDB1",
    "EDAR",
    "MFSD12",
    "OCA2",
    "HERC2",
    "SLC24A5",
    "TCHH",
]


def load_rows(csv_path: Path) -> List[Dict[str, str]]:
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def metric_value(row: Dict[str, str], metric: str) -> float:
    raw = row.get(metric, "")
    if raw is None or raw == "":
        return float("nan")
    return float(raw)


def sort_rows(rows: List[Dict[str, str]], metric: str) -> List[Dict[str, str]]:
    by_gene = {row["gene"]: row for row in rows}
    curated_rows = [by_gene[gene] for gene in CURATED_GENE_ORDER if gene in by_gene]
    curated_rows.sort(key=lambda row: metric_value(row, metric), reverse=True)
    curated_set = {row["gene"] for row in curated_rows}
    random_rows = sorted(
        (row for row in rows if row["gene"] not in curated_set),
        key=lambda row: metric_value(row, metric),
        reverse=True,
    )
    return curated_rows + random_rows


def plot_metric(rows: List[Dict[str, str]], metric: str, output_path: Path) -> None:
    ordered_rows = sort_rows(rows, metric)
    genes = [row["gene"] for row in ordered_rows]
    values = [metric_value(row, metric) for row in ordered_rows]
    colors = ["#2f6db0" if row["group"] == "curated11" else "#b8c2cc" for row in ordered_rows]

    width = max(16, len(genes) * 0.45)
    fig, ax = plt.subplots(figsize=(width, 7))
    bars = ax.bar(range(len(genes)), values, color=colors, edgecolor="black", linewidth=0.5)

    ontology = ordered_rows[0].get("ontology", "unknown") if ordered_rows else "unknown"
    ax.set_title(f"Single-Gene Pigmentation Results: {ontology} / {metric}")
    ax.set_ylabel(metric)
    ax.set_xlabel("Gene")
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90)
    ax.set_ylim(0.0, 1.0)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    curated_boundary = len([gene for gene in genes if gene in CURATED_GENE_ORDER])
    if curated_boundary:
        ax.axvline(curated_boundary - 0.5, color="black", linestyle=":", linewidth=1)
        ax.text(
            max(curated_boundary / 2 - 0.5, 0),
            0.98,
            "Curated 11",
            ha="center",
            va="top",
            fontsize=10,
            transform=ax.get_xaxis_transform(),
        )
        if curated_boundary < len(genes):
            ax.text(
                curated_boundary + (len(genes) - curated_boundary) / 2 - 0.5,
                0.98,
                "Random genes",
                ha="center",
                va="top",
                fontsize=10,
                transform=ax.get_xaxis_transform(),
            )

    for bar, value in zip(bars, values):
        if value == value:
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.01, f"{value:.2f}", ha="center", va="bottom", fontsize=7, rotation=90)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def sanitize_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value)


def group_rows_by_ontology(rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    grouped: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        ontology = row.get("ontology") or "unknown"
        grouped.setdefault(ontology, []).append(row)
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot single-gene pigmentation screening results per ontology")
    parser.add_argument("--results-csv", required=True, help="Path to results.csv")
    parser.add_argument(
        "--metric",
        default="test_accuracy",
        choices=["test_accuracy", "test_f1", "test_macro_f1", "val_accuracy", "best_val_accuracy"],
        help="Metric to plot",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output image path or directory. If multiple ontologies are present, one file is created per ontology.",
    )
    args = parser.parse_args()

    rows = load_rows(Path(args.results_csv))
    completed_rows = [row for row in rows if row.get("status") == "completed" and row.get(args.metric, "") not in ("", None)]
    if not completed_rows:
        raise ValueError(f"No completed rows with metric {args.metric} found in {args.results_csv}")

    grouped_rows = group_rows_by_ontology(completed_rows)
    output_path = Path(args.output)

    if len(grouped_rows) == 1:
        plot_metric(completed_rows, args.metric, output_path)
        print(f"Saved plot to: {output_path}")
        return

    output_path.mkdir(parents=True, exist_ok=True)
    for ontology, ontology_rows in sorted(grouped_rows.items()):
        ontology_output = output_path / f"{sanitize_name(ontology)}_{args.metric}.png"
        plot_metric(ontology_rows, args.metric, ontology_output)
        print(f"Saved plot to: {ontology_output}")


if __name__ == "__main__":
    main()
