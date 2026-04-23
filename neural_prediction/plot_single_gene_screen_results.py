#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def load_rows(results_csv: Path):
    with results_csv.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot neural_prediction single-gene screening results")
    parser.add_argument("--results-csv", required=True)
    parser.add_argument("--metric", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = load_rows(Path(args.results_csv))
    valid_rows = [row for row in rows if row.get(args.metric) not in ("", None)]
    for row in valid_rows:
        row[args.metric] = float(row[args.metric])

    curated = sorted([row for row in valid_rows if row["group"] == "curated11"], key=lambda row: row[args.metric], reverse=True)
    random_rows = sorted([row for row in valid_rows if row["group"] != "curated11"], key=lambda row: row[args.metric], reverse=True)
    ordered = curated + random_rows

    genes = [row["gene"] for row in ordered]
    values = [row[args.metric] for row in ordered]
    colors = ["#4c78a8" if row["group"] == "curated11" else "#f58518" for row in ordered]

    plt.figure(figsize=(max(12, len(genes) * 0.5), 6))
    plt.bar(range(len(genes)), values, color=colors)
    plt.xticks(range(len(genes)), genes, rotation=90)
    plt.ylabel(args.metric)
    plt.tight_layout()
    plt.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
