#!/usr/bin/env python3
"""CLI to materialize the dataset layout consumed by genotype_based_predictor."""

from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console

from genotype_based_predictor.dataset_layout import materialize_dataset

console = Console()


def main() -> int:
    parser = argparse.ArgumentParser(description="Materializa o dataset para o layout usado por genotype_based_predictor")
    parser.add_argument("source_dataset_dir", help="Diretorio do dataset de origem")
    parser.add_argument("target_dataset_dir", help="Diretorio do dataset de destino")
    args = parser.parse_args()

    source_dir = Path(args.source_dataset_dir).resolve()
    target_dir = Path(args.target_dataset_dir).resolve()
    materialize_dataset(source_dir, target_dir)
    console.print(f"[bold green]Concluido[/bold green]: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
