#!/usr/bin/env python3
"""CLI to materialize the dataset layout consumed by genotype_based_predictor."""

from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console

from genotype_based_predictor.data.layout import materialize_dataset, update_window_variant_sources

console = Console()


def main() -> int:
    parser = argparse.ArgumentParser(description="Materializa o dataset para o layout usado por genotype_based_predictor")
    parser.add_argument("source_dataset_dir", help="Diretorio do dataset de origem")
    parser.add_argument("target_dataset_dir", help="Diretorio do dataset de destino")
    parser.add_argument("--vcf-pattern", help="Padrao dos VCFs por cromossomo, com placeholder {chrom}")
    parser.add_argument("--vcf-root-dir", help="Diretorio com VCFs por cromossomo")
    parser.add_argument("--update-variant-sources", action="store_true", help="Atualiza apenas raw_variant_source nos metadados do dataset de destino")
    args = parser.parse_args()

    source_dir = Path(args.source_dataset_dir).resolve()
    target_dir = Path(args.target_dataset_dir).resolve()
    if not args.update_variant_sources:
        materialize_dataset(source_dir, target_dir)
    update_window_variant_sources(target_dir, vcf_pattern=args.vcf_pattern, vcf_root_dir=args.vcf_root_dir)
    console.print(f"[bold green]Concluido[/bold green]: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
