#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from hf_genomics_dataset.builder import build_hf_dataset_store
from hf_genomics_dataset.models import ConversionOptions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert existing filesystem genomics datasets into a normalized local HF dataset store"
    )
    parser.add_argument(
        "--source",
        action="append",
        required=True,
        help="Source dataset directory. Repeat the flag to merge multiple datasets.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for the local Hugging Face dataset store.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=0,
        help="Optional number of rows per temporary write shard. Use 0 to keep all rows in memory.",
    )
    parser.add_argument(
        "--no-sequences",
        action="store_true",
        help="Do not carry ref/H1/H2 sequence fields into gene_windows rows.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_dirs = [Path(value).resolve() for value in args.source]
    output_dir = Path(args.output).resolve()
    options = ConversionOptions(
        chunk_size=max(0, args.chunk_size),
        include_sequences=not args.no_sequences,
    )

    summary = build_hf_dataset_store(source_dirs=source_dirs, output_dir=output_dir, options=options)
    print(f"Converted dataset store at {output_dir}")
    print(f"Samples: {summary.total_samples}")
    print(f"Gene windows: {summary.total_gene_windows}")
    print(f"Genes: {len(summary.genes)}")


if __name__ == "__main__":
    main()
