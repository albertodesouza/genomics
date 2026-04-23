#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

from hf_genomics_dataset.validation import inspect_hf_dataset_store, validate_hf_dataset_store


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect and validate a local HF genomics dataset store")
    parser.add_argument("--dataset", required=True, help="Path to the built HF dataset store")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the inspection result as JSON",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset).resolve()

    validation = validate_hf_dataset_store(dataset_dir)
    inspection = inspect_hf_dataset_store(dataset_dir) if validation["valid"] else None

    result = {
        "validation": validation,
        "inspection": inspection,
    }

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    print(f"Dataset: {dataset_dir}")
    print(f"Valid: {validation['valid']}")
    if validation["errors"]:
        print("Errors:")
        for error in validation["errors"]:
            print(f"  - {error}")
    if validation["warnings"]:
        print("Warnings:")
        for warning in validation["warnings"]:
            print(f"  - {warning}")

    if inspection is None:
        return

    summary = inspection["summary"]
    print(f"Samples: {summary['samples']}")
    print(f"Gene windows: {summary['gene_windows']}")
    print(f"Genes: {summary['genes']}")
    print(f"Source datasets: {', '.join(summary['source_datasets'])}")
    print(f"Populations: {', '.join(summary['populations'])}")
    print(f"Superpopulations: {', '.join(summary['superpopulations'])}")
    print("Top genes:")
    for gene, count in summary["top_genes"]:
        print(f"  - {gene}: {count}")


if __name__ == "__main__":
    main()
