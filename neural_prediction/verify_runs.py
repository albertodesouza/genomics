#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify neural_prediction screening runs")
    parser.add_argument("--results-csv", required=True)
    args = parser.parse_args()

    results_csv = Path(args.results_csv)
    with results_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    missing = [row for row in rows if row["status"] == "completed" and not row.get("test_accuracy")]
    failed = [row for row in rows if row["status"] != "completed"]

    print(f"Total rows: {len(rows)}")
    print(f"Failed rows: {len(failed)}")
    print(f"Completed rows missing test_accuracy: {len(missing)}")
    for row in missing:
        print(f"MISSING {row['gene']} {row['group']}")
    for row in failed:
        print(f"FAILED {row['gene']} rc={row['returncode']}")


if __name__ == "__main__":
    main()
