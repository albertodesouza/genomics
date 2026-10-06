#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable

from genomics.core.metrics import public_classification_results


def _is_zeroish(value: Any) -> bool:
    try:
        return float(value) == 0.0
    except (TypeError, ValueError):
        return False


def iter_result_files(paths: Iterable[Path]) -> Iterable[Path]:
    for path in paths:
        if path.is_file():
            if path.name.endswith("results.json") or path.name == "best_summary.json":
                yield path
            elif path.name == "search_results.csv":
                yield path
            continue
        if path.is_dir():
            yield from path.rglob("*results.json")
            yield from path.rglob("best_summary.json")
            yield from path.rglob("search_results.csv")


def _convert_search_row(row: Dict[str, Any]) -> bool:
    aliases = {
        "val_accuracy": "val_weighted_accuracy",
        "val_precision": "val_weighted_precision",
        "val_recall": "val_weighted_recall",
        "val_f1": "val_weighted_f1_score",
    }
    changed = False
    for old_key, new_key in aliases.items():
        if old_key in row and (new_key not in row or _is_zeroish(row.get(new_key))):
            row[new_key] = row[old_key]
            changed = True
        if old_key in row:
            row.pop(old_key, None)
            changed = True
    return changed


def _convert_selection_metric(payload: Dict[str, Any]) -> bool:
    selection_metric = payload.get("selection_metric")
    aliases = {
        "val_accuracy": "val_weighted_accuracy",
        "val_precision": "val_weighted_precision",
        "val_recall": "val_weighted_recall",
        "val_f1": "val_weighted_f1_score",
    }
    if selection_metric in aliases:
        payload["selection_metric"] = aliases[selection_metric]
        return True
    return False


def _convert_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    required_weighted = {"weighted_accuracy", "weighted_precision", "weighted_recall", "weighted_f1_score"}
    required_legacy = {"accuracy", "precision", "recall", "f1"}
    required_recomputable = {"confusion_matrix", "per_class_metrics"}
    if required_legacy.issubset(payload) or required_recomputable.issubset(payload):
        return public_classification_results(payload)

    updated = dict(payload)
    changed = _convert_selection_metric(updated)
    best = updated.get("best")
    if isinstance(best, dict):
        changed = _convert_search_row(best) or changed
    rows = updated.get("results")
    if isinstance(rows, list):
        new_rows = []
        for row in rows:
            if isinstance(row, dict):
                row = dict(row)
                changed = _convert_search_row(row) or changed
            new_rows.append(row)
        updated["results"] = new_rows
    if _convert_search_row(updated):
        changed = True
    return updated if changed else payload


def update_file(path: Path, *, dry_run: bool) -> bool:
    if path.suffix == ".csv":
        return update_csv_file(path, dry_run=dry_run)
    with open(path, "r", encoding="utf-8") as f:
        payload: Dict[str, Any] = json.load(f)
    updated = _convert_payload(payload)
    if updated == payload:
        return False
    if not dry_run:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(updated, f, indent=2, ensure_ascii=False)
            f.write("\n")
    return True


def update_csv_file(path: Path, *, dry_run: bool) -> bool:
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = [dict(row) for row in reader]
        fieldnames = list(reader.fieldnames or [])
    aliases = {
        "val_accuracy": "val_weighted_accuracy",
        "val_precision": "val_weighted_precision",
        "val_recall": "val_weighted_recall",
        "val_f1": "val_weighted_f1_score",
    }
    if not any(old_key in fieldnames for old_key in aliases):
        return False
    for row in rows:
        for old_key, new_key in aliases.items():
            if old_key in row and new_key not in row:
                row[new_key] = row[old_key]
            row.pop(old_key, None)
    new_fieldnames = []
    for fieldname in fieldnames:
        if fieldname in aliases:
            mapped = aliases[fieldname]
            if mapped not in new_fieldnames:
                new_fieldnames.append(mapped)
        elif fieldname not in new_fieldnames:
            new_fieldnames.append(fieldname)
    if not dry_run:
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=new_fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description="Add explicit weighted metric aliases to result JSON files.")
    parser.add_argument("paths", nargs="+", type=Path, help="Result JSON files or directories to scan")
    parser.add_argument("--dry-run", action="store_true", help="Print files that would be updated without writing")
    args = parser.parse_args()

    updated = 0
    scanned = 0
    for path in sorted(set(iter_result_files(args.paths))):
        scanned += 1
        if update_file(path, dry_run=args.dry_run):
            updated += 1
            print(f"updated: {path}")
    print(f"scanned={scanned} updated={updated} dry_run={args.dry_run}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
