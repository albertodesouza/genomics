"""Free-form sample metadata for VCF imports: CSV, TSV, JSON or PLINK ``.fam``.

Nothing about the columns is assumed. The sample-ID column is detected as the one whose values
match the most VCF sample names; every other column becomes a sample field (facet, training
target, ...). Two optional roles map onto the canonical pedigree fields used elsewhere:
``family_id`` (keeps relatives in the same train/val/test split) and ``sex`` / ``sex_label``.
"""
from __future__ import annotations

import csv
import io
import json
import math
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

ID_COLUMN_NAMES = ("sample_id", "sample", "samples", "sampleid", "iid", "id", "individual", "individual_id", "subject", "subject_id", "name")
FAMILY_COLUMN_NAMES = ("family_id", "familyid", "family", "fid", "pedigree")
SEX_COLUMN_NAMES = ("sex", "gender", "sex_label")
FAM_COLUMNS = ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"]
MAX_CATEGORIES = 50
_NUMBER_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_SEX_LABELS = {"1": "Male", "2": "Female", "m": "Male", "male": "Male", "f": "Female", "female": "Female", "man": "Male", "woman": "Female"}


def clean_field_name(name: Any) -> str:
    text = re.sub(r"\s+", "_", str(name or "").strip().lstrip("﻿"))
    return re.sub(r"[^\w.\-]", "", text) or "field"


def coerce_value(value: Any) -> Any:
    """Strings that look numeric become int/float; blanks and NA markers become None."""
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return None if isinstance(value, float) and not math.isfinite(value) else value
    if isinstance(value, (list, dict)):
        return json.dumps(value, sort_keys=True)
    text = str(value).strip()
    if text == "" or text.lower() in ("na", "n/a", "nan", "none", "null"):
        return None
    if _NUMBER_RE.match(text):
        try:
            number = float(text)
        except ValueError:
            return text
        if number.is_integer() and "." not in text and "e" not in text.lower():
            return int(number)
        return number if math.isfinite(number) else None
    return text


def parse_table(text: str, filename: str = "") -> Tuple[List[str], List[Dict[str, Any]]]:
    """Parse a metadata file into ``(columns, rows)``; rows keep raw string values."""
    text = text.lstrip("﻿")
    stripped = text.strip()
    if not stripped:
        return [], []
    lower = filename.lower()
    if lower.endswith(".json") or stripped[0] in "[{":
        return _parse_json(json.loads(stripped))
    if lower.endswith(".fam"):
        rows = [dict(zip(FAM_COLUMNS, line.split())) for line in stripped.splitlines() if line.strip()]
        return FAM_COLUMNS, rows
    first = stripped.splitlines()[0]
    if "\t" in first:
        delimiter = "\t"
    elif first.count(";") > first.count(","):
        delimiter = ";"
    elif "," in first:
        delimiter = ","
    else:
        delimiter = None  # whitespace separated
    if delimiter is None:
        lines = [line.split() for line in stripped.splitlines() if line.strip()]
        header = [clean_field_name(h) for h in lines[0]]
        return header, [dict(zip(header, values)) for values in lines[1:]]
    reader = csv.reader(io.StringIO(stripped), delimiter=delimiter)
    lines = [row for row in reader if any(cell.strip() for cell in row)]
    header = _unique([clean_field_name(h) for h in lines[0]])
    return header, [dict(zip(header, row)) for row in lines[1:]]


def _parse_json(data: Any) -> Tuple[List[str], List[Dict[str, Any]]]:
    if isinstance(data, dict):
        for key in ("samples", "individuals", "records", "data"):
            if isinstance(data.get(key), (list, dict)):
                return _parse_json(data[key])
        if all(isinstance(v, dict) for v in data.values()):  # {sample: {field: value}}
            rows = [{"sample_id": sample, **{clean_field_name(k): v for k, v in fields.items()}} for sample, fields in data.items()]
            return _columns(rows), rows
        raise ValueError("JSON metadata must be a list of records or an object mapping sample ids to records")
    if isinstance(data, list):
        rows = [{clean_field_name(k): v for k, v in item.items()} for item in data if isinstance(item, dict)]
        if not rows:
            raise ValueError("JSON metadata list has no records")
        return _columns(rows), rows
    raise ValueError("Unsupported JSON metadata layout")


def _columns(rows: Iterable[Dict[str, Any]]) -> List[str]:
    names: List[str] = []
    for row in rows:
        for key in row:
            if key not in names:
                names.append(key)
    return names


def _unique(names: List[str]) -> List[str]:
    seen: Dict[str, int] = {}
    out = []
    for name in names:
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        out.append(name)
    return out


def guess_id_column(columns: Sequence[str], rows: Sequence[Dict[str, Any]], samples: Optional[Iterable[str]] = None) -> Optional[str]:
    """The column whose values match the most VCF samples (ties: conventional names first)."""
    if not columns:
        return None
    sample_set = set(samples or [])
    def rank(column: str) -> Tuple[int, int, int]:
        values = [str(r.get(column, "")).strip() for r in rows]
        matches = sum(1 for v in values if v in sample_set) if sample_set else 0
        unique = len(set(values)) == len(values)
        named = 1 if column.lower() in ID_COLUMN_NAMES else 0
        return (matches, named, int(unique))
    best = max(columns, key=lambda c: (rank(c), -columns.index(c)))
    return best


def guess_role(columns: Sequence[str], names: Sequence[str]) -> Optional[str]:
    lowered = {c.lower(): c for c in columns}
    for name in names:
        if name in lowered:
            return lowered[name]
    return None


def sex_label(value: Any) -> Optional[str]:
    if value is None:
        return None
    return _SEX_LABELS.get(str(value).strip().lower())


def records_by_sample(
    rows: Sequence[Dict[str, Any]],
    id_column: str,
    family_column: Optional[str] = None,
    sex_column: Optional[str] = None,
    drop_columns: Sequence[str] = (),
) -> Dict[str, Dict[str, Any]]:
    """``{sample: {field: value}}`` with typed values and the canonical role fields filled in."""
    out: Dict[str, Dict[str, Any]] = {}
    drop = set(drop_columns) | {id_column}
    for row in rows:
        sample = str(row.get(id_column) or "").strip()
        if not sample:
            continue
        record = {clean_field_name(k): coerce_value(v) for k, v in row.items() if k not in drop and k}
        if family_column and row.get(family_column) not in (None, ""):
            record["family_id"] = str(coerce_value(row.get(family_column)))
        if sex_column:
            label = sex_label(row.get(sex_column))
            if label:
                record.setdefault("sex_label", label)
                record["sex"] = 1 if label == "Male" else 2
        out[sample] = {k: v for k, v in record.items() if v is not None}
    return out


def describe_fields(records: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Per-field summary: categorical (with value counts), numeric or text."""
    names = _columns(records.values())
    fields = []
    for name in names:
        values = [r[name] for r in records.values() if r.get(name) is not None]
        counts: Dict[str, int] = {}
        for value in values:
            counts[str(value)] = counts.get(str(value), 0) + 1
        numeric = bool(values) and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values)
        categorical = 1 < len(counts) <= MAX_CATEGORIES and not (numeric and len(counts) > 12)
        fields.append({
            "name": name,
            "kind": "categorical" if categorical else ("numeric" if numeric else "text"),
            "distinct": len(counts),
            "missing": len(records) - len(values),
            "counts": sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:MAX_CATEGORIES] if categorical else [],
        })
    return fields
