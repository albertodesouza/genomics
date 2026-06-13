from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from genomics.predictors.snp_ancestry.pipeline import _build_statistics_path, load_config


def _safe_float(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def fst_like(frequencies: Sequence[float]) -> float:
    """Return the same multi-class Fst-like score used by the ancestry pipeline."""
    if not frequencies:
        return 0.0
    p_bar = sum(frequencies) / len(frequencies)
    if p_bar <= 0.0 or p_bar >= 1.0:
        return 0.0
    var_b = sum((p - p_bar) ** 2 for p in frequencies) / len(frequencies)
    return var_b / (p_bar * (1.0 - p_bar))


def marker_scores(rsid: str, stats: dict) -> dict:
    populations = stats["metadata"]["populations"]
    freqs = [_safe_float(v) for v in stats["allele_frequencies"][rsid]]
    chrom, pos = stats.get("snp_info", {}).get(rsid, ["", ""])
    overall = sum(freqs) / len(freqs) if freqs else 0.0
    maf = min(overall, 1.0 - overall)
    max_frequency = max(freqs) if freqs else 0.0
    min_frequency = min(freqs) if freqs else 0.0
    max_delta_frequency = max_frequency - min_frequency
    max_population = populations[freqs.index(max_frequency)] if freqs else ""
    min_population = populations[freqs.index(min_frequency)] if freqs else ""

    row = {
        "rsid": rsid,
        "chrom": chrom,
        "position": pos,
        "ref_allele": stats.get("ref_alleles", {}).get(rsid, ""),
        "maf": maf,
        "fst": fst_like(freqs),
        "max_delta_frequency": max_delta_frequency,
        "max_population": max_population,
        "min_population": min_population,
    }
    for pop, freq in zip(populations, freqs):
        row[f"freq_{pop}"] = freq
    return row


def load_statistics(config_path: Path, statistics_path: Optional[Path] = None) -> Tuple[dict, Path]:
    if statistics_path is None:
        config = load_config(config_path)
        statistics_path = Path(_build_statistics_path(config))
    if not statistics_path.exists():
        raise FileNotFoundError(
            f"Statistics file not found: {statistics_path}. Run `genomics snp-ancestry run --config {config_path}` first."
        )
    with open(statistics_path, "r", encoding="utf-8") as f:
        return json.load(f), statistics_path


def build_marker_table(
    stats: dict,
    score: str = "fst",
    top: Optional[int] = None,
    min_maf: Optional[float] = None,
    min_score: Optional[float] = None,
) -> List[dict]:
    valid_scores = {"fst", "maf", "max_delta_frequency"}
    if score not in valid_scores:
        raise ValueError(f"Unsupported score '{score}'. Choose one of: {', '.join(sorted(valid_scores))}")

    rows = [marker_scores(rsid, stats) for rsid in stats.get("allele_frequencies", {})]
    if min_maf is not None:
        rows = [row for row in rows if row["maf"] >= min_maf]
    if min_score is not None:
        rows = [row for row in rows if row[score] >= min_score]

    rows.sort(key=lambda row: (-_safe_float(row[score]), str(row["chrom"]), _position_sort_key(row["position"]), row["rsid"]))
    if top is not None:
        rows = rows[:top]
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        row["score"] = row[score]
        row["score_name"] = score
    return rows


def _position_sort_key(value: object) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _format_value(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return f"{value:.8g}"
    return value


def write_marker_table(rows: Sequence[dict], output: Path, populations: Iterable[str]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "rank",
        "rsid",
        "chrom",
        "position",
        "ref_allele",
        "score_name",
        "score",
        "fst",
        "maf",
        "max_delta_frequency",
        "max_population",
        "min_population",
    ]
    fieldnames.extend(f"freq_{pop}" for pop in populations)
    with open(output, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format_value(row.get(key, "")) for key in fieldnames})


def run_export(args: argparse.Namespace) -> int:
    stats, statistics_path = load_statistics(args.config, args.statistics)
    rows = build_marker_table(
        stats,
        score=args.score,
        top=args.top,
        min_maf=args.min_maf,
        min_score=args.min_score,
    )
    write_marker_table(rows, args.output, stats["metadata"]["populations"])
    print(f"Statistics: {statistics_path}")
    print(f"Markers: {len(rows)}")
    print(f"Output: {args.output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export ranked ancestry-informative markers from SNP ancestry statistics")
    parser.add_argument("--config", type=Path, required=True, help="SNP ancestry YAML config used to build statistics")
    parser.add_argument("--statistics", type=Path, default=None, help="Explicit statistics JSON path; defaults to the config-derived path")
    parser.add_argument("--output", type=Path, required=True, help="TSV output path")
    parser.add_argument("--score", choices=["fst", "maf", "max_delta_frequency"], default="fst")
    parser.add_argument("--top", type=int, default=None, help="Keep only the top N markers after filtering")
    parser.add_argument("--min-maf", type=float, default=None, help="Filter markers below this MAF")
    parser.add_argument("--min-score", type=float, default=None, help="Filter markers below the selected score")
    parser.set_defaults(func=run_export)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
