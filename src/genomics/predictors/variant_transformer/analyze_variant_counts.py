from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import numpy as np
import torch
from rich.console import Console
from rich.table import Table

console = Console()


def _load_json(path: Path):
    with open(path) as f:
        return json.load(f)


def _percentile(values: List[int], q: float) -> float:
    if not values:
        return 0.0
    return float(np.percentile(np.asarray(values, dtype=np.float64), q))


def _summary(values: List[int]) -> Dict[str, float]:
    if not values:
        return {"mean": 0.0, "median": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0}
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "p95": _percentile(values, 95),
        "p99": _percentile(values, 99),
        "max": float(arr.max()),
    }


def analyze_variant_counts(processed_dir: Path, central_window_size: int, output_dir: Path | None = None) -> Path:
    processed_dir = Path(processed_dir)
    output_dir = Path(output_dir) if output_dir else processed_dir / "reports"
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = _load_json(processed_dir / "metadata.json")
    sample_index = _load_json(processed_dir / "sample_index.json").get("samples", [])
    gene_vocab = _load_json(processed_dir / "gene_vocab.json")
    idx_to_gene = {int(idx): gene for gene, idx in gene_vocab.items()}
    regions = {region["gene_id"]: region for region in metadata.get("regions", [])}

    central_slices: Dict[str, tuple[int, int]] = {}
    for gene, region in regions.items():
        length = int(region["end"]) - int(region["start"]) + 1
        width = min(int(central_window_size), length)
        start = max((length - width) // 2, 0)
        central_slices[gene] = (start, start + width)

    all_counts_by_gene = {gene: 0 for gene in gene_vocab}
    central_counts_by_gene = {gene: 0 for gene in gene_vocab}
    all_per_sample_by_gene: Dict[str, List[int]] = {gene: [] for gene in gene_vocab}
    central_per_sample_by_gene: Dict[str, List[int]] = {gene: [] for gene in gene_vocab}
    total_tokens_per_sample: List[int] = []
    central_tokens_per_sample: List[int] = []
    sample_rows = []

    for row in sample_index:
        sample_path = processed_dir / row["path"]
        item = torch.load(sample_path, map_location="cpu", weights_only=True)
        genes = item["gene"].to(torch.long)
        positions = item["position_relative"].to(torch.long)
        sample_gene_counts = defaultdict(int)
        sample_central_gene_counts = defaultdict(int)

        for gene_idx, pos in zip(genes.tolist(), positions.tolist()):
            gene = idx_to_gene[int(gene_idx)]
            sample_gene_counts[gene] += 1
            all_counts_by_gene[gene] += 1
            central_start, central_end = central_slices[gene]
            if central_start <= int(pos) < central_end:
                sample_central_gene_counts[gene] += 1
                central_counts_by_gene[gene] += 1

        total_tokens = int(genes.numel())
        central_tokens = int(sum(sample_central_gene_counts.values()))
        total_tokens_per_sample.append(total_tokens)
        central_tokens_per_sample.append(central_tokens)
        sample_rows.append({
            "sample_id": row.get("sample_id"),
            "target": row.get("target"),
            "population": row.get("population"),
            "superpopulation": row.get("superpopulation"),
            "total_tokens": total_tokens,
            f"central_{central_window_size}_tokens": central_tokens,
        })

        for gene in gene_vocab:
            all_per_sample_by_gene[gene].append(int(sample_gene_counts[gene]))
            central_per_sample_by_gene[gene].append(int(sample_central_gene_counts[gene]))

    gene_rows = []
    for gene in sorted(gene_vocab, key=lambda g: gene_vocab[g]):
        region = regions.get(gene, {})
        region_len = int(region.get("end", 0)) - int(region.get("start", 0)) + 1 if region else 0
        all_summary = _summary(all_per_sample_by_gene[gene])
        central_summary = _summary(central_per_sample_by_gene[gene])
        total = int(all_counts_by_gene[gene])
        central = int(central_counts_by_gene[gene])
        gene_rows.append({
            "gene": gene,
            "chrom": region.get("chrom", ""),
            "region_start": region.get("start", ""),
            "region_end": region.get("end", ""),
            "region_len": region_len,
            "central_rel_start": central_slices[gene][0],
            "central_rel_end_exclusive": central_slices[gene][1],
            "all_total_tokens": total,
            "all_mean_per_sample": all_summary["mean"],
            "all_median_per_sample": all_summary["median"],
            "all_p95_per_sample": all_summary["p95"],
            "all_p99_per_sample": all_summary["p99"],
            "all_max_per_sample": all_summary["max"],
            f"central_{central_window_size}_total_tokens": central,
            f"central_{central_window_size}_mean_per_sample": central_summary["mean"],
            f"central_{central_window_size}_median_per_sample": central_summary["median"],
            f"central_{central_window_size}_p95_per_sample": central_summary["p95"],
            f"central_{central_window_size}_p99_per_sample": central_summary["p99"],
            f"central_{central_window_size}_max_per_sample": central_summary["max"],
            "central_fraction": float(central / total) if total else 0.0,
        })

    overall = {
        "processed_dir": str(processed_dir),
        "num_samples": len(sample_index),
        "num_genes": len(gene_vocab),
        "central_window_size": central_window_size,
        "all_total_tokens": int(sum(total_tokens_per_sample)),
        f"central_{central_window_size}_total_tokens": int(sum(central_tokens_per_sample)),
        "all_per_sample": _summary(total_tokens_per_sample),
        f"central_{central_window_size}_per_sample": _summary(central_tokens_per_sample),
        "samples_over_4096_all": int(sum(v > 4096 for v in total_tokens_per_sample)),
        f"samples_over_4096_central_{central_window_size}": int(sum(v > 4096 for v in central_tokens_per_sample)),
    }

    json_path = output_dir / f"variant_counts_central_{central_window_size}.json"
    csv_path = output_dir / f"variant_counts_by_gene_central_{central_window_size}.csv"
    sample_csv_path = output_dir / f"variant_counts_by_sample_central_{central_window_size}.csv"
    md_path = output_dir / f"variant_counts_report_central_{central_window_size}.md"

    with open(json_path, "w") as f:
        json.dump({"overall": overall, "genes": gene_rows, "samples": sample_rows}, f, indent=2)
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(gene_rows[0].keys()))
        writer.writeheader()
        writer.writerows(gene_rows)
    with open(sample_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(sample_rows[0].keys()))
        writer.writeheader()
        writer.writerows(sample_rows)

    with open(md_path, "w") as f:
        f.write(f"# Variant Count Report - central {central_window_size} bp\n\n")
        f.write(f"Dataset: `{processed_dir}`\n\n")
        f.write("## Overall\n\n")
        f.write(f"- Samples: {overall['num_samples']}\n")
        f.write(f"- Genes: {overall['num_genes']}\n")
        f.write(f"- Full-window total tokens: {overall['all_total_tokens']:,}\n")
        f.write(f"- Central {central_window_size} bp total tokens: {overall[f'central_{central_window_size}_total_tokens']:,}\n")
        f.write(f"- Samples > 4096 tokens full-window: {overall['samples_over_4096_all']}\n")
        f.write(f"- Samples > 4096 tokens central {central_window_size} bp: {overall[f'samples_over_4096_central_{central_window_size}']}\n\n")
        f.write("## Per-Gene Counts\n\n")
        f.write("| Gene | Full total | Full mean/sample | Full p99 | Full max | Central total | Central mean/sample | Central p99 | Central max | Central fraction |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in gene_rows:
            f.write(
                f"| {row['gene']} | {int(row['all_total_tokens']):,} | {row['all_mean_per_sample']:.2f} | {row['all_p99_per_sample']:.2f} | {row['all_max_per_sample']:.0f} | "
                f"{int(row[f'central_{central_window_size}_total_tokens']):,} | {row[f'central_{central_window_size}_mean_per_sample']:.2f} | "
                f"{row[f'central_{central_window_size}_p99_per_sample']:.2f} | {row[f'central_{central_window_size}_max_per_sample']:.0f} | {row['central_fraction']:.4f} |\n"
            )

    table = Table(title=f"Variant counts by gene - central {central_window_size} bp")
    table.add_column("Gene")
    table.add_column("Full total", justify="right")
    table.add_column("Full mean", justify="right")
    table.add_column("Central total", justify="right")
    table.add_column("Central mean", justify="right")
    table.add_column("Central frac", justify="right")
    for row in gene_rows:
        table.add_row(
            row["gene"],
            f"{int(row['all_total_tokens']):,}",
            f"{row['all_mean_per_sample']:.1f}",
            f"{int(row[f'central_{central_window_size}_total_tokens']):,}",
            f"{row[f'central_{central_window_size}_mean_per_sample']:.1f}",
            f"{row['central_fraction']:.3f}",
        )
    console.print(table)
    console.print(f"[green]Relatorio salvo:[/green] {md_path}")
    return md_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Conta tokens de variantes por gene e no recorte central")
    parser.add_argument("processed_dir", type=Path)
    parser.add_argument("--central-window-size", type=int, default=32768)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    analyze_variant_counts(args.processed_dir, args.central_window_size, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
