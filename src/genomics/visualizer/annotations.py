"""Gene models (genes, transcripts, exons, CDS) overlapping each dataset window.

Read from a GTF cache table (``gtf_cache.feather`` in the dataset directory, or any
``--gtf`` feather/parquet/csv table with GTF columns). Requires pandas (+ pyarrow for feather);
without it the browser simply shows no gene-model lane. Results are cached on disk per dataset.
"""
from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from genomics.visualizer.cache import stable_key
from genomics.visualizer.datasets import Dataset

FEATURES = ("gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon")


def find_gtf_table(dataset: Dataset, explicit: Optional[Path] = None) -> Optional[Path]:
    if explicit and Path(explicit).exists():
        return Path(explicit)
    for name in ("gtf_cache.feather", "gtf_cache.parquet"):
        candidate = dataset.path / name
        if candidate.exists():
            return candidate
    return None


class AnnotationService:
    def __init__(self, cache_dir: Optional[Path], gtf: Optional[Path] = None):
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.gtf = gtf
        self._memory: Dict[str, Dict[str, Any]] = {}
        self._lock = threading.Lock()

    def _cache_path(self, dataset: Dataset, table: Path) -> Optional[Path]:
        if self.cache_dir is None:
            return None
        stat = table.stat()
        key = stable_key({"v": 2, "dataset": dataset.fingerprint, "table": str(table), "mtime": stat.st_mtime_ns, "genes": dataset.genes})
        return self.cache_dir / "annotations" / f"{key}.json"

    def status(self, dataset: Dataset) -> Dict[str, Any]:
        table = find_gtf_table(dataset, self.gtf)
        return {"available": table is not None, "source": str(table) if table else None, "loaded": dataset.id in self._memory}

    def cached(self, dataset: Dataset) -> Optional[Dict[str, Any]]:
        hit = self._memory.get(dataset.id)
        if hit is not None:
            return hit
        table = find_gtf_table(dataset, self.gtf)
        if table is None:
            return {"source": None, "genes": {}}
        path = self._cache_path(dataset, table)
        if path is not None and path.exists():
            try:
                with open(path, "r", encoding="utf-8") as f:
                    payload = json.load(f)
                self._memory[dataset.id] = payload
                return payload
            except Exception:
                return None
        return None

    def build(self, dataset: Dataset, progress: Callable[[float, str], None]) -> Dict[str, Any]:
        with self._lock:
            hit = self.cached(dataset)
            if hit is not None:
                return hit
            table = find_gtf_table(dataset, self.gtf)
            if table is None:
                return {"source": None, "genes": {}}
            try:
                import pandas as pd
            except ImportError as exc:
                raise RuntimeError("pandas is required to read GTF annotations") from exc
            progress(0.05, f"Reading {table.name}")
            columns = ["Chromosome", "Feature", "Start", "End", "Strand", "gene_name", "gene_type", "transcript_id", "transcript_name", "transcript_type", "tag", "exon_number"]
            if table.suffix == ".feather":
                df = pd.read_feather(table, columns=columns)
            elif table.suffix == ".parquet":
                df = pd.read_parquet(table, columns=columns)
            else:
                df = pd.read_csv(table, usecols=columns, sep=None, engine="python")
            df = df[df["Feature"].isin(FEATURES)]
            progress(0.6, "Selecting features in gene windows")
            by_gene: Dict[str, List[Dict[str, Any]]] = {}
            for gene in dataset.genes:
                window = dataset.window(gene)
                if window.chromosome is None or window.start is None or window.end is None:
                    continue
                chrom = window.chromosome
                alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                sel = df[df["Chromosome"].astype(str).isin([chrom, alt_chrom]) & (df["End"] >= window.start) & (df["Start"] <= window.end)]
                by_gene[gene] = _features_to_models(sel, window.start)
            payload = {"source": str(table), "genes": by_gene}
            path = self._cache_path(dataset, table)
            if path is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
                with open(path, "w", encoding="utf-8") as f:
                    json.dump(payload, f, separators=(",", ":"))
            self._memory[dataset.id] = payload
            progress(1.0, "Annotations ready")
            return payload


def _features_to_models(df, window_start: int) -> List[Dict[str, Any]]:
    """Group GTF rows into gene -> transcripts -> exon/CDS blocks (window offsets, 0-based half-open).

    GTF cache tables (pyranges) are 0-based half-open while ``window_start`` is the 1-based
    position of window offset 0, so offset = Start - (window_start - 1).
    """
    genes: Dict[str, Dict[str, Any]] = {}
    transcripts: Dict[str, Dict[str, Any]] = {}
    for row in df.itertuples(index=False):
        feature = row.Feature
        start = int(row.Start) - window_start + 1
        end = int(row.End) - window_start + 1
        if feature == "gene":
            genes[str(row.gene_name)] = {"name": str(row.gene_name), "type": str(row.gene_type), "strand": str(row.Strand), "start": start, "end": end, "transcripts": []}
            continue
        tid = str(row.transcript_id) if row.transcript_id is not None else None
        if not tid or tid == "None":
            continue
        tx = transcripts.get(tid)
        if tx is None:
            tags = str(row.tag or "")
            tx = transcripts[tid] = {
                "id": tid,
                "name": str(row.transcript_name or tid),
                "gene": str(row.gene_name),
                "type": str(row.transcript_type or ""),
                "strand": str(row.Strand),
                "canonical": "Ensembl_canonical" in tags,
                "mane": "MANE_Select" in tags,
                "start": start,
                "end": end,
                "exons": [],
                "cds": [],
            }
        if feature == "transcript":
            tx["start"], tx["end"] = start, end
        elif feature == "exon":
            tx["exons"].append([start, end])
        elif feature == "CDS":
            tx["cds"].append([start, end])
    for tx in transcripts.values():
        gene = genes.setdefault(tx["gene"], {"name": tx["gene"], "type": tx["type"], "strand": tx["strand"], "start": tx["start"], "end": tx["end"], "transcripts": []})
        tx["exons"].sort()
        tx["cds"].sort()
        gene["transcripts"].append(tx)
    out = []
    for gene in genes.values():
        gene["transcripts"].sort(key=lambda t: (not t["mane"], not t["canonical"], t["start"]))
        out.append(gene)
    out.sort(key=lambda g: g["start"])
    return out
