"""Gene lookup for adding windows: every ``gene`` row of a GTF cache table, searchable by symbol or id.

The table is the dataset's ``gtf_cache.feather`` (or the 1000 Genomes one); loading keeps only
gene rows and a few columns, once per table.
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

COLUMNS = ["Chromosome", "Feature", "Start", "End", "Strand", "gene_name", "gene_id", "gene_type"]
# Searches rank protein-coding genes, then lncRNAs, before pseudogenes and small RNAs.
TYPE_RANK = {"protein_coding": 0, "lncRNA": 1}


class GeneIndex:
    def __init__(self) -> None:
        self._tables: Dict[str, List[Dict[str, Any]]] = {}
        self._lookups: Dict[str, Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]] = {}
        self._lock = threading.Lock()

    def genes(self, table: Path) -> List[Dict[str, Any]]:
        key = str(Path(table).resolve())
        with self._lock:
            cached = self._tables.get(key)
            if cached is not None:
                return cached
            import pandas as pd

            path = Path(table)
            if path.suffix in (".feather", ".parquet"):
                read = pd.read_feather if path.suffix == ".feather" else pd.read_parquet
                try:
                    df = read(path, columns=COLUMNS)
                except (KeyError, ValueError):  # no gene_type column
                    df = read(path, columns=COLUMNS[:-1])
            else:
                df = pd.read_csv(path, sep=None, engine="python")
            df = df[df["Feature"] == "gene"]
            rows = []
            for row in df.itertuples(index=False):
                chrom = str(row.Chromosome)
                rows.append({
                    "name": str(row.gene_name),
                    "id": str(row.gene_id).split(".")[0],
                    "type": str(getattr(row, "gene_type", "") or ""),
                    "chrom": chrom,
                    "start": int(row.Start) + 1,  # GTF tables from pyranges are 0-based
                    "end": int(row.End),
                    "strand": str(row.Strand),
                    "primary": "_" not in chrom,
                })
            self._tables[key] = rows
            return rows

    def lookup(self, table: Path) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:
        """(by gene name, by Ensembl gene id) over the table's genes, primary assembly first."""
        key = f"lookup:{Path(table).resolve()}"
        with self._lock:
            cached = self._lookups.get(key)
        if cached is not None:
            return cached
        by_name: Dict[str, Dict[str, Any]] = {}
        by_id: Dict[str, Dict[str, Any]] = {}
        for gene in sorted(self.genes(table), key=lambda g: not g["primary"]):
            by_name.setdefault(gene["name"].upper(), gene)
            by_id.setdefault(gene["id"].upper(), gene)
        with self._lock:
            self._lookups[key] = (by_name, by_id)
        return by_name, by_id

    def locate(self, table: Path, symbol: str, ensembl_id: str = "") -> Optional[Dict[str, Any]]:
        """GENCODE row of a gene, by Ensembl id first (symbols change between releases), then by name."""
        by_name, by_id = self.lookup(table)
        return (by_id.get(ensembl_id.upper()) if ensembl_id else None) or by_name.get(symbol.upper())

    def search(self, table: Path, query: str, limit: int = 25) -> List[Dict[str, Any]]:
        q = query.strip().upper()
        if not q:
            return []
        hits = []
        for gene in self.genes(table):
            name = gene["name"].upper()
            if name == q or gene["id"].upper() == q.split(".")[0]:
                rank = 0
            elif name.startswith(q):
                rank = 1
            elif q in name:
                rank = 2
            else:
                continue
            hits.append((rank, not gene["primary"], TYPE_RANK.get(gene["type"], 2), len(name), name, gene))
        hits.sort(key=lambda h: h[:5])
        return [h[-1] for h in hits[:limit]]


def find_table(dataset_path: Optional[Path], explicit: Optional[Path] = None) -> Optional[Path]:
    candidates = [Path(explicit)] if explicit else []
    if dataset_path is not None:
        candidates += [Path(dataset_path) / "gtf_cache.feather", Path(dataset_path) / "gtf_cache.parquet"]
    from genomics.visualizer.launch import default_gtf

    fallback = default_gtf()
    if fallback:
        candidates.append(Path(fallback))
    return next((c for c in candidates if c.exists()), None)
