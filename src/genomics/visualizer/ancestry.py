"""Genotype PCA over the dataset's windows, and cohort matching on the principal components.

PCA: biallelic SNVs with minor allele frequency >= ``min_maf`` over all samples, thinned to one site per
``spacing`` bp (a crude stand-in for LD pruning), from the cohort genotype matrices of the chosen windows
(``genotypes.CohortGenotypes``, built once per window and cached). Dosages are standardised per site,
``(g - 2p) / sqrt(2p(1 - p))``, and the sample scores come from the eigendecomposition of the n x n
genetic relationship matrix, as in EIGENSOFT. Only samples genotyped in every chosen window get PCs. By
default the windows covering the most samples are used, preferring the control windows (not listed in
the dataset metadata) so the pigmentation loci under selection shape the PCs only when nothing else
covers the cohort.

Matching: each sample of group A is paired with its nearest unused sample of group B in the space of the
first ``k`` PCs (greedy 1:1 nearest neighbour, A samples taken in order of their nearest distance, within a
caliper). Balance is reported as the standardised mean difference of each PC before and after matching.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from genomics.visualizer.cache import DiskArrayCache, LRUCache, stable_key
from genomics.visualizer.datasets import Dataset
from genomics.visualizer.jobs import JobCancelled, is_cancelled

MAX_COMPONENTS = 20


def window_coverage(dataset: Dataset, probes: int = 128) -> Dict[str, float]:
    """Estimated fraction of samples with a VCF per window (from an evenly spaced subset of samples)."""
    rows = dataset.samples
    if not rows:
        return {g: 0.0 for g in dataset.genes}
    step = max(1, len(rows) // probes)
    picked = [r["sample_id"] for r in rows[::step]]
    return {g: sum(dataset.sample_vcf_path(s, g) is not None for s in picked) / len(picked) for g in dataset.genes}


def default_windows(dataset: Dataset, coverage: Optional[Dict[str, float]] = None) -> List[str]:
    """Windows covering the most samples, preferring the control windows (not listed in the metadata),
    so pigmentation loci under selection shape the PCs only when nothing else covers the cohort."""
    coverage = coverage if coverage is not None else window_coverage(dataset)
    if not coverage:
        return list(dataset.genes)
    best = max(coverage.values())
    widest = [g for g in dataset.genes if coverage.get(g, 0.0) >= best - 1e-9]
    listed = {str(g) for g in dataset.metadata.get("genes") or []}
    controls = [g for g in widest if g not in listed]
    return controls or widest


def select_sites(geno, min_maf: float, spacing: int) -> List[int]:
    """Biallelic SNV sites with MAF >= min_maf, at least ``spacing`` bp apart (first come in position order)."""
    n_haps = 2 * len(geno.samples)
    counts = np.diff(geno.indptr)
    positions = geno.positions
    multi = np.zeros(geno.n_sites, bool)
    same = positions[1:] == positions[:-1]  # several ALT alleles at one position: not biallelic
    multi[1:] |= same
    multi[:-1] |= same
    keep: List[int] = []
    last = -(10 ** 12)
    for i in range(geno.n_sites):
        if multi[i] or len(geno.refs[i]) != 1 or len(geno.alts[i]) != 1:
            continue
        af = counts[i] / n_haps
        if min(af, 1.0 - af) < min_maf or positions[i] - last < spacing:
            continue
        keep.append(i)
        last = int(positions[i])
    return keep


def dosage_matrix(geno, sites: Sequence[int]) -> np.ndarray:
    """(n_samples, len(sites)) ALT dosages as float32."""
    n = len(geno.samples)
    out = np.zeros((n, len(sites)), np.float32)
    for j, site in enumerate(sites):
        haps = geno.carriers[geno.indptr[site]:geno.indptr[site + 1]]
        np.add.at(out[:, j], haps // 2, 1.0)
    return out


def pca_scores(dosages: np.ndarray, components: int) -> Dict[str, np.ndarray]:
    """Sample scores (n, k) and explained variance ratios (k,) of standardised dosages (n, m)."""
    n, m = dosages.shape
    if m < 2 or n < 3:
        raise ValueError(f"PCA needs at least 3 samples and 2 sites (got {n} samples, {m} sites)")
    p = dosages.mean(axis=0) / 2.0
    scale = np.sqrt(2.0 * p * (1.0 - p))
    ok = scale > 0
    x = ((dosages[:, ok] - 2.0 * p[ok]) / scale[ok]).astype(np.float64)
    grm = x @ x.T / x.shape[1]
    values, vectors = np.linalg.eigh(grm)
    order = np.argsort(values)[::-1][: min(components, n - 1)]
    values, vectors = np.clip(values[order], 0.0, None), vectors[:, order]
    # Deterministic signs: the sample with the largest |score| on each PC is positive.
    signs = np.sign(vectors[np.abs(vectors).argmax(axis=0), np.arange(vectors.shape[1])])
    signs[signs == 0] = 1.0
    scores = vectors * signs * np.sqrt(values)
    total = float(np.trace(grm))
    return {"scores": scores.astype(np.float32), "explained": (values / total if total > 0 else values).astype(np.float32), "n_sites": np.int64(x.shape[1])}


class AncestryService:
    def __init__(self, genotypes, disk: DiskArrayCache, cache_bytes: int = 64 << 20):
        self.genotypes = genotypes
        self.disk = disk
        self.memory = LRUCache(cache_bytes)

    @staticmethod
    def params(dataset: Dataset, genes: Optional[Sequence[str]], min_maf: float, spacing: int, components: int) -> Dict[str, Any]:
        genes = sorted(set(genes or default_windows(dataset)))
        unknown = [g for g in genes if g not in dataset.genes]
        if unknown:
            raise KeyError(f"Unknown window(s): {', '.join(unknown)}")
        if not 0.0 < min_maf < 0.5:
            raise ValueError("min_maf must be in (0, 0.5)")
        return {"genes": genes, "min_maf": round(float(min_maf), 4), "spacing": max(0, int(spacing)), "components": max(2, min(MAX_COMPONENTS, int(components)))}

    def key(self, dataset: Dataset, params: Dict[str, Any]) -> str:
        # The genotype keys version the inputs: a change in how windows are genotyped invalidates the PCA.
        inputs = [self.genotypes.key(dataset, g) for g in params["genes"]]
        return stable_key({"v": 1, "dataset": dataset.fingerprint, "path": str(dataset.path), "inputs": inputs, **params})

    def cached(self, dataset: Dataset, params: Dict[str, Any]) -> Optional[Dict[str, np.ndarray]]:
        key = self.key(dataset, params)
        hit = self.memory.get(key)
        if hit is None:
            hit = self.disk.load("ancestry_pca", key)
            if hit is not None:
                self.memory.put(key, hit)
        return hit

    def compute(self, dataset: Dataset, params: Dict[str, Any], progress=None) -> Dict[str, np.ndarray]:
        hit = self.cached(dataset, params)
        if hit is not None:
            return hit
        genes = params["genes"]
        matrices = []
        for i, gene in enumerate(genes):
            if is_cancelled(progress):
                raise JobCancelled()
            sub = progress.sub(0.9 * i / len(genes), 0.9 / len(genes), f"{gene} ({i + 1}/{len(genes)}): ") if progress is not None else None
            matrices.append(self.genotypes.genotypes(dataset, gene, sub))
        # Samples genotyped in every window (windows may cover different subsets of the cohort).
        common = set(matrices[0].samples)
        for geno in matrices[1:]:
            common &= set(geno.samples)
        samples = [s for s in matrices[0].samples if s in common]
        if len(samples) < 3:
            raise ValueError("Fewer than 3 samples are genotyped in all chosen windows")
        blocks, per_gene = [], []
        for geno in matrices:
            sites = select_sites(geno, params["min_maf"], params["spacing"])
            per_gene.append(len(sites))
            if sites:
                index = {s: k for k, s in enumerate(geno.samples)}
                blocks.append(dosage_matrix(geno, sites)[[index[s] for s in samples]])
        if progress is not None:
            progress(0.92, f"PCA over {sum(per_gene):,} sites")
        if not blocks:
            raise ValueError("No sites pass the filters")
        result = pca_scores(np.concatenate(blocks, axis=1), params["components"])
        arrays = {"samples": np.asarray(samples, dtype=str), "sites_per_window": np.asarray(per_gene, np.int64), **result}
        key = self.key(dataset, params)
        self.disk.save("ancestry_pca", key, arrays, compress=True)
        self.memory.put(key, arrays)
        return arrays

    @staticmethod
    def payload(result: Dict[str, np.ndarray], params: Dict[str, Any]) -> Dict[str, Any]:
        scores = result["scores"]
        return {
            **params,
            "samples": result["samples"].tolist(),
            "scores": np.round(scores.astype(np.float64), 6).tolist(),
            "explained": [float(v) for v in result["explained"]],
            "n_sites": int(result["n_sites"]),
            "sites_per_window": dict(zip(params["genes"], (int(v) for v in result["sites_per_window"]))),
        }


def smd(a: np.ndarray, b: np.ndarray) -> List[float]:
    """Standardised mean difference per column (pooled SD)."""
    if not len(a) or not len(b):
        return [float("nan")] * (a.shape[1] if a.ndim == 2 else 0)
    pooled = np.sqrt((a.var(axis=0, ddof=1 if len(a) > 1 else 0) + b.var(axis=0, ddof=1 if len(b) > 1 else 0)) / 2.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = (a.mean(axis=0) - b.mean(axis=0)) / pooled
    return [float(v) for v in d]


def match_groups(scores: np.ndarray, group_a: Sequence[int], group_b: Sequence[int], k: int, caliper: float) -> Dict[str, Any]:
    """Greedy 1:1 nearest-neighbour matching of rows ``group_a`` to ``group_b`` on the first ``k`` columns.

    ``caliper`` is in units of the pooled SD of the PC-space distance scale (sqrt of the summed variances of
    the k PCs over A and B); pairs farther apart are left unmatched. Returns index pairs and balance.
    """
    k = max(1, min(k, scores.shape[1]))
    a_idx, b_idx = np.asarray(group_a, int), np.asarray(group_b, int)
    if not a_idx.size or not b_idx.size:
        raise ValueError("Both groups need samples")
    xa, xb = scores[a_idx, :k].astype(np.float64), scores[b_idx, :k].astype(np.float64)
    pooled = np.concatenate([xa, xb])
    limit = caliper * math.sqrt(float(pooled.var(axis=0).sum())) if caliper > 0 else math.inf
    sq = (xa ** 2).sum(axis=1)[:, None] + (xb ** 2).sum(axis=1)[None, :] - 2.0 * xa @ xb.T
    dist = np.sqrt(np.maximum(sq, 0.0))  # (A, B); groups are at most a few thousand samples
    order = np.argsort(dist.min(axis=1), kind="stable")
    used = np.zeros(len(b_idx), bool)
    pairs = []
    for i in order:
        row = np.where(used, np.inf, dist[i])
        j = int(row.argmin())
        if not math.isfinite(row[j]) or row[j] > limit:
            continue
        used[j] = True
        pairs.append((int(a_idx[i]), int(b_idx[j]), float(row[j])))
    ma = np.asarray([p[0] for p in pairs], int)
    mb = np.asarray([p[1] for p in pairs], int)
    return {
        "pairs": pairs, "k": k, "caliper_distance": limit if math.isfinite(limit) else None,
        "balance_before": smd(scores[a_idx, :k], scores[b_idx, :k]),
        "balance_after": smd(scores[ma, :k], scores[mb, :k]) if len(pairs) else [float("nan")] * k,
    }
