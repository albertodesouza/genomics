"""Cohort genotypes per window and AlphaGenome effects of single variants (the Variant page).

``CohortGenotypes`` holds, for one window, every variant site carried by any sample of the dataset and
the haplotypes carrying its ALT allele (sparse: CSR over sites). It is built once per window from the
per-sample window VCFs and cached on disk, so the genotype of any site across the whole cohort is a
lookup. A sample whose VCF has no record at a site is homozygous reference there (the window VCFs list
non-reference calls only); multi-allelic records are split into one site per ALT allele.

``VariantEffect`` relates a site's genotypes to an AlphaGenome track summarised over a region (mean
signal per haplotype over ``[start, end)`` in reference coordinates), like an in-silico eQTL:

* diploid: the per-sample mean of H1 and H2 regressed on ALT dosage (0/1/2), with a t-test of the slope;
* haplotype: mean of ALT-carrying haplotypes vs REF-carrying haplotypes, as log2 fold change.

Both are associations across the cohort: other variants in linkage disequilibrium ride along with the
ALT allele, exactly as in an eQTL study, which is what makes the comparison with GTEx meaningful.
"""
from __future__ import annotations

import math
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.visualizer.cache import DiskArrayCache, LRUCache, stable_key
from genomics.visualizer.coords import parse_window_vcf
from genomics.visualizer.datasets import Dataset
from genomics.visualizer.jobs import JobCancelled, is_cancelled
from genomics.visualizer.signals import REFERENCE_SAMPLE

GENOTYPE_LABELS = ("0/0", "0/1", "1/1")


@dataclass
class CohortGenotypes:
    samples: List[str]          # columns: haplotype index = 2 * sample index + (0 for H1, 1 for H2)
    positions: np.ndarray       # int64, 1-based genomic POS per site (sorted)
    refs: List[str]
    alts: List[str]
    ids: List[str]
    indptr: np.ndarray          # int64 (n_sites + 1,)
    carriers: np.ndarray        # int32 haplotype indices carrying the site's ALT allele

    @property
    def n_sites(self) -> int:
        return int(self.positions.size)

    def find(self, pos: int, ref: str, alt: str) -> int:
        lo = int(np.searchsorted(self.positions, pos, side="left"))
        hi = int(np.searchsorted(self.positions, pos, side="right"))
        for i in range(lo, hi):
            if self.refs[i] == ref and self.alts[i] == alt:
                return i
        raise KeyError(f"No sample carries {pos} {ref}>{alt} in this window")

    def at(self, pos: int) -> List[int]:
        lo = int(np.searchsorted(self.positions, pos, side="left"))
        hi = int(np.searchsorted(self.positions, pos, side="right"))
        return list(range(lo, hi))

    def haplotype_alleles(self, site: int) -> np.ndarray:
        """(n_samples, 2) int8: 1 where H1/H2 carries the ALT allele."""
        out = np.zeros(2 * len(self.samples), np.int8)
        out[self.carriers[self.indptr[site]:self.indptr[site + 1]]] = 1
        return out.reshape(-1, 2)

    def dosage(self, site: int) -> np.ndarray:
        return self.haplotype_alleles(site).sum(axis=1).astype(np.int8)

    def allele_counts(self) -> np.ndarray:
        return np.diff(self.indptr).astype(np.int64)

    def to_arrays(self) -> Dict[str, np.ndarray]:
        return {
            "samples": np.asarray(self.samples, dtype=str),
            "positions": self.positions,
            "refs": np.asarray(self.refs, dtype=str),
            "alts": np.asarray(self.alts, dtype=str),
            "ids": np.asarray(self.ids, dtype=str),
            "indptr": self.indptr,
            "carriers": self.carriers,
        }

    @classmethod
    def from_arrays(cls, a: Dict[str, np.ndarray]) -> "CohortGenotypes":
        return cls(
            samples=a["samples"].tolist(), positions=a["positions"].astype(np.int64), refs=a["refs"].tolist(),
            alts=a["alts"].tolist(), ids=a["ids"].tolist(), indptr=a["indptr"].astype(np.int64), carriers=a["carriers"].astype(np.int32),
        )


def build_cohort_genotypes(dataset: Dataset, gene: str, samples: Sequence[str], progress=None, workers: int = 8) -> CohortGenotypes:
    window = dataset.window(gene)
    if window.start is None:
        raise ValueError(f"{gene}: the window has no genomic coordinates")
    samples = list(samples)
    calls: Dict[Tuple[int, str, str], List[int]] = {}
    ids: Dict[Tuple[int, str, str], str] = {}

    def parse(index: int):
        path = dataset.sample_vcf_path(samples[index], gene)
        if path is None or is_cancelled(progress):
            return index, None
        return index, parse_window_vcf(path, window.start)

    done = 0
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        for index, parsed in pool.map(parse, range(len(samples))):
            done += 1
            if progress is not None and (done % 64 == 0 or done == len(samples)):
                progress(0.95 * done / max(len(samples), 1), f"Reading genotypes {done}/{len(samples)}")
            if parsed is None:
                continue
            for i in range(parsed.positions.size):
                pos, ref = int(parsed.positions[i]), parsed.refs[i]
                for h in (0, 1):
                    allele = int(parsed.carried[i, h])
                    if allele <= 0 or allele > len(parsed.alt_alleles[i]):
                        continue
                    key = (pos, ref, parsed.alt_alleles[i][allele - 1])
                    calls.setdefault(key, []).append(2 * index + h)
                    ids.setdefault(key, parsed.ids[i])
    if is_cancelled(progress):
        raise JobCancelled()
    keys = sorted(calls)
    lengths = np.fromiter((len(calls[k]) for k in keys), dtype=np.int64, count=len(keys))
    indptr = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    carriers = np.fromiter((h for k in keys for h in sorted(calls[k])), dtype=np.int32, count=int(indptr[-1]))
    return CohortGenotypes(
        samples=samples, positions=np.asarray([k[0] for k in keys], np.int64), refs=[k[1] for k in keys],
        alts=[k[2] for k in keys], ids=[ids[k] for k in keys], indptr=indptr, carriers=carriers,
    )


# ----------------------------------------------------------------------------- statistics
def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction of the regularized incomplete beta (Numerical Recipes, Lentz's method)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        for aa in (m * (b - m) * x / ((qam + m2) * (a + m2)), -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))):
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + aa / c
            c = c if abs(c) > tiny else tiny
            h *= d * c
        if abs(d * c - 1.0) < 1e-12:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_test_p(t: float, df: float) -> float:
    """Two-sided p-value of Student's t with ``df`` degrees of freedom."""
    if not math.isfinite(t) or df <= 0:
        return float("nan")
    return _betainc(df / 2.0, 0.5, df / (df + t * t))


def regress_on_dosage(values: np.ndarray, dosage: np.ndarray) -> Dict[str, float]:
    """OLS of ``values`` on ALT dosage: slope (per ALT allele), its SE, t, two-sided p, r."""
    ok = np.isfinite(values)
    y = values[ok].astype(np.float64)
    x = dosage[ok].astype(np.float64)
    n = int(y.size)
    out = {"n": n, "slope": float("nan"), "se": float("nan"), "t": float("nan"), "p": float("nan"), "r": float("nan"), "intercept": float("nan")}
    if n < 3 or np.ptp(x) == 0:
        return out
    xm, ym = x.mean(), y.mean()
    sxx = float(((x - xm) ** 2).sum())
    sxy = float(((x - xm) * (y - ym)).sum())
    syy = float(((y - ym) ** 2).sum())
    slope = sxy / sxx
    intercept = ym - slope * xm
    rss = max(syy - slope * sxy, 0.0)
    df = n - 2
    se = math.sqrt(rss / df / sxx) if df > 0 else float("nan")
    t = slope / se if se > 0 else (math.copysign(float("inf"), slope) if slope else 0.0)
    p = t_test_p(t, df) if math.isfinite(t) else 0.0
    r = sxy / math.sqrt(sxx * syy) if syy > 0 else float("nan")
    out.update(slope=slope, se=se, t=t, p=p, r=r, intercept=intercept)
    return out


def summarize(values: np.ndarray) -> Dict[str, float]:
    v = values[np.isfinite(values)]
    if not v.size:
        return {"n": 0}
    q = np.quantile(v, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {"n": int(v.size), "mean": float(v.mean()), "sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "min": float(q[0]), "q1": float(q[1]), "median": float(q[2]), "q3": float(q[3]), "max": float(q[4])}


def nanmean_columns(values: np.ndarray) -> np.ndarray:
    """Column means ignoring NaN; NaN (without a warning) for all-NaN columns."""
    ok = np.isfinite(values)
    count = ok.sum(axis=0)
    total = np.where(ok, values, 0.0).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(count > 0, total / np.maximum(count, 1), np.nan).astype(np.float32)


def log2_fold_change(alt: float, ref: float) -> float:
    if not (math.isfinite(alt) and math.isfinite(ref)) or alt <= 0 or ref <= 0:
        return float("nan")
    return math.log2(alt / ref)


# ----------------------------------------------------------------------------- service
class GenotypeService:
    def __init__(self, signals, disk: DiskArrayCache, cache_bytes: int = 512 << 20):
        self.signals = signals
        self.disk = disk
        self.memory = LRUCache(cache_bytes)

    def key(self, dataset: Dataset, gene: str) -> str:
        return stable_key({"v": 1, "dataset": dataset.fingerprint, "path": str(dataset.path), "gene": gene})

    def cached(self, dataset: Dataset, gene: str) -> Optional[CohortGenotypes]:
        key = self.key(dataset, gene)
        hit = self.memory.get(("geno", key))
        if hit is not None:
            return hit
        stored = self.disk.load("cohort_genotypes", key)
        if stored is None:
            return None
        geno = CohortGenotypes.from_arrays(stored)
        self.memory.put(("geno", key), geno)
        return geno

    def genotypes(self, dataset: Dataset, gene: str, progress=None) -> CohortGenotypes:
        hit = self.cached(dataset, gene)
        if hit is not None:
            return hit
        samples = [row["sample_id"] for row in dataset.samples]
        workers = min(16, os.cpu_count() or 4)
        geno = build_cohort_genotypes(dataset, gene, samples, progress, workers=workers)
        key = self.key(dataset, gene)
        self.disk.save("cohort_genotypes", key, geno.to_arrays(), compress=True)
        self.memory.put(("geno", key), geno)
        return geno

    # -- region means per haplotype -------------------------------------------------
    def region_key(self, dataset: Dataset, gene: str, output: str, start: int, end: int) -> str:
        return stable_key({"v": 1, "dataset": dataset.fingerprint, "path": str(dataset.path), "gene": gene, "output": output, "start": start, "end": end})

    def cached_region_means(self, dataset: Dataset, gene: str, output: str, start: int, end: int) -> Optional[Dict[str, np.ndarray]]:
        key = self.region_key(dataset, gene, output, start, end)
        hit = self.memory.get(("region", key))
        if hit is None:
            hit = self.disk.load("region_means", key)
            if hit is not None:
                self.memory.put(("region", key), hit)
        return hit

    def region_means(self, dataset: Dataset, gene: str, output: str, start: int, end: int, progress=None) -> Dict[str, np.ndarray]:
        """``means``: (n_samples, 2, tracks) mean signal per haplotype over reference offsets [start, end)
        (NaN where a sample has no prediction), plus ``samples`` and ``reference`` (tracks,)."""
        hit = self.cached_region_means(dataset, gene, output, start, end)
        if hit is not None:
            return hit
        samples = [row["sample_id"] for row in dataset.samples]
        if progress is not None:
            self.signals.ensure_variants(dataset, gene, samples, progress.sub(0.0, 0.2, "Indel maps: ") if hasattr(progress, "sub") else progress)
        n_tracks = None
        jobs = [(i, s, h) for i, s in enumerate(samples) for h in (0, 1)]
        results: Dict[Tuple[int, int], np.ndarray] = {}

        def one(job):
            i, sample, h = job
            if is_cancelled(progress):
                return job, None
            path = dataset.prediction_path(sample, gene, ("H1", "H2")[h], output)
            if not path.exists():
                return job, None
            try:
                values = self.signals.haplotype_window(dataset, sample, gene, ("H1", "H2")[h], output, "reference", start, end, None, cache=False)
            except Exception:
                return job, None
            return job, nanmean_columns(values)

        done = 0
        with ThreadPoolExecutor(max_workers=self.signals.bulk_workers) as pool:
            for job, mean in pool.map(one, jobs):
                done += 1
                if progress is not None and (done % 128 == 0 or done == len(jobs)):
                    progress(0.2 + 0.8 * done / len(jobs), f"Region means {done}/{len(jobs)} haplotypes")
                if mean is not None:
                    results[(job[0], job[2])] = mean
                    n_tracks = mean.size
        if is_cancelled(progress):
            raise JobCancelled()
        if n_tracks is None:
            raise FileNotFoundError(f"No {output} predictions for {gene}")
        means = np.full((len(samples), 2, n_tracks), np.nan, np.float32)
        for (i, h), mean in results.items():
            means[i, h] = mean
        reference = np.full(n_tracks, np.nan, np.float32)
        try:
            ref = self.signals.haplotype_window(dataset, REFERENCE_SAMPLE, gene, "H1", output, "reference", start, end, None)
            reference = nanmean_columns(ref)
        except Exception:
            pass
        arrays = {"samples": np.asarray(samples, dtype=str), "means": means, "reference": reference}
        key = self.region_key(dataset, gene, output, start, end)
        self.disk.save("region_means", key, arrays, compress=True)
        self.memory.put(("region", key), arrays)
        return arrays

    # -- one site -------------------------------------------------------------------
    @staticmethod
    def site_record(geno: CohortGenotypes, site: int, chromosome: Optional[str]) -> Dict[str, object]:
        chrom = chromosome or ""
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        pos, ref, alt = int(geno.positions[site]), geno.refs[site], geno.alts[site]
        return {
            "pos": pos, "ref": ref, "alt": alt, "id": geno.ids[site], "chromosome": chrom,
            "variant_id": f"{chrom}_{pos}_{ref}_{alt}_b38",
            "type": "SNV" if len(ref) == 1 and len(alt) == 1 else ("insertion" if len(alt) > len(ref) else "deletion" if len(alt) < len(ref) else "MNV"),
        }

    def site_summary(self, dataset: Dataset, geno: CohortGenotypes, site: int, cohort: List[str], field: Optional[str]) -> Dict[str, object]:
        index = {s: i for i, s in enumerate(geno.samples)}
        rows = [index[s] for s in cohort if s in index]
        dose = geno.dosage(site)[rows]
        counts = {label: int((dose == k).sum()) for k, label in enumerate(GENOTYPE_LABELS)}
        n = len(rows)
        af = float(dose.sum()) / (2 * n) if n else float("nan")
        by_group = []
        if field:
            groups: Dict[str, List[int]] = {}
            for s, r in zip([s for s in cohort if s in index], rows):
                value = str(dataset.samples[dataset.sample_index[s]].get(field, "") or "")
                if value:
                    groups.setdefault(value, []).append(r)
            all_dose = geno.dosage(site)
            for value, members in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
                d = all_dose[members]
                by_group.append({"group": value, "n": len(members), "af": float(d.sum()) / (2 * len(members)),
                                 "counts": {label: int((d == k).sum()) for k, label in enumerate(GENOTYPE_LABELS)}})
        return {"cohort": n, "counts": counts, "af": af, "field": field, "by_group": by_group}

    def effect(self, geno: CohortGenotypes, site: int, region: Dict[str, np.ndarray], track: int, cohort: List[str], max_points: int = 4000) -> Dict[str, object]:
        index = {s: i for i, s in enumerate(region["samples"].tolist())}
        gindex = {s: i for i, s in enumerate(geno.samples)}
        members = [s for s in cohort if s in index and s in gindex]
        alleles = geno.haplotype_alleles(site)[[gindex[s] for s in members]]  # (n, 2)
        hap_values = region["means"][[index[s] for s in members], :, track].astype(np.float64)  # (n, 2)
        diploid = nanmean_columns(hap_values.T).astype(np.float64)
        dose = alleles.sum(axis=1)
        reg = regress_on_dosage(diploid, dose)
        groups = []
        for k, label in enumerate(GENOTYPE_LABELS):
            sel = dose == k
            vals = diploid[sel]
            ids = [s for s, m in zip(members, sel) if m]
            keep = np.isfinite(vals)
            vals, ids = vals[keep], [s for s, kk in zip(ids, keep) if kk]
            if vals.size > max_points:  # thin for the strip plot; summaries use every value
                pick = np.linspace(0, vals.size - 1, max_points).round().astype(int)
                shown, shown_ids = vals[pick], [ids[i] for i in pick]
            else:
                shown, shown_ids = vals, ids
            groups.append({"genotype": label, "summary": summarize(diploid[sel]), "values": shown.astype(np.float32).tolist(), "samples": shown_ids})
        ref_haps = hap_values[alleles == 0]
        alt_haps = hap_values[alleles == 1]
        ref_mean = float(nanmean_columns(ref_haps.reshape(-1, 1))[0])
        alt_mean = float(nanmean_columns(alt_haps.reshape(-1, 1))[0])
        reference = float(region["reference"][track]) if region["reference"].size > track else float("nan")
        cohort_means = region["means"][[index[s] for s in members]].reshape(-1, region["means"].shape[2])
        return {
            "track_means": [None if not np.isfinite(v) else float(v) for v in nanmean_columns(cohort_means)],
            "n_samples": len(members),
            "groups": groups,
            "regression": reg,
            "haplotypes": {"ref": summarize(ref_haps), "alt": summarize(alt_haps), "log2fc": log2_fold_change(alt_mean, ref_mean)},
            "reference_genome": reference,
            "direction": 0 if not math.isfinite(reg["slope"]) or reg["slope"] == 0 else (1 if reg["slope"] > 0 else -1),
        }
