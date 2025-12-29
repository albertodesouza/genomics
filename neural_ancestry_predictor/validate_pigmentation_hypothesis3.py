#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
validate_pigmentation_hypothesis3.py

Third, independent validation for the pigmentation hypothesis using external evidence:
- GWAS Catalog (EBI) associations/traits by rsID
- Open Targets Platform (GraphQL) variant/disease/phenotype associations (best-effort)
- Fallback: Ensembl Variation phenotypes by rsID

Input directory is fixed to:
  /home/lume2/genomics/neural_ancestry_predictor/top_regions_reports_central2

Outputs:
  pigmentation_validation3.md (English) in that same directory
  validation3_cache/ (per-rsID JSON cache)

Notes:
- This script is intentionally independent from validate_pigmentation_hypothesis.py(2).
- Requires: requests (pip install requests)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import threading
import time
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import requests


DEFAULT_INPUT_DIR = Path("/home/lume2/genomics/neural_ancestry_predictor/top_regions_reports_central2")
# Bump when changing how we fetch/parse external evidence so existing caches are invalidated.
CACHE_SCHEMA_VERSION = 3


RSID_RE = re.compile(r"\brs\d+\b")


PIGMENTATION_KEYWORDS = [
    # core
    "pigmentation",
    "skin color",
    "skin pigmentation",
    "melanin",
    "melanogenesis",
    "melanocyte",
    "tanning",
    "sunburn",
    "uv",
    "ultraviolet",
    # related phenotypes
    "hair color",
    "eye color",
    "iris",
    "freckles",
    "albinism",
    "oculocutaneous",
]


def _now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _stable_unique(seq: Iterable[str]) -> List[str]:
    seen: Dict[str, None] = {}
    out: List[str] = []
    for x in seq:
        if x not in seen:
            seen[x] = None
            out.append(x)
    return out


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if not r:
                continue
            rows.append({k: (v or "").strip() for k, v in r.items()})
    return rows


def _sleep_rate_limit(seconds: float) -> None:
    if seconds <= 0:
        return
    time.sleep(seconds)


def _request_json(
    method: str,
    url: str,
    *,
    headers: Optional[Dict[str, str]] = None,
    params: Optional[Dict[str, str]] = None,
    json_body: Any = None,
    timeout_s: int | Tuple[int, int] = 30,
    max_retries: int = 6,
    backoff_s: float = 1.0,
    rate_limit_s: float = 0.34,
) -> Tuple[Optional[dict], Optional[str]]:
    """
    Returns (json_dict, error_str).
    Handles 429/5xx with retry+backoff.
    """
    last_err: Optional[str] = None
    for attempt in range(max_retries):
        _sleep_rate_limit(rate_limit_s)
        try:
            r = requests.request(
                method,
                url,
                headers=headers,
                params=params,
                json=json_body,
                timeout=timeout_s,
            )
            if r.status_code in (429, 500, 502, 503, 504):
                last_err = f"HTTP {r.status_code}"
                time.sleep(backoff_s * (2 ** attempt))
                continue
            if not r.ok:
                text = (r.text or "")[:500]
                return None, f"HTTP {r.status_code}: {text}"
            try:
                return r.json(), None
            except Exception as e:  # noqa: BLE001
                return None, f"JSON decode error: {e}"
        except Exception as e:  # noqa: BLE001
            last_err = str(e)
            time.sleep(backoff_s * (2 ** attempt))
    return None, f"Request failed after retries: {last_err}"


@dataclass
class VariantMeta:
    rsid: str
    chrom: str
    pos_1based: int
    min_distance_to_center: int
    worst_impact: str
    most_severe_consequence: str
    gene_symbols: List[str]


@dataclass
class Evidence:
    rsid: str
    fetched_at: str
    cache_schema_version: int
    gwascatalog: Dict[str, Any]
    opentargets: Dict[str, Any]
    ensembl: Dict[str, Any]
    errors: Dict[str, str]


def _cache_path_for_rsid(rsid: str) -> Path:
    safe = rsid.replace("/", "_")
    return CACHE_DIR / f"{safe}.json"


def load_cached_evidence(rsid: str) -> Optional[Evidence]:
    p = _cache_path_for_rsid(rsid)
    if not p.exists():
        return None
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        if int(data.get("cache_schema_version", 0) or 0) != CACHE_SCHEMA_VERSION:
            return None
        return Evidence(
            rsid=data.get("rsid", rsid),
            fetched_at=data.get("fetched_at", ""),
            cache_schema_version=int(data.get("cache_schema_version", 0) or 0),
            gwascatalog=data.get("gwascatalog", {}) or {},
            opentargets=data.get("opentargets", {}) or {},
            ensembl=data.get("ensembl", {}) or {},
            errors=data.get("errors", {}) or {},
        )
    except Exception:
        return None


def save_cached_evidence(ev: Evidence) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    p = _cache_path_for_rsid(ev.rsid)
    # Atomic write to avoid corrupt cache on interruption.
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(asdict(ev), indent=2, sort_keys=True), encoding="utf-8")
    os.replace(tmp, p)


def fetch_ensembl_variation(rsid: str) -> Tuple[Dict[str, Any], Optional[str]]:
    # Ensembl Variation may omit phenotypes unless explicitly requested.
    url = f"https://rest.ensembl.org/variation/human/{rsid}"
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    data, err = _request_json("GET", url, headers=headers, params={"phenotypes": "1"}, timeout_s=(6, 20))
    if err or not (data or {}).get("phenotypes"):
        # Fallback to semicolon-style flags used elsewhere in the repo
        url2 = f"https://rest.ensembl.org/variation/human/{rsid}?pops=1;phenotypes=1"
        data2, err2 = _request_json("GET", url2, headers=headers, timeout_s=(6, 20))
        if err2:
            return {}, err or err2
        data = data2
    phenos = data.get("phenotypes") or []
    trait_list: List[str] = []
    for p in phenos:
        trait = (p.get("trait") or p.get("description") or "").strip()
        if trait:
            trait_list.append(trait)
    out = {
        "phenotypes": _stable_unique(trait_list),
        "clinical_significance": data.get("clinical_significance") or [],
        "MAF": data.get("MAF"),
        "minor_allele": data.get("minor_allele"),
    }
    return out, None


def fetch_gwas_catalog(rsid: str) -> Tuple[Dict[str, Any], Optional[str]]:
    """
    Best-effort GWAS Catalog fetch via the public REST API (HAL).
    We try:
      - SNP endpoint: /singleNucleotidePolymorphisms/{rsid}
      - Follow associations link (first pages only)
    """
    base = "https://www.ebi.ac.uk/gwas/rest/api"
    snp_url = f"{base}/singleNucleotidePolymorphisms/{rsid}"
    headers = {"Accept": "application/json"}
    snp, err = _request_json("GET", snp_url, headers=headers, timeout_s=(6, 30), rate_limit_s=0.34)
    if err:
        return {}, err

    assoc_traits: List[str] = []
    n_associations = 0

    links = (snp.get("_links") or {})
    assoc_href = ((links.get("associations") or {}).get("href") or "").strip()

    # Association pages can be large; cap to a few pages for runtime safety.
    next_url = assoc_href
    pages = 0
    max_pages = 3
    
    # IMPORTANT: GWAS Catalog often does NOT embed traits inside the association object.
    # Traits are usually reachable via association links:
    #   - _links.study.href -> study.diseaseTrait.trait
    #   - _links.efoTraits.href -> _embedded.efoTraits[].trait
    #
    # Following links for every association can be expensive; cap link calls per rsID.
    max_trait_link_calls = 10
    trait_link_calls = 0
    seen_link_hrefs: set[str] = set()

    while next_url and pages < max_pages:
        pages += 1
        assoc_page, err2 = _request_json("GET", next_url, headers=headers, timeout_s=(6, 40), rate_limit_s=0.34)
        if err2:
            break
        embedded = assoc_page.get("_embedded") or {}
        assocs = embedded.get("associations") or embedded.get("association") or []
        if isinstance(assocs, dict):
            assocs = [assocs]
        for a in assocs:
            n_associations += 1
            # Trait can appear in different places depending on API version
            trait = (
                (a.get("trait") or "")
                or (a.get("diseaseTrait") or "")
                or (a.get("reportedTrait") or "")
            )
            trait = str(trait).strip()
            if trait:
                assoc_traits.append(trait)
                continue

            # Follow links to recover traits (capped)
            a_links = a.get("_links") or {}

            # 1) Study -> diseaseTrait.trait
            if trait_link_calls < max_trait_link_calls:
                study_href = ((a_links.get("study") or {}).get("href") or "").strip()
                if study_href and study_href not in seen_link_hrefs:
                    seen_link_hrefs.add(study_href)
                    trait_link_calls += 1
                    study_page, errS = _request_json(
                        "GET",
                        study_href,
                        headers=headers,
                        timeout_s=(6, 40),
                        rate_limit_s=0.34,
                    )
                    if not errS and isinstance(study_page, dict):
                        dt = study_page.get("diseaseTrait") or {}
                        if isinstance(dt, dict):
                            st_trait = str(dt.get("trait") or "").strip()
                            if st_trait:
                                assoc_traits.append(st_trait)

            # 2) EFO traits -> embedded.efoTraits[].trait
            if trait_link_calls < max_trait_link_calls:
                efo_href = ((a_links.get("efoTraits") or {}).get("href") or "").strip()
                if efo_href and efo_href not in seen_link_hrefs:
                    seen_link_hrefs.add(efo_href)
                    trait_link_calls += 1
                    efo_page, errE = _request_json(
                        "GET",
                        efo_href,
                        headers=headers,
                        timeout_s=(6, 40),
                        rate_limit_s=0.34,
                    )
                    if not errE and isinstance(efo_page, dict):
                        efo_emb = efo_page.get("_embedded") or {}
                        efo_traits = efo_emb.get("efoTraits") or []
                        if isinstance(efo_traits, dict):
                            efo_traits = [efo_traits]
                        for t in efo_traits:
                            name = str((t or {}).get("trait") or "").strip()
                            if name:
                                assoc_traits.append(name)

        # HAL pagination: _links.next.href
        next_url = (((assoc_page.get("_links") or {}).get("next") or {}).get("href") or "").strip()

    out = {
        "snp_reported_genes": snp.get("reportedGene") or snp.get("genes") or None,
        "associations_count_capped": n_associations,
        "traits": _stable_unique([t for t in assoc_traits if t]),
        "source": "GWAS Catalog REST (capped pages)",
    }
    return out, None


def fetch_opentargets_best_effort(rsid: str) -> Tuple[Dict[str, Any], Optional[str]]:
    """
    Best-effort: Open Targets Platform GraphQL.
    We try to search the rsID and extract any trait/phenotype hints.
    This may fail depending on schema / availability; in that case we return {} with error.
    """
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    query = """
    query SearchVariant($q: String!) {
      search(queryString: $q, entityNames: ["variant"], page: { index: 0, size: 5 }) {
        hits {
          id
          entity
          name
          description
        }
      }
    }
    """

    payload = {"query": query, "variables": {"q": rsid}}
    data, err = _request_json("POST", url, headers=headers, json_body=payload, timeout_s=(6, 30), rate_limit_s=0.34)
    if err:
        return {}, err

    hits = (((data or {}).get("data") or {}).get("search") or {}).get("hits") or []
    hit_summaries = []
    for h in hits:
        if (h.get("entity") or "").lower() != "variant":
            continue
        hit_summaries.append(
            {
                "id": h.get("id"),
                "name": h.get("name"),
                "description": h.get("description"),
            }
        )

    # Keep it minimal: the presence of a variant hit is already a useful cross-ref.
    return {"variant_hits": hit_summaries, "source": "Open Targets Platform (search)"}, None


def _normalize_text(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip().lower())


def keyword_hits(texts: Sequence[str], keywords: Sequence[str]) -> Tuple[int, List[str]]:
    corpus = " | ".join([_normalize_text(t) for t in texts if t])
    hits: List[str] = []
    for kw in keywords:
        nkw = _normalize_text(kw)
        if nkw and nkw in corpus:
            hits.append(kw)
    return len(hits), hits


def compute_score(meta: VariantMeta, traits: Sequence[str]) -> Tuple[float, List[str]]:
    """
    Score rsID for pigmentation evidence:
      - keyword matches in traits (strongest)
      - proximity to center (smaller distance is better)
      - functional impact (HIGH/MODERATE)
    """
    n_hits, matched = keyword_hits(traits, PIGMENTATION_KEYWORDS)

    # Trait-based component
    score = 0.0
    score += 12.0 * n_hits
    score += 1.5 * math.log1p(len(traits))

    # Distance component (assume central window ±50bp typical here)
    # Map distance [0..50] to [1..0]
    d = max(0, meta.min_distance_to_center)
    dist_weight = max(0.0, 1.0 - min(d, 50) / 50.0)
    score += 6.0 * dist_weight

    # Functional impact component
    imp = (meta.worst_impact or "").upper()
    if imp == "HIGH":
        score += 4.0
    elif imp == "MODERATE":
        score += 2.0

    return score, matched


def parse_variant_meta() -> Dict[str, VariantMeta]:
    """
    Parse vep_annotated_variants.tsv and build per-rsID metadata.
    Note: One rsID may appear multiple times; we keep the 'best' (smallest distance).
    """
    raise RuntimeError("parse_variant_meta() must be called with file paths via main().")
    per_rsid: Dict[str, VariantMeta] = {}
    for r in rows:
        rsids = [x.strip() for x in (r.get("rsids") or "").split(",") if x.strip().startswith("rs")]
        if not rsids:
            continue

        chrom = r.get("chrom") or ""
        pos_s = r.get("pos_1based") or "0"
        try:
            pos = int(pos_s)
        except ValueError:
            pos = 0

        dist_s = r.get("min_distance_to_center") or "0"
        try:
            dist = int(dist_s)
        except ValueError:
            dist = 0

        worst_impact = r.get("worst_impact") or ""
        msc = r.get("most_severe_consequence") or ""
        genes = [g.strip() for g in (r.get("gene_symbols") or "").split(",") if g.strip()]

        for rsid in rsids:
            meta = VariantMeta(
                rsid=rsid,
                chrom=chrom,
                pos_1based=pos,
                min_distance_to_center=dist,
                worst_impact=worst_impact,
                most_severe_consequence=msc,
                gene_symbols=genes,
            )
            if rsid not in per_rsid:
                per_rsid[rsid] = meta
            else:
                # keep the closer-to-center instance as representative
                if meta.min_distance_to_center < per_rsid[rsid].min_distance_to_center:
                    per_rsid[rsid] = meta
    return per_rsid


def parse_gene_set() -> List[str]:
    raise RuntimeError("parse_gene_set() must be called with file paths via main().")
    genes = [r.get("gene") or "" for r in rows if (r.get("gene") or "").strip()]
    return sorted(_stable_unique([g for g in genes if g]))


def main() -> int:
    ap = argparse.ArgumentParser(description="Pigmentation validation v3 (external evidence cross-ref).")
    ap.add_argument(
        "--indir",
        type=str,
        default=str(DEFAULT_INPUT_DIR),
        help="Input directory containing annotate_deeplift_windows.py outputs (TSVs).",
    )
    ap.add_argument(
        "--out",
        type=str,
        default=None,
        help="Output markdown path. Default: <indir>/pigmentation_validation3.md",
    )
    ap.add_argument(
        "--cache-dir",
        type=str,
        default=None,
        help="Cache directory for per-rsID JSON. Default: <indir>/validation3_cache",
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=1,
        help="Print progress every N rsIDs while fetching evidence (default: 1). Use 0 to disable per-rsID logs.",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of concurrent workers for fetching rsID evidence (default: 1).",
    )
    ap.add_argument(
        "--disable-gwas",
        action="store_true",
        help="Disable GWAS Catalog calls (use if EBI is slow/unreachable).",
    )
    ap.add_argument(
        "--disable-opentargets",
        action="store_true",
        help="Disable Open Targets calls.",
    )
    ap.add_argument(
        "--disable-ensembl",
        action="store_true",
        help="Disable Ensembl Variation calls.",
    )
    ap.add_argument(
        "--connect-timeout",
        type=int,
        default=6,
        help="Connection timeout (seconds) for external APIs.",
    )
    ap.add_argument(
        "--read-timeout",
        type=int,
        default=30,
        help="Read timeout (seconds) for external APIs.",
    )
    args = ap.parse_args()

    input_dir = Path(args.indir).expanduser().resolve()
    vep_ann_tsv = input_dir / "vep_annotated_variants.tsv"
    occs_tsv = input_dir / "variant_occurrences.tsv"
    gene_sum_tsv = input_dir / "summary_by_gene.tsv"

    out_md = Path(args.out).expanduser().resolve() if args.out else (input_dir / "pigmentation_validation3.md")
    global OUT_MD
    OUT_MD = out_md

    cache_dir = Path(args.cache_dir).expanduser().resolve() if args.cache_dir else (input_dir / "validation3_cache")
    global CACHE_DIR
    CACHE_DIR = cache_dir

    print(f"[INFO] Input dir: {input_dir}", flush=True)
    print(f"[INFO] Output report: {OUT_MD}", flush=True)
    print(f"[INFO] Cache dir: {CACHE_DIR} (schema v{CACHE_SCHEMA_VERSION})", flush=True)

    for p in (vep_ann_tsv, occs_tsv, gene_sum_tsv):
        if not p.exists():
            raise SystemExit(f"Missing required input: {p}")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    def parse_variant_meta_from_tsv(vep_path: Path) -> Dict[str, VariantMeta]:
        """
        Parse vep_annotated_variants.tsv and build per-rsID metadata.
        Note: One rsID may appear multiple times; we keep the 'best' (smallest distance).
        """
        rows = _read_tsv(vep_path)
        per_rsid_local: Dict[str, VariantMeta] = {}
        for r in rows:
            rsids = [x.strip() for x in (r.get("rsids") or "").split(",") if x.strip().startswith("rs")]
            if not rsids:
                continue

            chrom = r.get("chrom") or ""
            pos_s = r.get("pos_1based") or "0"
            try:
                pos = int(pos_s)
            except ValueError:
                pos = 0

            dist_s = r.get("min_distance_to_center") or "0"
            try:
                dist = int(dist_s)
            except ValueError:
                dist = 0

            worst_impact = r.get("worst_impact") or ""
            msc = r.get("most_severe_consequence") or ""
            genes_local = [g.strip() for g in (r.get("gene_symbols") or "").split(",") if g.strip()]

            for rsid in rsids:
                meta = VariantMeta(
                    rsid=rsid,
                    chrom=chrom,
                    pos_1based=pos,
                    min_distance_to_center=dist,
                    worst_impact=worst_impact,
                    most_severe_consequence=msc,
                    gene_symbols=genes_local,
                )
                if rsid not in per_rsid_local:
                    per_rsid_local[rsid] = meta
                else:
                    if meta.min_distance_to_center < per_rsid_local[rsid].min_distance_to_center:
                        per_rsid_local[rsid] = meta
        return per_rsid_local

    def parse_gene_set_from_tsv(gene_sum_path: Path) -> List[str]:
        rows = _read_tsv(gene_sum_path)
        genes_local = [r.get("gene") or "" for r in rows if (r.get("gene") or "").strip()]
        return sorted(_stable_unique([g for g in genes_local if g]))

    print(f"[INFO] Reading TSVs...", flush=True)
    per_rsid = parse_variant_meta_from_tsv(vep_ann_tsv)
    genes = parse_gene_set_from_tsv(gene_sum_tsv)
    print(f"[INFO] Parsed rsIDs: {len(per_rsid)}", flush=True)
    print(f"[INFO] Parsed genes (from summary_by_gene.tsv): {len(genes)} -> {', '.join(genes)}", flush=True)

    rsids = sorted(per_rsid.keys(), key=lambda r: (per_rsid[r].min_distance_to_center, r))
    if not rsids:
        OUT_MD.write_text("# Pigmentation Hypothesis Validation (v3)\n\nNo rsIDs found.\n", encoding="utf-8")
        return 0

    # Fetch evidence with caching
    evidence_map: Dict[str, Evidence] = {}
    cache_hits = 0
    cache_misses = 0
    total = len(rsids)
    sources = []
    if not args.disable_gwas:
        sources.append("gwas")
    if not args.disable_opentargets:
        sources.append("opentargets")
    if not args.disable_ensembl:
        sources.append("ensembl")

    timeout_pair = (max(1, int(args.connect_timeout)), max(1, int(args.read_timeout)))
    print(
        f"[INFO] Fetching external evidence for {total} rsIDs (progress_every={args.progress_every}) "
        f"sources={sources} timeout={timeout_pair} ...",
        flush=True,
    )

    print_lock = threading.Lock()

    def log(msg: str) -> None:
        with print_lock:
            print(msg, flush=True)

    # Global rate limiters per source (shared across workers).
    # Keep overall request rate conservative to avoid 429s.
    rl_lock = threading.Lock()
    next_allowed: Dict[str, float] = {"gwas": 0.0, "ensembl": 0.0, "opentargets": 0.0}
    min_interval_s: Dict[str, float] = {"gwas": 0.5, "ensembl": 0.34, "opentargets": 0.34}

    def _rate_limit(source: str) -> None:
        with rl_lock:
            now = time.time()
            t = next_allowed.get(source, 0.0)
            if now < t:
                sleep_for = t - now
            else:
                sleep_for = 0.0
            next_allowed[source] = max(now, t) + min_interval_s.get(source, 0.34)
        if sleep_for > 0:
            time.sleep(sleep_for)

    def fetch_one(rsid: str, idx: int) -> Tuple[str, Optional[Evidence], bool]:
        """
        Returns: (rsid, evidence_or_none, was_cache_hit)
        """
        cached = load_cached_evidence(rsid)
        if cached:
            return rsid, cached, True

        errors: Dict[str, str] = {}

        # Progress log BEFORE doing network, so hangs are visible.
        if args.progress_every and (idx % args.progress_every == 0):
            log(f"[INFO] ({idx}/{total}) Fetching rsID {rsid} ...")

        gwas: Dict[str, Any] = {}
        if not args.disable_gwas:
            _rate_limit("gwas")
            gwas, err_g = fetch_gwas_catalog(rsid)
            if err_g:
                errors["gwas"] = err_g
                gwas = {}

        ot: Dict[str, Any] = {}
        if not args.disable_opentargets:
            _rate_limit("opentargets")
            ot, err_ot = fetch_opentargets_best_effort(rsid)
            if err_ot:
                errors["opentargets"] = err_ot
                ot = {}

        ens: Dict[str, Any] = {}
        if not args.disable_ensembl:
            _rate_limit("ensembl")
            ens, err_e = fetch_ensembl_variation(rsid)
            if err_e:
                errors["ensembl"] = err_e
                ens = {}

        ev = Evidence(
            rsid=rsid,
            fetched_at=_now_iso(),
            cache_schema_version=CACHE_SCHEMA_VERSION,
            gwascatalog=gwas,
            opentargets=ot,
            ensembl=ens,
            errors=errors,
        )
        save_cached_evidence(ev)

        if args.progress_every and (idx % args.progress_every == 0):
            gwas_n = len((gwas.get("traits") or [])) if gwas else 0
            ens_n = len((ens.get("phenotypes") or [])) if ens else 0
            ot_n = len((ot.get("variant_hits") or [])) if ot else 0
            err_keys = ",".join(sorted(errors.keys())) if errors else ""
            if err_keys:
                # Show both the list of failing sources and a compact reason per source.
                reasons = []
                for k in sorted(errors.keys()):
                    v = str(errors.get(k, "") or "").strip().replace("\n", " ")
                    if len(v) > 120:
                        v = v[:117] + "..."
                    reasons.append(f"{k}={v}" if v else k)
                reasons_str = "; ".join(reasons)
                suffix = f" errors_count={len(errors)} errors_sources=[{err_keys}] errors_reasons=({reasons_str})"
            else:
                suffix = ""
            log(f"[INFO] rsID {rsid}: GWAS_traits={gwas_n} Ensembl_phenotypes={ens_n} OT_hits={ot_n}{suffix}")

        return rsid, ev, False

    workers = max(1, int(args.workers))
    if workers > 1:
        log(f"[INFO] Using {workers} workers with global per-source rate limiting.")

    # Submit all rsIDs; caching prevents re-fetch.
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = []
        for idx, rsid in enumerate(rsids, start=1):
            futures.append(ex.submit(fetch_one, rsid, idx))

        for fut in as_completed(futures):
            rsid, ev, was_hit = fut.result()
            if ev:
                evidence_map[rsid] = ev
            if was_hit:
                cache_hits += 1
            else:
                cache_misses += 1

    print(f"[INFO] Cache hits: {cache_hits} | cache misses (fetched): {cache_misses}", flush=True)

    # Compute scores
    print("[INFO] Computing scores and generating report...", flush=True)
    scored: List[Tuple[str, float, List[str], List[str], Dict[str, Any]]] = []
    for rsid in rsids:
        meta = per_rsid[rsid]
        ev = evidence_map.get(rsid)
        traits: List[str] = []
        if ev:
            traits.extend(ev.gwascatalog.get("traits") or [])
            traits.extend(ev.ensembl.get("phenotypes") or [])
        traits = _stable_unique([t for t in traits if t])
        score, matched = compute_score(meta, traits)
        scored.append((rsid, score, matched, traits, {"meta": meta, "ev": ev}))

    scored.sort(key=lambda x: (-x[1], per_rsid[x[0]].min_distance_to_center, x[0]))

    # Summary counts
    def _has_any_external(ev: Evidence) -> bool:
        if not ev:
            return False
        if (ev.gwascatalog.get("traits") or []):
            return True
        if (ev.ensembl.get("phenotypes") or []):
            return True
        if (ev.opentargets.get("variant_hits") or []):
            return True
        return False

    n_total = len(rsids)
    n_with_any = sum(1 for r in rsids if _has_any_external(evidence_map.get(r)))
    n_with_kw = sum(1 for (r, _s, matched, traits, _ctx) in scored if matched and traits)

    top_n = 30
    top_rows = scored[:top_n]
    kw_hit_rows = [(r, s, matched, traits, ctx) for (r, s, matched, traits, ctx) in scored if matched]

    # Aggregate by gene (use the meta.gene_symbols list; keep any that are in the selected genes set when possible)
    gene_scores: Dict[str, Dict[str, Any]] = {g: {"score": 0.0, "rsids": []} for g in genes}
    for (rsid, s, matched, traits, ctx) in scored:
        meta: VariantMeta = ctx["meta"]
        for g in (meta.gene_symbols or []):
            if g not in gene_scores:
                gene_scores[g] = {"score": 0.0, "rsids": []}
            gene_scores[g]["score"] += s
            gene_scores[g]["rsids"].append(rsid)

    gene_rank = sorted(gene_scores.items(), key=lambda x: (-x[1]["score"], x[0]))

    # Hypothesis appraisal (simple, explicit heuristics)
    known_pigmentation_genes = {"SLC24A5", "OCA2", "HERC2", "TYR", "TYRP1", "MC1R", "SLC45A2"}
    top_genes_in_input = set(genes)
    known_in_input = sorted(known_pigmentation_genes.intersection(top_genes_in_input))
    kw_hit_rsids = [r for (r, _s, matched, _traits, _ctx) in scored if matched]
    kw_hit_genes = sorted({g for r in kw_hit_rsids for g in (per_rsid.get(r).gene_symbols if per_rsid.get(r) else [])})

    # Appraisal tiers:
    # - Strong support: >=3 keyword-hit rsIDs AND at least one hit maps to known pigmentation genes
    # - Moderate support: >=1 keyword-hit rsID OR >=2 rsIDs with "hair/eye/skin/pigmentation/melanin/albinism" in traits
    # - Weak/inconclusive: no keyword hits and low external evidence coverage
    if n_with_kw >= 3 and (set(kw_hit_genes) & known_pigmentation_genes):
        appraisal = "STRONG SUPPORT"
    elif n_with_kw >= 1:
        appraisal = "MODERATE SUPPORT"
    elif n_with_any >= max(5, int(0.05 * n_total)):
        appraisal = "WEAK / INCONCLUSIVE (external evidence exists, but no pigmentation-keyword matches)"
    else:
        appraisal = "INCONCLUSIVE (insufficient external evidence coverage)"

    # Write report (English)
    lines: List[str] = []
    lines.append("# Pigmentation Hypothesis Validation (v3)\n")
    lines.append("This report validates whether DeepLIFT-selected central variants show **external published evidence** related to pigmentation.\n")
    lines.append("**Data sources (best-effort):** GWAS Catalog (EBI), Open Targets (Platform), and Ensembl Variation.\n")
    lines.append("\n---\n")

    lines.append("## Inputs\n")
    lines.append(f"- `vep_annotated_variants.tsv`: `{vep_ann_tsv}`\n")
    lines.append(f"- `variant_occurrences.tsv`: `{occs_tsv}`\n")
    lines.append(f"- `summary_by_gene.tsv`: `{gene_sum_tsv}`\n")
    lines.append(f"- Cache directory: `{CACHE_DIR}`\n")
    lines.append("\n---\n")

    lines.append("## Summary\n")
    lines.append(f"- Total rsIDs in central window: **{n_total}**\n")
    lines.append(f"- rsIDs with any external annotations (traits/phenotypes/hits): **{n_with_any}**\n")
    lines.append(f"- rsIDs with pigmentation keyword hits: **{n_with_kw}**\n")
    lines.append("\n---\n")

    lines.append("## Keyword model\n")
    lines.append("Pigmentation-related keywords used (case-insensitive substring match):\n\n")
    for kw in PIGMENTATION_KEYWORDS:
        lines.append(f"- `{kw}`\n")
    lines.append("\n---\n")

    lines.append(f"## Top rsIDs by pigmentation-evidence score (top {top_n})\n\n")
    lines.append("| Rank | rsID | Score | Min dist to center (bp) | Impact | Consequence | Gene(s) | Matched keywords | Evidence traits (snippet) |\n")
    lines.append("|---:|---|---:|---:|---|---|---|---|---|\n")
    for idx, (rsid, score, matched, traits, ctx) in enumerate(top_rows, start=1):
        meta: VariantMeta = ctx["meta"]
        genes_str = ",".join(meta.gene_symbols) if meta.gene_symbols else ""
        matched_str = ",".join(matched) if matched else ""
        trait_snip = "; ".join(traits[:3]) + (f" (+{len(traits)-3} more)" if len(traits) > 3 else "")
        lines.append(
            f"| {idx} | {rsid} | {score:.3f} | {meta.min_distance_to_center} | {meta.worst_impact} | "
            f"{meta.most_severe_consequence} | {genes_str} | {matched_str} | {trait_snip} |\n"
        )
    lines.append("\n---\n")

    lines.append("## Pigmentation-specific hits (keyword matched)\n\n")
    if kw_hit_rows:
        lines.append("| rsID | Score | Min dist (bp) | Gene(s) | Matched keywords | Evidence traits (top) |\n")
        lines.append("|---|---:|---:|---|---|---|\n")
        for (rsid, score, matched, traits, ctx) in kw_hit_rows[:50]:
            meta: VariantMeta = ctx["meta"]
            genes_str = ",".join(meta.gene_symbols) if meta.gene_symbols else ""
            matched_str = ",".join(matched) if matched else ""
            trait_snip = "; ".join(traits[:5]) + (f" (+{len(traits)-5} more)" if len(traits) > 5 else "")
            lines.append(
                f"| {rsid} | {score:.3f} | {meta.min_distance_to_center} | {genes_str} | {matched_str} | {trait_snip} |\n"
            )
        lines.append("\n")
        lines.append(
            "Interpretation: the presence of pigmentation-related traits/phenotypes among DeepLIFT-centered variants "
            "supports that the pipeline is capturing biology connected to pigmentation.\n\n"
        )
    else:
        lines.append("No rsIDs matched the pigmentation keyword set in the retrieved external traits/phenotypes.\n\n")

    lines.append("## Gene-level aggregation (sum of rsID scores)\n\n")
    lines.append("| Gene | Aggregated score | #rsIDs |\n")
    lines.append("|---|---:|---:|\n")
    for gene, info in gene_rank:
        rs_gene = _stable_unique(info["rsids"])
        lines.append(f"| {gene} | {info['score']:.3f} | {len(rs_gene)} |\n")
    lines.append("\n---\n")

    lines.append("## Evidence details (top rsIDs)\n\n")
    for rsid, score, matched, traits, ctx in top_rows:
        meta: VariantMeta = ctx["meta"]
        ev: Evidence = ctx["ev"]
        lines.append(f"### {rsid}\n\n")
        lines.append(f"- Location: **{meta.chrom}:{meta.pos_1based}**, min distance to center: **{meta.min_distance_to_center}bp**\n")
        lines.append(f"- Impact/consequence: **{meta.worst_impact} / {meta.most_severe_consequence}**\n")
        if meta.gene_symbols:
            lines.append(f"- Gene symbols (VEP): **{', '.join(meta.gene_symbols)}**\n")
        if matched:
            lines.append(f"- Matched keywords: **{', '.join(matched)}**\n")
        lines.append("\n")

        # GWAS
        gwas_traits = (ev.gwascatalog.get("traits") or []) if ev else []
        if gwas_traits:
            lines.append("**GWAS Catalog traits (capped):**\n\n")
            for t in gwas_traits[:10]:
                lines.append(f"- {t}\n")
            if len(gwas_traits) > 10:
                lines.append(f"- ... (+{len(gwas_traits)-10} more)\n")
            lines.append("\n")
        else:
            lines.append("**GWAS Catalog traits:** (none found or fetch failed)\n\n")

        # Ensembl
        ens_traits = (ev.ensembl.get("phenotypes") or []) if ev else []
        if ens_traits:
            lines.append("**Ensembl Variation phenotypes:**\n\n")
            for t in ens_traits[:10]:
                lines.append(f"- {t}\n")
            if len(ens_traits) > 10:
                lines.append(f"- ... (+{len(ens_traits)-10} more)\n")
            lines.append("\n")
        else:
            lines.append("**Ensembl Variation phenotypes:** (none found or fetch failed)\n\n")

        # Open Targets
        ot_hits = (ev.opentargets.get("variant_hits") or []) if ev else []
        if ot_hits:
            lines.append("**Open Targets (variant search hits):**\n\n")
            for h in ot_hits[:5]:
                hid = h.get("id") or ""
                hname = h.get("name") or ""
                hdesc = h.get("description") or ""
                lines.append(f"- id={hid} name={hname} desc={hdesc}\n")
            lines.append("\n")
        else:
            lines.append("**Open Targets:** (no variant hits or fetch failed)\n\n")

        # Errors
        if ev and ev.errors:
            lines.append("**Fetch notes:**\n\n")
            for src, err in sorted(ev.errors.items()):
                lines.append(f"- {src}: {err}\n")
            lines.append("\n")

        lines.append("---\n\n")

    lines.append("## Robustness / Failure modes\n\n")
    lines.append("- External APIs can return incomplete data (rate limits, missing annotations, schema drift).\n")
    lines.append("- This report is **best-effort** and caches per-rsID responses for reproducibility.\n")
    lines.append("- If GWAS/Open Targets are empty but Ensembl has phenotypes, the score may still capture pigmentation terms.\n")
    lines.append("- The scoring is heuristic: it prioritizes central variants and known functional impacts.\n")
    lines.append("\n---\n\n")

    lines.append("## Conclusion / Hypothesis appraisal\n\n")
    lines.append(
        "Hypothesis: *\"It is possible to identify what in the DNA causes Africans to have higher skin pigmentation "
        "using AlphaGenome outputs processed by `neural_ancestry_predictor.py` and `annotate_deeplift_windows.py`.\"*\n\n"
    )
    lines.append(f"### Appraisal: **{appraisal}**\n\n")
    lines.append("This conclusion is based on the following measurable signals extracted from the pipeline outputs:\n\n")
    lines.append(f"- Total rsIDs analyzed (central-window variants): **{n_total}**\n")
    lines.append(f"- rsIDs with any external evidence (GWAS/Ensembl/Open Targets): **{n_with_any}**\n")
    lines.append(f"- rsIDs with pigmentation-keyword hits in external traits/phenotypes: **{n_with_kw}**\n")
    if known_in_input:
        lines.append(f"- Known pigmentation genes present in the analyzed gene set: **{', '.join(known_in_input)}**\n")
    else:
        lines.append("- Known pigmentation genes present in the analyzed gene set: **none detected**\n")
    if kw_hit_genes:
        lines.append(f"- Genes implicated by keyword-hit rsIDs (VEP mapping): **{', '.join(kw_hit_genes)}**\n")
    else:
        lines.append("- Genes implicated by keyword-hit rsIDs (VEP mapping): **none**\n")
    lines.append("\n")

    lines.append("### Interpretation\n\n")
    lines.append(
        "- **What this supports**: if DeepLIFT-centered variants overlap rsIDs whose external annotations mention pigmentation-related "
        "phenotypes (e.g., skin pigmentation, hair/eye color, melanin, albinism), then the DeepLIFT→variant-calling→annotation pipeline "
        "is capturing biology consistent with pigmentation differences.\n"
    )
    lines.append(
        "- **What this does NOT prove**: this report does not establish causal direction (\"higher\" pigmentation) nor quantify "
        "population-level effect sizes; it is evidence of biological relevance/consistency, not a full causal model.\n"
    )
    lines.append(
        "- **Next validation step (strongest)**: run the same pipeline on a contrasting target class (e.g., EUR vs AFR) and check whether "
        "the keyword-hit signal and implicated genes shift in the expected direction (e.g., `SLC24A5`, `HERC2/OCA2`, `TYR`).\n"
    )

    OUT_MD.write_text("".join(lines), encoding="utf-8")
    print(f"[SAVED] {OUT_MD}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


