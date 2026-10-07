"""What is known about the variants an individual carries: rsID, clinical significance, associated
traits and literature, with the identifiers needed to link each variant to its reference records.

* **Ensembl VEP** (``POST /vep/human/region``, batched): the co-located known variant of each
  GRCh38 allele, i.e. its dbSNP rsID, ClinVar significance, the PubMed ids that cite it, cross-references
  (ClinVar VCV, OMIM allelic variant, UniProt VAR, PharmGKB) and gnomAD / 1000 Genomes frequencies.
* **Ensembl Variation** (``GET /variation/human/<rsID>?phenotypes=1``), only for variants flagged as
  associated with a phenotype: the traits of the NHGRI-EBI GWAS Catalog and the ClinVar conditions.

Responses are cached on disk (``<cache>/variants``) and nothing is requested with ``--no-remote``.
Only positions and alleles are sent, the same as looking a variant up on dbSNP.
"""
from __future__ import annotations

import hashlib
import json
import re
import threading
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

from genomics.visualizer.remote import DAY, USER_AGENT, RemoteCache, RemoteError

VEP_REGION = "https://rest.ensembl.org/vep/human/region"
VARIATION = "https://rest.ensembl.org/variation/human/{rsid}"
BATCH = 150  # VEP accepts up to 200 variants per request
MAX_TRAITS = 12
Variant = Tuple[str, int, str, str]  # chromosome (no "chr"), POS, REF, ALT


def variant_key(pos: int, ref: str, alt: str) -> str:
    return f"{int(pos)}:{ref}:{alt}"


def _pheno_name(trait: str) -> str:
    """'Basal cell carcinoma PheCode 172.21' -> 'Basal cell carcinoma'."""
    return re.sub(r"\s+PheCode\s+[\d.]+$", "", str(trait or "")).strip()


class VariantKnowledge:
    def __init__(self, remote: RemoteCache, cache_dir: Optional[Path]):
        self.remote = remote
        self.cache_dir = Path(cache_dir) / "variants" if cache_dir else None
        self._lock = threading.Lock()

    # -- Ensembl VEP ---------------------------------------------------------------------------
    def _post(self, url: str, payload: Dict[str, Any]) -> Any:
        body = json.dumps(payload, sort_keys=True).encode("utf-8")
        key = hashlib.sha1(url.encode() + body).hexdigest()
        path = self.cache_dir / f"{key}.json" if self.cache_dir else None
        if path is not None and path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
        if self.remote.offline:
            raise RemoteError("Ensembl: not cached and remote lookups are disabled (--no-remote)")
        last: Optional[Exception] = None
        for attempt in range(4):
            request = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT, "Content-Type": "application/json", "Accept": "application/json"})
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    data = json.loads(response.read().decode("utf-8"))
                break
            except urllib.error.HTTPError as exc:  # Ensembl answers 429/503 under load
                last = exc
                if exc.code not in (429, 500, 502, 503, 504):
                    raise RemoteError(f"Ensembl VEP: HTTP {exc.code}")
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last = exc
            import time

            time.sleep(1.5 * (attempt + 1))
        else:
            raise RemoteError(f"Ensembl VEP unreachable: {last}")
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(data), encoding="utf-8")
        return data

    def colocated(self, variants: Sequence[Variant]) -> Dict[str, Dict[str, Any]]:
        """{"pos:ref:alt": known-variant record} for the variants Ensembl knows (others are absent)."""
        out: Dict[str, Dict[str, Any]] = {}
        clean = [v for v in variants if re.fullmatch(r"[ACGTN]+", v[2]) and re.fullmatch(r"[ACGTN]+", v[3])]
        for i in range(0, len(clean), BATCH):
            batch = clean[i:i + BATCH]
            lines = [f"{c} {p} . {r} {a} . . ." for c, p, r, a in batch]
            data = self._post(VEP_REGION, {"variants": lines, "pubmed": 1, "var_synonyms": 1})
            by_input = {str(item.get("input")): item for item in data or []}
            for (chrom, pos, ref, alt), line in zip(batch, lines):
                record = _known_record(by_input.get(line), ref, alt)
                if record is not None:
                    out[variant_key(pos, ref, alt)] = record
        return out

    # -- Ensembl Variation: phenotypes ------------------------------------------------------------
    def phenotypes(self, rsid: str) -> Dict[str, Any]:
        """Traits associated with ``rsid``: GWAS Catalog associations and ClinVar conditions."""
        try:
            data = self.remote.get_json(VARIATION.format(rsid=rsid), ttl=90 * DAY, params={"phenotypes": 1})
        except RemoteError as exc:
            return {"error": str(exc), "traits": []}
        traits: Dict[str, Dict[str, Any]] = {}
        for p in (data or {}).get("phenotypes") or []:
            name = _pheno_name(p.get("trait"))
            if not name or name.lower().startswith("clinvar: phenotype not specified") or name.lower() == "not provided":
                continue
            t = traits.setdefault(name, {"trait": name, "sources": set(), "studies": set(), "best_p": None, "risk_alleles": set()})
            t["sources"].add(p.get("source") or "")
            study = str(p.get("study") or "")
            if study:
                t["studies"].add(study)
            if p.get("risk_allele"):
                t["risk_alleles"].add(str(p["risk_allele"]))
            try:
                pv = float(p.get("pvalue")) if p.get("pvalue") not in (None, "") else None
            except (TypeError, ValueError):
                pv = None
            if pv is not None and (t["best_p"] is None or pv < t["best_p"]):
                t["best_p"] = pv
        rows = []
        for t in traits.values():
            rows.append({"trait": t["trait"], "sources": sorted(s for s in t["sources"] if s), "studies": sorted(t["studies"])[:8],
                         "best_p": t["best_p"], "risk_alleles": sorted(t["risk_alleles"])})
        # ClinVar conditions first, then GWAS traits by strength of evidence
        rows.sort(key=lambda r: ("ClinVar" not in r["sources"], r["best_p"] if r["best_p"] is not None else 1.0, r["trait"].lower()))
        return {"traits": rows, "clinical_significance": (data or {}).get("clinical_significance") or [], "total": len(rows)}

    def annotate(self, chromosome: str, variants: Sequence[Tuple[int, str, str]], with_phenotypes: bool = True) -> Dict[str, Any]:
        """Known-variant records (and their traits) for ``(pos, ref, alt)`` on ``chromosome``."""
        chrom = str(chromosome).replace("chr", "")
        try:
            known = self.colocated([(chrom, int(p), r, a) for p, r, a in variants])
        except RemoteError as exc:
            return {"available": False, "error": str(exc), "variants": {}}
        if with_phenotypes:
            for record in known.values():
                if record.get("phenotype_or_disease") and record.get("rsid"):
                    ph = self.phenotypes(record["rsid"])
                    record["traits"] = ph.get("traits", [])[:MAX_TRAITS]
                    record["trait_count"] = ph.get("total", 0)
        return {"available": True, "variants": known, "chromosome": chrom}


def _known_record(item: Optional[Dict[str, Any]], ref: str, alt: str) -> Optional[Dict[str, Any]]:
    """The dbSNP variant co-located with our allele in a VEP result (if any), as a compact record."""
    if not item:
        return None
    best = None
    for c in item.get("colocated_variants") or []:
        cid = str(c.get("id") or "")
        if not cid.startswith("rs"):
            continue
        alleles = str(c.get("allele_string") or "").split("/")
        # the same allele (VEP trims shared bases of indels: compare loosely for those)
        same = alt in alleles or (len(ref) != len(alt) and len(alleles) > 1)
        if best is None or same:
            best = (c, same)
            if same:
                break
    if best is None:
        return None
    c, same_allele = best
    synonyms = c.get("var_synonyms") or {}
    clinvar = [s for s in synonyms.get("ClinVar") or [] if str(s).startswith("VCV")]
    omim = [str(s) for s in synonyms.get("OMIM") or []]
    freq = (c.get("frequencies") or {}).get(alt) or {}
    clin = [str(s).replace("_", " ") for s in c.get("clin_sig") or []]
    pubmed = [int(p) for p in c.get("pubmed") or [] if str(p).isdigit()]
    return {
        "rsid": c["id"], "same_allele": bool(same_allele), "clinical_significance": clin,
        "phenotype_or_disease": bool(c.get("phenotype_or_disease")), "pubmed": pubmed[:400], "pubmed_count": len(pubmed),
        "clinvar": clinvar[:3], "omim": omim[:3], "uniprot": [str(s) for s in synonyms.get("UniProt") or []][:3],
        "pharmgkb": [str(s) for s in synonyms.get("PharmGKB") or []][:2],
        "frequency": {"gnomad_genomes": freq.get("gnomadg"), "gnomad_exomes": freq.get("gnomade"), "1000g": freq.get("af"),
                      "1000g_eur": freq.get("eur"), "1000g_afr": freq.get("afr"), "1000g_eas": freq.get("eas")},
        "known": bool(clin or pubmed or c.get("phenotype_or_disease")),
    }
