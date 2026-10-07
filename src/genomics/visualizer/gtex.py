"""GTEx portal API (v2) lookups for the Variant page: eQTLs of a variant, measured per tissue.

* ``variant``: resolves ``chr15_28120472_A_G_b38`` <-> rsID (GTEx v8 variant table, GRCh38).
* ``significant``: the variant's significant single-tissue eQTLs (every gene, every tissue).
* ``dynamic``: GTEx's association test for one gene x tissue, significant or not, with the donors'
  normalized expression by genotype (what the GTEx portal's violin plot shows).

GTEx's NES is the effect of the ALT allele on normalized expression; its sign is comparable with the
AlphaGenome dosage slope. Responses go through :class:`RemoteCache` (cached on disk, ``--no-remote``).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from genomics.visualizer.remote import DAY, RemoteCache, RemoteError

API = "https://gtexportal.org/api/v2"
DATASET = "gtex_v8"
GENCODE = "v26"
DEFAULT_TISSUES = ("Skin_Sun_Exposed_Lower_leg", "Skin_Not_Sun_Exposed_Suprapubic")
TTL = 90 * DAY  # GTEx v8 is a frozen release


def tissue_name_key(name: str) -> str:
    """'Thyroid gland' / 'Thyroid' -> 'thyroid'; 'Adrenal Gland' -> 'adrenal'."""
    return " ".join(w for w in str(name).lower().replace("_", " ").split() if w != "gland")


class GtexError(RuntimeError):
    pass


class GtexClient:
    def __init__(self, remote: RemoteCache):
        self.remote = remote

    def _get(self, path: str, **params: Any) -> Dict[str, Any]:
        try:
            data = self.remote.get_json(f"{API}/{path}", ttl=TTL, params={"datasetId": DATASET, **params})
        except RemoteError as exc:
            # GTEx answers unknown variants/genes with HTTP 400 and a JSON "detail".
            raise GtexError(str(exc)) from exc
        if not isinstance(data, dict):
            raise GtexError(f"{path}: unexpected response")
        return data

    def tissues(self) -> List[Dict[str, Any]]:
        data = self._get("dataset/tissueSiteDetail", itemsPerPage=250).get("data") or []
        return [
            {
                "id": t["tissueSiteDetailId"], "name": t.get("tissueSiteDetail") or t["tissueSiteDetailId"],
                "ontology": t.get("ontologyId"), "color": f"#{t['colorHex']}" if t.get("colorHex") else None,
                "samples": ((t.get("eqtlSampleSummary") or {}).get("totalCount")),
            }
            for t in data
        ]

    def variant(self, variant_id: Optional[str] = None, snp_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
        params = {"variantId": variant_id} if variant_id else {"snpId": snp_id}
        try:
            rows = self._get("dataset/variant", **params).get("data") or []
        except GtexError as exc:
            if _not_found(exc):
                return None
            raise
        return rows[0] if rows else None

    def gene(self, symbol: str) -> Optional[Dict[str, Any]]:
        try:
            rows = self._get("reference/gene", geneId=symbol, gencodeVersion=GENCODE, genomeBuild="GRCh38/hg38").get("data") or []
        except GtexError as exc:
            if _not_found(exc):
                return None
            raise
        exact = [r for r in rows if str(r.get("geneSymbol", "")).upper() == symbol.upper()]
        return (exact or rows or [None])[0]

    def median_expression(self, gencode_id: str, tissue: str) -> Optional[float]:
        """GTEx v8 median TPM of a gene in a tissue."""
        rows = self._get("expression/medianGeneExpression", gencodeId=gencode_id, tissueSiteDetailId=tissue).get("data") or []
        return float(rows[0]["median"]) if rows else None

    def median_transcripts(self, gencode_id: str, tissue: Optional[str]) -> Dict[str, float]:
        """GTEx v8 median TPM of each transcript of a gene in a tissue (``None``: mean over every tissue),
        by versionless Ensembl id (GENCODE v26)."""
        params = {"gencodeId": gencode_id, "itemsPerPage": 10000}
        if tissue:
            params["tissueSiteDetailId"] = tissue
        rows = self._get("expression/medianTranscriptExpression", **params).get("data") or []
        sums: Dict[str, List[float]] = {}
        for r in rows:
            if r.get("transcriptId"):
                sums.setdefault(str(r["transcriptId"]).split(".")[0], []).append(float(r.get("median") or 0.0))
        return {t: sum(v) / len(v) for t, v in sums.items()}

    def name_tissues(self) -> Dict[str, str]:
        """{lower-case name without "gland": GTEx tissue id} for single-site GTEx tissues ("Liver", "Adrenal Gland")."""
        out: Dict[str, List[str]] = {}
        for t in self.tissues():
            name = str(t.get("name") or "")
            if " - " in name or "(" in name:
                continue
            out.setdefault(tissue_name_key(name), []).append(t["id"])
        return {k: v[0] for k, v in out.items() if len(v) == 1}

    def ontology_tissues(self) -> Dict[str, str]:
        """{ontology term: GTEx tissue id} for GTEx tissues whose term is unique to them."""
        counts: Dict[str, List[str]] = {}
        for t in self.tissues():
            if t.get("ontology"):
                counts.setdefault(t["ontology"], []).append(t["id"])
        return {k: v[0] for k, v in counts.items() if len(v) == 1}

    def transcript_introns(self, gencode_id: str) -> Dict[str, List[tuple]]:
        """{versionless transcript id: intron chain [(last base of exon, first base of next exon)]} in GTEx v8's GENCODE v26."""
        rows = self._get("reference/exon", gencodeId=gencode_id).get("data") or []
        exons: Dict[str, List[tuple]] = {}
        for r in rows:
            exons.setdefault(str(r["transcriptId"]).split(".")[0], []).append((int(r["start"]), int(r["end"])))
        return {t: [(a[1], b[0]) for a, b in zip(sorted(e), sorted(e)[1:])] for t, e in exons.items()}

    def significant(self, variant_id: str) -> List[Dict[str, Any]]:
        try:
            rows = self._get("association/singleTissueEqtl", variantId=variant_id, itemsPerPage=250).get("data") or []
        except GtexError as exc:
            if _not_found(exc):
                return []
            raise
        out = [
            {"gene": r.get("geneSymbol"), "gencode_id": r.get("gencodeId"), "tissue": r.get("tissueSiteDetailId"), "nes": r.get("nes"), "p": r.get("pValue")}
            for r in rows
        ]
        return sorted(out, key=lambda r: (r["p"] if r["p"] is not None else 1.0))

    def dynamic(self, gencode_id: str, variant_id: str, tissue: str, max_points: int = 1000) -> Dict[str, Any]:
        try:
            d = self._get("association/dyneqtl", gencodeId=gencode_id, variantId=variant_id, tissueSiteDetailId=tissue)
        except GtexError as exc:
            return {"tissue": tissue, "error": "not tested in GTEx for this gene and tissue" if _not_found(exc) else str(exc)}
        values = d.get("data") or []
        genotypes = d.get("genotypes") or []
        by_genotype: Dict[int, List[float]] = {0: [], 1: [], 2: []}
        for v, g in zip(values, genotypes):
            if g in by_genotype and v is not None:
                by_genotype[g].append(float(v))
        for g, vals in by_genotype.items():
            if len(vals) > max_points:
                step = len(vals) / max_points
                by_genotype[g] = [vals[int(i * step)] for i in range(max_points)]
        return {
            "tissue": tissue, "nes": d.get("nes"), "p": d.get("pValue"), "t": d.get("tStatistic"), "maf": d.get("maf"),
            "p_threshold": d.get("pValueThreshold"),
            "counts": {"0/0": d.get("homoRefCount"), "0/1": d.get("hetCount"), "1/1": d.get("homoAltCount")},
            "values": {"0/0": by_genotype[0], "0/1": by_genotype[1], "1/1": by_genotype[2]},
        }

    def report(self, variant_id: str, gene: Optional[str], tissues: List[str]) -> Dict[str, Any]:
        """Everything the Variant page shows; raises :class:`GtexError` when GTEx cannot be reached."""
        info = self.variant(variant_id=variant_id)
        out: Dict[str, Any] = {"variant_id": variant_id, "in_gtex": info is not None, "rsid": (info or {}).get("snpId"), "dataset": DATASET}
        if info is None:
            return out
        out["significant"] = self.significant(variant_id)
        if gene:
            g = self.gene(gene)
            out["gene"] = {"symbol": gene, "gencode_id": (g or {}).get("gencodeId")}
            if g and g.get("gencodeId"):
                out["dynamic"] = [self.dynamic(g["gencodeId"], variant_id, t) for t in tissues]
        return out


def _not_found(exc: Exception) -> bool:
    """GTEx answers unknown variants, genes and untested gene x tissue pairs with HTTP 400."""
    return "HTTP Error 400" in str(exc)
