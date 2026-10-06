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
