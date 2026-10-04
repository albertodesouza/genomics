"""Gene and ontology knowledge from public databases, for gene pickers and info cards.

* **HGNC** complete set (approved symbols, names, aliases, previous symbols, gene groups and
  cross-references: Ensembl, NCBI Gene, UniProt, OMIM, RefSeq, MANE). Downloaded once
  (~17 MB) into the cache directory and refreshed monthly.
* **Gene Ontology** through QuickGO: term search, the human genes annotated to a term (and its
  descendants), and the GO annotations of a gene.
* **Ontology terms** (CL, UBERON, EFO, OBI, GO, ...) through the EBI Ontology Lookup Service
  (OLS4): label, definition, synonyms.

Every lookup goes through :class:`genomics.visualizer.remote.RemoteCache`, so it works offline
once fetched. Gene coordinates come from the dataset's GENCODE table (:mod:`gene_index`).
"""
from __future__ import annotations

import csv
import io
import re
import threading
from typing import Any, Dict, Iterable, List, Optional, Tuple

from genomics.visualizer.remote import DAY, RemoteCache, RemoteError

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
QUICKGO = "https://www.ebi.ac.uk/QuickGO/services"
OLS = "https://www.ebi.ac.uk/ols4/api"
HUMAN_TAXON = 9606
HGNC_KEYS = ("hgnc_id", "symbol", "name", "locus_group", "locus_type", "location", "alias_symbol", "prev_symbol",
             "gene_group", "gene_group_id", "entrez_id", "ensembl_gene_id", "uniprot_ids", "omim_id", "refseq_accession",
             "mane_select", "orphanet", "ucsc_id")
LIST_KEYS = ("alias_symbol", "prev_symbol", "gene_group", "gene_group_id", "uniprot_ids", "omim_id", "refseq_accession", "mane_select")
GO_ASPECTS = {"biological_process": "Biological process", "molecular_function": "Molecular function", "cellular_component": "Cellular component"}
CURIE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*:[A-Za-z0-9_.-]+$")
GO_RE = re.compile(r"^GO:\d{7}$")


def _split(value: str) -> List[str]:
    value = (value or "").strip().strip('"')
    return [v.strip() for v in value.split("|") if v.strip()] if value else []


class HgncTable:
    """Approved human gene symbols with aliases, previous symbols, groups and identifiers."""

    def __init__(self, text: str):
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        self.records: Dict[str, Dict[str, Any]] = {}
        self.by_alias: Dict[str, List[str]] = {}
        self.by_prev: Dict[str, List[str]] = {}
        self.by_ensembl: Dict[str, str] = {}
        self.by_hgnc: Dict[str, str] = {}
        self.groups: Dict[str, Dict[str, Any]] = {}
        for row in reader:
            if (row.get("status") or "Approved") != "Approved":
                continue
            rec: Dict[str, Any] = {}
            for key in HGNC_KEYS:
                raw = row.get(key) or ""
                rec[key] = _split(raw) if key in LIST_KEYS else raw.strip().strip('"')
            symbol = rec["symbol"]
            if not symbol:
                continue
            self.records[symbol] = rec
            self.by_hgnc[rec["hgnc_id"]] = symbol
            if rec["ensembl_gene_id"]:
                self.by_ensembl[rec["ensembl_gene_id"]] = symbol
            for alias in rec["alias_symbol"]:
                self.by_alias.setdefault(alias.upper(), []).append(symbol)
            for prev in rec["prev_symbol"]:
                self.by_prev.setdefault(prev.upper(), []).append(symbol)
            for name, gid in zip(rec["gene_group"], rec["gene_group_id"] + [""] * len(rec["gene_group"])):
                group = self.groups.setdefault(name, {"name": name, "id": gid, "genes": []})
                group["genes"].append(symbol)
        self.upper = {s.upper(): s for s in self.records}

    def __len__(self) -> int:
        return len(self.records)

    def get(self, symbol: str) -> Optional[Dict[str, Any]]:
        return self.records.get(self.resolve(symbol) or "")

    def resolve(self, text: str) -> Optional[str]:
        """Approved symbol for a symbol, alias, previous symbol, HGNC id or Ensembl gene id."""
        q = (text or "").strip()
        if not q:
            return None
        if q in self.records:
            return q
        up = q.upper()
        if up in self.upper:
            return self.upper[up]
        if up.startswith("HGNC:") and q in self.by_hgnc:
            return self.by_hgnc[q]
        if up.startswith("ENSG"):
            return self.by_ensembl.get(up.split(".")[0])
        for table in (self.by_prev, self.by_alias):
            hits = table.get(up)
            if hits and len(hits) == 1:
                return hits[0]
        return None

    def search(self, query: str, limit: int = 25) -> List[Tuple[Dict[str, Any], str]]:
        """(record, how it matched), best first: symbol, previous symbol, alias, then name."""
        q = query.strip()
        if not q:
            return []
        up = q.upper()
        ranked: Dict[str, Tuple[int, str]] = {}

        def add(symbol: str, rank: int, why: str) -> None:
            if symbol not in ranked or rank < ranked[symbol][0]:
                ranked[symbol] = (rank, why)

        if up in self.upper:
            add(self.upper[up], 0, "")
        if up.startswith("HGNC:") and q in self.by_hgnc:
            add(self.by_hgnc[q], 0, q)
        if up.startswith("ENSG") and up.split(".")[0] in self.by_ensembl:
            add(self.by_ensembl[up.split(".")[0]], 0, up.split(".")[0])
        for symbol in self.by_prev.get(up, []):
            add(symbol, 1, f"previous symbol {q}")
        for symbol in self.by_alias.get(up, []):
            add(symbol, 2, f"alias {q}")
        if len(ranked) < limit:
            for sym_up, symbol in self.upper.items():
                if sym_up.startswith(up):
                    add(symbol, 3, "")
            if len(up) >= 3:
                needle = q.lower()
                for symbol, rec in self.records.items():
                    if needle in rec["name"].lower():
                        add(symbol, 5, rec["name"])
                    elif len(up) >= 4 and any(a.upper().startswith(up) for a in rec["alias_symbol"]):
                        add(symbol, 4, "alias " + next(a for a in rec["alias_symbol"] if a.upper().startswith(up)))
        locus_rank = {"protein-coding gene": 0, "non-coding RNA": 1}
        order = sorted(ranked.items(), key=lambda kv: (kv[1][0], locus_rank.get(self.records[kv[0]]["locus_group"], 2), len(kv[0]), kv[0]))
        return [(self.records[s], why) for s, (_, why) in order[:limit]]

    def search_groups(self, query: str, limit: int = 25) -> List[Dict[str, Any]]:
        q = query.strip().lower()
        if not q:
            return []
        hits = [g for g in self.groups.values() if q in g["name"].lower()]
        hits.sort(key=lambda g: (not g["name"].lower().startswith(q), len(g["genes"]) > 400, g["name"].lower()))
        return [{"name": g["name"], "id": g["id"], "size": len(g["genes"])} for g in hits[:limit]]


def public_record(rec: Dict[str, Any]) -> Dict[str, Any]:
    """The fields the frontend shows and links (identifiers for HGNC, Ensembl, NCBI, UniProt, OMIM...)."""
    return {
        "symbol": rec["symbol"],
        "name": rec["name"],
        "hgnc_id": rec["hgnc_id"],
        "locus_group": rec["locus_group"],
        "locus_type": rec["locus_type"],
        "location": rec["location"],
        "aliases": rec["alias_symbol"],
        "previous": rec["prev_symbol"],
        "groups": [{"name": n, "id": i} for n, i in zip(rec["gene_group"], rec["gene_group_id"] + [""] * len(rec["gene_group"]))],
        "entrez_id": rec["entrez_id"],
        "ensembl_gene_id": rec["ensembl_gene_id"],
        "uniprot_ids": rec["uniprot_ids"],
        "omim_ids": rec["omim_id"],
        "refseq": rec["refseq_accession"],
        "mane_select": rec["mane_select"],
        "orphanet": rec["orphanet"],
    }


class KnowledgeBase:
    def __init__(self, remote: RemoteCache):
        self.remote = remote
        self._hgnc: Optional[HgncTable] = None
        self._hgnc_error: Optional[str] = None
        self._lock = threading.Lock()

    # -- HGNC ---------------------------------------------------------------------------
    def hgnc(self, required: bool = True) -> Optional[HgncTable]:
        with self._lock:
            if self._hgnc is None:
                try:
                    self._hgnc = HgncTable(self.remote.get_text(HGNC_URL, ttl=30 * DAY))
                    self._hgnc_error = None
                except RemoteError as exc:
                    self._hgnc_error = str(exc)
            if self._hgnc is None and required:
                raise RemoteError(f"HGNC gene table unavailable: {self._hgnc_error}")
            return self._hgnc

    def hgnc_status(self) -> Dict[str, Any]:
        return {"loaded": self._hgnc is not None, "genes": len(self._hgnc) if self._hgnc else 0, "error": self._hgnc_error, "source": HGNC_URL}

    def gene_names(self, symbols: Iterable[str]) -> Dict[str, Dict[str, str]]:
        """{symbol: {symbol (approved), name, hgnc_id, locus_type}} for the symbols HGNC knows."""
        table = self.hgnc()
        out = {}
        for symbol in symbols:
            rec = table.get(symbol)
            if rec:
                out[symbol] = {"symbol": rec["symbol"], "name": rec["name"], "hgnc_id": rec["hgnc_id"], "locus_type": rec["locus_type"]}
        return out

    def group_genes(self, name: str) -> List[str]:
        group = self.hgnc().groups.get(name)
        if group is None:
            raise KeyError(f"Unknown HGNC gene group: {name}")
        return sorted(group["genes"])

    # -- Gene Ontology (QuickGO) --------------------------------------------------------------
    def go_search(self, query: str, limit: int = 25) -> List[Dict[str, Any]]:
        q = query.strip()
        if not q:
            return []
        if GO_RE.match(q.upper()):
            term = self.go_term(q.upper())
            return [term] if term else []
        data = self.remote.get_json(f"{QUICKGO}/ontology/go/search", ttl=30 * DAY, params={"query": q, "limit": min(limit, 100)})
        return [{"id": r["id"], "name": r.get("name") or "", "aspect": r.get("aspect") or "", "obsolete": bool(r.get("isObsolete"))}
                for r in data.get("results") or [] if not r.get("isObsolete")]

    def go_term(self, go_id: str) -> Optional[Dict[str, Any]]:
        data = self.remote.get_json(f"{QUICKGO}/ontology/go/terms/{go_id}", ttl=30 * DAY)
        results = data.get("results") or []
        if not results:
            return None
        r = results[0]
        return {"id": r["id"], "name": r.get("name") or "", "aspect": r.get("aspect") or "", "definition": (r.get("definition") or {}).get("text") or "", "obsolete": bool(r.get("isObsolete"))}

    def go_genes(self, go_id: str, include_regulation: bool = False) -> Dict[str, Any]:
        """Human genes (reviewed UniProt entries) annotated to a GO term or its descendants."""
        go_id = go_id.strip().upper()
        if not GO_RE.match(go_id):
            raise ValueError(f"Not a GO id: {go_id}")
        relations = "is_a,part_of,occurs_in" + (",regulates,positively_regulates,negatively_regulates" if include_regulation else "")
        text = self.remote.get_text(f"{QUICKGO}/annotation/downloadSearch", ttl=7 * DAY, headers={"Accept": "text/tsv"}, params={
            "goId": go_id, "taxonId": HUMAN_TAXON, "goUsage": "descendants", "goUsageRelationships": relations,
            "geneProductType": "protein", "geneProductSubset": "Swiss-Prot", "selectedFields": "symbol,qualifier,goId,goName",
            "includeFields": "goName", "downloadLimit": 50000,
        })
        genes: Dict[str, Dict[str, Any]] = {}
        table = self.hgnc(required=False)
        for line in text.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) < 4 or parts[1].startswith("NOT"):
                continue
            symbol = parts[0]
            approved = table.resolve(symbol) if table else symbol
            if not approved:
                continue
            entry = genes.setdefault(approved, {"symbol": approved, "terms": {}})
            entry["terms"][parts[2]] = parts[3] if len(parts) > 3 else ""
        rows = sorted(genes.values(), key=lambda g: g["symbol"])
        for row in rows:
            row["terms"] = [{"id": k, "name": v} for k, v in sorted(row["terms"].items())]
        term = self.go_term(go_id)
        return {"term": term, "genes": rows, "relations": relations.split(",")}

    def gene_go(self, symbol: str) -> Dict[str, Any]:
        """GO annotations of a gene (its reviewed UniProt entries), grouped by aspect."""
        table = self.hgnc()
        rec = table.get(symbol)
        if rec is None:
            raise KeyError(f"{symbol} is not an approved HGNC symbol")
        accessions = rec["uniprot_ids"][:3]
        terms: Dict[str, Dict[str, Any]] = {}
        for accession in accessions:
            page, pages = 1, 1
            while page <= min(pages, 5):
                data = self.remote.get_json(f"{QUICKGO}/annotation/search", ttl=7 * DAY, params={
                    "geneProductId": f"UniProtKB:{accession}", "taxonId": HUMAN_TAXON, "includeFields": "goName", "limit": 100, "page": page})
                pages = int((data.get("pageInfo") or {}).get("total") or 1)
                for r in data.get("results") or []:
                    qualifier = r.get("qualifier") or ""
                    term = terms.setdefault(r["goId"], {"id": r["goId"], "name": r.get("goName") or "", "aspect": r.get("goAspect") or "", "qualifiers": set(), "evidence": set()})
                    term["qualifiers"].add(qualifier)
                    if r.get("goEvidence"):
                        term["evidence"].add(r["goEvidence"])
                page += 1
        grouped: Dict[str, List[Dict[str, Any]]] = {aspect: [] for aspect in GO_ASPECTS}
        for term in sorted(terms.values(), key=lambda t: t["name"].lower()):
            negated = all(q.startswith("NOT") for q in term["qualifiers"])
            grouped.setdefault(term["aspect"], []).append({"id": term["id"], "name": term["name"], "qualifiers": sorted(term["qualifiers"]), "evidence": sorted(term["evidence"]), "negated": negated})
        return {"symbol": rec["symbol"], "uniprot_ids": accessions, "aspects": [{"aspect": a, "label": GO_ASPECTS.get(a, a), "terms": grouped.get(a, [])} for a in GO_ASPECTS]}

    # -- ontology terms (OLS4) ------------------------------------------------------------------
    def ontology_term(self, curie: str) -> Dict[str, Any]:
        curie = curie.strip()
        if not CURIE_RE.match(curie):
            raise ValueError(f"Not a CURIE: {curie}")
        data = self.remote.get_json(f"{OLS}/terms", ttl=30 * DAY, params={"obo_id": curie, "size": 20})
        terms = (data.get("_embedded") or {}).get("terms") or []
        prefix = curie.split(":")[0].lower()
        # The defining ontology first (EFO, AISM and others import CL / UBERON terms).
        terms.sort(key=lambda t: (str(t.get("ontology_prefix") or t.get("ontology_name") or "").lower() != prefix, not t.get("is_defining_ontology", False)))
        if not terms:
            return {"curie": curie, "found": False}
        t = terms[0]
        description = t.get("description") or []
        return {
            "curie": curie,
            "found": True,
            "label": t.get("label") or "",
            "definition": description[0] if description else "",
            "synonyms": (t.get("synonyms") or [])[:12],
            "ontology": t.get("ontology_prefix") or t.get("ontology_name") or prefix,
            "iri": t.get("iri") or "",
            "obsolete": bool(t.get("is_obsolete")),
        }
