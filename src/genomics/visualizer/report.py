"""Gene report: what the reference genome makes of a gene, and what each of an individual's
haplotypes changes, transcript by transcript.

The pipeline, for one gene, one sample and one tissue:

1. **Reference baseline** (observed): the gene's level in the tissue (Human Protein Atlas
   single-cell nCPM for melanocytes, GTEx v8 median TPM for lymphoblastoid cells) split over its
   transcripts by GTEx v8's median transcript TPMs (bulk sun-exposed skin stands in for
   melanocytes, which GTEx does not profile). Each transcript gets a share and an absolute level.
2. **Each haplotype** (predicted): the transcript's level on that copy =
   ``reference level / 2 x gene fold x splicing shift x NMD``, where the gene fold is AlphaGenome's
   RNA-seq of the haplotype against the reference window, the splicing shift is the change in
   AlphaGenome usage of the transcript's weakest junction (shares renormalised over the gene), and
   an mRNA that gains a premature stop codon triggering nonsense-mediated decay keeps ``NMD_RESIDUAL``
   of its level (typical steady-state residue; the exact value varies by transcript).
3. **Products**: per transcript and haplotype, the protein change (combined effect of all the
   haplotype's variants), whether the mRNA is unstable (NMD), whether it carries a premature stop,
   and a protein similarity score against the reference protein (BLOSUM62 alignment, identity).
4. **Variants**: every variant the haplotype carries in the gene, each applied alone to the
   reference genome and classified on every coding transcript (synonymous, missense, nonsense,
   frameshift, splice, UTR, intronic), with AlphaMissense pathogenicity for missense changes of
   the base transcript (from the AlphaFold Database's per-protein table).
"""
from __future__ import annotations

import csv
import io
from typing import Any, Callable, Dict, List, Optional

import numpy as np

from genomics.visualizer.products import (HAPLOTYPES, REFERENCE, Junctions, _incomplete, build_product, coding_start_stop, introns_of,
                                          transcript_support, variant_effects)
from genomics.visualizer.remote import DAY
from genomics.visualizer.structures import AFDB_API, compare_proteins, protein_similarity

ProgressFn = Callable[[float, str], None]

NMD_RESIDUAL = 0.2  # share of an NMD-targeted mRNA left at steady state (typically 5-30%)
SUPPORT_FLOOR = 0.02  # junction usage below which a transcript's splicing shift is not estimated
MAX_SHIFT = 5.0
SHIFT_PSEUDOCOUNT = 0.05
AA_ONE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
          "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
AM_CLASSES = {"LBen": "likely benign", "Amb": "ambiguous", "LPath": "likely pathogenic"}
UNSTABLE = "unstable"


class ReportService:
    def __init__(self, app):
        self.app = app
        self._memory: Dict[tuple, Dict[str, Any]] = {}
        self._am: Dict[str, Dict[str, Any]] = {}

    def key(self, dataset, sample, gene, target, tissue: str) -> str:
        return f"report:{dataset.id}:{sample}:{gene}:{target or ''}:{tissue}"

    def cached(self, dataset, sample, gene, target, tissue) -> Optional[Dict[str, Any]]:
        return self._memory.get((dataset.fingerprint, sample, gene, target, tissue))

    # -- databases ---------------------------------------------------------------------------
    def _transcript_tpms(self, symbol: str, source: Optional[Dict[str, Any]], transcripts: List[Dict[str, Any]], window_start: int) -> Dict[str, Any]:
        """GTEx median TPM per current transcript: by Ensembl id, else by an identical intron chain.

        GTEx v8 uses GENCODE v26; transcripts annotated later (often the MANE Select one) carry new
        ids but usually the same intron chain as a v26 transcript, whose TPM they inherit.
        """
        if source is None:
            return {"values": {}, "error": "no transcript-level reference for this tissue"}
        try:
            gene = self.app.gtex.gene(symbol)
            if not gene or str(gene.get("geneSymbol", "")).upper() != symbol.upper():
                return {"values": {}, "error": f"{symbol} not found in GTEx"}
            tpm = self.app.gtex.median_transcripts(gene["gencodeId"], source["tissue"])
        except Exception as exc:  # offline or GTEx down: the baseline says so
            return {"values": {}, "error": str(exc)}
        values: Dict[str, Optional[float]] = {}
        matched_by: Dict[str, str] = {}
        claimed = set()
        for t in transcripts:
            key = t["id"].split(".")[0]
            if key in tpm:
                values[t["id"]] = tpm[key]
                claimed.add(key)
        missing = [t for t in transcripts if t["id"] not in values]
        if missing:
            try:
                chains = self.app.gtex.transcript_introns(gene["gencodeId"])
            except Exception:
                chains = {}
            by_chain: Dict[tuple, List[str]] = {}
            for tid, chain in chains.items():
                if chain and tid in tpm and tid not in claimed:
                    by_chain.setdefault(tuple(chain), []).append(tid)
            for t in sorted(missing, key=lambda t: (not t.get("mane"), not t.get("canonical"))):
                exons = sorted((int(a), int(b)) for a, b in t["exons"])
                chain = tuple((window_start + a[1] - 1, window_start + b[0]) for a, b in zip(exons, exons[1:]))
                candidates = [c for c in by_chain.get(chain, []) if c not in claimed]
                if chain and candidates:
                    values[t["id"]] = tpm[candidates[0]]
                    matched_by[t["id"]] = candidates[0]
                    claimed.add(candidates[0])
        return {"values": values, "matched_by": matched_by, "source": source["label"], "proxy": source["proxy"]}

    def alphamissense(self, transcript_id: str, protein: str) -> Dict[str, Any]:
        """{"R151C": (score, class)} for the transcript's protein when its UniProt sequence is identical."""
        hit = self._am.get(transcript_id)
        if hit is not None:
            return hit
        out: Dict[str, Any] = {"scores": {}, "accession": None}
        try:
            entry = self.app.structures.uniprot(transcript_id, protein)
            if entry is None or entry["sequence"] != protein:
                out["error"] = "no UniProt entry with exactly this protein"
            else:
                out["accession"] = entry["accession"]
                models = self.app.remote.get_json(AFDB_API.format(accession=entry["accession"]), ttl=90 * DAY)
                url = (models or [{}])[0].get("amAnnotationsUrl")
                if not url:
                    out["error"] = "AlphaMissense not available for this protein"
                else:
                    text = self.app.remote.get_text(url, ttl=180 * DAY)
                    for row in csv.DictReader(io.StringIO(text)):
                        out["scores"][row["protein_variant"]] = (float(row["am_pathogenicity"]), AM_CLASSES.get(row["am_class"], row["am_class"]))
        except Exception as exc:  # offline, no entry, service down: reported, never fatal
            out["error"] = str(exc)
        self._am[transcript_id] = out
        return out

    # -- the report ----------------------------------------------------------------------------
    def compute(self, dataset, sample: str, gene: str, models: List[Dict[str, Any]], target: Optional[str], tissue_info: Dict[str, Any], progress: ProgressFn) -> Dict[str, Any]:
        """``tissue_info``: a :func:`genomics.visualizer.expression.tissue_spec`."""
        app = self.app
        tissue = tissue_info["ontology"]
        progress(0.02, "Gene products")
        payload = app.products.products(dataset, sample, gene, models, target=target)
        if not payload.get("transcripts"):
            return {"available": False, "message": payload.get("message") or "No transcripts"}
        symbol = payload["target"]
        chosen = next((g for g in models if g["name"] == symbol), {})
        hit = app.expression.cached(dataset, sample, gene, target, [tissue_info])
        expression = hit if hit is not None else app.expression.compute(dataset, sample, gene, models, target, lambda f, m: progress(0.05 + 0.5 * f, m), tissues=[tissue_info])
        trow = next((t for t in expression.get("tissues", []) if t["ontology"] == tissue_info["ontology"]), {}) if expression.get("available") else {}
        progress(0.6, "GTEx transcript expression")
        tpms = self._transcript_tpms(symbol, tissue_info.get("transcripts"), payload["transcripts"], int(payload["window_start"] or 1))
        anchor = trow.get("anchor") or app.expression.references.level(tissue_info.get("anchor"), symbol, chosen.get("id"))
        level = anchor.get("value")
        unit = anchor.get("unit")
        transcripts = payload["transcripts"]
        splicing = payload.get("splicing") or {}
        column = next((t["index"] for t in splicing.get("tracks") or [] if t.get("ontology") == tissue), None)
        support = splicing.get("support") or {}
        junction_source = "stored" if column is not None else None
        if column is None:  # a tissue without stored junction predictions: use the ones predicted with its RNA-seq
            predicted = app.expression.junctions(dataset, sample, gene, tissue)
            if predicted and all(h in predicted for h in (REFERENCE,) + HAPLOTYPES):
                views_ = app.products.views(dataset, sample, gene)
                sets = {h: Junctions.from_arrays(arrays, meta, views_[h], payload["strand"]) for h, (arrays, meta) in predicted.items()}
                tracks = sets[REFERENCE].tracks
                cols = [i for i, r in enumerate(tracks) if r.get("ontology_curie") == tissue] or list(range(len(tracks)))
                full = transcript_support(transcripts, sets, payload["strand"])
                support = {tid: {h: [float(np.mean(v[cols]))] for h, v in per.items()} for tid, per in full.items()}
                column = 0
                junction_source = "on_demand"

        # 1. reference baseline
        values = tpms.get("values") or {}
        matched = {t["id"]: values.get(t["id"]) for t in transcripts}
        total = sum(v for v in matched.values() if v)
        baseline = []
        for t in transcripts:
            v = matched[t["id"]]
            share = (v / total) if (v is not None and total > 0) else None
            baseline.append({
                "id": t["id"], "name": t["name"], "type": t["type"], "mane": t["mane"], "coding": bool(t["products"]["ref"].get("coding")),
                "protein_length": len(t["products"]["ref"].get("protein") or "") or None, "gtex_tpm": v, "share": share,
                "level": level * share if (share is not None and level is not None) else None,
                "in_gtex": v is not None, "nmd_biotype": t["type"] == "nonsense_mediated_decay",
                "gtex_match": (tpms.get("matched_by") or {}).get(t["id"]),
            })

        # 2-3. haplotypes
        mrna = trow.get("mrna") or {}
        gene_fold = {h: (mrna.get("fold") or {}).get(h) for h in HAPLOTYPES}
        base = next(t for t in transcripts if t["id"] == payload["base"])
        am = self.alphamissense(base["id"], base["products"]["ref"].get("protein") or "") if base["products"]["ref"].get("protein") else {"scores": {}}
        progress(0.75, "Variant consequences")
        views = app.products.views(dataset, sample, gene)
        model_tx = {t["id"]: t for t in chosen.get("transcripts", [])}
        full = [model_tx[t["id"]] for t in transcripts if t["id"] in model_tx]  # with tags, CDS phase, strand
        reference_products = {}
        for t in full:
            start, stop = coding_start_stop(t)
            inc = _incomplete(t)
            reference_products[t["id"]] = build_product(views[REFERENCE], t["exons"], t["strand"], start, stop, inc["start"], with_sequence=True, phase=t.get("cds_phase", 0))
        variants = variant_effects(views[REFERENCE], full, model_tx.get(base["id"], {**base, "strand": payload["strand"]}), payload["strand"], reference_products,
                                   app.signals.variants(dataset, sample, gene), int(payload["window_start"] or 1))
        for v in variants:
            bc = v.get("base_change")
            if bc and bc.get("class") == "missense":
                key = _one_letter(bc.get("hgvs", ""))
                score = am["scores"].get(key) if key else None
                v["alphamissense"] = {"variant": key, "score": score[0], "class": score[1]} if score else None
        haplotypes = {}
        for h in HAPLOTYPES:
            rows = []
            shifted = {}
            for t, b in zip(transcripts, baseline):
                if b["share"] is None:
                    continue
                factor = 1.0
                s = support.get(t["id"])
                if s and column is not None and introns_of(t["exons"]):
                    ref_u, hap_u = s[REFERENCE][column], (s.get(h) or s[REFERENCE])[column]
                    if ref_u >= SUPPORT_FLOOR:
                        # pseudocount: a ratio of two small usages is mostly noise
                        factor = min(MAX_SHIFT, (hap_u + SHIFT_PSEUDOCOUNT) / (ref_u + SHIFT_PSEUDOCOUNT))
                shifted[t["id"]] = b["share"] * factor
            norm = sum(shifted.values())
            for t, b in zip(transcripts, baseline):
                ref_p = t["products"]["ref"]
                p = t["products"][h]
                change = p.get("change") or {}
                nmd_gained = bool((p.get("nmd") or {}).get("predicted")) and not bool((ref_p.get("nmd") or {}).get("predicted"))
                share_h = (shifted[t["id"]] / norm) if (t["id"] in shifted and norm > 0) else None
                copy_level = None
                fold = None
                if share_h is not None and gene_fold[h] is not None:
                    rel = gene_fold[h] * (share_h / b["share"] if b["share"] else 1.0) * (NMD_RESIDUAL if nmd_gained else 1.0)
                    fold = rel
                    copy_level = (level / 2) * b["share"] * rel if level is not None else None
                ptc = change.get("class") in ("stop_gained", "frameshift") and p.get("stop_found") and len(p.get("protein") or "") < len(ref_p.get("protein") or "")
                similarity = protein_similarity(ref_p.get("protein") or "", p.get("protein") or "") if ref_p.get("coding") else None
                subs = compare_proteins(ref_p.get("protein") or "", p.get("protein") or "").get("substitutions") if ref_p.get("coding") else None
                if subs and t["id"] == base["id"]:
                    for sub in subs:
                        score = am["scores"].get(f"{sub['ref']}{sub['pos']}{sub['alt']}")
                        sub["alphamissense"] = {"score": score[0], "class": score[1]} if score else None
                rows.append({
                    "id": t["id"], "name": t["name"], "share": share_h, "level": copy_level, "fold": fold,
                    "change": change, "nmd": p.get("nmd"), "nmd_gained": nmd_gained, "premature_stop": bool(ptc),
                    "protein_length": len(p.get("protein") or "") or None, "similarity": similarity, "substitutions": subs,
                    "start_lost": bool(p.get("start_lost")) and not bool(ref_p.get("start_lost")),
                })
            carried = [v for v in variants if h in v["haplotypes"]]
            counts: Dict[str, int] = {}
            for v in carried:
                counts[v["category"]] = counts.get(v["category"], 0) + 1
            base_row = next(r for r in rows if r["id"] == base["id"])

            def flagged(key: str) -> List[Dict[str, Any]]:
                # every affected transcript, with its share so the page can say how much it matters
                return [{"name": r["name"], "share": r["share"], "reference_share": next(b["share"] for b in baseline if b["id"] == r["id"]),
                         "change": r["change"], "base": r["id"] == base["id"]} for r in rows if r[key]]

            haplotypes[h] = {
                "gene_fold": gene_fold[h],
                "level": sum(r["level"] for r in rows if r["level"] is not None) if level is not None else None,
                "transcripts": rows,
                "variants": carried,
                "counts": counts,
                "flags": {
                    "unstable": flagged("nmd_gained"),
                    "premature_stop": flagged("premature_stop"),
                    "start_lost": flagged("start_lost"),
                    "missense": base_row["substitutions"] or [],
                    "base_change": base_row["change"],
                    "base_similarity": base_row["similarity"],
                },
            }
        progress(1.0, "Report ready")
        result = {
            "available": True, "sample": sample, "gene": gene, "target": symbol, "base": base["id"], "base_name": base["name"], "chromosome": payload.get("chromosome"),
            "tissue": {"ontology": tissue, "label": tissue_info["label"], "group": tissue_info["group"], "type": tissue_info.get("type"),
                       "gtex_tissue": tissue_info.get("gtex_tissue"), "junctions": junction_source,
                       "anchor": (tissue_info.get("anchor") or {}).get("label"), "anchor_matched_by_name": bool((tissue_info.get("anchor") or {}).get("matched_by_name"))},
            "level": {"value": level, "unit": unit, "source": anchor.get("source"), "error": anchor.get("error"), "expressed": (level or 0) >= 1.0 if level is not None else None},
            "transcript_source": {"label": tpms.get("source"), "proxy": tpms.get("proxy"), "error": tpms.get("error"),
                                  "unmatched": [b["name"] for b in baseline if not b["in_gtex"]]},
            "expression_source": trow.get("source"), "expression_problems": expression.get("problems") if expression.get("available") else [expression.get("message")],
            "baseline": baseline, "haplotypes": haplotypes,
            "alphamissense": {"accession": am.get("accession"), "error": am.get("error")},
            "notes": {
                "nmd_residual": NMD_RESIDUAL,
                "pipeline": __doc__.split("\n\n", 1)[1].strip() if __doc__ else "",
            },
        }
        self._memory[(dataset.fingerprint, sample, gene, target, tissue)] = result
        return result


def _one_letter(hgvs: str) -> Optional[str]:
    """``p.Arg151Cys`` -> ``R151C`` (single substitutions only)."""
    import re

    m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvs or "")
    if not m or m.group(1) not in AA_ONE or m.group(3) not in AA_ONE:
        return None
    return f"{AA_ONE[m.group(1)]}{m.group(2)}{AA_ONE[m.group(3)]}"
