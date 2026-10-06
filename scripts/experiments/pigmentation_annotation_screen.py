#!/usr/bin/env python3
"""The pigmentation annotation screen, run as an instrument instead of by hand.

WHY THIS EXISTS. The control draws are only interpretable if the screen that clears them
is fixed in advance and on record. For the first draw the screen was run by hand against
QuickGO and its result survives only as prose in
`scripts/experiments/random_gene_genotype_control.py`; for the second draw there is no
record at all, which is an open TODO in the paper's method. A screen nobody can re-run is
not a pre-registration. This script makes it re-runnable and dated.

THE SCREEN, and it is on these terms and no others, so it cannot be widened after seeing a
score:
  GO:0043473  pigmentation
  GO:0042438  melanin biosynthetic process
  GO:0030318  melanocyte differentiation
plus any GO term whose NAME contains pigment, melanin or melanocyte. Ancestors are not
walked: the rule as written is about annotations a gene carries, and inferring through the
ontology graph would be a different and wider screen than the one the draws were cleared
under.

WHAT IT DOES NOT DO. It does not decide membership of a draw. A draw is taken whole --
never substituted, reordered or trimmed -- so a gene that trips the screen is REPORTED,
and the decision about what to do with it is a matter for the paper's text, not for this
script. Nothing here writes a gene list.

THE PANEL IS THE POSITIVE CONTROL. A screen that clears everything is indistinguishable
from a screen that is broken, so the panel genes are screened too and are expected to
trip it. `--check` turns that expectation into an exit code.

Sources: gene symbol -> reviewed human UniProtKB accession (UniProt REST), accession ->
GO annotations (QuickGO annotation search, EBI), GO id -> term name (QuickGO ontology).
Network only; no API key, no cost, no local data.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path("/home/breno/I2CA/genomics")
OUT_DIR = REPO / "results" / "genotype_based_predictor" / "annotation_screen"

EXACT_TERMS = {
    "GO:0043473": "pigmentation",
    "GO:0042438": "melanin biosynthetic process",
    "GO:0030318": "melanocyte differentiation",
}
NAME_RE = re.compile(r"pigment|melanin|melanocyte", re.IGNORECASE)

PANEL = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR", "MFSD12", "OCA2", "HERC2",
         "SLC24A5", "TCHH"]
DRAWS = {
    # random_11_1 and random_11_2 are the two draws already in the paper; the third is the
    # earliest of the three (windows built 2025-12-29) and is unused so far. All three were
    # drawn for the unrelated non-longevous-dataset project, before this question existed.
    1: ["ECHDC3", "EIF1B", "FRA10AC1", "LACTB2", "LRRC36", "PPP1R3E", "PRSS55", "PSMC4",
        "SMCR8", "SPRED2", "TPM2"],
    2: ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D", "LYNX1", "OR11H12", "OR51S1",
        "TRHR", "TSPAN11"],
    3: ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6", "KIAA0319", "RIDA",
        "SEM1", "SFMBT2", "SUMF2"],
}
DRAW_DIR = {
    1: "non_longevous_results_genes_1000_random_11_1",
    2: "non_longevous_results_genes_1000_random_11_2",
    3: "non_longevous_results_genes_1000_random",
}


def get(url: str, tries: int = 4) -> dict:
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:  # noqa: BLE001 -- transient network, retried
            if k == tries - 1:
                raise
            print(f"    retry {k + 1} after {exc}", file=sys.stderr)
            time.sleep(2.0 * (k + 1))
    raise AssertionError("unreachable")


def accession(symbol: str) -> str | None:
    """Reviewed human accession for an exact gene-symbol match, or None."""
    q = urllib.parse.quote(f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true")
    d = get(f"https://rest.uniprot.org/uniprotkb/search?query={q}"
            f"&fields=accession,gene_names&format=json&size=10")
    res = d.get("results", [])
    if not res:
        return None
    # An exact-symbol query can still return more than one entry; prefer the one whose
    # primary gene name matches, rather than taking the first row on faith.
    for r in res:
        for g in r.get("genes", []):
            if g.get("geneName", {}).get("value", "").upper() == symbol.upper():
                return r["primaryAccession"]
    return res[0]["primaryAccession"]


def go_ids(acc: str) -> list[str]:
    """Distinct GO ids annotated to this accession in human, paged to exhaustion."""
    ids, page = set(), 1
    while True:
        d = get("https://www.ebi.ac.uk/QuickGO/services/annotation/search"
                f"?geneProductId={acc}&taxonId=9606&limit=200&page={page}")
        for r in d.get("results", []):
            ids.add(r["goId"])
        total = int(d.get("numberOfHits", 0))
        if page * 200 >= total or not d.get("results"):
            break
        page += 1
    return sorted(ids)


def term_names(ids: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for lo in range(0, len(ids), 100):
        chunk = ids[lo:lo + 100]
        d = get("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
                + ",".join(chunk))
        for r in d.get("results", []):
            out[r["id"]] = r.get("name", "")
    return out


def screen_gene(symbol: str) -> dict:
    acc = accession(symbol)
    if acc is None:
        return {"gene": symbol, "accession": None, "n_terms": 0, "hits": [],
                "status": "NO UNIPROT ENTRY"}
    ids = go_ids(acc)
    names = term_names(ids)
    hits = []
    for i in ids:
        nm = names.get(i, "")
        why = []
        if i in EXACT_TERMS:
            why.append("listed term")
        if NAME_RE.search(nm):
            why.append("name match")
        if why:
            hits.append({"go_id": i, "name": nm, "why": ", ".join(why)})
    return {"gene": symbol, "accession": acc, "n_terms": len(ids), "hits": hits,
            "status": "TRIPS SCREEN" if hits else "clear"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draw", type=int, action="append", choices=sorted(DRAWS),
                    help="repeatable; default: every draw")
    ap.add_argument("--panel", action="store_true", help="screen the panel as well")
    ap.add_argument("--genes", default=None, help="comma-separated symbols instead")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if a control trips, or if no panel gene trips at all")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    groups: dict[str, list[str]] = {}
    if a.genes:
        groups["ad hoc"] = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
    else:
        for d in (a.draw or sorted(DRAWS)):
            groups[f"control draw {d}"] = DRAWS[d]
        if a.panel or a.check:
            groups["panel"] = PANEL

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    record = {"screened_at_utc": stamp,
              "screen": {"exact_terms": EXACT_TERMS, "name_pattern": NAME_RE.pattern,
                         "ancestors_walked": False},
              "sources": {"symbol_to_accession": "UniProt REST (reviewed, taxon 9606)",
                          "annotations": "QuickGO annotation search (EBI)",
                          "term_names": "QuickGO ontology (EBI)"},
              "draw_directories": DRAW_DIR, "groups": {}}

    for label, genes in groups.items():
        print(f"\n=== {label} ===")
        rows = []
        for g in genes:
            r = screen_gene(g)
            rows.append(r)
            flag = "TRIPS" if r["hits"] else "clear"
            print(f"  {g:9s} {str(r['accession']):8s} {r['n_terms']:4d} terms  {flag}")
            for h in r["hits"]:
                print(f"      {h['go_id']} {h['name']}  [{h['why']}]")
        record["groups"][label] = rows

    out = a.out or (OUT_DIR / f"screen_{stamp[:10].replace('-', '')}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=1))
    print(f"\nwrote {out}")

    print("\n--- summary ---")
    for label, rows in record["groups"].items():
        n_trip = sum(1 for r in rows if r["hits"])
        lo = min((r["n_terms"] for r in rows), default=0)
        hi = max((r["n_terms"] for r in rows), default=0)
        print(f"  {label:18s} {n_trip}/{len(rows)} trip the screen; "
              f"{lo}-{hi} annotated terms per gene")

    # Two things are worth reporting and only one of them is a failure.
    #
    # A control that trips is a failure of the screen's premise and has to be seen.
    #
    # A PANEL gene that does not trip is expected, not broken: this screen is a reliable
    # instrument for establishing that a gene has no annotated pigmentation role and an
    # unreliable one for establishing that it has. SLC24A5 causes oculocutaneous albinism
    # type 6 and carries no pigmentation GO annotation at all. So the check requires only
    # that the panel trips SOMEWHERE -- a screen that clears every gene is indistinguishable
    # from a screen that is broken -- and prints the panel misses as the documented
    # limitation they are.
    #
    # A gene with zero annotations is cleared vacuously, which is not the same as being
    # cleared. It is reported separately, because the whole force of a clean draw is that
    # the absence is a real absence and not an empty record.
    if a.check:
        panel = record["groups"].get("panel", [])
        tripping_controls, vacuous = [], []
        for label, rows in record["groups"].items():
            for r in rows:
                if label.startswith("control"):
                    if r["hits"]:
                        tripping_controls.append(
                            f"{r['gene']}: {[h['go_id'] for h in r['hits']]}")
                    if r["n_terms"] == 0:
                        vacuous.append(f"{label} / {r['gene']}")
        panel_miss = [r["gene"] for r in panel if not r["hits"]]
        if panel:
            print(f"\npanel genes the screen does not detect ({len(panel_miss)} of "
                  f"{len(panel)}), expected and not a fault: {', '.join(panel_miss)}")
        if vacuous:
            print(f"cleared with zero annotations, so cleared vacuously -- flag in "
                  f"writing, do not drop: {', '.join(vacuous)}")
        if tripping_controls:
            print("\nCHECK FAILED (a report, not a licence to edit a draw):")
            for b in tripping_controls:
                print(f"  - control gene {b}")
            return 1
        if panel and not any(r["hits"] for r in panel):
            print("\nCHECK FAILED: no panel gene trips the screen, so the screen is "
                  "not measuring anything")
            return 1
        print("\ncheck passed: no control gene carries a pigmentation annotation"
              + (f", and {len(panel) - len(panel_miss)} of {len(panel)} panel genes trip"
                 if panel else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
