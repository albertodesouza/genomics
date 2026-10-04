"""Catalog of every AlphaGenome track: which tissues / cell types (ontology terms) each output has.

Built from ``DnaClient.output_metadata()`` (one cheap call on the hosted API or a self-hosted
server) and cached as JSON, so pickers can offer every ontology term with its track counts::

    python -m genomics.workflows.alphagenome.catalog [--output catalog.json] [--csv tracks.csv]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from genomics.workflows.alphagenome.outputs import OUTPUT_SPECS

CATALOG_VERSION = 1
# Track metadata columns kept per track (the rest is mostly provenance).
TRACK_KEYS = ("name", "strand", "ontology_curie", "biosample_name", "biosample_type", "biosample_life_stage", "gtex_tissue",
              "Assay title", "histone_mark", "transcription_factor", "data_source", "endedness", "genetically_modified")


def _clean(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            value = value.item()
        except (ValueError, AttributeError):
            pass
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def build_catalog(output_metadata: Any, source: str = "") -> Dict[str, Any]:
    """Summarise an ``OutputMetadata`` (any object with per-output DataFrame attributes)."""
    outputs: Dict[str, Dict[str, Any]] = {}
    terms: Dict[str, Dict[str, Any]] = {}
    tracks: Dict[str, List[Dict[str, Any]]] = {}
    for name, spec in OUTPUT_SPECS.items():
        frame = getattr(output_metadata, spec.attr, None)
        if frame is None:
            continue
        try:
            records = frame.to_dict(orient="records")
        except AttributeError:
            continue
        rows = [{k: _clean(r.get(k)) for k in TRACK_KEYS if _clean(r.get(k)) not in (None, "")} for r in records]
        tracks[name] = rows
        outputs[name] = {**spec.as_dict(), "tracks": len(rows), "ontologies": len({r.get("ontology_curie") for r in rows if r.get("ontology_curie")})}
        for row in rows:
            curie = row.get("ontology_curie")
            if not curie:
                continue
            term = terms.setdefault(curie, {"curie": curie, "name": row.get("biosample_name") or "", "type": row.get("biosample_type") or "", "outputs": {}, "marks": {}})
            if not term["name"] and row.get("biosample_name"):
                term["name"] = row["biosample_name"]
            term["outputs"][name] = term["outputs"].get(name, 0) + 1
            mark = row.get("transcription_factor") or row.get("histone_mark")
            if mark:
                term["marks"].setdefault(name, [])
                if mark not in term["marks"][name]:
                    term["marks"][name].append(mark)
    for term in terms.values():
        term["tracks"] = sum(term["outputs"].values())
        if not term["marks"]:
            term.pop("marks")
    return {
        "version": CATALOG_VERSION,
        "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "source": source,
        "outputs": outputs,
        "ontologies": sorted(terms.values(), key=lambda t: (t["name"].lower(), t["curie"])),
        "tracks": tracks,
    }


def fetch_catalog(client: Any, source: str = "") -> Dict[str, Any]:
    from alphagenome.models import dna_client

    return build_catalog(client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS), source)


def load_catalog(path: Path) -> Optional[Dict[str, Any]]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) and data.get("version") == CATALOG_VERSION and data.get("ontologies") else None


def save_catalog(catalog: Dict[str, Any], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(catalog, separators=(",", ":")), encoding="utf-8")
    os.replace(tmp, path)


def summary(catalog: Dict[str, Any]) -> Dict[str, Any]:
    """The catalog without per-track rows (what pickers need)."""
    return {k: v for k, v in catalog.items() if k != "tracks"}


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="genomics alphagenome catalog", description="List every AlphaGenome output track and ontology term")
    parser.add_argument("--output", type=Path, default=Path("alphagenome_catalog.json"), help="JSON catalog to write")
    parser.add_argument("--csv", type=Path, default=None, help="Also write one row per track")
    args = parser.parse_args(argv)
    from genomics.core.alphagenome_connection import create_dna_client

    catalog = fetch_catalog(create_dna_client(timeout=60.0), source=os.environ.get("ALPHAGENOME_ADDRESS") or "hosted API")
    save_catalog(catalog, args.output)
    if args.csv:
        import csv

        with open(args.csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["output_type", *TRACK_KEYS])
            writer.writeheader()
            for name, rows in catalog["tracks"].items():
                for row in rows:
                    writer.writerow({"output_type": name, **row})
    counts = ", ".join(f"{name.lower()} {o['tracks']}" for name, o in catalog["outputs"].items())
    print(f"{len(catalog['ontologies'])} ontology terms; tracks: {counts} -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
