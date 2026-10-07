#!/usr/bin/env python3
"""Check the visualizer's gene products (genomics.visualizer.products) against independent tools.

1. Reference proteins: every complete coding transcript translated from the reference window is
   compared with Ensembl's peptide for the same transcript (REST ``/sequence/id``). Needs network;
   ``--skip-ensembl`` turns it off. Transcripts whose Ensembl version differs from the GTF are
   reported separately (the sequence may have changed since the annotation release).
2. Haplotype changes: ``bcftools csq`` (haplotype-aware, ``-p a``) is run on the sample's window
   VCF with an Ensembl GFF3 of the same release as the GTF cache (GENCODE v46 = Ensembl 112), and
   its protein consequences per transcript and haplotype are compared with ours.

Usage::

    python3 scripts/diagnostics/validate_gene_products.py --sample HG00096 \\
        --gff3-dir /path/to/ensembl112 --reference /path/to/GRCh38.fa [--genes TYR,OCA2] [--json out.json]

The GFF3 directory holds ``Homo_sapiens.GRCh38.112.chromosome.<N>.gff3.gz`` files from
https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from genomics.core.data_registry import resolve_dataset
from genomics.visualizer.annotations import AnnotationService
from genomics.visualizer.cache import DiskArrayCache
from genomics.visualizer.datasets import DatasetCatalog
from genomics.visualizer.products import THREE_LETTER, ProductService
from genomics.visualizer.signals import SignalService

ONE_LETTER = {v: k for k, v in THREE_LETTER.items()}
PROTEIN_ALTERING = {"missense", "stop_gained", "frameshift", "inframe_deletion", "inframe_insertion", "inframe_altering", "start_lost", "stop_lost", "splice_donor", "splice_acceptor"}
CLASS_EQUIVALENTS = {
    "missense": {"missense"},
    "stop_gained": {"stop_gained"},
    "frameshift": {"frameshift"},
    "inframe_indel": {"inframe_deletion", "inframe_insertion", "inframe_altering"},
    "start_lost": {"start_lost"},
    "stop_lost": {"stop_lost"},
}


def _post(url: str, payload: Dict[str, Any], attempts: int = 5) -> Any:
    body = json.dumps(payload).encode("utf-8")
    for attempt in range(attempts):
        req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json", "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:  # Ensembl answers 429/503 under load
            if exc.code not in (429, 500, 502, 503, 504) or attempt == attempts - 1:
                raise
            time.sleep(2 ** attempt)


def ensembl_proteins(transcript_ids: List[str], batch_size: int = 10) -> Dict[str, Dict[str, str]]:
    """{versionless transcript id: {"seq", "version"}} from the Ensembl REST API."""
    out: Dict[str, Dict[str, str]] = {}
    for i in range(0, len(transcript_ids), batch_size):
        batch = transcript_ids[i:i + batch_size]
        for item in _post("https://rest.ensembl.org/sequence/id", {"ids": batch, "type": "protein"}):
            out[item["query"]] = {"seq": item["seq"]}
        for tid, item in _post("https://rest.ensembl.org/lookup/id", {"ids": batch}).items():
            if item and tid in out:
                out[tid]["version"] = str(item.get("version"))
    return out


def run_csq(vcf: Path, gff3: Path, reference: Path) -> Dict[Tuple[str, int], List[Dict[str, Any]]]:
    """{(versionless transcript, haplotype 1|2): [consequences]} from ``bcftools csq -Ot``."""
    cmd = ["bcftools", "csq", "-p", "a", "--unify-chr-names", "chr,-,chr", "-f", str(reference), "-g", str(gff3), "-Ot", str(vcf)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools csq failed: {proc.stderr.strip()[:500]}")
    out: Dict[Tuple[str, int], List[Dict[str, Any]]] = defaultdict(list)
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 6 or parts[0] != "CSQ" or parts[5].startswith("@"):
            continue
        fields = parts[5].split("|")
        if len(fields) < 5:
            continue
        classes = set(fields[0].lstrip("*").split("&"))
        tid = fields[2]
        out[(tid, int(parts[2]))].append({"classes": classes, "aa": fields[5] if len(fields) > 5 else "", "dna": fields[6] if len(fields) > 6 else "", "pos": int(parts[4])})
    return out


def csq_substitutions(records: List[Dict[str, Any]]) -> Set[Tuple[int, str, str]]:
    subs = set()
    for r in records:
        if "missense" not in r["classes"]:
            continue
        m = re.match(r"^(\d+)([A-Z*]+)>(\d+)([A-Z*]+)$", r["aa"])
        if m and len(m.group(2)) == len(m.group(4)):
            start = int(m.group(1))
            for k, (a, b) in enumerate(zip(m.group(2), m.group(4))):
                if a != b:
                    subs.add((start + k, a, b))
    return subs


def our_substitutions(ref: str, alt: str) -> Set[Tuple[int, str, str]]:
    return {(i + 1, a, b) for i, (a, b) in enumerate(zip(ref, alt)) if a != b} if len(ref) == len(alt) else set()


def compare(our_class: str, ref_protein: str, alt_protein: str, records: List[Dict[str, Any]]) -> Tuple[bool, str]:
    classes = set().union(*[r["classes"] for r in records]) if records else set()
    altering = classes & PROTEIN_ALTERING - {"splice_donor", "splice_acceptor"}
    if our_class in ("no_change", "utr_change", "noncoding_change"):
        return (not altering and "synonymous" not in classes), f"csq: {sorted(classes) or '-'}"
    if our_class == "synonymous":
        return (not altering and "synonymous" in classes), f"csq: {sorted(classes)}"
    if our_class == "missense":
        ours, theirs = our_substitutions(ref_protein, alt_protein), csq_substitutions(records)
        return ours == theirs and not (altering - {"missense"}), f"ours {sorted(ours)} csq {sorted(theirs)} {sorted(altering)}"
    wanted = CLASS_EQUIVALENTS.get(our_class, {our_class})
    return bool(altering & wanted), f"csq: {sorted(classes)}"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dataset-id", default="1kg_high_coverage")
    parser.add_argument("--dataset-dir", type=Path, help="Overrides --dataset-id")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genes", help="Comma-separated windows (default: all)")
    parser.add_argument("--gff3-dir", type=Path, help="Ensembl GFF3 per chromosome; without it the bcftools csq check is skipped")
    parser.add_argument("--reference", type=Path, help="Reference FASTA (chr-prefixed, indexed) for bcftools csq")
    parser.add_argument("--release", default="112")
    parser.add_argument("--skip-ensembl", action="store_true")
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--json", type=Path, help="Write the full report here")
    args = parser.parse_args(argv)

    path = args.dataset_dir or resolve_dataset(args.dataset_id).path
    dataset = DatasetCatalog().add(Path(path))
    signals = SignalService(1 << 30, DiskArrayCache(args.cache_dir))
    annotations = AnnotationService(args.cache_dir)
    models = annotations.build(dataset, lambda f, m: None)["genes"]
    service = ProductService(signals, annotations)
    genes = args.genes.split(",") if args.genes else list(dataset.genes)

    report: Dict[str, Any] = {"sample": args.sample, "genes": {}, "ensembl": {"match": 0, "mismatch": [], "version_differs": [], "missing": []},
                              "csq": {"agree": 0, "disagree": [], "skipped": []}}
    payloads = {}
    for gene in genes:
        if not dataset.gene_dir(args.sample, gene).is_dir():
            report["csq"]["skipped"].append(f"{gene}: no window for {args.sample}")
            continue
        payloads[gene] = service.products(dataset, args.sample, gene, models.get(gene, []))

    if not args.skip_ensembl:
        complete = {}
        for gene, p in payloads.items():
            for tx in p["transcripts"]:
                if tx["products"]["ref"].get("coding") and not tx["incomplete"]["start"] and not tx["incomplete"]["end"]:
                    complete[tx["id"].split(".")[0]] = (gene, tx)
        fetched = ensembl_proteins(sorted(complete))
        for tid, (gene, tx) in complete.items():
            got = fetched.get(tid)
            ours = tx["products"]["ref"].get("protein") or ""
            if got is None:
                report["ensembl"]["missing"].append(tx["id"])
            elif got["seq"] == ours:
                report["ensembl"]["match"] += 1
            elif got.get("version") and got["version"] != tx["id"].split(".")[-1]:
                report["ensembl"]["version_differs"].append(f"{tx['id']} (Ensembl now v{got['version']})")
            else:
                report["ensembl"]["mismatch"].append({"transcript": tx["id"], "name": tx["name"], "ours": len(ours), "ensembl": len(got["seq"])})
        e = report["ensembl"]
        print(f"[ensembl] reference proteins: {e['match']} identical, {len(e['mismatch'])} different, {len(e['version_differs'])} newer Ensembl version, {len(e['missing'])} not found")
        for m in e["mismatch"]:
            print(f"  MISMATCH {m}")

    if args.gff3_dir and args.reference:
        for gene, p in payloads.items():
            chrom = str(p["chromosome"]).replace("chr", "")
            gff3 = args.gff3_dir / f"Homo_sapiens.GRCh38.{args.release}.chromosome.{chrom}.gff3.gz"
            vcf = dataset.sample_vcf_path(args.sample, gene)
            if not gff3.exists() or vcf is None:
                report["csq"]["skipped"].append(f"{gene}: {'no GFF3 for chr' + chrom if not gff3.exists() else 'no VCF'}")
                continue
            csq = run_csq(vcf, gff3, args.reference)
            for tx in p["transcripts"]:
                if not tx["products"]["ref"].get("coding"):
                    continue
                tid = tx["id"].split(".")[0]
                for hap_index, hap in ((1, "H1"), (2, "H2")):
                    product = tx["products"][hap]
                    ok, detail = compare(product["change"]["class"], tx["products"]["ref"].get("protein") or "", product.get("protein") or "", csq.get((tid, hap_index), []))
                    if ok:
                        report["csq"]["agree"] += 1
                    else:
                        report["csq"]["disagree"].append({"gene": gene, "transcript": tx["name"], "haplotype": hap, "ours": product["change"], "detail": detail})
        c = report["csq"]
        print(f"[bcftools csq] transcript x haplotype calls: {c['agree']} agree, {len(c['disagree'])} disagree, {len(c['skipped'])} windows skipped")
        for d in c["disagree"]:
            print(f"  DISAGREE {d['gene']} {d['transcript']} {d['haplotype']}: ours {d['ours'].get('class')} {d['ours'].get('hgvs')} | {d['detail']}")
    else:
        print("[bcftools csq] skipped (needs --gff3-dir and --reference)")

    for gene, p in payloads.items():
        altered = []
        for tx in p["transcripts"]:
            for hap in ("H1", "H2"):
                ch = tx["products"][hap]["change"]
                if ch["class"] not in ("no_change", "utr_change", "synonymous", "noncoding_change"):
                    altered.append(f"{tx['name']} {hap} {ch['class']} {ch.get('hgvs', '')}")
        report["genes"][gene] = altered
    if args.json:
        args.json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 1 if report["ensembl"]["mismatch"] or report["csq"]["disagree"] else 0


if __name__ == "__main__":
    sys.exit(main())
