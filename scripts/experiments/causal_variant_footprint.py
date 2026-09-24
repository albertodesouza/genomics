#!/usr/bin/env python3
"""X-B: the transcriptional footprint of fine-mapped pigmentation variants.

WHAT THIS ANSWERS
-----------------
The knockdown probe finds that OCA2 and HERC2 do not move the classifier's
decision at all (|Delta| = 0.002 and 0.018), and the validity instrumentation
shows this is not a delivery failure for OCA2: it receives 1.5x MC1R's delivered
perturbation and moves the decision 3% as much. So the classifier does not use
OCA2. That is a statement about the *classifier*, and it leaves the more
interesting question open: is the locus unused because the classifier ignores it,
or because the frozen model never encoded anything there in the first place?

This script asks the second question directly, and it does not need the
classifier at all. For each fine-mapped pigmentation variant it toggles that
single base against the reference and measures how much the frozen model's
predicted RNA-seq changes. A variant the frozen model is blind to cannot reach
any downstream classifier no matter how good the classifier is.

THE PRE-REGISTRATION IS ALREADY IN THE REPOSITORY
-------------------------------------------------
`snps.csv` classifies each variant before any measurement:

  expected_rnaseq_signal = NULL_CODING   a coding change with no reason to move
                                         RNA-seq; the footprint should be small
  expected_rnaseq_signal = VISIBLE       a regulatory variant; the footprint
                                         should be large

and `control_class` in {NEG, POS, SUPP}. That column was written from the
literature, not from any prediction, so it is a genuine pre-registration and this
script is scored against it rather than interpreted after the fact.

WHAT THE OUTCOMES MEAN
----------------------
The panel's two most-cited variants make opposite predictions and both bear on
the paper's flat genes:

  rs1426654  (SLC24A5, A111T)   NULL_CODING -- a *missense* variant. If its
             footprint is small, then the single largest-effect pigmentation
             allele in the genome is invisible to the frozen model, and the
             classifier's strong SLC24A5 response is therefore NOT that variant.
  rs12913832 (HERC2 -> OCA2)    VISIBLE -- but it acts through a long-range
             enhancer loop, the regime `karollus2023current` identifies as the
             worst case for these models. If its footprint is small too, that is
             a mechanistic account of the OCA2/HERC2 nulls: the causal variant is
             invisible to the frozen model, so the classifier has nothing to use.

Either result converts "OCA2 is unused by the classifier" into a statement about
the frozen model's blind spot, which is a more useful finding than the one the
knockdown alone supports. A high footprint for the NULL_CODING variants would
instead say the frozen model responds to coding changes for reasons unrelated to
transcription, which would undercut reading its output as regulatory at all.

METHOD
------
Coordinates are resolved from Ensembl REST (GRCh38) at run time and cached, so
no coordinate is hardcoded here. For each variant:

  1. take the same 512 kbp reference window the pipeline uses, centred on the
     gene midpoint, so the crop geometry matches the classifier's exactly;
  2. predict the reference window;
  3. substitute the alternate allele at the variant position, predict again;
  4. report the change over the full window and, separately, over the central
     32,768 units the classifier actually reads -- the latter is directly
     comparable to |Delta_in| from the scramble knockdown.

Only two calls per variant, and the reference call is shared by every variant in
the same window, so the whole panel costs well under a hundred calls.

Usage:
  python3 scripts/experiments/causal_variant_footprint.py
  ... --panel-only          # only variants in the eleven classifier windows
  ... --rsids rs1426654,rs12913832
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
os.chdir(REPO_ROOT)

REF_FASTA = Path("/dados/GENOMICS_DATA/top3/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa")
SAMTOOLS = os.environ.get("SAMTOOLS_BIN", "/home/breno/miniforge3/envs/genomics/bin/samtools")
GTF_CACHE = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all/gtf_cache.feather")
SNPS_CSV = REPO_ROOT / "snps.csv"
WINDOW_SIZE = 524288
CROP = 32768
ONTOLOGIES = ["CL:1000458", "CL:0000346", "CL:2000092"]

OUT_PATH = REPO_ROOT / "results" / "genotype_based_predictor" / "causal_variant_footprint.json"
COORD_CACHE = REPO_ROOT / "results" / "genotype_based_predictor" / "rsid_coordinates_grch38.json"

# The eleven windows the classifier reads. A variant outside these cannot reach
# the classifier at all, which is worth separating from a variant the frozen
# model is blind to.
PANEL_GENES = {
    "SLC24A5", "TYR", "SLC45A2", "TYRP1", "DDB1", "MFSD12",
    "MC1R", "HERC2", "OCA2", "EDAR", "TCHH",
}

# Realised delivered perturbation from the 100 bp promoter scramble, for scale.
# The variant footprints below are directly comparable to these numbers: both are
# total absolute change over the same crop, same ontologies, same strands.
SCRAMBLE_DELTA_IN = {
    "TYR": 54296.6, "DDB1": 34150.3, "TYRP1": 26589.5, "SLC45A2": 23730.6,
    "MFSD12": 11868.3, "SLC24A5": 6537.0, "OCA2": 2109.0, "MC1R": 1391.0,
    "HERC2": 798.6, "EDAR": 59.0, "TCHH": 53.0,
}


def _log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc).isoformat()}] {msg}", flush=True)


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------

def _ensembl_variant(rsid: str, attempts: int = 4) -> dict | None:
    """GRCh38 mapping for an rsID. Returns None rather than raising, so one dead
    rsID does not abort a panel run."""
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "genomics-xb/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            break
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt == attempts:
                _log(f"  {rsid}: Ensembl lookup failed after {attempts} attempts -- {exc}")
                return None
            time.sleep(2 ** attempt)
    else:                                                          # pragma: no cover
        return None

    for m in payload.get("mappings", []):
        if m.get("assembly_name") != "GRCh38":
            continue
        chrom = str(m.get("seq_region_name", ""))
        # Skip patches and scaffolds; only primary chromosomes are in our reference.
        if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}:
            continue
        return {
            "chrom": f"chr{chrom}",
            "pos": int(m["start"]),                # 1-based
            "allele_string": m.get("allele_string"),
            "ancestral": payload.get("ancestral_allele") or m.get("ancestral_allele"),
            "most_severe_consequence": payload.get("most_severe_consequence"),
        }
    return None


def _load_coords(rsids: list[str]) -> dict:
    cache = {}
    if COORD_CACHE.exists():
        cache = json.loads(COORD_CACHE.read_text(encoding="utf-8"))
    missing = [r for r in rsids if r not in cache]
    if missing:
        _log(f"Resolving {len(missing)} rsID(s) against Ensembl GRCh38")
        for rsid in missing:
            info = _ensembl_variant(rsid)
            cache[rsid] = info if info else {"error": "not resolved"}
            if info:
                _log(f"  {rsid} -> {info['chrom']}:{info['pos']} ({info['allele_string']})")
        COORD_CACHE.parent.mkdir(parents=True, exist_ok=True)
        COORD_CACHE.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    return cache


# --------------------------------------------------------------------------
# Sequence and prediction
# --------------------------------------------------------------------------

def _gene_window(gtf, gene: str):
    """512 kbp window centred on the gene midpoint -- identical to the built dataset."""
    from alphagenome.data import gene_annotation

    interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene)
    interval = interval.resize(WINDOW_SIZE)
    chrom = interval.chromosome
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return chrom, int(interval.start), int(interval.end)          # 0-based half-open


def _reference_sequence(chrom: str, start: int, end: int) -> str:
    region = f"{chrom}:{start + 1}-{end}"
    proc = subprocess.run(
        [SAMTOOLS, "faidx", str(REF_FASTA), region],
        check=True, capture_output=True, text=True,
    )
    seq = "".join(l.strip() for l in proc.stdout.splitlines() if not l.startswith(">")).upper()
    if len(seq) < WINDOW_SIZE:
        seq = seq + "N" * (WINDOW_SIZE - len(seq))
    return seq[:WINDOW_SIZE]


def _alt_allele(allele_string: str | None, ref_base: str, effect_allele: str) -> str | None:
    """Pick the base to substitute in.

    Prefer the effect allele curated in snps.csv when it differs from the
    reference base; otherwise take the first allele in Ensembl's string that is a
    single base and differs from the reference. Returns None when no valid
    substitution exists (e.g. the effect allele *is* the reference)."""
    if effect_allele and len(effect_allele) == 1 and effect_allele != ref_base:
        return effect_allele
    for a in (allele_string or "").split("/"):
        if len(a) == 1 and a != ref_base and a in "ACGT":
            return a
    return None


def _predict(client, seq: str):
    from alphagenome.models import dna_client
    return client.predict_sequence(
        seq,
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
        ontology_terms=ONTOLOGIES,
    )


def _footprint(ref_values, alt_values) -> dict:
    """Change statistics, over the full window and over the classifier's crop.

    `crop_abs_delta` is the quantity directly comparable to the scramble's
    |Delta_in|: total absolute change over the same 32,768 units, same six
    channels."""
    import numpy as np

    ref = np.asarray(ref_values, dtype="float64")
    alt = np.asarray(alt_values, dtype="float64")
    diff = alt - ref

    centre = ref.shape[0] // 2
    lo = centre - CROP // 2
    ref_c, diff_c = ref[lo:lo + CROP, :], diff[lo:lo + CROP, :]

    ref_sum = float(np.abs(ref_c).sum())
    return {
        "window_abs_delta": float(np.abs(diff).sum()),
        "window_max_abs_delta": float(np.abs(diff).max()),
        "crop_abs_delta": float(np.abs(diff_c).sum()),
        "crop_max_abs_delta": float(np.abs(diff_c).max()),
        "crop_ref_signal": ref_sum,
        "crop_rel_delta": float(np.abs(diff_c).sum() / ref_sum) if ref_sum > 0 else None,
        "crop_changed_frac": float((np.abs(diff_c) > 1e-6).any(axis=1).mean()),
        "crop_signed_delta": float(diff_c.sum()),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rsids", help="Comma-separated subset (default: all of snps.csv)")
    ap.add_argument("--panel-only", action="store_true",
                    help="Only variants whose gene is one of the eleven classifier windows")
    ap.add_argument("--out", default=str(OUT_PATH))
    args = ap.parse_args()

    from dotenv import load_dotenv
    load_dotenv(Path.home() / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY not found in the environment or ~/.env.")

    import pandas as pd
    from alphagenome.models import dna_client

    variants = [r for r in csv.DictReader(SNPS_CSV.open(encoding="utf-8"))
                if r.get("rsid", "").startswith("rs")]
    if args.panel_only:
        variants = [v for v in variants if v["gene"] in PANEL_GENES]
    if args.rsids:
        keep = {s.strip() for s in args.rsids.split(",")}
        variants = [v for v in variants if v["rsid"] in keep]
    _log(f"{len(variants)} variant(s) to measure")

    coords = _load_coords([v["rsid"] for v in variants])

    _log(f"Loading GTF cache from {GTF_CACHE}")
    gtf = pd.read_feather(GTF_CACHE)
    client = dna_client.create(api_key)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    records = {}
    if out_path.exists():
        records = json.loads(out_path.read_text(encoding="utf-8")).get("per_variant", {})
        _log(f"Resuming: {len(records)} variant(s) already measured")

    # One reference prediction per gene window, reused by every variant in it.
    ref_cache: dict[str, tuple] = {}

    for idx, v in enumerate(variants, 1):
        rsid, gene = v["rsid"], v["gene"]
        tag = f"({idx}/{len(variants)}) {rsid} [{gene}]"
        if rsid in records and "error" not in records[rsid]:
            _log(f"{tag}: cached, skipping")
            continue

        info = coords.get(rsid) or {}
        if "error" in info or not info:
            _log(f"{tag}: no GRCh38 coordinate, skipping")
            records[rsid] = {"gene": gene, "error": "coordinate not resolved"}
            continue

        try:
            chrom, wstart, wend = _gene_window(gtf, gene)
        except Exception as exc:                                   # noqa: BLE001
            _log(f"{tag}: window lookup FAILED -- {exc}")
            records[rsid] = {"gene": gene, "error": f"window lookup: {exc}"}
            continue

        pos0 = info["pos"] - 1                                     # to 0-based
        if info["chrom"] != chrom or not (wstart <= pos0 < wend):
            _log(f"{tag}: {info['chrom']}:{info['pos']} is OUTSIDE the {gene} window "
                 f"{chrom}:{wstart}-{wend}; recording and skipping")
            records[rsid] = {
                "gene": gene, "chrom": info["chrom"], "pos": info["pos"],
                "error": "variant outside the gene's 512 kbp window",
            }
            continue

        offset = pos0 - wstart
        in_crop = abs(offset - WINDOW_SIZE // 2) < CROP // 2

        if gene not in ref_cache:
            _log(f"{tag}: predicting reference window {chrom}:{wstart}-{wend}")
            seq = _reference_sequence(chrom, wstart, wend)
            ref_cache[gene] = (seq, _predict(client, seq).rna_seq.values)
        ref_seq, ref_vals = ref_cache[gene]

        ref_base = ref_seq[offset]
        alt_base = _alt_allele(info.get("allele_string"), ref_base,
                               (v.get("effect_allele") or "").strip())
        if alt_base is None:
            _log(f"{tag}: no substitutable alternate allele "
                 f"(ref={ref_base}, alleles={info.get('allele_string')}); skipping")
            records[rsid] = {
                "gene": gene, "chrom": chrom, "pos": info["pos"], "ref_base": ref_base,
                "error": "no alternate allele distinct from reference",
            }
            continue

        alt_seq = ref_seq[:offset] + alt_base + ref_seq[offset + 1:]
        _log(f"{tag}: {ref_base}->{alt_base} at offset {offset} "
             f"({'inside' if in_crop else 'OUTSIDE'} the classifier crop)")
        try:
            alt_vals = _predict(client, alt_seq).rna_seq.values
        except Exception as exc:                                   # noqa: BLE001
            _log(f"{tag}: prediction FAILED -- {exc}")
            records[rsid] = {"gene": gene, "error": f"prediction: {exc}"}
            continue

        rec = {
            "gene": gene,
            "window": v.get("window"),
            "chrom": chrom,
            "pos": info["pos"],
            "ref_base": ref_base,
            "alt_base": alt_base,
            "window_offset": offset,
            "inside_classifier_crop": bool(in_crop),
            "consequence": v.get("consequence"),
            "expected_rnaseq_signal": v.get("expected_rnaseq_signal"),
            "control_class": v.get("control_class"),
            "scramble_delta_in": SCRAMBLE_DELTA_IN.get(gene),
        }
        rec.update(_footprint(ref_vals, alt_vals))
        if rec["scramble_delta_in"]:
            rec["fraction_of_scramble"] = rec["crop_abs_delta"] / rec["scramble_delta_in"]
        records[rsid] = rec
        _log(f"{tag}: crop |delta| = {rec['crop_abs_delta']:.1f} "
             f"({rec.get('fraction_of_scramble', float('nan')):.4f} of the scramble)")

        out_path.write_text(json.dumps({
            "generated": datetime.now(timezone.utc).isoformat(),
            "window_size": WINDOW_SIZE, "crop": CROP, "ontologies": ONTOLOGIES,
            "note": "Single-base substitution against the reference window. "
                    "crop_abs_delta is directly comparable to the promoter "
                    "scramble's |Delta_in| (same crop, channels and strands).",
            "per_variant": records,
        }, indent=2), encoding="utf-8")

    _summarise(records)
    _log(f"Wrote {out_path}")
    return 0


def _summarise(records: dict) -> None:
    ok = {k: v for k, v in records.items() if "error" not in v}
    if not ok:
        _log("No successful measurements to summarise.")
        return

    _log("")
    _log("=== footprint against the pre-registered expectation ===")
    header = f"{'rsid':<12} {'gene':<8} {'expected':<12} {'cls':<5} {'crop|d|':>10} {'/scramble':>10} {'crop':>6}"
    _log(header)
    for rsid, r in sorted(ok.items(), key=lambda kv: -kv[1]["crop_abs_delta"]):
        frac = r.get("fraction_of_scramble")
        _log(f"{rsid:<12} {r['gene']:<8} {str(r.get('expected_rnaseq_signal')):<12} "
             f"{str(r.get('control_class')):<5} {r['crop_abs_delta']:>10.1f} "
             f"{(f'{frac:.4f}' if frac is not None else '--'):>10} "
             f"{('in' if r['inside_classifier_crop'] else 'OUT'):>6}")

    import statistics
    from scipy.stats import mannwhitneyu

    def _group(label, inside_only=False):
        return [r["crop_abs_delta"] for r in ok.values()
                if r.get("expected_rnaseq_signal") == label
                and (r["inside_classifier_crop"] or not inside_only)]

    _log("")
    _log("=== the pre-registered test ===")
    _log("`snps.csv` labels each variant VISIBLE (regulatory; the frozen model "
         "should see it) or NULL_CODING (coding; it should not) from the "
         "literature, before any measurement. A separation licenses reading the "
         "frozen model's output as regulatory; an overlap does not.")

    for inside_only in (False, True):
        scope = "inside the classifier crop only" if inside_only else "all variants"
        vis, nul = _group("VISIBLE", inside_only), _group("NULL_CODING", inside_only)
        if len(vis) < 3 or len(nul) < 3:
            continue
        u, p = mannwhitneyu(vis, nul, alternative="greater")
        _log(f"[{scope}] VISIBLE n={len(vis)} median={statistics.median(vis):.1f}; "
             f"NULL_CODING n={len(nul)} median={statistics.median(nul):.1f}; "
             f"Mann-Whitney one-sided p={p:.3f} -> separation "
             f"{'PRESENT' if p < 0.05 else 'ABSENT'}")

    outside = [r for r in ok.values() if not r["inside_classifier_crop"]]
    pos = [r for r in ok.values() if r.get("control_class") == "POS"]
    pos_out = [r for r in pos if not r["inside_classifier_crop"]]
    _log("")
    _log(f"crop: {len(outside)}/{len(ok)} measured variants fall OUTSIDE the "
         f"32,768-unit window the classifier reads, INCLUDING {len(pos_out)} of "
         f"the {len(pos)} positive-control regulatory variants. That is the "
         f"asymmetry that matters: the crop removes the positive controls "
         f"specifically, because distal regulatory variants are distal.")


if __name__ == "__main__":
    raise SystemExit(main())
