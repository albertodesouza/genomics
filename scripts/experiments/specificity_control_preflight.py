#!/usr/bin/env python3
"""Pre-flight magnitude screen for the specificity-control gene panel.

The specificity control asks whether a gene with *no* phenotype role but a
*comparable delivered perturbation* fails to respond to the promoter scramble.
The second half of that sentence is the part that is easy to get wrong: EDAR and
TCHH in the published panel deliver only 59 and 53 units of absolute change to
the classifier input, against TYR's 54,297, so their flat responses carry no
information about reliance. A control gene picked purely for biological
irrelevance can land in the same regime and produce an uninterpretable result
after a full cohort run.

So before spending ~6,400 AlphaGenome calls, this script runs *one* reference
call per candidate gene and measures the quantity that predicts delivered
magnitude: the total predicted RNA-seq signal inside the 32,768-unit crop the
classifier actually reads, in the three retained ontologies, on both strands.

It reports the same statistic for the eleven panel genes, whose realised
|delta_in| is already known from results/.../readout_normalisation.json, so the
candidates can be ranked against a calibrated scale rather than in the abstract.

Usage:
  python3 scripts/experiments/specificity_control_preflight.py
  ... --candidates PSMC4,EIF1B,TPM2 --skip-panel          # narrow re-run
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
os.chdir(REPO_ROOT)

REF_FASTA = Path("/dados/GENOMICS_DATA/top3/refs/GRCh38_full_analysis_set_plus_decoy_hla.fa")
SAMTOOLS = os.environ.get("SAMTOOLS_BIN", "/home/breno/miniforge3/envs/genomics/bin/samtools")
GTF_CACHE = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all/gtf_cache.feather")
WINDOW_SIZE = 524288
CROP = 32768
ONTOLOGIES = ["CL:1000458", "CL:0000346", "CL:2000092"]

# The eleven windows of the published panel, for calibration.
PANEL_GENES = [
    "SLC24A5", "TYR", "SLC45A2", "TYRP1", "DDB1", "MFSD12",
    "MC1R", "HERC2", "OCA2", "EDAR", "TCHH",
]

# Candidates carry no Gene Ontology annotation under pigmentation, melanin
# biosynthesis or melanocyte differentiation; this is the same vetted list used
# for the dosage-based irrelevant-gene control in the identifiability audit.
DEFAULT_CANDIDATES = [
    "ECHDC3", "EIF1B", "FRA10AC1", "LACTB2", "LRRC36",
    "PPP1R3E", "PRSS55", "PSMC4", "SMCR8", "SPRED2", "TPM2",
]

# Realised delivered perturbation for the panel, from the published run.
DELTA_IN = {
    "TYR": 54296.6, "DDB1": 34150.3, "TYRP1": 26589.5, "SLC45A2": 23730.6,
    "MFSD12": 11868.3, "SLC24A5": 6537.0, "OCA2": 2109.0, "MC1R": 1391.0,
    "HERC2": 798.6, "EDAR": 59.0, "TCHH": 53.0,
}

OUT_PATH = REPO_ROOT / "results" / "genotype_based_predictor" / "specificity_control_preflight.json"


def _log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc).isoformat()}] {msg}", flush=True)


def _gene_window(gtf, gene: str):
    """512 kbp window centred on the gene midpoint, matching the built dataset."""
    from alphagenome.data import gene_annotation

    interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene)
    interval = interval.resize(WINDOW_SIZE)
    chrom = interval.chromosome
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return chrom, int(interval.start), int(interval.end)


def _reference_sequence(chrom: str, start: int, end: int) -> str:
    """Reference window via samtools faidx, 0-based half-open -> 1-based inclusive."""
    region = f"{chrom}:{start + 1}-{end}"
    proc = subprocess.run(
        [SAMTOOLS, "faidx", str(REF_FASTA), region],
        check=True, capture_output=True, text=True,
    )
    seq = "".join(line.strip() for line in proc.stdout.splitlines() if not line.startswith(">"))
    seq = seq.upper()
    if len(seq) < WINDOW_SIZE:
        seq = seq + "N" * (WINDOW_SIZE - len(seq))
    return seq[:WINDOW_SIZE]


def _crop_signal(values) -> dict:
    """Signal statistics over the central CROP units the classifier reads."""
    import numpy as np

    arr = np.asarray(values, dtype="float64")            # (524288, 6)
    centre = arr.shape[0] // 2
    lo = centre - CROP // 2
    crop = arr[lo:lo + CROP, :]
    # The classifier normalises by log1p over a per-track training maximum; the
    # ordering of loci is what matters here, so report the raw integral and the
    # log1p integral side by side.
    return {
        "crop_total_signal": float(np.abs(crop).sum()),
        "crop_total_log1p": float(np.log1p(np.abs(crop)).sum()),
        "crop_max": float(np.abs(crop).max()),
        "crop_nonzero_frac": float((np.abs(crop) > 1e-6).mean()),
        "window_total_signal": float(np.abs(arr).sum()),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", help="Comma-separated candidate gene symbols (default: vetted irrelevant panel)")
    ap.add_argument("--skip-panel", action="store_true", help="Do not re-measure the eleven published panel genes")
    ap.add_argument("--out", default=str(OUT_PATH))
    args = ap.parse_args()

    from dotenv import load_dotenv
    load_dotenv(Path.home() / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY not found in the environment or ~/.env.")

    import pandas as pd
    from alphagenome.models import dna_client

    _log(f"Loading GTF cache from {GTF_CACHE}")
    gtf = pd.read_feather(GTF_CACHE)
    client = dna_client.create(api_key)

    candidates = [g.strip() for g in args.candidates.split(",")] if args.candidates else list(DEFAULT_CANDIDATES)
    genes = ([] if args.skip_panel else list(PANEL_GENES)) + candidates

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    records = {}
    if out_path.exists():
        records = json.loads(out_path.read_text(encoding="utf-8")).get("per_gene", {})
        _log(f"Resuming: {len(records)} gene(s) already measured")

    for idx, gene in enumerate(genes, 1):
        if gene in records:
            _log(f"({idx}/{len(genes)}) {gene}: cached, skipping")
            continue
        try:
            chrom, start, end = _gene_window(gtf, gene)
        except Exception as exc:                                   # noqa: BLE001
            _log(f"({idx}/{len(genes)}) {gene}: window lookup FAILED -- {exc}")
            records[gene] = {"error": f"window lookup: {exc}"}
            continue

        _log(f"({idx}/{len(genes)}) {gene}: {chrom}:{start}-{end}, predicting reference window")
        try:
            seq = _reference_sequence(chrom, start, end)
            outputs = client.predict_sequence(
                seq,
                requested_outputs=[dna_client.OutputType.RNA_SEQ],
                ontology_terms=ONTOLOGIES,
            )
            stats = _crop_signal(outputs.rna_seq.values)
        except Exception as exc:                                   # noqa: BLE001
            _log(f"({idx}/{len(genes)}) {gene}: prediction FAILED -- {exc}")
            records[gene] = {"error": f"prediction: {exc}"}
            continue

        stats.update({
            "chromosome": chrom,
            "start": start,
            "end": end,
            "is_panel_gene": gene in PANEL_GENES,
            "published_delta_in": DELTA_IN.get(gene),
        })
        records[gene] = stats
        _log(f"    crop_total_signal={stats['crop_total_signal']:.1f} "
             f"log1p={stats['crop_total_log1p']:.1f} nonzero={stats['crop_nonzero_frac']:.3f}")

        payload = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "window_size": WINDOW_SIZE,
            "crop": CROP,
            "ontologies": ONTOLOGIES,
            "note": "Reference-only AlphaGenome call per gene. crop_total_signal is the "
                    "magnitude proxy used to match control genes to the responding panel "
                    "genes on delivered perturbation.",
            "per_gene": records,
        }
        out_path.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    _log(f"Wrote {out_path}")

    # Ranked summary, panel genes annotated with their realised delivered perturbation.
    ok = {g: r for g, r in records.items() if "error" not in r}
    ranked = sorted(ok.items(), key=lambda kv: -kv[1]["crop_total_signal"])
    print("\n{:<10} {:>16} {:>14} {:>12}  {}".format("gene", "crop_signal", "crop_log1p", "delta_in", "role"))
    for gene, rec in ranked:
        d = rec.get("published_delta_in")
        print("{:<10} {:>16.1f} {:>14.1f} {:>12}  {}".format(
            gene, rec["crop_total_signal"], rec["crop_total_log1p"],
            f"{d:.0f}" if d else "-", "panel" if rec["is_panel_gene"] else "CANDIDATE",
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
