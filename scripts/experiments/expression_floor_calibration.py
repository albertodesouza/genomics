#!/usr/bin/env python3
"""Calibrate the expression floor against random genomic windows.

WHY A FLOOR. A flat knockdown response is interpretable only where the perturbation
delivered to the classifier's input is comparable across genes; the paper withdraws
EDAR and TCHH on exactly this ground, but withdraws them *after* seeing they responded
little. A floor declared in advance turns that into an inclusion criterion, and it is
also what the pipeline should apply before spending 2,144 AlphaGenome calls building a
gene's windows for the whole cohort: LRRC36 (29) and PRSS55 (23) would never have been
built.

WHY CALIBRATE IT THIS WAY. The floor must not be set by looking at the panel-vs-control
outcome it will later be used to report -- that is threshold-shopping, and the
sensitivity analysis in single_gene_floor_sensitivity.py shows the p-value in this study
is non-monotonic in the threshold, so a shopped threshold would find whatever it wanted.
The anchor used here is external to the labels and to the panel: the distribution of
crop_total_signal over random genomic windows, in the same three ontologies. A gene
below that distribution is one the frozen model does not distinguish from an arbitrary
piece of genome, and no comparison involving it can be read.

This is the same device as the paper's random-window null, applied one level up: there,
a scramble at a random position; here, a whole window at a random position.

Cost: one AlphaGenome call per draw, no cohort involved.

Usage:
  python3 scripts/experiments/expression_floor_calibration.py --draws 40
  ... --percentile 95        # the floor is this percentile of the random-window draws
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
WINDOW_SIZE = 524288
CROP = 32768
ONTOLOGIES = ["CL:1000458", "CL:0000346", "CL:2000092"]
PREFLIGHT = REPO_ROOT / "results/genotype_based_predictor/specificity_control_preflight.json"
OUT = REPO_ROOT / "results/genotype_based_predictor/expression_floor_calibration.json"
SEED = 13


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def chrom_lengths():
    """Autosome lengths from the reference index; the panel spans autosomes only."""
    fai = REF_FASTA.with_suffix(REF_FASTA.suffix + ".fai")
    if not fai.exists():
        raise SystemExit(f"ABORT: {fai} missing; samtools faidx the reference first.")
    out = {}
    for line in fai.read_text().splitlines():
        f = line.split("\t")
        if f[0] in {f"chr{i}" for i in range(1, 23)}:
            out[f[0]] = int(f[1])
    return out


def ref_seq(chrom, start, end):
    p = subprocess.run([SAMTOOLS, "faidx", str(REF_FASTA), f"{chrom}:{start + 1}-{end}"],
                       check=True, capture_output=True, text=True)
    s = "".join(l.strip() for l in p.stdout.splitlines() if not l.startswith(">")).upper()
    return (s + "N" * WINDOW_SIZE)[:WINDOW_SIZE]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=40)
    ap.add_argument("--percentile", type=float, default=95.0)
    ap.add_argument("--max-n-frac", type=float, default=0.10,
                    help="reject a draw whose reference window is more than this fraction N")
    args = ap.parse_args()

    import numpy as np
    from dotenv import load_dotenv
    from alphagenome.models import dna_client

    load_dotenv(Path.home() / ".env")
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        raise SystemExit("ALPHAGENOME_API_KEY not found.")
    client = dna_client.create(key)

    lens = chrom_lengths()
    names = sorted(lens)
    weights = np.array([lens[c] for c in names], dtype=float)
    weights /= weights.sum()
    rng = np.random.default_rng(SEED)          # fixed: the floor must not move on a rerun

    prior = {}
    if OUT.exists():
        prior = json.loads(OUT.read_text()).get("draws", {})
        _log(f"resuming: {len(prior)} draws present")

    draws = dict(prior)
    attempts = 0
    while len(draws) < args.draws and attempts < args.draws * 6:
        attempts += 1
        c = names[int(rng.choice(len(names), p=weights))]
        start = int(rng.integers(0, max(1, lens[c] - WINDOW_SIZE)))
        key_ = f"{c}:{start}"
        if key_ in draws:
            continue
        try:
            seq = ref_seq(c, start, start + WINDOW_SIZE)
        except Exception as e:                                       # noqa: BLE001
            _log(f"  {key_}: faidx failed ({type(e).__name__}); skipped")
            continue
        n_frac = seq.count("N") / len(seq)
        if n_frac > args.max_n_frac:
            continue                    # centromeres and gaps are not "a piece of genome"
        try:
            o = client.predict_sequence(seq, requested_outputs=[dna_client.OutputType.RNA_SEQ],
                                        ontology_terms=ONTOLOGIES)
        except Exception as e:                                       # noqa: BLE001
            _log(f"  {key_}: prediction failed ({type(e).__name__}); skipped")
            continue
        a = np.abs(np.asarray(o.rna_seq.values, dtype="float64"))
        mid = a.shape[0] // 2
        crop = a[mid - CROP // 2: mid + CROP // 2, :]
        draws[key_] = {"crop_total_signal": float(crop.sum()), "n_frac": n_frac}
        if len(draws) % 5 == 0:
            _log(f"  {len(draws)}/{args.draws} draws")
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps({"generated_at": datetime.now(timezone.utc).isoformat(),
                                   "window_size": WINDOW_SIZE, "crop": CROP,
                                   "ontologies": ONTOLOGIES, "seed": SEED,
                                   "draws": draws}, indent=1))

    v = np.array([d["crop_total_signal"] for d in draws.values()])
    if not len(v):
        raise SystemExit("ABORT: no usable draws.")
    floor = float(np.percentile(v, args.percentile))
    _log(f"random windows (n={len(v)}): median={np.median(v):.1f} "
         f"p75={np.percentile(v,75):.1f} p90={np.percentile(v,90):.1f} "
         f"p95={np.percentile(v,95):.1f} max={v.max():.1f}")
    _log(f"FLOOR (p{args.percentile:g}) = {floor:.1f}")

    payload = json.loads(OUT.read_text())
    payload["floor"] = {"percentile": args.percentile, "value": floor, "n_draws": int(len(v))}

    if PREFLIGHT.exists():
        pre = json.loads(PREFLIGHT.read_text())["per_gene"]
        rows = [(g, r.get("crop_total_signal"), r.get("is_panel_gene"))
                for g, r in pre.items() if "crop_total_signal" in r]
        rows.sort(key=lambda t: -t[1])
        print(f"\n{'gene':10}{'class':9}{'crop_total_signal':>19}{'admissible':>12}")
        print("-" * 51)
        verdict = {}
        for g, s, isp in rows:
            ok = s >= floor
            verdict[g] = {"crop_total_signal": s, "admissible": bool(ok)}
            print(f"{g:10}{'panel' if isp else 'control':9}{s:19.1f}{'yes' if ok else 'NO':>12}")
        payload["gene_verdict"] = verdict
        ex = [g for g, v_ in verdict.items() if not v_["admissible"]]
        print(f"\nexcluded by the floor: {', '.join(ex) if ex else '(none)'}")
    OUT.write_text(json.dumps(payload, indent=1))
    _log(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
