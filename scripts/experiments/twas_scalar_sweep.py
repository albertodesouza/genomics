#!/usr/bin/env python3
"""Run the "poor man's PrediXcan" TWAS arm on all 42 published genes.

Reviewer ask (Q5): a TWAS-style baseline alongside the GWAS gene scores. The scalar
classifier already exists (genomics.predictors.genotype_based.analysis.
train_gene_expression_classifier) and was previously only exercised ad hoc in
notebooks/pigmentation_alphagenome_cnn_variant_scoring.ipynb Section 9. This script runs it
on every one of the 42 published per-gene arms (9 panel + 33 controls, same genes as
poolmax_final_table.csv), resolving each gene's config the same way
poolmax_knockdown_sweep.find_config does, so the scalar features come from the exact same
aligned tensor cache the CNN arm reads -- zero AlphaGenome calls, zero re-alignment.

Costs almost nothing: a LogisticRegression over a handful of scalar features per
individual, not a CNN. The slow part is materialising the scalar features (one
dataset[i] per individual, ~1072 individuals per gene), which is CPU/IO-bound, not GPU-bound
-- hence the higher default parallelism than the seed-robustness sweep.

Usage:
  python3 scripts/experiments/twas_scalar_sweep.py --dry-run
  python3 scripts/experiments/twas_scalar_sweep.py --jobs 6
  python3 scripts/experiments/twas_scalar_sweep.py --genes OCA2,TYR
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

REPO = Path("/home/breno/I2CA/genomics")
PY = "/home/breno/miniforge3/envs/genomics/bin/python3"
os.environ["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{os.environ.get('PATH', '')}"

LOG_DIR = REPO / "results/genotype_based_predictor/logs/twas_sweep"
OUT_TABLE = REPO / "results/genotype_based_predictor/twas_final_table.csv"

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
DRAW1 = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
         "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]
DRAW3 = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6",
         "KIAA0319", "RIDA", "SEM1", "SFMBT2", "SUMF2"]
CONTROL = DRAW1 + DRAW2 + DRAW3
GENES = PANEL + CONTROL

OUT_SUBDIR = "scalar_gene_expression_logreg_mean_separate"


def _log(m: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {m}", flush=True)


def resolve_config(gene: str) -> Path:
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {gene} resolves to {len(hits)} published configs, expected 1")
    return hits[0]


def out_dir_for(gene: str, cfg_path: Path) -> Path:
    import yaml
    cfg = yaml.safe_load(cfg_path.read_text())
    return REPO / cfg["dataset_input"]["results_dir"] / OUT_SUBDIR


def bal_acc(results_json: Path) -> dict:
    d = json.loads(results_json.read_text())
    pcm = d["per_class_metrics"]
    rec_strong = float(pcm["strong pigmentation"]["recall"])
    rec_weak = float(pcm["weak pigmentation"]["recall"])
    out = {"bal_acc": (rec_strong + rec_weak) / 2.0, "acc": float(d["weighted_accuracy"])}
    if "auc_vs_strong_pigmentation" in d:
        out["auc_vs_strong"] = float(d["auc_vs_strong_pigmentation"]["auc"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", default=None, help="comma-separated; default: all 42")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-done", action="store_true", default=True)
    args = ap.parse_args()

    genes = [g.strip().upper() for g in args.genes.split(",")] if args.genes else list(GENES)
    cfgs = {g: resolve_config(g) for g in genes}
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{env.get('PATH', '')}"

    todo = list(genes)
    running: list[tuple[str, subprocess.Popen, object]] = []
    done: dict[str, int] = {}
    t0 = time.monotonic()

    while todo or running:
        while todo and len(running) < args.jobs:
            g = todo.pop(0)
            od = out_dir_for(g, cfgs[g])
            if args.skip_done and (od / "test_results.json").exists():
                _log(f"{g:<9} SKIP: {od / 'test_results.json'} already present")
                done[g] = 0
                continue
            if args.dry_run:
                _log(f"{g:<9} DRY RUN: would run {cfgs[g].name} -> {od}")
                done[g] = 0
                continue
            log_path = LOG_DIR / f"{g.lower()}.log"
            fh = open(log_path, "w")
            cmd = [PY, "-m", "genomics.predictors.genotype_based.analysis.train_gene_expression_classifier",
                   str(cfgs[g])]
            p = subprocess.Popen(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, env=env)
            running.append((g, p, fh))
            _log(f"{g:<9} start (pid {p.pid}) [{len(running)}/{args.jobs} slots]")
        time.sleep(3)
        for item in list(running):
            g, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                done[g] = p.returncode
                mins = (time.monotonic() - t0) / 60
                _log(f"{g:<9} done rc={p.returncode} ({len(done)}/{len(genes)} total, {mins:.1f} min elapsed)")

    bad = {g: rc for g, rc in done.items() if rc != 0}
    if bad:
        _log(f"FAILED genes (rc != 0): {bad}")

    if args.dry_run:
        return 0

    import pandas as pd
    rows = []
    for g in genes:
        od = out_dir_for(g, cfgs[g])
        test_j, val_j = od / "test_results.json", od / "val_results.json"
        if not (test_j.exists() and val_j.exists()):
            _log(f"{g:<9} WARN: missing results json in {od}")
            continue
        t, v = bal_acc(test_j), bal_acc(val_j)
        rows.append({"gene": g, "bal_acc": t["bal_acc"], "bal_acc_val": v["bal_acc"],
                     "acc": t["acc"], "acc_val": v["acc"],
                     "auc_vs_strong": t.get("auc_vs_strong"),
                     "auc_vs_strong_val": v.get("auc_vs_strong")})

    out_df = pd.DataFrame(rows).sort_values("gene")
    out_df.to_csv(OUT_TABLE, index=False)
    _log(f"wrote {len(out_df)} rows to {OUT_TABLE}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
