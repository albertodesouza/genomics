#!/usr/bin/env python3
"""Retrain every published per-gene arm reading a non-melanocyte RNA-seq track.

Reviewer ask (Q2a): does performance survive if the classifier reads an unrelated
cell-type track instead of melanocyte? The only ontologies already predicted (and
therefore free -- zero new AlphaGenome calls) for this cohort/window set are the three
cached in every individual's rna_seq_metadata.json:

  CL:1000458  melanocyte of skin              (the published track)
  CL:0000346  hair follicle dermal papilla cell
  CL:2000092  hair follicular keratinocyte

Both alternates are skin/hair lineage, not a genuinely unrelated tissue (liver, blood,
...) -- that would require new AlphaGenome inference across all 42 genes x ~1072
individuals x 2 haplotypes, tracked separately. This sweep runs the free version: same
42 genes, same strand-per-gene restriction, only `ontology_terms` swapped.

DITA alignment itself does not depend on ontology (it is built from genotype/indel
structure), and alignment_cache already holds entries for all 42 published genes -- so
this reuses the expensive chain-building step and only re-embeds a different track's
values into it, via a fresh, isolated processed-tensor cache and results_dir (both
namespaced so nothing published is touched).

Usage:
  python3 scripts/experiments/wrong_celltype_sweep.py --write-configs
  python3 scripts/experiments/wrong_celltype_sweep.py --dry-run
  python3 scripts/experiments/wrong_celltype_sweep.py --jobs 3
  python3 scripts/experiments/wrong_celltype_sweep.py --ontology CL:2000092 --genes OCA2
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import yaml

REPO = Path("/home/breno/I2CA/genomics")
PY = "/home/breno/miniforge3/envs/genomics/bin/python3"
os.environ["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{os.environ.get('PATH', '')}"

ONTOLOGY_LABEL = {
    "CL:0000346": "hair follicle dermal papilla cell",
    "CL:2000092": "hair follicular keratinocyte",
}
DEFAULT_ONTOLOGY = "CL:0000346"

CONFIG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/wrong_celltype"
RUNS_ROOT = REPO / "results/genotype_based_predictor/runs_wrongcelltype"
LOG_DIR = REPO / "results/genotype_based_predictor/logs/wrongcelltype_sweep"

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
DRAW1 = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
         "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]
DRAW3 = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6",
         "KIAA0319", "RIDA", "SEM1", "SFMBT2", "SUMF2"]
CONTROL = DRAW1 + DRAW2 + DRAW3
GENES = PANEL + CONTROL

ALLOWED_DRIFT = {
    "dataset_input.ontology_terms", "dataset_input.results_dir",
    "dataset_input.processed_cache_dir", "wandb.run_name",
}


def _log(m: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {m}", flush=True)


def resolve_base_config(gene: str) -> Path:
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {gene} resolves to {len(hits)} published configs, expected 1")
    return hits[0]


def write_variant(gene: str, base_cfg_path: Path, ontology: str) -> Path:
    tag = ontology.replace(":", "").lower()
    cfg = yaml.safe_load(base_cfg_path.read_text())
    cfg["dataset_input"]["ontology_terms"] = [ontology]
    cfg["dataset_input"]["results_dir"] = f"results/genotype_based_predictor/runs_wrongcelltype/{gene.lower()}_{tag}"
    cfg["dataset_input"]["processed_cache_dir"] = f"results/cache/genotype_based_predictor/wrongcelltype_{tag}_{gene.lower()}"
    cfg.setdefault("wandb", {})["run_name"] = f"wrongcelltype-{tag}-{gene.lower()}"

    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    out = CONFIG_DIR / f"{gene.lower()}_{tag}.yaml"
    out.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return out


def verify_variant(base_cfg_path: Path, variant_path: Path, gene: str, ontology: str) -> None:
    base = yaml.safe_load(base_cfg_path.read_text())
    cfg = yaml.safe_load(variant_path.read_text())
    diffs = []

    def walk(a, b, pre=""):
        for k in sorted(set(a) | set(b)):
            p = f"{pre}.{k}" if pre else k
            va, vb = a.get(k, "<missing>"), b.get(k, "<missing>")
            if isinstance(va, dict) and isinstance(vb, dict):
                walk(va, vb, p)
            elif va != vb:
                diffs.append(p)

    walk(base, cfg)
    extra = set(diffs) - ALLOWED_DRIFT
    if extra:
        raise SystemExit(f"ABORT: {variant_path} drifts from its base in {sorted(extra)}")
    if cfg["dataset_input"]["ontology_terms"] != [ontology]:
        raise SystemExit(f"ABORT: {variant_path} ontology_terms is not [{ontology}]")
    if cfg["dataset_input"]["processed_cache_dir"] == base["dataset_input"]["processed_cache_dir"]:
        raise SystemExit(f"ABORT: {variant_path} shares the published tensor cache")


def find_run_dir(gene: str, tag: str) -> Path | None:
    hits = sorted(RUNS_ROOT.glob(f"{gene.lower()}_{tag}/*/test_best_accuracy_results.json"))
    return hits[0].parent if len(hits) == 1 else None


def bal_acc_from(run_dir: Path, split: str) -> float:
    d = json.loads((run_dir / f"{split}_best_accuracy_results.json").read_text())
    pcm = d["per_class_metrics"]
    return (float(pcm["strong pigmentation"]["recall"]) + float(pcm["weak pigmentation"]["recall"])) / 2.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", default=None, help="comma-separated; default: all 42")
    ap.add_argument("--ontology", default=DEFAULT_ONTOLOGY, choices=list(ONTOLOGY_LABEL))
    ap.add_argument("--jobs", type=int, default=3, help="parallel arms; full CNN training, GPU-bound (empirically caps at 3 on this box)")
    ap.add_argument("--write-configs", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-done", action="store_true", default=True)
    args = ap.parse_args()

    tag = args.ontology.replace(":", "").lower()
    genes = [g.strip().upper() for g in args.genes.split(",")] if args.genes else list(GENES)
    out_table = REPO / f"results/genotype_based_predictor/wrongcelltype_{tag}_final_table.csv"

    jobs_list: list[tuple[str, Path]] = []
    for g in genes:
        base_cfg = resolve_base_config(g)
        variant = write_variant(g, base_cfg, args.ontology)
        verify_variant(base_cfg, variant, g, args.ontology)
        jobs_list.append((g, variant))
    _log(f"{len(jobs_list)} arms verified against their published base config "
         f"(ontology={args.ontology}, {ONTOLOGY_LABEL[args.ontology]})")
    if args.write_configs:
        return 0

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{env.get('PATH', '')}"

    todo = list(jobs_list)
    running: list[tuple[str, subprocess.Popen, object]] = []
    done: dict[str, int] = {}
    t0 = time.monotonic()

    def launch(gene: str, cfg: Path):
        gtag = f"{gene.lower()}_{tag}"
        log_path = LOG_DIR / f"{gtag}.log"
        script = (f"set -e\n{PY} -m genomics genotype train {cfg}\n"
                  f"{PY} -m genomics genotype test {cfg}\n")
        fh = open(log_path, "w")
        p = subprocess.Popen(["bash", "-c", script], cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, env=env)
        return p, fh

    while todo or running:
        while todo and len(running) < args.jobs:
            g, cfg = todo.pop(0)
            if args.skip_done and find_run_dir(g, tag) is not None:
                _log(f"{g:<9} SKIP: already has test results")
                done[g] = 0
                continue
            if args.dry_run:
                _log(f"{g:<9} DRY RUN: would train+test {cfg.name}")
                done[g] = 0
                continue
            p, fh = launch(g, cfg)
            running.append((g, p, fh))
            _log(f"{g:<9} start (pid {p.pid}) [{len(running)}/{args.jobs} slots]")
        time.sleep(5)
        for item in list(running):
            g, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                done[g] = p.returncode
                mins = (time.monotonic() - t0) / 60
                _log(f"{g:<9} done rc={p.returncode} ({len(done)}/{len(jobs_list)} total, {mins:.1f} min elapsed)")

    bad = {g: rc for g, rc in done.items() if rc != 0}
    if bad:
        _log(f"FAILED arms (rc != 0): {bad}")
    if args.dry_run:
        return 0

    import pandas as pd
    rows = []
    for g, _cfg in jobs_list:
        run_dir = find_run_dir(g, tag)
        if run_dir is None:
            continue
        rows.append({"gene": g, "bal_acc": bal_acc_from(run_dir, "test"),
                     "bal_acc_val": bal_acc_from(run_dir, "val")})
    out_df = pd.DataFrame(rows).sort_values("gene")
    out_df.to_csv(out_table, index=False)
    _log(f"wrote {len(out_df)} rows to {out_table}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
