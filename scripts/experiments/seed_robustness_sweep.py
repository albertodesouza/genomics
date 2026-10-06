#!/usr/bin/env python3
"""Retrain every published per-gene arm under extra training seeds.

Answers reviewer Q3 (single seed per gene, no seed sweep, ranking stability unknown):
how much does balanced accuracy and the panel/control ranking move when only the CNN's
weight init and minibatch order change, everything else held fixed?

ONLY `training.random_seed` VARIES. `data_split.random_seed` and
`dataset_input.normalization_fit_random_seed` stay at 13, so the train/val/test
partition and the aligned-tensor cache are byte-identical to the published run --
`processed_cache_dir` is left untouched and the DITA/AlphaGenome cache is reused, not
rebuilt. This sweep costs zero AlphaGenome calls and zero re-alignment; it only pays for
CNN training+test.

Each of the 42 published arms (config resolved the same way
poolmax_knockdown_sweep.find_config does: the config.yaml sitting next to that gene's
checkpoint under runs_poolmax*/) is retrained at each seed in SEEDS, into an isolated
results_dir under runs_seedsweep/ so no published checkpoint is touched.

Usage:
  python3 scripts/experiments/seed_robustness_sweep.py --write-configs
  python3 scripts/experiments/seed_robustness_sweep.py --dry-run
  python3 scripts/experiments/seed_robustness_sweep.py --jobs 3
  python3 scripts/experiments/seed_robustness_sweep.py --genes OCA2,TYR --seeds 7
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

CONFIG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/seed_robustness"
RUNS_ROOT = REPO / "results/genotype_based_predictor/runs_seedsweep"
LOG_DIR = REPO / "results/genotype_based_predictor/logs/seedsweep"
OUT_TABLE = REPO / "results/genotype_based_predictor/seed_robustness_table.csv"
PUBLISHED_TABLE = REPO / "results/genotype_based_predictor/poolmax_final_table.csv"

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
DRAW1 = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
         "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]
DRAW3 = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6",
         "KIAA0319", "RIDA", "SEM1", "SFMBT2", "SUMF2"]
CONTROL = DRAW1 + DRAW2 + DRAW3
GENES = PANEL + CONTROL

SEEDS = [7, 21, 99]
PUBLISHED_SEED = 13

ALLOWED_DRIFT = {
    "training.random_seed",
    "dataset_input.results_dir",
    "wandb.run_name",
}


def _log(m: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {m}", flush=True)


def resolve_base_config(gene: str) -> Path:
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {gene} resolves to {len(hits)} published configs, expected 1")
    return hits[0]


def write_variant(gene: str, base_cfg_path: Path, seed: int) -> Path:
    cfg = yaml.safe_load(base_cfg_path.read_text())
    cfg.setdefault("training", {})["random_seed"] = seed
    cfg["dataset_input"]["results_dir"] = f"results/genotype_based_predictor/runs_seedsweep/{gene.lower()}_seed{seed}"
    cfg.setdefault("wandb", {})["run_name"] = f"seedsweep-{gene.lower()}-seed{seed}"

    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    out = CONFIG_DIR / f"{gene.lower()}_seed{seed}.yaml"
    out.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return out


def verify_variant(base_cfg_path: Path, variant_path: Path, gene: str, seed: int) -> None:
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
    if cfg["training"]["random_seed"] != seed:
        raise SystemExit(f"ABORT: {variant_path} training.random_seed != {seed}")
    if cfg["data_split"]["random_seed"] != PUBLISHED_SEED:
        raise SystemExit(f"ABORT: {variant_path} data_split.random_seed changed; cache would not be reused")
    if cfg["dataset_input"]["normalization_fit_random_seed"] != PUBLISHED_SEED:
        raise SystemExit(f"ABORT: {variant_path} normalization_fit_random_seed changed; cache would not be reused")
    if cfg["dataset_input"]["processed_cache_dir"] != base["dataset_input"]["processed_cache_dir"]:
        raise SystemExit(f"ABORT: {variant_path} processed_cache_dir changed; cache would not be reused")
    if f"{gene.lower()}_seed{seed}" not in cfg["dataset_input"]["results_dir"]:
        raise SystemExit(f"ABORT: {variant_path} has no isolated results_dir")


def find_run_dir(gene: str, seed: int) -> Path | None:
    hits = sorted(RUNS_ROOT.glob(f"{gene.lower()}_seed{seed}/*/test_best_accuracy_results.json"))
    return hits[0].parent if len(hits) == 1 else None


def bal_acc_from(run_dir: Path, split: str) -> dict:
    d = json.loads((run_dir / f"{split}_best_accuracy_results.json").read_text())
    pcm = d["per_class_metrics"]
    rec_strong = float(pcm["strong pigmentation"]["recall"])
    rec_weak = float(pcm["weak pigmentation"]["recall"])
    return {"bal_acc": (rec_strong + rec_weak) / 2.0}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", default=None, help="comma-separated; default: all 42")
    ap.add_argument("--seeds", default=None, help="comma-separated ints; default: 7,21,99")
    ap.add_argument("--jobs", type=int, default=3, help="parallel arms (memory: max ~3 workers/arm)")
    ap.add_argument("--write-configs", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-done", action="store_true", default=True)
    args = ap.parse_args()

    genes = [g.strip().upper() for g in args.genes.split(",")] if args.genes else list(GENES)
    seeds = [int(s) for s in args.seeds.split(",")] if args.seeds else list(SEEDS)

    jobs_list: list[tuple[str, int, Path]] = []
    for g in genes:
        base_cfg = resolve_base_config(g)
        for s in seeds:
            variant = write_variant(g, base_cfg, s)
            verify_variant(base_cfg, variant, g, s)
            jobs_list.append((g, s, variant))
    _log(f"{len(jobs_list)} (gene, seed) arms verified against their published base config")
    if args.write_configs:
        return 0

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{env.get('PATH', '')}"

    todo = list(jobs_list)
    running: list[tuple[str, int, subprocess.Popen, object]] = []
    done: dict[tuple[str, int], int] = {}
    t0 = time.monotonic()

    def launch(gene: str, seed: int, cfg: Path):
        tag = f"{gene.lower()}_seed{seed}"
        log_path = LOG_DIR / f"{tag}.log"
        script = (
            f"set -e\n"
            f"{PY} -m genomics genotype train {cfg}\n"
            f"{PY} -m genomics genotype test {cfg}\n"
        )
        fh = open(log_path, "w")
        p = subprocess.Popen(["bash", "-c", script], cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, env=env)
        return p, fh

    while todo or running:
        while todo and len(running) < args.jobs:
            g, s, cfg = todo.pop(0)
            run_dir_guess = RUNS_ROOT / f"{g.lower()}_seed{s}"
            existing = find_run_dir(g, s)
            if args.skip_done and existing is not None:
                _log(f"{g:<9} seed={s:<3} SKIP: already has test results at {existing}")
                done[(g, s)] = 0
                continue
            if args.dry_run:
                _log(f"{g:<9} seed={s:<3} DRY RUN: would train+test {cfg.name}")
                done[(g, s)] = 0
                continue
            p, fh = launch(g, s, cfg)
            running.append((g, s, p, fh))
            _log(f"{g:<9} seed={s:<3} start (pid {p.pid}) [{len(running)}/{args.jobs} slots]")
        time.sleep(5)
        for item in list(running):
            g, s, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                done[(g, s)] = p.returncode
                mins = (time.monotonic() - t0) / 60
                _log(f"{g:<9} seed={s:<3} done rc={p.returncode} "
                     f"({len(done)}/{len(jobs_list)} total, {mins:.1f} min elapsed)")

    bad = {k: v for k, v in done.items() if v != 0}
    if bad:
        _log(f"FAILED arms (rc != 0): {bad}")

    # Harvest a combined table: published seed 13 + every retrained seed.
    import pandas as pd
    rows = []
    pub = pd.read_csv(PUBLISHED_TABLE)
    for g in genes:
        r = pub[pub["gene"] == g]
        if len(r) == 1:
            rows.append({"gene": g, "seed": PUBLISHED_SEED,
                         "bal_acc": float(r["bal_acc"].iloc[0]),
                         "bal_acc_val": float(r["bal_acc_val"].iloc[0])})
    for g, s, _cfg in jobs_list:
        run_dir = find_run_dir(g, s)
        if run_dir is None:
            continue
        test_m = bal_acc_from(run_dir, "test")
        val_m = bal_acc_from(run_dir, "val")
        rows.append({"gene": g, "seed": s, "bal_acc": test_m["bal_acc"], "bal_acc_val": val_m["bal_acc"]})

    out_df = pd.DataFrame(rows).sort_values(["gene", "seed"])
    out_df.to_csv(OUT_TABLE, index=False)
    _log(f"wrote {len(out_df)} rows ({out_df['gene'].nunique()} genes x up to "
         f"{out_df.groupby('gene').size().max()} seeds) to {OUT_TABLE}")

    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
