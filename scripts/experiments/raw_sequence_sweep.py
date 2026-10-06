#!/usr/bin/env python3
"""Full sweep for the raw-sequence (no AlphaGenome) control, all 42 genes.

Two phases, run in order:
  1. export  -- CPU/IO only (DynamicIndelAligner + one-hot encoding), no GPU needed.
               Safe to run concurrently with GPU-bound sweeps.
  2. train   -- full CNN training+eval on the one-hot tensors, GPU-bound like the
               published arms; empirically caps at 3 concurrent jobs on this box.

Usage:
  python3 scripts/experiments/raw_sequence_sweep.py export --jobs 4
  python3 scripts/experiments/raw_sequence_sweep.py train --jobs 3
  python3 scripts/experiments/raw_sequence_sweep.py train --write-configs
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

CONFIG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/raw_sequence"
RUNS_ROOT = REPO / "results/genotype_based_predictor/runs_rawsequence"
ONEHOT_CACHE_DIR = REPO / "results/cache/genotype_based_predictor/raw_sequence_onehot"
LOG_DIR = REPO / "results/genotype_based_predictor/logs/rawsequence_sweep"
OUT_TABLE = REPO / "results/genotype_based_predictor/raw_sequence_final_table.csv"

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
    "dataset_input.results_dir", "wandb.run_name",
    "model.cnn2.kernel_stage1", "model.cnn2.stride_stage1",
    "dataset_input.tensor_layout",
}


_GENE_BY_UPPER = {g.upper(): g for g in GENES}


def _resolve_gene_casing(genes_arg: str) -> list[str]:
    """Map a --genes CLI value back to GENES' exact casing (C6orf52 is not all-caps;
    export_aligned_dna's window-directory lookup is case-sensitive on that gene)."""
    out = []
    for raw in genes_arg.split(","):
        g = raw.strip()
        out.append(_GENE_BY_UPPER.get(g.upper(), g))
    return out


def _log(m: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {m}", flush=True)


def resolve_base_config(gene: str) -> Path:
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {gene} resolves to {len(hits)} published configs, expected 1")
    return hits[0]


def write_train_variant(gene: str, base_cfg_path: Path) -> Path:
    cfg = yaml.safe_load(base_cfg_path.read_text())
    cfg["dataset_input"]["results_dir"] = f"results/genotype_based_predictor/runs_rawsequence/{gene.lower()}_rawseq"
    cfg["dataset_input"]["tensor_layout"] = "raw_center_crop"
    cfg.setdefault("wandb", {})["run_name"] = f"rawseq-{gene.lower()}"
    b = cfg["model"]["cnn2"]
    b["kernel_stage1"] = [4, b["kernel_stage1"][1]]
    b["stride_stage1"] = [4, b["stride_stage1"][1]]

    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    out = CONFIG_DIR / f"{gene.lower()}_rawseq.yaml"
    out.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return out


def verify_variant(base_cfg_path: Path, variant_path: Path, gene: str) -> None:
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
    if cfg["dataset_input"]["processed_cache_dir"] != base["dataset_input"]["processed_cache_dir"]:
        raise SystemExit(f"ABORT: {variant_path} must keep reading the published cache "
                          f"(only used here for labels/split, not features)")
    if cfg["model"]["cnn2"]["kernel_stage1"][0] != 4 or cfg["model"]["cnn2"]["stride_stage1"][0] != 4:
        raise SystemExit(f"ABORT: {variant_path} stage-1 kernel/stride height must be 4")


def find_run_dir(gene: str) -> Path | None:
    hits = sorted(RUNS_ROOT.glob(f"{gene.lower()}_rawseq/*/test_best_accuracy_results.json"))
    return hits[0].parent if len(hits) == 1 else None


def bal_acc_from(run_dir: Path, split: str) -> float:
    d = json.loads((run_dir / f"{split}_best_accuracy_results.json").read_text())
    pcm = d["per_class_metrics"]
    return (float(pcm["strong pigmentation"]["recall"]) + float(pcm["weak pigmentation"]["recall"])) / 2.0


def cmd_export(args) -> int:
    genes = _resolve_gene_casing(args.genes) if args.genes else list(GENES)
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
            if (ONEHOT_CACHE_DIR / f"{g}.npz").exists() and not args.force:
                _log(f"{g:<9} SKIP: one-hot cache already present")
                done[g] = 0
                continue
            cfg = resolve_base_config(g)
            log_path = LOG_DIR / f"export_{g.lower()}.log"
            fh = open(log_path, "w")
            cmd = [PY, str(REPO / "scripts/experiments/raw_sequence_export_cache.py"), str(cfg), "--gene", g]
            p = subprocess.Popen(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, env=env)
            running.append((g, p, fh))
            _log(f"{g:<9} export start (pid {p.pid}) [{len(running)}/{args.jobs} slots]")
        time.sleep(5)
        for item in list(running):
            g, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                done[g] = p.returncode
                mins = (time.monotonic() - t0) / 60
                _log(f"{g:<9} export done rc={p.returncode} ({len(done)}/{len(genes)} total, {mins:.1f} min elapsed)")
    bad = {g: rc for g, rc in done.items() if rc != 0}
    if bad:
        _log(f"FAILED exports (rc != 0): {bad}")
    return 1 if bad else 0


def cmd_train(args) -> int:
    genes = _resolve_gene_casing(args.genes) if args.genes else list(GENES)
    jobs_list = []
    for g in genes:
        base_cfg = resolve_base_config(g)
        variant = write_train_variant(g, base_cfg)
        verify_variant(base_cfg, variant, g)
        jobs_list.append((g, variant))
    _log(f"{len(jobs_list)} train arms verified against their published base config")
    if args.write_configs:
        return 0

    missing_onehot = [g for g, _ in jobs_list if not (ONEHOT_CACHE_DIR / f"{g.upper()}.npz").exists()]
    if missing_onehot:
        raise SystemExit(f"ABORT: run the export phase first for {missing_onehot}")

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{env.get('PATH', '')}"
    todo = list(jobs_list)
    running: list[tuple[str, subprocess.Popen, object]] = []
    done: dict[str, int] = {}
    t0 = time.monotonic()

    while todo or running:
        while todo and len(running) < args.jobs:
            g, cfg = todo.pop(0)
            if args.skip_done and find_run_dir(g) is not None:
                _log(f"{g:<9} SKIP: already has test results")
                done[g] = 0
                continue
            if args.dry_run:
                _log(f"{g:<9} DRY RUN: would train+eval {cfg.name}")
                done[g] = 0
                continue
            log_path = LOG_DIR / f"train_{g.lower()}.log"
            fh = open(log_path, "w")
            cmd = [PY, str(REPO / "scripts/experiments/train_raw_sequence.py"), str(cfg), "--gene", g]
            p = subprocess.Popen(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, env=env)
            running.append((g, p, fh))
            _log(f"{g:<9} train start (pid {p.pid}) [{len(running)}/{args.jobs} slots]")
        time.sleep(5)
        for item in list(running):
            g, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                done[g] = p.returncode
                mins = (time.monotonic() - t0) / 60
                _log(f"{g:<9} train done rc={p.returncode} ({len(done)}/{len(jobs_list)} total, {mins:.1f} min elapsed)")

    bad = {g: rc for g, rc in done.items() if rc != 0}
    if bad:
        _log(f"FAILED train arms (rc != 0): {bad}")
    if args.dry_run:
        return 0

    import pandas as pd
    rows = []
    for g, _cfg in jobs_list:
        run_dir = find_run_dir(g)
        if run_dir is None:
            continue
        rows.append({"gene": g, "bal_acc": bal_acc_from(run_dir, "test"),
                     "bal_acc_val": bal_acc_from(run_dir, "val")})
    out_df = pd.DataFrame(rows).sort_values("gene")
    out_df.to_csv(OUT_TABLE, index=False)
    _log(f"wrote {len(out_df)} rows to {OUT_TABLE}")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="phase", required=True)

    exp = sub.add_parser("export")
    exp.add_argument("--genes", default=None)
    exp.add_argument("--jobs", type=int, default=4)
    exp.add_argument("--force", action="store_true")
    exp.set_defaults(func=cmd_export)

    tr = sub.add_parser("train")
    tr.add_argument("--genes", default=None)
    tr.add_argument("--jobs", type=int, default=3)
    tr.add_argument("--write-configs", action="store_true")
    tr.add_argument("--dry-run", action="store_true")
    tr.add_argument("--skip-done", action="store_true", default=True)
    tr.set_defaults(func=cmd_train)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
