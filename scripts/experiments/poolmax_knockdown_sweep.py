#!/usr/bin/env python3
"""Replay the promoter knockdown on the 20 max-pool single-gene CNNs.

Each arm reads its own checkpoint and its own config from the run directory the poolmax sweep
wrote, so no arm can pick up another gene's model: the replay refuses to run if the config's
gene list is not exactly the arm's gene, and if the resolved run directory does not name it.

Costs nothing at the AlphaGenome API -- every scrambled prediction comes from the on-disk cache,
verified complete (6,480/6,480 entries) before this was launched.

One anchoring method only (biology_tss). The panel CSV carries three and the control sources
one or two, so an unfiltered replay would average a different number of scrambles per gene and
the panel/control contrast would partly measure that.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

REPO = Path("/home/breno/I2CA/genomics")
PY = "/home/breno/miniforge3/envs/genomics/bin/python"
OUT = REPO / "results/genotype_based_predictor/knockout_bulk/poolmax"
LOGS = REPO / "results/genotype_based_predictor/knockout_bulk/poolmax/logs"

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
CONTROL = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
           "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
# random_11_2, the second pre-existing random draw, taken whole. Both draws predate the
# specificity question and are independent of it, which is the only thing making the
# comparison interpretable; a draw is never reordered, substituted or trimmed.
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]
CONTROL = CONTROL + DRAW2



def find_config(gene: str) -> Path:
    """The config the poolmax sweep wrote next to this gene's checkpoint."""
    hits = sorted((REPO / "results/genotype_based_predictor").glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/config.yaml"))
    if not hits:
        raise SystemExit(f"ABORT: no poolmax config for {gene}")
    if len(hits) > 1:
        raise SystemExit(f"ABORT: {gene} has {len(hits)} poolmax configs, ambiguous:\n"
                         + "\n".join(f"  {h}" for h in hits))
    ck = hits[0].parent / "models" / "best_accuracy.pt"
    if not ck.exists():
        raise SystemExit(f"ABORT: {gene} config found but no checkpoint at {ck}")
    return hits[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4,
                    help="parallel arms. Each holds one dataset view and one 46k-param model; "
                         "unlike training there are no dataloader workers, so this is bounded "
                         "by RAM for the tensor views, not by cores.")
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    genes = a.genes or (PANEL + CONTROL)
    cfgs = {g: find_config(g) for g in genes}   # fail before launching anything
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{env.get('PATH','')}"

    print(f"{len(genes)} arms, {a.jobs} at a time", flush=True)
    t0 = time.monotonic()
    running: list[tuple[str, subprocess.Popen, object]] = []
    todo = list(genes)
    done: dict[str, int] = {}

    while todo or running:
        while todo and len(running) < a.jobs:
            g = todo.pop(0)
            out_csv = OUT / f"single_{g.lower()}.csv"
            if out_csv.exists():
                out_csv.unlink()            # the replay appends; a stale file would merge runs
            cmd = [PY, str(REPO / "scripts/experiments/single_gene_knockdown_replay.py"),
                   "--arm", g, "--config", str(cfgs[g]),
                   "--method", "biology_tss", "--out", str(out_csv)]
            if a.limit:
                cmd += ["--limit", str(a.limit)]
            fh = open(LOGS / f"{g.lower()}.log", "w")
            p = subprocess.Popen(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, env=env)
            running.append((g, p, fh))
            print(f"[{time.monotonic()-t0:7.1f}s] start {g} (pid {p.pid})", flush=True)
        time.sleep(2)
        for item in list(running):
            g, p, fh = item
            if p.poll() is not None:
                fh.close()
                running.remove(item)
                done[g] = p.returncode
                print(f"[{time.monotonic()-t0:7.1f}s] done  {g} rc={p.returncode} "
                      f"({len(done)}/{len(genes)})", flush=True)

    bad = {g: rc for g, rc in done.items() if rc != 0}
    print(f"\nsweep completo em {(time.monotonic()-t0)/60:.1f} min")
    if bad:
        print(f"FALHARAM (rc != 0): {bad}")
        print("rc=1 tambem significa 'entradas de cache faltando'; ver o log do braco.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
