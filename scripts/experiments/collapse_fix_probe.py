#!/usr/bin/env python3
"""Identifica a config que tira os bracos do plato da classe majoritaria.

Diagnostico que motiva isto (2026-09-11): 5 de 9 bracos do sweep
melanocyte+strand terminaram com matriz de confusao [[115,0],[47,0]] -- preditor
constante, acuracia balanceada 0.5. A informacao esta na entrada (logreg sobre
30 PCs da mesma entrada da 0.86-0.93 de cv-AUC), e nenhuma propriedade da
entrada prediz quem colapsa: MC1R e TYRP1 tem recuperabilidade identica
(0.9339 vs 0.9340) e desfechos opostos.

A causa medida e o orcamento de atualizacao. Com n_train=761 e batch_size=256 a
corrida inteira faz ceil(761/256)*400 = 1200 passos de otimizador. Os bracos que
escaparam levaram 111-615 passos; os cinco que nao escaparam precisavam de mais.
Replicacao da arquitetura confirma: sob LR constante o escape vem, sob
cosine T_max=400 e sob cosine_warm_restarts T_0=50 nao vem em 400 epocas.

Variantes, cada uma isolando um fator:
  A  batch 32, scheduler intacto        -> isola numero de passos (9600)
  B  batch 256, scheduler desligado     -> isola decaimento de LR (1200 passos)
  D  batch 32 + cosine T_max=num_epochs -> muitos passos E annealing pos-escape
                                           (a que eu mandaria para producao)
"""
import argparse, os, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO = Path("/home/breno/I2CA/genomics")
PY_BIN = "/home/breno/miniforge3/envs/genomics/bin/python3"
SRC_CFG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/single_gene_mel_strand_dita"
OUT_CFG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/collapse_fix"
LOG_DIR = REPO / "results/genotype_based_predictor/logs"

# O aligner bcftools_chain chama `bcftools` puro, entao o env conda tem de estar
# no PATH do subprocesso de treino -- nao basta lancar o python de la.
os.environ["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{os.environ.get('PATH', '')}"

VARIANTS = {
    "A": {"batch_size": 32, "sched": "keep"},
    "B": {"batch_size": 256, "sched": "off"},
    "D": {"batch_size": 32, "sched": "cosine"},
}


def log(msg, path=None):
    line = f"[{datetime.now(timezone.utc).isoformat()}] {msg}"
    print(line, flush=True)
    if path:
        with open(path, "a") as fh:
            fh.write(line + "\n")


def write_variant(gene: str, tag: str, spec: dict) -> Path:
    src = SRC_CFG_DIR / f"pigmentation_binary_single_{gene.lower()}_melstranddita.yaml"
    cfg = yaml.safe_load(src.read_text())
    di = cfg["dataset_input"]
    slug = f"{gene.lower()}_fix{tag}"
    # results_dir e cache isolados por variante: nenhum braco compartilha diretorio,
    # e caches separados evitam corrida de escrita entre variantes em paralelo.
    di["results_dir"] = f"results/genotype_based_predictor/runs_collapse_fix/{slug}"
    di["processed_cache_dir"] = f"results/cache/genotype_based_predictor/collapse_fix_{slug}"
    tr = cfg["training"]
    tr["batch_size"] = spec["batch_size"]
    if spec["sched"] == "off":
        tr["lr_scheduler"]["enabled"] = False
    elif spec["sched"] == "cosine":
        tr["lr_scheduler"]["enabled"] = True
        tr["lr_scheduler"]["type"] = "cosine"
        tr["lr_scheduler"]["T_max"] = tr["num_epochs"]
    cfg.setdefault("wandb", {})["run_name"] = f"collapse-fix-{tag}-{gene.lower()}"
    steps = -(-761 // spec["batch_size"]) * tr["num_epochs"]
    cfg["metadata"] = {"note": (
        f"COLLAPSE-FIX PROBE, variante {tag}, gene {gene}. batch_size={spec['batch_size']}, "
        f"scheduler={spec['sched']}, ~{steps} passos de otimizador em {tr['num_epochs']} epocas "
        f"(a config original fazia 1200). Gerado por scripts/experiments/collapse_fix_probe.py. "
        f"Identico ao braco melstranddita em split, seed, crop, layout e arquitetura."
    )}
    OUT_CFG_DIR.mkdir(parents=True, exist_ok=True)
    dst = OUT_CFG_DIR / f"{slug}.yaml"
    dst.write_text("# GERADO por scripts/experiments/collapse_fix_probe.py -- nao editar a mao\n"
                   + yaml.safe_dump(cfg, sort_keys=False))
    return dst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="MC1R")
    ap.add_argument("--variants", default="A,B,D")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    main_log = LOG_DIR / "collapse_fix_MAIN.log"
    tags = [t.strip() for t in a.variants.split(",") if t.strip()]
    log(f"collapse-fix probe: gene={a.gene} variantes={tags}", main_log)

    procs = []
    for tag in tags:
        spec = VARIANTS[tag]
        cfg = write_variant(a.gene, tag, spec)
        steps = -(-761 // spec["batch_size"]) * 400
        log(f"  variante {tag}: batch={spec['batch_size']} sched={spec['sched']} "
            f"~{steps} passos -> {cfg.name}", main_log)
        if a.dry_run:
            continue
        vlog = LOG_DIR / f"collapse_fix_{a.gene.lower()}_{tag}.log"
        fh = open(vlog, "w")
        p = subprocess.Popen([PY_BIN, "-m", "genomics", "genotype", "train", str(cfg)],
                             cwd=REPO, stdout=fh, stderr=subprocess.STDOUT)
        procs.append((tag, p, fh, vlog))
        log(f"  variante {tag} lancada, pid={p.pid}, log={vlog.name}", main_log)
        time.sleep(5)

    if a.dry_run:
        log("DRY RUN: nada lancado", main_log)
        return

    for tag, p, fh, vlog in procs:
        rc = p.wait()
        fh.close()
        log(f"  variante {tag} terminou rc={rc}", main_log)
    log("collapse-fix probe completo", main_log)


if __name__ == "__main__":
    main()
