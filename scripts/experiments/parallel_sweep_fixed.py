#!/usr/bin/env python3
"""Sweep de 22 bracos em paralelo sob a config corrigida, reportado em acuracia balanceada.

Contexto (2026-09-11): o sweep sequencial rodava um braco por vez a ~1 core de 20,
e 5 dos 9 bracos terminaram em preditor constante (bal_acc 0.5000). A causa medida e
o orcamento de atualizacao: n_train=761 com batch_size=256 da 1200 passos de otimizador
na corrida inteira, e os bracos que escaparam do plato levaram 111-615 passos.

Este script corrige as duas coisas de uma vez: a config (batch menor -> mais passos,
mais annealing pos-escape) e o uso da maquina (bracos em paralelo, 1 core cada).

O cache de alinhamento bcftools vive em dataset_dir e e compartilhado por todas as
configs -- indexado por conjunto de amostras e gene, nao por hiperparametro de treino.
Os 11 controles nao o tem, e construi-lo e o custo dominante deles. Esse custo nao se
perde se a config mudar depois.
"""
import argparse, os, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO = Path("/home/breno/I2CA/genomics")
PY_BIN = "/home/breno/miniforge3/envs/genomics/bin/python3"
SRC_CFG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/single_gene_mel_strand_dita"
OUT_CFG_DIR = REPO / "configs/predictors/genotype_based/pigmentation/fixed_sweep"
LOG_DIR = REPO / "results/genotype_based_predictor/logs"

os.environ["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{os.environ.get('PATH', '')}"

PANEL = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "MFSD12", "OCA2", "HERC2", "SLC24A5"]
NOT_PIGMENTATION = ["EDAR", "TCHH"]          # morfologia de pelo, fora da analise
CONTROL = ["TPM2", "SMCR8", "PSMC4", "PPP1R3E", "ECHDC3", "FRA10AC1",
           "SPRED2", "EIF1B", "LACTB2", "LRRC36", "PRSS55"]


def log(msg, path):
    line = f"[{datetime.now(timezone.utc).isoformat()}] {msg}"
    print(line, flush=True)
    with open(path, "a") as fh:
        fh.write(line + "\n")


def write_config(gene: str, batch_size: int, sched: str, epochs: int | None, tag: str,
                 pool: str | None = None) -> Path:
    src = SRC_CFG_DIR / f"pigmentation_binary_single_{gene.lower()}_melstranddita.yaml"
    if not src.exists():
        raise FileNotFoundError(src)
    cfg = yaml.safe_load(src.read_text())
    slug = f"{gene.lower()}_{tag}"
    di = cfg["dataset_input"]
    di["results_dir"] = f"results/genotype_based_predictor/runs_{tag}/{slug}"
    di["processed_cache_dir"] = f"results/cache/genotype_based_predictor/{tag}_{slug}"
    tr = cfg["training"]
    tr["batch_size"] = batch_size
    if epochs:
        tr["num_epochs"] = epochs
    if sched == "off":
        tr["lr_scheduler"]["enabled"] = False
    elif sched == "cosine":
        tr["lr_scheduler"]["enabled"] = True
        tr["lr_scheduler"]["type"] = "cosine"
        tr["lr_scheduler"]["T_max"] = tr["num_epochs"]
    if pool:
        # global_pool_type avg -> max resgata bracos cujo sinal nao sobrevive a media
        # da janela de 32.768 posicoes (DDB1 tem AUC 0.465 pos-pooling, abaixo do acaso).
        # Medido: levanta painel E controle, preservando a separacao entre os dois.
        cfg["model"]["cnn2"]["global_pool_type"] = pool
    steps = -(-761 // batch_size) * tr["num_epochs"]
    cfg.setdefault("wandb", {})["run_name"] = f"{tag}-{gene.lower()}"
    cfg["metadata"] = {"note": (
        f"SWEEP CORRIGIDO ({tag}), gene {gene}, melanocito CL:1000458 + strand do gene, "
        f"eixo DITA. batch_size={batch_size}, scheduler={sched}, num_epochs={tr['num_epochs']} "
        f"-> ~{steps} passos de otimizador (a config que colapsou fazia 1200), "
        f"global_pool_type={cfg['model']['cnn2']['global_pool_type']}. "
        f"Gerado por scripts/experiments/parallel_sweep_fixed.py. Identico ao braco "
        f"melstranddita em split, seed, crop, layout, ontologia e arquitetura."
    )}
    OUT_CFG_DIR.mkdir(parents=True, exist_ok=True)
    dst = OUT_CFG_DIR / f"{slug}.yaml"
    dst.write_text("# GERADO por scripts/experiments/parallel_sweep_fixed.py -- nao editar a mao\n"
                   + yaml.safe_dump(cfg, sort_keys=False))
    return dst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="all", help="'all', 'panel', 'control', ou lista por virgula")
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--sched", default="cosine", choices=["keep", "off", "cosine"])
    ap.add_argument("--epochs", type=int, default=None)
    ap.add_argument("--tag", default="fixed_sweep")
    ap.add_argument("--pool", default=None, choices=["avg", "max"])
    ap.add_argument("--max-parallel", type=int, default=6)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if a.genes == "all":
        genes = PANEL + NOT_PIGMENTATION + CONTROL
    elif a.genes == "panel":
        genes = PANEL
    elif a.genes == "control":
        genes = CONTROL
    else:
        genes = [g.strip().upper() for g in a.genes.split(",") if g.strip()]

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    main_log = LOG_DIR / f"{a.tag}_MAIN.log"
    steps = -(-761 // a.batch_size) * (a.epochs or 400)
    log(f"sweep {a.tag}: {len(genes)} bracos, batch={a.batch_size}, sched={a.sched}, "
        f"~{steps} passos/braco, max_parallel={a.max_parallel}", main_log)
    log(f"  genes: {', '.join(genes)}", main_log)

    configs = {}
    for g in genes:
        configs[g] = write_config(g, a.batch_size, a.sched, a.epochs, a.tag, a.pool)
    log(f"  {len(configs)} configs escritas em {OUT_CFG_DIR.relative_to(REPO)}", main_log)
    if a.dry_run:
        log("DRY RUN: nada lancado", main_log)
        return

    pending = list(genes)
    running = {}
    t0 = time.time()
    while pending or running:
        while pending and len(running) < a.max_parallel:
            g = pending.pop(0)
            vlog = LOG_DIR / f"{a.tag}_{g.lower()}.log"
            fh = open(vlog, "w")
            # train e depois test na mesma shell: sem o passo de test nao existe
            # test_best_accuracy_results.json, e sem ele nao ha matriz de confusao --
            # que e de onde sai a acuracia balanceada.
            cmd = (f"{PY_BIN} -m genomics genotype train {configs[g]} && "
                   f"{PY_BIN} -m genomics genotype test {configs[g]}")
            p = subprocess.Popen(["bash", "-c", cmd],
                                 cwd=REPO, stdout=fh, stderr=subprocess.STDOUT)
            running[g] = (p, fh, time.time())
            log(f"  + {g:<9} lancado pid={p.pid} ({len(running)} rodando, {len(pending)} na fila)",
                main_log)
            time.sleep(3)
        time.sleep(20)
        for g, (p, fh, ts) in list(running.items()):
            rc = p.poll()
            if rc is None:
                continue
            fh.close()
            del running[g]
            mins = (time.time() - ts) / 60
            log(f"  - {g:<9} rc={rc} em {mins:.1f} min ({len(running)} rodando, "
                f"{len(pending)} na fila)", main_log)
    log(f"sweep {a.tag} completo em {(time.time()-t0)/60:.1f} min", main_log)


if __name__ == "__main__":
    main()
