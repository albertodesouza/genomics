#!/usr/bin/env python3
"""Rankings and the head-to-head verdict, from gene_ranking_<arm>.json."""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from scipy import stats

R = Path("/home/breno/I2CA/genomics/results/genotype_based_predictor/gwas_ranking")
arm = sys.argv[1] if len(sys.argv) > 1 else "uncorrected"
tag = sys.argv[2] if len(sys.argv) > 2 else "window"
f = R / (arm if arm.endswith(".json") else f"gene_ranking_{arm}.json")
d = json.loads(f.read_text())
genes, isp = d["genes"], d["is_panel"]
gs, cnn = d["gene_scores"][tag], d["cnn"]
cmp_ = d["comparisons"][tag]

METHODS = [("CNN_bal_acc", "CNN bal_acc", lambda g: cnn[g]),
           ("GWAS_S_sidak", "Sidak/LiJi", lambda g: gs[g]["S_sidak"]),
           ("GWAS_S_gates", "GATES", lambda g: gs[g]["S_gates"]),
           ("GWAS_S_meanchi2", "mean chi2", lambda g: gs[g]["S_meanchi2"]),
           ("GWAS_S_min", "min p cru", lambda g: gs[g]["S_min"])]

print(f"### arm={arm}  janela={tag}  lambda_GC={d['lambda_gc_windows']:.2f}   "
      f"{cmp_['n_panel']} painel vs {cmp_['n_control']} controle\n")
print("ORDENACOES (rank 1 = melhor por aquele metodo)\n")
hdr = f"{'#':>2}  " + "".join(f"{lbl:<26}" for _, lbl, _ in METHODS)
print(hdr); print("-" * len(hdr))
order = {}
for key, lbl, fn in METHODS:
    v = [(g, fn(g)) for g in genes]
    v.sort(key=lambda t: -t[1])
    order[key] = v
for i in range(len(genes)):
    row = f"{i+1:>2}  "
    for key, lbl, _ in METHODS:
        g, s = order[key][i]
        mark = "P" if isp[genes.index(g)] else "c"
        row += f"{g:<9}{s:9.3f} [{mark}]   "
    print(row)

print(f"\nSEPARACAO PAINEL vs CONTROLE  (AUC = P(gene do painel > controle); "
      f"teste exato, {cmp_['n_assignments']} atribuicoes)\n")
print(f"{'metodo':<14}{'AUC':>8}{'rank medio painel':>19}{'rank medio ctrl':>17}{'p exato 1-cauda':>17}")
for key, lbl, _ in METHODS:
    m = cmp_["per_method"][key]
    print(f"{lbl:<14}{m['auc']:8.4f}{m['mean_rank_panel']:19.2f}"
          f"{m['mean_rank_control']:17.2f}{m['p_exact_one_sided']:17.4f}")

print("\nCNN CONTRA CADA ESTIMADOR DE GWAS  (delta AUC, teste exato pareado)\n")
print(f"{'GWAS':<14}{'AUC GWAS':>10}{'dAUC':>9}{'p exato 2-caudas':>18}"
      f"{'p exato 1-cauda':>17}{'DeLong z':>10}{'DeLong p':>10}{'rho vs CNN':>12}")
for key, lbl, _ in METHODS[1:]:
    c = cmp_["cnn_vs_gwas"][key]
    print(f"{lbl:<14}{cmp_['per_method'][key]['auc']:10.4f}{c['delta_auc_cnn_minus_gwas']:9.4f}"
          f"{c['p_exact_two_sided']:18.4f}{c['p_exact_one_sided_cnn_higher']:17.4f}"
          f"{c['delong_z']:10.2f}{c['delong_p_two_sided']:10.4f}{c['spearman_rho_vs_cnn']:12.3f}")

print("\nMEDIANAS")
for key, lbl, fn in METHODS:
    p = [fn(g) for g, b in zip(genes, isp) if b]
    c = [fn(g) for g, b in zip(genes, isp) if not b]
    print(f"  {lbl:<14} painel {np.median(p):9.3f}   controle {np.median(c):9.3f}   "
          f"sobreposicao: {sum(1 for x in c if x > min(p))}/{len(c)} controles acima do pior painel")
