#!/usr/bin/env bash
# Status dos bracos de treino em curso e dos que terminaram, em acuracia balanceada.
G=/home/breno/I2CA/genomics
R=$G/results/genotype_based_predictor
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo
echo "--- processos de treino ---"
ps -eo pid,etime,pcpu,rss,cmd --no-headers 2>/dev/null | grep "genomics.predictors.genotype_based.experiments.train" | grep -v grep \
  | awk '{printf "  pid=%-8s %8s  %5.1f%%  %5.1fGB  %s\n", $1, $2, $3, $4/1048576, $NF}' || echo "  (nenhum)"
echo "  carga: $(uptime | sed 's/.*load average: //')"
echo "  disco: $(df -h /dados 2>/dev/null | awk 'NR==2{print $4" livres"}')"
echo
echo "--- fase de cada braco ativo ---"
for f in $R/logs/collapse_fix_mc1r_*.log $R/logs/fixed_sweep_*.log; do
  [ -f "$f" ] || continue
  case "$f" in *_MAIN.log|*launcher.log) continue;; esac
  n=$(basename "$f" .log)
  last=$(grep -E "Computando normaliza|Epoch|epoca|Treinamento|Teste|Traceback" "$f" 2>/dev/null | tail -1 | cut -c1-70)
  printf "  %-28s %s\n" "$n" "${last:-iniciando}"
done
echo
echo "--- resultados (acuracia balanceada) ---"
for pat in "runs_collapse_fix/*" "runs_fixed_sweep/*" "runs_single_gene_*_melstranddita"; do
  $G/../genomics/../ 2>/dev/null
done
/home/breno/miniforge3/envs/genomics/bin/python3 $G/scripts/experiments/balanced_accuracy_report.py --glob "runs_collapse_fix/*" 2>/dev/null | sed -n '2,40p'
/home/breno/miniforge3/envs/genomics/bin/python3 $G/scripts/experiments/balanced_accuracy_report.py --glob "runs_fixed_sweep/*" 2>/dev/null | sed -n '2,40p'
