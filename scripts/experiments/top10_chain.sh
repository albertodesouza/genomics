#!/usr/bin/env bash
# Unattended continuation: wait for the top-10 train+test to finish, check the balanced
# accuracy, and only then run the knockdown and its report. Written so the whole thing
# completes with nobody watching.
#
# The balanced-accuracy gate is the same judgement I would apply by hand: a constant
# predictor scores exactly 0.5000 balanced (0.7099 raw), so anything at or near 0.5000
# means the run collapsed and perturbing it would produce numbers with no model behind
# them. On that outcome this stops and says so rather than filling a table.
set -uo pipefail
cd /home/breno/I2CA/genomics
export PATH=/home/breno/miniforge3/envs/genomics/bin:$PATH
PY=/home/breno/miniforge3/envs/genomics/bin/python
RUN=results/genotype_based_predictor/runs_top10
LOG=$RUN/TOP10_MAIN.log
say() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

say "aguardando o pipeline de treino/teste terminar"
while ! grep -qE "TOP10 PIPELINE rc=" "$LOG" 2>/dev/null; do
  if ! pgrep -f "experiments\.(train|evaluate_checkpoint)" >/dev/null 2>&1 \
     && ! grep -qE "TOP10 PIPELINE rc=" "$LOG" 2>/dev/null; then
    say "ABORTADO: o trainer sumiu sem marcador de conclusao"; exit 1
  fi
  sleep 30
done
RC=$(grep -oE "TOP10 PIPELINE rc=[0-9]+" "$LOG" | tail -1 | cut -d= -f2)
say "pipeline terminou com rc=$RC"
[ "$RC" != "0" ] && { say "ABORTADO: treino/teste falhou"; exit 1; }

say "lendo acuracia balanceada"
BAL=$($PY - <<'PYEOF'
import json, glob
f = glob.glob("results/genotype_based_predictor/runs_top10/top10_multigene/*/test_best_accuracy_results.json")
if not f:
    print("NONE"); raise SystemExit
d = json.load(open(f[0])); pcm = d["per_class_metrics"]
rs, rw = pcm["strong pigmentation"]["recall"], pcm["weak pigmentation"]["recall"]
cm = d["confusion_matrix"]
print(f"{(rs+rw)/2:.4f} {sum(cm[i][i] for i in range(2))/sum(map(sum,cm)):.4f} {rw:.4f} {rs:.4f}")
PYEOF
)
say "bal_acc acc rec_weak rec_strong = $BAL"
[ "$BAL" = "NONE" ] && { say "ABORTADO: sem resultado de teste"; exit 1; }
GATE=$($PY -c "print(1 if float('$(echo $BAL | cut -d' ' -f1)') > 0.55 else 0)")
if [ "$GATE" != "1" ]; then
  say "PARADO: acuracia balanceada indica colapso (preditor constante = 0.5000)."
  say "Nao rodei a perturbacao: ela mediria um modelo que nao aprendeu."
  exit 2
fi
say "modelo convergiu; seguindo para a perturbacao"

say "teste de fumaca do knockdown (2 individuos)"
rm -f /tmp/top10_smoke.csv
if ! $PY scripts/experiments/top10_knockdown.py --limit 2 --out /tmp/top10_smoke.csv \
     > $RUN/knockdown_smoke.log 2>&1; then
  say "ABORTADO: teste de fumaca falhou; ver $RUN/knockdown_smoke.log"; exit 1
fi
say "fumaca ok ($(( $(wc -l < /tmp/top10_smoke.csv) - 1 )) linhas)"

say "knockdown completo nos 10 genes"
rm -f results/genotype_based_predictor/knockout_bulk/top10/top10_knockdown.csv
$PY scripts/experiments/top10_knockdown.py > $RUN/knockdown_full.log 2>&1
say "knockdown rc=$? ; $(grep -c . results/genotype_based_predictor/knockout_bulk/top10/top10_knockdown.csv 2>/dev/null) linhas"

say "relatorio de regime"
$PY scripts/experiments/top10_regime_report.py > $RUN/REGIME_REPORT.txt 2>&1
cat $RUN/REGIME_REPORT.txt
say "CADEIA TOP10 COMPLETA"
