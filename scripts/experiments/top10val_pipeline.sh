#!/usr/bin/env bash
# Downstream of the top10val training: the promoter knockdown and its null band, both
# through the newly trained 10-gene model, then the comparison against the arm the paper
# currently reports.
#
# Costs nothing at the AlphaGenome API. Every scrambled prediction for these ten genes is
# already cached -- verified before launch: 324 promoter (biology_tss) and 324 null
# (random_windowd0 for the five panel genes, random_windowd1 for the five controls) per
# gene, i.e. 162 individuals x 2 haplotypes with no gap.
set -euo pipefail
PY=/home/breno/miniforge3/envs/genomics/bin/python
REPO=/home/breno/I2CA/genomics
cd "$REPO"

RUN=results/genotype_based_predictor/runs_top10val/top10val_multigene
CFG=$(ls "$RUN"/*/config.yaml)
OUT=results/genotype_based_predictor/knockout_bulk/top10val
mkdir -p "$OUT"

# Refuse to read a half-written run: the checkpoint and the test evaluation must both exist.
CK=$(dirname "$CFG")/models/best_accuracy.pt
[ -f "$CK" ] || { echo "ABORT: no checkpoint at $CK" >&2; exit 1; }
[ -f "$(dirname "$CFG")/test_best_accuracy_results.json" ] || {
    echo "ABORT: training did not finish its test evaluation" >&2; exit 1; }
echo "config: $CFG"

echo "== promoter knockdown (biology_tss)"
$PY scripts/experiments/top10_knockdown.py --config "$CFG" \
    --method biology_tss --out "$OUT/top10val_knockdown.csv"

echo "== null band (random window, one draw per individual/gene)"
$PY scripts/experiments/top10_knockdown.py --config "$CFG" \
    --method random_window,random_windowd1 --out "$OUT/top10val_null.csv"

echo "== done"; wc -l "$OUT"/*.csv
