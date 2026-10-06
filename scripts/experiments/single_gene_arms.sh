#!/usr/bin/env bash
# Single-gene arms: train the tracks-only classifier on ONE gene window, then knock that gene
# down on its own checkpoint. Two arms, SLC24A5 and TYR.
#
# WHAT THIS TESTS
# The panel knockdown measures RELIANCE, not information: a flat gene may carry nothing about
# the label, or carry something the classifier declined to use because a cheaper gene carries
# it too. Training on one gene removes the choice. Single-gene accuracy then measures the
# information the window carries; the single-gene knockdown measures the same 100 bp promoter
# scramble with reliance pinned at its maximum.
#
# WHY THESE TWO GENES
# Both respond in the eleven-gene panel, at opposite efficiency:
#   SLC24A5  |Delta| 0.503 at |Delta_in|  6,537  -> ratio 7.7e-5  (largest in the panel)
#   TYR      |Delta| 0.337 at |Delta_in| 54,297  -> ratio 6.2e-6  (12x less efficient)
#
# COSTS NOTHING AT THE ALPHAGENOME API. The 100 bp scrambles for both genes over all 162 test
# individuals are already on disk (notebooks/.cache/promoter_knockout_predictions, written by
# the published bulk run); the knockdown step is a replay through the new checkpoints.
#
# Usage:
#   scripts/experiments/single_gene_arms.sh              # both arms
#   scripts/experiments/single_gene_arms.sh slc24a5      # one arm
#   DRY_RUN=1 scripts/experiments/single_gene_arms.sh    # guards only
set -euo pipefail

REPO_ROOT="/home/breno/I2CA/genomics"
PY="/home/breno/miniforge3/envs/genomics/bin/python3"
DATASET="/dados/GENOMICS_DATA/v1/1kG_high_coverage"
LOG_DIR="results/genotype_based_predictor/logs"
MIN_FREE_GB=8          # a 1-gene signals_only view is ~1.6 GB (the 11-gene one is 18 GB)
EXPECTED_INDIVIDUALS=1072

cd "$REPO_ROOT"
export PATH="/home/breno/miniforge3/envs/genomics/bin:$PATH"
mkdir -p "$LOG_DIR"

ARMS="${*:-slc24a5 tyr}"

for arm in $ARMS; do
    case "$arm" in
        slc24a5) GENE=SLC24A5 ;;
        tyr)     GENE=TYR ;;
        *) echo "unknown arm: $arm (expected slc24a5|tyr)" >&2; exit 2 ;;
    esac
    CONFIG="configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment_single_${arm}.yaml"
    [ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 2; }

    # The run name is derived from architecture/layout only -- it does NOT include the gene list
    # -- so without an isolated results_dir this run resolves to the published no_alignment run's
    # directory and overwrites it. That happened on 2026-09-06 and destroyed a checkpoint.
    if ! grep -q "results_dir:.*runs_single_gene_${arm}" "$CONFIG"; then
        echo "ABORTED: $CONFIG does not set an isolated results_dir." >&2
        exit 1
    fi

    n=$(find "$DATASET/individuals" -maxdepth 3 -type d -name "$GENE" 2>/dev/null | wc -l)
    printf '=== single-gene arm %s (%s): %d / %d windows\n' "$arm" "$GENE" "$n" "$EXPECTED_INDIVIDUALS"
    [ "$n" -ge "$EXPECTED_INDIVIDUALS" ] || { echo "ABORTED: windows incomplete." >&2; exit 1; }

    avail=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
    echo "disk free: ${avail} GB (need >= ${MIN_FREE_GB})"
    [ "$avail" -ge "$MIN_FREE_GB" ] || { echo "ABORTED: too little disk." >&2; exit 1; }

    if [ "${DRY_RUN:-0}" != "0" ]; then
        echo "DRY_RUN: guards pass; would train $CONFIG"
        continue
    fi

    stamp=$(date +%Y%m%d_%H%M%S)
    log="$LOG_DIR/single_gene_${arm}_${stamp}.log"
    echo "training -> $log"
    "$PY" -m genomics genotype train "$CONFIG" > "$log" 2>&1
    echo "evaluating test split -> ${log%.log}.test.log"
    "$PY" -m genomics genotype test "$CONFIG" > "${log%.log}.test.log" 2>&1
    echo "done: results/genotype_based_predictor/runs_single_gene_${arm}"
done

echo "NEXT: scripts/experiments/single_gene_knockdown_replay.py (cached scrambles, no API cost)"
