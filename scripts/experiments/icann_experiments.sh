#!/usr/bin/env bash
set -Eeuo pipefail

# Sequential overnight runner for the ICANN genotype experiments.
# Each step gets its own log and a row in summary.tsv.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-results/experiments/icann/$TIMESTAMP}"
LOG_DIR="$RUN_ROOT/logs"
SUMMARY_FILE="$RUN_ROOT/summary.tsv"
MIN_FREE_GB="${MIN_FREE_GB:-30}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-0}"
RUN_VALIDATIONS="${RUN_VALIDATIONS:-1}"
RUN_PCA_VARIANCE="${RUN_PCA_VARIANCE:-1}"
RUN_RF="${RUN_RF:-1}"
RUN_XGBOOST="${RUN_XGBOOST:-1}"
RUN_CNN2="${RUN_CNN2:-1}"
RUN_SEARCH_RF_XGBOOST="${RUN_SEARCH_RF_XGBOOST:-1}"
RUN_SEARCH_CNN2="${RUN_SEARCH_CNN2:-1}"
RUN_Y_RANDOMIZATION="${RUN_Y_RANDOMIZATION:-1}"
GENOMICS_BIN="${GENOMICS_BIN:-genomics}"

export WANDB_MODE="${WANDB_MODE:-offline}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

mkdir -p "$LOG_DIR"
printf "step\tstatus\tstarted_at\tended_at\tlog\tcommand\n" > "$SUMMARY_FILE"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

free_gb() {
  df -BG "$ROOT_DIR" | awk 'NR == 2 {gsub("G", "", $4); print $4}'
}

check_disk() {
  local available
  available="$(free_gb)"
  if [ "$available" -lt "$MIN_FREE_GB" ]; then
    log "ERROR: only ${available}GB free at $ROOT_DIR; need at least ${MIN_FREE_GB}GB."
    return 1
  fi
}

run_step() {
  local name="$1"
  shift
  local started ended status log_file command_text
  started="$(date -Iseconds)"
  log_file="$LOG_DIR/${name}.log"
  command_text="$*"

  log "Starting $name"
  log "Command: $command_text"

  if ! check_disk; then
    ended="$(date -Iseconds)"
    printf "%s\tfailed:disk\t%s\t%s\t%s\t%s\n" "$name" "$started" "$ended" "$log_file" "$command_text" >> "$SUMMARY_FILE"
    if [ "$CONTINUE_ON_ERROR" = "1" ]; then
      return 0
    fi
    return 1
  fi

  set +e
  "$@" 2>&1 | tee "$log_file"
  status="${PIPESTATUS[0]}"
  set -e

  ended="$(date -Iseconds)"
  if [ "$status" -eq 0 ]; then
    log "Completed $name"
    printf "%s\tok\t%s\t%s\t%s\t%s\n" "$name" "$started" "$ended" "$log_file" "$command_text" >> "$SUMMARY_FILE"
    return 0
  fi

  log "FAILED $name with exit code $status"
  printf "%s\tfailed:%s\t%s\t%s\t%s\t%s\n" "$name" "$status" "$started" "$ended" "$log_file" "$command_text" >> "$SUMMARY_FILE"
  if [ "$CONTINUE_ON_ERROR" = "1" ]; then
    return 0
  fi
  return "$status"
}

run_if_enabled() {
  local enabled="$1"
  shift
  local name="$1"
  shift
  if [ "$enabled" = "1" ]; then
    run_step "$name" "$@"
  else
    log "Skipping $name"
  fi
}

log "ICANN genotype experiments"
log "Run root: $RUN_ROOT"
log "Logs: $LOG_DIR"
log "Summary: $SUMMARY_FILE"
log "Minimum free disk: ${MIN_FREE_GB}GB"
log "CONTINUE_ON_ERROR: $CONTINUE_ON_ERROR"
log "WANDB_MODE: $WANDB_MODE"

if [ "$RUN_VALIDATIONS" = "1" ]; then
  run_step audit_configs "$GENOMICS_BIN" audit-configs --fail-on-active-legacy
  run_step audit_data_1kg_high_coverage "$GENOMICS_BIN" audit-data --dataset-id 1kg_high_coverage --sample-limit 3 --fail-on-missing
fi

run_if_enabled "$RUN_PCA_VARIANCE" pca500_variance_randomized "$GENOMICS_BIN" genotype pca-variance configs/predictors/genotype_based/icann/search_rf_xgboost.yaml --output results/genotype_based_predictor/icann/pca_variance/pca500_variance_randomized.png --json-output results/genotype_based_predictor/icann/pca_variance/pca500_variance_randomized.json --max-components 500
run_if_enabled "$RUN_RF" rf_train "$GENOMICS_BIN" genotype train configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml
run_if_enabled "$RUN_XGBOOST" xgboost_train "$GENOMICS_BIN" genotype train configs/predictors/genotype_based/icann/genes_1000_all_xgboost.yaml
run_if_enabled "$RUN_CNN2" cnn2_train "$GENOMICS_BIN" genotype train configs/predictors/genotype_based/icann/genes_1000_all_cnn2.yaml
run_if_enabled "$RUN_SEARCH_RF_XGBOOST" rf_xgboost_search "$GENOMICS_BIN" genotype search configs/predictors/genotype_based/icann/search_rf_xgboost.yaml
run_if_enabled "$RUN_SEARCH_CNN2" cnn2_ablation_search "$GENOMICS_BIN" genotype search configs/predictors/genotype_based/icann/search_cnn2_ablation.yaml
run_if_enabled "$RUN_Y_RANDOMIZATION" cnn2_y_randomization_train "$GENOMICS_BIN" genotype train configs/predictors/genotype_based/icann/genes_1000_all_cnn2_y_randomization.yaml

log "Finished. Summary: $SUMMARY_FILE"
