#!/usr/bin/env bash
set -euo pipefail

# Sequential overnight runner for docs/IMPORTANT_EXPERIMENTS.md.
# Defaults run only canonical/ready experiments and skip legacy or blocked ones.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-results/overnight_important_experiments/$TIMESTAMP}"
LOG_DIR="$RUN_ROOT/logs"
CONFIG_DIR="$RUN_ROOT/configs"
SUMMARY_FILE="$RUN_ROOT/summary.tsv"
MIN_FREE_GB="${MIN_FREE_GB:-50}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-0}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OVERNIGHT_NUM_EPOCHS="${OVERNIGHT_NUM_EPOCHS:-1}"
OVERNIGHT_ALIGNMENT_AXIS_SPLITS="${OVERNIGHT_ALIGNMENT_AXIS_SPLITS:-train}"
RUN_CNN="${RUN_CNN:-1}"
RUN_RAW_BASELINE="${RUN_RAW_BASELINE:-1}"
RUN_SKLEARN="${RUN_SKLEARN:-1}"
RUN_XGBOOST="${RUN_XGBOOST:-1}"
RUN_MLC_LEGACY="${RUN_MLC_LEGACY:-0}"
RUN_VALIDATIONS="${RUN_VALIDATIONS:-1}"
WANDB_MODE="${WANDB_MODE:-offline}"
PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"
export PYTHON_BIN WANDB_MODE PYTORCH_CUDA_ALLOC_CONF
LAST_STEP_STATUS=0

mkdir -p "$LOG_DIR" "$CONFIG_DIR"
printf "step\tstatus\tstarted_at\tended_at\tlog\n" > "$SUMMARY_FILE"

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
    log "No files were deleted. See cleanup candidates printed at the end of this script output."
    return 1
  fi
}

run_step() {
  local name="$1"
  shift
  local started ended status log_file
  started="$(date -Iseconds)"
  log_file="$LOG_DIR/${name}.log"
  log "Starting $name"
  if ! check_disk; then
    LAST_STEP_STATUS=1
    ended="$(date -Iseconds)"
    printf "%s\tfailed:disk\t%s\t%s\t%s\n" "$name" "$started" "$ended" "$log_file" >> "$SUMMARY_FILE"
    print_cleanup_candidates
    if [ "$CONTINUE_ON_ERROR" = "1" ]; then
      return 0
    fi
    return 1
  fi
  set +e
  "$@" 2>&1 | tee "$log_file"
  status="${PIPESTATUS[0]}"
  LAST_STEP_STATUS="$status"
  set -e
  ended="$(date -Iseconds)"
  if [ "$status" -eq 0 ]; then
    log "Completed $name"
    printf "%s\tok\t%s\t%s\t%s\n" "$name" "$started" "$ended" "$log_file" >> "$SUMMARY_FILE"
    return 0
  fi
  log "FAILED $name with exit code $status"
  printf "%s\tfailed:%s\t%s\t%s\t%s\n" "$name" "$status" "$started" "$ended" "$log_file" >> "$SUMMARY_FILE"
  if [ "$CONTINUE_ON_ERROR" = "1" ]; then
    return 0
  fi
  return "$status"
}

train_and_eval() {
  local name="$1"
  local config="$2"
  local run_config
  run_config="$(config_for_run "$name" "$config")"
  run_step "${name}_train" "$PYTHON_BIN" -m genomics_cli genotype train "$run_config"
  if [ "$LAST_STEP_STATUS" -ne 0 ]; then
    log "Skipping ${name}_eval_test_best_accuracy because training failed."
    return 1
  fi
  run_step "${name}_eval_test_best_accuracy" "$PYTHON_BIN" -m genomics_cli genotype evaluate "$run_config" --checkpoint best_accuracy --split test
}

config_for_run() {
  local name="$1"
  local config="$2"
  "$PYTHON_BIN" - "$name" "$config" "$CONFIG_DIR" "$OVERNIGHT_NUM_EPOCHS" "$OVERNIGHT_ALIGNMENT_AXIS_SPLITS" <<'PY'
from pathlib import Path
import re
import sys
import yaml

name = sys.argv[1]
src = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
epochs_raw = sys.argv[4]
alignment_axis_splits_raw = sys.argv[5]
payload = yaml.safe_load(src.read_text()) or {}
if epochs_raw:
    training = payload.setdefault("training", {})
    training["num_epochs"] = int(epochs_raw)
    training["validation_frequency"] = 1
dataset_input = payload.setdefault("dataset_input", {})
if dataset_input.get("tensor_layout") == "haplotype_channels":
    splits = [item.strip() for item in alignment_axis_splits_raw.split(",") if item.strip()]
    if not splits:
        raise ValueError("OVERNIGHT_ALIGNMENT_AXIS_SPLITS nao pode ficar vazio para tensor_layout='haplotype_channels'")
    dataset_input["alignment_axis_splits"] = splits
safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
out = out_dir / f"{safe_name}_{src.name}"
out.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
print(out)
PY
}

module_available() {
  "$PYTHON_BIN" - "$1" <<'PY'
import importlib.util
import sys
raise SystemExit(0 if importlib.util.find_spec(sys.argv[1]) else 1)
PY
}

print_cleanup_candidates() {
  cat <<'EOF'

Cleanup candidates if disk is full. Review sizes first; this script does not delete anything.

Low-risk local cleanup:
  du -sh .pytest_cache results/genotype_based_predictor/runs_smoke results/variant_transformer_predictor/runs_smoke results/cache/genotype_based_predictor/smoke_train 2>/dev/null
  find /home/breno/I2CA/genomics -type d -name __pycache__ -prune -print
  find /home/breno/I2CA/genomics -type f -name '*.pyc' -print

Potentially useful cache cleanup after checking sizes:
  du -sh ~/.cache/pip ~/.cache/conda ~/.conda/pkgs 2>/dev/null

Large project/data directories to inspect before deleting anything:
  du -xh /home/breno/I2CA/genomics | sort -h | tail -50
  du -xh /dados/GENOMICS_DATA | sort -h | tail -50

Do not delete without explicit confirmation:
  /dados/GENOMICS_DATA/_deprecated/top3_20260601
EOF
}

log "Overnight important experiments"
log "Run root: $RUN_ROOT"
log "Logs: $LOG_DIR"
log "Temp configs: $CONFIG_DIR"
log "Summary: $SUMMARY_FILE"
log "Minimum free disk: ${MIN_FREE_GB}GB"
log "Python: $PYTHON_BIN"
log "OVERNIGHT_NUM_EPOCHS: ${OVERNIGHT_NUM_EPOCHS:-full config}"
log "OVERNIGHT_ALIGNMENT_AXIS_SPLITS: $OVERNIGHT_ALIGNMENT_AXIS_SPLITS"
log "WANDB_MODE: $WANDB_MODE"
log "PYTORCH_CUDA_ALLOC_CONF: $PYTORCH_CUDA_ALLOC_CONF"

if [ "$RUN_VALIDATIONS" = "1" ]; then
  run_step "audit_configs" "$PYTHON_BIN" -m genomics_cli audit-configs --fail-on-active-legacy
  run_step "audit_data_1kg_high_coverage" "$PYTHON_BIN" -m genomics_cli audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 3 --fail-on-missing
fi

if [ "$RUN_CNN" = "1" ]; then
  # Fastest DITA jobs first: fewer/no AlphaGenome signal channels before full signal+masks.
  train_and_eval "dita_masks_only" "configs/predictors/genotype_based/genes_1000_all_3ontologies_masks_only.yaml"
  train_and_eval "dita_masks_snp" "configs/predictors/genotype_based/genes_1000_all_3ontologies_masks_snp.yaml"
  train_and_eval "dita_signals_only" "configs/predictors/genotype_based/genes_1000_all_3ontologies_signals_only.yaml"
  train_and_eval "dita_variant_signal_mask" "configs/predictors/genotype_based/genes_1000_all_3ontologies_variant_signal_mask.yaml"
  train_and_eval "dita_signals_and_masks" "configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml"
fi

if [ "$RUN_RAW_BASELINE" = "1" ]; then
  train_and_eval "raw_center_crop_h1_h2" "configs/predictors/genotype_based/neural_legacy/genes_1000_all.yaml"
fi

if [ "$RUN_SKLEARN" = "1" ]; then
  if module_available sklearn && module_available scipy && module_available joblib; then
    run_step "pca400_rf_train_eval" "$PYTHON_BIN" -m genomics_cli genotype train "$(config_for_run pca400_rf_raw_center_crop configs/predictors/genotype_based/neural_legacy/genes_1000_all_pca400_rf.yaml)"
  else
    log "Skipping PCA+RF: sklearn/scipy/joblib not all available. Set RUN_SKLEARN=0 or install dependencies."
  fi
fi

if [ "$RUN_XGBOOST" = "1" ]; then
  if module_available sklearn && module_available scipy && module_available joblib && module_available xgboost; then
    run_step "pca400_xgboost_train_eval" "$PYTHON_BIN" -m genomics_cli genotype train "$(config_for_run pca400_xgboost_raw_center_crop configs/predictors/genotype_based/neural_legacy/genes_1000_all_pca400_xgboost.yaml)"
  else
    log "Skipping PCA+XGBoost: sklearn/scipy/joblib/xgboost not all available. Set RUN_XGBOOST=0 or install dependencies."
  fi
fi

if [ "$RUN_MLC_LEGACY" = "1" ]; then
  if [ -e /dados/GENOMICS_DATA/top3 ]; then
    run_step "mlc_snp_indel_dita_windows" "$PYTHON_BIN" -m genomics snp-ancestry run --config configs/predictors/snp_ancestry/genotype_windows_all_variants_mle.yaml
  else
    log "Skipping MLC legacy: /dados/GENOMICS_DATA/top3 is missing/quarantined. Do not recreate/delete without deciding the migration path."
  fi
fi

log "Finished. Summary: $SUMMARY_FILE"
print_cleanup_candidates
