#!/usr/bin/env bash
# Thirty-three-vs-nine specificity: build the third pre-existing random draw and put
# each gene through the same single-gene pipeline as the thirty-one already measured.
#
# WHY. Adding the second draw moved the variant-level gene score from AUC 0.667
# (p = 0.115) at 9v11 to AUC 0.737 (p = 0.021) at 9v22, and the accuracy ranking from
# 0.843 to 0.902 -- the earlier non-results were sample size, and this is the next step
# on the same curve. In the single-gene design adding controls costs no retrain of
# anything already measured: each gene carries its own classifier.
#
# SAME SETUP AS THE OTHER 31, WHICH IS THE WHOLE POINT. The windows are rebuilt here
# rather than reused from the draw's own directory: only 400 of the 1072 cohort
# individuals exist there, and AlphaGenome is not reproducible call to call (byte-
# identical input sequence, predictions differing at 944,619 of 3,145,728 positions).
# Reuse would have made these eleven arms non-comparable to the rest. See the
# 11c workflow config for the full reasoning.
#
# COST. 1072 individuals x 11 genes x 2 haplotypes = 23,584 AlphaGenome calls for the
# windows, plus 162 x 11 x 2 x 2 methods = 7,128 for the scrambles. ~30,700 calls,
# ~90 GB on disk, ~8 h of GPU for the eleven arms. The second batch of eight took
# 10.4 h wall for its build at 4 workers, so expect ~14 h here, then ~8 h of GPU:
# this does NOT finish overnight, it finishes late the following day.
#
# RESUMABLE at every stage: the build skips samples in its checkpoint, the prefetch
# skips cached keys, and the sweep skips arms whose knockdown CSV already exists.
# Rerunning this script after any failure picks up where it stopped and re-spends
# nothing.
set -uo pipefail
cd /home/breno/I2CA/genomics
PY=/home/breno/miniforge3/envs/genomics/bin/python3
# The builder shells out to bare `samtools` and `bcftools`, so the conda env must be on
# PATH for the workers, not merely used to launch python. Without this the reference
# window extraction fails per sample with FileNotFoundError and the build burns the
# whole cohort returning failures.
export PATH="/home/breno/miniforge3/envs/genomics/bin:$PATH"
LOG_DIR=results/genotype_based_predictor/logs
STAMP=$(date +%Y%m%d_%H%M%S)
GENES="ATP11B,BCL3,C6orf52,FBXO5,FOXN2,HERC6,KIAA0319,RIDA,SEM1,SFMBT2,SUMF2"
DS=/dados/GENOMICS_DATA/v1/1kG_high_coverage
NEED_WINDOWS=1072
# ~90 GB for this batch; 150 leaves room for the tensor caches the arms build and
# tear down (1.6 GB live per arm) plus the checkpoints.
MIN_FREE_GB=150

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%S%z)] $*"; }

# The window builder resolves ${ALPHAGENOME_API_KEY} from the process environment and,
# unlike the replay and prefetch scripts, never calls load_dotenv itself. Extract just
# that one variable rather than sourcing ~/.env, which is not guaranteed to be valid
# shell. The key is never echoed; only its length is reported.
export ALPHAGENOME_API_KEY="$(awk -F= '/^ALPHAGENOME_API_KEY=/{sub(/^[^=]*=/,"");gsub(/^["'"'"']|["'"'"']$/,"");print;exit}' "$HOME/.env")"
if [ -z "${ALPHAGENOME_API_KEY:-}" ]; then
  echo "ABORT: ALPHAGENOME_API_KEY not resolvable from ~/.env"; exit 1
fi
log "api key resolved (${#ALPHAGENOME_API_KEY} chars)"
free_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }

log "=== control expansion 11c: third draw -> 9 panel vs 33 control ==="
log "genes: $GENES"
log "cost: ~23,584 window calls + ~7,128 scramble calls, ~90 GB, ~8 h GPU for the arms"

if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
  log "ABORT: $(free_gb) GB free, need >= $MIN_FREE_GB"; exit 1
fi

# ---- stage 1: window build ------------------------------------------------------------
log "stage 1: window build starting (4 workers), $(free_gb) GB free"
$PY scripts/experiments/specificity_control_build_parallel.py \
    --config configs/workflows/non_longevous_dataset/specificity_control_genes_11c.yaml \
    --workers 4 > "$LOG_DIR/control11c_build_$STAMP.log" 2>&1
BUILD_RC=$?
log "stage 1: build exited rc=$BUILD_RC, $(free_gb) GB free"

# Completeness is judged from the windows on disk, not from the builder's exit code: a
# partial build that returns 0 would otherwise feed a silently truncated cohort forward.
# Counted exactly as single_gene_sweep.py counts, so the two cannot disagree.
INCOMPLETE=""
for g in ${GENES//,/ }; do
  n=$(ls -d "$DS"/individuals/*/windows/"$g" 2>/dev/null | wc -l)
  log "  $g: $n/$NEED_WINDOWS windows"
  [ "$n" -lt "$NEED_WINDOWS" ] && INCOMPLETE="$INCOMPLETE $g"
done
if [ -n "$INCOMPLETE" ]; then
  log "ABORT: incomplete windows for:$INCOMPLETE"
  log "The build is resumable -- rerun this script and it will pick up from the checkpoint."
  exit 1
fi
log "stage 1: all 11 genes complete at $NEED_WINDOWS individuals"

# ---- stage 2: scrambles for the eleven new genes ---------------------------------------
log "stage 2: scramble prefetch (knockdown + null)"
$PY scripts/experiments/scramble_prefetch.py --genes "$GENES" \
    > "$LOG_DIR/control11c_prefetch_$STAMP.log" 2>&1
PF_RC=$?
log "stage 2: prefetch exited rc=$PF_RC"
if [ "$PF_RC" -ne 0 ]; then
  log "ABORT: prefetch reported failures; the replay would emit an incomplete CSV."
  log "Rerun this script -- cached keys are skipped, only the failures are retried."
  exit 1
fi

# ---- stage 3: one classifier per control gene -------------------------------------------
if [ "$(free_gb)" -lt 40 ]; then
  log "ABORT before training: only $(free_gb) GB free, the tensor caches need room."
  exit 1
fi
log "stage 3: eleven single-gene arms, $(free_gb) GB free"
# parallel_sweep_fixed.py, with the exact hyperparameters of the 31 arms the paper
# reports (batch 256, scheduler off, 600 epochs, global max pool). NOT
# single_gene_sweep.py, which is the older gpavg architecture at 49,522 parameters.
# max-parallel 3: each arm spawns 6 dataloader workers, so 3 is the ceiling on this box.
$PY scripts/experiments/parallel_sweep_fixed.py --genes "$GENES" \
    --batch-size 256 --sched off --epochs 600 --pool max \
    --tag poolmax_ctrl4 --max-parallel 3 \
    > "$LOG_DIR/control11c_sweep_$STAMP.log" 2>&1
log "stage 3: sweep exited rc=$?"

# ---- stage 4: aggregate -----------------------------------------------------------------
# Writes to a PREVIEW path, not over the paper's table: the results are shown before
# anything is substituted.
log "stage 4: aggregate into the 42-arm preview table"
$PY scripts/experiments/poolmax_final_table.py \
    --csv-out results/genotype_based_predictor/poolmax_final_table_42arm.csv \
    --json-out results/genotype_based_predictor/poolmax_final_table_42arm.json \
    > "$LOG_DIR/control11c_report_$STAMP.log" 2>&1 || \
    log "stage 4: table needs the DRAW3 list added to poolmax_final_table.py first"
tail -40 "$LOG_DIR/control11c_report_$STAMP.log"
log "=== control expansion 11c finished; $(free_gb) GB free ==="
