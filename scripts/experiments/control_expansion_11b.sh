#!/usr/bin/env bash
# Twenty-two-vs-eleven specificity: build the random_11_2 draw and put each gene
# through the same single-gene pipeline as the twenty-two already measured.
#
# WHY. At 11v11 the comparison gave p = 0.106 and the leave-one-out was unstable in
# BOTH directions -- LRRC36 (control) alone gives 0.042, HERC2 (panel) alone gives
# 0.049. That is the sample size, not an outlier. Doubling the control arm is the
# cheapest fix available, because the single-gene design gives each gene its own
# classifier: nothing already measured has to be retrained.
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
GENES="CD47,COA1,EGF,FAM234B,FYB1,HSH2D,LYNX1,OR11H12,OR51S1,TRHR,TSPAN11"
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

log "=== control expansion 11b: random_11_2 -> 11 panel vs 22 control ==="
log "genes: $GENES"
log "cost: ~23,584 window calls + ~7,128 scramble calls, ~90 GB, ~8 h GPU for the arms"

if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
  log "ABORT: $(free_gb) GB free, need >= $MIN_FREE_GB"; exit 1
fi

# ---- stage 1: window build ------------------------------------------------------------
log "stage 1: window build starting (4 workers), $(free_gb) GB free"
$PY scripts/experiments/specificity_control_build_parallel.py \
    --config configs/workflows/non_longevous_dataset/specificity_control_genes_11b.yaml \
    --workers 4 > "$LOG_DIR/control11b_build_$STAMP.log" 2>&1
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
    > "$LOG_DIR/control11b_prefetch_$STAMP.log" 2>&1
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
# NOT --write-configs: that flag returns before training. Missing configs are
# written by the sweep itself, with the same drift verification.
$PY scripts/experiments/single_gene_sweep.py --genes "$GENES" \
    > "$LOG_DIR/control11b_sweep_$STAMP.log" 2>&1
log "stage 3: sweep exited rc=$?"

# ---- stage 4: aggregate -----------------------------------------------------------------
log "stage 4: aggregate report"
$PY scripts/experiments/single_gene_sweep_report.py > "$LOG_DIR/control11b_report_$STAMP.log" 2>&1
tail -60 "$LOG_DIR/control11b_report_$STAMP.log"
log "=== control expansion 11b finished; $(free_gb) GB free ==="
