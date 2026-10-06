#!/usr/bin/env bash
# Eleven-vs-eleven specificity: build the eight remaining control genes and put each
# through the same single-gene pipeline as the fourteen already measured.
#
# ORDERING. The window build is API-bound and touches no GPU; the running single-gene
# sweep is GPU-bound. They are started together on purpose -- serialising them would add
# the sweep's remaining hours to a job that is already the night's long pole. The eight
# new arms wait for BOTH.
#
# COST. 1072 individuals x 8 genes x 2 haplotypes = 17,152 AlphaGenome calls for the
# windows, plus 162 x 8 x 2 x 2 methods = 5,184 for the knockdown and null scrambles.
# ~22,300 calls, ~65 GB on disk, and ~6 h of GPU for the eight arms.
#
# RESUMABLE at every stage: the build skips samples in its checkpoint, the prefetch skips
# cached keys, and the sweep skips arms whose knockdown CSV exists.
set -uo pipefail
cd /home/breno/I2CA/genomics
PY=/home/breno/miniforge3/envs/genomics/bin/python3
# The builder shells out to bare `samtools` and `bcftools`, so the conda env must be on
# PATH for the workers, not merely used to launch python. Without this the reference
# window extraction fails per sample with FileNotFoundError and the build burns the
# cohort returning failures.
export PATH="/home/breno/miniforge3/envs/genomics/bin:$PATH"
LOG_DIR=results/genotype_based_predictor/logs
STAMP=$(date +%Y%m%d_%H%M%S)
GENES="PPP1R3E,ECHDC3,FRA10AC1,SPRED2,EIF1B,LACTB2,LRRC36,PRSS55"
DS=/dados/GENOMICS_DATA/v1/1kG_high_coverage
NEED_WINDOWS=1072
MIN_FREE_GB=120

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%S%z)] $*"; }

# The window builder resolves ${ALPHAGENOME_API_KEY} from the process environment and,
# unlike the replay and prefetch scripts, never calls load_dotenv itself. Extract just
# that one variable rather than sourcing ~/.env, which is not guaranteed to be valid
# shell. The key is never echoed; only its length is reported.
export ALPHAGENOME_API_KEY="$(awk -F= '/^ALPHAGENOME_API_KEY=/{sub(/^[^=]*=/,"");gsub(/^["'"'"']|["'"'"']$/,"");print;exit}' "$HOME/.env")"
if [ -z "${ALPHAGENOME_API_KEY:-}" ]; then
  echo "ABORT: ALPHAGENOME_API_KEY not resolvable from ~/.env"; exit 1
fi
free_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }

log "=== control expansion: 8 genes -> 11 panel vs 11 control ==="
log "genes: $GENES"
log "cost: ~17,152 window calls + ~5,184 scramble calls, ~65 GB, ~6 h GPU for the arms"

if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
  log "ABORT: $(free_gb) GB free, need >= $MIN_FREE_GB"; exit 1
fi

# ---- stage 1: window build (starts now, in parallel with the running sweep) ----------
log "stage 1: window build starting (4 workers)"
$PY scripts/experiments/specificity_control_build_parallel.py \
    --config configs/workflows/non_longevous_dataset/specificity_control_genes_8.yaml \
    --workers 4 > "$LOG_DIR/control8_build_$STAMP.log" 2>&1
BUILD_RC=$?
log "stage 1: build exited rc=$BUILD_RC"

# Completeness is judged from the windows on disk, not from the builder's exit code: a
# partial build that returns 0 would otherwise feed a silently truncated cohort forward.
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
log "stage 1: all 8 genes complete at $NEED_WINDOWS individuals"

# ---- stage 2: wait for the running single-gene sweep -----------------------------------
log "stage 2: waiting for the 14-gene sweep to finish"
while pgrep -f "scripts/experiments/single_gene_sweep.py" > /dev/null; do sleep 120; done
DONE14=$(ls results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout_single_*.csv 2>/dev/null | wc -l)
log "stage 2: sweep finished with $DONE14 arms measured"

# ---- stage 3: scrambles for the eight new genes ---------------------------------------
log "stage 3: scramble prefetch (knockdown + null)"
$PY scripts/experiments/scramble_prefetch.py --genes "$GENES" \
    > "$LOG_DIR/control8_prefetch_$STAMP.log" 2>&1
PF_RC=$?
log "stage 3: prefetch exited rc=$PF_RC"
if [ "$PF_RC" -ne 0 ]; then
  log "ABORT: prefetch reported failures; the replay would emit an incomplete CSV."
  log "Rerun this script -- cached keys are skipped, only the failures are retried."
  exit 1
fi

# ---- stage 4: one classifier per control gene -------------------------------------------
log "stage 4: eight single-gene arms"
# NOT --write-configs: that flag returns before training. Missing configs are
# written by the sweep itself, with the same drift verification.
$PY scripts/experiments/single_gene_sweep.py --genes "$GENES" \
    > "$LOG_DIR/control8_sweep_$STAMP.log" 2>&1
log "stage 4: sweep exited rc=$?"

# ---- stage 5: aggregate -----------------------------------------------------------------
log "stage 5: aggregate report"
$PY scripts/experiments/single_gene_sweep_report.py > "$LOG_DIR/control8_report_$STAMP.log" 2>&1
tail -40 "$LOG_DIR/control8_report_$STAMP.log"
log "=== control expansion finished; $(free_gb) GB free ==="
