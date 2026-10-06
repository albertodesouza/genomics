#!/bin/bash
# check_superpopulation_training.sh - Hourly safety-net health check for the unattended
# superpopulation CNN comparison launched by run_superpopulation_cnn_comparison.py.
#
# Invoked hourly via a Claude Code CronCreate job (see the job created alongside this script),
# not a system crontab entry -- so it only fires while that Claude Code session is alive
# (even if idle). It is deliberately conservative: it only ever restarts the orchestrator if
# the orchestrator process itself has died while work remains (self-healing against a crashed
# *orchestrator*, distinct from the orchestrator's own per-config retry logic for a crashed
# *training run*). It never kills a live training process and never deletes anything.
#
# Prints "ALL_DONE" as the last line once status.json reports all_done=true, so the calling
# Claude Code prompt knows to stop scheduling further checks for this job.

set -uo pipefail

REPO_ROOT="/home/breno/I2CA/genomics"
GENOMICS_PYTHON="/home/breno/miniforge3/envs/genomics/bin/python3"
RUN_DIR="$REPO_ROOT/results/genotype_based_predictor/superpopulation_comparison"
STATUS_FILE="$RUN_DIR/status.json"
LOCK_FILE="$RUN_DIR/orchestrator.pid"
ORCHESTRATOR="$REPO_ROOT/scripts/experiments/run_superpopulation_cnn_comparison.py"
HEALTH_LOG="$RUN_DIR/health_check.log"
LAUNCH_LOG="$RUN_DIR/orchestrator_launch.log"

mkdir -p "$RUN_DIR"

log() {
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" >> "$HEALTH_LOG"
}

if [ ! -f "$STATUS_FILE" ]; then
    log "no status.json yet at $STATUS_FILE -- orchestrator may not have started; nothing to check"
    exit 0
fi

ALL_DONE=$("$GENOMICS_PYTHON" -c "import json,sys; d=json.load(open('$STATUS_FILE')); print(d.get('all_done', False))" 2>>"$HEALTH_LOG")

if [ "$ALL_DONE" = "True" ]; then
    log "status.json reports all_done=true; comparison finished"
    echo "ALL_DONE"
    exit 0
fi

FREE_GB=$(df --output=avail -BG "$REPO_ROOT" | tail -1 | tr -dc '0-9')
log "free disk on repo filesystem: ${FREE_GB}GB"
if [ -n "$FREE_GB" ] && [ "$FREE_GB" -lt 30 ]; then
    log "WARNING: free disk below 30GB (${FREE_GB}GB) -- orchestrator's own guard should abort further configs, but flagging here too"
fi

PID=""
if [ -f "$LOCK_FILE" ]; then
    PID=$(cat "$LOCK_FILE" 2>/dev/null)
fi

if [ -n "$PID" ] && kill -0 "$PID" 2>/dev/null; then
    log "orchestrator alive (pid $PID)"
    CURRENT=$("$GENOMICS_PYTHON" -c "
import json
d = json.load(open('$STATUS_FILE'))
for name, c in d.get('configs', {}).items():
    if c.get('state') == 'running':
        print(name)
        break
" 2>>"$HEALTH_LOG")
    if [ -n "$CURRENT" ]; then
        LOG_FILE="$RUN_DIR/logs/${CURRENT}.log"
        if [ -f "$LOG_FILE" ]; then
            LAST_MTIME=$(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)
            NOW=$(date +%s)
            AGE=$((NOW - LAST_MTIME))
            log "current config: $CURRENT, log last updated ${AGE}s ago"
            if [ "$AGE" -gt 5400 ]; then
                log "WARNING: no log activity for ${CURRENT} in over 90 minutes -- possible hang (not auto-killing; flagging only)"
            fi
        fi
    fi
    exit 0
fi

log "orchestrator NOT running (pid file: ${PID:-<none>}) but all_done is not true -- restarting"
cd "$REPO_ROOT" || exit 1
setsid nohup "$GENOMICS_PYTHON" "$ORCHESTRATOR" >> "$LAUNCH_LOG" 2>&1 < /dev/null &
disown
log "restarted orchestrator, new pid $!"
