#!/usr/bin/env bash
# 14-gene specificity control: retrain the tracks-only classifier on the published
# eleven-gene pigmentation panel PLUS three magnitude-matched, pigmentation-
# irrelevant windows (TPM2, SMCR8, PSMC4), so the knockdown probe can be run on a
# gene that is unused BY CONSTRUCTION.
#
# WHAT THIS TESTS
# The probe's validity instrumentation establishes that a flat response indicates
# non-reliance only where delivered perturbation is comparable to a responding
# gene. The within-panel version of that control exists (OCA2 receives 1.5x
# MC1R's delivered perturbation and moves the decision 3% as much). The
# between-panel version does not: no gene with *no role in the phenotype* has
# ever been knocked down in this pipeline, because a control gene must be inside
# the classifier's input to have any Delta at all.
#
# WHY raw_center_crop AND NOT DITA
# The gene-level readout is a per-gene scalar and carries no positional index;
# the alignment section of the probe paper shows an unaligned classifier recovers
# the aligned model's ranking at Spearman rho = 0.83 with signs agreeing 10/11.
# No DITA axis exists for the control windows. The published unaligned tracks-only
# reference is 0.926/0.923 on pigmentation.
#
# PRE-REGISTERED READING -- fixed before the knockdown is run.
#   - all three controls flat, |Delta_in| confirmed in the responding range
#       -> specificity established; the probe's first positive claim.
#   - any control responding at panel-gene magnitude
#       -> the probe responds to delivered perturbation at any well-expressed
#          locus and the reliance reading of the ranking is in trouble.
#   - mixed -> report per-gene against delivered magnitude; never average.
# Report |Delta| and |Delta_in| together for all fourteen genes, always.
#
# LIKE-FOR-LIKE. The retrain changes the classifier, so the eleven PUBLISHED
# Delta values are not comparable to this run's. All fourteen are measured here.
#
# Usage:
#   scripts/experiments/specificity_14gene.sh            # train + test
#   DRY_RUN=1 scripts/experiments/specificity_14gene.sh  # check guards only
set -euo pipefail

REPO_ROOT="/home/breno/I2CA/genomics"
PY="/home/breno/miniforge3/envs/genomics/bin/python3"
DATASET="/dados/GENOMICS_DATA/v1/1kG_high_coverage"
CONFIG="configs/predictors/genotype_based/pigmentation/pigmentation_binary_no_alignment_14gene.yaml"
LOG_DIR="results/genotype_based_predictor/logs"
# A 14-gene signals_only view is ~23 GB (the 11-gene one is 18 GB). Filling the
# disk mid-run has already cost one experiment here (2026-08-31), and the
# training CLI has no free-space guard of its own.
MIN_FREE_GB=32
EXPECTED_INDIVIDUALS=1072
GENES="MC1R TYRP1 TYR SLC45A2 DDB1 EDAR MFSD12 OCA2 HERC2 SLC24A5 TCHH TPM2 SMCR8 PSMC4"

cd "$REPO_ROOT"
export PATH="/home/breno/miniforge3/envs/genomics/bin:$PATH"   # samtools/bcftools
mkdir -p "$LOG_DIR"

[ -f "$CONFIG" ] || { echo "missing config: $CONFIG" >&2; exit 2; }

# The run name is derived from architecture and layout only -- it does NOT
# include the gene list -- so this run would resolve to the same directory as the
# published 11-gene no_alignment run and overwrite it. That happened on
# 2026-09-06 with a 3-gene arm and destroyed the published checkpoints.
if ! grep -q 'results_dir:.*runs_specificity_14gene' "$CONFIG"; then
    echo "ABORTED: $CONFIG does not set an isolated results_dir." >&2
    echo "         Without it this run clobbers the published no_alignment run." >&2
    exit 1
fi

echo "=== 14-gene specificity control ==="
echo "window availability:"
missing=0
for gene in $GENES; do
    n=$(find "$DATASET/individuals" -maxdepth 3 -type d -name "$gene" 2>/dev/null | wc -l)
    printf '  %-8s %5d / %d\n' "$gene" "$n" "$EXPECTED_INDIVIDUALS"
    [ "$n" -ge "$EXPECTED_INDIVIDUALS" ] || missing=1
done
if [ "$missing" -ne 0 ]; then
    echo "ABORTED: windows incomplete; training now would fit on a partial cohort." >&2
    exit 1
fi

avail=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
echo "disk free: ${avail} GB (need >= ${MIN_FREE_GB})"
if [ "$avail" -lt "$MIN_FREE_GB" ]; then
    echo "ABORTED: too little disk for a ~23 GB tensor cache." >&2
    exit 1
fi

if [ "${DRY_RUN:-0}" != "0" ]; then
    echo "DRY_RUN: guards pass; would train $CONFIG"
    exit 0
fi

stamp=$(date +%Y%m%d_%H%M%S)
log="$LOG_DIR/specificity_14gene_${stamp}.log"
echo "training -> $log"
"$PY" -m genomics genotype train "$CONFIG" > "$log" 2>&1

# `train` writes only validation metrics into the run directory; the test split
# is evaluated by the separate `test` subcommand.
echo "evaluating test split -> ${log%.log}.test.log"
"$PY" -m genomics genotype test "$CONFIG" > "${log%.log}.test.log" 2>&1

echo "done. run dir: results/genotype_based_predictor/runs_specificity_14gene"
echo "NEXT: knockdown over all fourteen genes on this checkpoint, reporting"
echo "      |Delta| and |Delta_in| together. See scripts/experiments/bulk_knockout_pigmentation.py"
