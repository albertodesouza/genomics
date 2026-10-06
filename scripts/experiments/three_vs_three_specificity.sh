#!/usr/bin/env bash
# 3-vs-3 magnitude-matched specificity control: the tracks-side arms.
#
# WHAT THIS TESTS
# Audit test 2 of the paper -- the irrelevant-gene control -- has only ever been
# run against the dosage logistic regression, where eleven pigmentation-irrelevant
# windows reproduced the real panel exactly. It has never been run against the
# *tracks* pipeline, which is the pipeline whose phenotype specificity the
# predicted-transcriptome argument is defending. This script supplies the two
# tracks-only CNN arms of that test.
#
# THE PAIRING (optimal assignment, minimising total |log10(signal ratio)|; the
# proxy is reference-window crop signal from specificity_control_preflight.py,
# which tracks realised |Delta_in| over the published panel at rho = 0.855):
#
#     TPM2   48,071  <->  TYR      48,068   (1.00x)
#     SMCR8  23,951  <->  MFSD12   26,326   (0.91x)
#     PSMC4  18,620  <->  MC1R     14,995   (1.24x)
#     control total 90,642  vs  panel total 89,389   -- 1.4% apart
#
# THE GENOTYPE ROW IS ALREADY MEASURED (2026-09-06, CPU, seconds):
#     panel triple    dosage LR  1.0000 raw / 1.0000 PCA   (1222 variants)
#     control triple  dosage LR  0.9815 raw / 0.9938 PCA   ( 876 variants)
# Both triples' genotype content saturates the pigmentation proxy, so a gap that
# appears in the tracks row below cannot be explained by the control windows
# simply carrying less genotype information. That is what makes the 2x2
# interpretable, and it is why the LR arms are run first and reported alongside.
# Re-run them with:
#   python3 scripts/experiments/random_gene_genotype_control.py --task pigmentation \
#       --genes TPM2,SMCR8,PSMC4 --out results/.../three_vs_three_lr_control_pigmentation.json
#
# PRE-REGISTERED READING, fixed before either arm is trained:
#   - control arm approaches the panel arm  -> the tracks representation carries
#     ancestry as promiscuously as dosage does; the filter argument is refuted,
#     the same way the random-panel LR refuted the biological reading of Table I.
#   - control arm falls materially below    -> the first evidence of phenotype
#     specificity anywhere in this work.
#   - either way, report the full 2x2. Never the tracks row on its own.
#
# PREREQUISITE: the control windows must exist for all 1072 pigmentation
# individuals. The panel arm can run today; the control arm needs the
# specificity-control build to finish (see specificity_control_build_parallel.py).
# This script checks and refuses rather than training on a partial cohort.
#
# Usage:
#   scripts/experiments/three_vs_three_specificity.sh            # both arms
#   scripts/experiments/three_vs_three_specificity.sh panel      # one arm
#   DRY_RUN=1 scripts/experiments/three_vs_three_specificity.sh  # check guards only
set -euo pipefail

REPO_ROOT="/home/breno/I2CA/genomics"
PY="/home/breno/miniforge3/envs/genomics/bin/python3"
DATASET="/dados/GENOMICS_DATA/v1/1kG_high_coverage"
CONFIG_DIR="configs/predictors/genotype_based/pigmentation"
LOG_DIR="results/genotype_based_predictor/logs"
MIN_FREE_GB=25          # a 3-gene signals_only view is ~5 GB; refuse below this
EXPECTED_INDIVIDUALS=1072

cd "$REPO_ROOT"
export PATH="/home/breno/miniforge3/envs/genomics/bin:$PATH"   # samtools/bcftools
mkdir -p "$LOG_DIR"

arms=("${@:-panel control}")

free_gb() { df -BG --output=avail / | tail -1 | tr -dc '0-9'; }

check_windows() {
    # $1... gene symbols; every one must be present for the whole cohort.
    local missing=0
    for gene in "$@"; do
        local n
        n=$(find "$DATASET/individuals" -maxdepth 3 -type d -name "$gene" 2>/dev/null | wc -l)
        printf '    %-8s %5d / %d individuals\n' "$gene" "$n" "$EXPECTED_INDIVIDUALS"
        [ "$n" -ge "$EXPECTED_INDIVIDUALS" ] || missing=1
    done
    return $missing
}

for arm in $arms; do
    case "$arm" in
        panel)   genes="TYR MFSD12 MC1R"   ;;
        control) genes="TPM2 SMCR8 PSMC4"  ;;
        *) echo "unknown arm '$arm' (expected: panel, control)" >&2; exit 2 ;;
    esac
    config="$CONFIG_DIR/pigmentation_binary_no_alignment_triple_${arm}.yaml"
    [ -f "$config" ] || { echo "missing config: $config" >&2; exit 2; }

    echo "=== arm: $arm ($genes) ==="

    # The run name is derived from architecture and layout only -- it does NOT
    # include the gene list -- so a 3-gene raw_center_crop run resolves to the
    # same directory as the published 11-gene no_alignment run and overwrites its
    # provenance files. That happened on 2026-09-06. Refuse to launch unless the
    # config sends this arm to its own isolated tree.
    if ! grep -q 'results_dir:.*runs_three_vs_three/'"$arm" "$config"; then
        echo "  ABORTED: $config does not set an isolated results_dir." >&2
        echo "           Expected results_dir ending in runs_three_vs_three/$arm --" >&2
        echo "           without it this run would clobber the published no_alignment run." >&2
        exit 1
    fi
    echo "  window availability:"
    if ! check_windows $genes; then
        echo "  SKIPPED: windows incomplete for this arm; training now would fit on a" >&2
        echo "           partial cohort and would not be comparable to the other arm." >&2
        continue
    fi

    avail=$(free_gb)
    echo "  disk free: ${avail} GB (need >= ${MIN_FREE_GB})"
    if [ "$avail" -lt "$MIN_FREE_GB" ]; then
        echo "  ABORTED: too little disk. A training run materialises a new tensor cache" >&2
        echo "           and the CLI has no free-space guard -- filling the disk mid-run" >&2
        echo "           has already cost one experiment here (see paper TODO, 2026-08-31)." >&2
        exit 1
    fi

    if [ "${DRY_RUN:-0}" != "0" ]; then
        echo "  DRY_RUN: would train $config"
        continue
    fi

    stamp=$(date +%Y%m%d_%H%M%S)
    log="$LOG_DIR/three_vs_three_${arm}_${stamp}.log"
    echo "  training -> $log"
    "$PY" -m genomics genotype train "$config" > "$log" 2>&1

    # `train` writes only validation metrics into the run directory; the test
    # split is evaluated by the separate `test` subcommand. The report script
    # looks for exactly this output name, so keep it in sync with it.
    echo "  evaluating test split -> ${log%.log}.test.log"
    "$PY" -m genomics genotype test "$config" \
        --output-name "three_vs_three_${arm}_test_results" \
        > "${log%.log}.test.log" 2>&1
    echo "  done: $arm"
done

echo
echo "Assemble the 2x2 with: scripts/experiments/three_vs_three_report.py"
