#!/usr/bin/env bash
# Regenerate every paper figure and table from finished artefacts, then rebuild the PDF.
#
# Reads only what is already on disk: the gene-ranking JSONs, the poolmax final table,
# the knockdown replay CSVs and the processed tensor shards. Trains nothing, predicts
# nothing, and bills no AlphaGenome API call. Safe to run at any point -- an arm that
# has not finished is reported as skipped rather than faked.
set -euo pipefail

PY=/home/breno/miniforge3/envs/genomics/bin/python
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER=/home/breno/I2CA/paper-knockdown-probe
cd "$HERE"

echo "== cohort"                   && $PY fig_cohort.py
echo "== GWAS"                     && $PY fig_gwas.py
echo "== GWAS window Manhattan"    && $PY fig_gwas_window_manhattan.py
echo "== tables"                   && $PY tab_arms.py
echo "== accuracy violin"          && $PY fig_balacc_violin.py
echo "== rank comparison"          && $PY fig_rank_comparison.py
echo "== top-k admission"          && $PY fig_topk_admission.py
echo "== knockdown (single-gene)"  && $PY fig_delta_violin.py
echo "== knockdown (top-10)"       && $PY fig_top10_knockdown.py
echo "== appendix: class signal"   && $PY fig_app_class_mean_signal.py
echo "== appendix: pre/post"       && $PY fig_app_knockdown_prepost.py
echo "== appendix: per-individual arrows" && $PY fig_app_knockdown_individuals.py
echo "== appendix: expr vs response"  && $PY fig_app_expr_vs_response.py

echo "== LaTeX"
cd "$PAPER" && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex \
  | grep -E "Output written|Latexmk: All targets" || true
