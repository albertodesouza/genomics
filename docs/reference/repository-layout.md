# Repository Layout

## Canonical Layout

```text
src/genomics/
  cli.py
  workspace.py
  core/
  workflows/
  predictors/
  converters/

configs/
  genomes_analyzer/
  predictors/
  workflows/
  converters/
  legacy/

legacy/
  neural_ancestry_predictor_deprecated/
  neural_longevity_dataset/

scripts/
  env/
  ops/
  maintenance/
  diagnostics/
  experiments/
  dev/

native/
  genes_difference_count/

third_party/
  FROGAncestryCalc/
```

## Rules

- Active Python package code lives under `src/genomics/`.
- New commands should be exposed through `genomics ...`.
- Shared code belongs in `src/genomics/core/`.
- Canonical configs belong under `configs/`.
- Historical reproducibility code belongs under `legacy/`.
- External modified code belongs under `third_party/`.
- Native non-Python tools belong under `native/`.
- Root-level wrappers and old package directories should not be reintroduced.
