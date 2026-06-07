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
    genotype_based/
    variant_transformer/
    snp_ancestry/
  workflows/
    non_longevous_dataset/
    longevity_dataset/
    alphagenome/
  converters/
    vcf_to_23andme/

docs/
  components/
  concepts/
  getting-started/
  guides/
  historical/
    alphagenome/
    1000-genomes/
  operations/
    monster.md
    scripts.md
  reference/

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
- Legacy-only neural ancestry configs belong under `legacy/neural_ancestry_predictor_deprecated/configs/`, not under `configs/`.
- Historical reproducibility code belongs under `legacy/`.
- Historical notes and migration records belong under `docs/historical/`.
- External modified code belongs under `third_party/`.
- Native non-Python tools belong under `native/`.
- Root-level wrappers and old package directories should not be reintroduced.
