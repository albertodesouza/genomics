# Configuration

Canonical configs live under `configs/`.

| Area | Path |
|---|---|
| Genomes Analyzer | `configs/genomes_analyzer/` |
| Genotype predictor | `configs/predictors/genotype_based/` |
| Variant transformer | `configs/predictors/variant_transformer/` |
| SNP ancestry | `configs/predictors/snp_ancestry/` |
| Non-longevous dataset builder | `configs/workflows/non_longevous_dataset/` |
| AlphaGenome workflow | `configs/workflows/alphagenome/` |
| VCF to 23andMe converter | `configs/converters/vcf_to_23andme/` |
| Legacy neural ancestry code | `legacy/neural_ancestry_predictor_deprecated/configs/` |

## Auditing Configs

Use the config audit before long runs:

```bash
genomics audit-configs
genomics audit-configs --fail-on-active-legacy
```

`--fail-on-active-legacy` is the recommended gate for new experiments. It permits deprecated configs to exist but fails if an active config still depends on legacy data paths.

## Overrides

Several CLI commands accept path overrides such as:

- `--dataset-dir`
- `--processed-cache-dir`
- `--results-dir`
- `--processed-dir`
- `--consensus-dataset-dir`

Overrides are applied only when explicitly passed. This keeps config files reproducible while allowing smoke tests and temporary runs.
