# Configuration Layout

This directory contains the canonical configuration layout for new commands.
Historical config files remain in their original module directories as compatibility assets.

## Canonical Paths

| Domain | Path |
|---|---|
| Genomes Analyzer workflow | `configs/genomes_analyzer/` |
| Genotype-based predictor | `configs/predictors/genotype_based/` |
| Genotype raw-center-crop legacy baselines | `configs/predictors/genotype_based/neural_legacy/` |
| Variant Transformer predictor | `configs/predictors/variant_transformer/` |
| SNP ancestry predictor | `configs/predictors/snp_ancestry/` |
| Non-longevous dataset builder | `configs/workflows/non_longevous_dataset/` |
| Longevity dataset legacy workflow | `configs/workflows/longevity_dataset/` |
| AlphaGenome workflow | `configs/workflows/alphagenome/` |
| VCF to 23andMe converter | `configs/converters/vcf_to_23andme/` |
| Deprecated neural ancestry legacy configs | `configs/legacy/neural_ancestry_predictor_deprecated/` |

Prefer these paths in new commands and documentation. Existing module-local configs are kept until external scripts and historical notebooks are migrated.
