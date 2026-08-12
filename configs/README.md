# Configuration Layout

This directory contains the canonical configuration layout for new commands.
Historical config files remain in their original module directories as compatibility assets.

## Canonical Paths

| Domain | Path |
|---|---|
| Genomes Analyzer workflow | `configs/genomes_analyzer/` |
| Genotype-based predictor | `configs/predictors/genotype_based/` |
| Genotype ICANN paper configs | `configs/predictors/genotype_based/icann/` |
| Genotype raw-center-crop legacy baselines | `configs/predictors/genotype_based/neural_legacy/` |
| SNP ancestry ICANN paper configs | `configs/predictors/snp_ancestry/icann/` |
| Variant Transformer predictor | `configs/predictors/variant_transformer/` |
| SNP ancestry predictor | `configs/predictors/snp_ancestry/` |
| Non-longevous dataset builder | `configs/workflows/non_longevous_dataset/` |
| Longevity dataset legacy workflow | `configs/workflows/longevity_dataset/` |
| AlphaGenome workflow | `configs/workflows/alphagenome/` |
| VCF to 23andMe converter | `configs/converters/vcf_to_23andme/` |

Prefer these paths in new commands and documentation. Historical legacy-only configs live beside their deprecated code under `legacy/`.

Each canonical directory should contain a `default.yaml` with a portable, commented example. Use these defaults as the starting point for new runs, then copy or override values for named experiments.
