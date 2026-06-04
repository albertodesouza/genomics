# Components Overview

The repository is organized into active components, external tools, and legacy components.

| Component | CLI Prefix | Code |
|---|---|---|
| Genomes Analyzer | `genomics genomes-analyzer` | `src/genomics/workflows/genomes_analyzer/` |
| Genotype predictor | `genomics genotype` | `src/genomics/predictors/genotype_based/` |
| Variant transformer | `genomics variant` | `src/genomics/predictors/variant_transformer/` |
| SNP ancestry | `genomics snp-ancestry` | `src/genomics/predictors/snp_ancestry/` |
| VCF to 23andMe | `genomics convert vcf-to-23andme` | `src/genomics/converters/vcf_to_23andme/` |
| AlphaGenome | `genomics alphagenome` | `src/genomics/workflows/alphagenome/` |
| Dataset builders | `genomics dataset-builders` | `src/genomics/workflows/dataset_builders/` |
| Legacy neural ancestry | `genomics neural` | `legacy/neural_ancestry_predictor_deprecated/` |

All active Python code should import through `genomics.*` or relative imports inside the same package.
