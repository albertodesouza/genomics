# CLI Reference

The primary entrypoint is:

```bash
genomics
```

## Top-Level Commands

| Command | Purpose |
|---|---|
| `audit-configs` | Check configs for legacy paths and active/inactive status |
| `audit-data` | Validate registered dataset paths and expected artifacts |
| `config ...` | Describe, validate, and export typed config schemas |
| `completion bash` | Print Bash completion script |
| `convert vcf-to-23andme` | Convert VCF to 23andMe raw format |
| `snp-ancestry run` | Run SNP ancestry pipeline |
| `genomes-analyzer run` | Run FASTQ/BAM/CRAM/VCF operational workflow |
| `dataset-builders non-longevous ...` | Build derived 1000G/AlphaGenome datasets |
| `alphagenome ...` | Run AlphaGenome analysis/integration utilities |
| `genotype ...` | Dense/aligned genotype predictor workflows |
| `variant ...` | Sparse variant transformer workflows |
| `neural ...` | Deprecated neural ancestry reproducibility commands |

## Examples

```bash
genomics audit-configs --fail-on-active-legacy
genomics audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 3 --fail-on-missing
genomics config describe genotype
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics completion bash
```

```bash
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype split configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype test configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype search configs/predictors/genotype_based/icann/search_rf_xgboost.yaml
genomics genotype confidence-intervals configs/predictors/genotype_based/icann/search_rf_xgboost.yaml --experiment-dir results/genotype_based_predictor/icann/search/rf_xgboost_pca300/best --split test
genomics genotype evaluate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy --split test
```

`genomics genotype train` trains on the training split, validates during training, and reports final validation metrics for the `best_accuracy` checkpoint. It does not evaluate the test split. Use `genomics genotype test` after model/hyperparameter selection to evaluate `best_accuracy` on the held-out test split.

`genomics genotype split` materializes or validates the processed dataset cache, split metadata, and dataset report plots without training a model.

`genomics genotype search` runs validation-only hyperparameter search for configured sklearn baselines. Use `genomics genotype test` on the selected best directory after model selection.

`genomics genotype confidence-intervals` recomputes metrics and configured bootstrap confidence intervals for a saved model artifact/checkpoint without retraining.

```bash
genomics variant materialize --dataset-id 1kg_high_coverage --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
```

Before training variant transformer configs, verify that their processed datasets exist:

```bash
genomics audit-data --dataset-id variant_transformer_superpopulation --fail-on-missing
```

Use `--help` at any level:

```bash
genomics genotype train --help
genomics alphagenome analyze -- --help
```
