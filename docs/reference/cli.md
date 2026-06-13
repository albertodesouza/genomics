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
| `snp-ancestry markers` | Export ranked ancestry-informative markers from computed SNP ancestry statistics |
| `snp-ancestry prune` | Positionally prune ranked ancestry-informative markers |
| `snp-ancestry train-ml` | Train sklearn ancestry classifiers from exported AIMs |
| `snp-ancestry ablate` | Retrain sklearn baselines after removing top AIMs |
| `snp-ancestry plot` | Plot ML metrics, feature importance, and AIM-ablation curves |
| `genomes-analyzer run` | Run FASTQ/BAM/CRAM/VCF operational workflow |
| `dataset-builders non-longevous ...` | Build derived 1000G/AlphaGenome datasets |
| `alphagenome ...` | Run AlphaGenome analysis/integration utilities |
| `genotype ...` | Dense/aligned genotype predictor workflows |
| `variant ...` | Sparse variant transformer workflows |

## Examples

```bash
genomics audit-configs --fail-on-active-legacy
genomics audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 3 --fail-on-missing
genomics config describe genotype
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics completion bash
```

`genomics audit-configs` checks active genotype and variant transformer configs for legacy dataset paths and result/cache migration status. Use `--fail-on-active-legacy` as the normal gate for new runs.

`genomics audit-data` validates registered dataset IDs. Add `--check-bcftools-chain` before aligned genotype runs that use `alignment_mapping: bcftools_chain`.

## Genotype Commands

```bash
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype split configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype test configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype search configs/predictors/genotype_based/icann/search_rf_xgboost.yaml
genomics genotype search configs/predictors/genotype_based/icann/search_cnn2_ablation.yaml
genomics genotype stability configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml
genomics genotype confidence-intervals configs/predictors/genotype_based/icann/search_rf_xgboost.yaml --experiment-dir results/genotype_based_predictor/icann/search/rf_xgboost_pca300/best --split test
genomics genotype evaluate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy --split test
genomics genotype pca-variance configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml --output results/pca_variance.png --json-output results/pca_variance.json
genomics genotype workbench --host 127.0.0.1 --port 8780
genomics genotype sync-bcftools-artifacts --source-dir /path/to/consensus --target-dir /path/to/canonical --link-mode hardlink
genomics genotype single-gene-screen configs/predictors/genotype_based/neural_legacy/pigmentation_binary_single_gene_screen.yaml --dry-run
```

`genomics genotype train` trains on the training split, validates during training, and reports final validation metrics for the `best_accuracy` checkpoint. It does not evaluate the test split. Use `genomics genotype test` after model/hyperparameter selection to evaluate `best_accuracy` on the held-out test split.

`genomics genotype split` materializes or validates the processed dataset cache, split metadata, and dataset report plots without training a model.

`genomics genotype search` runs validation-only hyperparameter search for configured sklearn baselines or named PyTorch ablation candidates. Use `genomics genotype test` on the selected best directory after model selection.

`genomics genotype stability` evaluates sklearn model stability on the development split only. It keeps the original test split fixed, resamples `train+val` using `stability_analysis.strategy` (`repeated_random_split`, `randomized_split`, or `cross_validation`), and writes aggregate validation metrics. Use the held-out test split only after model selection. Set the same `stability_analysis.split_plan_path` across model configs to reuse exactly the same resampling plan by `sample_id`.

`genomics genotype confidence-intervals` recomputes metrics and configured bootstrap confidence intervals for a saved model artifact/checkpoint without retraining.

`genomics genotype pca-variance` computes and plots sklearn PCA explained variance for the selected processed dataset/config. Use `--force` to rebuild existing outputs.

`genomics genotype workbench` launches the local genotype workbench apps for inspecting datasets, aligned tensors, AlphaGenome tracks, and experiment outputs.

`genomics genotype sync-bcftools-artifacts` previews or applies hardlink/symlink/copy operations for consensus and chain artifacts required by the aligned `haplotype_channels` layout. Add `--apply` only after reviewing the preview.

`genomics genotype single-gene-screen` expands a base config into per-gene or per-ontology runs. Use `--dry-run` first to inspect generated commands and output paths.

## SNP Ancestry Commands

```bash
genomics snp-ancestry run configs/predictors/snp_ancestry/default.yaml
genomics snp-ancestry run configs/predictors/snp_ancestry/icann/gene_windows_h1_mlc.yaml
genomics snp-ancestry run configs/predictors/snp_ancestry/chr15_aims.yaml
genomics snp-ancestry markers --config configs/predictors/snp_ancestry/chr15_aims.yaml --top 500 --output results/snp_ancestry_predictor/chr15/aims_top500.tsv
genomics snp-ancestry prune --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --window-bp 50000 --output results/snp_ancestry_predictor/chr15/aims_top500_pruned_50kb.tsv
genomics snp-ancestry train-ml --config configs/predictors/snp_ancestry/chr15_aims.yaml --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --models logistic random_forest --output-dir results/snp_ancestry_predictor/chr15/ml
genomics snp-ancestry ablate --config configs/predictors/snp_ancestry/chr15_aims.yaml --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --remove-top 0 1 5 10 50 100 --output-dir results/snp_ancestry_predictor/chr15/ablation
genomics snp-ancestry plot --ml-dir results/snp_ancestry_predictor/chr15/ml --ablation-dir results/snp_ancestry_predictor/chr15/ablation --output-dir results/snp_ancestry_predictor/chr15/plots
```

`genomics snp-ancestry run` preserves the existing conversion, statistics, and prediction pipeline. The legacy `--config`/`-c` form is still accepted. `genomics snp-ancestry markers` is an optional post-processing command: it reads the statistics JSON produced by `run`, ranks markers by `fst`, `maf`, or `max_delta_frequency`, and writes an audit-friendly TSV with per-class allele frequencies. `genomics snp-ancestry prune` removes lower-ranked markers within a configured base-pair window of already kept markers, preserving the TSV columns and recalculating ranks. `genomics snp-ancestry train-ml` consumes an AIM TSV and the same per-individual 23andMe files to train sklearn `logistic` and/or `random_forest` baselines, writing metrics, predictions, feature importance, and `model.joblib` artifacts. `genomics snp-ancestry ablate` uses the same inputs, removes ranked marker prefixes such as the top 1, 5, or 10 AIMs, retrains the selected sklearn models, and writes `ablation.tsv` plus `summary.json` for measuring robustness to top-marker removal. `genomics snp-ancestry plot` turns those outputs into PNG confusion matrices, top-feature bar charts, and ablation curves.

## Variant Commands

```bash
genomics variant materialize --dataset-id 1kg_high_coverage --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
genomics variant evaluate configs/predictors/variant_transformer/repo_layout.example.yaml --checkpoint best_accuracy --split test
genomics variant analyze-counts /dados/GENOMICS_DATA/variant_transformer/superpopulation --central-window-size 32768
```

Before training variant transformer configs, verify that their processed datasets exist:

```bash
genomics audit-data --dataset-id variant_transformer_superpopulation --fail-on-missing
```

`genomics variant materialize` creates the sparse-token processed dataset from a canonical dataset, explicit VCF sources, or BED/sample metadata inputs. `genomics variant train` and `genomics variant evaluate` read that materialized dataset.

`genomics variant analyze-counts` summarizes token counts and central-window behavior for an existing processed dataset.

## AlphaGenome Commands

```bash
genomics alphagenome analyze -- -i sequence.fasta -k API_KEY -o results/
genomics alphagenome integrate -- --integrated --vcf vcf/sample.vcf.gz --ref refs/GRCh38.fa --api-key API_KEY --output integrated_analysis/
genomics alphagenome tracks --api-key API_KEY --output configs/workflows/alphagenome/tracks.json
```

Arguments after `--` are forwarded to the underlying AlphaGenome modules. `tracks` exports output and ontology metadata used when choosing `alphagenome_outputs` and `ontology_terms`.

## Dataset Builder Commands

```bash
genomics dataset-builders non-longevous build --config configs/workflows/non_longevous_dataset/default.yaml
genomics dataset-builders non-longevous build-window -- --help
genomics dataset-builders non-longevous visualize configs/workflows/non_longevous_dataset/default.yaml
```

`build-window` forwards arguments to the window builder module. Use `-- --help` to inspect the forwarded module's options.

Use `--help` at any level:

```bash
genomics genotype train --help
genomics alphagenome analyze -- --help
```
