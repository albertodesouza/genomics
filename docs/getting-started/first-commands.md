# First Commands

## Inspect The CLI

```bash
genomics --help
genomics genotype --help
genomics variant --help
```

## Audit Configuration And Data

```bash
genomics audit-configs
genomics audit-configs --fail-on-active-legacy
genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing
```

## Run The Operational Pipeline

```bash
genomics genomes-analyzer run --config configs/genomes_analyzer/config_human_30x_low_memory.yaml
```

## Train A Genotype Model

```bash
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
```

## Train A Variant Transformer

```bash
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
```

## Run SNP Ancestry

```bash
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml
```
