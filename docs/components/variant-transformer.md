# Variant Transformer

The variant transformer works on sparse variant tokens rather than dense aligned tensors.

## CLI

```bash
genomics variant materialize --dataset-id 1kg_high_coverage --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
genomics variant evaluate configs/predictors/variant_transformer/repo_layout.example.yaml --checkpoint best_accuracy --split test
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/predictors/variant_transformer/` |
| Configs | `configs/predictors/variant_transformer/` |

## Subcommands

| Command | Purpose |
|---|---|
| `materialize` | Build a processed sparse-variant dataset |
| `train` | Train the transformer model |
| `evaluate` | Evaluate a checkpoint |
| `analyze-counts` | Analyze variant token counts and central windows |

## When To Use

Use this pipeline when sparse variant representation is preferable to dense sequence/tensor inputs, especially for experiments focused on tokenized variants, indel sizes, and central-window sequence neighborhoods.

## Implementation Structure

| Module | Responsibility |
|---|---|
| `materialize_dataset.py` | Converts source VCF/dataset inputs into processed sparse variant records |
| `variant_schema.py` | Defines record fields and schema expectations for processed variants |
| `allele_codec.py` | Encodes REF/ALT allele information into model-friendly IDs/features |
| `dataset.py` | Loads processed records and builds PyTorch datasets/batches |
| `model.py` | Transformer model implementation |
| `rope.py` | Rotary positional encoding support |
| `training.py` | Training loop integration |
| `evaluation.py`, `evaluate_checkpoint.py` | Checkpoint loading and metrics generation |
| `config.py` | YAML config parsing and validation |
| `analyze_variant_counts.py` | Dataset statistics and central-window analysis |

## Materialization Flow

The transformer separates materialization from training. Materialization reads source sample metadata, regions, VCFs, targets, and split settings, then writes a processed dataset. Training then reads that processed dataset without repeatedly scanning large VCF files.

```text
registered dataset or explicit VCF paths
  -> region/gene filtering
  -> allele and indel encoding
  -> sample-level target assignment
  -> train/val/test split
  -> processed sparse records
  -> transformer training
```

## Model Inputs

The model works with sparse variant tokens instead of dense nucleotide windows. Tokens encode chromosome/position context, allele information, indel size constraints, and target labels. Configs control limits such as `l_max`, `max_indel_size`, central window size, class set, and split behavior.

## Output Artifacts

Materialized datasets should include enough metadata to trace source data, regions, class mappings, split parameters, and sample membership. Training outputs follow the shared run layout with manifests, copied configs, checkpoints, and evaluation JSONs.

## Registered Derived Datasets

The data registry knows these variant transformer outputs:

| Dataset ID | Default Path | Config |
|---|---|---|
| `variant_transformer_superpopulation` | `/dados/GENOMICS_DATA/variant_transformer/superpopulation` | `configs/predictors/variant_transformer/superpopulation.yaml` |
| `variant_transformer_superpopulation_32k` | `/dados/GENOMICS_DATA/variant_transformer/superpopulation_32k` | use `--central-window-size 32768` |
| `variant_transformer_pigmentation_binary` | `/dados/GENOMICS_DATA/variant_transformer/pigmentation_binary` | `configs/predictors/variant_transformer/pigmentation_binary.yaml` |

Materialize them before training:

```bash
genomics variant materialize \
  --dataset-id 1kg_high_coverage \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation \
  --target superpopulation \
  --classes AFR AMR EAS EUR SAS

genomics variant materialize \
  --dataset-id 1kg_high_coverage \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation_32k \
  --target superpopulation \
  --classes AFR AMR EAS EUR SAS \
  --central-window-size 32768

genomics variant materialize \
  --dataset-id 1kg_high_coverage \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/pigmentation_binary \
  --target pigmentation_binary \
  --classes non_pigmentation pigmentation
```

Check materialization status with:

```bash
genomics audit-data --dataset-id variant_transformer_superpopulation --fail-on-missing
genomics audit-data --dataset-id variant_transformer_superpopulation_32k --fail-on-missing
genomics audit-data --dataset-id variant_transformer_pigmentation_binary --fail-on-missing
```
