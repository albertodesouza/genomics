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

## End-To-End Inputs

The variant transformer does not read AlphaGenome dense tensors. Its source data is a combination of sample metadata, genomic regions, and phased VCF genotypes.

| Input | Source | Used For |
|---|---|---|
| Canonical dataset directory | `--dataset-id`, `--dataset-dir`, or config `dataset.source_dataset_id` / `source_dataset_dir` | Discovers samples, metadata, window regions, and VCF source hints |
| `dataset_metadata.json` | Canonical dataset root | Sample/pedigree metadata, target values, window catalog, and raw VCF source configuration |
| `references/windows/<gene>/window_metadata.json` | Canonical dataset tree | Region coordinates when regions are loaded from the dataset |
| `--regions-bed` | Optional BED with `chrom start end gene_id` | Explicit region set when not using the dataset window catalog |
| `--samples-metadata` | Optional TSV/JSON | Explicit sample metadata and target labels when not using canonical dataset metadata |
| Phased VCF files | Dataset VCF sources, `--vcf-pattern`, `--vcf-root-dir`, `KG1000_VCF_PATTERN`, or `KG1000_VCF_ROOT_DIR` | Raw REF/ALT/genotype rows used to create sparse variant tokens |
| `--target` / config `dataset.target` | CLI/config | Metadata field to predict, for example `superpopulation` or `pigmentation_binary` |
| `--classes` / config `dataset.classes` | CLI/config | Ordered list of labels retained in the dataset and mapped to target indices |

VCFs must contain genotype (`GT`) fields for the requested samples. The standard path assumes phased genotypes such as `0|1`, because haplotype identity is part of each token. With `unphased_policy: skip`, unphased calls such as `0/1` are ignored. With `unphased_policy: error`, materialization fails on the first unphased call.

## Transformations During Materialization

Materialization converts chromosome-level VCF rows into one `.pt` token record per sample:

```text
source regions + samples + phased VCFs
  -> optional gene filtering
  -> optional central-window cropping
  -> batched VCF query per region and sample batch
  -> phased genotype parsing
  -> one token per non-reference haplotype allele
  -> allele/type/length/position encoding
  -> per-sample token sorting
  -> target index assignment
  -> train/val/test split metadata
  -> processed_dir/{samples/*.pt, metadata.json, sample_index.json, splits.json, gene_vocab.json}
```

### Region Selection

Regions come from one of two places:

| Mode | Behavior |
|---|---|
| Dataset windows | `load_regions_from_dataset` reads the canonical window catalog and creates a region for each selected gene/window |
| BED file | `load_regions` reads `chrom`, 0-based BED `start`, half-open `end`, and `gene_id`, then converts the start to 1-based coordinates for VCF querying |

`--genes` filters the region list before extraction. `--central-window-size` crops each region around its center before VCF extraction. For example, a 32 KiB central window uses only the middle 32768 bp of each configured gene/window, reducing token counts and memory use.

### Sample And Target Selection

Samples are loaded from canonical metadata or an explicit metadata file. The materializer keeps only samples whose target value is included in `classes`. The class order is fixed by `classes`, so target indices are stable:

```text
classes: [AFR, AMR, EAS, EUR, SAS]
AFR -> 0
AMR -> 1
EAS -> 2
EUR -> 3
SAS -> 4
```

Splits are created after token extraction and written to `splits.json`. With `family_split_mode: family_aware`, related samples stay in the same split.

### VCF Extraction And Genotype Parsing

For each region and sample batch, materialization queries these VCF fields:

```text
CHROM, POS, REF, ALT, GT for each requested sample
```

When `bcftools` is available, the query uses `bcftools query -s <sample-list> -r <region>`. If `bcftools` is not on `PATH`, the code falls back to a slower gzip reader. Multi-allelic ALT fields are split by comma. Symbolic ALT values such as `<DEL>` are ignored.

For each sample genotype, the parser keeps only phased calls. It then emits one token for each haplotype whose allele is non-reference and resolvable in the ALT list:

| Genotype | Tokens Emitted |
|---|---|
| `0|0` | none |
| `0|1` | one token for `H2` and ALT allele 1 |
| `1|0` | one token for `H1` and ALT allele 1 |
| `1|2` | one token for `H1` ALT 1 and one token for `H2` ALT 2 |
| `./.` or unphased with `skip` | none |

### Token Fields

Each emitted token starts as a `VariantToken` with biological fields:

| Field | Meaning |
|---|---|
| `chrom` | Chromosome from the VCF row |
| `position` | 1-based genomic position from the VCF row |
| `position_relative` | `position - region.start`; used by RoPE during attention |
| `gene_id` | Region/gene identifier |
| `haplotype` | `H1` or `H2`, based on the phased GT side |
| `reference_allele` | VCF REF allele |
| `alternate_allele` | Selected ALT allele for that haplotype |
| `variant_type` | `SNP`, `INS`, or `DEL` from REF/ALT length comparison |
| `length` | `len(ALT) - len(REF)`, preserving insertion/deletion direction |

The tokens are sorted by chromosome, position, gene, and haplotype before saving, so sequence order is deterministic.

### Tensor Encoding

The saved `.pt` file for each sample contains tensors derived from those token fields:

| Tensor Key | Shape Before Batching | Encoding |
|---|---|---|
| `variant_type` | `(tokens,)` | `SNP -> 0`, `INS -> 1`, `DEL -> 2` |
| `haplotype` | `(tokens,)` | `H1 -> 0`, `H2 -> 1` |
| `gene` | `(tokens,)` | Integer ID from `gene_vocab.json` |
| `length_norm` | `(tokens, 1)` | `length / max_indel_size`, clipped to `[-1, 1]` |
| `position_relative` | `(tokens,)` | Offset within the selected region |
| `position` | `(tokens,)` | Absolute genomic position retained for analysis/debugging |
| `ref_allele` | `(tokens, l_max)` | Base IDs for REF, truncated/padded to `l_max` |
| `alt_allele` | `(tokens, l_max)` | Base IDs for ALT, truncated/padded to `l_max` |
| `target` | scalar | Class index |
| `num_tokens` | scalar/int metadata | Number of emitted variant tokens |

Base IDs are `A -> 0`, `C -> 1`, `G -> 2`, `T -> 3`, `N -> 4`, and `PAD -> 5`. Alleles longer than `l_max` are truncated; shorter alleles are padded with `PAD`.

### Materialized Dataset Files

A processed dataset has this structure:

```text
<processed_dir>/
  metadata.json
  gene_vocab.json
  sample_index.json
  splits.json
  samples/
    <sample_id>.pt
```

`metadata.json` records the target, class order, region list, VCF source mapping, `l_max`, `max_indel_size`, optional `central_window_size`, split sizes, and split strategy details. `sample_index.json` records each sample path, target, token count, family, population, and superpopulation. `splits.json` stores the sample IDs in each split.

## Model Inputs

The model works with sparse variant tokens instead of dense nucleotide windows. Tokens encode chromosome/position context, allele information, indel size constraints, haplotype identity, gene identity, and target labels. Configs control limits such as `l_max`, `max_indel_size`, central window size, class set, split behavior, and maximum sequence length.

At training time, `VariantTokenDataset` reads one `.pt` file per sample. `collate_variant_tokens` pads variable-length samples in a batch and adds a classifier token:

| Batch Field | Shape | Meaning |
|---|---|---|
| `variant_type` | `(batch, max_tokens)` | Padded token type IDs |
| `haplotype` | `(batch, max_tokens)` | Padded haplotype IDs |
| `gene` | `(batch, max_tokens)` | Padded gene IDs |
| `length_norm` | `(batch, max_tokens, 1)` | Padded normalized indel lengths |
| `position_relative` | `(batch, max_tokens)` | Token positions used for rotary attention |
| `position` | `(batch, max_tokens)` | Absolute positions retained in the batch |
| `ref_allele` | `(batch, max_tokens, l_max)` | Padded REF base IDs |
| `alt_allele` | `(batch, max_tokens, l_max)` | Padded ALT base IDs |
| `attention_mask` | `(batch, max_tokens + 1)` | Marks the prepended classifier token and real variant tokens as valid |
| `targets` | `(batch,)` | Target class indices |
| `lengths` | `(batch,)` | Original token count per sample |
| `sample_ids` | Python list | Sample IDs for diagnostics |

If `dataset.max_sequence_length` is set, samples longer than that limit are handled by `dataset.truncate_policy`:

| Policy | Behavior |
|---|---|
| `error` | Fail when a sample has more tokens than allowed |
| `keep_first` | Keep the earliest tokens after sorting |
| `keep_last` | Keep the latest tokens after sorting |

`dataset.loading_strategy` controls IO behavior. `lazy` loads each sample file on demand. `preload` loads all selected split samples into memory when the dataset object is created.

## Transformer Architecture

The available active model for this pipeline is `VariantTransformerClassifier` in `model.py`. It is configured by the `model` block in the YAML.

```text
padded sparse token tensors
  -> token-type, haplotype, gene, length, and allele encoders
  -> projection to d_model
  -> prepend learned CLS token
  -> RoPE multi-head self-attention blocks
  -> CLS representation
  -> classifier MLP
  -> class logits
```

### Token Encoder

Each token is embedded from several sources and then projected to `d_model`:

| Component | Config Dimension | Transformation |
|---|---|---|
| Variant type | `d_type` | Embedding over `SNP`, `INS`, `DEL` |
| Haplotype | `d_hap` | Embedding over `H1`, `H2` |
| Gene | `d_gene` | Embedding over `gene_vocab.json` IDs |
| Normalized length | `d_len` | Linear projection from one scalar |
| REF/ALT alleles | `d_base`, `d_allele` | Base embeddings for fixed-length REF/ALT strings, flattened and passed through an MLP |
| Combined token | `d_model` | Linear projection of concatenated component embeddings |

### Attention Blocks

The transformer uses `layers` repeated blocks. Each block applies layer normalization, multi-head self-attention, residual connections, and an MLP with hidden size `d_model * mlp_ratio`. Rotary positional encoding is applied to queries and keys using `position_relative`, so attention receives information about each token's position within its source region.

Important constraints and settings:

| Field | Meaning |
|---|---|
| `d_model` | Width of token embeddings and transformer hidden states |
| `heads` | Number of attention heads; `d_model` must be divisible by `heads` |
| `d_model / heads` | Per-head dimension; must be even for RoPE |
| `layers` | Number of transformer blocks |
| `mlp_ratio` | Expansion factor in each block MLP |
| `dropout` | Dropout in allele encoder, attention, block MLP, and classifier head |
| `rope_base` | Base frequency parameter for rotary positional encoding |

The classifier reads only the final CLS token representation, applies `LayerNorm`, then a small MLP to produce class logits. Training uses the configured class order from the processed dataset metadata.

## Available Models

The variant pipeline currently exposes one active architecture:

| Model | File | Purpose |
|---|---|---|
| `VariantTransformerClassifier` | `src/genomics/predictors/variant_transformer/model.py` | Sparse-token transformer classifier for phased VCF-derived variant sequences |

Model scale is controlled through YAML rather than through multiple model names. For smoke tests or constrained hardware, reduce `d_model`, `layers`, `heads`, `batch_size`, and/or `dataset.max_sequence_length`. For larger experiments, increase these values while keeping the RoPE constraints above valid.

## Output Artifacts

Materialized datasets should include enough metadata to trace source data, regions, class mappings, split parameters, and sample membership. Training outputs follow the shared run layout with manifests, copied configs, checkpoints, and evaluation JSONs.

Typical training outputs are written under:

```text
<dataset.results_dir>/<experiment_name>/
```

The experiment name includes the target and main architecture scale, for example:

```text
variant_transformer_superpopulation_d256_l6_h8_<hash>/
```

The run directory contains copied/resolved config files, `manifest.json`, checkpoints such as `models/best_accuracy.pt` and `models/best_loss.pt`, and evaluation JSON/plot artifacts produced by the shared metrics utilities.

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
