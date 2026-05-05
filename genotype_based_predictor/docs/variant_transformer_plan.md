# Variant Transformer Plan

## Goal

Implement the sparse variant-token proposal as a new module, without mixing it into the dense CNN pipeline already living in `genotype_based_predictor/`.

The new approach should:

- use the same canonical dataset in `/dados/GENOMICS_DATA/v1/1kG_high_coverage`
- preserve logical views and family-aware split support
- tokenize SNPs/INDELs per haplotype
- enrich tokens with AlphaGenome-derived functional context
- train a transformer on variable-length token sequences

The implementation should avoid large-scale duplication by extracting only the infrastructure that is truly representation-agnostic.

## Recommended Package Layout

Create two siblings:

- `genomics_common/`
- `variant_transformer_predictor/`

Keep `genotype_based_predictor/` focused on the current dense aligned-tensor pipeline.

### Proposed Tree

```text
genomics_common/
  __init__.py
  dataset_layout.py
  genomic_dataset.py
  data_splitting.py
  runtime.py
  view_resolution.py

variant_transformer_predictor/
  __init__.py
  config.py
  dataset.py
  tokenizer.py
  variant_context.py
  collate.py
  data_pipeline.py
  train.py
  training.py
  experiment.py
  models/
    __init__.py
    transformer_model.py
  configs/
  docs/
```

## What Should Be Shared

Only move code that is independent of the input representation.

### Move Entirely

These files already look like shared infrastructure and should be moved out of `genotype_based_predictor/`.

1. `genotype_based_predictor/genomic_dataset.py`
2. `genotype_based_predictor/dataset_layout.py`
3. `genotype_based_predictor/data_splitting.py`

New locations:

- `genomics_common/genomic_dataset.py`
- `genomics_common/dataset_layout.py`
- `genomics_common/data_splitting.py`

### Extract Partially

These files contain a mix of generic and dense-pipeline-specific code.

#### From `genotype_based_predictor/utils.py`

Move only:

- `set_random_seeds(...)`
- `worker_init_fn(...)`

Keep in `genotype_based_predictor/utils.py`:

- `block_reduce_2d(...)`
- `taint_sample(...)`
- any dense/debug-specific helpers

Suggested destination:

- `genomics_common/runtime.py`

#### From `genotype_based_predictor/data_pipeline.py`

Extract only the view/dataset resolution helpers:

- `_build_view_definition(...)`
- `_build_resolved_view_definition(...)`
- `_resolve_runtime_dataset_dir(...)`
- `_load_base_dataset(...)`

Suggested destination:

- `genomics_common/view_resolution.py`

The cache-building logic should stay local to each pipeline, because dense tensor caches and sparse token caches will diverge heavily.

#### From `genotype_based_predictor/experiment.py`

Potentially extract a tiny generic core later if needed:

- interrupt flag handling
- experiment directory setup pattern

For the first iteration, duplication here is acceptable. The current file is small and still imports dense config naming.

## What Should Stay Dense-Specific

Do not move these. They are tightly coupled to the aligned dense representation.

1. `genotype_based_predictor/dataset.py`
2. `genotype_based_predictor/dynamic_indel_alignment.py`
3. `genotype_based_predictor/indel_tensor_builder.py`
4. `genotype_based_predictor/normalization.py`
5. `genotype_based_predictor/models/nn_model.py`
6. `genotype_based_predictor/models/cnn_model.py`
7. `genotype_based_predictor/models/cnn2_model.py`
8. `genotype_based_predictor/interpretability.py`
9. `genotype_based_predictor/synthetic_indel_demo.py`

## Why a New Module Is Better

The transformer proposal is not just another INDEL handling mode. It changes:

- the sample representation
- the dataset contract
- the collate strategy
- the model family
- the cache format

Keeping this inside `genotype_based_predictor/` would force branching logic across:

- `config.py`
- `dataset.py`
- `data_pipeline.py`
- `train.py`
- `models/`

That would make both paradigms harder to evolve and debug.

## Implementation Phases

## Phase 1. Extract Shared Infrastructure

Goal: create `genomics_common/` without changing behavior.

### Tasks

1. Create `genomics_common/` package.
2. Move:
   - `genomic_dataset.py`
   - `dataset_layout.py`
   - `data_splitting.py`
3. Extract from `utils.py`:
   - `set_random_seeds(...)`
   - `worker_init_fn(...)`
4. Create `view_resolution.py` for shared view loading/resolution.
5. Update imports in `genotype_based_predictor/`.

### Success Criteria

- `genotype_based_predictor/` still runs unchanged
- shared code no longer lives in a dense-specific namespace

## Phase 2. Create the New Module Skeleton

Goal: create a clean sparse pipeline shell.

### Files to Create

- `variant_transformer_predictor/config.py`
- `variant_transformer_predictor/dataset.py`
- `variant_transformer_predictor/tokenizer.py`
- `variant_transformer_predictor/variant_context.py`
- `variant_transformer_predictor/collate.py`
- `variant_transformer_predictor/data_pipeline.py`
- `variant_transformer_predictor/models/transformer_model.py`
- `variant_transformer_predictor/train.py`

### Shared Imports

The new module should directly reuse:

- `GenomicDataset`
- dataset layout helpers
- family-aware split helpers
- seed/runtime helpers
- view resolution helpers

## Phase 3. Define the Sparse Config Contract

Goal: define a config that matches the token pipeline instead of the dense tensor pipeline.

### New Config Sections

`dataset_input`

- `view_path`
- `processed_cache_dir`
- `cache_tokenized_samples`
- `alphagenome_outputs`
- `selected_track_index`
- `max_tokens_per_sample`
- `token_context_mode`
- `allele_pad_length`

`model`

- transformer depth
- `d_model`
- `nhead`
- `ffn_dim`
- dropout
- pooling strategy (`cls` by default)

`data_split`

- preserve `family_split_mode`
- preserve `train/val/test` fractions

`training`

- standard optimizer/scheduler fields

### Deliberate Non-Goals

Do not carry over dense-only config fields such as:

- `tensor_layout`
- `indel_include_valid_mask`
- `normalization_method` tied to dense tracks

## Phase 4. Implement Variant Tokenization

Goal: map VCF variants to canonical per-haplotype tokens.

### New Component

- `variant_transformer_predictor/tokenizer.py`

### Responsibilities

1. Query the VCF for variants intersecting each gene window.
2. Separate events by haplotype using phased genotype.
3. Normalize events into a canonical token representation.
4. Handle SNP/INS/DEL first.
5. Count and optionally skip unsupported complex cases in a controlled way.

### Token Schema

Start with a plain Python dict contract:

```python
{
    "sample_id": "HG00096",
    "gene": "OCA2",
    "haplotype": "H1",
    "chrom": "chr15",
    "pos_ref": 12345678,
    "pos_rel": 12483,
    "variant_type": "SNP",
    "ref": "A",
    "alt": "G",
    "event_len": 0,
}
```

### Scope Control

For v1:

- accept phased, simple SNP/INS/DEL events
- explicitly map `<DEL>` if needed
- log unsupported events instead of over-designing full VCF normalization

## Phase 5. Add AlphaGenome Functional Context

Goal: enrich tokens with local functional information.

### New Component

- `variant_transformer_predictor/variant_context.py`

### Recommended First Implementation

Start with the cheapest version:

- point value at the variant position

Optional next step:

- individual minus reference point value

Do not start with local windows around each variant. That adds a large amount of complexity and cost.

### Output Shape

Each token gets a fixed-size numeric vector:

```python
token["functional_context"] = [x1, x2, ..., xk]
```

where `k` depends on the number of selected tracks/features.

## Phase 6. Build the Dataset Contract

Goal: make `Dataset.__getitem__` return token fields, not dense tensors.

### `variant_transformer_predictor/dataset.py`

For each sample:

1. load metadata and windows from the canonical dataset
2. tokenize variants across selected genes and haplotypes
3. attach functional context
4. encode the target label
5. return a structured sample payload

### Recommended Return Format

Return separate fields instead of prematurely concatenating everything into one giant vector:

```python
{
    "variant_type_ids": ...,
    "hap_ids": ...,
    "gene_ids": ...,
    "positions": ...,
    "ref_base_ids": ...,
    "alt_base_ids": ...,
    "event_len": ...,
    "functional_context": ...,
    "target": ...,
}
```

This keeps the model cleaner and makes ablations easier.

## Phase 7. Implement Variable-Length Collation

Goal: batch samples with different token counts.

### `variant_transformer_predictor/collate.py`

Responsibilities:

1. pad each field to the max token length in the batch
2. create `attention_mask`
3. optionally truncate to `max_tokens_per_sample`

### Batch Contract

```python
{
    "variant_type_ids": Tensor[B, T],
    "hap_ids": Tensor[B, T],
    "gene_ids": Tensor[B, T],
    "positions": Tensor[B, T],
    "ref_base_ids": Tensor[B, T, A],
    "alt_base_ids": Tensor[B, T, A],
    "event_len": Tensor[B, T],
    "functional_context": Tensor[B, T, C],
    "attention_mask": Tensor[B, T],
    "targets": Tensor[B],
}
```

## Phase 8. Implement the Transformer Model

Goal: create a first baseline model with minimal novelty.

### `variant_transformer_predictor/models/transformer_model.py`

Recommended structure:

1. embedding tables for:
   - variant type
   - haplotype
   - gene
2. a positional encoding path for relative position
3. a tiny allele encoder for `ref` and `alt`
4. a small projection for event length
5. a linear projection for functional context
6. fuse all token parts into `d_model`
7. prepend `[CLS]`
8. run `TransformerEncoder`
9. classify from `[CLS]`

### Keep v1 Simple

Do not start with:

- cross-attention between haplotypes
- hybrid CNN-transformer designs
- sparse attention kernels

Those are follow-up experiments.

## Phase 9. Add Training and Evaluation

Goal: make the sparse pipeline runnable end-to-end.

### Reuse Strategy

Reuse concepts, not dense-specific modules.

Good candidates for later extraction:

- generic optimizer setup
- generic scheduler setup
- generic checkpoint helpers

Avoid over-sharing training code until the sparse batch contract stabilizes.

## Phase 10. Minimum Viable Experiment

Goal: produce a first trustworthy baseline before optimizing.

### Recommended MVP Scope

1. same 11 genes
2. same `superpopulation` target
3. family-aware split enabled
4. simple SNP/INS/DEL support
5. point functional context only
6. moderate transformer encoder
7. `max_tokens_per_sample` safeguard

### Success Criteria

- token pipeline runs end-to-end
- batches are valid
- model trains without shape/path issues
- results are comparable against the dense CNN baseline

## Code Movement Checklist

### Files to Move

- `genotype_based_predictor/genomic_dataset.py`
- `genotype_based_predictor/dataset_layout.py`
- `genotype_based_predictor/data_splitting.py`

### Functions to Extract

From `genotype_based_predictor/utils.py`:

- `set_random_seeds`
- `worker_init_fn`

From `genotype_based_predictor/data_pipeline.py`:

- `_build_view_definition`
- `_build_resolved_view_definition`
- `_resolve_runtime_dataset_dir`
- `_load_base_dataset`

### Files to Leave Untouched in Dense Pipeline

- `genotype_based_predictor/dataset.py`
- `genotype_based_predictor/dynamic_indel_alignment.py`
- `genotype_based_predictor/indel_tensor_builder.py`
- `genotype_based_predictor/normalization.py`
- `genotype_based_predictor/models/*`

## Risks and Scope Control

### Main Risk

Trying to implement all of this at once will explode scope.

The riskiest early temptations are:

- full support for every VCF edge case
- sophisticated local AlphaGenome context windows
- generic shared cache abstractions
- fancy transformer variants before a baseline exists

### Scope Discipline

Build the first version in this order:

1. shared infrastructure extraction
2. sparse tokenization without advanced context
3. collate + baseline transformer
4. only then richer functional context
5. only then caching and performance work

## Immediate Next Step

The first concrete implementation step should be:

1. create `genomics_common/`
2. move the three representation-agnostic files there
3. extract runtime/view helpers
4. adjust imports in `genotype_based_predictor/`
5. verify the dense pipeline still runs

Only after that should `variant_transformer_predictor/` start being created.
