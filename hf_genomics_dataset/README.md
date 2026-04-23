## hf_genomics_dataset

Phase 1 of the Hugging Face dataset migration lives here.

This module converts the existing filesystem dataset layout into a normalized,
local Hugging Face dataset store with two tables:

- `samples`: one row per individual
- `gene_windows`: one row per `(sample_id, gene)`

The output is intended to be saved locally with `datasets.Dataset.save_to_disk()`.

Initial scope:

- read existing dataset folders produced by `build_non_longevous_dataset`
- merge them into one canonical local store
- preserve provenance and source manifests
- keep the conversion code organized and explicit

This phase does not yet modify training code.

The current CLI is a conversion/import tool. It reads already-built filesystem
datasets and writes a normalized local Hugging Face dataset store. A native
build-from-scratch pipeline will be added separately.

### Convert

Install the module requirements:

```bash
pip install -r hf_genomics_dataset/requirements.txt
```

Convert one or more existing filesystem datasets into a canonical local HF store:

```bash
python3 hf_genomics_dataset/convert_dataset.py \
  --source /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000 \
  --source /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_random \
  --source /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_random_11_1 \
  --source /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_random_11_2 \
  --output /dados/GENOMICS_DATA/hf/genomics_master
```

Lower-memory conversion is supported with temporary shard writes:

```bash
python3 hf_genomics_dataset/convert_dataset.py \
  --source /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000 \
  --output /dados/GENOMICS_DATA/hf/genomics_master \
  --chunk-size 256
```

`--chunk-size` now streams rows sample-by-sample and flushes temporary shards as
the buffers fill, instead of accumulating a full source dataset in memory first.
It still does not make conversion fully constant-memory because duplicate tracking
and final shard merging keep some state, but it substantially reduces peak row
buffering.

By default, `gene_windows` rows carry the reference and haplotype sequences when
those files are present. If you want a smaller store, you can disable that:

```bash
python3 hf_genomics_dataset/convert_dataset.py \
  --source /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000 \
  --output /dados/GENOMICS_DATA/hf/genomics_master \
  --no-sequences
```

### Inspect

Validate and inspect an already-built store:

```bash
python3 hf_genomics_dataset/inspect_dataset.py \
  --dataset /dados/GENOMICS_DATA/hf/genomics_master
```

JSON output is also supported:

```bash
python3 hf_genomics_dataset/inspect_dataset.py \
  --dataset /dados/GENOMICS_DATA/hf/genomics_master \
  --json
```

### Duplicate handling

The builder treats duplicates conservatively:

- identical duplicate `sample_id` rows are skipped
- identical duplicate `(sample_id, gene)` rows are skipped
- conflicting duplicates fail the build

This keeps the canonical store deterministic and forces source conflicts to be resolved explicitly.

### Provenance

Each row carries source-level provenance fields in addition to `source_dataset`.

For `samples`:

- `source_path`
- `source_sample_path`

For `gene_windows`:

- `source_path`
- `source_sample_path`
- `source_window_path`
- `source_metadata_path`
