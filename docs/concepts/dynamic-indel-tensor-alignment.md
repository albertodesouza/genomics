# Dynamic Indel Tensor Alignment

Dynamic Indel Tensor Alignment, or DITA, is the aligned tensor representation used by genotype configs with `tensor_layout: haplotype_channels` and `alignment_mapping: bcftools_chain`.

The goal is to compare AlphaGenome-derived signal tracks across samples even when samples contain insertions, deletions, and SNPs inside the same reference window.

## Why Alignment Is Needed

AlphaGenome predictions are dense arrays indexed along a sequence. If one sample haplotype has an insertion or deletion relative to the reference, the same array index no longer corresponds to the same biological position across samples.

The legacy `raw_center_crop` layout avoids this problem by using a fixed center crop and treating the result as a baseline representation. That is useful for historical comparisons, but it does not explicitly repair coordinate shifts caused by indels.

DITA builds a shared coordinate axis for each window and maps sample haplotype predictions onto that axis. This lets the model see both the AlphaGenome signal and where sequence variation changed the coordinate structure.

## Inputs

For each sample/window, DITA reads the canonical dataset files plus BCFtools consensus artifacts:

```text
individuals/<sample>/windows/<target>/
  predictions_H1/
  predictions_H2/
  <sample>.H1.window.raw.fa
  <sample>.H2.window.raw.fa
  <sample>.window.consensus_ready.vcf.gz
  <sample>.window.consensus_ready.vcf.gz.tbi
  <sample>.window.vcf.gz
  <sample>.window.vcf.gz.tbi
references/windows/<target>/
  ref.window.fa
  window_metadata.json
```

Validate availability with:

```bash
genomics audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --fail-on-missing
```

## Coordinate Model

Each target window has a reference coordinate interval from `window_metadata.json`. For each sample haplotype, the pipeline uses the consensus/chain mapping to connect three coordinate systems:

| Coordinate System | Meaning |
|---|---|
| Reference window | The canonical genomic interval for the target. |
| Sample haplotype sequence | The per-sample sequence after applying variants. |
| Expanded aligned axis | A shared axis that can represent reference positions plus inserted sequence slots. |

The expanded aligned axis preserves comparable reference positions while adding room for insertions. Deletions are represented as positions where a sample lacks sequence or valid signal.

## Feature Rows

DITA tensors combine signal rows and mask rows. The exact row count depends on config settings such as `genes_to_use`, `haplotype_mode`, `alphagenome_outputs`, `ontology_terms`, and `feature_mode`.

| Row Type | Meaning |
|---|---|
| AlphaGenome signal | Selected prediction tracks mapped onto the aligned axis. |
| Insertion mask (`MI`) | Marks aligned positions introduced by sample insertions. |
| Deletion mask (`MD`) | Marks reference positions deleted in the sample haplotype. |
| Validity mask (`MV`) | Marks positions where mapped signal is valid for that sample/haplotype. |
| SNP mask (`MS`) | Optional mask for SNP positions when `indel_include_snp_mask` is enabled. |

`feature_mode` controls which rows are used:

| `feature_mode` | Rows |
|---|---|
| `signals_only` | AlphaGenome signal rows only. |
| `masks_only` | Indel/SNP/validity mask rows only. |
| `signals_and_masks` | Signals plus masks. |

## Signal Transformations

The aligned signal can be transformed before normalization:

| Setting | Meaning |
|---|---|
| `alphagenome_signal_transform: absolute` | Use stored AlphaGenome predictions directly. |
| `alphagenome_signal_transform: delta_reference` | Subtract reference-only predictions for the same track/window. |
| `alphagenome_signal_variant_mask: true` | Zero signal positions that do not overlap variation in the selected sample universe. |

`delta_reference` needs a reference-only prediction dataset. Keep configs that point to historical reference-only paths inactive until the reference-only dataset is materialized under a canonical path.

## Operational Checks

Before running a DITA genotype experiment:

1. Validate config schema with `genomics config validate <config.yaml>`.
2. Validate the source dataset with `genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing`.
3. Validate chain artifacts with `genomics audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --fail-on-missing`.
4. Run `genomics genotype split <config.yaml>` to materialize or validate the processed cache without training.

If chain artifacts are missing but available in a separate tree, sync them with `genomics genotype sync-bcftools-artifacts` before rebuilding caches.
