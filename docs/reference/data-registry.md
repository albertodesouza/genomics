# Data Registry

The data registry maps stable dataset IDs to concrete filesystem paths. Active configs should prefer dataset IDs over absolute paths when the component supports them.

The implementation lives in `src/genomics/core/data_registry.py`; default roots are defined in `src/genomics/workspace.py`.

## Environment Roots

| Variable | Purpose |
|---|---|
| `GENOMICS_DATA_ROOT` | Overrides the default root used for registered datasets. |
| `GENOMICS_RESULTS_ROOT` | Overrides the default root for run outputs and caches. |

If `GENOMICS_DATA_ROOT` is not set, dataset paths are resolved relative to the workspace default data root.

## Common Dataset IDs

| Dataset ID | Role |
|---|---|
| `1kg_high_coverage` | Canonical high-coverage 1000 Genomes dataset used by active genotype and variant workflows. |
| `legacy_top3_1kg_high_coverage` | Historical `top3` dataset retained for provenance and reproducibility only. |
| `variant_transformer_superpopulation` | Materialized sparse-token variant transformer dataset for superpopulation labels. |
| `variant_transformer_superpopulation_32k` | Materialized variant transformer dataset using 32 KiB central windows. |
| `variant_transformer_pigmentation_binary` | Materialized variant transformer dataset for the pigmentation target. |

Use `genomics audit-data` to list the registered datasets that are available in the installed package:

```bash
genomics audit-data
```

## Canonical Dataset Layout

The default `1kg_high_coverage` dataset is expected to follow this layout for predictor and dataset-builder compatibility:

```text
<dataset>/
  dataset_metadata.json
  layout_metadata.json
  references/
    windows/
      <target>/
        ref.window.fa
        window_metadata.json
  individuals/
    <sample>/
      individual_metadata.json
      windows/
        <target>/
          predictions_H1/
          predictions_H2/
          <sample>.H1.window.fixed.fa
          <sample>.H2.window.fixed.fa
```

`dataset_metadata.json` is the dataset-level index. It records sample IDs, target metadata, pedigree or family data when available, window catalogs, and raw VCF source hints.

`layout_metadata.json` marks the dataset as following the active canonical layout. It is used as a compatibility signal by audits and readers.

`references/windows/<target>/window_metadata.json` records the reference coordinate context for a target window, including chromosome, start/end coordinates, target name, and source VCF metadata when available.

`individuals/<sample>/individual_metadata.json` stores sample-level labels such as population, superpopulation, family, sex, and configured derived targets.

## BCFtools Chain Artifacts

The genotype predictor's `haplotype_channels` layout requires consensus and chain-derived artifacts for each selected sample/window when `alignment_mapping: bcftools_chain` is used.

For each sample and target window, the audit expects these files:

```text
individuals/<sample>/windows/<target>/
  <sample>.H1.window.raw.fa
  <sample>.H2.window.raw.fa
  <sample>.window.consensus_ready.vcf.gz
  <sample>.window.consensus_ready.vcf.gz.tbi
  <sample>.window.vcf.gz
  <sample>.window.vcf.gz.tbi
```

Check a small subset before training aligned genotype models:

```bash
genomics audit-data \
  --dataset-id 1kg_high_coverage \
  --check-bcftools-chain \
  --sample-limit 3 \
  --fail-on-missing
```

Use `--genes` to restrict the check to selected windows:

```bash
genomics audit-data \
  --dataset-id 1kg_high_coverage \
  --check-bcftools-chain \
  --genes HBB BRCA1 \
  --fail-on-missing
```

## Syncing Chain Artifacts

If consensus/chain artifacts exist in a separate dataset tree, sync them into the canonical dataset with:

```bash
genomics genotype sync-bcftools-artifacts \
  --source-dir /path/to/consensus_dataset \
  --target-dir /path/to/canonical_dataset \
  --link-mode hardlink \
  --sample-limit 3
```

The command previews actions by default. Add `--apply` to modify the target dataset:

```bash
genomics genotype sync-bcftools-artifacts \
  --source-dir /path/to/consensus_dataset \
  --target-dir /path/to/canonical_dataset \
  --link-mode hardlink \
  --apply
```

`--link-mode hardlink` is preferred when source and target are on the same filesystem. Use `symlink` or `copy` when hardlinks are not possible.

## Config Usage

Genotype configs use `dataset_input.dataset_id`:

```yaml
dataset_input:
  dataset_id: "1kg_high_coverage"
```

Variant transformer configs use `dataset.source_dataset_id` for the source canonical dataset and `dataset.processed_dir` for the materialized sparse-token output:

```yaml
dataset:
  source_dataset_id: "1kg_high_coverage"
  processed_dir: "/path/to/materialized/variant_transformer"
```

## Audit Behavior

`genomics audit-data` checks that registered dataset paths exist and that required top-level metadata and layout directories are present. With `--check-bcftools-chain`, it also checks the aligned genotype artifacts listed above.

`--fail-on-missing` returns a non-zero exit code if any checked dataset has missing required artifacts. This is suitable for CI or pre-run validation.
