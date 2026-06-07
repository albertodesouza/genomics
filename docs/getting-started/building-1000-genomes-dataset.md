# Building 1000 Genomes Dataset

This step builds the canonical dataset used by the neural genotype predictor. It starts from 1000 Genomes metadata, phased VCFs, a reference genome, and AlphaGenome track generation.

Run the commands below in order. Use the links at the end of the page when you need conceptual or component-level details.

## 1. Inspect Data Configuration

Start with the default non-longevous dataset builder config. Confirm local paths such as `metadata_csv`, `reference.fasta`, `vcf_pattern`, and output settings before running a large job.

```bash
genomics dataset-builders non-longevous build --help
```

For larger or more specific builds, use one of the configs under:

```text
configs/workflows/non_longevous_dataset/
```

The config points to 1000 Genomes metadata, reference FASTA, VCF patterns, sample selection rules, AlphaGenome output settings, and output directories.

## 2. Discover AlphaGenome Tracks

Before large runs, export track metadata so you know which output families and ontology terms are available:

```bash
genomics alphagenome tracks \
  --api-key API_KEY \
  --output configs/workflows/alphagenome/tracks.json
```

For sequence sizing and output selection, see [AlphaGenome Sequences And Outputs](../guides/alphagenome-sequences.md).

## 3. Build The Dataset

Run the builder:

```bash
genomics dataset-builders non-longevous build \
  --config configs/workflows/non_longevous_dataset/default.yaml
```

The resulting canonical layout contains dataset metadata, reference windows, per-individual metadata, haplotype sequences, and AlphaGenome prediction directories.

## 4. Validate Dataset Availability

When the dataset is registered locally, validate it with:

```bash
genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing
```

If you plan to train aligned DITA models, also check BCFtools chain artifacts:

```bash
genomics audit-data \
  --dataset-id 1kg_high_coverage \
  --check-bcftools-chain \
  --sample-limit 3 \
  --fail-on-missing
```

## More Details

| Topic | Read More |
|---|---|
| Data roots and dataset IDs | [Data And Results](../concepts/data-and-results.md) |
| Canonical dataset layout and audit checks | [Data Registry](../reference/data-registry.md) |
| AlphaGenome workflow | [AlphaGenome](../components/alphagenome.md) |
| Dataset builder internals | [Dataset Builders](../components/dataset-builders.md) |
| Configuration conventions | [Configuration](../concepts/configuration.md) |
