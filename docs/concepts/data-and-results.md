# Data And Results

The workspace separates source data, derived datasets, and run outputs.

## Data Roots

Canonical data paths are resolved by `genomics.workspace` and `genomics.core.data_registry`.

Important environment variables:

| Variable | Purpose |
|---|---|
| `GENOMICS_DATA_ROOT` | Override the default root for datasets |
| `GENOMICS_RESULTS_ROOT` | Override the default root for results |

## Registered Dataset IDs

Use dataset IDs in configs when possible instead of hard-coded absolute paths.

Common IDs include:

| Dataset ID | Purpose |
|---|---|
| `1kg_high_coverage` | Canonical high-coverage 1000 Genomes dataset |
| `legacy_top3_1kg_high_coverage` | Historical `top3` dataset, kept for provenance only |
| `variant_transformer_superpopulation` | Materialized variant transformer dataset |
| `variant_transformer_superpopulation_32k` | 32k-window variant transformer dataset |
| `variant_transformer_pigmentation_binary` | Variant transformer pigmentation target dataset |

Validate data availability with:

```bash
genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing
```

See [Data Registry](../reference/data-registry.md) for the canonical dataset layout, registered dataset IDs, and `bcftools_chain` artifact checks.

## Results

Runs generally write under:

```text
results/<pipeline>/runs/<run_name>/
```

Run directories usually contain:

| File or directory | Purpose |
|---|---|
| `manifest.json` | Run metadata and status |
| `config.yaml` | Copy of the config used |
| `resolved_config.json` | Expanded/validated config |
| `models/` | Checkpoints and model artifacts |
| `plots/` | Figures and diagnostics |
| `reports/` | Reports and metrics |
