# Config Schemas

Large YAML configs are easier to maintain when their fields are discoverable from the command line.

The `genomics config` command exposes typed schemas for config families that already have Pydantic models.

Currently supported kinds:

| Kind | Component |
|---|---|
| `genotype` | Genotype-based predictor |
| `variant` | Variant transformer predictor |

## Describe Fields

Print every known field, including nested sections, types, defaults, required status, and available descriptions:

```bash
genomics config describe genotype
genomics config describe variant
```

Machine-readable output is available with:

```bash
genomics config describe genotype --json
```

## Validate YAML

Validate a config file against its typed schema:

```bash
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics config validate configs/predictors/variant_transformer/repo_layout.example.yaml
```

The command infers the kind from canonical config paths. If the file lives outside the canonical tree, pass `--kind` explicitly:

```bash
genomics config validate /tmp/experiment.yaml --kind genotype
```

Use JSON output for automation:

```bash
genomics config validate configs/predictors/variant_transformer/repo_layout.example.yaml --json
```

## JSON Schema

Export JSON Schema for editor integration or external tools:

```bash
genomics config schema genotype > genotype.schema.json
genomics config schema variant > variant.schema.json
```

## Genotype Preprocessing Fit Subsampling

Genotype configs can reduce preprocessing cost by fitting normalization and sklearn PCA on stratified subsets while still applying the fitted transform to every sample.

Normalization parameters are controlled under `dataset_input`:

```yaml
dataset_input:
  normalization_fit_splits: ["train"]      # stricter; ["train", "val"] preserves historical behavior
  normalization_fit_sample_fraction: 0.5
  normalization_fit_min_samples: 700
  normalization_fit_min_samples_per_class: 75
  normalization_fit_random_seed: 13
  normalization_fit_stratify: true
```

Use `normalization_fit_splits: ["train"]` when validation should remain unseen by fitted preprocessing parameters. The fitted normalization is still applied to train, validation, and test. `log` and `zscore` usually tolerate moderate stratified subsampling. `minmax_keep_zero` is more sensitive to missed extremes; prefer larger fractions such as `0.8` or keep `1.0`.

Sklearn scaler/PCA fitting is controlled under `model.sklearn`:

```yaml
model:
  sklearn:
    pca_fit_sample_fraction: 0.5
    pca_fit_min_samples: 700
    pca_fit_min_samples_per_class: 75
    pca_fit_random_seed: 13
    pca_fit_stratify: true
```

Only the scaler/PCA fit uses the subset. The PCA transform, classifier training, validation, and test evaluation still use the full split arrays. Changing either normalization-fit or PCA-fit parameters invalidates the corresponding cache so stale transforms are not reused. See [Genotype Predictor: Dataset Processing Flow](../components/genotype-predictor.md#dataset-processing-flow) for the full processing sequence and the intermediate artifact paths to inspect.

See [Preprocessing And Leakage](../concepts/preprocessing-and-leakage.md) for the evaluation rationale behind these fields.

## Current Scope

This first schema layer covers `genotype` and `variant` configs because those components already use typed Pydantic configuration models.

Configs for SNP ancestry, Genomes Analyzer, dataset builders, AlphaGenome, and converters still use more free-form YAML loading in parts of the codebase. Add typed models for those components before exposing them through `genomics config`.
