# Preprocessing And Leakage

Genomics ML experiments often include preprocessing steps that learn parameters from data before model training. Examples include normalization, sklearn scaling, PCA, feature selection, and repeated resampling for stability analysis.

The core rule is simple: do not fit preprocessing parameters on the held-out test split before final evaluation.

## Split Roles

| Split | Role |
|---|---|
| Train | Fit model parameters and standard preprocessing parameters. |
| Validation | Select models, hyperparameters, checkpoints, and experiment variants. |
| Test | Final held-out evaluation after model selection. |

Family-aware splitting keeps related samples in the same split when family metadata is available. This reduces leakage through relatives and is the preferred mode for 1000 Genomes derived datasets.

## Genotype Normalization

Genotype configs control normalization fitting under `dataset_input`:

```yaml
dataset_input:
  normalization_fit_splits: ["train"]
  normalization_fit_sample_fraction: 0.3
  normalization_fit_min_samples: 700
  normalization_fit_min_samples_per_class: 75
  normalization_fit_random_seed: 13
  normalization_fit_stratify: true
```

Only the fit subset is used to learn normalization parameters. The learned transform is then applied to all train, validation, and test samples when the processed cache is materialized.

Use `normalization_fit_splits: ["train"]` for standard evaluation. Including validation can reproduce older behavior but makes validation less independent for model selection. Do not include `test` unless the run is explicitly leakage-prone or exploratory.

## Sklearn Scaler And PCA

Sklearn baselines flatten processed genotype tensors, fit a scaler/PCA bundle, and train a classifier. PCA fit sampling is controlled under `model.sklearn`:

```yaml
model:
  sklearn:
    pca_fit_sample_fraction: 0.5
    pca_fit_min_samples: 700
    pca_fit_min_samples_per_class: 75
    pca_fit_random_seed: 13
    pca_fit_stratify: true
```

The scaler/PCA fit can use a subset to reduce cost. The PCA transform, classifier training, validation, and test evaluation still use the full selected split data.

Changing normalization-fit or PCA-fit settings invalidates the corresponding cache so stale preprocessing parameters are not silently reused.

## Stability Analysis

`genomics genotype stability` estimates how sensitive sklearn results are to development-set resampling. It keeps the original test split fixed and resamples only the development data (`train+val`).

Supported strategies are configured with `stability_analysis.strategy`:

| Strategy | Purpose |
|---|---|
| `repeated_random_split` | Repeated train/validation draws from the development set. |
| `randomized_split` | A single randomized development split. |
| `cross_validation` | Cross-validation over the development set. |

Use the same `stability_analysis.split_plan_path` across comparable model configs to reuse the exact same sample-level plan. This makes model comparisons less sensitive to different random resampling plans.

Run final `genomics genotype test` only after model family and hyperparameters are selected.

## Confidence Intervals

`genomics genotype confidence-intervals` reloads a saved model artifact or checkpoint and recomputes metrics for a selected split. When configured, metrics can include bootstrap confidence intervals.

Use confidence intervals on the final test split only for the selected model. For exploratory comparisons, compute intervals on validation or stability outputs first.

## Practical Workflow

```bash
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype split configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype evaluate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy --split val
genomics genotype test configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy
```

For sklearn searches:

```bash
genomics genotype search configs/predictors/genotype_based/icann/search_rf_xgboost.yaml
genomics genotype stability configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml
genomics genotype test configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml --experiment-dir <selected-run>
```
