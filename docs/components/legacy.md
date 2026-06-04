# Legacy Components

Legacy code is kept under `legacy/` for reproducibility. New development should not import legacy modules unless the work is explicitly about reproducing historical results.

## Legacy Neural Ancestry

Path: `legacy/neural_ancestry_predictor_deprecated/`

CLI wrapper:

```bash
genomics neural train legacy/neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml
genomics neural test legacy/neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml
genomics neural pca-cache legacy/neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml
```

For new ancestry/model experiments, prefer:

- `genomics genotype ...`
- `genomics variant ...`
- `genomics snp-ancestry ...`

### Implementation Structure

| Module | Purpose |
|---|---|
| `neural_ancestry_predictor_deprecated.py` | Historical monolithic CLI entrypoint now delegated by `genomics neural` |
| `config.py` | Legacy YAML parsing and dataset path resolution |
| `data_pipeline.py` | Legacy processed-cache creation and loading |
| `models.py` | Historical NN/CNN model definitions |
| `training.py` | Training loop and checkpoint handling |
| `evaluation.py` | Test/evaluation helpers |
| `sklearn_baselines.py` | SVM/RF/XGBoost baseline support |
| `tune_sklearn_baselines.py` | Hyperparameter search for sklearn baselines |
| `interpretability.py`, `annotate_deeplift_windows.py` | DeepLIFT/interpretability utilities |
| `verify_processed_dataset.py` | Dataset consistency checks |

The legacy package now imports shared infrastructure from `genomics.core` where possible, but the model/data assumptions remain historical.

### Why It Remains

This code produced historical benchmarks and plots. It remains available so old results can be reproduced or compared, but it is not the recommended implementation for new experiments.

Legacy-only configs live beside this code under `legacy/neural_ancestry_predictor_deprecated/configs/`. Configs for reproducing the old raw-center-crop representation on the active genotype predictor live under `configs/predictors/genotype_based/neural_legacy/`.

## Legacy Longevity Dataset

Path: `legacy/neural_longevity_dataset/`

This code remains coupled to historical `top3` layout assumptions and is not the recommended path for new dataset builds.

### Implementation Structure

| Module | Purpose |
|---|---|
| `neural_longevity_dataset.py` | Historical dataset construction pipeline |
| `longevity_train.py` | Historical training script for longevity experiments |

Use the active dataset builder under `src/genomics/workflows/dataset_builders/non_longevous/` for new dataset creation.
