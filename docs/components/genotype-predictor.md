# Genotype Predictor

The genotype predictor trains dense/aligned models from processed genotype and AlphaGenome-derived tensors.

For the reproduction commands, read [Training And Running Predictions](../getting-started/training-and-running-predictions.md). DITA background is covered in [Dynamic Indel Tensor Alignment](../concepts/dynamic-indel-tensor-alignment.md).

## CLI

```bash
genomics genotype prepare-cache configs/predictors/genotype_based/repo_layout.example.yaml
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype evaluate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy --split test
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/predictors/genotype_based/` |
| Configs | `configs/predictors/genotype_based/` |
| Legacy baseline configs | `configs/predictors/genotype_based/neural_legacy/` |

## Tensor Layouts

| Layout | Use |
|---|---|
| `haplotype_channels` | Canonical aligned layout for new experiments |
| `raw_center_crop` | Legacy baseline layout, mostly for historical comparisons |

See [Dynamic Indel Tensor Alignment](../concepts/dynamic-indel-tensor-alignment.md) for the coordinate model, mask rows, and required `bcftools_chain` artifacts behind the aligned layout.

## Implementation Structure

| Area | Modules | Responsibility |
|---|---|---|
| Config | `config.py` | Parses predictor YAML into a structured pipeline config |
| Data pipeline | `data/pipeline.py`, `data/processed_dataset.py`, `data/genomic_dataset.py` | Loads source datasets, prepares processed caches, returns PyTorch datasets |
| Splitting | `data/splitting.py`, `genomics.core.splitting` | Builds deterministic train/val/test splits, including family-aware splits |
| Tensor layout | `data/layout.py`, `data/normalization.py` | Defines tensor channel conventions and normalization behavior |
| Alignment | `alignment/*.py` | Builds dynamic indel-aligned tensors and bcftools chain-based mappings |
| Models | `models/cnn_model.py`, `models/cnn2_model.py`, `models/nn_model.py`, `models/sklearn_models.py` | CNN/MLP/sklearn baselines |
| Training | `experiments/train.py`, `experiments/training.py`, `genomics.core.training_utils` | Training loops, schedulers, checkpoints, metrics |
| Evaluation | `experiments/evaluate_checkpoint.py`, `experiments/evaluation.py` | Loads checkpoints and writes split metrics |
| Analysis | `analysis/*.py` | PCA plots, UMAP/shard plots, interpretability helpers |
| Apps | `apps/*.py` | Workbench and lightweight viewers |

## Data Representation

The genotype predictor consumes per-sample, per-window tensors. The current canonical direction is `haplotype_channels`, where haplotypes are represented as channels and indel-aware alignment keeps AlphaGenome signals, masks, and sequence-derived features on a comparable coordinate axis.

The older `raw_center_crop` layout is preserved for baselines and for reproducing older results. It crops fixed windows around target positions and avoids dynamic indel alignment.

Feature channels commonly include:

- AlphaGenome signal tracks;
- insertion masks;
- deletion masks;
- validity masks;
- SNP masks;
- optional transformed or reference-normalized signals.

## End-To-End Inputs

The pipeline starts from the canonical dataset layout registered in `genomics.core.data_registry`. A config can point to it with `dataset_input.dataset_id` or directly with `dataset_input.dataset_dir`. For the default 1000 Genomes high-coverage dataset, the reader expects this information to be present in the dataset tree:

| Input | Where It Is Read | Used For |
|---|---|---|
| `dataset_metadata.json` | Dataset root | Sample list, pedigree/population metadata, target distributions, raw VCF source hints, and window catalog fallback |
| `layout_metadata.json` | Dataset root | Confirms that the directory follows the canonical active layout |
| `references/windows/<gene>/ref.window.fa` | Reference window directory | Reference sequence and coordinate context for each gene/window |
| `references/windows/<gene>/window_metadata.json` | Reference window directory | Chromosome, start/end coordinates, target gene/window name, and raw VCF source metadata |
| `individuals/<sample>/individual_metadata.json` | Sample directory | Sample-level metadata such as family, population, superpopulation, and optional extra targets |
| `individuals/<sample>/windows/<gene>/predictions_H1/` | Sample window directory | AlphaGenome prediction tracks for haplotype 1 |
| `individuals/<sample>/windows/<gene>/predictions_H2/` | Sample window directory | AlphaGenome prediction tracks for haplotype 2 |
| `individuals/<sample>/windows/<gene>/<sample>.H1.window.fixed.fa` | Sample window directory | Haplotype 1 sequence, mainly for legacy/raw sequence-compatible paths |
| `individuals/<sample>/windows/<gene>/<sample>.H2.window.fixed.fa` | Sample window directory | Haplotype 2 sequence, mainly for legacy/raw sequence-compatible paths |
| Consensus/chain artifacts | Canonical dataset plus `consensus_dataset_dir` | Required by `tensor_layout: haplotype_channels` to align AlphaGenome signal coordinates across indels |

The config narrows this raw dataset into a logical view before any tensor is built. Important selectors are:

| Config Field | Effect |
|---|---|
| `genes_to_use` | Selects the gene/window directories that become rows in the final tensor |
| `alphagenome_outputs` | Selects prediction output families such as RNA-seq or ATAC-style tracks |
| `ontology_terms` | Filters multi-track AlphaGenome outputs to a specific ontology/cell-type subset when metadata is available |
| `sample_ids` / `sample_ids_path` | Limits the sample universe to an explicit list |
| `superpopulations_to_use` / `populations_to_use` | Limits samples by metadata-driven ancestry groups |
| `haplotype_mode` | Uses only `H1`, only `H2`, or both haplotypes as model input |
| `output.prediction_target` | Chooses the supervised label, for example `superpopulation`, `population`, `frog_likelihood`, or a configured derived target |
| `output.known_classes` | Fixes class order and avoids discovering classes by scanning the dataset |
| `output.derived_targets` | Defines labels by mapping values from another metadata field into new classes |

## Transformations Before Modeling

The data path can be read as a sequence of deterministic transformations:

```text
canonical dataset + YAML config
  -> logical view selection
  -> target/class mapping
  -> family-aware train/val/test split
  -> per-sample, per-gene window extraction
  -> tensor layout construction
  -> signal transformation and masking
  -> normalization parameter fit
  -> normalized processed tensor cache
  -> DataLoaders or sklearn PCA arrays
  -> model training/evaluation
```

### Logical View And Split

`prepare_data` first resolves the dataset path, loads `GenomicDataset`, creates the target mapping, removes samples without valid targets, and builds the train/validation/test split. With `family_split_mode: family_aware`, all members of a family are kept in the same split to reduce leakage through relatives. The resulting sample IDs are saved in `split_index.json` and the resolved view is saved in `view_definition.json`.

See [Preprocessing And Leakage](../concepts/preprocessing-and-leakage.md) for split roles, normalization fitting, PCA fitting, stability analysis, and confidence interval guidance.

### Window Extraction

For each selected sample and gene, the base reader loads the configured AlphaGenome predictions from `predictions_H1` and `predictions_H2`. The pipeline then selects the configured output tracks, optionally filters ontology terms, takes the configured central window (`window_center_size`), and optionally applies `downsample_factor`.

In `raw_center_crop`, this step is intentionally simple: it crops the fixed window and stacks the selected tracks in the historical row order. This representation is useful for legacy comparisons but does not explicitly repair coordinate shifts caused by indels.

In `haplotype_channels`, the pipeline builds an indel-aware coordinate axis for each gene. The bcftools chain mapper reads consensus/chain artifacts for each sample and haplotype, maps AlphaGenome prediction positions back into a shared expanded axis, and creates comparable signal rows across samples even when insertions or deletions are present.

### Feature Channels

The final row count is driven by the selected genes, haplotypes, ontology tracks, and mask options. Conceptually, each gene contributes rows like this:

| Setting | Rows Per Gene |
|---|---|
| `feature_mode: signals_only` | AlphaGenome signal rows for each selected haplotype and ontology/output track |
| `feature_mode: masks_only` | Indel/SNP/validity mask rows only |
| `feature_mode: signals_and_masks` | Signal rows plus mask rows |
| `haplotype_mode: H1+H2` | Doubles haplotype-specific rows compared with a single-haplotype mode |
| `indel_include_valid_mask: true` | Adds validity mask rows indicating aligned positions with real mapped signal |
| `indel_include_snp_mask: true` | Adds SNP mask rows in addition to insertion/deletion masks |

Mask rows identify where the aligned axis contains insertions, deletions, valid mapped signal, and optionally SNP positions. Signal rows contain AlphaGenome values after any configured signal transform.

### Signal Transform And Masking

`alphagenome_signal_transform` controls the value before normalization:

| Transform | Meaning |
|---|---|
| `absolute` | Use the AlphaGenome prediction value as stored in the dataset |
| `delta_reference` | Subtract the configured reference-only prediction for the same output/track/window, producing a sample-vs-reference signal delta |

`alphagenome_signal_variant_mask: true` zeros signal values at positions that do not overlap a SNP/indel in the selected sample universe. This creates a representation focused on variable positions rather than carrying all dense AlphaGenome signal positions.

### Normalization

Normalization parameters are fitted once, then applied to every split. The fit subset is controlled separately from the train/validation/test tensors by `normalization_fit_splits`, `normalization_fit_sample_fraction`, minimum sample settings, and stratification settings.

| Method | Fitted Parameter | Applied Transformation |
|---|---|---|
| `zscore` | Per-track mean and standard deviation | `(x - mean) / std` |
| `minmax_keep_zero` | Per-track positive maximum | `x / max`, preserving exact zeros |
| `log` | Per-track `log1p(max_positive)` | `log1p(x) / log1p(max_positive)` |

The test split is not needed to fit normalization by default. If you include `test` in a fit set through future code changes, treat that as a deliberate leakage-prone experiment rather than a standard evaluation setup.

### Processed Sample Shape

Neural models receive tensors as either `(rows, length)` or `(haplotypes, rows, length)` before batching. The model implementations flatten the haplotype dimension into rows when needed, so training commonly sees batch tensors equivalent to `(batch, rows, length)`.

Sklearn baselines flatten each processed tensor to one vector per sample, fit/transform it with `StandardScaler` and PCA, then train a classifier on the resulting component matrix.

## Cache And Run Artifacts

`prepare-cache` materializes processed data so training does not repeatedly parse heavy genomic inputs. Training creates experiment directories with configs, manifests, checkpoints, and metrics through shared `genomics.core` infrastructure.

Typical artifacts:

| Artifact | Purpose |
|---|---|
| processed cache directory | Reusable tensors and split metadata |
| `split_index.json` | Train/val/test sample split |
| `manifest.json` | Run metadata |
| `models/best_accuracy.pt` | Best accuracy checkpoint |
| `models/best_loss.pt` | Best loss checkpoint |
| sklearn `.joblib` files | Serialized sklearn baselines and PCA/scaler bundles |

## Dataset Processing Flow

The genotype data pipeline is implemented mainly in `src/genomics/predictors/genotype_based/data/pipeline.py` and `src/genomics/predictors/genotype_based/data/processed_dataset.py`. It has two distinct preprocessing moments:

- Fit preprocessing parameters, such as normalization statistics and sklearn scaler/PCA, optionally on a subset.
- Apply the fitted transforms to the full train/validation/test splits and materialize reusable artifacts.

For sklearn baselines, this distinction is important: `normalization_fit_sample_fraction` and `pca_fit_sample_fraction` reduce only the fitting set. The processed tensor cache, PCA transform arrays, classifier training, validation, and test evaluation still use the full split data.

### Resolve The Concrete Paths

Use the config helpers to print the exact cache and run paths for a genotype YAML:

```bash
python3 - <<'PY'
from pathlib import Path
from genomics.predictors.genotype_based.config import load_config, get_dataset_cache_dir, generate_experiment_name

cfg = load_config(Path("configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml"))
experiment_name = generate_experiment_name(cfg)

print("experiment_name:", experiment_name)
print("dataset_cache_dir:", get_dataset_cache_dir(cfg))
print("experiment_dir:", Path(cfg.dataset_input.results_dir) / experiment_name)
PY
```

The processed dataset cache has this shape:

```text
<processed_cache_dir>/datasets/<dataset_name>/
```

For example, with `processed_cache_dir: results/cache/genotype_based_predictor`:

```text
results/cache/genotype_based_predictor/datasets/<dataset_name>/
```

The experiment directory has this shape:

```text
<results_dir>/<experiment_name>/
```

For example, with `results_dir: results/genotype_based_predictor/icann/runs`:

```text
results/genotype_based_predictor/icann/runs/<experiment_name>/
```

The sklearn PCA cache is stored next to the processed dataset caches:

```text
<processed_cache_dir>/pca_cache/<dataset_name>_pca<components>[_<backend>]/
```

For a randomized streaming PCA with 300 components:

```text
results/cache/genotype_based_predictor/pca_cache/<dataset_name>_pca300_randomized_streaming/
```

### Experiment Directory

`genomics genotype train ...` first creates an experiment directory through `setup_experiment_run`. Check these files to confirm the resolved run configuration:

| Path | Meaning |
|---|---|
| `<experiment_dir>/config.yaml` | Copy of the YAML passed on the command line |
| `<experiment_dir>/resolved_config.json` | Fully validated/resolved config, including defaults and resolved dataset paths |
| `<experiment_dir>/manifest.json` | Run metadata, status, artifact pointers, environment, and later model/PCA paths |
| `<experiment_dir>/models/` | Checkpoints, sklearn artifacts, and copied normalization parameters |
| `<experiment_dir>/plots/` | Classification plots and diagnostics |
| `<experiment_dir>/reports/` | Reserved reports directory |

For sklearn baselines, the trained classifier artifact is saved as:

```text
<experiment_dir>/models/sklearn_baseline.joblib
```

Validation metrics are saved as:

```text
<experiment_dir>/val_results.json
```

### Processed Dataset Cache Validation

Before rebuilding data, `prepare_data` checks whether the processed cache is complete and compatible with the current config. A valid processed cache contains:

| Path | Meaning |
|---|---|
| `<dataset_cache_dir>/.cache_complete` | Completion flag; absent means an interrupted/incomplete cache |
| `<dataset_cache_dir>/metadata.json` | Cache version, sample counts, input shape, processing parameters, split sizes |
| `<dataset_cache_dir>/normalization_params.json` | Fitted normalization parameters used for all splits |
| `<dataset_cache_dir>/split_index.json` | Train/validation/test sample IDs |
| `<dataset_cache_dir>/view_definition.json` | Requested and resolved dataset view definition |
| `<dataset_cache_dir>/shards_index.json` | Shard layout when `cache_processed_tensors: true` |

`metadata.json` is the quickest way to verify which preprocessing parameters produced the cache. Useful fields include:

```json
{
  "processing_params": {
    "normalization_method": "log",
    "normalization_fit_splits": ["train"],
    "normalization_fit_sample_fraction": 0.3,
    "normalization_fit_min_samples": 700,
    "normalization_fit_min_samples_per_class": 75,
    "normalization_fit_random_seed": 13,
    "normalization_fit_stratify": true,
    "cache_processed_tensors": true,
    "processed_shard_size": 64
  },
  "splits": {
    "train_size": 0,
    "val_size": 0,
    "test_size": 0,
    "random_seed": 13
  }
}
```

Changing normalization-fit parameters invalidates the processed cache so stale normalization parameters are not reused.

### Split And Normalization Fit

When the processed cache is missing or invalid, the pipeline rebuilds it from the source dataset:

1. Resolve `dataset_input.dataset_id` or `dataset_input.dataset_dir` to a concrete dataset directory.
2. Build valid sample indices from the source dataset.
3. Build deterministic train/validation/test groups, including family-aware grouping when configured.
4. Select the base indices used only to fit normalization parameters.
5. Fit normalization on that subset.
6. Build a full `ProcessedGenomicDataset` using the fitted normalization parameters.
7. Materialize the full train/validation/test processed cache.

The normalization subset is controlled by:

```yaml
dataset_input:
  normalization_fit_splits: ["train"]
  normalization_fit_sample_fraction: 0.3
  normalization_fit_min_samples: 700
  normalization_fit_min_samples_per_class: 75
  normalization_fit_random_seed: 13
  normalization_fit_stratify: true
```

The console logs the selected fit size, for example:

```text
Normalization fit subset: 700/<train-size> amostras (splits=['train'], fraction=0.3, stratify=True)
```

For `normalization_method: log` with `normalization_value: 0.0`, the pipeline scans the selected samples and records per-track positive maxima as `log_max` values. The result is written to:

```text
<dataset_cache_dir>/normalization_params.json
<experiment_dir>/models/normalization_params.json
```

The file has this general structure:

```json
{
  "method": "log",
  "per_track": true,
  "num_tracks": 0,
  "track_params": [
    {"log_max": 0.0}
  ]
}
```

The values above are illustrative; inspect the actual file for the real track count and fitted maxima.

### Processed Tensor Materialization

After fitting normalization parameters, the pipeline applies normalization while saving the full processed dataset cache. This is separate from the normalization fit subset. With `cache_processed_tensors: true`, each split is saved as `.pt` shards.

During an in-progress rebuild, files are written under a temporary directory:

```text
<dataset_cache_dir>_tmp_<pid>/
```

Only after a successful build is the temporary directory renamed to:

```text
<dataset_cache_dir>/
```

Shard files are named like this:

```text
<dataset_cache_dir>/train_data_shard_00000.pt
<dataset_cache_dir>/train_data_shard_00001.pt
<dataset_cache_dir>/val_data_shard_00000.pt
<dataset_cache_dir>/test_data_shard_00000.pt
```

`shards_index.json` records which files belong to each split:

```json
{
  "train": {
    "shard_size": 64,
    "num_samples": 0,
    "shards": ["train_data_shard_00000.pt"]
  },
  "val": {
    "shard_size": 64,
    "num_samples": 0,
    "shards": ["val_data_shard_00000.pt"]
  },
  "test": {
    "shard_size": 64,
    "num_samples": 0,
    "shards": ["test_data_shard_00000.pt"]
  }
}
```

Each shard is a `torch.save(...)` file containing a Python list of `(features_tensor, target_tensor)` pairs. Reading one item from a lazy sharded cache loads the containing shard with `torch.load(...)`, so random subset access can still touch many shards even when the requested sample count is small.

The cache build log includes split-level profiling lines such as:

```text
[profile train 256/...] elapsed=... | base_fetch=... | window_process=... | normalize=...
```

Use those counters to distinguish where rebuild time is spent:

- `base_fetch`: reading source dataset entries.
- `window_process`: extracting/cropping/aligning windows into model tensors.
- `normalize`: applying already fitted normalization parameters while saving tensors.

### Loading The Processed Cache

When the processed cache is valid, `load_processed_dataset` creates `CachedProcessedDataset` instances for train, validation, and test. In lazy mode, the dataset reads shard files on demand according to `shards_index.json`.

For sklearn baselines, DataLoaders use `num_workers=0`. This keeps loading deterministic and avoids multiplying memory use, but it also means shard deserialization and flattening happen serially in the main process.

### Sklearn PCA Cache

For sklearn baselines with `model.sklearn.use_pca_cache: true`, PCA artifacts are cached separately from processed tensors.

The PCA fit subset is controlled by:

```yaml
model:
  sklearn:
    pca_fit_sample_fraction: 0.3
    pca_fit_min_samples: 700
    pca_fit_min_samples_per_class: 75
    pca_fit_random_seed: 13
    pca_fit_stratify: true
```

The PCA cache directory contains:

| Path | Meaning |
|---|---|
| `<pca_cache_dir>/.pca_cache_complete` | Completion flag |
| `<pca_cache_dir>/pca_metadata.json` | Fit size, feature count, PCA backend, randomized PCA settings, split sizes |
| `<pca_cache_dir>/scaler_pca.joblib` | Fitted scaler and PCA object |
| `<pca_cache_dir>/X_train.npy` | Full train split transformed to PCA components |
| `<pca_cache_dir>/y_train.npy` | Train labels |
| `<pca_cache_dir>/X_val.npy` | Full validation split transformed to PCA components |
| `<pca_cache_dir>/y_val.npy` | Validation labels |
| `<pca_cache_dir>/X_test.npy` | Full test split transformed to PCA components |
| `<pca_cache_dir>/y_test.npy` | Test labels |
| `<pca_cache_dir>/randomized_pca_components.<dtype>.memmap` | Streaming randomized PCA components, when using `randomized_streaming` |

`pca_metadata.json` is the best place to confirm the PCA fitting subset:

```json
{
  "pca_components_requested": 300,
  "pca_n_components_effective": 300,
  "n_features_original": 0,
  "n_train": 0,
  "n_fit": 700,
  "n_val": 0,
  "n_test": 0,
  "pca_backend": "randomized_streaming",
  "pca_fit_sample_fraction": 0.3,
  "pca_fit_min_samples": 700,
  "pca_fit_min_samples_per_class": 75,
  "randomized_pca_n_iter": 2,
  "randomized_pca_feature_chunk_size": 8192,
  "randomized_pca_dtype": "float64"
}
```

The PCA fit subset is used for fitting the scaler/PCA only. `X_train.npy`, `X_val.npy`, and `X_test.npy` are still full transformed splits.

### Stability Analysis Outputs

When running `genomics genotype stability`, the initial processed cache is prepared in the same way, but each stability split fits its own scaler/PCA and classifier. The output root comes from `stability_analysis.output_dir`.

Typical files are:

```text
<stability_output_dir>/prepare_cache/
<stability_output_dir>/pca/<run_name>/
<stability_output_dir>/runs/<run_name>/val_results.json
<stability_output_dir>/stability_results.json
<stability_output_dir>/stability_results.csv
<stability_output_dir>/split_plan_indices.json
```

Use the per-run PCA directories to inspect intermediate PCA artifacts for each repeated split.

## BCFtools Chain Alignment

The `alignment/bcftools_chain_mapper.py` path rebuilds haplotype consensus sequences and reads BCFtools chain files to map AlphaGenome predictions back to a global aligned coordinate system. `sync-bcftools-artifacts` copies or links the required consensus and chain artifacts into the canonical dataset tree.

## Available Models

The active model is selected with `model.type`. All model types consume the same processed dataset cache; the difference is how they transform the cached tensor into predictions.

| `model.type` | Implementation | Input Handling | Main Use |
|---|---|---|---|
| `NN` | `models/nn_model.py` | Flattens the processed tensor into one dense vector, then applies configurable fully connected layers | Simple neural baseline for small/medium feature spaces |
| `CNN` | `models/cnn_model.py` | Treats rows x positions as a 2D image with one channel, applies one Conv2D block, optional pooling, then an MLP head | Lightweight convolutional baseline |
| `CNN2` | `models/cnn2_model.py` | Uses a multi-stage Conv2D design: stage 1 groups tracks across each gene, stages 2/3 scan along the genomic axis, then global pooling and a classifier head | Main dense neural architecture for gene x track x position tensors |
| `SVM` | `models/sklearn_models.py` | Flattens tensors, applies `StandardScaler`, PCA, then `LinearSVC` | Linear margin baseline after dimensionality reduction |
| `LOGREG` | `models/sklearn_models.py` | Flattens tensors, applies `StandardScaler`, PCA, then multinomial/binary logistic regression | Interpretable linear baseline and calibration-friendly comparisons |
| `RF` | `models/sklearn_models.py` | Flattens tensors, applies `StandardScaler`, PCA, then `RandomForestClassifier` | Nonlinear tree baseline on PCA components |
| `XGBOOST` | `models/sklearn_models.py` | Flattens tensors, applies `StandardScaler`, PCA, then `xgboost.XGBClassifier` | Gradient-boosted tree baseline when `xgboost` is installed |

### Neural Model Inputs

The neural models are trained through the shared PyTorch loop in `experiments/training.py`. The processed cache provides `(features_tensor, target_tensor)` pairs. During batching, `features_tensor` becomes a batch tensor and `target_tensor` becomes the class index or regression target.

The common neural hyperparameters are:

| Field | Applies To | Meaning |
|---|---|---|
| `model.hidden_layers` | `NN`, `CNN` | Hidden layer sizes in the MLP/head |
| `model.activation` | `NN`, `CNN` | Activation for hidden layers: `relu`, `tanh`, or `sigmoid` |
| `model.dropout_rate` | `NN`, `CNN`, `CNN2` | Dropout after hidden/head layers |
| `training.optimizer` | Neural models | `adam`, `adamw`, or `sgd` |
| `training.learning_rate` | Neural models | Optimizer learning rate |
| `training.batch_size` | Neural models | Samples per optimization batch |
| `lr_scheduler.*` | Neural models | Optional learning-rate schedule |

`CNN2` has dedicated geometry settings under `model.cnn2`. `kernel_stage1` and `stride_stage1` define how rows are grouped at the first stage. For gene-organized tensors, the vertical kernel is usually chosen to span all rows for one gene, so the first convolution sees the complete per-gene feature block. `kernel_stages23`, `stride_stages23`, and `padding_stages23` then process along the genomic coordinate axis. `global_pool_type` collapses the remaining width before the final classifier.

### Sklearn Model Inputs

The sklearn baselines are not trained through gradient descent. Their data path is:

```text
processed tensor cache
  -> flatten each sample to one feature vector
  -> StandardScaler fit on train/PCA-fit subset
  -> PCA fit on train/PCA-fit subset
  -> transform full train, val, and test splits
  -> train classifier on full transformed train split
  -> evaluate validation or test split
```

PCA settings live under `model.sklearn`:

| Field | Meaning |
|---|---|
| `pca_components` | Requested number of PCA components; the effective value is capped by sample/feature limits |
| `use_pca_cache` | Reuses transformed `X_train.npy`, `X_val.npy`, and `X_test.npy` when compatible |
| `pca_backend` | `incremental` or `randomized_streaming` |
| `pca_fit_sample_fraction` | Fraction of the training split used only to fit scaler/PCA |
| `pca_fit_min_samples` | Minimum total samples for scaler/PCA fit subset |
| `pca_fit_min_samples_per_class` | Minimum per-class samples for stratified fit subset |
| `pca_fit_stratify` | Preserves class balance while selecting the fit subset |

Classifier-specific settings are grouped under `model.sklearn.svm`, `model.sklearn.logistic_regression`, `model.sklearn.random_forest`, and `model.sklearn.xgboost`.

### Model Selection Pattern

Use the validation split, stability analysis, and confidence intervals for model selection. Keep the held-out test split for the final selected configuration. The intended sequence for sklearn experiments is usually:

```bash
genomics genotype search configs/predictors/genotype_based/icann/search_rf_xgboost.yaml
genomics genotype stability configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml
genomics genotype test configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml
```

## Useful Subcommands

| Command | Purpose |
|---|---|
| `prepare-cache` | Materialize processed cache before training |
| `train` | Train a model from config |
| `evaluate` | Evaluate a checkpoint on a split |
| `pca-variance` | Analyze PCA variance for sklearn baselines |
| `workbench` | Launch interactive workbench/dashboard |
| `sync-bcftools-artifacts` | Materialize chain/consensus artifacts |
| `single-gene-screen` | Run single-gene ablation/screen experiments |

## Adding A New Model

Add the model implementation under `models/`, expose it in the model factory/config path, and reuse the existing `experiments/training.py` and `genomics.core.training_utils` loop unless the model needs a genuinely different training contract.

## Related Docs

- [Research Pipelines](../historical/RESEARCH_PIPELINES.md)
- [Important Experiments](../historical/IMPORTANT_EXPERIMENTS.md)
