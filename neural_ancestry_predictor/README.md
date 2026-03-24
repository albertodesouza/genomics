# Neural Ancestry Predictor

> **🧬 Ancestry Prediction using Neural Networks and AlphaGenome Data**

This module implements a YAML-configurable neural network that predicts ancestry (superpopulation, population, or FROG likelihood) from AlphaGenome predictions stored in a PyTorch dataset.

## 📑 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Architecture](#architecture)
- [Data Processing](#data-processing)
- [Processed Dataset Cache](#processed-dataset-cache)
- [Dataset Verification](#dataset-verification)
- [Training](#training)
- [Sklearn Baselines and PCA Cache](#sklearn-baselines-and-pca-cache)
- [Debug and Visualization](#debug-and-visualization)
- [Testing and Evaluation](#testing-and-evaluation)
- [Model Interpretability](#model-interpretability)
- [Variant Annotation](#variant-annotation-deeplift-post-processing)
- [Hypothesis Validation](#hypothesis-validation)
- [Weights & Biases](#weights--biases)
- [Hyperparameter Tuning](#hyperparameter-tuning)
- [FAQ](#faq)

---

## Overview

### What does this module do?

The **Neural Ancestry Predictor** trains a neural network to predict the genetic ancestry of individuals using:

- **Input**: AlphaGenome predictions (e.g., ATAC-seq, RNA-seq) from genomic windows
- **Output**: Superpopulation (AFR, AMR, EAS, EUR, SAS), Population (26 classes), or FROG likelihood (150 values)

### Features

- ✅ **Fully configurable via YAML**
- ✅ **Supports multiple targets** (superpopulation, population, FROG likelihood)
- ✅ **Flexible processing** of windows, haplotypes, and AlphaGenome outputs
- ✅ **Weights & Biases integration** for tracking and visualization
- ✅ **Automatic checkpointing** to save trained models
- ✅ **Detailed metrics** (accuracy, precision, recall, F1, confusion matrix)
- ✅ **Automatic normalization** with parameter caching
- ✅ **Interactive debug visualization** to inspect inputs and predictions sample-by-sample

---

## Installation

### 1. Dependencies

```bash
# Navigate to directory
cd genomics/neural_ancestry_predictor

# Install Python dependencies
pip3 install torch numpy pandas pyyaml scikit-learn rich

# Optional: Weights & Biases (for tracking)
pip3 install wandb
wandb login  # Authenticate with your W&B account
```

### 2. Dataset

This module requires a PyTorch dataset created by `build_non_longevous_dataset`:

```bash
# Example: dataset at /dados/GENOMICS_DATA/top3/non_longevous_results
# Must contain:
#   - individuals/
#   - dataset_metadata.json
```

See [build_non_longevous_dataset/docs/PYTORCH_DATASET.md](../build_non_longevous_dataset/docs/PYTORCH_DATASET.md) for more information.

---

## Quick Start

### Train Model

```bash
cd neural_ancestry_predictor
python3 neural_ancestry_predictor.py --config configs/default.yaml
```

### Test Model

```bash
python3 neural_ancestry_predictor.py --config configs/default.yaml --mode test
```

### Example Output

```
╭─────────────────────────────────────────╮
│ 🧬 Genomics                             │
│                                         │
│ Neural Ancestry Predictor               │
│ Mode: train                             │
│ Target: superpopulation                 │
│ Config: configs/default.yaml            │
╰─────────────────────────────────────────╯

Device: cuda
[INFO] GenomicLongevityDataset initialized:
  • Dataset: non_longevous_1000g
  • Individuals: 78
  • Load predictions: True
  • Load sequences: False

Computing normalization parameters...
Normalization: mean=0.123456, std=0.654321

Dataset split:
  • Train: 54 samples
  • Validation: 12 samples
  • Test: 12 samples

Model created:
  • Input size: 11000
  • Hidden layers: [128, 64]
  • Output size: 5
  • Activation: relu
  • Dropout: 0.2
  • Total parameters: 1,415,237

╭─────────────────────────────────────────╮
│ Starting Training                       │
│                                         │
│ Epochs: 100                             │
│ Batch size: 16                          │
│ Learning rate: 0.001                    │
╰─────────────────────────────────────────╯

Validation - Epoch 5: Loss=0.8234, Accuracy=0.7500
✓ Checkpoint saved: models/best_accuracy.pt

...

✓ Training completed!
```

---

## Configuration

### YAML File Structure

The `configs/default.yaml` file contains **8 main sections**:

#### A) Dataset Input Parameters

Controls **what** to load from the dataset and **how** to process it:

```yaml
dataset_input:
  dataset_dir: "/path/to/dataset"           # Path to PyTorch dataset
  alphagenome_outputs: ["ATAC"]             # Which outputs to use (RNA_SEQ, ATAC, CAGE, etc.)
  haplotype_mode: "H1+H2"                   # "H1", "H2" or "H1+H2"
  window_center_size: 100                   # Size of central region (bases)
  downsample_factor: 1                      # Downsampling factor (1 = none)
  processed_cache_dir: null                 # Cache directory (null = disabled)
```

**Processed Dataset Cache:**

The `processed_cache_dir` parameter enables caching of processed data to save time in subsequent runs:

- **`null`** (default): Cache disabled, always reprocess data
- **"processed_datasets/my_cache"**: Enables cache with specified directory

When cache is enabled:
1. **First run**: Processes data (normalization + splits) and saves to cache (~30-60 seconds for 78 samples)
2. **Subsequent runs**: Loads from cache if compatible (~2-5 seconds)
3. **Auto-invalidation**: Cache is rebuilt if parameters change (window_center_size, haplotype_mode, etc.)

**Benefits:**
- ⚡ **10-20x faster** startup for subsequent runs
- 🎯 **Reproducibility**: Same splits guaranteed across runs
- 💾 **Disk space**: ~10-50 MB per cache (depends on dataset size)

**Force Reprocessing:**
- Set to `null` in config, OR
- Delete cache directory manually

**Impact on Dimensionality:**

- Each window has ~1M bases per output
- `window_center_size=100` → extracts 100 bases from center
- `downsample_factor=2` → uses 1 in every 2 bases (reduces to 50)
- `haplotype_mode="H1+H2"` → doubles size (2 haplotypes)
- Final dimension = `n_windows × n_outputs × n_haplotypes × (window_center_size / downsample_factor)`

**Example:**
- 55 windows (SNPs)
- 1 output (ATAC)
- 2 haplotypes (H1+H2)
- window_center_size=100, downsample_factor=1
- **Dimension = 55 × 1 × 2 × 100 = 11,000 features**

#### B) Output Parameters

Defines **what** the network should predict:

```yaml
output:
  prediction_target: "superpopulation"  # "superpopulation", "population" or "frog_likelihood"
```

| Target | Type | Classes | Difficulty |
|--------|------|---------|------------|
| `superpopulation` | Classification | 5 (AFR, AMR, EAS, EUR, SAS) | Easy ⭐ |
| `population` | Classification | 26 | Medium ⭐⭐ |
| `frog_likelihood` | Regression | 150 values | Hard ⭐⭐⭐ |

#### C) Model Architecture

Defines the network **architecture**:

```yaml
model:
  hidden_layers: [128, 64]    # List of neurons per hidden layer
  activation: "relu"          # "relu", "tanh" or "sigmoid"
  dropout_rate: 0.2           # Dropout rate (0.0 to 1.0)
```

**Resulting architecture:**

```
Input (11000) → Dense(128) → ReLU → Dropout(0.2) →
Dense(64) → ReLU → Dropout(0.2) →
Dense(5) → Softmax
```

#### D) Training Parameters

Controls the **training process**:

```yaml
training:
  optimizer: "adam"              # "adam", "adamw" or "sgd"
  learning_rate: 0.001           # Learning rate
  loss_function: "cross_entropy" # "cross_entropy" or "mse"
  batch_size: 16                 # Number of samples per batch
  num_epochs: 100                # Number of epochs
  validation_frequency: 5        # Validate every N epochs
```

#### E) Data Split

Defines **dataset division**:

```yaml
data_split:
  train_split: 0.7      # 70% for training
  val_split: 0.15       # 15% for validation
  test_split: 0.15      # 15% for testing
  random_seed: 42       # Seed for reproducibility
```

#### F) Weights & Biases

Configuration for **tracking and visualization**:

```yaml
wandb:
  use_wandb: false                          # Enable W&B
  project_name: "neural-ancestry-predictor" # Project name
  run_name: null                            # Run name (auto if null)
  log_frequency: 10                         # Log every N batches
```

#### G) Checkpointing

Controls **model saving**:

```yaml
checkpointing:
  checkpoint_dir: "models"       # Directory for checkpoints
  save_frequency: 10             # Save every N epochs
  load_checkpoint: null          # Path to existing checkpoint
```

#### H) Mode

Defines **operation mode**:

```yaml
mode: "train"  # "train" or "test"
```

---

## Architecture

### Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    NEURAL ANCESTRY PREDICTOR                 │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Input: AlphaGenome Predictions                             │
│  ┌──────────────────────────────────────────┐               │
│  │ Window 1 (SNP/Gene)                      │               │
│  │  ├─ H1: [ATAC: 100 bases]                │               │
│  │  └─ H2: [ATAC: 100 bases]                │               │
│  │ Window 2                                  │               │
│  │  ├─ H1: [ATAC: 100 bases]                │               │
│  │  └─ H2: [ATAC: 100 bases]                │               │
│  │ ...                                       │               │
│  │ Window 55                                 │               │
│  │  ├─ H1: [ATAC: 100 bases]                │               │
│  │  └─ H2: [ATAC: 100 bases]                │               │
│  └──────────────────────────────────────────┘               │
│           ↓ Concatenation + Normalization                   │
│  ┌──────────────────────────────────────────┐               │
│  │ Feature Vector [11000 elements]          │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  ┌──────────────────────────────────────────┐               │
│  │ Dense Layer (128 neurons)                │               │
│  │ ReLU Activation                          │               │
│  │ Dropout (0.2)                            │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  ┌──────────────────────────────────────────┐               │
│  │ Dense Layer (64 neurons)                 │               │
│  │ ReLU Activation                          │               │
│  │ Dropout (0.2)                            │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  ┌──────────────────────────────────────────┐               │
│  │ Output Layer (5 neurons)                 │               │
│  │ Softmax                                  │               │
│  └──────────────────────────────────────────┘               │
│           ↓                                                  │
│  Output: [AFR, AMR, EAS, EUR, SAS] probabilities           │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Components

1. **ProcessedGenomicDataset**: Dataset wrapper that:
   - Loads data from `GenomicLongevityDataset`
   - Extracts central region from windows
   - Applies downsampling
   - Combines haplotypes
   - Normalizes (z-score)
   - Caches normalization parameters

2. **AncestryPredictor**: PyTorch model that:
   - Dynamic construction based on config
   - Fully connected layers (Dense)
   - Dropout for regularization
   - Softmax for classification or linear for regression

3. **Trainer**: Manages training:
   - Training loop with progress bars
   - Periodic validation
   - Automatic checkpointing
   - W&B logging

4. **Tester**: Manages testing:
   - Inference on test set
   - Detailed metrics
   - Confusion matrix
   - Classification report

---

## Data Processing

### Processing Pipeline

```
Original Dataset → Center     → Downsampling → Haplotype  → Normalization → Tensor
                   Extraction                   Combination
```

#### 1. Central Region Extraction

Each window has ~1M bases. We extract the central region:

```python
# window_center_size = 100
# Original array: [1000000 elements]
center_idx = 500000
start = center_idx - 50  # 499950
end = center_idx + 50    # 500050
extracted = array[start:end]  # [100 elements]
```

**Why center?** Assumes the central region is more relevant (closer to gene/SNP).

#### 2. Downsampling

Further reduces dimensionality:

```python
# downsample_factor = 2
downsampled = extracted[::2]  # [50 elements]
```

#### 3. Haplotype Combination

```python
# haplotype_mode = "H1+H2"
features = concatenate([H1_features, H2_features])

# haplotype_mode = "H1"
features = H1_features
```

#### 4. Normalization

Z-score normalization using the entire dataset:

```python
# Computed once at the beginning
mean = mean(all_training_data)
std = std(all_training_data)

# Applied to each sample
normalized = (features - mean) / std
```

Parameters saved in `models/normalization_params.json` for reuse.

---

## Processed Dataset Cache

### Overview

The Neural Ancestry Predictor supports **caching of processed datasets** to dramatically speed up subsequent runs. The cache stores:
- Normalized feature vectors (z-score normalized)
- Train/validation/test splits
- Normalization parameters
- Metadata for validation

### Enabling Cache

Edit `configs/default.yaml`:

```yaml
dataset_input:
  processed_cache_dir: "processed_datasets/my_experiment"  # Enable cache
```

### Cache Workflow

**First Execution (cache miss):**
```
1. Load raw dataset                    [5s]
2. Process windows (extract + downsample) [10s]
3. Compute normalization parameters    [30s]
4. Split dataset                       [1s]
5. Save to cache                       [10s]
6. Start training                      ---
   Total preprocessing: ~56s
```

**Subsequent Executions (cache hit):**
```
1. Validate cache compatibility        [1s]
2. Load processed data from cache      [3s]
3. Start training                      ---
   Total preprocessing: ~4s
   
   ⚡ 14x faster!
```

### Cache Validation

The cache is automatically validated and rebuilt if any of these parameters change:
- `dataset_dir`
- `alphagenome_outputs`
- `haplotype_mode`
- `window_center_size`
- `downsample_factor`
- `random_seed` (for splits)

Example validation output:
```
[yellow]Parâmetro window_center_size diferente: cache=100, atual=200[/yellow]
[bold yellow]Cache Inválido ou Incompatível[/bold yellow]
Parâmetros mudaram. Reprocessando dataset...
```

### Cache Structure

```
processed_datasets/my_experiment/
├── metadata.json              # Creation date, config hash, parameters
├── normalization_params.json  # Mean and std for z-score
├── splits.json                # Train/val/test indices
├── train_data.pt              # Processed training data (~20 MB)
├── val_data.pt                # Processed validation data (~5 MB)
└── test_data.pt               # Processed test data (~5 MB)
```

### Use Cases

**1. Hyperparameter Tuning**

Enable cache once, then tune model hyperparameters without reprocessing data:

```yaml
# First run
dataset_input:
  processed_cache_dir: "processed_datasets/base_cache"
model:
  hidden_layers: [128, 64]

# Subsequent runs (instant data loading)
model:
  hidden_layers: [256, 128, 64]  # Changed
  dropout_rate: 0.3              # Changed
```

**2. Multiple Experiments**

Use different cache directories for different data configurations:

```yaml
# Experiment 1: ATAC only
dataset_input:
  alphagenome_outputs: ["ATAC"]
  processed_cache_dir: "processed_datasets/atac_only"

# Experiment 2: ATAC + RNA_SEQ
dataset_input:
  alphagenome_outputs: ["ATAC", "RNA_SEQ"]
  processed_cache_dir: "processed_datasets/atac_rna"
```

**3. Debugging**

Keep cache enabled during debugging to focus on model code:

```yaml
dataset_input:
  processed_cache_dir: "processed_datasets/debug_cache"
```

### Disabling Cache

Set to `null` to always reprocess:

```yaml
dataset_input:
  processed_cache_dir: null  # Cache disabled
```

Or delete the cache directory:

```bash
rm -rf processed_datasets/my_experiment
```

### Performance Comparison

| Dataset Size | First Run | Cached Run | Speedup |
|--------------|-----------|------------|---------|
| 78 samples   | ~56 sec   | ~4 sec     | 14x     |
| 200 samples  | ~2 min    | ~8 sec     | 15x     |
| 500 samples  | ~5 min    | ~15 sec    | 20x     |

*Times measured on standard workstation with 78 samples, 55 windows, ATAC output*

### Disk Space

| Dataset Size | Cache Size |
|--------------|------------|
| 78 samples   | ~30 MB     |
| 200 samples  | ~80 MB     |
| 500 samples  | ~200 MB    |

### Troubleshooting

**Q: Cache not loading?**

Check validation output. If parameters changed, cache is auto-rebuilt.

**Q: Want to force reprocessing?**

```bash
# Option 1: Delete cache
rm -rf processed_datasets/my_experiment

# Option 2: Disable in YAML
processed_cache_dir: null
```

**Q: Different experiments sharing cache?**

Use unique cache directories per experiment:
```yaml
processed_cache_dir: "processed_datasets/exp_001"
```

**Q: Out of disk space?**

Clean old caches:
```bash
rm -rf processed_datasets/old_*
```

---

## Dataset Verification

### Overview

The `verify_processed_dataset.py` tool is essential for **validating the quality and consistency** of processed genomic datasets before training. It compares data from different stages of the pipeline to detect bugs, inconsistencies, or processing errors.

### What does it do?

- 🔍 **Compares datasets** from different processing stages
- 📊 **Visualizes RNA-seq tracks** for individual genes and samples
- 📈 **Computes metrics** (MAE, correlation) to quantify differences
- 🎮 **Interactive navigation** to explore multiple samples
- 👥 **Compare individuals** side-by-side to see genetic variations

### Comparison Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `alphagenome_ref_x_dataset_dir` | AlphaGenome (reference) vs processed dataset | Validate reference genome processing |
| `alphagenome_ind_x_dataset_dir` | AlphaGenome (individual) vs processed dataset | Validate individual variant processing |
| `dataset_dir_x_cache_dir` | Original dataset vs processed cache | Validate normalization and cache |
| `alphagenome_x_alphagenome_ref` | API methods comparison | Validate AlphaGenome API consistency |

### Interactive Modes

**Single Mode** (default):
- Navigate through samples one at a time
- Use ← → to move between individuals
- View 6 RNA-seq tracks per gene (3 cell types × 2 strands)

**Comparison Mode**:
- Compare two individuals simultaneously
- Navigate both together (← →) or independently (A/D for second individual)
- Switch genes with W/Z keys
- Visualize genetic variations between individuals

### Quick Start

```bash
cd neural_ancestry_predictor

# Verify cache vs dataset
python3 verify_processed_dataset.py --config configs/verify_processed_dataset.yaml

# Verify specific gene
python3 verify_processed_dataset.py --config configs/verify_tyr_only.yaml

# Compare two individuals interactively
# Edit config: interactive_comparison_mode: "comparison"
python3 verify_processed_dataset.py --config configs/verify_comparison.yaml
```

### When to Use

- ✅ **Before training**: Ensure processed data is correct
- ✅ **After dataset creation**: Validate `build_non_longevous_dataset` output
- ✅ **After normalization**: Check if cache matches original data
- ✅ **When debugging**: Visualize what the network actually sees
- ✅ **Comparing individuals**: Understand genetic variations

### Key Features

- **11 ancestry-related genes**: SLC24A5, SLC45A2, OCA2, HERC2, MC1R, EDAR, MFSD12, DDB1, TCHH, TYR, TYRP1
- **3 cell type ontologies**: Melanocyte (CL:1000458), Dermal Papilla (CL:0000346), Keratinocyte (CL:2000092)
- **Strand-specific**: Tracks for both + and - DNA strands
- **Metrics**: MAE (Mean Absolute Error) and Pearson correlation
- **Flexible filtering**: View all genes or specific ones

### Example Output

```
════════════════════════════════════════════════════════════
       VERIFICAÇÃO DE DATASET PROCESSADO                   
════════════════════════════════════════════════════════════

✓ Sample: HG00120 (índice 0, global 0)

═══════════════════════════════════════════════════════
           MÉTRICAS DE COMPARAÇÃO                      
═══════════════════════════════════════════════════════
┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━┓
┃ Métrica               ┃             Valor ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━┩
│ MAE Média (global)    │          0.004782 │
│ Correlação Média      │          0.987512 │
└───────────────────────┴───────────────────┘
```

### Interpreting Results

**MAE (Mean Absolute Error):**
- < 0.01: ✅ Excellent match
- 0.01-0.05: ⚠️ Good match, minor differences
- \> 0.05: ❌ Possible pipeline bug

**Pearson Correlation:**
- \> 0.9: ✅ Excellent correlation
- 0.7-0.9: ✅ Good correlation
- < 0.7: ❌ Investigate discrepancies

### Full Documentation

For complete documentation including:
- Configuration options
- All comparison modes
- Troubleshooting guide
- Advanced usage examples

See **[docs/README_VERIFY.md](docs/README_VERIFY.md)**

---

## Training

### Run Training

```bash
python3 neural_ancestry_predictor.py --config configs/default.yaml
```

### During Training

The program will:

1. **Load dataset** and print summary
2. **Compute normalization** (may take a few minutes)
3. **Split data** into train/validation/test
4. **Create model** and print architecture
5. **Train** with progress bars per epoch
6. **Validate** every N epochs
7. **Save checkpoints**:
   - `best_loss.pt`: Best validation loss
   - `best_accuracy.pt`: Best validation accuracy
   - `epoch_N.pt`: Periodic checkpoints

### Monitoring

**Terminal:**
- Progress bars per epoch
- Validation loss and accuracy
- Warnings and errors

**Files:**
- `models/training_history.json`: Complete history
- `models/normalization_params.json`: Normalization parameters

**Weights & Biases** (if enabled):
- Real-time loss graphs
- Validation metrics
- Gradient histograms
- Run comparisons

### Resume Training

To continue from a checkpoint:

```yaml
# configs/default.yaml
checkpointing:
  load_checkpoint: "models/epoch_50.pt"
```

---

## Sklearn Baselines and PCA Cache

### Overview

In addition to `NN`/`CNN`/`CNN2`, the project now supports sklearn baselines:
- `SVM` (LinearSVC, optional probability calibration)
- `RF` (Random Forest)
- `XGBOOST` (if `xgboost` is installed)

These baselines use a preprocessing pipeline with:
1. `StandardScaler` fit on train split
2. `IncrementalPCA` fit on train split
3. Classifier fit on PCA-reduced train matrix
4. Evaluation on train/val/test

### Main config keys

```yaml
model:
  type: "SVM"  # "SVM", "RF", "XGBOOST", "NN", "CNN", "CNN2"
  sklearn:
    pca_components: 300
    use_pca_cache: true
    pca_align_n_train: false
    svm:
      C: 1.0
      max_iter: 20000
      calibrate_probabilities: true
      calibration_cv: 3
```

### PCA cache for repeated experiments

With `model.sklearn.use_pca_cache: true`, reduced matrices are persisted in:
- `processed_cache_dir/pca_cache/<dataset_tag>_pca<requested_k>/`

Saved artifacts include:
- `X_train.npy`, `y_train.npy`, `X_val.npy`, `y_val.npy`, `X_test.npy`, `y_test.npy`
- `scaler_pca.joblib` (fitted scaler + PCA)
- `pca_metadata.json`

Benefits:
- Avoids recomputing PCA when only classifier/hyperparameters change
- Reuses the same cache across SVM/RF/XGBoost runs with matching dataset + PCA settings

### Important PCA details

- `pca_components` can be automatically reduced to a LAPACK-safe value for very large `D` (feature dimension), preventing int32 overflow in large intermediate matrices.
- `pca_align_n_train: true` can reduce effective `k` to a divisor of `n_train`, avoiding tail padding in final `IncrementalPCA.partial_fit`.

### Progress bars

PCA steps now use Rich progress bars when running with console output:
- StandardScaler fit (train)
- IncrementalPCA fit (train)
- PCA transforms (train/val/test in cache build; train in no-cache path)

### Utility scripts added

```bash
# Build/rebuild PCA cache from config
python3 sklearn_pca_cache.py --config configs/genes_1000.yaml
python3 sklearn_pca_cache.py --config configs/genes_1000.yaml --force

# Plot cumulative explained variance for IncrementalPCA
python3 plot_sklearn_pca_variance.py --config configs/genes_1000.yaml --output runs/pca_variance.png

# Plot confusion matrix PNGs from train/val/test_results.json
python3 plot_sklearn_confusion_matrices.py --experiment-dirs /path/to/svm_run /path/to/rf_run /path/to/xgboost_run

# Hyperparameter search on cached PCA matrices (validation metrics + plots)
python3 tune_sklearn_baselines.py --config configs/genes_1000.yaml --output-dir runs/sklearn_tune_001
python3 tune_sklearn_baselines.py --config configs/genes_1000.yaml --output-dir runs/sklearn_tune_custom \\
    --grid-json configs/sklearn_tune_grid.example.json --models SVM,RF

# Large grids: reduce chart clutter
python3 tune_sklearn_baselines.py --config configs/genes_1000.yaml --output-dir runs/sklearn_tune_big \\
    --top-k-bars 8 --label-max-len 32
```

Outputs:
- `tuning_results.csv` / `tuning_results.json` with per-trial metrics
- `tuning_best_per_model_val_metrics.png` (best trial per model, selected by validation F1)
- `tuning_top_trials_val_f1.png` (top-k bars per model; controlled by `--top-k-bars` and `--label-max-len`)
- `tuning_ranked_val_f1_curves.png` (ranked curves, better for large search spaces)
- `tuning_train_vs_val_f1.png` (overfitting diagnostic)
- `tuning_<model>_val_f1_heatmap.png` when exactly two hyperparameters vary for that model
- optional `tuning_<model>_<param>_vs_val_f1.png` when exactly one numeric hyperparameter varies

---

## Debug and Visualization

### Interactive Data Inspection

The module includes a **debug visualization mode** that lets you inspect input data and network predictions sample-by-sample during training. This is invaluable for debugging when the network isn't learning.

### Enable Visualization

```yaml
# configs/default.yaml
debug:
  enable_visualization: true
  max_samples_per_epoch: 10  # Optional: limit to first N samples
```

**Important:** When enabled, `batch_size` is automatically forced to 1.

### What You'll See

For each training sample, a visualization window will display:

**Upper Panel - Input Features:**
- Line plot of all input feature values (x=feature index, y=feature value)
- Feature statistics (min, max, mean, std)
- Total number of features

**Lower Panel - Network Output:**
- Bar chart of predicted probabilities for each class
- **Green bar**: True class (target)
- **Red border**: Predicted class
- Correctness indicator (✓ CORRETO or ✗ ERRADO)
- Class names and probabilities

### Usage

1. **Enable in config:**
   ```yaml
   debug:
     enable_visualization: true
     max_samples_per_epoch: 20  # Visualize only first 20 samples
   ```

2. **Run training:**
   ```bash
   python3 neural_ancestry_predictor.py --config configs/default.yaml
   ```

3. **Inspect each sample:**
   - The visualization window will appear
   - Press **any key** in the graph window to proceed to next sample
   - Analyze patterns in features and predictions

### Use Cases

- **Network not learning:** Check if input features look reasonable (not all zeros/NaN)
- **Always wrong:** See if predictions are stuck on one class
- **Data quality:** Verify normalization is working correctly
- **Class imbalance:** Observe distribution of true labels
- **Feature patterns:** Look for distinguishable patterns between classes

### Example Workflow

```bash
# 1. Enable visualization for debugging
vim configs/default.yaml  # Set enable_visualization: true

# 2. Limit to first 10 samples per epoch
#    (saves time while debugging)
#    max_samples_per_epoch: 10

# 3. Run with small number of epochs
#    num_epochs: 2

# 4. Inspect visualizations and press any key in graph window to advance

# 5. Once issue identified, disable and train normally
#    enable_visualization: false
#    num_epochs: 100
```

---

## Testing and Evaluation

### Run Test

```bash
python3 neural_ancestry_predictor.py --config configs/default.yaml --mode test
```

Or configure in YAML:

```yaml
mode: "test"
checkpointing:
  load_checkpoint: "models/best_accuracy.pt"
```

### Results

The test generates:

**1. Overall Metrics:**

```
╔═══════════════════════════════════════╗
║        Performance Metrics            ║
╠═══════════════════════════════════════╣
║ Accuracy           │ 0.9167           ║
║ Precision (weighted)│ 0.9250          ║
║ Recall (weighted)  │ 0.9167           ║
║ F1-Score (weighted)│ 0.9183           ║
╚═══════════════════════════════════════╝
```

**2. Classification Report:**

```
              precision    recall  f1-score   support

         AFR       0.92      0.95      0.93        20
         AMR       0.88      0.85      0.86        13
         EAS       0.95      0.93      0.94        15
         EUR       0.90      0.93      0.91        15
         SAS       0.93      0.90      0.91        15

    accuracy                           0.92        78
   macro avg       0.92      0.91      0.91        78
weighted avg       0.92      0.92      0.92        78
```

**3. Confusion Matrix:**

```
╔════════════════════════════════════════════════╗
║              Confusion Matrix                  ║
╠════════════════════════════════════════════════╣
║ True \ Pred │  AFR  │  AMR  │  EAS  │  EUR  │  SAS  ║
║ AFR         │   19  │    1  │    0  │    0  │    0  ║
║ AMR         │    1  │   11  │    0  │    1  │    0  ║
║ EAS         │    0  │    0  │   14  │    1  │    0  ║
║ EUR         │    0  │    1  │    0  │   14  │    0  ║
║ SAS         │    0  │    0  │    1  │    0  │   14  ║
╚════════════════════════════════════════════════╝
```

### Interpretation

- **Accuracy**: % of correct predictions
- **Precision**: % of correct positive predictions
- **Recall**: % of positive cases identified
- **F1-Score**: Harmonic mean of precision and recall
- **Confusion Matrix**: Where the model makes mistakes

**Example analysis:**
- AFR: High recall (0.95) → identifies Africans well
- AMR: Lower precision (0.88) → sometimes confuses with others
- Strong diagonal → well-calibrated model

---

## Model Interpretability

### Overview

The Neural Ancestry Predictor includes **DeepLIFT** (Deep Learning Important Features) for understanding which genomic regions contribute most to ancestry predictions. This is essential for:

- 🔬 **Validating** that the model learns biologically meaningful patterns
- 📊 **Identifying** ancestry-informative genomic markers
- 🧬 **Extracting** DNA sequences for further analysis (BLAT/BLAST)
- 📝 **Publishing** interpretable results

### Quick Start

```yaml
# configs/genes_interp.yaml
debug:
  enable_visualization: true
  interpretability:
    enabled: true
    method: "deeplift"
    save_images: true
    output_dir: "interpretability_results"
    deeplift:
      baseline: "mean"      # Use dataset mean as reference
      target_class: "AFR"   # Analyze AFR class patterns

mode: "test"
```

```bash
python3 neural_ancestry_predictor.py --config configs/genes_interp.yaml
```

### Key Features

| Feature | Description |
|---------|-------------|
| **Class Mean Mode** | Average attributions across all samples of a class for robust patterns |
| **Top Regions** | Automatically identifies the 5 most important genomic regions |
| **Individual Search** | Finds the individual with maximum attribution in each region |
| **DNA Extraction** | Extracts 1000bp sequences (H1 + H2) for BLAT/BLAST analysis |
| **Idempotency** | Skips already-generated outputs for efficient re-runs |

### Example Output

The system generates:
1. **Visualization PNG**: Heatmaps of input data and attribution maps
2. **Top Regions Report**: Detailed analysis with DNA sequences

```
Top 5 Regiões Mais Ativas (DeepLIFT):
  1. DDB1: valor = 0.041460, chr11: 61,082,789
  2. OCA2: valor = 0.039910, chr15: 28,000,123
  3. HERC2: valor = 0.032900, chr15: 28,356,789
  ...
```

### Full Documentation

For complete documentation including theoretical background, configuration options, visualization details, and interpretation guide:

📚 **[docs/DEEPLIFT.md](docs/DEEPLIFT.md)**

---

## Variant Annotation (DeepLIFT Post-Processing)

### Overview

After running DeepLIFT interpretation, you can use **annotate_deeplift_windows.py** to analyze the genetic variants in the identified genomic regions. This tool:

- 🔬 **Calls variants** by comparing individual DNA sequences against the hg38 reference genome
- 🧬 **Annotates variants** using Ensembl VEP (consequence, impact, rsIDs)
- 📊 **Fetches gene information** from Ensembl, NCBI, and UniProt APIs
- 📝 **Generates reports** with HIGH impact variant analysis and phenotype validation

### Quick Start

```bash
# Using YAML configuration (recommended)
python3 annotate_deeplift_windows.py --config configs/annotate_deeplift.yaml

# Override YAML parameters from command line
python3 annotate_deeplift_windows.py --config configs/annotate_deeplift.yaml --central-window 50

# Legacy mode (without YAML)
python3 annotate_deeplift_windows.py \
    top_regions_class_mean_AFR_250samples_deeplift.txt \
    --outdir variant_analysis \
    --central-window 50
```

### Key Features

| Feature | Description |
|---------|-------------|
| **YAML Configuration** | All parameters in a single config file for reproducibility |
| **Central Window Filter** | Focus on variants near DeepLIFT center (90-95% reduction) |
| **Haplotype Selection** | Analyze H1, H2, or both |
| **HIGH Impact Analysis** | Detailed info on stop_gained, splice variants, frameshifts |
| **Expression Strand** | Compare with DeepLIFT visualization |
| **Multi-API Integration** | UCSC, VEP, Ensembl, NCBI, UniProt |

### Full Documentation

📚 **[docs/ANNOTATE_DEEPLIFT_WINDOWS.md](docs/ANNOTATE_DEEPLIFT_WINDOWS.md)**

---

## Hypothesis Validation

### Overview

After running the DeepLIFT → VEP pipeline, you can use **validate_pigmentation_hypothesis.py** to systematically validate whether the identified genomic regions correspond to known phenotype-associated genes. This is essential for:

- 🔬 **Validating** the pipeline before applying to new phenotypes (e.g., longevity)
- 📊 **Comparing** results with published GWAS studies
- 🧬 **Analyzing** biological mechanisms of detected variants
- 📝 **Generating** comprehensive validation reports

### Quick Start

```bash
# Basic validation
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2

# With population frequencies from gnomAD
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 --gnomad

# Custom output file
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 -o my_report.md
```

### Three-Level Validation

The script performs validation at three levels:

| Level | What it validates | Key metrics |
|-------|-------------------|-------------|
| **Level 1: Genes** | Are detected genes known pigmentation genes? | % match with literature, Crawford 2017 overlap |
| **Level 2: Variants** | Are known pigmentation rsIDs present? | Known SNPs found, HIGH/MODERATE impact counts |
| **Level 3: Mechanism** | Do variants make biological sense? | stop_gained, splice variants, missense analysis |

### Example Output

```
=== Validation Result ===
  Status: VALIDATED (++)
  Score: 76.0/100
[SAVED] Validation report: top_regions_reports_central2/pigmentation_validation.md
```

### Validation Scores

| Score | Status |
|-------|--------|
| ≥ 80 | STRONGLY VALIDATED (+++) |
| 60-79 | VALIDATED (++) |
| 40-59 | PARTIALLY VALIDATED (+) |
| < 40 | NOT VALIDATED (-) |

### Application to Longevity

Once validated with pigmentation (where results are easy to verify), the same pipeline can be applied to longevity studies:

1. Train with longevity data (longevous vs non-longevous individuals)
2. Run DeepLIFT to identify important genomic regions
3. Annotate variants with `annotate_deeplift_windows.py`
4. Create a `validate_longevity_hypothesis.py` with longevity genes (FOXO3, APOE, CETP, etc.)

### Full Documentation

📚 **[docs/VALIDATE_PIGMENTATION_HYPOTHESIS.md](docs/VALIDATE_PIGMENTATION_HYPOTHESIS.md)**

---

## Weights & Biases

### Configure W&B

```bash
# Install
pip3 install wandb

# Authenticate
wandb login

# Enable in config
```

```yaml
wandb:
  use_wandb: true
  project_name: "neural-ancestry-predictor"
  run_name: "experiment-atac-h1h2-100bases"  # Optional
```

### Available Visualizations

1. **Loss Curves**: Train vs Validation loss
2. **Accuracy**: Evolution of accuracy
3. **Confusion Matrix**: Interactive matrix
4. **Gradients**: Gradient histograms
5. **Parameters**: Weight distributions
6. **System Metrics**: GPU, CPU, RAM

### Compare Experiments

In the W&B dashboard, you can:
- Overlay graphs from multiple runs
- Filter by hyperparameters
- Generate comparison tables
- Export graphs for papers (PNG, SVG, PDF)

---

## Hyperparameter Tuning

### Dimensionality too high?

**Problem**: Training too slow, insufficient memory

**Solutions**:
```yaml
dataset_input:
  window_center_size: 50        # Reduce from 100 to 50
  downsample_factor: 2          # Use 1 in every 2 bases
  haplotype_mode: "H1"          # Use only one haplotype
  alphagenome_outputs: ["ATAC"] # Use only 1 output
```

### Underfitting (high loss)?

**Problem**: Model can't learn patterns

**Solutions**:
```yaml
model:
  hidden_layers: [256, 128, 64]  # More layers/neurons
  dropout_rate: 0.0               # Remove regularization

training:
  num_epochs: 200                 # Train longer
```

### Overfitting (val_loss > train_loss)?

**Problem**: Model memorizes training but doesn't generalize

**Solutions**:
```yaml
model:
  dropout_rate: 0.5               # Increase dropout
  hidden_layers: [64, 32]         # Reduce capacity

training:
  learning_rate: 0.0001           # Lower learning rate
```

### Unstable training (loss oscillates)?

**Problem**: Unstable gradients

**Solutions**:
```yaml
training:
  learning_rate: 0.0001           # Lower LR
  optimizer: "adamw"              # Try another optimizer
  batch_size: 32                  # Larger batch
```

### Slow convergence?

**Problem**: Training takes too long

**Solutions**:
```yaml
training:
  learning_rate: 0.01             # Higher LR (careful!)
  batch_size: 64                  # Larger batch
  optimizer: "adam"               # Adam usually faster than SGD
```

---

## FAQ

### Q: How do I speed up repeated runs?

**A**: Enable **processed dataset cache**!

```yaml
dataset_input:
  processed_cache_dir: "processed_datasets/my_cache"
```

This saves processed data (normalized features + splits) to disk. Subsequent runs load from cache in ~4 seconds instead of reprocessing for ~56 seconds.

**Benefits:**
- 10-20x faster startup
- Perfect for hyperparameter tuning
- Automatic invalidation when data parameters change

See [Processed Dataset Cache](#processed-dataset-cache) section for details.

### Q: How long does training take?

**A**: Depends on:
- **Dataset size**: 78 samples → minutes; 1000 samples → hours
- **Dimensionality**: 100 bases → fast; 10000 bases → slow
- **Hardware**: GPU → 10-100x faster than CPU
- **Epochs**: 100 epochs → proportional

**Estimates** (78 samples, 100 bases, GPU):
- Normalization: ~30 seconds
- Epoch: ~5 seconds
- 100 epochs: ~8 minutes

### Q: Which target should I use?

**A**: Recommendations:
- **Beginner**: `superpopulation` (5 classes, easier)
- **Intermediate**: `population` (26 classes)
- **Advanced**: `frog_likelihood` (regression, 150 values)

### Q: Do I need a GPU?

**A**: Not mandatory, but **highly recommended**:
- CPU: Works, but ~50-100x slower
- GPU: Nvidia with CUDA (RTX 3060 or higher ideal)

To install PyTorch with GPU:
```bash
# CUDA 11.8
pip3 install torch --index-url https://download.pytorch.org/whl/cu118
```

### Q: How to interpret accuracy?

**A**: Depends on baseline:
- **Random guessing** (5 classes): 20%
- **Good model**: 70-80%
- **Excellent model**: 85-95%
- **Perfect**: 100% (beware of overfitting!)

Always compare with both validation AND test.

### Q: Can I use multiple AlphaGenome outputs?

**A**: Yes! Increases dimensionality but may improve performance:

```yaml
dataset_input:
  alphagenome_outputs:
    - "ATAC"
    - "RNA_SEQ"
    - "CAGE"
```

Dimension grows linearly: 1 output → 11k features; 3 outputs → 33k features.

### Q: How to add more genes/SNPs to the dataset?

**A**: Recreate the dataset with `build_non_longevous_dataset`:

```yaml
# build_non_longevous_dataset/configs/custom.yaml
build_window_params:
  mode: "snp"
  snp:
    snp_list_file: "path/to/your_snps.txt"  # Add more SNPs
```

More windows → higher dimensionality → may improve or worsen (curse of dimensionality).

### Q: "CUDA out of memory" error?

**A**: Reduce memory usage:

```yaml
training:
  batch_size: 4          # Reduce batch
  
dataset_input:
  window_center_size: 50  # Reduce dimension
  downsample_factor: 2    # Increase downsampling
```

Or use CPU:
```bash
# Force CPU
export CUDA_VISIBLE_DEVICES=""
python3 neural_ancestry_predictor.py --config configs/default.yaml
```

### Q: How to export graphs for papers?

**A**: 

**Option 1: Weights & Biases**
- In dashboard, click on graph → "Export" → PNG/SVG
- High quality, ideal for publication

**Option 2: Programmatically**
```python
import matplotlib.pyplot as plt
import json

# Load history
with open('models/training_history.json') as f:
    history = json.load(f)

# Plot
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(history['epoch'], history['train_loss'], label='Train')
plt.plot(history['epoch'], history['val_loss'], label='Validation')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.savefig('loss_curve.png', dpi=300, bbox_inches='tight')
```

---

## File Structure

```
neural_ancestry_predictor/
├── neural_ancestry_predictor.py      # Main training/inference program
├── sklearn_pca_cache.py              # Shared StandardScaler+IncrementalPCA disk cache for sklearn baselines
├── plot_sklearn_pca_variance.py      # Cumulative explained variance plot (IncrementalPCA)
├── plot_sklearn_confusion_matrices.py # Confusion matrix plotting from sklearn result JSONs
├── tune_sklearn_baselines.py          # Hyperparameter tuning on cached PCA (val metrics + plots)
├── annotate_deeplift_windows.py      # DeepLIFT output variant annotation
├── validate_pigmentation_hypothesis.py  # Hypothesis validation script
├── verify_processed_dataset.py       # Dataset verification tool
├── configs/
│   ├── default.yaml                  # Default configuration
│   ├── genes_interp.yaml             # Interpretability configuration
│   └── ...
├── docs/
│   ├── DEEPLIFT.md                   # DeepLIFT documentation
│   ├── ANNOTATE_DEEPLIFT_WINDOWS.md  # Variant annotation documentation
│   ├── VALIDATE_PIGMENTATION_HYPOTHESIS.md  # Hypothesis validation documentation
│   └── ...
├── models/                           # Checkpoints (created automatically)
│   ├── best_loss.pt
│   ├── best_accuracy.pt
│   ├── normalization_params.json
│   └── training_history.json
└── README.md                         # This documentation
```

---

## References

- **PyTorch**: https://pytorch.org/
- **Weights & Biases**: https://wandb.ai/
- **AlphaGenome**: https://alphagenome.ai/
- **1000 Genomes**: http://www.internationalgenome.org/
- **FROGAncestryCalc**: Ancestry inference via AISNPs

---

## Support

For issues or questions:

1. Check this README
2. Consult [build_non_longevous_dataset/docs/PYTORCH_DATASET.md](../build_non_longevous_dataset/docs/PYTORCH_DATASET.md)
3. Run with debug mode: add prints to the code
4. Check W&B logs (if enabled)

**Author**: Alberto F. De Souza (via ChatGPT)  
**Date**: 2025-11-14  
**Version**: 1.0

---

## Changelog

### v1.3 (2026-03-22)
- ✨ **NEW: sklearn baselines** (`SVM`, `RF`, `XGBOOST`) integrated in main train/test flow
- ✨ **NEW: shared PCA cache** (`sklearn_pca_cache.py`) with metadata validation and force rebuild option
- ✨ **NEW: robust IncrementalPCA guardrails**
  - automatic LAPACK-safe cap for very high-dimensional inputs
  - optional `pca_align_n_train` to avoid tail padding
- ✨ **NEW: progress bars for PCA/scaler stages** in cache and no-cache paths
- ✨ **NEW: diagnostics utilities**
  - `plot_sklearn_pca_variance.py` for cumulative explained variance
  - `plot_sklearn_confusion_matrices.py` for train/val/test confusion-matrix PNGs
  - `tune_sklearn_baselines.py` for grid search on cached PCA with validation metric plots
- ✨ **NEW: large-grid friendly tuning plots**
  - ranked validation curves (`tuning_ranked_val_f1_curves.png`)
  - automatic 2-parameter validation heatmaps per model
  - anti-clutter knobs: `--top-k-bars`, `--label-max-len`
- 🔧 Matplotlib backend handling improved for headless execution (`TkAgg` fallback to `Agg`)

### v1.2 (2025-12-23)
- ✨ **NEW: Hypothesis Validation Script** (`validate_pigmentation_hypothesis.py`)
- ✨ Three-level validation (genes, variants, mechanism)
- ✨ Built-in database of 14 known pigmentation genes with GWAS references
- ✨ Comparison with Crawford et al., 2017 (Science)
- ✨ Optional gnomAD population frequency integration
- ✨ Automatic validation score calculation (0-100)
- 📚 Comprehensive documentation

### v1.1 (2025-11-14)
- ✨ **NEW: Processed dataset cache** for 10-20x faster repeated runs
- ✨ Automatic cache validation and invalidation
- 📚 Comprehensive cache documentation with examples

### v1.0 (2025-11-14)
- ✨ Initial implementation
- ✨ Support for superpopulation, population and FROG likelihood
- ✨ Weights & Biases integration
- ✨ Configurable window and haplotype processing
- ✨ Automatic checkpointing and normalization
- ✨ Detailed metrics and visualizations
