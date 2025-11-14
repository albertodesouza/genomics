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
- [Training](#training)
- [Testing and Evaluation](#testing-and-evaluation)
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

---

## Installation

### 1. Dependencies

```bash
# Navigate to directory
cd genomics/neural_ancestry_predictor

# Install Python dependencies
pip install torch numpy pandas pyyaml scikit-learn rich

# Optional: Weights & Biases (for tracking)
pip install wandb
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

See `build_non_longevous_dataset/docs/PYTORCH_DATASET.md` for more information.

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
```

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

## Weights & Biases

### Configure W&B

```bash
# Install
pip install wandb

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
pip install torch --index-url https://download.pytorch.org/whl/cu118
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
├── neural_ancestry_predictor.py    # Main program
├── configs/
│   └── default.yaml                 # Default configuration
├── models/                          # Checkpoints (created automatically)
│   ├── best_loss.pt
│   ├── best_accuracy.pt
│   ├── epoch_10.pt
│   ├── epoch_20.pt
│   ├── normalization_params.json
│   └── training_history.json
└── README.md                        # This documentation
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
2. Consult `build_non_longevous_dataset/docs/PYTORCH_DATASET.md`
3. Run with debug mode: add prints to the code
4. Check W&B logs (if enabled)

**Author**: Alberto F. De Souza (via ChatGPT)  
**Date**: 2025-11-14  
**Version**: 1.0

---

## Changelog

### v1.0 (2025-11-14)
- ✨ Initial implementation
- ✨ Support for superpopulation, population and FROG likelihood
- ✨ Weights & Biases integration
- ✨ Configurable window and haplotype processing
- ✨ Automatic checkpointing and normalization
- ✨ Detailed metrics and visualizations
