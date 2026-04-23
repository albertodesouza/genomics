## neural_prediction

`neural_prediction` is a refactored replacement for `neural_ancestry_predictor`.

Current scope:

- local Hugging Face dataset backend only
- classification tasks with categorical or derived targets
- family-aware train/val/test split
- processed tensor dataset for AlphaGenome predictions
- neural models with a clean training loop (`MLP`, `CNN`, `CNN2`)
- typed config objects
- optional checkpoint save/load
- `train` and `test` modes
- optional Weights & Biases logging
- resumable single-gene screening workflow
- screening result plotting and run verification utilities
- interpretability helpers for `CNN`/`CNN2` and general attribution via DeepLIFT-style gradients
- sklearn baselines with shared PCA cache support
- processed-dataset verification CLI
- confusion-matrix and PCA-variance plotting helpers
- sample visualization utility
- interrupted checkpoint save on Ctrl-C

This module intentionally starts smaller than the legacy codebase so it can stay
maintainable while the new architecture stabilizes.

### Train

```bash
python3 -m neural_prediction.main \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --output-dir /tmp/neural_prediction_train
```

### Test

```bash
python3 -m neural_prediction.main \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --output-dir /tmp/neural_prediction_test
```

For test mode, set these fields in the config first:

```yaml
mode: "test"
checkpointing:
  load_checkpoint: "/path/to/models/best_accuracy.pt"
```

### Verify Processed Dataset

```bash
python3 -m neural_prediction.verify_processed_dataset \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml
```

### Single-Gene Screening

```bash
python3 -m neural_prediction.run_single_gene_screen \
  --base-config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --output-root /tmp/neural_prediction_single_gene
```

Plot screening results:

```bash
python3 -m neural_prediction.plot_single_gene_screen_results \
  --results-csv /tmp/neural_prediction_single_gene/results.csv \
  --metric test_accuracy \
  --output /tmp/neural_prediction_single_gene/test_accuracy.png
```

Verify screening runs:

```bash
python3 -m neural_prediction.verify_runs \
  --results-csv /tmp/neural_prediction_single_gene/results.csv
```

### Attribution

GradCAM:

```bash
python3 -m neural_prediction.run_attribution \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --checkpoint /path/to/models/best_accuracy.pt \
  --method gradcam \
  --sample-index 0 \
  --output /tmp/gradcam_sample0.json
```

DeepLIFT-style attribution:

```bash
python3 -m neural_prediction.run_attribution \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --checkpoint /path/to/models/best_accuracy.pt \
  --method deeplift \
  --baseline mean \
  --sample-index 0 \
  --output /tmp/deeplift_sample0.json
```

### Sklearn Baselines

```bash
python3 -m neural_prediction.run_sklearn_baselines \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --output-dir /tmp/neural_prediction_sklearn
```

Plot PCA variance from the generated cache:

```bash
python3 -m neural_prediction.plot_pca_variance \
  --pca-cache-dir /dados/GENOMICS_DATA/hf/cache/pca_cache/<cache_name> \
  --output /tmp/pca_variance.png
```

### Plot Confusion Matrix

```bash
python3 -m neural_prediction.plot_confusion_matrix \
  --results /tmp/neural_prediction_train/results.json \
  --output /tmp/neural_prediction_train/test_confusion.png
```

### Visualize One Processed Sample

```bash
python3 -m neural_prediction.visualize_sample \
  --config neural_prediction/configs/pigmentation_binary_hf.yaml \
  --sample-index 0 \
  --output /tmp/sample0.png
```
