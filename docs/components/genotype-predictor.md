# Genotype Predictor

The genotype predictor trains dense/aligned models from processed genotype and AlphaGenome-derived tensors.

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

## BCFtools Chain Alignment

The `alignment/bcftools_chain_mapper.py` path rebuilds haplotype consensus sequences and reads BCFtools chain files to map AlphaGenome predictions back to a global aligned coordinate system. `sync-bcftools-artifacts` copies or links the required consensus and chain artifacts into the canonical dataset tree.

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

- [Research Pipelines](../RESEARCH_PIPELINES.md)
- [Important Experiments](../IMPORTANT_EXPERIMENTS.md)
