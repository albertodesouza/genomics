# Dataset Builders

Dataset builders create derived machine-learning datasets from 1000 Genomes and AlphaGenome outputs.

## Non-Longevous Builder

```bash
genomics dataset-builders non-longevous build --config configs/workflows/non_longevous_dataset/default.yaml
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/workflows/dataset_builders/non_longevous/` |
| Configs | `configs/workflows/non_longevous_dataset/` |

## Related Commands

| Command | Purpose |
|---|---|
| `build` | Run the configured dataset builder |
| `build-window` | Forward arguments to the window builder |
| `visualize` | Visualize generated AlphaGenome outputs |

## Implementation Structure

| Module | Responsibility |
|---|---|
| `build_non_longevous_dataset.py` | Main CLI/config-driven builder orchestration |
| `dataset_builder.py` | Lower-level builder operations and reusable pipeline logic |
| `build_window_and_predict.py` | Builds haplotype windows and runs AlphaGenome prediction for each window |
| `genomic_dataset.py` | PyTorch dataset wrapper for generated outputs |
| `frog_ancestry_parser.py` | Parses FROG ancestry likelihood outputs when available |
| `alphagenome_output_visualization.py` | Visualizes generated AlphaGenome arrays and metadata |

## Build Flow

```text
1000 Genomes metadata
  -> sample filtering by superpopulation/population/sex/config rules
  -> VCF/window selection
  -> haplotype consensus extraction
  -> AlphaGenome prediction
  -> per-individual output directories
  -> PyTorch dataset metadata and splits
```

The builder is designed to be idempotent. Existing per-individual outputs can be reused, and failed/missing samples can be inspected and rerun.

## Dataset Contents

Generated datasets typically contain:

- individual metadata;
- split metadata;
- AlphaGenome output arrays;
- sequence/window metadata;
- target labels such as superpopulation, population, sex, or derived targets;
- optional FROG ancestry likelihood vectors.

## Relationship To Predictors

The genotype predictor and legacy neural ancestry predictor can consume datasets produced by this builder. New experiments should prefer canonical configs under `configs/workflows/non_longevous_dataset/` and predictor configs under `configs/predictors/`.

## Legacy Builder

The longevity dataset builder is kept under `legacy/neural_longevity_dataset/` for historical reproducibility and still depends on old `top3` assumptions.
