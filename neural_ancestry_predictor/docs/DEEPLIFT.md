# DeepLIFT Interpretability

> **Understanding What Your Model Has Learned**

This document describes the DeepLIFT (Deep Learning Important Features) implementation in the Neural Ancestry Predictor, which provides interpretable attribution maps to understand which genomic regions contribute most to ancestry predictions.

## Table of Contents

- [Overview](#overview)
- [Theoretical Background](#theoretical-background)
- [Configuration](#configuration)
- [Operation Modes](#operation-modes)
- [Top Regions Analysis](#top-regions-analysis)
- [DNA Sequence Extraction](#dna-sequence-extraction)
- [Visualization](#visualization)
- [Output Files](#output-files)
- [Idempotency](#idempotency)
- [Example Workflow](#example-workflow)
- [Interpreting Results](#interpreting-results)

---

## Overview

DeepLIFT is a feature attribution method that explains neural network predictions by comparing neuron activations to a reference (baseline) and propagating contribution scores back to the input layer. This implementation:

- ✅ Computes signed attribution scores for every input element
- ✅ Supports both individual sample and class-mean analysis
- ✅ Identifies top contributing genomic regions per gene
- ✅ Finds individuals with maximum attribution in each region
- ✅ Extracts DNA sequences (H1 and H2 haplotypes) for further analysis
- ✅ Generates publication-ready visualizations

---

## Theoretical Background

### Rescale Rule Approximation

This implementation uses the **Rescale Rule** approximation of DeepLIFT, which is mathematically equivalent to:

```
attribution = gradient × (input - baseline)
```

Where:
- **gradient**: Backpropagated gradients from the target class output
- **input**: The actual input being analyzed
- **baseline**: A reference input (zeros or dataset mean)

This approximation is valid for networks with ReLU activations and provides the same results as the full DeepLIFT algorithm with significantly simpler implementation.

### Reference

Shrikumar et al., 2017 - *"Learning Important Features Through Propagating Activation Differences"*

---

## Configuration

DeepLIFT is configured in the `debug.interpretability` section of your YAML config:

```yaml
debug:
  enable_visualization: true
  visualization:
    width: 600
    height: 660
    downsample_aggregation: "max"  # "max" or "mean" for visualization
  interpretability:
    enabled: true
    method: "deeplift"             # "deeplift" or "gradcam"
    save_images: true
    output_dir: "interpretability_results"
    deeplift:
      baseline: "mean"             # "zeros" or "mean"
      target_class: "AFR"          # Class name or "predicted"

mode: "test"
test_dataset: "train"              # Which split to analyze
```

### Configuration Parameters

| Parameter | Options | Description |
|-----------|---------|-------------|
| `method` | `"deeplift"`, `"gradcam"` | Interpretability method to use |
| `save_images` | `true`, `false` | Save visualization PNGs and reports |
| `output_dir` | string | Directory for output files |
| `baseline` | `"zeros"`, `"mean"` | Reference input for DeepLIFT |
| `target_class` | class name or `"predicted"` | Which class to compute attributions for |

### Baseline Options

#### `baseline: "zeros"`
- Uses a zero tensor as reference
- Fast computation (no preprocessing needed)
- Good for understanding what makes the input different from "nothing"

#### `baseline: "mean"`
- Uses the mean of up to 10,000 training samples as reference
- Computed once and cached for efficiency
- **Recommended**: Shows what makes the input different from the "average" input
- Prevents data leakage (only uses training set)

---

## Operation Modes

### Individual Mode (`target_class: "predicted"`)

Computes attributions for each sample individually, using the predicted class as the target:

```yaml
deeplift:
  target_class: "predicted"
```

- Generates one visualization per sample
- Shows which regions contributed to that specific prediction
- Useful for debugging individual predictions

### Class Mean Mode (`target_class: "AFR"`, etc.)

Computes the mean attribution across all samples of the specified class:

```yaml
deeplift:
  target_class: "AFR"  # Or "AMR", "EAS", "EUR", "SAS"
```

- Averages attributions over all samples with that target label
- Shows **consistent** patterns that the model uses for that class
- More robust to noise than individual samples
- **Recommended for publication and analysis**

When in class mean mode:
1. The input visualization shows the mean input for that class
2. The DeepLIFT map shows the mean attribution
3. The title indicates the number of samples averaged

---

## Top Regions Analysis

After computing attributions, the system identifies the **Top 5 Most Active Regions**:

### How Regions are Identified

1. For each gene, find the maximum positive attribution value across all 6 tracks
2. Record the column index (position within the gene window)
3. Rank genes by their peak attribution value
4. Select the top 5 genes

### Individual Search

For each top region, the system finds the **individual with the highest attribution** in that specific position:

1. Iterates through all samples of the target class
2. Computes individual DeepLIFT attributions
3. Extracts the value at the specific (gene, column) position
4. Records the sample with maximum value

The output includes:
- Sample ID
- Superpopulation and population
- Attribution value
- Values at all 5 top regions for comparison

---

## DNA Sequence Extraction

For each individual identified in the top regions analysis, the system extracts **1000 base pairs** of DNA sequence centered on the peak attribution position.

### Extraction Process

1. Locates the FASTA file: `{dataset_dir}/individuals/{sample_id}/windows/{gene}/{sample_id}.H1.window.fixed.fa`
2. Reads the complete sequence
3. Extracts 500 bases before and 500 bases after the center
4. Repeats for H2 haplotype

### Purpose

The extracted sequences can be used for:
- **BLAT/BLAST searches** to identify proteins
- **Variant analysis** to find mutations
- **Motif discovery** to identify regulatory elements
- **Comparison** between individuals with different ancestries

---

## Visualization

### Plot Components

When running in test mode with visualization enabled, you'll see a multi-panel figure:

#### Panel 1: Input Data
- Grayscale heatmap of the (normalized) input
- Y-axis: Genes (11 genes × 6 tracks = 66 rows)
- X-axis: Gene position (0 to window_center_size)

#### Panel 2: DeepLIFT Attribution Map
- Diverging colormap: Blue (negative) ← Black (zero) → Red (positive)
- Positive values indicate features that **support** the target class
- Negative values indicate features that **oppose** the target class
- Green circles mark the top 5 most active regions
- X-axis shows top regions with their attribution values

#### Panel 3: Output Probabilities
- Bar chart of predicted class probabilities
- Green bar: True class
- Red border: Predicted class

### Colormap

The custom colormap uses:
- **Blue (0.2, 0.4, 0.9)**: Maximum negative attribution
- **Black (0.0, 0.0, 0.0)**: Zero attribution
- **Red (0.9, 0.3, 0.2)**: Maximum positive attribution

### Scale

- **Individual mode**: Fixed scale [-0.10, 0.10] for comparison across samples
- **Class mean mode**: Dynamic scale based on data range

---

## Output Files

When `save_images: true`, the following files are generated:

### Visualization PNG

**Class mean mode:**
```
{output_dir}/class_mean_{class}_{N}samples_deeplift.png
```

**Individual mode:**
```
{output_dir}/{sample_id}_{predicted_class}_{correct|wrong}_deeplift.png
```

### Top Regions Report

**Filename:**
```
{output_dir}/top_regions_class_mean_{class}_{N}samples_deeplift.txt
```

**Format:**
```
Top 5 Regiões Mais Ativas (DeepLIFT)
================================================================================
Classe: AFR (250 amostras)
================================================================================

1. DDB1: valor_medio = 0.041460, chr11: 61,082,789
   Indivíduo com maior valor: HG02635 (AFR/YRI, valor = 0.089234)
   Valores nas 5 regiões: DDB1=0.08923, OCA2=0.04512, HERC2=0.03245, TYR=0.02834, SLC24A5=0.02156
   DNA H1 (1000bp centradas em chr11:61,082,789):
   >HG02635_H1_DDB1_center_61082789
   ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
   ...
   DNA H2 (1000bp centradas em chr11:61,082,789):
   >HG02635_H2_DDB1_center_61082789
   ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
   ...

2. OCA2: valor_medio = 0.039910, chr15: 28,000,123
   ...
```

---

## Idempotency

The visualization system is **idempotent**: if output files already exist, they are skipped:

```
[dim]⏭ Skipping visualization (already exists): .../class_mean_AFR_250samples_deeplift.png[/dim]
```

This allows you to:
- Re-run the analysis without recomputing everything
- Add new samples without regenerating existing results
- Resume interrupted runs

To force regeneration, delete the existing output files.

---

## Example Workflow

### 1. Configure for Class Mean Analysis

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
      baseline: "mean"
      target_class: "AFR"

mode: "test"
test_dataset: "train"
```

### 2. Run Analysis

```bash
python3 neural_ancestry_predictor.py --config configs/genes_interp.yaml
```

### 3. Review Console Output

```
Top 5 Regiões Mais Ativas (DeepLIFT):
  1. DDB1: valor = 0.041460, chr11: 61,082,789
  2. OCA2: valor = 0.039910, chr15: 28,000,123
  3. HERC2: valor = 0.032900, chr15: 28,356,789
  4. TYR: valor = 0.030990, chr11: 89,178,456
  5. SLC24A5: valor = 0.023640, chr15: 48,426,123

Buscando indivíduos com maior DeepLIFT em cada região...
  → 250 amostras da classe alvo
  ✓ DDB1: HG02635 (AFR/YRI, valor = 0.089234)
  ✓ OCA2: NA19143 (AFR/ESN, valor = 0.078562)
  ...
```

### 4. Examine Output Files

- View `interpretability_results/class_mean_AFR_250samples_deeplift.png`
- Analyze `interpretability_results/top_regions_class_mean_AFR_250samples_deeplift.txt`
- Use extracted DNA sequences for BLAT/BLAST analysis

---

## Interpreting Results

### What Do Positive Attributions Mean?

Positive attribution values indicate input features that **increase** the model's confidence in the target class. For ancestry prediction:

- High positive attribution in a gene region suggests that gene contains **ancestry-informative variation**
- The pattern of attribution across tracks (cell types) shows **which biological processes** are most informative

### Known Pigmentation Genes

The model typically highlights genes known to be associated with pigmentation and ancestry:

| Gene | Chromosome | Known Association |
|------|------------|-------------------|
| SLC24A5 | chr15 | Skin pigmentation, European ancestry |
| SLC45A2 | chr5 | Skin/hair pigmentation |
| OCA2 | chr15 | Eye/skin color, oculocutaneous albinism |
| HERC2 | chr15 | Blue/brown eye color |
| MC1R | chr16 | Red hair, fair skin |
| TYR | chr11 | Tyrosinase, melanin synthesis |
| TYRP1 | chr9 | Tyrosinase-related protein |

### Validation

To validate that the model has learned biologically meaningful patterns:

1. Compare top regions with known GWAS hits for pigmentation/ancestry
2. Check if identified variants are in functional regions
3. Verify that different ancestries show distinct patterns
4. Cross-reference with published ancestry-informative markers (AIMs)

---

## Troubleshooting

### All attributions are zero

- Check that the model is trained and loaded correctly
- Verify that `enable_visualization: true` is set
- Ensure `mode: "test"` is configured

### Class not found

If you specify a target class that doesn't exist:
```
[yellow]⚠ Classe 'XYZ' não encontrada, usando predita[/yellow]
```
Use one of: `AFR`, `AMR`, `EAS`, `EUR`, `SAS` (for superpopulation)

### DNA extraction fails

```
[yellow]⚠ Arquivo FASTA não encontrado: .../HG00096.H1.window.fixed.fa[/yellow]
```
Verify that `dataset_dir` points to a complete dataset with individual FASTA files.

### Memory issues

For large datasets, the individual search can be memory-intensive. Consider:
- Reducing the dataset size
- Using a subset of samples
- Running on a machine with more RAM

---

## References

- Shrikumar, A., Greenside, P., & Kundaje, A. (2017). Learning Important Features Through Propagating Activation Differences. *ICML*.
- Ancona, M., Ceolini, E., Öztireli, C., & Gross, M. (2018). Towards better understanding of gradient-based attribution methods for Deep Neural Networks. *ICLR*.

---

**Author**: Neural Ancestry Predictor Team  
**Last Updated**: 2025-12-22  
**Version**: 1.0

