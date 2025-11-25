# verify_processed_dataset.py - Processed Dataset Verification

## 📋 Overview

A program to verify and compare processed genomic sequence data against AlphaGenome predictions. Visualizes tracks (genes × RNA-seq ontologies) in overlaid line charts to detect potential bugs in the processing pipeline and validate data consistency.

## 🎯 Objective

Detect inconsistencies and validate comparisons between:
- **Processed cache**: Normalized data stored in `.pt` files for training
- **Original dataset**: Processed data in `.npz` files from individual genomes
- **AlphaGenome (reference)**: Raw API predictions using reference genome
- **AlphaGenome (individual)**: Raw API predictions using individual's genome

## ✨ Features

### 🔄 Comparison Modes

The program supports 4 distinct comparison modes:

1. **`alphagenome_ref_x_dataset_dir`**: 
   - Compares AlphaGenome predictions (reference genome) vs. processed dataset data
   - Useful for validating whether the dataset correctly reflects the reference genome

2. **`alphagenome_ind_x_dataset_dir`**: 
   - Compares AlphaGenome predictions (individual genome with variants) vs. processed dataset data
   - Validates whether individual variant processing is correct
   - Uses `build_window_and_predict.py` as a library

3. **`dataset_dir_x_cache_dir`**: 
   - Compares original dataset data (`.npz`) vs. processed cache (`.pt`)
   - Validates normalization and transformation pipeline

4. **`alphagenome_x_alphagenome_ref`**: 
   - Compares two ways of calling the AlphaGenome API
   - `predict_interval` (without FASTA) vs. `predict_sequence` (with extracted FASTA)
   - Validates API consistency

### 🎮 Interactive Navigation Modes

#### "single" Mode (default)
Traditional navigation through one individual at a time:
- **← (left arrow)**: Go to previous sample
- **→ (right arrow)**: Go to next sample
- **Q**: Exit program

#### "comparison" Mode (new!)
Interactive comparison between two individuals:
- **← → (arrows)**: Navigate both individuals simultaneously forward/backward
- **A**: Navigate only the second individual backward
- **D**: Navigate only the second individual forward
- **W**: Advance to next gene (both individuals)
- **Z**: Go back to previous gene (both individuals)
- **Q**: Exit program

In "comparison" mode:
- Displays 6 tracks stacked vertically (3 ontologies × 2 strands)
- First individual in solid blue
- Second individual in dashed red
- Legend shows: `sample_id (population/superpopulation)`
- Current gene displayed in main title

### 🧬 Gene Filtering
- View all 11 genes (66 tracks)
- View only a specific gene (6 tracks)
- View multiple selected genes

### ⚙️ YAML Configuration
- All settings in YAML file
- Multiple predefined configurations
- Easy customization

## 🚀 Usage

### Basic Usage (interactive mode, all genes)

```bash
cd neural_ancestry_predictor

python3 verify_processed_dataset.py --config configs/verify_processed_dataset.yaml
```

### Verify Only a Specific Gene (e.g., MC1R)

```bash
python3 verify_processed_dataset.py --config configs/verify_tyr_only.yaml
```

### Two-Individual Comparison Mode

```yaml
# configs/verify_comparison.yaml
interactive_mode: true
interactive_comparison_mode: "comparison"  # Enable comparison mode
comparison_mode: "dataset_dir_x_cache_dir"
gene_filter: "MC1R"  # Recommended: one gene at a time
```

```bash
python3 verify_processed_dataset.py --config configs/verify_comparison.yaml
# Use ← → to navigate both
# Use A D to navigate only the second individual
# Use W Z to change genes
```

### Create Custom Configuration

Create a `.yaml` file in the `configs/` directory:

```yaml
# configs/my_verification.yaml
cache_dir: "/path/to/cache"
dataset_dir: "/path/to/dataset"
split: "test"
index: 0
gene_filter: "TYR"  # null for all, "GENE" for one, ["G1", "G2"] for multiple

# Comparison mode
comparison_mode: "alphagenome_ref_x_dataset_dir"

# Navigation mode
interactive_mode: true
interactive_comparison_mode: "single"  # or "comparison"

show_navigation_help: true
verbose_metrics: true
show_stats_in_plot: true
save_plots: false
output_dir: null
```

Then execute:

```bash
python3 verify_processed_dataset.py --config configs/my_verification.yaml
```

## 📝 Available Settings

### YAML Configuration File

```yaml
# ═══════════════════════════════════════════════════════════════
# DIRECTORIES (required)
# ═══════════════════════════════════════════════════════════════

cache_dir: "/path/to/cache"           # Processed cache (.pt, metadata.json)
dataset_dir: "/path/to/dataset"       # Original dataset (.npz)

# ═══════════════════════════════════════════════════════════════
# SAMPLE SELECTION
# ═══════════════════════════════════════════════════════════════

split: "test"                         # train, val, or test
index: 0                              # Initial index (changes with ← →)
sample_id: null                       # Specific ID (optional, overrides split/index)

# ═══════════════════════════════════════════════════════════════
# COMPARISON MODE
# ═══════════════════════════════════════════════════════════════

comparison_mode: "dataset_dir_x_cache_dir"
# Options:
#   - alphagenome_ref_x_dataset_dir: AlphaGenome (ref) vs dataset
#   - alphagenome_ind_x_dataset_dir: AlphaGenome (ind) vs dataset
#   - dataset_dir_x_cache_dir: Dataset vs cache (default)
#   - alphagenome_x_alphagenome_ref: API interval vs sequence

# ═══════════════════════════════════════════════════════════════
# GENE FILTER
# ═══════════════════════════════════════════════════════════════

gene_filter: null                     # Options:
                                      # null: all 11 genes (66 tracks)
                                      # "TYR": only TYR gene (6 tracks)
                                      # ["TYR", "TYRP1"]: multiple genes (12 tracks)

# ═══════════════════════════════════════════════════════════════
# INTERACTIVE NAVIGATION
# ═══════════════════════════════════════════════════════════════

interactive_mode: true                # true: keyboard navigation
                                      # false: view only one sample

interactive_comparison_mode: "single" # "single": one individual at a time (← →)
                                      # "comparison": two individuals (← → A D W Z)

show_navigation_help: true            # Show instructions in plot

# ═══════════════════════════════════════════════════════════════
# ALPHAGENOME API (OPTIONAL)
# ═══════════════════════════════════════════════════════════════

alphagenome_api:
  enabled: false                      # Enable API calls
  api_key: null                       # Uses ALPHAGENOME_API_KEY from environment
  rate_limit_delay: 0.5               # Delay between calls (seconds)
  ontology_terms: ["CL:1000458", "CL:0000346", "CL:2000092"]

# ═══════════════════════════════════════════════════════════════
# RAW MODE (OPTIONAL)
# ═══════════════════════════════════════════════════════════════

raw_mode:
  enabled: false                      # Raw mode (AlphaGenome only, no cache)
  source: "files"                     # "files" (.npz) or "api" (API call)
  window_size_key: "SEQUENCE_LENGTH_16KB"  
  # Options: SEQUENCE_LENGTH_2KB, SEQUENCE_LENGTH_16KB, 
  #         SEQUENCE_LENGTH_100KB (128 KiB), SEQUENCE_LENGTH_500KB (512 KiB),
  #         SEQUENCE_LENGTH_1MB

# ═══════════════════════════════════════════════════════════════
# VISUALIZATION
# ═══════════════════════════════════════════════════════════════

verbose_metrics: true                 # Display detailed metrics in console
show_stats_in_plot: true              # Show statistics in plot

# ═══════════════════════════════════════════════════════════════
# SAVING (optional)
# ═══════════════════════════════════════════════════════════════

save_plots: false                     # Save plots automatically
output_dir: null                      # Directory to save (creates if doesn't exist)
output_prefix: "verify"               # File prefix (e.g., verify_HG00120.png)
```

## 🌐 AlphaGenome API Mode

### What is it?

API mode allows verifying the dataset by comparing cache data directly with **real-time** predictions from the AlphaGenome API, instead of using pre-computed `.npz` files.

### When to use?

- ✅ Verify if original `.npz` files were generated correctly
- ✅ Test with custom sequences (individual variants)
- ✅ Complete end-to-end pipeline validation
- ✅ Debug inconsistencies in processed data

### How to enable?

1. Set up API key:
```bash
export ALPHAGENOME_API_KEY="your_api_key_here"
```

2. Create/edit YAML configuration:
```yaml
# Enable API mode
alphagenome_api:
  enabled: true  # Activate API calls
  api_key: null  # Uses ALPHAGENOME_API_KEY from environment
  rate_limit_delay: 0.5  # Delay between calls (seconds)
  ontology_terms: ["CL:1000458", "CL:0000346", "CL:2000092"]

# Recommended: test one gene at a time
gene_filter: "MC1R"
```

3. Run:
```bash
python3 verify_processed_dataset.py --config configs/verify_api_test.yaml
```

### Requirements

- 📦 `alphagenome` package installed (`pip install alphagenome`)
- 🔑 Valid AlphaGenome API key
- 🌐 Internet connection
- 📁 `.fa` sequence files in `dataset_dir/individuals/{sample}/windows/{gene}/`

### Important

⚠️ **API Quota**: Each gene consumes 1 API call. Use `gene_filter` to save quota!  
⚠️ **Rate limiting**: Configure `rate_limit_delay` appropriately (default: 0.5s)  
⚠️ **Ontologies**: Must match **exactly** those used in original dataset creation

### AlphaGenome Constants

AlphaGenome uses specific constants for window sizes:

| Constant | Size (bp) | Size (KiB) | Usage |
|----------|-----------|------------|-------|
| `SEQUENCE_LENGTH_2KB` | 2048 | 2 KiB | Quick tests |
| `SEQUENCE_LENGTH_16KB` | 16384 | 16 KiB | Small genes |
| `SEQUENCE_LENGTH_100KB` | 131072 | 128 KiB | Medium genes |
| `SEQUENCE_LENGTH_500KB` | 524288 | 512 KiB | Large genes |
| `SEQUENCE_LENGTH_1MB` | 1048576 | 1 MiB | Extended regions |

⚠️ **Note**: Constant names use KB (base 1000), but sizes are in powers of 2 (KiB).

## 📊 Output

### Console ("single" Mode)

```
════════════════════════════════════════════════════════════
       PROCESSED DATASET VERIFICATION                   
════════════════════════════════════════════════════════════

Interactive mode activated:
  • Split: test
  • Total samples: 13
  • Initial index: 0
  • Use ← → to navigate, 'q' to exit

✓ Sample: NA19472 (index 0, global 65)

═══════════════════════════════════════════════════════
           COMPARISON METRICS                      
═══════════════════════════════════════════════════════
┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┓
┃ Metric                ┃              Value ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━┩
│ Mean MAE (global)     │           0.004782 │
│ Maximum MAE           │           0.048390 │
│ Minimum MAE           │           0.000001 │
│ Track with max MAE    │ 45 (gene 7, ont 3) │
│                       │                    │
│ Mean Correlation      │           0.337512 │
│ Minimum Correlation   │          -0.043882 │
│ Maximum Correlation   │           0.988911 │
│ Track with min corr.  │ 55 (gene 9, ont 1) │
└───────────────────────┴────────────────────┘
═══════════════════════════════════════════════════════
```

### Console ("comparison" Mode)

```
════════════════════════════════════════════════════════════
       PROCESSED DATASET VERIFICATION                   
════════════════════════════════════════════════════════════

Interactive comparison mode activated:
  • Split: train
  • Total samples: 54
  • Total genes: 11
  • Use ← → (both), A D (ind2), W Z (genes), Q (exit)

✓ Showing: MC1R
  Individual 1: HG00120 (GBR/EUR)
  Individual 2: HG00120 (GBR/EUR)
```

### Plot ("single" Mode)

The plot shows:
- **Solid blue**: Data from first set (cache, dataset, or AlphaGenome)
- **Dashed red**: Data from second set (AlphaGenome, cache, or reference)
- **Gray dotted lines**: Separation between genes
- **Y-axis labels**: Ontology and cell type in two lines
  - Line 1: `CL:XXXXX (±)` (ontology code and strand)
  - Line 2: Cell type name (e.g., `Melanocyte`)
- **Statistics box**: MAE and correlation per track
- **Navigation instructions**: Bottom center (if interactive mode)

### Plot ("comparison" Mode)

The plot shows:
- **6 subplots stacked vertically** (3 ontologies × 2 strands)
- **Solid blue**: First individual's tracks
- **Dashed red**: Second individual's tracks
- **Main title**: Current gene and information for both individuals
- **Y-axis labels**: 
  - `CL:1000458 (+)` / `Melanocyte`
  - `CL:0000346 (+)` / `Dermal Papilla`
  - `CL:2000092 (+)` / `Keratinocyte`
  - (and (-) versions for negative strand)
- **Legend**: `sample_id (population/superpopulation)` for each individual

### RNA-seq Ontologies

The dataset uses 3 cell type ontologies, each with 2 strands (+/-):

| CL Code | Cell Type | Description |
|---------|-----------|-------------|
| CL:1000458 | Melanocyte | Melanin-producing cells |
| CL:0000346 | Dermal Papilla | Dermal papilla cells |
| CL:2000092 | Keratinocyte | Epidermal cells |

Each gene has **6 tracks**: 3 ontologies × 2 strands (+ and -)

## 🧬 Available Genes

1. **SLC24A5** - Solute carrier family 24 member 5
2. **SLC45A2** - Solute carrier family 45 member 2
3. **OCA2** - OCA2 melanosomal transmembrane protein
4. **HERC2** - HECT and RLD domain containing E3 ubiquitin protein ligase 2
5. **MC1R** - Melanocortin 1 receptor
6. **EDAR** - Ectodysplasin A receptor
7. **MFSD12** - Major facilitator superfamily domain containing 12
8. **DDB1** - Damage specific DNA binding protein 1
9. **TCHH** - Trichohyalin
10. **TYR** - Tyrosinase
11. **TYRP1** - Tyrosinase related protein 1

Each gene has 6 tracks corresponding to 3 ontologies × 2 strands.

## 🔍 Interpreting Results

### MAE (Mean Absolute Error)
- **< 0.01**: ✅ Excellent match
- **0.01 - 0.05**: ⚠️ Good match, minor differences
- **> 0.05**: ❌ Possible pipeline bug

### Pearson Correlation
- **> 0.9**: ✅ Excellent correlation
- **0.7 - 0.9**: ✅ Good correlation
- **0.5 - 0.7**: ⚠️ Moderate correlation
- **< 0.5**: ❌ Weak correlation - investigate

### In the Plot
1. **Overlap**: Blue and red lines should be close
2. **Patterns**: Curve shapes should be similar
3. **Outliers**: Tracks with large divergence need investigation

### "comparison" Mode
- In comparison mode, there are no MAE/correlation metrics
- Focus is on **visual comparison** of tracks from two individuals
- Differences between individuals may indicate real genetic variations

## 📚 Usage Examples

### Example 1: Browse All Test Samples

```yaml
# configs/scan_test_all.yaml
cache_dir: "/dados/GENOMICS_DATA/top3/non_longevous_results_runs_genes/datasets/rna_seq_H1_100000_ds1_log_split0.7-0.15-0.15_seed5"
dataset_dir: "/dados/GENOMICS_DATA/top3/non_longevous_results_genes"
split: "test"
index: 0
gene_filter: null  # All genes
interactive_mode: true
interactive_comparison_mode: "single"
comparison_mode: "dataset_dir_x_cache_dir"
```

```bash
python3 verify_processed_dataset.py --config configs/scan_test_all.yaml
# Press → to navigate through 13 test samples
# Press 'q' to exit
```

### Example 2: Verify Specific Gene Across All Samples

```yaml
# configs/verify_tyr_all.yaml
gene_filter: "TYR"
interactive_mode: true
interactive_comparison_mode: "single"
split: "train"  # 54 samples
comparison_mode: "dataset_dir_x_cache_dir"
```

```bash
python3 verify_processed_dataset.py --config configs/verify_tyr_all.yaml
# Navigate with → to see TYR in all training samples
```

### Example 3: Compare Two Individuals Interactively

```yaml
# configs/compare_individuals.yaml
interactive_mode: true
interactive_comparison_mode: "comparison"
comparison_mode: "dataset_dir_x_cache_dir"
gene_filter: "MC1R"  # Recommended: one gene at a time
split: "train"
```

```bash
python3 verify_processed_dataset.py --config configs/compare_individuals.yaml
# ← → : navigate both individuals
# A D : navigate only the second individual
# W Z : change the displayed gene
# Q   : exit
```

### Example 4: Validate AlphaGenome vs Dataset (reference)

```yaml
# configs/validate_alphagenome_ref.yaml
comparison_mode: "alphagenome_ref_x_dataset_dir"
gene_filter: "MC1R"
interactive_mode: true
interactive_comparison_mode: "single"
sample_id: "HG02445"
```

```bash
python3 verify_processed_dataset.py --config configs/validate_alphagenome_ref.yaml
# Compare AlphaGenome (reference genome) with processed dataset
```

### Example 5: Validate AlphaGenome vs Dataset (individual)

```yaml
# configs/validate_alphagenome_ind.yaml
comparison_mode: "alphagenome_ind_x_dataset_dir"
gene_filter: "MC1R"
interactive_mode: false
sample_id: "HG02445"
```

```bash
python3 verify_processed_dataset.py --config configs/validate_alphagenome_ind.yaml
# Compare AlphaGenome (individual genome) with processed dataset
# Uses build_window_and_predict.py to generate individual predictions
```

### Example 6: Save Plots of All Samples

```yaml
# configs/save_all_plots.yaml
interactive_mode: false
save_plots: true
output_dir: "verification_plots"
output_prefix: "verify"
gene_filter: "TYR"
comparison_mode: "dataset_dir_x_cache_dir"
```

```bash
# Create script to process all samples
for i in {0..12}; do
    sed "s/index: 0/index: $i/" configs/save_all_plots.yaml > /tmp/config_$i.yaml
    python3 verify_processed_dataset.py --config /tmp/config_$i.yaml
done
```

## 🛠️ Recommended Workflow

### 1. Initial Check (Overview)
```bash
# Check a few samples with all genes
python3 verify_processed_dataset.py --config configs/verify_processed_dataset.yaml
# Use → to see 3-5 different samples
# If MAE > 0.05 or Corr < 0.5 → investigate
```

### 2. Gene-Specific Investigation
```bash
# If problem detected, isolate problematic gene
python3 verify_processed_dataset.py --config configs/verify_tyr_only.yaml
# Check if problem is gene-specific or general
```

### 3. AlphaGenome Validation
```bash
# Validate data against AlphaGenome (reference)
# Edit config for comparison_mode: "alphagenome_ref_x_dataset_dir"
python3 verify_processed_dataset.py --config configs/verify_alphagenome_ref.yaml

# Validate data against AlphaGenome (individual)
# Edit config for comparison_mode: "alphagenome_ind_x_dataset_dir"
python3 verify_processed_dataset.py --config configs/verify_alphagenome_ind.yaml
```

### 4. Individual Comparison
```bash
# Compare expression patterns between individuals
# Edit config for interactive_comparison_mode: "comparison"
python3 verify_processed_dataset.py --config configs/compare_individuals.yaml
# Use A D to navigate the second individual independently
# Use W Z to change genes
```

### 5. Systematic Scan
```bash
# Verify all splits
python3 verify_processed_dataset.py --config configs/verify_train.yaml
python3 verify_processed_dataset.py --config configs/verify_val.yaml  
python3 verify_processed_dataset.py --config configs/verify_test.yaml
```

### 6. Documentation
```bash
# Save plots for report
# Configure save_plots: true and output_dir
python3 verify_processed_dataset.py --config configs/save_for_report.yaml
```

## 🐛 Troubleshooting

### Error: "Configuration file does not exist"
- Check YAML file path
- Verify you're in the correct directory

### Error: "Cache dir does not exist"
- Check `cache_dir` in YAML
- Verify read permissions

### Error: "Gene not found"
- Check gene name spelling
- Use complete list: SLC24A5, SLC45A2, OCA2, HERC2, MC1R, EDAR, MFSD12, DDB1, TCHH, TYR, TYRP1

### Error: "Invalid window_size_key: SEQUENCE_LENGTH_512KB"
- AlphaGenome uses `SEQUENCE_LENGTH_500KB` for 524288 bp (512 KiB)
- AlphaGenome uses `SEQUENCE_LENGTH_100KB` for 131072 bp (128 KiB)
- Check correct constants in table above

### Error: "Sequence length X not supported by the model"
- AlphaGenome only supports: 2048, 16384, 131072, 524288, 1048576 bp
- Check if dataset was generated with one of these sizes
- Use corresponding `window_size_key` for dataset size

### MAE too high (> 0.05)
1. Check that cache was generated with same parameters
2. Verify normalization method
3. Check if AlphaGenome data is correct
4. Use `alphagenome_ref_x_dataset_dir` mode to validate

### Window not responding to keys
- Ensure matplotlib window has focus
- On some systems, click on window before pressing keys

### Plot too "crowded"
- Use `gene_filter` to visualize fewer tracks
- Example: `gene_filter: "TYR"` shows only 6 tracks
- In "comparison" mode, use one gene at a time

### alphagenome_ind_x_dataset_dir mode slow
- This mode calls `build_window_and_predict.py` to generate predictions
- Uses temporary files in `/tmp/GENOMICS_DATA/top3`
- Requires access to reference FASTA and VCF
- Automatically cleans up temporary files at the end

## 🔧 Technical Details

### Window Size and Context

The program implements sophisticated logic to handle window sizes:

1. **Automatic detection**: Reads dataset and detects its original size (full_length)
2. **AlphaGenome prediction**: Calls API with same full_length to maintain context
3. **Visualization**: Extracts smaller central window (viz_length) from both for comparison
4. **Centralization**: Uses consistent `extract_center_window()` function

Example:
- Dataset generated with 524288 bp (512 KiB)
- AlphaGenome called with 524288 bp (same context)
- Visualization shows only 16384 bp (center)
- Both extracted from same center → perfect alignment

### Boundary Conditions

Genes near chromosome start/end are handled correctly:

1. **Clipped coordinates**: `samtools faidx` may return smaller sequences
2. **'N' padding**: Sequences are padded at start/end as needed
3. **Alignment maintained**: Relative positions preserved even with padding

### Logic Centralization

`extract_center_window()` function used consistently:
- Removes code duplication
- Ensures same centralization logic
- Fixes off-by-one bug from previous versions

## 📦 Dependencies

- Python 3.10+
- PyTorch
- NumPy
- Matplotlib
- SciPy
- PyYAML
- Rich
- Pandas
- alphagenome (optional, for API mode)
- samtools (optional, for alphagenome_ind mode)
- bcftools (optional, for alphagenome_ind mode)

All already installed in `genomics` environment.

## 📄 Example Files

The repository includes:
- `configs/verify_processed_dataset.yaml` - Default configuration (all genes, single mode)
- `configs/verify_tyr_only.yaml` - MC1R gene only (filter example)
- `configs/verify_raw_test.yaml` - Raw mode (AlphaGenome only)

## 🔄 Version History

| Date | Version | Changes |
|------|---------|---------|
| 2025-11-23 | 1.0 | Initial version with interactive navigation |
| 2025-11-24 | 1.1 | Added gene filter and API mode |
| 2025-11-25 | 2.0 | Added comparison modes and interactive_comparison_mode |
| 2025-11-25 | 2.1 | Fixed window_size, boundary conditions, AlphaGenome constants |
| 2025-11-25 | 2.2 | Two-line ontology labels, refined "comparison" mode |

## 👥 Author

ChatGPT (for Alberto)  
Created: 2025-11-23  
Updated: 2025-11-25

## 📝 Final Notes

This program is an essential tool to ensure quality and consistency of the genomic data processing pipeline. Use it regularly during development and before training models to avoid training with incorrect data.

For technical questions or bugs, consult the code documentation or get in touch.
