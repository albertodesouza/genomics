# Non-Longevous Dataset Builder

> **📁 Location**: This module is in `build_non_longevous_dataset/`

Pipeline for building datasets from non-longevous individuals from the 1000 Genomes Project.

> **🚀 Quick Start**: New to this module? Check the **[Quick Start Guide](docs/QUICKSTART.md)** for a step-by-step introduction!

## 📑 Table of Contents

- [📋 Description](#-description)
  - [build_window_and_predict.py](#build_window_and_predictpy)
  - [Window Modes: Gene vs SNP](#window-modes-gene-vs-snp)
- [🔥 NEW: PyTorch Dataset](#-new-pytorch-dataset)
- [🎨 Visualization Tool](#-visualization-tool)
- [🔧 Requirements](#-requirements)
- [📊 CSV Format](#-csv-format)
- [🚀 Basic Usage](#-basic-usage)
  - [Step 1: Analyze Metadata](#step-1-analyze-metadata)
  - [Step 2: Configure Sample Selection](#step-2-configure-sample-selection)
  - [Step 3: Run Complete Pipeline](#step-3-run-complete-pipeline)
- [🧬 Window Modes](#-window-modes)
  - [Gene Mode (default)](#gene-mode-default)
  - [SNP Mode (AISNP analysis)](#snp-mode-aisnp-analysis)
  - [Gene List Mode](#gene-list-mode)
- [📁 Output Structure](#-output-structure)
- [🔄 Idempotence](#-idempotence)
- [⚙️ Advanced Options](#️-advanced-options)
- [📊 Output Example (Step 1)](#-output-example-step-1)
- [🧬 1000 Genomes Superpopulations](#-1000-genomes-superpopulations)
- [💡 Tips](#-tips)
- [🔍 Troubleshooting](#-troubleshooting)
- [📚 See Also](#-see-also)

---

## 📋 Description

This program analyzes a CSV file with metadata from 1000 Genomes Project individuals, allows sample selection based on custom criteria, and executes genomic analyses using `build_window_and_predict.py` for each selected individual.

### build_window_and_predict.py

The module includes `build_window_and_predict.py`, a script for:
- Extracting 1 Mb genomic windows centered on specific genes or SNPs
- Applying individual variants from 1000 Genomes to the reference genome
- [Generating consensus sequences per haplotype (H1 and H2)](docs/HAPLOTYPES.md)
- Running AlphaGenome predictions for functional analysis (RNA-seq, ATAC-seq, etc.)

📚 **AlphaGenome Documentation**:
- **[AlphaGenome Predictions Guide](docs/ALPHAGENOME_PREDICTIONS.md)** - Complete guide on running predictions, analyzing results, performance optimization, and common use cases
- **[Tissues/Cells Guide](docs/ALPHAGENOME_TISSUES.md)** - How to discover and use tissue/cell ontologies (UBERON, CL) with AlphaGenome

### Window Modes: Gene vs SNP

The pipeline supports two operating modes:

1. **Gene Mode** (default): Creates windows centered on gene bodies
   - Ideal for functional genomics studies
   - Requires GTF annotation
   - Supports single gene or gene list

2. **SNP Mode**: Creates windows centered on specific SNP positions
   - Ideal for ancestry-informative SNP (AISNP) analysis
   - Integrates with FROGAncestryCalc
   - No GTF required
   - Uses pre-defined SNP coordinates

---

## 🔥 NEW: PyTorch Dataset

> **✨ New**: `build_non_longevous_dataset` now automatically generates a complete PyTorch Dataset!
> 
> **📚 Full Documentation**: See the **[Complete PyTorch Dataset Guide](docs/PYTORCH_DATASET.md)** for API reference, advanced examples, and detailed usage instructions.

### Features

The pipeline now organizes data in a structure compatible with PyTorch Dataset/DataLoader:

- ✅ **`GenomicLongevityDataset` class** ready to use
- ✅ **Complete metadata** per individual and global
- ✅ **Ancestry data** from FROGAncestryCalc included
- ✅ **Lazy loading** (on-demand) to save memory
- ✅ **AlphaGenome predictions** (.npz) for each haplotype
- ✅ **DataLoader support** with custom collate function

### Output Structure

```
non_longevous_dataset/
├── dataset_metadata.json              # Global metadata
└── individuals/
    ├── HG01879/
    │   ├── individual_metadata.json   # Metadata + FROG likelihood
    │   └── windows/
    │       ├── CYP2B6/               # Genomic window
    │       │   ├── ref.window.fa
    │       │   ├── HG01879.H1.window.fixed.fa
    │       │   ├── HG01879.H2.window.fixed.fa
    │       │   ├── predictions_H1/
    │       │   │   ├── rna_seq.npz
    │       │   │   └── atac.npz
    │       │   └── predictions_H2/
    │       │       ├── rna_seq.npz
    │       │       └── atac.npz
    │       └── rs10497191/           # Another window
    └── HG01880/
        └── ...
```

### Usage Example

```python
from genomic_dataset import GenomicLongevityDataset
from torch.utils.data import DataLoader

# Load dataset
dataset = GenomicLongevityDataset(
    dataset_dir='non_longevous_results',
    load_predictions=True,
    load_sequences=False
)

print(f"Total: {len(dataset)} individuals")

# Access a sample
input_data, output_data = dataset[0]

print(f"Sample: {output_data['sample_id']}")
print(f"Population: {output_data['population']}")
print(f"Windows: {list(input_data['windows'].keys())}")

# Use with DataLoader
from genomic_dataset import collate_fn_variable_windows

dataloader = DataLoader(
    dataset,
    batch_size=4,
    shuffle=True,
    collate_fn=collate_fn_variable_windows
)

for batch_input, batch_output in dataloader:
    # Train neural network
    pass
```

### Data Format

Each sample returns `(input_data, output_data)`:

**Input Data** (for neural network):
```python
{
    'windows': {
        'CYP2B6': {
            'h1_sequence': str (optional),
            'h2_sequence': str (optional),
            'predictions_h1': {
                'rna_seq': np.ndarray,  # (1000000,)
                'atac': np.ndarray      # (1000000,)
            },
            'predictions_h2': {...}
        }
    }
}
```

**Output Data** (labels):
```python
{
    'sample_id': str,
    'longevity': 0,  # 0 for non-longevous, 1 for longevous
    'sex': 1,  # 1=Male, 2=Female
    'population': 'ACB',
    'superpopulation': 'AFR',
    'frog_likelihood': np.ndarray  # (150,) likelihood per population
}
```

### Complete Documentation

📚 **[Complete PyTorch Dataset Documentation](docs/PYTORCH_DATASET.md)**

Includes:
- Complete API
- Advanced examples
- Custom transformations
- Collate functions
- Compression and distribution
- FAQ

### Complete Example

Run the example script:

```bash
cd build_non_longevous_dataset
python3 examples/load_dataset_example.py --dataset-dir non_longevous_results
```

### Pipeline Configuration

To generate dataset metadata, enable in config YAML:

```yaml
pipeline:
  steps:
    run_predictions: true
    generate_dataset_metadata: true  # New step!
```

---

## 🎨 Visualization Tool

> **✨ New**: Interactive visualizer for exploring AlphaGenome predictions!
> 
> **📚 Full Documentation**: See the **[Visualization Guide](docs/VISUALIZATION.md)** for complete keyboard controls, usage examples, and tips.

### Quick Start

Visualize AlphaGenome predictions interactively:

```bash
cd build_non_longevous_dataset
python3 alphagenome_output_visualization.py configs/small.yaml
```

### Features

- ⌨️ **Keyboard navigation** through thousands of .npz files
- 🎭 **Overlay mode** to compare multiple predictions side-by-side
- 📊 **Group mode** to view ATAC + RNA-seq together
- 🎨 **High-contrast colors** (blue, red, green, orange, etc.)
- 📈 **Real-time statistics** (mean, standard deviation)
- 🚀 **Smart caching** for improved performance

### Keyboard Controls

| Key | Action |
|-----|--------|
| `Space` / `→` | Next file |
| `←` | Previous file |
| `a` / `A` | Activate/deactivate overlay mode |
| `d` / `D` | Switch output type (atac ↔ rna_seq) |
| `g` / `G` | Activate/deactivate group mode |
| `q` / `ESC` | Exit |

### Use Cases

1. **Quality Control**: Verify predictions look reasonable
2. **Compare Haplotypes**: Overlay H1 and H2 predictions
3. **Population Analysis**: Compare functional variation across samples
4. **Multi-omics View**: See ATAC and RNA-seq together

**→ [Complete Visualization Documentation](docs/VISUALIZATION.md)**

---

## 🔧 Requirements

- Python 3.8+
- Python packages:
  - pandas
  - pyyaml
  - numpy
  - alphagenome (for predictions)
  - **torch** (to use PyTorch Dataset - optional but recommended)
- Tools:
  - samtools
  - bcftools
- Files:
  - `build_window_and_predict.py` (included in this module)
  - GRCh38 reference genome (.fa + .fai)
  - 1000 Genomes VCFs (filtered and phased)

## 📊 CSV Format

The CSV file must contain the following columns:

```
FamilyID,SampleID,FatherID,MotherID,Sex,Population,Superpopulation
```

Where:
- **SampleID**: Unique individual identifier (e.g., HG00096)
- **Sex**: 1 = Male, 2 = Female
- **Population**: Population (e.g., ACB, GBR, CHB)
- **Superpopulation**: Superpopulation (AFR, EUR, EAS, SAS, AMR)

Example:
```csv
BB01,HG01879,0,0,1,ACB,AFR
BB01,HG01880,0,0,2,ACB,AFR
Y001,HG00096,0,0,1,GBR,EUR
```

## 🚀 Basic Usage

### Step 1: Analyze Metadata

First, analyze the CSV file to see statistics about available samples:

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml

# OR, from project root:
python3 build_non_longevous_dataset/build_non_longevous_dataset.py --config build_non_longevous_dataset/configs/default.yaml
```

This will print:
- Total number of samples
- How many superpopulations exist
- How many people in each superpopulation
- How many populations in each superpopulation
- Sex distribution in each population

### Step 2: Configure Sample Selection

Edit the `configs/default.yaml` file to configure:

1. **CSV path**:
```yaml
data_sources:
  metadata_csv: "../../docs/1000_genomes_metadata.csv"  # Relative to configs/
```

2. **Selection criteria**:
```yaml
sample_selection:
  level: "superpopulation"  # or "population"
  samples_per_group: 2       # how many samples per group
  sex_filter: "all"          # "all", "male", or "female"
```

3. **Window mode and target parameters**:
```yaml
build_window_params:
  mode: "gene"              # "gene" or "snp"
  
  gene:
    symbol: "CYP2B6"        # gene of interest
    # OR use gene_list_file for multiple genes
  
  window_size: 1000000      # 1 Mb window
  predict: true             # run AlphaGenome predictions
  outputs: "RNA_SEQ,ATAC"   # output types
  ontology: "UBERON:0002107,UBERON:0000955"  # tissues (liver, brain)
```

4. **Enable additional steps**:
```yaml
pipeline:
  steps:
    analyze_metadata: true    # Step 1: analyze CSV
    select_samples: true      # Step 2: select samples
    validate_vcfs: false      # Step 3: validate VCFs (optional)
    run_predictions: true     # Step 4: run predictions
    generate_report: true     # Step 5: generate report
```

### Step 3: Run Complete Pipeline

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

## 🧬 Window Modes

### Gene Mode (default)

Creates 1 Mb windows centered on gene bodies. Best for functional genomics.

**Single Gene:**
```yaml
build_window_params:
  mode: "gene"
  gene:
    symbol: "CYP2B6"      # HGNC gene symbol
    # OR
    id: "ENSG00000197894"  # ENSEMBL gene ID
```

**Output window:** `individuals/HG00096/windows/CYP2B6/`

### SNP Mode (AISNP analysis)

Creates 1 Mb windows centered on SNP positions. Integrates with FROGAncestryCalc for ancestry analysis.

> 🧬 **Complete AISNP Guide**: For detailed information about Ancestry-Informative SNPs analysis, see the **[AISNP Mode Documentation](docs/AISNP_MODE.md)**.

**Configuration:**
```yaml
build_window_params:
  mode: "snp"
  snp:
    snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"
```

**Features:**
- Processes all 55 AISNPs automatically
- No GTF annotation required
- Creates separate windows for each SNP
- Output windows organized hierarchically: `individuals/HG00096/windows/rs10497191/`, `individuals/HG00096/windows/rs1079597/`, etc.

**SNP File Format:**
Tab-delimited with header:
```
ALFRED_UID	dbSNP_rsnumber	chrom	chrom_pos	alleles
SI047925B	rs10497191	2	157810705	C/T
SI000148N	rs1079597	11	113425564	C/T
```

**Integration with FROGAncestryCalc:**
1. Extract AISNP genotypes using FROGAncestryCalc tools
2. Run ancestry inference to identify population-specific patterns
3. Use this pipeline to extract functional predictions around each AISNP
4. Compare epigenetic/functional profiles between ancestries

### Gene List Mode

Process multiple genes in a single run.

**Configuration:**
```yaml
build_window_params:
  mode: "gene"
  gene:
    gene_list_file: "my_genes.txt"
```

**Gene list file format** (one gene per line):
```
CYP2B6
APOE
TP53
BRCA1
# Comments start with #
ENSG00000197894
```

**Output:** Creates one window directory per gene: `individuals/HG00096/windows/CYP2B6/`, `individuals/HG00096/windows/APOE/`, etc.

## 📁 Output Structure

The pipeline generates a well-organized directory structure with all results:

```
non_longevous_results/
├── metadata_statistics.json                      # CSV statistics
├── selected_samples.csv                          # Selected samples
├── dataset_metadata.json                         # Global dataset metadata (if enabled)
├── non_longevous_dataset_checkpoint.json         # Checkpoint (idempotence)
├── processing_summary.txt                        # Final report
└── individuals/                                  # Individual results
    ├── HG00096/
    │   ├── individual_metadata.json              # Metadata + ancestry data
    │   └── windows/
    │       ├── CYP2B6/                          # Gene mode window
    │       │   ├── ref.window.fa
    │       │   ├── HG00096.H1.window.fixed.fa
    │       │   ├── HG00096.H2.window.fixed.fa
    │       │   ├── predictions_H1/
    │       │   │   ├── rna_seq.npz
    │       │   │   ├── rna_seq_metadata.json
    │       │   │   ├── atac.npz
    │       │   │   └── atac_metadata.json
    │       │   └── predictions_H2/
    │       │       ├── rna_seq.npz
    │       │       ├── rna_seq_metadata.json
    │       │       ├── atac.npz
    │       │       └── atac_metadata.json
    │       └── rs10497191/                      # SNP mode window
    │           ├── ref.window.fa
    │           ├── HG00096.H1.window.fixed.fa
    │           ├── HG00096.H2.window.fixed.fa
    │           ├── predictions_H1/
    │           │   ├── rna_seq.npz
    │           │   └── atac.npz
    │           └── predictions_H2/
    │               ├── rna_seq.npz
    │               └── atac.npz
    └── HG00097/
        └── ...
```

**Key Features:**

- **Hierarchical organization**: `individuals/` → `{SampleID}/` → `windows/` → `{Target}/`
- **Multiple windows**: Each sample can have multiple windows (genes or SNPs)
- **Complete haplotypes**: Reference sequence + H1 and H2 haplotypes per window
- **AlphaGenome predictions**: Separate directories for H1 and H2 predictions
- **Metadata files**: JSON metadata alongside each .npz prediction file
- **PyTorch-ready**: Structure compatible with `GenomicLongevityDataset`

**For SNP mode with 55 AISNPs and 78 samples:**
- Total sample directories: 78 (`individuals/{SampleID}/`)
- Total window directories: 4,290 (78 samples × 55 SNPs)
- Each window contains: 1 reference + 2 haplotypes + 2 prediction directories
- Total .npz files (2 outputs): ~17,160 files (78 × 55 × 2 haplotypes × 2 outputs)

> 📖 **Detailed Structure Documentation**: See the **[Complete Structure Guide](docs/STRUCTURE.md)** for a comprehensive explanation of all output files, their formats, and how they integrate with the PyTorch Dataset.

## 🔄 Idempotence

The program is idempotent and maintains a checkpoint file. If execution is interrupted:

1. Already processed samples will **not** be reprocessed
2. The pipeline will continue from where it stopped
3. To reprocess everything, delete the checkpoint file:
   ```bash
   rm non_longevous_results/non_longevous_dataset_checkpoint.json
   ```

## ⚙️ Advanced Options

### Select Only Some Populations

```yaml
sample_selection:
  level: "population"
  samples_per_group: 5
  include_groups: ["GBR", "CHB", "YRI"]  # only these populations
```

### Exclude Populations

```yaml
sample_selection:
  exclude_groups: ["ACB", "ASW"]  # exclude these
```

### Filter by Sex

```yaml
sample_selection:
  sex_filter: "male"  # males only
```

### Disable AlphaGenome Predictions (faster)

```yaml
build_window_params:
  predict: false  # only extract sequences
```

### Process Only Haplotype 1 (faster)

```yaml
build_window_params:
  skip_h2: true  # skip H2, only build H1
```

### Configure API Rate Limiting

To respect AlphaGenome API usage policies, the pipeline applies a configurable delay **after each individual API call** (not just between samples):

```yaml
build_window_params:
  api_rate_limit_delay: 0.5  # seconds between API calls
```

**Important:** Each sample requires multiple API calls:
- **SNP mode (55 AISNPs)**: 55 SNPs × 2 haplotypes = **110 API calls per sample**
- **Gene mode**: Number of genes × 2 haplotypes

**Time estimation example** (78 samples, 55 SNPs, 0.5s delay):
- Total API calls: 78 × 55 × 2 = **8,580 calls**
- Delay time alone: ~71 minutes
- API processing time (~7s/call): ~16.7 hours
- **Total estimated time: ~18-20 hours**

Recommended values:
- `0.5` - Moderate rate limiting (recommended for most cases)
- `1.0` - Conservative rate limiting (very respectful)
- `0.0` - No rate limiting (use with caution, may violate API policies)

## 📊 Output Example (Step 1)

```
================================================================================
DATASET STATISTICS - 1000 GENOMES PROJECT
================================================================================

📊 TOTAL SAMPLES: 56

🌍 SUPERPOPULATIONS: 5
--------------------------------------------------------------------------------

  AFR:
    • Total individuals: 16
    • Male: 8
    • Female: 8
    • Number of populations: 2
    • Populations: ACB, ASW

  AMR:
    • Total individuals: 10
    • Male: 5
    • Female: 5
    • Number of populations: 2
    • Populations: MXL, PUR

  EAS:
    • Total individuals: 10
    • Male: 5
    • Female: 5
    • Number of populations: 2
    • Populations: CHB, CHS

  EUR:
    • Total individuals: 10
    • Male: 5
    • Female: 5
    • Number of populations: 2
    • Populations: GBR, TSI

  SAS:
    • Total individuals: 10
    • Male: 5
    • Female: 5
    • Number of populations: 2
    • Populations: GIH, ITU

🏘️  POPULATIONS: 10
--------------------------------------------------------------------------------

  AFR:
    ACB: 10 individuals (♂ 5, ♀ 5)
    ASW: 6 individuals (♂ 3, ♀ 3)

  AMR:
    MXL: 4 individuals (♂ 2, ♀ 2)
    PUR: 6 individuals (♂ 3, ♀ 3)

  ...
```

## 🧬 1000 Genomes Superpopulations

- **AFR**: African
- **AMR**: Ad Mixed American
- **EAS**: East Asian
- **EUR**: European
- **SAS**: South Asian

## 💡 Tips

1. **Start with analysis**: Run only the `analyze_metadata` step first to understand your data
2. **Test with few samples**: Use `samples_per_group: 1` or `2` for quick tests
3. **Use checkpoint**: The system saves progress automatically
4. **Multiple ontologies**: Separate by comma: `"UBERON:0002107,UBERON:0000955,CL:0002601"`
5. **VCF per chromosome**: Make sure the VCF contains the chromosome of your gene/SNP
6. **SNP mode for ancestry**: Use SNP mode with FROGAncestryCalc's AISNP list for ancestry-specific functional analysis
7. **Gene list for batch**: Use gene list mode to process multiple genes in parallel

## 🔍 Troubleshooting

### Error: CSV not found
```
[ERROR] CSV file not found: ...
```
**Solution**: Check the path in `data_sources.metadata_csv` in YAML

### Error: VCF pattern contains {chrom}
```
[WARN] VCF pattern contains {chrom}, but chromosome was not determined.
```
**Solution**: Provide the full VCF path or specify the chromosome

### Error: API key not found
```
RuntimeError: AlphaGenome API key not provided
```
**Solution**: 
```bash
export ALPHAGENOME_API_KEY="your_key_here"
```

### Error: SNP file not found (SNP mode)
```
[ERROR] SNP mode requires snp.snp_list_file in config
```
**Solution**: Check that the SNP file path is correct and relative to config directory:
```yaml
snp:
  snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"
```

### Error: No gene specified (Gene mode)
```
[ERROR] Gene mode requires either gene.symbol, gene.id, or gene.gene_list_file
```
**Solution**: Specify at least one gene parameter in config

## 📚 See Also

- **[PyTorch Dataset Guide](docs/PYTORCH_DATASET.md)** - Complete Dataset documentation
- **[Visualization Guide](docs/VISUALIZATION.md)** - Interactive visualization tool
- [FROGAncestryCalc](../FROGAncestryCalc/README.md) - Ancestry inference using AISNPs
- [AlphaGenome Predictions Guide](docs/ALPHAGENOME_PREDICTIONS.md)
- [AlphaGenome Tissues Guide](docs/ALPHAGENOME_TISSUES.md)
- [Dataset Examples](examples/load_dataset_example.py) - Complete usage examples

## 📝 Author

Alberto F. De Souza  
Last updated: 2025-11-09
