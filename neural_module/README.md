# Neural Module - DNA Analysis with AlphaGenome

> **📁 Location**: This module is in `neural_module/`

## 📑 Index

- [📋 Description](#-description)
- [🚀 Features](#-features)
- [🔗 Integration with Genomic Pipelines](#-integration-with-genomic-pipelines)
- [📦 Installation](#-installation)
- [🎯 Usage](#-usage)
- [📊 Available Output Types](#-available-output-types)
- [📝 Input Format (FASTA)](#-input-format-fasta)
- [📁 Output Structure](#-output-structure)
- [⚙️ Command Line Options](#️-command-line-options)
- [💡 Advanced Examples](#-advanced-examples)
- [🔧 Troubleshooting](#-troubleshooting)
- [📚 Resources](#-resources)
- [🔬 Use Cases](#-use-cases)
- [⚠️ Limitations](#️-limitations)
- [🧬 Haplotype Processing Details](#-haplotype-processing-details)
- [📄 License](#-license)
- [🤝 Contributions](#-contributions)
- [📧 Support](#-support)

---

## 📋 Description

`neural_module.py` is a Python tool that integrates the [AlphaGenome](https://github.com/google-deepmind/alphagenome) API from Google DeepMind to perform advanced DNA sequence analyses.

## 🚀 Features

- **Sequence Analysis**: Prediction of multiple functional aspects of DNA sequences
- **Gene Expression Prediction**: RNA-seq, CAGE
- **Chromatin Features**: ATAC-seq, H3K27AC, H3K4ME3, H3K27ME3, H3K36ME3, H3K9ME3, CTCF
- **Variant Analysis**: Prediction of SNP variant effects
- **Visualizations**: Automatic generation of plots in multiple formats (PNG, PDF, SVG)
- **FASTA Support**: Reading of standard FASTA files

## 🔗 Integration with Genomic Pipelines

The Neural Module includes **neural_integration.py**, a powerful bridge tool that connects traditional variant calling pipelines with AI-powered functional analysis.

### Key Capabilities:
- 🔄 **Automated Workflow**: VCF → Sequence Extraction → Neural Analysis → Results Correlation
- 📊 **Multiple Input Formats**: VCF (variants), BED (regions), GTF (genes)
- 🎯 **4 Operation Modes**: Integrated analysis, VCF extraction, BED extraction, Gene extraction
- 🧬 **Smart Extraction**: Automatically extracts ±5kb regions around variants
- 📈 **Results Correlation**: Links variant calls with functional predictions

### Quick Example:
```bash
# Extract variants and analyze with AlphaGenome in one command
cd neural_module
python neural_integration.py \
  --integrated \
  --vcf ../vcf/sample.vcf.gz \
  --ref ../refs/GRCh38.fa \
  --api-key YOUR_API_KEY \
  --output integrated_results/
```

📖 **Complete Integration Guide**: [docs/INTEGRATION.md](docs/INTEGRATION.md)

---

## 📦 Installation

### 1. Install AlphaGenome

Run the provided installation script:

```bash
bash install_alphagenome.sh
```

Or manually:

```bash
git clone https://github.com/google-deepmind/alphagenome.git
pip install ./alphagenome
```

### 2. Install Additional Dependencies

```bash
pip install rich matplotlib
```

### 3. Obtain API Key

Visit [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/) and request your free API key for non-commercial use.

## 🎯 Usage

### Basic Example

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/
```

### Analysis with Specific Outputs

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --outputs RNA_SEQ ATAC H3K27AC
```

### Variant Analysis

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --variant 1000 A C
```

This command analyzes the effect of an A→C variant at position 1000 (relative to the start of the sequence).

### Multiple Output Formats

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --formats png pdf svg --dpi 600
```

### Analysis Only (No Plots)

```bash
python neural_module.py -i sequence.fasta -k YOUR_API_KEY -o results/ \
    --no-plots
```

## 📊 Available Output Types

| Output | Description |
|--------|-------------|
| `RNA_SEQ` | Gene expression prediction via RNA-seq |
| `CAGE` | Cap Analysis of Gene Expression |
| `ATAC` | Chromatin accessibility (ATAC-seq) |
| `H3K27AC` | Marker of active regulatory elements |
| `H3K4ME3` | Marker of active promoters |
| `H3K27ME3` | Marker of gene repression |
| `H3K36ME3` | Marker of active gene bodies |
| `H3K9ME3` | Marker of heterochromatin |
| `CTCF` | Insulator binding factor |

## 📝 Input Format (FASTA)

The FASTA file must follow the standard format:

```
>sequence_id_1 description
ATCGATCGATCGATCG...
>sequence_id_2 description
GCTAGCTAGCTAGCTA...
```

### Requirements:

- Sequences with 100 bp to 1,000,000 bp (1 Mbp)
- Valid characters: A, C, G, T, N (and IUPAC ambiguity codes)

## 📁 Output Structure

```
results/
├── sequence_id_1_RNA_SEQ.png
├── sequence_id_1_ATAC.png
├── sequence_id_1_H3K27AC.png
├── sequence_id_2_RNA_SEQ.png
├── ...
└── analysis_report.json
```

### JSON Report

The `analysis_report.json` file contains:

```json
{
  "timestamp": "2025-10-16T10:30:00",
  "total_sequences": 2,
  "successful_analyses": 2,
  "sequences": [
    {
      "id": "sequence_id_1",
      "length": 50000,
      "status": "success",
      "outputs": ["RNA_SEQ", "ATAC", "H3K27AC"]
    }
  ]
}
```

## ⚙️ Command Line Options

### Required

- `-i, --input`: Input FASTA file
- `-k, --api-key`: AlphaGenome API key
- `-o, --output`: Output directory

### Optional

- `--outputs`: Desired output types (default: RNA_SEQ CAGE ATAC H3K27AC H3K4ME3)
- `--chromosome`: Reference chromosome (default: chr1)
- `--start`: Reference start position (default: 1000000)
- `--variant POS REF ALT`: Analyze variant at position POS with bases REF>ALT
- `--formats`: Plot formats (png, pdf, svg) (default: png)
- `--dpi`: Plot resolution (default: 300)
- `--no-plots`: Don't generate plots (analysis only)

## 💡 Advanced Examples

### 1. Complete Analysis of a Genomic Region

```bash
python neural_module.py \
    -i chr1_region.fasta \
    -k YOUR_API_KEY \
    -o chr1_analysis/ \
    --chromosome chr1 \
    --start 1000000 \
    --outputs RNA_SEQ CAGE ATAC H3K27AC H3K4ME3 H3K27ME3 CTCF \
    --formats png pdf \
    --dpi 600
```

### 2. Quick Analysis of Multiple Sequences

```bash
python neural_module.py \
    -i multiple_sequences.fasta \
    -k YOUR_API_KEY \
    -o quick_analysis/ \
    --outputs RNA_SEQ ATAC \
    --formats png \
    --dpi 150
```

### 3. Pathogenic Variant Analysis

```bash
python neural_module.py \
    -i disease_gene.fasta \
    -k YOUR_API_KEY \
    -o variant_effect/ \
    --variant 5000 G T \
    --formats png pdf svg
```

## 🔧 Troubleshooting

### Error: "AlphaGenome is not installed"

```bash
pip install git+https://github.com/google-deepmind/alphagenome.git
```

### Error: "Invalid API key"

Check that your API key is correct and active at [alphagenomedocs.com](https://www.alphagenomedocs.com/)

### Error: "Sequence too long"

AlphaGenome supports up to 1 Mbp. Split longer sequences into smaller chunks.

### Problem with Matplotlib

```bash
pip install --upgrade matplotlib seaborn
```

## 📚 Resources

- **AlphaGenome Documentation**: [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/)
- **GitHub Repository**: [https://github.com/google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome)
- **Paper**: Avsec et al. 2025 - "AlphaGenome: advancing regulatory variant effect prediction"

## 🔬 Use Cases

1. **Regulatory Variant Analysis**: Identify the impact of SNPs in regulatory regions
2. **Functional Element Prediction**: Identify promoters, enhancers, insulators
3. **Gene Expression Studies**: Predict expression levels in different tissues
4. **Chromatin Analysis**: Study accessibility and histone modifications
5. **Functional Genomics**: Characterize unknown DNA sequences

## ⚠️ Limitations

- Requires internet connection (online API)
- Limited query rate (check terms of use)
- Free use only for non-commercial research
- Sequences from 100 bp to 1 Mbp
- Processing time varies with sequence size
- **Haplotype processing**: AlphaGenome processes each [haplotype](../build_non_longevous_dataset/docs/HAPLOTYPES.md) independently (see [Haplotype Processing Details](#-haplotype-processing-details) below)

## 🧬 Haplotype Processing Details

### How AlphaGenome Handles Haplotypes

AlphaGenome **does NOT process both [haplotypes](../build_non_longevous_dataset/docs/HAPLOTYPES.md) simultaneously**. When analyzing genomic sequences from diploid organisms (like humans), each haplotype is treated as a completely independent sequence.

### Separate Processing Model

When you have phased genomic data with separate consensus sequences for H1 (haplotype 1) and H2 (haplotype 2):

```python
# Each haplotype is submitted separately
h1_seq = load_sequence("consensus_H1.fasta")  # 1 Mb from haplotype 1
h2_seq = load_sequence("consensus_H2.fasta")  # 1 Mb from haplotype 2

# AlphaGenome processes them independently
predictions_h1 = client.predict_sequence(h1_seq, ...)  # First call: H1 only
predictions_h2 = client.predict_sequence(h2_seq, ...)  # Second call: H2 only
```

### Key Implications

**What happens:**
- ✅ Each haplotype gets its own independent functional prediction
- ✅ You can identify allele-specific differences by comparing H1 vs H2 results
- ✅ Useful for detecting regulatory variants affecting one allele differently

**What doesn't happen:**
- ❌ AlphaGenome doesn't "know" it's analyzing two haplotypes from the same individual
- ❌ No modeling of interactions between haplotypes
- ❌ No representation of the true diploid cellular state
- ❌ Each prediction assumes the sequence exists in isolation

### Practical Example

Suppose you have a regulatory variant present only on H1:

```
Position:        chr2:1,000,500
Variant:         C→T (only on H1, not on H2)

AlphaGenome Results:
├── H1 prediction: Shows increased ATAC signal at position 500,500
│                  (T allele creates new TF binding site)
└── H2 prediction: Shows normal ATAC signal at position 500,500
                   (C allele maintains reference state)

Real cell state: Intermediate effect, with allele-specific binding
```

### Working with Results

To approximate diploid cellular behavior, you can post-process the predictions:

```python
# Load predictions from both haplotypes
h1_atac = np.load("predictions_H1/atac.npz")['values']
h2_atac = np.load("predictions_H2/atac.npz")['values']

# Calculate allele-specific differences
difference = h1_atac - h2_atac
significant_regions = np.where(np.abs(difference) > threshold)

# Approximate diploid state (simple average)
diploid_prediction = (h1_atac + h2_atac) / 2
```

### Architecture Limitation

This is a fundamental limitation of the AlphaGenome model architecture, which was trained on:
- **Input**: Single DNA sequence of up to 1,000,000 bases (A, C, G, T)
- **Output**: Functional predictions for that sequence alone

The model has no built-in mechanism to:
- Accept two sequences simultaneously
- Model *cis* interactions between homologous chromosomes
- Represent diploid cellular states

### When This Matters Most

**Critical for analyses involving:**
- Allele-specific expression (ASE)
- Imprinted genes
- Regulatory variants with dominant/recessive effects
- Parent-of-origin effects

**Less critical for:**
- Homozygous variants (both haplotypes identical)
- Analyses focused on structural features
- Initial variant screening

### Integration with Other Modules

For pipelines that generate haplotype-resolved consensus sequences (like `build_non_longevous_dataset`), this module will:
1. Process H1 and H2 as separate jobs
2. Save predictions with `_H1` and `_H2` suffixes
3. Allow downstream comparative analysis

For detailed information about haplotype generation and phasing, see:
- [Understanding Haplotypes and Consensus Sequences](../build_non_longevous_dataset/docs/HAPLOTYPES.md)
- [AISNP Haplotype Analysis Mode](../build_non_longevous_dataset/docs/AISNP_MODE.md)

## 📄 License

This module is distributed under the Apache 2.0 license, compatible with AlphaGenome.

## 🤝 Contributions

Contributions are welcome! Please open issues or pull requests in the project repository.

## 📧 Support

For AlphaGenome questions: alphagenome@google.com
For questions about this module: open an issue in the repository

---

**Developed for integration with genomes_analyzer.py**
