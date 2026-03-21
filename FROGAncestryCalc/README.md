# FROGAncestryCalc - Modified Version

FROG-kb (Forensic Resource/Reference On Genetics - Knowledge base) Ancestry Inference Batch Likelihood Computation Tool - Modified to use pipe delimiters.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Index
- [Quick Start](#-quick-start)
- [Configuration](#-configuration)
  - [Available AI Panels](#available-ai-panels)
- [Input File Format](#-input-file-format)
  - [Format Specifications](#format-specifications)
  - [Preparing Your Input File](#preparing-your-input-file)
- [Extracting SNPs from Genomic Data](#-extracting-snps-from-genomic-data)
  - [Available Tools](#available-tools)
  - [Option 1: From 1000 Genomes Project](#option-1-from-1000-genomes-project)
  - [Option 2: From Your Own VCF File](#option-2-from-your-own-vcf-file)
  - [VCF Conversion Notes](#vcf-conversion-notes)
  - [Option 3: From Whole Genome Sequencing](#option-3-from-whole-genome-sequencing)
  - [Option 4: From a BAM File](#option-4-from-a-bam-file)
  - [SNP List Reference](#snp-list-reference)
  - [Installation of Required Tools](#installation-of-required-tools)
  - [Notes on Genome Builds](#notes-on-genome-builds)
- [Output Files](#-output-files)
  - [Visualization](#visualization)
- [Admixture Estimation](#-admixture-estimation)
  - [MLE Mode](#mle-mode)
  - [ADMIXTURE Mode](#admixture-mode)
  - [Panel Selection](#panel-selection)
  - [Caveats and Limitations](#caveats-and-limitations)
- [Project Structure](#️-project-structure)
- [Requirements](#️-requirements)
- [Modifications from Original](#-modifications-from-original)
- [Error Handling](#-error-handling)
- [Maintenance](#️-maintenance)
- [Population Coverage](#-population-coverage)
- [Troubleshooting](#-troubleshooting)
- [License](#-license)
- [Acknowledgments](#-acknowledgments)
- [Support](#-support)
- [Appendix: bcftools Commands for 1000 Genomes VCF Extraction](#-appendix-bcftools-commands-for-1000-genomes-vcf-extraction)
- [Appendix: How FROGAncestryCalc Handles Haplotypes](#-appendix-how-frogancestrycalc-handles-haplotypes)

---

## 🚀 Quick Start

### Run Analysis
```bash
./run.sh
```

### Recompile Code
```bash
./recompile.sh
```

## 📋 Configuration

Edit the `FROGAncestryCalc.properties` file:

```properties
homePath=.
inputFilename=55_aisnp_1000_Genome.txt
panelInfo=55AI
```

**⚠️ IMPORTANT:** Update the properties file with the appropriate input file name and AI panel name before starting a new job.

### Available AI Panels

| Panel Code | Description | SNP Count |
|-----------|-------------|-----------|
| `55AI` | KiddLab - Set of 55 AISNPs | 55 |
| `128AI` | Seldin's list of 128 AISNPs | 128 |
| `34plex` | SNPforID 34-plex | 34 |
| `combined` | Combined panel (Kiddlab-55 + Seldin's-128 + SNPforID34-plex) | 192 |
| `precision` | Precision ID Ancestry Panel | 165 |

## 📂 Input File Format

Place your input files in the `input/` directory with the following format:

```
Individual|rs10497191|rs1079597|rs11652805|...|rs9522149
HG02561_GWD|NN|CC|CC|CC|...|TT
HG02562_GWD|TT|CT|CC|CC|...|TT
```

### Format Specifications

- ✅ **Delimiter:** pipe `|`
- ✅ **Line endings:** Unix (LF)
- ✅ **Encoding:** UTF-8
- ✅ **First line:** Header with "Individual" + ordered list of SNP IDs
- ✅ **Following lines:** Individual ID + genotypes
- ✅ **SNP order:** Must match the order in the corresponding sample file
- ✅ **Individual IDs:** Must be unique

### Preparing Your Input File

1. Follow the SNP order given in the sample files for your chosen AI panel (see `sampleInputFiles/`)
2. SNP labels and genotypes must be ordered by rs number (alphanumeric)
3. Use the sorting function in Excel or similar tools (ascending order)
4. Ensure all Individual Identifiers are unique
5. Consult the appropriate file in `SNPInfo/` to find valid alleles for each SNP
6. Use accepted genotype notations:
   - Two-allele format: `AA`, `TT`, `GG`, `CC`, `AT`, `AG`, etc.
   - Missing data: `NN`

## 🧬 Extracting SNPs from Genomic Data

The `tools/` directory contains scripts to extract the required SNPs from various genomic data sources:

### Available Tools

| Tool | Description |
|------|-------------|
| `vcf_to_frog.py` | Convert VCF files to FROGAncestryCalc format |
| `plot_frog_ancestry_pie.py` | Generate a pie chart and summary TSV from likelihood output |
| `admixture_estimate.py` | Estimate ancestry admixture proportions (MLE and/or ADMIXTURE) |
| `run_frog_from_bam.sh` | End-to-end pipeline: BAM → variant calling → FROG analysis |
| `extract_snps_from_1000genomes.sh` | Download and extract SNPs from 1000 Genomes Project |
| `extract_snps_from_wgs.sh` | Extract SNPs from whole genome sequencing data (FASTQ/BAM/VCF) |
| `aisnps_55_list.txt` | List of the 55 AISNP rs IDs |
| `aisnps_128_list.txt` | List of the 128 AISNP rs IDs |
| `aisnps_55_grch37.bed` | BED coordinates for 55 AISNPs (GRCh37) |
| `aisnps_128_grch37.bed` | BED coordinates for 128 AISNPs (GRCh37) |

### Option 1: From 1000 Genomes Project

Download and extract data from the 1000 Genomes Project. Supports both GRCh37 Phase 3 (default) and GRCh38 High Coverage:

**Using GRCh37 Phase 3 (default, ~20GB):**
```bash
# Extract all samples (auto-detects existing VCFs)
./tools/extract_snps_from_1000genomes.sh

# Extract specific samples only
echo -e "HG02561\nHG02562\nHG03055" > my_samples.txt
./tools/extract_snps_from_1000genomes.sh -s my_samples.txt -o input/my_samples.txt
```

**Using GRCh38 High Coverage (~35GB, requires coordinate conversion):**
```bash
# First, generate GRCh38 coordinates (only once)
python3 tools/convert_grch37_to_grch38.py

# Then extract data
./tools/extract_snps_from_1000genomes.sh -b grch38 -s my_samples.txt -o input/my_samples.txt
```

**Requirements:**
- `bcftools` (install via: `conda install -c bioconda bcftools`)
- `wget`
- Python 3 with `requests` (for GRCh38)
- Disk space: ~20 GB (GRCh37) or ~35 GB (GRCh38)

### Option 2: From Your Own VCF File

If you already have a VCF file (from sequencing, microarray, or other sources):

```bash
# Recommended for typical GRCh37/hg19 variant-only VCFs
python3 tools/vcf_to_frog.py \
    your_samples.vcf.gz \
    tools/aisnps_55_list.txt \
    input/your_data.txt \
    SNPInfo/55_aisnps_alleles.txt \
    --missing-mode nn
```

```bash
# GRCh38/hg38 VCFs
python3 tools/vcf_to_frog.py \
    your_samples.vcf.gz \
    tools/aisnps_55_list.txt \
    input/your_data.txt \
    SNPInfo/55_aisnps_alleles_grch38.txt \
    --missing-mode nn
```

```bash
# Backward-compatible behavior: missing SNPs become REF/REF
python3 tools/vcf_to_frog.py \
    your_samples.vcf.gz \
    tools/aisnps_55_list.txt \
    input/your_data.txt \
    SNPInfo/55_aisnps_alleles.txt
```

### VCF Conversion Notes

- `vcf_to_frog.py` accepts an optional `alleles_file` argument that is used to:
  - map genomic coordinates (`chrom:pos`) to `rsID`
  - optionally fill missing SNPs as `REF/REF`
- Many VCFs store `.` in the `ID` column instead of `rsIDs`. In that case, providing the correct `SNPInfo/*.txt` file is strongly recommended, otherwise the converter may not recognize the target AISNPs.
- Missing SNP handling is controlled by `--missing-mode`:
  - `--missing-mode ref` (default): missing SNPs become `REF/REF` when reference alleles are available
  - `--missing-mode nn`: missing SNPs always become `NN`
- For **variant-only VCFs** (common in sequencing and consumer pipelines), `--missing-mode nn` is recommended because missing sites may mean:
  - homozygous reference
  - filtered out
  - not confidently called
- Use `--missing-mode ref` only when absence of a variant confidently implies homozygous reference, such as in an all-sites VCF or another well-characterized pipeline.

### Option 3: From Whole Genome Sequencing

Extract SNPs from raw sequencing data:

```bash
# From VCF
./tools/extract_snps_from_wgs.sh \
    -i sample.vcf.gz \
    -t vcf \
    -o input/sample.txt

# From BAM (aligned reads)
./tools/extract_snps_from_wgs.sh \
    -i sample.bam \
    -t bam \
    -r GRCh38.fa \
    -o input/sample.txt

# From FASTQ (paired-end)
./tools/extract_snps_from_wgs.sh \
    -i sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -t fastq \
    -r GRCh38.fa \
    -o input/sample.txt
```

**Requirements for WGS:**
- `bcftools`, `samtools` (for all types)
- `bwa` (for FASTQ alignment)
- Reference genome (GRCh37/hg19 or GRCh38/hg38)

### Option 4: From a BAM File

If you have aligned reads in BAM format, the `run_frog_from_bam.sh` script handles the full pipeline automatically:

```bash
# Basic usage (auto-detects genome build and reference)
bash tools/run_frog_from_bam.sh -i sample.bam

# With pie chart generation
bash tools/run_frog_from_bam.sh -i sample.bam -p

# With explicit reference genome
bash tools/run_frog_from_bam.sh -i sample.bam -r /path/to/reference.fasta -p
```

The script performs the following steps:
1. Detects the BAM genome build (GRCh37 or GRCh38) and chromosome naming convention
2. Locates or downloads a matching reference genome
3. Calls genotypes at the 55 AISNP positions using `bcftools mpileup` + `bcftools call`
4. Converts the resulting VCF to FROGAncestryCalc format
5. Runs FROGAncestryCalc and optionally generates visualizations

For the 128 AISNP panel, manual extraction is required (see [Admixture Estimation](#-admixture-estimation) for details).

### SNP List Reference

The 55 AISNPs used in this panel are:

<details>
<summary>Click to expand SNP list</summary>

```
rs10497191, rs1079597, rs11652805, rs1229984, rs12439433, rs12498138,
rs12913832, rs1426654, rs1462906, rs1572018, rs16891982, rs174570,
rs17642714, rs1800414, rs1834619, rs1871534, rs1876482, rs192655,
rs200354, rs2024566, rs2042762, rs2166624, rs2196051, rs2238151,
rs2593595, rs260690, rs2814778, rs310644, rs3737576, rs3811801,
rs3814134, rs3823159, rs3827760, rs3916235, rs4411548, rs4471745,
rs459920, rs4833103, rs4891825, rs4918664, rs671, rs6754311,
rs6990312, rs7226659, rs7251928, rs7326934, rs735480, rs7554936,
rs7657799, rs7722456, rs798443, rs7997709, rs870347, rs917115,
rs9522149
```

Full list available in: `tools/aisnps_55_list.txt`

</details>

### Installation of Required Tools

The easiest way to install bioinformatics tools is via conda:

```bash
# Create a new conda environment with all tools
conda create -n genomics \
    bcftools samtools bwa gatk4 \
    python=3.9 -c bioconda -c conda-forge

# Activate the environment
conda activate genomics

# Verify installation
bcftools --version
samtools --version
bwa
```

### Notes on Genome Builds

FROGAncestryCalc now supports both major genome builds:

- **GRCh37/hg19 (default)**:
  - 1000 Genomes Phase 3 data
  - Most legacy sequencing data
  - Coordinates in `SNPInfo/55_aisnps_alleles.txt`
  
- **GRCh38/hg38**:
  - 1000 Genomes High Coverage data
  - Modern sequencing data
  - Coordinates in `SNPInfo/55_aisnps_alleles_grch38.txt` (generated by conversion script)

**To use GRCh38:**
```bash
# Generate GRCh38 coordinates (only once)
python3 tools/convert_grch37_to_grch38.py

# Then use -b grch38 flag
./tools/extract_snps_from_1000genomes.sh -b grch38
```

**Check your data's build:**
```bash
# Check VCF header
bcftools view -h your_file.vcf.gz | grep "##reference"

# Check BAM header
samtools view -H your_file.bam | grep "@SQ"
```

## 📊 Output Files

Generated in the `output/` directory:

| File | Description |
|------|-------------|
| `*_likelihood.txt` | Likelihood values for ancestral population for each individual across 155 populations |
| `*_orderOfMag.txt` | Order of magnitude of the likelihoods |
| `*_rankOrder.txt` | Population rankings by likelihood for each individual |

All output files are tab-delimited and can be opened in Excel.

**Note:** Output files from previous jobs (including any `errFile.txt`) are deleted at the start of a new job.

### Visualization

You can generate a quick visualization from the likelihood output:

```bash
# Auto-detect a single *_likelihood.txt file in output/
python3 tools/plot_frog_ancestry_pie.py output

# Or point directly to a specific likelihood file
python3 tools/plot_frog_ancestry_pie.py output/your_sample_likelihood.txt

# Show more populations explicitly
python3 tools/plot_frog_ancestry_pie.py output --top-n 10
```

The visualization script creates:

| File | Description |
|------|-------------|
| `*_ancestry_pie.png` | Pie chart of normalized relative likelihoods |
| `*_ancestry_summary.tsv` | Full population ranking with relative percentages and metadata |

**Important:** The displayed percentages are **normalized relative likelihoods**, not exact biological admixture fractions. They are useful for visualization and comparison, but should not be interpreted as precise ancestry percentages.

## 🧪 Admixture Estimation

FROGAncestryCalc's native analysis computes the likelihood that an individual belongs to each reference population independently, ranking populations by best fit. However, it does **not** estimate **admixture proportions** -- i.e., which fraction of an individual's ancestry derives from each ancestral group.

The `admixture_estimate.py` script fills this gap by providing two complementary approaches for estimating ancestry admixture:

### MLE Mode

Maximum Likelihood Estimation using FROG's own reference allele frequencies (`referenceData_*.txt`). Estimates the proportions of ancestry from each superpopulation that best explain the individual's observed genotypes.

```bash
# MLE with the default 55 AISNP panel
python3 tools/admixture_estimate.py mle input/sample_55aisnps.txt

# MLE with the 128 AISNP panel
python3 tools/admixture_estimate.py mle input/sample_128aisnps.txt --panel 128AI

# Control bootstrap iterations (default: 100; set to 0 to skip)
python3 tools/admixture_estimate.py mle input/sample.txt --panel 128AI --n-bootstrap 50
```

Output files:
- `*_admixture_mle.tsv` -- proportions with optional confidence intervals
- `*_admixture_mle_pie.png` -- pie chart visualization

### ADMIXTURE Mode

Runs the official [ADMIXTURE](https://dalexander.github.io/admixture/) tool in supervised mode with a 1000 Genomes Phase 3 reference panel. Requires a FROG-formatted reference file with 1000 Genomes genotypes.

```bash
# ADMIXTURE supervised mode with 1000G reference
python3 tools/admixture_estimate.py admixture input/sample.txt \
    --reference input/whole_1000genomes_55aisnps.txt

# Run both MLE and ADMIXTURE and generate a comparison chart
python3 tools/admixture_estimate.py both input/sample.txt \
    --reference input/whole_1000genomes_55aisnps.txt
```

**Requirements:** `admixture` and `plink` must be installed and available in `$PATH`.

### Panel Selection

The `--panel` argument selects which reference data and alleles file to use:

| Panel | Reference Data | Alleles File | Populations | SNPs |
|-------|---------------|--------------|-------------|------|
| `55AI` (default) | `referenceData_55.txt` | `55_aisnps_alleles.txt` | 158 | 55 |
| `128AI` | `referenceData_128.txt` | `128_aisnps_alleles.txt` | 128 | 128 |

When `--panel` is specified, the `--reference-data` and `--alleles` arguments are set automatically. They can still be overridden explicitly if needed.

### Extracting 128 AISNPs from a BAM File

To use the 128 AISNP panel, extract genotypes from your BAM at the 128 positions:

```bash
# 1. Create target positions (GRCh37; use appropriate coords for GRCh38)
awk 'BEGIN{OFS="\t"} NR>1{print $3,$4}' SNPInfo/128_aisnps_alleles.txt | \
    LC_ALL=C sort -k1,1V -k2,2n > /tmp/128_targets.tsv
bgzip -f /tmp/128_targets.tsv && tabix -f -s1 -b2 -e2 /tmp/128_targets.tsv.gz

# 2. Call genotypes
bcftools mpileup -Ou -f reference.fasta -T /tmp/128_targets.tsv.gz -q20 -Q20 sample.bam | \
    bcftools call -m -Oz -o output/sample_128aisnps.vcf.gz
bcftools index output/sample_128aisnps.vcf.gz

# 3. Convert to FROG format
python3 tools/vcf_to_frog.py \
    output/sample_128aisnps.vcf.gz \
    tools/aisnps_128_list.txt \
    input/sample_128aisnps.txt \
    SNPInfo/128_aisnps_alleles.txt \
    --missing-mode nn

# 4. Run admixture estimation
python3 tools/admixture_estimate.py mle input/sample_128aisnps.txt --panel 128AI
```

### Caveats and Limitations

1. **FROGAncestryCalc vs. admixture estimation.** FROGAncestryCalc answers "which *single* population does this individual most resemble?", not "what mixture of populations produced this individual?". For admixed individuals (e.g., most people of Latin American, Caribbean, or African-American ancestry), the likelihood rankings may be misleading because the best-matching single population is often a proxy rather than a true ancestral source. Use `admixture_estimate.py` for admixture proportions.

2. **Small marker panels.** Both 55 and 128 AISNPs are very small panels compared to commercial ancestry services (which use hundreds of thousands of SNPs). Results have wide confidence intervals and should be treated as rough estimates, not precise measurements. The 128 panel provides better resolution than the 55 panel but still cannot distinguish closely related populations.

3. **Superpopulation-level estimation only.** The MLE admixture model groups populations into broad superpopulations (African, European, East Asian, South Asian, American, Middle Eastern & North African, Central Asian, Oceanian). Estimating proportions for all 128+ individual populations from 128 SNPs is statistically underdetermined (as many parameters as data points) and would produce unreliable results.

4. **"American" ancestry reflects Indigenous/Amerindian heritage.** In the 128 AISNP panel, the "American" superpopulation includes reference populations such as Maya, Karitiana, Surui, Quechua, Nahuas, and Pima. A high "American" component in the MLE results indicates Indigenous American ancestry. The 55 panel lacks sufficient Amerindian reference populations, so Indigenous ancestry may instead appear as "East Asian" (due to shared deep ancestry) or be absorbed by other categories.

5. **Panel-dependent results.** The 55 and 128 panels use largely different sets of SNPs (only 13 overlap) and different reference populations (158 vs. 128). Results from the two panels are not directly comparable and may assign substantially different proportions to the same superpopulations. This reflects differences in the informative content and reference data of each panel, not errors in the method.

6. **ADMIXTURE mode requires a 1000 Genomes reference file.** The supervised ADMIXTURE mode needs a FROG-formatted file containing genotypes for 1000 Genomes individuals at the target SNPs (`whole_1000genomes_*aisnps.txt`). This file is currently available only for the 55 AISNP panel. For the 128 panel, use MLE mode.

7. **No mitochondrial or Y-chromosome markers.** All AISNPs are autosomal. The analysis captures genome-wide admixture proportions but provides no information about specific maternal (mtDNA) or paternal (Y-chromosome) lineages.

## 🗂️ Project Structure

```
FROGAncestryCalc/
├── src/                        # Modified source code
│   ├── bean/                   # Data classes
│   ├── dv/                     # Validation (modified for pipes)
│   ├── main/                   # Main application class
│   ├── read/                   # File reading (modified)
│   ├── sub/                    # Helper classes
│   └── write/                  # Output writing
├── bin/                        # Compiled classes
├── input/                      # Input files directory
│   ├── ind/                    # Working directory (do not delete)
│   └── indGenotype/            # Working directory (do not delete)
├── output/                     # Results directory
├── SNPInfo/                    # SNP information for each panel
│   ├── 55_aisnps_alleles.txt
│   ├── 128_aisnps_alleles.txt
│   ├── 34_plex_alleles.txt
│   ├── combined_alleles.txt
│   └── precision_alleles.txt
├── sampleInputFiles/           # Sample input files
│   ├── 55_aisnps_sample.txt
│   ├── 128_aisnps_sample.txt
│   ├── 34_plex_sample.txt
│   ├── combined_sample.txt
│   └── precision_sample.txt
├── log/                        # Execution logs
│   └── workingLog.txt
├── tools/                      # Data extraction and analysis tools
│   ├── vcf_to_frog.py          # VCF converter
│   ├── plot_frog_ancestry_pie.py  # Likelihood visualization
│   ├── admixture_estimate.py   # Admixture proportion estimation (MLE / ADMIXTURE)
│   ├── run_frog_from_bam.sh    # BAM-to-FROG end-to-end pipeline
│   ├── extract_snps_from_1000genomes.sh  # 1000G extractor
│   ├── extract_snps_from_wgs.sh          # WGS extractor
│   ├── aisnps_55_list.txt      # SNP list (55 panel)
│   ├── aisnps_128_list.txt     # SNP list (128 panel)
│   ├── aisnps_55_grch37.bed    # BED coordinates (55 panel, GRCh37)
│   └── aisnps_128_grch37.bed   # BED coordinates (128 panel, GRCh37)
├── obsolete/                   # Original files with bugs
├── run.sh                      # Execution script
├── recompile.sh                # Recompilation script
├── FROGAncestryCalc.properties # Configuration file
└── MODIFICACOES.md             # Technical modification details
```

## ⚙️ Requirements

- **Java:** 17+ (OpenJDK recommended)
- **Shell:** Bash
- **OS:** Linux/Unix
- **Python 3** (for tools): `numpy`, `scipy`, `matplotlib`
- **Optional** (for ADMIXTURE mode): `admixture`, `plink` (install via conda/bioconda)

## 🔄 Modifications from Original

This modified version includes the following improvements:

### 1. **Pipe Delimiter Support**
- Changed from comma (`,`) to pipe (`|`) delimiter
- Modified validation and parsing logic
- Updated error messages

### 2. **Locale Fix**
- Added `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8` to prevent number formatting issues
- Resolves `NumberFormatException` with scientific notation (e.g., `5.652E-62`)

### 3. **Linux Compatibility**
- Unix line endings (LF)
- Proper path handling
- Shell script optimizations

### Modified Files:
- `src/dv/ValidateFileHeader.java` - Validation logic
- `src/read/ReadTxtFiles.java` - File parsing

For complete technical details, see [`MODIFICATIONS.md`](MODIFICATIONS.md)

## 📝 Error Handling

### Error File
If validation errors occur, check `output/errFile.txt` for details.

### Working Log
View `log/workingLog.txt` for:
- Processing information for all jobs
- Copy of error messages
- Timestamps and status updates

**Note:** The log file accumulates across jobs until manually deleted.

## 🛠️ Maintenance

### Recompiling After Code Changes

```bash
./recompile.sh
```

Or manually:

```bash
cd /path/to/FROGAncestryCalc
rm -rf bin
mkdir bin
javac -d bin -sourcepath src $(find src -name "*.java")
cp -r src/read/data bin/read/
```

### Cleaning Up

```bash
# Clean output files
rm -f output/*.txt

# Clean working directories
rm -f input/ind/* input/indGenotype/*

# Clean logs (optional)
rm -f log/workingLog.txt
```

## 📚 Population Coverage

The tool calculates ancestry likelihoods across multiple reference panels:

- **55AI panel:** 158 populations (22 from 1000 Genomes + 136 from ALFRED/HGDP/other)
- **128AI panel:** 128 populations (all from ALFRED/HGDP/other)
- **combined panel:** 192 SNPs across the union of populations from both panels

Major population groups include:
- African populations (Yoruba, Mbuti, Biaka, Luhya, etc.)
- European populations (Danes, Finns, British, Toscani, etc.)
- East Asian populations (Han Chinese variants, Japanese, Korean, Vietnamese, etc.)
- South Asian populations (Gujarati, Punjabi, Bengali, Telugu, etc.)
- American populations (Maya, Pima, Karitiana, Peruvian, etc.)
- Middle Eastern populations (Druze, Palestinian Arabs, Saudi, etc.)
- And many more regional populations...

### Population Mapping Files

Two comprehensive mapping files are available to understand the population structure:

#### 1. **`population_mapping_complete.csv`** - All Populations

Lists all populations in FROGAncestryCalc across all panels and maps each to a superpopulation category (African, European, East Asian, South Asian, American, Middle Eastern & North African, Central Asian, Oceanian). Includes:
- 22 populations from **1000 Genomes Project Phase 3** (with full metadata)
- All populations from **ALFRED, HGDP, and other sources** used across 55AI, 128AI, and other panels

#### 2. **`population_mapping_1000genomes.csv`** - 1000 Genomes Only

Contains only the **22 populations from 1000 Genomes Phase 3** that relates:

1. **FROGAncestryCalc population names** → **1000 Genomes population codes**
2. **1000 Genomes populations** → **Superpopulations**

#### Population-Superpopulation Structure in 1000 Genomes Project

The **1000 Genomes Project Phase 3 / High Coverage** has **26 populations** organized into 5 major superpopulations:

| Superpopulation | 1000G Total | In FROGAncestryCalc | Populations in FROGAncestryCalc | Missing from FROG |
|----------------|-------------|---------------------|----------------------------------|-------------------|
| African (AFR) | 7 | 6 | LWK, ASW, YRI, GWD, MSL, ESN | **ACB** |
| European (EUR) | 5 | 5 | CEU, GBR, FIN, IBS, TSI | none |
| East Asian (EAS) | 5 | 5 | CHB, CHS, JPT, CDX, KHV | none |
| South Asian (SAS) | 5 | 5 | GIH, PJL, BEB, STU, ITU | none |
| Admixed American (AMR) | 4 | 1 | PEL | **CLM, MXL, PUR** |
| **TOTAL** | **26** | **22** | | **4 missing** |

**⚠️ Important Notes:**

1. **FROGAncestryCalc has 22 of 26 populations** from 1000 Genomes. The 4 missing populations are:
   - **ACB** (African Caribbean in Barbados) - AFR
   - **CLM** (Colombian in Medellin) - AMR
   - **MXL** (Mexican Ancestry in Los Angeles) - AMR
   - **PUR** (Puerto Rican in Puerto Rico) - AMR

2. **Use High Coverage data when possible**: The 1000 Genomes Project has two main releases:
   - **Phase 3** (2015): ~2,500 individuals, low coverage (~7.4x)
   - **High Coverage** (2022): ~3,200 individuals, high coverage (~30x) - **RECOMMENDED**
   
   The High Coverage data provides better quality and includes 602 complete trios. Use GRCh38 coordinates for High Coverage data.

#### Using the Mapping Files

**View all 160 populations:**
```bash
# See all populations with their source database
cat population_mapping_complete.csv

# Count populations by source
grep "1000 Genomes" population_mapping_complete.csv | wc -l  # 22 populations
grep "ALFRED/HGDP" population_mapping_complete.csv | wc -l   # 138 populations
```

**Work with 1000 Genomes populations only:**
```bash
# View only 1000 Genomes populations
cat population_mapping_1000genomes.csv

# Filter by superpopulation (e.g., African)
grep ",AFR," population_mapping_1000genomes.csv

# Filter by specific population code
grep ",GBR," population_mapping_1000genomes.csv

# Or from the complete file
grep "1000 Genomes" population_mapping_complete.csv | grep ",AFR,"
```

#### Example Mapping Entries

**From `population_mapping_complete.csv` (includes all populations):**
```csv
FROGAncestryCalc_Population,Source_Database,Pop_Code_1000G,Full_Name_1000G,Superpopulation_1000G,Superpopulation_Full_Name
Yoruba(YRI),1000 Genomes Phase 3,YRI,Yoruba in Ibadan Nigeria,AFR,African
British(GBR),1000 Genomes Phase 3,GBR,British in England and Scotland,EUR,European
Mbuti,ALFRED/HGDP/Other,,,,
Maya-Yucatec,ALFRED/HGDP/Other,,,,
Druze,ALFRED/HGDP/Other,,,,
```

**From `population_mapping_1000genomes.csv` (1000 Genomes only):**
```csv
FROGAncestryCalc_Population,Pop_Code_1000G,Full_Name_1000G,Superpopulation_1000G,Superpopulation_Full_Name
Yoruba(YRI),YRI,Yoruba in Ibadan Nigeria,AFR,African
British(GBR),GBR,British in England and Scotland,EUR,European
Han Chinese(CHB),CHB,Han Chinese in Beijing China,EAS,East Asian
Gujarati(GIH),GIH,Gujarati Indian in Houston TX,SAS,South Asian
Peruvian(PEL),PEL,Peruvian in Lima Peru,AMR,Admixed American
```

#### References

1. **The 1000 Genomes Project Consortium** (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74. https://doi.org/10.1038/nature15393

2. **Byrska-Bishop, M., et al.** (2022). High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell*, 185(18), 3426-3440.e19. https://doi.org/10.1016/j.cell.2022.08.004

3. **1000 Genomes Project Official Website**. Population descriptions and metadata. https://www.internationalgenome.org/

4. **Kidd, K.K., et al.** (2014). Progress toward an efficient panel of SNPs for ancestry inference. *Forensic Science International: Genetics*, 10, 23-32. https://doi.org/10.1016/j.fsigen.2014.01.002

5. **ALFRED (ALlele FREquency Database)**. The ALlele FREquency Database. https://alfred.med.yale.edu/

6. **FROG-kb (Forensic Resource/Reference On Genetics - Knowledge base)**. https://frog.med.yale.edu/

**Important Notes:**
- The **1000 Genomes Project** has **26 populations total** (Phase 3 / High Coverage data)
- FROGAncestryCalc includes only **22 of those 26** (missing: ACB, CLM, MXL, PUR)
- Of the **160 total populations** in FROGAncestryCalc:
  - **22 are from 1000 Genomes** (13.75%)
  - **138 are from ALFRED, HGDP, and regional studies** (86.25%)
- The `population_mapping_complete.csv` file lists all 160 populations and identifies their source
- The `population_mapping_1000genomes.csv` file contains only the 22 populations from 1000 Genomes with full metadata

## 🐛 Troubleshooting

### Common Issues

1. **"Your input file is not pipe delimited"**
   - Ensure file uses pipe `|` as delimiter, not comma
   - Check for proper line endings (Unix LF, not Windows CRLF)

2. **"NumberFormatException"**
   - Make sure to run with proper locale settings
   - Use `./run.sh` which handles this automatically

3. **"Missing SNPs" or "Wrong SNP count"**
   - Verify SNP order matches the sample file
   - Check that all required SNPs are present
   - Ensure no extra columns or missing data

4. **"No SNPs found in VCF"**
   - Many VCFs do not store `rsIDs` in the `ID` column
   - Provide `SNPInfo/55_aisnps_alleles.txt` for GRCh37 or `SNPInfo/55_aisnps_alleles_grch38.txt` for GRCh38 so the converter can map `chrom:pos` to `rsID`
   - Verify that your genome build matches the alleles file you supplied

5. **"Ancestry result looks overly concentrated"**
   - This is common when plotting normalized likelihoods as percentages
   - For variant-only VCFs, avoid assuming that missing sites are homozygous reference unless you are sure that is correct
   - Re-run `vcf_to_frog.py` with `--missing-mode nn` for a more conservative input file
   - Prefer interpreting `*_rankOrder.txt` and `*_ancestry_summary.tsv` as rankings/relative affinity, not exact admixture fractions

6. **Java version issues**
   - Requires Java 17 or higher
   - Check version: `java -version`

## 📄 License

MIT License

Copyright (c) 2019 haseenaR

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

## 🙏 Acknowledgments

- Original FROG-kb tool by haseenaR
- Reference populations from various genomics databases
- Modified for improved usability and Linux compatibility

## 📞 Support

For issues related to:
- **Original tool:** Refer to FROG-kb documentation
- **This modified version:** Check `MODIFICACOES.md` for technical details

---

## 📚 Appendix: bcftools Commands for 1000 Genomes VCF Extraction

This appendix provides useful `bcftools` commands for working directly with 1000 Genomes VCF files.

### 1. Extract Specific Individual from a Chromosome

```bash
# Example: extract individual HG02561 from chromosome 1
bcftools view -s HG02561 \
  /path/to/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  -Oz -o HG02561_chr1.vcf.gz
```

### 2. Extract Multiple Individuals

```bash
# Extract HG02561, HG02562 and HG03055 from chromosome 1
bcftools view -s HG02561,HG02562,HG03055 \
  /path/to/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  -Oz -o three_samples_chr1.vcf.gz
```

### 3. Extract Only SNPs (no INDELs or structural variants)

```bash
# Extract only biallelic SNPs from HG02561
bcftools view -s HG02561 -v snps -m2 -M2 \
  /path/to/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  -Oz -o HG02561_chr1_snps.vcf.gz
```

### 4. Extract Specific Genomic Regions (by coordinates)

```bash
# Extract specific region from chromosome 2 (positions 158667216-158667217)
bcftools view -s HG02561 -r chr2:158667216-158667217 \
  /path/to/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

### 5. List Available Samples in a VCF

```bash
# View which samples are available
bcftools query -l \
  /path/to/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | head -20
```

### 6. Use BED File for Specific Regions

```bash
# Extract specific positions using BED file
bcftools view -s HG02561,HG02562,HG03055 \
  -R aisnps_55_grch38.bed \
  /path/to/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

### 7. Extract and View Genotypes in Readable Format

```bash
# Query specific format: chromosome, position, REF, ALT, and genotypes
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  -s HG02561,HG02562,HG03055 \
  -r chr1:1000000-2000000 \
  /path/to/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

### 8. Concatenate Multiple Chromosomes

```bash
# Concatenate all chromosomes into a single VCF
bcftools concat --threads 16 \
  /path/to/1kGP_high_coverage_Illumina.chr*.vcf.gz \
  -Oz -o all_chromosomes.vcf.gz

# Index the concatenated file
bcftools index all_chromosomes.vcf.gz
```

### Useful bcftools view Options

| Option | Description |
|--------|-------------|
| `-s` | Select samples (comma-separated) |
| `-r` | Specific region(s), format `chr:start-end` |
| `-R` | BED file with regions |
| `-v snps` | Only SNPs |
| `-v indels` | Only INDELs |
| `-m2 -M2` | Only biallelic variants |
| `-Oz` | Output compressed with bgzip |
| `-o` | Output file |
| `--threads N` | Use N threads for processing |

### Important Note on Genome Builds

⚠️ **Genome Build Compatibility:** The 1000 Genomes High Coverage VCFs use **GRCh38/hg38** coordinates. If you're working with older references or annotation files that use **GRCh37/hg19** coordinates, you will need to:

1. Use liftOver or CrossMap to convert coordinates between builds, OR
2. Use the Phase 3 VCFs (which use GRCh37/hg19) instead of High Coverage VCFs

For FROGAncestryCalc's 55 AISNPs, the coordinates in `SNPInfo/55_aisnps_alleles.txt` are in **GRCh37/hg19** format.

---

## 📚 Appendix: How FROGAncestryCalc Handles Haplotypes

### Overview

FROGAncestryCalc does **not** work with [haplotypes](../build_non_longevous_dataset/docs/HAPLOTYPES.md) separately. Instead, it processes genotypes as unphased allele pairs, which is sufficient for ancestry inference based on allele frequencies.

### Genotype Processing

The module uses a two-allele representation for each SNP:

```
Individual|rs10497191|rs1079597|rs11652805|...
HG02561|NN|CC|CC|...
HG02562|TT|CT|CC|...
```

Where:
- `CC`, `TT`, `GG`, `AA` = Homozygous (both haplotypes have the same allele)
- `CT`, `AG`, etc. = Heterozygous (haplotypes have different alleles)
- `NN` = Missing data

### Phase Information is Not Preserved

When extracting genotypes from phased VCFs (which contain haplotype information), the `vcf_to_frog.py` script:

1. **Accepts phased genotypes**: Can read `0|1` (phased) or `0/1` (unphased)
2. **Sorts alleles alphabetically**: Always outputs alleles in alphabetical order
3. **Loses phase information**: `0|1` and `1|0` both become the same output (e.g., `AG`)

**Example:**
```
VCF Genotype    →    FROGAncestryCalc Output
0|1 (A|G)       →    AG
1|0 (G|A)       →    AG  (same as above!)
0|0 (A|A)       →    AA
1|1 (G|G)       →    GG
```

### Why This Approach Works

Ancestry inference algorithms used by FROGAncestryCalc:

1. **Calculate allele frequencies** in reference populations
2. **Compute likelihood** of observing the genotype set in each population
3. **Compare likelihoods** across populations

For this allele frequency-based analysis, **it doesn't matter which allele is on which haplotype**. What matters is:
- How many copies of each allele the individual has (0, 1, or 2)
- The frequencies of those alleles in reference populations

### Implications

**Advantages:**
- ✅ Simpler input format
- ✅ Works with unphased data
- ✅ Sufficient for population ancestry inference
- ✅ Compatible with various genotyping platforms

**Limitations:**
- ❌ Cannot detect haplotype-specific patterns
- ❌ Cannot analyze phase-dependent effects
- ❌ Cannot track which parental lineage contributed which allele

### When Haplotype Information Matters

For analyses that require haplotype information (such as functional genomics or regulatory studies), see the `build_non_longevous_dataset` module, which:
- Extracts separate consensus sequences for H1 and H2
- Preserves phase information from phased VCFs
- Enables haplotype-specific functional predictions

### Related Documentation

- [Understanding Haplotypes and Consensus Sequences](../build_non_longevous_dataset/docs/HAPLOTYPES.md) - Detailed explanation of haplotype concepts
- [AISNP Mode](../build_non_longevous_dataset/docs/AISNP_MODE.md) - Functional analysis of AISNPs with haplotype resolution

---

**Note:** The original JAR and Windows batch files have been moved to `obsolete/` as they contained bugs incompatible with this input format.

