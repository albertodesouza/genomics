# VCF to 23andMe Converter

Converts DNA data files from **VCF** (Variant Call Format) to the **23andMe** raw data format (V3 or V5). Automatically annotates VCF files with dbSNP rsIDs when needed.

---

## Table of Contents

- [Description](#description)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
  - [YAML Fields Reference](#yaml-fields-reference)
- [Usage](#usage)
- [dbSNP Annotation Pipeline](#dbsnp-annotation-pipeline)
  - [When Annotation Triggers](#when-annotation-triggers)
  - [How It Works](#how-it-works)
  - [Chromosome Name Handling](#chromosome-name-handling)
- [SNP Panel Filtering](#snp-panel-filtering)
- [Genome Build Auto-Detection](#genome-build-auto-detection)
- [Supported Formats](#supported-formats)
  - [Input: VCF](#input-vcf)
  - [Output: 23andMe](#output-23andme)
  - [V5 vs V3](#v5-vs-v3)
- [Examples](#examples)
- [Limitations](#limitations)

---

## Description

This converter reads a VCF file (gzipped or plain text) containing genomic variants and produces a file in the 23andMe raw data format, widely used by personal genomics services and ancestry analysis tools such as GEDmatch, Promethease, and DNA.Land.

Key features:

- **Automatic rsID annotation**: VCFs from variant callers (e.g. GATK HaplotypeCaller) often lack rsIDs. The converter detects this and annotates the VCF with dbSNP rsIDs via `bcftools` before conversion.
- **Genome build auto-detection**: automatically identifies GRCh37 or GRCh38 from VCF headers.
- **SNP panel filtering**: optionally restrict output to a specific set of SNPs (e.g. 23andMe V5 chip panel).
- **Multi-sample VCF support**: extract a specific sample from multi-sample VCFs.
- **V3 and V5 output formats**: generate 23andMe files in either format version.

---

## Prerequisites

- **Conda `genomics` environment** (see `scripts/start_genomics_universal.sh` in the repository root)
- **Python 3.10+**
- **PyYAML**
- **bcftools** (required for dbSNP annotation; included in the genomics conda environment)
- **wget** or **curl** (required for downloading dbSNP files)

---

## Installation

The conda `genomics` environment already includes all required dependencies. If you need to install them manually:

```bash
pip install pyyaml
conda install -c bioconda bcftools
```

---

## Configuration

The converter is configured through a YAML file in the `configs/` subdirectory. A template is provided at `configs/default.yaml`.

To get started, copy the template and edit it:

```bash
cp configs/default.yaml configs/my_conversion.yaml
nano configs/my_conversion.yaml
```

### YAML Fields Reference

| Section | Field | Type | Required | Default | Description |
|---------|-------|------|:--------:|---------|-------------|
| `input` | `vcf_file` | string | Yes | -- | Path to input VCF file (`.vcf` or `.vcf.gz`) |
| `input` | `sample` | string | No | first sample | Sample name to extract (for multi-sample VCFs) |
| `output` | `output_file` | string | Yes | -- | Path to the output 23andMe file |
| `output` | `format_version` | string | No | `"V5"` | Output format: `"V5"` or `"V3"` |
| `output` | `skip_no_rsid` | bool | No | `false` | If `true`, skip variants without rsID; if `false`, assign positional IDs (`i{chrom}_{pos}`) |
| `output` | `snp_panel` | string | No | -- | Path to a custom SNP panel file (one rsID per line); takes priority over `filter_by_chip_panel` |
| `output` | `filter_by_chip_panel` | bool | No | `false` | Auto-download and filter by the Illumina chip panel matching `format_version` (V5: GSA, V3: OmniExpress) |
| `annotation` | `enabled` | bool | No | `true` | Annotate VCF with dbSNP rsIDs if missing |
| `annotation` | `ref_dir` | string | Yes* | -- | Directory for dbSNP and chip panel files (*required when `enabled: true` or `filter_by_chip_panel: true`) |
| -- | `reference_build` | string | No | `"GRCh37"` | Fallback genome build: `"GRCh37"` or `"GRCh38"` |

---

## Usage

Before running the converter, activate the genomics conda environment:

```bash
source ../scripts/start_genomics_universal.sh
```

Then run the converter with a config file:

```bash
python vcf_to_23andme.py --config configs/my_conversion.yaml
```

To see command-line help:

```bash
python vcf_to_23andme.py --help
```

---

## dbSNP Annotation Pipeline

### When Annotation Triggers

Many VCFs produced by variant callers (GATK HaplotypeCaller, DeepVariant, etc.) do not include rsIDs -- the ID column contains `.` for all variants. The converter automatically detects this by sampling the first 500 variant lines: if more than 80% have a missing ID, annotation is triggered.

If the VCF already has rsIDs, this step is skipped entirely.

### How It Works

When annotation is needed, the converter:

1. **Downloads the dbSNP common variants file** to the directory specified in `annotation.ref_dir`. The file (~1.2 GB for GRCh37) is downloaded only on the first run and reused for subsequent conversions.
   - GRCh37: dbSNP b151 common variants
   - GRCh38: dbSNP b151 common variants

2. **Checks chromosome name compatibility** between the input VCF and the dbSNP file, and generates a renamed copy of the dbSNP file if needed (see [Chromosome Name Handling](#chromosome-name-handling)).

3. **Ensures the input VCF is indexed** (`bcftools` requires a tabix index). If the VCF is regular-gzipped (not bgzipped), the converter creates a bgzipped copy with a `.bgz.vcf.gz` suffix for indexing. **The original input file is never modified.**

4. **Runs `bcftools annotate -c ID`** to produce a new annotated VCF with rsIDs from dbSNP.

5. **Saves the annotated VCF** in the same directory as the input file, with a `.rsid.vcf.gz` suffix. For example:
   - Input: `A56102411139365.vcf.gz`
   - bgzip copy (if needed): `A56102411139365.bgz.vcf.gz`
   - Annotated output: `A56102411139365.rsid.vcf.gz`

6. The annotated VCF is preserved for future use. On subsequent runs, if the `.rsid.vcf.gz` file already exists, annotation is skipped.

### Chromosome Name Handling

VCF files use different chromosome naming conventions:

| Style | Example | Typical source |
|-------|---------|----------------|
| Numeric | `1`, `2`, ..., `X`, `Y` | GRCh37 VCFs, 1000 Genomes |
| chr-prefixed | `chr1`, `chr2`, ..., `chrX` | GRCh38 VCFs, UCSC |
| RefSeq | `NC_000001.11`, ... | NCBI dbSNP GRCh38 |

The converter automatically detects the naming styles of both the input VCF and the dbSNP file. If they differ, it generates a renamed copy of the dbSNP file with matching chromosome names (stored in `ref_dir` for reuse).

---

## SNP Panel Filtering

By default, the converter includes **all variants** from the VCF in the output. This is the most useful behavior for downstream tools, which select the SNPs they need.

There are two ways to filter by a SNP panel:

### Automatic chip panel download

Set `output.filter_by_chip_panel: true` to automatically download and use the Illumina chip manifest matching the selected `format_version`:

| Format | Illumina Chip | Approx. SNPs |
|--------|--------------|:------------:|
| V5 | Infinium Global Screening Array v1.0 | ~640K |
| V3 | HumanOmniExpress-12 v1.1 | ~960K |

The manifest is downloaded once, rsIDs are extracted, and the resulting panel file is cached in `annotation.ref_dir` (e.g. `panel_23andme_V5.txt`). On subsequent runs, the cached panel is reused instantly.

**Note**: 23andMe customizes these base Illumina chips (adding and removing some probes), so the Illumina manifests are a very close approximation (~99%) but not an exact match to the actual 23andMe chip content.

### Custom panel file

Alternatively, provide a custom panel file via `output.snp_panel` (one rsID per line). This takes priority over `filter_by_chip_panel`:

```
rs12345
rs67890
rs11111
```

When any panel is active, the converter reports how many of the panel SNPs were found in the VCF (panel coverage).

---

## Genome Build Auto-Detection

The converter tries to detect whether the VCF uses GRCh37 or GRCh38, using the following heuristics (in priority order):

1. **`##reference=` header line**: looks for keywords like `GRCh37`, `hg19`, `b37`, `GRCh38`, or `hg38`.
2. **`##contig=` header lines**: compares declared contig lengths against known chromosome sizes for each build (e.g. chr1 is 249,250,621 bp in GRCh37 and 248,956,422 bp in GRCh38).
3. **Chromosome naming convention**: the `chr` prefix in chromosome names suggests GRCh38; absence suggests GRCh37 (weakest heuristic).
4. **Fallback**: if none of the above is conclusive, uses the `reference_build` value from the YAML config (default: `GRCh37`).

---

## Supported Formats

### Input: VCF

The VCF (Variant Call Format) is the standard for representing genomic variants. The converter supports:

- `.vcf` (plain text) and `.vcf.gz` (gzip compressed) files
- Single-sample and multi-sample VCFs
- Phased (`0|1`) and unphased (`0/1`) genotypes
- Variants with or without rsIDs (automatic annotation when missing)

### Output: 23andMe

The output file uses tab-separated columns:

```
rsid    chromosome    position    genotype
```

Genotype encoding:

| Type | Example | Description |
|------|---------|-------------|
| Homozygous SNP | `AA`, `GG` | Both alleles identical |
| Heterozygous SNP | `AG`, `CT` | Different alleles |
| Hemizygous (male X, Y, MT) | `A`, `G` | Single allele |
| No-call | `--` | Genotype not determined |
| Homozygous insertion | `II` | Insertion on both alleles |
| Homozygous deletion | `DD` | Deletion on both alleles |
| Heterozygous indel | `DI` | One reference allele, one indel |

### V5 vs V3

- **V5**: current format with a detailed header including data quality disclaimers and reference assembly links.
- **V3**: older format with a shorter header.

Both formats use the same tab-separated data structure (4 columns).

---

## Examples

### Basic conversion with automatic annotation

```yaml
input:
  vcf_file: "/data/genomes/my_genome.vcf.gz"

output:
  output_file: "/data/genomes/my_genome_23andme.txt"
  format_version: "V5"

annotation:
  enabled: true
  ref_dir: "/data/refs"

reference_build: "GRCh37"
```

```bash
source ../scripts/start_genomics_universal.sh
python vcf_to_23andme.py --config configs/my_conversion.yaml
```

### Multi-sample VCF (extracting a specific sample)

```yaml
input:
  vcf_file: "/data/1000genomes/ALL.chr1.vcf.gz"
  sample: "NA12878"

output:
  output_file: "/data/genomes/NA12878_23andme.txt"
  format_version: "V5"

annotation:
  enabled: false

reference_build: "GRCh38"
```

### Filtering by chip panel (automatic download)

```yaml
input:
  vcf_file: "/data/genomes/my_genome.vcf.gz"

output:
  output_file: "/data/genomes/my_genome_23andme_v5.txt"
  format_version: "V5"
  filter_by_chip_panel: true

annotation:
  enabled: true
  ref_dir: "/data/refs"

reference_build: "GRCh37"
```

### Filtering by custom SNP panel

```yaml
input:
  vcf_file: "/data/genomes/my_genome.vcf.gz"

output:
  output_file: "/data/genomes/my_genome_23andme_panel.txt"
  format_version: "V5"
  snp_panel: "/data/panels/my_custom_panel.txt"

annotation:
  enabled: true
  ref_dir: "/data/refs"

reference_build: "GRCh37"
```

### Conversion without annotation (VCF already has rsIDs)

```yaml
input:
  vcf_file: "/data/genomes/annotated_genome.vcf.gz"

output:
  output_file: "/data/genomes/annotated_genome_23andme.txt"
  format_version: "V3"

annotation:
  enabled: false

reference_build: "GRCh37"
```

---

## Limitations

- **Variants without rsID**: when the VCF ID column is `.` and annotation is disabled (or dbSNP does not cover the variant), a positional identifier `i{chrom}_{pos}` is assigned. This can be changed with `skip_no_rsid: true`.
- **dbSNP coverage**: the converter uses the dbSNP common variants file, which covers the most frequently observed variants. Rare or novel variants may not receive rsIDs. For maximum coverage, you can replace the dbSNP file in `ref_dir` with the full `All_*.vcf.gz` file from NCBI (~15 GB for GRCh37).
- **Multi-allelic variants**: genotypes are determined from VCF allele indices, supporting all multi-allelic combinations.
- **Indels**: represented as `II`/`DD`/`DI` following 23andMe conventions.
- **Sex chromosomes**: the converter does not auto-detect sample sex; haploid genotypes on X/Y/MT appear as single letters only when the VCF genotype is already haploid.
- **No liftover**: the converter does not perform coordinate conversion between builds. Output positions match the input VCF coordinates.
- **Chip panel approximation**: the auto-downloaded Illumina manifests (GSA for V5, OmniExpress for V3) are the base chips used by 23andMe, but 23andMe customizes them by adding/removing probes. The manifests cover ~99% of the actual 23andMe chip content.
