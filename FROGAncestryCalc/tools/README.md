# Genomic Data Extraction Tools

This directory contains scripts to extract the 55 Ancestry Informative SNPs (AISNPs) from various genomic data sources.

## 📁 Files

- **`aisnps_55_list.txt`** - List of 55 SNP rsIDs (one per line)
- **`vcf_to_frog.py`** - Python script to convert VCF to FROGAncestryCalc format
- **`extract_snps_from_1000genomes.sh`** - Download and extract from 1000 Genomes
- **`extract_snps_from_wgs.sh`** - Extract from WGS data (FASTQ/BAM/VCF)

## 🚀 Quick Start

### From 1000 Genomes Project

```bash
# Basic usage - extracts all samples (auto-detects existing VCFs)
./tools/extract_snps_from_1000genomes.sh

# Extract specific samples
echo -e "HG02561\nHG02562\nHG03055" > samples.txt
./tools/extract_snps_from_1000genomes.sh -s samples.txt -o input/my_data.txt
```

### From Your Own VCF

```bash
python3 tools/vcf_to_frog.py \
    your_data.vcf.gz \
    tools/aisnps_55_list.txt \
    input/output.txt
```

### From Sequencing Data

```bash
# From VCF
./tools/extract_snps_from_wgs.sh -i sample.vcf.gz -t vcf -o input/sample.txt

# From BAM
./tools/extract_snps_from_wgs.sh -i sample.bam -t bam -r genome.fa -o input/sample.txt

# From FASTQ
./tools/extract_snps_from_wgs.sh -i R1.fq.gz -2 R2.fq.gz -t fastq -r genome.fa -o input/sample.txt
```

## 📖 Detailed Documentation

### vcf_to_frog.py

Converts VCF files to FROGAncestryCalc pipe-delimited format.

**Usage:**
```bash
python3 vcf_to_frog.py <input.vcf.gz> <snp_list.txt> <output.txt>
```

**Arguments:**
- `input.vcf.gz` - VCF file (can be gzipped or plain text)
- `snp_list.txt` - File with SNP IDs (one rsID per line)
- `output.txt` - Output file in FROGAncestryCalc format

**Features:**
- Handles both phased (`0|1`) and unphased (`0/1`) genotypes
- Automatically handles gzipped files
- Converts numeric genotypes to allele notation (e.g., `0/1` → `AG`)
- Reports missing SNPs
- Handles multi-allelic sites (uses first alternate allele)

**Example:**
```bash
python3 tools/vcf_to_frog.py \
    1000genomes.vcf.gz \
    tools/aisnps_55_list.txt \
    input/1000g_samples.txt
```

**Output Format:**
```
Individual|rs10497191|rs1079597|rs11652805|...
Sample1|CC|CT|AA|...
Sample2|TT|CC|AG|...
```

---

### extract_snps_from_1000genomes.sh

Downloads 1000 Genomes data (Phase 3 GRCh37 or High Coverage GRCh38) and extracts the 55 AISNPs.
Auto-detects existing VCF files and skips download if already present.

**Usage:**
```bash
./tools/extract_snps_from_1000genomes.sh [options]
```

**Options:**
- `-b BUILD` - Genome build: `grch37` or `grch38` (default: `grch37`)
- `-d DIR` - VCF directory (default: auto-set based on build)
  - GRCh37: `/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes_GRCh37`
  - GRCh38: `/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes`
  - Script auto-detects existing VCFs and skips download if present
- `-o FILE` - Output file (default: `input/1000genomes_55aisnps.txt`)
- `-s FILE` - Sample list file (one sample ID per line)
- `-h` - Show help message

**Requirements:**
- `bcftools` (install via conda: `conda install -c bioconda bcftools`)
- `wget`
- Python 3 with `requests` library
- Disk space:
  - GRCh37 Phase 3: ~20 GB
  - GRCh38 High Coverage: ~35 GB
- Internet connection (only if downloading)
- **For GRCh38:** Run `tools/convert_grch37_to_grch38.py` first to generate GRCh38 coordinates

**What it does:**
1. Checks for existing VCF files in specified directory
2. Downloads missing VCF files for chromosomes 1-22 and X (if needed)
   - Uses robust download with automatic retry (up to 5 attempts)
   - Detects stalled downloads with 5-minute timeout
   - Supports resume for interrupted downloads
   - Exponential backoff between retries (2s, 4s, 8s, 16s)
3. Checks for cached extracted SNPs (reuses if found)
4. Concatenates all chromosomes (if cache not found)
5. Extracts only the 55 target SNPs and caches result
6. Optionally filters for specific samples
7. Converts to FROGAncestryCalc format

**Examples:**
```bash
# Extract all samples using GRCh37 Phase 3 (default)
./tools/extract_snps_from_1000genomes.sh

# Extract using GRCh38 High Coverage (requires coordinate conversion first)
python3 tools/convert_grch37_to_grch38.py
./tools/extract_snps_from_1000genomes.sh -b grch38

# Extract specific samples with GRCh37
./tools/extract_snps_from_1000genomes.sh \
    -b grch37 \
    -s african_samples.txt \
    -o input/african_1000g.txt

# Use custom VCF directory with GRCh38
./tools/extract_snps_from_1000genomes.sh \
    -b grch38 \
    -d /path/to/vcf/directory \
    -o input/custom.txt
```

**Data Sources:**
- **GRCh37 Phase 3:**
  - FTP: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
  - Genome build: GRCh37/hg19
  - Version: Phase 3 (low coverage)
- **GRCh38 High Coverage:**
  - FTP: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
  - Genome build: GRCh38/hg38
  - Version: High Coverage (30x)

---

### extract_snps_from_wgs.sh

Processes whole genome sequencing data and extracts the 55 AISNPs.

**Usage:**
```bash
./tools/extract_snps_from_wgs.sh -i INPUT -t TYPE -o OUTPUT [options]
```

**Required Arguments:**
- `-i FILE` - Input file (FASTQ, BAM, or VCF)
- `-t TYPE` - Input type: `fastq`, `bam`, or `vcf`
- `-o FILE` - Output file for FROGAncestryCalc

**Optional Arguments:**
- `-r FILE` - Reference genome FASTA (required for `fastq` and `bam`)
- `-n NAME` - Sample name (default: derived from input filename)
- `-2 FILE` - Second FASTQ file for paired-end reads
- `-k` - Keep intermediate files
- `-h` - Show help message

**Requirements:**
- `bcftools`, `samtools` - For variant calling and manipulation
- `bwa` - For FASTQ alignment (only needed for `-t fastq`)
- Reference genome - GRCh37/hg19 or GRCh38/hg38 FASTA file

**Processing Pipeline:**

**For VCF:**
1. Extract 55 SNPs from VCF
2. Convert to FROGAncestryCalc format

**For BAM:**
1. Call variants with bcftools mpileup
2. Extract 55 SNPs
3. Convert to FROGAncestryCalc format

**For FASTQ:**
1. Align reads with BWA-MEM
2. Sort and index BAM
3. Call variants with bcftools
4. Extract 55 SNPs
5. Convert to FROGAncestryCalc format

**Examples:**

```bash
# From VCF (fastest)
./tools/extract_snps_from_wgs.sh \
    -i patient001.vcf.gz \
    -t vcf \
    -o input/patient001.txt

# From BAM
./tools/extract_snps_from_wgs.sh \
    -i patient001.bam \
    -t bam \
    -r /path/to/GRCh38.fa \
    -o input/patient001.txt \
    -n Patient001

# From FASTQ (paired-end)
./tools/extract_snps_from_wgs.sh \
    -i patient001_R1.fastq.gz \
    -2 patient001_R2.fastq.gz \
    -t fastq \
    -r /path/to/GRCh38.fa \
    -o input/patient001.txt

# Keep intermediate files for inspection
./tools/extract_snps_from_wgs.sh \
    -i sample.bam \
    -t bam \
    -r genome.fa \
    -o input/sample.txt \
    -k
```

---

## 🛠️ Installation

### Using Conda (Recommended)

```bash
# Create environment with all tools
conda create -n genomics \
    python=3.9 \
    bcftools samtools bwa gatk4 \
    -c bioconda -c conda-forge

# Activate environment
conda activate genomics

# Verify installation
bcftools --version
samtools --version
bwa
python3 --version
```

### Manual Installation

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install bcftools samtools bwa python3 wget
```

**macOS (with Homebrew):**
```bash
brew install bcftools samtools bwa python3 wget
```

## 📝 Notes

### Genome Build Compatibility

The FROGAncestryCalc tools support both major genome builds:

- **GRCh37/hg19 (default)**:
  - 1000 Genomes Phase 3 data
  - Coordinates in `SNPInfo/55_aisnps_alleles.txt`
  - Most legacy sequencing data
  
- **GRCh38/hg38**:
  - 1000 Genomes High Coverage data
  - Coordinates in `SNPInfo/55_aisnps_alleles_grch38.txt` (generated via conversion script)
  - Modern sequencing data

**Converting Coordinates:**

To use GRCh38 data, first generate the GRCh38 coordinate file:

```bash
python3 tools/convert_grch37_to_grch38.py
```

This script:
- Queries the Ensembl REST API for each SNP
- Converts GRCh37 coordinates to GRCh38
- Creates `SNPInfo/55_aisnps_alleles_grch38.txt`
- Takes ~5 minutes due to API rate limiting

**Checking Your Data's Build:**

```bash
# Check VCF header
bcftools view -h your_file.vcf.gz | grep "##reference"

# Check BAM header
samtools view -H your_file.bam | grep "@SQ"
```

Common reference indicators:
- GRCh37/hg19: "b37", "hg19", "GRCh37"
- GRCh38/hg38: "b38", "hg38", "GRCh38"

### SNP Coordinates

The 55 AISNPs span multiple chromosomes. Missing SNPs in output may indicate:
- SNP not covered by sequencing
- Low quality/filtered out
- Wrong genome build (use `-b` flag to match your data)
- SNP not in VCF (microarray data may miss some SNPs)

### Performance Tips

- The script automatically caches extracted SNPs for faster subsequent runs
- VCF files are preserved for future use
- Downloads can be resumed if interrupted (the script detects partial files)
- For large datasets, consider parallel processing:
  ```bash
  parallel -j 4 './tools/extract_snps_from_wgs.sh -i {} -t vcf -o input/{/.}.txt' ::: *.vcf.gz
  ```

### Download Issues

If downloads are stalling or failing:

1. **The script will automatically retry** up to 5 times with exponential backoff
2. **Downloads can be resumed** - just run the script again, it will skip completed files
3. **Check your connection:**
   ```bash
   ping ftp.1000genomes.ebi.ac.uk
   ```
4. **Verify disk space:**
   ```bash
   df -h /dados/GENOMICS_DATA/
   ```
5. **Try a different mirror** (if available) by manually downloading and placing files in the target directory
6. **Use curl instead of wget** - the script automatically falls back to curl if wget fails

### Quality Filtering

For critical applications, consider adding quality filters:

```bash
# Filter VCF before conversion
bcftools view -i 'QUAL>=30 && DP>=10' input.vcf.gz | \
bcftools view -i "ID=@tools/aisnps_55_list.txt" -Oz -o filtered.vcf.gz

python3 tools/vcf_to_frog.py filtered.vcf.gz tools/aisnps_55_list.txt input/sample.txt
```

## 🐛 Troubleshooting

**Error: "bcftools: command not found"**
- Install bcftools: `conda install -c bioconda bcftools`

**Error: "No SNPs found in VCF"**
- Check if SNP IDs use rsID format (not chr:pos)
- Verify genome build matches
- Ensure VCF contains these specific SNPs

**Error: "Reference genome not found"**
- Download reference genome (GRCh37 or GRCh38)
- Provide full path with `-r` option

**Warning: "Missing X SNPs"**
- Normal if using targeted sequencing or microarray
- Check coverage of your data
- Some SNPs may need imputation

## 📚 Additional Resources

- [1000 Genomes Project](https://www.internationalgenome.org/)
- [bcftools documentation](http://samtools.github.io/bcftools/)
- [VCF format specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- [UCSC Genome Browser](https://genome.ucsc.edu/)

## 📄 License

These tools are released under the MIT License, consistent with the FROGAncestryCalc project.

