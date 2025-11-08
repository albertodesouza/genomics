# 🔗 Neural Module + Genomes Analyzer Integration

## 📖 Overview

**Neural Integration** (`neural_integration.py`) is a bridge tool that connects the traditional genomic analysis pipeline (`genomes_analyzer.py`) with AI-based neural analysis (`neural_module.py` + AlphaGenome).

It automates the workflow:
**VCF/BED/Genes → Sequence Extraction → Neural Analysis → Results Correlation**

---

## 🎯 What Does Neural Integration Do?

### 1. **Intelligent Sequence Extraction**
- ✅ Extract genomic regions from VCF (variants)
- ✅ Extract regions from BED files (regions of interest)
- ✅ Extract specific genes from GTF (with flanking regions)
- ✅ Convert everything to FASTA ready for AlphaGenome

### 2. **Automated Neural Analysis**
- ✅ Automatically run `neural_module.py`
- ✅ Configure appropriate parameters
- ✅ Manage AlphaGenome outputs

### 3. **Results Correlation**
- ✅ Correlate variants with neural predictions
- ✅ Generate integration reports
- ✅ Create combined visualizations

---

## 🚀 Operation Modes

`neural_integration.py` operates in 4 different modes:

### Mode 1: **Complete Integrated Analysis** 🌟

Automatically executes the entire workflow: extraction + analysis + correlation.

```bash
cd neural_module
python neural_integration.py \
  --integrated \
  --vcf ../vcf/NA12878.vcf.gz \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output integrated_results/
```

**When to use**: After running `genomes_analyzer.py` and having a variant VCF.

**What happens**:
1. Extracts regions around each variant (±5kb)
2. Converts to FASTA
3. Runs neural analysis with AlphaGenome
4. Correlates variants with predictions
5. Generates integrated reports

---

### Mode 2: **Extract Sequences from VCF**

Only extracts sequences without running neural analysis.

```bash
python neural_integration.py \
  --extract-vcf \
  --vcf ../vcf/NA12878.vcf.gz \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --output variants_sequences.fasta
```

**When to use**: When you only want to prepare sequences for later analysis.

**Output**: FASTA file with ±5kb regions around each variant.

---

### Mode 3: **Extract Sequences from BED**

Extracts sequences from regions specified in BED file.

```bash
python neural_integration.py \
  --extract-bed \
  --bed regions_of_interest.bed \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --output regions_sequences.fasta
```

**When to use**: When you have specific regions of interest (enhancers, promoters, etc.).

**Expected BED format**:
```
chr1    1000000    1002048    region_1
chr2    5000000    5002048    region_2
```

---

### Mode 4: **Extract Specific Genes**

Extracts gene sequences by name, with flanking regions.

```bash
python neural_integration.py \
  --extract-genes \
  --genes BRCA1 TP53 HBB CFTR \
  --gtf ../refs/gencode.v38.annotation.gtf.gz \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --output genes_sequences.fasta \
  --flank 10000
```

**When to use**: To analyze specific genes of interest with regulatory context.

**Parameters**:
- `--flank`: Flanking bases (default: 10kb before and after gene)

---

## 📊 Practical Use Cases

### Case 1: Analyze High-Impact Variants

After identifying high-impact variants in VCF:

```bash
# Step 1: Filter high-impact variants (example)
bcftools view -i 'INFO/ANN~"HIGH"' ../vcf/sample.vcf.gz > high_impact.vcf

# Step 2: Integrated analysis
cd neural_module
python neural_integration.py \
  --integrated \
  --vcf high_impact.vcf \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output neural_high_impact/ \
  --outputs RNA_SEQ ATAC CHIP_HISTONE
```

**Result**: Neural predictions for regions with high-impact variants.

---

### Case 2: Analyze Candidate Genes

You've identified disease candidate genes:

```bash
cd neural_module
python neural_integration.py \
  --extract-genes \
  --genes BRCA1 BRCA2 TP53 PTEN \
  --gtf ../refs/gencode.v38.annotation.gtf.gz \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --output candidate_genes.fasta \
  --flank 20000

# Then analyze with neural_module
python neural_module.py \
  -i candidate_genes.fasta \
  -k YOUR_API_KEY \
  -o candidate_genes_neural/
```

**Result**: Complete functional analysis of candidate genes.

---

### Case 3: Non-Coding Regulatory Regions

You have regulatory regions of interest in BED:

```bash
# enhancers.bed contains enhancer regions
cd neural_module
python neural_integration.py \
  --extract-bed \
  --bed enhancers.bed \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --output enhancers.fasta

# Analyze accessibility and epigenetic markers
python neural_module.py \
  -i enhancers.fasta \
  -k YOUR_API_KEY \
  -o enhancers_analysis/ \
  --outputs ATAC DNASE CHIP_HISTONE CHIP_TF
```

**Result**: Regulatory activity predictions.

---

### Case 4: Trio Analysis (De Novo Variants)

After identifying de novo variants in trio:

```bash
# Step 1: Get de novo variants from pipeline
# (genomes_analyzer.py already generates trio/denovo_candidates.vcf)

# Step 2: Neural analysis of de novo variants
cd neural_module
python neural_integration.py \
  --integrated \
  --vcf ../trio/denovo_candidates.vcf \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output neural_denovo/ \
  --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF
```

**Result**: Predicted functional impact of de novo variants.

---

## 🔄 Complete Integrated Workflow

### Complete Pipeline: DNA → Variants → Neural Predictions

```bash
# ════════════════════════════════════════════════════════════════
# STEP 1: Traditional Genomic Analysis
# ════════════════════════════════════════════════════════════════

conda activate genomics
python genomes_analyzer.py --config config_human_30x.yaml

# Outputs:
# - vcf/NA12878.vcf.gz (variants)
# - trio/denovo_candidates.vcf (de novo)
# - bam/*.bam (alignments)

# ════════════════════════════════════════════════════════════════
# STEP 2: Filter Variants of Interest (Optional)
# ════════════════════════════════════════════════════════════════

# Example: Exonic high-impact variants
bcftools view -i 'INFO/ANN~"HIGH|MODERATE" && INFO/ANN~"exonic"' \
  vcf/NA12878.vcf.gz > variants_of_interest.vcf

# ════════════════════════════════════════════════════════════════
# STEP 3: Integrated Neural Analysis
# ════════════════════════════════════════════════════════════════

cd neural_module
python neural_integration.py \
  --integrated \
  --vcf ../variants_of_interest.vcf \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_ALPHAGENOME_KEY \
  --output integrated_analysis/ \
  --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF

# Outputs:
# - integrated_analysis/neural_results/ (AlphaGenome predictions)
# - integrated_analysis/correlation_report.json (correlation)

# ════════════════════════════════════════════════════════════════
# STEP 4: Interpret Results
# ════════════════════════════════════════════════════════════════

# View correlation report
cat integrated_analysis/correlation_report.json | jq .

# View neural visualizations
ls integrated_analysis/neural_results/*.png

# View ontology metadata
cat integrated_analysis/neural_results/*_metadata.csv
```

---

## 📁 Output Structure

### Integrated Mode (`--integrated`)

```
integrated_analysis/
├── variants_sequences.fasta          # Sequences extracted from VCF
├── neural_results/                   # neural_module results
│   ├── variant_1_*_RNA_SEQ.png
│   ├── variant_1_*_RNA_SEQ_enhanced.png
│   ├── variant_1_*_RNA_SEQ_heatmap.png
│   ├── variant_1_*_RNA_SEQ_metadata.csv
│   ├── variant_1_*_RNA_SEQ_metadata.json
│   ├── ... (other outputs)
│   ├── variant_1_*_comparison.png
│   ├── variant_1_*_dashboard.png
│   └── analysis_report.json
└── correlation_report.json           # Variants × predictions correlation
```

### Extraction Mode (`--extract-*`)

```
sequences.fasta                       # Extracted sequences ready for use
```

---

## 💡 Advanced Examples

### Example 1: Chromosome-Specific Analysis

```bash
# Analyze only chr11 variants (HBB gene)
cd neural_module
python neural_integration.py \
  --integrated \
  --vcf <(bcftools view -r chr11 ../vcf/sample.vcf.gz) \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output chr11_neural/
```

### Example 2: Metabolic Pathway Genes

```bash
# Extract all genes from a pathway (example: DNA repair)
cd neural_module
python neural_integration.py \
  --extract-genes \
  --genes BRCA1 BRCA2 ATM CHEK2 TP53 PALB2 RAD51 \
  --gtf ../refs/gencode.v38.annotation.gtf.gz \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --output dna_repair_genes.fasta \
  --flank 15000
```

### Example 3: Variant Prioritization

```bash
#!/bin/bash
# Script to prioritize variants with neural analysis

# 1. Rare high-impact variants
bcftools view -i 'INFO/AF<0.01 && INFO/ANN~"HIGH"' \
  vcf/sample.vcf.gz > rare_high_impact.vcf

# 2. Neural analysis
cd neural_module
python neural_integration.py \
  --integrated \
  --vcf rare_high_impact.vcf \
  --ref ../refs/GRCh38.d1.vd1.fa \
  --api-key $ALPHAGENOME_KEY \
  --output prioritized_variants/ \
  --outputs RNA_SEQ ATAC CHIP_HISTONE

# 3. View results
python -c "
import json
with open('prioritized_variants/correlation_report.json') as f:
    data = json.load(f)
    print(f'Variants analyzed: {data[\"summary\"][\"total_sequences\"]}')
    print(f'Successful predictions: {data[\"summary\"][\"successful_predictions\"]}')
"
```

---

## ⚙️ Advanced Configuration

### Customize Window Size

By default, extracts ±5kb around variants. To change, edit `neural_integration.py`:

```python
# Line ~76
start = max(1, pos - 10000)  # Was 5000, now 10kb
end = pos + 10000
```

### Add Quality Filters

Filter variants before neural analysis:

```bash
# Only PASS variants with DP ≥ 20 and GQ ≥ 30
bcftools view -f PASS -i 'FORMAT/DP>=20 && FORMAT/GQ>=30' \
  vcf/sample.vcf.gz | \
python neural_module/neural_integration.py \
  --integrated \
  --vcf /dev/stdin \
  --ref refs/GRCh38.d1.vd1.fa \
  --api-key YOUR_API_KEY \
  --output high_quality_neural/
```

---

## 🔍 Interpreting Results

### Correlation Report

The `correlation_report.json` file contains:

```json
{
  "vcf_source": "vcf/sample.vcf.gz",
  "neural_results": {
    "timestamp": "2025-10-18T...",
    "total_sequences": 42,
    "successful_analyses": 40,
    "sequences": [
      {
        "id": "variant_1_chr11_5227002",
        "length": 10000,
        "status": "success",
        "outputs": ["RNA_SEQ", "CAGE", "ATAC", "CHIP_HISTONE", "CHIP_TF"]
      },
      ...
    ]
  },
  "summary": {
    "total_sequences": 42,
    "successful_predictions": 40
  }
}
```

### Results Analysis

1. **Individual Visualizations**: See `neural_results/*_enhanced.png` for each variant
2. **Heatmaps**: Compare multiple tracks in `*_heatmap.png`
3. **Dashboard**: Statistical summary in `*_dashboard.png`
4. **Metadata**: Tissue/cell information in `*_metadata.csv`

---

## ❓ FAQ

### Q: Does neural_integration require genomes_analyzer installed?

**A**: No! `neural_integration.py` is independent. It only requires:
- `bcftools` (for VCF)
- `bedtools` (for BED)
- `samtools` (for sequence extraction)
- `neural_module.py` (for neural analysis)

All are already in the `genomics` environment.

### Q: Can I use with VCFs from other pipelines?

**A**: Yes! Works with any standard VCF, not just those generated by `genomes_analyzer.py`.

### Q: How much does it cost to use with AlphaGenome?

**A**: AlphaGenome is free for non-commercial use. See https://www.alphagenomedocs.com/

### Q: Can I analyze only specific variants?

**A**: Yes! Use `bcftools view` to filter the VCF first:

```bash
# Only variants at specific positions
bcftools view -t chr11:5227002 vcf/sample.vcf.gz | \
python neural_module/neural_integration.py --integrated --vcf /dev/stdin ...
```

### Q: How to add variant analysis (REF vs ALT)?

**A**: Use `neural_module.py` directly after extraction:

```bash
# 1. Extract sequences
cd neural_module
python neural_integration.py \
  --extract-vcf \
  --vcf ../variants.vcf \
  --ref ../genome.fa \
  --output sequences.fasta

# 2. Analyze each variant with --variant
# (requires additional script to parse VCF)
```

---

## 🔗 Related Resources

- **[Main Neural Module](NEURAL_MODULE.md)** - Complete documentation
- **[Usage Guide](USAGE.md)** - How to use neural_module.py
- **[Download Sequences](../../docs/DOWNLOAD_SEQUENCES.md)** - Download real genomes
- **[Results Interpretation](RESULTS.md)** - Understanding predictions

---

## 🚀 Next Steps

After mastering basic integration:

1. **Automate**: Create scripts for batch analysis
2. **Prioritize**: Combine variant scores with neural predictions
3. **Validate**: Compare predictions with experimental data (if available)
4. **Publish**: Include neural analyses in your reports

---

**Created**: October 2025  
**Translated**: November 2025  
**Part of**: Neural Module Documentation
