# AISNP Mode Documentation

## Overview

**AISNP Mode** is a specialized operating mode of `build_window_and_predict.py` designed for analyzing ancestry-informative SNPs (AISNPs) from the FROGAncestryCalc module. Instead of centering genomic windows on gene bodies, AISNP mode centers 1 Mb windows on specific SNP positions.

## What are AISNPs?

**Ancestry Informative SNPs (AISNPs)** are genetic variants that show significant frequency differences between populations from different ancestral backgrounds. These SNPs are particularly useful for:

- Inferring individual ancestry composition
- Population stratification correction
- Forensic identification
- Understanding population-specific genetic adaptations

The FROGAncestryCalc module uses a panel of **55 validated AISNPs** from the Kidd Lab (ALFRED database) to infer ancestry across 155 worldwide populations.

## Why Use AISNP Mode?

AISNP mode enables you to:

1. **Analyze functional context** around ancestry markers
2. **Compare epigenetic profiles** between ancestries at AISNP loci
3. **Identify regulatory elements** that may contribute to population differences
4. **Study gene expression differences** near ancestry-informative positions
5. **Generate AlphaGenome predictions** for tissue-specific functional effects

## Integration with FROGAncestryCalc

This mode is designed to work seamlessly with FROGAncestryCalc:

```
FROGAncestryCalc Module                     build_non_longevous_dataset Module
┌─────────────────────────┐                ┌────────────────────────────────┐
│ 1. Extract 55 AISNPs    │                │ 4. Extract 1 Mb windows        │
│    from 1000G VCFs      │───────────────>│    around each AISNP           │
│                         │                │                                │
│ 2. Calculate ancestry   │                │ 5. Apply individual variants   │
│    likelihoods          │                │                                │
│                         │                │ 6. Run AlphaGenome predictions │
│ 3. Identify population  │                │                                │
│    clusters             │<───────────────│ 7. Compare functional profiles │
│                         │                │    between ancestries          │
└─────────────────────────┘                └────────────────────────────────┘
```

## Configuration

### Basic AISNP Mode Configuration

```yaml
build_window_params:
  # Enable SNP mode
  mode: "snp"
  
  # Specify SNP list file (GRCh38 coordinates)
  snp:
    snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"
  
  # Window size (centered on SNP position)
  window_size: 1000000     # 1 Mb (500kb upstream + 500kb downstream)
  
  # AlphaGenome predictions
  predict: true
  outputs: "ATAC,CHIP_HISTONE,RNA_SEQ"
  ontology: "UBERON:0002107"  # liver (or other tissues)
  
  # Haplotype options
  skip_h2: false              # Keep both haplotypes
  also_iupac: false
```

### Sample Selection for Ancestry Analysis

Select samples from different superpopulations to compare:

```yaml
sample_selection:
  level: "superpopulation"
  samples_per_group: 10        # 10 samples per ancestry
  include_groups: ["AFR", "EUR", "EAS", "SAS", "AMR"]
  sex_filter: "all"
```

## SNP File Format

The SNP list file must be tab-delimited with the following format:

```
ALFRED_UID	dbSNP_rsnumber	chrom	chrom_pos	alleles
SI047925B	rs10497191	2	157810705	C/T
SI000148N	rs1079597	11	113425564	C/T
SI014486X	rs11652805	17	64991033	C/T
```

**Columns:**
- `ALFRED_UID`: Unique identifier from ALFRED database
- `dbSNP_rsnumber`: rsID (e.g., rs10497191)
- `chrom`: Chromosome number (1-22, X, Y)
- `chrom_pos`: Position on chromosome (1-based, GRCh38)
- `alleles`: Reference/alternate alleles (e.g., C/T)

**Important:** The SNP file must use **GRCh38 coordinates** to match the 1000 Genomes High Coverage data.

## Output Structure

For each sample and each AISNP, a separate directory is created:

```
non_longevous_results/
├── HG00096__rs10497191/          # Sample HG00096, SNP rs10497191
│   ├── ref.window.fa             # Reference sequence (1 Mb)
│   ├── HG00096.H1.window.fixed.fa  # Haplotype 1 with variants applied
│   ├── HG00096.H2.window.fixed.fa  # Haplotype 2 with variants applied
│   ├── HG00096.window.vcf.gz      # Variants in this region
│   ├── predictions_H1/
│   │   ├── atac.npz              # ATAC-seq predictions
│   │   ├── atac_metadata.json
│   │   ├── chip_histone.npz      # Histone ChIP-seq predictions
│   │   └── chip_histone_metadata.json
│   └── predictions_H2/
│       └── ...
│
├── HG00096__rs1079597/           # Next SNP for same sample
│   └── ...
│
├── HG00097__rs10497191/          # Next sample, first SNP
│   └── ...
└── ...
```

**Total directories** for N samples and 55 AISNPs: N × 55

**Example:** 100 samples = 5,500 directories

## Use Cases

### 1. Compare Chromatin Accessibility Between Ancestries

**Goal:** Determine if ATAC-seq profiles differ between African and European ancestries at AISNP loci.

**Configuration:**
```yaml
sample_selection:
  include_groups: ["AFR", "EUR"]
  samples_per_group: 20

build_window_params:
  mode: "snp"
  snp:
    snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"
  predict: true
  outputs: "ATAC"
  ontology: "UBERON:0002107"  # liver
```

**Analysis:**
1. Extract ATAC predictions for each sample at each AISNP
2. Group by ancestry
3. Compare mean accessibility profiles
4. Identify AISNPs with significant differences

### 2. Study Regulatory Elements Near AISNPs

**Goal:** Identify histone modifications and transcription factor binding sites near ancestry markers.

**Configuration:**
```yaml
build_window_params:
  mode: "snp"
  predict: true
  outputs: "CHIP_HISTONE,CHIP_TF"
  ontology: "CL:0002601"  # senescent cell (longevity-related)
```

**Analysis:**
1. Extract ChIP-seq predictions
2. Identify peaks near AISNP positions
3. Annotate regulatory regions
4. Compare between ancestries

### 3. RNA Expression Differences Near AISNPs

**Goal:** Determine if genes near AISNPs show differential expression between ancestries.

**Configuration:**
```yaml
build_window_params:
  mode: "snp"
  predict: true
  outputs: "RNA_SEQ,CAGE"
  ontology: "UBERON:0000955,UBERON:0002107"  # brain, liver
```

### 4. Focused AISNP Subset Analysis

If you want to analyze only specific AISNPs (not all 55), create a custom SNP file:

**custom_aisnps.txt:**
```
ALFRED_UID	dbSNP_rsnumber	chrom	chrom_pos	alleles
SI047925B	rs10497191	2	157810705	C/T
SI000148N	rs1079597	11	113425564	C/T
SI014486X	rs11652805	17	64991033	C/T
```

**Configuration:**
```yaml
snp:
  snp_list_file: "custom_aisnps.txt"
```

## Comparison: AISNP Mode vs Gene Mode

| Feature | Gene Mode | AISNP Mode |
|---------|-----------|------------|
| **Target** | Gene bodies | SNP positions |
| **Window centering** | Gene TSS or midpoint | SNP position (1-based) |
| **GTF required** | Yes | No |
| **Number of targets** | 1-100s | Typically 55 (AISNPs) |
| **Best for** | Functional genomics | Ancestry analysis |
| **Output naming** | `SAMPLE__GENE` | `SAMPLE__rsID` |
| **Integration** | General | FROGAncestryCalc |

## Workflow Example

### Complete AISNP Analysis Workflow

```bash
# ═══════════════════════════════════════════════════════════════
# Step 1: Prepare FROGAncestryCalc SNP coordinates (GRCh38)
# ═══════════════════════════════════════════════════════════════
cd FROGAncestryCalc

# Generate GRCh38 coordinates if not already done
python3 tools/convert_grch37_to_grch38.py

# Verify file exists
ls -lh SNPInfo/55_aisnps_alleles_grch38.txt

# ═══════════════════════════════════════════════════════════════
# Step 2: Extract AISNP genotypes from 1000 Genomes
# ═══════════════════════════════════════════════════════════════
./tools/extract_snps_from_1000genomes.sh -b grch38

# ═══════════════════════════════════════════════════════════════
# Step 3: Run ancestry inference
# ═══════════════════════════════════════════════════════════════
./run.sh

# Review ancestry results
cat output/1000genomes_55aisnps_likelihood.txt

# ═══════════════════════════════════════════════════════════════
# Step 4: Configure build_non_longevous_dataset for AISNP mode
# ═══════════════════════════════════════════════════════════════
cd ../build_non_longevous_dataset

# Edit configs/default.yaml
nano configs/default.yaml

# Set:
#   mode: "snp"
#   snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"

# ═══════════════════════════════════════════════════════════════
# Step 5: Run AISNP window extraction and predictions
# ═══════════════════════════════════════════════════════════════
python3 build_non_longevous_dataset.py --config configs/default.yaml

# ═══════════════════════════════════════════════════════════════
# Step 6: Analyze results
# ═══════════════════════════════════════════════════════════════
# Load predictions and compare between ancestries
python3 analyze_aisnp_predictions.py
```

## Performance Considerations

### Disk Space

- **Per sample, per AISNP**: ~50-500 MB (depending on predictions)
- **55 AISNPs × 100 samples**: ~275 GB - 2.75 TB
- **Recommendation**: Start with fewer samples or disable predictions for initial tests

### Processing Time

- **Without predictions**: ~30 seconds per AISNP per sample
- **With predictions** (1 tissue): ~2-3 minutes per AISNP per sample
- **55 AISNPs × 100 samples × 2 min**: ~183 hours (7.6 days)
- **Recommendation**: Use parallelization (`n_workers: 8`) and checkpoint system

### Optimization Tips

1. **Disable H2**: `skip_h2: true` → 2x faster
2. **Fewer tissues**: Use 1-2 specific ontologies instead of all
3. **Subset AISNPs**: Analyze 10-20 most informative AISNPs first
4. **Parallel processing**: Increase `n_workers` based on available CPUs
5. **Checkpoint system**: Pipeline resumes automatically if interrupted

## Troubleshooting

### Issue: SNP file not found

```
FileNotFoundError: [Errno 2] No such file or directory: '55_aisnps_alleles_grch38.txt'
```

**Solution:**
```bash
cd FROGAncestryCalc
python3 tools/convert_grch37_to_grch38.py
```

### Issue: VCF doesn't contain chromosome

```
[ERROR] No variants found in region chr2:157310705-158310705
```

**Solution:** Ensure VCF pattern includes all chromosomes where AISNPs are located. AISNPs span multiple chromosomes (2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 17, 19, 22, X).

### Issue: Too many output directories

```
OSError: [Errno 28] No space left on device
```

**Solution:**
- Reduce `samples_per_group`
- Disable predictions: `predict: false`
- Use `skip_h2: true`
- Analyze subset of AISNPs

## Advanced Topics

### Custom AISNP Panels

You can use your own AISNP panel by creating a compatible SNP file:

```python
# create_custom_aisnps.py
import pandas as pd

# Your custom SNPs (GRCh38 coordinates)
snps = [
    {"ALFRED_UID": "CUSTOM001", "dbSNP_rsnumber": "rs123456", 
     "chrom": "1", "chrom_pos": 12345678, "alleles": "A/G"},
    # ... more SNPs
]

df = pd.DataFrame(snps)
df.to_csv("custom_aisnps.txt", sep="\t", index=False)
```

### Downstream Analysis

After generating predictions, you can:

1. **Load predictions into Python:**
```python
import numpy as np
import json

# Load ATAC predictions
data = np.load("HG00096__rs10497191/predictions_H1/atac.npz")
atac_values = data['values']  # Shape: (num_tracks, 1000000)

# Load metadata
with open("HG00096__rs10497191/predictions_H1/atac_metadata.json") as f:
    metadata = json.load(f)
    track_names = metadata['track_names']
```

2. **Compare between ancestries:**
```python
import matplotlib.pyplot as plt

# Plot mean ATAC profile around AISNP
afr_samples = ["HG01879", "HG01880", ...]  # African
eur_samples = ["HG00096", "HG00097", ...]  # European

aisnp = "rs10497191"
afr_profiles = [load_atac(s, aisnp) for s in afr_samples]
eur_profiles = [load_atac(s, aisnp) for s in eur_samples]

plt.plot(np.mean(afr_profiles, axis=0), label="AFR")
plt.plot(np.mean(eur_profiles, axis=0), label="EUR")
plt.legend()
plt.title(f"ATAC profile at {aisnp}")
plt.show()
```

## References

- **FROGAncestryCalc**: See `../FROGAncestryCalc/README.md`
- **55 AISNPs Panel**: Kidd Lab ALFRED database
- **AlphaGenome**: DeepMind functional genomics predictions
- **1000 Genomes High Coverage**: GRCh38-aligned whole genome sequences

## See Also

- [Main README](../README.md)
- [Quick Start Guide](../QUICKSTART.md)
- [AlphaGenome Predictions Guide](ALPHAGENOME_PREDICTIONS.md)
- [FROGAncestryCalc Documentation](../../FROGAncestryCalc/README.md)

