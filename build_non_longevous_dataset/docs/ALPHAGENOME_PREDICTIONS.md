# AlphaGenome Predictions - Usage Guide

## 📋 Overview

The `build_window_and_predict.py` script now saves complete AlphaGenome predictions as NumPy arrays, enabling detailed analysis of prediction data (ATAC-seq, RNA-seq, etc.) for each nucleotide in the sequence.

## 🚀 Running Predictions

### Basic example:

```bash
# From build_non_longevous_dataset directory
python3 build_window_and_predict.py \
  --sample HG00096 \
  --gene CYP2B6 \
  --ref-fasta ../refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --vcf ../longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --outdir ./alphagenome_output \
  --predict \
  --outputs "ATAC" \
  --ontology "UBERON:0002107"
```

### View available outputs:

```bash
python3 build_window_and_predict.py --list-outputs
```

Example output:
```
Available OutputType attributes in AlphaGenome:
  ATAC
  CAGE
  DNASE
  H3K27AC
  H3K27ME3
  H3K36ME3
  H3K4ME1
  H3K4ME3
  H3K9ME3
  RNA
  ...
```

## 📁 Output Structure

After execution, files are organized as follows:

```
alphagenome/HG00096__CYP2B6/
├── gtf_cache.feather                         # GTF cache (reused)
├── ref.window.fa                             # Reference sequence (1 Mb)
├── HG00096.window.vcf.gz                     # Sample variants in region
├── HG00096.window.consensus_ready.vcf.gz     # Filtered variants
├── HG00096.H1.window.raw.fa                  # H1 consensus (before adjustment)
├── HG00096.H1.window.fixed.fa                # H1 consensus (exactly 1 million bases)
├── HG00096.H2.window.raw.fa                  # H2 consensus (before adjustment)
├── HG00096.H2.window.fixed.fa                # H2 consensus (exactly 1 million bases)
├── predictions_H1/                           # ⭐ AlphaGenome predictions for H1
│   ├── atac.npz                              #    NumPy arrays (1M values)
│   └── atac_metadata.json                    #    Track metadata
├── predictions_H2/                           # ⭐ AlphaGenome predictions for H2
│   ├── atac.npz
│   └── atac_metadata.json
├── prediction_H1.ok.txt                      # H1 completion marker
└── prediction_H2.ok.txt                      # H2 completion marker
```

## 📊 Analyzing the Results

### Included analysis script:

```bash
# Basic file analysis
python3 ~/genomics/read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz

# Generate plot for a region
python3 ~/genomics/read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz \
  --plot --start 0 --end 10000 --output atac_plot.png

# Compare H1 vs H2 haplotypes
python3 ~/genomics/read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz \
  --compare alphagenome/HG00096__CYP2B6/predictions_H2/atac.npz
```

### Python example:

```python
import numpy as np
import json
from pathlib import Path

# Load H1 predictions
data_h1 = np.load('alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz')

# View available tracks
print(f"Tracks: {data_h1.files}")  # Ex: ['track_0', 'track_1', ...]

# Access specific track
track_0 = data_h1['track_0']  # Array with ~1 million values

# Basic statistics
print(f"Shape: {track_0.shape}")
print(f"Mean:  {track_0.mean():.6f}")
print(f"Std:   {track_0.std():.6f}")
print(f"Min:   {track_0.min():.6f}")
print(f"Max:   {track_0.max():.6f}")

# Load metadata
with open('alphagenome/HG00096__CYP2B6/predictions_H1/atac_metadata.json') as f:
    metadata = json.load(f)
    
print(f"Metadata: {metadata}")

# Analyze specific region (e.g., first 1000 nucleotides)
region = track_0[0:1000]
print(f"Mean in region 0-1000: {region.mean():.6f}")

# Compare H1 vs H2
data_h2 = np.load('alphagenome/HG00096__CYP2B6/predictions_H2/atac.npz')
track_h2 = data_h2['track_0']

# Absolute difference
diff = np.abs(track_0 - track_h2)
print(f"Mean difference between H1 and H2: {diff.mean():.6f}")
print(f"Positions with difference > 0.1: {(diff > 0.1).sum()}")

# Save processed results
np.save('difference_h1_h2.npy', diff)
```

## 🧬 Common Tissue CURIEs

For use with `--tissue`:

| CURIE | Tissue/Cell |
|-------|---------------|
| `UBERON:0002107` | Liver |
| `UBERON:0000955` | Brain |
| `UBERON:0000948` | Heart |
| `UBERON:0002048` | Lung |
| `UBERON:0001264` | Pancreas |
| `CL:0000182` | Hepatocyte |
| `CL:0000540` | Neuron |
| `CL:0000746` | Cardiomyocyte |

If you don't specify `--tissue` or use an invalid value, AlphaGenome returns predictions for **all** available tissues/cells.

## 🔄 Idempotence

The script is completely idempotent:

- ✅ GTF cache (reused across all genes)
- ✅ FASTA sequences (skipped if already exist)
- ✅ Processed VCFs (skipped if already exist)
- ✅ Predictions (skipped if `.ok.txt` markers exist)

You can run the same command multiple times and only incomplete steps will be executed.

## ⚡ Performance

### First execution (without cache):
- GTF download: ~10-30 seconds
- Reference extraction: ~1-2 seconds
- VCF subset: ~2-5 seconds
- Consensus (H1+H2): ~3-5 seconds
- AlphaGenome predictions: ~30-60 seconds (depends on API)
- **Total: ~1-2 minutes**

### Subsequent executions (with cache):
- GTF loading: ~0.5 seconds
- Skip all completed steps
- **Total: ~1 second** (if everything already exists)

## 💾 Disk Space

Per case (sample + gene):

- FASTA sequences: ~3-5 MB
- VCFs: ~0.5-2 MB (depends on number of variants)
- NPZ predictions (compressed): ~8-20 MB per output type per haplotype
- **Estimated total**: ~20-50 MB per case

The GTF cache (~50-100 MB) is shared across all cases.

## 🎯 Use Cases

### 1. Functional impact analysis
Compare predictions between haplotypes to identify variants with functional effects:

```python
diff = np.abs(h1_predictions - h2_predictions)
high_impact_positions = np.where(diff > threshold)[0]
```

### 2. Gene epigenetic profile
Analyze chromatin profile (ATAC, H3K27ac, etc.) around a gene of interest.

### 3. Tissue-specific effects
Compare predictions using different `--tissue` to see if variants have tissue-specific effects.

### 4. Population analysis
Run for multiple samples (e.g., 1000 Genomes) and compare profiles across populations.

## 📚 References

- [AlphaGenome Documentation](https://alphafold.com/alphagenome)
- [UBERON Ontology Browser](https://www.ebi.ac.uk/ols/ontologies/uberon)
- [Cell Ontology (CL)](https://www.ebi.ac.uk/ols/ontologies/cl)

## 🐛 Troubleshooting

### Error: "Invalid ontology_curie"
Use CURIEs in the format `TYPE:ID` (e.g., `UBERON:0002107`), not free text.

### Error: "Output type not found"
Use `--list-outputs` to see valid names. Use the exact name (e.g., `ATAC`, not `ATAC-seq`).

### Empty arrays or None
Some output/tissue combinations may not have data. Check warnings in the log.

### Out of memory
Predictions use ~2-4 GB of RAM. For multiple samples, process sequentially.

