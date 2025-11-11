# Guide: How to Discover Available Tissues/Cells in AlphaGenome

## 🔍 New Feature: `--list-tissues`

The `build_window_and_predict.py` script now includes functionality to list all available tissue/cell ontologies in AlphaGenome!

## 📋 Available Commands

### 1. List ALL tissues/cells

```bash
python3 ~/genomics/build_window_and_predict.py --list-tissues
```

**Expected output:**
```
[INFO] Loading tissue metadata from AlphaGenome (this may take a few seconds)...

================================================================================
Available tissues/cells in AlphaGenome (XXX total)
================================================================================

CURIE                     Biosample Name                                Type           
------------------------- --------------------------------------------- ---------------
CL:0000182                hepatocyte                                    primary cell   
CL:0000540                neuron                                        primary cell   
...
UBERON:0000955            brain                                         tissue         
UBERON:0002048            lung                                          tissue         
UBERON:0002107            liver                                         tissue         
...

================================================================================
Usage: --tissue CURIE (e.g., --tissue UBERON:0002107)
================================================================================
```

### 2. Filter by name (RECOMMENDED)

Search only for tissues/cells containing "brain":

```bash
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue brain
```

Other useful examples:
```bash
# Search for liver
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue liver

# Search for heart
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue heart

# Search for lung
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue lung

# Search for T cells
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue "T cell"

# Search for neuron
python3 ~/genomics/build_window_and_predict.py --list-tissues --filter-tissue neuron
```

### 3. List available output types

```bash
python3 ~/genomics/build_window_and_predict.py --list-outputs
```

**Output:**
```
Available OutputType attributes in AlphaGenome:
  ATAC
  CAGE
  CHIP_HISTONE
  CHIP_TF
  CONTACT_MAPS
  DNASE
  PROCAP
  RNA_SEQ
  SPLICE_JUNCTIONS
  SPLICE_SITES
  SPLICE_SITE_USAGE
```

## 🎯 Using CURIEs in Predictions

After finding the desired CURIE, use it in the prediction command:

```bash
python3 ~/genomics/build_window_and_predict.py \
  --sample HG00096 \
  --gene CYP2B6 \
  --ref-fasta refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --vcf longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --outdir ./alphagenome \
  --predict \
  --outputs "CAGE" \
  --tissue "UBERON:0002107"
```

## 📚 Ontology Types

AlphaGenome uses standardized ontologies:

| Prefix | Name | Description | Examples |
|---------|------|-----------|----------|
| **UBERON** | Uber-anatomy ontology | Anatomy and tissues | UBERON:0002107 (liver), UBERON:0000955 (brain) |
| **CL** | Cell Ontology | Cell types | CL:0000182 (hepatocyte), CL:0000540 (neuron) |
| **CLO** | Cell Line Ontology | Cell lines | CLO:0000001 (HeLa) |
| **EFO** | Experimental Factor Ontology | Experimental factors | EFO:0000001 |
| **NTR** | New Term Requested | New requested terms | NTR:XXX |

## 🌐 Online Ontology Browsers

If you prefer to browse visually:

### UBERON (Tissues/Anatomy)
- **URL**: https://www.ebi.ac.uk/ols/ontologies/uberon
- **Use**: Search for organs and tissues
- **Common examples**:
  - Liver: UBERON:0002107
  - Brain: UBERON:0000955
  - Heart: UBERON:0000948
  - Lung: UBERON:0002048
  - Kidney: UBERON:0002113
  - Pancreas: UBERON:0001264
  - Spleen: UBERON:0002106
  - Blood: UBERON:0000178

### Cell Ontology (Cell Types)
- **URL**: https://www.ebi.ac.uk/ols/ontologies/cl
- **Use**: Search for specific cell types
- **Common examples**:
  - Hepatocyte: CL:0000182
  - Neuron: CL:0000540
  - Cardiomyocyte: CL:0000746
  - T cell: CL:0000084
  - B cell: CL:0000236
  - Macrophage: CL:0000235
  - Fibroblast: CL:0000057

## 💡 Usage Tips

### Tip 1: Always filter by name first
Avoid listing all tissues at once (there are hundreds!). Use `--filter-tissue`:

```bash
# ❌ Bad: Lists everything (hundreds of lines)
python3 build_window_and_predict.py --list-tissues

# ✅ Good: Lists only what matters
python3 build_window_and_predict.py --list-tissues --filter-tissue brain
```

### Tip 2: Save complete list for reference
```bash
python3 ~/genomics/build_window_and_predict.py --list-tissues > tissues_complete.txt
```

### Tip 3: Combine with grep for advanced search
```bash
# Search for tissues related to nervous system
python3 ~/genomics/build_window_and_predict.py --list-tissues | grep -i nerve

# Search for immune system cells
python3 ~/genomics/build_window_and_predict.py --list-tissues | grep -i "immune\|lymph\|T cell"
```

### Tip 4: For multi-tissue analysis
If you need to compare multiple tissues, run the script multiple times with different `--tissue`:

```bash
# Liver
python3 build_window_and_predict.py ... --tissue "UBERON:0002107" --outdir ./results_liver

# Brain
python3 build_window_and_predict.py ... --tissue "UBERON:0000955" --outdir ./results_brain

# Heart
python3 build_window_and_predict.py ... --tissue "UBERON:0000948" --outdir ./results_heart
```

## ⚠️ Important: API Key Required

The `--list-tissues` option requires a valid AlphaGenome API key because it needs to connect to the server to fetch metadata.

**Ways to provide the API key:**

1. Via environment variable (recommended):
   ```bash
   export ALPHAGENOME_API_KEY="your-key-here"
   python3 build_window_and_predict.py --list-tissues
   ```

2. Via argument:
   ```bash
   python3 build_window_and_predict.py --list-tissues --api-key "your-key-here"
   ```

## 📊 Statistics by Output Type

Number of biosamples and tracks per output type:

| Output Type | Unique Biosamples | Total Tracks |
|-------------|-------------------|-----------------|
| RNA_SEQ | 285 | 667 |
| CAGE | 264 | 546 |
| DNASE | 305 | 305 |
| ATAC | 167 | 167 |
| CHIP_HISTONE | 219 | 1116 |
| PROCAP | 6 | 12 |

## 🚀 Complete Workflow

```bash
# 1. View available outputs
python3 build_window_and_predict.py --list-outputs

# 2. Search for tissue of interest
python3 build_window_and_predict.py --list-tissues --filter-tissue liver

# 3. Copy the desired CURIE (e.g., UBERON:0002107)

# 4. Run prediction
python3 build_window_and_predict.py \
  --sample HG00096 \
  --gene CYP2B6 \
  --ref-fasta refs/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --vcf longevity_dataset/vcf_chromosomes/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  --outdir ./alphagenome \
  --predict \
  --outputs "CAGE" \
  --tissue "UBERON:0002107"

# 5. Analyze results
python3 read_alphagenome_predictions.py \
  alphagenome/HG00096__CYP2B6/predictions_H1/cage.npz
```

## 📝 References

- [AlphaGenome Documentation](https://alphafold.com/alphagenome)
- [UBERON Ontology Browser](https://www.ebi.ac.uk/ols/ontologies/uberon)
- [Cell Ontology Browser](https://www.ebi.ac.uk/ols/ontologies/cl)
- [Ontology Lookup Service (OLS)](https://www.ebi.ac.uk/ols/)

