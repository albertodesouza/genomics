# 📊 Available Outputs in AlphaGenome

## Complete Output List

AlphaGenome offers **11 different output types** for DNA analysis:

### 🧬 Gene Expression

#### 1. `RNA_SEQ`
- **Description**: Gene expression level prediction via RNA-seq
- **Use**: Identify active genes and their expression levels
- **Resolution**: Single base

#### 2. `CAGE`
- **Description**: Cap Analysis of Gene Expression
- **Use**: Identify transcription start sites (promoters)
- **Resolution**: Single base

#### 3. `PROCAP`
- **Description**: Precision Run-On sequencing (PRO-cap)
- **Use**: Map transcription start sites with high precision
- **Resolution**: Single base

### 🔬 Chromatin Accessibility

#### 4. `ATAC`
- **Description**: Assay for Transposase-Accessible Chromatin
- **Use**: Identify open/accessible chromatin regions
- **Resolution**: Single base
- **Application**: Find active regulatory elements

#### 5. `DNASE`
- **Description**: DNase I hypersensitivity
- **Use**: Identify accessible chromatin regions
- **Resolution**: Single base
- **Similar to**: ATAC-seq

### ⚛️ Chromatin Modifications

#### 6. `CHIP_HISTONE`
- **Description**: ChIP-seq for histone markers
- **Includes**:
  - H3K27AC - Active enhancers
  - H3K4ME3 - Active promoters
  - H3K27ME3 - Repression (Polycomb)
  - H3K36ME3 - Active gene bodies
  - H3K9ME3 - Heterochromatin
  - And other markers
- **Use**: Identify epigenetic state of chromatin
- **Resolution**: Single base

#### 7. `CHIP_TF`
- **Description**: ChIP-seq for transcription factors
- **Includes**:
  - CTCF - Insulators and loops
  - And other transcription factors
- **Use**: Identify TF binding sites
- **Resolution**: Single base

### 🧩 3D Structure and Splicing

#### 8. `CONTACT_MAPS`
- **Description**: Chromatin contact maps (Hi-C like)
- **Use**: Predict 3D chromatin interactions
- **Application**: Understand nuclear organization

#### 9. `SPLICE_JUNCTIONS`
- **Description**: Splicing junction prediction
- **Use**: Identify exons and introns
- **Application**: Study alternative splicing

#### 10. `SPLICE_SITES`
- **Description**: Splicing site prediction (5' and 3')
- **Use**: Identify splice donor/acceptor sites
- **Resolution**: Single base

#### 11. `SPLICE_SITE_USAGE`
- **Description**: Splicing site usage prediction
- **Use**: Quantify alternative site usage
- **Application**: Study isoforms

---

## 💡 Usage Examples

### Basic Analysis (Expression and Accessibility)
```bash
python neural_module/neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ CAGE ATAC
```

### Complete Epigenetic Analysis
```bash
python neural_module/neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs ATAC CHIP_HISTONE CHIP_TF
```

### Splicing Analysis
```bash
python neural_module/neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ SPLICE_JUNCTIONS SPLICE_SITES SPLICE_SITE_USAGE
```

### 3D Analysis
```bash
python neural_module/neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs CONTACT_MAPS CHIP_TF
```

### Complete Analysis
```bash
python neural_module/neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF DNASE
```

---

## 🔍 How to Check Available Outputs

Run the verification script:
```bash
python scripts/check_alphagenome_outputs.py
```

This script will list all outputs that your AlphaGenome installation supports.

---

## 📋 Comparison with Old Nomenclature

If you've seen documentation with different names, here's the correspondence:

| Old Name | Correct Name | Type |
|----------|--------------|------|
| H3K27AC | CHIP_HISTONE | Histone marker |
| H3K4ME3 | CHIP_HISTONE | Histone marker |
| H3K27ME3 | CHIP_HISTONE | Histone marker |
| CTCF | CHIP_TF | Transcription factor |

**Note**: `CHIP_HISTONE` and `CHIP_TF` return predictions for **multiple** markers/factors, not just one.

---

## 🎯 Recommendations by Use Case

### 1. Regulatory Variant Analysis
```bash
--outputs RNA_SEQ CAGE ATAC CHIP_HISTONE
```
Identifies impact on expression, promoters, and enhancers.

### 2. Gene Expression Studies
```bash
--outputs RNA_SEQ CAGE PROCAP
```
Characterizes expression and transcription start sites.

### 3. Regulatory Element Analysis
```bash
--outputs ATAC DNASE CHIP_HISTONE CHIP_TF
```
Identifies enhancers, promoters, and binding sites.

### 4. Splicing Studies
```bash
--outputs RNA_SEQ SPLICE_JUNCTIONS SPLICE_SITES SPLICE_SITE_USAGE
```
Analyzes splicing patterns and isoforms.

### 5. 3D Genomics
```bash
--outputs CONTACT_MAPS CHIP_TF
```
Studies nuclear organization and loops.

---

## ⚠️ Important Notes

1. **Processing Time**: More outputs will take longer to analyze
2. **File Sizes**: Some outputs generate large files (especially CONTACT_MAPS)
3. **Resolution**: Most offer predictions at single-base resolution
4. **Combinations**: You can combine as many outputs as you want

---

## 🆘 Troubleshooting

### Error: "Output 'X' not available"
**Cause**: Output doesn't exist or incorrect name  
**Solution**: Run `python scripts/check_alphagenome_outputs.py` to see correct list

### Error: "No valid output specified"
**Cause**: All requested outputs are invalid  
**Solution**: Use outputs from the list above

### Process too slow
**Cause**: Too many outputs or sequence too long  
**Solution**: Reduce number of outputs or split sequence

---

**For more information**: https://www.alphagenomedocs.com/

*Last updated: October 2025*

