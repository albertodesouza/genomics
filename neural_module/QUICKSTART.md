# 🚀 Neural Module - Quick Start Guide

## ⚡ Get Started in 5 Minutes

### 1. Check Requirements

```bash
bash scripts/check_requirements.sh
```

### 2. Install AlphaGenome

```bash
bash ../install_alphagenome.sh
```

### 3. Obtain API Key

Visit: [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/)

### 4. Run Analysis

```bash
cd neural_module
python neural_module.py \
    -i ../example_sequence.fasta \
    -k YOUR_API_KEY \
    -o results/
```

---

## 📋 Useful Commands

### Basic Analysis

```bash
# Complete analysis with default outputs
python neural_module.py -i sequence.fasta -k API_KEY -o results/
```

### Custom Analysis

```bash
# Choose specific outputs
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --outputs RNA_SEQ ATAC H3K27AC

# High resolution
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --dpi 600 \
    --formats png pdf svg
```

### Variant Analysis

```bash
# Analyze effect of A→C variant at position 1000
python neural_module.py \
    -i sequence.fasta \
    -k API_KEY \
    -o results/ \
    --variant 1000 A C
```

---

## 🔗 Integration with genomes_analyzer.py

### Complete Workflow

```bash
# 1. Run genomic pipeline
python genomes_analyzer.py config.yaml

# 2. Extract sequences of interest
cd neural_module
python neural_integration.py \
    --extract-vcf \
    --vcf ../results/variants.vcf \
    --ref ../reference.fa \
    --output extracted.fasta

# 3. Analyze with AlphaGenome
python neural_module.py \
    -i extracted.fasta \
    -k API_KEY \
    -o neural_results/

# 4. OR do everything at once:
python neural_integration.py \
    --integrated \
    --vcf ../results/variants.vcf \
    --ref ../reference.fa \
    --api-key API_KEY \
    --output integrated_results/
```

### Extract Specific Genes

```bash
cd neural_module
python neural_integration.py \
    --extract-genes \
    --genes BRCA1 TP53 EGFR \
    --gtf ../annotations.gtf \
    --ref ../reference.fa \
    --output genes.fasta
```

---

## 📊 Output Types

| Output | Description | Typical Use |
|--------|-------------|-------------|
| `RNA_SEQ` | Gene expression | Identify active genes |
| `CAGE` | Cap analysis | Identify promoters |
| `ATAC` | Accessibility | Open regulatory regions |
| `H3K27AC` | Histone marker | Active enhancers |
| `H3K4ME3` | Histone marker | Active promoters |
| `H3K27ME3` | Histone marker | Polycomb repression |
| `CTCF` | Transcription factor | Insulators/loops |

---

## 💡 Usage Examples

### Example 1: Pathogenic Variant Research

```bash
# Analyze variant in disease gene
python neural_module.py \
    -i disease_gene.fasta \
    -k API_KEY \
    -o pathogenic_analysis/ \
    --variant 5234 G A \
    --outputs RNA_SEQ CAGE ATAC H3K27AC \
    --formats png pdf
```

### Example 2: Genomic Region Characterization

```bash
# Complete region analysis
python neural_module.py \
    -i genomic_region.fasta \
    -k API_KEY \
    -o region_analysis/ \
    --outputs RNA_SEQ CAGE ATAC H3K27AC H3K4ME3 H3K27ME3 CTCF \
    --dpi 600
```

### Example 3: Multiple Variant Comparison

```bash
# Create FASTAs with different variants
# Then analyze each one
for variant in variant1.fasta variant2.fasta variant3.fasta; do
    python neural_module.py \
        -i $variant \
        -k API_KEY \
        -o comparison/$(basename $variant .fasta)/
done
```

---

## 🧪 Tests

### Quick Test

```bash
# Test with example sequence
python neural_module.py \
    -i ../example_sequence.fasta \
    -k API_KEY \
    -o test_results/
```

### Complete Test Suite

```bash
cd scripts
bash test.sh YOUR_API_KEY
```

### Programmatic Examples

```bash
# Run all examples
cd examples
python neural_example.py -k API_KEY

# Run specific example
python neural_example.py -k API_KEY -e 3
```

---

## 📁 File Structure

```
neural_module/
├── neural_module.py              # Main module
├── neural_integration.py         # Integration with genomes_analyzer
├── neural_visualizations_advanced.py
├── README.md                     # Main documentation
├── QUICKSTART.md                 # This guide
├── CHANGELOG.md
├── ESTRUTURA.md
├── configs/
│   └── default.yaml              # Configuration
├── docs/
│   ├── NEURAL_MODULE.md          # Complete documentation
│   ├── INTEGRATION.md
│   ├── USAGE.md
│   ├── RESULTS.md
│   └── ...
├── examples/
│   └── neural_example.py         # Usage examples
└── scripts/
    ├── demo.sh                   # Demonstration
    ├── test.sh                   # Test suite
    ├── check_requirements.sh     # Requirements check
    └── show_summary.sh
```

---

## ❓ Troubleshooting

### Problem: "AlphaGenome is not installed"

```bash
bash ../install_alphagenome.sh
```

### Problem: "Invalid API key"

- Check if you copied the key correctly
- Make sure it's active at alphagenomedocs.com

### Problem: "Sequence too long"

- AlphaGenome supports up to 1 Mbp
- Split longer sequences into chunks

### Problem: Plots are not generated

```bash
pip install --upgrade matplotlib seaborn
```

### Problem: ImportError

```bash
# Reinstall dependencies
pip install --force-reinstall rich matplotlib numpy
```

---

## 🔍 Verify Installation

```bash
# From neural_module directory
cd neural_module

# Check neural_module
python -c "from neural_module import AlphaGenomeAnalyzer; print('✓ OK')"

# Check AlphaGenome
python -c "from alphagenome.models import dna_client; print('✓ OK')"

# Check dependencies
python -c "import rich, matplotlib, numpy; print('✓ OK')"
```

---

## 📚 Additional Resources

- **Complete Documentation**: `README.md`
- **Demonstration**: `bash scripts/demo.sh`
- **AlphaGenome API**: [https://www.alphagenomedocs.com/](https://www.alphagenomedocs.com/)
- **AlphaGenome GitHub**: [https://github.com/google-deepmind/alphagenome](https://github.com/google-deepmind/alphagenome)

---

## 💬 Support

### Questions about AlphaGenome
- Email: alphagenome@google.com
- Documentation: https://www.alphagenomedocs.com/

### Questions about neural_module
- Open an issue in the repository
- Consult the complete documentation

---

## ⚠️ Important Notes

1. **Non-Commercial Use**: Free API only for research
2. **Limited Rate**: Suitable for ~1000s of predictions
3. **Internet Required**: API requires online connection
4. **Security**: Never share your API key

---

## 🎯 Next Steps

After running the basic test:

1. ✅ Explore different available outputs
2. ✅ Test variant analysis
3. ✅ Integrate with your existing pipeline
4. ✅ Experiment with your own sequences
5. ✅ Use programmatic examples for custom workflows

**Happy analyzing! 🧬🚀**
