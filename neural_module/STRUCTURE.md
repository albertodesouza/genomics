# 📁 Neural Module Structure

This directory contains the **Neural Module** - DNA analysis tools using AlphaGenome (Google DeepMind).

## 📂 File Structure

```
neural_module/
├── neural_module.py              # Main program - sequence analysis
├── neural_integration.py         # Genomic pipeline integration
├── neural_visualizations_advanced.py  # Advanced visualizations
├── README.md                     # Main documentation
├── QUICKSTART.md                 # Quick start guide
├── CHANGELOG.md                  # Change history
├── STRUCTURE.md                  # This file
├── configs/
│   └── default.yaml              # Default configuration
├── docs/
│   ├── NEURAL_MODULE.md          # Complete technical documentation
│   ├── INDEX.md                  # Documentation index
│   ├── INTEGRATION.md            # Integration guide
│   ├── USAGE.md                  # Detailed usage guide
│   ├── RESULTS.md                # Results interpretation
│   ├── INSTALL.md                # Installation guide
│   ├── LEIA-ME.md                # Portuguese documentation
│   └── README.html               # HTML documentation
├── examples/
│   └── neural_example.py         # Usage examples
└── scripts/
    ├── demo.sh                   # Usage demonstration
    ├── test.sh                   # Test script
    ├── show_summary.sh           # Display results summary
    └── check_requirements.sh     # Check dependencies
```

## 🚀 How to Use

### From the module directory:
```bash
cd neural_module
python neural_module.py -i ../example_sequence.fasta -k YOUR_API_KEY -o results/
```

### From the project root:
```bash
python neural_module/neural_module.py \
  -i example_sequence.fasta \
  -k YOUR_API_KEY \
  -o results/
```

### Using scripts:
```bash
cd neural_module/scripts
bash demo.sh      # View usage examples
bash test.sh      # Run tests
```

## 📝 Main Programs

### neural_module.py
Main DNA sequence analysis with AlphaGenome:
- Gene expression predictions (RNA-seq, CAGE)
- Chromatin accessibility (ATAC-seq, DNase-seq)
- Histone markers (H3K27AC, H3K4ME3, etc.)
- Transcription factors (CTCF)
- Variant analysis (SNP effects)
- High-quality visualizations

### neural_integration.py
Integration with genomic pipelines:
- Sequence extraction from VCF
- Conversion to FASTA
- AlphaGenome analysis
- Results correlation
- 4 operation modes: integrated, vcf, bed, gene

### neural_visualizations_advanced.py
Advanced visualizations:
- Multi-output heatmaps
- Comparative plots
- Interactive dashboards
- Multi-format export

## 📚 Documentation

- **README.md**: Overview and quick start
- **QUICKSTART.md**: 5-minute guide
- **docs/NEURAL_MODULE.md**: Complete technical documentation
- **docs/INTEGRATION.md**: Integration with other pipelines
- **docs/USAGE.md**: Advanced usage guide
- **docs/RESULTS.md**: Interpreting results
- **docs/INSTALL.md**: Detailed installation

## 🔗 Integration

This module can be used:
1. **Standalone**: Direct FASTA sequence analysis
2. **Integrated**: With genomes_analyzer.py via neural_integration.py
3. **Programmatic**: Import in Python scripts
4. **Pipeline**: In automated workflows

## ✅ Tested and Working

✓ FASTA sequence analysis
✓ AlphaGenome predictions
✓ Advanced visualizations
✓ VCF integration
✓ Variant analysis
✓ Self-contained module

---

**Location**: `neural_module/` (organized on 2025-11-04)

