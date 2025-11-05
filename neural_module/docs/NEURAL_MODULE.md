# 🧬 Neural Module - DNA Analysis with AlphaGenome AI

## 📖 Overview

**Neural Module** is a complete implementation for DNA analysis using Google DeepMind's [AlphaGenome](https://github.com/google-deepmind/alphagenome) API. This module allows you to predict DNA functional characteristics through artificial intelligence, including:

- 🧬 **Gene Expression** (RNA-seq, CAGE, PRO-cap)
- 🔬 **Chromatin Accessibility** (ATAC-seq, DNase-seq)
- ⚛️ **Epigenetic Markers** (H3K27AC, H3K4ME3, H3K27ME3, etc.)
- 🔗 **Transcription Factors** (CTCF and others)
- 🧩 **3D Structure** (Contact Maps)
- ✂️ **Splicing** (Junction sites, site usage)

## 🎯 Key Features

✅ **11 analysis types** supported by AlphaGenome  
✅ **Advanced visualizations** (heatmaps, dashboards, comparisons)  
✅ **Variant analysis** with functional effect prediction  
✅ **Complete and intuitive** command-line interface  
✅ **Programmatic use** as a Python library  
✅ **Integration** with existing genomic pipelines  
✅ **Complete documentation** in Portuguese and English  

---

## 📚 Documentation

### 🚀 Quick Start
- **[Installation Guide](INSTALL.md)** - How to install and configure
- **[Download Sequences](../../DOWNLOAD_SEQUENCES.md)** - How to download real genomic sequences
- **[Usage Guide](USAGE.md)** - How to run analyses
- **[Interpreting Results](RESULTS.md)** - How to interpret visualizations

### 📖 Detailed Documentation
- **[Complete README](../README.md)** - Complete technical documentation
- **[Quick Start](../QUICKSTART.md)** - Get started in 5 minutes
- **[Available Outputs](../../OUTPUTS_DISPONIVEIS.md)** - List of all analysis types
- **[Supported Sizes](../../TAMANHOS_SUPORTADOS.md)** - Sequence size restrictions

### 🔧 Advanced Features
- **[Advanced Visualizations](../../VISUALIZACOES_AVANCADAS.md)** - Visualization guide
- **[Programmatic Usage](../examples/neural_example.py)** - Python code examples
- **[Pipeline Integration](INTEGRATION.md)** - Bridge with genomes_analyzer

### 🐛 Troubleshooting
- **[Applied Fixes](../../CORRECOES_APLICADAS.md)** - Resolved issues
- **[FAQ](../README.md#troubleshooting)** - Frequently asked questions

---

## 💡 Usage Example

### Basic Analysis
```bash
cd neural_module
python neural_module.py \
    -i ../example_sequence.fasta \
    -k YOUR_API_KEY \
    -o results/
```

### Variant Analysis (Sickle Cell Anemia)
```bash
python neural_module.py \
    -i ../example_sickle_cell.fasta \
    -k YOUR_API_KEY \
    -o sickle_cell_analysis/ \
    --variant 1024 A T
```

### Complete Analysis with Advanced Visualizations
```bash
python neural_module.py \
    -i gene_region.fasta \
    -k YOUR_API_KEY \
    -o comprehensive_analysis/ \
    --outputs RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF \
    --dpi 600 \
    --formats png pdf
```

---

## 🧪 Included Example: Sickle Cell Anemia

The `example_sequence.fasta` file contains the **HBB** gene region (Beta-globin) which, when mutated, causes sickle cell anemia. This is one of the most studied genetic diseases and serves as an excellent example of how a single mutation can affect gene function.

**Mutation**: Position 1024: A→T (GAG→GTG, Glu→Val)

To analyze:
```bash
cd neural_module
python neural_module.py \
    -i ../example_sequence.fasta \
    -k YOUR_API_KEY \
    -o sickle_cell/ \
    --variant 1024 A T
```

---

## 📊 Generated Outputs

### Visualizations (Default Mode: Advanced)
- **Enhanced tracks** - Multiple subplots with metadata
- **Heatmaps** - Comparison of multiple tracks
- **Dashboard** - Complete statistical summary
- **Multi-output comparison** - All outputs in one plot
- **Variant analysis** - 3 panels (overlay, difference, zoom)

### Reports
- **analysis_report.json** - JSON summary
- **Ontology metadata** - Tissue/cell information

---

## 🗂️ File Structure

```
neural_module/
├── neural_module.py                    # 🌟 Main module
├── neural_integration.py               # 🔗 Pipeline integration
├── neural_visualizations_advanced.py   # 🎨 Advanced visualizations
├── README.md                           # 📖 Main documentation
├── QUICKSTART.md                       # 🚀 Quick start guide
├── CHANGELOG.md                        # 📝 Change history
├── STRUCTURE.md                        # 📁 Structure guide
├── configs/
│   └── default.yaml                    # ⚙️ Configuration
├── docs/
│   ├── NEURAL_MODULE.md                # 📖 This file
│   ├── INSTALL.md                      # 🚀 Installation guide
│   ├── USAGE.md                        # 💡 Usage guide
│   ├── RESULTS.md                      # 📊 Interpreting results
│   └── ...
├── examples/
│   └── neural_example.py               # 📝 Usage examples
└── scripts/
    ├── demo.sh                         # 🎬 Demonstration
    ├── test.sh                         # 🧪 Tests
    └── ...
```

---

## 🔗 Useful Links

- **AlphaGenome Documentation**: https://www.alphagenomedocs.com/
- **AlphaGenome GitHub**: https://github.com/google-deepmind/alphagenome
- **Get API Key**: https://www.alphagenomedocs.com/
- **Paper**: Avsec et al. 2025 - "AlphaGenome: advancing regulatory variant effect prediction"

---

## 🤝 Contributing

Neural Module is part of the genomic analysis project. To contribute:
1. Report bugs via issues
2. Suggest improvements
3. Share use cases

---

## 📄 License

Compatible with Apache 2.0 (AlphaGenome license)  
Free use for non-commercial research

---

## 📞 Support

- **AlphaGenome**: alphagenome@google.com
- **Documentation**: See the guides listed above
- **Issues**: Open an issue in the repository

---

## ✨ Project Status

🟢 **Stable and Ready to Use**

- ✅ Fully tested
- ✅ Complete documentation
- ✅ Working examples
- ✅ Professional visualizations
- ✅ Pipeline integration

---

**Developed with ❤️ for advanced genomic analysis**

*Last update: November 2025*
