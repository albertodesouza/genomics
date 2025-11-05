# 📁 Neural Longevity Dataset Module Structure

This directory contains the **Neural Longevity Dataset Builder** - a tool for creating machine learning datasets from longevity-associated genomic markers.

## 📂 File Structure

\`\`\`
neural_longevity_dataset/
├── neural_longevity_dataset.py  # Main program - dataset builder
├── longevity_train.py           # Model training script
├── README.md                    # Main documentation
├── QUICKSTART.md                # Quick start guide
├── STRUCTURE.md                 # This file
├── configs/
│   ├── default.yaml             # Default configuration
│   ├── top3.yaml                # Top 3 samples config
│   └── train.yaml               # Training configuration
├── docs/
│   └── PROJECT.md               # Project documentation
├── examples/
└── scripts/
    └── test_dryrun.sh           # Dry-run test script
\`\`\`

## 🚀 How to Use

### From the module directory:
\`\`\`bash
cd neural_longevity_dataset
python neural_longevity_dataset.py --config configs/default.yaml
\`\`\`

### From the project root:
\`\`\`bash
python neural_longevity_dataset/neural_longevity_dataset.py \\
  --config neural_longevity_dataset/configs/default.yaml
\`\`\`

### Using scripts:
\`\`\`bash
cd neural_longevity_dataset/scripts
bash test_dryrun.sh
\`\`\`

## 📝 Main Programs

### neural_longevity_dataset.py
Main dataset builder for longevity research:
- Downloads 1000 Genomes samples
- Calls variants with bcftools
- Selects central points (variants)
- Extracts FASTA windows centered on variants
- Processes with AlphaGenome for features
- Consolidates balanced PyTorch dataset (train/val/test)

### longevity_train.py
Model training script:
- Loads PyTorch dataset from neural_longevity_dataset.py
- Trains neural network on longevity markers
- Validates and evaluates performance
- Saves trained models

## 📚 Documentation

- **README.md**: Overview and quick start
- **QUICKSTART.md**: 5-minute guide
- **docs/PROJECT.md**: Complete project documentation

## 🔗 Integration

This module can be used:
1. **Standalone**: Build datasets from 1000 Genomes data
2. **With neural_module**: Use AlphaGenome for feature extraction
3. **With training**: Use longevity_train.py for model training
4. **Pipeline**: In automated workflows

## ✅ Tested and Working

✓ 1000 Genomes High Coverage VCF download
✓ Variant calling with bcftools
✓ Central point selection
✓ FASTA extraction with ALT allele application
✓ AlphaGenome feature extraction
✓ PyTorch dataset generation
✓ Self-contained module

---

**Location**: \`neural_longevity_dataset/\` (organized on 2025-11-05)
