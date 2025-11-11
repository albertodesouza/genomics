# 📁 Structure of build_non_longevous_dataset Module

This directory contains the **Non-Longevous Dataset Builder** module completely organized.

## 📂 File Structure

```
build_non_longevous_dataset/
├── build_non_longevous_dataset.py    # Main program
├── build_window_and_predict.py       # Script to extract genomic windows and predictions
├── frog_ancestry_parser.py           # FROGAncestryCalc likelihood parser
├── dataset_builder.py                # Individual and global metadata builders
├── genomic_dataset.py                # PyTorch Dataset implementation
├── __init__.py                       # Module initialization
├── README.md                         # Complete documentation
├── configs/
│   ├── default.yaml                  # Default configuration
│   └── small.yaml                    # Small dataset configuration
├── docs/
│   ├── AISNP_MODE.md                 # AISNP mode documentation
│   ├── ALPHAGENOME_PREDICTIONS.md    # AlphaGenome predictions guide
│   ├── ALPHAGENOME_TISSUES.md        # AlphaGenome tissues/cells guide
│   ├── HAPLOTYPES.md                 # Haplotype generation guide
│   ├── IMPLEMENTATION.md             # Technical implementation details
│   ├── PYTORCH_DATASET.md            # PyTorch Dataset complete documentation
│   ├── PYTORCH_DATASET_SUMMARY.md    # Implementation summary
│   ├── QUICKSTART.md                 # Quick start guide
│   └── STRUCTURE.md                  # This file
├── examples/
│   └── load_dataset_example.py       # PyTorch Dataset usage examples
└── scripts/
    └── test.sh                       # Test script
```

## 🚀 How to Use

### From the module directory:
```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

### From the project root:
```bash
python3 build_non_longevous_dataset/build_non_longevous_dataset.py \
  --config build_non_longevous_dataset/configs/default.yaml
```

### Using the test script:
```bash
cd build_non_longevous_dataset/scripts
bash test.sh
```

## 📝 Relative Paths

Paths in the `configs/default.yaml` file are relative to the `configs/` directory:
- `../../doc/file.csv` → `/path/to/genomics/doc/file.csv`
- `../../refs/genome.fa` → `/path/to/genomics/refs/genome.fa`

## 📚 Documentation

- **README.md** (root): Complete module documentation
- **docs/QUICKSTART.md**: Quick guide to get started
- **docs/IMPLEMENTATION.md**: Technical implementation details
- **docs/PYTORCH_DATASET.md**: Complete PyTorch Dataset documentation
- **docs/PYTORCH_DATASET_SUMMARY.md**: Implementation summary
- **docs/STRUCTURE.md**: This file - module structure

## ✅ Tested and Working

✓ Execution from module directory
✓ Execution from project root
✓ Functional test script
✓ Correct relative path resolution
✓ Integration with build_window_and_predict.py (included in the module)

## 🔗 build_window_and_predict.py

This module includes `build_window_and_predict.py`, which is responsible for:
- Extracting 1 Mb genomic windows around specific genes
- Applying variants from 1000 Genomes samples
- Generating consensus sequences per haplotype
- Running AlphaGenome predictions (optional)

📚 Additional documentation:
- [AlphaGenome Predictions Guide](ALPHAGENOME_PREDICTIONS.md)
- [Tissues/Cells Guide](ALPHAGENOME_TISSUES.md)

