# 📁 Structure of build_non_longevous_dataset Module

This directory contains the **Non-Longevous Dataset Builder** module completely organized.

## 📂 File Structure

```
build_non_longevous_dataset/
├── build_non_longevous_dataset.py    # Main program
├── build_window_and_predict.py       # Script to extract genomic windows and predictions
├── README.md                         # Complete documentation
├── QUICKSTART.md                     # Quick start guide
├── IMPLEMENTATION.md                 # Technical implementation details
├── STRUCTURE.md                      # This file
├── configs/
│   └── default.yaml                  # Default configuration
├── docs/
│   ├── ALPHAGENOME_PREDICTIONS.md    # AlphaGenome predictions guide
│   └── ALPHAGENOME_TISSUES.md        # AlphaGenome tissues/cells guide
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

- **README.md**: Complete module documentation
- **QUICKSTART.md**: Quick guide to get started
- **IMPLEMENTATION.md**: Technical implementation details

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
- [AlphaGenome Predictions Guide](docs/ALPHAGENOME_PREDICTIONS.md)
- [Tissues/Cells Guide](docs/ALPHAGENOME_TISSUES.md)

