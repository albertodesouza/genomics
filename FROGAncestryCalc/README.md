# FROGAncestryCalc - Modified Version

FROG-kb (Forensic Resource/Reference On Genetics - Knowledge base) Ancestry Inference Batch Likelihood Computation Tool - Modified to use pipe delimiters.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🚀 Quick Start

### Run Analysis
```bash
./run.sh
```

### Recompile Code
```bash
./recompile.sh
```

## 📋 Configuration

Edit the `FROGAncestryCalc.properties` file:

```properties
homePath=/home/lume2/genomics/frog/FROGAncestryCalc
inputFilename=55_aisnp_1000_Genome.txt
panelInfo=55AI
```

**⚠️ IMPORTANT:** Update the properties file with the appropriate input file name and AI panel name before starting a new job.

### Available AI Panels

| Panel Code | Description | SNP Count |
|-----------|-------------|-----------|
| `55AI` | KiddLab - Set of 55 AISNPs | 55 |
| `128AI` | Seldin's list of 128 AISNPs | 128 |
| `34plex` | SNPforID 34-plex | 34 |
| `combined` | Combined panel (Kiddlab-55 + Seldin's-128 + SNPforID34-plex) | 192 |
| `precision` | Precision ID Ancestry Panel | 165 |

## 📂 Input File Format

Place your input files in the `input/` directory with the following format:

```
Individual|rs10497191|rs1079597|rs11652805|...|rs9522149
HG02561_GWD|NN|CC|CC|CC|...|TT
HG02562_GWD|TT|CT|CC|CC|...|TT
```

### Format Specifications

- ✅ **Delimiter:** pipe `|`
- ✅ **Line endings:** Unix (LF)
- ✅ **Encoding:** UTF-8
- ✅ **First line:** Header with "Individual" + ordered list of SNP IDs
- ✅ **Following lines:** Individual ID + genotypes
- ✅ **SNP order:** Must match the order in the corresponding sample file
- ✅ **Individual IDs:** Must be unique

### Preparing Your Input File

1. Follow the SNP order given in the sample files for your chosen AI panel (see `sampleInputFiles/`)
2. SNP labels and genotypes must be ordered by rs number (alphanumeric)
3. Use the sorting function in Excel or similar tools (ascending order)
4. Ensure all Individual Identifiers are unique
5. Consult the appropriate file in `SNPInfo/` to find valid alleles for each SNP
6. Use accepted genotype notations:
   - Two-allele format: `AA`, `TT`, `GG`, `CC`, `AT`, `AG`, etc.
   - Missing data: `NN`

## 📊 Output Files

Generated in the `output/` directory:

| File | Description |
|------|-------------|
| `*_likelihood.txt` | Likelihood values for ancestral population for each individual across 155 populations |
| `*_orderOfMag.txt` | Order of magnitude of the likelihoods |
| `*_rankOrder.txt` | Population rankings by likelihood for each individual |

All output files are tab-delimited and can be opened in Excel.

**Note:** Output files from previous jobs (including any `errFile.txt`) are deleted at the start of a new job.

## 🗂️ Project Structure

```
FROGAncestryCalc/
├── src/                        # Modified source code
│   ├── bean/                   # Data classes
│   ├── dv/                     # Validation (modified for pipes)
│   ├── main/                   # Main application class
│   ├── read/                   # File reading (modified)
│   ├── sub/                    # Helper classes
│   └── write/                  # Output writing
├── bin/                        # Compiled classes
├── input/                      # Input files directory
│   ├── ind/                    # Working directory (do not delete)
│   └── indGenotype/            # Working directory (do not delete)
├── output/                     # Results directory
├── SNPInfo/                    # SNP information for each panel
│   ├── 55_aisnps_alleles.txt
│   ├── 128_aisnps_alleles.txt
│   ├── 34_plex_alleles.txt
│   ├── combined_alleles.txt
│   └── precision_alleles.txt
├── sampleInputFiles/           # Sample input files
│   ├── 55_aisnps_sample.txt
│   ├── 128_aisnps_sample.txt
│   ├── 34_plex_sample.txt
│   ├── combined_sample.txt
│   └── precision_sample.txt
├── log/                        # Execution logs
│   └── workingLog.txt
├── obsolete/                   # Original files with bugs
├── run.sh          # Execution script
├── recompile.sh                # Recompilation script
├── FROGAncestryCalc.properties # Configuration file
└── MODIFICACOES.md             # Technical modification details
```

## ⚙️ Requirements

- **Java:** 17+ (OpenJDK recommended)
- **Shell:** Bash
- **OS:** Linux/Unix

## 🔄 Modifications from Original

This modified version includes the following improvements:

### 1. **Pipe Delimiter Support**
- Changed from comma (`,`) to pipe (`|`) delimiter
- Modified validation and parsing logic
- Updated error messages

### 2. **Locale Fix**
- Added `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8` to prevent number formatting issues
- Resolves `NumberFormatException` with scientific notation (e.g., `5.652E-62`)

### 3. **Linux Compatibility**
- Unix line endings (LF)
- Proper path handling
- Shell script optimizations

### Modified Files:
- `src/dv/ValidateFileHeader.java` - Validation logic
- `src/read/ReadTxtFiles.java` - File parsing

For complete technical details, see [`MODIFICACOES.md`](MODIFICACOES.md)

## 📝 Error Handling

### Error File
If validation errors occur, check `output/errFile.txt` for details.

### Working Log
View `log/workingLog.txt` for:
- Processing information for all jobs
- Copy of error messages
- Timestamps and status updates

**Note:** The log file accumulates across jobs until manually deleted.

## 🛠️ Maintenance

### Recompiling After Code Changes

```bash
./recompile.sh
```

Or manually:

```bash
cd /path/to/FROGAncestryCalc
rm -rf bin
mkdir bin
javac -d bin -sourcepath src $(find src -name "*.java")
cp -r src/read/data bin/read/
```

### Cleaning Up

```bash
# Clean output files
rm -f output/*.txt

# Clean working directories
rm -f input/ind/* input/indGenotype/*

# Clean logs (optional)
rm -f log/workingLog.txt
```

## 📚 Population Coverage

The tool calculates ancestry likelihoods for **155 populations** including:

- African populations (Yoruba, Mbuti, Biaka, etc.)
- European populations (Danes, Finns, British, etc.)
- Asian populations (Han Chinese, Japanese, Korean, etc.)
- American populations (Maya, Pima, Karitiana, etc.)
- Middle Eastern populations (Druze, Bedouin, Palestinian, etc.)
- And many more...

## 🐛 Troubleshooting

### Common Issues

1. **"Your input file is not pipe delimited"**
   - Ensure file uses pipe `|` as delimiter, not comma
   - Check for proper line endings (Unix LF, not Windows CRLF)

2. **"NumberFormatException"**
   - Make sure to run with proper locale settings
   - Use `./run.sh` which handles this automatically

3. **"Missing SNPs" or "Wrong SNP count"**
   - Verify SNP order matches the sample file
   - Check that all required SNPs are present
   - Ensure no extra columns or missing data

4. **Java version issues**
   - Requires Java 17 or higher
   - Check version: `java -version`

## 📄 License

MIT License

Copyright (c) 2019 haseenaR

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

## 🙏 Acknowledgments

- Original FROG-kb tool by haseenaR
- Reference populations from various genomics databases
- Modified for improved usability and Linux compatibility

## 📞 Support

For issues related to:
- **Original tool:** Refer to FROG-kb documentation
- **This modified version:** Check `MODIFICACOES.md` for technical details

---

**Note:** The original JAR and Windows batch files have been moved to `obsolete/` as they contained bugs incompatible with this input format.

