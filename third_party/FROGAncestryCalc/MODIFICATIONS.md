# Modifications to FROGAncestryCalc

## Summary

The source code was modified to accept input files delimited by **pipe (|)** instead of comma (,).

## Modified Files

### 1. `src/dv/ValidateFileHeader.java`
- Changed `split("\\,")` → `split("\\|")`
- Changed `indexOf(',')` → `indexOf('|')`
- Changed `countOccurrences(record, ',')` → `countOccurrences(record, '|')`
- Updated error messages: "comma delimited" → "pipe delimited"

### 2. `src/read/ReadTxtFiles.java`
- Changed `split("\\,")` → `split("\\|")` in genotype reading

## Input File Format

Input files must have the following format:

```
Individual|rs10497191|rs1079597|rs11652805|...(other SNPs)
HG02561_GWD|NN|CC|CC|...(genotypes)
HG02562_GWD|TT|CT|CC|...(genotypes)
```

**Important:**
- Delimiter: pipe `|` (not comma)
- Line endings: Unix (LF)
- Encoding: UTF-8

## How to Run

### From compiled sources:

```bash
./run_from_source.sh
```

Or manually:

```bash
LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 java -cp bin main.ComputeBatchAnalysis
```

### Note on Locale

It's essential to use `LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8` to avoid problems with number formatting in scientific notation (e.g., `5.652E-62` instead of `5,652E-62`).

## Recompilation

If you need to recompile after modifications:

```bash
cd /home/lume2/genomics/FROGAncestryCalc
rm -rf bin
mkdir bin
javac -d bin -sourcepath src $(find src -name "*.java")
cp -r src/read/data bin/read/
```

## Output Files

The application generates 3 files in `output/`:
- `*_likelihood.txt` - Ancestry probabilities for each population
- `*_orderOfMag.txt` - Order of magnitude of probabilities
- `*_rankOrder.txt` - Ranking of populations by probability

## Backup

Original files were saved with `.bak` extension:
- `src/dv/ValidateFileHeader.java.bak`

