# Genes Difference Count — Pairwise Genetic Comparison Tool

A high-performance C++ tool for comparing gene sequences between family members (trio analysis: father, mother, child) using parallel processing and optimized algorithms.

## 📋 Overview

This tool performs comprehensive pairwise comparison of gene sequences from three FASTA files, calculating genetic differences with support for IUPAC nucleotide codes. It generates detailed statistics on matches, substitutions, and indels for each gene across all three family member pairs.

## ✨ Key Features

- **🚀 Parallel Processing**: OpenMP-based parallelization for maximum performance
- **🧬 IUPAC Compatibility**: Full support for IUPAC nucleotide ambiguity codes
- **⚡ Smart Alignment**: Hirschberg algorithm for short sequences (<200K bases), fast estimation for longer ones
- **📊 Comprehensive Statistics**: Detailed per-gene metrics including matches, substitutions, indels, and ambiguous bases
- **💾 CSV Output**: Structured results in CSV format for easy analysis
- **🔄 Idempotent**: Reuses existing CSV output to avoid recomputation
- **🎯 Optimized**: Aggressive compiler optimizations for ultra-fast execution

## 🔧 Configuration

Before compiling, you need to configure the input and output file paths in the source code. Edit `genes_difference_count.cpp` and modify the constants at the top of the file (lines 16-24):

```cpp
// Configuration constants
const std::string IN_PAI = "fasta/NA12891.genes.consensus.fa";      // Father's FASTA
const std::string IN_MAE = "fasta/NA12892.genes.consensus.fa";      // Mother's FASTA
const std::string IN_FILHA = "fasta/NA12878.genes.consensus.fa";    // Child's FASTA
const std::string OUT_CSV = "fasta/family_pairwise_differences.csv"; // Output CSV
const int MAX_ALIGN_LEN = 200000;    // Maximum length for full alignment
const int GAP_PENALTY = -1;          // Gap penalty for alignment
const int MATCH_SCORE = 1;           // Match score for alignment
const int MISMATCH_SCORE = -1;       // Mismatch penalty for alignment
```

### Input FASTA Format

The tool expects FASTA files with headers in this format:
```
>SampleID|ENSG_ID::chr:start-end(strand)
SEQUENCE_DATA
```

Example:
```
>NA12878|ENSG00000223972::chr1:11869-14409(+)
ATCGATCGATCG...
```

## 🔨 Compilation

### Prerequisites

- **C++ Compiler**: g++ with C++17 support
- **OpenMP**: For parallel processing (usually included with g++)
- **Build Tools**: make

Install on Ubuntu/Debian:
```bash
sudo apt install build-essential
```

### Build Commands

The Makefile provides several targets:

#### 1. Build Optimized Version (Recommended)
```bash
cd genes_difference_count
make
```

This creates the `genes_difference_count` executable with maximum optimizations:
- `-O3`: Highest optimization level
- `-march=native`: CPU-specific optimizations
- `-mtune=native`: Tuning for current CPU
- `-flto`: Link-time optimization
- `-funroll-loops`: Loop unrolling
- `-ffast-math`: Fast floating-point math
- `-fopenmp`: OpenMP parallelization

#### 2. Build Debug Version
```bash
make debug
```

Creates `genes_difference_count_debug` with debugging symbols and no optimizations.

#### 3. Clean Build Artifacts
```bash
make clean
```

#### 4. Check Dependencies
```bash
make deps
```

#### 5. Show Compilation Info
```bash
make info
```

#### 6. Build and Test
```bash
make test
```

## 🚀 Usage

### Basic Usage

After configuration and compilation:

```bash
cd genes_difference_count
./genes_difference_count
```

The program will:
1. Load all three FASTA files
2. Find genes common to all three samples
3. Perform pairwise comparisons:
   - Father vs Mother
   - Father vs Child
   - Mother vs Child
4. Output detailed statistics
5. Save results to CSV

### Example Output

```
🧬 Iniciando análise de diferenças genéticas (versão C++ PARALELA)...
🚀 Usando 16 threads para processamento paralelo
📁 Carregando arquivos FASTA...
   ✓ Pai: 19,969 genes carregados
   ✓ Mãe: 19,969 genes carregados
   ✓ Filha: 19,969 genes carregados

📊 Encontrados 19,969 genes em comum para análise
🔬 Cada gene será comparado em 3 pares: pai×mãe, pai×filha, mãe×filha
⏱️  Total de comparações: 59,907

🚀 Iniciando processamento paralelo...

🧬 Gene 1/19969: ENSG00000223972
   🔬 pai × mãe (2540 vs 2540 bases) - alinhamento completo ➤ 127 diferenças (5.00000%) em 0.012s
   🔬 pai × filha (2540 vs 2540 bases) - alinhamento completo ➤ 64 diferenças (2.51968%) em 0.010s
   🔬 mãe × filha (2540 vs 2540 bases) - alinhamento completo ➤ 63 diferenças (2.48031%) em 0.011s
   📊 Progresso: 3/59907 (0.0%) - ETA: 12m 34s

...

✅ ANÁLISE CONCLUÍDA!
⏱️  Tempo total: 8m 23s
📊 19,969 genes analisados em 59,907 comparações
⚡ Velocidade média: 119.2 comparações/segundo
```

### Idempotent Behavior

If the output CSV file already exists, the program will:
- Read and parse the existing CSV
- Display summary statistics
- Skip recomputation

To force recomputation, delete the CSV file:
```bash
rm fasta/family_pairwise_differences.csv
```

## 📊 Output Format

The tool generates a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier (ENSG ID) |
| `pair` | Comparison pair (pai_vs_mãe, pai_vs_filha, mãe_vs_filha) |
| `aligned_cols` | Total aligned columns (including gaps) |
| `informative_cols` | Informative columns (excluding N's) |
| `matches` | Number of matching/compatible bases |
| `subs` | Number of substitutions (mismatches) |
| `indels` | Number of insertions/deletions |
| `ambiguous_N` | Number of ambiguous bases (N) ignored |
| `differences_total` | Total differences (subs + indels) |

### Example CSV Content

```csv
gene,pair,aligned_cols,informative_cols,matches,subs,indels,ambiguous_N,differences_total
ENSG00000223972,pai_vs_mãe,2540,2540,2413,127,0,0,127
ENSG00000223972,pai_vs_filha,2540,2540,2476,64,0,0,64
ENSG00000223972,mãe_vs_filha,2540,2540,2477,63,0,0,63
ENSG00000227232,pai_vs_mãe,1351,1351,1298,53,0,0,53
...
```

## 📈 Performance

### Optimization Techniques

1. **Parallel Processing**: OpenMP parallelization across genes
2. **Cache-Based IUPAC Lookup**: O(1) compatibility checks using pre-computed cache
3. **Adaptive Algorithm**: Full alignment for short sequences, fast estimation for long ones
4. **Memory Pre-allocation**: Reduced allocations during processing
5. **Vectorized Operations**: Chunk-based processing for improved cache locality
6. **Thread-Safe Accumulation**: Lock-free counters with atomic operations

### Expected Performance

On a system with 16 cores:
- **Speed**: ~100-150 comparisons/second
- **Typical Runtime**: 8-12 minutes for ~60,000 comparisons (20,000 genes × 3 pairs)
- **Memory Usage**: Moderate (~1-2 GB depending on FASTA sizes)

### Performance Tips

1. **More Cores = Faster**: The tool scales well with available CPU cores
2. **SSD Storage**: Faster I/O for loading FASTA files
3. **Sufficient RAM**: Ensure adequate memory for loading all sequences
4. **Native Compilation**: The `-march=native` flag optimizes for your specific CPU

## 🧬 Algorithm Details

### IUPAC Nucleotide Compatibility

The tool supports all IUPAC nucleotide codes:

| Code | Bases | Meaning |
|------|-------|---------|
| A | A | Adenine |
| C | C | Cytosine |
| G | G | Guanine |
| T | T | Thymine |
| R | A,G | Purine |
| Y | C,T | Pyrimidine |
| S | G,C | Strong |
| W | A,T | Weak |
| K | G,T | Keto |
| M | A,C | Amino |
| B | C,G,T | Not A |
| D | A,G,T | Not C |
| H | A,C,T | Not G |
| V | A,C,G | Not T |
| N | A,C,G,T | Any base |

Two bases are considered compatible if their IUPAC code sets have any intersection.

### Alignment Strategy

1. **Short Sequences** (≤200,000 bases):
   - Uses Hirschberg algorithm for space-efficient alignment
   - Dynamic programming with O(n×m) time, O(min(n,m)) space

2. **Long Sequences** (>200,000 bases):
   - Fast estimation without full alignment
   - Position-by-position comparison
   - Handles length differences as indels

## 🐛 Troubleshooting

### Compilation Errors

**Error**: `g++: command not found`
```bash
sudo apt install build-essential
```

**Error**: OpenMP not found
```bash
sudo apt install libomp-dev
```

### Runtime Errors

**Error**: Cannot open input file
- Check that file paths in `genes_difference_count.cpp` are correct
- Ensure FASTA files exist at specified locations

**Error**: Out of memory
- Reduce the number of genes
- Use a system with more RAM
- Process files in batches

### Performance Issues

**Slow execution**:
- Ensure you compiled with optimizations (`make`, not `make debug`)
- Check available CPU cores: `nproc`
- Monitor CPU usage: `htop`

## 📚 Technical Notes

### File Structure

```
genes_difference_count/
├── genes_difference_count.cpp   # Main source code
├── Makefile                      # Build configuration
└── README.md                     # This file
```

### Dependencies

- C++17 standard library
- OpenMP (included with g++)
- No external libraries required

### Thread Safety

- All statistics accumulation is thread-safe
- Uses mutexes for shared data structures
- Atomic counters for progress tracking

## 📖 Related Tools

This tool is part of the Genomes Analyzer pipeline. For more information:
- **Main Pipeline**: [../README.md](../README.md)
- **Other Modules**: See Appendix 2 in main README

## 🤝 Contributing

When modifying the code:
1. Maintain the existing code style
2. Test with debug build first (`make debug`)
3. Ensure thread safety for parallel sections
4. Update this README if adding features

## 📄 License

Part of the Genomes Analyzer project. See [../LICENSE](../LICENSE) for details.

---

**Last Updated**: November 2025

