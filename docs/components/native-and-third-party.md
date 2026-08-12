# Native And Third-Party Tools

## Native: Genes Difference Count

Path: `native/genes_difference_count/`

This is a C++/OpenMP tool for pairwise gene sequence comparison. It is not part of the Python package, but it is kept in the repository because it supports operational genomic comparisons.

### Implementation

| File | Purpose |
|---|---|
| `genes_difference_count.cpp` | Main C++ implementation |
| `Makefile` | Compiler flags and build target |
| `scripts/generate_gene_fastas.sh` | Generates FASTA inputs from configured genomic sources |
| `scripts/run_cpp_analysis.sh` | Example runner for compiled analysis |

The tool compares gene FASTA sequences across pairs, handles IUPAC ambiguity codes, and uses OpenMP to parallelize comparisons. It is useful when the Python workflow has already produced gene FASTA inputs and a fast pairwise comparison pass is needed.

Typical flow:

```bash
cd native/genes_difference_count
make
./genes_difference_count
```

## Third-Party: FROGAncestryCalc

Path: `third_party/FROGAncestryCalc/`

This is a modified external ancestry calculator focused on AISNP panels. It includes tools to extract SNPs from VCF and other genomic data sources.

### Implementation

| Path | Purpose |
|---|---|
| `src/` | Java source for likelihood computation |
| `SNPInfo/` | AISNP panel definitions and allele metadata |
| `tools/vcf_to_frog.py` | Converts VCF genotypes into FROG input format |
| `tools/extract_snps_from_1000genomes.sh` | Downloads/extracts AISNPs from 1000 Genomes resources |
| `run.sh` | Runs the Java analysis with configured properties |

The repository version has local modifications for pipe-delimited genotypes and additional extraction workflows. Treat this component as third-party: avoid importing it from active `genomics` Python code.

Typical flow:

```bash
cd third_party/FROGAncestryCalc
python3 tools/vcf_to_frog.py sample.vcf.gz tools/aisnps_55_list.txt input/sample_data.txt
./run.sh
```

Use the SNP ancestry pipeline for integrated `genomics` CLI runs. Use FROGAncestryCalc when reproducing or comparing against AISNP-panel likelihood methods.

In the reproduction flow, FROGAncestryCalc is part of the [classical SNP superpopulation path](../getting-started/training-and-running-predictions.md#classical-snp-superpopulation-path), alongside the integrated SNP ancestry pipeline.
