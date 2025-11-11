# 🚀 Quick Start: Non-Longevous Dataset Builder

## ⚡ Quick Test (5 minutes)

### Step 1: Analyze Metadata

Run the program with the configured CSV:

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

You'll see statistics about:
- 5 superpopulations (AFR, AMR, EAS, EUR, SAS)
- 10 populations
- 56 total individuals
- Sex distribution in each population

### Step 2: Configure For Your Project

1. **Prepare your CSV with 1000 Genomes metadata**:
   - Download from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/
   - Or use the provided example for testing

2. **Edit `configs/default.yaml`**:

```yaml
data_sources:
  # Your CSV with metadata (path relative to configs/ directory)
  metadata_csv: "../../doc/your_file.csv"
  
  # GRCh38 reference (path relative to configs/ directory)
  reference:
    fasta: "../../path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  
  # VCF pattern (per chromosome)
  vcf_pattern: "/path/to/vcfs/1kGP_high_coverage.{chrom}.vcf.gz"

sample_selection:
  # Choose: "superpopulation" or "population"
  level: "superpopulation"
  
  # How many samples per group
  samples_per_group: 2
  
  # Sex filter: "all", "male", or "female"
  sex_filter: "all"

build_window_params:
  # Window mode: "gene" or "snp"
  mode: "gene"
  
  # Gene mode settings
  gene:
    symbol: "CYP2B6"         # Gene of interest
    # OR use gene_list_file: "my_genes.txt" for multiple genes
  
  # SNP mode settings (alternative to gene mode)
  # snp:
  #   snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"
  
  # Run AlphaGenome predictions
  predict: true
  
  # Output types
  outputs: "RNA_SEQ,ATAC"
  
  # Tissues/cells (CURIEs)
  ontology: "UBERON:0002107,UBERON:0000955"  # liver, brain

pipeline:
  steps:
    analyze_metadata: true    # ✓ Analysis
    select_samples: true      # ✓ Selection
    run_predictions: true     # ✓ Run build_window_and_predict.py
    generate_report: true     # ✓ Final report
```

3. **Configure AlphaGenome API Key** (if using predictions):

```bash
export ALPHAGENOME_API_KEY="your_key_here"
```

### Step 3: Run Complete Pipeline

```bash
cd build_non_longevous_dataset
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

## 📊 Expected Outputs

```
non_longevous_results/
├── metadata_statistics.json              # CSV statistics
├── selected_samples.csv                  # Selected samples
├── non_longevous_dataset_checkpoint.json # Checkpoint (idempotence)
├── processing_summary.txt                # Final report
│
├── HG00096__CYP2B6/                      # Results per sample
│   ├── ref.window.fa                     # Reference
│   ├── HG00096.H1.window.fixed.fa        # Haplotype 1
│   ├── HG00096.H2.window.fixed.fa        # Haplotype 2
│   ├── predictions_H1/                   # H1 predictions
│   │   ├── rna_seq.npz
│   │   ├── rna_seq_metadata.json
│   │   ├── atac.npz
│   │   └── atac_metadata.json
│   └── predictions_H2/                   # H2 predictions
│       └── ...
│
└── HG00097__CYP2B6/                      # Next sample
    └── ...
```

## 🎯 Common Use Cases

### Case 1: Compare African and European Populations

```yaml
sample_selection:
  level: "superpopulation"
  samples_per_group: 10
  include_groups: ["AFR", "EUR"]
```

### Case 2: Only Females from Specific Populations

```yaml
sample_selection:
  level: "population"
  samples_per_group: 5
  include_groups: ["GBR", "CHB", "YRI"]
  sex_filter: "female"
```

### Case 3: All Superpopulations (Balanced)

```yaml
sample_selection:
  level: "superpopulation"
  samples_per_group: 20
  include_groups: []  # all
  sex_filter: "all"
```

### Case 4: Only Extract Sequences (No Predictions)

```yaml
build_window_params:
  mode: "gene"
  gene:
    symbol: "BRCA1"
  predict: false       # Disable AlphaGenome
  skip_h2: false       # Keep both haplotypes
```

### Case 5: AISNP Ancestry Analysis (SNP Mode)

```yaml
build_window_params:
  mode: "snp"
  snp:
    snp_list_file: "../FROGAncestryCalc/SNPInfo/55_aisnps_alleles_grch38.txt"
  predict: true
  outputs: "ATAC,CHIP_HISTONE"
  ontology: "UBERON:0002107"  # liver
```

**Result**: Creates 55 windows (one per AISNP) for each sample, ideal for comparing epigenetic profiles across ancestries.

### Case 6: Multiple Genes in Batch

```yaml
build_window_params:
  mode: "gene"
  gene:
    gene_list_file: "cancer_genes.txt"  # One gene per line
  predict: true
```

**cancer_genes.txt**:
```
BRCA1
BRCA2
TP53
PTEN
```

## 🔧 Troubleshooting

### Error: CSV not found

```
[ERROR] CSV file not found
```

**Solution**: Check the path in `data_sources.metadata_csv`

### Error: VCF pattern with {chrom}

```
[WARN] VCF pattern contains {chrom}, but chromosome was not determined
```

**Solution**: Currently, the program assumes you know which chromosome contains your gene. For CYP2B6 gene (chromosome 19), configure:

```yaml
vcf_pattern: "/path/to/1kGP_high_coverage.chr19.vcf.gz"
```

Alternatively, let `build_window_and_predict.py` determine it automatically by downloading the GTF.

### Error: AlphaGenome API Key

```
RuntimeError: AlphaGenome API key not provided
```

**Solution**:
```bash
export ALPHAGENOME_API_KEY="your_key"
# Or add to ~/.bashrc
```

### Error: SNP file not found

```
FileNotFoundError: [Errno 2] No such file or directory: '55_aisnps_alleles_grch38.txt'
```

**Solution**: Make sure the SNP file exists. If using SNP mode with FROGAncestryCalc, generate the GRCh38 file:
```bash
cd FROGAncestryCalc
python3 tools/convert_grch37_to_grch38.py
```

## 💡 Performance Tips

1. **Start small**: Use `samples_per_group: 1` or `2` for testing
2. **Disable H2**: Use `skip_h2: true` to process 2x faster
3. **Fewer ontologies**: Use 1-3 specific tissues instead of all
4. **Parallelization**: Configure `n_workers: 8` (adjust for your CPU)
5. **Checkpoint**: If interrupted, the program continues from where it stopped
6. **SNP mode**: Much faster when you don't need full genes, only specific positions

## 📚 Complete Documentation

See `README.md` for detailed documentation.

## 🧬 Common Genes for Analysis

- **CYP2B6** (chr19): Drug metabolism
- **APOE** (chr19): Alzheimer's risk
- **BRCA1** (chr17): Breast cancer
- **BRCA2** (chr13): Breast/ovarian cancer
- **FOXO3** (chr6): Longevity
- **TP53** (chr17): Tumor suppressor

## 🎓 Complete Workflow Example

```bash
# 1. Enter module directory
cd build_non_longevous_dataset

# 2. Analyze available data
python3 build_non_longevous_dataset.py --config configs/default.yaml

# 3. Edit configuration based on statistics
nano configs/default.yaml

# 4. Enable additional steps (edit YAML)
# select_samples: true
# run_predictions: true
# generate_report: true

# 5. Run complete pipeline
python3 build_non_longevous_dataset.py --config configs/default.yaml

# 6. Check results
ls -lh ../non_longevous_results/
cat ../non_longevous_results/processing_summary.txt

# 7. If interrupted, simply run again
# The checkpoint will ensure already processed samples are skipped
python3 build_non_longevous_dataset.py --config configs/default.yaml
```

## ✅ Preparation Checklist

Before running the complete pipeline:

- [ ] CSV with metadata prepared
- [ ] GRCh38 reference genome (.fa + .fai) available
- [ ] 1000 Genomes VCFs downloaded
- [ ] VCFs indexed (.tbi)
- [ ] AlphaGenome API key configured (if using predictions)
- [ ] Sufficient disk space (~500MB-2GB per sample, more for SNP mode)
- [ ] YAML configuration reviewed and adjusted
- [ ] For SNP mode: SNP list file exists and is in GRCh38 coordinates

---

**Ready!** You're prepared to build your non-longevous dataset! 🎉

For SNP/AISNP analysis, also see: [AISNP Mode Documentation](AISNP_MODE.md)
