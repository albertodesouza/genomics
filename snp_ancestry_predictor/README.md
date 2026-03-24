# SNP Ancestry Predictor

> Ancestry prediction from genome-wide SNP genotypes using allele frequency-based methods.

This module extracts SNP genotypes from multi-sample 1000 Genomes VCF files, converts them to the 23andMe raw data format, computes per-population allele frequency statistics, and predicts ancestry at the superpopulation or population level.

---

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1 — Generate 23andMe Files](#step-1--generate-23andme-files)
  - [Step 2 — Compute Statistics](#step-2--compute-statistics)
  - [Step 3 — Predict Ancestry](#step-3--predict-ancestry)
- [Configuration Reference](#configuration-reference)
- [Ancestry Estimation Methods](#ancestry-estimation-methods)
  - [Maximum Likelihood Classification (MLC)](#maximum-likelihood-classification-mlc)
  - [Admixture MLE](#admixture-mle)
- [Output Format](#output-format)
- [Idempotency](#idempotency)
- [Performance and Scalability](#performance-and-scalability)
- [FAQ](#faq)

---

## Overview

### What does this module do?

The **SNP Ancestry Predictor** is a three-step pipeline that:

1. **Extracts SNP genotypes** from multi-sample 1000 Genomes VCF files and writes per-individual files in the widely-used 23andMe raw data format.
2. **Computes reference allele frequency statistics** from a configurable subset of individuals (train, validation, and/or test splits).
3. **Predicts ancestry** for evaluation individuals using either Maximum Likelihood Classification (single population assignment) or Admixture MLE (mixture proportion estimation).

### Key Features

- Fully configurable via a single YAML file
- Batch VCF processing via `bcftools` for high throughput
- Supports both superpopulation (5 classes: AFR, AMR, EAS, EUR, SAS) and population (26 classes) prediction levels
- Two estimation methods: Maximum Likelihood Classification and Admixture MLE
- Idempotent execution: safely resume after interruption
- MAF filtering and Fst-based SNP selection for optimal ancestry discrimination
- Detailed evaluation metrics: accuracy, precision, recall, F1, confusion matrix
- Rich console output with progress bars

---

## Prerequisites

- **Conda `genomics` environment** (see `scripts/start_genomics_universal.sh` in the repository root)
- **Python 3.10+**
- **bcftools** (included in the genomics conda environment)
- **Python packages**: `pyyaml`, `numpy`, `scipy`, `rich`

### Data Requirements

- Multi-sample VCF files from the 1000 Genomes Project (per-chromosome, phased)
- A `splits_metadata.json` file defining train/val/test splits (produced by `neural_ancestry_predictor`)
- An individuals directory with per-sample subdirectories (produced by `build_non_longevous_dataset`)

---

## Installation

The conda `genomics` environment already includes all required dependencies. Activate it before running:

```bash
cd genomics/snp_ancestry_predictor
source ../scripts/start_genomics_universal.sh
```

If you need to install Python dependencies manually:

```bash
pip install pyyaml numpy scipy rich
```

---

## Quick Start

### 1. Activate the environment

```bash
cd genomics/snp_ancestry_predictor
source ../scripts/start_genomics_universal.sh
```

### 2. Run the full pipeline

```bash
python3 snp_ancestry_predictor.py --config configs/default.yaml
```

This will execute all three steps sequentially:
1. Generate 23andMe files for all individuals
2. Compute allele frequency statistics from the training set
3. Predict ancestry for the test set and report metrics

### 3. Run individual steps

Disable steps you don't need in the YAML:

```yaml
pipeline:
  steps:
    generate_23andme: false    # skip — files already exist
    compute_statistics: false  # skip — statistics already computed
    predict_ancestry: true     # run only prediction
```

---

## Pipeline Steps

### Step 1 — Generate 23andMe Files

**Purpose:** Extract per-individual SNP genotypes from multi-sample VCF files and write them in 23andMe format.

**How it works:**

1. Reads `splits_metadata.json` to obtain the list of all sample IDs.
2. Checks which individuals already have a 23andMe file (idempotency).
3. Queries sample names from the first VCF to verify membership.
4. For each chromosome VCF, uses `bcftools query` to extract genotypes for all pending individuals in a single pass.
5. Filters to biallelic SNPs with rsIDs (configurable); optionally restricts to a custom SNP panel.
6. Writes one `<sample_id>_23andme.txt` file per individual in their directory.

**Output:** One file per individual at `<individuals_dir>/<sample_id>/<sample_id>_23andme.txt`.

**Automatic batching:** If the operating system file-descriptor limit is lower than the number of individuals, the program automatically batches processing and re-reads VCFs for each batch. In most cases the limit is raised transparently.

### Step 2 — Compute Statistics

**Purpose:** Compute per-population allele frequencies from a reference subset.

**How it works:**

1. Selects individuals belonging to the configured reference subsets (default: train only).
2. For each chromosome VCF, uses `bcftools query` to extract genotypes for reference individuals only.
3. For each biallelic SNP with an rsID, counts the reference allele occurrences per population (or superpopulation).
4. Computes allele frequencies and applies a minimum MAF filter.
5. Optionally selects the top *N* SNPs by Fst (inter-population frequency variance).
6. Saves all statistics to a single JSON file.

**Output:** A JSON file containing:
- `metadata` — reference subsets, level, population sizes, filter parameters
- `ref_alleles` — the VCF REF allele for each retained SNP
- `snp_info` — chromosome and position for each SNP
- `allele_frequencies` — per-population REF allele frequencies

**Fst computation:**

The simplified Fst used for SNP selection is:

```
Fst(i) = Var_between(p_i) / (p_bar_i * (1 - p_bar_i))
```

where `p_bar_i` is the mean allele frequency across populations and `Var_between` is the variance of per-population frequencies. Higher Fst indicates greater differentiation across populations and therefore more ancestry-informative SNPs.

### Step 3 — Predict Ancestry

**Purpose:** Predict ancestry labels for evaluation individuals and report classification metrics.

**How it works:**

1. Loads the statistics JSON from Step 2.
2. For each evaluation individual, reads their 23andMe file (loading only the SNPs present in the statistics).
3. Applies the configured method (`mle` or `admixture_mle`) to predict ancestry.
4. Computes accuracy, precision, recall, F1-score, and confusion matrix.
5. Prints results to the console and saves them as JSON.

**Output:** A JSON file at `<results_dir>/predictions_<method>_<level>.json` containing metrics and per-individual predictions.

---

## Configuration Reference

The YAML configuration file has five sections. All paths can be absolute or relative to the working directory.

### `input`

| Field | Type | Description |
|-------|------|-------------|
| `splits_metadata` | string | Path to `splits_metadata.json` |
| `individuals_dir` | string | Directory with per-individual subdirectories |
| `vcf_pattern` | string | VCF file pattern; `{chrom}` is replaced by chromosome name |
| `chromosomes` | list | Chromosome names to process (default: chr1–chr22 + chrX) |

### `conversion`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `format_version` | string | `"V5"` | 23andMe output format: `"V5"` or `"V3"` |
| `skip_no_rsid` | bool | `true` | Skip variants without rsID |
| `snp_panel` | string/null | `null` | Path to SNP panel file (one rsID per line) |
| `output_filename` | string | `"{sample_id}_23andme.txt"` | Filename template for output files |

### `statistics`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_file` | string | — | Path for the statistics JSON output |
| `reference_subsets` | list | `["train"]` | Split subsets to use as reference |
| `min_maf` | float | `0.01` | Minimum minor allele frequency |
| `max_snps` | int/null | `500000` | Maximum SNPs (top by Fst); null = keep all |
| `snp_panel` | string/null | `null` | Optional SNP panel filter for statistics |

### `prediction`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `level` | string | `"superpopulation"` | `"superpopulation"` (5 classes) or `"population"` (26 classes) |
| `method` | string | `"mle"` | `"mle"` or `"admixture_mle"` |
| `evaluation_subsets` | list | `["test"]` | Split subsets to evaluate |
| `results_dir` | string | `"results"` | Directory for prediction output |
| `admixture_mle.n_restarts` | int | `20` | Random restarts for the admixture optimiser |

### `pipeline`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `steps.generate_23andme` | bool | `true` | Enable/disable Step 1 |
| `steps.compute_statistics` | bool | `true` | Enable/disable Step 2 |
| `steps.predict_ancestry` | bool | `true` | Enable/disable Step 3 |

---

## Ancestry Estimation Methods

### Maximum Likelihood Classification (MLC)

MLC assigns each individual to the single population that maximises the likelihood of their observed genotypes.

**Model:** Under Hardy-Weinberg equilibrium, the probability of observing genotype *g* at SNP *i* given population allele frequency *p* is:

```
P(g_i = 2 | p) = p^2                    (homozygous reference)
P(g_i = 1 | p) = 2 * p * (1 - p)        (heterozygous)
P(g_i = 0 | p) = (1 - p)^2              (homozygous alternate)
```

where *g_i* is the count of the reference allele (0, 1, or 2), and *p* is the reference allele frequency in the candidate population.

**Classification rule:** The log-likelihood for population *k* across all *N* usable SNPs is:

```
L(k) = sum_{i=1}^{N} log P(g_i | p_{i,k})
```

The predicted population is:

```
k* = argmax_k L(k)
```

**Properties:**
- Fast: vectorised computation over all SNPs and populations simultaneously.
- Produces a single hard classification.
- Works well when populations are clearly separated and individual ancestry is predominantly from one source.

### Admixture MLE

Admixture MLE estimates the proportion of ancestry from each reference population, modelling the individual as a mixture.

**Model:** Let alpha = (alpha_1, ..., alpha_K) be the admixture proportions with alpha_k >= 0 and sum_k alpha_k = 1. The mixed allele frequency at SNP *i* is:

```
p_mixed_i = sum_{k=1}^{K} alpha_k * p_{i,k}
```

The genotype probability under the mixture model is:

```
P(g_i | alpha) = P(g_i | p_mixed_i)
```

using the same Hardy-Weinberg formulas as MLC but with the mixed frequency `p_mixed_i`.

**Estimation:** The proportions are found by minimising the negative log-likelihood:

```
NLL(alpha) = - sum_{i=1}^{N} log P(g_i | p_mixed_i(alpha))
```

subject to:
- `sum_k alpha_k = 1`
- `alpha_k >= 0` for all *k*

This constrained optimisation is solved with Sequential Least-Squares Programming (SLSQP) from `scipy.optimize.minimize`, using multiple random restarts (default: 20) to mitigate local minima. The initial restart uses uniform proportions; subsequent restarts draw from a symmetric Dirichlet distribution.

**Classification:** The predicted ancestry is the population with the highest estimated proportion:

```
k* = argmax_k alpha_k
```

**Properties:**
- Provides a full vector of ancestry proportions per individual.
- More appropriate for admixed individuals.
- Slower than MLC due to iterative optimisation (mitigated by vectorised likelihood computation).
- Based on the implementation in `FROGAncestryCalc/tools/admixture_estimate.py`.

---

## Output Format

### Statistics file (Step 2)

```json
{
  "metadata": {
    "reference_subsets": ["train"],
    "level": "superpopulation",
    "populations": ["AFR", "AMR", "EAS", "EUR", "SAS"],
    "pop_sizes": {"AFR": 200, "AMR": 85, ...},
    "n_individuals": 910,
    "n_snps": 500000,
    "min_maf": 0.01,
    "max_snps": 500000,
    "created_at": "2026-03-23T..."
  },
  "ref_alleles": {"rs12345": "A", ...},
  "snp_info": {"rs12345": ["1", "12345"], ...},
  "allele_frequencies": {
    "rs12345": [0.45, 0.52, 0.61, 0.48, 0.55]
  }
}
```

The frequency arrays follow the order defined in `metadata.populations`.

### Predictions file (Step 3)

```json
{
  "metadata": {
    "method": "mle",
    "level": "superpopulation",
    "n_snps": 500000,
    "n_individuals": 195
  },
  "metrics": {
    "accuracy": 0.98,
    "weighted_precision": 0.98,
    "weighted_recall": 0.98,
    "weighted_f1": 0.98
  },
  "per_class_metrics": {
    "AFR": {"precision": 1.0, "recall": 0.98, "f1": 0.99, "support": 40}
  },
  "confusion_matrix": {
    "AFR": {"AFR": 39, "AMR": 1, "EAS": 0, "EUR": 0, "SAS": 0}
  },
  "predictions": [
    {"sample_id": "HG02577", "true": "AFR", "predicted": "AFR", "correct": true, "snps_used": 498321}
  ]
}
```

---

## Idempotency

Each step is designed to be safely re-run without recomputing completed work:

| Step | Condition to skip |
|------|-------------------|
| Step 1 | The 23andMe file for an individual already exists |
| Step 2 | The statistics JSON file already exists |
| Step 3 | Always runs (prediction is fast) |

To force recomputation:
- **Step 1:** Delete the individual's `*_23andme.txt` file(s).
- **Step 2:** Delete the statistics JSON file.

---

## Performance and Scalability

### Estimated runtimes

| Step | Typical duration | Bottleneck |
|------|-----------------|------------|
| Step 1 (1300 individuals, 23 chromosomes) | 2–4 hours | VCF I/O + bcftools decompression |
| Step 2 (910 reference individuals, 23 chromosomes) | 1–2 hours | VCF I/O |
| Step 3 (195 test individuals, 500K SNPs) | 10–30 minutes | 23andMe file reading |

### Memory considerations

- **Step 1:** Memory usage is minimal; genotypes are streamed and written to files.
- **Step 2:** Allele frequency counters are kept in Python dictionaries. For all biallelic SNPs with rsIDs (~35M), this requires approximately 15–25 GB. Setting `max_snps` has no effect on peak memory during accumulation, but setting `statistics.snp_panel` does.
- **Step 3:** Only the SNPs in the statistics file are loaded per individual (~500K entries). Memory usage is negligible.

### Reducing runtime

- Use a **SNP panel** (`conversion.snp_panel` or `statistics.snp_panel`) to restrict to a predefined set of ancestry-informative SNPs.
- Set `max_snps` to 100,000–500,000 (sufficient for accurate superpopulation classification).
- Reduce the number of chromosomes in `input.chromosomes` for quick experiments.

---

## FAQ

### How many SNPs do I need for accurate ancestry prediction?

For superpopulation-level classification (5 classes), as few as 1,000 highly informative SNPs can achieve >95% accuracy. The default `max_snps: 500000` provides excellent accuracy with manageable file sizes.

### Can I use a custom set of ancestry-informative SNPs (AISNPs)?

Yes. Create a text file with one rsID per line and set both `conversion.snp_panel` and `statistics.snp_panel` to its path. For example, the 55 or 128 AISNPs from FROGAncestryCalc.

### What is the difference between `mle` and `admixture_mle`?

- **`mle`** assigns each individual to one population — the one with the highest genotype likelihood. Fast and simple.
- **`admixture_mle`** estimates the proportion of ancestry from each population. Better for admixed individuals, but slower.

### Can I predict at the population level (26 classes)?

Yes. Set `prediction.level: "population"` in the YAML. Note that population-level prediction is harder due to smaller reference group sizes and closer genetic distances between some populations.

### How do I add more reference individuals?

Modify `statistics.reference_subsets` to include more splits:

```yaml
statistics:
  reference_subsets: ["train", "val"]
```

### What if `bcftools` is not found?

Activate the genomics conda environment before running:

```bash
source ../scripts/start_genomics_universal.sh
```

### What does "batched mode" mean in Step 1?

If your OS limits the number of simultaneously open files (typically 1024), the program automatically processes individuals in batches. This is slower because VCFs are re-read for each batch, but it runs correctly without manual intervention.

### How do I force Step 2 to recompute?

Delete the statistics file:

```bash
rm /dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000/snp_ancestry_statistics.json
```

---

## File Structure

```
snp_ancestry_predictor/
├── snp_ancestry_predictor.py      # Main pipeline script
├── configs/
│   └── default.yaml               # Default configuration
└── README.md                      # This documentation
```

---

## References

- 1000 Genomes Project: http://www.internationalgenome.org/
- bcftools: https://samtools.github.io/bcftools/
- Hardy-Weinberg principle: Fisher, R.A. (1918)
- Admixture estimation: Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.
- Ancestry-informative SNPs: Kidd, K.K. et al. (2014). Progress toward an efficient panel of SNPs for ancestry inference. *Forensic Science International: Genetics*, 10, 23–32.

---

**Author:** Alberto F. De Souza  
**Date:** 2026-03-23  
**Version:** 1.0
