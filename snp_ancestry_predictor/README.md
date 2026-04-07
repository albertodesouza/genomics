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
  - [Admixture EM](#admixture-em)
  - [Admixture MLE](#admixture-mle)
- [Output Format](#output-format)
  - [Per-individual result files](#per-individual-result-files)
- [Plotting Ancestry Pie Charts](#plotting-ancestry-pie-charts)
- [Computing Gene-Window SNP Panels](#computing-gene-window-snp-panels)
- [Idempotency](#idempotency)
- [Performance and Scalability](#performance-and-scalability)
- [FAQ](#faq)

---

## Overview

### What does this module do?

The **SNP Ancestry Predictor** is a three-step pipeline that:

1. **Extracts SNP genotypes** from multi-sample 1000 Genomes VCF files and writes per-individual files in the widely-used 23andMe raw data format.
2. **Computes reference allele frequency statistics** from a configurable subset of individuals (train, validation, and/or test splits).
3. **Predicts ancestry** for evaluation individuals using Maximum Likelihood Classification (single population assignment), Admixture EM (fast mixture proportion estimation), or Admixture MLE (numerical optimisation-based mixture estimation).

### Key Features

- Fully configurable via a single YAML file
- Batch VCF processing via `bcftools` for high throughput
- Supports both superpopulation (5 classes: AFR, AMR, EAS, EUR, SAS) and population (26 classes) prediction levels
- Three estimation methods: Maximum Likelihood Classification, Admixture EM (fast), and Admixture MLE
- Idempotent execution: safely resume after interruption
- MAF filtering and Fst-based SNP selection for optimal ancestry discrimination
- Detailed evaluation metrics: accuracy, precision, recall, F1, confusion matrix
- Per-individual result files with ancestry proportions (JSON)
- Pie chart visualization via `plot_ancestry_pie.py`
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

### 4. Predict ancestry for a single individual

Once the statistics have been computed (Step 2), you can predict ancestry for any individual by passing their 23andMe file directly:

```bash
python3 snp_ancestry_predictor.py \
    --config configs/default.yaml \
    --individual /path/to/sample_23andme.txt
```

This mode:
- Skips Steps 1, 2, and the batch Step 3 entirely
- Loads pre-computed statistics from `statistics.output_dir` (the filename is auto-generated from the active parameters)
- Runs the configured method (`mle`, `admixture_em`, or `admixture_mle`) on the provided file
- Displays a table of ancestry proportions and the predicted population
- Saves the result as a JSON file in the current working directory
- The individual does not need to be present in `splits_metadata.json`

You can then generate a pie chart from the saved result:

```bash
python3 plot_ancestry_pie.py \
    /dados/.../ancestry_results/individuals/sample_admixture_em_superpopulation.json
```

---

## Pipeline Steps

### Step 1 — Generate 23andMe Files

**Purpose:** Extract per-individual SNP genotypes from multi-sample VCF files and write them in 23andMe format.

**How it works:**

1. Reads `splits_metadata.json` to obtain the list of all sample IDs.
2. Checks which individuals already have a 23andMe file (idempotency).
3. Queries sample names from the first VCF to verify membership.
4. When `annotate_dbsnp: true`, downloads the dbSNP common variants VCF (once, cached in `ref_dir`) and pipes each chromosome VCF through `bcftools annotate -c ID` to add rsIDs on-the-fly, without creating intermediate files on disk.
5. For each chromosome VCF, uses `bcftools query` to extract genotypes for all pending individuals in a single pass. Filters to biallelic SNPs with proper rsIDs (`ID~"^rs"`). Optionally restricts to a chip panel or custom SNP panel file.
6. Writes one `<sample_id>_23andme.txt` file per individual in their directory.
7. Reports the number of panel SNPs found per individual and aborts if any individual has zero SNPs.

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

$$F_{ST}(i) = \frac{\mathrm{Var}_{\text{between}}(p_i)}{\bar{p}_i \,(1 - \bar{p}_i)}$$

where $\bar{p}_i$ is the mean allele frequency across populations and $\mathrm{Var}_{\text{between}}$ is the variance of per-population frequencies. Higher $F_{ST}$ indicates greater differentiation across populations and therefore more ancestry-informative SNPs.

### Step 3 — Predict Ancestry

**Purpose:** Predict ancestry labels for evaluation individuals and report classification metrics.

**How it works:**

1. Loads the statistics JSON from Step 2.
2. For each evaluation individual, reads their 23andMe file (loading only the SNPs present in the statistics).
3. Applies the configured method (`mle`, `admixture_em`, or `admixture_mle`) to predict ancestry.
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
| `skip_no_rsid` | bool | `true` | Skip variants without rsID (only keep variants with `rs` IDs) |
| `annotate_dbsnp` | bool | `false` | Annotate VCFs with rsIDs from dbSNP on-the-fly (required when VCFs lack rsIDs) |
| `genome_build` | string | `"GRCh38"` | Genome build for dbSNP download (`"GRCh37"` or `"GRCh38"`) |
| `filter_by_chip_panel` | bool | `false` | Auto-download and apply the Illumina chip panel matching `format_version` |
| `ref_dir` | string | `"refs"` | Directory where downloaded reference files (dbSNP, chip panel) are cached |
| `snp_panel` | string/null | `null` | Path to custom SNP panel file (one rsID per line); overrides `filter_by_chip_panel` |
| `output_filename` | string | `"{sample_id}_23andme.txt"` | Filename template for output files |

### `statistics`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `output_dir` | string | — | Directory for statistics JSON files. The filename is auto-generated from `reference_subsets`, `level`, `min_maf`, `max_snps`, and `snp_panel`, so different configurations produce distinct files that are reused or recomputed automatically. |
| `reference_subsets` | list | `["train"]` | Split subsets to use as reference |
| `min_maf` | float | `0.01` | Minimum minor allele frequency |
| `max_snps` | int/null | `500000` | Maximum SNPs (top by Fst); null = keep all |
| `snp_panel` | string/null | `null` | Optional SNP panel filter for statistics |

### `prediction`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `level` | string | `"superpopulation"` | `"superpopulation"` (5 classes) or `"population"` (26 classes) |
| `haplotype_mode` | string | `"H1+H2"` | `"H1+H2"` (diploid, both haplotypes), `"H1"` (first haplotype only), or `"H2"` (second haplotype only). Affects Steps 2 and 3. Use `"H1"` to match `neural_ancestry_predictor`'s default for a fair comparison. |
| `method` | string | `"admixture_em"` | `"mle"`, `"admixture_em"`, or `"admixture_mle"` |
| `evaluation_subsets` | list | `["test"]` | Split subsets to evaluate |
| `results_dir` | string | `"results"` | Directory for prediction output |
| `admixture_em.max_iter` | int | `1000` | Maximum EM iterations per individual |
| `admixture_em.tol` | float | `1e-7` | Convergence tolerance on proportions |
| `admixture_mle.n_restarts` | int | `20` | Random restarts for the L-BFGS-B optimiser |

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

**Model:** Under Hardy-Weinberg equilibrium, the probability of observing genotype $g_i$ at SNP $i$ given population allele frequency $p$ is:

$$P(g_i = 2 \mid p) = p^2 \qquad \text{(homozygous reference)}$$

$$P(g_i = 1 \mid p) = 2\,p\,(1 - p) \qquad \text{(heterozygous)}$$

$$P(g_i = 0 \mid p) = (1 - p)^2 \qquad \text{(homozygous alternate)}$$

where $g_i \in \{0, 1, 2\}$ is the count of the tracked allele, and $p$ is its frequency in the candidate population.

**Classification rule:** The log-likelihood for population $k$ across all $N$ usable SNPs is:

$$\mathcal{L}(k) = \sum_{i=1}^{N} \log P(g_i \mid p_{i,k})$$

The predicted population is:

$$k^* = \arg\max_k \; \mathcal{L}(k)$$

**Properties:**
- Fast: vectorised computation over all SNPs and populations simultaneously.
- Produces a single hard classification.
- Works well when populations are clearly separated and individual ancestry is predominantly from one source.

### Admixture EM

Admixture EM estimates ancestry proportions using the classical Expectation-Maximization algorithm from FRAPPE (Tang et al., 2005). This is the **recommended** method for admixture estimation — it is much faster than `admixture_mle` and produces equivalent results.

**Model:** The same mixture model as Admixture MLE (below). Let $\boldsymbol{\alpha} = (\alpha_1, \dots, \alpha_K)$ be the admixture proportions with $\alpha_k \ge 0$ and $\sum_k \alpha_k = 1$. The mixed allele frequency at SNP $j$ is:

$$p_j^{\text{mix}} = \sum_{k=1}^{K} \alpha_k \, f_{kj}$$

**Algorithm:** Each allele copy is treated as belonging latently to one of the $K$ populations. The EM algorithm iterates:

**E-step** — compute the posterior responsibility that allele copies at SNP $j$ came from population $k$:

$$r_{jk}^{\text{ref}} = \frac{\alpha_k \, f_{kj}}{p_j^{\text{mix}}} \qquad r_{jk}^{\text{alt}} = \frac{\alpha_k \, (1 - f_{kj})}{1 - p_j^{\text{mix}}}$$

**M-step** — update the proportions from the expected allele counts:

$$\alpha_k^{\text{new}} = \frac{1}{2J} \sum_{j=1}^{J} \left[ g_j \cdot r_{jk}^{\text{ref}} + (2 - g_j) \cdot r_{jk}^{\text{alt}} \right]$$

where $g_j \in \{0, 1, 2\}$ is the genotype dose and $J$ is the number of SNPs.

Iteration continues until $\max_k |\alpha_k^{\text{new}} - \alpha_k| < \text{tol}$ or `max_iter` is reached.

**Classification:** The predicted ancestry is $k^* = \arg\max_k \; \alpha_k$.

**Properties:**
- **Fast:** Each iteration is a few numpy matrix–vector multiplies with no scipy overhead. Typically converges in 50–200 iterations.
- **Monotonic:** The likelihood is guaranteed to increase (or stay the same) at every step.
- **No restarts needed:** Initialises at uniform proportions; convergence to the global optimum is typical for this problem.
- Provides the same full vector of ancestry proportions as `admixture_mle`.

### Admixture MLE

Admixture MLE estimates the proportion of ancestry from each reference population using numerical optimisation. It uses the same mixture model as Admixture EM but solves it with L-BFGS-B via `scipy.optimize.minimize`.

**Model:** Let $\boldsymbol{\alpha} = (\alpha_1, \dots, \alpha_K)$ be the admixture proportions with $\alpha_k \ge 0$ and $\sum_k \alpha_k = 1$. The mixed allele frequency at SNP $i$ is:

$$p_i^{\text{mix}} = \sum_{k=1}^{K} \alpha_k \, p_{i,k}$$

The genotype probability under the mixture model is:

$$P(g_i \mid \boldsymbol{\alpha}) = P\!\left(g_i \mid p_i^{\text{mix}}\right)$$

using the same Hardy-Weinberg formulas as MLC but with the mixed frequency $p_i^{\text{mix}}$.

**Estimation:** The proportions are found by minimising the negative log-likelihood:

$$\mathrm{NLL}(\boldsymbol{\alpha}) = -\sum_{i=1}^{N} \log P\!\left(g_i \mid p_i^{\text{mix}}(\boldsymbol{\alpha})\right)$$

The simplex constraint is enforced via softmax reparametrisation ($\alpha = \text{softmax}(\theta)$), allowing unconstrained L-BFGS-B optimisation with analytical gradient. Multiple random restarts (default: 20) mitigate local minima, with early stopping after 3 consecutive non-improving restarts.

**Classification:** The predicted ancestry is the population with the highest estimated proportion:

$$k^* = \arg\max_k \; \alpha_k$$

**Properties:**
- Provides a full vector of ancestry proportions per individual.
- More appropriate for admixed individuals.
- **Significantly slower** than Admixture EM due to scipy overhead and multiple restarts. Use `admixture_em` instead unless you have a specific reason to prefer numerical optimisation.

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
    {
      "sample_id": "HG02577", "true": "AFR", "predicted": "AFR",
      "correct": true, "snps_used": 498321,
      "proportions": {"AFR": 0.9834, "AMR": 0.0052, "EAS": 0.0011, "EUR": 0.0068, "SAS": 0.0035}
    }
  ]
}
```

All three methods include a `proportions` field. For `mle`, proportions are computed via softmax normalization of the log-likelihoods. For `admixture_em` and `admixture_mle`, they are the directly estimated admixture fractions.

### Per-individual result files

In addition to the aggregate predictions file, Step 3 saves a separate JSON file for each individual in the `individuals/` subdirectory of `results_dir`:

```
<results_dir>/individuals/<sample_id>_<method>_<level>.json
```

Example:

```json
{
  "sample_id": "HG02577",
  "true": "AFR",
  "predicted": "AFR",
  "correct": true,
  "snps_used": 498321,
  "proportions": {"AFR": 0.9834, "AMR": 0.0052, "EAS": 0.0011, "EUR": 0.0068, "SAS": 0.0035},
  "method": "mle",
  "level": "superpopulation"
}
```

These files can be used directly with `plot_ancestry_pie.py` to generate pie chart visualizations (see [Plotting Ancestry Pie Charts](#plotting-ancestry-pie-charts)).

---

## Plotting Ancestry Pie Charts

The `plot_ancestry_pie.py` script generates a pie chart from a per-individual result JSON file produced by Step 3.

### Usage

```bash
source ../scripts/start_genomics_universal.sh
python3 plot_ancestry_pie.py <result_json> [-o output.png] [--title "..."]
```

### Examples

Plot ancestry for a specific individual:

```bash
python3 plot_ancestry_pie.py \
    /dados/.../ancestry_results/individuals/HG02577_mle_superpopulation.json
```

Specify a custom output path and title:

```bash
python3 plot_ancestry_pie.py \
    /dados/.../ancestry_results/individuals/HG02577_admixture_mle_superpopulation.json \
    -o /home/user/HG02577_ancestry.png \
    --title "HG02577 Admixture Ancestry"
```

### Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `result_json` | *(required)* | Path to the per-individual JSON file |
| `-o`, `--output` | `<input>.png` | Output PNG file path |
| `--title` | Auto-generated | Custom chart title |
| `--min-percent` | `0.5` | Populations below this % are grouped into "Other" |
| `--dpi` | `150` | Output resolution |

### Colour Scheme

At the superpopulation level, the chart uses the canonical 1000 Genomes Project colours:

| Population | Colour |
|-----------|--------|
| AFR | Orange `#E8832A` |
| AMR | Red `#ED1E24` |
| EAS | Green `#108C44` |
| EUR | Blue `#2D59A4` |
| SAS | Purple `#6F3198` |

At the population level, colours are assigned automatically from the `tab20` palette.

---

## Computing Gene-Window SNP Panels

The `compute_snp_panel.py` script generates a **subset** of an existing SNP panel (e.g. the 23andMe V5 panel) by retaining only SNPs whose genomic positions fall within the central windows of genes specified in a `neural_ancestry_predictor` configuration.

This is useful for restricting the ancestry predictor to the same genomic regions used by the neural model, enabling a direct comparison of both approaches on the same data.

### Usage

```bash
python3 compute_snp_panel.py --config configs/compute_snp_panel.yaml
```

### How it works

1. Reads the `neural_ancestry_predictor` config to obtain the gene list and `window_center_size`
2. Loads per-gene genomic coordinates from `individual_metadata.json` in the dataset
3. Computes the central window of `window_center_size` bp for each gene
4. Maps rsIDs from the input panel to genomic positions using an existing 23andMe file
5. Filters the panel to SNPs within gene windows and prints per-gene counts

### Configuration

| Field | Type | Description |
|-------|------|-------------|
| `neural_config` | string | Path to the `neural_ancestry_predictor` YAML config (defines `genes_to_use` and `window_center_size`). Relative paths are resolved from the config file's directory. |
| `dataset_dir` | string | Dataset directory containing `individuals/` with `individual_metadata.json` files (for gene coordinates) |
| `individuals_dir` | string | Directory with per-individual 23andMe files from Step 1 (defaults to `dataset_dir/individuals`) |
| `input_panel` | string | Path to the input SNP panel to filter (one rsID per line) |
| `output_panel` | string | Output path for the filtered SNP panel |

### Example output

```
Neural config: ../../neural_ancestry_predictor/configs/genes_1000.yaml
  Genes: 11
  Window center size: 32,768 bp

Gene windows (32,768 bp center):
  DDB1        chr11:61,298,900-61,331,667
  EDAR        chr2:108,929,778-108,962,545
  ...

SNPs per gene:
  Gene        Chrom   Window                         SNPs
  ----------  ------  ---------------------------  ------
  DDB1        chr11    61,298,900-  61,331,667         42
  EDAR        chr2    108,929,778- 108,962,545         38
  ...
                                        TOTAL         350

Output panel saved: /dados/.../panel_23andme_V5_genes_1000_w32768.txt
  350 SNPs from 11 genes
```

After generating the panel, point `statistics.snp_panel` in the main `default.yaml` to the output file. Step 1 does **not** need to be re-run because the new panel is a subset of the V5 panel already used in Step 1.

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
| Step 3 — mle (195 test individuals, 500K SNPs) | 2–5 minutes | 23andMe file reading |
| Step 3 — admixture_em (195 test individuals, 500K SNPs) | 5–15 minutes | EM iterations |
| Step 3 — admixture_mle (195 test individuals, 500K SNPs) | 30–90 minutes | L-BFGS-B optimisation |

### Memory considerations

- **Step 1:** Memory usage is minimal; genotypes are streamed and written to files.
- **Step 2:** Allele frequency counters are kept in Python dictionaries. For all biallelic SNPs with rsIDs (~35M), this requires approximately 15–25 GB. Setting `max_snps` has no effect on peak memory during accumulation, but setting `statistics.snp_panel` does.
- **Step 3:** Only the SNPs in the statistics file are loaded per individual (~500K entries). Memory usage is negligible.

### Reducing runtime

- Enable **chip panel filtering** (`conversion.filter_by_chip_panel: true`) to restrict 23andMe files to the ~640K SNPs on the Illumina chip panel, dramatically reducing file size and disk usage.
- Use a custom **SNP panel** (`conversion.snp_panel` or `statistics.snp_panel`) to further restrict to a specific set of ancestry-informative SNPs.
- Set `max_snps` to 100,000–500,000 (sufficient for accurate superpopulation classification).
- Reduce the number of chromosomes in `input.chromosomes` for quick experiments.

---

## FAQ

### How many SNPs do I need for accurate ancestry prediction?

For superpopulation-level classification (5 classes), as few as 1,000 highly informative SNPs can achieve >95% accuracy. The default `max_snps: 500000` provides excellent accuracy with manageable file sizes.

### Why do I need `annotate_dbsnp`?

Some VCF datasets (e.g., 1000 Genomes NYGC high-coverage) use positional IDs like `1:10399:C:A` instead of standard rsIDs. Since the pipeline filters by rsID (`ID~"^rs"`), these VCFs produce zero matches. When `annotate_dbsnp: true`, each chromosome VCF is piped through `bcftools annotate -a dbsnp.vcf.gz -c ID` to add rsIDs from the NCBI dbSNP common variants database on-the-fly. The dbSNP VCF (~3 GB) is downloaded once, chromosome names are remapped if needed, and the result is cached in `ref_dir`. No intermediate annotated VCF files are created on disk — the annotation is streamed via a pipe.

### What is chip panel filtering?

When `conversion.filter_by_chip_panel` is `true`, the program automatically downloads the Illumina chip manifest matching `format_version` (e.g., Infinium Global Screening Array v1.0 for V5, ~640K rsIDs), extracts the rsIDs, and caches the panel file in `ref_dir`. Only SNPs present in the chip panel are written to the 23andMe files, producing files of ~20 MB per individual instead of hundreds of MB.

### Can I use a custom set of ancestry-informative SNPs (AISNPs)?

Yes. Create a text file with one rsID per line and set `conversion.snp_panel` to its path (this overrides `filter_by_chip_panel`). You can also set `statistics.snp_panel` for Step 2. For example, the 55 or 128 AISNPs from FROGAncestryCalc.

### What is the difference between `mle`, `admixture_em`, and `admixture_mle`?

- **`mle`** assigns each individual to one population — the one with the highest genotype likelihood. Fast and simple, but proportions are not meaningful (essentially 100%/0% with many SNPs).
- **`admixture_em`** estimates ancestry proportions via the EM algorithm. **Recommended** — fast and produces meaningful proportions for pie charts.
- **`admixture_mle`** estimates ancestry proportions via L-BFGS-B numerical optimisation. Same model as `admixture_em` but significantly slower due to scipy overhead and multiple restarts.

### How do I compare results with `neural_ancestry_predictor`?

Set `prediction.haplotype_mode: "H1"` so that the SNP predictor uses only the first haplotype, matching the default `haplotype_mode: "H1"` in the neural config. You can also use `compute_snp_panel.py` to restrict the SNP panel to the same gene windows used by the neural model.

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
├── plot_ancestry_pie.py           # Pie chart generator for individual results
├── compute_snp_panel.py           # Gene-window SNP panel generator
├── configs/
│   ├── default.yaml               # Default configuration
│   └── compute_snp_panel.yaml     # Configuration for compute_snp_panel.py
└── README.md                      # This documentation
```

---

## References

- 1000 Genomes Project: http://www.internationalgenome.org/
- bcftools: https://samtools.github.io/bcftools/
- Hardy-Weinberg principle: Fisher, R.A. (1918)
- FRAPPE (EM algorithm): Tang, H., Peng, J., Wang, P. & Risch, N.J. (2005). Estimation of individual admixture: analytical and study design considerations. *Genetic Epidemiology*, 28(4), 289–301.
- Admixture estimation: Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.
- Ohana: Cheng, J.Y., Mailund, T. & Nielsen, R. (2017). Fast admixture analysis and population tree estimation for SNP and NGS data. *Bioinformatics*, 33(14), 2148–2155.
- Ancestry-informative SNPs: Kidd, K.K. et al. (2014). Progress toward an efficient panel of SNPs for ancestry inference. *Forensic Science International: Genetics*, 10, 23–32.

---

**Author:** Alberto F. De Souza  
**Date:** 2026-03-23  
**Version:** 1.0
