# Validate Pigmentation Hypothesis

> **🧬 Systematic Validation of the Pigmentation Hypothesis**

This document describes the `validate_pigmentation_hypothesis.py` script, which validates whether the DeepLIFT → VEP pipeline correctly identifies genetic variants associated with pigmentation in African populations.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Line](#command-line)
- [Validation Levels](#validation-levels)
- [Gene Database](#gene-database)
- [Known Pigmentation SNPs](#known-pigmentation-snps)
- [Scoring System](#scoring-system)
- [Output Report](#output-report)
- [gnomAD Integration](#gnomad-integration)
- [Application to Longevity](#application-to-longevity)
- [Examples](#examples)
- [References](#references)

---

## Overview

### What does this script do?

The `validate_pigmentation_hypothesis.py` validates the hypothesis that:

> *"Using AlphaGenome output processed by `neural_ancestry_predictor.py` and `annotate_deeplift_windows.py`, it is possible to identify what in the DNA of Africans causes them to have greater skin pigmentation."*

### Three-Level Validation

```
┌─────────────────────────────────────────────────────────────┐
│                    HYPOTHESIS VALIDATION                     │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  LEVEL 1: GENES                                              │
│  ├─ Do detected genes correspond to pigmentation genes?      │
│  ├─ How many are primary genes? Regulatory? UV response?     │
│  └─ What is the overlap with Crawford et al., 2017 (Science)?│
│                                                              │
│  LEVEL 2: VARIANTS                                           │
│  ├─ Which known pigmentation rsIDs were found?               │
│  ├─ How many HIGH/MODERATE impact variants exist?            │
│  └─ How many missense variants were detected?                │
│                                                              │
│  LEVEL 3: MECHANISM                                          │
│  ├─ Do the variants make biological sense?                   │
│  ├─ stop_gained → truncated protein                          │
│  ├─ splice_variant → abnormal splicing                       │
│  └─ missense → amino acid change                             │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## Installation

### Dependencies

```bash
# Required dependency (already included in Python)
# No additional installation needed for basic usage

# Optional: for gnomAD population frequencies
pip install requests
```

### Required Files

The script requires output files from `annotate_deeplift_windows.py`:

```
input_dir/
├── summary_by_gene.tsv           # REQUIRED
├── vep_annotated_variants.tsv    # Recommended
└── ... (other files)
```

---

## Quick Start

### Basic Usage

```bash
cd neural_ancestry_predictor

# Validate annotate_deeplift_windows.py results
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2
```

### With Population Frequencies (gnomAD)

```bash
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 --gnomad
```

### Specifying Output File

```bash
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 --output my_report.md
```

---

## Command Line

### Syntax

```bash
python3 validate_pigmentation_hypothesis.py <input_dir> [options]
```

### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `input_dir` | Required | Directory with `annotate_deeplift_windows.py` outputs |
| `--output`, `-o` | Optional | Output report path (default: `<input_dir>/pigmentation_validation.md`) |
| `--gnomad` | Flag | Fetch population frequencies from gnomAD API |

### Examples

```bash
# Basic usage
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2

# With gnomAD (slower, requires internet)
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 --gnomad

# Custom output
python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 -o results/validation.md
```

---

## Validation Levels

### Level 1: Gene Validation

Verifies whether genes detected by DeepLIFT correspond to known pigmentation genes in the literature.

**Calculated metrics:**
- % of genes that are primary pigmentation genes
- % of genes that are regulatory pigmentation genes
- % of genes associated with UV response
- % overlap with Crawford et al., 2017 (Science)

**Example output:**

```
=== Level 1: Gene Validation ===
  Matched pigmentation genes: ['SLC24A5', 'HERC2', 'DDB1', 'OCA2', 'TYR']
  Unmatched genes: []
  Crawford et al. 2017 overlap: 4/5 (80.0%)
```

### Level 2: Variant Validation

Verifies whether found rsIDs include known pigmentation SNPs and analyzes the functional impact of variants.

**Calculated metrics:**
- Total rsIDs found
- Known pigmentation rsIDs (rs1426654, rs16891982, etc.)
- HIGH impact variant count
- MODERATE impact variant count
- Missense variant count

**Example output:**

```
=== Level 2: Variant Validation ===
  Total rsIDs found: 67
  Known pigmentation rsIDs: ['rs1426654']
  HIGH impact variants: 12
  MODERATE impact variants: 70
  Missense variants: 71
```

### Level 3: Mechanism Validation

Analyzes whether high-impact variants make biological sense in the context of pigmentation.

**Types of variants analyzed:**
- `stop_gained` → Truncated protein, possible loss of function
- `splice_donor_variant` / `splice_acceptor_variant` → Abnormal splicing
- `frameshift_variant` → Altered reading frame
- `missense_variant` → Amino acid change

**Example output:**

```
=== Level 3: Mechanism Validation ===
  SLC24A5: 4 stop_gained (truncated protein); 1 splice variants (abnormal splicing); 
           50 missense variants (amino acid changes). 
           Function: Calcium antiporter in melanosomes, critical for melanin synthesis...
  HERC2: 3 stop_gained (truncated protein); 2 splice variants (abnormal splicing); 
         20 missense variants (amino acid changes). 
         Function: E3 ubiquitin ligase, regulates OCA2 expression...
```

---

## Gene Database

The script includes an internal database with 14 known pigmentation genes:

### Primary Pigmentation Genes

| Gene | Function | OMIM | GWAS Reference |
|------|----------|------|----------------|
| **SLC24A5** | Ca²⁺ antiporter in melanosomes | 609802 | Lamason 2005, Crawford 2017 |
| **SLC45A2** | Melanosomal transport | 606202 | Graf 2005, Stokowski 2007 |
| **OCA2** | Modulates pH for tyrosinase activity | 611409 | Sturm 2008, Eiberg 2008 |
| **TYR** | Rate-limiting enzyme in melanin biosynthesis | 606933 | Sulem 2007, Han 2008 |
| **TYRP1** | Melanin biosynthesis stabilization | 115501 | Sulem 2007 |
| **MC1R** | Melanocortin 1 receptor | 155555 | Valverde 1995, Box 1997 |
| **MFSD12** | African pigmentation variant | 618170 | Crawford 2017 |

### Regulatory Genes

| Gene | Function | OMIM | GWAS Reference |
|------|----------|------|----------------|
| **HERC2** | Regulates OCA2 expression | 605837 | Sturm 2008, Eiberg 2008 |
| **ASIP** | MC1R antagonist | 600201 | Kanetsky 2002 |
| **IRF4** | Regulates TYR expression | 601900 | Han 2008, Praetorius 2013 |
| **BNC2** | Transcription factor in melanocytes | 608669 | Jacobs 2013 |

### UV Response Genes

| Gene | Function | OMIM | GWAS Reference |
|------|----------|------|----------------|
| **DDB1** | UV DNA damage repair | 600045 | Crawford 2017 |

### Other Related Genes

| Gene | Function | OMIM | GWAS Reference |
|------|----------|------|----------------|
| **EDAR** | Hair morphology (Asians) | 604095 | Fujimoto 2008 |
| **KITLG** | Melanocyte development | 184745 | Miller 2007, Sulem 2007 |

---

## Known Pigmentation SNPs

The script checks for the presence of these well-characterized SNPs:

| rsID | Gene | Population Effect | Consequence |
|------|------|-------------------|-------------|
| **rs1426654** | SLC24A5 | Light skin in Europeans (A allele: ~100% EUR, ~0% AFR) | missense (Ala111Thr) |
| **rs16891982** | SLC45A2 | Light skin/hair in Europeans | missense (Leu374Phe) |
| **rs12913832** | HERC2 | Eye color (A: brown, G: blue) | intron (regulatory) |
| **rs1800407** | OCA2 | Eye color modifier | missense (Arg419Gln) |
| **rs1042602** | TYR | Skin/hair color | missense (Ser192Tyr) |
| **rs1805007** | MC1R | Red hair, fair skin | missense (Arg151Cys) |
| **rs1805008** | MC1R | Red hair, fair skin | missense (Arg160Trp) |
| **rs11230664** | DDB1 | Darker skin in Africans | intron |
| **rs10424065** | MFSD12 | African-specific pigmentation | regulatory |
| **rs3827760** | EDAR | Hair thickness (Asians) | missense (Val370Ala) |
| **rs12203592** | IRF4 | Freckles, sun sensitivity | intron (regulatory) |

---

## Scoring System

### Score Calculation

The validation score (0-100) is calculated with weights:

| Component | Weight | Description |
|-----------|--------|-------------|
| **Gene match** | 40% | % of detected genes that are known pigmentation genes |
| **Known rsIDs** | 20% | Presence of known pigmentation rsIDs (binary) |
| **Crawford overlap** | 20% | % overlap with Crawford et al., 2017 genes |
| **Mechanism** | 20% | Based on presence of functional variants |

### Mechanism Criteria

| Condition | Points |
|-----------|--------|
| HIGH + MODERATE variants | 20 |
| Only HIGH or only MODERATE | 15 |
| Only missense | 10 |
| Only LOW/MODIFIER | 5 |

### Validation Status

| Score | Status |
|-------|--------|
| ≥ 80 | **STRONGLY VALIDATED (+++)** |
| 60-79 | **VALIDATED (++)** |
| 40-59 | **PARTIALLY VALIDATED (+)** |
| < 40 | **NOT VALIDATED (-)** |

---

## Output Report

The script generates a Markdown report (`pigmentation_validation.md`) with the following sections:

### 1. Executive Summary

```markdown
## Executive Summary

This report validates whether the DeepLIFT-identified genomic regions
correspond to known pigmentation genes...

### Validation Status: VALIDATED (++)
```

### 2. Validation Scores

Table with validation metrics:
- Total genes detected
- Primary pigmentation genes
- Regulatory genes
- Unknown genes

### 3. Gene Classifications

Detailed table of each detected gene:
- Category (primary, regulatory, uv_response)
- OMIM ID
- Function
- HIGH/MODERATE variant counts

### 4. Detailed Gene Analysis

For each gene:
- Category and function
- OMIM link
- Relevant GWAS studies
- Population effect
- Variant statistics

### 5. Comparison with Published GWAS

Comparison with Crawford et al., 2017 (Science):
- Which genes from the study were detected
- Overlap percentage

### 6. Biological Mechanism Analysis

Analysis of high-impact variants for each gene.

### 7. Conclusion

Final interpretation of validation and implications for longevity.

---

## gnomAD Integration

### Usage

```bash
python3 validate_pigmentation_hypothesis.py input_dir --gnomad
```

### What does gnomAD provide?

The gnomAD (Genome Aggregation Database) provides allele frequencies by population:

| Population | Code |
|------------|------|
| African | AFR |
| European | EUR |
| East Asian | EAS |
| American | AMR |
| South Asian | SAS |

### Interpretation

Typical pigmentation variants show very different frequencies between populations:

```
| rsID       | AFR    | EUR    | EAS    |
|------------|--------|--------|--------|
| rs1426654  | 0.0012 | 0.9987 | 0.0003 |  ← Strong differentiation
| rs16891982 | 0.0008 | 0.9523 | 0.0002 |  ← Strong differentiation
```

If detected variants show similar frequencies across populations, they are probably **not** functional pigmentation variants.

### Limitations

- gnomAD API can be slow
- Rate limiting may occur
- Only the first 10 rsIDs are queried
- Not all rsIDs have data in gnomAD

---

## Application to Longevity

### Why validate with pigmentation first?

Pigmentation is an "obvious" phenotype to validate the pipeline because:

1. **Well-known genes**: SLC24A5, OCA2, TYR are extensively studied
2. **Clear difference**: Africans have clearly different pigmentation from Europeans
3. **Available data**: 1000 Genomes has superpopulation for each individual
4. **Easy verification**: We can confirm if the pipeline identifies the correct genes

### Transition to longevity

Once the pipeline is validated with pigmentation, the same approach can be applied to longevity:

1. **Train** `neural_ancestry_predictor.py` with longevous vs non-longevous individuals
2. **Run** DeepLIFT to identify important genomic regions
3. **Annotate** with `annotate_deeplift_windows.py` to find variants
4. **Analyze** variants in known longevity genes (FOXO3, APOE, CETP, etc.)

### Known Longevity Genes

For future reference, genes associated with longevity include:

| Gene | Function | Reference |
|------|----------|-----------|
| FOXO3 | Transcription factor, stress response | Willcox 2008 |
| APOE | Lipid metabolism | Schachter 1994 |
| CETP | Cholesterol ester transfer | Barzilai 2003 |
| IGF1R | IGF-1 receptor | Suh 2008 |
| TERT | Telomerase | Atzmon 2010 |
| LMNA | Nuclear envelope | Eriksson 2003 |

---

## Examples

### Example 1: Basic Validation

```bash
$ python3 validate_pigmentation_hypothesis.py top_regions_reports_central2

[INFO] Reading data from: top_regions_reports_central2
[INFO] Detected genes: SLC24A5, HERC2, DDB1, OCA2, TYR
[INFO] Loaded 373 VEP-annotated variants

=== Level 1: Gene Validation ===
  Matched pigmentation genes: ['SLC24A5', 'HERC2', 'DDB1', 'OCA2', 'TYR']
  Unmatched genes: []
  Crawford et al. 2017 overlap: 4/5 (80.0%)

=== Level 2: Variant Validation ===
  Total rsIDs found: 67
  Known pigmentation rsIDs: None
  HIGH impact variants: 12
  MODERATE impact variants: 70
  Missense variants: 71

=== Level 3: Mechanism Validation ===
  SLC24A5: 4 stop_gained (truncated protein); 1 splice variants...
  HERC2: 3 stop_gained (truncated protein); 2 splice variants...
  DDB1: Regulatory/intronic variants only...
  OCA2: Regulatory/intronic variants only...
  TYR: Regulatory/intronic variants only...

=== Validation Result ===
  Status: VALIDATED (++)
  Score: 76.0/100
[SAVED] Validation report: top_regions_reports_central2/pigmentation_validation.md
```

### Example 2: With gnomAD

```bash
$ python3 validate_pigmentation_hypothesis.py top_regions_reports_central2 --gnomad

[INFO] Reading data from: top_regions_reports_central2
...
=== Fetching gnomAD Frequencies ===
  Fetched frequencies for 5 variants

=== Validation Result ===
  Status: VALIDATED (++)
  Score: 76.0/100
[SAVED] Validation report: top_regions_reports_central2/pigmentation_validation.md
```

---

## References

### Pigmentation Studies

- **Lamason RL et al. (2005)** - SLC24A5, a putative cation exchanger, affects pigmentation in zebrafish and humans. *Science*
- **Crawford NG et al. (2017)** - Loci associated with skin pigmentation identified in African populations. *Science*
- **Sturm RA et al. (2008)** - A single SNP in an evolutionary conserved region within intron 86 of the HERC2 gene determines human blue-brown eye color. *Am J Hum Genet*
- **Eiberg H et al. (2008)** - Blue eye color in humans may be caused by a perfectly associated founder mutation in a regulatory element. *Hum Genet*
- **Sulem P et al. (2007)** - Genetic determinants of hair, eye and skin pigmentation in Europeans. *Nat Genet*

### Databases

- **OMIM**: https://omim.org/
- **gnomAD**: https://gnomad.broadinstitute.org/
- **Ensembl VEP**: https://www.ensembl.org/vep
- **1000 Genomes**: https://www.internationalgenome.org/

### Longevity Studies

- **Willcox BJ et al. (2008)** - FOXO3A genotype is strongly associated with human longevity. *PNAS*
- **Schachter F et al. (1994)** - Genetic associations with human longevity at the APOE and ACE loci. *Nat Genet*
- **Barzilai N et al. (2003)** - Unique lipoprotein phenotype and genotype associated with exceptional longevity. *JAMA*

---

**Author**: Alberto F. De Souza  
**Date**: 2025-12-23  
**Version**: 1.0
