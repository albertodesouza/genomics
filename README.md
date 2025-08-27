# Genomes Analyzer

_A technical-scientific guide to `genomes_analyzer.py`_

---

## Abstract

High-throughput DNA sequencing has transformed biological and clinical research, but turning raw reads into actionable variant calls still requires a reliable, explainable, and resource-efficient pipeline. **Genomes Analyzer** is a Python-driven workflow designed to take one or more whole-genome or exome samples from raw FASTQ (or pre-aligned BAM/CRAM) to compressed, indexed VCF. It emphasizes clear provenance, conservative defaults, and transparent logging, while remaining pragmatic about compute and memory on commodity Linux workstations. The pipeline integrates widely used open-source tools—`FastQC` for read quality control, `cutadapt` for adapter/quality trimming, `bwa-mem2` for alignment, `samtools` for sorting/indexing/duplication marking, and a selectable variant-calling backend: either **GATK** (HaplotypeCaller → GenotypeGVCFs) or **BCFtools** (`mpileup` → `call` → `norm`). To keep runtimes reasonable, Genomes Analyzer splits the reference into **shards** (BED intervals) and executes them in parallel with safe concatenation and indexing, while preserving chromosome order. It prints human-readable “progress dashboards” (including the shard BED preview and per-shard heartbeats) so researchers can monitor long runs without guesswork.

This document introduces the pipeline to readers **without assuming prior genomics expertise**. We define essential terms (e.g., FASTQ, BAM/CRAM, VCF, read mapping, duplicate marking, genotyping, normalization), explain each processing step and its rationale, and show representative snippets of intermediate and final outputs. We also document configuration via YAML, including paths, sample declarations, quality/coverage thresholds, and parallelization parameters. Finally, we provide an appendix detailing the external genomics tools used, typical command patterns, and commonly tuned parameters. The goal is to enable reproducible analyses with sensible defaults, while making it straightforward to adapt the pipeline to different datasets, hardware budgets, and scientific objectives.

---

## Introduction

**Genomes Analyzer** is a modular, script-based pipeline for small teams or single investigators who need trustworthy variant calls from short-read sequencing. It is particularly useful when:

- You want a **transparent** pipeline using standard tools with well-understood behavior.
- You need to **choose between** GATK and BCFtools callers without rewriting your workflow.
- You want robust **parallel sharding** of the reference genome to utilize multi-core CPUs.
- You value **clear logging**—what is running, on which intervals, and how progress looks over time.

### Key inputs and outputs

- **Inputs**: paired FASTQ files per sample (or a pre-aligned BAM/CRAM), a reference FASTA (`refs/reference.fa`) with index (`.fai`) and dictionary (`.dict`).
- **Primary outputs**:
  - Per-sample aligned, duplicate-marked **BAM** (`bam/<sample>.mkdup.bam`) + index (`.bai`).
  - Per-sample compressed, indexed **VCF** (`vcf/<sample>.vcf.gz`) + index (`.tbi`).
  - Per-sample **sharded VCFs** for parallel calling (`vcf/shards/<sample>/<sample>.part_XX.vcf.gz` + `.tbi`).
  - Quality control artifacts (FastQC reports, optional `bcftools stats`, etc.).

---

## Genomes Analyzer Pipeline

Below is the end-to-end flow. Terms are defined as they appear.

```
FASTQ → (QC) → (Trimming) → Alignment → Sort/Index → Mark Duplicates
   → [Optional: BQSR] → Shard Reference → Variant Calling (GATK or BCFtools)
   → Normalize/Index → Concatenate Shards → Final VCF (+ QC)
```

### 1) Read QC (FastQC)

**What it is**: Rapid quality assessment of raw reads (per-cycle Phred scores, over-represented sequences, adapter content).

**Why it matters**: Detects library/adaptor issues that can propagate to false positives later.

**Representative snippet** (from `fastqc_data.txt`):
```
>>Per base sequence quality  pass
#Base    Mean    Median  Lower   Upper
1        33.8    34      31      36
...
>>Overrepresented sequences  warn
```

**Outputs**: `qc/<sample>_R1_fastqc.html`, `qc/<sample>_R2_fastqc.html` (+ `.zip`).

---

### 2) Adapter/quality trimming (cutadapt)

**What it is**: Removes leftover adapters and optionally low-quality trailing bases.

**Why it matters**: Adapters reduce mapping accuracy and inflate false variant calls.

**Representative snippet** (from `cutadapt` log):
```
=== Summary ===
Total read pairs processed:    374,102,311
Pairs written (passing filters): 372,918,421 (99.7%)
Total basepairs processed:     112.2 Gbp
Quality-trimmed:               1.6 Gbp (1.4%)
```

**Outputs**: `trimmed/<sample>_R1.fastq.gz`, `trimmed/<sample>_R2.fastq.gz`.

---

### 3) Alignment (bwa-mem2)

**What it is**: Aligns reads to the reference genome, recording positions where reads best match.

**Why it matters**: Accurate mapping is prerequisite to reliable variant calling.

**Representative SAM header**:
```
@SQ SN:chr1 LN:248956422
@SQ SN:chr2 LN:242193529
@PG ID:bwa-mem2 PN:bwa-mem2 VN:2.2.1 CL:bwa-mem2 mem -t 16 ...
```

**Outputs**: `aligned/<sample>.sam` (intermediate; often piped onward).

---

### 4) Sorting & indexing (samtools)

**What it is**: Sorts alignments by genomic coordinate; builds BAM index for random access.

**Why it matters**: Many downstream tools assume coordinate-sorted, indexed BAM.

**Outputs**: `bam/<sample>.sorted.bam`, `bam/<sample>.sorted.bam.bai`.

---

### 5) Duplicate marking (samtools markdup)

**What it is**: Marks PCR/optical duplicates (reads that are copies of the same original fragment).

**Why it matters**: Duplicates inflate apparent depth and can cause false positives.

**Outputs**: `bam/<sample>.mkdup.bam`, `bam/<sample>.mkdup.bam.bai`.

---

### 6) [Optional] Base Quality Score Recalibration (BQSR)

**What it is**: Adjusts base quality scores using known variant sites.

**When to use**: If you have reliable known-sites for your organism and sequencing chemistry benefits from BQSR. It is off by default in many lightweight pipelines.

---

### 7) Sharding the reference (BED intervals)

**What it is**: Splits the reference genome into **balanced BED intervals** (“shards”), often by chromosome or chromosome groups, to parallelize variant calling.

**Why it matters**: Shortens wall-clock time on multi-core machines.

**Illustrative preview table** (printed before calling):

```
┏━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━━━━━┓
┃ part ┃ contigs (preview)                                                              ┃ #regions ┃      total ┃ BED         ┃
┡━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━━━━━┩
│   01 │ chr1                                                                           │        1 │ 248.96 Mbp │ part_01.bed │
│   02 │ chr2                                                                           │        1 │ 242.19 Mbp │ part_02.bed │
│   …  │ …                                                                              │      …   │        …   │ part_NN.bed │
└──────┴────────────────────────────────────────────────────────────────────────────────┴──────────┴────────────┴─────────────┘
```

**Outputs**: `vcf/shards/<sample>/part_XX.bed`.

---

### 8) Variant calling

Two interchangeable backends:

#### A) GATK HaplotypeCaller → GenotypeGVCFs (Java)

- **HaplotypeCaller** produces a per-sample **gVCF** that encodes evidence at variant and non-variant sites.
- **GenotypeGVCFs** converts gVCF to a standard per-sample **VCF**.
- Supports **scatter/gather** across shards.

Representative HaplotypeCaller command (simplified):
```
gatk --java-options "-Xmx24g" HaplotypeCaller \
  -R refs/reference.fa -I bam/<sample>.mkdup.bam \
  -O vcf/<sample>.g.vcf.gz \
  [-L vcf/shards/<sample>/part_01.bed]
```

#### B) BCFtools mpileup → call → norm (C)

- `mpileup`: summarizes per-position evidence.
- `call -m`: multiallelic caller producing SNPs/indels.
- `norm`: left-aligns and splits multiallelics for consistency.

Representative per-shard pipeline:
```
bcftools mpileup -f refs/reference.fa -Ou --threads 2 \
  -q 20 -Q 20 --max-depth 250 \
  -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
  -R vcf/shards/<sample>/part_03.bed bam/<sample>.mkdup.bam |
bcftools call -m -v -Ou |
bcftools norm -f refs/reference.fa --threads 2 --multiallelics -both \
  -Oz -o vcf/shards/<sample>/<sample>.part_03.tmp.vcf.gz
```

**Per-shard heartbeat** (printed during runs):
```
[<sample>] shard 03… 2m30s • out 3.6MB • I/O +44.4MB read, +608KB write
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 250
Lines total/split/joined/...: 307319/2973/0/45665/0/0/0
```

**Outputs**:
- GATK path: `vcf/<sample>.g.vcf.gz`, then `vcf/<sample>.vcf.gz` (+ `.tbi`).
- BCFtools path: `vcf/shards/<sample>/<sample>.part_XX.vcf.gz` (+ `.tbi`).

---

### 9) Concatenation of shards (safe, ordered)

**What it is**: Joins shard VCFs into a single, coordinate-ordered per-sample VCF.

**Robustness features**:
- Deterministic ordering by shard index.
- Uses `bcftools concat` with `-a` (allow non-contiguity) and **dedup strategy `-D none`** to avoid creating spurious duplicate rows while still accepting non-contiguous blocks (e.g., when chromosomes are grouped in shards).
- Final index with `tabix`.

Representative command (generated by the script):
```
bcftools concat -a -D none -Oz -o vcf/<sample>.tmp.vcf.gz \
  vcf/shards/<sample>/<sample>.part_01.vcf.gz \
  ... \
  vcf/shards/<sample>/<sample>.part_16.vcf.gz
tabix -p vcf vcf/<sample>.tmp.vcf.gz
```

**Outputs**: `vcf/<sample>.vcf.gz`, `vcf/<sample>.vcf.gz.tbi`.

---

### 10) Optional VCF QC

- `bcftools stats vcf/<sample>.vcf.gz > qc/<sample>.bcftools.stats.txt`
- `plot-vcfstats -p qc/<sample>_plots qc/<sample>.bcftools.stats.txt`

**Representative VCF record**:
```
#CHROM POS      ID REF ALT QUAL FILTER INFO                 FORMAT  <sample>
chr21  33033448 .  G   A   162  PASS   DP=34;MQ=60;QD=18.0  GT:AD:DP:GQ:PL  0/1:16,18:34:99:...
```

---

## How to Use the Genomes Analyzer

### Installation

**Requirements (Linux/macOS recommended)**

- Python ≥ 3.10
- Conda or mamba (recommended)
- ~16 CPU cores & 32–64 GB RAM for whole-genome speed (but runs on less)
- Disk space: input FASTQs + ~2–3× for intermediates

**Create an environment**

```bash
# with mamba (faster) or conda
mamba create -n genomes python=3.11 fastqc cutadapt bwa-mem2 samtools \
  bcftools htslib tabix -c bioconda -c conda-forge
# GATK (optional; requires Java):
mamba install -n genomes gatk4 openjdk -c bioconda -c conda-forge

mamba activate genomes
```

**Project layout (example)**
```
.
├── genomes_analyzer.py
├── config.yaml
├── refs/
│   ├── reference.fa
│   ├── reference.fa.fai
│   └── reference.dict
└── fastq/
    ├── NA12878_R1.fastq.gz
    └── NA12878_R2.fastq.gz
```

### Configuration (YAML)

Below is a representative `config.yaml` with comments explaining each field.

```yaml
general:
  threads: 16             # total logical cores to utilize
  mem_gb: 64              # available RAM (soft hint for tools that accept it)
  limit_to_canonical: false   # if true, restricts to canonical chromosomes only (no decoy/ALT)
  work_dir: "."           # base working directory (optional)

paths:
  reference: refs/reference.fa
  known_sites: []         # optional VCF(s) for BQSR or hard filters
  adapters: []            # optional adapter sequences for cutadapt
  fastq_dir: fastq        # where FASTQs live (if using FASTQ input)

params:
  # Choose the variant caller backend
  variant_caller: "bcftools"          # "bcftools" or "gatk"

  # BCFtools path (mpileup → call → norm)
  bcf_scatter_parts: 16               # how many BED shards to generate
  bcf_max_parallel: 8                 # max shards to run simultaneously
  bcf_threads_io: 2                   # htslib I/O threads per shard
  bcf_mapq: 20                        # minimum mapping quality (mpileup -q)
  bcf_baseq: 20                       # minimum base quality (mpileup -Q)
  bcf_max_depth: 250                  # cap per-position read depth (mpileup --max-depth)
  keep_alt_decoy: true                # include ALT/decoy contigs in shards
  prefer_bam_for_calling: true        # prefer BAM over CRAM when both exist
  force_recall: false                 # ignore cache, recompute shards/VCFs

  # GATK path (HaplotypeCaller → GenotypeGVCFs)
  hc_scatter_parts: 16                # number of shards for HaplotypeCaller (if used)
  hc_java_mem_gb: 24                  # HaplotypeCaller Java heap
  hc_threads_native: 4                # native threads for HaplotypeCaller
  hc_pcr_free: true                   # adjust model for PCR-free libraries (if applicable)
  hc_extra_args: []                   # additional raw args, e.g. ["--dont-use-soft-clipped-bases"]

  gg_java_mem_gb: 12                  # GenotypeGVCFs Java heap
  gg_extra_args: []                   # additional raw args for GenotypeGVCFs

samples:
  - id: "NA12878"
    # Option A: from FASTQ (paired-end)
    r1: fastq/NA12878_R1.fastq.gz
    r2: fastq/NA12878_R2.fastq.gz
    # Option B: from pre-aligned data (comment FASTQs and provide BAM/CRAM)
    # bam: bam/NA12878.mkdup.bam

  - id: "NA12891"
    r1: fastq/NA12891_R1.fastq.gz
    r2: fastq/NA12891_R2.fastq.gz
```

> **Tip (parallelization)**: On a 16-core workstation, starting points that work well are `bcf_scatter_parts: 16`, `bcf_max_parallel: 8–12`, `bcf_threads_io: 1–2`. If you see only one core busy, increase `bcf_max_parallel` until you saturate CPU without thrashing disk I/O.

### Running the pipeline

```bash
python genomes_analyzer.py --config config.yaml
```

**Common flags** (if supported by your copy of the script; optional):
- `--dry-run` : show what would be done without executing.
- `--resume`  : reuse cached outputs when possible (default behavior).
- `--force`   : ignore caches.

### What to expect on the console

- A **shard BED preview** table describing how the reference is split.
- Per-shard **heartbeat** lines showing runtime, bytes read/written, and the growing `.tmp.vcf.gz`.
- Clear **cache notices** (e.g., “shard pronto (cache) → …”) to avoid redundant work.
- Safe concatenation summary with the exact `bcftools concat` command.

---

## Conclusion

Genomes Analyzer provides a clear, reproducible path from raw reads to variant calls using a small, dependable stack of open-source tools. By letting you switch between GATK and BCFtools and by exposing explicit sharding and parallelization controls, it adapts to many datasets and machines. Its console output favors transparency—what is running, on which genome segments, and with which performance—so you can trust both the process and the product. This guide should enable you to install, configure, and operate the pipeline confidently, even if you are new to genomics workflows.

---

## Appendix 1 — Tools & typical usage

> Commands are illustrative and omit some options for brevity. Consult each tool’s manual for full details.

### FastQC

- **Purpose**: Per-sample read quality report.
- **Run**:
  ```bash
  fastqc -t 8 fastq/*_R{1,2}.fastq.gz -o qc/
  ```

### cutadapt

- **Purpose**: Adapter & quality trimming.
- **Common options**:
  - `-a/-A`: adapter sequences (R1/R2).
  - `-q 20,20`: trim low-quality ends (Phred ≥20).
  - `-m 30`: minimum length after trimming.
- **Run**:
  ```bash
  cutadapt -j 8 -q 20,20 -m 30 \
    -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -o trimmed/S_R1.fastq.gz -p trimmed/S_R2.fastq.gz \
    fastq/S_R1.fastq.gz fastq/S_R2.fastq.gz
  ```

### BWA-MEM2

- **Purpose**: Fast alignment of short paired-end reads.
- **Common options**:
  - `-t`: threads
  - `-R`: read group (ID, sample, platform)
- **Run**:
  ```bash
  bwa-mem2 index refs/reference.fa
  bwa-mem2 mem -t 16 -R "@RG\tID:S\tSM:S\tPL:ILLUMINA" \
    refs/reference.fa trimmed/S_R1.fastq.gz trimmed/S_R2.fastq.gz \
    | samtools view -b -o bam/S.unsorted.bam
  ```

### samtools

- **Purpose**: BAM/CRAM processing.
- **Sort & index**:
  ```bash
  samtools sort -@ 8 -o bam/S.sorted.bam bam/S.unsorted.bam
  samtools index bam/S.sorted.bam
  ```
- **Mark duplicates**:
  ```bash
  samtools markdup -@ 8 -s bam/S.sorted.bam bam/S.mkdup.bam
  samtools index bam/S.mkdup.bam
  ```
- **Flagstat**:
  ```bash
  samtools flagstat -@ 8 bam/S.mkdup.bam > qc/S.flagstat.txt
  ```

### BCFtools (C backend)

- **mpileup** (summary of evidence per position)
  - `-q`: min mapping quality
  - `-Q`: min base quality
  - `--max-depth`: cap for per-position reads (performance/FP control)
  - `-a`: which annotations to include (e.g., `FORMAT/AD,FORMAT/DP`)
- **call -m**: multiallelic model (SNPs/indels)
- **norm**: left-align, split multiallelics
- **concat**: concatenate shards
- **stats**: VCF QC summary
- **Examples**:
  ```bash
  bcftools mpileup -f refs/reference.fa -q 20 -Q 20 --max-depth 250 \
    -a FORMAT/AD,FORMAT/DP -Ou -R part_03.bed bam/S.mkdup.bam |
  bcftools call -m -v -Ou |
  bcftools norm -f refs/reference.fa --multiallelics -both -Oz -o S.part_03.vcf.gz
  tabix -p vcf S.part_03.vcf.gz

  bcftools concat -a -D none -Oz -o S.vcf.gz S.part_*.vcf.gz
  tabix -p vcf S.vcf.gz

  bcftools stats S.vcf.gz > qc/S.bcftools.stats.txt
  ```

### GATK (Java backend)

- **HaplotypeCaller** (gVCF mode):
  ```bash
  gatk --java-options "-Xmx24g -Dsamjdk.compression_level=2" HaplotypeCaller \
    -R refs/reference.fa -I bam/S.mkdup.bam \
    -O vcf/S.g.vcf.gz \
    [-L vcf/shards/S/part_01.bed] --native-pair-hmm-threads 4
  ```
- **GenotypeGVCFs** (gVCF → VCF):
  ```bash
  gatk --java-options "-Xmx12g -Dsamjdk.compression_level=2" GenotypeGVCFs \
    -R refs/reference.fa -V vcf/S.g.vcf.gz -O vcf/S.vcf.gz
  ```

---

### Frequently Asked Questions

**Q: Which backend should I choose?**  
A: If your lab standardizes on GATK, use it for consistency. If you prefer a fast, C-based stack with simpler dependencies, the BCFtools path is robust and efficient. The pipeline lets you switch via `params.variant_caller`.

**Q: Why split into shards?**  
A: Short-read callers parallelize poorly within a single region; sharding gives coarse-grained parallelism across chromosomes/blocks, improving wall-clock time.

**Q: I see only one core busy.**  
A: Increase `bcf_max_parallel` (number of shards running concurrently). Keep `bcf_threads_io` small (1–2) to avoid I/O bottlenecks.

---
