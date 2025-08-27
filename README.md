# Genomes Analyzer

_A technical-scientific guide to `genomes_analyzer.py`_

---

## Abstract

High-throughput DNA sequencing has transformed biological and clinical research, but turning raw reads into actionable variant calls still requires a reliable, explainable, and resource-efficient pipeline. **Genomes Analyzer** is a Python-driven workflow designed to take one or more whole-genome or exome samples from raw FASTQ (or pre-aligned BAM/CRAM) to compressed, indexed VCF. It emphasizes clear provenance, conservative defaults, and transparent logging while remaining pragmatic about compute and memory on commodity Linux workstations. The pipeline integrates widely used open-source tools—`FastQC` for read quality control, `cutadapt` for adapter/quality trimming, `bwa-mem2` for alignment, `samtools` for sorting/indexing/duplicate marking, and a selectable variant-calling backend: either **GATK** (HaplotypeCaller → GenotypeGVCFs) or **BCFtools** (`mpileup` → `call` → `norm`). To keep runtimes reasonable, Genomes Analyzer splits the reference into **shards** (BED intervals) and executes them in parallel with safe concatenation and indexing, while preserving chromosome order. It prints human-readable progress dashboards (including the shard BED preview and per-shard heartbeats) so researchers can monitor long runs without guesswork.

This document introduces the pipeline to readers **without assuming prior genomics expertise**. We define essential terms (e.g., FASTQ, BAM/CRAM, VCF, read mapping, duplicate marking, genotyping, normalization), explain each processing step and its rationale, and show representative snippets of intermediate and final outputs. We also document configuration via YAML, including paths, sample declarations, quality/coverage thresholds, and parallelization parameters. Finally, we provide an appendix detailing the external genomics tools used, typical command patterns, and commonly tuned parameters. The goal is to enable reproducible analyses with sensible defaults while making it straightforward to adapt the pipeline to different datasets, hardware budgets, and scientific objectives.

---

## Introduction

**Genomes Analyzer** is a modular, script-based pipeline for small teams or single investigators who need trustworthy variant calls from short-read sequencing. The script `genomes_analyzer.py` orchestrates a complete workflow from raw sequencing reads to an indexed Variant Call Format (VCF) file. The tool is particularly useful when:

- You want a **transparent** pipeline using standard tools with well-understood behavior.
- You need to **choose between** GATK and BCFtools callers without rewriting your workflow.
- You want robust **parallel sharding** of the reference genome to utilize multi-core CPUs.
- You value **clear logging**—what is running, on which intervals, and how progress looks over time.

### Key inputs and outputs

| Type | Description | Examples |
|------|-------------|----------|
| **Inputs** | Paired FASTQ files per sample (or a pre-aligned BAM/CRAM) and a reference FASTA with its index and dictionary. | `fastq/NA12878_R1.fastq.gz`, `fastq/NA12878_R2.fastq.gz`, `refs/GRCh38.d1.vd1.fa` |
| **Primary outputs** | Aligned BAM/BAI, per-shard VCFs, final VCF/VCF index, and quality-control reports. | `bam/NA12878.mkdup.bam`, `vcf/NA12878.vcf.gz`, `qc/NA12878_fastqc.html` |

---

## Genomes Analyzer Pipeline

The workflow transforms raw reads into variant calls through a series of well-defined stages. Each acronym is introduced before use so that readers new to genomics can follow along.

```
FASTQ
  ├── Quality control (FastQC)
  ├── Adapter trimming (cutadapt)
  ├── Alignment (BWA‑MEM2)
  ├── Sort & index (samtools)
  ├── Mark duplicates (samtools markdup)
  ├── [Optional] Base Quality Score Recalibration (GATK)
  ├── Shard reference genome (BED intervals)
  ├── Variant calling
  │     ├── GATK: HaplotypeCaller → GenotypeGVCFs
  │     └── BCFtools: mpileup → call → norm
  ├── Concatenate shards (preserve order)
  └── Final VCF quality control
```

### 1. Read Quality Control — FastQC

**FASTQ** files store sequencing reads and their per-base quality scores. `FastQC` scans these files and produces HTML reports summarizing quality metrics such as Phred scores, nucleotide composition, and overrepresented sequences.

*Why it matters*: Early detection of poor-quality cycles or adapter contamination prevents misleading alignments and spares compute time.

*Representative snippet* (`fastqc_data.txt`):

```
>>Per base sequence quality  pass
#Base    Mean    Median  Lower   Upper
1        33.8    34      31      36
...
>>Overrepresented sequences  warn
```

*Outputs*: `qc/<sample>_R1_fastqc.html`, `qc/<sample>_R2_fastqc.html` plus compressed archives containing raw metrics.

### 2. Adapter and Quality Trimming — cutadapt

Sequencing libraries often carry leftover **adapter** sequences and low-quality ends. `cutadapt` removes these artifacts and can filter short reads.

*Why it matters*: Adapter sequences and low-quality bases reduce mapping accuracy and inflate false variant calls.

*Representative snippet* (log extract):

```
=== Summary ===
Total read pairs processed:     374,102,311
Pairs written (passing filters):372,918,421 (99.7%)
Total basepairs processed:      112.2 Gbp
Quality-trimmed:                1.6 Gbp (1.4%)
```

*Outputs*: `trimmed/<sample>_R1.fastq.gz`, `trimmed/<sample>_R2.fastq.gz`.

### 3. Alignment — BWA‑MEM2

`bwa-mem2` maps each read pair to the reference genome. The output is a **SAM** (Sequence Alignment/Map) stream that records candidate genomic coordinates, alignment scores, and flags. Alignments are typically converted on the fly to the binary **BAM** format.

*Why it matters*: Accurate mapping is a prerequisite to reliable variant detection. Each read receives a mapping quality score indicating confidence in its genomic position.

*Representative SAM header*:

```
@SQ SN:chr1 LN:248956422
@SQ SN:chr2 LN:242193529
@PG ID:bwa-mem2 PN:bwa-mem2 VN:2.2.1 CL:bwa-mem2 mem -t 16 ...
```

*Outputs*: `aligned/<sample>.sam` (usually streamed downstream).

### 4. Sorting & Indexing — samtools

`samtools sort` arranges BAM records by genomic coordinate and writes an index (`.bai`) to enable random access.

*Why it matters*: Variant callers expect coordinate-sorted and indexed BAM files to quickly fetch reads overlapping regions of interest.

*Outputs*: `bam/<sample>.sorted.bam`, `bam/<sample>.sorted.bam.bai`.

### 5. Duplicate Marking — samtools markdup

PCR amplification and optical artifacts can yield **duplicate reads**—multiple observations of the same DNA fragment. `samtools markdup` flags these duplicates so that callers can ignore them.

*Why it matters*: Treating duplicates as independent evidence biases allele counts and may cause false positives.

*Outputs*: `bam/<sample>.mkdup.bam`, `bam/<sample>.mkdup.bam.bai`.

### 6. Optional Base Quality Score Recalibration — GATK BQSR

Sequencing machines sometimes misestimate base quality scores. **Base Quality Score Recalibration (BQSR)** uses known variant sites to adjust these scores.

*Why it matters*: More accurate base qualities improve probabilistic models in downstream callers. This step is optional and can be skipped to save time.

*Outputs*: recalibrated BAM (if applied) plus BQSR reports.

### 7. Sharding the Reference Genome

Large genomes are divided into manageable **BED** intervals ("shards"). Each shard is processed independently to leverage parallelism.

*Why it matters*: Short-read variant callers parallelize poorly within a single region; sharding provides coarse-grained parallelism across chromosomes/blocks and reduces wall-clock time.

*Outputs*: `vcf/shards/<sample>/<sample>.part_XX.vcf.gz` and corresponding `.tbi` indexes.

### 8. Variant Calling

At this stage, aligned reads are converted into genomic variants. Two backends are available:

#### A. GATK HaplotypeCaller → GenotypeGVCFs

1. **HaplotypeCaller** analyzes each shard, performing local re-assembly to produce per-sample genomic VCFs (**gVCFs**).
2. **GenotypeGVCFs** merges gVCFs and emits a final VCF with genotypes.

*Outputs*: `vcf/<sample>.g.vcf.gz`, `vcf/<sample>.vcf.gz` and their indexes.

#### B. BCFtools mpileup → call → norm

1. **mpileup** summarizes read evidence at each position, applying mapping-quality (MAPQ) and base-quality (BASEQ) filters.
2. **call** infers SNPs and indels using a multiallelic model.
3. **norm** left-aligns indels and splits multiallelic records for compatibility.

*Outputs*: sharded VCFs (`part_XX.vcf.gz`) and final merged VCF (`vcf/<sample>.vcf.gz`).

### 9. Concatenation of Shards

Per-shard VCFs are concatenated in chromosome order and re-indexed so that downstream tools see a seamless, genome-wide VCF.

*Outputs*: `vcf/<sample>.vcf.gz`, `vcf/<sample>.vcf.gz.tbi`.

### 10. Optional VCF Quality Control

`bcftools stats` and related utilities generate summary metrics and plots to assess callset quality.

*Outputs*: `qc/<sample>.bcftools.stats.txt` and optional graphical summaries via `plot-vcfstats` or `MultiQC`.

### Output overview

| Path | Type | Description |
|------|------|-------------|
| `bam/<sample>.mkdup.bam` | BAM | Coordinate-sorted, duplicate-marked alignments |
| `vcf/shards/<sample>/part_XX.vcf.gz` | VCF | Per-shard variant calls |
| `vcf/<sample>.vcf.gz` | VCF | Final, genome-wide variant calls |
| `qc/<sample>_R1_fastqc.html` | HTML | Read quality report |

---

## How to Use the Genomes Analyzer

### Installation

1. Install a Python environment (e.g., via [conda](https://docs.conda.io/)).
2. Optional helper script: `install_genomics_env.sh` installs common genomics packages.
3. Activate the environment: `source start_genomics.sh`.

### YAML Configuration: `config_human_30x_low_memory.yaml`

`genomes_analyzer.py` is driven by a YAML file. The example `config_human_30x_low_memory.yaml` targets a human trio at ~30× coverage using memory‑efficient settings. Important sections are summarized below.

#### project

| Field | Description | Example |
|-------|-------------|---------|
| `name` | Project identifier used in output paths. | `human_30x_trio_demo` |
| `organism` | Latin name of the species. | `homo_sapiens` |
| `reference.name` | Identifier for the reference genome build. | `GRCh38.d1.vd1` |
| `reference.fasta_url` | URL to a gzipped FASTA that will be downloaded and unpacked. | `https://api.gdc.cancer.gov/...834` |
| `reference.bwa_index_url` | Pre-built BWA index (saves RAM during alignment). | `https://api.gdc.cancer.gov/...7225` |

#### general

| Field | Purpose | Typical value |
|-------|---------|---------------|
| `force_indexes` | Rebuild alignment indexes even if present. | `false` |
| `sort_mem_mb` | Memory (MiB) per thread for `samtools sort`. | `512` |
| `bwa_batch_k` | Reads per batch for BWA; smaller uses less RAM. | `20000000` |
| `aln_threads` | Threads for alignment. | `16` |
| `gene_presence_min_mean_cov` | Minimum average coverage to consider a gene present. | `5.0` |

#### params

Alignment and variant-calling options.

| Field | Description | Example |
|-------|-------------|---------|
| `aligner` | Choose `bwa` or `bwa-mem2`. | `bwa` |
| `variant_caller` | `gatk` or `bcftools`. | `bcftools` |
| `bcf_mapq` / `bcf_baseq` | Minimum mapping/base quality in `mpileup`. | `20` |
| `bcf_scatter_parts` | Number of BED shards for parallel calling. | `16` |
| `hc_java_mem_gb` | Heap size for GATK HaplotypeCaller. | `24` |

#### storage

| Field | Description |
|-------|-------------|
| `base_dir` | Root directory for results. |
| `work_dir` | Working directory for temporary files. |
| `temp_dir` | High‑speed disk for intermediates. |

#### download

Defines how sequencing data are retrieved from the Sequence Read Archive (SRA) or European Nucleotide Archive (ENA).

| Field | Description | Example |
|-------|-------------|---------|
| `tool` | Download utility (`sra_toolkit`). | `sra_toolkit` |
| `use_ascp` | Use Aspera for faster transfers if available. | `false` |
| `threads` | Threads for download/conversion. | `8` |

#### execution

Run‑time behaviour and resilience.

| Field | Description | Example |
|-------|-------------|---------|
| `verbose` | Verbose logging. | `true` |
| `resume` | Skip completed steps. | `true` |
| `progress_interval_sec` | Interval for progress reports. | `60` |

#### size_control

Controls downsampling or disk usage.

| Field | Description | Example |
|-------|-------------|---------|
| `downsample.enabled` | Whether to subsample reads. | `true` |
| `downsample.fraction` | Fraction of reads to retain (e.g., `0.25` ≈7.5×). | `0.25` |

#### samples

Defines the biological samples to process.

| Field | Description |
|-------|-------------|
| `sample_id` | Unique identifier (e.g., `NA12878`). |
| `study` | Accession of the sequencing study. |
| `runs` | List of SRA/ENA run IDs containing the reads. |

#### steps

Ordered list of pipeline actions. Typical values include `fetch_fastqs`, `qc_reads`, `align_and_sort`, `mark_duplicates`, `bqsr`, `call_genes`, and `summarize`.

### Running the pipeline

Execute the workflow by pointing the script to your YAML file:

```bash
source start_genomics.sh
./genomes_analyzer.py --config config_human_30x_low_memory.yaml
```

The console prints progress panels, including per-shard heartbeats when variant calling is parallelized.

---

## Conclusion

Genomes Analyzer provides a clear, reproducible path from raw reads to variant calls using a dependable stack of open‑source tools. By letting users switch between GATK and BCFtools and by exposing explicit sharding and parallelization controls, it adapts to many datasets and machines. The detailed documentation and YAML configuration aim to demystify genomics processing for newcomers while remaining efficient for experienced practitioners.

---

## Appendix 1 — Tools & typical usage

The pipeline wraps several established command‑line programs. Table 1 summarizes their roles and highlights common parameters.

| Tool | Role | Typical command | Key parameters |
|------|------|----------------|----------------|
| FastQC | Read quality assessment | `fastqc -t 8 fastq/*_R{1,2}.fastq.gz -o qc/` | `-t` threads |
| cutadapt | Adapter and quality trimming | `cutadapt -j 8 -q 20,20 -m 30 -a ADAPT1 -A ADAPT2 -o out_R1.fastq.gz -p out_R2.fastq.gz in_R1.fastq.gz in_R2.fastq.gz` | `-q` quality cutoff, `-m` min length, `-a/-A` adapters |
| BWA‑MEM2 | Short‑read alignment | `bwa-mem2 mem -t 16 -R "@RG\tID:S\tSM:S" refs.fa R1.fq.gz R2.fq.gz` | `-t` threads, `-R` read group |
| samtools | BAM/CRAM processing | `samtools sort -@8 -o out.bam in.bam` | `-@` threads |
| BCFtools | Variant calling and manipulation | `bcftools mpileup -f refs.fa -q20 -Q20 | bcftools call -m -v -Oz -o out.vcf.gz` | `-q/-Q` quality filters, `-m` multiallelic model |
| GATK | Variant calling (Java) | `gatk --java-options "-Xmx24g" HaplotypeCaller -R refs.fa -I in.bam -O out.g.vcf.gz` | `--java-options` memory, `-L` intervals |

Each tool offers many more options; consult the official manuals for advanced usage. The examples above match how `genomes_analyzer.py` invokes them under default settings.

