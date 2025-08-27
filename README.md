# Genomes Analyzer

_A technical-scientific guide to `genomes_analyzer.py`_

## Index
- [Abstract](#abstract)
- [Introduction](#introduction)
- [Genomes Analyzer Pipeline](#genomes-analyzer-pipeline)
- [How to Use the Genomes Analyzer](#how-to-use-the-genomes-analyzer)
- [Conclusion](#conclusion)
- [Appendix 1 — Tools & typical usage](#appendix-1--tools--typical-usage)

---

## Abstract

High-throughput DNA sequencing has transformed biological and clinical research, but turning raw reads into actionable variant calls still requires a reliable, explainable, and resource-efficient pipeline. **Genomes Analyzer** is a Python-driven workflow designed to take one or more whole-genome or exome samples from raw [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) (or pre-aligned [BAM/CRAM](https://en.wikipedia.org/wiki/SAM_(file_format))) to compressed, indexed [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf). It emphasizes clear provenance, conservative defaults, and transparent logging while remaining pragmatic about compute and memory on commodity Linux workstations. The pipeline integrates widely used open-source tools—[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for read quality control, [cutadapt](https://cutadapt.readthedocs.io/) for adapter/quality trimming, [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) for alignment, [samtools](http://www.htslib.org/) for sorting/indexing/duplicate marking, and a selectable variant-calling backend: either [GATK](https://gatk.broadinstitute.org/hc/en-us) (HaplotypeCaller → GenotypeGVCFs) or [BCFtools](https://samtools.github.io/bcftools/) (`mpileup` → `call` → `norm`). To keep runtimes reasonable, Genomes Analyzer splits the reference into **shards** ([BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) intervals) and executes them in parallel with safe concatenation and indexing, while preserving chromosome order. It prints human-readable progress dashboards (including the shard BED preview and per-shard heartbeats) so researchers can monitor long runs without guesswork.

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

Before running the workflow, make sure the environment is prepared as described in [How to Use the Genomes Analyzer](#how-to-use-the-genomes-analyzer). The workflow transforms raw [sequencing reads](https://en.wikipedia.org/wiki/Read_(biology))—unaltered sequences produced by high-throughput instruments—into *variant calls*: structured records that capture single‑nucleotide changes, small insertions/deletions, and other deviations from a reference genome, typically stored in the text-based [Variant Call Format (VCF)](https://en.wikipedia.org/wiki/Variant_Call_Format). These stages are well-defined, and each acronym is introduced before use so that readers new to genomics can follow along.

```
FASTQ
  ├── Quality control (FastQC)
  ├── Adapter trimming (cutadapt)
  ├── Alignment (BWA‑MEM2)
  ├── Sort & index (samtools)
  ├── Mark duplicates (samtools markdup)
  ├── [Optional] Base Quality Score Recalibration (GATK)
  ├── CRAM conversion & coverage (mosdepth)
  ├── Shard reference genome (BED intervals)
  ├── Variant calling
  │     ├── GATK: HaplotypeCaller → GenotypeGVCFs
  │     └── BCFtools: mpileup → call → norm
  ├── Concatenate shards (preserve order)
  ├── Final VCF quality control
  ├── Gene list & per-sample coverage
  ├── Pairwise & trio comparisons
  └── [Optional] RNA-seq module
```

### 1. Read Quality Control — [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

**FASTQ** files store sequencing reads and their per-base quality scores. `FastQC` scans these files and produces HTML reports summarizing quality metrics such as [Phred scores](https://en.wikipedia.org/wiki/Phred_quality_score), nucleotide composition, and overrepresented sequences.

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

### 2. Adapter and Quality Trimming — [cutadapt](https://cutadapt.readthedocs.io/)

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

### 3. Alignment — [BWA‑MEM2](https://github.com/bwa-mem2/bwa-mem2)

`bwa-mem2` maps each read pair to the reference genome. The output is a **[SAM](https://en.wikipedia.org/wiki/SAM_(file_format))** (Sequence Alignment/Map) stream that records candidate genomic coordinates, alignment scores, and flags. Alignments are typically converted on the fly to the binary **BAM** format.

*Why it matters*: Accurate mapping is a prerequisite to reliable variant detection. Each read receives a mapping quality score indicating confidence in its genomic position.

*Representative SAM header*:

```
@SQ SN:chr1 LN:248956422
@SQ SN:chr2 LN:242193529
@PG ID:bwa-mem2 PN:bwa-mem2 VN:2.2.1 CL:bwa-mem2 mem -t 16 ...
```

*Outputs*: `aligned/<sample>.sam` (usually streamed downstream).

### 4. Sorting & Indexing — [samtools](http://www.htslib.org/)

`samtools sort` arranges BAM records by genomic coordinate and writes an index (`.bai`) to enable random access.

*Why it matters*: Variant callers expect coordinate-sorted and indexed BAM files to quickly fetch reads overlapping regions of interest.

*Outputs*: `bam/<sample>.sorted.bam`, `bam/<sample>.sorted.bam.bai`.

### 5. Duplicate Marking — [samtools markdup](http://www.htslib.org/)

PCR amplification ([PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction)) and optical artifacts can yield **duplicate reads**—multiple observations of the same DNA fragment captured more than once. `samtools markdup` flags these duplicates so that callers can ignore them.

*Why it matters*: Treating duplicates as independent evidence biases allele counts and may cause false positives.

*Outputs*: `bam/<sample>.mkdup.bam`, `bam/<sample>.mkdup.bam.bai`.

### 6. Optional Base Quality Score Recalibration — [GATK](https://gatk.broadinstitute.org/hc/en-us) BQSR

Sequencing machines sometimes misestimate base quality scores. **Base Quality Score Recalibration (BQSR)** uses known variant sites to adjust these scores.

*Why it matters*: More accurate base qualities improve probabilistic models in downstream callers. This step is optional and can be skipped to save time.

*Outputs*: recalibrated BAM (if applied) plus BQSR reports.

### 7. CRAM Conversion & Coverage — samtools + [mosdepth](https://github.com/brentp/mosdepth)

After duplicate marking (and optional BQSR) the pipeline compresses BAM files to CRAM and runs `mosdepth` to summarize per-base coverage.

*Why it matters*: CRAM greatly reduces disk usage while coverage metrics reveal under‑covered regions and overall depth.

*Outputs*: `bam/<sample>.mkdup.cram`, `bam/<sample>.mkdup.cram.crai`, `bam/<sample>.mosdepth.summary.txt`.

### 8. Sharding the Reference Genome

Large genomes are divided into manageable **[BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)** intervals ("shards"). Each shard is processed independently to leverage parallelism.

*Why it matters*: Short-read variant callers parallelize poorly within a single region; sharding provides coarse-grained parallelism across chromosomes/blocks and reduces wall-clock time.

*Outputs*: `vcf/shards/<sample>/<sample>.part_XX.vcf.gz` and corresponding `.tbi` indexes.

### 9. Variant Calling

At this stage, aligned reads are converted into genomic variants. Two backends are available:

#### A. [GATK](https://gatk.broadinstitute.org/hc/en-us) HaplotypeCaller → GenotypeGVCFs

1. **HaplotypeCaller** analyzes each shard, performing local re-assembly to produce per-sample genomic VCFs (**gVCFs**).
2. **GenotypeGVCFs** merges gVCFs and emits a final VCF with genotypes.

*Outputs*: `vcf/<sample>.g.vcf.gz`, `vcf/<sample>.vcf.gz` and their indexes.

#### B. [BCFtools](https://samtools.github.io/bcftools/) mpileup → call → norm

1. **mpileup** summarizes read evidence at each position, applying mapping-quality (MAPQ) and base-quality (BASEQ) filters.
2. **call** infers SNPs and indels using a multiallelic model.
3. **norm** left-aligns indels and splits multiallelic records for compatibility.

*Outputs*: sharded VCFs (`part_XX.vcf.gz`) and final merged VCF (`vcf/<sample>.vcf.gz`).

### 10. Concatenation of Shards

Per-shard VCFs are concatenated in chromosome order and re-indexed so that downstream tools see a seamless, genome-wide VCF.

*Outputs*: `vcf/<sample>.vcf.gz`, `vcf/<sample>.vcf.gz.tbi`.

### 11. Optional VCF Quality Control

[`bcftools stats`](https://samtools.github.io/bcftools/howtos/stats.html) and related utilities generate summary metrics and plots to assess callset quality.

*Outputs*: `qc/<sample>.bcftools.stats.txt` and optional graphical summaries via `plot-vcfstats` or `MultiQC`.

### 12. Gene List & Coverage Reports

`genomes_analyzer.py` extracts gene coordinates from the reference GTF and, using `mosdepth`, calculates breadth and depth per gene for each sample.

*Why it matters*: Gene-level summaries highlight targets with insufficient coverage and facilitate downstream presence/absence analyses.

*Outputs*: `genes/genes.bed`, `genes/<sample>_gene_presence.tsv`.

### 13. Pairwise & Trio Comparisons

For cohorts or family trios, the pipeline can compare VCFs pairwise and flag candidate de novo variants absent in the parents.

*Why it matters*: Automated comparison streamlines interpretation and quality control across samples.

*Outputs*: reports in `comparisons/` and `trio/` directories.

### 14. Optional RNA-seq Module

If RNA-seq samples are defined in the YAML, a lightweight expression pipeline (HISAT2 → StringTie → gffcompare) is executed after DNA analysis.

*Outputs*: transcript assemblies and expression tables under `rna/`.

### Output overview

| Path | Type | Description |
|------|------|-------------|
| `bam/<sample>.mkdup.bam` | BAM | Coordinate-sorted, duplicate-marked alignments |
| `bam/<sample>.mkdup.cram` | CRAM | Compressed alignments with index |
| `bam/<sample>.mosdepth.summary.txt` | TXT | Coverage summary (mosdepth) |
| `vcf/shards/<sample>/part_XX.vcf.gz` | VCF | Per-shard variant calls |
| `vcf/<sample>.vcf.gz` | VCF | Final, genome-wide variant calls |
| `genes/<sample>_gene_presence.tsv` | TSV | Per-gene coverage report |
| `qc/<sample>_R1_fastqc.html` | HTML | Read quality report |

---

## How to Use the Genomes Analyzer

### Uninstallation (if desired)

```bash
CONDA_BASE="/home/lume2/miniforge3"
if ! command -v mamba >/dev/null 2>&1; then
  conda install -n base -c conda-forge -y mamba
fi

conda activate
conda env remove --name genomics
conda deactivate
```

### Installation

```bash
CONDA_BASE="/home/lume2/miniforge3"
if ! command -v mamba >/dev/null 2>&1; then
  conda install -n base -c conda-forge -y mamba
fi

conda activate
./install_genomics_env.sh
source vep_install.sh
```

### Starting the environment

Leave any active conda environment and initialize the session:

```bash
conda deactivate
source start_genomics.sh
```

### Running the pipeline

Execute the workflow by pointing the script to your YAML file:

```bash
conda deactivate
source start_genomics.sh
./genomes_analyzer.py --config config_human_30x_low_memory.yaml
```

The console prints progress panels, including per-shard heartbeats when variant calling is parallelized.

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

---

## Conclusion

Genomes Analyzer provides a clear, reproducible path from raw reads to variant calls using a dependable stack of open‑source tools. By letting users switch between GATK and BCFtools and by exposing explicit sharding and parallelization controls, it adapts to many datasets and machines. The detailed documentation and YAML configuration aim to demystify genomics processing for newcomers while remaining efficient for experienced practitioners.

---

## Appendix 1 — Tools & typical usage

The pipeline wraps several established command‑line programs. Table 1 summarizes their roles and highlights common parameters.

| Tool | Role | Typical command | Key parameters |
|------|------|----------------|----------------|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Read quality assessment | `fastqc -t 8 fastq/*_R{1,2}.fastq.gz -o qc/` | `-t` threads |
| [cutadapt](https://cutadapt.readthedocs.io/) | Adapter and quality trimming | `cutadapt -j 8 -q 20,20 -m 30 -a ADAPT1 -A ADAPT2 -o out_R1.fastq.gz -p out_R2.fastq.gz in_R1.fastq.gz in_R2.fastq.gz` | `-q` quality cutoff, `-m` min length, `-a/-A` adapters |
| [BWA‑MEM2](https://github.com/bwa-mem2/bwa-mem2) | Short‑read alignment | `bwa-mem2 mem -t 16 -R "@RG\tID:S\tSM:S" refs.fa R1.fq.gz R2.fq.gz` | `-t` threads, `-R` read group |
| [samtools](http://www.htslib.org/) | BAM/CRAM processing | `samtools sort -@8 -o out.bam in.bam` | `-@` threads |
| [BCFtools](https://samtools.github.io/bcftools/) | Variant calling and manipulation | `bcftools mpileup -f refs.fa -q20 -Q20 | bcftools call -m -v -Oz -o out.vcf.gz` | `-q/-Q` quality filters, `-m` multiallelic model |
| [GATK](https://gatk.broadinstitute.org/hc/en-us) | Variant calling (Java) | `gatk --java-options "-Xmx24g" HaplotypeCaller -R refs.fa -I in.bam -O out.g.vcf.gz` | `--java-options` memory, `-L` intervals |

Each tool offers many more options; consult the official manuals for advanced usage. The examples above match how `genomes_analyzer.py` invokes them under default settings.

### FastQC

`FastQC` evaluates base-call quality, sequence content and other metrics for raw reads. The example command above scans paired FASTQ files and writes an HTML report and a compressed archive to `qc/`. Successful runs produce files named `*_fastqc.html` and `*_fastqc.zip`.

### cutadapt

`cutadapt` removes adapter sequences and trims low-quality bases. Running the sample command yields trimmed FASTQ files (`out_R1.fastq.gz`, `out_R2.fastq.gz`) and a log describing how many reads were modified or discarded.

### BWA‑MEM2

`bwa-mem2` aligns reads to a reference genome. The output is typically piped to `samtools` to create BAM files. Expect alignment statistics on `stderr` and an alignment stream on `stdout`.

### samtools

`samtools` manipulates SAM/BAM/CRAM files. The sorting example writes a coordinate-sorted BAM (`out.bam`) and reports progress percentages while running. Additional subcommands handle indexing, duplicate marking and more.

### BCFtools

`BCFtools` calls and filters variants. The example pipeline (`mpileup` → `call`) emits a compressed VCF (`out.vcf.gz`) and writes call statistics to `stderr`. The `norm` step (not shown) normalizes indels and splits multiallelic records.

### GATK

The Java-based `GATK` suite provides sophisticated variant calling. `HaplotypeCaller` generates per-sample gVCFs, and tools like `GenotypeGVCFs` merge them. Outputs include `out.g.vcf.gz` and corresponding indexes, with detailed progress logs in the console.

