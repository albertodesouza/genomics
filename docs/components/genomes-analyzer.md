# Genomes Analyzer

The Genomes Analyzer workflow is an operational sequencing pipeline for taking configured DNA samples from raw reads to aligned reads, coverage summaries, variant calls, and downstream relationship/gene reports. It is designed for large human short-read datasets such as 1000 Genomes high-coverage samples, but the implementation also has optional paths for local FASTQs, long-read alignment, supervised ancestry, and RNA-seq.

## CLI

```bash
genomics genomes-analyzer run --config configs/genomes_analyzer/config_human_30x_low_memory.yaml
```

The command loads the YAML config, normalizes it to the internal schema, creates the configured `storage.base_dir`, changes into that directory, and writes all runtime artifacts below it.

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/workflows/genomes_analyzer/` |
| Configs | `configs/genomes_analyzer/` |
| Operations scripts | `scripts/ops/` |
| Diagnostics | `scripts/diagnostics/` |

## Pipeline Overview

The default DNA path is:

```text
YAML config
  -> reference FASTA/GTF download and indexes
  -> SRA/ENA or local FASTQ staging
  -> optional FASTQ downsampling
  -> FastQC, MultiQC, cutadapt trimming
  -> BWA/BWA-MEM2 alignment and coordinate-sorted BAM
  -> Picard duplicate marking and BAM/CRAM indexing
  -> mosdepth genome coverage
  -> BCFtools or GATK variant calls, sharded by reference intervals
  -> optional VEP annotation
  -> gene list, per-gene coverage, pairwise/trio/paternity/ancestry reports
```

Each stage checks for its expected outputs and reuses them when they are current. This matters because complete genome runs can take hours or days. Most generated files live in stable subdirectories such as `refs/`, `raw/`, `fastq/`, `fastq_ds/`, `qc/`, `trimmed/`, `bam/`, `vcf/`, `genes/`, `comparisons/`, `rnaseq/`, and `logs/` under `storage.base_dir`.

## Step-By-Step Data Flow

| Step | Main inputs and origin | Format | Main outputs | Format |
|---|---|---|---|---|
| 1. References and indexes | `project.reference.fasta_url` and `project.reference.gtf_url` from the config. Optional `bwa_index_url` can provide a prebuilt aligner index. | Reference genome as FASTA, usually compressed upstream as `.fa.gz`; gene annotation as GTF, usually `.gtf.gz`; optional BWA/BWA-MEM2 index archive. | `refs/reference.fa`, `refs/genes.gtf`, `refs/reference.fa.fai`, sequence dictionary, BWA/BWA-MEM2 index files, optional HISAT2 index for RNA-seq. | FASTA, GTF, FAI text index, `.dict`, BWA files such as `.amb`, `.ann`, `.bwt`, `.pac`, `.sa`, or BWA-MEM2 files such as `.bwt.2bit.64`/`.0123`. |
| 2. FASTQ acquisition | `samples[].runs` accessions in the config, or local FASTQ paths in legacy-style configs. SRA data comes from SRA Toolkit `prefetch`; optional ENA fallback queries the ENA read-run API for FASTQ URLs. | SRA container (`.sra`) or paired/single FASTQ files (`.fastq.gz`/`.fq.gz`). | `raw/<accession>/...` for SRA cache, `fastq/<accession>_1.fastq.gz`, `fastq/<accession>_2.fastq.gz`, or copied local FASTQs. | `.sra`, optional `.vdbcache`, gzipped FASTQ. |
| 3. Space estimate | FASTQs already staged in `fastq/`. | `.fastq.gz`/`.fq.gz`. | Console warning/summary only. | No persistent data product beyond logs. |
| 4. Optional downsample | FASTQs in `fastq/`, controlled by `size_control.downsample`. | Gzipped FASTQ. | Downsampled reads in `fastq_ds/`, typically named with `.ds.fastq.gz`. | Gzipped FASTQ. |
| 5. QC and trimming | FASTQs from `fastq_ds/` when downsampling produced files; otherwise from `fastq/`. Adapter sequences come from normalized config defaults or config overrides. | Gzipped FASTQ. | FastQC reports in `qc/`, MultiQC report in `qc/`, trimmed reads in `trimmed/`. | FastQC `.html` and `.zip`, MultiQC `.html`, trimmed `.trim.fq.gz`. |
| 6. Alignment, sort, duplicates | Trimmed FASTQs from `trimmed/` when present; otherwise downsampled/raw FASTQs. Reference and aligner indexes from `refs/`. | FASTQ plus FASTA/index files. | `bam/<sample>.sorted.bam`, `bam/<sample>.sorted.bam.bai`, `bam/<sample>.mkdup.bam`, `bam/<sample>.mkdup.bam.bai`, duplicate metrics. Sorted BAMs may be removed depending on cleanup settings. | BAM, BAI, Picard metrics text. |
| 7. CRAM and whole-genome coverage | Duplicate-marked BAM or CRAM in `bam/`; `refs/reference.fa` is required when writing or reading CRAM. | BAM/CRAM plus FASTA. | `bam/<sample>.mkdup.cram`, `bam/<sample>.mkdup.cram.crai`, mosdepth summaries such as `bam/<sample>.mosdepth.summary.txt` and `bam/<sample>.mosdepth.global.dist.txt`. | CRAM, CRAI, mosdepth text summaries and distribution files. |
| 8. Variant calling | Duplicate-marked BAM/CRAM, `refs/reference.fa`, `refs/reference.fa.fai`, and generated interval BED shards. Caller is selected by `params.variant_caller`. | BAM/CRAM, FASTA, FAI, BED. | BCFtools path: `vcf/shards/<sample>/<sample>.part_XX.vcf.gz` plus indexes and final `vcf/<sample>.vcf.gz`. GATK path: optional sharded `*.g.vcf.gz`, final `vcf/<sample>.g.vcf.gz`, and genotyped `vcf/<sample>.vcf.gz`. | BGZF-compressed VCF/gVCF (`.vcf.gz`, `.g.vcf.gz`) plus Tabix indexes (`.tbi`). |
| 9. Optional VEP annotation | Final per-sample VCFs in `vcf/`, local VEP installation/cache configured by `params.vep_*`. | Indexed `.vcf.gz`. | Annotated VCFs under `vep/` when `params.annotate_with_vep` and normalized annotation settings enable this path. | VCF or compressed VCF depending on VEP settings. |
| 10. Gene list | `refs/genes.gtf` from the configured GTF URL. | GTF. | `genes/gene_list.txt`. | Plain text list, one gene name per line. |
| 11. Gene presence and coverage | Duplicate-marked BAM/CRAM, `refs/genes.gtf`, generated `genes/genes.bed`, and `refs/reference.fa` for CRAM. | BAM/CRAM, GTF-derived BED. | `genes/genes.bed`, `genes/<sample>.genes.regions.bed.gz`, `genes/<sample>.genes.thresholds.bed.gz`, `genes/<sample>_gene_presence.tsv`. | BED, gzipped mosdepth BED-like outputs, TSV. |
| 12. Pairwise comparisons | Final per-sample VCFs in `vcf/`. | Indexed `.vcf.gz`. | `comparisons/<A>_vs_<B>.merge.vcf.gz`, `comparisons/<A>_vs_<B>.metrics.tsv`, `comparisons/<A>_vs_<B>.summary.md`, `comparisons/summary_pairwise.csv`. | Compressed VCF, Tabix index, TSV, Markdown, CSV. |
| 13. Trio de novo | Final VCFs and trio IDs from `general.trio_child_id` and `general.trio_parent_ids`, or fallback heuristics. | Indexed `.vcf.gz`. | Trio candidate report files under the workflow output directories. | Usually TSV/Markdown summaries derived from VCF genotypes. |
| 14. Paternity | Final VCFs and trio configuration. | Indexed `.vcf.gz`. | Paternity likelihood/summary reports under the workflow output directories. | TSV/Markdown-style reports. |
| 15. Optional ancestry | Final VCF-derived genotype data plus configured ancestry assets. The code uses PLINK/ADMIXTURE-style tooling when enabled. | VCF and PLINK-style genotype/reference assets, depending on config. | Supervised ADMIXTURE outputs and ancestry summaries. | PLINK files such as `.bed/.bim/.fam`, ADMIXTURE outputs such as `.Q/.P`, and text summaries. |
| 16. Optional RNA-seq | `rna_samples` from legacy-style config, `refs/genes.gtf`, and HISAT2 index under `refs/`. | FASTQ, GTF, HISAT2 index. | `rnaseq/<run>.rnaseq.bam`, BAM index, `rnaseq/<run>.transcripts.gtf`, `rnaseq/cmp.*` from gffcompare. | BAM/BAI, GTF, gffcompare text outputs. |

## Detailed Stage Notes

### 1. References And Indexes

The reference FASTA and GTF are declared in the YAML under `project.reference`. The default example uses GENCODE GRCh38 URLs. The workflow downloads them into `refs/`, normalizes the FASTA name to `refs/reference.fa`, and normalizes the annotation to `refs/genes.gtf`.

The FASTA drives every coordinate-based downstream step: alignment, CRAM conversion, mosdepth with CRAM, BCFtools/GATK variant calling, VCF normalization, and interval ordering. The GTF drives gene list generation and per-gene coverage. Keep the GTF assembly consistent with the FASTA; mixing GRCh37 and GRCh38, or using a FASTA with different contig names, will produce empty or invalid downstream results.

Index outputs depend on the selected aligner:

| Tool | Input | Output |
|---|---|---|
| `samtools faidx` | `refs/reference.fa` | `refs/reference.fa.fai` |
| GATK/Picard dictionary tooling | `refs/reference.fa` | sequence dictionary such as `refs/reference.dict` or equivalent |
| `bwa index` | `refs/reference.fa` | `refs/reference.fa.amb`, `.ann`, `.bwt`, `.pac`, `.sa` |
| `bwa-mem2 index` | `refs/reference.fa` | BWA-MEM2 index files such as `.bwt.2bit.64` and related files |
| `hisat2-build` when RNA-seq is enabled | `refs/reference.fa` | `refs/<assembly_name>.*.ht2` index files |

If `limit_to_canonical` is enabled, the reference can be restricted to canonical contigs before indexing. Do not combine that option with a prebuilt index unless the prebuilt index was created from the exact same restricted FASTA.

### 2. Read Acquisition

The normalized user-facing config uses `samples` like this:

```yaml
samples:
  - sample_id: "NA12878"
    study: "PRJEB31736"
    runs:
      - "ERR3239334"
```

For each run accession, the workflow can:

| Source path | How it is obtained | Result |
|---|---|---|
| SRA Toolkit | `prefetch <accession>` stores an `.sra` object under `raw/`; `fasterq-dump --split-files` converts it. | `fastq/<accession>_1.fastq.gz` and, for paired-end runs, `fastq/<accession>_2.fastq.gz`. |
| ENA fallback | The ENA filereport API returns `fastq_http`/`fastq_ftp` URLs; the workflow downloads matching `_1.fastq.gz` and `_2.fastq.gz` files and can validate gzip/MD5 when metadata is available. | `fastq/<accession>_1.fastq.gz`, `fastq/<accession>_2.fastq.gz`. |
| Local FASTQ | Legacy-style configs may provide `fastq1`/`fastq2`; the files are copied or staged into the run directory. | Files under `fastq/` with the original FASTQ basename. |

FASTQ files are four-line-per-read text records, usually gzip-compressed in this workflow. Paired-end data is detected from filename patterns such as `_1.fastq.gz` and `_2.fastq.gz`.

### 3. QC, Downsampling, And Trimming

Downsampling is controlled by `size_control.downsample.enabled`, `fraction`, and `seed`. It creates a reduced read set under `fastq_ds/`. This is useful for smoke tests or low-memory runs; it should be disabled or set to `fraction: 1.0` for production-depth analysis.

QC and trimming use:

| Tool | Input | Output |
|---|---|---|
| FastQC | Raw or downsampled FASTQ | `qc/<read>_fastqc.html`, `qc/<read>_fastqc.zip` |
| MultiQC | The `qc/` directory | `qc/multiqc_report.html` |
| cutadapt | Raw/downsampled FASTQ plus adapter sequences | `trimmed/<run>_1.trim.fq.gz`, `trimmed/<run>_2.trim.fq.gz` |

After trimming, FastQC/MultiQC are run again for trimmed FASTQs when reports are missing or stale.

### 4. Alignment And Duplicate Marking

For short reads, the workflow aligns with `bwa-mem2 mem` when a BWA-MEM2 index is available, otherwise it can fall back to `bwa mem`. The aligner emits SAM on stdout; the pipeline streams that into `samtools view` and `samtools sort`, producing a coordinate-sorted BAM.

The main file transformations are:

```text
trimmed FASTQ
  -> streamed SAM
  -> unsaved BAM stream
  -> bam/<sample>.sorted.bam
  -> bam/<sample>.mkdup.bam
```

`picard MarkDuplicates` creates `bam/<sample>.mkdup.bam` and duplicate metrics. The duplicate-marked file is the preferred input for coverage and variant calling. Long-read configs can use the minimap2 path, but the common human 30x configs use short-read BWA/BWA-MEM2.

### 5. CRAM And Coverage

When `general.use_cram` is true, `samtools view -C -T refs/reference.fa` converts duplicate-marked BAM to CRAM:

```text
bam/<sample>.mkdup.bam
  -> bam/<sample>.mkdup.cram
  -> bam/<sample>.mkdup.cram.crai
```

CRAM is smaller than BAM but is reference-dependent. Any later tool reading CRAM must use the same `refs/reference.fa` used to create it.

The same stage runs mosdepth for whole-genome coverage. Typical outputs include:

| Output | Meaning |
|---|---|
| `bam/<sample>.mosdepth.summary.txt` | Per-contig and total coverage summary. |
| `bam/<sample>.mosdepth.global.dist.txt` | Global depth distribution. |
| Other `bam/<sample>.mosdepth.*` files | Additional mosdepth distribution/region files depending on mosdepth behavior and options. |

### 6. Variant Calling

Variant calling starts from `bam/<sample>.mkdup.bam` or `bam/<sample>.mkdup.cram`. The reference index `refs/reference.fa.fai` is used to derive contig order and balanced interval shards. Each shard is represented as a BED file.

The BCFtools path uses:

```text
aligned BAM/CRAM + refs/reference.fa + shard BED
  -> bcftools mpileup
  -> bcftools call
  -> bcftools norm
  -> vcf/shards/<sample>/<sample>.part_XX.vcf.gz
  -> vcf/<sample>.vcf.gz
```

The GATK path uses:

```text
aligned BAM/CRAM + refs/reference.fa + optional shard BED
  -> HaplotypeCaller
  -> vcf/<sample>.g.vcf.gz or sharded gVCFs
  -> GenotypeGVCFs
  -> vcf/<sample>.vcf.gz
```

All final VCFs are BGZF-compressed and indexed with Tabix (`.tbi`). That indexing is required for downstream BCFtools queries, pairwise merges, trio analysis, and annotation.

### 7. Annotation With VEP

When VEP annotation is enabled, the workflow reads final indexed VCFs from `vcf/` and calls Ensembl VEP using `params.vep_species`, `params.vep_assembly`, `params.vep_dir_cache`, `params.vep_fork`, and `params.vep_extra`. VEP requires a local VEP installation and compatible cache. The input VCF assembly must match the VEP cache assembly.

### 8. Gene Reports

Gene reporting has two layers:

| Layer | Input | Output |
|---|---|---|
| Gene list | `refs/genes.gtf` | `genes/gene_list.txt` containing unique gene names. |
| Gene coverage/presence | `refs/genes.gtf`, duplicate-marked BAM/CRAM | `genes/genes.bed`, mosdepth per-gene files, `genes/<sample>_gene_presence.tsv`. |

`genes/genes.bed` is derived from GTF `gene` features. Its rows contain chromosome, zero-based start, end, gene ID, gene name, gene type/biotype, and description. mosdepth then runs with `--by genes/genes.bed --thresholds 1,5,10`. The final presence TSV contains mean coverage, breadth at 1x, and a boolean `present` call based on `general.gene_presence_min_mean_cov` and `general.gene_presence_min_breadth_1x`.

### 9. Relationship Reports

Pairwise comparison reads per-sample VCFs, merges each pair with `bcftools merge`, extracts genotypes with `bcftools query`, and writes metrics such as exact genotype concordance, allele sharing, and IBS0 fraction.

The trio and paternity stages also consume final VCFs. Trio IDs should be configured explicitly with `general.trio_child_id` and `general.trio_parent_ids` when possible. If they are not configured, the code falls back to heuristics such as using `NA12878` as a child when present, or the first three configured sample IDs.

### 10. Optional Ancestry

The ancestry step is disabled in the default example. When enabled and configured with the required assets, it uses PLINK/ADMIXTURE-style processing. Input starts from VCF-derived genotypes and external reference/population assets, and output includes PLINK binary files, ADMIXTURE `.Q/.P` matrices, and text summaries.

### 11. Optional RNA-seq

RNA-seq is also disabled in the normalized default config. Legacy-style configs can provide `rna_samples`. For each RNA run, FASTQs are staged like DNA reads, aligned with HISAT2 against `refs/<assembly_name>`, sorted with samtools, assembled with StringTie against `refs/genes.gtf`, and compared with gffcompare.

The main RNA-seq artifacts are:

| Artifact | Format |
|---|---|
| `rnaseq/<run>.rnaseq.bam` and index | BAM/BAI |
| `rnaseq/<run>.transcripts.gtf` | GTF |
| `rnaseq/cmp.*` | gffcompare text/stat outputs |

## Important Config Sections

| Section | Purpose |
|---|---|
| `project` | Project name, organism label, reference metadata. |
| `project.reference` | FASTA, GTF, and optional aligner index URLs. |
| `storage` | Base directory, temporary directory, and intermediate retention. |
| `download` | SRA Toolkit download behavior. |
| `execution` | Threads, resume behavior, ENA fallback, aligner choice, progress/stall settings. |
| `size_control` | Optional read downsampling. |
| `general` | Operational flags such as CRAM usage, reference/index forcing, trio IDs, and gene presence thresholds. |
| `params` | Variant caller, BCFtools/GATK sharding, quality thresholds, memory, and VEP settings. |
| `samples` | DNA sample IDs and SRA/ENA run accessions. |
| `steps` | Enables or disables individual pipeline stages. |
| `ancestry` | Optional supervised ADMIXTURE settings. |

## Output Directory Map

The exact set of files depends on enabled steps, caller choice, and cleanup settings, but a typical run under `storage.base_dir` looks like this:

```text
refs/           reference FASTA, GTF, FASTA indexes, aligner indexes
raw/            SRA prefetch cache and related download artifacts
fastq/          staged raw FASTQs
fastq_ds/       optional downsampled FASTQs
qc/             FastQC and MultiQC reports
trimmed/        cutadapt-trimmed FASTQs
bam/            sorted/mkdup BAMs, CRAMs, indexes, mosdepth coverage
vcf/            final VCFs plus shard directories
vep/            optional VEP-annotated VCFs
genes/          gene list, gene BED, per-gene mosdepth, presence TSVs
comparisons/    pairwise merged VCFs and relationship metrics
rnaseq/         optional RNA-seq BAM/GTF/gffcompare outputs
logs/           tool logs for downloads, conversions, and variant shards
intervals/      generated BED interval files for variant calling
```

## Resume And Failure Recovery

The workflow is designed around resumability:

- output existence and modification time are checked before expensive recomputation;
- state markers record completed steps where applicable;
- variant shards are written independently, so failed shards can be inspected via `vcf/shards/<sample>/` and `logs/bcftools_<sample>_part_XX.log`;
- SRA/FASTQ download logs are kept in `logs/`;
- diagnostic scripts under `scripts/diagnostics/` help investigate common failures such as empty BCFtools outputs, stalled downloads, and reference mismatches.

When debugging, first verify that the reference, BAM/CRAM, and VCF contig names agree. Most downstream failures in this workflow are caused by mismatched assemblies, missing indexes, incomplete FASTQ downloads, or CRAM reads being evaluated without the correct FASTA.

## External Tools

The workflow invokes external command-line tools as subprocesses. The environment must provide the tools needed by the enabled stages. Common DNA tools include `prefetch`, `fasterq-dump`, `wget`/`curl`, `fastqc`, `multiqc`, `cutadapt`, `bwa` or `bwa-mem2`, `samtools`, `picard`, `mosdepth`, `bcftools`, `tabix`, and optionally `gatk` and `vep`. Optional ancestry and RNA-seq stages need tools such as `plink`, `admixture`, `hisat2`, `stringtie`, and `gffcompare`.

Keep `genomics --help` lightweight: optional scientific and ML dependencies should remain lazy and should not be required just to inspect CLI help.
