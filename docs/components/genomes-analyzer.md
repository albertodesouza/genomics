# Genomes Analyzer

The Genomes Analyzer workflow processes sequencing data into aligned reads, variant calls, QC reports, ancestry summaries, and optional RNA-seq outputs.

## CLI

```bash
genomics genomes-analyzer run --config configs/genomes_analyzer/config_human_30x_low_memory.yaml
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/workflows/genomes_analyzer/` |
| Configs | `configs/genomes_analyzer/` |
| Operations scripts | `scripts/ops/` |
| Diagnostics | `scripts/diagnostics/` |

## Pipeline Stages

1. Prepare references and indexes.
2. Fetch FASTQs when configured.
3. Run read QC and trimming.
4. Align and sort reads.
5. Mark duplicates.
6. Convert to CRAM and compute coverage.
7. Call variants with BCFtools or GATK.
8. Concatenate and index VCF shards.
9. Build gene lists and per-gene coverage reports.
10. Run pairwise, trio, paternity, and ancestry analyses when enabled.
11. Run optional RNA-seq steps.

## Implementation Structure

The workflow is split into focused modules instead of one large script. The main orchestration lives in `pipeline.py` and the command-line entrypoint lives in `cli.py`.

| Module | Responsibility |
|---|---|
| `cli.py` | Parses workflow arguments and starts the run |
| `config.py` | Loads YAML, normalizes legacy/current config shapes, exposes typed accessors |
| `pipeline.py` | Orchestrates the end-to-end step graph and idempotent execution |
| `state.py` | Tracks completed steps and resumability markers |
| `references.py` | Downloads/prepares reference FASTA, indexes, dictionaries, GTF resources |
| `fastq.py` | Fetches FASTQ inputs from configured sources and manages raw read files |
| `qc.py` | Runs FastQC, cutadapt, MultiQC-style checks, and quality summaries |
| `alignment.py` | Runs BWA/BWA-MEM2, sorting, indexing, duplicate marking, CRAM conversion |
| `variants.py` | Handles sharding, BCFtools/GATK calls, concatenation, VEP annotation |
| `reports.py` | Builds final reports, gene coverage summaries, paternity/trio outputs |
| `ancestry.py` | Runs supervised ADMIXTURE/PLINK ancestry workflows |
| `rnaseq.py` | Optional HISAT2/StringTie/gffcompare RNA-seq flow |
| `space.py` | Disk/memory estimation helpers for large runs |
| `utils.py` | Shared shell execution, logging, path, and small helper utilities |

## Execution Model

The implementation favors explicit files and idempotent checkpoints. Long-running stages check for expected outputs before recomputing. This is important because genome processing jobs can run for hours or days and may be restarted after interruption.

The variant-calling path scatters the reference into BED shards, runs shard jobs independently, then concatenates shards in chromosome order. The BCFtools path uses `mpileup`, `call`, and `norm`; the GATK path uses HaplotypeCaller and GenotypeGVCFs where configured.

External tools are invoked as subprocesses. The workflow assumes the Conda environment supplies tools such as `bwa-mem2`, `samtools`, `bcftools`, `gatk`, `mosdepth`, `plink`, `admixture`, `fastqc`, `cutadapt`, and optionally VEP.

## Data Flow

```text
configured samples
  -> FASTQ download or existing BAM/CRAM
  -> QC/trimming
  -> aligned BAM
  -> duplicate-marked BAM/CRAM
  -> per-shard VCFs
  -> final indexed VCF
  -> gene, paternity, ancestry, trio, QC reports
```

## Resume And Failure Recovery

The workflow is designed around resumability:

- output existence is checked before expensive recomputation;
- state files mark step completion;
- shard outputs can be inspected independently;
- diagnostic scripts under `scripts/diagnostics/` help investigate common failures such as empty BCFtools outputs or reference mismatches.

## Important Config Sections

| Section | Purpose |
|---|---|
| `project` | Project name, organism, reference metadata |
| `general` | Runtime options, trio IDs, paternity options |
| `params` | Aligner, caller, sharding, VEP settings |
| `storage` | Base, work, and temp directories |
| `download` | SRA/ENA download behavior |
| `steps` | Enabled pipeline stages |
| `samples` | Sample IDs and run accessions |
| `ancestry` | Supervised ADMIXTURE settings |

## Outputs

Typical outputs include:

- BAM/CRAM files and indexes;
- per-shard and final VCFs;
- FastQC/MultiQC reports;
- mosdepth coverage summaries;
- gene coverage TSVs;
- paternity and ancestry reports.

## Extension Points

Add new operational stages by keeping the stage logic in a focused module and wiring it through `pipeline.py`. Prefer writing outputs to a stable subdirectory and make the step idempotent by checking for a small set of required output files.
