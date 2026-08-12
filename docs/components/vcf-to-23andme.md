# VCF To 23andMe Converter

Converts VCF files into 23andMe raw data format.

## CLI

```bash
genomics convert vcf-to-23andme --config configs/converters/vcf_to_23andme/default.yaml
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/converters/vcf_to_23andme/` |
| Configs | `configs/converters/vcf_to_23andme/` |

## Features

- automatic rsID annotation with dbSNP through BCFtools;
- GRCh37/GRCh38 auto-detection from VCF metadata;
- V3 and V5 output formats;
- optional chip-panel filtering;
- multi-sample VCF support;
- read-only handling of the input VCF.

## Implementation Structure

The converter is implemented primarily in `converter.py` and exposed through `__main__.py` for CLI delegation.

Key responsibilities in `converter.py`:

| Area | Implementation detail |
|---|---|
| Build detection | Compares contig metadata and chromosome lengths against GRCh37/GRCh38 constants |
| dbSNP support | Downloads/caches dbSNP VCFs and indexes for GRCh37 or GRCh38 |
| Annotation | Uses BCFtools to annotate missing rsIDs when configured |
| Chromosome normalization | Handles `chr1`, `1`, RefSeq accessions, `X/Y/MT` variants |
| Chip panels | Downloads/caches Illumina panel manifests and extracts rsIDs |
| Output formatting | Emits 23andMe-like V3/V5 headers and tab-separated SNP records |

## Conversion Flow

```text
input VCF
  -> detect sample and genome build
  -> optionally annotate rsIDs with dbSNP
  -> optionally filter to chip panel or custom SNP panel
  -> normalize chromosome names and genotypes
  -> write 23andMe raw text file
```

## Config Sections

| Section | Purpose |
|---|---|
| `input` | Input VCF path and optional sample name |
| `output` | Output path, format version, rsID/chip-panel behavior |
| `annotation` | dbSNP annotation and reference cache directory |
| `reference_build` | Fallback build when auto-detection is inconclusive |

## External Requirements

dbSNP annotation requires `bcftools` on `PATH`. Use `source scripts/env/start_genomics_universal.sh` before running conversion in the project environment.
