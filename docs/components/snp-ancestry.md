# SNP Ancestry

The SNP ancestry predictor uses allele frequencies and SNP panels to infer ancestry or configured derived targets.

## CLI

```bash
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/predictors/snp_ancestry/` |
| Configs | `configs/predictors/snp_ancestry/` |
| Default refs location | `configs/predictors/snp_ancestry/refs/` |

## Pipeline Steps

1. Convert multi-sample VCF data into 23andMe-like per-individual files.
2. Compute reference allele-frequency statistics.
3. Predict ancestry with maximum-likelihood or admixture-style methods.

## Implementation Structure

| Module | Responsibility |
|---|---|
| `pipeline.py` | End-to-end conversion, frequency computation, prediction, reporting |
| `generate_regions.py` | Utility for generating region/panel files from configured windows |
| `__main__.py` | Module entrypoint used by the `genomics` CLI delegation |

The implementation reuses the VCF-to-23andMe converter for normalization, chip-panel handling, dbSNP annotation, chromosome naming, and output formatting.

## Config Sections

| Section | Purpose |
|---|---|
| `input` | Splits, individuals directory, multi-sample VCF pattern, chromosomes |
| `conversion` | 23andMe output settings, dbSNP annotation, chip-panel filtering, parallelism |
| `statistics` | Reference subsets, allele-frequency output, SNP panel, MAF/max SNP filters |
| `prediction` | Target level, derived targets, haplotype mode, prediction method, results dir |

## Prediction Methods

The maximum-likelihood classifier scores an individual's observed genotypes against population/reference allele frequencies. The admixture path estimates mixture proportions with constrained optimization when the required optional scientific stack is available.

Haplotype handling is controlled by config. Some experiments use one haplotype; others combine H1 and H2 to approximate diploid genotypes.

## Intermediate Artifacts

Typical outputs include:

- per-individual 23andMe-like files;
- allele-frequency tables by population or target class;
- selected SNP panels or filtered SNP sets;
- prediction JSON/TSV summaries;
- confusion matrices and classification metrics when labels are available.

## Notes

Some historical configs may still refer to legacy `/dados/GENOMICS_DATA/top3` paths. New configs should use canonical datasets and registered dataset IDs when possible.
