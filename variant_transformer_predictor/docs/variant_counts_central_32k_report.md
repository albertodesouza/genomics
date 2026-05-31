# Variant Count Report - Central 32 kb

Dataset analyzed:

```text
/dados/GENOMICS_DATA/variant_transformer/superpopulation
```

This report compares the full 524,288 bp windows currently materialized for each gene against the central 32,768 bp region of each gene/window.

## Overall Summary

| Metric | Value |
|---|---:|
| Samples | 3,202 |
| Genes/windows | 11 |
| Full-window total tokens | 40,870,647 |
| Central 32 kb total tokens | 2,066,495 |
| Full-window mean tokens/sample | 12,764.10 |
| Full-window median tokens/sample | 12,496.50 |
| Full-window p95 tokens/sample | 15,016.85 |
| Full-window p99 tokens/sample | 15,402.95 |
| Full-window max tokens/sample | 16,112 |
| Central 32 kb mean tokens/sample | 645.38 |
| Central 32 kb median tokens/sample | 616 |
| Central 32 kb p95 tokens/sample | 928 |
| Central 32 kb p99 tokens/sample | 992 |
| Central 32 kb max tokens/sample | 1,078 |
| Samples > 4,096 tokens, full-window | 3,202 |
| Samples > 4,096 tokens, central 32 kb | 0 |
| Samples > 1,024 tokens, central 32 kb | 14 |

## Interpretation

The full 524 kb windows are too large for the initial Transformer configuration: every sample exceeds 4,096 tokens.

Using only the central 32 kb of each gene/window reduces the average sequence length from 12,764 tokens to 645 tokens per sample. The maximum central 32 kb sequence length is 1,078 tokens.

The model configs now use:

```yaml
max_sequence_length: 1024
truncate_policy: "keep_first"
```

This means 14 out of 3,202 samples will be truncated by at most 54 tokens in the central 32 kb dataset. If strict non-truncating training is preferred, set `max_sequence_length: 1088` or `max_sequence_length: 1152` instead.

## Per-Gene Counts

| Gene | Full total | Full mean/sample | Full p99 | Full max | Central 32 kb total | Central 32 kb mean/sample | Central 32 kb p99 | Central 32 kb max | Central fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DDB1 | 2,836,991 | 886.01 | 1,240.00 | 1,295 | 106,135 | 33.15 | 77.00 | 93 | 0.0374 |
| EDAR | 4,068,878 | 1,270.73 | 1,524.98 | 1,684 | 160,354 | 50.08 | 120.00 | 194 | 0.0394 |
| HERC2 | 2,532,868 | 791.03 | 1,272.97 | 1,397 | 121,357 | 37.90 | 114.99 | 124 | 0.0479 |
| MC1R | 4,780,875 | 1,493.09 | 2,192.97 | 2,374 | 279,110 | 87.17 | 174.00 | 183 | 0.0584 |
| MFSD12 | 4,550,288 | 1,421.08 | 1,864.95 | 1,981 | 253,438 | 79.15 | 127.99 | 163 | 0.0557 |
| OCA2 | 4,137,653 | 1,292.21 | 1,791.00 | 1,975 | 385,947 | 120.53 | 207.99 | 228 | 0.0933 |
| SLC24A5 | 2,754,528 | 860.25 | 1,182.99 | 1,479 | 72,628 | 22.68 | 65.99 | 77 | 0.0264 |
| SLC45A2 | 3,174,598 | 991.44 | 1,366.00 | 1,475 | 265,428 | 82.89 | 116.00 | 125 | 0.0836 |
| TCHH | 4,172,628 | 1,303.13 | 1,640.00 | 1,988 | 55,740 | 17.41 | 59.00 | 97 | 0.0134 |
| TYR | 2,850,784 | 890.31 | 1,557.00 | 1,674 | 219,238 | 68.47 | 169.99 | 188 | 0.0769 |
| TYRP1 | 5,010,556 | 1,564.82 | 1,894.98 | 2,091 | 147,120 | 45.95 | 121.00 | 163 | 0.0294 |

## Generated Source Files

The runtime analysis generated these files under `/dados`:

```text
/dados/GENOMICS_DATA/variant_transformer/superpopulation/reports/variant_counts_report_central_32768.md
/dados/GENOMICS_DATA/variant_transformer/superpopulation/reports/variant_counts_by_gene_central_32768.csv
/dados/GENOMICS_DATA/variant_transformer/superpopulation/reports/variant_counts_by_sample_central_32768.csv
/dados/GENOMICS_DATA/variant_transformer/superpopulation/reports/variant_counts_central_32768.json
```

To regenerate the analysis:

```bash
python3 -m variant_transformer_predictor.analyze_variant_counts \
  /dados/GENOMICS_DATA/variant_transformer/superpopulation \
  --central-window-size 32768
```

To materialize a central 32 kb dataset for training:

```bash
python3 -m variant_transformer_predictor.materialize_dataset \
  --dataset-dir /dados/GENOMICS_DATA/v1/1kG_high_coverage \
  --vcf-root-dir /dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes \
  --target superpopulation \
  --classes AFR AMR EAS EUR SAS \
  --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation_32k \
  --sample-batch-size 128 \
  --central-window-size 32768
```
