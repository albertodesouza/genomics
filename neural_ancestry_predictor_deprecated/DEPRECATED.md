# Deprecated

`neural_ancestry_predictor_deprecated` is deprecated and kept only for historical reproducibility. The old active package name `neural_ancestry_predictor` should not be used for new work.

Use `genotype_based_predictor` for new runs:

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/neural_legacy/genes_1000_all.yaml
python3 -m genomics_cli genotype train genotype_based_predictor/configs/neural_legacy/pigmentation_binary.yaml
```

Migration status:

| Historical config | Replacement |
| --- | --- |
| `neural_ancestry_predictor_deprecated/configs/default.yaml` | `genotype_based_predictor/configs/neural_legacy/default.yaml` |
| `neural_ancestry_predictor_deprecated/configs/default_genes.yaml` | `genotype_based_predictor/configs/neural_legacy/default_genes.yaml` |
| `neural_ancestry_predictor_deprecated/configs/genes_1000.yaml` | `genotype_based_predictor/configs/neural_legacy/genes_1000.yaml` |
| `neural_ancestry_predictor_deprecated/configs/genes_1000_all.yaml` | `genotype_based_predictor/configs/neural_legacy/genes_1000_all.yaml` |
| `neural_ancestry_predictor_deprecated/configs/pigmentation_binary.yaml` | `genotype_based_predictor/configs/neural_legacy/pigmentation_binary.yaml` |
| `neural_ancestry_predictor_deprecated/configs/pigmentation_binary_single_gene_screen.yaml` | `genotype_based_predictor/configs/neural_legacy/pigmentation_binary_single_gene_screen.yaml` |

The replacements use `tensor_layout: raw_center_crop` to reproduce the historical unaligned AlphaGenome center-crop representation. Canonical aligned experiments should use `genotype_based_predictor` configs with `tensor_layout: haplotype_channels` and `alignment_mapping: bcftools_chain`.

The migrated path has been smoke-tested with a short pigmentation run through:

```bash
python3 -m genomics_cli genotype train genotype_based_predictor/configs/neural_legacy/pigmentation_binary.yaml
python3 -m genomics_cli genotype evaluate genotype_based_predictor/configs/neural_legacy/pigmentation_binary.yaml --checkpoint best_accuracy --split test
```

Do not add new experiments here. Add them under `genotype_based_predictor/configs/` instead.

Single-gene screen replacement:

```bash
python3 -m genomics_cli genotype single-gene-screen \
  genotype_based_predictor/configs/neural_legacy/pigmentation_binary_single_gene_screen.yaml \
  --dry-run
```
