# Training And Running Predictions

After the 1000 Genomes dataset is available, you can train neural genotype models or run classical SNP-based superpopulation prediction.

Run the commands below in order for the neural path. Use the classical SNP section when you want a VCF/SNP baseline instead of AlphaGenome neural tensors.

## Neural Genotype Models

Prepare or validate the processed tensor cache:

```bash
genomics genotype split configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
```

Train a model:

```bash
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
```

Evaluate a saved checkpoint:

```bash
genomics genotype evaluate \
  configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml \
  --checkpoint best_accuracy \
  --split val
```

Run final held-out test evaluation after model selection:

```bash
genomics genotype test \
  configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml \
  --checkpoint best_accuracy
```

## Enabling DITA

DITA is the aligned tensor representation for indel-aware neural models. Enable it when you want to address coordinate shifts caused by indels in dense AlphaGenome tracks.

Set these fields in the genotype config:

```yaml
dataset_input:
  tensor_layout: "haplotype_channels"
  alignment_mapping: "bcftools_chain"
  feature_mode: "signals_and_masks"
  indel_include_valid_mask: true
  indel_include_snp_mask: true
```

Validate required artifacts before training:

```bash
genomics audit-data \
  --dataset-id 1kg_high_coverage \
  --check-bcftools-chain \
  --sample-limit 3 \
  --fail-on-missing
```

If artifacts exist in a separate tree, preview syncing them with:

```bash
genomics genotype sync-bcftools-artifacts \
  --source-dir /path/to/consensus_dataset \
  --target-dir /path/to/canonical_dataset \
  --link-mode hardlink
```

Add `--apply` only after reviewing the preview.

## Classical SNP Superpopulation Path

If you do not need AlphaGenome neural tensors, run the integrated SNP ancestry pipeline:

```bash
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml
```

For FROGAncestryCalc likelihood comparisons, use the third-party tool directly:

```bash
cd third_party/FROGAncestryCalc
python3 tools/vcf_to_frog.py sample.vcf.gz tools/aisnps_55_list.txt input/sample_data.txt
./run.sh
```

## More Details

| Topic | Read More |
|---|---|
| Genotype model inputs, caches, and outputs | [Genotype Predictor](../components/genotype-predictor.md) |
| Why DITA exists and what masks mean | [Dynamic Indel Tensor Alignment](../concepts/dynamic-indel-tensor-alignment.md) |
| Split policy, preprocessing, stability, and leakage | [Preprocessing And Leakage](../concepts/preprocessing-and-leakage.md) |
| SNP ancestry math and outputs | [SNP Ancestry](../components/snp-ancestry.md) |
| FROGAncestryCalc and native/third-party tools | [Native And Third-Party](../components/native-and-third-party.md) |
| All CLI flags | [CLI Reference](../reference/cli.md) |
