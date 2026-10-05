# CLI Reference

The primary entrypoint is:

```bash
genomics
```

## Visualizer

```bash
genomics visualize                                     # default dataset and runs root, http://127.0.0.1:8780
genomics visualize --dataset /path/to/dataset --open   # any canonical-layout dataset
genomics visualize --dataset-id 1kg_high_coverage --annotations phenotypes.tsv --memory-mb 8192
```

Options: `--dataset DIR` / `--dataset-id ID` (repeatable), `--annotations TABLE`, `--runs-root DIR` (repeatable), `--gtf TABLE`, `--consensus-dataset-dir DIR`, `--cache-dir DIR`, `--no-disk-cache`, `--memory-mb N`, `--workers N`, `--model-window N`, `--jobs-dir DIR`, `--no-jobs`, `--pigmentation-config PATH`, `--alphagenome-address URL`, `--alphagenome-ca-cert PEM`, `--alphagenome-server-dir DIR`, `--alphagenome-server-python PYTHON`, `--alphagenome-server-port PORT`, `--host`, `--port`, `--open`, `--no-add-datasets`, `--no-remote`, `--verbose`.

Without `--port` the visualizer uses 8780, or the next free port if 8780 is taken. If a visualizer already running on that port has every requested dataset (and no `--annotations`, `--runs-root`, `--gtf` or `--consensus-dataset-dir` is given), the command prints its URL, opens it with `--open`, and exits; an explicit `--port` that is busy is an error. `--open` waits until the server is listening and, on a machine without a graphical display, prints the URL instead. In an SSH session the command prints the `ssh -L` tunnel to reach it from your computer.

## Top-Level Commands

| Command | Purpose |
|---|---|
| `visualize` | Interactive visualizer for datasets, AlphaGenome tracks, sequences and experiments (see [Visualizer](../components/visualizer.md)) |
| `doctor` | Which features this machine can run (Python packages, `bcftools`/`samtools`, AlphaGenome backend and local server, PyTorch, GPU, free disk) and the command that enables each missing piece; `--json` for scripts. See [Requirements](../getting-started/requirements.md) |
| `audit-configs` | Check configs for legacy paths and active/inactive status |
| `audit-data` | Validate registered dataset paths and expected artifacts |
| `config ...` | Describe, validate, and export typed config schemas |
| `completion bash` | Print Bash completion script |
| `data ensure-1kg-vcf` | Download 1000 Genomes high-coverage chromosome VCFs into the canonical dataset layout |
| `references ensure-grch38` | Download/register the canonical full GRCh38 FASTA |
| `convert vcf-to-23andme` | Convert VCF to 23andMe raw format |
| `snp-ancestry run` | Run SNP ancestry pipeline |
| `snp-ancestry markers` | Export ranked ancestry-informative markers from computed SNP ancestry statistics |
| `snp-ancestry prune` | Positionally prune ranked ancestry-informative markers |
| `snp-ancestry train-ml` | Train sklearn ancestry classifiers from exported AIMs |
| `snp-ancestry ablate` | Retrain sklearn baselines after removing top AIMs |
| `snp-ancestry plot` | Plot ML metrics, feature importance, and AIM-ablation curves |
| `genomes-analyzer run` | Run FASTQ/BAM/CRAM/VCF operational workflow |
| `dataset-builders non-longevous ...` | Build derived 1000G/AlphaGenome datasets |
| `dataset-builders vcf-import --spec SPEC` | Build a canonical dataset from any phased VCF plus free-form sample metadata |
| `alphagenome predict-dataset DIR` | AlphaGenome predictions (any of the 11 outputs, any tissues) for every haplotype window of a canonical dataset (resumable) |
| `alphagenome catalog` | List every AlphaGenome output track and ontology term (JSON, optional per-track CSV) |
| `alphagenome server setup\|start\|check` | Run AlphaGenome on this machine's NVIDIA GPU behind the hosted API's gRPC interface (see [AlphaGenome commands](#alphagenome-commands)) |
| `alphagenome ...` | Run AlphaGenome analysis/integration utilities |
| `genotype ...` | Dense/aligned genotype predictor workflows |
| `variant ...` | Sparse variant transformer workflows |

## Examples

```bash
genomics audit-configs --fail-on-active-legacy
genomics audit-data --dataset-id 1kg_high_coverage --check-bcftools-chain --sample-limit 3 --fail-on-missing
genomics config describe genotype
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics references ensure-grch38
genomics data ensure-1kg-vcf --chrom chr15
genomics completion bash
```

`genomics audit-configs` checks active genotype and variant transformer configs for legacy dataset paths and result/cache migration status. Use `--fail-on-active-legacy` as the normal gate for new runs.

`genomics audit-data` validates registered dataset IDs. Add `--check-bcftools-chain` before aligned genotype runs that use `alignment_mapping: bcftools_chain`.

## Genotype Commands

```bash
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype split configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype test configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics genotype search configs/predictors/genotype_based/icann/search_rf_xgboost.yaml
genomics genotype search configs/predictors/genotype_based/icann/search_cnn2_ablation.yaml
genomics genotype stability configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml
genomics genotype confidence-intervals configs/predictors/genotype_based/icann/search_rf_xgboost.yaml --experiment-dir results/genotype_based_predictor/icann/search/rf_xgboost_pca300/best --split test
genomics genotype evaluate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --checkpoint best_accuracy --split test
genomics genotype pca-variance configs/predictors/genotype_based/icann/genes_1000_all_rf.yaml --output results/pca_variance.png --json-output results/pca_variance.json
genomics genotype compare-aligned-signals configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml --max-samples 100 --max-pairs 500 --output-dir results/genotype_based_predictor/analysis/aligned_signal_similarity
genomics genotype workbench --host 127.0.0.1 --port 8780
genomics genotype sync-bcftools-artifacts --source-dir /path/to/consensus --target-dir /path/to/canonical --link-mode hardlink
genomics genotype single-gene-screen configs/predictors/genotype_based/neural_legacy/pigmentation_binary_single_gene_screen.yaml --dry-run
```

`genomics genotype train` trains on the training split, validates during training, and reports final validation metrics for the `best_accuracy` checkpoint. It does not evaluate the test split. Use `genomics genotype test` after model/hyperparameter selection to evaluate `best_accuracy` on the held-out test split.

`genomics genotype split` materializes or validates the processed dataset cache, split metadata, and dataset report plots without training a model.

`genomics genotype search` runs validation-only hyperparameter search for configured sklearn baselines or named PyTorch ablation candidates. Use `genomics genotype test` on the selected best directory after model selection.

`genomics genotype stability` evaluates sklearn model stability on the development split only. It keeps the original test split fixed, resamples `train+val` using `stability_analysis.strategy` (`repeated_random_split`, `randomized_split`, or `cross_validation`), and writes aggregate validation metrics. Use the held-out test split only after model selection. Set the same `stability_analysis.split_plan_path` across model configs to reuse exactly the same resampling plan by `sample_id`.

`genomics genotype confidence-intervals` recomputes metrics and configured bootstrap confidence intervals for a saved model artifact/checkpoint without retraining.

`genomics genotype pca-variance` computes and plots sklearn PCA explained variance for the selected processed dataset/config. Use `--force` to rebuild existing outputs.

`genomics genotype compare-aligned-signals` reads the processed aligned tensor cache and compares AlphaGenome signal channels between pairs of individuals using only positions where both individuals have `valid_mask=1`. By default it analyzes the `train` split only; pass `--splits train val test` to include other splits deliberately. It writes global pairwise similarity, top absolute differences, per-position superpopulation effects (`eta_squared`, group mean delta, standardized delta), a sparse top-effect pairwise summary, and `summary.json`. Use `--max-samples` and `--max-pairs` for a fast pilot run; add `--permutations 1000` to test the global between-vs-within superpopulation MAD difference.

`genomics genotype workbench` opens the unified visualizer (same as `genomics visualize --dataset <dataset-dir> --runs-root <runs-root>`); see [Visualizer](../components/visualizer.md). The Perturbation Lab (in-silico overwrite/scramble/revert edits of a haplotype, re-predicted by AlphaGenome and re-scored by any trained model) now runs inside the visualizer; `--pigmentation-config` only chooses the run it selects first, and `--pigmentation-lab-port` is ignored. `--legacy` starts the previous multi-process workbench instead.

`genomics genotype sync-bcftools-artifacts` previews or applies hardlink/symlink/copy operations for consensus and chain artifacts required by the aligned `haplotype_channels` layout. Add `--apply` only after reviewing the preview.

`genomics genotype single-gene-screen` expands a base config into per-gene or per-ontology runs. Use `--dry-run` first to inspect generated commands and output paths.

## SNP Ancestry Commands

```bash
genomics snp-ancestry run configs/predictors/snp_ancestry/default.yaml
genomics snp-ancestry run configs/predictors/snp_ancestry/icann/gene_windows_h1_mlc.yaml
genomics snp-ancestry run configs/predictors/snp_ancestry/chr15_aims.yaml
genomics snp-ancestry markers --config configs/predictors/snp_ancestry/chr15_aims.yaml --top 500 --output results/snp_ancestry_predictor/chr15/aims_top500.tsv
genomics snp-ancestry prune --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --window-bp 50000 --output results/snp_ancestry_predictor/chr15/aims_top500_pruned_50kb.tsv
genomics snp-ancestry train-ml --config configs/predictors/snp_ancestry/chr15_aims.yaml --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --models logistic random_forest --output-dir results/snp_ancestry_predictor/chr15/ml
genomics snp-ancestry ablate --config configs/predictors/snp_ancestry/chr15_aims.yaml --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --remove-top 0 1 5 10 50 100 --output-dir results/snp_ancestry_predictor/chr15/ablation
genomics snp-ancestry plot --ml-dir results/snp_ancestry_predictor/chr15/ml --ablation-dir results/snp_ancestry_predictor/chr15/ablation --output-dir results/snp_ancestry_predictor/chr15/plots
```

`genomics snp-ancestry run` preserves the existing conversion, statistics, and prediction pipeline. The legacy `--config`/`-c` form is still accepted. `genomics snp-ancestry markers` is an optional post-processing command: it reads the statistics JSON produced by `run`, ranks markers by `fst`, `maf`, or `max_delta_frequency`, and writes an audit-friendly TSV with per-class allele frequencies. `genomics snp-ancestry prune` removes lower-ranked markers within a configured base-pair window of already kept markers, preserving the TSV columns and recalculating ranks. `genomics snp-ancestry train-ml` consumes an AIM TSV and the same per-individual 23andMe files to train sklearn `logistic` and/or `random_forest` baselines, writing metrics, predictions, feature importance, and `model.joblib` artifacts. `genomics snp-ancestry ablate` uses the same inputs, removes ranked marker prefixes such as the top 1, 5, or 10 AIMs, retrains the selected sklearn models, and writes `ablation.tsv` plus `summary.json` for measuring robustness to top-marker removal. `genomics snp-ancestry plot` turns those outputs into PNG confusion matrices, top-feature bar charts, and ablation curves.

## Variant Commands

```bash
genomics variant materialize --dataset-id 1kg_high_coverage --output-dir /dados/GENOMICS_DATA/variant_transformer/superpopulation
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
genomics variant evaluate configs/predictors/variant_transformer/repo_layout.example.yaml --checkpoint best_accuracy --split test
genomics variant analyze-counts /dados/GENOMICS_DATA/variant_transformer/superpopulation --central-window-size 32768
```

Before training variant transformer configs, verify that their processed datasets exist:

```bash
genomics audit-data --dataset-id variant_transformer_superpopulation --fail-on-missing
```

`genomics variant materialize` creates the sparse-token processed dataset from a canonical dataset, explicit VCF sources, or BED/sample metadata inputs. `genomics variant train` and `genomics variant evaluate` read that materialized dataset.

`genomics variant analyze-counts` summarizes token counts and central-window behavior for an existing processed dataset.

## AlphaGenome Commands

```bash
genomics alphagenome analyze -- -i sequence.fasta -k API_KEY -o results/
genomics alphagenome integrate -- --integrated --vcf vcf/sample.vcf.gz --ref refs/GRCh38.fa --api-key API_KEY --output integrated_analysis/
genomics alphagenome tracks --api-key API_KEY --output configs/workflows/alphagenome/tracks.json
genomics alphagenome chr15-local --config configs/workflows/alphagenome/chr15_local.yaml --sample HG00096 --outputs RNA_SEQ --max-windows 4 --haplotype H1 --strand plus --batch-size 4
```

Local AlphaGenome server (an `alphagenome_research` checkout in its own environment with CUDA `jax`):

```bash
genomics alphagenome server setup                          # clone, create env "alphagenome", install jax[cuda12], check weights
genomics alphagenome server setup --download-weights       # also download the gated Hugging Face weights (~700 MB)
genomics alphagenome server setup --python /path/to/env/bin/python --jax cuda13 --dry-run
genomics alphagenome server start                          # 127.0.0.1:50051, foreground; --host 0.0.0.0 to share, --port N
genomics alphagenome server check --predict                # checkout, env, weights, connection and a 16 kb prediction
```

`setup` takes `--dir` (checkout, default `$ALPHAGENOME_SERVER_DIR`, `../alphagenome_research` or
`~/.local/share/genomics/alphagenome_research`), `--conda-env NAME`, `--python PATH`, `--jax {cuda12,cuda13,cpu}`,
`--update` (git pull), `--reinstall`, `--dry-run`. `start` takes `--host`, `--port`, `--model-version`
(`all_folds` or `fold_0`…`fold_3`), `--checkpoint DIR` (local, e.g. Kaggle), `--plaintext`, `--allow-cpu`.
Clients reach it with `ALPHAGENOME_ADDRESS=grpc://127.0.0.1:50051` (every command that builds its client
through `genomics.core.alphagenome_connection.create_dna_client`, including `predict-dataset`, `catalog`,
`tracks`, `analyze` and the non-longevous builder's `--predict`) or from the visualizer's AlphaGenome page.

Arguments after `--` are forwarded to the underlying AlphaGenome modules. `tracks` exports output and ontology metadata used when choosing `alphagenome_outputs` and `ontology_terms`.

`chr15-local` runs the local `alphagenome_research` JAX implementation over phased 1000 Genomes haplotypes for chr15. The initial implementation uses 1,048,576 bp windows with 524,288 bp stride, applies phased SNVs only, predicts both plus and minus strands by default, and stores per-window outputs with a `window_plan.json` that records the nearest-center stitching intervals. `variants.include_indels` is present in the config but intentionally raises until coordinate remapping for indels is implemented. By default the chr15 VCF is resolved from `variants.dataset_id` under `raw_variants/vcf_chromosomes`. Use `--batch-size` for true JAX batching through the internal `model._predict` path, and increase it gradually until GPU memory is well utilized. Use `--shard-index` and `--num-shards` to split samples across workers, or `--window-shard-index` and `--num-window-shards` to split chr15 windows across workers for the same sample set. The `--ref-fasta`, `--vcf`, `--output-dir`, and `--outputs` flags override the YAML without editing it.

Before running `chr15-local` without `--ref-fasta`, register the canonical full reference FASTA:

```bash
genomics references ensure-grch38
```

This writes `${GENOMICS_DATA_ROOT:-/dados/GENOMICS_DATA}/references/GRCh38_full_analysis_set_plus_decoy_hla.fa` and indexes it with `samtools faidx` unless `--skip-index` is used.

If the canonical 1000 Genomes dataset exists but lacks raw chromosome VCFs, download chr15 into the registered layout with:

```bash
genomics data ensure-1kg-vcf --chrom chr15
```

This writes to `${GENOMICS_DATA_ROOT:-/dados/GENOMICS_DATA}/v1/1kG_high_coverage/raw_variants/vcf_chromosomes/`. Use `--chrom all` for chromosomes 1-22 and X.

Batching note: `--batch-size > 1` stacks sequences into `[B, 1048576, 4]` and calls the JAX-jitted AlphaGenome predictor once per batch. Start with `--batch-size 2`, then try `4`, `8`, etc. If you hit OOM, reduce the value. Process/window sharding remains useful across multiple GPUs or for independent resume/retry.

## Dataset Builder Commands

```bash
genomics dataset-builders non-longevous build --config configs/workflows/non_longevous_dataset/default.yaml
genomics dataset-builders non-longevous build-window -- --help
genomics dataset-builders non-longevous visualize configs/workflows/non_longevous_dataset/default.yaml
```

`build-window` forwards arguments to the window builder module. Use `-- --help` to inspect the forwarded module's options.

```bash
genomics dataset-builders vcf-import --spec import.yaml            # build (resumable; re-run to add samples or genes)
genomics dataset-builders vcf-import --spec import.yaml --inspect  # samples, contigs and phasing of the VCF
genomics alphagenome predict-dataset /data/my_cohort --outputs RNA_SEQ,CAGE --ontology CL:1000458,UBERON:0002107
genomics alphagenome predict-dataset /data/my_cohort --outputs all --ontology CL:1000458   # every output type
genomics alphagenome predict-dataset /data/my_cohort --outputs RNA_SEQ,CAGE --ontology CL:1000458 --haplotypes ref   # the reference window of each gene only
genomics alphagenome catalog --output catalog.json --csv tracks.csv                         # every track / tissue
```

`vcf-import` builds the canonical layout (reference windows, per-sample window VCFs, H1/H2 consensus FASTAs, `dataset_metadata.json`) from any phased, bgzipped VCF (or a `{chrom}` pattern) and a reference FASTA, without assuming 1000 Genomes metadata. `vcf` and `reference_fasta` may be URLs (`https://`, `ftp://`, `s3://`, `gs://`): only each window's region is fetched through the remote `.tbi`/`.csi`/`.fai` index, which is cached in `~/.cache/genomics/remote_index`; `vcf_overrides` (`{chrom: file}`) covers chromosomes whose file does not follow the `{chrom}` pattern (e.g. `...chrX...v2.vcf.gz`). The spec (JSON or YAML) takes `name`, `output_dir`, `vcf`, `reference_fasta`, `window_size` (an AlphaGenome length: 16384, 131072, 524288 or 1048576), `genes` (symbols or Ensembl ids, resolved with `gtf`, default the 1000 Genomes dataset's `gtf_cache.feather`) and/or `regions` (`{name, chrom, start, end}`), optional `samples`, and sample metadata as `metadata_file` (CSV/TSV/JSON/PLINK `.fam`; `id_column`, `family_column`, `sex_column` optional, the id column is detected otherwise) or inline `metadata` (`{sample: {field: value}}`). Every metadata column becomes a sample field (facet, group, training target). Windows are centred like `build-window`, and the outputs are byte-identical to the 1000 Genomes builder's for the same sample and window. A spec with `"extend": true` adds windows to an existing dataset (an import or the 1000 Genomes dataset) and changes only its window lists. `predict-dataset` writes `predictions_<H>/<output>.npz` + `<output>_metadata.json` with the hosted API or the server in `ALPHAGENOME_ADDRESS`, skips windows that already have the outputs, and updates the dataset metadata. Every AlphaGenome output is supported: per-base tracks (RNA-seq, CAGE, PRO-cap, DNase, ATAC, splice sites and splice-site usage), 128 bp ChIP-seq tracks (histone marks, TFs; `values` has one row per bin and the file stores `resolution`), splice junctions (`starts`/`ends`/`strands`/`values`) and 2048 bp contact maps (`values` is bins × bins × tracks). `--ontology` takes any CURIE in the AlphaGenome catalog (`alphagenome catalog`); `--all-tissues` predicts every track. `--haplotypes` takes `H1`, `H2` and/or `ref`: `ref` predicts each window's reference sequence (`references/windows/<gene>/ref.window.fa`) into `references/windows/<gene>/predictions_ref/` (same files; recorded as `reference_outputs` in the window metadata), the baseline the visualizer draws next to individuals. Both are also started from the visualizer (Overview, Jobs, Import).

Use `--help` at any level:

```bash
genomics genotype train --help
genomics alphagenome analyze -- --help
```
