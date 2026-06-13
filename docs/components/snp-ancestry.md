# SNP Ancestry

The SNP ancestry predictor implements an allele-frequency ancestry model over a configured SNP or variant panel. It first normalizes multi-sample VCF genotypes into per-individual 23andMe-like files, then estimates class-specific reference allele frequencies, and finally predicts each individual's ancestry either as a single maximum-likelihood class or as a vector of admixture proportions.

For the reproduction branch that uses classical VCF/SNP classification instead of AlphaGenome neural tensors, start with [Training And Running Predictions](../getting-started/training-and-running-predictions.md#classical-snp-superpopulation-path).

## CLI

```bash
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml
genomics snp-ancestry markers --config configs/predictors/snp_ancestry/chr15_aims.yaml --top 500 --output results/snp_ancestry_predictor/chr15/aims_top500.tsv
genomics snp-ancestry prune --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --window-bp 50000 --output results/snp_ancestry_predictor/chr15/aims_top500_pruned_50kb.tsv
genomics snp-ancestry train-ml --config configs/predictors/snp_ancestry/chr15_aims.yaml --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --output-dir results/snp_ancestry_predictor/chr15/ml
genomics snp-ancestry ablate --config configs/predictors/snp_ancestry/chr15_aims.yaml --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv --remove-top 0 1 5 10 50 100 --output-dir results/snp_ancestry_predictor/chr15/ablation
genomics snp-ancestry plot --ml-dir results/snp_ancestry_predictor/chr15/ml --ablation-dir results/snp_ancestry_predictor/chr15/ablation --output-dir results/snp_ancestry_predictor/chr15/plots
```

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/predictors/snp_ancestry/` |
| Configs | `configs/predictors/snp_ancestry/` |
| Default refs location | `results/snp_ancestry_predictor/refs/` |

## Pipeline Steps

1. Convert multi-sample VCF data into 23andMe-like per-individual files.
2. Compute reference allele-frequency statistics.
3. Predict ancestry with maximum-likelihood or admixture-style methods.

## Implementation Structure

| Module | Responsibility |
|---|---|
| `pipeline.py` | End-to-end conversion, frequency computation, prediction, reporting |
| `markers.py` | Optional AIM ranking/export from computed statistics |
| `prune.py` | Optional positional pruning for ranked AIM TSVs |
| `train_ml.py` | Optional sklearn baselines over exported AIMs |
| `ablate.py` | Optional robustness study that removes top AIMs and retrains sklearn baselines |
| `plots.py` | Optional PNG plots for ML metrics, feature importance, and ablation curves |
| `generate_regions.py` | Utility for generating region/panel files from configured windows |
| `__main__.py` | Module entrypoint used by the `genomics` CLI delegation |

The implementation reuses the VCF-to-23andMe converter for normalization, chip-panel handling, dbSNP annotation, chromosome naming, and output formatting.

## AIM Export

The existing `run` command is unchanged and remains the supported path for configs such as `configs/predictors/snp_ancestry/icann/gene_windows_h1_mlc.yaml` and `configs/predictors/snp_ancestry/icann/illumina_gsa_diploid_mlc.yaml`. Marker analysis is an opt-in post-processing step over the statistics JSON produced by `run`:

```bash
genomics snp-ancestry run --config configs/predictors/snp_ancestry/chr15_aims.yaml
genomics snp-ancestry markers \
  --config configs/predictors/snp_ancestry/chr15_aims.yaml \
  --score fst \
  --top 500 \
  --output results/snp_ancestry_predictor/chr15/aims_top500.tsv
```

The TSV contains `rsid`, genomic position, tracked allele, `maf`, the same multi-class Fst-like score used during statistics top-N selection, maximum allele-frequency delta across classes, and one `freq_<class>` column per ancestry class. This makes the AIM set auditable before downstream ML, ablation, or functional annotation.

## Positional Pruning

Ranked AIMs can contain nearby redundant SNPs. Use positional pruning to keep the highest-ranked marker in each local window:

```bash
genomics snp-ancestry prune \
  --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv \
  --window-bp 50000 \
  --output results/snp_ancestry_predictor/chr15/aims_top500_pruned_50kb.tsv
```

The command reads the AIM TSV in rank order. For each chromosome, it keeps a marker if it is more than `--window-bp` base pairs away from every already kept marker on that chromosome. Lower-ranked markers inside the window are removed, original TSV columns are preserved, and `rank` is recalculated. A small JSON summary is written next to the output TSV by default.

## ML Baselines

After exporting AIMs, train sklearn classifiers with:

```bash
genomics snp-ancestry train-ml \
  --config configs/predictors/snp_ancestry/chr15_aims.yaml \
  --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv \
  --models logistic random_forest \
  --output-dir results/snp_ancestry_predictor/chr15/ml
```

`train-ml` reuses the same split metadata, 23andMe files, tracked alleles, and `haplotype_mode` as the allele-frequency model. Missing marker doses are imputed with train-set marker means. Each model writes `metrics.json`, `predictions_<split>.tsv`, `feature_importance.tsv`, and `model.joblib` under its model-specific output directory.

## AIM Ablation

To test whether performance depends on a small number of top-ranked AIMs, remove ranked prefixes from the marker TSV and retrain the same sklearn baselines:

```bash
genomics snp-ancestry ablate \
  --config configs/predictors/snp_ancestry/chr15_aims.yaml \
  --markers results/snp_ancestry_predictor/chr15/aims_top500.tsv \
  --models logistic random_forest \
  --remove-top 0 1 5 10 50 100 \
  --output-dir results/snp_ancestry_predictor/chr15/ablation
```

The `0` condition is the baseline with no marker removal. Larger `--remove-top` values drop that many rows from the ranked AIM TSV before retraining. The command writes `ablation.tsv` with one row per model/removal condition and `summary.json` with full split metrics. By default, the headline metric comes from the test split when present, then validation, then train.

## Plots

Generate visual summaries from `train-ml` and `ablate` outputs with:

```bash
genomics snp-ancestry plot \
  --ml-dir results/snp_ancestry_predictor/chr15/ml \
  --ablation-dir results/snp_ancestry_predictor/chr15/ablation \
  --splits test val \
  --top-features 50 \
  --output-dir results/snp_ancestry_predictor/chr15/plots
```

When `--ml-dir` is supplied, the command reads each model's `metrics.json` and `feature_importance.tsv`, then writes `confusion_matrix_<split>.png` and `feature_importance_topN.png`. When `--ablation-dir` is supplied, it reads `ablation.tsv` and writes `ablation_curves.png` for the selected metrics. Plotting is deliberately separate from training so non-graphical runs do not require importing matplotlib.

## Config Sections

| Section | Purpose |
|---|---|
| `input` | Splits, individuals directory, multi-sample VCF pattern, chromosomes |
| `conversion` | 23andMe output settings, dbSNP annotation, chip-panel filtering, parallelism |
| `statistics` | Reference subsets, allele-frequency output, SNP panel, MAF/max SNP filters |
| `prediction` | Target level, derived targets, haplotype mode, prediction method, results dir |

## Mathematical Implementation

The implementation in `src/genomics/predictors/snp_ancestry/pipeline.py` uses the same statistical representation in Steps 2 and 3. Let:

- \(K\) be the number of ancestry classes at `prediction.level`.
- \(C = \{c_1, \ldots, c_K\}\) be the ordered class list stored in the statistics JSON.
- \(S\) be the set of retained rsIDs after panel, region, variant-type, MAF, and optional Fst filters.
- \(g_{i,s}\) be the observed genotype string for individual \(i\) at marker \(s\).
- \(a_s\) be the tracked allele for marker \(s\); this is the first observed allele in `alleles` encoding, or `0` in numeric `gt` encoding.
- \(d_{i,s}\) be the dose, meaning the number of copies of tracked allele \(a_s\) observed in \(g_{i,s}\).
- \(m\) be the ploidy used by the model: \(m = 1\) for `H1` or `H2`, and \(m = 2\) for `H1+H2`.

For haplotype mode `H1`, only the first allele contributes, so \(d_{i,s} \in \{0, 1\}\). For `H2`, only the second allele contributes. For `H1+H2`, both alleles contribute, so \(d_{i,s} \in \{0, 1, 2\}\). Missing calls (`--`, `.`, or malformed genotypes) are skipped.

### Step 1: VCF Normalization

Step 1 queries each configured chromosome with `bcftools`. The include expression keeps biallelic variants only:

```text
(TYPE="snp" and/or TYPE="indel") and N_ALT=1
```

When `skip_no_rsid` is enabled, the expression also requires an rsID-like identifier:

```text
ID~"^rs"
```

If `annotate_dbsnp` is enabled, `bcftools annotate -c ID` is piped into `bcftools query` so that dbSNP rsIDs are added before filtering and output. Optional chip-panel, custom SNP-panel, and region filters further restrict the emitted rows.

The code supports two genotype encodings:

- `alleles`: biallelic SNP genotypes are converted from VCF `GT` values into allele strings such as `AA`, `AG`, or `GG`. For a VCF row with `REF = r` and `ALT = b`, genotype `0/0` becomes `rr`, `0/1` or `1/0` becomes `rb` or `br` preserving allele order when available, and `1/1` becomes `bb`.
- `gt`: the stored genotype remains a biallelic code such as `0/0`, `0/1`, `1/0`, or `1/1`. This is selected by default when non-SNP variants such as indels are included.

Each individual receives a file with rows of the form:

```text
rsid    chromosome    position    genotype
```

### Step 2: Reference Frequency Estimation

Step 2 loads the configured split metadata and keeps reference individuals whose split is listed in `statistics.reference_subsets` and whose label exists at `prediction.level`. For class `c`, define the reference set:

\[
R_c = \{i : \operatorname{split}(i) \in \operatorname{reference\_subsets} \land \operatorname{label}(i) = c\}
\]

For each retained marker `s` and class `c`, only individuals with a usable genotype at `s` contribute. Define:

\[
R_{s,c} = \{i \in R_c : g_{i,s}\ \text{is usable}\}
\]

The pipeline accumulates tracked-allele counts and allele totals:

\[
x_{s,c} = \sum_{i \in R_{s,c}} d_{i,s}
\]

\[
n_{s,c} = \sum_{i \in R_{s,c}} m
\]

The estimated class-specific tracked-allele frequency is:

\[
p_{s,c} = \frac{x_{s,c}}{n_{s,c}}
\]

If a class has no usable allele observations at a marker, the implementation uses \(p_{s,c} = 0.5\) for that class at that marker.

The global tracked-allele frequency and minor allele frequency are:

\[
p_s = \frac{\sum_c x_{s,c}}{\sum_c n_{s,c}}
\]

\[
\operatorname{MAF}_s = \min(p_s, 1 - p_s)
\]

Markers with \(\operatorname{MAF}_s < \texttt{statistics.min\_maf}\) are discarded.

When `statistics.max_snps` is set and more markers remain, markers are ranked by the implemented Fst-like separation score:

\[
\bar{p}_s = \frac{1}{K} \sum_c p_{s,c}
\]

\[
F_{ST,s} = \frac{\operatorname{Var}_c(p_{s,c})}{\bar{p}_s(1 - \bar{p}_s)}
\]

where:

\[
\operatorname{Var}_c(p_{s,c}) = \frac{1}{K} \sum_c (p_{s,c} - \bar{p}_s)^2
\]

If \(\bar{p}_s\) is 0 or 1, \(F_{ST,s}\) is set to 0. The top `statistics.max_snps` markers by this score are retained.

The statistics JSON stores:

- `metadata.populations`: the class order `C`.
- `ref_alleles`: the tracked allele \(a_s\) per marker.
- `snp_info`: chromosome and position per marker.
- `allele_frequencies`: the vector \([p_{s,c_1}, \ldots, p_{s,c_K}]\) per marker.

### Step 3a: Maximum-Likelihood Classification

With `prediction.method: "mle"`, the pipeline assumes markers are conditionally independent given an ancestry class. For individual \(i\), marker \(s\), and class \(c\), the Bernoulli or binomial genotype probability is:

For single-haplotype modes (\(m = 1\)):

\[
P(d_{i,s} \mid c) = p_{s,c}^{d_{i,s}}(1 - p_{s,c})^{1 - d_{i,s}}
\]

For diploid mode (\(m = 2\)):

\[
P(d_{i,s}=2 \mid c) = p_{s,c}^2
\]

\[
P(d_{i,s}=1 \mid c) = 2p_{s,c}(1 - p_{s,c})
\]

\[
P(d_{i,s}=0 \mid c) = (1 - p_{s,c})^2
\]

Frequencies are clipped to `[1e-6, 1 - 1e-6]` before taking logarithms. The log-likelihood for class `c` is:

\[
L_i(c) = \sum_{s \in S_i} \log P(d_{i,s} \mid c)
\]

where `S_i` is the subset of retained markers with usable genotypes for individual `i`. The predicted class is:

\[
\hat{c}_i = \arg\max_c L_i(c)
\]

The reported `proportions` for MLE are not fitted admixture fractions. They are a softmax-normalized version of class log-likelihoods:

\[
q_i(c) = \frac{\exp(L_i(c) - \max_j L_i(c_j))}{\sum_j \exp(L_i(c_j) - \max_l L_i(c_l))}
\]

### Step 3b: Admixture Proportions With EM

With `prediction.method: "admixture_em"`, the pipeline estimates a simplex vector:

\[
\alpha_i = (\alpha_{i,1}, \ldots, \alpha_{i,K})
\]

\[
\sum_c \alpha_{i,c} = 1, \qquad \alpha_{i,c} \ge 0
\]

The mixed tracked-allele frequency at marker `s` is:

\[
\pi_{i,s} = \sum_c \alpha_{i,c}p_{s,c}
\]

The likelihood uses the same Bernoulli/binomial observation model as MLE, replacing \(p_{s,c}\) with \(\pi_{i,s}\). EM starts from a uniform vector \(\alpha_c = 1 / K\). At each iteration, the implementation computes:

\[
q_{i,s} = 1 - \pi_{i,s}
\]

\[
r_s = \frac{d_{i,s}}{\pi_{i,s}}
\]

\[
t_s = \frac{m - d_{i,s}}{q_{i,s}}
\]

Then for each class `c`:

\[
u_c = \sum_s p_{s,c}r_s + \sum_s (1 - p_{s,c})t_s
\]

\[
\alpha_{c}^{new} = \alpha_c u_c
\]

\[
\alpha_{c}^{new} = \frac{\alpha_{c}^{new}}{\sum_j \alpha_{j}^{new}}
\]

The loop stops when:

\[
\max_c |\alpha_{c}^{new} - \alpha_c| < \texttt{prediction.admixture\_em.tol}
\]

or when `prediction.admixture_em.max_iter` is reached. The predicted class for metrics is \(\arg\max_c \alpha_{i,c}\).

### Step 3c: Admixture Proportions With L-BFGS-B

With `prediction.method: "admixture_mle"`, the same admixture likelihood is optimized as a negative log-likelihood. To enforce simplex constraints, the implementation uses unconstrained parameters `theta` and a softmax transform:

\[
\alpha_c(\theta) = \frac{\exp(\theta_c - \max_j \theta_j)}{\sum_j \exp(\theta_j - \max_l \theta_l)}
\]

The minimized objective is:

\[
\operatorname{NLL}(\theta) = -\sum_{s \in S_i} \log P(d_{i,s} \mid \pi_{i,s}(\alpha(\theta)))
\]

The code supplies the analytical gradient to SciPy `minimize(..., method="L-BFGS-B", jac=True)`. It runs from a uniform initialization plus random Dirichlet restarts controlled by `prediction.admixture_mle.n_restarts`, and keeps the solution with the lowest negative log-likelihood.

### Metrics And Outputs

For batch evaluation, each prediction stores `sample_id`, split, true label, predicted label, correctness, number of SNPs used, and proportions. The summary metrics are computed from the confusion counts:

\[
TP_c = \operatorname{count}(true = c \land predicted = c)
\]

\[
FP_c = \operatorname{count}(true \ne c \land predicted = c)
\]

\[
FN_c = \operatorname{count}(true = c \land predicted \ne c)
\]

\[
\operatorname{precision}_c = \frac{TP_c}{TP_c + FP_c}
\]

\[
\operatorname{recall}_c = \frac{TP_c}{TP_c + FN_c}
\]

\[
F1_c = \frac{2\operatorname{precision}_c\operatorname{recall}_c}{\operatorname{precision}_c + \operatorname{recall}_c}
\]

\[
\operatorname{accuracy} = \frac{\operatorname{count}(correct)}{\operatorname{count}(total)}
\]

Weighted precision, recall, and F1 are support-weighted averages over classes.

## Prediction Methods

The configured method is selected by `prediction.method`:

| Method | Output interpretation | Implementation |
|---|---|---|
| `mle` | Single ancestry assignment plus likelihood-softmax scores | Independent-marker maximum likelihood over class allele frequencies |
| `admixture_em` | Estimated ancestry proportions | EM updates on the mixture likelihood; fast and does not require SciPy |
| `admixture_mle` | Estimated ancestry proportions | Softmax-parameterized numerical maximum likelihood with L-BFGS-B |

Haplotype handling is controlled by `prediction.haplotype_mode`. Some experiments use one haplotype (`H1` or `H2`); others combine both (`H1+H2`) to model diploid genotypes.

## Intermediate Artifacts

Typical outputs include:

- per-individual 23andMe-like files;
- allele-frequency tables by population or target class;
- selected SNP panels or filtered SNP sets;
- prediction JSON/TSV summaries;
- confusion matrices and classification metrics when labels are available.

## Notes

Some historical configs may still refer to legacy `/dados/GENOMICS_DATA/top3` paths. New configs should use canonical datasets and registered dataset IDs when possible.
