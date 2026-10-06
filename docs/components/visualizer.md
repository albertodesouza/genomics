# Visualizer

`genomics visualize` starts a single local web app for exploring any canonical-layout dataset:
cohort metadata, AlphaGenome prediction tracks, haplotype sequences against the reference, and
experiment runs. It also imports new datasets from a VCF, runs AlphaGenome predictions, trains and
evaluates models as background jobs, and hosts the Perturbation Lab. It replaces the former multi-process workbench (six servers behind a proxy, shown
in iframes); `genomics genotype workbench` now opens the same app.

```bash
genomics visualize                                   # default dataset (1kg_high_coverage) and runs root
genomics visualize --dataset /path/to/dataset --open # any dataset directory
genomics visualize --dataset-id 1kg_high_coverage --dataset /other/dataset --annotations phenotypes.tsv
```

**Starting.** The port is 8780, or the next free one when 8780 is taken (`--port` fixes it;
a busy explicit port is an error). Running `genomics visualize` again while a visualizer is up
reuses it when it already has the requested datasets: the URL is printed (and opened with
`--open`) and the command exits, so `genomics visualize --open` is also the way to get back to a
running instance. `--open` opens the browser once the port is listening; without a graphical
display (no `DISPLAY` / `WAYLAND_DISPLAY`, e.g. over SSH) it prints the URL instead of starting a
text browser. In an SSH session the startup message includes the
`ssh -N -L 8780:localhost:8780 user@host` tunnel to run on your computer. With `--host 0.0.0.0`
the printed URL uses `localhost`.

Code: `src/genomics/visualizer/` (Python stdlib HTTP server + JSON API, static single-page app in
`static/`, no build step and no new runtime dependencies).
`tests/test_visualizer_frontend.py` opens every page in headless Chromium against a synthetic dataset
and fails on any JavaScript error; it also round-trips the figure export. It needs the `test-ui` extra
and `python3 -m playwright install chromium` (skipped otherwise), and runs in the *Visualizer tests*
GitHub workflow.

## Pages

| Page | What it does |
|---|---|
| Overview | Dataset summary, cohort composition by any sample field, gene windows, outputs/tracks, open other datasets |
| Samples | Faceted cohort builder (counts update across facets), virtualized table, sample details, pin individuals, CSV export, save the cohort as a training `.view.json` (former View Builder). *Region scalars* turn a track's signal over a region into a sample field; the *Ancestry PCA* view shows genotype PCs and matches two groups on them; see [Region scalars](#region-scalars) and [Ancestry PCA and matching](#ancestry-pca-and-matching) |
| Tracks | Canvas genome browser: overview strip, ruler, gene models, a sequence lane and one panel per AlphaGenome track. Tracks of several outputs (RNA-seq, CAGE, DNase, ChIP, …) can be shown together, up to 16 at once; drag a panel's grip (⋮⋮, or focus it and press ↑ / ↓) to reorder them (the order is remembered, new tracks are added at the bottom). Shows pinned **individuals**, cohort **group means** ± SD (or each group's difference from the cohort mean) by any facet, or a **population heatmap** with every cohort sample as a row. *Compare with* adds AlphaGenome's prediction of the **reference genome** (dashed line) and the **observed** ENCODE / FANTOM5 signal behind each track (a lane under it, with an agreement badge); see [Reference and observed data](#reference-and-observed-data). *Region scalar…* defines a per-sample scalar from the view. Track names open a track card, gene models and the gene button a gene card (see [Links to databases](#links-to-ontologies-and-databases)). The sequence lane shows the letter frequency per base over the pinned haplotypes (or the reference genome): letters scaled by frequency when zoomed in, stacked base composition per bin when zoomed out |
| Sequence | Pinned haplotypes against the reference: gene models, bases when zoomed in (mismatches coloured, matches as dots, deletions and insertions marked), mismatch/indel density when zoomed out, variant lane and genotype table |
| Variant | One site across the cohort: genotype counts and ALT frequency by any sample field, AlphaGenome's predicted signal by genotype over a region (an in-silico eQTL) and GTEx's measured eQTL of the same variant, with whether the directions agree; see [Variant page](#variant-page) |
| Perturbation Lab | Edit an individual's haplotypes (scramble, overwrite, revert to reference, custom sequence), re-predict them with AlphaGenome and re-score them with any trained model; gene models, original vs edited tracks, class means and class probabilities. A *saturation scan* slides one edit across a range and shows how much each window moves a class probability |
| Experiments | Runs table with any numeric metric (`weighted_*` included), training curves, run comparison, confusion matrix, per-class metrics, config, plots. Start a training run, evaluate a checkpoint on a split, or open a run in the Perturbation Lab |
| Jobs | Background jobs (imports, AlphaGenome predictions, training, evaluation): progress, live log, cancel. They keep running when the tab closes or the visualizer stops |
| AlphaGenome | Chooses the AlphaGenome backend (hosted API, remote server, or a server started on this machine) used by the lab and by prediction jobs |
| System | What this machine can run, feature by feature (packages, `bcftools`/`samtools`, AlphaGenome backend, PyTorch), hardware, and free disk space where the visualizer writes; the same report as `genomics doctor`, with the command that enables each missing piece |

Navigation: drag to pan, Ctrl/⌘+scroll to zoom, Shift+drag to zoom to a region, double-click to
zoom in, ←/→ and +/− on the focused plot. The locus box accepts `chr:start-end`, `start-end` or a
single position. Pinned samples, cohort filters and the locus are shared between pages and kept
per dataset in the browser; the URL is shareable.

**Figure export.** *Export* on the Tracks, Sequence and Perturbation Lab pages downloads the current
view as a **PNG** (2× or 4× pixels) or an **SVG**. The SVG is vector: the plots are re-drawn into SVG
paths and text (not embedded bitmaps; only the population heatmap stays a bitmap), so it can be
edited in Inkscape or Illustrator. The figure has a title and caption with the gene, locus, outputs,
what is shown (individuals, group means by a field, heatmap; reference / observed lanes), the cohort
filters (or the Perturbation Lab's model, individual and edits), the dataset and the date. *Light
background* (on by default) draws the figure with the light theme whatever theme the page uses.
The overview strip, buttons and hints are left out.

## Reference and observed data

**Reference genome.** AlphaGenome can predict each window's reference sequence (the haplotype
`ref` of a prediction job: *Reference genome* in the prediction form, *Predict the reference…* on
the Tracks page, or `genomics alphagenome predict-dataset DIR --haplotypes ref`). It is stored in
`references/windows/<gene>/predictions_ref/` and drawn as a dashed dark line in individuals and
group-mean views: the baseline every haplotype deviates from. Its tracks are matched to the
dataset's by ontology term, strand and assay, so it may be predicted with more tissues; tracks it
lacks are left empty. It is identical in genomic and haplotype coordinates and follows the
training axis (insertion columns are gaps). On TYRP1 the reference prediction correlates 0.98–0.999
per track with an individual's in genomic coordinates.

**Observed data.** *Observed data (ENCODE / FANTOM5)* draws, under each track, the experimental
signal of the experiments AlphaGenome was trained on, found from the track's metadata
(`genomics.visualizer.observed`):

| Output | Source | Matching | Signal |
|---|---|---|---|
| CAGE | FANTOM5 hg38 data hub | libraries whose sample class `derives_from` the track's CL / UBERON term in the FANTOM5 ontology; cell lines (EFO) by name | CTSS TPM (or read counts), per strand |
| RNA-seq, DNase, ATAC, histone / TF ChIP-seq, PRO-cap | ENCODE portal | experiments with the same biosample term, assay title and target | one GRCh38 bigWig per experiment: the strand's signal of unique reads (RNA), read-depth normalized signal (DNase), fold change over control (ChIP / ATAC) |
| GTEx RNA-seq, splicing, contact maps | – | reported as unavailable (no per-base public coverage) | – |

At most 4 libraries / experiments are averaged (untreated samples first). Values are read over the
window with a small bigWig reader that uses HTTP range requests (a 524 kb window costs a few
requests, ~5–10 s the first time), then cached per file and window under `<cache>/observed`.
Observed values are in the source's units, not AlphaGenome's: compare shapes and peaks (each lane
has its own scale). Coverage of the AlphaGenome catalog: 524 of 546 CAGE tracks have FANTOM5
libraries. On TYRP1 the observed FANTOM5 melanocyte TSS peaks (+0, +67, +70 bp from the MANE TSS)
are exactly where AlphaGenome predicts them. Hovering a lane shows its value; clicking its name
lists the experiments with links.

**Agreement badge.** Each observed lane carries a badge comparing it with the predicted track over the
bins in view: Pearson r of the values and of log(1 + x), the share of the observed peak bins (top
decile) that are also predicted peak bins, and the ratio of the means (predicted / observed). Units
differ between sources, so the ratio is only meaningful, and only flagged, for FANTOM5 CAGE in TPM
(beyond 4× either way); the badge is also flagged when log r < 0.3. Hover it for the details. On
TYRP1, melanocyte RNA-seq agrees well (r 0.84, 96% of the observed peak bins shared), while melanocyte
CAGE is predicted about 4.5× higher than FANTOM5 measures. The statistics depend on the view: zoom to
a gene or exon to judge it locally.

## Region scalars

A region scalar is one number per sample: the mean of one AlphaGenome track over a region, averaged
over H1 and H2, optionally summed over the region's length and log2(x + 1) transformed. Define one
with *Region scalar…* on the Tracks page (prefilled with the view's window, output, first track and
locus) or *New…* under *Region scalars* in the Samples sidebar. The region is a gene's TSS ± a
flank (500 bp by default), a gene body, or a custom genomic range. The scalar becomes a numeric sample field (a
Samples column, exported with the CSV) and, with bins, a categorical field `<name>_bin` with quantile
bins `Q1` (lowest) … `Qk`, so it can filter the cohort and split group means on Tracks. The Samples
sidebar shows each scalar's histogram and bin edges.

Values come from the per-haplotype region means of the Variant page (computed for all tracks of an
output at once and cached under `<cache>/region_means`); for melanocyte RNA-seq over the TYRP1 gene
body that takes about 20 s for the 3,202 samples of the 1000 Genomes dataset, and other tracks of the
same output and region are then instant. Definitions are kept per dataset path in
`~/.config/genomics/visualizer_scalars.json` and re-applied from the cache whenever the samples are
listed. A scalar needs every sample's prediction of that output: per-sample CAGE predictions, for
example, must be made before CAGE scalars are useful.

## Ancestry PCA and matching

*Ancestry PCA* on the Samples page (the *Table / Ancestry PCA* switch) shows a principal component
analysis of the cohort's genotypes, as a scatter of any two PCs coloured by any sample field (the
current cohort in front, the other samples dimmed; click a point to pin or unpin the sample).

**How it is computed** (`genomics.visualizer.ancestry`). Sites are biallelic SNVs with minor allele
frequency ≥ 5% over all samples, thinned to one per 2,000 bp (a crude stand-in for LD pruning), from
the cohort genotype matrices of the chosen windows (see [Variant page](#variant-page)). Dosages are
standardised per site and the scores come from the eigendecomposition of the samples' genetic
relationship matrix, as in EIGENSOFT; 10 components by default (up to 20). Only samples genotyped in
every chosen window get PCs. By default the windows covering the most samples are used, preferring
control windows (windows not listed in the dataset metadata's genes). *Windows & filters…* changes
the windows, MAF, spacing and number of components. Results are cached under `<cache>/ancestry_pca`.
*Add PC fields* adds `pc1` … `pcK` as numeric sample fields.

**Matching.** *Match two groups* pairs each sample of group A with its nearest unused sample of group
B on the first k PCs (greedy 1:1 nearest neighbour, among the samples of the current cohort, within a
caliper in pooled SDs of the PC space; 0 = no limit). The table shows the standardised mean difference
of each PC before and after matching. The matched samples become a categorical sample field (named by
you; group A's and B's labels, empty for unmatched samples) for *Use as cohort* or *Group means in
Tracks →*, so a comparison of the groups is not driven by ancestry differences between them.

In the 1000 Genomes dataset the 33 control windows have VCFs for only 1,072 of the 3,202 samples (AFR
and EUR), so the default PCA uses the 11 pigmentation windows, which cover every sample. On it, the
pigmentation labels' *weak* (EUR) and *strong* (AFR) groups differ by 10.4 pooled SDs on PC1 and no
pair falls within the default caliper: the label is population membership, and no ancestry-matched
comparison exists within this cohort. AMR vs EUR, by contrast, gives 79 pairs with every PC balanced
(|SMD| < 0.12).

## Links to ontologies and databases

Technical terms link to their databases (all links open in a new tab):

- **Genes** (gene button on Tracks, Sequence and Perturbation Lab, gene models, Overview, sample
  windows, gene pickers) open a gene card: HGNC approved name, locus type, cytoband, aliases,
  previous symbols and gene groups; links to HGNC, Ensembl, NCBI Gene, UCSC, UniProt, QuickGO and
  AmiGO, GTEx, Human Protein Atlas, gnomAD, OMIM, Open Targets, ClinVar and GeneCards; the
  dataset windows of the gene; and its **Gene Ontology** annotations (QuickGO) by aspect with
  evidence codes, each term linked to QuickGO and AmiGO.
- **Tracks** (track names on Tracks, rows of the Overview track table) open a track card: the
  biosample term (CL / UBERON / EFO) with its OLS definition and synonyms and links to OLS,
  Ontobee, CELLxGENE CellGuide and the ENCODE biosample page; the assay's OBI / EFO term (e.g.
  CAGE → OBI:0001674, DNase-seq → OBI:0001853); the ChIP target (ENCODE target page, gene card for
  TFs); the data source; and the observed experiments (ENCODE experiments and files, FANTOM5
  SSTAR sample pages).
- **Ontology CURIEs** in the tissue picker and track tables link to OLS; 1000 Genomes sample ids
  and population codes to the IGSR data portal; rsIDs of regulatory presets to dbSNP.

**Gene lists come from HGNC and the Gene Ontology.** Gene searches (prediction form, Import)
match HGNC approved symbols, previous symbols, aliases, names, HGNC and Ensembl ids (`OCA3` finds
TYRP1) and take coordinates from the GENCODE table; gene pickers label windows with HGNC names.
*Genes of a Gene Ontology term or HGNC gene group* lists the human genes annotated to a GO term or
its descendants (reviewed UniProt entries, QuickGO; optionally genes that regulate it) or the
members of an HGNC gene group, to add as windows.

Lookups (HGNC complete set, ~17 MB once; QuickGO; OLS; ENCODE; FANTOM5) are cached under
`<cache>/remote` and still served when a service is unreachable. `--no-remote` disables them
(only cached answers are used; gene searches fall back to GENCODE names).

## Coordinate systems

AlphaGenome predicts on each haplotype's own consensus sequence, so position *i* of a prediction is
position *i* of that haplotype, not of the reference. The Tracks and Sequence pages offer:

| Coordinates | Meaning |
|---|---|
| Genomic (default) | Reference positions shared by every sample. Each haplotype is remapped through the indels in its phased window VCF, using `bcftools consensus`' overlap rules; deleted reference bases are gaps. Needs only `ref.window.fa` and the per-sample window VCF (no bcftools). Validated against the bcftools_chain training alignment: identical on 2,230 haplotypes across 40+ genes |
| Haplotype | Raw prediction index (positions drift after indels) |
| Training axis | The bcftools_chain expanded alignment used to build CNN tensors (shared insertion columns). Uses the same caches as training; uncached entries are built with bcftools on first use |

The CNN model window is drawn as a blue band in every view.

## Performance

- Prediction arrays, FASTA sequences, VCF indel events and alignment entries are kept in byte-bounded
  LRU caches (`--memory-mb`, default 4096). `.npz` members are inflated with one GIL-free zlib call
  and BGZF VCF blocks are inflated directly, so loading parallelises across threads.
- Coordinate remapping and binning (NaN-aware mean/min/max per pixel) are vectorised; arrays travel
  as base64 float32, and responses are gzip-compressed.
- Cohort-wide work (group means, population heatmaps, gene annotations, training axis) runs as
  background jobs with progress and a Cancel button. Group means are computed once per
  gene/output/tracks/haplotypes/cohort and cached on disk (`results/cache/visualizer`, or
  `--cache-dir`), as are per-gene indel indexes, so later views are instant.
- Frontend: canvas rendering, requests debounced and cancelled when superseded, data fetched with
  margins so panning redraws immediately.

On the 1000 Genomes dataset (3,202 samples, 524 kb windows), 12 samples × 2 haplotypes over the
32 kb model window: the legacy track viewer needed 0.4–0.5 s per request for one track on every
pan/zoom; the visualizer returns all 6 tracks in 0.1–0.3 s on first load and ~25 ms afterwards.
Group means over the whole cohort (6,404 haplotypes) take ~50 s the first time and ~0.1 s afterwards
(also across restarts).

## Using other datasets

Any directory with `dataset_metadata.json` and the canonical layout works
(`references/windows/<gene>/ref.window.fa`, `individuals/<sample>/windows/<gene>/...`,
`predictions_<H>/<output>.npz` with optional `<output>_metadata.json`). Genes are discovered from
`references/windows/`, outputs, haplotypes and tracks from the prediction files, and sample facets
from whatever fields `individuals_pedigree` carries. Add more facets with
`--annotations table.tsv` (first column or `sample_id` = sample). Gene models are read from
`<dataset>/gtf_cache.feather` (or `--gtf`) when pandas/pyarrow are installed. Datasets can also be
opened from the Overview page while the server runs (disable with `--no-add-datasets`).

## Importing a dataset

**Import dataset** (Overview → Datasets, or Jobs) builds a canonical-layout dataset from any
phased, bgzipped VCF, so a cohort that has nothing to do with 1000 Genomes can use every page:

1. **Variants**: a VCF (or a per-chromosome pattern with `{chrom}`) and the reference FASTA it was
   called against. *Read samples* lists the VCF's samples and contigs and warns about unphased
   genotypes (H1/H2 then follow allele order). `chr15` vs `15` naming is handled.
2. **Sample metadata** (optional, any columns): upload a CSV/TSV/JSON/PLINK `.fam` file, point to
   one on the server, paste a table, or type it. The sample-id column is detected as the one that
   matches the most VCF samples; a family column (keeps relatives in one split) and a sex column
   can be mapped. Unmatched ids are listed both ways. The editable grid lets you add columns, fill a
   column by pasting from a spreadsheet, set a value for selected rows, or derive a column from the
   sample ids with a regular expression. Choose to import every VCF sample or only those with
   metadata. Every column becomes a sample field: a facet, a group for group means, a training target.
3. **Windows**: gene symbols / Ensembl ids (resolved with the GENCODE table, default the 1000
   Genomes dataset's `gtf_cache.feather`) and/or `NAME=chr:start[-end]` regions; one AlphaGenome
   window (16 kb to 1 Mb, default 524 kb) is centred on each, exactly like the 1000 Genomes builder.
4. **Dataset** name and output directory; optionally **AlphaGenome** predictions right afterwards.

The import runs as a background job (`genomics dataset-builders vcf-import --spec`), writes the
same files as the 1000 Genomes builder (byte-identical for the same sample and window), keeps each
window's cohort VCF in `references/windows/<gene>/cohort.window.vcf.gz` (used by the training
aligner) and is resumable: re-running adds samples or genes. When it finishes the dataset is
opened in the visualizer and remembered (`~/.config/genomics/visualizer_datasets.json`; *Remove*
on the Overview forgets it, files are kept).

### Quickstart with a public cohort

The Import page opens with a **Quickstart** card for two public phased GRCh38 cohorts, so the whole
workflow (import → AlphaGenome → training) can be tried without preparing any file:

| Cohort | Samples | Source |
|---|---|---|
| 1000 Genomes high coverage (Byrska-Bishop et al. 2022) | 3,202 · 26 populations | EBI FTP (or the local panel recorded by the canonical dataset, when present) |
| gnomAD HGDP + 1000 Genomes, SHAPEIT5 (Koenig et al. 2023) | 4,091 · 78 populations (925 HGDP · 52) | gnomAD's public Google Cloud bucket |

Choose how many unrelated samples to take per population (or every sample) and, for gnomAD,
*HGDP only* (no overlap with 1000 Genomes) or both projects. **Import now** starts the import with the
pigmentation genes of the canonical dataset (OCA2, HERC2, SLC24A5, SLC45A2, TYR, MC1R); **Review in
the form** fills in steps 1–4 instead, to change samples, genes or add AlphaGenome predictions.
Sample metadata is downloaded once and flattened (population, region, sex, family, coordinates;
cached under the visualizer cache's `quickstart/`). Variants are not downloaded in bulk: bcftools reads
each window's region from the remote file. 1 sample per population × 6 windows takes about a minute
and ~0.7 GB.

Remembered datasets whose directory has gone (moved, deleted, or a cleaned-up temporary directory)
are listed as *missing* in Overview → Datasets, with a button to forget them.

## Background jobs

Imports, AlphaGenome predictions, training and evaluation run as separate processes started in
their own session by `python -m genomics.visualizer.task_runner`, so they survive closing the
browser and stopping the visualizer; a restarted visualizer lists them again (a job whose process
disappeared, e.g. after a reboot, shows as *lost*). Each job has a folder under `--jobs-dir`
(default `results/visualizer/jobs`) with `task.json` (the exact commands, so any job can be re-run
from a shell), `state.json` and `log.txt`. Training/evaluation jobs share a `gpu` queue and
prediction jobs an `alphagenome` queue: one at a time each, in order. `--no-jobs` disables starting
jobs from the UI. The top bar shows running jobs; a toast reports when one finishes.

**AlphaGenome predictions** (Overview, Jobs): every AlphaGenome output (RNA-seq, CAGE, PRO-cap,
DNase, ATAC, histone and TF ChIP-seq in 128 bp bins, splice sites, splice-site usage, splice
junctions and contact maps), tissues/cell types, samples (all, the Samples-page cohort, or pinned),
windows and haplotypes. The tissue picker lists every ontology term of the AlphaGenome track
catalog (read once from the selected backend and cached as `alphagenome_catalog.json` in the cache
directory) with its track counts per output, searchable by name, CURIE, TF or histone mark and
filterable by biosample type. Only terms with tracks for the chosen outputs are listed, and
*Select all* takes every listed term (search and type filters apply). Windows can be any
existing window, any GENCODE gene (searched in the dataset's `gtf_cache.feather`), custom
`NAME=chr:start-end` regions, or curated regulatory elements that are not genes (pigmentation
enhancers near OCA2, KITLG, IRF4, BNC2 and TMEM138/DDB1, plus classic enhancers such as the LCT,
FTO/IRX3, MYC 8q24, BCL11A, SORT1 and 9p21 loci; see
`genomics.workflows.alphagenome.regulatory_regions`). New windows are first built for every
sample from the dataset's source VCF (`vcf-import` with `"extend": true`), then predicted. The form
estimates tracks and disk use. Windows that already have the outputs are skipped unless
*overwrite* is set. Uses the selected AlphaGenome backend (`genomics alphagenome predict-dataset`).
ChIP-seq tracks are drawn in the track browser as 128 bp steps; contact maps and junctions are
stored and listed on the Overview but not drawn. Training uses per-base outputs only.

**Training** (Experiments, Overview, Jobs): start from any config under
`configs/predictors/genotype_based` or an existing run, then pick the sample field to predict and
the classes (values can be grouped into classes, which becomes a derived target), the samples, the
AlphaGenome output and tissue tracks, the genes, window, features, model, epochs and splits. The
CNN kernel heights follow the number of tracks per gene. *Preview config* shows the YAML; runs are
written to `<runs root>/<run name>/` and appear on the Experiments page with their curves; the best
checkpoint is evaluated on the test split afterwards. **Evaluate…** on a run evaluates any
checkpoint on any split (`<split>_<checkpoint>_results.json`).

## Variant page

Open a variant from the Sequence page (*Variant →* in the variants table, or Shift+click a variant
tick) or type it on the Variant page: an rsID (resolved through GTEx), `chr15:28120472 A>G`,
`chr15_28120472_A_G_b38` or a position (lists the sites there). Without a variant the page lists the
window's common sites (ALT frequency ≥ 5%).

| Card | Shows |
|---|---|
| Header | Position, type, rsID (dbSNP), links to gnomAD and the GTEx variant page; genotype counts and ALT frequency in the cohort and by a sample field (superpopulation by default) |
| Direction: AlphaGenome vs GTEx | Per GTEx tissue: GTEx NES, the AlphaGenome slope and a verdict (*agree* / *disagree* when both have p < 0.05, else which side is not significant) |
| AlphaGenome by genotype | Each sample's predicted signal (mean of H1 and H2 over the region) by genotype as box and points; the OLS slope per ALT allele with SE, t and p; ALT vs REF haplotypes as log2 fold change; AlphaGenome on the reference genome (dashed) |
| GTEx v8 eQTL | For the target gene and chosen tissues (sun-exposed and not sun-exposed skin by default): GTEx's association (NES, p) with the donors' normalized expression by genotype, and the variant's significant eQTLs in every tissue |

**Choices.** *Target gene*: an annotated gene of the window (the window's gene, else the one
containing the variant, else the nearest). *Region*: the target's gene body, its promoter (TSS ± 1 kb),
the variant ± 1 kb or a custom range. *Track*: by default the track with the highest cohort mean over
the region on the target gene's strand. *Tracks by genotype →* opens the Tracks page with group means
by genotype at the site (group field `variant:<pos>:<ref>:<alt>`).

**How it is computed.** Genotypes of every sample at every site of a window come from the samples'
phased window VCFs (a sample with a VCF but no record at a site is homozygous REF; a sample without a
VCF for the window is not genotyped and left out of the counts; multi-allelic records are split per ALT
allele). They are read once per window as a background job, in parallel worker processes (3–8 s per
window for the 3,202 samples of the 1000 Genomes dataset), and cached on disk (`<cache>/cohort_genotypes`). The mean
signal of every haplotype over a region is computed for all tracks of an output at once (about 80 s
for 6,404 haplotypes) and cached (`<cache>/region_means`), so switching track, cohort filters or site
afterwards is instant. The AlphaGenome slope is an association across the cohort, like an eQTL: it
includes the effect of variants in linkage disequilibrium with this one. GTEx lookups use the GTEx
portal API v2 (dataset `gtex_v8`: GRCh38, GENCODE v26; only variants with MAF ≥ 1% among GTEx donors
were tested) and are cached under `<cache>/remote` (`--no-remote` uses cached answers only). GTEx's NES
is the effect of the ALT allele on normalized expression: compare signs, not magnitudes.

On rs12913832 (HERC2 intron, the OCA2 enhancer; G is the blue-eye allele) GTEx skin shows lower OCA2
with G (NES −0.13 sun-exposed, −0.17 not sun-exposed), while AlphaGenome's melanocyte RNA-seq over
the OCA2 gene body shows no association (slope p = 0.89 over 3,202 samples).

## Perturbation Lab

Pick a trained run (NN, CNN or CNN2; any target, layout and feature mode) and load it; the lab
scores individuals through the run's own window processing and saved normalization, so the
baseline equals what the model saw in training. Choose an individual (true class and split are
shown) and a gene, drag on the sequence lane to select a region, add edits for H1, H2 or both:

| Edit | Effect |
|---|---|
| Scramble | Shuffle the region's bases (seeded) |
| Overwrite | Every base becomes A, C, G, T or N |
| Revert | Bases aligned to the reference are reverted (removes SNVs; indels stay) |
| Sequence | Replace the region with your own sequence of the same length |

Selections are on the genomic axis and are applied to each haplotype through its own indels; all
edits keep the haplotype length, so edited predictions line up with the originals.
**Run AlphaGenome + model** re-predicts each edited haplotype (the model's input output plus the
displayed output, with the tissues of the stored tracks, reordered into the stored track order),
rebuilds the model input and shows every class probability before and after. Tracks show original
and edited signal (or their difference), optionally with each class's mean over the cohort; the
sequence lane shows reference, original and edited bases. AlphaGenome results are cached by
sequence (`<cache>/perturb`). Re-predicting an unedited haplotype reproduces its stored tracks.

**Saturation scan.** The *Saturation scan* section slides one edit (scramble, overwrite or revert to
the reference) across a range (the model window, the view or the selection) in windows of a chosen
size every *step* bp (step = window: tiled; smaller: overlapping), on H1, H2 or both, and re-scores
the model after each window. It is a background job (`POST /api/perturb/scan`, at most 256 windows);
each window costs one AlphaGenome call per edited haplotype, and windows where the edit changes no
base (reverting a stretch without variants) are scored as unchanged without a call. The result is a
lane above the tracks with the change in one class probability per window (the individual's true
class by default; green up, red down, grey ticks for unchanged windows; hover for every class) and a
list of the windows with the largest changes: *Show* zooms to one, *+* adds it as an edit so *Run*
shows its effect on the tracks. *Download TSV* saves every window's genomic range, changed bases,
probabilities and changes.

## AlphaGenome backend

The Perturbation Lab and prediction jobs call AlphaGenome. The **AlphaGenome** page selects where
those calls go; the choice is saved to `~/.config/genomics/visualizer_alphagenome.json`, used
in-process by the lab and passed to jobs as `ALPHAGENOME_ADDRESS` / `ALPHAGENOME_TLS_CA_CERT`. A
job keeps the backend it was started with. If the visualizer exits while prediction jobs use a
server it started on this machine, that server is left running for them.

| Mode | Uses |
|---|---|
| Hosted API | Google's AlphaGenome API with `ALPHAGENOME_API_KEY` (env or `~/.env`) |
| Remote server | A self-hosted `alphagenome_research` `server.py`, e.g. `grpc://10.0.0.5:50051` (plaintext), `grpcs://host:50051` (TLS) or `host:port` (TLS detected). For a self-signed TLS server set the CA certificate (`certs/ca.crt` from `scripts/generate_certs.sh`). No API key is needed |
| This machine | **Start server** serves the model on this machine's NVIDIA GPU from an `alphagenome_research` checkout in its own Python environment, shows its state (loading model → ready), packages and log, and stops it when the visualizer exits (unless background prediction jobs still use it). It listens on `127.0.0.1:50051` (`$ALPHAGENOME_SERVER_HOST`, `$ALPHAGENOME_SERVER_PORT` or `--alphagenome-server-port`) and uses TLS when `certs/server.crt` and `certs/server.key` exist in the checkout. A server already listening there (e.g. `genomics alphagenome server start` in a terminal) is detected and used |

**Test connection** connects and calls the service (`GetMetadata`); **Test prediction** also runs a
16 kb `predict_sequence`, which proves the model runs (the first prediction on a fresh server
includes JIT compilation).

**Setting up "This machine".** Run once:

```bash
genomics alphagenome server setup        # --dry-run shows the commands; --download-weights fetches the weights
```

It clones [FeLiPeOLi7/alphagenome_research](https://github.com/FeLiPeOLi7/alphagenome_research) (a fork of
Google DeepMind's research code whose `server.py` speaks the hosted API's gRPC protocol) next to this
repository, or under `~/.local/share/genomics`, creates the conda env `alphagenome` (a venv in the
checkout without conda), installs it with a CUDA-enabled `jax` whose plugin matches `jaxlib`, and checks
the Hugging Face weights (gated: accept the terms at
<https://huggingface.co/google/alphagenome-all-folds> and `hf auth login`). `genomics alphagenome server
check [--predict]` reports the checkout, environment, weights and a running server. The page lists what
is missing ("Cannot start: …") and the setup command.

The server process is `genomics/workflows/alphagenome/model_server.py`, run by the server environment's
interpreter (it does not need `genomics` installed there). It reuses the checkout's `AlphaGenomeServer`
servicer and adds what running `server.py` directly lacks: `GetMetadata` returns the model's track
metadata (`server.py` returns an empty message, which left the track catalog and tissue picker empty), the
bind address and port are options (`server.py` always listens on `0.0.0.0:50051`), cached weights load
without network or an interactive Hugging Face login, XLA does not preallocate 75% of GPU memory, and
CPU-only `jax` is refused up front. Predictions match the hosted API's (Pearson ≥ 0.9998 per track on an
OCA2 haplotype). Memory and timings: [Requirements](../getting-started/requirements.md#alphagenome-on-your-own-gpu).

The checkout is found with `--alphagenome-server-dir` (or `$ALPHAGENOME_SERVER_DIR`, then
`../alphagenome_research` next to this repository, then `~/.local/share/genomics/alphagenome_research`)
and runs with `--alphagenome-server-python` (or `$ALPHAGENOME_SERVER_PYTHON`; by default the checkout's
`.venv`, the conda env `alphagenome`, or the first other conda env that has `alphagenome_research`,
`alphagenome>=0.7` and a CUDA-enabled `jax`). Environments are inspected without importing JAX.
`--alphagenome-address URL` / `--alphagenome-ca-cert PEM` select a remote server from the command
line, overriding the saved setting.

Outside the visualizer, `genomics.core.alphagenome_connection.create_dna_client()` honours the same
environment variables, so any code that builds its client through it can use a self-hosted server.

## Legacy apps

The previous standalone apps under `src/genomics/predictors/genotype_based/apps/` are still
importable and runnable with `python -m ...`; `genomics genotype workbench --legacy` starts the old
workbench. The standalone Pigmentation Sequence Lab (`python -m genomics.predictors.genotype_based.apps.pigmentation_sequence_lab`) still works; inside the visualizer it is replaced by the Perturbation Lab.
