# Visualizer

`genomics visualize` starts a single local web app for exploring any canonical-layout dataset:
cohort metadata, AlphaGenome prediction tracks, haplotype sequences against the reference, and
experiment runs. It replaces the former multi-process workbench (six servers behind a proxy, shown
in iframes); `genomics genotype workbench` now opens the same app.

```bash
genomics visualize                                   # default dataset (1kg_high_coverage) and runs root
genomics visualize --dataset /path/to/dataset --open # any dataset directory
genomics visualize --dataset-id 1kg_high_coverage --dataset /other/dataset --annotations phenotypes.tsv
```

Code: `src/genomics/visualizer/` (Python stdlib HTTP server + JSON API, static single-page app in
`static/`, no build step and no new runtime dependencies).

## Pages

| Page | What it does |
|---|---|
| Overview | Dataset summary, cohort composition by any sample field, gene windows, outputs/tracks, open other datasets |
| Samples | Faceted cohort builder (counts update across facets), virtualized table, sample details, pin individuals, CSV export, save the cohort as a training `.view.json` (former View Builder) |
| Tracks | Canvas genome browser: overview strip, ruler, gene models, one panel per AlphaGenome track. Shows pinned **individuals**, cohort **group means** ± SD (or each group's difference from the cohort mean) by any facet, or a **population heatmap** with every cohort sample as a row |
| Sequence | Pinned haplotypes against the reference: bases when zoomed in (mismatches coloured, matches as dots, deletions and insertions marked), mismatch/indel density when zoomed out, variant lane and genotype table |
| Experiments | Runs table with any numeric metric (`weighted_*` included), training curves, run comparison, confusion matrix, per-class metrics, config, plots |
| Labs | Launches the Pigmentation Sequence Lab on demand (needs a trained checkpoint and `ALPHAGENOME_API_KEY`) |

Navigation: drag to pan, Ctrl/⌘+scroll to zoom, Shift+drag to zoom to a region, double-click to
zoom in, ←/→ and +/− on the focused plot. The locus box accepts `chr:start-end`, `start-end` or a
single position. Pinned samples, cohort filters and the locus are shared between pages and kept
per dataset in the browser; the URL is shareable.

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

## Legacy apps

The previous standalone apps under `src/genomics/predictors/genotype_based/apps/` are still
importable and runnable with `python -m ...`; `genomics genotype workbench --legacy` starts the old
workbench. The Pigmentation Sequence Lab keeps running as its own process (launched from Labs).
