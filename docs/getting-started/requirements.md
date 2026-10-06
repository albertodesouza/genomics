# Requirements

What you need to run the visualizer and each of its features, and how much disk, memory and compute
datasets take. To see what **this** machine already has, run:

```bash
genomics doctor
```

It checks every feature below (packages, `bcftools`/`samtools`, the AlphaGenome backend, the GPU and
free disk space) and prints the command that enables whatever is missing. The visualizer shows the same
report on its **System** page.

## Features at a glance

| Feature | Python packages | Other software | Hardware | Network |
|---|---|---|---|---|
| **Browse** datasets, tracks, sequences, experiments, jobs (`genomics visualize`) | base install (`numpy`, `PyYAML`) | – | any machine; 4 GB RAM for the default cache (`--memory-mb`) | none (gene and ontology cards look up public databases when online; `--no-remote` turns that off) |
| **Gene models** on Tracks/Sequence, gene search | `[visualizer]` (`pandas`, `pyarrow`) | – | – | – |
| **Start a training run** from the visualizer (building and validating its config) | `[visualizer]` (`pydantic`) | – | – | – |
| **Import a dataset** from a phased VCF | `[visualizer]` | `bcftools`, `samtools` (htslib) | – | only when the VCF/FASTA are URLs (regions are streamed) |
| **AlphaGenome predictions** (prediction jobs, reference tracks, Perturbation Lab, track catalog) | `[visualizer]` (AlphaGenome client, Python ≥ 3.10) | one backend: hosted API key, a remote server, or a [local server](#alphagenome-on-your-own-gpu) | none for the hosted API; an NVIDIA GPU for a local server | hosted API or the server's address |
| **Train / evaluate** models, Perturbation Lab scoring | `[genotype]` (`torch`, `scikit-learn`, `scipy`, `pydantic`, …) | `bcftools` for the training-axis alignment | an NVIDIA GPU is strongly recommended; 16 GB+ RAM | – |

Install commands for each column are in [Installation](installation.md). In short:

```bash
scripts/env/install.sh                        # conda env with everything the visualizer can use
scripts/env/install.sh --training             # + training and Perturbation Lab scoring
scripts/env/install.sh --alphagenome-server   # + AlphaGenome on this machine's GPU
```

## AlphaGenome backends

Predictions need exactly one of these; the visualizer's **AlphaGenome** page switches between them.

| Backend | What you need | Speed (524 kb window, 3 outputs, 1 tissue) | Notes |
|---|---|---|---|
| **Hosted API** (Google) | `ALPHAGENOME_API_KEY` in the environment or in `~/.env` (free for non-commercial use: <https://deepmind.google.com/science/alphagenome>) | ~4 s per call | subject to the API's quotas; sequences are sent to Google |
| **Remote server** | the address of a machine running the local server below, e.g. `grpc://10.0.0.5:50051` | same as local, plus network | no API key needed |
| **This machine** | an NVIDIA GPU, the server environment and the weights (`genomics alphagenome server setup`) | ~0.8 s per call once warm | private; no quotas |

A local server returns the same tracks as the hosted API: on an OCA2 haplotype (HepG2 RNA-seq, CAGE,
DNase) every track correlates ≥ 0.9998 with the hosted prediction (differences come from bfloat16
arithmetic).

### AlphaGenome on your own GPU

The model runs from [`alphagenome_research`](https://github.com/FeLiPeOLi7/alphagenome_research) (a fork of
Google DeepMind's research code with a gRPC server that speaks the hosted API's protocol) in its **own**
Python environment, because it needs JAX with CUDA, TensorFlow and `alphagenome>=0.7`.

| Need | Detail |
|---|---|
| GPU | NVIDIA with a CUDA 12 driver (≥ 525). GPU memory held by the server (measured on a DGX Spark, GB10): ~2 GB idle, ~4 GB after 16 kb windows, ~8 GB after 131 kb, ~25 GB after 524 kb, ~58 GB after 1 Mb; it is kept until the server stops. On unified-memory machines (DGX Spark, GH200) that is RAM the visualizer and training also use. Google recommends an H100. CPU-only JAX is refused (`--allow-cpu` forces it; minutes per prediction) |
| Disk | ~7.5 GB for the server environment (JAX, CUDA libraries, TensorFlow) + 0.7 GB of weights (`~/.cache/huggingface`) |
| Weights | gated on Hugging Face: accept the terms at <https://huggingface.co/google/alphagenome-all-folds>, then `hf auth login` once (or a Kaggle checkpoint with `--checkpoint DIR`) |
| Network at start-up | the model reads GENCODE annotation tables and the hg38 FASTA index from Google Cloud Storage each time it starts |
| Time | loading the model: ~1.5 min. The first prediction for each new window length / output set compiles the model: ~15 s (16 kb), ~25 s (131 kb), ~1 min (524 kb), ~2 min (1 Mb). Afterwards, per call: 0.05 s (16 kb), 0.2 s (131 kb), 0.8 s (524 kb), 1.9 s (1 Mb) |

```bash
genomics alphagenome server setup     # clone, create the env (conda "alphagenome" or a venv), CUDA jax, check weights
genomics alphagenome server start     # serve on 127.0.0.1:50051 (Ctrl+C stops it)
genomics alphagenome server check --predict
```

Or click **Start server** on the visualizer's AlphaGenome page (*This machine*). The server listens on
127.0.0.1 only; `--host 0.0.0.0` (or `ALPHAGENOME_SERVER_HOST=0.0.0.0` for the visualizer) shares it with
other machines, without authentication. See [Visualizer → AlphaGenome backend](../components/visualizer.md#alphagenome-backend).

## Storage

Datasets use the [canonical layout](../concepts/data-and-results.md): one directory per sample and window.
Measured on the 1000 Genomes dataset (524 kb windows):

| Item | Size |
|---|---|
| Sequences of one sample-window (H1 + H2, raw and fixed FASTA) | 2.1 MB |
| Variants of one sample-window (window VCF + consensus-ready VCF) | 2.2 MB |
| One 1 bp-resolution track (RNA-seq, CAGE, DNase, ATAC, PRO-cap, splice sites), one haplotype | 0.3–0.9 MB (compressed float32; sparse signals are smaller) |
| One 128 bp ChIP-seq track (histone, TF), one haplotype | ~6 KB |
| Contact maps / splice junctions, one haplotype | ~0.1 MB each |
| Every output for one cell type (~575 tracks, 539 of them ChIP-TF), one haplotype | ~13 MB |

**Rule of thumb** for 524 kb windows (scale linearly with the window length):

```text
dataset size ≈ samples × windows × (4.3 MB + 2 haplotypes × Σ tracks × ~0.35 MB)
```

| Example | Size |
|---|---|
| Quickstart import: 1 sample per population × 6 genes, RNA-seq for 1 tissue | ~0.7 GB |
| 1000 Genomes: 3,202 samples × 44 genes, RNA-seq (6 tracks) | ~1.2 TB (380 MB per sample) |
| + training-axis alignment cache (`<dataset>/alignment_cache`, built on first use) | ~430 GB |
| + genotype predictor training caches (`results/cache/genotype_based_predictor`) | ~250 GB |
| Visualizer cache (cohort aggregates, indel indexes, observed data; `results/cache/visualizer`) | ~1 GB |

Where things are written:

| Location | Default | Change with |
|---|---|---|
| Registered datasets (`1kg_high_coverage`, …) | `/dados/GENOMICS_DATA` | `GENOMICS_DATA_ROOT` |
| Imported datasets | the output directory chosen in the Import form | – |
| Runs, caches | `results/` in the repository | `GENOMICS_RESULTS_ROOT`, `--cache-dir` |
| Background jobs (logs, state) | `results/visualizer/jobs` | `--jobs-dir` |
| AlphaGenome weights | `~/.cache/huggingface` | `HF_HOME` |

`genomics doctor` lists free space for each of these; the System page also shows it for every open dataset.

## Compute

| Task | Typical cost |
|---|---|
| Opening a 524 kb window for 12 samples × 2 haplotypes | 0.1–0.3 s first time, ~25 ms afterwards |
| Group means over 6,404 haplotypes (whole 1000 Genomes cohort) | ~50 s the first time, cached afterwards |
| Importing 1 sample per population × 6 windows from remote VCFs | about a minute |
| AlphaGenome, 1 haplotype-window (524 kb) | ~4 s hosted, ~0.8 s local GPU (after warm-up) |
| AlphaGenome for a whole cohort | samples × 2 haplotypes × windows calls: 3,202 × 2 × 44 ≈ 282 k calls ≈ 13 days hosted, ~3 days on one local GPU (DGX Spark) |
| Training a CNN on one gene set | minutes to hours on a GPU, depending on samples, windows and tracks |

Prediction jobs run one at a time per backend (`alphagenome` queue) and training/evaluation jobs one at
a time on the GPU (`gpu` queue); both are resumable, so a large cohort can be predicted in several sessions.
