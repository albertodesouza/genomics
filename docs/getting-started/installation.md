# Installation

**What do I need?** See [Requirements](requirements.md) for what each feature needs (packages, tools,
GPU, disk). After installing, `genomics doctor` tells you what works on this machine and how to enable
the rest.

## Recommended: one command

From the repository root, on Linux or macOS with [Miniforge](https://github.com/conda-forge/miniforge)
(or any conda):

```bash
scripts/env/install.sh
conda activate genomics
genomics visualize --open
```

`install.sh` creates (or reuses) the conda environment `genomics` with Python 3.11, `bcftools` and
`samtools`, installs this repository with the `[visualizer]` extra, and ends with `genomics doctor`.
It is safe to re-run. Options:

| Option | Adds |
|---|---|
| `--training` | PyTorch, scikit-learn, SciPy, pydantic (`[genotype]`): training, evaluation and Perturbation Lab scoring. On NVIDIA ARM machines (DGX Spark) PyTorch comes from its CUDA wheel index; set `TORCH_INDEX` to choose another |
| `--alphagenome-server` | AlphaGenome on this machine's NVIDIA GPU, in a separate environment (`genomics alphagenome server setup`, below) |
| `--env NAME` / `--prefix DIR` | another environment name or location |
| `--no-conda` | a plain virtual environment in `.venv`; install `bcftools`/`samtools` yourself (`apt install bcftools samtools`) |

No conda yet? `scripts/env/install_conda_universal.sh` installs Miniforge.

## With pip

```bash
python3 -m pip install -e .                        # CLI + visualizer (numpy, PyYAML)
python3 -m pip install -e ".[visualizer]"          # + gene models, VCF import and AlphaGenome client (recommended)
python3 -m pip install -e ".[visualizer,genotype]" # + training and Perturbation Lab scoring
```

| Extra | Packages | For |
|---|---|---|
| *(base)* | `numpy`, `PyYAML` | `genomics` CLI, browsing datasets in the visualizer |
| `visualizer` | `pandas`, `pyarrow`, `alphagenome` (Python ≥ 3.10) | gene models, gene search, VCF import, AlphaGenome predictions |
| `genotype` | `torch`, `scikit-learn`, `scipy`, `pydantic`, `captum`, `matplotlib`, … | genotype-based predictor: training, evaluation, Perturbation Lab scoring |
| `variant` | `torch`, `scikit-learn`, … | variant transformer predictor |
| `snp-ancestry` | `scikit-learn`, `scipy`, … | SNP ancestry predictor |
| `alphagenome` | AlphaGenome client + plotting | the standalone `genomics alphagenome analyze/integrate` workflow |
| `test`, `docs` | `pytest`, `mkdocs-material` | development |
| `all` | everything above | – |

Dataset import also needs `bcftools` and `samtools` on `PATH` (or in the same environment's `bin`).
The package supports Python ≥ 3.8; the AlphaGenome client needs ≥ 3.10, and 3.10–3.11 are what the
project is tested with.

## AlphaGenome

Predictions use one backend, chosen on the visualizer's **AlphaGenome** page:

- **Hosted API**: put `ALPHAGENOME_API_KEY=...` in `~/.env` (or export it). Keys are free for
  non-commercial use at <https://deepmind.google.com/science/alphagenome>.
- **This machine** (NVIDIA GPU): install the model server once, then click *Start server* in the
  visualizer (or run it yourself):

  ```bash
  genomics alphagenome server setup                     # clone alphagenome_research, create its env with CUDA jax, check weights
  genomics alphagenome server setup --download-weights  # after accepting the terms on Hugging Face and `hf auth login`
  genomics alphagenome server start                     # or the visualizer's Start server button
  genomics alphagenome server check --predict
  ```

  `setup` clones <https://github.com/FeLiPeOLi7/alphagenome_research> next to this repository (or under
  `~/.local/share/genomics`), creates the conda env `alphagenome` (a venv in the checkout without conda)
  and installs it with `jax[cuda12]`. Use `--jax cuda13` for CUDA 13 drivers, `--python PATH` for an
  existing environment, `--dry-run` to see the commands. Hardware and timing: [Requirements](requirements.md#alphagenome-on-your-own-gpu).
- **Remote server**: the address of a machine running `genomics alphagenome server start --host 0.0.0.0`.

## Checking the installation

```bash
genomics doctor          # per-feature report with fixes; --json for scripts
python3 -m pytest tests  # with the [test] extra
```

## Full bioinformatics toolchain

The FASTQ → VCF pipeline (`genomes-analyzer`) and other workflows use many more command-line tools. The
original environment scripts install them:

```bash
scripts/env/install_genomics_env.sh                     # x86_64: bcftools, samtools, bwa(-mem2), GATK, plink, admixture, ...
scripts/env/install_genomics_env_linux_aarch64_nvidia.sh # DGX Spark / ARM + NVIDIA (no GATK/plink/admixture builds)
```

The package metadata lives in `pyproject.toml`. The root `setup.py` is only a compatibility shim for
older editable installs.

## VEP

Install Ensembl VEP with:

```bash
source scripts/maintenance/vep_install.sh
```

Optional environment variables:

| Variable | Purpose |
|---|---|
| `VEP_BRANCH` | Ensembl VEP branch or tag |
| `VEP_SPECIES` | Species name, usually `homo_sapiens` |
| `VEP_ASSEMBLY` | Genome assembly, usually `GRCh38` |
| `VEP_CACHE_DIR` | VEP cache directory |
| `VEP_DIR` | VEP installation directory |

## Documentation Tools

```bash
python3 -m pip install -e ".[docs]"
mkdocs serve
```
