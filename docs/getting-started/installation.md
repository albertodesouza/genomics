# Installation

## Python Package

Install from the repository root:

```bash
python3 -m pip install -e .
```

This installs the `genomics` console script:

```bash
genomics --help
```

The package metadata lives in `pyproject.toml`. The root `setup.py` is only a compatibility shim for older editable installs.

## Conda Environment

The project includes scripts for creating a Conda environment with genomics command-line tools:

```bash
scripts/env/install_genomics_env.sh
```

The environment script installs tools such as `bcftools`, `samtools`, `bwa-mem2`, `gatk4`, `plink`, `plink2`, `admixture`, `fastqc`, `multiqc`, `cutadapt`, and Python support libraries.

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

Material for MkDocs is optional for development, but recommended for browsing docs locally:

```bash
python3 -m pip install -e ".[docs]"
mkdocs serve
```
