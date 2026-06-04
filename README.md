# Genomics

Multi-pipeline genomics workspace with one primary command-line interface: `genomics`.

The repository contains an operational genome processing workflow, dataset builders, AlphaGenome integration, ancestry/model predictors, shared ML infrastructure, a native C++ tool, a modified third-party ancestry calculator, and legacy reproducibility code.

## Quick Start

Install the package in editable mode from the repository root:

```bash
python3 -m pip install -e .
```

Use the CLI:

```bash
genomics --help
genomics genomes-analyzer run --config configs/genomes_analyzer/config_human_30x_low_memory.yaml
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml
```

Activate the Conda environment and Bash completion:

```bash
source scripts/env/start_genomics_universal.sh
```

The activation script loads Bash completion automatically when `genomics` is installed.

## Documentation

Public documentation is published with GitHub Pages at:

https://albertodesouza.github.io/genomics/

The full documentation is organized under `docs/` and can be served with Material for MkDocs:

```bash
python3 -m pip install -e ".[docs]"
mkdocs serve
```

Start with:

- [Documentation Home](docs/index.md)
- [Installation](docs/getting-started/installation.md)
- [CLI Reference](docs/reference/cli.md)
- [Repository Layout](docs/reference/repository-layout.md)
- [Configuration Layout](configs/README.md)

## Main Components

| Component | Documentation |
|---|---|
| Genomes Analyzer workflow | [docs/components/genomes-analyzer.md](docs/components/genomes-analyzer.md) |
| Genotype-based predictor | [docs/components/genotype-predictor.md](docs/components/genotype-predictor.md) |
| Variant transformer predictor | [docs/components/variant-transformer.md](docs/components/variant-transformer.md) |
| SNP ancestry predictor | [docs/components/snp-ancestry.md](docs/components/snp-ancestry.md) |
| VCF to 23andMe converter | [docs/components/vcf-to-23andme.md](docs/components/vcf-to-23andme.md) |
| AlphaGenome workflow | [docs/components/alphagenome.md](docs/components/alphagenome.md) |
| Dataset builders | [docs/components/dataset-builders.md](docs/components/dataset-builders.md) |
| Native and third-party tools | [docs/components/native-and-third-party.md](docs/components/native-and-third-party.md) |
| Legacy code | [docs/components/legacy.md](docs/components/legacy.md) |

## Development Checks

Run the test suite:

```bash
python3 -m pytest tests
```

Run documentation build when MkDocs is installed:

```bash
mkdocs build --strict
```
