# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Snapshot

Multi-pipeline genomics workspace exposed through a single console-script entrypoint, `genomics` (also runnable as `python3 -m genomics`). It bundles an operational genome-processing workflow, dataset builders, AlphaGenome integration, ancestry/model predictors, shared ML infrastructure, a native C++ tool, a modified third-party ancestry calculator, and legacy reproducibility code.

Do not reintroduce old root-level package wrappers or legacy module entrypoints (e.g. `genomics_cli.py`, `genomics_pipeline/`, `genotype_based_predictor/`, `variant_transformer_predictor/`, `snp_ancestry_predictor/`, `vcf_to_23andme/`, `genomes_analyzer_pipeline/`, `build_non_longevous_dataset/`, `neural_module/`, `neural_ancestry_predictor_deprecated/`, `neural_longevity_dataset/`). Some of these still exist as top-level directories from before the unification — treat them as historical, not as places for new active code.

## Setup

```bash
python3 -m pip install -e .                # base install
python3 -m pip install -e ".[test]"         # + pytest
python3 -m pip install -e ".[genotype]"     # + torch/sklearn/scipy for genotype predictor
python3 -m pip install -e ".[variant]"      # + torch/sklearn for variant transformer
python3 -m pip install -e ".[snp-ancestry]" # + sklearn/scipy for SNP ancestry
python3 -m pip install -e ".[alphagenome]"  # + alphagenome client
python3 -m pip install -e ".[docs]"         # + mkdocs-material
python3 -m pip install -e ".[all]"          # everything
source scripts/env/start_genomics_universal.sh   # activates conda env + Bash completion for `genomics`
```

Package supports Python `>=3.8` — avoid APIs from newer versions; use `importlib.resources.read_text(...)` style APIs already present in the repo for package resources.

## Commands

```bash
genomics --help
genomics genomes-analyzer run --config configs/genomes_analyzer/config_human_30x_low_memory.yaml
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml

# config schema tools
genomics config describe genotype
genomics config describe variant
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics config schema genotype

# repo/dataset audits
genomics audit-configs --fail-on-active-legacy
genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing
```

### Tests

```bash
python3 -m pytest tests                                                        # full suite
python3 -m pytest tests/test_repository_layout.py tests/test_config_compatibility.py
python3 -m pytest tests/test_genomics_cli.py tests/test_genomics_namespace.py
python3 -m pytest tests/test_canonical_dataset_layout.py
python3 -m pytest tests/test_some_file.py::test_name                           # single test
python3 -m compileall -q src legacy tests
```

Other validation:

```bash
mkdocs build --strict                       # docs, if mkdocs-material installed
python3 -m build --sdist --wheel            # package build, if `build` installed
```

If a command can't run because an optional dependency/tool is missing, say so explicitly rather than silently skipping it.

## Architecture

All active Python code lives under `src/genomics/`, organized by responsibility rather than by historical module:

- `genomics.core` — shared infrastructure used across predictors/workflows: config I/O and typed schemas (`config_schema.py`), dataset registry/resolution (`data_registry.py`), reference genome registry (`reference_registry.py`), dataset metadata, splitting, metrics, checkpointing, training utils, torch/sklearn helpers, wandb integration.
- `genomics.workflows` — operational pipelines and dataset builders: `genomes_analyzer/` (FASTQ → alignment → variants → reports pipeline), `alphagenome/` (AlphaGenome model integration), `dataset_builders/non_longevous/` (dataset construction).
- `genomics.predictors` — model/predictor pipelines, each a self-contained subpackage: `genotype_based/` (CNN models over aligned genomic windows, with `alignment/`, `data/`, `models/`, `experiments/`, `analysis/`, `apps/`, `tools/` subfolders), `variant_transformer/` (transformer over variant calls), `snp_ancestry/` (classical ML ancestry classifier).
- `genomics.converters` — format converters, e.g. `vcf_to_23andme/`.
- `src/genomics/cli.py` is the single authoritative entrypoint wiring all subcommands (`genomics <subcommand> ...`); it forwards to package modules via `_run_module` when a module already has a CLI-compatible `main`. Keep `genomics --help` and per-subcommand `--help` working without importing heavy optional dependencies (torch, sklearn, scipy, alphagenome, plotting libs) — use lazy imports inside command handlers.
- `src/genomics/workspace.py` defines repository and default data path locations used across the CLI and predictors.

Other top-level locations: `configs/` (canonical YAML configs, see below), `docs/` (MkDocs site, `docs/components/` per-component pages, `docs/reference/` for CLI/layout reference), `scripts/{env,ops,maintenance,diagnostics,experiments,dev}/`, `legacy/` (historical reproducibility code only — do not add new active behavior here), `native/` (non-Python tools, e.g. `genes_difference_count/`), `third_party/` (modified third-party code, e.g. `FROGAncestryCalc/` — be conservative, preserve historical context).

### Configs

Canonical config directories mirror the predictor/workflow structure, documented in `configs/README.md`:

- `configs/genomes_analyzer/`
- `configs/predictors/genotype_based/` (plus `icann/` and `neural_legacy/` subfolders)
- `configs/predictors/variant_transformer/`
- `configs/predictors/snp_ancestry/` (plus `icann/`)
- `configs/workflows/non_longevous_dataset/`, `configs/workflows/longevity_dataset/`, `configs/workflows/alphagenome/`
- `configs/converters/vcf_to_23andme/`

Each canonical directory should have a `default.yaml` as a portable, commented starting point. Don't add new YAMLs directly under `configs/` root; don't hardcode machine-specific paths in active configs. Legacy-only neural ancestry configs live under `legacy/neural_ancestry_predictor_deprecated/configs/`, not under `configs/`. Add new config families to `src/genomics/core/config_schema.py` only once they have a typed schema and loader. Run `genomics audit-configs --fail-on-active-legacy` after config changes.

### Dataset layout

Default registered dataset ID: `1kg_high_coverage`. Resolution is centralized in `src/genomics/core/data_registry.py`; workspace defaults in `src/genomics/workspace.py`. Prefer `dataset_id` + `resolve_dataset(...)` over embedding raw filesystem paths (e.g. `/dados/...`).

Canonical dataset directory layout (genotype/dataset-builder compatible):

```text
dataset_metadata.json
layout_metadata.json
references/windows/<target>/ref.window.fa
references/windows/<target>/window_metadata.json
individuals/<sample>/windows/<target>/...
```

Preserve this layout, and `dataset_metadata.json`/`layout_metadata.json` generation/read compatibility, when touching dataset builders/readers/predictors. Validate with `genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing` when the dataset is available locally. Some variant-transformer derived datasets need materialization before that audit can pass — see `docs/components/variant-transformer.md`.

## Coding conventions

- Imports rooted at `genomics.*` (or relative within a package); keep optional heavy imports (torch, sklearn, scipy, alphagenome, plotting) lazy unless needed at import time for metadata/help text.
- Use `pathlib.Path` for filesystem paths; be explicit about encoding on YAML/JSON I/O.
- Add tests when changing layout, CLI behavior, config compatibility, or dataset compatibility.
- Shell scripts: keep them in their categorized `scripts/` subdirectory, quote paths/variables, avoid hardcoded local paths unless the script is explicitly machine-specific and documented as such.
- When adding a CLI command, update `docs/reference/cli.md` and the Bash completion at `src/genomics/completions/genomics.bash`.
- Don't commit generated caches/artifacts: `__pycache__/`, `*.pyc`, `*.egg-info/`, `build/`, `dist/`, `results/`, `site/`, `wandb/`, `runs/`. Large genomic datasets stay outside Git — reference them via dataset IDs/env vars/documented external paths.

## Git practices

- Never use destructive commands (`git reset --hard`, `git checkout --`, etc.) unless explicitly requested and confirmed.
- If unrelated uncommitted changes exist in the worktree, leave them untouched and mention them if they affect validation.
