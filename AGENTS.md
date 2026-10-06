# AGENTS.md

Guidance for AI agents and automation working in this repository.

## Project Snapshot

This repository is a Python package named `genomics` plus related scripts, configs, native code, third-party code, documentation, and legacy reproducibility modules.

The supported user-facing entrypoint is:

```bash
genomics ...
```

`python3 -m genomics` is also supported. Do not reintroduce old root-level package wrappers or legacy module entrypoints.

## First Checks

Start by inspecting the relevant files before editing. Useful anchors:

- `README.md`: current quick start and high-level component list.
- `docs/reference/repository-layout.md`: canonical repository layout.
- `docs/reference/cli.md`: supported CLI commands and examples.
- `configs/README.md`: canonical config locations.
- `pyproject.toml`: package metadata, optional dependency groups, console script.
- `src/genomics/cli.py`: authoritative CLI wiring.
- `src/genomics/workspace.py`: repository and data path defaults.
- `src/genomics/core/data_registry.py`: registered datasets and dataset path resolution.

Check worktree state before larger edits:

```bash
git status --short
```

Never revert or overwrite unrelated user changes.

## Canonical Layout Rules

Active Python code belongs under:

```text
src/genomics/
```

Use these namespaces for new active code:

- `genomics.core`: shared infrastructure.
- `genomics.workflows`: operational workflows and dataset builders.
- `genomics.predictors`: model/predictor pipelines.
- `genomics.converters`: converters.

Other canonical top-level locations:

- `configs/`: canonical configs, organized by domain.
- `docs/`: MkDocs documentation.
- `scripts/env/`: environment activation and installation scripts.
- `scripts/ops/`: operational run/monitor/setup scripts.
- `scripts/maintenance/`: maintenance and repair scripts.
- `scripts/diagnostics/`: diagnostic scripts.
- `scripts/experiments/`: experiment orchestration scripts.
- `scripts/dev/`: development demos/utilities.
- `legacy/`: historical reproducibility code only.
- `native/`: native non-Python tools.
- `third_party/`: modified third-party code.

Do not recreate removed root-level packages or wrappers such as:

- `genomics_cli.py`
- `genomics_pipeline/`
- `genotype_based_predictor/`
- `variant_transformer_predictor/`
- `snp_ancestry_predictor/`
- `vcf_to_23andme/`
- `genomes_analyzer_pipeline/`
- `build_non_longevous_dataset/`
- `neural_module/`
- `neural_ancestry_predictor_deprecated/`
- `neural_longevity_dataset/`

## CLI Rules

Expose new user-facing behavior through `src/genomics/cli.py` under `genomics ...`.

Preferred examples:

```bash
genomics --help
genomics audit-configs --fail-on-active-legacy
genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing
genomics genomes-analyzer run --config configs/genomes_analyzer/config_human_30x_low_memory.yaml
genomics genotype train configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics variant train configs/predictors/variant_transformer/repo_layout.example.yaml
genomics snp-ancestry run --config configs/predictors/snp_ancestry/default.yaml
```

When adding CLI commands:

- Keep `genomics --help` and subcommand `--help` working without heavy optional dependencies.
- Use lazy imports inside command handlers for optional/heavy dependencies such as `torch`, `sklearn`, `scipy`, `alphagenome`, and plotting libraries.
- Prefer forwarding to package modules via `_run_module` when a module already has a CLI-compatible main.
- Update `docs/reference/cli.md` when adding or changing supported commands.
- Update Bash completion in `src/genomics/completions/genomics.bash` when adding top-level or common subcommands.

## Config Schema Tools

Typed config discovery is exposed through:

```bash
genomics config describe genotype
genomics config describe variant
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics config validate configs/predictors/variant_transformer/repo_layout.example.yaml
genomics config schema genotype
```

Use these commands when changing config models or large YAML files. Add new config families to `src/genomics/core/config_schema.py` only after they have a typed schema and loader.

## Config Rules

Canonical config paths are documented in `configs/README.md` and currently include:

- `configs/genomes_analyzer/`
- `configs/predictors/genotype_based/`
- `configs/predictors/genotype_based/neural_legacy/`
- `configs/predictors/variant_transformer/`
- `configs/predictors/snp_ancestry/`
- `configs/workflows/non_longevous_dataset/`
- `configs/workflows/longevity_dataset/`
- `configs/workflows/alphagenome/`
- `configs/converters/vcf_to_23andme/`

For new configs and docs:

- Use canonical paths under `configs/`; do not add duplicate YAMLs directly under `configs/` root.
- Prefer dataset IDs over absolute dataset directories when supported.
- Avoid hardcoded machine-specific paths in active configs.
- Keep historical paths only in clearly historical docs or third-party history files.
- Keep legacy-only neural ancestry configs under `legacy/neural_ancestry_predictor_deprecated/configs/`, not under `configs/`.
- Run `genomics audit-configs --fail-on-active-legacy` after config changes.

## Dataset Rules

The default registered dataset ID is:

```text
1kg_high_coverage
```

Dataset resolution is centralized in `src/genomics/core/data_registry.py` and workspace defaults in `src/genomics/workspace.py`.

The canonical dataset layout for genotype/dataset-builder compatibility includes:

```text
dataset_metadata.json
layout_metadata.json
references/windows/<target>/ref.window.fa
references/windows/<target>/window_metadata.json
individuals/<sample>/windows/<target>/...
```

When changing dataset builders, readers, or predictors:

- Preserve the canonical layout above.
- Keep `dataset_metadata.json` and `layout_metadata.json` generation/read behavior compatible.
- Use `dataset_id` and `resolve_dataset(...)` where possible instead of embedding `/dados/...` paths.
- Validate with `genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing` when the local dataset is available.
- Use `--check-bcftools-chain` only when the task involves bcftools chain artifacts.

Known operational caveat: some variant transformer derived datasets may need materialization before `audit-data --fail-on-missing` can pass for their dataset IDs. See `docs/components/variant-transformer.md`.

## Dependencies

The package supports Python `>=3.8`.

Optional dependency groups are declared in `pyproject.toml`:

- `alphagenome`
- `genotype`
- `variant`
- `snp-ancestry`
- `legacy`
- `test`
- `docs`
- `all`

Install examples:

```bash
python3 -m pip install -e .
python3 -m pip install -e ".[test]"
python3 -m pip install -e ".[docs]"
python3 -m pip install -e ".[all]"
```

Do not make base package import or `genomics --help` require optional ML/scientific dependencies unless there is a deliberate project decision.

For Python 3.8 compatibility, avoid APIs introduced only in newer Python versions. For package resources, prefer APIs already used in the repository, such as `importlib.resources.read_text(...)`.

## Tests And Validation

Run the full test suite for broad changes:

```bash
python3 -m pytest tests
```

Useful targeted checks:

```bash
python3 -m pytest tests/test_repository_layout.py tests/test_config_compatibility.py
python3 -m pytest tests/test_genomics_cli.py tests/test_genomics_namespace.py
python3 -m pytest tests/test_canonical_dataset_layout.py
python3 -m compileall -q src legacy tests
python3 -m genomics audit-configs --fail-on-active-legacy
python3 -m genomics audit-data --dataset-id 1kg_high_coverage --fail-on-missing
```

Documentation validation when MkDocs is installed:

```bash
mkdocs build --strict
```

Package build validation when `build` is installed:

```bash
python3 -m build --sdist --wheel
```

If a command cannot run because an optional tool or dependency is missing, report that explicitly instead of silently skipping it.

## Documentation

The public docs are built with Material for MkDocs from `docs/` using `mkdocs.yml`.

When changing user-facing behavior:

- Update the relevant page under `docs/components/`, `docs/reference/`, or `docs/getting-started/`.
- Keep `README.md` short and introductory; put detailed instructions in `docs/`.
- If adding a new documentation page, add it to `mkdocs.yml` navigation or explain why it should remain unlisted.
- Do not commit the generated `site/` directory.

## Scripts

Use categorized script paths under `scripts/`; do not add root-level script wrappers.

Environment activation scripts:

```bash
source scripts/env/start_genomics_universal.sh
source scripts/env/start_genomics.sh
```

These scripts load Bash completion automatically when `genomics` is installed.

When moving or renaming scripts, update docs and tests that assert the categorized layout.

## Legacy, Native, And Third-Party Code

Legacy code under `legacy/` is for reproducibility. Keep new active behavior in `src/genomics/` instead.

Modified third-party code belongs under `third_party/`. Be conservative with edits there and preserve historical context unless the task explicitly requires cleanup.

Native tools belong under `native/`. Avoid mixing native build products into the Python package unless explicitly needed.

## Generated Files And Large Outputs

Do not commit generated caches, results, or build artifacts. Common generated locations include:

- `__pycache__/`
- `*.pyc`
- `*.egg-info/`
- `build/`
- `dist/`
- `results/`
- `site/`
- `wandb/`
- `runs/`

Large genomic datasets should remain outside Git and be referenced through dataset IDs, environment variables, or documented external paths.

## Coding Style

Prefer minimal, focused changes.

For Python changes:

- Use package imports rooted at `genomics.*` or relative imports inside a package.
- Keep optional imports lazy when they are not needed for import-time metadata or help text.
- Prefer `pathlib.Path` for filesystem paths.
- Keep YAML/JSON file IO explicit about encoding.
- Add tests for layout, CLI behavior, config compatibility, or dataset compatibility when changing those areas.

For shell scripts:

- Preserve current categorized locations.
- Quote paths and variables where practical.
- Avoid hardcoded local paths in new scripts unless the script is explicitly machine-specific and documented as such.

## Git Practices

- Do not commit unless explicitly asked.
- Before committing, inspect `git status --short` and review the diff.
- Stage only intended files.
- Never use destructive commands such as `git reset --hard` or `git checkout --` unless explicitly requested and confirmed.
- If unrelated user changes exist, leave them untouched and mention them if they affect validation.
