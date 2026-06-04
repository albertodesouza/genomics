# Experiments

Experiment documentation is split by purpose:

| Document | Purpose |
|---|---|
| [Research Pipelines](../RESEARCH_PIPELINES.md) | Active ML pipeline conventions and smoke status |
| [Important Experiments](../IMPORTANT_EXPERIMENTS.md) | Consolidated benchmark and reproduction commands |
| [Architecture](../ARCHITECTURE.md) | Package layout and dependency rules |

## Recommended Workflow

1. Pick a canonical config under `configs/`.
2. Run `genomics audit-configs --fail-on-active-legacy`.
3. Run `genomics audit-data --dataset-id <id> --fail-on-missing`.
4. Start with a short smoke run when possible.
5. Record important results in `docs/IMPORTANT_EXPERIMENTS.md`.

## Reproducibility Artifacts

Each integrated run should preserve:

- the original config;
- a resolved config;
- a manifest;
- metrics and model artifacts;
- enough path and dataset metadata to rerun later.
