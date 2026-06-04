# Config Schemas

Large YAML configs are easier to maintain when their fields are discoverable from the command line.

The `genomics config` command exposes typed schemas for config families that already have Pydantic models.

Currently supported kinds:

| Kind | Component |
|---|---|
| `genotype` | Genotype-based predictor |
| `variant` | Variant transformer predictor |

## Describe Fields

Print every known field, including nested sections, types, defaults, required status, and available descriptions:

```bash
genomics config describe genotype
genomics config describe variant
```

Machine-readable output is available with:

```bash
genomics config describe genotype --json
```

## Validate YAML

Validate a config file against its typed schema:

```bash
genomics config validate configs/predictors/genotype_based/genes_1000_all_3ontologies.yaml
genomics config validate configs/predictors/variant_transformer/repo_layout.example.yaml
```

The command infers the kind from canonical config paths. If the file lives outside the canonical tree, pass `--kind` explicitly:

```bash
genomics config validate /tmp/experiment.yaml --kind genotype
```

Use JSON output for automation:

```bash
genomics config validate configs/predictors/variant_transformer/repo_layout.example.yaml --json
```

## JSON Schema

Export JSON Schema for editor integration or external tools:

```bash
genomics config schema genotype > genotype.schema.json
genomics config schema variant > variant.schema.json
```

## Current Scope

This first schema layer covers `genotype` and `variant` configs because those components already use typed Pydantic configuration models.

Configs for SNP ancestry, Genomes Analyzer, dataset builders, AlphaGenome, and converters still use more free-form YAML loading in parts of the codebase. Add typed models for those components before exposing them through `genomics config`.
