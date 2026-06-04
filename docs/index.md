# Genomics Workspace

`genomics` is a multi-pipeline workspace for genome processing and genomics ML experiments. The project is organized around one installed command, `genomics`, and a canonical Python package under `src/genomics/`.

## What Is Included

| Area | Purpose |
|---|---|
| Operational workflow | FASTQ/BAM/CRAM/VCF processing through `genomics genomes-analyzer ...` |
| ML predictors | Dense genotype models and sparse variant transformer models |
| SNP ancestry | Allele-frequency ancestry prediction and related conversion steps |
| AlphaGenome | Sequence and variant-effect analysis helpers |
| Dataset builders | 1000 Genomes derived dataset construction workflows |
| Converters | VCF to 23andMe raw format conversion |
| Shared core | Config IO, dataset registry, experiment manifests, metrics, splitting, checkpoints |
| Native/third-party | C++ gene comparison tool and modified FROGAncestryCalc |
| Legacy | Historical reproducibility modules isolated under `legacy/` |

## Recommended Reading Order

1. [Installation](getting-started/installation.md)
2. [Environment](getting-started/environment.md)
3. [First Commands](getting-started/first-commands.md)
4. [Components Overview](components/index.md)
5. [CLI Reference](reference/cli.md)
6. [Repository Layout](reference/repository-layout.md)

## Primary Command

```bash
genomics --help
```

All active workflows should be invoked through `genomics ...`, not through old module paths or removed root scripts.
