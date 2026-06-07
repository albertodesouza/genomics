# AlphaGenome Workflow

The AlphaGenome workflow wraps sequence analysis, variant-effect analysis, integration helpers, and output track discovery.

## CLI

```bash
genomics alphagenome analyze -- -i sequence.fasta -k API_KEY -o results/
genomics alphagenome integrate -- --integrated --vcf vcf/sample.vcf.gz --ref refs/GRCh38.fa --api-key API_KEY --output integrated_analysis/
genomics alphagenome tracks --api-key API_KEY --output configs/workflows/alphagenome/tracks.json
```

Arguments after `--` are forwarded to the underlying workflow modules.

## Code And Configs

| Kind | Path |
|---|---|
| Code | `src/genomics/workflows/alphagenome/` |
| Configs | `configs/workflows/alphagenome/` |
| Sequence and output guide | `docs/guides/alphagenome-sequences.md` |
| Track discovery script | `scripts/diagnostics/check_alphagenome_outputs.py` |

## Typical Uses

- predict functional genomic features from FASTA sequences;
- compare variant effects;
- export output ontology metadata;
- integrate VCF/BED/gene-list regions with AlphaGenome predictions.

## Implementation Structure

| Module | Responsibility |
|---|---|
| `neural_module.py` | Main AlphaGenome sequence/variant analysis workflow |
| `neural_integration.py` | Extracts regions from VCF/BED/gene lists and invokes analysis |
| `alphagenome_output_tracks.py` | Queries/export available AlphaGenome tracks and ontology metadata |
| `neural_visualizations_advanced.py` | Plotting and comparison helpers |

Imports of the external `alphagenome` package are intentionally lazy. This keeps `genomics --help` and basic package imports working in environments where AlphaGenome is not installed.

## Analyze Flow

```text
FASTA input
  -> sequence validation and size checks
  -> AlphaGenome model request
  -> selected output tracks
  -> arrays, metadata, plots, and optional variant comparison
```

The workflow supports the sequence sizes accepted by AlphaGenome. See [AlphaGenome Sequences And Outputs](../guides/alphagenome-sequences.md) for practical extraction, sizing, and output selection examples.

## Integration Flow

The integration command can derive sequences from variant or region inputs and run AlphaGenome analysis in one pass.

Supported sources include:

- VCF regions;
- BED intervals;
- gene lists;
- manually supplied FASTA files.

Use `genomics alphagenome tracks` to export the currently available AlphaGenome output and ontology metadata:

```bash
genomics alphagenome tracks --api-key API_KEY --output configs/workflows/alphagenome/tracks.json
```

That metadata is useful when selecting `alphagenome_outputs` and `ontology_terms` for dataset builders and genotype predictor configs.

## Outputs

Outputs depend on selected tracks and options, but commonly include:

- NumPy arrays for predicted tracks;
- CSV/JSON metadata;
- plots comparing reference and alternate sequence predictions;
- summaries of available output ontology labels.
