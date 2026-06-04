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
| Sequence guide | `docs/DOWNLOAD_SEQUENCES.md` |
| Output types | `docs/AVAILABLE_OUTPUTS.md` |
| Supported sizes | `docs/SUPPORTED_SIZES.md` |

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

The workflow supports the sequence sizes accepted by AlphaGenome. See [Supported Sizes](../SUPPORTED_SIZES.md) for practical resizing examples.

## Integration Flow

The integration command can derive sequences from variant or region inputs and run AlphaGenome analysis in one pass.

Supported sources include:

- VCF regions;
- BED intervals;
- gene lists;
- manually supplied FASTA files.

## Outputs

Outputs depend on selected tracks and options, but commonly include:

- NumPy arrays for predicted tracks;
- CSV/JSON metadata;
- plots comparing reference and alternate sequence predictions;
- summaries of available output ontology labels.
