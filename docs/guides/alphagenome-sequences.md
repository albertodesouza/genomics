# AlphaGenome Sequences And Outputs

This guide covers the practical AlphaGenome inputs used by `genomics alphagenome ...` and by dataset builders that call AlphaGenome internally.

## Commands

```bash
genomics alphagenome analyze -- -i sequence.fasta -k API_KEY -o results/
genomics alphagenome integrate -- --integrated --vcf vcf/sample.vcf.gz --ref refs/GRCh38.fa --api-key API_KEY --output integrated_analysis/
genomics alphagenome tracks --api-key API_KEY --output configs/workflows/alphagenome/tracks.json
```

Arguments after `--` are forwarded to the underlying AlphaGenome modules.

## Supported Sequence Sizes

AlphaGenome accepts a fixed set of sequence lengths:

| Length | Typical Use |
|---:|---|
| 2,048 bp | Small windows and quick tests. |
| 16,384 bp | Individual genes or local regulatory regions. |
| 131,072 bp | Larger genes, gene clusters, and broad regulatory regions. |
| 524,288 bp | Large loci and regional context. |
| 1,048,576 bp | Megabase-scale context and contact-map-style analyses. |

Check FASTA length before submitting:

```bash
grep -v "^>" sequence.fasta | tr -d '\n' | wc -c
```

If a sequence is shorter than the target length, pad with `N`. If it is longer, crop around the variant, gene, or region of interest. For variant analysis, center the variant whenever possible.

## Getting Sequences

You can extract real genomic sequence from Ensembl, UCSC, NCBI, or a local indexed reference. Keep the assembly consistent with downstream VCFs and annotations.

Example Ensembl REST request for a 2,048 bp GRCh38 interval:

```bash
curl 'https://rest.ensembl.org/sequence/region/human/11:5227002..5229049?coord_system_version=GRCh38' \
  -H 'Content-type:text/x-fasta' > HBB_2048bp.fasta
```

Example local reference extraction:

```bash
samtools faidx GRCh38.fa chr11:5227002-5229049 > HBB_2048bp.fasta
```

Before use, confirm that the FASTA contains only expected bases (`A`, `C`, `G`, `T`, and optionally `N`) and that the length is one of the supported sizes.

## Output Families

AlphaGenome output names are families, not individual assays. Common families include:

| Output | Meaning |
|---|---|
| `RNA_SEQ` | Gene expression signal. |
| `CAGE` | Transcription start site activity. |
| `PROCAP` | Precision run-on transcription initiation signal. |
| `ATAC` | Chromatin accessibility. |
| `DNASE` | DNase hypersensitivity/accessibility. |
| `CHIP_HISTONE` | Histone modification tracks such as H3K27ac and H3K4me3. |
| `CHIP_TF` | Transcription factor tracks such as CTCF. |
| `CONTACT_MAPS` | Chromatin contact predictions. |
| `SPLICE_JUNCTIONS` | Splice junction predictions. |
| `SPLICE_SITES` | Donor/acceptor splice site predictions. |
| `SPLICE_SITE_USAGE` | Splice site usage predictions. |

Use family names in commands:

```bash
genomics alphagenome analyze -- \
  -i sequence.fasta \
  -k API_KEY \
  -o results/ \
  --outputs RNA_SEQ ATAC CHIP_HISTONE
```

Do not request individual names such as `H3K27AC` directly. Use `CHIP_HISTONE` or `CHIP_TF` and select ontology/track metadata downstream when supported.

## Track And Ontology Discovery

Export available output and ontology metadata with:

```bash
genomics alphagenome tracks --api-key API_KEY --output configs/workflows/alphagenome/tracks.json
```

For diagnostics outside the CLI, use:

```bash
python scripts/diagnostics/check_alphagenome_outputs.py
python scripts/diagnostics/check_dna_client_methods.py
```

The exported metadata helps choose genotype predictor config fields such as:

```yaml
dataset_input:
  alphagenome_outputs: ["RNA_SEQ"]
  ontology_terms: ["UBERON:0000955"]
```

## Recommended Output Sets

| Use Case | Outputs |
|---|---|
| Quick smoke test | `RNA_SEQ` |
| Expression and promoters | `RNA_SEQ CAGE PROCAP` |
| Regulatory variant analysis | `RNA_SEQ CAGE ATAC CHIP_HISTONE` |
| Chromatin accessibility | `ATAC DNASE` |
| Epigenetic state | `CHIP_HISTONE CHIP_TF` |
| Splicing | `RNA_SEQ SPLICE_JUNCTIONS SPLICE_SITES SPLICE_SITE_USAGE` |
| 3D context | `CONTACT_MAPS CHIP_TF` |

More outputs increase runtime and storage. For large sequence lengths, start with a narrow output set and expand only after the workflow is validated.

## Relationship To Legacy Guides

Older notes remain in the repository for provenance:

- `docs/historical/alphagenome/DOWNLOAD_SEQUENCES.md`
- `docs/historical/alphagenome/AVAILABLE_OUTPUTS.md`
- `docs/historical/alphagenome/SUPPORTED_SIZES.md`
- `docs/historical/APPLIED_FIXES.md`

Use this page as the active guide for new runs.
