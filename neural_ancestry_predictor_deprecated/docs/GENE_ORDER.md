# Gene Order in Neural Ancestry Predictor

## Overview

The **gene order** is a critical piece of metadata that defines the sequence in which genes appear in the processed cache data. This order directly impacts how the neural network models interpret the input data matrix.

## Why Gene Order Matters

### Data Structure

The input to the neural networks is a 2D matrix with shape `[num_rows, window_center_size]` where:
- `num_rows = num_genes × num_tracks_per_gene`
- For RNA-seq data: `num_tracks_per_gene = 6` (6 ontology tracks)
- For 11 genes: `num_rows = 66` (11 genes × 6 tracks)

Each gene occupies 6 consecutive rows in the matrix, corresponding to its 6 RNA-seq ontology tracks.

### Example

For a dataset with genes in order `["MC1R", "TYRP1", "TYR", ...]`:
- Rows 0-5: MC1R tracks (6 ontologies)
- Rows 6-11: TYRP1 tracks (6 ontologies)
- Rows 12-17: TYR tracks (6 ontologies)
- ... and so on

## How Gene Order is Determined

### Source of Truth

The gene order comes from the **base dataset** (`GenomicLongevityDataset`) which reads it from:
```
dataset_dir/individuals/{sample_id}/individual_metadata.json
```

The `"windows"` field in this file contains the list of genes in the order they were processed.

### Important Notes

1. **Not Alphabetical**: The gene order is **NOT** alphabetical. It depends on how the base dataset was built.

2. **Consistent Across Samples**: All samples in the same base dataset have the same gene order.

3. **Stored in Cache**: When a processed cache is created, the gene order is automatically saved in the cache metadata.

## Gene Order in Cache Metadata

### New Caches

When you create a new processed cache, the system automatically:
1. Reads gene order from the first sample of the base dataset
2. Saves it in `cache_dir/metadata.json` under the key `"gene_order"`
3. Also saves `"tracks_per_gene": 6`

Example metadata:
```json
{
  "gene_order": ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", 
                 "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"],
  "tracks_per_gene": 6,
  ...
}
```

### Old Caches (Backward Compatibility)

For caches created before this feature was added:
1. On first load, the system detects missing `gene_order`
2. Temporarily loads the base dataset to extract the gene order
3. Updates the cache metadata.json (does **NOT** regenerate cache data)
4. Subsequent loads read from the updated metadata

This ensures **no forced cache regeneration** while maintaining correctness.

## Using Gene Order

### In YAML Configuration

The `genes_to_use` parameter in the YAML config allows you to select specific genes:

```yaml
dataset_input:
  genes_to_use:
    - "MC1R"
    - "TYR"
    - "OCA2"
    - "SLC24A5"
    - "SLC45A2"
```

**Important**: 
- Gene names must match exactly those in the gene_order
- Genes are used in the dataset's order, not the order you list them
- Invalid gene names will cause an error with available options displayed

### In Model Code

The neural network models automatically receive the gene order through the config:

```python
# Models read gene_order from config (populated from cache metadata)
GENE_ORDER = config['dataset_input'].get('gene_order', 
    DEFAULT_FALLBACK_ORDER)
```

This makes the code:
- **Portable**: Works with different datasets automatically
- **Robust**: No hard-coded assumptions about gene order
- **Verifiable**: Gene order is explicitly documented in cache metadata

## Validation

When using `genes_to_use`, the system validates that:
1. Each gene exists in the dataset's gene_order
2. Reports available genes if validation fails
3. Selects the correct rows based on gene position

Example error message:
```
ValueError: Gene inválido: INVALID_GENE. 
Opções: ['MC1R', 'TYRP1', 'TYR', 'SLC45A2', 'DDB1', 
         'EDAR', 'MFSD12', 'OCA2', 'HERC2', 'SLC24A5', 'TCHH']
```

## Best Practices

1. **Always check cache metadata**: Verify gene_order in `metadata.json` before analysis

2. **Document experiments**: When reporting results, include the gene_order used

3. **Consistent base datasets**: When comparing experiments, ensure they use the same base dataset (and thus same gene order)

4. **Update old caches**: The first run with an old cache will update metadata automatically - this is expected behavior

## Troubleshooting

### "Gene order not found" Warning

If you see:
```
⚠ Cache não contém ordem dos genes. Reconstruindo...
```

This is **normal** for old caches. The system will:
- Load the base dataset temporarily  
- Extract gene order
- Update metadata
- Continue normally

No action needed - this happens once per old cache.

### Gene Order Mismatch

If gene order changes between cache creation and usage, you'll need to regenerate the cache. Signs of mismatch:
- Unexpected model performance
- Genes appear in wrong order during inspection

**Solution**: Delete cache and regenerate with current base dataset.

## Technical Details

### Where Gene Order is Used

1. **save_processed_dataset()**: Reads from base dataset, saves to metadata
2. **load_processed_dataset()**: Reads from metadata, passes to config
3. **prepare_data()**: Passes gene_order to config
4. **Model __init__()**: Uses config['dataset_input']['gene_order']
5. **Model forward()**: Selects correct rows based on gene positions

### Fallback Behavior

If gene_order cannot be determined:
- Models use a default fallback order
- Warning is displayed
- May cause incorrect behavior if actual order differs

**Always ensure gene_order is properly saved in cache metadata.**

## See Also

- [Dynamic Gene Loading](DYNAMIC_GENE_LOADING.md) - How genes are loaded dynamically
- [Dataset Verification](README_VERIFY.md) - Tools to verify gene order correctness
- Main [README](../README.md) - General usage and configuration

