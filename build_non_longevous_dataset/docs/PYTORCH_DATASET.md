# PyTorch Dataset for Longevity Genomics

> **📦 New**: `build_non_longevous_dataset` now generates a complete PyTorch Dataset for training neural networks for longevity inference!

## 📑 Table of Contents

- [Overview](#overview)
- [Dataset Structure](#dataset-structure)
- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Complete API](#complete-api)
- [Data Format](#data-format)
- [Advanced Examples](#advanced-examples)
- [Compression and Distribution](#compression-and-distribution)
- [FAQ](#faq)

---

## Overview

The `build_non_longevous_dataset` module now generates a data structure compatible with PyTorch Dataset/DataLoader, enabling efficient neural network training.

### Features

- ✅ **PyTorch Compatible** Dataset and DataLoader
- ✅ **Lazy loading** (on-demand) to save memory
- ✅ **Complete metadata** per individual and global
- ✅ **Ancestry data** from FROGAncestryCalc included
- ✅ **Multiple genomic windows** per individual
- ✅ **AlphaGenome predictions** (.npz) for each haplotype
- ✅ **Optional FASTA sequences**
- ✅ **Customizable transformations**

---

## Dataset Structure

### Directory Structure

```
non_longevous_dataset/
├── dataset_metadata.json              # Metadados globais
├── individuals/
│   ├── HG01879/
│   │   ├── individual_metadata.json   # Metadados do indivíduo
│   │   └── windows/
│   │       ├── CYP2B6/               # Janela genômica (gene/SNP)
│   │       │   ├── ref.window.fa
│   │       │   ├── HG01879.H1.window.fixed.fa
│   │       │   ├── HG01879.H2.window.fixed.fa
│   │       │   ├── predictions_H1/
│   │       │   │   ├── rna_seq.npz   # Predições RNA-seq
│   │       │   │   └── atac.npz      # Predições ATAC-seq
│   │       │   └── predictions_H2/
│   │       │       ├── rna_seq.npz
│   │       │       └── atac.npz
│   │       └── rs10497191/           # Outra janela (SNP)
│   │           └── ...
│   ├── HG01880/
│   │   └── ...
│   └── ...
└── population_mapping.json            # Mapeamento FROG → 1000G (opcional)
```

### Global Metadata (`dataset_metadata.json`)

```json
{
  "dataset_name": "non_longevous_1000g",
  "creation_date": "2025-11-11T10:30:00",
  "total_individuals": 78,
  "window_size": 1000000,
  "individuals": ["HG01879", "HG01880", ...],
  "population_distribution": {
    "ACB": 3,
    "GBR": 3,
    ...
  },
  "superpopulation_distribution": {
    "AFR": 21,
    "EUR": 15,
    ...
  },
  "alphagenome_outputs": ["RNA_SEQ", "ATAC"],
  "ontologies": ["UBERON:0002107", "UBERON:0000955"],
  "frog_population_count": 150
}
```

### Individual Metadata (`individual_metadata.json`)

```json
{
  "sample_id": "HG01879",
  "sex": 1,
  "population": "ACB",
  "superpopulation": "AFR",
  "longevity": false,
  "frog_likelihood": [1.798e-56, 1.213e-57, ...],
  "frog_population_names": ["Taiwanese Han", "Ami", ...],
  "windows": ["CYP2B6", "rs10497191"],
  "window_metadata": {
    "CYP2B6": {
      "type": "gene",
      "chromosome": "chr19",
      "start": 41497445,
      "end": 42497445,
      "window_size": 1000000,
      "outputs": ["RNA_SEQ", "ATAC"],
      "ontologies": ["UBERON:0002107"]
    }
  }
}
```

---

## Installation

### Requirements

```bash
# PyTorch
pip install torch

# Pipeline dependencies (if not already installed)
pip install pandas numpy pyyaml
```

### Verification

```python
import torch
from genomic_dataset import GenomicLongevityDataset

print("✓ Imports working correctly!")
```

---

## Basic Usage

### 1. Load Dataset

```python
from pathlib import Path
from genomic_dataset import GenomicLongevityDataset

# Load dataset
dataset = GenomicLongevityDataset(
    dataset_dir=Path('non_longevous_results'),
    load_predictions=True,   # Load AlphaGenome predictions
    load_sequences=False,    # Don't load FASTA sequences (faster)
    cache_sequences=False    # Don't cache sequences in memory
)

print(f"Total individuals: {len(dataset)}")
dataset.print_summary()
```

### 2. Access a Sample

```python
# Access by index
input_data, output_data = dataset[0]

# Or by sample ID
input_data, output_data = dataset.get_sample_by_id('HG01879')

# Examine output data (labels)
print(f"Sample: {output_data['sample_id']}")
print(f"Population: {output_data['population']}")
print(f"Longevity: {output_data['longevity']}")  # 0 for non-longevous
print(f"FROG likelihood shape: {output_data['frog_likelihood'].shape}")

# Examine input data
for window_name, window_data in input_data['windows'].items():
    print(f"\nWindow: {window_name}")
    
    # H1 predictions
    for output_type, array in window_data['predictions_h1'].items():
        print(f"  H1 {output_type}: shape={array.shape}")
    
    # H2 predictions
    for output_type, array in window_data['predictions_h2'].items():
        print(f"  H2 {output_type}: shape={array.shape}")
```

### 3. Use with DataLoader

```python
from torch.utils.data import DataLoader
from genomic_dataset import collate_fn_variable_windows

# Create DataLoader
dataloader = DataLoader(
    dataset,
    batch_size=4,
    shuffle=True,
    num_workers=2,
    collate_fn=collate_fn_variable_windows  # Required for variable windows
)

# Iterate
for batch_input, batch_output in dataloader:
    # batch_input: list of dicts (one per individual)
    # batch_output: dict with lists/tensors
    
    sample_ids = batch_output['sample_id']      # List of IDs
    longevity = batch_output['longevity']        # Tensor: (batch_size,)
    populations = batch_output['population']     # List of strings
    likelihoods = batch_output['frog_likelihood'] # Tensor: (batch_size, 150)
    
    # Process batch...
    break
```

---

## Complete API

### GenomicLongevityDataset

```python
class GenomicLongevityDataset(torch.utils.data.Dataset):
    def __init__(
        self,
        dataset_dir: Path,
        load_predictions: bool = True,
        load_sequences: bool = True,
        transform: Optional[Callable] = None,
        cache_sequences: bool = False
    ):
        """
        Args:
            dataset_dir: Directory containing individuals/ and dataset_metadata.json
            load_predictions: If True, loads .npz predictions
            load_sequences: If True, loads FASTA sequences
            transform: Transformation function (input_data, output_data) -> (input_data, output_data)
            cache_sequences: If True, keeps sequences in memory (use with caution)
        """
```

#### Main Methods

- `__len__()` → int: Number of individuals
- `__getitem__(idx)` → Tuple[Dict, Dict]: Returns (input_data, output_data)
- `get_sample_by_id(sample_id)` → Tuple[Dict, Dict]: Access by ID
- `get_summary()` → Dict: Returns dataset statistics
- `print_summary()`: Prints formatted summary
- `clear_cache()`: Clears sequence cache

### Returned Data Format

#### Input Data (neural network input)

```python
input_data = {
    'windows': {
        'CYP2B6': {
            'ref_sequence': str,      # Reference sequence (if load_sequences=True)
            'h1_sequence': str,       # Haplotype 1 (if load_sequences=True)
            'h2_sequence': str,       # Haplotype 2 (if load_sequences=True)
            'predictions_h1': {
                'rna_seq': np.ndarray,  # Shape: (1000000,)
                'atac': np.ndarray      # Shape: (1000000,)
            },
            'predictions_h2': {
                'rna_seq': np.ndarray,
                'atac': np.ndarray
            }
        },
        # Other windows...
    }
}
```

#### Output Data (labels/metadados)

```python
output_data = {
    'sample_id': str,                    # E.g., 'HG01879'
    'longevity': int,                    # 0 (non-longevous) or 1 (longevous)
    'sex': int,                          # 1 (Male) or 2 (Female)
    'population': str,                   # E.g., 'ACB'
    'superpopulation': str,              # E.g., 'AFR'
    'frog_likelihood': np.ndarray,       # Shape: (150,) - likelihood per population
    'frog_population_names': List[str]   # Corresponding population names
}
```

### Collate Function

```python
def collate_fn_variable_windows(batch):
    """
    Collates batch of samples with variable number of windows.
    
    Returns:
        batch_input: List of dicts (not stacked due to variable sizes)
        batch_output: Dict with stacked fields when possible
    """
```

---

## Advanced Examples

### Custom Transformation

```python
def normalize_predictions(input_data, output_data):
    """Normalizes predictions to mean=0, std=1."""
    for window_name, window_data in input_data['windows'].items():
        for hap in ['predictions_h1', 'predictions_h2']:
            if hap in window_data:
                for output_type, array in window_data[hap].items():
                    mean, std = array.mean(), array.std()
                    if std > 0:
                        window_data[hap][output_type] = (array - mean) / std
    return input_data, output_data

# Use transformation
dataset = GenomicLongevityDataset(
    dataset_dir='non_longevous_results',
    transform=normalize_predictions
)
```

### Filter by Population

```python
class FilteredDataset(GenomicLongevityDataset):
    """Dataset filtered by population."""
    
    def __init__(self, dataset_dir, target_population, **kwargs):
        super().__init__(dataset_dir, **kwargs)
        
        # Filter individuals
        self.individuals = [
            ind for ind in self.individuals
            if self._get_population(ind) == target_population
        ]
    
    def _get_population(self, sample_id):
        metadata = self._load_individual_metadata(sample_id)
        return metadata.get('population', '')

# Use
dataset_gbr = FilteredDataset('non_longevous_results', target_population='GBR')
```

### H1 vs H2 Differences Analysis

```python
def analyze_haplotype_differences(dataset, window_name='CYP2B6'):
    """Analyzes differences between haplotypes."""
    differences = []
    
    for idx in range(len(dataset)):
        input_data, output_data = dataset[idx]
        
        if window_name in input_data['windows']:
            window = input_data['windows'][window_name]
            
            for output_type in window['predictions_h1'].keys():
                h1 = window['predictions_h1'][output_type]
                h2 = window['predictions_h2'][output_type]
                
                diff = np.abs(h1 - h2).mean()
                differences.append({
                    'sample_id': output_data['sample_id'],
                    'population': output_data['population'],
                    'output_type': output_type,
                    'mean_diff': diff
                })
    
    return pd.DataFrame(differences)

# Use
df_diffs = analyze_haplotype_differences(dataset)
print(df_diffs.groupby('population')['mean_diff'].mean())
```

### Custom Collate to Tensors

```python
def collate_to_tensors(batch):
    """
    Collate that converts predictions to stacked tensors.
    WARNING: Assumes all have same number and size of windows!
    """
    batch_input = []
    batch_output = {
        'sample_id': [],
        'longevity': [],
        'frog_likelihood': []
    }
    
    for input_data, output_data in batch:
        # Stack predictions from all windows
        all_predictions = []
        for window_data in input_data['windows'].values():
            for pred in window_data['predictions_h1'].values():
                all_predictions.append(pred)
            for pred in window_data['predictions_h2'].values():
                all_predictions.append(pred)
        
        batch_input.append(torch.tensor(np.stack(all_predictions)))
        batch_output['sample_id'].append(output_data['sample_id'])
        batch_output['longevity'].append(output_data['longevity'])
        batch_output['frog_likelihood'].append(output_data['frog_likelihood'])
    
    # Stack
    batch_input = torch.stack(batch_input)
    batch_output['longevity'] = torch.tensor(batch_output['longevity'])
    batch_output['frog_likelihood'] = torch.tensor(np.array(batch_output['frog_likelihood']))
    
    return batch_input, batch_output
```

---

## Compression and Distribution

### Compress Dataset

```bash
# Compress entire dataset
cd non_longevous_results
tar -czf ../non_longevous_dataset_v1.tar.gz individuals/ dataset_metadata.json

# Check size
ls -lh ../non_longevous_dataset_v1.tar.gz
```

### Decompress and Use

```bash
# Decompress
tar -xzf non_longevous_dataset_v1.tar.gz -C dataset_directory/

# Use
python3 -c "
from genomic_dataset import GenomicLongevityDataset
dataset = GenomicLongevityDataset('dataset_directory')
print(f'Loaded: {len(dataset)} individuals')
"
```

### Upload to Google Drive

```python
# Example with pydrive2
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive

# Authenticate
gauth = GoogleAuth()
gauth.LocalWebserverAuth()
drive = GoogleDrive(gauth)

# Upload
file = drive.CreateFile({'title': 'non_longevous_dataset_v1.tar.gz'})
file.SetContentFile('non_longevous_dataset_v1.tar.gz')
file.Upload()

print(f"Uploaded! ID: {file['id']}")
```

### Automatic Download (Future)

```python
# Planned for future version
class AutoDownloadDataset(GenomicLongevityDataset):
    """Dataset that automatically downloads from Google Drive if needed."""
    
    def __init__(self, dataset_dir, drive_file_id=None, **kwargs):
        if not dataset_dir.exists() and drive_file_id:
            self._download_from_drive(drive_file_id, dataset_dir)
        
        super().__init__(dataset_dir, **kwargs)
```

---

## FAQ

### Q: How much disk space does the dataset occupy?

**A:** Depends on number of individuals and windows:
- Per individual: ~20-50 MB
- 78 individuals (small config): ~1.5-4 GB
- Compressed: ~30-50% smaller

### Q: Can I load only predictions without sequences?

**A:** Yes! Use `load_sequences=False`:

```python
dataset = GenomicLongevityDataset(
    dataset_dir='...',
    load_sequences=False  # Faster and uses less memory
)
```

### Q: How to combine longevous and non-longevous datasets?

**A:** Use `torch.utils.data.ConcatDataset`:

```python
from torch.utils.data import ConcatDataset

dataset_non_long = GenomicLongevityDataset('non_longevous_results')
dataset_long = GenomicLongevityDataset('longevous_results')

combined = ConcatDataset([dataset_non_long, dataset_long])
```

### Q: What if individuals have different numbers of windows?

**A:** Use `collate_fn_variable_windows` which returns a list of dicts instead of stacked tensors. You'll need to process each sample individually in the batch or create a custom padding strategy.

### Q: How to handle limited memory?

**A:** 
1. Use `load_sequences=False` if you don't need sequences
2. Don't use `cache_sequences=True`
3. Use `num_workers=0` in DataLoader if you have issues
4. Process in smaller batches

### Q: Can I modify metadata after creation?

**A:** Yes, just edit the JSON files and reload the dataset:

```python
# Edit individual_metadata.json manually
# Then reload
dataset = GenomicLongevityDataset(...)
```

### Q: How to validate dataset integrity?

**A:** Use validation methods:

```python
from dataset_builder import IndividualDatasetBuilder

builder = IndividualDatasetBuilder(...)
validation = builder.validate_window_files('CYP2B6')
print(validation)  # Shows which files exist
```

---

## References

- **PyTorch Dataset Tutorial**: https://pytorch.org/tutorials/beginner/basics/data_tutorial.html
- **AlphaGenome**: Gene expression and epigenetic predictions
- **FROGAncestryCalc**: Ancestry inference via AISNPs
- **1000 Genomes Project**: http://www.internationalgenome.org/

---

## Support

For issues or questions:
1. Check this document
2. Run `python3 examples/load_dataset_example.py`
3. Consult the source code in `genomic_dataset.py`

**Author**: Alberto F. De Souza  
**Last updated**: 2025-11-11

