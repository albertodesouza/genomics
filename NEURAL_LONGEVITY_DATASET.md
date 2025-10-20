# Neural Longevity Dataset Builder

Comprehensive usage guide for `neural_longevity_dataset.py`, the CLI pipeline that assembles training-ready datasets for longevity studies by combining public genomes, AlphaGenome predictions, and PyTorch-friendly serialization.

## 1. Overview

`neural_longevity_dataset.py` automates the following stages:

1. Discover and download high-coverage whole-genome runs (e.g., 1000 Genomes Project) from the ENA portal.
2. Call individual variants against GRCh38 using `bcftools` (CRAM → VCF).
3. Select central variant positions according to configurable filters and strategies.
4. Extract FASTA windows centred on each variant while applying the alternate allele sequence.
5. Run AlphaGenome on every FASTA and collect the predicted feature vectors.
6. Build balanced PyTorch dataset splits (`train`, `val`, `test`) with cached tensors and metadata.

All persistent data **must** live under `/dados/GENOMICS_DATA/top3/`, and the script must be executed from that directory. Source code and configuration files stay in the repository at `~/genomics`.

## 2. Prerequisites

- Linux workstation with Python 3.10+.
- System tools: `wget`, `samtools`, `bcftools`, `tabix`, `pigz`.
- AlphaGenome CLI installed and accessible in `$PATH`.
- PyTorch, NumPy, pandas, cyvcf2, and pyfaidx installed in the active environment.
- ENA access over HTTPS (no authentication required for public runs).

## 3. Directory Layout Requirement

The dataset root is fixed to `/dados/GENOMICS_DATA/top3/`. The script enforces this at runtime. Example layout after the pipeline runs:

```
/dados/GENOMICS_DATA/top3/
├── longevity_dataset/              # project.output_dir from the YAML
│   ├── checkpoint.json
│   ├── cram/
│   ├── vcf/
│   ├── windows/
│   ├── alphagenome/
│   └── torch_dataset/
├── alphagenome_cache/              # if alphagenome.cache_dir is relative
└── ...
```

> **Important:** If the directory does not exist, the script will create it (subject to permissions). Always run the CLI from inside `/dados/GENOMICS_DATA/top3/` so relative paths resolve correctly.

## 4. Configuration File

The YAML configuration stays in the repository (e.g., `~/genomics/longevity_config.yaml`). Below is a minimal template illustrating the required sections:

```yaml
project:
  name: "Longevity Dataset"
  description: "Prototype build using 1000G high coverage"
  output_dir: "longevity_dataset"   # created inside /dados/GENOMICS_DATA/top3/

data_sources:
  reference:
    fasta: "refs/reference.fa"      # relative to the YAML file directory
  longevous:
    ena_project: "PRJEB31736"
    sample_range: [0, 10]
  non_longevous:
    ena_project: "PRJEB31736"
    sample_range: [10, 20]

sequence_extraction:
  window_size: 1000
  use_alternate_allele: true
  center_on_variant: true

alphagenome:
  api_key: "${ALPHAGENOME_API_KEY}"
  outputs: ["functional", "epigenetic"]
  cache_dir: "alphagenome_cache"    # stored under /dados/GENOMICS_DATA/top3/
  cache_results: true

dataset:
  random_seed: 1234
  splits: {train: 0.7, val: 0.15, test: 0.15}
  balance_classes: true
  features:
    sequence: {one_hot: true}
    position: {normalize: true}
    alphagenome_predictions:
      statistics: ["mean", "std", "max", "min"]
    metadata:
      fields: ["sample_id", "variant_type", "quality"]

pipeline:
  steps: [
    "download_samples",
    "select_central_points",
    "extract_sequences",
    "run_alphagenome",
    "build_dataset"
  ]

debug:
  dry_run: false
```

### Path Resolution Rules

- Relative paths in the YAML are resolved against the directory that contains the YAML file (inside `~/genomics`).
- `project.output_dir` **must** remain inside `/dados/GENOMICS_DATA/top3/`. The script aborts otherwise.
- `alphagenome.cache_dir` is resolved inside `/dados/GENOMICS_DATA/top3/` when relative.

## 5. Quickstart

1. Ensure `/dados/GENOMICS_DATA/top3/` exists and has the necessary permissions.
2. Change directory:
   ```bash
   cd /dados/GENOMICS_DATA/top3/
   ```
3. Execute the builder, pointing to the YAML stored in the repo:
   ```bash
   python3 ~/genomics/neural_longevity_dataset.py \
     --config ~/genomics/longevity_config.yaml
   ```
4. Monitor the Rich progress panels. On completion, review the outputs under `/dados/GENOMICS_DATA/top3/<output_dir>/`.

## 6. Running Specific Steps

Use `--steps` to run a subset of the pipeline:

```bash
python3 ~/genomics/neural_longevity_dataset.py \
  --config ~/genomics/longevity_config.yaml \
  --steps download_samples extract_sequences
```

Use `--dry-run` to preview actions without creating/modifying files.

## 7. Output Artifacts

| Directory | Contents |
|-----------|----------|
| `cram/` | Downloaded CRAM/CRAI pairs from ENA. |
| `vcf/` | Variant calls per sample (`.vcf.gz` + `.tbi`). |
| `windows/` | FASTA windows, metadata CSVs, and sequence index. |
| `alphagenome/` | Cached AlphaGenome predictions (`.pkl`, JSON). |
| `torch_dataset/` | `train.pkl`, `val.pkl`, `test.pkl`, plus `samples.csv` summary. |

Additional logs: `checkpoint.json`, per-step runtime metrics, and state caches for resumable execution.

## 8. Troubleshooting

- **RuntimeError: must be executed from /dados/GENOMICS_DATA/top3** – change directory to that path before launching the script.
- **Reference not found** – verify `data_sources.reference.fasta` resolves correctly relative to the YAML file.
- **AlphaGenome authentication errors** – confirm the CLI is installed and `alphagenome.api_key` is valid.
- **Permissions denied when creating directories** – adjust file-system permissions on `/dados/GENOMICS_DATA/top3/` or run with sufficient privileges.

## 9. Updating the Dataset

To refresh downloads or rerun specific stages:

1. Delete the relevant subdirectory inside the project output folder (e.g., remove `vcf/` to force variant recalling).
2. Rerun the script with the desired steps; the checkpoint will skip completed phases unless you remove or edit `checkpoint.json`.

## 10. Integration with PyTorch

After execution, load `torch_dataset/train.pkl` (and `val.pkl`, `test.pkl`) using `LongevityDataset` defined in `neural_longevity_dataset.py`. These pickles contain:

- `sequences`: one-hot encoded tensors (NumPy arrays).
- `positions`: normalized genomic coordinates.
- `alphagenome_features`: aggregated AlphaGenome statistics per sample.
- `labels`: binary class labels (longevous vs. control).
- `metadata`: dictionaries with provenance (sample ID, variant info, QC metrics).

Example snippet:

```python
from neural_longevity_dataset import LongevityDataset

dataset = LongevityDataset("/dados/GENOMICS_DATA/top3/longevity_dataset/torch_dataset/train.pkl")
print(dataset.get_class_distribution())
```

## 11. Maintenance Notes

- Keep the script and configs in `~/genomics` under version control.
- Document environment versions (`bcftools`, `samtools`, AlphaGenome) to ensure reproducibility.
- Review `checkpoint.json` and log files after failures to resume safely.

Happy dataset building!
