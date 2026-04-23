from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from datasets import load_from_disk


class HFGenomicDataset:
    def __init__(self, dataset_path: Path, genes_to_use: List[str] | None = None):
        self.dataset_path = Path(dataset_path)
        self.samples = load_from_disk(str(self.dataset_path / "samples"))
        self.gene_windows = load_from_disk(str(self.dataset_path / "gene_windows"))
        self.genes_to_use = set(genes_to_use or [])

        self.sample_rows = [self.samples[i] for i in range(len(self.samples))]
        self.sample_ids = [row["sample_id"] for row in self.sample_rows]
        self.sample_index = {sample_id: idx for idx, sample_id in enumerate(self.sample_ids)}
        self.windows_by_sample = self._index_gene_windows()
        self.dataset_metadata = self._build_dataset_metadata()

    def _index_gene_windows(self) -> Dict[str, Dict[str, Dict]]:
        windows_by_sample: Dict[str, Dict[str, Dict]] = defaultdict(dict)
        for idx in range(len(self.gene_windows)):
            row = self.gene_windows[idx]
            gene = row["gene"]
            if self.genes_to_use and gene not in self.genes_to_use:
                continue
            windows_by_sample[row["sample_id"]][gene] = row
        return dict(windows_by_sample)

    def _build_dataset_metadata(self) -> Dict:
        individuals_pedigree = {}
        population_distribution: Dict[str, int] = defaultdict(int)
        superpopulation_distribution: Dict[str, int] = defaultdict(int)

        for row in self.sample_rows:
            sample_id = row["sample_id"]
            individuals_pedigree[sample_id] = {
                "sample_id": sample_id,
                "family_id": row.get("family_id") or sample_id,
                "population": row.get("population") or "",
                "superpopulation": row.get("superpopulation") or "",
            }
            if row.get("population"):
                population_distribution[row["population"]] += 1
            if row.get("superpopulation"):
                superpopulation_distribution[row["superpopulation"]] += 1

        return {
            "individuals": list(self.sample_ids),
            "individuals_pedigree": individuals_pedigree,
            "population_distribution": dict(population_distribution),
            "superpopulation_distribution": dict(superpopulation_distribution),
        }

    def __len__(self) -> int:
        return len(self.sample_rows)

    def __getitem__(self, idx: int) -> Tuple[Dict, Dict]:
        sample_row = self.sample_rows[idx]
        sample_id = sample_row["sample_id"]
        window_rows = self.windows_by_sample.get(sample_id, {})

        input_data = {
            "windows": {
                gene: {
                    "ref_sequence": row.get("ref_sequence") or None,
                    "h1_sequence": row.get("h1_sequence") or None,
                    "h2_sequence": row.get("h2_sequence") or None,
                    "predictions_h1": _decode_predictions(row.get("predictions_h1", {})),
                    "predictions_h2": _decode_predictions(row.get("predictions_h2", {})),
                }
                for gene, row in window_rows.items()
            }
        }

        output_data = {
            "sample_id": sample_id,
            "family_id": sample_row.get("family_id") or sample_id,
            "longevity": 1 if sample_row.get("longevity") else 0,
            "sex": sample_row.get("sex", -1),
            "population": sample_row.get("population") or "",
            "superpopulation": sample_row.get("superpopulation") or "",
            "frog_likelihood": np.asarray(sample_row.get("frog_likelihood", []), dtype=np.float32),
            "frog_population_names": sample_row.get("frog_population_names", []),
        }
        return input_data, output_data


def _decode_predictions(predictions: Dict) -> Dict[str, np.ndarray]:
    decoded = {}
    for output_name, matrix in predictions.items():
        if not matrix:
            continue
        array = np.asarray(matrix, dtype=np.float32)
        if array.ndim == 2 and array.shape[1] == 1:
            array = array[:, 0]
        decoded[output_name] = array
    return decoded
