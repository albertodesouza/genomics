"""Dataset reader for the on-disk layout consumed by this package."""

from __future__ import annotations

import json
import importlib.util
from pathlib import Path
from typing import Any, Dict


def _load_genomic_longevity_dataset_class():
    repo_root = Path(__file__).resolve().parent.parent
    genomic_dataset_path = repo_root / "build_non_longevous_dataset" / "genomic_dataset.py"
    if not genomic_dataset_path.exists():
        raise FileNotFoundError(f"genomic_dataset.py não encontrado: {genomic_dataset_path}")

    spec = importlib.util.spec_from_file_location("build_non_longevous_dataset.genomic_dataset", genomic_dataset_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Não foi possível carregar módulo: {genomic_dataset_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.GenomicLongevityDataset


GenomicLongevityDataset = _load_genomic_longevity_dataset_class()


class GenomicDataset(GenomicLongevityDataset):
    """Compatibility reader for the package dataset layout."""

    def _load_window_data(self, sample_id: str, window_name: str) -> Dict:
        window_dir = self.dataset_dir / "individuals" / sample_id / "windows" / window_name
        ref_window_dir = self.dataset_dir / "references" / "windows" / window_name

        window_data: Dict[str, Any] = {}

        if self.load_sequences:
            ref_fasta = ref_window_dir / "ref.window.fa"
            h1_fasta = window_dir / f"{sample_id}.H1.window.fixed.fa"
            h2_fasta = window_dir / f"{sample_id}.H2.window.fixed.fa"

            try:
                window_data["ref_sequence"] = self._load_fasta_sequence(ref_fasta)
            except FileNotFoundError:
                window_data["ref_sequence"] = None

            try:
                window_data["h1_sequence"] = self._load_fasta_sequence(h1_fasta)
            except FileNotFoundError:
                window_data["h1_sequence"] = None

            try:
                window_data["h2_sequence"] = self._load_fasta_sequence(h2_fasta)
            except FileNotFoundError:
                window_data["h2_sequence"] = None

        if self.load_predictions:
            window_data["predictions_h1"] = self._load_predictions(window_dir / "predictions_H1")
            window_data["predictions_h2"] = self._load_predictions(window_dir / "predictions_H2")

        meta_path = ref_window_dir / "window_metadata.json"
        if meta_path.exists():
            with open(meta_path) as f:
                window_data["window_metadata"] = json.load(f)

        return window_data
