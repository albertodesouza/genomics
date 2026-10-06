"""Owns the one config + runtime dataset + loaded CNN2 checkpoint an interactive app needs to
score an individual's pigmentation prediction, optionally with edited track data substituted in
for one or more (gene, haplotype) pairs.

Promotes notebooks/genotype_cnn_alignment_deeplift.ipynb's ``build_individual_tensor``/
``run_cnn_logits`` (Section 9) into a reusable class that reads its gene list from the loaded
config instead of a notebook global, and is meant to be constructed exactly once per server
process (loading the checkpoint and building the runtime dataset are both non-trivial-cost
operations -- see ``experiments/evaluate_checkpoint.py`` for the equivalent one-shot CLI pattern).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


class PigmentationModelContext:
    def __init__(self, config_path: Union[str, Path], checkpoint: str = "best_accuracy", device: Optional[str] = None):
        import torch

        from genomics.predictors.genotype_based.config import (
            generate_experiment_name,
            get_dataset_cache_dir,
            get_experiment_runs_dir,
            load_config,
        )
        from genomics.predictors.genotype_based.data.normalization import apply_normalization
        from genomics.predictors.genotype_based.data.pipeline import _make_runtime_processed_datasets, _resolve_runtime_dataset_dir
        from genomics.predictors.genotype_based.models import CNN2AncestryPredictor

        self._torch = torch
        self._apply_normalization = apply_normalization

        self.config = load_config(Path(config_path))
        if self.config.model.type.upper() != "CNN2":
            raise ValueError(
                f"PigmentationModelContext so suporta model.type='CNN2', config tem {self.config.model.type!r}"
            )
        if self.config.dataset_input.tensor_layout != "haplotype_channels":
            raise ValueError("PigmentationModelContext requer dataset_input.tensor_layout='haplotype_channels'")

        runtime_dataset_dir = _resolve_runtime_dataset_dir(self.config)
        cache_dir = get_dataset_cache_dir(self.config)
        self.full_ds, _train_ds, _val_ds, _test_ds = _make_runtime_processed_datasets(
            runtime_dataset_dir, cache_dir, self.config
        )
        self.normalization_params = self.full_ds.normalization_params

        self.device = torch.device(device) if device else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        experiment_dir = get_experiment_runs_dir(self.config) / generate_experiment_name(self.config)
        checkpoint_path = self._resolve_checkpoint(experiment_dir, checkpoint)
        if not checkpoint_path.exists():
            raise FileNotFoundError(f"Checkpoint nao encontrado: {checkpoint_path}")

        self.model = CNN2AncestryPredictor(
            self.config, self.full_ds.get_input_shape(), self.full_ds.get_num_classes()
        ).to(self.device)
        state = torch.load(checkpoint_path, map_location=self.device)
        self.model.load_state_dict(state.get("model_state_dict", state))
        self.model.eval()

        self.class_names: List[str] = self.full_ds.get_class_names()
        if "strong pigmentation" not in self.class_names:
            raise ValueError(f"Classe 'strong pigmentation' nao encontrada em {self.class_names}")
        self.strong_pigmentation_idx = self.class_names.index("strong pigmentation")

        self._baseline_tensor_cache: Dict[str, "torch.Tensor"] = {}

    @staticmethod
    def _resolve_checkpoint(experiment_dir: Path, checkpoint: str) -> Path:
        path = Path(checkpoint)
        if path.exists():
            return path.resolve()
        if path.suffix != ".pt":
            path = path.with_suffix(".pt")
        return experiment_dir / "models" / path.name

    @property
    def genes(self) -> List[str]:
        return list(self.config.dataset_input.genes_to_use or [])

    @property
    def rna_seq_ontology_terms(self) -> List[str]:
        return list(self.config.dataset_input.ontology_terms or [])

    @property
    def window_center_size(self) -> int:
        return int(self.config.dataset_input.window_center_size)

    def load_raw_prediction(self, sample_id: str, gene: str, haplotype: str, output: str = "rna_seq") -> Tuple[np.ndarray, list]:
        dataset_dir = Path(self.config.dataset_input.dataset_dir)
        pred_dir = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
        values = np.load(pred_dir / f"{output}.npz")["values"]
        import json

        metadata = json.loads((pred_dir / f"{output}_metadata.json").read_text())["metadata"]
        return values, metadata

    def build_individual_tensor(
        self,
        sample_id: str,
        overrides: Optional[Dict[Tuple[str, str], Tuple[np.ndarray, list]]] = None,
    ) -> "torch.Tensor":
        """Rebuild the CNN2 input tensor for ``sample_id`` across all of the config's genes,
        substituting ``overrides[(gene, haplotype)] = (raw_array, track_metadata)`` in place of
        the on-disk cached RNA-seq prediction for that gene/haplotype -- every other gene/haplotype
        uses the individual's real, unmodified prediction."""
        overrides = overrides or {}
        h1_rows, h2_rows = [], []
        for gene in self.genes:
            for haplotype, rows in (("H1", h1_rows), ("H2", h2_rows)):
                if (gene, haplotype) in overrides:
                    array, meta = overrides[(gene, haplotype)]
                else:
                    array, meta = self.load_raw_prediction(sample_id, gene, haplotype)
                result = self.full_ds._process_window_haplotype_channels(
                    sample_id, gene, haplotype, {"rna_seq": array}, {"rna_seq": meta}
                )
                if result is None:
                    raise RuntimeError(f"nao foi possivel construir tensor para {sample_id}/{gene}/{haplotype}")
                signals, _masks = result
                rows.append(signals)
        stacked = np.stack([np.concatenate(h1_rows, axis=0), np.concatenate(h2_rows, axis=0)], axis=0)
        return self._apply_normalization(self._torch.from_numpy(stacked), self.normalization_params)

    def baseline_tensor(self, sample_id: str) -> "torch.Tensor":
        if sample_id not in self._baseline_tensor_cache:
            self._baseline_tensor_cache[sample_id] = self.build_individual_tensor(sample_id)
        return self._baseline_tensor_cache[sample_id]

    def run_cnn_logits(self, tensor: "torch.Tensor") -> "torch.Tensor":
        with self._torch.no_grad():
            return self.model(tensor.unsqueeze(0).float().to(self.device))[0]

    def predict_proba(self, tensor: "torch.Tensor") -> "torch.Tensor":
        with self._torch.no_grad():
            return self._torch.softmax(self.run_cnn_logits(tensor), dim=0)

    def strong_pigmentation_probability(self, tensor: "torch.Tensor") -> float:
        return float(self.predict_proba(tensor)[self.strong_pigmentation_idx].item())
