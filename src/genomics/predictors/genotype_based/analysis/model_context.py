"""A trained genotype_based run (config + checkpoint) loaded once for interactive scoring.

Scores an individual through exactly the path used to build training tensors -- the run's window
processing (alignment mapping, feature mode, masks) and saved normalization -- optionally with
edited AlphaGenome predictions substituted for some ``(gene, haplotype)`` windows. Unlike
:class:`~genomics.predictors.genotype_based.analysis.pigmentation_model_context.PigmentationModelContext`
it is not tied to one target, model type or tensor layout: any torch model trained by
``genomics genotype train`` works, and every class probability is reported.
"""
from __future__ import annotations

import threading
from collections import OrderedDict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

Overrides = Dict[Tuple[str, str], Dict[str, Tuple[np.ndarray, Optional[list]]]]
ProgressFn = Callable[[float, str], None]


def _absolute(path: Optional[str], base: Path) -> Optional[str]:
    if not path:
        return path
    candidate = Path(path).expanduser()
    return str(candidate if candidate.is_absolute() else (base / candidate).resolve())


def resolve_checkpoint(experiment_dir: Path, checkpoint: str) -> Path:
    path = Path(checkpoint)
    if path.is_absolute() and path.exists():
        return path
    if path.suffix != ".pt":
        path = path.with_suffix(".pt")
    return Path(experiment_dir) / "models" / path.name


class GenotypeModelContext:
    def __init__(
        self,
        config_path: Union[str, Path],
        checkpoint: str = "best_accuracy",
        experiment_dir: Optional[Union[str, Path]] = None,
        device: Optional[str] = None,
        progress: Optional[ProgressFn] = None,
        base_dir: Optional[Path] = None,
    ):
        import torch

        from genomics.predictors.genotype_based.config import (
            generate_experiment_name,
            get_dataset_cache_dir,
            get_experiment_runs_dir,
            load_config,
        )
        from genomics.predictors.genotype_based.data.pipeline import (
            _load_split_index,
            _make_runtime_processed_datasets,
            _resolve_runtime_dataset_dir,
        )
        from genomics.predictors.genotype_based.experiments.train import _build_model
        from genomics.predictors.genotype_based.models import SKLEARN_BASELINE_TYPES

        progress = progress or (lambda *_: None)
        self._torch = torch
        self._lock = threading.RLock()
        progress(0.02, "Reading the run configuration")
        self.config_path = Path(config_path)
        self.config = load_config(self.config_path)
        if self.config.model.type.upper() in SKLEARN_BASELINE_TYPES:
            raise ValueError(f"{self.config.model.type} is a scikit-learn baseline; the lab needs a torch model (NN, CNN or CNN2)")
        if self.config.output.prediction_target == "frog_likelihood":
            raise ValueError("Regression targets (frog_likelihood) are not supported")
        # Relative paths in run configs are relative to the repository root (where training ran).
        if base_dir is None:
            from genomics.workspace import repo_root

            base_dir = repo_root()
        di = self.config.dataset_input
        di.processed_cache_dir = _absolute(di.processed_cache_dir, base_dir)
        di.results_dir = _absolute(di.results_dir, base_dir)
        if getattr(di, "consensus_dataset_dir", None):
            di.consensus_dataset_dir = _absolute(di.consensus_dataset_dir, base_dir)

        self.dataset_dir = _resolve_runtime_dataset_dir(self.config)
        cache_dir = get_dataset_cache_dir(self.config)
        if not (cache_dir / "normalization_params.json").exists():
            raise FileNotFoundError(f"No processed cache for this run at {cache_dir} (normalization parameters are needed to score samples)")
        progress(0.1, "Building the runtime dataset and alignment (first use can take minutes)")
        self.full_ds, _train, _val, _test = _make_runtime_processed_datasets(self.dataset_dir, cache_dir, self.config)
        try:
            self.split_index: Dict[str, List[str]] = _load_split_index(cache_dir)
        except (OSError, ValueError):
            self.split_index = {}
        self._split_of = {str(s): split for split, ids in self.split_index.items() if isinstance(ids, list) for s in ids}

        progress(0.85, "Loading the checkpoint")
        self.experiment_dir = Path(experiment_dir) if experiment_dir else get_experiment_runs_dir(self.config) / generate_experiment_name(self.config)
        self.checkpoint_path = resolve_checkpoint(self.experiment_dir, checkpoint)
        if not self.checkpoint_path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {self.checkpoint_path}")
        self.device = torch.device(device) if device else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model = _build_model(self.config, self.full_ds).to(self.device)
        state = torch.load(self.checkpoint_path, map_location=self.device)
        self.model.load_state_dict(state.get("model_state_dict", state) if isinstance(state, dict) else state)
        self.model.eval()
        self.class_names: List[str] = [str(c) for c in self.full_ds.get_class_names()]
        self._windows: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
        self._baseline: "OrderedDict[str, np.ndarray]" = OrderedDict()
        progress(1.0, "Model ready")

    # -- description ----------------------------------------------------------------------------
    @property
    def genes(self) -> List[str]:
        return list(self.config.dataset_input.genes_to_use or [])

    @property
    def outputs(self) -> List[str]:
        return [str(o).lower() for o in self.config.dataset_input.alphagenome_outputs or []]

    @property
    def window_center_size(self) -> int:
        return int(self.config.dataset_input.window_center_size)

    @property
    def target(self) -> str:
        return str(self.config.output.prediction_target)

    def label(self, sample_id: str) -> Optional[str]:
        pedigree = (self.full_ds.dataset_metadata.get("individuals_pedigree") or {}).get(sample_id) or {}
        try:
            return self.full_ds._get_target_value(pedigree)
        except ValueError:
            return None

    def split_of(self, sample_id: str) -> Optional[str]:
        return self._split_of.get(sample_id)

    def labels(self) -> Dict[str, Optional[str]]:
        individuals = [str(s) for s in self.full_ds.dataset_metadata.get("individuals") or []]
        return {s: self.label(s) for s in individuals}

    # -- scoring --------------------------------------------------------------------------------
    def _sample_windows(self, sample_id: str) -> Dict[str, Any]:
        hit = self._windows.get(sample_id)
        if hit is not None:
            self._windows.move_to_end(sample_id)
            return hit
        base = self.full_ds.base_dataset
        windows = {}
        for gene in self.genes:
            data = base._load_window_data(sample_id, gene)
            data.setdefault("window_metadata", {}).setdefault("sample_id", sample_id)
            windows[gene] = data
        self._windows[sample_id] = windows
        while len(self._windows) > 3:
            self._windows.popitem(last=False)
        return windows

    def tensor(self, sample_id: str, overrides: Optional[Overrides] = None):
        """Normalized model input for ``sample_id``; ``overrides[(gene, "H1")] = {output: (array, track_metadata)}``."""
        windows = self._sample_windows(sample_id)
        if overrides:
            windows = {gene: dict(data) for gene, data in windows.items()}
            for (gene, haplotype), outputs in overrides.items():
                if gene not in windows:
                    continue
                suffix = haplotype.lower()
                preds = dict(windows[gene].get(f"predictions_{suffix}") or {})
                metas = dict(windows[gene].get(f"prediction_metadata_{suffix}") or {})
                for output, (array, meta) in outputs.items():
                    preds[output] = array
                    if meta is not None:
                        metas[output] = meta
                windows[gene][f"predictions_{suffix}"] = preds
                windows[gene][f"prediction_metadata_{suffix}"] = metas
        features = self.full_ds._process_windows(windows, sample_id=sample_id)
        if features.size == 0:
            raise RuntimeError(f"No model input could be built for {sample_id} (missing windows for {', '.join(self.genes)})")
        return self.full_ds._normalize_features_tensor(self._torch.FloatTensor(features))

    def predict_proba(self, tensor) -> np.ndarray:
        with self._lock, self._torch.no_grad():
            logits = self.model(tensor.unsqueeze(0).float().to(self.device))[0]
            return self._torch.softmax(logits, dim=0).cpu().numpy().astype(np.float64)

    def baseline(self, sample_id: str) -> np.ndarray:
        hit = self._baseline.get(sample_id)
        if hit is None:
            hit = self.predict_proba(self.tensor(sample_id))
            self._baseline[sample_id] = hit
            while len(self._baseline) > 256:
                self._baseline.popitem(last=False)
        return hit

    def score(self, sample_id: str, overrides: Optional[Overrides] = None) -> np.ndarray:
        if not overrides:
            return self.baseline(sample_id)
        return self.predict_proba(self.tensor(sample_id, overrides))

    def raw_prediction(self, sample_id: str, gene: str, haplotype: str, output: str) -> Tuple[np.ndarray, Optional[list]]:
        windows = self._sample_windows(sample_id)
        data = windows.get(gene) or {}
        suffix = haplotype.lower()
        array = (data.get(f"predictions_{suffix}") or {}).get(output)
        meta = (data.get(f"prediction_metadata_{suffix}") or {}).get(output)
        if array is None:
            raise FileNotFoundError(f"No {output} prediction for {sample_id}/{gene}/{haplotype}")
        return array, meta
