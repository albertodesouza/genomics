from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset, Subset

from .config import PredictionConfig
from .hf_dataset import HFGenomicDataset


def zscore_normalize(tensor: torch.Tensor, mean: float, std: float) -> torch.Tensor:
    std = 1.0 if std < 1e-8 else std
    return (tensor - mean) / std


def log_normalize(tensor: torch.Tensor, log_max: float) -> torch.Tensor:
    log_max = 1.0 if log_max <= 0 else log_max
    return torch.log1p(torch.clamp(tensor, min=0.0)) / log_max


def minmax_keep_zero(tensor: torch.Tensor, xmax: float) -> torch.Tensor:
    xmax = 1.0 if xmax <= 0 else xmax
    return torch.where(tensor == 0, tensor, tensor / xmax)


@dataclass
class DataBundle:
    dataset: "ProcessedPredictionDataset"
    train_loader: DataLoader
    val_loader: DataLoader
    test_loader: DataLoader


class ProcessedPredictionDataset(Dataset):
    def __init__(self, base_dataset: HFGenomicDataset, config: PredictionConfig, normalization_params: Optional[Dict] = None):
        self.base_dataset = base_dataset
        self.config = config
        self.alphagenome_outputs = config.dataset_input.alphagenome_outputs
        self.haplotype_mode = config.dataset_input.haplotype_mode
        self.window_center_size = int(config.dataset_input.window_center_size)
        self.downsample_factor = int(config.dataset_input.downsample_factor)
        self.prediction_target = config.output.prediction_target
        self.derived_targets = config.output.derived_targets
        self.normalization_method = config.dataset_input.normalization_method
        self.target_to_idx: Dict[str, int] = {}
        self.idx_to_target: Dict[int, str] = {}

        self._create_target_mappings()
        self.valid_sample_indices = self._discover_valid_sample_indices()
        self.normalization_params = normalization_params or self._compute_normalization_params(self.valid_sample_indices)

    def _create_target_mappings(self) -> None:
        known_classes = self.config.output.known_classes or []
        if known_classes:
            ordered = list(known_classes)
        else:
            values = []
            for idx in range(len(self.base_dataset)):
                _, output_data = self.base_dataset[idx]
                value = self._get_target_value(output_data)
                if value is not None:
                    values.append(value)
            ordered = sorted(set(values))
        self.target_to_idx = {value: idx for idx, value in enumerate(ordered)}
        self.idx_to_target = {idx: value for value, idx in self.target_to_idx.items()}

    def _discover_valid_sample_indices(self) -> List[int]:
        valid = []
        for idx in range(len(self.base_dataset)):
            _, output_data = self.base_dataset[idx]
            value = self._get_target_value(output_data)
            if value in self.target_to_idx:
                valid.append(idx)
        return valid

    def _get_target_value(self, output_data: Dict) -> Optional[str]:
        if self.prediction_target == "frog_likelihood":
            return None
        if self.prediction_target == "population":
            return output_data.get("population")
        if self.prediction_target == "superpopulation":
            return output_data.get("superpopulation")
        if self.prediction_target in self.derived_targets:
            derived = self.derived_targets[self.prediction_target]
            source_value = output_data.get(derived.source_field)
            if source_value is None:
                return None
            for class_name, source_values in derived.class_map.items():
                if source_value in source_values:
                    return class_name
            if derived.exclude_unmapped:
                return None
            return source_value
        raise ValueError(f"Unsupported prediction_target: {self.prediction_target}")

    def _process_windows(self, windows: Dict) -> np.ndarray:
        rows = []
        for gene in self.config.dataset_input.genes_to_use:
            window_data = windows.get(gene)
            if not window_data:
                continue
            if self.haplotype_mode in {"H1", "H1+H2"}:
                h1_rows = self._process_haplotype(window_data.get("predictions_h1", {}))
                if h1_rows is not None:
                    rows.append(h1_rows if h1_rows.ndim == 2 else h1_rows.reshape(1, -1))
            if self.haplotype_mode in {"H2", "H1+H2"}:
                h2_rows = self._process_haplotype(window_data.get("predictions_h2", {}))
                if h2_rows is not None:
                    rows.append(h2_rows if h2_rows.ndim == 2 else h2_rows.reshape(1, -1))
        if not rows:
            return np.zeros((0, 0), dtype=np.float32)
        return np.vstack(rows)

    def _process_haplotype(self, predictions: Dict[str, np.ndarray]) -> Optional[np.ndarray]:
        rows = []
        for output_name in self.alphagenome_outputs:
            if output_name not in predictions:
                continue
            array = predictions[output_name]
            if array.ndim == 2 and array.shape[1] > 1:
                for track_idx in range(array.shape[1]):
                    rows.append(self._downsample(self._extract_center(array[:, track_idx])))
            else:
                flattened = array.flatten() if array.ndim > 1 else array
                rows.append(self._downsample(self._extract_center(flattened)))
        if not rows:
            return None
        return rows[0] if len(rows) == 1 else np.vstack(rows)

    def _extract_center(self, array: np.ndarray) -> np.ndarray:
        if self.window_center_size <= 0 or self.window_center_size >= len(array):
            return array
        center = len(array) // 2
        half = self.window_center_size // 2
        return array[max(0, center - half): min(len(array), center + half)]

    def _downsample(self, array: np.ndarray) -> np.ndarray:
        if self.downsample_factor <= 1:
            return array
        return array[:: self.downsample_factor]

    def _compute_normalization_params(self, indices: List[int]) -> Dict:
        track_values = None
        for base_idx in indices:
            input_data, _ = self.base_dataset[base_idx]
            processed = self._process_windows(input_data["windows"])
            if processed.size == 0:
                continue
            if track_values is None:
                track_values = [[] for _ in range(processed.shape[0])]
            for track_idx in range(processed.shape[0]):
                track_values[track_idx].append(processed[track_idx])

        if not track_values:
            return {"method": self.normalization_method, "per_track": False, "mean": 0.0, "std": 1.0}

        track_params = []
        for values in track_values:
            merged = np.concatenate(values)
            if self.normalization_method == "log":
                nonzero = merged[merged > 0]
                log_max = float(np.log1p(nonzero.max())) if len(nonzero) else 1.0
                track_params.append({"log_max": log_max})
            elif self.normalization_method == "minmax_keep_zero":
                nonzero = merged[merged > 0]
                xmax = float(nonzero.max()) if len(nonzero) else 1.0
                track_params.append({"max": xmax})
            else:
                mean = float(np.mean(merged))
                std = float(np.std(merged))
                track_params.append({"mean": mean, "std": 1.0 if std < 1e-8 else std})
        return {
            "method": self.normalization_method,
            "per_track": True,
            "track_params": track_params,
        }

    def __len__(self) -> int:
        return len(self.valid_sample_indices)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        base_idx = self.valid_sample_indices[idx]
        input_data, output_data = self.base_dataset[base_idx]
        features = torch.tensor(self._process_windows(input_data["windows"]), dtype=torch.float32)
        if self.normalization_params.get("per_track"):
            rows = []
            for track_idx in range(features.shape[0]):
                track = features[track_idx: track_idx + 1, :]
                params = self.normalization_params["track_params"][track_idx]
                if self.normalization_method == "log":
                    rows.append(log_normalize(track, params["log_max"]))
                elif self.normalization_method == "minmax_keep_zero":
                    rows.append(minmax_keep_zero(track, params["max"]))
                else:
                    rows.append(zscore_normalize(track, params["mean"], params["std"]))
            features = torch.cat(rows, dim=0)

        if self.prediction_target == "frog_likelihood":
            target = torch.tensor(output_data.get("frog_likelihood", []), dtype=torch.float32)
        else:
            target_value = self._get_target_value(output_data)
            target = torch.tensor(self.target_to_idx[target_value], dtype=torch.long)
        return features, target

    def get_input_shape(self) -> Tuple[int, int]:
        for idx in range(len(self)):
            features, _ = self[idx]
            if features.numel() > 0:
                return tuple(features.shape)
        raise ValueError("No valid features found in dataset")

    def get_num_classes(self) -> int:
        if self.prediction_target == "frog_likelihood":
            _, output_data = self.base_dataset[self.valid_sample_indices[0]]
            return int(len(output_data.get("frog_likelihood", [])))
        return len(self.target_to_idx)


def get_dataset_cache_dir(config: PredictionConfig) -> Optional[Path]:
    if not config.dataset_input.processed_cache_dir:
        return None
    return Path(config.dataset_input.processed_cache_dir) / generate_dataset_name(config)


def generate_dataset_name(config: PredictionConfig) -> str:
    fingerprint = {
        "dataset": config.dataset_input.hf_dataset_path,
        "outputs": config.dataset_input.alphagenome_outputs,
        "haplotype_mode": config.dataset_input.haplotype_mode,
        "window_center_size": config.dataset_input.window_center_size,
        "downsample_factor": config.dataset_input.downsample_factor,
        "genes_to_use": config.dataset_input.genes_to_use,
        "normalization_method": config.dataset_input.normalization_method,
        "prediction_target": config.output.prediction_target,
    }
    digest = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    return f"dataset_{digest}"


def save_processed_cache(cache_dir: Path, dataset: ProcessedPredictionDataset, train_indices: List[int], val_indices: List[int], test_indices: List[int]) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "normalization_params": dataset.normalization_params,
        "train_indices": train_indices,
        "val_indices": val_indices,
        "test_indices": test_indices,
        "valid_sample_indices": dataset.valid_sample_indices,
    }
    with (cache_dir / "cache.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def load_processed_cache(cache_dir: Path) -> Dict:
    with (cache_dir / "cache.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_family_splits(base_dataset: HFGenomicDataset, processed_dataset: ProcessedPredictionDataset, config: PredictionConfig) -> Tuple[List[int], List[int], List[int]]:
    pedigree = base_dataset.dataset_metadata.get("individuals_pedigree", {})
    groups: Dict[str, List[int]] = {}
    valid_set = set(processed_dataset.valid_sample_indices)
    for base_idx, sample_id in enumerate(base_dataset.sample_ids):
        if base_idx not in valid_set:
            continue
        family_id = pedigree.get(sample_id, {}).get("family_id") or sample_id
        groups.setdefault(family_id, []).append(base_idx)

    family_groups = list(groups.values())
    rng = np.random.default_rng(config.data_split.random_seed)
    rng.shuffle(family_groups)

    total_groups = len(family_groups)
    train_count = int(total_groups * config.data_split.train_split)
    val_count = int(total_groups * config.data_split.val_split)
    train_groups = family_groups[:train_count]
    val_groups = family_groups[train_count: train_count + val_count]
    test_groups = family_groups[train_count + val_count:]

    processed_idx_by_base_idx = {
        base_idx: processed_idx for processed_idx, base_idx in enumerate(processed_dataset.valid_sample_indices)
    }

    def flatten(group_list: List[List[int]]) -> List[int]:
        return [processed_idx_by_base_idx[base_idx] for group in group_list for base_idx in group]

    return flatten(train_groups), flatten(val_groups), flatten(test_groups)


def create_dataloaders(config: PredictionConfig) -> DataBundle:
    base_dataset = HFGenomicDataset(
        dataset_path=config.dataset_input.hf_dataset_path,
        genes_to_use=config.dataset_input.genes_to_use,
    )
    cache_dir = get_dataset_cache_dir(config)
    cached = load_processed_cache(cache_dir) if cache_dir is not None and (cache_dir / "cache.json").exists() else None

    processed_dataset = ProcessedPredictionDataset(
        base_dataset=base_dataset,
        config=config,
        normalization_params=cached.get("normalization_params") if cached else None,
    )
    if cached:
        train_indices = cached["train_indices"]
        val_indices = cached["val_indices"]
        test_indices = cached["test_indices"]
    else:
        train_indices, val_indices, test_indices = build_family_splits(base_dataset, processed_dataset, config)
        if cache_dir is not None:
            save_processed_cache(cache_dir, processed_dataset, train_indices, val_indices, test_indices)

    train_dataset = Subset(processed_dataset, train_indices)
    val_dataset = Subset(processed_dataset, val_indices)
    test_dataset = Subset(processed_dataset, test_indices)

    num_workers = 0 if config.dataset_input.loading_strategy == "lazy" else 2
    train_loader = DataLoader(train_dataset, batch_size=config.training.batch_size, shuffle=True, num_workers=num_workers)
    val_loader = DataLoader(val_dataset, batch_size=config.training.batch_size, shuffle=False, num_workers=num_workers)
    test_loader = DataLoader(test_dataset, batch_size=config.training.batch_size, shuffle=False, num_workers=num_workers)

    return DataBundle(
        dataset=processed_dataset,
        train_loader=train_loader,
        val_loader=val_loader,
        test_loader=test_loader,
    )
