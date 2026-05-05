# -*- coding: utf-8 -*-
"""
dataset.py — ProcessedGenomicDataset e CachedProcessedDataset.
"""

import json
import math
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn
from torch.utils.data import Dataset

from genotype_based_predictor.config import PipelineConfig
from genotype_based_predictor.dynamic_indel_alignment import DynamicIndelAligner
from genotype_based_predictor.indel_tensor_builder import build_aligned_haplotype_tensor
from genotype_based_predictor.normalization import (
    apply_normalization,
    log_normalize,
    minmax_keep_zero,
    zscore_normalize,
)
from genotype_based_predictor.utils import taint_sample

console = Console()


# ===========================================================================
# ProcessedGenomicDataset
# ===========================================================================

class ProcessedGenomicDataset(Dataset):
    """
    Dataset que processa dados brutos on-the-fly.

    Aplica extração de centro de janela, downsampling, combinação de
    haplótipos e normalização. Os parâmetros de normalização são
    computados apenas sobre treino+validação para evitar data leakage.

    Parameters
    ----------
    base_dataset : GenomicLongevityDataset
        Dataset base.
    config : PipelineConfig
        Configuração validada do experimento.
    normalization_params : Optional[Dict]
        Parâmetros pré-computados. Se None e ``compute_normalization=True``,
        serão calculados varrendo o dataset.
    compute_normalization : bool
        Se True, computa parâmetros varrendo todo o dataset.
    """

    def __init__(
        self,
        base_dataset: Any,
        config: PipelineConfig,
        normalization_params: Optional[Dict] = None,
        compute_normalization: bool = True,
    ):
        self.base_dataset = base_dataset
        self.config = config

        di = config.dataset_input
        self.alphagenome_outputs = di.alphagenome_outputs
        self.ontology_terms = di.ontology_terms
        self.haplotype_mode = di.haplotype_mode
        self.window_center_size = di.window_center_size
        self.downsample_factor = di.downsample_factor
        self.normalization_method = di.normalization_method
        self.normalization_value = di.normalization_value
        self.selected_track_index = di.selected_track_index
        self.genes_to_use = set(di.genes_to_use or [])
        self.indel_neutral_value = di.indel_neutral_value
        self.indel_include_valid_mask = di.indel_include_valid_mask
        out = config.output
        self.prediction_target = out.prediction_target
        self.derived_targets = out.derived_targets
        self._derived_target_config = self.derived_targets.get(self.prediction_target)
        self.dataset_metadata = getattr(self.base_dataset, "dataset_metadata", {}) or {}
        self.individuals = self.dataset_metadata.get("individuals", [])
        self.selected_sample_ids = self._resolve_selected_sample_ids()
        self.dynamic_indel_aligner = DynamicIndelAligner(
            Path(di.dataset_dir),
            selected_sample_ids=self.selected_sample_ids,
        )
        self._base_item_cache: OrderedDict[int, Tuple[Any, Any]] = OrderedDict()
        self._base_item_cache_limit = 4
        self._processed_item_cache: OrderedDict[int, Tuple[torch.Tensor, torch.Tensor]] = OrderedDict()
        self._processed_item_cache_limit = 8
        self.profile_stats = {
            "base_fetch_s": 0.0,
            "window_process_s": 0.0,
            "normalization_s": 0.0,
            "getitem_calls": 0,
        }

        self.valid_sample_indices: List[int] = list(range(len(self.base_dataset)))
        self._valid_sample_index_set = set(self.valid_sample_indices)
        self.target_to_idx: Dict[str, int] = {}
        self.idx_to_target: Dict[int, str] = {}

        if normalization_params is not None:
            self.normalization_params = normalization_params
        elif compute_normalization:
            self.normalization_params = self._compute_normalization_params()
        else:
            self.normalization_params = {"method": self.normalization_method, "mean": 0.0, "std": 1.0}

        self._create_target_mappings()

    def _resolve_selected_sample_ids(self) -> Optional[set[str]]:
        di = self.config.dataset_input
        selected = set(di.sample_ids or [])

        if di.sample_ids_path:
            sample_ids_path = Path(di.sample_ids_path)
            with open(sample_ids_path) as f:
                if sample_ids_path.suffix.lower() == ".json":
                    payload = json.load(f)
                    if isinstance(payload, list):
                        selected.update(str(item) for item in payload)
                    else:
                        selected.update(str(item) for item in payload.get("sample_ids", []))
                else:
                    selected.update(line.strip() for line in f if line.strip())

        if not selected and not di.superpopulations_to_use and not di.populations_to_use:
            return None
        return selected

    # ------------------------------------------------------------------
    # Normalização
    # ------------------------------------------------------------------

    def _compute_normalization_params(self) -> Dict:
        if self.normalization_method in {"log", "minmax_keep_zero"} and self.normalization_value not in (0, 0.0, None):
            num_genes = len(self.config.dataset_input.genes_to_use or []) or 1
            num_ontologies = len(self.ontology_terms) if self.ontology_terms else 1
            num_tracks = (2 * num_ontologies + 6) * num_genes
            if self.normalization_method == "log":
                track_params = [{"log_max": float(self.normalization_value)} for _ in range(num_tracks)]
            else:
                track_params = [{"max": float(self.normalization_value)} for _ in range(num_tracks)]

            console.print(
                f"[green]✓ Normalização rápida: usando normalization_value={self.normalization_value} para {num_tracks} tracks[/green]"
            )
            return {
                "method": self.normalization_method,
                "per_track": True,
                "num_tracks": len(track_params),
                "track_params": track_params,
            }

        track_max: Optional[np.ndarray] = None
        track_count: Optional[np.ndarray] = None
        track_mean: Optional[np.ndarray] = None
        track_M2: Optional[np.ndarray] = None
        num_processed = 0

        console.print(f"[cyan]Computando normalização ({self.normalization_method})...[/cyan]")

        with Progress(
            SpinnerColumn(), TextColumn("{task.description}"),
            BarColumn(), TextColumn("{task.percentage:>3.0f}%"),
            console=console,
        ) as progress:
            task = progress.add_task("Coletando valores...", total=len(self.base_dataset))

            for idx in range(len(self.base_dataset)):
                try:
                    input_data, _ = self.base_dataset[idx]
                    processed = self._process_windows(input_data["windows"])
                    if processed.size > 0:
                        num_tracks = processed.shape[0]
                        if self.normalization_method == "zscore":
                            if track_count is None:
                                track_count = np.zeros(num_tracks, dtype=np.float64)
                                track_mean = np.zeros(num_tracks, dtype=np.float64)
                                track_M2 = np.zeros(num_tracks, dtype=np.float64)
                            for ti in range(num_tracks):
                                row = processed[ti].astype(np.float64)
                                n = row.size
                                if n == 0:
                                    continue
                                batch_mean = row.mean()
                                batch_var = row.var()
                                old_count = track_count[ti]
                                new_count = old_count + n
                                delta = batch_mean - track_mean[ti]
                                track_mean[ti] += delta * n / new_count
                                track_M2[ti] += batch_var * n + delta * delta * old_count * n / new_count
                                track_count[ti] = new_count
                        else:
                            if track_max is None:
                                track_max = np.zeros(num_tracks, dtype=np.float64)
                            for ti in range(num_tracks):
                                row = processed[ti]
                                nz = row[row > 0]
                                if len(nz) > 0:
                                    track_max[ti] = max(track_max[ti], float(nz.max()))
                        num_processed += 1
                except Exception:
                    pass
                progress.update(task, advance=1)

        if track_max is None and track_count is None:
            return {"method": self.normalization_method, "per_track": False, "mean": 0.0, "std": 1.0}

        track_params: List[Dict] = []
        num_tracks = len(track_count) if track_count is not None else len(track_max)
        for ti in range(num_tracks):
            if self.normalization_method == "zscore":
                mean = float(track_mean[ti])
                variance = float(track_M2[ti] / track_count[ti]) if track_count[ti] > 0 else 0.0
                std = float(max(np.sqrt(variance), 1e-8))
                track_params.append({"mean": mean, "std": std})
            elif self.normalization_method == "minmax_keep_zero":
                xmax = float(track_max[ti]) if track_max[ti] > 0 else 1.0
                track_params.append({"max": xmax})
            elif self.normalization_method == "log":
                xmax = float(track_max[ti]) if track_max[ti] > 0 else 1.0
                track_params.append({"log_max": float(np.log1p(xmax))})
            else:
                mean = float(track_mean[ti]) if track_mean is not None else 0.0
                variance = float(track_M2[ti] / track_count[ti]) if track_count is not None and track_count[ti] > 0 else 0.0
                std = float(max(np.sqrt(variance), 1e-8))
                track_params.append({"mean": mean, "std": std})

        console.print(f"[green]✓ Normalização: {num_processed} amostras, {len(track_params)} tracks[/green]")
        return {
            "method": self.normalization_method,
            "per_track": True,
            "num_tracks": len(track_params),
            "track_params": track_params,
        }

    # ------------------------------------------------------------------
    # Target mappings
    # ------------------------------------------------------------------

    def _create_target_mappings(self):
        known = self.config.output.known_classes
        if known:
            self.target_to_idx = {t: i for i, t in enumerate(known)}
            self.idx_to_target = {i: t for t, i in self.target_to_idx.items()}
            self._discover_valid_sample_indices()
            return

        # Tentar carregar do metadata.json do dataset
        dataset_metadata = getattr(self.base_dataset, "dataset_metadata", None)
        if dataset_metadata is None:
            meta_file = Path(self.config.dataset_input.dataset_dir) / "dataset_metadata.json"
            if meta_file.exists():
                with open(meta_file) as f:
                    dataset_metadata = json.load(f)

        classes_from_metadata: Optional[List[str]] = None
        if dataset_metadata:
            pt = self.prediction_target
            if pt == "superpopulation":
                dist = dataset_metadata.get("superpopulation_distribution", {})
                classes_from_metadata = sorted(dist.keys()) if dist else None
            elif pt == "population":
                dist = dataset_metadata.get("population_distribution", {})
                classes_from_metadata = sorted(dist.keys()) if dist else None
            elif pt in self.derived_targets:
                pedigree = dataset_metadata.get("individuals_pedigree", {})
                classes_from_metadata = sorted({
                    t for t in (self._get_target_value(p) for p in pedigree.values())
                    if t is not None
                })

        if classes_from_metadata:
            self.target_to_idx = {t: i for i, t in enumerate(classes_from_metadata)}
            self.idx_to_target = {i: t for t, i in self.target_to_idx.items()}
            self._discover_valid_sample_indices()
            console.print(f"[green]✓ Classes: {classes_from_metadata}[/green]")
            return

        # Fallback: varrer dataset
        console.print("[yellow]Varrendo dataset para descobrir classes...[/yellow]")
        unique_targets: set = set()
        for idx in range(len(self.base_dataset)):
            try:
                _, out = self.base_dataset[idx]
                t = self._get_target_value(out)
                if t is not None:
                    unique_targets.add(t)
            except Exception:
                pass
        sorted_t = sorted(unique_targets)
        self.target_to_idx = {t: i for i, t in enumerate(sorted_t)}
        self.idx_to_target = {i: t for t, i in self.target_to_idx.items()}
        self._discover_valid_sample_indices()
        console.print(f"[green]✓ Classes encontradas: {sorted_t}[/green]")

    def _discover_valid_sample_indices(self):
        if self.prediction_target == "frog_likelihood":
            self.valid_sample_indices = list(range(len(self.base_dataset)))
            self._valid_sample_index_set = set(self.valid_sample_indices)
            return

        dataset_metadata = getattr(self.base_dataset, "dataset_metadata", None)
        if dataset_metadata is None:
            meta_file = Path(self.config.dataset_input.dataset_dir) / "dataset_metadata.json"
            if meta_file.exists():
                with open(meta_file) as f:
                    dataset_metadata = json.load(f)

        if dataset_metadata:
            individuals = dataset_metadata.get("individuals", [])
            pedigree = dataset_metadata.get("individuals_pedigree", {})
            if individuals and pedigree:
                valid = []
                with Progress(
                    SpinnerColumn(), TextColumn("{task.description}"),
                    BarColumn(), TextColumn("{task.percentage:>3.0f}%"),
                    TimeElapsedColumn(), console=console,
                ) as prog:
                    task = prog.add_task(
                        f"Filtrando para '{self.prediction_target}'...",
                        total=len(individuals),
                    )
                    for idx, sid in enumerate(individuals):
                        if not self._sample_matches_logical_view(sid, pedigree.get(sid, {})):
                            prog.update(task, advance=1)
                            continue
                        t = self._get_target_value(pedigree.get(sid, {}))
                        if t in self.target_to_idx:
                            valid.append(idx)
                        prog.update(task, advance=1)
                self.valid_sample_indices = valid
                self._valid_sample_index_set = set(valid)
                return

        # Fallback
        valid = []
        for idx in range(len(self.base_dataset)):
            try:
                _, out = self.base_dataset[idx]
                sample_id = self.individuals[idx] if idx < len(self.individuals) else None
                if sample_id and not self._sample_matches_logical_view(sample_id, out):
                    continue
                if self._get_target_value(out) in self.target_to_idx:
                    valid.append(idx)
            except Exception:
                pass
        self.valid_sample_indices = valid
        self._valid_sample_index_set = set(valid)

    def _sample_matches_logical_view(self, sample_id: str, metadata: Dict) -> bool:
        di = self.config.dataset_input
        if self.selected_sample_ids is not None and sample_id not in self.selected_sample_ids:
            return False
        if di.superpopulations_to_use and metadata.get("superpopulation") not in set(di.superpopulations_to_use):
            return False
        if di.populations_to_use and metadata.get("population") not in set(di.populations_to_use):
            return False
        return True

    def _get_target_value(self, output_data: Dict) -> Optional[str]:
        pt = self.prediction_target
        if pt == "superpopulation":
            return output_data.get("superpopulation")
        elif pt == "population":
            return output_data.get("population")
        elif pt in self.derived_targets:
            return self._get_derived_target_value(output_data)
        elif pt == "frog_likelihood":
            return None
        raise ValueError(f"prediction_target inválido: {pt}")

    def _get_derived_target_value(self, output_data: Dict) -> Optional[str]:
        dc = self._derived_target_config
        if dc is None:
            return None
        source_value = output_data.get(dc.source_field)
        if source_value is None:
            return None
        for cls, values in dc.class_map.items():
            if source_value in values:
                return cls
        return None if dc.exclude_unmapped else source_value

    # ------------------------------------------------------------------
    # Processamento de janelas
    # ------------------------------------------------------------------

    def _process_window_haplotype_channels(
        self,
        sample_id: str,
        window_name: str,
        haplotype: str,
        predictions: Dict,
        prediction_metadata: Optional[Dict] = None,
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        if len(self.alphagenome_outputs) != 1:
            raise ValueError("tensor_layout='haplotype_channels' atualmente requer exatamente um alphagenome_output")

        output_type = self.alphagenome_outputs[0]
        if output_type not in predictions:
            return None

        array = predictions[output_type]
        track_meta = prediction_metadata.get(output_type) if prediction_metadata else None
        if self.downsample_factor != 1:
            raise ValueError("tensor_layout='haplotype_channels' nao suporta downsample_factor diferente de 1")
        if self.window_center_size <= 0:
            raise ValueError("tensor_layout='haplotype_channels' requer window_center_size > 0")
        entry = self.dynamic_indel_aligner.get_haplotype_entry(window_name, sample_id, haplotype)
        if entry is None:
            return None

        expanded_length = self.dynamic_indel_aligner.get_expanded_length(window_name)
        signal_rows = []
        shared_masks = None
        if array.ndim == 2:
            if track_meta:
                track_indices = self._filter_track_indices(output_type, track_meta)
            else:
                track_indices = [self.selected_track_index]
            for track_index in track_indices:
                if not (0 <= track_index < array.shape[1]):
                    raise ValueError(f"selected_track_index={track_index} fora do limite para {output_type} com shape={array.shape}")
                row = np.asarray(array[:, track_index], dtype=np.float32)
                aligned = build_aligned_haplotype_tensor(
                    row=row,
                    entry=entry,
                    expanded_length=expanded_length,
                    neutral_value=self.indel_neutral_value,
                    include_valid_mask=True,
                )
                signal_rows.append(aligned[0:1, :])
                if shared_masks is None:
                    shared_masks = aligned[1:, :]
        else:
            row = np.asarray(array.flatten() if array.ndim > 1 else array, dtype=np.float32)
            aligned = build_aligned_haplotype_tensor(
                row=row,
                entry=entry,
                expanded_length=expanded_length,
                neutral_value=self.indel_neutral_value,
                include_valid_mask=True,
            )
            signal_rows.append(aligned[0:1, :])
            shared_masks = aligned[1:, :]

        if not signal_rows or shared_masks is None:
            return None

        signals = np.concatenate(signal_rows, axis=0)
        center = signals.shape[1] // 2
        half = self.window_center_size // 2
        start = max(0, center - half)
        end = min(signals.shape[1], start + self.window_center_size)
        return signals[:, start:end], shared_masks[:, start:end]

    def _filter_track_indices(self, output_type: str, track_metadata: List[Dict]) -> List[int]:
        if not self.ontology_terms:
            return list(range(len(track_metadata)))
        requested = set(self.ontology_terms)
        indices = [i for i, m in enumerate(track_metadata) if m.get("ontology_curie") in requested]
        if not indices:
            available = sorted({m.get("ontology_curie", "") for m in track_metadata})
            raise ValueError(f"Nenhuma track para {sorted(requested)} em {output_type}. Disponíveis: {available}")
        return indices

    def _process_windows(self, windows: Dict, sample_id: Optional[str] = None) -> np.ndarray:
        if self.config.dataset_input.tensor_layout != "haplotype_channels":
            raise ValueError("O pipeline atual requer tensor_layout='haplotype_channels'")
        return self._process_windows_haplotype_channels(windows, sample_id=sample_id)

    def _process_windows_haplotype_channels(self, windows: Dict, sample_id: Optional[str] = None) -> np.ndarray:
        if self.haplotype_mode != "H1+H2":
            raise ValueError("tensor_layout='haplotype_channels' requer haplotype_mode='H1+H2'")
        gene_names = list(self.config.dataset_input.genes_to_use or windows.keys())
        if not gene_names:
            return np.array([]).reshape(2, 4, 0)

        h1_rows = []
        h2_rows = []

        for gene_name in gene_names:
            window_data = windows.get(gene_name)
            if window_data is None:
                continue

            window_sample_id = sample_id or self._infer_sample_id_from_window(window_data)
            if window_sample_id is None:
                raise ValueError("Nao foi possivel inferir sample_id para tensor_layout='haplotype_channels'")

            h1 = self._process_window_haplotype_channels(
                window_sample_id,
                gene_name,
                "H1",
                window_data.get("predictions_h1", {}),
                window_data.get("prediction_metadata_h1"),
            )
            h2 = self._process_window_haplotype_channels(
                window_sample_id,
                gene_name,
                "H2",
                window_data.get("predictions_h2", {}),
                window_data.get("prediction_metadata_h2"),
            )
            if h1 is None or h2 is None:
                continue

            h1_signals, h1_masks = h1
            h2_signals, h2_masks = h2
            h1_rows.append(np.concatenate([h1_signals, h1_masks], axis=0))
            h2_rows.append(np.concatenate([h2_signals, h2_masks], axis=0))

        if not h1_rows or not h2_rows:
            return np.array([]).reshape(2, 4, 0)

        return np.stack([
            np.concatenate(h1_rows, axis=0),
            np.concatenate(h2_rows, axis=0),
        ], axis=0)

    def _infer_sample_id_from_window(self, window_data: Dict) -> Optional[str]:
        window_meta = window_data.get("window_metadata") or {}
        sample_id = window_meta.get("sample_id") or window_meta.get("individual_id")
        if sample_id:
            return sample_id
        return None

    # ------------------------------------------------------------------
    # Dataset interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self.valid_sample_indices)

    def _load_base_item(self, base_idx: int) -> Tuple[Any, Any]:
        base_cached = self._base_item_cache.get(base_idx)
        if base_cached is not None:
            self._base_item_cache.move_to_end(base_idx)
            return base_cached

        t0 = time.perf_counter()
        input_data, output_data = self.base_dataset[base_idx]
        self.profile_stats["base_fetch_s"] += time.perf_counter() - t0
        self._base_item_cache[base_idx] = (input_data, output_data)
        while len(self._base_item_cache) > self._base_item_cache_limit:
            self._base_item_cache.popitem(last=False)
        return input_data, output_data

    def _normalize_features_tensor(self, features_tensor: torch.Tensor) -> torch.Tensor:
        method = self.normalization_params.get("method", "zscore")
        per_track = self.normalization_params.get("per_track", False)

        if features_tensor.ndim != 3:
            raise ValueError(f"Esperado tensor 3D (2,4,L), recebido shape={tuple(features_tensor.shape)}")
        if per_track:
            t0 = time.perf_counter()
            flat = features_tensor.view(-1, features_tensor.shape[-1])
            track_params = self.normalization_params["track_params"]
            rows = []
            for ti in range(flat.shape[0]):
                p = track_params[ti]
                row = flat[ti: ti + 1, :]
                if method == "zscore":
                    row = zscore_normalize(row, p["mean"], p["std"])
                elif method == "minmax_keep_zero":
                    row = minmax_keep_zero(row, p["max"])
                elif method == "log":
                    row = log_normalize(row, p["log_max"])
                rows.append(row)
            features_tensor = torch.cat(rows, dim=0).view_as(features_tensor)
            self.profile_stats["normalization_s"] += time.perf_counter() - t0
            return features_tensor

        t0 = time.perf_counter()
        features_tensor = apply_normalization(features_tensor, self.normalization_params)
        self.profile_stats["normalization_s"] += time.perf_counter() - t0
        return features_tensor

    def _build_target_tensor(self, output_data: Dict, features_tensor: torch.Tensor) -> torch.Tensor:
        if self.prediction_target == "frog_likelihood":
            return torch.FloatTensor(output_data.get("frog_likelihood", np.zeros(150)))

        target_value = self._get_target_value(output_data)
        if target_value in self.target_to_idx:
            target_idx = self.target_to_idx[target_value]
            target_tensor = torch.tensor(target_idx, dtype=torch.long)
            if self.config.debug.taint_at_runtime:
                d = self.config.debug
                features_tensor = taint_sample(
                    features_tensor, target_idx, len(self.target_to_idx),
                    taint_type=d.taint_type,
                    taint_value=d.taint_value,
                    taint_horizontal_size=d.taint_horizontal_size,
                    taint_vertical_size=d.taint_vertical_size,
                    taint_horizontal_step=d.taint_horizontal_step,
                    taint_vertical_step=d.taint_vertical_step,
                )
            return target_tensor

        return torch.tensor(-1, dtype=torch.long)

    def _process_single_index(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        base_idx = self.valid_sample_indices[idx]
        input_data, output_data = self._load_base_item(base_idx)

        sample_id = self._infer_sample_id_from_input(input_data)
        if sample_id is None and base_idx < len(self.individuals):
            sample_id = self.individuals[base_idx]
            self._inject_sample_id_into_windows(input_data, sample_id)

        t0 = time.perf_counter()
        features = self._process_windows(input_data["windows"], sample_id=sample_id)
        self.profile_stats["window_process_s"] += time.perf_counter() - t0
        features_tensor = self._normalize_features_tensor(torch.FloatTensor(features))
        target_tensor = self._build_target_tensor(output_data, features_tensor)
        self.profile_stats["getitem_calls"] += 1
        return features_tensor, target_tensor

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        cached = self._processed_item_cache.get(idx)
        if cached is not None:
            self._processed_item_cache.move_to_end(idx)
            return cached

        features_tensor, target_tensor = self._process_single_index(idx)

        self._processed_item_cache[idx] = (features_tensor, target_tensor)
        while len(self._processed_item_cache) > self._processed_item_cache_limit:
            self._processed_item_cache.popitem(last=False)

        return features_tensor, target_tensor

    def get_items_bulk(self, indices: List[int]) -> List[Tuple[torch.Tensor, torch.Tensor]]:
        items: List[Tuple[torch.Tensor, torch.Tensor]] = []
        for idx in indices:
            cached = self._processed_item_cache.get(idx)
            if cached is not None:
                self._processed_item_cache.move_to_end(idx)
                items.append(cached)
                continue

            item = self._process_single_index(idx)
            self._processed_item_cache[idx] = item
            while len(self._processed_item_cache) > self._processed_item_cache_limit:
                self._processed_item_cache.popitem(last=False)
            items.append(item)
        return items

    def _infer_sample_id_from_input(self, input_data: Dict) -> Optional[str]:
        windows = input_data.get("windows", {})
        for window_data in windows.values():
            sample_id = self._infer_sample_id_from_window(window_data)
            if sample_id:
                return sample_id
        return None

    def _inject_sample_id_into_windows(self, input_data: Dict, sample_id: str) -> None:
        for window_data in input_data.get("windows", {}).values():
            window_meta = window_data.setdefault("window_metadata", {})
            window_meta.setdefault("sample_id", sample_id)

    def get_num_classes(self) -> int:
        return len(self.target_to_idx)

    def get_class_names(self) -> List[str]:
        return [self.idx_to_target[i] for i in range(len(self.idx_to_target))]

    def get_input_shape(self) -> Tuple[int, int]:
        if self.config.dataset_input.tensor_layout == "haplotype_channels":
            effective = self.window_center_size // max(self.downsample_factor, 1)
            num_genes = len(self.config.dataset_input.genes_to_use or []) or 1
            num_ontologies = len(self.ontology_terms) if self.ontology_terms else 1
            return ((2 * num_ontologies + 6) * num_genes, effective)

        for idx in range(min(64, len(self))):
            try:
                f, _ = self[idx]
                s = f.shape
                if len(s) != 3:
                    raise ValueError(f"Esperado shape 3D canônico (2,4,L), recebido {tuple(s)}")
                return (s[0] * s[1], s[2])
            except Exception:
                continue
        effective = self.window_center_size // max(self.downsample_factor, 1)
        num_genes = len(self.config.dataset_input.genes_to_use or []) or 1
        num_ontologies = len(self.ontology_terms) if self.ontology_terms else 1
        return ((2 * num_ontologies + 6) * num_genes, effective)

    def get_input_size(self) -> int:
        r, c = self.get_input_shape()
        return r * c


# ===========================================================================
# CachedProcessedDataset
# ===========================================================================

class CachedProcessedDataset(Dataset):
    """
    Dataset que carrega dados pré-processados do cache em disco (``.pt``).

    Modos de carregamento
    ---------------------
    ``preload``
        Carrega tudo em RAM na inicialização (rápido, usa mais memória).
    ``lazy``
        Carrega sob demanda com LRU cache em memória.

    Parameters
    ----------
    data_file : Path
        Arquivo ``.pt`` com lista de tuplas ``(features_tensor, target_tensor)``.
    target_to_idx : Dict[str, int]
    idx_to_target : Dict[int, str]
    config : PipelineConfig
    split_name : Optional[str]
    length_hint : Optional[int]
    """

    def __init__(
        self,
        data_file: Path,
        target_to_idx: Dict[str, int],
        idx_to_target: Dict[int, str],
        config: PipelineConfig,
        split_name: Optional[str] = None,
        length_hint: Optional[int] = None,
    ):
        self.data_file = data_file
        self.target_to_idx = target_to_idx
        self.idx_to_target = idx_to_target
        self.config = config
        self.split_name = (split_name or self._infer_split_name()).lower()
        self._length_hint = length_hint
        self._sample_metadata = self._load_sample_metadata()
        self._shard_index = self._load_shard_index()

        self.loading_strategy = config.data_loading.loading_strategy
        self.cache_size = config.data_loading.cache_size

        console.print(f"[cyan]Preparando {data_file.name}...[/cyan]")
        if self.loading_strategy == "preload" and self._shard_index is None:
            self.data = torch.load(data_file)
            self._length = len(self.data)
            self._data_loaded = True
            console.print(f"[green]✓ {self._length} samples (preload)[/green]")
        elif self.loading_strategy == "preload" and self._shard_index is not None:
            self._length = self._determine_length()
            self._data_loaded = False
            self.data = None
            self._cache = {}
            self._cache_order = []
            self._loaded_shard_path = None
            console.print(f"[green]✓ {self._length} samples (sharded preload-lazy)[/green]")
        elif self.loading_strategy == "lazy":
            self._length = self._determine_length()
            self._data_loaded = False
            self.data = None
            self._cache: Dict[int, Tuple] = {}
            self._cache_order: List[int] = []
            self._loaded_shard_path = None
            console.print(f"[green]✓ {self._length} samples (lazy)[/green]")
        else:
            raise ValueError(f"loading_strategy inválida: {self.loading_strategy}")

    def _infer_split_name(self) -> str:
        stem = self.data_file.stem.lower()
        for key, val in {"train": "train", "val": "val", "validation": "val", "test": "test"}.items():
            if key in stem:
                return val
        return "train"

    def _determine_length(self) -> int:
        if self._length_hint is not None:
            return self._length_hint
        meta_path = self.data_file.parent / "metadata.json"
        if meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
            sz = meta.get("splits", {}).get(f"{self.split_name}_size")
            if sz:
                return sz
        raise RuntimeError(f"Não foi possível determinar tamanho de {self.data_file.name}")

    def _load_sample_metadata(self) -> Optional[List[Dict]]:
        try:
            split_index_file = self.data_file.parent / "split_index.json"
            if split_index_file.exists():
                with open(split_index_file) as f:
                    split_index = json.load(f)
                sample_ids = split_index.get(self.split_name, []) or []
                if sample_ids:
                    return [{"sample_id": sid} for sid in sample_ids]

            splits_file = self.data_file.parent / "splits_metadata.json"
            if not splits_file.exists():
                return None
            with open(splits_file) as f:
                splits = json.load(f)
            if splits.get("format_version", 1) >= 2:
                return splits.get(self.split_name, []) or None
        except Exception:
            return None
        return None

    def _load_shard_index(self) -> Optional[Dict]:
        shard_index_file = self.data_file.parent / "shards_index.json"
        if not shard_index_file.exists():
            return None
        with open(shard_index_file) as f:
            index = json.load(f)
        return index.get(self.split_name)

    def get_sample_id(self, idx: int) -> str:
        if self._sample_metadata and idx < len(self._sample_metadata):
            return self._sample_metadata[idx].get("sample_id", f"#{idx + 1}")
        return f"#{idx + 1}"

    def get_sample_metadata(self, idx: int) -> Dict:
        if self._sample_metadata and idx < len(self._sample_metadata):
            return self._sample_metadata[idx]
        return {"sample_id": f"#{idx + 1}", "superpopulation": "UNK", "population": "UNK", "sex": 0}

    def __len__(self) -> int:
        return self._length

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, int]:
        if self.loading_strategy == "preload":
            if self._shard_index is None:
                features, target = self.data[idx]
            else:
                features, target = self._get_item_from_shards(idx)
        else:
            if idx in self._cache:
                features, target = self._cache[idx]
            else:
                if self._shard_index is None:
                    if not self._data_loaded or self.data is None:
                        self.data = torch.load(self.data_file)
                        self._data_loaded = True
                    features, target = self.data[idx]
                else:
                    features, target = self._get_item_from_shards(idx)
                if self.cache_size > 0:
                    if len(self._cache) >= self.cache_size:
                        oldest = self._cache_order.pop(0)
                        del self._cache[oldest]
                    self._cache[idx] = (features, target)
                    self._cache_order.append(idx)

        if (
            self.config.debug.taint_at_runtime
            and self.config.output.prediction_target != "frog_likelihood"
        ):
            d = self.config.debug
            target_class = target.item() if target.ndim == 0 else target[0].item()
            features = taint_sample(
                features, target_class, len(self.target_to_idx),
                taint_type=d.taint_type,
                taint_value=d.taint_value,
                taint_horizontal_size=d.taint_horizontal_size,
                taint_vertical_size=d.taint_vertical_size,
                taint_horizontal_step=d.taint_horizontal_step,
                taint_vertical_step=d.taint_vertical_step,
            )

        return features, target, idx

    def unload_data(self):
        if self.loading_strategy == "lazy" and self.data is not None:
            import gc
            del self.data
            self.data = None
            self._data_loaded = False
            self._loaded_shard_path = None
            gc.collect()

    def get_num_classes(self) -> int:
        return len(self.target_to_idx)

    def get_class_names(self) -> List[str]:
        return [self.idx_to_target[i] for i in range(len(self.idx_to_target))]

    def get_input_shape(self) -> Tuple[int, int]:
        if self._shard_index is not None and self._length > 0:
            features, _target, _idx = self[0]
            s = features.shape
            if len(s) == 3:
                return (s[0] * s[1], s[2])
            if len(s) == 2:
                return (s[0], s[1])
            return (1, s[0])

        if self.loading_strategy == "lazy" and not self._data_loaded:
            self.data = torch.load(self.data_file)
            self._data_loaded = True
        if self.data and len(self.data) > 0:
            s = self.data[0][0].shape
            if len(s) == 3:
                return (s[0] * s[1], s[2])
            if len(s) == 2:
                return (s[0], s[1])
            return (1, s[0])
        return (0, 0)

    def _get_item_from_shards(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        if self._shard_index is None:
            raise RuntimeError("Shard index ausente")

        shard_size = self._shard_index["shard_size"]
        shard_id = idx // shard_size
        item_offset = idx % shard_size
        shard_name = self._shard_index["shards"][shard_id]
        shard_path = self.data_file.parent / shard_name

        if self._loaded_shard_path != shard_path or self.data is None:
            self.data = torch.load(shard_path)
            self._data_loaded = True
            self._loaded_shard_path = shard_path

        return self.data[item_offset]

    def get_input_size(self) -> int:
        r, c = self.get_input_shape()
        return r * c
