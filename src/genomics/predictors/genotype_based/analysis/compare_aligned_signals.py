from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch

from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config


@dataclass(frozen=True)
class SampleItem:
    split: str
    split_index: int
    sample_id: str
    superpopulation: str
    population: str
    features: torch.Tensor


@dataclass(frozen=True)
class ChannelLayout:
    genes: List[str]
    signal_channels_per_gene: int
    mask_channels_per_gene: int
    channels_per_gene: int
    valid_mask_offset: Optional[int]
    signal_names: List[str]


def _read_json(path: Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def _open_csv(path: Path, fieldnames: Sequence[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    f = open(path, "w", newline="", encoding="utf-8")
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    return f, writer


def _trim_top_differences(rows: List[Dict], top_k: int) -> None:
    if top_k <= 0:
        rows.clear()
        return
    limit = max(top_k * 10, top_k)
    if len(rows) <= limit:
        return
    rows.sort(key=lambda row: float(row["abs_diff"]), reverse=True)
    del rows[top_k:]


def _safe_float(value: float) -> Optional[float]:
    if value is None or not math.isfinite(float(value)):
        return None
    return float(value)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2:
        return float("nan")
    x_std = float(np.std(x))
    y_std = float(np.std(y))
    if x_std == 0.0 or y_std == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _cosine(x: np.ndarray, y: np.ndarray) -> float:
    denom = float(np.linalg.norm(x) * np.linalg.norm(y))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(x, y) / denom)


def _metadata_for_samples(dataset_dir: Path) -> Dict[str, Dict]:
    meta_path = dataset_dir / "dataset_metadata.json"
    if meta_path.exists():
        meta = _read_json(meta_path)
        pedigree = meta.get("individuals_pedigree", {}) or {}
        return {str(k): dict(v or {}) for k, v in pedigree.items()}
    return {}


def _class_name(config, target_idx: int) -> str:
    known = config.output.known_classes
    if known and 0 <= target_idx < len(known):
        return str(known[target_idx])
    return str(target_idx)


def _split_sample_ids(cache_dir: Path, split: str) -> List[str]:
    split_index = _read_json(cache_dir / "split_index.json")
    return [str(sample_id) for sample_id in split_index.get(split, [])]


def _validate_processed_cache(cache_dir: Path) -> None:
    missing = []
    if not cache_dir.exists():
        missing.append(str(cache_dir))
    for file_name in ("split_index.json", "shards_index.json"):
        path = cache_dir / file_name
        if not path.exists():
            missing.append(str(path))
    if not missing:
        return
    details = "\n".join(f"  - {path}" for path in missing)
    raise SystemExit(
        "Cache processado ausente ou incompleto para compare-aligned-signals.\n"
        f"Arquivos esperados nao encontrados:\n{details}\n\n"
        "Materialize o cache antes de comparar sinais, por exemplo:\n"
        "  genomics genotype prepare-cache <mesmo-config-yaml>\n\n"
        "Se voce ja tem um cache em outro local, passe explicitamente:\n"
        "  genomics genotype compare-aligned-signals ... --cache-dir /caminho/do/cache"
    )


def _iter_split_items(cache_dir: Path, split: str) -> Iterable[Tuple[int, Tuple[torch.Tensor, torch.Tensor]]]:
    shard_index_path = cache_dir / "shards_index.json"
    if shard_index_path.exists():
        split_index = _read_json(shard_index_path).get(split)
        if not split_index:
            return
        shard_size = int(split_index.get("shard_size", 0) or 0)
        for shard_id, shard_name in enumerate(split_index.get("shards", [])):
            shard_items = torch.load(cache_dir / shard_name, map_location="cpu")
            for offset, item in enumerate(shard_items):
                idx = shard_id * shard_size + offset if shard_size > 0 else offset
                yield idx, item
        return

    data_file = cache_dir / f"{split}_data.pt"
    if not data_file.exists():
        return
    items = torch.load(data_file, map_location="cpu")
    for idx, item in enumerate(items):
        yield idx, item


def _load_samples(
    *,
    cache_dir: Path,
    dataset_dir: Path,
    splits: Sequence[str],
    max_samples: Optional[int],
    sample_ids_filter: Optional[set],
    config,
) -> List[SampleItem]:
    pedigree = _metadata_for_samples(dataset_dir)
    loaded: List[SampleItem] = []
    for split in splits:
        split_sample_ids = _split_sample_ids(cache_dir, split)
        for split_idx, (item_idx, item) in enumerate(_iter_split_items(cache_dir, split)):
            if item_idx >= len(split_sample_ids):
                continue
            sample_id = split_sample_ids[item_idx]
            if sample_ids_filter and sample_id not in sample_ids_filter:
                continue
            features, target = item[:2]
            if not isinstance(features, torch.Tensor):
                features = torch.as_tensor(features)
            features = features.detach().cpu().float()
            ped = pedigree.get(sample_id, {})
            superpopulation = str(ped.get("superpopulation") or "UNK")
            if superpopulation == "UNK" and torch.is_tensor(target) and target.ndim == 0:
                superpopulation = _class_name(config, int(target.item()))
            loaded.append(
                SampleItem(
                    split=split,
                    split_index=item_idx,
                    sample_id=sample_id,
                    superpopulation=superpopulation,
                    population=str(ped.get("population") or "UNK"),
                    features=features,
                )
            )
            if max_samples is not None and len(loaded) >= max_samples:
                return loaded
    return loaded


def _build_layout(config, cache_dir: Path) -> ChannelLayout:
    di = config.dataset_input
    if di.tensor_layout != "haplotype_channels":
        raise ValueError("compare-aligned-signals requer tensor_layout='haplotype_channels'")
    if di.feature_mode == "masks_only":
        raise ValueError("compare-aligned-signals requer canais de sinal; feature_mode='masks_only' nao e suportado")

    genes = list(di.genes_to_use or [])
    view_path = cache_dir / "view_definition.json"
    if not genes and view_path.exists():
        view = _read_json(view_path)
        resolved = view.get("resolved_view", {}) if isinstance(view, dict) else {}
        genes = [str(g) for g in resolved.get("resolved_genes", [])]
    if not genes:
        raise ValueError("Nao foi possivel inferir genes_to_use para decodificar canais")

    signal_channels = len(di.ontology_terms or []) or 1
    mask_channels = 0 if di.feature_mode == "signals_only" else 2 + int(di.indel_include_valid_mask) + int(di.indel_include_snp_mask)
    valid_mask_offset = None
    if di.feature_mode != "signals_only" and di.indel_include_valid_mask:
        valid_mask_offset = signal_channels + 2
    if valid_mask_offset is None:
        raise ValueError(
            "compare-aligned-signals requer indel_include_valid_mask=true e feature_mode com mascaras "
            "para usar apenas posicoes validas nos dois individuos."
        )
    channels_per_gene = signal_channels + mask_channels
    output_name = di.alphagenome_outputs[0] if di.alphagenome_outputs else "signal"
    if di.ontology_terms:
        signal_names = [f"{output_name}:{term}" for term in di.ontology_terms]
    else:
        signal_names = [output_name]
    return ChannelLayout(
        genes=genes,
        signal_channels_per_gene=signal_channels,
        mask_channels_per_gene=mask_channels,
        channels_per_gene=channels_per_gene,
        valid_mask_offset=valid_mask_offset,
        signal_names=signal_names,
    )


def _valid_mask(features: np.ndarray, hap_idx: int, gene_idx: int, layout: ChannelLayout) -> np.ndarray:
    start = gene_idx * layout.channels_per_gene
    if layout.valid_mask_offset is None:
        return np.ones(features.shape[-1], dtype=bool)
    return features[hap_idx, start + layout.valid_mask_offset, :] > 0.5


def _signal_row(features: np.ndarray, hap_idx: int, gene_idx: int, signal_idx: int, layout: ChannelLayout) -> np.ndarray:
    start = gene_idx * layout.channels_per_gene
    return features[hap_idx, start + signal_idx, :]


def _pair_rows(
    a: SampleItem,
    b: SampleItem,
    layout: ChannelLayout,
    min_valid_positions: int,
    top_differences: List[Dict],
    top_k: int,
) -> Iterable[Dict]:
    fa = a.features.numpy()
    fb = b.features.numpy()
    if fa.shape != fb.shape or fa.ndim != 3:
        return
    if fa.shape[0] < 2:
        return

    for hap_idx, hap_name in enumerate(("H1", "H2")):
        if hap_idx >= fa.shape[0]:
            continue
        for gene_idx, gene in enumerate(layout.genes):
            valid = _valid_mask(fa, hap_idx, gene_idx, layout) & _valid_mask(fb, hap_idx, gene_idx, layout)
            n_valid = int(valid.sum())
            if n_valid < min_valid_positions:
                continue
            valid_positions = np.where(valid)[0]
            for signal_idx, signal_name in enumerate(layout.signal_names):
                x = _signal_row(fa, hap_idx, gene_idx, signal_idx, layout)[valid]
                y = _signal_row(fb, hap_idx, gene_idx, signal_idx, layout)[valid]
                diff = np.abs(x - y)
                if diff.size:
                    best_local = int(np.argmax(diff))
                    best_position = int(valid_positions[best_local])
                    best_diff = float(diff[best_local])
                    if top_k > 0:
                        top_differences.append({
                            "split_a": a.split,
                            "split_b": b.split,
                            "sample_a": a.sample_id,
                            "sample_b": b.sample_id,
                            "superpopulation_a": a.superpopulation,
                            "superpopulation_b": b.superpopulation,
                            "same_superpopulation": a.superpopulation == b.superpopulation,
                            "haplotype": hap_name,
                            "gene": gene,
                            "signal": signal_name,
                            "position": best_position,
                            "abs_diff": best_diff,
                            "value_a": float(x[best_local]),
                            "value_b": float(y[best_local]),
                            "n_valid_positions": n_valid,
                        })
                        _trim_top_differences(top_differences, top_k)
                yield {
                    "split_a": a.split,
                    "split_b": b.split,
                    "sample_a": a.sample_id,
                    "sample_b": b.sample_id,
                    "superpopulation_a": a.superpopulation,
                    "superpopulation_b": b.superpopulation,
                    "same_superpopulation": a.superpopulation == b.superpopulation,
                    "haplotype": hap_name,
                    "gene": gene,
                    "signal": signal_name,
                    "n_valid_positions": n_valid,
                    "pearson": _safe_float(_pearson(x, y)),
                    "cosine": _safe_float(_cosine(x, y)),
                    "mad": _safe_float(float(np.mean(diff))) if diff.size else None,
                    "rmse": _safe_float(float(np.sqrt(np.mean((x - y) ** 2)))) if diff.size else None,
                    "max_abs_diff": _safe_float(best_diff if diff.size else float("nan")),
                    "fraction_equal_tol": _safe_float(float(np.mean(diff <= 1e-6))) if diff.size else None,
                }


def _mean(values: List[float]) -> Optional[float]:
    cleaned = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    if not cleaned:
        return None
    return float(sum(cleaned) / len(cleaned))


def _summarize_pairwise(pairwise_path: Path) -> Dict:
    by_same = {"true": {"pearson": [], "mad": [], "rmse": [], "cosine": []}, "false": {"pearson": [], "mad": [], "rmse": [], "cosine": []}}
    by_superpop_pair: Dict[str, Dict[str, List[float]]] = {}
    count = 0
    with open(pairwise_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            count += 1
            same_key = "true" if row.get("same_superpopulation") == "True" else "false"
            pair_key = "|".join(sorted([row.get("superpopulation_a", "UNK"), row.get("superpopulation_b", "UNK")]))
            by_superpop_pair.setdefault(pair_key, {"pearson": [], "mad": [], "rmse": [], "cosine": []})
            for metric in ("pearson", "mad", "rmse", "cosine"):
                value = row.get(metric)
                if value:
                    try:
                        parsed = float(value)
                    except ValueError:
                        continue
                    if math.isfinite(parsed):
                        by_same[same_key][metric].append(parsed)
                        by_superpop_pair[pair_key][metric].append(parsed)
    return {
        "pairwise_rows": count,
        "same_superpopulation": {key: {metric: _mean(vals) for metric, vals in metrics.items()} for key, metrics in by_same.items()},
        "superpopulation_pairs": {key: {metric: _mean(vals) for metric, vals in metrics.items()} for key, metrics in sorted(by_superpop_pair.items())},
    }


def _effect_row(
    *,
    values: np.ndarray,
    labels: np.ndarray,
    superpopulations: Sequence[str],
    min_samples: int,
    min_groups: int,
) -> Optional[Dict]:
    valid = np.isfinite(values)
    n_total = int(valid.sum())
    if n_total < min_samples:
        return None
    x = values[valid].astype(np.float64, copy=False)
    y = labels[valid]
    total_mean = float(np.mean(x))
    total_ss = float(np.sum((x - total_mean) ** 2))
    if total_ss <= 0.0:
        return None

    group_stats = {}
    between_ss = 0.0
    within_ss = 0.0
    eligible_groups = 0
    for group in superpopulations:
        group_values = x[y == group]
        n_group = int(group_values.size)
        if n_group == 0:
            group_stats[group] = {"n": 0, "mean": None}
            continue
        group_mean = float(np.mean(group_values))
        group_stats[group] = {"n": n_group, "mean": group_mean}
        eligible_groups += 1
        between_ss += n_group * (group_mean - total_mean) ** 2
        within_ss += float(np.sum((group_values - group_mean) ** 2))
    if eligible_groups < min_groups:
        return None

    means = [(group, stats["mean"]) for group, stats in group_stats.items() if stats["mean"] is not None]
    if len(means) < min_groups:
        return None
    top_group, top_mean = max(means, key=lambda item: item[1])
    bottom_group, bottom_mean = min(means, key=lambda item: item[1])
    group_delta = float(top_mean - bottom_mean)
    pooled_denom = max(n_total - eligible_groups, 1)
    pooled_std = math.sqrt(max(within_ss / pooled_denom, 0.0))
    standardized_delta = group_delta / pooled_std if pooled_std > 0.0 else float("nan")
    return {
        "n_samples": n_total,
        "n_superpopulations": eligible_groups,
        "total_variance": float(total_ss / n_total),
        "within_superpopulation_variance": float(within_ss / n_total),
        "between_superpopulation_variance": float(between_ss / n_total),
        "eta_squared": float(between_ss / total_ss),
        "max_group_mean_delta": group_delta,
        "standardized_group_delta": _safe_float(standardized_delta),
        "top_superpopulation": top_group,
        "bottom_superpopulation": bottom_group,
        "top_superpopulation_mean": float(top_mean),
        "bottom_superpopulation_mean": float(bottom_mean),
        "group_stats": group_stats,
    }


def _trim_top_rows(rows: List[Dict], key: str, top_k: int) -> None:
    if top_k <= 0:
        rows.clear()
        return
    limit = max(top_k * 10, top_k)
    if len(rows) <= limit:
        return
    rows.sort(key=lambda row: float(row.get(key) or 0.0), reverse=True)
    del rows[top_k:]


def _write_position_effects(
    *,
    samples: Sequence[SampleItem],
    layout: ChannelLayout,
    output_dir: Path,
    top_k: int,
    min_samples: int,
    min_groups: int,
    write_all: bool,
) -> Dict:
    labels = np.asarray([sample.superpopulation for sample in samples], dtype=object)
    superpopulations = sorted({str(label) for label in labels if str(label) and str(label) != "UNK"})
    if len(superpopulations) < min_groups:
        return {"enabled": True, "reason": "fewer_than_min_groups", "superpopulations": superpopulations}

    feature_arrays = [sample.features.numpy() for sample in samples]
    top_eta_rows: List[Dict] = []
    top_delta_rows: List[Dict] = []
    all_count = 0
    all_path = output_dir / "position_superpopulation_effects.csv"
    top_eta_path = output_dir / "top_positions_by_eta_squared.csv"
    top_delta_path = output_dir / "top_positions_by_group_delta.csv"
    fieldnames = [
        "haplotype", "gene", "signal", "position", "n_samples", "n_superpopulations",
        "total_variance", "within_superpopulation_variance", "between_superpopulation_variance",
        "eta_squared", "max_group_mean_delta", "standardized_group_delta",
        "top_superpopulation", "bottom_superpopulation", "top_superpopulation_mean", "bottom_superpopulation_mean",
    ]
    for group in superpopulations:
        fieldnames.extend([f"n_{group}", f"mean_{group}"])

    all_file = None
    all_writer = None
    if write_all:
        all_file, all_writer = _open_csv(all_path, fieldnames)
    try:
        for hap_idx, hap_name in enumerate(("H1", "H2")):
            for gene_idx, gene in enumerate(layout.genes):
                valid_stack = np.stack([_valid_mask(arr, hap_idx, gene_idx, layout) for arr in feature_arrays], axis=0)
                for signal_idx, signal_name in enumerate(layout.signal_names):
                    values = np.stack(
                        [_signal_row(arr, hap_idx, gene_idx, signal_idx, layout) for arr in feature_arrays],
                        axis=0,
                    ).astype(np.float32, copy=False)
                    values = np.where(valid_stack, values, np.nan)
                    for position in range(values.shape[1]):
                        effect = _effect_row(
                            values=values[:, position],
                            labels=labels,
                            superpopulations=superpopulations,
                            min_samples=min_samples,
                            min_groups=min_groups,
                        )
                        if effect is None:
                            continue
                        row = {
                            "haplotype": hap_name,
                            "gene": gene,
                            "signal": signal_name,
                            "position": int(position),
                            **{k: effect[k] for k in fieldnames if k in effect},
                        }
                        for group in superpopulations:
                            stats = effect["group_stats"].get(group, {})
                            row[f"n_{group}"] = int(stats.get("n", 0) or 0)
                            row[f"mean_{group}"] = stats.get("mean")
                        all_count += 1
                        if all_writer is not None:
                            all_writer.writerow(row)
                        top_eta_rows.append(row)
                        top_delta_rows.append(row)
                        _trim_top_rows(top_eta_rows, "eta_squared", top_k)
                        _trim_top_rows(top_delta_rows, "max_group_mean_delta", top_k)
    finally:
        if all_file is not None:
            all_file.close()

    top_eta_rows.sort(key=lambda row: float(row.get("eta_squared") or 0.0), reverse=True)
    top_delta_rows.sort(key=lambda row: float(row.get("max_group_mean_delta") or 0.0), reverse=True)
    for path, rows in ((top_eta_path, top_eta_rows[:top_k]), (top_delta_path, top_delta_rows[:top_k])):
        f, writer = _open_csv(path, fieldnames)
        try:
            for row in rows:
                writer.writerow(row)
        finally:
            f.close()

    return {
        "enabled": True,
        "position_rows_considered": all_count,
        "superpopulations": superpopulations,
        "top_eta_squared": top_eta_rows[:top_k],
        "outputs": {
            "top_positions_by_eta_squared": str(top_eta_path),
            "top_positions_by_group_delta": str(top_delta_path),
            "position_superpopulation_effects": str(all_path) if write_all else None,
        },
    }


def _sparse_pairwise_summary(
    *,
    samples: Sequence[SampleItem],
    layout: ChannelLayout,
    top_effect_rows: Sequence[Dict],
    max_pairs: Optional[int],
) -> Dict:
    if not top_effect_rows:
        return {"enabled": False, "reason": "no_top_effect_rows"}
    gene_to_idx = {gene: idx for idx, gene in enumerate(layout.genes)}
    signal_to_idx = {signal: idx for idx, signal in enumerate(layout.signal_names)}
    hap_to_idx = {"H1": 0, "H2": 1}
    coordinates = []
    for row in top_effect_rows:
        try:
            coordinates.append((
                hap_to_idx[str(row["haplotype"])],
                gene_to_idx[str(row["gene"])],
                signal_to_idx[str(row["signal"])],
                int(row["position"]),
            ))
        except (KeyError, ValueError):
            continue
    if not coordinates:
        return {"enabled": False, "reason": "no_valid_coordinates"}

    feature_arrays = [sample.features.numpy() for sample in samples]
    by_same = {"true": {"mad": [], "rmse": [], "pearson": [], "cosine": []}, "false": {"mad": [], "rmse": [], "pearson": [], "cosine": []}}
    pair_count = 0
    for i, j in itertools.combinations(range(len(samples)), 2):
        if max_pairs is not None and pair_count >= max_pairs:
            break
        xi = []
        xj = []
        ai = feature_arrays[i]
        aj = feature_arrays[j]
        for hap_idx, gene_idx, signal_idx, position in coordinates:
            if not (_valid_mask(ai, hap_idx, gene_idx, layout)[position] and _valid_mask(aj, hap_idx, gene_idx, layout)[position]):
                continue
            xi.append(float(_signal_row(ai, hap_idx, gene_idx, signal_idx, layout)[position]))
            xj.append(float(_signal_row(aj, hap_idx, gene_idx, signal_idx, layout)[position]))
        if len(xi) < 2:
            continue
        x = np.asarray(xi, dtype=np.float64)
        y = np.asarray(xj, dtype=np.float64)
        diff = np.abs(x - y)
        same_key = "true" if samples[i].superpopulation == samples[j].superpopulation else "false"
        by_same[same_key]["mad"].append(float(np.mean(diff)))
        by_same[same_key]["rmse"].append(float(np.sqrt(np.mean((x - y) ** 2))))
        by_same[same_key]["pearson"].append(_pearson(x, y))
        by_same[same_key]["cosine"].append(_cosine(x, y))
        pair_count += 1
    summary = {key: {metric: _mean(vals) for metric, vals in metrics.items()} for key, metrics in by_same.items()}
    within = summary.get("true", {}).get("mad")
    between = summary.get("false", {}).get("mad")
    return {
        "enabled": True,
        "top_positions": len(coordinates),
        "num_pairs_compared": pair_count,
        "same_superpopulation": summary,
        "between_minus_within_mad": None if within is None or between is None else float(between - within),
        "between_within_mad_ratio": None if not within or between is None else float(between / within),
    }


def _permutation_test(pairwise_path: Path, samples: Sequence[SampleItem], permutations: int, seed: int) -> Dict:
    if permutations <= 0:
        return {"enabled": False}
    labels_by_sample = {sample.sample_id: sample.superpopulation for sample in samples}
    sample_ids = [sample.sample_id for sample in samples]
    labels = [labels_by_sample[sample_id] for sample_id in sample_ids]
    rows = []
    with open(pairwise_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                mad = float(row["mad"])
            except (TypeError, ValueError):
                continue
            if math.isfinite(mad):
                rows.append((row["sample_a"], row["sample_b"], mad))
    if not rows:
        return {"enabled": True, "reason": "no_pairwise_rows"}

    def delta(label_map: Dict[str, str]) -> float:
        within = []
        between = []
        for a, b, mad in rows:
            if label_map.get(a) == label_map.get(b):
                within.append(mad)
            else:
                between.append(mad)
        if not within or not between:
            return float("nan")
        return float(np.mean(between) - np.mean(within))

    observed = delta(labels_by_sample)
    rng = np.random.default_rng(seed)
    permuted = []
    for _ in range(permutations):
        shuffled = list(labels)
        rng.shuffle(shuffled)
        value = delta(dict(zip(sample_ids, shuffled)))
        if math.isfinite(value):
            permuted.append(value)
    if not permuted or not math.isfinite(observed):
        return {"enabled": True, "observed_between_minus_within_mad": _safe_float(observed), "reason": "insufficient_permutations"}
    arr = np.asarray(permuted, dtype=np.float64)
    p_value = float((np.sum(arr >= observed) + 1) / (arr.size + 1))
    z_score = float((observed - float(np.mean(arr))) / max(float(np.std(arr)), 1e-12))
    return {
        "enabled": True,
        "permutations": int(arr.size),
        "observed_between_minus_within_mad": float(observed),
        "permuted_mean": float(np.mean(arr)),
        "permuted_std": float(np.std(arr)),
        "p_value_greater_equal": p_value,
        "z_score": z_score,
    }


def run(args: argparse.Namespace) -> int:
    config = load_config(Path(args.config))
    cache_dir = Path(args.cache_dir) if args.cache_dir else get_dataset_cache_dir(config)
    dataset_dir = Path(args.dataset_dir) if args.dataset_dir else Path(config.dataset_input.dataset_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _validate_processed_cache(cache_dir)

    layout = _build_layout(config, cache_dir)
    sample_filter = None
    if args.sample_ids:
        sample_filter = {str(item) for item in args.sample_ids}
    samples = _load_samples(
        cache_dir=cache_dir,
        dataset_dir=dataset_dir,
        splits=args.splits,
        max_samples=args.max_samples,
        sample_ids_filter=sample_filter,
        config=config,
    )
    if len(samples) < 2:
        raise SystemExit("Menos de 2 individuos carregados para comparacao")

    pairwise_fields = [
        "split_a", "split_b", "sample_a", "sample_b", "superpopulation_a", "superpopulation_b",
        "same_superpopulation", "haplotype", "gene", "signal", "n_valid_positions",
        "pearson", "cosine", "mad", "rmse", "max_abs_diff", "fraction_equal_tol",
    ]
    top_fields = [
        "split_a", "split_b", "sample_a", "sample_b", "superpopulation_a", "superpopulation_b",
        "same_superpopulation", "haplotype", "gene", "signal", "position", "abs_diff",
        "value_a", "value_b", "n_valid_positions",
    ]
    pairwise_path = output_dir / "pairwise_valid_signal_similarity.csv"
    top_path = output_dir / "top_valid_signal_differences.csv"
    top_differences: List[Dict] = []

    max_pairs = args.max_pairs
    pair_count = 0
    pair_file, pair_writer = _open_csv(pairwise_path, pairwise_fields)
    try:
        for a, b in itertools.combinations(samples, 2):
            if max_pairs is not None and pair_count >= max_pairs:
                break
            wrote_pair = False
            for row in _pair_rows(a, b, layout, args.min_valid_positions, top_differences, args.top_k):
                pair_writer.writerow(row)
                wrote_pair = True
            if wrote_pair:
                pair_count += 1
    finally:
        pair_file.close()

    top_differences.sort(key=lambda row: float(row["abs_diff"]), reverse=True)
    top_file, top_writer = _open_csv(top_path, top_fields)
    try:
        for row in top_differences[: args.top_k]:
            top_writer.writerow(row)
    finally:
        top_file.close()

    position_effects = _write_position_effects(
        samples=samples,
        layout=layout,
        output_dir=output_dir,
        top_k=args.effect_top_k,
        min_samples=args.effect_min_samples,
        min_groups=args.effect_min_groups,
        write_all=args.write_all_position_effects,
    )
    top_effect_rows = position_effects.get("top_eta_squared", []) if position_effects.get("enabled") else []
    sparse_summary = _sparse_pairwise_summary(
        samples=samples,
        layout=layout,
        top_effect_rows=top_effect_rows,
        max_pairs=max_pairs,
    )
    permutation_summary = _permutation_test(pairwise_path, samples, args.permutations, args.permutation_seed)

    summary = {
        "config": str(Path(args.config).resolve()),
        "cache_dir": str(cache_dir.resolve()),
        "dataset_dir": str(dataset_dir.resolve()),
        "splits": list(args.splits),
        "num_samples_loaded": len(samples),
        "num_pairs_compared": pair_count,
        "layout": {
            "genes": layout.genes,
            "signal_channels_per_gene": layout.signal_channels_per_gene,
            "mask_channels_per_gene": layout.mask_channels_per_gene,
            "channels_per_gene": layout.channels_per_gene,
            "valid_mask_offset": layout.valid_mask_offset,
            "signal_names": layout.signal_names,
        },
        "outputs": {
            "pairwise_valid_signal_similarity": str(pairwise_path),
            "top_valid_signal_differences": str(top_path),
            **position_effects.get("outputs", {}),
        },
        "pairwise_summary": _summarize_pairwise(pairwise_path),
        "position_effects_summary": {
            key: value for key, value in position_effects.items() if key not in {"top_eta_squared", "outputs"}
        },
        "sparse_top_eta_pairwise_summary": sparse_summary,
        "permutation_test": permutation_summary,
    }
    _write_json(output_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare aligned AlphaGenome signal similarity using only positions valid in both individuals."
    )
    parser.add_argument("config", type=Path)
    parser.add_argument("--cache-dir", type=Path, default=None, help="Processed dataset cache dir. Defaults to config-derived cache.")
    parser.add_argument("--dataset-dir", type=Path, default=None, help="Dataset dir for sample metadata. Defaults to config dataset_dir.")
    parser.add_argument("--splits", nargs="+", choices=["train", "val", "test"], default=["train", "val", "test"])
    parser.add_argument("--sample-ids", nargs="*", default=None)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--max-pairs", type=int, default=None)
    parser.add_argument("--min-valid-positions", type=int, default=10)
    parser.add_argument("--top-k", type=int, default=1000)
    parser.add_argument("--effect-top-k", type=int, default=1000, help="Number of top superpopulation-associated positions to rank.")
    parser.add_argument("--effect-min-samples", type=int, default=10, help="Minimum valid samples required for per-position effect metrics.")
    parser.add_argument("--effect-min-groups", type=int, default=2, help="Minimum superpopulation groups required for per-position effect metrics.")
    parser.add_argument("--write-all-position-effects", action="store_true", help="Write all per-position effects, not only top rankings.")
    parser.add_argument("--permutations", type=int, default=0, help="Permutation count for global between-vs-within MAD test.")
    parser.add_argument("--permutation-seed", type=int, default=13)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
