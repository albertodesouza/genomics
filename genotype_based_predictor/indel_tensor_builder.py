from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from genotype_based_predictor.dynamic_indel_alignment import DynamicIndelAligner


@dataclass(frozen=True)
class CenterWindowSlice:
    center_ref_start: int
    center_ref_end: int
    center_expanded_start: int
    center_expanded_end: int

    @property
    def center_ref_length(self) -> int:
        return self.center_ref_end - self.center_ref_start

    @property
    def center_expanded_length(self) -> int:
        return self.center_expanded_end - self.center_expanded_start


def compute_center_window_slice(
    gene_mapping: Dict[str, object],
    ref_length: int,
    center_window_size: int,
    haplotypes: Sequence[str] = ("H1", "H2"),
) -> CenterWindowSlice:
    center_half = center_window_size // 2
    center_ref_start = max(0, (ref_length // 2) - center_half)
    center_ref_end = min(ref_length, center_ref_start + center_window_size)

    expanded_indices_covering_center: List[int] = []
    for sample_payload in gene_mapping.get("samples", {}).values():
        for haplotype in haplotypes:
            hap_entry = sample_payload.get(haplotype)
            if hap_entry is None:
                continue
            for source_idx, target_idx in zip(hap_entry.get("copy_from_indices", []), hap_entry.get("expanded_indices", [])):
                if center_ref_start <= int(source_idx) < center_ref_end:
                    expanded_indices_covering_center.append(int(target_idx))
            for target_idx in hap_entry.get("deletion_indices", []):
                expanded_indices_covering_center.append(int(target_idx))

    if not expanded_indices_covering_center:
        raise RuntimeError("Nao foi possivel localizar indices expandidos para a janela central.")

    return CenterWindowSlice(
        center_ref_start=center_ref_start,
        center_ref_end=center_ref_end,
        center_expanded_start=min(expanded_indices_covering_center),
        center_expanded_end=max(expanded_indices_covering_center) + 1,
    )


def build_aligned_haplotype_tensor(
    row: np.ndarray,
    entry: Dict[str, object],
    expanded_length: int,
    neutral_value: float = 0.0,
    include_valid_mask: bool = True,
    expanded_slice: Optional[Tuple[int, int]] = None,
) -> np.ndarray:
    if expanded_slice is None:
        slice_start = 0
        slice_end = expanded_length
    else:
        slice_start, slice_end = expanded_slice
    local_length = max(int(slice_end) - int(slice_start), 0)

    values = np.full(local_length, float(neutral_value), dtype=np.float32)
    valid_mask = np.zeros(local_length, dtype=np.float32)
    ins_mask = np.zeros(local_length, dtype=np.float32)
    del_mask = np.zeros(local_length, dtype=np.float32)

    copy_from = entry.get("copy_from_indices", [])
    copy_to = entry.get("expanded_indices", [])
    if len(copy_from) != len(copy_to):
        raise ValueError("Mapeamento dinamico invalido: copy_from_indices e expanded_indices com comprimentos diferentes.")

    row_len = len(row)
    for source_idx, target_idx in zip(copy_from, copy_to):
        source_idx = int(source_idx)
        target_idx = int(target_idx)
        if 0 <= source_idx < row_len and slice_start <= target_idx < slice_end:
            local_idx = target_idx - slice_start
            values[local_idx] = row[source_idx]
            valid_mask[local_idx] = 1.0

    for target_idx in entry.get("insertion_indices", []):
        target_idx = int(target_idx)
        if slice_start <= target_idx < slice_end:
            ins_mask[target_idx - slice_start] = 1.0

    for target_idx in entry.get("deletion_indices", []):
        target_idx = int(target_idx)
        if slice_start <= target_idx < slice_end:
            del_mask[target_idx - slice_start] = 1.0

    output_rows = [values, ins_mask, del_mask]
    if include_valid_mask:
        output_rows.append(valid_mask)
    return np.vstack(output_rows)


def build_global_aligned_tensor(
    *,
    dataset,
    aligner: DynamicIndelAligner,
    gene_name: str,
    sample_ids: Sequence[str],
    haplotypes: Sequence[str],
    output_name: str,
    track_index: int,
    neutral_value: float,
    include_valid_mask: bool,
    center_window_size: Optional[int] = None,
) -> Tuple[np.ndarray, Dict[str, object], List[Dict[str, object]]]:
    gene_mapping = aligner._load_gene_mapping(gene_name)
    example_window = dataset._load_window_data(sample_ids[0], gene_name)
    example_prediction = example_window[f"predictions_{haplotypes[0].lower()}"][output_name]
    ref_length = int(example_prediction.shape[0])
    expanded_length = aligner.get_expanded_length(gene_name)

    expanded_slice = None
    center_slice_meta = None
    local_length = expanded_length
    if center_window_size is not None:
        center_slice = compute_center_window_slice(
            gene_mapping=gene_mapping,
            ref_length=ref_length,
            center_window_size=int(center_window_size),
            haplotypes=haplotypes,
        )
        expanded_slice = (center_slice.center_expanded_start, center_slice.center_expanded_end)
        local_length = center_slice.center_expanded_length
        center_slice_meta = {
            "center_ref_start": center_slice.center_ref_start,
            "center_ref_end": center_slice.center_ref_end,
            "center_expanded_start": center_slice.center_expanded_start,
            "center_expanded_end": center_slice.center_expanded_end,
            "center_ref_length": center_slice.center_ref_length,
            "center_expanded_length": center_slice.center_expanded_length,
        }

    channel_count = 4 if include_valid_mask else 3
    tensor = np.zeros((len(sample_ids), len(haplotypes), channel_count, local_length), dtype=np.float32)
    summary_rows: List[Dict[str, object]] = []

    sample_id_to_index = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    haplotype_to_index = {haplotype: idx for idx, haplotype in enumerate(haplotypes)}

    for sample_id in sample_ids:
        window_data = dataset._load_window_data(sample_id, gene_name)
        for haplotype in haplotypes:
            entry = aligner.get_haplotype_entry(gene_name, sample_id, haplotype)
            if entry is None:
                continue
            prediction = window_data[f"predictions_{haplotype.lower()}"][output_name]
            if prediction.ndim == 1:
                row = np.asarray(prediction, dtype=np.float32)
            else:
                row = np.asarray(prediction[:, track_index], dtype=np.float32)

            aligned = build_aligned_haplotype_tensor(
                row=row,
                entry=entry,
                expanded_length=expanded_length,
                neutral_value=neutral_value,
                include_valid_mask=include_valid_mask,
                expanded_slice=expanded_slice,
            )
            tensor[sample_id_to_index[sample_id], haplotype_to_index[haplotype]] = aligned

            summary_rows.append({
                "sample_id": sample_id,
                "haplotype": haplotype,
                "copied_positions": int(aligned[3].sum()) if include_valid_mask else int(np.count_nonzero(aligned[0] != neutral_value)),
                "insertion_positions": int(aligned[1].sum()),
                "deletion_positions": int(aligned[2].sum()),
                "valid_positions": int(aligned[3].sum()) if include_valid_mask else None,
                "tensor_shape": tuple(aligned.shape),
            })

    metadata = {
        "expanded_length": expanded_length,
        "local_length": local_length,
        "include_valid_mask": include_valid_mask,
        "center_slice": center_slice_meta,
        "sample_ids": list(sample_ids),
        "haplotypes": list(haplotypes),
    }
    return tensor, metadata, summary_rows
