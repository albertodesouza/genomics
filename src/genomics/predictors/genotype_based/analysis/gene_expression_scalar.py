"""PrediXcan-style scalar gene-expression features.

Reduces each individual's per-position AlphaGenome-derived tracks (the same processed tensor
consumed by CNN2AncestryPredictor/NNAncestryPredictor/the sklearn baselines) to one scalar per
(gene, ontology term, strand), by summing (or averaging) the already log-normalized per-track
values over a gene's window. This mirrors classical TWAS methods (PrediXcan/S-PrediXcan
"predicted expression -> phenotype" features) far more literally than feeding a full positional
tensor to a model, and is meant for training a plain linear/logistic classifier on top -- see
notebooks/pigmentation_alphagenome_cnn_variant_scoring.ipynb Section 9.

Reuses the same processed-tensor pipeline (`ProcessedGenomicDataset` via
`_make_runtime_processed_datasets`) so results are directly comparable to the CNN/NN/sklearn
arms trained on the full tensor from the same config, rather than a separately-computed feature
set.

Haplotypes (H1/H2) can be kept as two separate features per track (`haplotype_mode="separate"`,
the default -- matching how the full-tensor CNN/NN arms see them, as two distinct rows) or
merged into one "diploid total" feature per track (`haplotype_mode="sum"`). Merging always sums
rather than averages the two haplotypes regardless of `reduction`: `reduction` controls how a
single haplotype's signal is aggregated *within* a gene window (total vs. average signal across
position), while haplotype merging represents each haplotype's independent contribution to total
expression adding up -- a real biological assumption (roughly: total transcript count = maternal
allele's contribution + paternal allele's contribution) that need not hold for genes with strong
allele-specific or dominant/recessive regulatory effects, hence being opt-in rather than default.
"""

from __future__ import annotations

from typing import Any, List, Tuple

import numpy as np
from torch.utils.data import Subset

from genomics.predictors.genotype_based.analysis.interpretability import build_haplotype_track_layout

_HAPLOTYPE_MODES = ("separate", "sum")


def scalar_feature_names(config: Any, haplotype_mode: str = "separate") -> List[str]:
    """One name per scalar feature, in the same order `reduce_tensor_to_scalars` returns them."""
    if haplotype_mode not in _HAPLOTYPE_MODES:
        raise ValueError(f"Unsupported haplotype_mode '{haplotype_mode}'; choose one of {_HAPLOTYPE_MODES}")
    layout = build_haplotype_track_layout(config)
    base_names = [
        f"{row['gene_name']}__{row.get('ontology_curie') or row.get('biosample_name') or 'na'}__{row['strand']}"
        for row in layout
        if row["track_type"] == "signal"
    ]
    if haplotype_mode == "sum":
        return base_names
    return [f"H1__{name}" for name in base_names] + [f"H2__{name}" for name in base_names]


def reduce_tensor_to_scalars(features: np.ndarray, config: Any, reduction: str = "sum", haplotype_mode: str = "separate") -> np.ndarray:
    """Reduce one individual's `(2, rows, L)` haplotype-major tensor to a scalar feature vector,
    aggregating over the position axis (`reduction`: "sum" = total signal, "mean" = average
    signal per position) and either keeping both haplotypes as separate features
    (`haplotype_mode="separate"`, default -- returns `2 * n_signal_rows` values, H1 block then
    H2 block) or summing them into one "diploid total" feature per track
    (`haplotype_mode="sum"` -- returns `n_signal_rows` values; the sum is always a sum,
    independent of `reduction`, since it represents each haplotype's contribution adding up to
    the individual's total).

    `features` is expected already log-normalized, matching what
    `ProcessedGenomicDataset`/`CachedProcessedDataset.__getitem__` returns -- i.e. this sums
    per-track log1p-scaled values, not raw AlphaGenome signal. Only rows with
    `track_type == "signal"` are reduced; INDEL mask rows (present only when
    `feature_mode != "signals_only"`) are dropped, since a mask isn't an expression-like
    quantity to sum.
    """
    if reduction not in ("sum", "mean"):
        raise ValueError(f"Unsupported reduction '{reduction}'; choose 'sum' or 'mean'")
    if haplotype_mode not in _HAPLOTYPE_MODES:
        raise ValueError(f"Unsupported haplotype_mode '{haplotype_mode}'; choose one of {_HAPLOTYPE_MODES}")
    layout = build_haplotype_track_layout(config)
    signal_row_idx = [row["row_index"] for row in layout if row["track_type"] == "signal"]
    reduce_fn = np.sum if reduction == "sum" else np.mean
    per_haplotype = reduce_fn(features[:, signal_row_idx, :], axis=-1)  # (2, n_signal_rows)
    if haplotype_mode == "separate":
        return per_haplotype.reshape(-1)  # H1 block (n_signal_rows,) then H2 block (n_signal_rows,)
    return per_haplotype.sum(axis=0)  # (n_signal_rows,) -- H1 + H2, always summed


def _sample_id_for_subset_item(subset: Subset, position: int) -> str:
    processed_ds = subset.dataset
    proc_idx = subset.indices[position]
    base_idx = processed_ds.valid_sample_indices[proc_idx]
    return processed_ds._sample_id_for_base_index(base_idx)


def materialize_scalar_split(
    dataset: Subset,
    config: Any,
    reduction: str = "sum",
    haplotype_mode: str = "separate",
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Iterate one train/val/test split (a `torch.utils.data.Subset` over a
    `ProcessedGenomicDataset`, as returned per-split by `_make_runtime_processed_datasets`) and
    reduce every individual's tensor to the scalar feature vector from
    `reduce_tensor_to_scalars`.

    Returns `(X, y, sample_ids)`: `X` is `(n_samples, n_features)`, `y` is the integer target
    class index (matching the model's `idx_to_target`), and `sample_ids[i]` is the 1000
    Genomes sample ID for row `i`.
    """
    names = scalar_feature_names(config, haplotype_mode=haplotype_mode)
    n = len(dataset)
    X = np.zeros((n, len(names)), dtype=np.float64)
    y = np.zeros(n, dtype=np.int64)
    sample_ids: List[str] = []
    for i in range(n):
        features_tensor, target_tensor = dataset[i]
        X[i] = reduce_tensor_to_scalars(features_tensor.numpy(), config, reduction=reduction, haplotype_mode=haplotype_mode)
        y[i] = int(target_tensor.item())
        sample_ids.append(_sample_id_for_subset_item(dataset, i))
    return X, y, sample_ids
