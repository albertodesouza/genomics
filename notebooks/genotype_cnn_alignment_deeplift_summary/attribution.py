"""DeepLIFT class-mean attribution and per-gene track-curve extraction."""

import numpy as np
import torch

from genomics.predictors.genotype_based.analysis.interpretability import (
    DeepLIFT,
    build_haplotype_track_layout,
)

_class_input_stack_cache = {}


def compute_class_mean_deeplift(model, dataset, class_indices, class_names, model_name="model"):
    """DeepLIFT (Gradient x (input - baseline), baseline_type="mean") class-mean attribution for
    each class in `class_indices`. Using the dataset-mean input as baseline (rather than an
    all-zeros tensor) makes delta = input - baseline reflect deviation from a typical sample,
    instead of collapsing to gradient(x) * x -- dominated by raw input magnitude rather than
    class-discriminative signal.

    Returns {class_idx: {"class_name":..., "mean_attr":..., "mean_input":..., "n":...}}.
    """
    deeplift = DeepLIFT(model)
    results = {}
    for class_idx in class_indices:
        mean_attr, mean_input, n = deeplift.generate_class_mean(class_idx, dataset=dataset, baseline_type="mean")
        results[class_idx] = {
            "class_name": class_names[class_idx], "mean_attr": mean_attr, "mean_input": mean_input, "n": n,
        }
        print(f"{model_name} / {class_names[class_idx]}: DeepLIFT computed over {n} test samples")
    return results


def compute_population_mean_input(deeplift_results, class_indices):
    """Population-wide mean input as the n-weighted combination of the per-class means -- no
    extra data pass needed."""
    n_by_class = {c: deeplift_results[c]["n"] for c in class_indices}
    n_total = sum(n_by_class.values())
    return sum(deeplift_results[c]["mean_input"] * n_by_class[c] for c in n_by_class) / n_total


def compute_logodds_attributions(model, dataset, target_idx, other_idx, deeplift_results, class_names):
    """Per-population attribution toward log-odds(target/other) = logit_target - logit_other.

    `deeplift_results[c]["mean_attr"]` (from `compute_class_mean_deeplift`) only explains samples
    of class c toward class c's own output logit. DeepLIFT's Rescale multipliers at each
    nonlinearity depend only on the forward activations for a given input/baseline pair, not on
    which output neuron backprop started from -- so attribution is linear in the output-side
    target: attribution(A - B) = attribution(A) - attribution(B) exactly, no new interpretation
    needed. That combination needs each population's attribution toward *both* logits, and
    `compute_class_mean_deeplift`'s cache only has each population's attribution toward its *own*
    logit, so this takes two new DeepLIFT passes for the missing half via
    `DeepLIFT.generate_class_mean_cross`.

    Returns (logodds_attr_target_pop, logodds_attr_other_pop).
    """
    cross_deeplift = DeepLIFT(model)
    cross_attr_other_toward_target, _, n_other_cross = cross_deeplift.generate_class_mean_cross(
        other_idx, target_idx, dataset=dataset, baseline_type="mean",
    )
    cross_attr_target_toward_other, _, n_target_cross = cross_deeplift.generate_class_mean_cross(
        target_idx, other_idx, dataset=dataset, baseline_type="mean",
    )
    print(
        f"cross attribution computed (intermediate, combined into log-odds below): "
        f"{class_names[other_idx]} samples -> {class_names[target_idx]} logit "
        f"(n={n_other_cross}); {class_names[target_idx]} samples -> "
        f"{class_names[other_idx]} logit (n={n_target_cross})"
    )

    logodds_attr_target_pop = deeplift_results[target_idx]["mean_attr"] - cross_attr_target_toward_other
    logodds_attr_other_pop = cross_attr_other_toward_target - deeplift_results[other_idx]["mean_attr"]
    return logodds_attr_target_pop, logodds_attr_other_pop


def gene_track_curves(tensor, config, gene):
    '''Per-(ontology_curie, strand) curve for one gene from a (haplotype, row, L)-shaped tensor
    (either a signed DeepLIFT mean_attr or a mean_input), averaged over both haplotypes --
    {(ontology_curie, strand): 1D curve}, same length as the window. For haplotype_channels
    models `layout` has one row per (gene, track) and `tensor` carries an explicit haplotype axis
    (size 2), so the haplotype average happens on that axis. For no_alignment (raw_center_crop),
    `layout` instead carries two rows per (gene, track) -- one tagged haplotype="H1", one "H2" --
    since H1/H2 are stacked as extra rows rather than a tensor axis (see
    `build_haplotype_track_layout`'s docstring); those two rows are averaged together below so
    both haplotypes' signal (not just H1's) contributes to the curve.'''
    layout = build_haplotype_track_layout(config)
    rows = [r for r in layout if r["gene_name"] == gene and r["track_type"] == "signal"]
    t = tensor.detach().cpu()
    if t.ndim == 2:
        t = t.unsqueeze(0)
    grouped = {}
    for r in rows:
        key = (r["ontology_curie"], r["strand"])
        grouped.setdefault(key, []).append(t[:, r["row_index"], :].mean(dim=0))
    return {key: torch.stack(curves).mean(dim=0).numpy() for key, curves in grouped.items()}


def smooth(x, window=513):
    if window <= 1:
        return x
    kernel = np.ones(window) / window
    return np.convolve(x, kernel, mode="same")


def collect_class_input_stack(loaders, name, class_idx):
    """Every test-split sample's raw input tensor (no gradients) for one model/class, stacked
    on a new leading sample axis -- the same per-class population `compute_class_mean_deeplift`
    already averages down to `mean_input`, kept here unaveraged so per-position std across
    individuals can be computed for hatched-band plots."""
    key = (name, class_idx)
    if key not in _class_input_stack_cache:
        dataset = loaders[name]["test"].dataset
        samples = [dataset[i][0] for i in range(len(dataset)) if int(dataset[i][1]) == class_idx]
        _class_input_stack_cache[key] = torch.stack(samples)
    return _class_input_stack_cache[key]


def class_track_curves_mean_std(loaders, name, class_idx, gene, config):
    """Per-(ontology_curie, strand) (mean_curve, std_curve) pair across the individuals of one
    class, haplotype-averaged per individual first (matching `gene_track_curves`'s convention)
    then reduced over the sample axis -- {(ontology_curie, strand): (mean_curve, std_curve)}.

    For haplotype_channels models `stack` carries an explicit haplotype axis (dim 1, size 2) that
    is averaged away up front. For no_alignment (raw_center_crop) there's no such axis -- H1/H2
    are stacked as extra rows instead (`layout` has two rows per (gene, track), tagged
    haplotype="H1"/"H2") -- so those per-individual row pairs are averaged together per key below
    instead, keeping both haplotypes' signal rather than only H1's.
    """
    stack = collect_class_input_stack(loaders, name, class_idx)
    if stack.ndim == 3:  # (n, row, L): no_alignment, no explicit haplotype axis
        stack = stack.unsqueeze(1)
    arr = stack.numpy()  # (n, hap_axis, row, L); hap_axis is a stand-in size-1 axis for no_alignment
    layout = build_haplotype_track_layout(config)
    rows = [r for r in layout if r["gene_name"] == gene and r["track_type"] == "signal"]
    grouped = {}
    for r in rows:
        key = (r["ontology_curie"], r["strand"])
        per_individual = arr[:, :, r["row_index"], :].mean(axis=1)  # (n, L)
        grouped.setdefault(key, []).append(per_individual)
    result = {}
    for key, curves in grouped.items():
        combined = np.stack(curves).mean(axis=0)  # (n, L): haplotype-averaged per individual
        result[key] = (combined.mean(axis=0), combined.std(axis=0))
    return result
