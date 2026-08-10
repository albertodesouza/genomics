"""Tests PigmentationModelContext.build_individual_tensor/run_cnn_logits/predict_proba against a
fake ProcessedGenomicDataset and a tiny torch model, bypassing the real (expensive, bcftools- and
checkpoint-dependent) __init__ -- this environment doesn't have `bcftools` installed, so building a
genuine on-disk runtime dataset via `_make_runtime_processed_datasets` isn't possible here. Instead
these tests target the context's own wiring: gene/haplotype iteration, override substitution, and
tensor stacking/normalization/inference, using collaborators shaped exactly like the real ones.
"""
from types import SimpleNamespace

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from genomics.predictors.genotype_based.analysis.pigmentation_model_context import PigmentationModelContext
from genomics.predictors.genotype_based.data.normalization import apply_normalization


class _FakeProcessedDataset:
    """Stands in for ProcessedGenomicDataset._process_window_haplotype_channels: in this fake, the
    "raw prediction" the caller already resolved (real cache or an edited override) is treated as
    the finished per-haplotype signal tensor, so this test targets PigmentationModelContext's own
    gene/haplotype/override wiring rather than the real bcftools_chain alignment math."""

    def _process_window_haplotype_channels(self, sample_id, gene, haplotype, predictions, prediction_metadata=None):
        return predictions["rna_seq"], None


class _TinyModel(torch.nn.Module):
    def __init__(self, in_features: int, num_classes: int):
        super().__init__()
        torch.manual_seed(0)
        self.linear = torch.nn.Linear(in_features, num_classes)

    def forward(self, x):
        return self.linear(x.reshape(x.shape[0], -1))


GENES = ["GENE1", "GENE2"]
CHANNELS_PER_GENE = 2
WINDOW = 8


def _build_context(canned_predictions):
    ctx = object.__new__(PigmentationModelContext)
    ctx._torch = torch
    ctx._apply_normalization = apply_normalization
    ctx.config = SimpleNamespace(
        dataset_input=SimpleNamespace(
            genes_to_use=GENES,
            ontology_terms=["CL:1000458", "CL:0000346"],
            window_center_size=WINDOW,
        )
    )
    ctx.full_ds = _FakeProcessedDataset()
    ctx.normalization_params = {"method": "zscore", "mean": 0.0, "std": 1.0}  # identity
    ctx.device = torch.device("cpu")
    in_features = 2 * CHANNELS_PER_GENE * len(GENES) * WINDOW
    ctx.model = _TinyModel(in_features, num_classes=2)
    ctx.class_names = ["weak pigmentation", "strong pigmentation"]
    ctx.strong_pigmentation_idx = 1
    ctx._baseline_tensor_cache = {}
    ctx.load_raw_prediction = lambda sample_id, gene, haplotype, output="rna_seq": canned_predictions[(gene, haplotype)]
    return ctx


def _canned(value: float):
    return np.full((CHANNELS_PER_GENE, WINDOW), value, dtype=np.float32), []


def test_build_individual_tensor_shape_and_gene_ordering():
    canned = {
        ("GENE1", "H1"): _canned(1.0),
        ("GENE1", "H2"): _canned(2.0),
        ("GENE2", "H1"): _canned(3.0),
        ("GENE2", "H2"): _canned(4.0),
    }
    ctx = _build_context(canned)

    tensor = ctx.build_individual_tensor("SAMPLE1")

    assert tensor.shape == (2, CHANNELS_PER_GENE * len(GENES), WINDOW)
    assert torch.allclose(tensor[0, 0:2, :], torch.full((2, WINDOW), 1.0))
    assert torch.allclose(tensor[0, 2:4, :], torch.full((2, WINDOW), 3.0))
    assert torch.allclose(tensor[1, 0:2, :], torch.full((2, WINDOW), 2.0))
    assert torch.allclose(tensor[1, 2:4, :], torch.full((2, WINDOW), 4.0))


def test_build_individual_tensor_overrides_only_targeted_gene_haplotype():
    canned = {
        ("GENE1", "H1"): _canned(1.0),
        ("GENE1", "H2"): _canned(2.0),
        ("GENE2", "H1"): _canned(3.0),
        ("GENE2", "H2"): _canned(4.0),
    }
    ctx = _build_context(canned)
    baseline = ctx.build_individual_tensor("SAMPLE1")

    override_array, override_meta = _canned(99.0)
    edited = ctx.build_individual_tensor("SAMPLE1", overrides={("GENE1", "H1"): (override_array, override_meta)})

    assert torch.allclose(edited[0, 0:2, :], torch.full((2, WINDOW), 99.0))
    assert torch.allclose(edited[0, 2:4, :], baseline[0, 2:4, :])
    assert torch.allclose(edited[1, :, :], baseline[1, :, :])


def test_run_cnn_logits_and_predict_proba_shapes():
    canned = {
        ("GENE1", "H1"): _canned(1.0),
        ("GENE1", "H2"): _canned(2.0),
        ("GENE2", "H1"): _canned(3.0),
        ("GENE2", "H2"): _canned(4.0),
    }
    ctx = _build_context(canned)
    tensor = ctx.build_individual_tensor("SAMPLE1")

    logits = ctx.run_cnn_logits(tensor)
    assert logits.shape == (2,)
    assert torch.isfinite(logits).all()

    proba = ctx.predict_proba(tensor)
    assert proba.shape == (2,)
    assert abs(float(proba.sum()) - 1.0) < 1e-5

    prob_strong = ctx.strong_pigmentation_probability(tensor)
    assert 0.0 <= prob_strong <= 1.0
    assert prob_strong == pytest.approx(float(proba[ctx.strong_pigmentation_idx]))


def test_editing_a_gene_changes_the_prediction():
    canned = {
        ("GENE1", "H1"): _canned(1.0),
        ("GENE1", "H2"): _canned(2.0),
        ("GENE2", "H1"): _canned(3.0),
        ("GENE2", "H2"): _canned(4.0),
    }
    ctx = _build_context(canned)
    baseline_tensor = ctx.build_individual_tensor("SAMPLE1")
    baseline_prob = ctx.strong_pigmentation_probability(baseline_tensor)

    override_array, override_meta = _canned(-50.0)
    edited_tensor = ctx.build_individual_tensor("SAMPLE1", overrides={("GENE2", "H2"): (override_array, override_meta)})
    edited_prob = ctx.strong_pigmentation_probability(edited_tensor)

    assert baseline_prob != pytest.approx(edited_prob)


def test_baseline_tensor_is_cached_per_sample():
    canned = {
        ("GENE1", "H1"): _canned(1.0),
        ("GENE1", "H2"): _canned(2.0),
        ("GENE2", "H1"): _canned(3.0),
        ("GENE2", "H2"): _canned(4.0),
    }
    ctx = _build_context(canned)
    assert ctx._baseline_tensor_cache == {}

    first = ctx.baseline_tensor("SAMPLE1")
    assert "SAMPLE1" in ctx._baseline_tensor_cache
    second = ctx.baseline_tensor("SAMPLE1")
    assert torch.equal(first, second)
