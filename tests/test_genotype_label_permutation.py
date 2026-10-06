from collections import Counter

import pytest
import torch

from genomics.predictors.genotype_based.config import PipelineConfig
from genomics.predictors.genotype_based.data.processed_dataset import ProcessedGenomicDataset


class _DummyBaseDataset:
    def __init__(self):
        self.dataset_metadata = {
            "individuals": ["S1", "S2", "S3", "S4", "S5"],
            "individuals_pedigree": {
                "S1": {"superpopulation": "AFR"},
                "S2": {"superpopulation": "EUR"},
                "S3": {"superpopulation": "EAS"},
                "S4": {"superpopulation": "EUR"},
                "S5": {"superpopulation": "AFR"},
            },
        }

    def __len__(self):
        return len(self.dataset_metadata["individuals"])

    def __getitem__(self, idx):
        sample_id = self.dataset_metadata["individuals"][idx]
        return {}, self.dataset_metadata["individuals_pedigree"][sample_id]


def _config(tmp_path):
    return PipelineConfig.model_validate(
        {
            "dataset_input": {
                "dataset_dir": str(tmp_path),
                "alphagenome_outputs": ["rna_seq"],
                "tensor_layout": "raw_center_crop",
            },
            "output": {
                "prediction_target": "superpopulation",
                "known_classes": ["AFR", "EAS", "EUR"],
            },
            "label_permutation": {
                "enabled": True,
                "random_seed": 13,
            },
        }
    )


def test_label_permutation_is_reproducible_and_preserves_distribution(tmp_path):
    ds1 = ProcessedGenomicDataset(
        _DummyBaseDataset(),
        _config(tmp_path),
        normalization_params={"mean": 0.0, "std": 1.0},
        compute_normalization=False,
    )
    ds2 = ProcessedGenomicDataset(
        _DummyBaseDataset(),
        _config(tmp_path),
        normalization_params={"mean": 0.0, "std": 1.0},
        compute_normalization=False,
    )

    original = [
        ds1.dataset_metadata["individuals_pedigree"][sample_id]["superpopulation"]
        for sample_id in ds1.dataset_metadata["individuals"]
    ]
    permuted = [
        ds1.permuted_targets_by_sample_id[sample_id]
        for sample_id in ds1.dataset_metadata["individuals"]
    ]

    assert ds1.permuted_targets_by_sample_id == ds2.permuted_targets_by_sample_id
    assert Counter(permuted) == Counter(original)

    target = ds1._build_target_tensor(
        {"superpopulation": "AFR"},
        torch.zeros(1),
        sample_id="S1",
    )
    assert int(target.item()) == ds1.target_to_idx[ds1.permuted_targets_by_sample_id["S1"]]


class _StratifiedBaseDataset(_DummyBaseDataset):
    """Two superpopulations whose pigmentation class is perfectly determined by ancestry."""

    def __init__(self):
        super().__init__()
        people = {
            "S1": ("AFR", "strong"), "S2": ("AFR", "strong"), "S3": ("AFR", "strong"),
            "S4": ("EUR", "weak"), "S5": ("EUR", "weak"), "S6": ("EUR", "strong"),
        }
        self.dataset_metadata = {
            "individuals": sorted(people),
            "individuals_pedigree": {s: {"superpopulation": sp, "pigmentation": p} for s, (sp, p) in people.items()},
        }


def _stratified_config(tmp_path, stratify_field, seed=13):
    config = _config(tmp_path)
    config.output.prediction_target = "pigmentation"
    config.output.known_classes = ["strong", "weak"]
    config.label_permutation.stratify_field = stratify_field
    config.label_permutation.random_seed = seed
    return config


def _permute(tmp_path, stratify_field, seed=13):
    ds = ProcessedGenomicDataset(
        _StratifiedBaseDataset(), _stratified_config(tmp_path, stratify_field, seed),
        normalization_params={"mean": 0.0, "std": 1.0}, compute_normalization=False,
    )
    return ds, ds.permuted_targets_by_sample_id


def test_permutation_within_a_field_keeps_each_group_class_distribution(tmp_path):
    ds, permuted = _permute(tmp_path, "superpopulation")
    pedigree = ds.dataset_metadata["individuals_pedigree"]
    assert set(permuted) == set(pedigree)
    assert Counter(permuted.values()) == Counter(r["pigmentation"] for r in pedigree.values())
    # Labels only move inside a superpopulation, so each group keeps its own class counts:
    # the group -> class association (the confound) survives, individual-level signal does not.
    for group in ("AFR", "EUR"):
        members = [s for s, r in pedigree.items() if r["superpopulation"] == group]
        assert Counter(permuted[s] for s in members) == Counter(pedigree[s]["pigmentation"] for s in members)
    assert _permute(tmp_path, "superpopulation")[1] == permuted  # reproducible

    # AFR is all one class, so permuting inside it cannot change a label; EUR is mixed and does.
    assert all(permuted[s] == "strong" for s in ("S1", "S2", "S3"))
    assert {permuted[s] for s in ("S4", "S5", "S6")} == {"weak", "strong"}

    # A global permutation is free to move a "weak" label into AFR, whatever the seed; a
    # stratified one never can, which is the whole point of the control.
    afr = ("S1", "S2", "S3")
    seeds = range(6)
    assert not any(_permute(tmp_path, "superpopulation", k)[1][s] == "weak" for k in seeds for s in afr)
    assert any(_permute(tmp_path, None, k)[1][s] == "weak" for k in seeds for s in afr)
    assert Counter(_permute(tmp_path, None)[1].values()) == Counter(permuted.values())


def test_permutation_within_a_missing_field_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="stratify_field"):
        _permute(tmp_path, "no_such_field")
    with pytest.raises(ValueError, match="stratify_field"):
        PipelineConfig.model_validate({
            "dataset_input": {"dataset_dir": str(tmp_path), "alphagenome_outputs": ["rna_seq"], "tensor_layout": "raw_center_crop"},
            "output": {"prediction_target": "superpopulation"},
            "label_permutation": {"enabled": True, "stratify_field": "   "},
        })
