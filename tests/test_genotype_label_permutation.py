from collections import Counter

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
