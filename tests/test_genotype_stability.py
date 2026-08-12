from genomics.predictors.genotype_based.experiments.stability import _plan_from_sample_ids, _plan_with_sample_ids
from genomics.predictors.genotype_based.config import PipelineConfig


def _minimal_config():
    return PipelineConfig.model_validate(
        {
            "dataset_input": {
                "dataset_dir": "/tmp/dataset",
                "alphagenome_outputs": ["rna_seq"],
            },
            "output": {"prediction_target": "superpopulation"},
            "model": {"type": "RF"},
            "stability_analysis": {"enabled": True, "random_seed": 13},
        }
    )


def test_shared_stability_plan_roundtrips_sample_ids():
    config = _minimal_config()
    dev_sample_ids = ["A", "B", "C", "D"]
    test_sample_ids = ["T1"]
    plan = [
        {
            "name": "repeat_001",
            "kind": "repeated_random_split",
            "seed": 13,
            "train_idx": [0, 2],
            "val_idx": [1, 3],
        }
    ]

    payload = _plan_with_sample_ids(config, plan, dev_sample_ids, test_sample_ids)
    restored = _plan_from_sample_ids(payload, dev_sample_ids)

    assert payload["runs"][0]["train_sample_ids"] == ["A", "C"]
    assert payload["runs"][0]["val_sample_ids"] == ["B", "D"]
    assert payload["test_sample_ids"] == ["T1"]
    assert restored == plan
