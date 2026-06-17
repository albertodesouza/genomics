import csv
import json
from argparse import Namespace

import torch

from genomics.predictors.genotype_based.analysis.compare_aligned_signals import run


def _write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def test_compare_aligned_signals_uses_pairwise_valid_positions(tmp_path):
    cache_dir = tmp_path / "cache"
    dataset_dir = tmp_path / "dataset"
    output_dir = tmp_path / "out"

    _write_json(cache_dir / "split_index.json", {"train": ["S1", "S2"], "val": [], "test": []})
    _write_json(
        cache_dir / "shards_index.json",
        {"train": {"shard_size": 64, "num_samples": 2, "shards": ["train_data_shard_00000.pt"]}},
    )
    _write_json(
        dataset_dir / "dataset_metadata.json",
        {
            "individuals_pedigree": {
                "S1": {"superpopulation": "AFR", "population": "YRI"},
                "S2": {"superpopulation": "EUR", "population": "GBR"},
            }
        },
    )

    # Shape: (haplotypes=2, channels=signal+ins+del+valid, length=4).
    a = torch.zeros((2, 4, 4), dtype=torch.float32)
    b = torch.zeros((2, 4, 4), dtype=torch.float32)
    a[:, 3, :] = 1.0
    b[:, 3, :] = 1.0
    b[:, 3, 2] = 0.0
    a[:, 0, :] = torch.tensor([1.0, 2.0, 1000.0, 4.0])
    b[:, 0, :] = torch.tensor([1.0, 3.0, -1000.0, 7.0])
    torch.save([(a, torch.tensor(0)), (b, torch.tensor(1))], cache_dir / "train_data_shard_00000.pt")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "dataset_input:",
                f"  dataset_dir: {dataset_dir}",
                "  alphagenome_outputs: [rna_seq]",
                "  genes_to_use: [MC1R]",
                "  tensor_layout: haplotype_channels",
                "  feature_mode: signals_and_masks",
                "  indel_include_valid_mask: true",
                "  consensus_dataset_dir: /tmp/unused",
                "output:",
                "  prediction_target: superpopulation",
            ]
        )
    )

    rc = run(
        Namespace(
            config=config_path,
            cache_dir=cache_dir,
            dataset_dir=dataset_dir,
            splits=["train"],
            sample_ids=None,
            max_samples=None,
            max_pairs=None,
            min_valid_positions=1,
            top_k=10,
            effect_top_k=10,
            effect_min_samples=1,
            effect_min_groups=2,
            write_all_position_effects=True,
            permutations=5,
            permutation_seed=13,
            output_dir=output_dir,
        )
    )

    assert rc == 0
    with open(output_dir / "pairwise_valid_signal_similarity.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 2
    assert {row["haplotype"] for row in rows} == {"H1", "H2"}
    assert {int(row["n_valid_positions"]) for row in rows} == {3}
    assert {float(row["max_abs_diff"]) for row in rows} == {3.0}

    with open(output_dir / "top_valid_signal_differences.csv") as f:
        top_rows = list(csv.DictReader(f))
    assert {int(row["position"]) for row in top_rows} == {3}

    with open(output_dir / "top_positions_by_eta_squared.csv") as f:
        effect_rows = list(csv.DictReader(f))
    assert effect_rows
    assert int(effect_rows[0]["position"]) in {1, 3}
    assert float(effect_rows[0]["eta_squared"]) == 1.0

    summary = json.loads((output_dir / "summary.json").read_text())
    assert summary["sparse_top_eta_pairwise_summary"]["enabled"] is True
    assert summary["permutation_test"]["enabled"] is True
