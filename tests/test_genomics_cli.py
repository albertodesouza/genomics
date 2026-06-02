import subprocess
import sys
from pathlib import Path

import yaml

import genomics_cli


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def test_genotype_config_overrides_preserve_yaml_defaults(monkeypatch, tmp_path):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path / "data"))
    config_path = tmp_path / "config.yaml"
    _write_yaml(
        config_path,
        {
            "dataset_input": {
                "dataset_id": "1kg_high_coverage",
                "results_dir": "custom/runs",
                "processed_cache_dir": "custom/cache",
                "consensus_dataset_dir": "custom/consensus",
            }
        },
    )

    args = genomics_cli.build_parser().parse_args(["genotype", "train", str(config_path)])
    temp_config = genomics_cli._genotype_config_with_overrides(args)
    payload = yaml.safe_load(temp_config.read_text(encoding="utf-8"))
    dataset_input = payload["dataset_input"]

    assert dataset_input["results_dir"] == "custom/runs"
    assert dataset_input["processed_cache_dir"] == "custom/cache"
    assert dataset_input["consensus_dataset_dir"] == "custom/consensus"
    assert dataset_input["dataset_dir"] == str((tmp_path / "data" / "v1" / "1kG_high_coverage").resolve())


def test_genotype_config_explicit_overrides_replace_yaml(monkeypatch, tmp_path):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path / "data"))
    config_path = tmp_path / "config.yaml"
    _write_yaml(
        config_path,
        {
            "dataset_input": {
                "dataset_id": "1kg_high_coverage",
                "results_dir": "custom/runs",
                "processed_cache_dir": "custom/cache",
                "consensus_dataset_dir": "custom/consensus",
            }
        },
    )

    args = genomics_cli.build_parser().parse_args(
        [
            "genotype",
            "train",
            str(config_path),
            "--results-dir",
            str(tmp_path / "override-runs"),
            "--processed-cache-dir",
            str(tmp_path / "override-cache"),
            "--consensus-dataset-dir",
            str(tmp_path / "override-consensus"),
        ]
    )
    temp_config = genomics_cli._genotype_config_with_overrides(args)
    payload = yaml.safe_load(temp_config.read_text(encoding="utf-8"))
    dataset_input = payload["dataset_input"]

    assert dataset_input["results_dir"] == str((tmp_path / "override-runs").resolve())
    assert dataset_input["processed_cache_dir"] == str((tmp_path / "override-cache").resolve())
    assert dataset_input["consensus_dataset_dir"] == str((tmp_path / "override-consensus").resolve())


def test_variant_config_overrides_preserve_yaml_defaults(tmp_path):
    config_path = tmp_path / "config.yaml"
    _write_yaml(
        config_path,
        {
            "dataset": {
                "processed_dir": "custom/processed",
                "results_dir": "custom/runs",
            }
        },
    )

    args = genomics_cli.build_parser().parse_args(["variant", "train", str(config_path)])
    temp_config = genomics_cli._variant_config_with_overrides(args)
    payload = yaml.safe_load(temp_config.read_text(encoding="utf-8"))

    assert payload["dataset"]["processed_dir"] == "custom/processed"
    assert payload["dataset"]["results_dir"] == "custom/runs"


def test_audit_configs_has_no_active_legacy_references():
    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, "-m", "genomics_cli", "audit-configs", "--fail-on-active-legacy"],
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr
