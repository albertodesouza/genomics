import subprocess
import sys
import json
from pathlib import Path

import yaml

from genomics import cli as genomics_cli


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
        [sys.executable, "-m", "genomics", "audit-configs", "--fail-on-active-legacy"],
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_audit_data_reports_dataset_status(monkeypatch, tmp_path, capsys):
    data_root = tmp_path / "data"
    dataset_dir = data_root / "v1" / "1kG_high_coverage"
    dataset_dir.mkdir(parents=True)
    (dataset_dir / "dataset_metadata.json").write_text(
        json.dumps({"individuals": ["HG00096"], "genes": ["DDB1"], "window_catalog": {"DDB1": {}}}),
        encoding="utf-8",
    )
    (dataset_dir / "layout_metadata.json").write_text(json.dumps({"layout_version": 1}), encoding="utf-8")
    (dataset_dir / "individuals" / "HG00096" / "windows" / "DDB1").mkdir(parents=True)
    ref_window_dir = dataset_dir / "references" / "windows" / "DDB1"
    ref_window_dir.mkdir(parents=True)
    (ref_window_dir / "ref.window.fa").write_text(">DDB1\nACGT\n", encoding="utf-8")
    (ref_window_dir / "window_metadata.json").write_text(json.dumps({"chromosome": "chr11"}), encoding="utf-8")
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(data_root))

    rc = genomics_cli.main(["audit-data", "--dataset-id", "1kg_high_coverage", "--json", "--fail-on-missing"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert rc == 0
    assert payload[0]["status"] == "ok"
    assert payload[0]["individual_count"] == 1
    assert payload[0]["gene_count"] == 1


def test_audit_data_detects_missing_bcftools_chain_artifacts(monkeypatch, tmp_path, capsys):
    data_root = tmp_path / "data"
    dataset_dir = data_root / "v1" / "1kG_high_coverage"
    window_dir = dataset_dir / "individuals" / "HG00096" / "windows" / "DDB1"
    window_dir.mkdir(parents=True)
    ref_window_dir = dataset_dir / "references" / "windows" / "DDB1"
    ref_window_dir.mkdir(parents=True)
    (ref_window_dir / "ref.window.fa").write_text(">DDB1\nACGT\n", encoding="utf-8")
    (ref_window_dir / "window_metadata.json").write_text(json.dumps({"chromosome": "chr11"}), encoding="utf-8")
    (dataset_dir / "dataset_metadata.json").write_text(
        json.dumps({"individuals": ["HG00096"], "genes": ["DDB1"], "window_catalog": {"DDB1": {}}}),
        encoding="utf-8",
    )
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(data_root))

    rc = genomics_cli.main([
        "audit-data",
        "--dataset-id",
        "1kg_high_coverage",
        "--check-bcftools-chain",
        "--sample-limit",
        "1",
        "--json",
        "--fail-on-missing",
    ])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert rc == 2
    assert payload[0]["status"] == "error"
    assert payload[0]["bcftools_chain"]["missing"] == len(genomics_cli.BCFTOOLS_CHAIN_REQUIRED_TEMPLATES)
    assert str(window_dir) in payload[0]["bcftools_chain"]["examples"][0]


def test_completion_bash_outputs_completion_script(capsys):
    rc = genomics_cli.main(["completion", "bash"])
    captured = capsys.readouterr()

    assert rc == 0
    assert "complete -F _genomics_completion genomics" in captured.out
    assert "genotype" in captured.out
    assert "test" in captured.out


def test_genotype_test_command_parses():
    args = genomics_cli.build_parser().parse_args([
        "genotype",
        "test",
        "configs/predictors/genotype_based/default.yaml",
    ])

    assert args.genotype_command == "test"
    assert args.checkpoint == "best_accuracy"


def test_config_describe_genotype_lists_known_fields(capsys):
    rc = genomics_cli.main(["config", "describe", "genotype"])
    captured = capsys.readouterr()

    assert rc == 0
    assert "Genotype-Based Predictor config schema" in captured.out
    assert "dataset_input.dataset_id" in captured.out
    assert "training.batch_size" in captured.out


def test_config_validate_infers_genotype_kind(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("GENOMICS_DATA_ROOT", str(tmp_path / "data"))
    config_dir = tmp_path / "configs" / "predictors" / "genotype_based"
    config_dir.mkdir(parents=True)
    config_path = config_dir / "minimal.yaml"
    _write_yaml(
        config_path,
        {
            "dataset_input": {
                "dataset_id": "1kg_high_coverage",
                "alphagenome_outputs": ["rna_seq"],
            },
            "output": {"prediction_target": "superpopulation"},
        },
    )

    rc = genomics_cli.main(["config", "validate", str(config_path)])
    captured = capsys.readouterr()

    assert rc == 0
    assert "ok:" in captured.out
    assert "(genotype)" in captured.out


def test_config_validate_reports_errors(tmp_path, capsys):
    config_path = tmp_path / "bad.yaml"
    _write_yaml(config_path, {"dataset": {"processed_dir": "processed"}, "model": {"d_model": 255, "heads": 8}})

    rc = genomics_cli.main(["config", "validate", str(config_path), "--kind", "variant"])
    captured = capsys.readouterr()

    assert rc == 2
    assert "config validation failed" in captured.err
