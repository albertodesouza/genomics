import json

from genomics_pipeline.experiment import setup_experiment_run, update_manifest


def test_setup_experiment_run_creates_standard_layout(tmp_path):
    config_path = tmp_path / "source.yaml"
    config_path.write_text("answer: 42\n", encoding="utf-8")

    run = setup_experiment_run(
        pipeline="unit_pipeline",
        run_name="unit_run",
        runs_root=tmp_path / "runs",
        config_path=config_path,
        resolved_config={"answer": 42},
        extra={"purpose": "test"},
    )

    assert run.run_dir == tmp_path / "runs" / "unit_run"
    assert run.models_dir.is_dir()
    assert run.logs_dir.is_dir()
    assert run.plots_dir.is_dir()
    assert run.reports_dir.is_dir()
    assert (run.run_dir / "config.yaml").read_text(encoding="utf-8") == "answer: 42\n"

    resolved = json.loads((run.run_dir / "resolved_config.json").read_text(encoding="utf-8"))
    assert resolved == {"answer": 42}

    manifest = json.loads(run.manifest_path.read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 1
    assert manifest["pipeline"] == "unit_pipeline"
    assert manifest["run_name"] == "unit_run"
    assert manifest["status"] == "initialized"
    assert manifest["extra"] == {"purpose": "test"}
    assert manifest["artifacts"]["models_dir"] == "models"
    assert "python" in manifest["environment"]


def test_update_manifest_preserves_existing_fields(tmp_path):
    config_path = tmp_path / "source.yaml"
    config_path.write_text("answer: 42\n", encoding="utf-8")
    run = setup_experiment_run(
        pipeline="unit_pipeline",
        run_name="unit_run",
        runs_root=tmp_path / "runs",
        config_path=config_path,
        resolved_config={"answer": 42},
    )

    update_manifest(run.run_dir, status="completed", best_accuracy=0.5)

    manifest = json.loads(run.manifest_path.read_text(encoding="utf-8"))
    assert manifest["pipeline"] == "unit_pipeline"
    assert manifest["status"] == "completed"
    assert manifest["best_accuracy"] == 0.5
    assert "updated_at" in manifest
