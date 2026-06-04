from __future__ import annotations

import json
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class ExperimentRun:
    pipeline: str
    name: str
    root_dir: Path
    run_dir: Path
    models_dir: Path
    logs_dir: Path
    plots_dir: Path
    reports_dir: Path
    manifest_path: Path


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _git_commit() -> Optional[str]:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parents[1],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return None


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, default=_json_default)


def _read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    return payload if isinstance(payload, dict) else {}


def setup_experiment_run(
    *,
    pipeline: str,
    run_name: str,
    runs_root: Path,
    config_path: Path,
    resolved_config: Dict[str, Any],
    extra: Optional[Dict[str, Any]] = None,
) -> ExperimentRun:
    runs_root = Path(runs_root)
    run_dir = runs_root / run_name
    models_dir = run_dir / "models"
    logs_dir = run_dir / "logs"
    plots_dir = run_dir / "plots"
    reports_dir = run_dir / "reports"
    for path in (models_dir, logs_dir, plots_dir, reports_dir):
        path.mkdir(parents=True, exist_ok=True)

    config_path = Path(config_path).resolve()
    shutil.copyfile(config_path, run_dir / "config.yaml")
    _write_json(run_dir / "resolved_config.json", resolved_config)

    manifest = {
        "schema_version": 1,
        "pipeline": pipeline,
        "run_name": run_name,
        "run_dir": str(run_dir.resolve()),
        "created_at": _utc_now(),
        "config_path": str(config_path),
        "artifacts": {
            "config": "config.yaml",
            "resolved_config": "resolved_config.json",
            "models_dir": "models",
            "logs_dir": "logs",
            "plots_dir": "plots",
            "reports_dir": "reports",
        },
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "git_commit": _git_commit(),
        },
        "status": "initialized",
    }
    if extra:
        manifest["extra"] = extra
    manifest_path = run_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    return ExperimentRun(
        pipeline=pipeline,
        name=run_name,
        root_dir=runs_root,
        run_dir=run_dir,
        models_dir=models_dir,
        logs_dir=logs_dir,
        plots_dir=plots_dir,
        reports_dir=reports_dir,
        manifest_path=manifest_path,
    )


def update_manifest(run_dir: Path, **updates: Any) -> None:
    manifest_path = Path(run_dir) / "manifest.json"
    manifest = _read_json(manifest_path)
    manifest.update(updates)
    manifest["updated_at"] = _utc_now()
    _write_json(manifest_path, manifest)
