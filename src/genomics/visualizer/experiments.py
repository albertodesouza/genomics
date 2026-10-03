"""Experiment run discovery and summaries for the visualizer's Experiments page.

A run is any directory under a runs root containing ``manifest.json``, ``config.yaml`` or
``*_results.json``. Metrics are read generically: every numeric top-level value of each results
file (``weighted_accuracy``, ``accuracy``, ``f1``...) is exposed as ``<file stem>.<key>``.
Summaries are cached by file modification times so listing hundreds of runs stays fast.
"""
from __future__ import annotations

import json
import math
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

MAX_JSON_BYTES = 4_000_000
MAX_TEXT_BYTES = 256_000
IMAGE_TYPES = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".svg": "image/svg+xml", ".gif": "image/gif"}
TEXT_TYPES = {".yaml", ".yml", ".json", ".txt", ".log", ".md", ".csv", ".tsv"}


def _load_json(path: Path) -> Any:
    if path.stat().st_size > MAX_JSON_BYTES:
        raise ValueError(f"{path.name} is too large to preview")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _finite(value: Any) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _numeric_metrics(data: Any) -> Dict[str, float]:
    if not isinstance(data, dict):
        return {}
    out = {}
    for key, value in data.items():
        number = _finite(value)
        if number is not None:
            out[str(key)] = number
    return out


class ExperimentService:
    def __init__(self, roots: List[Path]):
        self.roots = [Path(r).expanduser().resolve() for r in roots]
        self._summaries: Dict[str, Tuple[Tuple, Dict[str, Any]]] = {}
        self._lock = threading.Lock()

    def _is_run(self, path: Path) -> bool:
        return (path / "manifest.json").exists() or (path / "config.yaml").exists() or any(path.glob("*_results.json"))

    def run_dirs(self) -> Dict[str, Path]:
        runs: Dict[str, Path] = {}
        for root in self.roots:
            if not root.is_dir():
                continue
            for child in sorted(root.iterdir()):
                if not child.is_dir():
                    continue
                if self._is_run(child):
                    runs[self._run_id(root, child)] = child
                else:
                    for grandchild in sorted(child.iterdir()):
                        if grandchild.is_dir() and self._is_run(grandchild):
                            runs[self._run_id(root, grandchild)] = grandchild
        return runs

    def _run_id(self, root: Path, run: Path) -> str:
        rel = run.relative_to(root).as_posix()
        return rel if len(self.roots) == 1 else f"{root.name}/{rel}"

    def _resolve(self, run_id: str) -> Path:
        runs = self.run_dirs()
        path = runs.get(run_id)
        if path is None:
            raise KeyError(f"Unknown run: {run_id}")
        return path

    def _signature(self, path: Path) -> Tuple:
        items = [path.stat().st_mtime_ns]
        for name in ("manifest.json", "config.yaml"):
            p = path / name
            items.append(p.stat().st_mtime_ns if p.exists() else 0)
        items.extend((p.name, p.stat().st_mtime_ns) for p in sorted(path.glob("*_results.json")))
        return tuple(items)

    def summary(self, run_id: str, path: Path) -> Dict[str, Any]:
        signature = self._signature(path)
        with self._lock:
            cached = self._summaries.get(run_id)
        if cached is not None and cached[0] == signature:
            return cached[1]
        manifest: Dict[str, Any] = {}
        if (path / "manifest.json").exists():
            try:
                manifest = _load_json(path / "manifest.json")
            except Exception:
                manifest = {}
        metrics: Dict[str, float] = {}
        results = []
        for result_path in sorted(path.glob("*_results.json")):
            stem = result_path.name[: -len("_results.json")]
            try:
                data = _load_json(result_path)
                file_metrics = _numeric_metrics(data)
            except Exception as exc:
                results.append({"file": result_path.name, "error": str(exc)})
                continue
            results.append({"file": result_path.name, "stem": stem, "metrics": file_metrics})
            for key, value in file_metrics.items():
                metrics[f"{stem}.{key}"] = value
        for key in ("best_val_accuracy", "best_val_loss", "last_epoch"):
            number = _finite(manifest.get(key))
            if number is not None:
                metrics[key] = number
        models_dir = path / "models"
        checkpoints = sorted(p.name for p in models_dir.glob("*.pt")) if models_dir.is_dir() else []
        summary = {
            "id": run_id,
            "name": path.name,
            "group": path.parent.name if path.parent not in self.roots else None,
            "path": str(path),
            "status": manifest.get("status"),
            "pipeline": manifest.get("pipeline"),
            "created_at": manifest.get("created_at"),
            "updated_at": manifest.get("updated_at") or manifest.get("created_at"),
            "mtime": path.stat().st_mtime,
            "metrics": metrics,
            "results": results,
            "checkpoints": checkpoints,
            "has_history": (models_dir / "training_history.json").exists(),
        }
        with self._lock:
            self._summaries[run_id] = (signature, summary)
        return summary

    def list(self) -> Dict[str, Any]:
        runs = [self.summary(run_id, path) for run_id, path in self.run_dirs().items()]
        metric_names: Dict[str, int] = {}
        for run in runs:
            for key in run["metrics"]:
                metric_names[key] = metric_names.get(key, 0) + 1
        return {
            "roots": [str(r) for r in self.roots],
            "runs": runs,
            "metrics": sorted(metric_names, key=lambda k: (-metric_names[k], k)),
        }

    def detail(self, run_id: str) -> Dict[str, Any]:
        path = self._resolve(run_id)
        summary = dict(self.summary(run_id, path))
        results = []
        for result_path in sorted(path.glob("*_results.json")):
            try:
                data = _load_json(result_path)
            except Exception as exc:
                results.append({"file": result_path.name, "error": str(exc)})
                continue
            if not isinstance(data, dict):
                continue
            results.append(
                {
                    "file": result_path.name,
                    "metrics": _numeric_metrics(data),
                    "confusion_matrix": data.get("confusion_matrix"),
                    "class_names": list((data.get("per_class_metrics") or {}).keys()) or data.get("class_names"),
                    "per_class_metrics": data.get("per_class_metrics"),
                    "classification_report": data.get("classification_report"),
                }
            )
        summary["results_detail"] = results
        history_path = path / "models" / "training_history.json"
        if history_path.exists():
            try:
                history = _load_json(history_path)
                summary["history"] = {
                    k: [(_finite(x) if not isinstance(x, bool) else None) for x in v]
                    for k, v in history.items()
                    if isinstance(v, list) and v and all(isinstance(x, (int, float)) for x in v)
                }
            except Exception as exc:
                summary["history_error"] = str(exc)
        config_path = path / "config.yaml"
        if config_path.exists():
            summary["config"] = config_path.read_text(encoding="utf-8", errors="replace")[:MAX_TEXT_BYTES]
        manifest_path = path / "manifest.json"
        if manifest_path.exists():
            try:
                summary["manifest"] = _load_json(manifest_path)
            except Exception:
                pass
        files = []
        for sub in ("plots", "reports", "."):
            directory = path / sub
            if not directory.is_dir():
                continue
            for file in sorted(directory.iterdir()):
                if file.is_file() and (file.suffix.lower() in IMAGE_TYPES or (sub != "." and file.suffix.lower() in TEXT_TYPES)):
                    files.append({"path": file.relative_to(path).as_posix(), "size": file.stat().st_size, "image": file.suffix.lower() in IMAGE_TYPES})
        summary["files"] = files
        return summary

    def file(self, run_id: str, relative: str) -> Tuple[bytes, str]:
        run = self._resolve(run_id)
        target = (run / relative).resolve()
        if run not in target.parents or not target.is_file():
            raise KeyError("File not found")
        suffix = target.suffix.lower()
        if suffix in IMAGE_TYPES:
            return target.read_bytes(), IMAGE_TYPES[suffix]
        if suffix in TEXT_TYPES:
            return target.read_bytes()[:MAX_TEXT_BYTES], "text/plain; charset=utf-8"
        raise KeyError("Unsupported file type")
