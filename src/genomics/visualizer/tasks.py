"""Detached tasks: dataset imports, AlphaGenome predictions, training and evaluation.

Unlike :mod:`genomics.visualizer.jobs` (threads inside the server for interactive work), a task is
a separate process tree started through :mod:`genomics.visualizer.task_runner` in its own session:
it keeps running when the browser closes or the visualizer exits, and a restarted visualizer
finds it again. Everything about a task lives in ``<tasks_dir>/<task id>/``::

    task.json    what to run (steps, working directory, parameters shown in the UI)
    state.json   written by the runner: status, progress, message, exit code, result
    log.txt      combined output of every step
    runner.log   the runner's own stderr (normally empty)

Secrets are never written there: steps inherit the server's environment (API key, AlphaGenome
server address) when they are started.
"""
from __future__ import annotations

import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from collections import deque
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from genomics.visualizer.task_runner import process_start_time

ACTIVE = ("starting", "queued", "running")
FINISHED = ("done", "failed", "cancelled", "lost")
STARTUP_GRACE_SECONDS = 15.0
MAX_LOG_LINES = 2000

Hook = Callable[[Dict[str, Any]], None]


def _read_json(path: Path) -> Dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


class TaskManager:
    def __init__(self, root: Path, python: Optional[str] = None):
        self.root = Path(root).expanduser().resolve()
        self.python = python or sys.executable
        self._lock = threading.Lock()
        self._hooks: Dict[str, List[Hook]] = {}
        self._handled: set = set()

    # -- hooks ---------------------------------------------------------------------------------
    def on_finished(self, kind: str, hook: Hook) -> None:
        """Call ``hook(task)`` once when a task of ``kind`` stops (recorded in the task dir, so
        a restarted visualizer does not repeat it)."""
        self._hooks.setdefault(kind, []).append(hook)

    def _run_hooks(self, task: Dict[str, Any]) -> None:
        if task["status"] not in FINISHED or task["id"] in self._handled:
            return
        self._handled.add(task["id"])
        marker = Path(task["dir"]) / ".hooks_done"
        if marker.exists() or not self._hooks.get(task["kind"]):
            return
        for hook in self._hooks.get(task["kind"], []):
            try:
                hook(task)
            except Exception as exc:  # a hook must not break listing
                print(f"visualizer: task hook for {task['id']} failed: {exc}", file=sys.stderr)
        try:
            marker.touch()
        except OSError:
            pass

    # -- creation ------------------------------------------------------------------------------
    def create(
        self,
        kind: str,
        title: str,
        steps: List[Dict[str, Any]],
        params: Optional[Dict[str, Any]] = None,
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
        resource: Optional[str] = None,
        files: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """Write the task and start its runner. ``files`` are extra text files for the task dir;
        steps can refer to them as ``{task_dir}/<name>``."""
        self.root.mkdir(parents=True, exist_ok=True)
        task_id = f"{time.strftime('%Y%m%d-%H%M%S')}-{kind}-{uuid.uuid4().hex[:6]}"
        task_dir = self.root / task_id
        task_dir.mkdir()
        for name, text in (files or {}).items():
            (task_dir / name).write_text(text, encoding="utf-8")
        resolved_steps = []
        for step in steps:
            step = dict(step)
            step["command"] = [str(c).replace("{task_dir}", str(task_dir)) for c in step["command"]]
            if step.get("history_glob"):
                step["history_glob"] = str(step["history_glob"]).replace("{task_dir}", str(task_dir))
            resolved_steps.append(step)
        spec = {
            "id": task_id,
            "kind": kind,
            "title": title,
            "created": time.time(),
            "steps": resolved_steps,
            "cwd": str(cwd) if cwd else None,
            "params": params or {},
            "resource": resource,
            "locks_dir": str(self.root / "locks"),
        }
        (task_dir / "task.json").write_text(json.dumps(spec, indent=2, default=str), encoding="utf-8")
        (task_dir / "state.json").write_text(json.dumps({"status": "starting", "progress": 0.0, "message": "Starting"}), encoding="utf-8")
        child_env = dict(os.environ if env is None else env)
        child_env.setdefault("PYTHONUNBUFFERED", "1")
        child_env.setdefault("COLUMNS", "160")
        # Tools installed in the interpreter's environment (bcftools, samtools) even when the
        # visualizer was started without activating it.
        child_env["PATH"] = os.pathsep.join([str(Path(self.python).parent), child_env.get("PATH", "")])
        with open(task_dir / "runner.log", "w", encoding="utf-8") as runner_log:
            subprocess.Popen(
                [self.python, "-m", "genomics.visualizer.task_runner", str(task_dir)],
                cwd=str(cwd) if cwd else None,
                env=child_env,
                stdin=subprocess.DEVNULL,
                stdout=runner_log,
                stderr=subprocess.STDOUT,
                start_new_session=True,  # survives the visualizer and its terminal
                close_fds=True,
            )
        return self.get(task_id)

    # -- reading -------------------------------------------------------------------------------
    def _dir(self, task_id: str) -> Path:
        path = (self.root / task_id).resolve()
        if path.parent != self.root or not (path / "task.json").exists():
            raise KeyError(f"Unknown task: {task_id}")
        return path

    def _runner_alive(self, state: Dict[str, Any]) -> bool:
        pid = state.get("runner_pid")
        if not pid:
            return False
        started = process_start_time(int(pid))
        if started is None:
            return False
        expected = state.get("runner_start")
        return started == -1 or expected in (None, -1) or started == expected

    def _describe(self, task_dir: Path) -> Dict[str, Any]:
        spec = _read_json(task_dir / "task.json")
        state = _read_json(task_dir / "state.json")
        status = state.get("status", "starting")
        if status in ACTIVE and not self._runner_alive(state):
            # The runner died without recording an outcome (killed, machine rebooted).
            if status != "starting" or time.time() - float(spec.get("created") or 0) > STARTUP_GRACE_SECONDS:
                status = "lost"
                state.update(status="lost", message=state.get("message") or "The task process is gone", finished=state.get("finished") or time.time())
                try:
                    (task_dir / "state.json").write_text(json.dumps(state, indent=2, default=str), encoding="utf-8")
                except OSError:
                    pass
        started = state.get("started") or spec.get("created")
        finished = state.get("finished")
        return {
            "id": spec.get("id", task_dir.name),
            "kind": spec.get("kind"),
            "title": spec.get("title"),
            "params": spec.get("params") or {},
            "created": spec.get("created"),
            "started": state.get("started"),
            "finished": finished,
            "elapsed": round((finished or time.time()) - float(started or time.time()), 1),
            "status": status,
            "progress": float(state.get("progress") or 0.0),
            "message": state.get("message") or "",
            "last_line": state.get("last_line") or "",
            "step": state.get("step"),
            "steps": len(spec.get("steps") or []),
            "step_title": state.get("step_title"),
            "exit_code": state.get("exit_code"),
            "result": state.get("result") or {},
            "resource": spec.get("resource"),
            "dir": str(task_dir),
        }

    def get(self, task_id: str, log_lines: int = 0) -> Dict[str, Any]:
        task_dir = self._dir(task_id)
        task = self._describe(task_dir)
        if log_lines:
            task["log"] = self.log_tail(task_id, log_lines)
            spec = _read_json(task_dir / "task.json")
            task["commands"] = [{"title": s.get("title"), "command": s.get("command")} for s in spec.get("steps") or []]
        self._run_hooks(task)
        return task

    def list(self, kind: Optional[str] = None, limit: int = 200) -> List[Dict[str, Any]]:
        if not self.root.is_dir():
            return []
        tasks = []
        for task_dir in sorted((p for p in self.root.iterdir() if (p / "task.json").exists()), key=lambda p: p.name, reverse=True)[:limit]:
            task = self._describe(task_dir)
            if kind and task["kind"] != kind:
                continue
            self._run_hooks(task)
            tasks.append(task)
        return tasks

    def active(self) -> List[Dict[str, Any]]:
        return [t for t in self.list() if t["status"] in ACTIVE]

    def log_tail(self, task_id: str, lines: int = 200) -> List[str]:
        path = self._dir(task_id) / "log.txt"
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                return [line.rstrip("\n") for line in deque(handle, maxlen=max(1, min(lines, MAX_LOG_LINES)))]
        except OSError:
            return []

    def log_path(self, task_id: str) -> Path:
        return self._dir(task_id) / "log.txt"

    # -- control -------------------------------------------------------------------------------
    def cancel(self, task_id: str) -> Dict[str, Any]:
        task_dir = self._dir(task_id)
        state = _read_json(task_dir / "state.json")
        if state.get("status") in ACTIVE and self._runner_alive(state):
            try:
                os.killpg(int(state["runner_pid"]), signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                pass
        return self.get(task_id)

    def delete(self, task_id: str) -> None:
        task_dir = self._dir(task_id)
        task = self._describe(task_dir)
        if task["status"] in ACTIVE:
            raise RuntimeError("Cancel the task before removing it")
        shutil.rmtree(task_dir)
