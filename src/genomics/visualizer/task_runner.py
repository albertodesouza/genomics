"""Runs one visualizer task detached from the server: ``python -m genomics.visualizer.task_runner DIR``.

The server starts this in its own session, so a task keeps running when the browser tab closes or
the visualizer itself exits. The runner executes the steps of ``DIR/task.json`` in order (stopping
at the first failure), appends their output to ``DIR/log.txt`` and keeps ``DIR/state.json`` up to
date: status, progress, message and exit code. Steps report progress by printing
``@@progress <fraction> <message>`` lines and structured results with ``@@result <json object>``;
a training step can instead name its ``training_history.json`` so progress follows the epochs.

A task with a ``resource`` waits for an exclusive lock on ``<locks_dir>/<resource>.lock`` first
(status ``queued``), so for example only one GPU job runs at a time. The lock is released by the
kernel when the runner exits, whatever the reason.
"""
from __future__ import annotations

import fcntl
import glob
import json
import os
import shlex
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, Optional

PROGRESS_PREFIX = "@@progress "
RESULT_PREFIX = "@@result "


def write_json_atomic(path: Path, payload: Any) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def process_start_time(pid: int) -> Optional[int]:
    """Kernel start time of ``pid`` (clock ticks since boot), to detect PID reuse; None if gone."""
    try:
        with open(f"/proc/{pid}/stat", "rb") as handle:
            stat = handle.read()
    except OSError:
        try:
            os.kill(pid, 0)
        except (ProcessLookupError, PermissionError):
            return None
        return -1  # alive, start time unknown (non-Linux)
    # The command name (field 2) may contain spaces; fields after it are space separated.
    fields = stat[stat.rfind(b")") + 2:].split()
    return int(fields[19]) if len(fields) > 19 else -1


class Runner:
    def __init__(self, task_dir: Path):
        self.dir = Path(task_dir)
        self.spec: Dict[str, Any] = json.loads((self.dir / "task.json").read_text(encoding="utf-8"))
        self.state_path = self.dir / "state.json"
        try:
            self.state: Dict[str, Any] = json.loads(self.state_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            self.state = {}
        self.child: Optional[subprocess.Popen] = None
        self.cancelled = False
        self._lock = threading.Lock()
        self._last_save = 0.0
        self._step_base = 0.0
        self._step_share = 1.0

    # -- state ---------------------------------------------------------------------------
    def save(self, force: bool = False) -> None:
        with self._lock:
            now = time.time()
            if not force and now - self._last_save < 0.5:
                return
            self.state["updated"] = now
            write_json_atomic(self.state_path, self.state)
            self._last_save = now

    def progress(self, fraction: float, message: Optional[str] = None) -> None:
        fraction = max(0.0, min(1.0, float(fraction)))
        self.state["progress"] = round(self._step_base + self._step_share * fraction, 4)
        if message is not None:
            self.state["message"] = message[:300]
        self.save()

    def _on_signal(self, _signum, _frame) -> None:
        self.cancelled = True
        self.state["message"] = "Cancelling…"
        child = self.child
        if child is not None and child.poll() is None:
            try:
                child.terminate()
            except ProcessLookupError:
                pass

    # -- resource lock -----------------------------------------------------------------------
    def _older_waiting(self, resource: str) -> bool:
        """Whether an older task for the same resource is still waiting (keeps the queue FIFO)."""
        for sibling in self.dir.parent.iterdir():
            if sibling.name >= self.dir.name or not (sibling / "task.json").exists():
                continue
            try:
                spec = json.loads((sibling / "task.json").read_text(encoding="utf-8"))
                state = json.loads((sibling / "state.json").read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if spec.get("resource") != resource or state.get("status") not in ("starting", "queued"):
                continue
            pid = state.get("runner_pid")
            if pid and process_start_time(int(pid)) is not None:
                return True
        return False

    def _acquire(self, resource: str):
        locks_dir = Path(self.spec.get("locks_dir") or self.dir.parent / "locks")
        locks_dir.mkdir(parents=True, exist_ok=True)
        handle = open(locks_dir / f"{resource}.lock", "a+")
        announced = False
        while not self.cancelled:
            if not self._older_waiting(resource):
                try:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    return handle
                except BlockingIOError:
                    pass
            if not announced:
                self.state.update(status="queued", message=f"Waiting for another {resource} job to finish")
                self.save(force=True)
                announced = True
            time.sleep(1.0)
        handle.close()
        return None

    # -- output parsing ------------------------------------------------------------------------
    def _handle_line(self, line: str) -> None:
        if line.startswith(PROGRESS_PREFIX):
            parts = line[len(PROGRESS_PREFIX):].strip().split(" ", 1)
            try:
                self.progress(float(parts[0]), parts[1] if len(parts) > 1 else None)
            except ValueError:
                pass
        elif line.startswith(RESULT_PREFIX):
            try:
                payload = json.loads(line[len(RESULT_PREFIX):])
            except ValueError:
                return
            if isinstance(payload, dict):
                self.state.setdefault("result", {}).update(payload)
                self.save(force=True)
        elif line.strip():
            self.state["last_line"] = line.strip()[:300]
            self.save()

    def _pump(self, stream, log) -> None:
        """Copy child output to the log; carriage-return progress bars keep only their last frame."""
        buffer = b""
        while True:
            chunk = os.read(stream.fileno(), 65536)
            if not chunk:
                break
            buffer += chunk
            *lines, buffer = buffer.split(b"\n")
            for raw in lines:
                text = raw.rsplit(b"\r", 1)[-1] if b"\r" in raw.rstrip(b"\r") else raw.rstrip(b"\r")
                line = text.decode("utf-8", "replace")
                if not line.startswith("@@"):
                    log.write(line + "\n")
                self._handle_line(line)
            log.flush()
        if buffer:
            line = buffer.rsplit(b"\r", 1)[-1].decode("utf-8", "replace")
            if not line.startswith("@@"):
                log.write(line + "\n")
            self._handle_line(line)
            log.flush()

    def _watch_history(self, step: Dict[str, Any], stop: threading.Event) -> None:
        """Training progress from the epochs recorded in training_history.json."""
        pattern = step.get("history_glob")
        epochs = int(step.get("epochs") or 0)
        if not pattern or epochs <= 0:
            return
        while not stop.wait(5.0):
            for path in sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)[:1]:
                try:
                    history = json.loads(Path(path).read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    continue
                done = len(history.get("epoch") or history.get("train_loss") or [])
                if done:
                    message = f"Epoch {done}/{epochs}"
                    for key in ("val_accuracy", "train_loss"):
                        values = [v for v in (history.get(key) or []) if isinstance(v, (int, float))]
                        if values:
                            message += f" · {key} {values[-1]:.4g}"
                    self.progress(done / epochs, message)

    # -- main ----------------------------------------------------------------------------------
    def run(self) -> int:
        signal.signal(signal.SIGTERM, self._on_signal)
        signal.signal(signal.SIGINT, self._on_signal)
        signal.signal(signal.SIGHUP, signal.SIG_IGN)
        self.state.update(
            status="starting",
            runner_pid=os.getpid(),
            runner_start=process_start_time(os.getpid()),
            progress=0.0,
            message="Starting",
        )
        self.save(force=True)
        lock = None
        resource = self.spec.get("resource")
        if resource:
            lock = self._acquire(str(resource))
        steps = self.spec.get("steps") or []
        total_weight = sum(float(s.get("weight") or 1.0) for s in steps) or 1.0
        code = 0
        if not self.cancelled:
            self.state.update(status="running", started=time.time(), message="Running")
            self.save(force=True)
            with open(self.dir / "log.txt", "a", encoding="utf-8") as log:
                done_weight = 0.0
                for index, step in enumerate(steps):
                    if self.cancelled:
                        break
                    weight = float(step.get("weight") or 1.0)
                    self._step_base = done_weight / total_weight
                    self._step_share = weight / total_weight
                    title = step.get("title") or f"Step {index + 1}"
                    self.state.update(step=index, steps=len(steps), step_title=title)
                    self.progress(0.0, title)
                    command = [str(c) for c in step["command"]]
                    log.write(f"\n=== {title}\n$ {' '.join(shlex.quote(c) for c in command)}\n")
                    log.flush()
                    try:
                        self.child = subprocess.Popen(
                            command,
                            cwd=self.spec.get("cwd") or None,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL,
                        )
                    except OSError as exc:
                        log.write(f"Could not start: {exc}\n")
                        code = 127
                        break
                    stop = threading.Event()
                    watcher = threading.Thread(target=self._watch_history, args=(step, stop), daemon=True)
                    watcher.start()
                    try:
                        self._pump(self.child.stdout, log)
                    finally:
                        code = self.child.wait()
                        stop.set()
                    log.write(f"=== exit code {code}\n")
                    log.flush()
                    if code != 0:
                        break
                    done_weight += weight
                    self.progress(1.0)
        if lock is not None:
            lock.close()
        if self.cancelled:
            status, message = "cancelled", "Cancelled"
        elif code == 0:
            status, message = "done", "Done"
            self.state["progress"] = 1.0
        else:
            status, message = "failed", f"Failed (exit code {code}): {self.state.get('last_line', '')}"[:300]
        self.state.update(status=status, message=message, exit_code=None if self.cancelled else code, finished=time.time())
        self.save(force=True)
        return 0 if status == "done" else 1


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        print("usage: python -m genomics.visualizer.task_runner TASK_DIR", file=sys.stderr)
        return 2
    return Runner(Path(args[0])).run()


if __name__ == "__main__":
    raise SystemExit(main())
