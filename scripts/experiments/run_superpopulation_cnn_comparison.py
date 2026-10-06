#!/usr/bin/env python3
"""Sequentially trains the 5 superpopulation-classification CNN2 configs that mirror the
pigmentation alignment/mask comparison (see notebooks/genotype_cnn_alignment_deeplift_summary.ipynb
and configs/predictors/genotype_based/superpopulation/*.yaml).

Meant to run unattended (launched via nohup/setsid, survives the launching shell exiting).
Designed for safety over speed:
  - Runs configs strictly one at a time (no parallel GPU/RAM contention).
  - Checks free disk space before each run; aborts the remaining queue instead of risking a
    disk-full mid-write if space gets critically low.
  - On a crashed run (non-zero exit), retries a bounded number of times, resuming from the
    latest checkpoint found under that experiment's models/ dir (epoch_N.pt > best_accuracy.pt >
    best_loss.pt) rather than restarting from epoch 1, by writing a patched copy of the config
    with checkpointing.load_checkpoint set.
  - Writes a JSON status file after every state change so an external health check (e.g. an
    hourly cron) can see progress without parsing logs.
"""
from __future__ import annotations

import copy
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path("/home/breno/I2CA/genomics")
GENOMICS_ENV_BIN = "/home/breno/miniforge3/envs/genomics/bin"
GENOMICS_BIN = f"{GENOMICS_ENV_BIN}/genomics"
GENOMICS_PYTHON = f"{GENOMICS_ENV_BIN}/python3"


def _subprocess_env() -> dict:
    # Launching genomics by absolute path puts the right python on argv[0], but does not
    # activate the conda env -- PATH stays whatever the parent shell had, so bcftools (which
    # lives in this same env's bin/, needed for alignment_mapping="bcftools_chain" consensus
    # building) is invisible to the training subprocess even though it's right next to the
    # genomics executable. Prepend the env's bin dir explicitly rather than relying on
    # inherited PATH.
    env = os.environ.copy()
    env["PATH"] = GENOMICS_ENV_BIN + os.pathsep + env.get("PATH", "")
    return env

sys.path.insert(0, str(REPO_ROOT / "src"))

RUN_DIR = REPO_ROOT / "results" / "genotype_based_predictor" / "superpopulation_comparison"
LOG_DIR = RUN_DIR / "logs"
SCRATCH_DIR = RUN_DIR / "resume_configs"
STATUS_PATH = RUN_DIR / "status.json"
LOCK_PATH = RUN_DIR / "orchestrator.pid"

CONFIG_DIR = REPO_ROOT / "configs" / "predictors" / "genotype_based" / "superpopulation"
CONFIGS = [
    ("no_alignment", CONFIG_DIR / "superpopulation_no_alignment.yaml"),
    ("dita_masks_1", CONFIG_DIR / "superpopulation_all_masks.yaml"),
    ("dita_masks_0.5", CONFIG_DIR / "superpopulation_dita_masks_0_5.yaml"),
    ("dita_masks_0.25", CONFIG_DIR / "superpopulation_dita_masks_0_25.yaml"),
    ("dita_no_masks", CONFIG_DIR / "superpopulation.yaml"),
]

MAX_RETRIES = 3
RETRY_BACKOFF_SECONDS = 60
MIN_FREE_DISK_GB = 30


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _log(msg: str) -> None:
    line = f"[{_now()}] {msg}"
    print(line, flush=True)


def _write_status(status: dict) -> None:
    status["updated_at"] = _now()
    tmp = STATUS_PATH.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(status, indent=2))
    tmp.replace(STATUS_PATH)


def _free_disk_gb(path: Path) -> float:
    usage = shutil.disk_usage(path)
    return usage.free / (1024**3)


def _experiment_dir_for(config_path: Path):
    from genomics.predictors.genotype_based.config import (
        generate_experiment_name,
        get_experiment_runs_dir,
        load_config,
    )

    config = load_config(config_path)
    return get_experiment_runs_dir(config) / generate_experiment_name(config)


def _latest_checkpoint(models_dir: Path) -> str | None:
    if not models_dir.is_dir():
        return None
    epoch_ckpts = []
    for p in models_dir.glob("epoch_*.pt"):
        m = re.match(r"epoch_(\d+)\.pt$", p.name)
        if m:
            epoch_ckpts.append((int(m.group(1)), p.name))
    if epoch_ckpts:
        epoch_ckpts.sort()
        return epoch_ckpts[-1][1]
    if (models_dir / "best_accuracy.pt").exists():
        return "best_accuracy.pt"
    if (models_dir / "best_loss.pt").exists():
        return "best_loss.pt"
    return None


def _make_resume_config(original_config_path: Path, checkpoint_name: str, attempt: int) -> Path:
    payload = yaml.safe_load(original_config_path.read_text())
    payload = copy.deepcopy(payload)
    payload.setdefault("checkpointing", {})["load_checkpoint"] = checkpoint_name
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    out_path = SCRATCH_DIR / f"{original_config_path.stem}_resume_attempt{attempt}.yaml"
    out_path.write_text(yaml.safe_dump(payload, sort_keys=False, allow_unicode=True))
    return out_path


def _run_once(config_path: Path, log_path: Path) -> int:
    with open(log_path, "a", encoding="utf-8") as log_fh:
        log_fh.write(f"\n===== {_now()} :: launching {GENOMICS_BIN} genotype train {config_path} =====\n")
        log_fh.flush()
        proc = subprocess.run(
            [GENOMICS_BIN, "genotype", "train", str(config_path)],
            cwd=str(REPO_ROOT),
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            env=_subprocess_env(),
        )
        log_fh.write(f"\n===== {_now()} :: exited with code {proc.returncode} =====\n")
    return proc.returncode


def run_config(name: str, config_path: Path, status: dict) -> str:
    log_path = LOG_DIR / f"{name}.log"
    status["configs"][name]["state"] = "running"
    status["configs"][name]["started_at"] = _now()
    _write_status(status)

    current_config = config_path
    for attempt in range(1, MAX_RETRIES + 1):
        free_gb = _free_disk_gb(REPO_ROOT)
        if free_gb < MIN_FREE_DISK_GB:
            _log(f"[{name}] ABORTING QUEUE: only {free_gb:.1f} GB free (< {MIN_FREE_DISK_GB} GB floor)")
            status["configs"][name]["state"] = "aborted_low_disk"
            status["configs"][name]["free_disk_gb_at_abort"] = free_gb
            status["aborted_low_disk"] = True
            _write_status(status)
            return "aborted_low_disk"

        _log(f"[{name}] attempt {attempt}/{MAX_RETRIES} using {current_config.name} (free disk {free_gb:.1f} GB)")
        status["configs"][name]["attempt"] = attempt
        status["configs"][name]["config_used"] = str(current_config)
        _write_status(status)

        rc = _run_once(current_config, log_path)

        if rc == 0:
            _log(f"[{name}] SUCCESS on attempt {attempt}")
            status["configs"][name]["state"] = "success"
            status["configs"][name]["finished_at"] = _now()
            _write_status(status)
            return "success"

        _log(f"[{name}] FAILED (exit {rc}) on attempt {attempt}")
        if attempt == MAX_RETRIES:
            break

        try:
            experiment_dir = _experiment_dir_for(config_path)
            checkpoint = _latest_checkpoint(experiment_dir / "models")
        except Exception as exc:  # noqa: BLE001 - best-effort resume lookup, never fatal
            _log(f"[{name}] could not resolve experiment dir for resume lookup: {exc!r}")
            checkpoint = None

        if checkpoint:
            _log(f"[{name}] resuming from checkpoint {checkpoint!r}")
            current_config = _make_resume_config(config_path, checkpoint, attempt + 1)
        else:
            _log(f"[{name}] no checkpoint found yet; retrying from scratch")
            current_config = config_path

        _log(f"[{name}] backing off {RETRY_BACKOFF_SECONDS}s before retry")
        time.sleep(RETRY_BACKOFF_SECONDS)

    _log(f"[{name}] giving up after {MAX_RETRIES} attempts")
    status["configs"][name]["state"] = "failed"
    status["configs"][name]["finished_at"] = _now()
    _write_status(status)
    return "failed"


def _previously_succeeded() -> set[str]:
    if not STATUS_PATH.exists():
        return set()
    try:
        previous = json.loads(STATUS_PATH.read_text())
    except Exception:
        return set()
    return {
        name
        for name, c in previous.get("configs", {}).items()
        if c.get("state") == "success"
    }


def main() -> int:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    if LOCK_PATH.exists():
        old_pid = int(LOCK_PATH.read_text().strip() or 0)
        if old_pid and _pid_alive(old_pid):
            _log(f"Another orchestrator instance is already running (pid {old_pid}). Exiting.")
            return 1
        _log(f"Removing stale lock file for dead pid {old_pid}")
        LOCK_PATH.unlink()

    # If a previous orchestrator run (e.g. before a bugfix) already got some configs to
    # "success", don't burn hours redoing them -- only (re)run what isn't already done.
    already_succeeded = _previously_succeeded()

    LOCK_PATH.write_text(str(os.getpid()))

    status = {
        "run_started_at": _now(),
        "pid": os.getpid(),
        "aborted_low_disk": False,
        "configs": {
            name: {"state": "pending", "config_path": str(path)} for name, path in CONFIGS
        },
    }
    _write_status(status)
    _log(f"Superpopulation CNN comparison orchestrator starting (pid={os.getpid()})")
    _log(f"Configs: {[name for name, _ in CONFIGS]}")
    if already_succeeded:
        _log(f"Skipping already-succeeded configs from a prior run: {sorted(already_succeeded)}")

    outcomes = {}
    for name, config_path in CONFIGS:
        if name in already_succeeded:
            status["configs"][name]["state"] = "success"
            status["configs"][name]["note"] = "carried forward from a prior orchestrator run"
            _write_status(status)
            outcomes[name] = "success"
            continue
        outcome = run_config(name, config_path, status)
        outcomes[name] = outcome
        if outcome == "aborted_low_disk":
            _log("Stopping remaining queue due to low disk space.")
            break

    status["run_finished_at"] = _now()
    status["all_done"] = True
    status["outcomes"] = outcomes
    _write_status(status)

    _log("Run complete. Outcomes: " + json.dumps(outcomes))
    if LOCK_PATH.exists():
        LOCK_PATH.unlink()
    return 0 if all(v == "success" for v in outcomes.values()) else 1


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


if __name__ == "__main__":
    raise SystemExit(main())
