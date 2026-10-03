"""Optional heavyweight tools launched on demand from the visualizer.

The Pigmentation Sequence Lab needs a trained CNN checkpoint (torch) and an AlphaGenome backend
(API key or a self-hosted server, see :mod:`genomics.visualizer.alphagenome`), so instead of
starting it with the visualizer it is spawned as a separate process only when the user asks for
it from the Labs page.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional

from genomics.visualizer.alphagenome import AlphaGenomeBackend, port_open

PIGMENTATION_LAB_MODULE = "genomics.predictors.genotype_based.apps.pigmentation_sequence_lab"


class LabsService:
    def __init__(
        self,
        dataset_dir: Optional[Path],
        consensus_dir: Optional[Path] = None,
        pigmentation_config: Optional[Path] = None,
        port: int = 8781,
        log_dir: Optional[Path] = None,
        alphagenome: Optional[AlphaGenomeBackend] = None,
    ):
        self.dataset_dir = Path(dataset_dir) if dataset_dir else None
        self.consensus_dir = Path(consensus_dir) if consensus_dir else self.dataset_dir
        self.pigmentation_config = Path(pigmentation_config) if pigmentation_config else None
        self.port = port
        self.log_dir = Path(log_dir) if log_dir else Path(tempfile.gettempdir()) / "genomics_visualizer_logs"
        self._process: Optional[subprocess.Popen] = None
        self._lock = threading.Lock()
        self._reasons: Optional[List[str]] = None
        self.alphagenome = alphagenome
        self._backend_label: Optional[str] = None

    def _config(self) -> Path:
        if self.pigmentation_config:
            return self.pigmentation_config
        from genomics.workspace import repo_root

        return repo_root() / "configs" / "predictors" / "genotype_based" / "pigmentation" / "pigmentation_binary.yaml"

    def _pigmentation_reasons(self) -> List[str]:
        return self._static_reasons() + self._backend_reasons()

    def _backend_reasons(self) -> List[str]:
        if self.alphagenome is not None:
            return self.alphagenome.reasons()
        from genomics.core.alphagenome_connection import backend_configured

        return [] if backend_configured() else ["ALPHAGENOME_API_KEY not set (env or ~/.env)"]

    def _static_reasons(self) -> List[str]:
        if self._reasons is not None:
            return self._reasons
        reasons = []
        if self.dataset_dir is None:
            reasons.append("no dataset loaded")
        config = self._config()
        if not config.exists():
            reasons.append(f"config not found: {config}")
        else:
            try:
                from genomics.predictors.genotype_based.apps.genomics_workbench import _pigmentation_checkpoint_exists

                if not _pigmentation_checkpoint_exists(config):
                    reasons.append("trained 'best_accuracy' checkpoint not found for the config")
            except Exception as exc:
                reasons.append(f"could not inspect config: {exc}")
        self._reasons = reasons
        return reasons

    @property
    def running(self) -> bool:
        return self._process is not None and self._process.poll() is None

    def list(self) -> List[Dict[str, Any]]:
        reasons = self._pigmentation_reasons()
        running = self.running
        current = self.alphagenome.label if self.alphagenome is not None else None
        return [
            {
                "key": "pigmentation",
                "title": "Pigmentation Sequence Lab",
                "description": "Edit a haplotype in silico (overwrite / scramble a region), re-predict with AlphaGenome and "
                "see the effect on the trained CNN2 pigmentation classifier. Requires a checkpoint and an AlphaGenome backend.",
                "available": not reasons,
                "reasons": [] if running else reasons,
                "running": running,
                "backend": self._backend_label if running else current,
                "backend_changed": bool(running and current is not None and self._backend_label != current),
                "url": f"http://127.0.0.1:{self.port}/",
                "log": str(self.log_dir / "pigmentation_lab.log"),
            }
        ]

    def start(self, key: str) -> Dict[str, Any]:
        if key != "pigmentation":
            raise KeyError(f"Unknown lab: {key}")
        with self._lock:
            if self.running:
                return self.list()[0]
            reasons = self._pigmentation_reasons()
            if reasons:
                raise RuntimeError("; ".join(reasons))
            if port_open(self.port):
                raise RuntimeError(f"Port {self.port} is already in use")
            self.log_dir.mkdir(parents=True, exist_ok=True)
            log = open(self.log_dir / "pigmentation_lab.log", "w", encoding="utf-8")
            cmd = [
                sys.executable, "-m", PIGMENTATION_LAB_MODULE, str(Path(self.dataset_dir).resolve()),
                "--host", "127.0.0.1", "--port", str(self.port),
                "--consensus-dataset-dir", str(Path(self.consensus_dir).resolve()),
                "--config", str(self._config().resolve()),
            ]
            env = self.alphagenome.child_env() if self.alphagenome is not None else None
            # Config paths (results_dir, caches) are relative to the repository root, so run the
            # lab from there regardless of where the visualizer was started.
            from genomics.workspace import repo_root

            self._process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=str(repo_root()))
            self._backend_label = self.alphagenome.label if self.alphagenome is not None else None
        return self.list()[0]

    def stop(self, key: str) -> Dict[str, Any]:
        if key != "pigmentation":
            raise KeyError(f"Unknown lab: {key}")
        with self._lock:
            self._stop_process()
        return self.list()[0]

    def _stop_process(self) -> None:
        if self.running and self._process is not None:
            self._process.terminate()
            try:
                self._process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self._process.kill()
        self._process = None

    def stop_all(self) -> None:
        self._stop_process()
        if self.alphagenome is not None:
            self.alphagenome.local.stop()
