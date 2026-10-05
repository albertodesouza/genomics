"""AlphaGenome backend used by the visualizer: hosted API, a remote server, or a server on this machine.

The choice is saved to ``~/.config/genomics/visualizer_alphagenome.json``. The Perturbation Lab
calls it in-process (:meth:`AlphaGenomeBackend.create_client`); background prediction jobs get it
through ``ALPHAGENOME_ADDRESS`` / ``ALPHAGENOME_TLS_CA_CERT`` (read by
:func:`genomics.core.alphagenome_connection.create_dna_client`).

"This machine" serves the model from an ``alphagenome_research`` checkout in its own Python
environment (JAX + weights) through :mod:`genomics.workflows.alphagenome.model_server`; discovery,
probing and installation live in :mod:`genomics.workflows.alphagenome.local_server`
(``genomics alphagenome server setup`` prepares everything). The server listens on 127.0.0.1:50051
by default (``$ALPHAGENOME_SERVER_HOST`` / ``$ALPHAGENOME_SERVER_PORT``) and uses TLS when
``certs/server.crt`` and ``certs/server.key`` exist in the checkout.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import tempfile
import threading
from collections import deque
from pathlib import Path
from typing import Any, Dict, List, Optional

from genomics.core.alphagenome_connection import (
    ADDRESS_ENV,
    CA_CERT_ENV,
    api_key_available,
    parse_address,
)
from genomics.workflows.alphagenome.local_server import (  # noqa: F401 (re-exported for callers/tests)
    MODEL_SERVER,
    PROBE_CPU,
    PROBE_GPU,
    PROBE_MISSING,
    SERVER_DIR_ENV,
    SERVER_PYTHON_ENV,
    default_host,
    default_port,
    default_server_dir,
    default_server_python,
    port_open,
    probe,
    server_command,
    server_env,
    tls_files,
)

MODES = ("cloud", "remote", "local")
MODE_LABELS = {"cloud": "AlphaGenome API (hosted)", "remote": "Remote server", "local": "Server on this machine"}
SETUP_COMMAND = "genomics alphagenome server setup"


def default_settings_path() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base) / "genomics" / "visualizer_alphagenome.json"


class LocalAlphaGenomeServer:
    """Starts/stops the model server (``model_server.py`` in the server environment) as a child process."""

    def __init__(
        self,
        server_dir: Optional[Path],
        python: Optional[Path],
        log_dir: Path,
        port: Optional[int] = None,
        host: Optional[str] = None,
        launcher: Path = MODEL_SERVER,
    ):
        self.server_dir = Path(server_dir).expanduser().resolve() if server_dir else None
        self.python = Path(python).expanduser() if python else None
        self.log_path = Path(log_dir) / "alphagenome_server.log"
        self.port = port or default_port()
        self.host = host or default_host()
        self.launcher = Path(launcher)
        self._process: Optional[subprocess.Popen] = None
        self._lock = threading.Lock()

    # -- configuration -------------------------------------------------------------------
    @property
    def connect_host(self) -> str:
        return "127.0.0.1" if self.host in ("0.0.0.0", "::", "") else self.host

    @property
    def tls(self) -> bool:
        if self.server_dir is None:
            return False
        files = tls_files(self.server_dir)
        return files["cert"].exists() and files["key"].exists()

    @property
    def address(self) -> str:
        return f"{'grpcs' if self.tls else 'grpc'}://{self.connect_host}:{self.port}"

    @property
    def ca_cert(self) -> Optional[str]:
        if not self.tls:
            return None
        ca = tls_files(self.server_dir)["ca"]  # type: ignore[arg-type]
        return str(ca) if ca.exists() else None

    def reasons(self) -> List[str]:
        reasons = []
        if self.server_dir is None:
            reasons.append(f"alphagenome_research checkout not found: run `{SETUP_COMMAND}` (or set {SERVER_DIR_ENV} / --alphagenome-server-dir)")
        elif not (self.server_dir / "server.py").exists():
            reasons.append(f"server.py not found in {self.server_dir}")
        if self.python is None or not self.python.exists():
            reasons.append(f"Python interpreter not found: {self.python}: run `{SETUP_COMMAND}` (or set {SERVER_PYTHON_ENV} / --alphagenome-server-python)")
        else:
            state = probe(self.python)
            if state.code == PROBE_MISSING:
                reasons.append(f"{self.python}: {'; '.join(state.problems)}; run `{SETUP_COMMAND}`")
            elif state.code == PROBE_CPU:
                reasons.append(f"{self.python}: {'; '.join(state.problems)}, and the server needs a GPU; run `{SETUP_COMMAND}`")
        if self.tls and self.ca_cert is None:
            reasons.append("server uses TLS but certs/ca.crt is missing, so clients cannot verify it")
        return reasons

    # -- lifecycle -----------------------------------------------------------------------
    @property
    def owned_running(self) -> bool:
        return self._process is not None and self._process.poll() is None

    def state(self) -> str:
        listening = port_open(self.port, self.connect_host)
        if self.owned_running:
            return "ready" if listening else "starting"
        if listening:
            return "external"
        if self._process is not None:
            return "exited"
        return "stopped"

    def log_tail(self, lines: int = 40) -> List[str]:
        try:
            with open(self.log_path, "r", encoding="utf-8", errors="replace") as handle:
                return [line.rstrip("\n") for line in deque(handle, maxlen=lines)]
        except OSError:
            return []

    def info(self) -> Dict[str, Any]:
        state = self.state()
        return {
            "state": state,
            "pid": self._process.pid if self.owned_running else None,
            "exit_code": self._process.returncode if state == "exited" and self._process is not None else None,
            "address": self.address,
            "host": self.host,
            "port": self.port,
            "tls": self.tls,
            "ca_cert": self.ca_cert,
            "server_dir": str(self.server_dir) if self.server_dir else None,
            "python": str(self.python) if self.python else None,
            "reasons": self.reasons(),
            "environment": probe(self.python).describe() if self.python and self.python.exists() else None,
            "setup_command": SETUP_COMMAND,
            "log": str(self.log_path),
            "log_tail": self.log_tail() if state in ("starting", "ready", "exited") else [],
        }

    def start(self) -> Dict[str, Any]:
        with self._lock:
            if self.owned_running:
                return self.info()
            if port_open(self.port, self.connect_host):
                raise RuntimeError(f"port {self.port} is already in use (an AlphaGenome server may already be running)")
            reasons = self.reasons()
            if reasons:
                raise RuntimeError("; ".join(reasons))
            assert self.server_dir is not None and self.python is not None
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.log_path, "w", encoding="utf-8") as log:
                self._process = subprocess.Popen(
                    server_command(self.python, self.server_dir, self.host, self.port, launcher=self.launcher),
                    cwd=str(self.server_dir),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    env=server_env(self.python),
                    start_new_session=True,  # Ctrl-C in the visualizer's terminal is handled by stop()
                )
        return self.info()

    def stop(self) -> Dict[str, Any]:
        with self._lock:
            process = self._process
            if process is not None and process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait(timeout=5)
                except ProcessLookupError:
                    pass
            self._process = None
        return self.info()


class AlphaGenomeBackend:
    """The AlphaGenome backend selected in the UI, persisted across visualizer restarts."""

    def __init__(
        self,
        local: LocalAlphaGenomeServer,
        settings_path: Optional[Path] = None,
        address: Optional[str] = None,
        ca_cert: Optional[str] = None,
    ):
        self.local = local
        self.settings_path = Path(settings_path) if settings_path else default_settings_path()
        self._lock = threading.Lock()
        self.settings = self._load()
        if address:  # command line wins over saved settings
            self.settings.update(mode="remote", address=address, ca_cert=ca_cert or self.settings.get("ca_cert") or "")

    def _load(self) -> Dict[str, str]:
        settings = {"mode": "cloud", "address": "", "ca_cert": ""}
        if os.environ.get(ADDRESS_ENV):
            settings.update(mode="remote", address=os.environ[ADDRESS_ENV], ca_cert=os.environ.get(CA_CERT_ENV, ""))
        try:
            saved = json.loads(self.settings_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return settings
        if saved.get("mode") in MODES:
            settings["mode"] = saved["mode"]
        for key in ("address", "ca_cert"):
            if isinstance(saved.get(key), str):
                settings[key] = saved[key]
        return settings

    def update(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        mode = payload.get("mode", self.settings["mode"])
        if mode not in MODES:
            raise ValueError(f"mode must be one of {', '.join(MODES)}")
        address = str(payload.get("address", self.settings["address"]) or "").strip()
        ca_cert = str(payload.get("ca_cert", self.settings["ca_cert"]) or "").strip()
        if mode == "remote":
            parse_address(address, ca_cert or None)  # raises ValueError on a malformed address
            if ca_cert and not Path(ca_cert).expanduser().is_file():
                raise ValueError(f"CA certificate not found: {ca_cert}")
        with self._lock:
            self.settings = {"mode": mode, "address": address, "ca_cert": ca_cert}
            self.settings_path.parent.mkdir(parents=True, exist_ok=True)
            self.settings_path.write_text(json.dumps(self.settings, indent=2) + "\n", encoding="utf-8")
        return self.describe()

    # -- what lab processes should use ------------------------------------------------------
    def _endpoint(self, settings: Optional[Dict[str, str]] = None):
        settings = settings or self.settings
        mode = settings["mode"]
        if mode == "remote":
            return settings["address"] or None, settings["ca_cert"] or None
        if mode == "local":
            return self.local.address, self.local.ca_cert
        return None, None

    @property
    def label(self) -> str:
        address, _ = self._endpoint()
        return f"{MODE_LABELS[self.settings['mode']]}" + (f" ({address})" if address else "")

    def reasons(self) -> List[str]:
        mode = self.settings["mode"]
        if mode == "cloud":
            return [] if api_key_available() else ["ALPHAGENOME_API_KEY not set (env or ~/.env); or use a self-hosted AlphaGenome server"]
        if mode == "remote":
            return [] if self.settings["address"] else ["no AlphaGenome server address set"]
        state = self.local.state()
        if state in ("ready", "external"):
            return []
        if state == "starting":
            return ["the local AlphaGenome server is still loading the model"]
        return ["the local AlphaGenome server is not running (start it above)"]

    def child_env(self, base: Optional[Dict[str, str]] = None) -> Dict[str, str]:
        env = dict(os.environ if base is None else base)
        env.pop(ADDRESS_ENV, None)
        env.pop(CA_CERT_ENV, None)
        address, ca_cert = self._endpoint()
        if address:
            env[ADDRESS_ENV] = address
            if ca_cert:
                env[CA_CERT_ENV] = str(Path(ca_cert).expanduser())
        return env

    def create_client(self, timeout: float = 60.0):
        """A ``DnaClient`` for the selected backend, for calls made inside the visualizer."""
        from genomics.core.alphagenome_connection import create_dna_client, resolve_api_key

        address, ca_cert = self._endpoint()
        if address:
            return create_dna_client(address=address, ca_cert=ca_cert, timeout=timeout)
        from alphagenome.models import dna_client

        return dna_client.create(api_key=resolve_api_key(), timeout=timeout)

    def describe(self) -> Dict[str, Any]:
        return {
            "settings": dict(self.settings),
            "modes": [{"value": m, "label": MODE_LABELS[m]} for m in MODES],
            "label": self.label,
            "reasons": self.reasons(),
            "api_key_available": api_key_available(),
            "local": self.local.info(),
            "settings_path": str(self.settings_path),
        }

    def test(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        """Check the backend described by ``payload`` (unsaved form values) or the saved one."""
        settings = {**self.settings, **{k: str(v or "").strip() for k, v in payload.items() if k in ("mode", "address", "ca_cert")}}
        if settings["mode"] not in MODES:
            return {"ok": False, "message": f"unknown mode {settings['mode']!r}"}
        if settings["mode"] == "remote" and not settings["address"]:
            return {"ok": False, "message": "enter a server address first"}
        if settings["mode"] == "local" and self.local.state() not in ("ready", "external"):
            return {"ok": False, "message": f"the local server is {self.local.state()}"}
        address, ca_cert = self._endpoint(settings)
        try:
            from genomics.core.alphagenome_connection import check_connection

            return check_connection(address, ca_cert, timeout=float(payload.get("timeout") or 10), predict=bool(payload.get("predict")))
        except ImportError as exc:
            return {"ok": False, "message": f"{exc}; install the client with: pip install -e '.[alphagenome]'"}


def create_backend(args: Any, log_dir: Optional[Path]) -> AlphaGenomeBackend:
    local = LocalAlphaGenomeServer(
        server_dir=getattr(args, "alphagenome_server_dir", None) or default_server_dir(),
        python=getattr(args, "alphagenome_server_python", None) or default_server_python(),
        port=getattr(args, "alphagenome_server_port", None),
        log_dir=log_dir or Path(tempfile.gettempdir()) / "genomics_visualizer_logs",
    )
    return AlphaGenomeBackend(
        local,
        address=getattr(args, "alphagenome_address", None),
        ca_cert=str(args.alphagenome_ca_cert) if getattr(args, "alphagenome_ca_cert", None) else None,
    )
