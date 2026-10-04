"""AlphaGenome backend used by the visualizer: hosted API, a remote server, or a server on this machine.

The choice is saved to ``~/.config/genomics/visualizer_alphagenome.json``. The Perturbation Lab
calls it in-process (:meth:`AlphaGenomeBackend.create_client`); background prediction jobs get it
through ``ALPHAGENOME_ADDRESS`` / ``ALPHAGENOME_TLS_CA_CERT`` (read by
:func:`genomics.core.alphagenome_connection.create_dna_client`).

"This machine" runs ``server.py`` from an ``alphagenome_research`` checkout (JAX + model weights)
in its own Python environment, by default the conda env named ``alphagenome``. That server always
listens on 0.0.0.0:50051 and uses TLS when ``certs/server.crt`` and ``certs/server.key`` exist.
"""
from __future__ import annotations

import json
import os
import signal
import socket
import subprocess
import sys
import tempfile
import threading
from collections import deque
from pathlib import Path
from typing import Any, Dict, List, Optional

from genomics.core.alphagenome_connection import (
    ADDRESS_ENV,
    CA_CERT_ENV,
    DEFAULT_SERVER_PORT,
    api_key_available,
    parse_address,
)

MODES = ("cloud", "remote", "local")
MODE_LABELS = {"cloud": "AlphaGenome API (hosted)", "remote": "Remote server", "local": "Server on this machine"}
SERVER_DIR_ENV = "ALPHAGENOME_SERVER_DIR"
SERVER_PYTHON_ENV = "ALPHAGENOME_SERVER_PYTHON"
SERVER_CONDA_ENV = "alphagenome"
# Exit 0: alphagenome_research + jax with a CUDA plugin; 2: CPU-only jax; 1: missing packages.
# alphagenome_research needs a recent alphagenome SDK (alphagenome.io); older SDKs lack it.
_PROBE_CODE = (
    "import importlib.util, sys\n"
    "def has(m):\n"
    "    try: return importlib.util.find_spec(m) is not None\n"
    "    except ImportError: return False  # parent package missing\n"
    "if not all(has(m) for m in ('alphagenome_research', 'alphagenome.io', 'jax')): sys.exit(1)\n"
    "sys.exit(0 if any(has(m) for m in ('jax_cuda12_plugin', 'jax_cuda13_plugin', 'jax_plugins.xla_cuda12', 'jax_plugins.xla_cuda13')) else 2)\n"
)
_PROBE_CACHE: Dict[str, int] = {}


def probe_python(python: Path) -> int:
    """0 = can run the server on GPU, 2 = CPU-only jax, 1 = missing packages or not runnable."""
    key = str(python)
    if key not in _PROBE_CACHE:
        try:
            _PROBE_CACHE[key] = subprocess.run([key, "-c", _PROBE_CODE], timeout=60, capture_output=True).returncode
        except (OSError, subprocess.TimeoutExpired):
            _PROBE_CACHE[key] = 1
    return _PROBE_CACHE[key]


def port_open(port: int, host: str = "127.0.0.1") -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(0.2)
        return sock.connect_ex((host, port)) == 0


def default_settings_path() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base) / "genomics" / "visualizer_alphagenome.json"


def default_server_dir() -> Optional[Path]:
    if os.environ.get(SERVER_DIR_ENV):
        return Path(os.environ[SERVER_DIR_ENV]).expanduser()
    from genomics.workspace import repo_root

    candidate = repo_root().parent / "alphagenome_research"
    return candidate if (candidate / "server.py").exists() else None


def default_server_python() -> Path:
    """``$ALPHAGENOME_SERVER_PYTHON``; else the first conda env (``alphagenome`` first) that has
    ``alphagenome_research`` and a CUDA-enabled jax, else one with CPU-only jax, else this one."""
    if os.environ.get(SERVER_PYTHON_ENV):
        return Path(os.environ[SERVER_PYTHON_ENV]).expanduser()
    bases = []
    if os.environ.get("CONDA_EXE"):
        bases.append(Path(os.environ["CONDA_EXE"]).resolve().parents[1])
    prefix = Path(sys.prefix)
    bases.append(prefix.parent.parent if prefix.parent.name == "envs" else prefix)
    candidates: List[Path] = []
    for base in dict.fromkeys(bases):
        preferred = base / "envs" / SERVER_CONDA_ENV / "bin" / "python"
        others = sorted((base / "envs").glob("*/bin/python")) if (base / "envs").is_dir() else []
        candidates += [c for c in [preferred, *others] if c.exists() and c not in candidates]
    candidates.append(Path(sys.executable))
    codes = [probe_python(c) for c in candidates]
    for wanted in (0, 2):
        if wanted in codes:
            return candidates[codes.index(wanted)]
    return candidates[0]


class LocalAlphaGenomeServer:
    """Starts/stops ``alphagenome_research/server.py`` as a child process and reports its state."""

    def __init__(self, server_dir: Optional[Path], python: Optional[Path], log_dir: Path, port: int = DEFAULT_SERVER_PORT):
        self.server_dir = Path(server_dir).expanduser().resolve() if server_dir else None
        self.python = Path(python).expanduser() if python else None
        self.log_path = Path(log_dir) / "alphagenome_server.log"
        self.port = port
        self._process: Optional[subprocess.Popen] = None
        self._lock = threading.Lock()

    # -- configuration -------------------------------------------------------------------
    def _cert_paths(self):
        assert self.server_dir is not None
        cert = self.server_dir / os.environ.get("ALPHAGENOME_TLS_CERT", "certs/server.crt")
        key = self.server_dir / os.environ.get("ALPHAGENOME_TLS_KEY", "certs/server.key")
        ca = self.server_dir / "certs" / "ca.crt"
        return cert, key, ca

    @property
    def tls(self) -> bool:
        if self.server_dir is None:
            return False
        cert, key, _ = self._cert_paths()
        return cert.exists() and key.exists()

    @property
    def address(self) -> str:
        return f"{'grpcs' if self.tls else 'grpc'}://127.0.0.1:{self.port}"

    @property
    def ca_cert(self) -> Optional[str]:
        if not self.tls:
            return None
        ca = self._cert_paths()[2]
        return str(ca) if ca.exists() else None

    def reasons(self) -> List[str]:
        reasons = []
        if self.server_dir is None:
            reasons.append(f"alphagenome_research checkout not found (set {SERVER_DIR_ENV} or --alphagenome-server-dir)")
        elif not (self.server_dir / "server.py").exists():
            reasons.append(f"server.py not found in {self.server_dir}")
        if self.python is None or not self.python.exists():
            reasons.append(f"Python interpreter not found: {self.python} (set {SERVER_PYTHON_ENV} or --alphagenome-server-python)")
        else:
            code = probe_python(self.python)
            if code == 1:
                reasons.append(f"{self.python} lacks alphagenome_research, a recent alphagenome SDK or jax (pip install -e {self.server_dir or 'alphagenome_research'})")
            elif code == 2:
                reasons.append(f"jax in {self.python} has no CUDA support, and server.py needs a GPU (pip install -U 'jax[cuda12]', or set {SERVER_PYTHON_ENV})")
        if self.tls and self.ca_cert is None:
            reasons.append("server uses TLS but certs/ca.crt is missing, so clients cannot verify it")
        return reasons

    # -- lifecycle -----------------------------------------------------------------------
    @property
    def owned_running(self) -> bool:
        return self._process is not None and self._process.poll() is None

    def state(self) -> str:
        listening = port_open(self.port)
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
            "tls": self.tls,
            "ca_cert": self.ca_cert,
            "server_dir": str(self.server_dir) if self.server_dir else None,
            "python": str(self.python) if self.python else None,
            "reasons": self.reasons(),
            "log": str(self.log_path),
            "log_tail": self.log_tail() if state in ("starting", "ready", "exited") else [],
        }

    def start(self) -> Dict[str, Any]:
        with self._lock:
            if self.owned_running:
                return self.info()
            if port_open(self.port):
                raise RuntimeError(f"port {self.port} is already in use (an AlphaGenome server may already be running)")
            reasons = self.reasons()
            if reasons:
                raise RuntimeError("; ".join(reasons))
            assert self.server_dir is not None and self.python is not None
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
            env = os.environ.copy()
            env_bin = str(self.python.parent)
            env["PATH"] = env_bin + os.pathsep + env.get("PATH", "")
            env["CONDA_PREFIX"] = str(self.python.parent.parent)
            env["PYTHONUNBUFFERED"] = "1"
            with open(self.log_path, "w", encoding="utf-8") as log:
                self._process = subprocess.Popen(
                    [str(self.python), "server.py"],
                    cwd=str(self.server_dir),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    env=env,
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
        log_dir=log_dir or Path(tempfile.gettempdir()) / "genomics_visualizer_logs",
    )
    return AlphaGenomeBackend(
        local,
        address=getattr(args, "alphagenome_address", None),
        ca_cert=str(args.alphagenome_ca_cert) if getattr(args, "alphagenome_ca_cert", None) else None,
    )
