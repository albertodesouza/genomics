"""Startup helpers for ``genomics visualize``: binding the port, finding a visualizer already
running on it, the URL to show and open, and whether a browser can be opened here (stdlib only).
"""
from __future__ import annotations

import errno
import getpass
import http.client
import json
import os
import socket
import sys
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple, TypeVar

DEFAULT_PORT = 8780
PORT_ATTEMPTS = 20
APP_NAME = "genomics-visualizer"
LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}
WILDCARD_HOSTS = {"", "0.0.0.0", "::"}

S = TypeVar("S")


class StartupError(RuntimeError):
    """The server cannot start (reported to the user without a traceback)."""


def browser_host(host: str) -> str:
    """Host to put in URLs: a wildcard bind address is reachable here as localhost."""
    return "localhost" if host in WILDCARD_HOSTS else host


def server_url(host: str, port: int) -> str:
    name = browser_host(host)
    if ":" in name:
        name = f"[{name}]"
    return f"http://{name}:{port}/"


def probe_visualizer(host: str, port: int, timeout: float = 1.5) -> Optional[Dict[str, Any]]:
    """``/api/status`` of a visualizer listening on ``host:port``, or None if something else is there."""
    conn = http.client.HTTPConnection(browser_host(host), port, timeout=timeout)
    try:
        conn.request("GET", "/api/status")
        response = conn.getresponse()
        if response.status != 200:
            return None
        status = json.loads(response.read().decode("utf-8"))
    except (OSError, ValueError, http.client.HTTPException):
        return None
    finally:
        conn.close()
    if not isinstance(status, dict):
        return None
    # Instances started before the "app" field existed are recognised by their status fields.
    if status.get("app") == APP_NAME or {"version", "uptime", "datasets", "caches"} <= set(status):
        return status
    return None


def _requested_paths(args: Any) -> Tuple[set, set]:
    paths = {str(Path(p).expanduser().resolve()) for p in getattr(args, "dataset", None) or []}
    ids = set(getattr(args, "dataset_id", None) or [])
    return paths, ids


def can_reuse(args: Any, status: Dict[str, Any]) -> bool:
    """A running visualizer serves this invocation when it already has every requested dataset and
    no other data option (annotations, runs roots, consensus dataset, GTF) was given."""
    if any(getattr(args, name, None) for name in ("annotations", "runs_root", "consensus_dataset_dir", "gtf")):
        return False
    paths, ids = _requested_paths(args)
    datasets = [d for d in status.get("datasets") or [] if isinstance(d, dict)]
    open_paths = {str(d.get("path")) for d in datasets}
    open_ids = {str(d.get("id")) for d in datasets}
    return paths <= open_paths and ids <= open_ids


def _in_use(exc: OSError) -> bool:
    return exc.errno in (errno.EADDRINUSE, getattr(errno, "WSAEADDRINUSE", errno.EADDRINUSE))


def bind_first_free(factory: Callable[[int], S], port: int, attempts: int) -> Tuple[S, int]:
    """Bind ``port`` or, when it is taken, the next free one among ``attempts`` ports."""
    last: Optional[OSError] = None
    for candidate in range(port, port + max(1, attempts)):
        try:
            return factory(candidate), candidate
        except OSError as exc:
            if not _in_use(exc):
                raise
            last = exc
    assert last is not None
    raise last


def display_available() -> bool:
    """Whether ``webbrowser`` would open a graphical browser (not a text browser in this terminal)."""
    if sys.platform.startswith("linux") or "bsd" in sys.platform:
        return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    return True


def in_ssh_session() -> bool:
    return bool(os.environ.get("SSH_CONNECTION") or os.environ.get("SSH_TTY"))


def ssh_tunnel_hint(host: str, port: int) -> Optional[str]:
    """How to reach a loopback-bound server from the machine this SSH session comes from."""
    if not in_ssh_session() or host not in LOOPBACK_HOSTS:
        return None
    try:
        user = getpass.getuser()
    except Exception:
        user = "user"
    return f"ssh -N -L {port}:localhost:{port} {user}@{socket.gethostname()}"
