"""Install, check and start a self-hosted AlphaGenome server on this machine.

The model runs from an ``alphagenome_research`` checkout (https://github.com/FeLiPeOLi7/alphagenome_research,
a fork of google-deepmind/alphagenome_research with a gRPC ``server.py`` compatible with the hosted
AlphaGenome API) in **its own Python environment**: it needs JAX with CUDA, TensorFlow and the
``alphagenome>=0.7`` SDK, which conflict with the ``genomics`` environment (torch). This module runs
in the ``genomics`` environment and drives that one through subprocesses:

* :func:`setup` clones the checkout, creates the environment (conda env ``alphagenome``, else a venv
  inside the checkout), installs it with a CUDA-enabled jax and checks the weights;
* :func:`probe` inspects an interpreter (packages, versions, CUDA plugin) without importing JAX;
* :func:`server_command` is the command that serves the model (:mod:`.model_server`), used by
  ``genomics alphagenome server start`` and by the visualizer's *Start server*.

Only the standard library is imported here, so ``genomics --help`` stays light.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

REPO_URL = "https://github.com/FeLiPeOLi7/alphagenome_research"
SERVER_DIR_ENV = "ALPHAGENOME_SERVER_DIR"
SERVER_PYTHON_ENV = "ALPHAGENOME_SERVER_PYTHON"
SERVER_HOST_ENV = "ALPHAGENOME_SERVER_HOST"
SERVER_PORT_ENV = "ALPHAGENOME_SERVER_PORT"
SERVER_CONDA_ENV = "alphagenome"
DEFAULT_PORT = 50051
DEFAULT_HOST = "127.0.0.1"
DEFAULT_MODEL_VERSION = "all_folds"
SERVER_PYTHON_VERSION = "3.11"
MODEL_SERVER = Path(__file__).resolve().with_name("model_server.py")

# Probe exit codes, kept stable for callers: 0 = ready for GPU, 2 = CPU-only jax, 1 = missing packages.
PROBE_GPU, PROBE_MISSING, PROBE_CPU = 0, 1, 2
CUDA_PLUGINS = ("jax-cuda12-plugin", "jax-cuda13-plugin")
_PROBE_CODE = r"""
import importlib.util, json, sys
try:
    from importlib import metadata
except ImportError:
    metadata = None
def has(m):
    try: return importlib.util.find_spec(m) is not None
    except (ImportError, ValueError): return False
def version(dist):
    try: return metadata.version(dist)
    except Exception: return None
info = {"python": "%d.%d.%d" % sys.version_info[:3], "prefix": sys.prefix,
        "modules": {m: has(m) for m in ("alphagenome", "alphagenome.io", "alphagenome_research", "jax", "grpc", "huggingface_hub")},
        "versions": {d: version(d) for d in ("alphagenome", "alphagenome_research", "jax", "jaxlib", "jax-cuda12-plugin", "jax-cuda13-plugin", "grpcio")}}
print(json.dumps(info))
"""
_PROBE_CACHE: Dict[str, "Probe"] = {}


@dataclass
class Probe:
    """What an interpreter can do for the server. ``code`` is one of the ``PROBE_*`` constants."""

    python: str
    code: int
    problems: List[str] = field(default_factory=list)
    info: Dict = field(default_factory=dict)

    @property
    def ready(self) -> bool:
        return self.code == PROBE_GPU

    def describe(self) -> str:
        versions = self.info.get("versions", {})
        parts = [f"Python {self.info.get('python', '?')}"]
        for dist in ("alphagenome", "alphagenome_research", "jax", *CUDA_PLUGINS):
            if versions.get(dist):
                parts.append(f"{dist} {versions[dist]}")
        return ", ".join(parts)


def _analyse(python: str, info: Dict) -> Probe:
    modules, versions = info.get("modules", {}), info.get("versions", {})
    problems: List[str] = []
    if not modules.get("alphagenome_research"):
        problems.append("alphagenome_research is not installed")
    if not modules.get("alphagenome.io"):
        found = versions.get("alphagenome")
        problems.append(f"the alphagenome SDK is {'version ' + found if found else 'missing'}; alphagenome_research needs alphagenome>=0.7")
    if not modules.get("jax"):
        problems.append("jax is not installed")
    if problems:
        return Probe(python, PROBE_MISSING, problems, info)
    jaxlib = versions.get("jaxlib")
    plugins = {d: versions[d] for d in CUDA_PLUGINS if versions.get(d)}
    if not plugins:
        return Probe(python, PROBE_CPU, [f"jax {versions.get('jax')} has no CUDA plugin (CPU only)"], info)
    if jaxlib and not any(v == jaxlib for v in plugins.values()):
        found = ", ".join(f"{d} {v}" for d, v in plugins.items())
        return Probe(python, PROBE_CPU, [f"{found} does not match jaxlib {jaxlib}, so jax falls back to CPU"], info)
    return Probe(python, PROBE_GPU, [], info)


def probe(python: Path, refresh: bool = False) -> Probe:
    """Inspect ``python`` (cached per path; ``refresh`` re-runs it, e.g. after an install)."""
    key = str(python)
    if refresh or key not in _PROBE_CACHE:
        try:
            done = subprocess.run([key, "-c", _PROBE_CODE], timeout=60, capture_output=True, text=True)
            info = json.loads(done.stdout.strip().splitlines()[-1]) if done.returncode == 0 and done.stdout.strip() else None
        except (OSError, subprocess.TimeoutExpired, ValueError):
            info = None
        if info is None:
            _PROBE_CACHE[key] = Probe(key, PROBE_MISSING, [f"{key} is not a runnable Python interpreter"])
        else:
            _PROBE_CACHE[key] = _analyse(key, info)
    return _PROBE_CACHE[key]


def probe_python(python: Path) -> int:
    """``PROBE_GPU`` (0), ``PROBE_CPU`` (2) or ``PROBE_MISSING`` (1)."""
    return probe(python).code


def port_open(port: int, host: str = "127.0.0.1") -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(0.2)
        return sock.connect_ex((host, port)) == 0


# -- locations ------------------------------------------------------------------------------------

def _repo_checkout_parent() -> Optional[Path]:
    """The directory holding this repository, when ``genomics`` runs from a source checkout."""
    from genomics.workspace import repo_root

    root = repo_root()
    return root.parent if (root / "pyproject.toml").is_file() else None


def data_home() -> Path:
    return Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share") / "genomics"


def candidate_server_dirs() -> List[Path]:
    dirs = []
    if os.environ.get(SERVER_DIR_ENV):
        dirs.append(Path(os.environ[SERVER_DIR_ENV]).expanduser())
    parent = _repo_checkout_parent()
    if parent is not None:
        dirs.append(parent / "alphagenome_research")
    dirs.append(data_home() / "alphagenome_research")
    return list(dict.fromkeys(dirs))


def default_server_dir() -> Optional[Path]:
    """``$ALPHAGENOME_SERVER_DIR``, ``../alphagenome_research`` next to this repository, or the copy
    ``setup`` makes under ``~/.local/share/genomics``: the first with ``server.py``."""
    for candidate in candidate_server_dirs():
        if (candidate / "server.py").is_file():
            return candidate
    return None


def default_install_dir() -> Path:
    """Where ``setup`` clones the checkout when none exists."""
    return candidate_server_dirs()[0]


def conda_executable() -> Optional[str]:
    for name in ("CONDA_EXE", "MAMBA_EXE"):
        if os.environ.get(name) and Path(os.environ[name]).exists():
            return os.environ[name]
    return shutil.which("conda") or shutil.which("mamba")


def conda_base() -> Optional[Path]:
    exe = conda_executable()
    if exe:
        return Path(exe).resolve().parents[1]
    prefix = Path(sys.prefix)
    if prefix.parent.name == "envs":
        return prefix.parent.parent
    return None


def candidate_pythons(server_dir: Optional[Path] = None) -> List[Path]:
    """Interpreters that may run the server: the checkout's ``.venv``, the conda env ``alphagenome``,
    every other conda env, then this interpreter."""
    candidates: List[Path] = []
    if server_dir is not None:
        candidates.append(Path(server_dir) / ".venv" / "bin" / "python")
    bases = [b for b in (conda_base(),) if b is not None]
    prefix = Path(sys.prefix)
    bases.append(prefix.parent.parent if prefix.parent.name == "envs" else prefix)
    for base in dict.fromkeys(bases):
        candidates.append(base / "envs" / SERVER_CONDA_ENV / "bin" / "python")
        if (base / "envs").is_dir():
            candidates.extend(sorted((base / "envs").glob("*/bin/python")))
    found = [c for c in dict.fromkeys(candidates) if c.exists()]
    return found + [Path(sys.executable)]


def default_server_python(server_dir: Optional[Path] = None) -> Path:
    """``$ALPHAGENOME_SERVER_PYTHON``; else the first candidate that can serve on GPU, else one with
    CPU-only jax (reported as such), else the first candidate."""
    if os.environ.get(SERVER_PYTHON_ENV):
        return Path(os.environ[SERVER_PYTHON_ENV]).expanduser()
    candidates = candidate_pythons(server_dir if server_dir is not None else default_server_dir())
    codes = [probe_python(c) for c in candidates]
    for wanted in (PROBE_GPU, PROBE_CPU):
        if wanted in codes:
            return candidates[codes.index(wanted)]
    return candidates[0]


def default_host() -> str:
    return os.environ.get(SERVER_HOST_ENV) or DEFAULT_HOST


def default_port() -> int:
    try:
        return int(os.environ.get(SERVER_PORT_ENV) or DEFAULT_PORT)
    except ValueError:
        return DEFAULT_PORT


def server_command(
    python: Path,
    server_dir: Path,
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    extra: Sequence[str] = (),
    launcher: Path = MODEL_SERVER,
) -> List[str]:
    return [str(python), str(launcher), "--server-dir", str(server_dir), "--host", host, "--port", str(port), *extra]


def server_env(python: Path, base: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Environment for the server process: its interpreter's ``bin`` first on PATH, as ``conda activate`` would."""
    env = dict(os.environ if base is None else base)
    env_bin = Path(python).parent
    env["PATH"] = str(env_bin) + os.pathsep + env.get("PATH", "")
    if (env_bin.parent / "conda-meta").is_dir():
        env["CONDA_PREFIX"] = str(env_bin.parent)
    else:
        env.pop("CONDA_PREFIX", None)
    env["PYTHONUNBUFFERED"] = "1"
    env.pop("PYTHONPATH", None)  # the genomics env's paths must not leak into the server env
    return env


def tls_files(server_dir: Path) -> Dict[str, Path]:
    return {
        "cert": server_dir / os.environ.get("ALPHAGENOME_TLS_CERT", "certs/server.crt"),
        "key": server_dir / os.environ.get("ALPHAGENOME_TLS_KEY", "certs/server.key"),
        "ca": server_dir / "certs" / "ca.crt",
    }


_WEIGHTS_CODE = r"""
import sys
try:
    import huggingface_hub
    print(huggingface_hub.snapshot_download(repo_id=sys.argv[1], local_files_only=True))
except Exception as exc:
    print("MISSING " + type(exc).__name__)
    sys.exit(3)
"""


def hf_repo_id(model_version: str = DEFAULT_MODEL_VERSION) -> str:
    return f"google/alphagenome-{model_version.replace('_', '-').lower()}"


def cached_weights(python: Path, model_version: str = DEFAULT_MODEL_VERSION) -> Optional[str]:
    """Path of the cached Hugging Face checkpoint as seen by ``python``, or None when not downloaded."""
    try:
        done = subprocess.run([str(python), "-c", _WEIGHTS_CODE, hf_repo_id(model_version)], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    lines = done.stdout.strip().splitlines()
    return lines[-1] if done.returncode == 0 and lines else None


# -- setup ----------------------------------------------------------------------------------------

def _run(cmd: Sequence[str], dry_run: bool, env: Optional[Dict[str, str]] = None, cwd: Optional[Path] = None) -> None:
    print("$ " + " ".join(cmd), flush=True)
    if dry_run:
        return
    done = subprocess.run(list(cmd), env=env, cwd=str(cwd) if cwd else None)
    if done.returncode != 0:
        raise SystemExit(f"command failed with exit code {done.returncode}: {' '.join(cmd)}")


def _jax_extra(jax: str) -> Optional[str]:
    return {"cuda12": "cuda12", "cuda13": "cuda13", "cuda": "cuda12", "cpu": None}[jax]


def setup(
    server_dir: Optional[Path] = None,
    conda_env: Optional[str] = None,
    python: Optional[Path] = None,
    jax: str = "cuda12",
    download_weights: bool = False,
    update: bool = False,
    reinstall: bool = False,
    dry_run: bool = False,
) -> int:
    """Clone, create the environment, install with a CUDA-enabled jax, and check the weights."""
    server_dir = Path(server_dir).expanduser().resolve() if server_dir else (default_server_dir() or default_install_dir())
    print(f"[1/4] alphagenome_research checkout: {server_dir}")
    if (server_dir / "server.py").is_file():
        if update and (server_dir / ".git").exists():
            _run(["git", "-C", str(server_dir), "pull", "--ff-only"], dry_run)
        else:
            print("      found (use --update to git pull)")
    elif server_dir.exists() and any(server_dir.iterdir()):
        raise SystemExit(f"{server_dir} exists but has no server.py; pass another --dir")
    else:
        if not shutil.which("git"):
            raise SystemExit(f"git is required to clone {REPO_URL}")
        server_dir.parent.mkdir(parents=True, exist_ok=True)
        _run(["git", "clone", REPO_URL, str(server_dir)], dry_run)

    print("[2/4] Python environment for the server")
    if python is not None:
        python = Path(python).expanduser()
        print(f"      using {python}")
    else:
        base = conda_base()
        conda = conda_executable()
        if conda and base is not None and not (conda_env is None and (server_dir / ".venv").exists()):
            name = conda_env or SERVER_CONDA_ENV
            python = base / "envs" / name / "bin" / "python"
            if python.exists():
                print(f"      conda env {name!r} exists ({python})")
            else:
                _run([conda, "create", "-y", "-n", name, f"python={SERVER_PYTHON_VERSION}", "pip"], dry_run)
        else:
            python = server_dir / ".venv" / "bin" / "python"
            if python.exists():
                print(f"      venv exists ({python})")
            else:
                if sys.version_info < (3, 10):
                    raise SystemExit("no conda found and this Python is older than 3.10; install Miniforge or pass --python")
                _run([sys.executable, "-m", "venv", str(server_dir / ".venv")], dry_run)

    print("[3/4] packages")
    env = server_env(python)
    state = probe(python, refresh=True) if python.exists() else Probe(str(python), PROBE_MISSING, ["not created yet"])
    if state.code == PROBE_MISSING or reinstall:
        _run([str(python), "-m", "pip", "install", "--upgrade", "pip"], dry_run, env)
        _run([str(python), "-m", "pip", "install", "-e", str(server_dir)], dry_run, env)
        state = probe(python, refresh=True) if not dry_run else state
    extra = _jax_extra(jax)
    if extra and (state.code != PROBE_GPU or reinstall):
        jax_version = state.info.get("versions", {}).get("jax")
        spec = f"jax[{extra}]=={jax_version}" if jax_version else f"jax[{extra}]"
        _run([str(python), "-m", "pip", "install", spec], dry_run, env)
    if dry_run:
        print("[4/4] weights: (dry run)")
        return 0
    state = probe(python, refresh=True)
    print(f"      {state.describe()}")
    for problem in state.problems:
        print(f"      problem: {problem}")

    print(f"[4/4] model weights ({hf_repo_id()})")
    weights = cached_weights(python)
    if weights:
        print(f"      cached: {weights}")
    elif download_weights:
        code = f"import huggingface_hub; print(huggingface_hub.snapshot_download(repo_id={hf_repo_id()!r}))"
        _run([str(python), "-c", code], False, env)
        weights = cached_weights(python)
    else:
        print(weights_help())
    ok = state.ready or (jax == "cpu" and state.code == PROBE_CPU)
    print()
    if ok:
        print("Ready. Start the server with:\n  genomics alphagenome server start\n"
              "or open the visualizer's AlphaGenome page and choose 'Server on this machine'.")
        if os.environ.get(SERVER_PYTHON_ENV) is None and default_server_python(server_dir) != python:
            print(f"(Another environment is picked by default; set {SERVER_PYTHON_ENV}={python} to use this one.)")
    return 0 if ok and weights else 1


def weights_help() -> str:
    repo = hf_repo_id()
    return (
        f"      not downloaded. The weights (~700 MB) are gated on Hugging Face:\n"
        f"        1. accept the terms at https://huggingface.co/{repo}\n"
        f"        2. log in once in the server environment: hf auth login   (or export HF_TOKEN=...)\n"
        f"        3. re-run with --download-weights (or let the first server start download them)"
    )


# -- check / start --------------------------------------------------------------------------------

def check(server_dir: Optional[Path], python: Optional[Path], address: Optional[str], predict: bool, port: int) -> int:
    server_dir = Path(server_dir).expanduser() if server_dir else default_server_dir()
    python = Path(python).expanduser() if python else default_server_python(server_dir)
    ok = True
    print(f"checkout : {server_dir or 'not found'}" + ("" if server_dir else f" (run: genomics alphagenome server setup; looked in {', '.join(map(str, candidate_server_dirs()))})"))
    ok &= server_dir is not None
    state = probe(python)
    print(f"python   : {python}\n           {state.describe() if state.info else ''}")
    for problem in state.problems:
        print(f"  problem: {problem}")
    ok &= state.ready
    weights = cached_weights(python) if state.info else None
    print(f"weights  : {weights or 'not downloaded'}")
    if not weights:
        print(weights_help())
    gpu = shutil.which("nvidia-smi")
    print(f"gpu      : {'nvidia-smi found' if gpu else 'nvidia-smi not found (no NVIDIA driver?)'}")
    if address is None and port_open(port):
        address = f"127.0.0.1:{port}"
        print(f"server   : something listens on port {port}")
    if address:
        try:
            from genomics.core.alphagenome_connection import check_connection
        except ImportError as exc:  # pragma: no cover - grpc missing
            print(f"connect  : {exc}")
            return 1
        result = check_connection(address, timeout=15 if not predict else 600, predict=predict)
        print(f"connect  : {'OK' if result['ok'] else 'FAILED'}: {result['message']}")
        ok &= bool(result["ok"])
    return 0 if ok else 1


def start(server_dir: Optional[Path], python: Optional[Path], host: str, port: int, extra: Sequence[str]) -> int:
    server_dir = Path(server_dir).expanduser() if server_dir else default_server_dir()
    if server_dir is None:
        raise SystemExit("alphagenome_research checkout not found; run: genomics alphagenome server setup")
    python = Path(python).expanduser() if python else default_server_python(server_dir)
    state = probe(python)
    if state.code == PROBE_MISSING:
        raise SystemExit(f"{python} cannot run the server: {'; '.join(state.problems)}. Run: genomics alphagenome server setup")
    if state.code == PROBE_CPU and "--allow-cpu" not in extra:
        raise SystemExit(f"{python}: {'; '.join(state.problems)}. Run: genomics alphagenome server setup (or pass --allow-cpu)")
    if port_open(port):
        raise SystemExit(f"port {port} is already in use (is a server already running? check with: genomics alphagenome server check)")
    cmd = server_command(python, server_dir, host, port, extra)
    print("$ " + " ".join(cmd), flush=True)
    try:
        return subprocess.call(cmd, env=server_env(python), cwd=str(server_dir))
    except KeyboardInterrupt:
        return 130


DESCRIPTION = "Run AlphaGenome on this machine's GPU behind the hosted API's gRPC interface"


def build_parser() -> argparse.ArgumentParser:
    return add_arguments(argparse.ArgumentParser(prog="genomics alphagenome server", description=DESCRIPTION))


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    sub = parser.add_subparsers(dest="action", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--dir", type=Path, default=None, help=f"alphagenome_research checkout (default: ${SERVER_DIR_ENV}, ../alphagenome_research, or ~/.local/share/genomics/alphagenome_research)")
    common.add_argument("--python", type=Path, default=None, help=f"Interpreter of the server environment (default: ${SERVER_PYTHON_ENV}, else the conda env '{SERVER_CONDA_ENV}' or another env that has a CUDA jax)")

    s = sub.add_parser("setup", parents=[common], help="Clone alphagenome_research, create its environment with CUDA jax, check the weights")
    s.add_argument("--conda-env", default=None, help=f"Conda environment to create/use (default: {SERVER_CONDA_ENV}); without conda a venv is made in the checkout")
    s.add_argument("--jax", choices=["cuda12", "cuda13", "cpu"], default="cuda12", help="jax build to install (default: cuda12, works with NVIDIA drivers >= 525)")
    s.add_argument("--download-weights", action="store_true", help="Download the model weights now (needs a Hugging Face login with the terms accepted)")
    s.add_argument("--update", action="store_true", help="git pull the checkout")
    s.add_argument("--reinstall", action="store_true", help="Reinstall the packages even when they look fine")
    s.add_argument("--dry-run", action="store_true", help="Print the commands without running them")

    st = sub.add_parser("start", parents=[common], help="Serve the model in the foreground (Ctrl+C stops it)")
    st.add_argument("--host", default=default_host(), help=f"Bind address (default: {DEFAULT_HOST}, this machine only; 0.0.0.0 for the network, unauthenticated)")
    st.add_argument("--port", type=int, default=default_port(), help=f"Port (default: {DEFAULT_PORT})")
    st.add_argument("--model-version", default=DEFAULT_MODEL_VERSION, help="all_folds (default) or fold_0 ... fold_3")
    st.add_argument("--checkpoint", default=None, metavar="DIR", help="Local checkpoint directory instead of Hugging Face")
    st.add_argument("--plaintext", action="store_true", help="Ignore certs/server.{crt,key} in the checkout")
    st.add_argument("--allow-cpu", action="store_true", help="Run on CPU-only jax (very slow)")

    c = sub.add_parser("check", parents=[common], help="Report the checkout, environment, weights and (when running) the server")
    c.add_argument("--address", default=None, help="Server to test (default: 127.0.0.1:<port> when something listens there)")
    c.add_argument("--port", type=int, default=default_port())
    c.add_argument("--predict", action="store_true", help="Also run a 16 kb test prediction")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


def run(args: argparse.Namespace) -> int:
    if args.action == "setup":
        return setup(args.dir, args.conda_env, args.python, args.jax, args.download_weights, args.update, args.reinstall, args.dry_run)
    if args.action == "start":
        extra: List[str] = ["--model-version", args.model_version]
        if args.checkpoint:
            extra += ["--checkpoint", args.checkpoint]
        if args.plaintext:
            extra.append("--plaintext")
        if args.allow_cpu:
            extra.append("--allow-cpu")
        return start(args.dir, args.python, args.host, args.port, extra)
    return check(args.dir, args.python, args.address, args.predict, args.port)


if __name__ == "__main__":
    sys.exit(main())
