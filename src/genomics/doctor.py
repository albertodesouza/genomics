"""``genomics doctor``: what this machine can do, feature by feature, and how to enable the rest.

Only the standard library is used (plus light ``genomics`` modules), and optional packages are
detected without importing them, so this runs quickly even on a half-installed environment. The
same report feeds the visualizer's System page (``/api/system``).
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from importlib import util as importlib_util
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from importlib import metadata as importlib_metadata
except ImportError:  # pragma: no cover - Python 3.7
    importlib_metadata = None  # type: ignore[assignment]

OK, WARN, MISSING, INFO = "ok", "warn", "missing", "info"
SYMBOLS = {OK: "✓", WARN: "!", MISSING: "✗", INFO: "·"}
GB = 1024 ** 3


@dataclass
class Check:
    name: str
    status: str
    detail: str = ""
    fix: str = ""


@dataclass
class Feature:
    key: str
    title: str
    purpose: str
    checks: List[Check] = field(default_factory=list)
    required: bool = False

    @property
    def status(self) -> str:
        statuses = {c.status for c in self.checks}
        if MISSING in statuses:
            return MISSING
        if WARN in statuses:
            return WARN
        return OK


# -- probes ---------------------------------------------------------------------------------------

def dist_version(dist: str) -> Optional[str]:
    if importlib_metadata is None:
        return None
    try:
        return importlib_metadata.version(dist)
    except Exception:
        return None


def has_module(name: str) -> bool:
    try:
        return importlib_util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def package_check(module: str, dist: Optional[str] = None, fix: str = "", why: str = "", status_if_missing: str = MISSING) -> Check:
    dist = dist or module
    if has_module(module):
        version = dist_version(dist)
        return Check(dist, OK, version or "installed")
    return Check(dist, status_if_missing, f"not installed{f' ({why})' if why else ''}", fix)


def find_tool(tool: str) -> Optional[str]:
    """``tool`` on PATH, or in this interpreter's ``bin`` (an environment that is not activated),
    where the visualizer's background jobs also look."""
    return shutil.which(tool) or shutil.which(tool, path=str(Path(sys.executable).parent))


def tool_version(tool: str) -> Optional[str]:
    path = find_tool(tool)
    if not path:
        return None
    try:
        done = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=10)
        first = (done.stdout or done.stderr).strip().splitlines()
        return first[0] if first else path
    except (OSError, subprocess.TimeoutExpired):
        return path


def gpus() -> List[Dict[str, str]]:
    if not shutil.which("nvidia-smi"):
        return []
    try:
        done = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"],
                              capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return []
    rows = []
    for line in done.stdout.strip().splitlines():
        parts = [p.strip() for p in line.split(",")]
        if parts and parts[0]:
            memory = parts[1] if len(parts) > 1 and not parts[1].startswith("[") else ""  # "[N/A]" on unified-memory GPUs (GB10)
            rows.append({"name": parts[0], "memory": memory, "driver": parts[2] if len(parts) > 2 else ""})
    return rows


def memory_gb() -> Optional[float]:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) / 1024 / 1024
    except (OSError, ValueError, IndexError):
        pass
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / GB
    except (ValueError, OSError, AttributeError):
        return None


def free_space(path: Path) -> Optional[Dict[str, Any]]:
    """Free space on the filesystem holding ``path`` (or its nearest existing parent)."""
    probe = Path(path)
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    try:
        usage = shutil.disk_usage(str(probe))
    except OSError:
        return None
    return {"path": str(path), "exists": Path(path).exists(), "free_gb": usage.free / GB, "total_gb": usage.total / GB}


# -- features -------------------------------------------------------------------------------------

def core_feature() -> Feature:
    feature = Feature("core", "Visualizer and CLI", "genomics visualize: browse datasets, tracks, sequences, experiments and jobs", required=True)
    version = "%d.%d.%d" % sys.version_info[:3]
    feature.checks.append(Check("Python", OK if sys.version_info >= (3, 8) else MISSING, f"{version} ({sys.executable})", "Python 3.10 or newer is recommended"))
    feature.checks.append(package_check("numpy", fix="pip install -e ."))
    feature.checks.append(package_check("yaml", "PyYAML", fix="pip install -e ."))
    return feature


def import_feature() -> Feature:
    feature = Feature("import", "Gene models and dataset import", "gene models on Tracks/Sequence, gene search, importing a dataset from a VCF")
    feature.checks.append(package_check("pandas", fix="pip install -e '.[visualizer]'", why="gene models, sample tables"))
    feature.checks.append(package_check("pyarrow", fix="pip install -e '.[visualizer]'", why="reads gtf_cache.feather"))
    for tool in ("bcftools", "samtools"):
        found = tool_version(tool)
        feature.checks.append(Check(tool, OK if found else MISSING, found or "not on PATH (VCF import, new prediction windows, training-axis alignment)",
                                    "conda install -c conda-forge -c bioconda bcftools samtools   (or apt install bcftools samtools)"))
    return feature


def alphagenome_feature(include_local: bool = True, backend: Optional[Dict[str, Any]] = None) -> Feature:
    """``backend`` (``{"label", "reasons"}``) is the backend a running visualizer has selected; it
    replaces the environment checks."""
    feature = Feature("alphagenome", "AlphaGenome predictions", "prediction jobs, reference-genome tracks, Perturbation Lab, track catalog")
    if sys.version_info < (3, 10):
        feature.checks.append(Check("alphagenome", MISSING, "the AlphaGenome client needs Python >= 3.10", "use a Python 3.10+ environment"))
    else:
        feature.checks.append(package_check("alphagenome", fix="pip install -e '.[visualizer]'", why="AlphaGenome client"))
    if backend is not None:
        reasons = backend.get("reasons") or []
        feature.checks.append(Check("backend", MISSING if reasons else OK, f"{backend.get('label')}{': ' + '; '.join(reasons) if reasons else ''}",
                                    "choose or start one on the AlphaGenome page" if reasons else ""))
        return feature
    from genomics.core.alphagenome_connection import ADDRESS_ENV, api_key_available

    address = os.environ.get(ADDRESS_ENV)
    key = api_key_available()
    backends = []
    if key:
        backends.append(Check("hosted API key", OK, "ALPHAGENOME_API_KEY found (environment or ~/.env)"))
    if address:
        backends.append(Check("self-hosted server", OK, f"ALPHAGENOME_ADDRESS={address}"))
    local_ok = False
    if include_local:
        local = local_server_check()
        local_ok = local.status == OK
        backends.append(local)
    if not key and not address and not local_ok:
        backends.insert(0, Check("backend", MISSING, "no hosted API key, no remote server and no local server ready",
                                 "get a key at https://deepmind.google.com/science/alphagenome and put ALPHAGENOME_API_KEY=... in ~/.env, "
                                 "or run a local server: genomics alphagenome server setup"))
    elif not key and include_local:
        backends.append(Check("hosted API key", INFO, "not set (not needed with a self-hosted or local server)"))
    feature.checks.extend(backends)
    return feature


def local_server_check() -> Check:
    """The local AlphaGenome server (GPU) as one check: checkout, environment, weights."""
    from genomics.workflows.alphagenome import local_server

    server_dir = local_server.default_server_dir()
    if server_dir is None:
        return Check("local server", INFO, "not set up (optional: runs AlphaGenome on this machine's NVIDIA GPU)", "genomics alphagenome server setup")
    python = local_server.default_server_python(server_dir)
    state = local_server.probe(python)
    if not state.ready:
        return Check("local server", WARN, f"checkout {server_dir}; {python}: {'; '.join(state.problems)}", "genomics alphagenome server setup")
    weights = local_server.cached_weights(python)
    if not weights:
        return Check("local server", WARN, f"environment ready ({python}) but the model weights are not downloaded",
                     "genomics alphagenome server setup --download-weights (after accepting the terms on Hugging Face)")
    running = local_server.port_open(local_server.default_port())
    return Check("local server", OK, f"ready ({state.describe()}){'; listening on port %d' % local_server.default_port() if running else ''}",
                 "" if running else "start it on the visualizer's AlphaGenome page, or: genomics alphagenome server start")


def training_feature() -> Feature:
    feature = Feature("training", "Training, evaluation and Perturbation Lab scoring", "train/evaluate genotype-based models; re-score edited haplotypes")
    fix = "pip install -e '.[genotype]'"
    torch_version = dist_version("torch") if has_module("torch") else None
    if torch_version:
        cuda = "+cu" in torch_version or "+rocm" in torch_version
        feature.checks.append(Check("torch", OK if cuda or not gpus() else WARN, torch_version + ("" if cuda else " (no CUDA tag: CPU build, or a build whose CUDA support is unknown)"),
                                    "" if cuda else "install a CUDA build of PyTorch (https://pytorch.org/get-started/locally/)"))
    else:
        feature.checks.append(Check("torch", MISSING, "not installed", fix))
    for module, dist in (("sklearn", "scikit-learn"), ("scipy", "scipy"), ("pydantic", "pydantic"), ("pandas", "pandas"), ("matplotlib", "matplotlib")):
        feature.checks.append(package_check(module, dist, fix=fix))
    return feature


def hardware_section() -> Dict[str, Any]:
    gpu_rows = gpus()
    return {
        "platform": f"{platform.system()} {platform.machine()}",
        "cpus": os.cpu_count(),
        "memory_gb": memory_gb(),
        "gpus": gpu_rows,
    }


def storage_section() -> Dict[str, Any]:
    from genomics import workspace

    rows = {
        "data root (GENOMICS_DATA_ROOT)": workspace.data_root(),
        "results root (GENOMICS_RESULTS_ROOT)": workspace.results_root(),
        "visualizer cache": workspace.cache_path("visualizer"),
    }
    spaces = {label: free_space(path) for label, path in rows.items()}
    try:
        from genomics.core.data_registry import resolve_dataset

        dataset = resolve_dataset(workspace.DEFAULT_DATASET_ID)
        default_dataset = {"id": workspace.DEFAULT_DATASET_ID, "path": str(dataset.path), "exists": (Path(dataset.path) / "dataset_metadata.json").is_file()}
    except Exception as exc:  # registry problems should not break the report
        default_dataset = {"id": workspace.DEFAULT_DATASET_ID, "path": None, "exists": False, "error": str(exc)}
    return {"locations": spaces, "default_dataset": default_dataset}


def report(include_local: bool = True, backend: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    features = [core_feature(), import_feature(), alphagenome_feature(include_local, backend), training_feature()]
    return {
        "genomics": dist_version("genomics"),
        "features": [{**asdict(f), "status": f.status} for f in features],
        "hardware": hardware_section(),
        "storage": storage_section(),
    }


# -- output ---------------------------------------------------------------------------------------

def _fmt_gb(value: Optional[float]) -> str:
    if value is None:
        return "?"
    return f"{value / 1024:.1f} TB" if value >= 1024 else f"{value:.0f} GB"


def print_report(data: Dict[str, Any]) -> None:
    print(f"genomics {data.get('genomics') or '(not installed as a package)'}\n")
    for feature in data["features"]:
        print(f"{SYMBOLS[feature['status']]} {feature['title']}: {feature['purpose']}")
        for check in feature["checks"]:
            print(f"    {SYMBOLS[check['status']]} {check['name']}: {check['detail']}")
            if check["fix"] and check["status"] in (MISSING, WARN, INFO):
                print(f"        → {check['fix']}")
        print()
    hw = data["hardware"]
    gpu_text = "; ".join(f"{g['name']} {g['memory']}".strip() for g in hw["gpus"]) or "no NVIDIA GPU detected"
    print(f"· Hardware: {hw['platform']}, {hw['cpus']} CPUs, {_fmt_gb(hw['memory_gb'])} RAM, {gpu_text}")
    print("· Storage:")
    for label, space in data["storage"]["locations"].items():
        if space is None:
            continue
        note = "" if space["exists"] else " (does not exist yet)"
        print(f"    {label}: {space['path']}{note}, {_fmt_gb(space['free_gb'])} free")
    ds = data["storage"]["default_dataset"]
    print(f"    default dataset {ds['id']}: {ds['path']} ({'found' if ds['exists'] else 'not found here; import or open another dataset'})")
    print("\nSizes and hardware needed per feature: docs/getting-started/requirements.md")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genomics doctor", description="Check which genomics features this machine can run, and how to enable the rest")
    add_arguments(parser)
    return parser


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument("--json", action="store_true", help="Print the report as JSON")
    parser.add_argument("--no-local-server", action="store_true", help="Skip probing the local AlphaGenome server environment")
    return parser


def run(args: argparse.Namespace) -> int:
    data = report(include_local=not args.no_local_server)
    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print_report(data)
    core = next(f for f in data["features"] if f["key"] == "core")
    return 1 if core["status"] == MISSING else 0


def main(argv: Optional[List[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
