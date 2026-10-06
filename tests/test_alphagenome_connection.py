"""Self-hosted AlphaGenome connections and the visualizer's AlphaGenome backend settings."""
from __future__ import annotations

import http.client
import json
import shutil
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from genomics.core.alphagenome_connection import (
    ADDRESS_ENV,
    CA_CERT_ENV,
    backend_configured,
    endpoint_from_env,
    parse_address,
    speaks_tls,
)
from genomics.visualizer.alphagenome import AlphaGenomeBackend, LocalAlphaGenomeServer
from genomics.workflows.alphagenome import local_server as server_module


def _free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch, tmp_path):
    for name in (ADDRESS_ENV, CA_CERT_ENV, "ALPHAGENOME_API_KEY"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))  # no ~/.env key


# -- address parsing -----------------------------------------------------------------------------

@pytest.mark.parametrize(
    "address, target, transport",
    [
        ("grpc://10.0.0.5:50051", "10.0.0.5:50051", "insecure"),
        ("grpcs://gpu-box", "gpu-box:50051", "tls"),
        ("https://example.org/", "example.org:443", "tls"),
        ("10.0.0.5", "10.0.0.5:50051", "auto"),
        ("[::1]:6000", "[::1]:6000", "auto"),
    ],
)
def test_parse_address(address, target, transport):
    endpoint = parse_address(address)
    assert (endpoint.target, endpoint.transport) == (target, transport)


@pytest.mark.parametrize("address", ["", "ftp://host", "host:0", "host:99999", "host:port"])
def test_parse_address_rejects_bad_input(address):
    with pytest.raises(ValueError):
        parse_address(address)


def test_env_selects_self_hosted_server(monkeypatch):
    assert endpoint_from_env() is None and not backend_configured()
    monkeypatch.setenv(ADDRESS_ENV, "grpcs://gpu-box:6000")
    monkeypatch.setenv(CA_CERT_ENV, "/certs/ca.crt")
    endpoint = endpoint_from_env()
    assert endpoint.target == "gpu-box:6000" and endpoint.transport == "tls" and endpoint.ca_cert == "/certs/ca.crt"
    assert backend_configured()  # no API key needed for a self-hosted server


# -- real gRPC round trips against a fake AlphaGenome service ---------------------------------------

def _fake_server(port, credentials=None):
    grpc = pytest.importorskip("grpc")
    pytest.importorskip("alphagenome")
    from concurrent import futures

    from alphagenome.protos import dna_model_service_pb2, dna_model_service_pb2_grpc

    class Servicer(dna_model_service_pb2_grpc.DnaModelServiceServicer):
        def GetMetadata(self, request, context):
            yield dna_model_service_pb2.MetadataResponse()

    server = grpc.server(futures.ThreadPoolExecutor(max_workers=2))
    dna_model_service_pb2_grpc.add_DnaModelServiceServicer_to_server(Servicer(), server)
    if credentials is None:
        server.add_insecure_port(f"127.0.0.1:{port}")
    else:
        server.add_secure_port(f"127.0.0.1:{port}", credentials)
    server.start()
    return server


def _self_signed_certs(directory: Path):
    if shutil.which("openssl") is None:
        pytest.skip("openssl not available")
    run = lambda *args: subprocess.run(["openssl", *args], check=True, capture_output=True)  # noqa: E731
    run("req", "-x509", "-newkey", "rsa:2048", "-nodes", "-keyout", str(directory / "ca.key"), "-out", str(directory / "ca.crt"), "-days", "2", "-subj", "/CN=Test CA")
    run("req", "-newkey", "rsa:2048", "-nodes", "-keyout", str(directory / "server.key"), "-out", str(directory / "server.csr"), "-subj", "/CN=127.0.0.1")
    ext = directory / "san.ext"
    ext.write_text("subjectAltName=IP:127.0.0.1,DNS:localhost\n", encoding="utf-8")
    run("x509", "-req", "-in", str(directory / "server.csr"), "-CA", str(directory / "ca.crt"), "-CAkey", str(directory / "ca.key"),
        "-CAcreateserial", "-out", str(directory / "server.crt"), "-days", "2", "-extfile", str(ext))
    return directory / "ca.crt", directory / "server.crt", directory / "server.key"


def test_check_connection_plaintext_with_autodetect():
    from genomics.core.alphagenome_connection import check_connection, create_dna_client

    port = _free_port()
    server = _fake_server(port)
    try:
        assert not speaks_tls("127.0.0.1", port, timeout=2)
        result = check_connection(f"127.0.0.1:{port}", timeout=5)
        assert result["ok"], result
        assert result["transport"] == "insecure" and result["target"] == f"grpc://127.0.0.1:{port}"
        client = create_dna_client(address=f"grpc://127.0.0.1:{port}", timeout=5)
        assert type(client).__name__ == "DnaClient"
    finally:
        server.stop(None)


def test_check_connection_tls_needs_the_ca(tmp_path):
    grpc = pytest.importorskip("grpc")
    from genomics.core.alphagenome_connection import check_connection

    ca, cert, key = _self_signed_certs(tmp_path)
    port = _free_port()
    server = _fake_server(port, grpc.ssl_server_credentials(((key.read_bytes(), cert.read_bytes()),)))
    try:
        assert speaks_tls("127.0.0.1", port, timeout=2)
        untrusted = check_connection(f"127.0.0.1:{port}", timeout=2)
        assert not untrusted["ok"] and "CA certificate" in untrusted["message"]
        trusted = check_connection(f"127.0.0.1:{port}", ca_cert=str(ca), timeout=5)
        assert trusted["ok"] and trusted["transport"] == "tls", trusted
    finally:
        server.stop(None)


def test_check_connection_reports_unreachable_and_missing_ca(tmp_path):
    pytest.importorskip("grpc")
    from genomics.core.alphagenome_connection import check_connection

    result = check_connection(f"grpc://127.0.0.1:{_free_port()}", timeout=1)
    assert not result["ok"] and "could not connect" in result["message"]
    assert result["target"].startswith("grpc://127.0.0.1:")
    result = check_connection("grpcs://127.0.0.1:1", ca_cert=str(tmp_path / "missing.crt"))
    assert not result["ok"] and "not found" in result["message"]


# -- visualizer backend ------------------------------------------------------------------------------

def _local(tmp_path, port=None, launcher_code=None):
    server_dir = tmp_path / "alphagenome_research"
    server_dir.mkdir(exist_ok=True)
    (server_dir / "server.py").write_text("print('hi')\n", encoding="utf-8")
    launcher = tmp_path / "fake_model_server.py"
    launcher.write_text(launcher_code or "print('hi')\n", encoding="utf-8")
    # pretend alphagenome_research + CUDA jax are installed
    server_module._PROBE_CACHE[sys.executable] = server_module.Probe(sys.executable, server_module.PROBE_GPU)
    return LocalAlphaGenomeServer(server_dir, Path(sys.executable), tmp_path / "logs", port=port or _free_port(), launcher=launcher)


def test_backend_settings_persist_and_drive_lab_env(tmp_path):
    settings_path = tmp_path / "settings.json"
    backend = AlphaGenomeBackend(_local(tmp_path), settings_path=settings_path)
    assert backend.settings["mode"] == "cloud"
    assert backend.reasons() and "ALPHAGENOME_API_KEY" in backend.reasons()[0]

    ca = tmp_path / "ca.crt"
    ca.write_text("pem", encoding="utf-8")
    backend.update({"mode": "remote", "address": "grpcs://gpu-box:50051", "ca_cert": str(ca)})
    assert backend.reasons() == []
    env = backend.child_env({ADDRESS_ENV: "stale", "PATH": "/bin"})
    assert env[ADDRESS_ENV] == "grpcs://gpu-box:50051" and env[CA_CERT_ENV] == str(ca) and env["PATH"] == "/bin"

    reloaded = AlphaGenomeBackend(_local(tmp_path), settings_path=settings_path)
    assert reloaded.settings == backend.settings

    with pytest.raises(ValueError):
        backend.update({"mode": "remote", "address": "ftp://nope"})
    with pytest.raises(ValueError):
        backend.update({"mode": "remote", "address": "host:50051", "ca_cert": str(tmp_path / "missing.crt")})

    backend.update({"mode": "cloud"})
    assert ADDRESS_ENV not in backend.child_env({ADDRESS_ENV: "stale"})

    override = AlphaGenomeBackend(_local(tmp_path), settings_path=settings_path, address="grpc://cli-host:1234")
    assert override.settings["mode"] == "remote" and override.settings["address"] == "grpc://cli-host:1234"


def test_local_server_lifecycle_and_tls_detection(tmp_path, monkeypatch):
    port = _free_port()
    # Stands in for model_server.py: must receive --server-dir, --host and --port.
    code = (
        "import argparse, socket, time\n"
        "p = argparse.ArgumentParser(); p.add_argument('--server-dir'); p.add_argument('--host'); p.add_argument('--port', type=int)\n"
        "a = p.parse_args()\n"
        "s = socket.socket(); s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)\n"
        "time.sleep(0.3); s.bind((a.host, a.port)); s.listen()\n"
        "print('listening', a.server_dir, flush=True)\n"
        "time.sleep(60)\n"
    )
    local = _local(tmp_path, port=port, launcher_code=code)
    backend = AlphaGenomeBackend(local, settings_path=tmp_path / "settings.json")
    backend.update({"mode": "local"})
    assert local.state() == "stopped" and local.address == f"grpc://127.0.0.1:{port}"
    assert "not running" in backend.reasons()[0]

    local.start()
    try:
        deadline = time.time() + 10
        while local.state() != "ready" and time.time() < deadline:
            time.sleep(0.1)
        assert local.state() == "ready" and backend.reasons() == []
        assert backend.child_env({})[ADDRESS_ENV] == f"grpc://127.0.0.1:{port}"
        with pytest.raises(RuntimeError):
            LocalAlphaGenomeServer(local.server_dir, Path(sys.executable), tmp_path, port=port, launcher=local.launcher).start()  # port taken
    finally:
        local.stop()
    assert local.state() == "stopped"
    assert any(line.startswith("listening") and str(local.server_dir) in line for line in local.log_tail())

    certs = local.server_dir / "certs"
    certs.mkdir()
    for name in ("server.crt", "server.key", "ca.crt"):
        (certs / name).write_text("pem", encoding="utf-8")
    assert local.tls and local.address.startswith("grpcs://") and local.ca_cert == str(certs / "ca.crt")


def test_local_server_reports_missing_checkout(tmp_path):
    local = LocalAlphaGenomeServer(None, tmp_path / "nope" / "python", tmp_path)
    reasons = local.reasons()
    assert any("alphagenome_research" in r for r in reasons) and any("Python interpreter" in r for r in reasons)
    with pytest.raises(RuntimeError):
        local.start()


def test_default_server_python_prefers_a_gpu_env(tmp_path, monkeypatch):
    base = tmp_path / "conda"
    pythons = {}
    for env in ("alphagenome", "cpu-only", "gpu"):
        python = base / "envs" / env / "bin" / "python"
        python.parent.mkdir(parents=True)
        python.write_text("", encoding="utf-8")
        pythons[env] = python
    codes = {str(pythons["alphagenome"]): 2, str(pythons["cpu-only"]): 2, str(pythons["gpu"]): 0, sys.executable: 1}
    monkeypatch.setattr(server_module, "probe_python", lambda python: codes.get(str(python), 1))
    monkeypatch.delenv(server_module.SERVER_PYTHON_ENV, raising=False)
    (base / "bin").mkdir()
    (base / "bin" / "conda").write_text("", encoding="utf-8")
    monkeypatch.setenv("CONDA_EXE", str(base / "bin" / "conda"))
    assert server_module.default_server_python(tmp_path) == pythons["gpu"]
    codes[str(pythons["gpu"])] = 2
    assert server_module.default_server_python(tmp_path) == pythons["alphagenome"]
    monkeypatch.setenv(server_module.SERVER_PYTHON_ENV, "/opt/py")
    assert server_module.default_server_python() == Path("/opt/py")


def test_probe_classifies_environments():
    def info(**versions):
        modules = {"alphagenome": True, "alphagenome.io": True, "alphagenome_research": True, "jax": True}
        return {"python": "3.11.0", "modules": modules, "versions": {"alphagenome": "0.7.0", "jax": "0.10.2", "jaxlib": "0.10.2", **versions}}

    assert server_module._analyse("py", info(**{"jax-cuda12-plugin": "0.10.2"})).code == server_module.PROBE_GPU
    cpu = server_module._analyse("py", info())
    assert cpu.code == server_module.PROBE_CPU and "no CUDA plugin" in cpu.problems[0]
    mismatch = server_module._analyse("py", info(**{"jax-cuda12-plugin": "0.10.1"}))
    assert mismatch.code == server_module.PROBE_CPU and "does not match jaxlib" in mismatch.problems[0]
    old_sdk = info(**{"jax-cuda12-plugin": "0.10.2", "alphagenome": "0.6.1"})
    old_sdk["modules"]["alphagenome.io"] = False
    missing = server_module._analyse("py", old_sdk)
    assert missing.code == server_module.PROBE_MISSING and "alphagenome>=0.7" in missing.problems[0]
    # A real interpreter is probed without importing jax; this one has no alphagenome_research.
    probed = server_module.probe(Path(sys.executable), refresh=True)
    assert probed.info["python"].startswith("%d.%d" % sys.version_info[:2])


def test_local_server_reasons_point_to_setup(tmp_path):
    server_dir = tmp_path / "alphagenome_research"
    server_dir.mkdir()
    (server_dir / "server.py").write_text("", encoding="utf-8")
    python = tmp_path / "py"
    server_module._PROBE_CACHE[str(python)] = server_module.Probe(str(python), server_module.PROBE_CPU, ["jax 0.10.2 has no CUDA plugin (CPU only)"])
    python.write_text("", encoding="utf-8")
    local = LocalAlphaGenomeServer(server_dir, python, tmp_path)
    reasons = local.reasons()
    assert len(reasons) == 1 and "needs a GPU" in reasons[0] and "genomics alphagenome server setup" in reasons[0]
    assert local.address == "grpc://127.0.0.1:50051" and local.info()["host"] == "127.0.0.1"
    shared = LocalAlphaGenomeServer(server_dir, python, tmp_path, host="0.0.0.0", port=50070)
    assert shared.address == "grpc://127.0.0.1:50070"


def test_server_command_and_env(tmp_path, monkeypatch):
    cmd = server_module.server_command(Path("/envs/ag/bin/python"), tmp_path, "127.0.0.1", 50052, ["--allow-cpu"])
    assert cmd[:2] == ["/envs/ag/bin/python", str(server_module.MODEL_SERVER)]
    assert cmd[2:] == ["--server-dir", str(tmp_path), "--host", "127.0.0.1", "--port", "50052", "--allow-cpu"]
    env = server_module.server_env(Path("/envs/ag/bin/python"), {"PATH": "/usr/bin", "PYTHONPATH": "/genomics/src", "CONDA_PREFIX": "/envs/genomics"})
    assert env["PATH"].startswith("/envs/ag/bin") and "PYTHONPATH" not in env and "CONDA_PREFIX" not in env


def test_model_server_is_standalone():
    """model_server.py runs in the server environment, which has no genomics package."""
    source = server_module.MODEL_SERVER.read_text(encoding="utf-8")
    assert "genomics." not in source.split('"""', 2)[2].replace("genomics.core.alphagenome_connection", "")
    done = subprocess.run([sys.executable, "-S", str(server_module.MODEL_SERVER), "--help"], capture_output=True, text=True, timeout=30)
    assert done.returncode == 0 and "--server-dir" in done.stdout
    # XLA must be told not to preallocate before jax is first imported (check_device), or the server
    # takes 75% of GPU memory (~90 GB of a DGX Spark's unified memory).
    main = source[source.index("def main("):]
    assert main.index('"XLA_PYTHON_CLIENT_PREALLOCATE", "false"') < main.index("check_device(args.allow_cpu)")


def test_model_server_metadata_response():
    pytest.importorskip("alphagenome")
    pd = pytest.importorskip("pandas")
    import importlib.util

    spec = importlib.util.spec_from_file_location("model_server", server_module.MODEL_SERVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from alphagenome.models import dna_client, dna_output

    rna = pd.DataFrame({"name": ["UBERON:0002107 total RNA-seq"], "strand": ["+"], "ontology_curie": ["UBERON:0002107"]})
    junctions = pd.DataFrame({"name": ["UBERON:0002107"], "ontology_curie": ["UBERON:0002107"]})

    class FakeModel:
        def output_metadata(self, organism):
            return dna_output.OutputMetadata(rna_seq=rna, splice_junctions=junctions)

    response = module.metadata_response(FakeModel(), dna_client.Organism.HOMO_SAPIENS)
    parsed = dna_client.construct_output_metadata(iter([response]))
    assert list(parsed.rna_seq["name"]) == ["UBERON:0002107 total RNA-seq"]
    assert list(parsed.splice_junctions["name"]) == ["UBERON:0002107"] and parsed.cage is None


def test_backend_client_uses_selected_endpoint(tmp_path, monkeypatch):
    from genomics.core import alphagenome_connection

    backend = AlphaGenomeBackend(_local(tmp_path), settings_path=tmp_path / "settings.json")
    backend.update({"mode": "remote", "address": "grpc://gpu-box:50051"})
    calls = []
    monkeypatch.setattr(alphagenome_connection, "create_dna_client", lambda **kw: calls.append(kw) or "client")
    assert backend.create_client(timeout=5) == "client"
    assert calls[-1]["address"] == "grpc://gpu-box:50051" and calls[-1]["timeout"] == 5


def test_http_alphagenome_routes(tmp_path):
    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.server import Handler, Server, VisualizerApp

    backend = AlphaGenomeBackend(_local(tmp_path), settings_path=tmp_path / "settings.json")
    app = VisualizerApp(DatasetCatalog(), cache_dir=None, memory_bytes=64 << 20, workers=1, runs_roots=[], alphagenome=backend)
    port = _free_port()
    server = Server(("127.0.0.1", port), type("H", (Handler,), {"app": app}))
    threading.Thread(target=server.serve_forever, daemon=True).start()

    def call(path, method="GET", body=None):
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
        conn.request(method, path, body=json.dumps(body) if body is not None else None, headers={"Content-Type": "application/json"})
        res = conn.getresponse()
        data = json.loads(res.read())
        conn.close()
        return res.status, data

    try:
        status, data = call("/api/alphagenome")
        assert status == 200 and data["settings"]["mode"] == "cloud" and data["local"]["state"] == "stopped"
        status, data = call("/api/alphagenome/settings", "POST", {"mode": "remote", "address": "grpc://127.0.0.1:1"})
        assert status == 200 and data["label"] == "Remote server (grpc://127.0.0.1:1)"
        assert call("/api/alphagenome/settings", "POST", {"mode": "bogus"})[0] == 400
        status, data = call("/api/alphagenome/test", "POST", {"mode": "remote", "address": ""})
        assert status == 200 and not data["ok"]
        status, data = call("/api/alphagenome/test", "POST", {"mode": "local"})
        assert not data["ok"] and "stopped" in data["message"]
        status, data = call("/api/perturb/models")
        assert status == 200 and data["runs"] == [] and data["backend"]["label"].startswith("Remote server")
        assert call("/api/perturb/score?sample=S1")[0] == 400  # no model loaded
    finally:
        server.shutdown()
        server.server_close()
        app.shutdown()


def test_checkpoint_check_resolves_relative_results_dir_against_repo_root(tmp_path, monkeypatch):
    from genomics.predictors.genotype_based import config as config_module
    from genomics.predictors.genotype_based.apps import genomics_workbench

    checkpoint = tmp_path / "repo" / "results" / "runs" / "exp" / "models" / "best_accuracy.pt"
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_bytes(b"")
    monkeypatch.setattr(config_module, "load_config", lambda path: object())
    monkeypatch.setattr(config_module, "get_experiment_runs_dir", lambda config: Path("results/runs"))
    monkeypatch.setattr(config_module, "generate_experiment_name", lambda config: "exp")
    monkeypatch.setattr(genomics_workbench, "repo_root", lambda: tmp_path / "repo")
    monkeypatch.chdir(tmp_path)  # not the repo root, as when `genomics visualize` runs from $HOME
    assert genomics_workbench._pigmentation_checkpoint_exists(Path("cfg.yaml"))


def test_read_dotenv_needs_no_python_dotenv(tmp_path):
    from genomics.core.alphagenome_connection import read_dotenv

    env = tmp_path / ".env"
    env.write_text("# comment\nexport ALPHAGENOME_API_KEY='abc 123'\nOTHER=x # trailing\nEMPTY=\nbroken line\n", encoding="utf-8")
    values = read_dotenv(env)
    assert values == {"ALPHAGENOME_API_KEY": "abc 123", "OTHER": "x", "EMPTY": ""}
    assert read_dotenv(tmp_path / "missing") == {}
