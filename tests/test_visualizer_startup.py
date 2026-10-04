"""Startup of ``genomics visualize``: port fallback, reusing a running visualizer, URLs, browser."""
from __future__ import annotations

import json
import socket
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

from genomics import cli as genomics_cli
from genomics.visualizer import server as visualizer_server
from genomics.visualizer import startup
from genomics.visualizer.cli_args import build_arg_parser


def _listening_socket() -> socket.socket:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("127.0.0.1", 0))
    sock.listen(1)
    return sock


class _StubServer:
    """An HTTP server answering ``/api/status`` with ``status`` (None: 404 on every path)."""

    def __init__(self, status):
        body = json.dumps(status).encode("utf-8") if status is not None else None

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):  # noqa: N802
                if body is None or self.path != "/api/status":
                    self.send_response(404)
                    self.end_headers()
                    return
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *args):
                pass

        self.httpd = HTTPServer(("127.0.0.1", 0), Handler)
        self.port = self.httpd.server_address[1]
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *exc):
        self.httpd.shutdown()
        self.httpd.server_close()


def _status(*datasets):
    return {"app": startup.APP_NAME, "version": "1.0", "uptime": 1.0, "caches": {}, "datasets": list(datasets)}


def _no_create_app(_args):
    raise AssertionError("the datasets must not be loaded")


def test_server_url_uses_a_reachable_host():
    assert startup.server_url("127.0.0.1", 8780) == "http://127.0.0.1:8780/"
    assert startup.server_url("0.0.0.0", 8780) == "http://localhost:8780/"
    assert startup.server_url("::", 8781) == "http://localhost:8781/"
    assert startup.server_url("::1", 8781) == "http://[::1]:8781/"


def test_bind_first_free_skips_ports_in_use():
    taken = _listening_socket()
    port = taken.getsockname()[1]
    try:
        server, chosen = startup.bind_first_free(lambda p: visualizer_server.Server(("127.0.0.1", p), visualizer_server.Handler), port, 20)
        try:
            assert chosen > port
            assert server.server_address[1] == chosen
        finally:
            server.server_close()
    finally:
        taken.close()


def test_bind_server_falls_back_only_without_explicit_port(monkeypatch):
    taken = _listening_socket()
    port = taken.getsockname()[1]
    try:
        monkeypatch.setattr(visualizer_server, "DEFAULT_PORT", port)
        server = visualizer_server.bind_server("127.0.0.1", None)
        try:
            assert server.server_address[1] > port
        finally:
            server.server_close()
        with pytest.raises(startup.StartupError, match="in use by another program"):
            visualizer_server.bind_server("127.0.0.1", port)
    finally:
        taken.close()


def test_probe_recognises_only_a_visualizer():
    with _StubServer(_status()) as stub:
        assert startup.probe_visualizer("127.0.0.1", stub.port)["app"] == startup.APP_NAME
    with _StubServer({"version": "1.0", "uptime": 3.0, "datasets": [], "caches": {}}) as stub:
        assert startup.probe_visualizer("127.0.0.1", stub.port) is not None  # started before the "app" field
    with _StubServer({"status": "ok"}) as stub:
        assert startup.probe_visualizer("127.0.0.1", stub.port) is None
    with _StubServer(None) as stub:
        assert startup.probe_visualizer("127.0.0.1", stub.port) is None
    free = _listening_socket()
    port = free.getsockname()[1]
    free.close()
    assert startup.probe_visualizer("127.0.0.1", port) is None


def test_can_reuse_needs_every_requested_dataset(tmp_path):
    status = _status({"id": "1kg_high_coverage", "path": str(tmp_path.resolve())})
    parse = build_arg_parser().parse_args
    assert startup.can_reuse(parse([]), status)
    assert startup.can_reuse(parse(["--dataset", str(tmp_path)]), status)
    assert startup.can_reuse(parse(["--dataset-id", "1kg_high_coverage"]), status)
    assert not startup.can_reuse(parse(["--dataset", str(tmp_path / "other")]), status)
    assert not startup.can_reuse(parse(["--dataset-id", "other"]), status)
    assert not startup.can_reuse(parse(["--annotations", "extra.tsv"]), status)


def test_main_reuses_a_running_visualizer(monkeypatch, tmp_path, capsys):
    opened = []
    monkeypatch.setattr(visualizer_server, "create_app", _no_create_app)
    monkeypatch.setattr(visualizer_server.webbrowser, "open", lambda url: opened.append(url))
    monkeypatch.setattr(startup, "display_available", lambda: True)
    with _StubServer(_status({"id": "tiny", "path": str(tmp_path.resolve())})) as stub:
        code = visualizer_server.main(["--port", str(stub.port), "--dataset", str(tmp_path), "--open"])
    assert code == 0
    assert opened == [f"http://127.0.0.1:{stub.port}/"]
    assert "already running" in capsys.readouterr().out


def test_main_reports_a_busy_explicit_port(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(visualizer_server, "create_app", _no_create_app)
    with _StubServer(_status()) as stub:
        assert visualizer_server.main(["--port", str(stub.port), "--dataset", str(tmp_path)]) == 1
    assert "without the requested datasets" in capsys.readouterr().err
    with _StubServer(None) as stub:
        assert visualizer_server.main(["--port", str(stub.port)]) == 1
    assert "in use by another program" in capsys.readouterr().err


def test_open_without_display_prints_the_url(monkeypatch, capsys):
    monkeypatch.setattr(startup, "display_available", lambda: False)
    monkeypatch.setattr(visualizer_server.webbrowser, "open", lambda url: pytest.fail("no browser without a display"))
    visualizer_server.open_in_browser("http://127.0.0.1:8780/", wait=True)
    assert "open http://127.0.0.1:8780/" in capsys.readouterr().err


def test_display_and_ssh_detection(monkeypatch):
    for name in ("DISPLAY", "WAYLAND_DISPLAY", "SSH_CONNECTION", "SSH_TTY"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(startup.sys, "platform", "linux")
    assert not startup.display_available()
    monkeypatch.setenv("WAYLAND_DISPLAY", "wayland-0")
    assert startup.display_available()
    assert startup.ssh_tunnel_hint("127.0.0.1", 8780) is None
    monkeypatch.setenv("SSH_CONNECTION", "10.0.0.2 51000 10.0.0.5 22")
    hint = startup.ssh_tunnel_hint("127.0.0.1", 8781)
    assert hint is not None and hint.startswith("ssh -N -L 8781:localhost:8781 ")
    assert startup.ssh_tunnel_hint("0.0.0.0", 8781) is None


def test_genomics_visualize_forwards_every_option(tmp_path):
    argv = [
        "--dataset", str(tmp_path / "a"), "--dataset-id", "1kg_high_coverage", "--annotations", "x.tsv",
        "--consensus-dataset-dir", "c", "--runs-root", "r", "--gtf", "g.feather", "--cache-dir", "cache",
        "--no-disk-cache", "--memory-mb", "512", "--workers", "3", "--model-window", "16384",
        "--jobs-dir", "jobs", "--no-jobs", "--pigmentation-config", "p.yaml",
        "--alphagenome-address", "grpc://h:50051", "--alphagenome-ca-cert", "ca.crt",
        "--alphagenome-server-dir", "srv", "--alphagenome-server-python", "py",
        "--host", "0.0.0.0", "--port", "9000", "--open", "--no-add-datasets", "--no-remote", "--verbose",
    ]
    direct = vars(build_arg_parser().parse_args(argv))
    args = genomics_cli.build_parser().parse_args(["visualize", *argv])
    forwarded = vars(build_arg_parser().parse_args([str(a) for a in genomics_cli._visualizer_args(args)]))
    assert forwarded == direct
    defaults = genomics_cli.build_parser().parse_args(["visualize"])
    assert vars(build_arg_parser().parse_args([str(a) for a in genomics_cli._visualizer_args(defaults)]))["port"] is None
