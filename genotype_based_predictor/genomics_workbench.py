from __future__ import annotations

import argparse
import http.client
import json
import os
import signal
import socket
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urlparse


DEFAULT_DATASET_DIR = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
DEFAULT_RUNS_ROOT = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage_runs")
DEFAULT_ALIGNED_TSV_ROOT = Path("genotype_based_predictor/aligned_dna_genes_1000_all")


@dataclass
class ViewerSpec:
    key: str
    title: str
    description: str
    module: str
    port: int
    args: List[str]
    enabled: bool = True
    disabled_reason: str = ""


@dataclass
class ViewerProcess:
    spec: ViewerSpec
    process: Optional[subprocess.Popen]
    log_path: Path

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.spec.port}"

    @property
    def route_url(self) -> str:
        return f"/apps/{self.spec.key}/"

    @property
    def running(self) -> bool:
        return self.process is not None and self.process.poll() is None

    @property
    def returncode(self) -> Optional[int]:
        return None if self.process is None else self.process.poll()


def _port_available(host: str, port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(0.2)
        return sock.connect_ex((host, port)) != 0


def _has_tsvs(path: Path) -> bool:
    if path.is_file() and path.suffix == ".tsv":
        return True
    return path.is_dir() and any(path.glob("*.tsv"))


def _build_specs(args: argparse.Namespace) -> List[ViewerSpec]:
    dataset_dir = Path(args.dataset_dir).resolve()
    runs_root = Path(args.runs_root).resolve()
    aligned_tsv_root = Path(args.aligned_tsv_root).resolve()

    specs = [
        ViewerSpec(
            key="dataset",
            title="Dataset Browser",
            description="Explore individuos, genes, populacoes, superpopulacoes e metadados do dataset.",
            module="genotype_based_predictor.dataset_browser",
            port=args.dataset_port,
            args=[str(dataset_dir), "--host", "127.0.0.1", "--port", str(args.dataset_port)],
            enabled=dataset_dir.exists(),
            disabled_reason=f"Dataset nao encontrado: {dataset_dir}",
        ),
        ViewerSpec(
            key="view-builder",
            title="View Builder",
            description="Construa arquivos .view.json via formulario, sem editar JSON manualmente.",
            module="genotype_based_predictor.view_builder",
            port=args.view_builder_port,
            args=[str(dataset_dir), "--host", "127.0.0.1", "--port", str(args.view_builder_port)],
            enabled=dataset_dir.exists(),
            disabled_reason=f"Dataset nao encontrado: {dataset_dir}",
        ),
        ViewerSpec(
            key="experiments",
            title="Experiment Dashboard",
            description="Compare runs, metricas, configs, checkpoints e matrizes de confusao.",
            module="genotype_based_predictor.experiment_dashboard",
            port=args.experiment_port,
            args=[str(runs_root), "--host", "127.0.0.1", "--port", str(args.experiment_port)],
            enabled=runs_root.exists(),
            disabled_reason=f"Diretorio de runs nao encontrado: {runs_root}",
        ),
        ViewerSpec(
            key="tracks",
            title="AlphaGenome Track Viewer",
            description="Visualize tracks AlphaGenome por individuo, gene, haplotipo e output.",
            module="genotype_based_predictor.alphagenome_track_viewer",
            port=args.track_port,
            args=[str(dataset_dir), "--host", "127.0.0.1", "--port", str(args.track_port)],
            enabled=dataset_dir.exists(),
            disabled_reason=f"Dataset nao encontrado: {dataset_dir}",
        ),
        ViewerSpec(
            key="alignment",
            title="Aligned DNA Viewer",
            description="Compare sequencias alinhadas, mutacoes e X/gaps a partir de TSVs por gene.",
            module="genotype_based_predictor.aligned_dna_viewer",
            port=args.alignment_port,
            args=[str(aligned_tsv_root), "--host", "127.0.0.1", "--port", str(args.alignment_port), "--dataset-dir", str(dataset_dir)],
            enabled=_has_tsvs(aligned_tsv_root),
            disabled_reason=f"Nenhum .tsv encontrado em: {aligned_tsv_root}",
        ),
    ]
    return specs


class Workbench:
    def __init__(self, specs: List[ViewerSpec], log_dir: Path):
        self.specs = specs
        self.log_dir = log_dir
        self.processes: Dict[str, ViewerProcess] = {}

    def start(self) -> None:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        for spec in self.specs:
            log_path = self.log_dir / f"{spec.key}.log"
            if not spec.enabled:
                self.processes[spec.key] = ViewerProcess(spec, None, log_path)
                continue
            if not _port_available("127.0.0.1", spec.port):
                spec.enabled = False
                spec.disabled_reason = f"Porta em uso: {spec.port}"
                self.processes[spec.key] = ViewerProcess(spec, None, log_path)
                continue
            log_file = open(log_path, "w", encoding="utf-8")
            cmd = [sys.executable, "-m", spec.module] + spec.args
            process = subprocess.Popen(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                cwd=Path.cwd(),
                env=os.environ.copy(),
                text=True,
            )
            self.processes[spec.key] = ViewerProcess(spec, process, log_path)
        time.sleep(0.8)

    def stop(self) -> None:
        for viewer in self.processes.values():
            if viewer.process is None or viewer.process.poll() is not None:
                continue
            viewer.process.terminate()
        deadline = time.time() + 5
        for viewer in self.processes.values():
            proc = viewer.process
            if proc is None:
                continue
            remaining = max(deadline - time.time(), 0.1)
            try:
                proc.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=2)

    def payload(self) -> Dict[str, object]:
        return {
            "viewers": [
                {
                    "key": viewer.spec.key,
                    "title": viewer.spec.title,
                    "description": viewer.spec.description,
                    "url": viewer.route_url,
                    "internal_url": viewer.url,
                    "port": viewer.spec.port,
                    "enabled": viewer.spec.enabled,
                    "running": viewer.running,
                    "returncode": viewer.returncode,
                    "disabled_reason": viewer.spec.disabled_reason,
                    "log_path": str(viewer.log_path),
                }
                for viewer in self.processes.values()
            ]
        }


class WorkbenchHandler(BaseHTTPRequestHandler):
    workbench: Workbench

    def log_message(self, fmt: str, *args):
        return

    def _send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self, text: str, status: int = 200) -> None:
        body = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self._send_html(INDEX_HTML)
        elif parsed.path == "/api/status":
            self._send_json(self.workbench.payload())
        elif parsed.path.startswith("/apps/"):
            self._proxy_to_viewer("GET")
        else:
            self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path.startswith("/apps/"):
            self._proxy_to_viewer("POST")
        else:
            self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)

    def _proxy_to_viewer(self, method: str) -> None:
        parsed = urlparse(self.path)
        parts = parsed.path.split("/")
        if len(parts) < 3 or parts[1] != "apps":
            self._send_json({"error": "invalid app route"}, status=HTTPStatus.NOT_FOUND)
            return

        key = parts[2]
        viewer = self.workbench.processes.get(key)
        if viewer is None:
            self._send_json({"error": f"unknown app: {key}"}, status=HTTPStatus.NOT_FOUND)
            return
        if not viewer.running:
            self._send_json(
                {"error": f"app is not running: {key}", "log_path": str(viewer.log_path)},
                status=HTTPStatus.BAD_GATEWAY,
            )
            return

        upstream_path = "/" + "/".join(parts[3:])
        if upstream_path == "/":
            upstream_path = "/"
        if parsed.query:
            upstream_path = f"{upstream_path}?{parsed.query}"

        body = None
        if method in {"POST", "PUT", "PATCH"}:
            length = int(self.headers.get("Content-Length", "0") or "0")
            body = self.rfile.read(length) if length > 0 else None

        headers = {
            key: value
            for key, value in self.headers.items()
            if key.lower() not in {"host", "connection", "content-length", "accept-encoding"}
        }
        if body is not None:
            headers["Content-Length"] = str(len(body))

        conn = http.client.HTTPConnection("127.0.0.1", viewer.spec.port, timeout=30)
        try:
            conn.request(method, upstream_path, body=body, headers=headers)
            response = conn.getresponse()
            payload = response.read()
            content_type = response.getheader("Content-Type", "")

            if "text/html" in content_type:
                text = payload.decode("utf-8", errors="replace")
                text = _rewrite_viewer_html(text, key)
                payload = text.encode("utf-8")

            self.send_response(response.status, response.reason)
            excluded = {"connection", "transfer-encoding", "content-length", "content-encoding"}
            for header, value in response.getheaders():
                if header.lower() not in excluded:
                    self.send_header(header, value)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
        except Exception as exc:
            self._send_json({"error": f"proxy error for {key}: {exc}"}, status=HTTPStatus.BAD_GATEWAY)
        finally:
            conn.close()


def _rewrite_viewer_html(text: str, app_key: str) -> str:
    """Make absolute in-app API calls work under /apps/<key>/ routing."""
    prefix = f"/apps/{app_key}"
    shim = f"""
<script>
(function() {{
  const prefix = {json.dumps(prefix)};
  const originalFetch = window.fetch.bind(window);
  window.fetch = function(resource, init) {{
    if (typeof resource === 'string' && resource.startsWith('/api')) {{
      resource = prefix + resource;
    }} else if (resource instanceof Request) {{
      const url = new URL(resource.url);
      if (url.origin === window.location.origin && url.pathname.startsWith('/api')) {{
        resource = new Request(prefix + url.pathname + url.search, resource);
      }}
    }}
    return originalFetch(resource, init);
  }};
}})();
</script>
"""
    if "</head>" in text:
        text = text.replace("</head>", shim + "\n</head>", 1)
    else:
        text = shim + text
    replacements = {
        "fetch('/api/": f"fetch('{prefix}/api/",
        'fetch("/api/': f'fetch("{prefix}/api/',
        "fetch(`/api/": f"fetch(`{prefix}/api/",
        "fetch('/api": f"fetch('{prefix}/api",
        'fetch("/api': f'fetch("{prefix}/api',
        "fetch(`/api": f"fetch(`{prefix}/api",
        "href=\"/api/": f"href=\"{prefix}/api/",
        "href='/api/": f"href='{prefix}/api/",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


INDEX_HTML = r"""
<!doctype html>
<html lang="pt-BR">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Genomics Workbench</title>
  <style>
    :root { color-scheme: dark; --bg:#0b1020; --panel:#111827; --card:#172033; --line:#2b3648; --text:#e5edf7; --muted:#aebbd0; --blue:#78c7ff; --green:#7ee787; --red:#ff6b6b; --yellow:#ffd166; }
    * { box-sizing: border-box; }
    body { margin: 0; background: var(--bg); color: var(--text); font-family: Inter, ui-sans-serif, system-ui, sans-serif; }
    header { padding: 16px 20px; background: linear-gradient(135deg, #080c18, #101b34); border-bottom: 1px solid var(--line); display:flex; align-items:center; justify-content:space-between; gap:16px; }
    h1 { margin:0; font-size:21px; letter-spacing:.2px; }
    .subtitle { margin:4px 0 0; color:var(--muted); font-size:13px; }
    button { border:1px solid #3d8bfd; background:#1f6feb; color:white; border-radius:8px; padding:8px 12px; cursor:pointer; font-weight:700; }
    main { display:grid; grid-template-columns: 360px 1fr; height: calc(100vh - 73px); }
    aside { border-right:1px solid var(--line); background:var(--panel); overflow:auto; padding:14px; }
    section { min-width:0; background:#050914; }
    .card { width:100%; text-align:left; margin:0 0 10px; padding:13px; border:1px solid var(--line); border-radius:14px; background:linear-gradient(180deg, #1a2438, #141d2e); color:var(--text); cursor:pointer; display:block; }
    .card:hover { border-color:#4c607c; }
    .card.active { border-color:var(--blue); box-shadow:0 0 0 1px rgba(120,199,255,.2) inset, 0 8px 28px rgba(0,0,0,.25); }
    .card.disabled { opacity:.55; cursor:not-allowed; }
    .card h2 { margin:0 0 5px; font-size:15px; }
    .card p { margin:0; color:var(--muted); font-size:12px; line-height:1.35; }
    .status { display:flex; align-items:center; gap:7px; margin-top:9px; font-size:12px; color:#cbd6e6; }
    .dot { width:9px; height:9px; border-radius:50%; background:var(--red); display:inline-block; }
    .dot.ok { background:var(--green); }
    .dot.warn { background:var(--yellow); }
    iframe { width:100%; height:100%; border:0; background:#0d1117; }
    .empty { padding:28px; color:var(--muted); }
    .topbar { height:44px; border-bottom:1px solid var(--line); background:#0d1324; display:flex; align-items:center; justify-content:space-between; padding:0 12px; color:var(--muted); font-size:13px; }
    .topbar a { color:var(--blue); text-decoration:none; }
    #frameWrap { height:calc(100% - 42px); }
    code { color:var(--yellow); font-size:11px; }
    @media (max-width: 900px) { main { grid-template-columns: 1fr; height:auto; } aside { border-right:0; border-bottom:1px solid var(--line); } section { height:75vh; } }
  </style>
</head>
<body>
  <header>
    <div>
      <h1>Genomics Workbench</h1>
      <p class="subtitle">Visualizacoes locais com roteamento unico, leitura sob demanda e execucao isolada por aplicativo.</p>
    </div>
    <button id="refresh">Atualizar status</button>
  </header>
  <main>
    <aside id="nav"></aside>
    <section>
      <div class="topbar"><span id="current">Selecione uma visualizacao</span><a id="open" target="_blank" href="#" style="display:none">abrir em nova aba</a></div>
      <div id="frameWrap"><div class="empty">As visualizacoes sao executadas como servidores locais independentes e exibidas aqui. Se uma delas aparecer desabilitada, veja o motivo no card e o log indicado.</div></div>
    </section>
  </main>
  <script>
    const nav = document.getElementById('nav');
    const frameWrap = document.getElementById('frameWrap');
    const current = document.getElementById('current');
    const openLink = document.getElementById('open');
    let activeKey = null;

    async function status() {
      const res = await fetch('/api/status');
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || res.statusText);
      return data.viewers || [];
    }

    function render(viewers) {
      nav.innerHTML = '';
      for (const v of viewers) {
        const card = document.createElement('div');
        card.className = 'card' + (v.key === activeKey ? ' active' : '') + (!v.enabled || !v.running ? ' disabled' : '');
        const dotClass = v.running ? 'ok' : (v.enabled ? 'warn' : '');
        const reason = v.running ? `ativo em ${v.url}` : (v.disabled_reason || `parado (returncode=${v.returncode})`);
        const internal = v.running ? `<p style="margin-top:6px">interno: <code>${v.internal_url}</code></p>` : '';
        card.innerHTML = `<h2>${v.title}</h2><p>${v.description}</p><div class="status"><span class="dot ${dotClass}"></span><span>${reason}</span></div>${internal}<p style="margin-top:6px">log: <code>${v.log_path}</code></p>`;
        if (v.enabled && v.running) {
          card.onclick = () => openViewer(v);
        }
        nav.appendChild(card);
      }
      if (!activeKey) {
        const first = viewers.find(v => v.enabled && v.running);
        if (first) openViewer(first);
      }
    }

    function openViewer(v) {
      activeKey = v.key;
      current.textContent = v.title;
      openLink.href = v.url;
      openLink.style.display = 'inline';
      frameWrap.innerHTML = '';
      const iframe = document.createElement('iframe');
      iframe.src = v.url;
      frameWrap.appendChild(iframe);
      status().then(render).catch(err => console.error(err));
    }

    async function refresh() {
      try { render(await status()); }
      catch (err) { nav.innerHTML = `<div class="empty">${err.message}</div>`; }
    }

    document.getElementById('refresh').onclick = refresh;
    refresh();
    setInterval(refresh, 10000);
  </script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch all genotype_based_predictor web visualizations")
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument("--runs-root", type=Path, default=DEFAULT_RUNS_ROOT)
    parser.add_argument("--aligned-tsv-root", type=Path, default=DEFAULT_ALIGNED_TSV_ROOT)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8780)
    parser.add_argument("--dataset-port", type=int, default=8770)
    parser.add_argument("--view-builder-port", type=int, default=8771)
    parser.add_argument("--experiment-port", type=int, default=8772)
    parser.add_argument("--track-port", type=int, default=8774)
    parser.add_argument("--alignment-port", type=int, default=8765)
    parser.add_argument("--log-dir", type=Path, default=None)
    args = parser.parse_args()

    if not _port_available(args.host, args.port):
        raise SystemExit(f"Porta principal em uso: {args.host}:{args.port}")

    log_dir = args.log_dir or Path(tempfile.gettempdir()) / "genomics_workbench_logs"
    workbench = Workbench(_build_specs(args), log_dir.resolve())
    workbench.start()

    class Handler(WorkbenchHandler):
        pass

    Handler.workbench = workbench
    server = ThreadingHTTPServer((args.host, args.port), Handler)

    def shutdown(_signum=None, _frame=None):
        workbench.stop()
        server.server_close()
        raise SystemExit(0)

    signal.signal(signal.SIGINT, shutdown)
    signal.signal(signal.SIGTERM, shutdown)

    print(f"Genomics Workbench: http://{args.host}:{args.port}")
    print(f"Logs: {log_dir.resolve()}")
    for item in workbench.payload()["viewers"]:
        state = "ON" if item["running"] else "OFF"
        reason = f" ({item['disabled_reason']})" if item["disabled_reason"] and not item["running"] else ""
        print(f"[{state}] {item['title']}: {item['url']}{reason}")

    try:
        server.serve_forever()
    finally:
        workbench.stop()


if __name__ == "__main__":
    main()
