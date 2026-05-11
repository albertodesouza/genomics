from __future__ import annotations

import argparse
import html
import json
import sys
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, quote, urlparse


METRIC_KEYS = ("accuracy", "precision", "recall", "f1", "mse", "mae")
CHECKPOINT_FILES = ("best_accuracy.pt", "best_loss.pt", "final.pt")
DEFAULT_MAX_RESULT_BYTES = 2_000_000
DEFAULT_MAX_CONFIG_BYTES = 128_000


def _json_safe_float(value: Any) -> Optional[float]:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _read_text_preview(path: Path, max_bytes: int) -> Dict[str, Any]:
    size = path.stat().st_size
    with open(path, "rb") as f:
        raw = f.read(max_bytes + 1)
    truncated = len(raw) > max_bytes or size > max_bytes
    if truncated:
        raw = raw[:max_bytes]
    return {
        "path": str(path),
        "size_bytes": size,
        "truncated": truncated,
        "text": raw.decode("utf-8", errors="replace"),
    }


def _load_json_limited(path: Path, max_bytes: int) -> Any:
    size = path.stat().st_size
    if size > max_bytes:
        raise ValueError(f"result JSON is too large to load ({size} bytes > {max_bytes} bytes): {path.name}")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _extract_metrics(data: Any) -> Dict[str, float]:
    if not isinstance(data, dict):
        return {}
    metrics: Dict[str, float] = {}
    for key in METRIC_KEYS:
        value = _json_safe_float(data.get(key))
        if value is not None:
            metrics[key] = value
    return metrics


def _has_confusion_matrix(data: Any) -> bool:
    cm = data.get("confusion_matrix") if isinstance(data, dict) else None
    return isinstance(cm, list) and bool(cm) and all(isinstance(row, list) for row in cm)


def _sort_value(experiment: Dict[str, Any], sort_key: str) -> Tuple[int, float, str]:
    metrics = experiment.get("metrics", {})
    value = metrics.get(sort_key)
    if isinstance(value, (int, float)):
        # Larger is better for score metrics, smaller is better for error metrics.
        signed = float(value) if sort_key not in {"mse", "mae"} else -float(value)
        return (0, -signed, experiment["name"])
    return (1, 0.0, experiment["name"])


class ExperimentRepository:
    def __init__(self, runs_root: Path, max_result_bytes: int = DEFAULT_MAX_RESULT_BYTES):
        self.runs_root = Path(runs_root).expanduser().resolve()
        self.max_result_bytes = max_result_bytes
        if not self.runs_root.exists():
            raise FileNotFoundError(f"runs root not found: {self.runs_root}")
        if not self.runs_root.is_dir():
            raise NotADirectoryError(f"runs root is not a directory: {self.runs_root}")

    def _experiment_dirs(self) -> List[Path]:
        return sorted((p for p in self.runs_root.iterdir() if p.is_dir()), key=lambda p: p.name.lower())

    def _experiment_dir(self, name: str) -> Path:
        if not name or "/" in name or "\\" in name or name in {".", ".."}:
            raise ValueError("invalid experiment name")
        path = (self.runs_root / name).resolve()
        if path.parent != self.runs_root or not path.is_dir():
            raise FileNotFoundError(f"experiment not found: {name}")
        return path

    def _result_path(self, exp_dir: Path, filename: str) -> Path:
        if not filename or Path(filename).name != filename or not filename.endswith("_results.json"):
            raise ValueError("invalid result filename")
        path = (exp_dir / filename).resolve()
        if path.parent != exp_dir or not path.is_file():
            raise FileNotFoundError(f"result file not found: {filename}")
        return path

    def _summarize_result_file(self, path: Path, include_confusion: bool = False) -> Dict[str, Any]:
        item: Dict[str, Any] = {
            "file": path.name,
            "size_bytes": path.stat().st_size,
            "metrics": {},
            "has_confusion_matrix": False,
        }
        try:
            data = _load_json_limited(path, self.max_result_bytes)
            item["metrics"] = _extract_metrics(data)
            item["has_confusion_matrix"] = _has_confusion_matrix(data)
            if include_confusion and item["has_confusion_matrix"]:
                item["confusion_matrix"] = data.get("confusion_matrix")
        except Exception as exc:
            item["error"] = str(exc)
        return item

    def _summarize_experiment(self, exp_dir: Path, include_results: bool = False) -> Dict[str, Any]:
        config_path = exp_dir / "config.yaml"
        models_dir = exp_dir / "models"
        model_files = {name: (models_dir / name).is_file() for name in CHECKPOINT_FILES}
        result_paths = sorted(exp_dir.glob("*_results.json"), key=lambda p: p.name.lower())

        result_items = [self._summarize_result_file(path, include_confusion=include_results) for path in result_paths]
        metrics: Dict[str, float] = {}
        metric_sources: Dict[str, str] = {}
        for item in result_items:
            for key, value in item.get("metrics", {}).items():
                current = metrics.get(key)
                better = current is None or (value < current if key in {"mse", "mae"} else value > current)
                if better:
                    metrics[key] = value
                    metric_sources[key] = item["file"]

        summary: Dict[str, Any] = {
            "name": exp_dir.name,
            "path": str(exp_dir),
            "mtime": exp_dir.stat().st_mtime,
            "config_yaml": config_path.is_file(),
            "models": model_files,
            "result_files": [item["file"] for item in result_items],
            "result_count": len(result_items),
            "metrics": metrics,
            "metric_sources": metric_sources,
            "has_confusion_matrix": any(item.get("has_confusion_matrix") for item in result_items),
        }
        if include_results:
            summary["results"] = result_items
            if config_path.is_file():
                summary["config_preview"] = _read_text_preview(config_path, DEFAULT_MAX_CONFIG_BYTES)
        return summary

    def list_experiments(self, query: str = "", sort_key: str = "accuracy") -> List[Dict[str, Any]]:
        query = query.strip().lower()
        experiments = [self._summarize_experiment(path) for path in self._experiment_dirs()]
        if query:
            experiments = [
                exp for exp in experiments
                if query in exp["name"].lower()
                or any(query in filename.lower() for filename in exp.get("result_files", []))
            ]
        if sort_key in {*METRIC_KEYS, "name", "mtime", "result_count"}:
            if sort_key == "name":
                experiments.sort(key=lambda exp: exp["name"].lower())
            elif sort_key == "mtime":
                experiments.sort(key=lambda exp: exp["mtime"], reverse=True)
            elif sort_key == "result_count":
                experiments.sort(key=lambda exp: (exp["result_count"], exp["name"]), reverse=True)
            else:
                experiments.sort(key=lambda exp: _sort_value(exp, sort_key))
        return experiments

    def get_experiment(self, name: str) -> Dict[str, Any]:
        return self._summarize_experiment(self._experiment_dir(name), include_results=True)

    def get_result(self, name: str, filename: str) -> Dict[str, Any]:
        exp_dir = self._experiment_dir(name)
        path = self._result_path(exp_dir, filename)
        data = _load_json_limited(path, self.max_result_bytes)
        return {
            "name": name,
            "file": filename,
            "size_bytes": path.stat().st_size,
            "metrics": _extract_metrics(data),
            "has_confusion_matrix": _has_confusion_matrix(data),
            "data": data,
        }


class ExperimentDashboardHandler(BaseHTTPRequestHandler):
    repository: ExperimentRepository

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
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
        try:
            if parsed.path == "/":
                self._send_html(INDEX_HTML)
            elif parsed.path == "/api/experiments":
                self._handle_experiments(parsed.query)
            elif parsed.path == "/api/experiment":
                self._handle_experiment(parsed.query)
            elif parsed.path == "/api/result":
                self._handle_result(parsed.query)
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _handle_experiments(self, query: str) -> None:
        qs = parse_qs(query)
        q = qs.get("q", [""])[0]
        sort = qs.get("sort", ["accuracy"])[0]
        self._send_json({"experiments": self.repository.list_experiments(q, sort), "sort": sort, "q": q})

    def _handle_experiment(self, query: str) -> None:
        qs = parse_qs(query)
        name = qs.get("name", [""])[0]
        self._send_json(self.repository.get_experiment(name))

    def _handle_result(self, query: str) -> None:
        qs = parse_qs(query)
        name = qs.get("name", [""])[0]
        filename = qs.get("file", [""])[0]
        self._send_json(self.repository.get_result(name, filename))


INDEX_HTML = r"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Genomics Experiment Dashboard</title>
  <style>
    :root { color-scheme: dark; --bg: #0b1120; --panel: #111827; --panel2: #172033; --text: #e5e7eb; --muted: #94a3b8; --line: #283548; --blue: #60a5fa; --green: #34d399; --yellow: #fbbf24; --red: #fb7185; }
    * { box-sizing: border-box; }
    body { margin: 0; background: var(--bg); color: var(--text); font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
    body::before { content:""; position:fixed; inset:0; pointer-events:none; background: linear-gradient(180deg, rgba(96,165,250,.10), transparent 260px); z-index:-1; }
    header { padding: 1.25rem clamp(1rem, 3vw, 2rem); border-bottom: 1px solid var(--line); background: #0b1120; position: sticky; top: 0; z-index: 4; }
    h1 { margin: 0 0 .25rem; font-size: clamp(1.35rem, 3vw, 2.2rem); letter-spacing: -.03em; }
    .sub { color: var(--muted); font-size: .95rem; }
    main { display: grid; grid-template-columns: minmax(0, 1.35fr) minmax(320px, .8fr); gap: 1rem; padding: 1rem; max-width: 1500px; margin: 0 auto; }
    .card { background: #111827; border: 1px solid var(--line); border-radius: 16px; overflow: hidden; box-shadow: 0 20px 50px rgba(0, 0, 0, .22); }
    .toolbar { display: flex; flex-wrap: wrap; gap: .75rem; align-items: center; padding: .9rem; border-bottom: 1px solid var(--line); background: #172033; }
    input, select, button { background: #0f172a; color: var(--text); border: 1px solid #334155; border-radius: 10px; padding: .55rem .7rem; font: inherit; }
    input { min-width: min(26rem, 100%); flex: 1; }
    button { cursor: pointer; }
    button:hover, th:hover { border-color: var(--blue); color: white; }
    .status { color: var(--muted); font-size: .9rem; margin-left: auto; }
    .error { color: var(--red); }
    .table-wrap { overflow: auto; max-height: calc(100vh - 12rem); }
    table { width: 100%; border-collapse: collapse; min-width: 920px; }
    th, td { padding: .62rem .7rem; border-bottom: 1px solid var(--line); text-align: left; white-space: nowrap; }
    th { position: sticky; top: 0; background: #101827; z-index: 2; color: #cbd5e1; font-size: .8rem; text-transform: uppercase; letter-spacing: .04em; cursor: pointer; }
    tr { cursor: pointer; }
    tbody tr:hover, tbody tr.active { background: rgba(96, 165, 250, .12); }
    .pill { display: inline-block; min-width: 1.8rem; text-align: center; border-radius: 999px; padding: .12rem .45rem; font-size: .78rem; border: 1px solid var(--line); color: var(--muted); }
    .yes { color: var(--green); border-color: rgba(52, 211, 153, .45); }
    .no { color: var(--red); border-color: rgba(251, 113, 133, .38); }
    .metric { font-variant-numeric: tabular-nums; }
    .details { padding: 1rem; display: grid; gap: .9rem; }
    .details h2 { margin: 0; font-size: 1.05rem; }
    .grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: .55rem; }
    .kv { background: #0f172a; border: 1px solid var(--line); border-radius: 12px; padding: .65rem; }
    .kv b { display: block; font-size: .78rem; color: var(--muted); margin-bottom: .2rem; text-transform: uppercase; letter-spacing: .04em; }
    pre { margin: 0; padding: .8rem; overflow: auto; max-height: 25rem; background: #020617; border: 1px solid var(--line); border-radius: 12px; color: #dbeafe; font-size: .82rem; line-height: 1.4; }
    .cm { overflow: auto; }
    .cm table { min-width: auto; width: auto; }
    .cm td { text-align: right; font-variant-numeric: tabular-nums; }
    .empty { color: var(--muted); padding: 1rem; }
    @media (max-width: 980px) { main { grid-template-columns: 1fr; } .table-wrap { max-height: none; } }
  </style>
</head>
<body>
  <header>
    <h1>Genomics Experiment Dashboard</h1>
    <div class="sub">Local MVP for scanning immediate run directories, metrics JSON, config presence, and checkpoint presence without loading .pt files.</div>
  </header>
  <main>
    <section class="card">
      <div class="toolbar">
        <input id="q" placeholder="Filter experiments or result files" autocomplete="off">
        <select id="sort">
          <option value="accuracy">accuracy</option><option value="f1">f1</option><option value="precision">precision</option><option value="recall">recall</option><option value="mse">mse</option><option value="mae">mae</option><option value="mtime">mtime</option><option value="result_count">result_count</option><option value="name">name</option>
        </select>
        <button id="refresh">Refresh</button>
        <span id="status" class="status">Loading...</span>
      </div>
      <div class="table-wrap">
        <table>
          <thead><tr>
            <th data-sort="name">Experiment</th><th data-sort="accuracy">Acc</th><th data-sort="f1">F1</th><th data-sort="precision">Precision</th><th data-sort="recall">Recall</th><th data-sort="mse">MSE</th><th data-sort="mae">MAE</th><th>Config</th><th>best_acc</th><th>best_loss</th><th>final</th><th data-sort="result_count">Results</th>
          </tr></thead>
          <tbody id="rows"></tbody>
        </table>
      </div>
    </section>
    <aside class="card details" id="details"><div class="empty">Select an experiment to inspect result JSON previews and confusion matrices.</div></aside>
  </main>
<script>
const $ = (id) => document.getElementById(id);
let activeName = '';

function esc(value) {
  return String(value ?? '').replace(/[&<>"']/g, (c) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}
function fmt(value) {
  return typeof value === 'number' && Number.isFinite(value) ? value.toFixed(4) : '—';
}
function pill(value) {
  return `<span class="pill ${value ? 'yes' : 'no'}">${value ? 'yes' : 'no'}</span>`;
}
async function getJson(url) {
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok || data.error) throw new Error(data.error || res.statusText);
  return data;
}
async function loadExperiments() {
  const q = $('q').value;
  const sort = $('sort').value;
  $('status').textContent = 'Loading...';
  $('status').className = 'status';
  try {
    const data = await getJson(`/api/experiments?q=${encodeURIComponent(q)}&sort=${encodeURIComponent(sort)}`);
    $('rows').innerHTML = data.experiments.map((exp) => {
      const m = exp.metrics || {}, models = exp.models || {};
      return `<tr class="${exp.name === activeName ? 'active' : ''}" data-name="${esc(exp.name)}" onclick="loadExperiment(this.dataset.name)">
        <td title="${esc(exp.path)}">${esc(exp.name)}</td><td class="metric">${fmt(m.accuracy)}</td><td class="metric">${fmt(m.f1)}</td><td class="metric">${fmt(m.precision)}</td><td class="metric">${fmt(m.recall)}</td><td class="metric">${fmt(m.mse)}</td><td class="metric">${fmt(m.mae)}</td><td>${pill(exp.config_yaml)}</td><td>${pill(models['best_accuracy.pt'])}</td><td>${pill(models['best_loss.pt'])}</td><td>${pill(models['final.pt'])}</td><td>${esc(exp.result_count)}</td>
      </tr>`;
    }).join('') || '<tr><td colspan="12" class="empty">No experiments found.</td></tr>';
    $('status').textContent = `${data.experiments.length} experiments`;
  } catch (err) {
    $('status').textContent = err.message;
    $('status').className = 'status error';
  }
}
async function loadExperiment(name) {
  activeName = name;
  loadExperiments();
  $('details').innerHTML = '<div class="empty">Loading details...</div>';
  try {
    const exp = await getJson(`/api/experiment?name=${encodeURIComponent(name)}`);
    const metrics = Object.entries(exp.metrics || {}).map(([k, v]) => `<div class="kv"><b>${esc(k)}</b><span class="metric">${fmt(v)}</span><br><small>${esc((exp.metric_sources || {})[k] || '')}</small></div>`).join('') || '<div class="empty">No metrics found.</div>';
    const files = (exp.results || []).map((r) => `<option value="${esc(r.file)}">${esc(r.file)}${r.has_confusion_matrix ? ' [CM]' : ''}</option>`).join('');
    $('details').innerHTML = `<h2>${esc(exp.name)}</h2>
      <div class="sub">${esc(exp.path)}</div>
      <div class="grid">${metrics}</div>
      <div class="grid"><div class="kv"><b>config.yaml</b>${pill(exp.config_yaml)}</div><div class="kv"><b>result files</b>${esc(exp.result_count)}</div><div class="kv"><b>best_accuracy.pt</b>${pill(exp.models['best_accuracy.pt'])}</div><div class="kv"><b>best_loss.pt</b>${pill(exp.models['best_loss.pt'])}</div><div class="kv"><b>final.pt</b>${pill(exp.models['final.pt'])}</div><div class="kv"><b>confusion matrix</b>${pill(exp.has_confusion_matrix)}</div></div>
      <div><b>Result JSON preview</b><br><select id="resultFile">${files}</select> <button onclick="loadResult()">Open</button></div>
      <div id="resultPreview" class="empty">Choose a result file.</div>
      <div><b>config.yaml preview</b><pre>${esc(exp.config_preview ? exp.config_preview.text : 'config.yaml not found')}</pre></div>`;
    if (files) loadResult();
  } catch (err) {
    $('details').innerHTML = `<div class="empty error">${esc(err.message)}</div>`;
  }
}
async function loadResult() {
  const fileEl = $('resultFile');
  if (!fileEl || !fileEl.value) return;
  $('resultPreview').innerHTML = '<div class="empty">Loading result...</div>';
  try {
    const data = await getJson(`/api/result?name=${encodeURIComponent(activeName)}&file=${encodeURIComponent(fileEl.value)}`);
    const cm = data.has_confusion_matrix ? renderConfusionMatrix(data.data.confusion_matrix) : '<div class="empty">No confusion matrix in this result.</div>';
    $('resultPreview').innerHTML = `<div class="cm">${cm}</div><pre>${esc(JSON.stringify(data.data, null, 2))}</pre>`;
  } catch (err) {
    $('resultPreview').innerHTML = `<div class="empty error">${esc(err.message)}</div>`;
  }
}
function renderConfusionMatrix(cm) {
  const max = Math.max(1, ...cm.flat().map((x) => Number(x) || 0));
  const rows = cm.map((row, i) => `<tr><th>${i}</th>${row.map((value) => {
    const n = Number(value) || 0, alpha = 0.12 + 0.55 * (n / max);
    return `<td style="background: rgba(96,165,250,${alpha})">${esc(value)}</td>`;
  }).join('')}</tr>`).join('');
  return `<table><thead><tr><th></th>${cm.map((_, i) => `<th>${i}</th>`).join('')}</tr></thead><tbody>${rows}</tbody></table>`;
}
$('refresh').onclick = loadExperiments;
$('q').addEventListener('input', () => { clearTimeout(window.__qTimer); window.__qTimer = setTimeout(loadExperiments, 180); });
$('sort').onchange = loadExperiments;
document.querySelectorAll('th[data-sort]').forEach((th) => th.addEventListener('click', () => { $('sort').value = th.dataset.sort; loadExperiments(); }));
loadExperiments();
</script>
</body>
</html>
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Local dashboard for genotype experiment runs.")
    parser.add_argument("runs_root", type=Path, help="Directory containing immediate child experiment run directories.")
    parser.add_argument("--host", default="127.0.0.1", help="Host interface to bind (default: 127.0.0.1).")
    parser.add_argument("--port", type=int, default=8772, help="Port to bind (default: 8772).")
    parser.add_argument("--max-result-bytes", type=int, default=DEFAULT_MAX_RESULT_BYTES, help="Maximum result JSON size to load.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        repository = ExperimentRepository(args.runs_root, max_result_bytes=args.max_result_bytes)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    class Handler(ExperimentDashboardHandler):
        pass

    Handler.repository = repository
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    url = f"http://{args.host}:{args.port}/"
    print(f"Serving experiment dashboard at {url}")
    print(f"Runs root: {html.escape(str(repository.runs_root))}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping experiment dashboard.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
