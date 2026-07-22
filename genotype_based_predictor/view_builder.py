from __future__ import annotations

import argparse
import json
import re
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlparse


DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8771
DEFAULT_WINDOW_CENTER_SIZE = 32768
DEFAULT_DOWNSAMPLE_FACTOR = 1
DEFAULT_NORMALIZATION_METHOD = "log"
NORMALIZATION_METHODS = ["zscore", "minmax_keep_zero", "log"]


def _json_response(handler: BaseHTTPRequestHandler, payload: object, status: int = 200) -> None:
    body = json.dumps(payload, indent=2).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def _text_response(
    handler: BaseHTTPRequestHandler,
    text: str,
    content_type: str = "text/html; charset=utf-8",
    status: int = 200,
) -> None:
    body = text.encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", content_type)
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def _read_json_body(handler: BaseHTTPRequestHandler) -> Dict[str, Any]:
    length = int(handler.headers.get("Content-Length", "0") or "0")
    if length <= 0:
        return {}
    raw = handler.rfile.read(length).decode("utf-8")
    payload = json.loads(raw)
    if not isinstance(payload, dict):
        raise ValueError("Payload JSON deve ser um objeto")
    return payload


def _split_csv(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        parts = value
    else:
        parts = str(value).split(",")
    return [str(part).strip() for part in parts if str(part).strip()]


def _safe_view_name(name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip()).strip("._-")
    if not clean:
        raise ValueError("Informe um nome de view valido")
    return clean


def _default_output_path(name: str) -> Path:
    return Path(__file__).resolve().parent / "views" / f"{_safe_view_name(name)}.view.json"


class MetadataRepository:
    def __init__(self, dataset_dir: Path):
        self.dataset_dir = Path(dataset_dir).expanduser().resolve()
        self.metadata_path = self.dataset_dir / "dataset_metadata.json"
        self.metadata = self._load_metadata()
        self.genes = [str(g) for g in self.metadata.get("genes", [])]
        self.individual_ids = [str(i) for i in self.metadata.get("individuals", [])]
        self.pedigree = self.metadata.get("individuals_pedigree", {}) or {}
        self.population_distribution = self.metadata.get("population_distribution", {}) or {}
        self.superpopulation_distribution = self.metadata.get("superpopulation_distribution", {}) or {}

    def _load_metadata(self) -> Dict[str, Any]:
        if not self.metadata_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json nao encontrado em {self.dataset_dir}")
        with open(self.metadata_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            raise ValueError(f"Metadados invalidos: {self.metadata_path}")
        return data

    def individuals_payload(self) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        for sample_id in self.individual_ids:
            meta = self.pedigree.get(sample_id, {}) or {}
            rows.append(
                {
                    "sample_id": sample_id,
                    "population": meta.get("population"),
                    "superpopulation": meta.get("superpopulation"),
                    "sex": meta.get("sex"),
                    "sex_label": meta.get("sex_label"),
                }
            )
        return rows

    def options(self) -> Dict[str, Any]:
        outputs = self.metadata.get("alphagenome_outputs", []) or []
        ontologies = self.metadata.get("ontologies", []) or []
        return {
            "dataset_dir": str(self.dataset_dir),
            "metadata_path": str(self.metadata_path),
            "dataset_name": self.metadata.get("dataset_name"),
            "total_individuals": len(self.individual_ids),
            "genes": self.genes,
            "individuals": self.individuals_payload(),
            "populations": sorted(self.population_distribution),
            "superpopulations": sorted(self.superpopulation_distribution),
            "population_distribution": self.population_distribution,
            "superpopulation_distribution": self.superpopulation_distribution,
            "alphagenome_outputs": [str(v).lower() for v in outputs] or ["rna_seq"],
            "ontologies": [str(v) for v in ontologies],
            "ontology_details": self.metadata.get("ontology_details", {}) or {},
            "normalization_methods": NORMALIZATION_METHODS,
            "defaults": {
                "name": "new_view",
                "description": "",
                "alphagenome_outputs": "rna_seq",
                "window_center_size": DEFAULT_WINDOW_CENTER_SIZE,
                "downsample_factor": DEFAULT_DOWNSAMPLE_FACTOR,
                "normalization_method": DEFAULT_NORMALIZATION_METHOD,
                "output_path": str(_default_output_path("new_view")),
            },
        }

    def build_preview(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        name = _safe_view_name(str(payload.get("name") or "new_view"))
        selected_genes = self._selected_genes(payload)
        selected_individuals = self._selected_individuals(payload)
        alphagenome_outputs = _split_csv(payload.get("alphagenome_outputs")) or ["rna_seq"]
        ontology_terms = _split_csv(payload.get("ontology_terms"))

        window_center_size = int(payload.get("window_center_size") or DEFAULT_WINDOW_CENTER_SIZE)
        downsample_factor = int(payload.get("downsample_factor") or DEFAULT_DOWNSAMPLE_FACTOR)
        normalization_method = str(payload.get("normalization_method") or DEFAULT_NORMALIZATION_METHOD)

        if window_center_size < 1:
            raise ValueError("window_center_size deve ser >= 1")
        if downsample_factor < 1:
            raise ValueError("downsample_factor deve ser >= 1")
        if normalization_method not in NORMALIZATION_METHODS:
            raise ValueError(f"normalization_method deve ser um de: {', '.join(NORMALIZATION_METHODS)}")

        output_path = str(payload.get("output_path") or _default_output_path(name))
        view = {
            "name": name,
            "description": str(payload.get("description") or ""),
            "dataset_dir": str(self.dataset_dir),
            "alphagenome_outputs": alphagenome_outputs,
            "haplotype_mode": "H1+H2",
            "tensor_layout": "haplotype_channels",
            "window_center_size": window_center_size,
            "downsample_factor": downsample_factor,
            "genes_to_use": selected_genes,
            "sample_ids": selected_individuals,
            "sample_ids_path": None,
            "superpopulations_to_use": _split_csv(payload.get("superpopulations")) or None,
            "populations_to_use": _split_csv(payload.get("populations")) or None,
            "normalization_method": normalization_method,
            "ontology_terms": ontology_terms or None,
        }

        yaml_snippet = f"dataset_input:\n  view_path: {Path(output_path).expanduser()}"
        return {
            "view": view,
            "counts": {
                "genes": len(selected_genes),
                "individuals": len(selected_individuals),
                "alphagenome_outputs": len(alphagenome_outputs),
                "ontology_terms": len(ontology_terms),
            },
            "filters": {
                "query": str(payload.get("individual_query") or "").strip(),
                "populations": _split_csv(payload.get("populations")),
                "superpopulations": _split_csv(payload.get("superpopulations")),
            },
            "output_path": output_path,
            "yaml_snippet": yaml_snippet,
        }

    def _selected_genes(self, payload: Dict[str, Any]) -> List[str]:
        selected = _split_csv(payload.get("genes"))
        if not selected:
            selected = list(self.genes)
        unknown = sorted(set(selected) - set(self.genes))
        if unknown:
            raise ValueError(f"Genes nao encontrados no dataset: {', '.join(unknown[:10])}")
        return selected

    def _selected_individuals(self, payload: Dict[str, Any]) -> List[str]:
        explicit = _split_csv(payload.get("sample_ids"))
        allowed = set(self.individual_ids)
        if payload.get("explicit_sample_ids"):
            unknown = sorted(set(explicit) - allowed)
            if unknown:
                raise ValueError(f"Individuos nao encontrados no dataset: {', '.join(unknown[:10])}")
            base = explicit
        else:
            base = list(self.individual_ids)

        populations = set(_split_csv(payload.get("populations")))
        superpopulations = set(_split_csv(payload.get("superpopulations")))
        query = str(payload.get("individual_query") or "").strip().lower()

        selected: List[str] = []
        for sample_id in base:
            meta = self.pedigree.get(sample_id, {}) or {}
            population = str(meta.get("population") or "")
            superpopulation = str(meta.get("superpopulation") or "")
            if populations and population not in populations:
                continue
            if superpopulations and superpopulation not in superpopulations:
                continue
            if query:
                haystack = f"{sample_id} {population} {superpopulation} {meta.get('sex_label') or ''}".lower()
                if query not in haystack:
                    continue
            selected.append(sample_id)
        return selected


class ViewBuilderHandler(BaseHTTPRequestHandler):
    repository: MetadataRepository

    def log_message(self, fmt: str, *args: object) -> None:
        return

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/":
                _text_response(self, INDEX_HTML)
            elif parsed.path == "/api/options":
                _json_response(self, self.repository.options())
            else:
                _json_response(self, {"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            _json_response(self, {"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        try:
            payload = _read_json_body(self)
            if parsed.path == "/api/preview":
                _json_response(self, self.repository.build_preview(payload))
            elif parsed.path == "/api/save":
                self._handle_save(payload)
            else:
                _json_response(self, {"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            _json_response(self, {"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _handle_save(self, payload: Dict[str, Any]) -> None:
        preview = self.repository.build_preview(payload)
        output_path = Path(preview["output_path"]).expanduser()
        if not output_path.is_absolute():
            output_path = Path.cwd() / output_path
        output_path = output_path.resolve()

        overwrite = bool(payload.get("overwrite"))
        if output_path.exists() and not overwrite:
            _json_response(
                self,
                {"error": f"Arquivo ja existe: {output_path}. Envie overwrite=true para sobrescrever."},
                status=HTTPStatus.CONFLICT,
            )
            return

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(preview["view"], f, indent=2)
            f.write("\n")

        yaml_snippet = f"dataset_input:\n  view_path: {output_path}"
        _json_response(
            self,
            {
                "saved": True,
                "path": str(output_path),
                "counts": preview["counts"],
                "yaml_snippet": yaml_snippet,
                "view": preview["view"],
            },
        )


INDEX_HTML = r"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Genomics View Builder</title>
  <style>
    :root { color-scheme: dark; --bg:#10151f; --panel:#172131; --muted:#93a4bd; --text:#eef4ff; --line:#26374f; --blue:#5aa9ff; --green:#4dd688; --yellow:#ffd166; --red:#ff6b6b; }
    * { box-sizing: border-box; }
    body { margin: 0; font-family: system-ui, -apple-system, Segoe UI, sans-serif; background: radial-gradient(circle at top left, #20324b 0, var(--bg) 42rem); color: var(--text); }
    main { width: min(1920px, calc(100vw - 28px)); margin: 0 auto; padding: 18px 14px 28px; }
    h1 { margin: 0 0 6px; font-size: clamp(28px, 5vw, 46px); letter-spacing: -0.04em; }
    p { color: var(--muted); }
    .grid { display: grid; grid-template-columns: minmax(720px, 1.25fr) minmax(520px, .95fr); gap: 18px; align-items: start; }
    .card { background: color-mix(in srgb, var(--panel) 94%, transparent); border: 1px solid var(--line); border-radius: 18px; padding: 18px; box-shadow: 0 18px 70px rgba(0,0,0,.28); }
    label { display: block; margin: 12px 0 6px; color: #c9d7ea; font-weight: 650; }
    input, textarea, select, button { width: 100%; border-radius: 12px; border: 1px solid var(--line); background: #0e1521; color: var(--text); padding: 10px 12px; font: inherit; }
    textarea { min-height: 74px; resize: vertical; }
    select[multiple] { min-height: 150px; }
    button { cursor: pointer; font-weight: 750; background: linear-gradient(135deg, #2773ca, #42b883); border: 0; margin-top: 12px; min-width: 0; white-space: nowrap; }
    button.secondary { background: #22324a; border: 1px solid var(--line); }
    button.danger { background: #5c2630; border: 1px solid #7d3440; }
    .row { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px; }
    .three { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px; }
    .hint { font-size: 13px; color: var(--muted); margin-top: 6px; }
    .status { padding: 10px 12px; border-radius: 12px; margin: 14px 0; border: 1px solid var(--line); color: var(--muted); white-space: pre-wrap; }
    .ok { color: var(--green); border-color: color-mix(in srgb, var(--green) 45%, var(--line)); }
    .err { color: var(--red); border-color: color-mix(in srgb, var(--red) 50%, var(--line)); }
    .chips { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
    .chip { border: 1px solid var(--line); background: #0d1724; color: #dce8fb; padding: 5px 8px; border-radius: 999px; font-size: 12px; }
    pre { overflow: auto; background: #07101b; border: 1px solid var(--line); border-radius: 14px; padding: 14px; line-height: 1.45; max-height: 620px; }
    .individuals { max-height: 220px; overflow: auto; border: 1px solid var(--line); border-radius: 12px; padding: 8px; background: #0e1521; }
    .individuals label { display: flex; gap: 8px; align-items: center; margin: 4px 0; font-weight: 500; font-size: 13px; color: #d7e2f2; }
    .individuals input { width: auto; }
    .filter-grid { display: grid; grid-template-columns: repeat(3, minmax(260px, 1fr)); gap: 12px; }
    .filter-panel { border: 1px solid var(--line); border-radius: 14px; background: #0e1521; padding: 10px; }
    .filter-panel h3 { margin: 0 0 8px; font-size: 13px; color: #dce8fb; display: flex; justify-content: space-between; gap: 8px; }
    .check-list { max-height: 220px; overflow: auto; margin-top: 8px; }
    .check-list label { display: flex; gap: 8px; align-items: center; margin: 4px 0; font-weight: 500; font-size: 13px; color: #d7e2f2; }
    .check-list input { width: auto; }
    .split-actions { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 6px; }
    .split-actions button { margin-top: 0; padding: 7px 6px; font-size: 12px; border-radius: 9px; overflow: hidden; text-overflow: ellipsis; }
    @media (max-width: 1100px) { .filter-grid { grid-template-columns: 1fr; } }
    @media (max-width: 900px) { main { padding: 14px; } .grid, .row, .three { grid-template-columns: 1fr; } }
  </style>
</head>
<body>
<main>
  <h1>Genomics View Builder</h1>
  <p>Create `.view.json` files from `dataset_metadata.json` without editing JSON by hand.</p>
  <div id="status" class="status">Loading dataset options...</div>
  <div class="grid">
    <section class="card">
      <div class="row">
        <div><label>View name</label><input id="name" value="new_view"></div>
        <div><label>Output path</label><input id="output_path"></div>
      </div>
      <label>Description</label><textarea id="description"></textarea>

      <div class="filter-panel">
        <h3><span>Genes</span><span id="geneCount"></span></h3>
        <input id="geneSearch" placeholder="Search gene">
        <div class="split-actions">
          <button class="secondary" onclick="setChecks('genes', true)">All</button>
          <button class="secondary" onclick="setChecks('genes', false)">None</button>
          <button class="secondary" onclick="preview()">Preview</button>
        </div>
        <div id="genes" class="check-list"></div>
      </div>
      <div class="hint">If no genes are selected, the API previews all genes.</div>

      <div class="filter-grid">
        <div class="filter-panel">
          <h3><span>Superpopulations</span><span id="superpopCount"></span></h3>
          <input id="superpopSearch" placeholder="Search superpopulation">
          <div class="split-actions">
            <button class="secondary" onclick="setChecks('superpopulations', true)">All</button>
            <button class="secondary" onclick="setChecks('superpopulations', false)">None</button>
            <button class="secondary" onclick="renderPopulations(); preview()">Apply</button>
          </div>
          <div id="superpopulations" class="check-list"></div>
        </div>
        <div class="filter-panel">
          <h3><span>Populations</span><span id="populationCount"></span></h3>
          <input id="populationSearch" placeholder="Search population">
          <div class="split-actions">
            <button class="secondary" onclick="setChecks('populations', true)">All</button>
            <button class="secondary" onclick="setChecks('populations', false)">None</button>
            <button class="secondary" onclick="renderIndividuals(); preview()">Apply</button>
          </div>
          <div id="populations" class="check-list"></div>
        </div>
        <div class="filter-panel">
          <h3><span>Individuals</span><span id="individualCount"></span></h3>
          <input id="individual_query" placeholder="Search sample, population, superpopulation, sex">
          <div class="split-actions">
            <button class="secondary" onclick="selectVisibleIndividuals(true)">Visible</button>
            <button class="secondary" onclick="selectVisibleIndividuals(false)">None</button>
            <button class="secondary" onclick="renderIndividuals(); preview()">Apply</button>
          </div>
          <div id="individuals" class="check-list"></div>
          <div class="hint" id="individualNote"></div>
        </div>
      </div>

      <div class="three">
        <div><label>AlphaGenome outputs</label><input id="alphagenome_outputs" value="rna_seq"></div>
        <div><label>window_center_size</label><input id="window_center_size" type="number" min="1" value="32768"></div>
        <div><label>downsample_factor</label><input id="downsample_factor" type="number" min="1" value="1"></div>
      </div>
      <div class="row">
        <div><label>normalization_method</label><select id="normalization_method"></select></div>
        <div><label>ontology_terms</label><input id="ontology_terms" placeholder="CL:0000346, CL:1000458"></div>
      </div>
      <label><input id="overwrite" type="checkbox" style="width:auto"> overwrite existing output file</label>
      <button onclick="saveView()">Save view JSON</button>
    </section>
    <aside class="card">
      <h2>Preview</h2>
      <div id="summary" class="chips"></div>
      <pre id="preview">{}</pre>
      <h2>YAML</h2>
      <pre id="yaml">dataset_input:\n  view_path: ...</pre>
    </aside>
  </div>
</main>
<script>
let options = null;
let visibleIndividuals = [];

async function fetchJson(url, init) {
  const res = await fetch(url, init);
  const data = await res.json();
  if (!res.ok || data.error) throw new Error(data.error || `HTTP ${res.status}`);
  return data;
}
function setStatus(text, cls='') {
  const el = document.getElementById('status');
  el.className = `status ${cls}`;
  el.textContent = text;
}
function checkedValues(id) { return Array.from(document.querySelectorAll(`#${id} input:checked`)).map(o => o.value); }
function csvValue(id) { return checkedValues(id).join(','); }
function input(id) { return document.getElementById(id).value; }
function selectAll(id) { Array.from(document.getElementById(id).options).forEach(o => o.selected = true); preview(); }
function clearSelect(id) { Array.from(document.getElementById(id).options).forEach(o => o.selected = false); preview(); }
function fillSelect(id, values, labels) {
  const el = document.getElementById(id);
  el.innerHTML = '';
  values.forEach(v => {
    const opt = document.createElement('option');
    opt.value = v;
    opt.textContent = labels && labels[v] != null ? `${v} (${labels[v]})` : v;
    el.appendChild(opt);
  });
}
function fillCheckList(id, rows, labels, searchId, countId) {
  const previous = new Set(checkedValues(id));
  const hadPrevious = document.querySelectorAll(`#${id} input`).length > 0;
  const needle = input(searchId).trim().toLowerCase();
  const filtered = rows.filter(v => !needle || String(v).toLowerCase().includes(needle));
  const el = document.getElementById(id);
  el.innerHTML = '';
  filtered.forEach(v => {
    const label = document.createElement('label');
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.value = v;
    cb.checked = hadPrevious ? previous.has(v) : true;
    label.appendChild(cb);
    label.append(labels && labels[v] != null ? `${v} (${labels[v]})` : v);
    el.appendChild(label);
  });
  document.getElementById(countId).textContent = `${checkedValues(id).length}/${filtered.length}`;
  document.querySelectorAll(`#${id} input`).forEach(cb => cb.addEventListener('change', () => {
    document.getElementById(countId).textContent = `${checkedValues(id).length}/${filtered.length}`;
    if (id === 'superpopulations') renderPopulations();
    else renderIndividuals();
  }));
}

function renderPopulations() {
  const previous = new Set(checkedValues('populations'));
  const hadPrevious = document.querySelectorAll('#populations input').length > 0;
  const supers = new Set(checkedValues('superpopulations'));
  const counts = new Map();
  for (const row of options.individuals) {
    if (supers.size && !supers.has(row.superpopulation || '')) continue;
    const pop = row.population || '';
    if (!pop) continue;
    counts.set(pop, (counts.get(pop) || 0) + 1);
  }
  const rows = Array.from(counts.entries()).sort((a, b) => a[0].localeCompare(b[0])).map(([id, count]) => ({id, count}));
  const needle = input('populationSearch').trim().toLowerCase();
  const filtered = rows.filter(row => !needle || row.id.toLowerCase().includes(needle));
  const el = document.getElementById('populations');
  el.innerHTML = '';
  filtered.forEach(row => {
    const label = document.createElement('label');
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.value = row.id;
    cb.checked = !hadPrevious || previous.has(row.id);
    label.appendChild(cb);
    label.append(`${row.id} (${row.count})`);
    el.appendChild(label);
  });
  document.getElementById('populationCount').textContent = `${checkedValues('populations').length}/${filtered.length}`;
  document.querySelectorAll('#populations input').forEach(cb => cb.addEventListener('change', renderIndividuals));
  renderIndividuals();
}
function setChecks(id, checked) {
  document.querySelectorAll(`#${id} input`).forEach(el => { el.checked = checked; });
  if (id === 'genes') document.getElementById('geneCount').textContent = `${checkedValues(id).length}/${document.querySelectorAll(`#${id} input`).length}`;
  if (id === 'superpopulations') document.getElementById('superpopCount').textContent = `${checkedValues(id).length}/${document.querySelectorAll(`#${id} input`).length}`;
  if (id === 'populations') document.getElementById('populationCount').textContent = `${checkedValues(id).length}/${document.querySelectorAll(`#${id} input`).length}`;
  if (id === 'superpopulations') renderPopulations();
  else renderIndividuals();
  preview();
}
function currentPayload() {
  const rendered = Array.from(document.querySelectorAll('#individuals input'));
  const checked = rendered.filter(el => el.checked).map(el => el.value);
  const explicitSampleIds = rendered.length > 0 && checked.length !== rendered.length;
  return {
    name: input('name'),
    description: input('description'),
    genes: csvValue('genes'),
    superpopulations: csvValue('superpopulations'),
    populations: csvValue('populations'),
    individual_query: input('individual_query'),
    sample_ids: checked.join(','),
    explicit_sample_ids: explicitSampleIds,
    alphagenome_outputs: input('alphagenome_outputs'),
    window_center_size: Number(input('window_center_size')),
    downsample_factor: Number(input('downsample_factor')),
    normalization_method: input('normalization_method'),
    ontology_terms: input('ontology_terms'),
    output_path: input('output_path'),
    overwrite: document.getElementById('overwrite').checked
  };
}
function matchesIndividual(row) {
  const pops = new Set(checkedValues('populations'));
  const supers = new Set(checkedValues('superpopulations'));
  const q = input('individual_query').trim().toLowerCase();
  if (pops.size && !pops.has(row.population || '')) return false;
  if (supers.size && !supers.has(row.superpopulation || '')) return false;
  if (q) {
    const text = `${row.sample_id} ${row.population || ''} ${row.superpopulation || ''} ${row.sex_label || ''}`.toLowerCase();
    if (!text.includes(q)) return false;
  }
  return true;
}
function renderIndividuals() {
  const wrap = document.getElementById('individuals');
  const rendered = Array.from(document.querySelectorAll('#individuals input'));
  const previous = new Set(rendered.filter(el => el.checked).map(el => el.value));
  const hadRendered = rendered.length > 0;
  visibleIndividuals = options.individuals.filter(matchesIndividual);
  wrap.innerHTML = '';
  visibleIndividuals.slice(0, 500).forEach(row => {
    const label = document.createElement('label');
    const cb = document.createElement('input');
    cb.type = 'checkbox'; cb.value = row.sample_id; cb.checked = hadRendered ? previous.has(row.sample_id) : true;
    label.appendChild(cb);
    label.append(`${row.sample_id} | ${row.superpopulation || '-'} / ${row.population || '-'} | ${row.sex_label || '-'}`);
    wrap.appendChild(label);
  });
  const note = document.createElement('div');
  document.getElementById('individualCount').textContent = `${Array.from(document.querySelectorAll('#individuals input:checked')).length}/${Math.min(visibleIndividuals.length, 500)}`;
  document.getElementById('individualNote').textContent = `Showing ${Math.min(visibleIndividuals.length, 500)} of ${visibleIndividuals.length} filtered individuals. Preview/save uses checked visible IDs; narrow search to select beyond first 500.`;
}
function selectVisibleIndividuals(checked) {
  document.querySelectorAll('#individuals input').forEach(el => el.checked = checked);
  preview();
}
async function preview() {
  try {
    renderIndividuals();
    const data = await fetchJson('/api/preview', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(currentPayload())});
    document.getElementById('preview').textContent = JSON.stringify(data.view, null, 2);
    document.getElementById('yaml').textContent = data.yaml_snippet;
    document.getElementById('summary').innerHTML = Object.entries(data.counts).map(([k,v]) => `<span class="chip">${k}: ${v}</span>`).join('');
    setStatus(`Preview ready for ${data.counts.genes} genes and ${data.counts.individuals} individuals.`, 'ok');
  } catch (err) { setStatus(err.message, 'err'); }
}
async function saveView() {
  try {
    const data = await fetchJson('/api/save', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(currentPayload())});
    document.getElementById('preview').textContent = JSON.stringify(data.view, null, 2);
    document.getElementById('yaml').textContent = data.yaml_snippet;
    setStatus(`Saved ${data.path}\n${data.yaml_snippet}`, 'ok');
  } catch (err) { setStatus(err.message, 'err'); }
}
async function init() {
  try {
    options = await fetchJson('/api/options');
    fillCheckList('genes', options.genes, null, 'geneSearch', 'geneCount');
    fillCheckList('superpopulations', options.superpopulations, options.superpopulation_distribution, 'superpopSearch', 'superpopCount');
    renderPopulations();
    fillSelect('normalization_method', options.normalization_methods);
    document.getElementById('normalization_method').value = options.defaults.normalization_method;
    document.getElementById('output_path').value = options.defaults.output_path;
    document.getElementById('alphagenome_outputs').value = options.defaults.alphagenome_outputs;
    renderIndividuals();
    setStatus(`Loaded ${options.genes.length} genes and ${options.total_individuals} individuals from ${options.metadata_path}.`, 'ok');
    preview();
  } catch (err) { setStatus(err.message, 'err'); }
}
['name','description','individual_query','alphagenome_outputs','window_center_size','downsample_factor','normalization_method','ontology_terms','output_path'].forEach(id => {
  window.addEventListener('load', () => document.getElementById(id).addEventListener('change', preview));
});
window.addEventListener('load', () => {
  document.getElementById('geneSearch').addEventListener('input', () => fillCheckList('genes', options.genes, null, 'geneSearch', 'geneCount'));
  document.getElementById('superpopSearch').addEventListener('input', () => fillCheckList('superpopulations', options.superpopulations, options.superpopulation_distribution, 'superpopSearch', 'superpopCount'));
  document.getElementById('populationSearch').addEventListener('input', renderPopulations);
  document.getElementById('individual_query').addEventListener('input', () => { renderIndividuals(); preview(); });
});
window.addEventListener('load', init);
</script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Local web builder for genotype predictor view JSON files.")
    parser.add_argument("dataset_dir", help="Dataset directory containing dataset_metadata.json")
    parser.add_argument("--host", default=DEFAULT_HOST, help=f"Bind host (default: {DEFAULT_HOST})")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help=f"Bind port (default: {DEFAULT_PORT})")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repository = MetadataRepository(Path(args.dataset_dir))
    handler = type("ConfiguredViewBuilderHandler", (ViewBuilderHandler,), {"repository": repository})
    server = ThreadingHTTPServer((args.host, args.port), handler)
    url = f"http://{args.host}:{args.port}/"
    print(f"View Builder serving {repository.dataset_dir}")
    print(f"Open {url}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping View Builder")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
