from __future__ import annotations

import argparse
import json
import math
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse


DEFAULT_PAGE_SIZE = 50
MAX_PAGE_SIZE = 500


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError(f"JSON raiz deve ser um objeto: {path}")
    return payload


def _safe_int(value: str, default: int, minimum: int, maximum: int) -> int:
    if value == "":
        return default
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"Valor inteiro invalido: {value}") from exc
    return max(minimum, min(parsed, maximum))


def _split_multi(values: List[str]) -> List[str]:
    result: List[str] = []
    for value in values:
        for part in str(value).split(","):
            part = part.strip()
            if part:
                result.append(part)
    return result


def _json_sanitize(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): _json_sanitize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_sanitize(item) for item in value]
    if isinstance(value, tuple):
        return [_json_sanitize(item) for item in value]
    return value


def _sort_distribution(distribution: Dict[str, Any]) -> List[Dict[str, Any]]:
    items: List[Tuple[str, int]] = []
    for key, value in distribution.items():
        try:
            count = int(value)
        except (TypeError, ValueError):
            count = 0
        items.append((str(key), count))
    items.sort(key=lambda item: (-item[1], item[0]))
    return [{"id": key, "count": count} for key, count in items]


class DatasetRepository:
    def __init__(self, dataset_dir: Path):
        self.dataset_dir = Path(dataset_dir).resolve()
        self.metadata_path = self.dataset_dir / "dataset_metadata.json"
        if not self.metadata_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json nao encontrado em {self.dataset_dir}")
        self.metadata = _load_json(self.metadata_path)
        self.individuals = self._build_individual_index()
        self.by_sample_id = {row["sample_id"]: row for row in self.individuals}

    def _build_individual_index(self) -> List[Dict[str, Any]]:
        raw_individuals = self.metadata.get("individuals", [])
        if not isinstance(raw_individuals, list):
            raise ValueError("Campo 'individuals' deve ser uma lista em dataset_metadata.json")

        pedigree = self.metadata.get("individuals_pedigree", {})
        if not isinstance(pedigree, dict):
            pedigree = {}

        rows: List[Dict[str, Any]] = []
        for idx, raw_sample_id in enumerate(raw_individuals):
            sample_id = str(raw_sample_id)
            ped = pedigree.get(sample_id, {})
            if not isinstance(ped, dict):
                ped = {}
            sex = ped.get("sex_label", ped.get("sex", ""))
            rows.append(
                {
                    "sample_id": sample_id,
                    "population": str(ped.get("population", "")),
                    "superpopulation": str(ped.get("superpopulation", "")),
                    "sex": sex,
                    "index": idx,
                }
            )
        return rows

    def summary(self) -> Dict[str, Any]:
        population_distribution = self.metadata.get("population_distribution", {})
        superpopulation_distribution = self.metadata.get("superpopulation_distribution", {})
        if not isinstance(population_distribution, dict):
            population_distribution = {}
        if not isinstance(superpopulation_distribution, dict):
            superpopulation_distribution = {}

        return {
            "dataset_dir": str(self.dataset_dir),
            "dataset_name": self.metadata.get("dataset_name"),
            "metadata_path": str(self.metadata_path),
            "total_individuals": len(self.individuals),
            "total_genes": len(self.metadata.get("genes", []) or []),
            "total_populations": len(population_distribution),
            "total_superpopulations": len(superpopulation_distribution),
            "population_distribution": _sort_distribution(population_distribution),
            "superpopulation_distribution": _sort_distribution(superpopulation_distribution),
            "sex_distribution": self.metadata.get("sex_distribution", {}),
            "window_size": self.metadata.get("window_size"),
            "window_types": self.metadata.get("window_types", []),
            "alphagenome_outputs": self.metadata.get("alphagenome_outputs", []),
            "ontologies": self.metadata.get("ontologies", []),
            "created_at": self.metadata.get("creation_date"),
            "last_updated": self.metadata.get("last_updated"),
        }

    def genes_payload(self) -> Dict[str, Any]:
        raw_genes = self.metadata.get("genes", [])
        if not isinstance(raw_genes, list):
            raw_genes = []
        window_catalog = self.metadata.get("window_catalog", {})
        if not isinstance(window_catalog, dict):
            window_catalog = {}
        gene_strands = self.metadata.get("gene_strands", {})
        if not isinstance(gene_strands, dict):
            gene_strands = {}

        genes: List[Dict[str, Any]] = []
        for raw_gene in raw_genes:
            gene = str(raw_gene)
            catalog = window_catalog.get(gene, {})
            if not isinstance(catalog, dict):
                catalog = {}
            genes.append(
                {
                    "gene": gene,
                    "type": catalog.get("type"),
                    "chromosome": catalog.get("chromosome"),
                    "start": catalog.get("start"),
                    "end": catalog.get("end"),
                    "window_size": catalog.get("window_size", self.metadata.get("window_size")),
                    "strand": gene_strands.get(gene),
                    "outputs": catalog.get("outputs", self.metadata.get("alphagenome_outputs", [])),
                    "ontologies": catalog.get("ontologies", self.metadata.get("ontologies", [])),
                }
            )
        return {"total": len(genes), "genes": genes}

    def filter_options_payload(self) -> Dict[str, Any]:
        return {
            "populations": _sort_distribution(self.metadata.get("population_distribution", {}) or {}),
            "superpopulations": _sort_distribution(self.metadata.get("superpopulation_distribution", {}) or {}),
            "individuals": self.individuals,
        }

    def individuals_payload(
        self,
        page: int,
        page_size: int,
        q: str,
        populations: List[str],
        superpopulations: List[str],
        samples: List[str],
    ) -> Dict[str, Any]:
        q_lower = q.strip().lower()
        population_set = {p for p in populations if p}
        superpopulation_set = {s for s in superpopulations if s}
        sample_set = {s for s in samples if s}

        rows = self.individuals
        if q_lower:
            rows = [row for row in rows if q_lower in row["sample_id"].lower()]
        if population_set:
            rows = [row for row in rows if row.get("population") in population_set]
        if superpopulation_set:
            rows = [row for row in rows if row.get("superpopulation") in superpopulation_set]
        if sample_set:
            rows = [row for row in rows if row.get("sample_id") in sample_set]

        total = len(rows)
        start = (page - 1) * page_size
        end = start + page_size
        return {
            "page": page,
            "page_size": page_size,
            "total": total,
            "total_pages": (total + page_size - 1) // page_size if total else 0,
            "individuals": rows[start:end],
        }

    def individual_payload(self, sample_id: str) -> Dict[str, Any]:
        sample_id = sample_id.strip()
        if not sample_id:
            raise ValueError("Parametro id e obrigatorio")
        base = self.by_sample_id.get(sample_id)
        if base is None:
            raise ValueError(f"Individuo nao encontrado: {sample_id}")

        pedigree = self.metadata.get("individuals_pedigree", {})
        if not isinstance(pedigree, dict):
            pedigree = {}
        pedigree_meta = pedigree.get(sample_id, {})
        if not isinstance(pedigree_meta, dict):
            pedigree_meta = {}

        metadata_path = self.dataset_dir / "individuals" / sample_id / "individual_metadata.json"
        individual_metadata: Dict[str, Any] = {}
        metadata_error: Optional[str] = None
        if metadata_path.exists():
            try:
                individual_metadata = _load_json(metadata_path)
            except Exception as exc:  # keep detail endpoint useful even with one bad file
                metadata_error = str(exc)

        windows = individual_metadata.get("windows", []) if isinstance(individual_metadata, dict) else []
        if not isinstance(windows, list):
            windows = []

        return {
            "sample_id": sample_id,
            "summary": base,
            "pedigree": pedigree_meta,
            "individual_metadata": individual_metadata,
            "metadata_path": str(metadata_path) if metadata_path.exists() else None,
            "metadata_error": metadata_error,
            "windows": [str(window) for window in windows],
            "window_count": len(windows),
        }


class DatasetBrowserHandler(BaseHTTPRequestHandler):
    repository: DatasetRepository

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(_json_sanitize(payload), ensure_ascii=False, indent=2, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, text: str, content_type: str = "text/html; charset=utf-8", status: int = 200) -> None:
        body = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/":
                self._send_text(INDEX_HTML)
            elif parsed.path == "/api/summary":
                self._send_json(self.repository.summary())
            elif parsed.path == "/api/genes":
                self._send_json(self.repository.genes_payload())
            elif parsed.path == "/api/filter-options":
                self._send_json(self.repository.filter_options_payload())
            elif parsed.path == "/api/individuals":
                self._handle_individuals(parsed.query)
            elif parsed.path == "/api/individual":
                self._handle_individual(parsed.query)
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _handle_individuals(self, query: str) -> None:
        qs = parse_qs(query)
        page = _safe_int(qs.get("page", ["1"])[0], 1, 1, 1_000_000)
        page_size = _safe_int(qs.get("page_size", [str(DEFAULT_PAGE_SIZE)])[0], DEFAULT_PAGE_SIZE, 1, MAX_PAGE_SIZE)
        q = qs.get("q", [""])[0]
        populations = _split_multi(qs.get("population", []))
        superpopulations = _split_multi(qs.get("superpopulation", []))
        samples = _split_multi(qs.get("samples", []))
        self._send_json(self.repository.individuals_payload(page, page_size, q, populations, superpopulations, samples))

    def _handle_individual(self, query: str) -> None:
        qs = parse_qs(query)
        sample_id = qs.get("id", [""])[0]
        self._send_json(self.repository.individual_payload(sample_id))


INDEX_HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Genomics Dataset Browser</title>
  <style>
    :root { color-scheme: dark; --bg:#081016; --panel:#111b25; --panel2:#162333; --line:#284157; --text:#e7f0f7; --muted:#8ca2b3; --blue:#5db4ff; --green:#6ee7a8; --yellow:#f6c85f; --red:#ff6b6b; }
    * { box-sizing: border-box; }
    body { margin:0; background:radial-gradient(circle at top left,#122438,#081016 45%); color:var(--text); font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
    header { padding:24px; border-bottom:1px solid var(--line); background:rgba(8,16,22,.88); position:sticky; top:0; z-index:3; backdrop-filter:blur(8px); }
    h1 { margin:0 0 6px; font-size:24px; letter-spacing:.02em; }
    .sub { color:var(--muted); font-size:13px; }
    main { padding:20px; display:grid; gap:18px; grid-template-columns:minmax(0,1.4fr) minmax(320px,.8fr); }
    section, .card { background:rgba(17,27,37,.94); border:1px solid var(--line); border-radius:16px; box-shadow:0 18px 50px rgba(0,0,0,.25); }
    section { padding:16px; }
    .cards { display:grid; grid-template-columns:repeat(4,minmax(120px,1fr)); gap:12px; grid-column:1/-1; }
    .card { padding:14px; }
    .label { color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.08em; }
    .value { font-size:28px; font-weight:750; margin-top:4px; color:var(--green); }
    .filters { display:grid; grid-template-columns:1fr; gap:12px; margin:12px 0; }
    .filter-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:12px; }
    .filter-panel { border:1px solid var(--line); border-radius:14px; background:#0d1823; padding:10px; min-height:0; }
    .filter-panel h3 { margin:0 0 8px; font-size:13px; color:#cfe0ef; display:flex; justify-content:space-between; gap:8px; }
    .check-list { max-height:210px; overflow:auto; margin-top:8px; padding-right:3px; }
    .check-list label { display:flex; gap:8px; align-items:center; padding:4px 2px; color:#d8e7f3; font-size:13px; }
    .check-list input { width:auto; padding:0; }
    .filter-actions { display:flex; gap:8px; flex-wrap:wrap; margin-top:8px; }
    .filter-note { color:var(--muted); font-size:12px; margin-top:6px; }
    input, select, button { background:var(--panel2); color:var(--text); border:1px solid var(--line); border-radius:10px; padding:10px 12px; font:inherit; }
    button { cursor:pointer; color:white; background:linear-gradient(135deg,#1e6fb0,#264db7); border-color:#367ac2; }
    button.secondary { background:#172536; color:var(--blue); }
    table { width:100%; border-collapse:collapse; overflow:hidden; border-radius:12px; }
    th, td { text-align:left; padding:10px 11px; border-bottom:1px solid #203447; font-size:13px; }
    th { color:#bad0df; background:#132132; position:sticky; top:82px; }
    tr:hover td { background:#14263a; }
    a { color:var(--blue); text-decoration:none; }
    .row-actions { display:flex; gap:8px; align-items:center; justify-content:space-between; margin-top:12px; color:var(--muted); }
    .genes { display:grid; grid-template-columns:repeat(auto-fill,minmax(120px,1fr)); gap:8px; max-height:330px; overflow:auto; }
    .gene { padding:9px 10px; border:1px solid var(--line); border-radius:10px; background:#0d1823; }
    .gene b { display:block; color:var(--yellow); }
    pre { white-space:pre-wrap; overflow:auto; max-height:430px; background:#071018; border:1px solid var(--line); border-radius:12px; padding:12px; color:#d9e7f2; }
    .status { min-height:20px; color:var(--muted); margin:8px 0; }
    .error { color:var(--red); }
    .pill { display:inline-block; border:1px solid var(--line); border-radius:999px; padding:3px 8px; margin:2px; color:#c8d7e3; background:#102032; font-size:12px; }
    @media (max-width: 1100px) { .filter-grid { grid-template-columns:1fr; } }
    @media (max-width: 900px) { main { grid-template-columns:1fr; padding:12px; } .cards { grid-template-columns:repeat(2,1fr); } th { position:static; } }
  </style>
</head>
<body>
  <header>
    <h1>Genomics Dataset Browser</h1>
    <div class="sub" id="datasetPath">Loading dataset metadata...</div>
  </header>
  <main>
    <div class="cards">
      <div class="card"><div class="label">Individuals</div><div class="value" id="cardIndividuals">-</div></div>
      <div class="card"><div class="label">Genes</div><div class="value" id="cardGenes">-</div></div>
      <div class="card"><div class="label">Populations</div><div class="value" id="cardPopulations">-</div></div>
      <div class="card"><div class="label">Superpops</div><div class="value" id="cardSuperpops">-</div></div>
    </div>

    <section>
      <h2>Individuals</h2>
      <div class="filters">
        <input id="q" placeholder="Free search in selected individuals, e.g. HG00096">
        <div class="filter-grid">
          <div class="filter-panel">
            <h3><span>Superpopulations</span><span id="superpopCount"></span></h3>
            <input id="superpopSearch" placeholder="Search superpopulation">
            <div class="filter-actions"><button class="secondary" id="superpopAll">All</button><button class="secondary" id="superpopNone">None</button></div>
            <div class="check-list" id="superpopChecks"></div>
          </div>
          <div class="filter-panel">
            <h3><span>Populations</span><span id="populationCount"></span></h3>
            <input id="populationSearch" placeholder="Search population">
            <div class="filter-actions"><button class="secondary" id="populationAll">All</button><button class="secondary" id="populationNone">None</button></div>
            <div class="check-list" id="populationChecks"></div>
          </div>
          <div class="filter-panel">
            <h3><span>Individuals</span><span id="sampleCount"></span></h3>
            <input id="sampleSearch" placeholder="Search sample id">
            <div class="filter-actions"><button class="secondary" id="sampleAllVisible">All visible</button><button class="secondary" id="sampleNone">None</button></div>
            <div class="check-list" id="sampleChecks"></div>
            <div class="filter-note" id="sampleNote"></div>
          </div>
        </div>
        <button id="applyFilters">Apply selected filters</button>
      </div>
      <div id="status" class="status"></div>
      <table>
        <thead><tr><th>Sample</th><th>Population</th><th>Superpopulation</th><th>Sex</th><th></th></tr></thead>
        <tbody id="individualRows"></tbody>
      </table>
      <div class="row-actions">
        <button class="secondary" id="prevPage">Previous</button>
        <span id="pageInfo">Page -</span>
        <button class="secondary" id="nextPage">Next</button>
      </div>
    </section>

    <section>
      <h2>Selected Individual</h2>
      <div id="detailStatus" class="status">Select a row to inspect metadata and windows.</div>
      <div id="windows"></div>
      <pre id="details">{}</pre>
    </section>

    <section>
      <h2>Genes</h2>
      <div class="genes" id="genes"></div>
    </section>

    <section>
      <h2>Distributions</h2>
      <div class="label">Population</div>
      <div id="populationDist"></div>
      <br>
      <div class="label">Superpopulation</div>
      <div id="superpopulationDist"></div>
    </section>
  </main>
  <script>
    let state = { page: 1, pageSize: 50, totalPages: 0 };
    let allIndividuals = [];
    let populationRows = [];
    let superpopulationRows = [];

    function apiUrl(url) {
      if (!url.startsWith('/api')) return url;
      const parts = window.location.pathname.split('/').filter(Boolean);
      if (parts[0] === 'apps' && parts[1]) return `/apps/${parts[1]}${url}`;
      return url;
    }

    async function fetchJson(url) {
      const response = await fetch(apiUrl(url));
      const data = await response.json();
      if (!response.ok || data.error) throw new Error(data.error || response.statusText);
      return data;
    }

    function setStatus(message, isError=false) {
      const node = document.getElementById('status');
      node.textContent = message;
      node.className = isError ? 'status error' : 'status';
    }

    function renderDist(id, rows) {
      document.getElementById(id).innerHTML = rows.map(r => `<span class="pill">${r.id}: ${r.count}</span>`).join('');
    }

    function selectedValues(containerId) {
      return Array.from(document.querySelectorAll(`#${containerId} input:checked`)).map(el => el.value);
    }

    function setChecks(containerId, checked) {
      document.querySelectorAll(`#${containerId} input`).forEach(el => { el.checked = checked; });
    }

    function renderCheckList(containerId, rows, searchId, countId, formatter, checkedByDefault=true) {
      const previous = new Set(selectedValues(containerId));
      const hadPrevious = document.querySelectorAll(`#${containerId} input`).length > 0;
      const needle = document.getElementById(searchId).value.trim().toLowerCase();
      const filtered = rows.filter(row => !needle || String(row.id).toLowerCase().includes(needle));
      document.getElementById(containerId).innerHTML = filtered.map(row => {
        const checked = (!hadPrevious ? checkedByDefault : previous.has(row.id)) ? 'checked' : '';
        return `<label><input type="checkbox" value="${row.id}" ${checked}>${formatter(row)}</label>`;
      }).join('');
      document.getElementById(countId).textContent = `${selectedValues(containerId).length}/${filtered.length}`;
      document.querySelectorAll(`#${containerId} input`).forEach(el => el.addEventListener('change', () => {
        document.getElementById(countId).textContent = `${selectedValues(containerId).length}/${filtered.length}`;
      }));
    }

    function renderSampleChecks() {
      const previous = new Set(selectedValues('sampleChecks'));
      const hadPrevious = document.querySelectorAll('#sampleChecks input').length > 0;
      const needle = document.getElementById('sampleSearch').value.trim().toLowerCase();
      const selectedSuperpops = new Set(selectedValues('superpopChecks'));
      const selectedPops = new Set(selectedValues('populationChecks'));
      let rows = allIndividuals;
      if (selectedSuperpops.size) rows = rows.filter(row => selectedSuperpops.has(row.superpopulation));
      if (selectedPops.size) rows = rows.filter(row => selectedPops.has(row.population));
      if (needle) rows = rows.filter(row => row.sample_id.toLowerCase().includes(needle));
      const visible = rows.slice(0, 500);
      document.getElementById('sampleChecks').innerHTML = visible.map(row => {
        const checked = (!hadPrevious ? false : previous.has(row.sample_id)) ? 'checked' : '';
        return `<label><input type="checkbox" value="${row.sample_id}" ${checked}>${row.sample_id} | ${row.superpopulation || '-'} / ${row.population || '-'}</label>`;
      }).join('');
      document.getElementById('sampleCount').textContent = `${selectedValues('sampleChecks').length}/${visible.length}`;
      document.getElementById('sampleNote').textContent = `Showing ${visible.length} of ${rows.length} matching individuals. Leave all unchecked to use population/superpopulation filters only.`;
      document.querySelectorAll('#sampleChecks input').forEach(el => el.addEventListener('change', () => {
        document.getElementById('sampleCount').textContent = `${selectedValues('sampleChecks').length}/${visible.length}`;
      }));
    }

    async function loadSummary() {
      const s = await fetchJson('/api/summary');
      document.getElementById('datasetPath').textContent = `${s.dataset_name || 'dataset'} - ${s.dataset_dir}`;
      document.getElementById('cardIndividuals').textContent = s.total_individuals;
      document.getElementById('cardGenes').textContent = s.total_genes;
      document.getElementById('cardPopulations').textContent = s.total_populations;
      document.getElementById('cardSuperpops').textContent = s.total_superpopulations;
      renderDist('populationDist', s.population_distribution);
      renderDist('superpopulationDist', s.superpopulation_distribution);
    }

    async function loadGenes() {
      const data = await fetchJson('/api/genes');
      document.getElementById('genes').innerHTML = data.genes.map(g => `<div class="gene"><b>${g.gene}</b><span>${g.chromosome || ''}${g.start ? ':' + g.start + '-' + g.end : ''}</span></div>`).join('');
    }

    async function loadIndividuals() {
      const params = new URLSearchParams({
        page: state.page,
        page_size: state.pageSize,
        q: document.getElementById('q').value,
      });
      const populations = selectedValues('populationChecks');
      const superpopulations = selectedValues('superpopChecks');
      const samples = selectedValues('sampleChecks');
      if (populations.length) params.set('population', populations.join(','));
      if (superpopulations.length) params.set('superpopulation', superpopulations.join(','));
      if (samples.length) params.set('samples', samples.join(','));
      const data = await fetchJson('/api/individuals?' + params.toString());
      state.totalPages = data.total_pages;
      document.getElementById('individualRows').innerHTML = data.individuals.map(row => `
        <tr>
          <td><a href="#" data-id="${row.sample_id}">${row.sample_id}</a></td>
          <td>${row.population || ''}</td>
          <td>${row.superpopulation || ''}</td>
          <td>${row.sex || ''}</td>
          <td><button class="secondary" data-id="${row.sample_id}">Details</button></td>
        </tr>`).join('');
      document.getElementById('pageInfo').textContent = `Page ${data.page} / ${data.total_pages || 1} - ${data.total} matches`;
      setStatus(`Showing ${data.individuals.length} individuals. Page size is capped at ${state.pageSize}.`);
      document.querySelectorAll('[data-id]').forEach(el => el.addEventListener('click', ev => { ev.preventDefault(); loadIndividual(el.dataset.id); }));
    }

    async function loadAllIndividualsForFilters() {
      const data = await fetchJson('/api/filter-options');
      allIndividuals = data.individuals || [];
      populationRows = data.populations || [];
      superpopulationRows = data.superpopulations || [];
      renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`, true);
      renderCheckList('populationChecks', populationRows, 'populationSearch', 'populationCount', row => `${row.id} (${row.count})`, true);
      renderSampleChecks();
    }

    async function loadIndividual(sampleId) {
      try {
        document.getElementById('detailStatus').textContent = 'Loading ' + sampleId + '...';
        const data = await fetchJson('/api/individual?id=' + encodeURIComponent(sampleId));
        document.getElementById('detailStatus').textContent = `${sampleId}: ${data.window_count} windows`;
        document.getElementById('windows').innerHTML = data.windows.map(w => `<span class="pill">${w}</span>`).join('');
        document.getElementById('details').textContent = JSON.stringify(data, null, 2);
      } catch (err) {
        document.getElementById('detailStatus').textContent = `Error loading ${sampleId}: ${err.message || err}`;
        document.getElementById('detailStatus').className = 'status error';
      }
    }

    document.getElementById('applyFilters').addEventListener('click', () => { state.page = 1; loadIndividuals().catch(e => setStatus(e.message, true)); });
    document.getElementById('prevPage').addEventListener('click', () => { if (state.page > 1) { state.page--; loadIndividuals().catch(e => setStatus(e.message, true)); } });
    document.getElementById('nextPage').addEventListener('click', () => { if (!state.totalPages || state.page < state.totalPages) { state.page++; loadIndividuals().catch(e => setStatus(e.message, true)); } });
    document.getElementById('q').addEventListener('keydown', ev => { if (ev.key === 'Enter') { state.page = 1; loadIndividuals().catch(e => setStatus(e.message, true)); } });
    document.getElementById('superpopSearch').addEventListener('input', () => renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`, true));
    document.getElementById('populationSearch').addEventListener('input', () => renderCheckList('populationChecks', populationRows, 'populationSearch', 'populationCount', row => `${row.id} (${row.count})`, true));
    document.getElementById('sampleSearch').addEventListener('input', renderSampleChecks);
    document.getElementById('superpopAll').addEventListener('click', () => { setChecks('superpopChecks', true); renderSampleChecks(); });
    document.getElementById('superpopNone').addEventListener('click', () => { setChecks('superpopChecks', false); renderSampleChecks(); });
    document.getElementById('populationAll').addEventListener('click', () => { setChecks('populationChecks', true); renderSampleChecks(); });
    document.getElementById('populationNone').addEventListener('click', () => { setChecks('populationChecks', false); renderSampleChecks(); });
    document.getElementById('sampleAllVisible').addEventListener('click', () => setChecks('sampleChecks', true));
    document.getElementById('sampleNone').addEventListener('click', () => setChecks('sampleChecks', false));
    document.getElementById('superpopChecks').addEventListener('change', renderSampleChecks);
    document.getElementById('populationChecks').addEventListener('change', renderSampleChecks);

    Promise.all([loadSummary(), loadGenes(), loadAllIndividualsForFilters()]).then(loadIndividuals).catch(e => setStatus(e.message, true));
  </script>
</body>
</html>
"""


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Local web browser for genomics dataset metadata.")
    parser.add_argument("dataset_dir", type=Path, help="Dataset directory containing dataset_metadata.json")
    parser.add_argument("--host", default="127.0.0.1", help="Host/interface to bind")
    parser.add_argument("--port", type=int, default=8770, help="Port to bind")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    repository = DatasetRepository(args.dataset_dir)
    DatasetBrowserHandler.repository = repository
    server = ThreadingHTTPServer((args.host, args.port), DatasetBrowserHandler)
    print(f"Dataset browser: http://{args.host}:{args.port}/")
    print(f"Dataset: {repository.dataset_dir}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping dataset browser.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
