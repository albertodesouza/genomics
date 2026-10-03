// Global app state shared by pages: active dataset, cohort filters, pinned individuals, locus.
import { api } from './api.js';

const listeners = new Set();
const MAX_PINNED = 32;

export const state = {
  status: null,
  datasets: [],
  datasetId: null,
  summary: null,
  samples: null, // { fields, columns, rows, index: Map(id -> row), col: {name -> idx} }
  filters: {}, // { field: [values] } — the cohort (Samples page)
  search: '',
  pinned: [], // sample ids used by Tracks / Sequence
  locus: {}, // per dataset: { gene, start, end, coords }
  genes: new Map(), // gene info cache
};

export function subscribe(fn) { listeners.add(fn); return () => listeners.delete(fn); }
export function emit(topic, detail) { for (const fn of listeners) fn(topic, detail); }

function storageKey(name) { return `gv.${state.datasetId || 'none'}.${name}`; }

function load(name, fallback) {
  try {
    const raw = localStorage.getItem(storageKey(name));
    return raw ? JSON.parse(raw) : fallback;
  } catch (e) { return fallback; }
}

function save(name, value) {
  try { localStorage.setItem(storageKey(name), JSON.stringify(value)); } catch (e) { /* private mode */ }
}

export const ds = () => `/api/d/${encodeURIComponent(state.datasetId)}`;

export async function loadStatus() {
  state.status = await api('/api/status');
  state.datasets = state.status.datasets;
  return state.status;
}

export async function setDataset(id) {
  state.datasetId = id;
  try { localStorage.setItem('gv.dataset', id); } catch (e) { /* ignore */ }
  state.summary = null;
  state.samples = null;
  state.genes = new Map();
  const [summary, samples] = await Promise.all([api(`${ds()}/summary`), api(`${ds()}/samples`)]);
  state.summary = summary;
  const col = {};
  samples.columns.forEach((name, i) => { col[name] = i; });
  const index = new Map();
  samples.rows.forEach((row) => index.set(row[col.sample_id], row));
  state.samples = { ...samples, col, index };
  state.filters = sanitizeFilters(load('filters', {}));
  state.search = '';
  state.pinned = (load('pinned', []) || []).filter((id) => index.has(id)).slice(0, MAX_PINNED);
  if (!state.pinned.length) state.pinned = samples.rows.slice(0, 3).map((r) => r[col.sample_id]);
  state.locus = load('locus', {}) || {};
  emit('dataset');
}

function sanitizeFilters(filters) {
  const out = {};
  const fields = new Set(state.samples.fields.filter((f) => f.kind === 'categorical').map((f) => f.name));
  for (const [k, v] of Object.entries(filters || {})) if (fields.has(k) && Array.isArray(v) && v.length) out[k] = v;
  return out;
}

export function setFilters(filters) {
  state.filters = sanitizeFilters(filters);
  save('filters', state.filters);
  emit('cohort');
}

export function setPinned(ids) {
  const seen = new Set();
  state.pinned = ids.filter((id) => state.samples.index.has(id) && !seen.has(id) && seen.add(id)).slice(0, MAX_PINNED);
  save('pinned', state.pinned);
  emit('pinned');
}

export function togglePinned(id) {
  if (state.pinned.includes(id)) setPinned(state.pinned.filter((x) => x !== id));
  else setPinned([...state.pinned, id]);
}

export function setLocus(patch) {
  state.locus = { ...state.locus, ...patch };
  save('locus', state.locus);
}

export function sampleValue(row, field) {
  const i = state.samples.col[field];
  return i === undefined ? undefined : row[i];
}

/** Sample rows matching filters (optionally ignoring one field, for facet counts). */
export function filterRows(filters = state.filters, search = state.search, ignoreField = null) {
  const { rows, col } = state.samples;
  const active = Object.entries(filters).filter(([k, v]) => k !== ignoreField && v && v.length).map(([k, v]) => [col[k], new Set(v.map(String))]);
  const needle = (search || '').trim().toLowerCase();
  return rows.filter((row) => {
    for (const [i, set] of active) if (!set.has(String(row[i] ?? ''))) return false;
    if (needle) {
      let hit = false;
      for (const v of row) if (v !== null && String(v).toLowerCase().includes(needle)) { hit = true; break; }
      if (!hit) return false;
    }
    return true;
  });
}

export function cohortSize() { return state.samples ? filterRows(state.filters, '').length : 0; }

export async function geneInfo(gene) {
  if (state.genes.has(gene)) return state.genes.get(gene);
  const info = await api(`${ds()}/genes/${encodeURIComponent(gene)}`);
  state.genes.set(gene, info);
  return info;
}

export function categoricalFields() {
  return state.samples ? state.samples.fields.filter((f) => f.kind === 'categorical') : [];
}

export const MAX_SERIES_COLORS = 8;
export function seriesColor(i) {
  const css = getComputedStyle(document.documentElement);
  return css.getPropertyValue(`--series-${(i % MAX_SERIES_COLORS) + 1}`).trim() || '#2a78d6';
}
