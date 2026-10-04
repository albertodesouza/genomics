// Pickers for AlphaGenome predictions (prediction form and import page): every output type, every
// tissue / cell type of the AlphaGenome track catalog, and the windows to predict (existing ones,
// any HGNC / GENCODE gene, genes of a Gene Ontology term or HGNC gene group, explicit regions and
// curated regulatory elements).
import { api } from './api.js';
import { h, clear, toast, checkbox, select, fmtInt, debounce } from './ui.js';
import { curieLink, dbsnpUrl, goTermLinks, hgncGroupUrl, xref } from './links.js';
import { geneChip } from './cards.js';

// Offered while the catalog is loading or unavailable (any CURIE can still be typed).
export const COMMON_ONTOLOGIES = [
  { curie: 'CL:1000458', name: 'melanocyte of skin' },
  { curie: 'CL:0000346', name: 'hair follicle dermal papilla cell' },
  { curie: 'CL:2000092', name: 'hair follicular keratinocyte' },
  { curie: 'CL:0000236', name: 'B cell' },
  { curie: 'CL:0000182', name: 'hepatocyte' },
  { curie: 'CL:0000540', name: 'neuron' },
  { curie: 'UBERON:0002107', name: 'liver' },
  { curie: 'UBERON:0000955', name: 'brain' },
  { curie: 'UBERON:0000948', name: 'heart' },
  { curie: 'UBERON:0002048', name: 'lung' },
];
const TERM_TYPES = [
  { value: '', label: 'All types' },
  { value: 'tissue', label: 'Tissues' },
  { value: 'primary_cell', label: 'Primary cells' },
  { value: 'cell_line', label: 'Cell lines' },
  { value: 'in_vitro_differentiated_cells', label: 'Differentiated cells' },
  { value: 'organoid', label: 'Organoids' },
];
const MAX_ROWS = 400;

export const fmtBytes = (n) => (n >= 1e12 ? `${(n / 1e12).toFixed(1)} TB` : n >= 1e9 ? `${(n / 1e9).toFixed(1)} GB` : n >= 1e6 ? `${(n / 1e6).toFixed(0)} MB` : `${Math.max(1, Math.round(n / 1e3))} kB`);

/** Labelled block that is not a <label> (pickers contain several inputs). */
export function section(label, control, extra = null) {
  return h('div', { class: 'field' }, h('div', { style: { display: 'flex', alignItems: 'baseline', gap: '8px' } }, h('span', { class: 'label' }, label), extra), control);
}

let catalogPromise = null;
/** The AlphaGenome track catalog (fetched once by the server from the selected backend). */
export function loadCatalog(refresh = false) {
  if (!catalogPromise || refresh) {
    catalogPromise = api('/api/alphagenome/catalog', { params: refresh ? { refresh: 1 } : {} }).catch((err) => { catalogPromise = null; throw err; });
  }
  return catalogPromise;
}

/** Uncompressed bytes of one haplotype window's output (mirrors outputs.stored_bytes). */
export function storedBytes(spec, length, tracks) {
  if (spec.kind === 'junctions') return Math.floor(length / 64) * (tracks * 4 + 17);
  const bins = Math.max(1, Math.floor(length / (spec.resolution || 1)));
  return (spec.kind === 'contact_map' ? bins * bins : bins) * tracks * 4;
}

function resolutionNote(spec) {
  if (spec.kind === 'junctions') return 'junction list';
  if (spec.kind === 'contact_map') return `2-D, ${spec.resolution} bp bins`;
  return spec.resolution > 1 ? `${spec.resolution} bp bins` : '';
}

// ---------------------------------------------------------------------------------- outputs
/** Output checkboxes grouped by assay family, with resolution and track counts. */
export function outputPicker(specs, selected, { existing = {}, onChange } = {}) {
  const chosen = new Set(selected);
  const counts = {};
  const host = h('div', { style: { display: 'grid', gap: '8px' } });
  const groups = {};
  for (const spec of Object.values(specs)) (groups[spec.group] = groups[spec.group] || []).push(spec);
  const changed = () => { render(); if (onChange) onChange(); };
  const render = () => {
    clear(host);
    for (const [group, list] of Object.entries(groups)) {
      host.appendChild(h('div', { style: { display: 'grid', gridTemplateColumns: '110px 1fr', gap: '8px', alignItems: 'start' } },
        h('span', { class: 'muted', style: { fontSize: '12px', paddingTop: '2px' } }, group),
        h('div', { class: 'chips', style: { gap: '4px 14px' } }, list.map((spec) => {
          const note = resolutionNote(spec);
          const disk = existing[spec.attr];
          const n = counts[spec.name];
          return checkbox(h('span', { title: spec.description },
            spec.label,
            note ? h('span', { class: 'muted' }, ` · ${note}`) : null,
            n !== undefined ? h('span', { class: 'muted' }, ` · ${fmtInt(n)} tracks`) : null,
            disk !== undefined ? h('span', { class: 'muted' }, ` (${disk} on disk)`) : null),
          chosen.has(spec.name), (v) => { if (v) chosen.add(spec.name); else chosen.delete(spec.name); changed(); });
        }))));
    }
  };
  render();
  const all = h('button', { class: 'btn small ghost', type: 'button', onclick: () => { Object.keys(specs).forEach((n) => chosen.add(n)); changed(); } }, 'All');
  const none = h('button', { class: 'btn small ghost', type: 'button', onclick: () => { chosen.clear(); changed(); } }, 'None');
  const el = section('Outputs', host, h('span', { style: { display: 'flex', gap: '4px' } }, all, none));
  el.values = () => Object.keys(specs).filter((n) => chosen.has(n));
  el.setCounts = (next) => { for (const k of Object.keys(counts)) delete counts[k]; Object.assign(counts, next || {}); render(); };
  return el;
}

// ------------------------------------------------------------------------- tissues / cells
/**
 * Tissue / cell-type picker over the full AlphaGenome catalog (falls back to the dataset's terms
 * and a few common ones while the catalog loads or when the backend is unavailable).
 * ``outputs()`` returns the chosen outputs: only terms with tracks for them are listed, with
 * per-output track counts.
 */
export function ontologyPicker(known, selected, { outputs = () => [], onChange } = {}) {
  const chosen = new Set(selected || []);
  const terms = new Map();
  const addTerm = (t) => { if (!terms.has(t.curie)) terms.set(t.curie, { outputs: {}, type: '', ...t }); else if (!terms.get(t.curie).name && t.name) terms.get(t.curie).name = t.name; };
  for (const t of [...(known || []), ...COMMON_ONTOLOGIES]) addTerm(t);
  let catalog = null;
  let catalogError = null;
  let query = '';
  let type = '';

  const status = h('div', { class: 'muted', style: { fontSize: '12px', display: 'flex', gap: '8px', alignItems: 'center', flexWrap: 'wrap' } });
  const summary = h('div', { class: 'chips' });
  const list = h('div', { class: 'checklist', style: { maxHeight: '240px', border: '1px solid var(--line)', borderRadius: '8px', padding: '4px' } });
  const search = h('input', { class: 'input', placeholder: 'Search tissues, cell types, CURIEs, TFs or marks…', style: { flex: '1', minWidth: '180px' } });
  search.addEventListener('input', debounce(() => { query = search.value.trim().toLowerCase(); renderList(); }, 120));
  const typeSelect = select(TERM_TYPES, type, (v) => { type = v; renderList(); });
  let shown = [];
  const selectAll = h('button', { class: 'btn small', type: 'button', title: 'Select every tissue / cell type listed (search and type filters apply)', onclick: () => { shown.forEach((t) => chosen.add(t.curie)); renderAll(); } }, 'Select all');
  const extra = h('input', { class: 'input mono', placeholder: 'Add CURIEs, e.g. UBERON:0002107, CL:0000236', style: { flex: '1' } });
  const addCuries = () => {
    for (const curie of extra.value.split(/[\s,;]+/).map((x) => x.trim()).filter(Boolean)) {
      if (!/^[A-Za-z]+:\d+$/.test(curie)) { toast(`Not a CURIE: ${curie}`, 'error'); continue; }
      if (catalog && !terms.has(curie)) toast(`${curie} has no AlphaGenome tracks`, 'error', 6000);
      addTerm({ curie, name: '' });
      chosen.add(curie);
    }
    extra.value = '';
    renderAll();
  };
  extra.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); addCuries(); } });

  const chosenOutputs = () => outputs() || [];
  const tissueOutputs = () => chosenOutputs().filter((o) => !catalog || !catalog.outputs[o] || catalog.outputs[o].tissue_specific !== false);
  const countsFor = (term) => {
    const outs = chosenOutputs();
    return outs.map((o) => [o, (term.outputs || {})[o] || 0]).filter(([, n]) => n > 0);
  };
  const matches = (term) => {
    if (type && term.type !== type) return false;
    if (catalog && tissueOutputs().length && !tissueOutputs().some((o) => (term.outputs || {})[o])) return false;
    if (!query) return true;
    const marks = Object.values(term.marks || {}).flat().join(' ');
    return `${term.name} ${term.curie} ${term.type} ${marks}`.toLowerCase().includes(query);
  };

  function renderStatus() {
    clear(status);
    if (catalog) {
      status.append(`AlphaGenome catalog: ${fmtInt(catalog.ontologies.length)} tissues / cell types, ${fmtInt(Object.values(catalog.outputs).reduce((a, o) => a + o.tracks, 0))} tracks`,
        h('button', { class: 'btn small ghost', type: 'button', title: `Fetched ${catalog.fetched_at} from ${catalog.source}`, onclick: () => fetchCatalog(true) }, 'Refresh'));
    } else if (catalogError) {
      status.append(h('span', null, `Full catalog unavailable (${catalogError}); showing known terms, CURIEs can be typed.`), h('button', { class: 'btn small ghost', type: 'button', onclick: () => fetchCatalog(true) }, 'Retry'));
    } else {
      status.append('Loading the AlphaGenome track catalog…');
    }
  }

  function renderSummary() {
    clear(summary);
    if (!chosen.size) { summary.appendChild(h('span', { class: 'muted', style: { fontSize: '12px' } }, 'None selected')); return; }
    for (const curie of chosen) {
      const t = terms.get(curie) || { curie, name: '' };
      summary.appendChild(h('span', { class: 'pill removable', title: curie }, t.name || curie,
        h('button', { type: 'button', title: 'Remove', onclick: () => { chosen.delete(curie); renderAll(); } }, '×')));
    }
    if (chosen.size > 1) summary.appendChild(h('button', { class: 'btn small ghost', type: 'button', onclick: () => { chosen.clear(); renderAll(); } }, `Clear ${fmtInt(chosen.size)}`));
  }

  function renderList() {
    clear(list);
    shown = [...terms.values()].filter(matches).sort((a, b) => (chosen.has(b.curie) - chosen.has(a.curie)) || (a.name || a.curie).localeCompare(b.name || b.curie));
    selectAll.disabled = !shown.length || shown.every((t) => chosen.has(t.curie));
    selectAll.textContent = shown.length > 1 ? `Select all ${fmtInt(shown.length)}` : 'Select all';
    for (const term of shown.slice(0, MAX_ROWS)) {
      const box = h('input', { type: 'checkbox' });
      box.checked = chosen.has(term.curie);
      box.addEventListener('change', () => { if (box.checked) chosen.add(term.curie); else chosen.delete(term.curie); renderSummary(); selectAll.disabled = shown.every((t) => chosen.has(t.curie)); if (onChange) onChange(); });
      const counts = countsFor(term).map(([o, n]) => `${o.toLowerCase()} ${n}`).join(' · ');
      list.appendChild(h('label', { title: Object.entries(term.marks || {}).map(([o, m]) => `${o}: ${m.join(', ')}`).join('\n') || term.curie },
        box, h('span', null, term.name || term.curie),
        term.type ? h('span', { class: 'muted', style: { fontSize: '11px' } }, term.type.replace(/_/g, ' ')) : null,
        h('span', { class: 'meta mono' }, counts ? `${counts} · ` : '', curieLink(term.curie))));
    }
    const outs = tissueOutputs();
    const scope = catalog && outs.length ? ` with ${outs.map((o) => o.toLowerCase()).join(' / ')} tracks` : '';
    const footer = shown.length > MAX_ROWS ? `Showing ${MAX_ROWS} of ${fmtInt(shown.length)}${scope}; refine the search (Select all still takes all ${fmtInt(shown.length)}).` : `${fmtInt(shown.length)} tissues / cell types${scope}`;
    list.appendChild(h('div', { class: 'muted', style: { fontSize: '11.5px', padding: '4px 6px' } }, footer));
  }

  function renderAll() { renderStatus(); renderSummary(); renderList(); if (onChange) onChange(); }

  async function fetchCatalog(refresh) {
    catalogError = null;
    renderStatus();
    try {
      catalog = await loadCatalog(refresh);
      for (const t of catalog.ontologies) {
        const existing = terms.get(t.curie);
        terms.set(t.curie, { ...t, name: t.name || (existing && existing.name) || '' });
      }
    } catch (err) { catalogError = err.message; }
    renderAll();
  }

  const el = h('div', { style: { display: 'grid', gap: '6px' } },
    status, summary,
    h('div', { style: { display: 'flex', gap: '6px', flexWrap: 'wrap', alignItems: 'center' } }, search, typeSelect, selectAll),
    list,
    h('div', { style: { display: 'flex', gap: '6px' } }, extra, h('button', { class: 'btn', type: 'button', onclick: addCuries }, 'Add')));
  el.values = () => [...chosen];
  el.refresh = () => { renderList(); };
  el.catalog = () => catalog;
  /** Tracks per chosen output for the current selection (null without the catalog). */
  el.trackCounts = () => {
    if (!catalog) return null;
    const out = {};
    for (const o of chosenOutputs()) {
      const info = catalog.outputs[o];
      if (!info) continue;
      if (info.tissue_specific === false) { out[o] = info.tracks; continue; }
      out[o] = [...chosen].reduce((sum, c) => sum + (((terms.get(c) || {}).outputs || {})[o] || 0), 0);
    }
    return out;
  };
  renderAll();
  fetchCatalog(false);
  return el;
}

// ---------------------------------------------------------------------------------- windows
/** Gene symbol / Ensembl id search box (GET /api/genes/search) with a suggestion list. */
export function geneSearch(datasetId, onPick, { placeholder = 'Add any gene: HGNC symbol, alias, previous symbol, name or Ensembl id (e.g. KITLG, OCA3, tyrosinase)' } = {}) {
  const input = h('input', { class: 'input', placeholder, autocomplete: 'off', style: { width: '100%' } });
  const menu = h('div', { class: 'suggest', hidden: true });
  let hits = [];
  let active = -1;
  let source = '';
  const pick = (gene) => { onPick(gene); input.value = ''; menu.hidden = true; hits = []; };
  const render = () => {
    clear(menu);
    menu.hidden = !hits.length;
    hits.forEach((g, i) => menu.appendChild(h('div', {
      class: `suggest-item${i === active ? ' on' : ''}`,
      title: g.hgnc_id ? `${g.hgnc_id} · ${g.locus_type || ''}` : 'GENCODE gene (no HGNC record)',
      onmousedown: (e) => { e.preventDefault(); pick(g); },
    }, h('b', null, g.symbol || g.name), g.symbol && g.symbol !== g.name ? h('span', { class: 'muted' }, ` (GENCODE ${g.name})`) : null,
    h('span', { class: 'muted' }, ` ${g.full_name || g.type.replace(/_/g, ' ')}`),
    g.matched && !g.matched.startsWith('ENSG') && g.matched !== g.full_name ? h('span', { class: 'pill' }, g.matched) : null,
    h('span', { class: 'meta mono' }, `${g.location ? `${g.location} · ` : ''}${g.chrom}:${fmtInt(g.start)}-${fmtInt(g.end)} ${g.strand} · ${g.id}`))));
    if (source) menu.appendChild(h('div', { class: 'muted', style: { fontSize: '11px', padding: '3px 8px' } }, `Source: ${source}`));
  };
  const lookup = debounce(async () => {
    const q = input.value.trim();
    if (q.length < 2) { hits = []; render(); return; }
    try {
      const res = await api('/api/genes/search', { params: { q, ds: datasetId || '' } });
      if (input.value.trim() !== q) return;
      hits = res.genes;
      source = res.source === 'HGNC + GENCODE' ? 'HGNC approved symbols, aliases and names · GENCODE coordinates' : 'GENCODE gene names (HGNC unavailable)';
      active = hits.length ? 0 : -1;
      render();
    } catch (err) { toast(err.message, 'error'); }
  }, 150);
  input.addEventListener('input', lookup);
  input.addEventListener('blur', () => { setTimeout(() => { menu.hidden = true; }, 100); });
  input.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
      e.preventDefault();
      if (hits.length) { active = (active + (e.key === 'ArrowDown' ? 1 : hits.length - 1)) % hits.length; render(); }
    } else if (e.key === 'Enter') {
      e.preventDefault();
      if (hits[active]) pick(hits[active]);
    } else if (e.key === 'Escape') { menu.hidden = true; }
  });
  return h('div', { style: { position: 'relative' } }, input, menu);
}

/**
 * Genes from a database: search Gene Ontology terms (QuickGO; genes annotated to the term or its
 * descendants) or HGNC gene groups, review the genes and add the chosen ones. ``onAdd`` receives
 * GENCODE gene records (as ``geneSearch`` hits); genes without GENCODE coordinates are listed but
 * cannot be added.
 */
export function geneSetPicker(datasetId, onAdd, { addLabel = 'Add' } = {}) {
  const input = h('input', { class: 'input', placeholder: 'Gene Ontology term or HGNC gene group, e.g. melanin biosynthetic process, GO:0042438, melanosome, Tyrosinase family', autocomplete: 'off', style: { width: '100%' } });
  const results = h('div', { class: 'geneset-results', hidden: true });
  const panel = h('div', { style: { display: 'grid', gap: '6px' } });
  let current = null;
  let regulation = false;
  const search = debounce(async () => {
    const q = input.value.trim();
    if (q.length < 3) { results.hidden = true; return; }
    try {
      const res = await api('/api/genesets/search', { params: { q, limit: 15 } });
      if (input.value.trim() !== q) return;
      clear(results);
      for (const t of res.go) {
        results.appendChild(h('div', { class: 'suggest-item', onclick: () => load('go', t.id) },
          h('span', { class: 'pill' }, 'GO'), h('b', null, t.name), h('span', { class: 'muted' }, ` ${String(t.aspect || '').replace(/_/g, ' ')}`),
          h('span', { class: 'meta' }, xref({ ...goTermLinks(t.id)[0], label: t.id }, { compact: true }))));
      }
      for (const g of res.hgnc_groups) {
        results.appendChild(h('div', { class: 'suggest-item', onclick: () => load('hgnc', g.name, g.id) },
          h('span', { class: 'pill' }, 'HGNC group'), h('b', null, g.name), h('span', { class: 'muted' }, ` ${fmtInt(g.size)} genes`),
          h('span', { class: 'meta' }, g.id ? xref({ label: `group ${g.id}`, url: hgncGroupUrl(g.id), title: 'HGNC gene group page' }, { compact: true }) : null)));
      }
      if (!res.go.length && !res.hgnc_groups.length) results.appendChild(h('div', { class: 'muted', style: { padding: '6px', fontSize: '12px' } }, 'No Gene Ontology term or HGNC group matches.'));
      for (const err of res.errors || []) results.appendChild(h('div', { class: 'muted', style: { padding: '6px', fontSize: '12px' } }, err));
      results.hidden = false;
    } catch (err) { toast(err.message, 'error'); }
  }, 250);
  input.addEventListener('input', search);

  async function load(source, id, groupId = '') {
    current = { source, id };
    results.hidden = true;
    clear(panel).appendChild(h('div', { class: 'muted', style: { fontSize: '12px' } }, h('span', { class: 'spinner', style: { marginRight: '6px' } }), source === 'go' ? 'Loading human genes annotated to the term (QuickGO)…' : 'Loading the HGNC gene group…'));
    let res;
    try { res = await api('/api/genesets/genes', { params: { source, id, ds: datasetId || '', regulation: regulation ? 1 : '' } }); } catch (err) { clear(panel).appendChild(h('div', { class: 'error-box' }, err.message)); return; }
    if (current.id !== id) return;
    const usable = res.genes.filter((g) => g.gencode);
    const chosen = new Set(usable.filter((g) => !g.in_dataset).map((g) => g.symbol));
    const list = h('div', { class: 'checklist', style: { maxHeight: '220px' } });
    const addBtn = h('button', { class: 'btn small primary', type: 'button' });
    const refresh = () => { addBtn.textContent = `${addLabel} ${chosen.size} gene${chosen.size === 1 ? '' : 's'}`; addBtn.disabled = !chosen.size; };
    const boxes = [];
    for (const g of res.genes) {
      const box = h('input', { type: 'checkbox', disabled: !g.gencode });
      box.checked = chosen.has(g.symbol);
      box.addEventListener('change', () => { if (box.checked) chosen.add(g.symbol); else chosen.delete(g.symbol); refresh(); });
      boxes.push([box, g]);
      const loc = g.gencode ? `${g.gencode.chrom}:${fmtInt(g.gencode.start)}` : 'no GENCODE locus';
      list.appendChild(h('label', { title: (g.terms || []).map((t) => `${t.id} ${t.name}`).join('\n') || g.hgnc_id },
        box, geneChip(g.symbol), h('span', { class: 'muted', style: { fontSize: '11.5px' } }, g.full_name),
        g.in_dataset ? h('span', { class: 'pill' }, 'in dataset') : null,
        h('span', { class: 'meta mono' }, loc)));
    }
    const setAll = (on) => { for (const [box, g] of boxes) { if (!g.gencode) continue; box.checked = on; if (on) chosen.add(g.symbol); else chosen.delete(g.symbol); } refresh(); };
    addBtn.addEventListener('click', () => {
      const genes = usable.filter((g) => chosen.has(g.symbol)).map((g) => ({ ...g.gencode, symbol: g.symbol, full_name: g.full_name, hgnc_id: g.hgnc_id }));
      onAdd(genes);
      toast(`Added ${genes.length} gene${genes.length === 1 ? '' : 's'} from ${res.title}`);
    });
    refresh();
    const term = res.term;
    clear(panel).append(
      h('div', { style: { display: 'flex', gap: '8px', alignItems: 'baseline', flexWrap: 'wrap' } },
        h('b', null, res.title), source === 'go' ? xref({ ...goTermLinks(id)[0], label: id }) : groupId ? xref({ label: `HGNC group ${groupId}`, url: hgncGroupUrl(groupId), title: 'HGNC gene group page' }) : null,
        h('span', { class: 'muted', style: { fontSize: '12px' } }, `${fmtInt(res.genes.length)} human genes${source === 'go' ? ' (reviewed UniProt entries annotated to the term or its descendants)' : ''}`)),
      term && term.definition ? h('div', { class: 'muted', style: { fontSize: '12px' } }, term.definition) : null,
      source === 'go' ? checkbox('Include genes that regulate it (regulates / positively / negatively regulates)', regulation, (v) => { regulation = v; load('go', id); }) : null,
      list,
      h('div', { style: { display: 'flex', gap: '6px', alignItems: 'center' } },
        h('button', { class: 'btn small ghost', type: 'button', onclick: () => setAll(true) }, 'All'),
        h('button', { class: 'btn small ghost', type: 'button', onclick: () => setAll(false) }, 'None'),
        addBtn,
        res.genes.length > usable.length ? h('span', { class: 'muted', style: { fontSize: '11.5px' } }, `${res.genes.length - usable.length} without a GENCODE locus in this gene table`) : null));
  }

  return h('div', { style: { display: 'grid', gap: '8px', paddingTop: '6px' } }, input, results, panel,
    h('p', { class: 'help', style: { margin: 0 } }, 'Gene Ontology annotations from QuickGO (EBI), gene groups from HGNC. Genes already in the dataset start unchecked.'));
}

/** Curated regulatory elements (non-gene windows) as a grouped checklist. */
export function presetChecklist(presets, isChosen, onToggle, { existing = new Set() } = {}) {
  const groups = {};
  for (const p of presets) (groups[p.category] = groups[p.category] || []).push(p);
  return h('div', { class: 'checklist', style: { maxHeight: '220px' } }, Object.entries(groups).map(([category, items]) => [
    h('div', { class: 'muted', style: { fontSize: '11px', fontWeight: 600, textTransform: 'uppercase', letterSpacing: '.04em', padding: '6px 6px 2px' } }, category),
    items.map((p) => {
      const box = h('input', { type: 'checkbox' });
      box.checked = isChosen(p);
      box.addEventListener('change', () => onToggle(p, box.checked));
      return h('label', { title: `${p.description}\n${p.reference}`, style: { alignItems: 'flex-start' } }, box,
        h('span', { style: { display: 'grid', gap: '1px' } },
          h('span', null, h('b', null, p.name), existing.has(p.name) ? h('span', { class: 'muted' }, ' · in dataset') : null),
          h('span', { class: 'muted', style: { fontSize: '11.5px' } }, p.description)),
        h('span', { class: 'meta mono' }, p.rsid ? xref({ label: p.rsid, url: dbsnpUrl(p.rsid), title: 'dbSNP' }, { compact: true }) : '', ` · ${p.chrom}:${fmtInt(p.position)}`));
    })]));
}

/**
 * Windows to predict: the dataset's windows plus new ones (any gene, NAME=chr:start-end regions,
 * regulatory presets) that are built from the dataset's VCF before predicting.
 */
export function windowPicker(options, datasetId, { onChange, selected = null } = {}) {
  const existing = new Set(options.genes);
  const chosen = new Set(selected && selected.length ? selected.filter((g) => existing.has(g)) : options.genes);
  const newGenes = new Map(); // name -> gene record
  const presets = new Map(); // name -> preset (new windows)
  let regionsText = '';
  let query = '';
  const changed = () => { if (onChange) onChange(); };
  const canExtend = options.extension && options.extension.available;

  const list = h('div', { class: 'checklist', style: { maxHeight: '170px' } });
  const renderList = () => {
    clear(list);
    for (const w of options.windows) {
      if (query && !w.name.toLowerCase().includes(query)) continue;
      const box = h('input', { type: 'checkbox' });
      box.checked = chosen.has(w.name);
      box.addEventListener('change', () => { if (box.checked) chosen.add(w.name); else chosen.delete(w.name); renderCount(); changed(); });
      list.appendChild(h('label', null, box, h('span', null, w.name), w.type !== 'gene' ? h('span', { class: 'muted', style: { fontSize: '11px' } }, w.type) : null,
        h('span', { class: 'meta mono' }, w.chromosome ? `${w.chromosome}:${fmtInt(w.start)}` : '')));
    }
  };
  const filter = h('input', { class: 'input', placeholder: `Filter ${options.windows.length} windows`, style: { flex: '1' } });
  filter.addEventListener('input', () => { query = filter.value.trim().toLowerCase(); renderList(); });
  const setAll = (on) => { for (const w of options.windows) if (!query || w.name.toLowerCase().includes(query)) { if (on) chosen.add(w.name); else chosen.delete(w.name); } renderList(); renderCount(); changed(); };
  const count = h('span', { class: 'muted', style: { fontSize: '12px' } });

  const added = h('div', { class: 'chips' });
  const renderAdded = () => {
    clear(added);
    for (const [name, g] of newGenes) {
      added.appendChild(h('span', { class: 'pill removable', title: g.chrom ? `${g.chrom}:${g.start}-${g.end} ${g.type || ''}` : name }, name, h('span', { class: 'muted' }, ' new'),
        h('button', { type: 'button', onclick: () => { newGenes.delete(name); renderAdded(); renderCount(); changed(); } }, '×')));
    }
  };
  const addGene = (g) => {
    if (existing.has(g.name)) { chosen.add(g.name); renderList(); toast(`${g.name} is already a window; selected it`); }
    else newGenes.set(g.name, g);
    renderAdded(); renderCount(); changed();
  };
  const regions = h('textarea', { class: 'input mono', rows: 2, style: { width: '100%' }, placeholder: 'NAME=chr:start[-end], one per line, e.g.\nmy_enhancer=chr2:1200000-1203000' });
  regions.addEventListener('input', () => { regionsText = regions.value; renderCount(); changed(); });
  const presetList = presetChecklist(options.region_presets || [], (p) => (existing.has(p.name) ? chosen.has(p.name) : presets.has(p.name)), (p, on) => {
    if (existing.has(p.name)) { if (on) chosen.add(p.name); else chosen.delete(p.name); renderList(); }
    else if (on) presets.set(p.name, p); else presets.delete(p.name);
    renderCount(); changed();
  }, { existing });

  const regionLines = () => regionsText.split('\n').map((x) => x.trim()).filter(Boolean);
  const total = () => chosen.size + newGenes.size + presets.size + regionLines().length;
  const newCount = () => newGenes.size + presets.size + regionLines().length;
  function renderCount() {
    const n = newCount();
    count.textContent = `${fmtInt(total())} selected${n ? ` · ${n} new (built first)` : ''}`;
  }

  const addHost = canExtend
    ? h('div', { style: { display: 'grid', gap: '8px' } },
      geneSearch(datasetId, addGene), added,
      h('details', { class: 'collapsible' }, h('summary', null, 'Genes of a Gene Ontology term or HGNC gene group'),
        geneSetPicker(datasetId, (genes) => { genes.forEach((g) => addGene(g)); })),
      h('details', { class: 'collapsible' }, h('summary', null, `Regulatory elements (not genes) · ${(options.region_presets || []).length}`),
        h('p', { class: 'muted', style: { fontSize: '12px', margin: '6px 0' } }, 'Windows centred on a causal or lead variant of well-studied enhancers, so distal elements sit in the middle of the window instead of near its edge.'),
        presetList),
      h('details', { class: 'collapsible' }, h('summary', null, 'Custom regions'), regions),
      h('div', { class: 'muted', style: { fontSize: '12px' } }, `New windows are built for every sample from ${options.extension.vcf} (${fmtInt(options.extension.window_size)} bp, centred on the gene or region), then predicted.`))
    : h('div', { class: 'notice' }, `New windows cannot be built for this dataset: ${(options.extension && options.extension.problems || ['unknown source']).join('; ')}. Existing windows can be predicted.`);

  renderList();
  renderAdded();
  renderCount();
  const el = h('div', { style: { display: 'grid', gap: '10px' } },
    section('Windows in the dataset', h('div', { style: { display: 'grid', gap: '6px' } },
      h('div', { style: { display: 'flex', gap: '6px', alignItems: 'center' } }, filter,
        h('button', { class: 'btn small ghost', type: 'button', onclick: () => setAll(true) }, 'All'),
        h('button', { class: 'btn small ghost', type: 'button', onclick: () => setAll(false) }, 'None')),
      list), count),
    section('Add windows', addHost));
  el.values = () => ({
    genes: [...chosen],
    new_genes: [...newGenes.keys()],
    new_regions: [...[...presets.values()].map((p) => ({ name: p.name, chrom: p.chrom, start: p.start, end: p.end })), ...regionLines()],
  });
  el.count = total;
  el.newCount = newCount;
  return el;
}
