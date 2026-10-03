// Samples: faceted cohort builder, virtualized table, sample details, pinning and view export.
import { api } from '../api.js';
import { navigate } from '../app.js';
import { state, ds, filterRows, setFilters, setPinned, togglePinned, categoricalFields } from '../state.js';
import { h, clear, icon, fmtInt, drawer, modal, toast, downloadText, debounce, field, select, errorBox, escapeHtml } from '../ui.js';

const ROW_H = 30;

export async function mount(root) {
  const layout = h('div', { class: 'samples-layout' });
  root.appendChild(layout);
  const facetsEl = h('aside', { class: 'facets', 'aria-label': 'Filters' });
  const main = h('section', { class: 'samples-main' });
  layout.append(facetsEl, main);

  const fields = state.samples.fields;
  const columns = fields.map((f) => f.name);
  const col = state.samples.col;
  let sort = { field: null, dir: 1 };
  let rows = [];
  let active = null;

  // ----- header & actions -----
  const countEl = h('div', null);
  const search = h('input', { class: 'input', type: 'search', placeholder: 'Search samples or any field…', value: state.search, style: { width: '280px' } });
  search.addEventListener('input', debounce(() => { state.search = search.value; refresh(); }, 120));
  const head = h('div', { style: { padding: '14px 16px', display: 'grid', gap: '10px', borderBottom: '1px solid var(--line)', background: 'var(--surface)' } },
    h('div', { style: { display: 'flex', alignItems: 'center', gap: '12px', flexWrap: 'wrap' } },
      h('h1', null, 'Samples'), countEl, h('div', { style: { flex: '1' } }), search),
    h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap', alignItems: 'center' } },
      h('button', { class: 'btn', title: 'Pin the first 5 matching samples', onclick: () => pinFirst(5) }, icon('pin', 14), 'Pin first 5'),
      h('button', { class: 'btn', title: 'Pin one sample from each value of a field', onclick: pinPerGroup }, icon('pin', 14), 'Pin one per group…'),
      h('button', { class: 'btn ghost', onclick: () => setPinned([]) }, 'Clear pins'),
      h('div', { style: { flex: '1' } }),
      h('button', { class: 'btn', onclick: exportCsv }, icon('download', 14), 'Export CSV'),
      h('button', { class: 'btn', onclick: () => viewBuilder(rows) }, 'Save as training view…'),
      h('button', { class: 'btn primary', onclick: () => navigate('tracks') }, 'Compare pinned in Tracks')),
    h('div', { class: 'chips', id: 'pinnedChips' }));
  const pinnedChips = head.querySelector('#pinnedChips');

  // ----- virtual table -----
  const widths = columns.map((c) => (c === 'sample_id' ? 130 : 150));
  const template = `40px ${widths.map((w) => `minmax(${w}px, 1fr)`).join(' ')}`;
  const headRow = h('div', { class: 'vtable-head', style: { gridTemplateColumns: template } },
    h('div', { title: 'Pinned' }, icon('pin', 13)),
    fields.map((f) => h('div', { onclick: () => { sort = { field: f.name, dir: sort.field === f.name ? -sort.dir : 1 }; refresh(); } }, f.label, h('span', { class: 'sort-mark' }))));
  const body = h('div', { class: 'vtable-body' });
  const vtable = h('div', { class: 'vtable', tabindex: '0' }, headRow, body);
  vtable.addEventListener('scroll', () => renderRows());
  main.append(head, vtable);

  function renderRows() {
    const top = vtable.scrollTop;
    const height = vtable.clientHeight || 600;
    const first = Math.max(0, Math.floor(top / ROW_H) - 5);
    const last = Math.min(rows.length, Math.ceil((top + height) / ROW_H) + 5);
    body.style.height = `${rows.length * ROW_H}px`;
    clear(body);
    const pinned = new Set(state.pinned);
    const frag = document.createDocumentFragment();
    for (let i = first; i < last; i++) {
      const row = rows[i];
      const id = row[col.sample_id];
      const box = h('input', { type: 'checkbox', 'aria-label': `Pin ${id}`, onclick: (e) => { e.stopPropagation(); togglePinned(id); } });
      box.checked = pinned.has(id);
      const el = h('div', { class: `vtable-row${active === id ? ' active' : ''}`, style: { top: `${i * ROW_H}px`, gridTemplateColumns: template }, onclick: () => openSample(id) },
        h('div', null, box),
        columns.map((c) => h('div', { title: String(row[col[c]] ?? '') }, c === 'sample_id' ? h('b', null, row[col[c]]) : (row[col[c]] ?? '–'))));
      frag.appendChild(el);
    }
    body.appendChild(frag);
  }

  function renderPinned() {
    clear(pinnedChips);
    if (!state.pinned.length) { pinnedChips.appendChild(h('span', { class: 'muted', style: { fontSize: '12px' } }, 'No pinned individuals. Tick rows to pin them for Tracks and Sequence.')); return; }
    pinnedChips.appendChild(h('span', { class: 'label', style: { alignSelf: 'center' } }, 'Pinned'));
    for (const id of state.pinned) pinnedChips.appendChild(h('span', { class: 'pill removable' }, id, h('button', { title: `Unpin ${id}`, onclick: () => togglePinned(id) }, '×')));
  }

  // ----- facets -----
  const facetSearch = {};
  function renderFacets() {
    clear(facetsEl);
    const cats = categoricalFields();
    facetsEl.appendChild(h('div', { class: 'facet' },
      h('div', { class: 'facet-head' }, h('h3', null, 'Cohort filters'), Object.keys(state.filters).length ? h('button', { class: 'btn small ghost', onclick: () => setFilters({}) }, 'Reset') : null),
      h('p', { class: 'muted', style: { margin: 0, fontSize: '12px' } }, 'The cohort drives group means and population heatmaps on the Tracks page.')));
    for (const f of cats) {
      const selected = new Set(state.filters[f.name] || []);
      const base = filterRows(state.filters, state.search, f.name);
      const counts = new Map();
      for (const row of base) { const v = String(row[col[f.name]] ?? ''); if (v) counts.set(v, (counts.get(v) || 0) + 1); }
      const values = f.counts.map((c) => c.value);
      const max = Math.max(1, ...counts.values());
      const needle = (facetSearch[f.name] || '').toLowerCase();
      const list = h('div', { class: 'facet-values' });
      for (const value of values) {
        if (needle && !value.toLowerCase().includes(needle)) continue;
        const n = counts.get(value) || 0;
        const box = h('input', { type: 'checkbox' });
        box.checked = selected.has(value);
        box.addEventListener('change', () => {
          const next = new Set(state.filters[f.name] || []);
          if (box.checked) next.add(value); else next.delete(value);
          setFilters({ ...state.filters, [f.name]: [...next] });
        });
        list.appendChild(h('label', { class: `facet-value${n ? '' : ' zero'}` }, box, h('span', { class: 'fv-name', title: value }, value), h('span', { class: 'fv-count' }, fmtInt(n)), h('span', { class: 'fv-bar', style: { width: `${(n / max) * 100}%` } })));
      }
      const searchBox = values.length > 10 ? h('input', { class: 'input', type: 'search', placeholder: `Filter ${f.label.toLowerCase()}…`, value: facetSearch[f.name] || '', style: { height: '26px' } }) : null;
      if (searchBox) searchBox.addEventListener('input', debounce(() => { facetSearch[f.name] = searchBox.value; renderFacets(); const el = facetsEl.querySelector(`[data-facet="${f.name}"] input[type=search]`); if (el) { el.focus(); el.setSelectionRange(el.value.length, el.value.length); } }, 150));
      facetsEl.appendChild(h('div', { class: 'facet', dataset: { facet: f.name } },
        h('div', { class: 'facet-head' }, h('h3', null, f.label), selected.size ? h('button', { class: 'btn small ghost', onclick: () => { const next = { ...state.filters }; delete next[f.name]; setFilters(next); } }, `Clear (${selected.size})`) : h('span', { class: 'muted', style: { fontSize: '11.5px' } }, `${values.length} values`)),
        searchBox, list));
    }
  }

  function refresh() {
    rows = filterRows(state.filters, state.search);
    if (sort.field) {
      const i = col[sort.field];
      const numeric = fields.find((f) => f.name === sort.field)?.kind === 'numeric';
      rows = rows.slice().sort((a, b) => {
        const x = a[i]; const y = b[i];
        if (x === y) return 0;
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        return (numeric ? x - y : String(x).localeCompare(String(y), undefined, { numeric: true })) * sort.dir;
      });
    }
    headRow.querySelectorAll('.sort-mark').forEach((el, idx) => { el.textContent = fields[idx].name === sort.field ? (sort.dir > 0 ? ' ↑' : ' ↓') : ''; });
    clear(countEl).append(h('span', { class: 'pill' }, `${fmtInt(rows.length)} of ${fmtInt(state.samples.rows.length)}`));
    renderFacets();
    renderRows();
  }

  function pinFirst(n) { setPinned([...state.pinned, ...rows.slice(0, n).map((r) => r[col.sample_id])]); }

  function pinPerGroup() {
    const cats = categoricalFields();
    if (!cats.length) return;
    let fieldName = cats[0].name;
    const m = modal('Pin one sample per group', h('div', { style: { display: 'grid', gap: '10px' } },
      h('p', { class: 'secondary', style: { margin: 0 } }, 'Pins the first matching sample of every value (up to 8) among the currently filtered rows — handy for comparing populations side by side.'),
      field('Group by', select(cats.map((f) => ({ value: f.name, label: f.label })), fieldName, (v) => { fieldName = v; }))),
      [h('button', { class: 'btn', onclick: () => m.close() }, 'Cancel'), h('button', { class: 'btn primary', onclick: () => {
        const seen = new Map();
        for (const row of rows) { const v = row[col[fieldName]]; if (v && !seen.has(v)) seen.set(v, row[col.sample_id]); if (seen.size >= 8) break; }
        setPinned([...seen.values()]);
        m.close();
        toast(`Pinned ${seen.size} samples`);
      } }, 'Pin')]);
  }

  function exportCsv() {
    const esc = (v) => { const s = v === null || v === undefined ? '' : String(v); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
    const text = [columns.join(','), ...rows.map((r) => columns.map((c) => esc(r[col[c]])).join(','))].join('\n');
    downloadText(`${state.datasetId}_samples_${rows.length}.csv`, text, 'text/csv');
  }

  async function openSample(id) {
    active = id;
    renderRows();
    const body = h('div', { style: { display: 'grid', gap: '16px' } }, h('div', { class: 'muted' }, 'Loading…'));
    const isPinned = state.pinned.includes(id);
    drawer(id, body, {
      onClose: () => { active = null; renderRows(); },
      actions: h('button', { class: 'btn small', onclick: () => { togglePinned(id); } }, icon('pin', 13), isPinned ? 'Unpin' : 'Pin'),
    });
    try {
      const detail = await api(`${ds()}/samples/${encodeURIComponent(id)}`);
      clear(body);
      const sample = detail.sample;
      body.appendChild(h('dl', { class: 'kv' }, Object.entries(sample).flatMap(([k, v]) => [h('dt', null, k), h('dd', null, v === null || v === undefined || v === '' ? '–' : String(v))])));
      const meta = detail.individual_metadata || {};
      const extra = Object.entries(meta).filter(([k, v]) => !(k in sample) && k !== 'windows' && (typeof v !== 'object' || v === null));
      if (extra.length) body.appendChild(h('div', null, h('div', { class: 'label', style: { marginBottom: '6px' } }, 'Individual metadata'), h('dl', { class: 'kv' }, extra.flatMap(([k, v]) => [h('dt', null, k), h('dd', null, String(v))]))));
      body.appendChild(h('div', null,
        h('div', { class: 'label', style: { marginBottom: '6px' } }, `Gene windows (${detail.windows.length})`),
        h('div', { class: 'table-wrap', style: { maxHeight: '360px' } }, h('table', { class: 'table' },
          h('thead', null, h('tr', null, h('th', null, 'Gene'), h('th', null, 'Predictions'), h('th', null, 'VCF'), h('th', null, ''))),
          h('tbody', null, detail.windows.map((w) => h('tr', null,
            h('td', null, h('b', null, w.gene)),
            h('td', null, w.outputs.length ? `${w.outputs.join(', ')} · ${w.haplotypes.join('/')}` : h('span', { class: 'muted' }, 'none')),
            h('td', null, w.has_vcf ? 'yes' : h('span', { class: 'muted' }, 'no')),
            h('td', null, h('div', { style: { display: 'flex', gap: '4px', justifyContent: 'flex-end' } },
              h('button', { class: 'btn small', disabled: !w.outputs.length, onclick: () => { setPinned([id, ...state.pinned.filter((x) => x !== id)]); navigate('tracks', { gene: w.gene }); } }, 'Tracks'),
              h('button', { class: 'btn small', onclick: () => { setPinned([id, ...state.pinned.filter((x) => x !== id)]); navigate('sequence', { gene: w.gene }); } }, 'Sequence'))))))))));
    } catch (err) { clear(body).appendChild(errorBox(err)); }
  }

  renderPinned();
  refresh();
  const onResize = () => renderRows();
  window.addEventListener('resize', onResize);
  return {
    onEvent(topic) {
      if (topic === 'cohort') refresh();
      if (topic === 'pinned') { renderPinned(); renderRows(); }
    },
    unmount() { window.removeEventListener('resize', onResize); },
  };
}

function viewBuilder(rows) {
  const col = state.samples.col;
  const sampleIds = rows.map((r) => r[col.sample_id]);
  const genes = state.summary.genes.map((g) => g.gene);
  const defaults = state.summary.genes.filter((g) => g.in_metadata).map((g) => g.gene);
  const name = h('input', { class: 'input', value: 'new_view' });
  const description = h('input', { class: 'input', placeholder: 'What is this cohort for?' });
  const outputs = h('input', { class: 'input mono', value: (state.summary.outputs[0] || 'rna_seq') });
  const window = h('input', { class: 'input', type: 'number', value: '32768', min: '1' });
  const downsample = h('input', { class: 'input', type: 'number', value: '1', min: '1' });
  const normalization = select(['log', 'zscore', 'minmax_keep_zero'], 'log', () => {});
  const outputPath = h('input', { class: 'input mono', placeholder: 'default: apps/views/<name>.view.json' });
  const geneBoxes = genes.map((g) => { const b = h('input', { type: 'checkbox', value: g }); b.checked = defaults.length ? defaults.includes(g) : true; return h('label', { class: 'check' }, b, g); });
  const preview = h('pre', { class: 'code-block' }, '…');
  const payload = (overwrite = false) => ({
    name: name.value, description: description.value, alphagenome_outputs: outputs.value,
    window_center_size: Number(window.value), downsample_factor: Number(downsample.value), normalization_method: normalization.value,
    genes: geneBoxes.map((l) => l.querySelector('input')).filter((b) => b.checked).map((b) => b.value),
    sample_ids: sampleIds, output_path: outputPath.value || null, overwrite,
  });
  const update = debounce(async () => {
    try {
      const res = await api(`${ds()}/views/preview`, { method: 'POST', body: payload() });
      preview.innerHTML = `<span class="muted"># ${escapeHtml(res.output_path)}\n# ${res.counts.samples} samples · ${res.counts.genes} genes</span>\n${escapeHtml(res.yaml_snippet)}`;
    } catch (err) { preview.textContent = err.message; }
  }, 200);
  const form = h('div', { style: { display: 'grid', gap: '12px' } },
    h('p', { class: 'secondary', style: { margin: 0 } }, `Exports the ${fmtInt(sampleIds.length)} samples currently matching your filters as a .view.json usable as dataset_input in genotype configs.`),
    h('div', { class: 'grid cols-2', style: { gap: '10px' } }, field('Name', name), field('Description', description), field('AlphaGenome outputs', outputs), field('Normalization', normalization), field('Window center size', window), field('Downsample factor', downsample)),
    field('Output path', outputPath),
    h('div', null, h('div', { class: 'label', style: { marginBottom: '6px' } }, 'Genes'), h('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(110px, 1fr))', gap: '4px' } }, geneBoxes)),
    h('div', null, h('div', { class: 'label', style: { marginBottom: '6px' } }, 'Preview'), preview));
  form.addEventListener('input', update);
  form.addEventListener('change', update);
  const save = async (overwrite) => {
    try {
      const res = await api(`${ds()}/views/save`, { method: 'POST', body: payload(overwrite) });
      toast(`Saved ${res.saved_path}`);
      m.close();
    } catch (err) {
      if (err.status === 409 && confirm(`${err.message}\n\nOverwrite?`)) save(true);
      else if (err.status !== 409) toast(err.message, 'error');
    }
  };
  const m = modal('Save cohort as training view', form, [h('button', { class: 'btn', onclick: () => m.close() }, 'Cancel'), h('button', { class: 'btn primary', onclick: () => save(false) }, 'Save view')]);
  update();
}
