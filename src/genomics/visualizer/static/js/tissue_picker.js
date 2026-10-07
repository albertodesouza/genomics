// Tissue / cell-type picker for the Gene products page: every ontology term with AlphaGenome RNA-seq
// (GET /api/products/tissues), searchable, each labelled with where its observed baseline comes from.
import { api } from './api.js';
import { h, clear, debounce } from './ui.js';

const BASELINE = {
  curated: { label: 'curated', title: 'Default tissue with a curated baseline', tone: 'good' },
  gtex: { label: 'GTEx', title: 'GTEx tissue: observed gene and transcript levels (TPM)', tone: 'good' },
  hpa_single_cell: { label: 'HPA cell type?', title: 'Absolute level from a Human Protein Atlas single-cell type of the same name, if there is one; transcript shares averaged over GTEx tissues', tone: '' },
  hpa_cell_line: { label: 'HPA cell line?', title: 'Absolute level from a Human Protein Atlas cell line of the same name, if there is one; transcript shares averaged over GTEx tissues', tone: '' },
  none: { label: 'relative only', title: 'No observed reference: fold changes only; transcript shares averaged over GTEx tissues', tone: 'quiet' },
};
const TYPES = [
  { value: '', label: 'All types' }, { value: 'tissue', label: 'Tissues' }, { value: 'primary_cell', label: 'Primary cells' },
  { value: 'cell_line', label: 'Cell lines' }, { value: 'in_vitro_differentiated_cells', label: 'Differentiated cells' },
];
const FALLBACK = [
  { curie: 'CL:1000458', name: 'melanocyte of skin', default: true, baseline: 'curated' },
  { curie: 'CL:2000045', name: 'foreskin melanocyte', default: true, baseline: 'curated' },
  { curie: 'EFO:0000572', name: 'lymphoblastoid cells (GTEx)', default: true, baseline: 'curated' },
  { curie: 'EFO:0002784', name: 'GM12878 (LCL)', default: true, baseline: 'curated' },
];
const MAX_ROWS = 200;
let cache = null;

function loadTissues() {
  if (!cache) cache = api('/api/products/tissues').catch((err) => { cache = null; throw err; });
  return cache;
}

function badge(kind) {
  const b = BASELINE[kind] || BASELINE.none;
  return h('span', { class: `pill cat-chip ${b.tone}`, title: b.title }, b.label);
}

/** A button showing the chosen tissue; clicking opens a searchable list. ``onChange(curie, row)``. */
export function tissuePicker(selected, onChange) {
  let rows = FALLBACK;
  let catalog = false;
  let current = selected;
  let query = '';
  let type = '';
  const label = h('span', { class: 'tissue-current' });
  const button = h('button', { class: 'btn tissue-button', type: 'button', onclick: () => toggle() }, label, h('span', { class: 'muted' }, '▾'));
  const search = h('input', { class: 'input', placeholder: 'Search 285 tissues and cell types, or a CURIE…', spellcheck: 'false' });
  const typeSel = h('select', { class: 'select', onchange: (e) => { type = e.target.value; renderList(); } }, TYPES.map((t) => h('option', { value: t.value }, t.label)));
  const list = h('div', { class: 'tissue-list' });
  const note = h('div', { class: 'muted small tissue-note' });
  const panel = h('div', { class: 'tissue-panel', hidden: true },
    h('div', { class: 'tissue-panel-head' }, search, typeSel), list, note);
  const el = h('div', { class: 'tissue-picker' }, button, panel);
  search.addEventListener('input', debounce(() => { query = search.value.trim().toLowerCase(); renderList(); }, 100));
  search.addEventListener('keydown', (e) => { if (e.key === 'Escape') close(); });

  function renderLabel() {
    const row = rows.find((r) => r.curie === current) || { curie: current, name: current };
    clear(label);
    label.append(row.name || row.curie, ' ', h('span', { class: 'muted mono small' }, row.curie));
  }

  function renderList() {
    clear(list);
    const matches = rows.filter((r) => (!type || r.type === type) && (!query || `${r.name} ${r.curie} ${r.type} ${r.gtex_tissue || ''}`.toLowerCase().includes(query)));
    const defaults = matches.filter((r) => r.default);
    const others = matches.filter((r) => !r.default);
    const section = (title, items) => {
      if (!items.length) return;
      list.appendChild(h('div', { class: 'tissue-section' }, title));
      for (const r of items.slice(0, MAX_ROWS)) {
        list.appendChild(h('button', { type: 'button', class: `tissue-row${r.curie === current ? ' on' : ''}`, onclick: () => pick(r) },
          h('span', { class: 'tissue-name' }, r.name || r.curie),
          h('span', { class: 'muted small' }, (r.type || '').replace(/_/g, ' ')),
          badge(r.baseline),
          r.junctions === 0 ? h('span', { class: 'pill cat-chip quiet', title: 'AlphaGenome has no splice-junction track for this term: no splicing shift' }, 'no junctions') : null,
          h('span', { class: 'mono muted small tissue-curie' }, r.curie)));
      }
    };
    section('Defaults (stored predictions for this dataset where available)', defaults);
    section(`All AlphaGenome tissues and cell types with RNA-seq (${others.length})`, others);
    if (!matches.length) {
      const curie = search.value.trim();
      list.appendChild(/^[A-Za-z]+:\d+$/.test(curie)
        ? h('button', { type: 'button', class: 'tissue-row', onclick: () => pick({ curie, name: curie }) }, `Use ${curie}`)
        : h('div', { class: 'empty' }, 'No match.'));
    }
    note.textContent = catalog
      ? 'A tissue not stored in the dataset is predicted on demand (three AlphaGenome calls, cached). GTEx = observed gene and transcript levels; HPA = an absolute gene level when the Human Protein Atlas has a cell type or line of that name (checked when the report runs).'
      : 'The AlphaGenome track catalog is unavailable (no backend?): only the default tissues are listed; any CURIE can be typed.';
  }

  function pick(r) {
    current = r.curie;
    renderLabel();
    close();
    onChange(r.curie, r);
  }

  function onDoc(e) { if (!el.contains(e.target)) close(); }
  function close() { panel.hidden = true; document.removeEventListener('mousedown', onDoc); }
  function toggle() {
    if (!panel.hidden) { close(); return; }
    panel.hidden = false;
    renderList();
    search.focus();
    document.addEventListener('mousedown', onDoc);
  }

  renderLabel();
  loadTissues().then((data) => {
    if (data.tissues && data.tissues.length) { rows = data.tissues; catalog = !!data.catalog; }
    renderLabel();
    if (!panel.hidden) renderList();
  }).catch(() => {});
  el.value = () => current;
  return el;
}
