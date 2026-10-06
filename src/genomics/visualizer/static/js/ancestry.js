// Ancestry panel (Samples page): genotype PCA scatter coloured by a sample field, PCs as sample fields,
// and matching two groups on the PCs (stored as a categorical sample field).
import { api, apiJob, isAbort, Latest } from './api.js';
import { navigate } from './app.js';
import { state, ds, filterRows, setFilters, categoricalFields, cohortDescription, reloadSamples, seriesColor, togglePinned } from './state.js';
import { h, clear, select, field, fmtInt, fmtNum, fmtPct, toast, errorBox, jobOverlay, modal, showTooltip, hideTooltip, escapeHtml, checkbox } from './ui.js';
import { scatterPlot, onResize } from './plot.js';
import { exportMenu, slug } from './figure.js';

const STORE_KEY = 'gv.ancestry';

function loadPrefs() {
  try { return JSON.parse(localStorage.getItem(`${STORE_KEY}.${state.datasetId}`) || '{}') || {}; } catch (e) { return {}; }
}
function savePrefs(prefs) {
  try { localStorage.setItem(`${STORE_KEY}.${state.datasetId}`, JSON.stringify(prefs)); } catch (e) { /* private mode */ }
}

export function mountAncestry(host) {
  const prefs = loadPrefs();
  const cats = categoricalFields();
  const cfg = {
    genes: prefs.genes || null, min_maf: prefs.min_maf || 0.05, spacing: prefs.spacing ?? 2000, components: prefs.components || 10,
    x: prefs.x || 0, y: prefs.y || 1,
    color: cats.some((f) => f.name === prefs.color) ? prefs.color : ((cats.find((f) => f.name === 'superpopulation') || cats[0] || {}).name || ''),
  };
  let windows = null; let pca = null; let plot = null; let derived = [];
  const latest = new Latest();
  const disposers = [];

  const scroll = h('div', { class: 'ancestry', style: { flex: '1', minHeight: '0', overflow: 'auto', padding: '16px', display: 'grid', gap: '16px', alignContent: 'start' } });
  host.appendChild(scroll);
  const pcaBody = h('div', { class: 'card-body', style: { position: 'relative', minHeight: '120px' } });
  const pcaHead = h('div', { class: 'card-head' });
  const pcaCard = h('section', { class: 'card' }, pcaHead, pcaBody);
  const matchBody = h('div', { class: 'card-body' });
  const matchCard = h('section', { class: 'card' }, h('div', { class: 'card-head' }, h('h2', null, 'Match two groups on ancestry')), matchBody);
  scroll.append(pcaCard, matchCard);

  const persist = () => savePrefs({ genes: cfg.genes, min_maf: cfg.min_maf, spacing: cfg.spacing, components: cfg.components, x: cfg.x, y: cfg.y, color: cfg.color });
  const query = () => ({ genes: (cfg.genes || windows.default).join(','), min_maf: cfg.min_maf, spacing: cfg.spacing, components: cfg.components });

  // ------------------------------------------------------------------ PCA
  async function init() {
    try {
      windows = await api(`${ds()}/ancestry/windows`);
    } catch (err) { pcaBody.appendChild(errorBox(err)); return; }
    derived = windows.derived || [];
    if (cfg.genes) cfg.genes = cfg.genes.filter((g) => windows.windows.some((w) => w.gene === g));
    if (!cfg.genes || !cfg.genes.length) cfg.genes = null;
    renderHead();
    const chosen = cfg.genes || windows.default;
    const allCached = chosen.every((g) => (windows.windows.find((w) => w.gene === g) || {}).cached);
    if (allCached) { compute(); return; }
    try {
      const hit = await api(`${ds()}/ancestry/pca`, { params: { ...query(), cached_only: 1 }, signal: latest.next() });
      if (hit && hit.cached === false) renderIdle();
      else { pca = hit; renderPca(); renderMatch(); }
    } catch (err) { if (!isAbort(err)) renderIdle(); }
  }

  function windowSummary() {
    const chosen = cfg.genes || windows.default;
    const byGene = new Map(windows.windows.map((w) => [w.gene, w]));
    const coverage = Math.min(...chosen.map((g) => (byGene.get(g) || { coverage: 0 }).coverage));
    const listed = chosen.filter((g) => (byGene.get(g) || {}).listed).length;
    const kind = listed === chosen.length ? 'all listed in the dataset metadata (pigmentation loci)' : listed ? `${listed} listed in the metadata` : 'control windows';
    return `${chosen.length} window${chosen.length === 1 ? '' : 's'} (${kind}) · about ${fmtInt(coverage * windows.samples)} samples genotyped in all`;
  }

  function renderHead() {
    clear(pcaHead).append(
      h('h2', null, 'Genotype PCA'),
      h('span', { class: 'muted', style: { fontSize: '12px', flex: '1' } }, windowSummary()),
      h('button', { class: 'btn small', onclick: chooseWindows }, 'Windows & filters…'),
      exportMenu(() => figureOptions()));
  }

  function renderIdle() {
    clear(pcaBody).append(h('div', { class: 'empty', style: { display: 'grid', gap: '10px', justifyItems: 'center' } },
      h('p', { style: { margin: 0 } }, 'PCA of common SNVs (MAF ≥ ', fmtPct(cfg.min_maf, 0), `, one per ${fmtInt(cfg.spacing)} bp) from the window VCFs. The first run reads every sample's VCFs (seconds per window) and is cached.`),
      h('button', { class: 'btn primary', onclick: () => compute() }, 'Compute PCA')));
  }

  async function compute() {
    const signal = latest.next();
    clear(pcaBody);
    const overlay = jobOverlay(pcaBody, 'Genotype PCA');
    try {
      pca = await apiJob(`${ds()}/ancestry/pca`, { params: query(), signal, onProgress: (job) => overlay.update(job) });
      renderPca();
      renderMatch();
    } catch (err) {
      if (!isAbort(err)) { pca = null; clear(pcaBody).appendChild(errorBox(err)); }
    } finally { overlay.remove(); }
  }

  function chooseWindows() {
    const chosen = new Set(cfg.genes || windows.default);
    const boxes = windows.windows.map((w) => {
      const box = h('input', { type: 'checkbox', value: w.gene });
      box.checked = chosen.has(w.gene);
      return { w, box, el: h('label', { class: 'check', title: `${fmtPct(w.coverage, 0)} of samples have a VCF${w.listed ? ' · listed in the dataset metadata' : ''}` }, box, w.gene, h('span', { class: 'muted', style: { fontSize: '11px' } }, ` ${fmtPct(w.coverage, 0)}${w.listed ? ' ●' : ''}`)) };
    });
    const preset = (pick) => boxes.forEach((b) => { b.box.checked = pick(b.w); });
    const minMaf = h('input', { class: 'input', type: 'number', min: '0.01', max: '0.49', step: '0.01', value: String(cfg.min_maf) });
    const spacing = h('input', { class: 'input', type: 'number', min: '0', step: '500', value: String(cfg.spacing) });
    const components = h('input', { class: 'input', type: 'number', min: '2', max: '20', value: String(cfg.components) });
    const body = h('div', { style: { display: 'grid', gap: '12px' } },
      h('p', { class: 'secondary', style: { margin: 0 } }, 'Only samples genotyped in every chosen window get PCs. ● marks windows listed in the dataset metadata (the pigmentation loci, under selection); the others are control windows. Percentages are the share of samples with a VCF.'),
      h('div', { style: { display: 'flex', gap: '6px', flexWrap: 'wrap' } },
        h('button', { class: 'btn small', onclick: () => preset((w) => windows.default.includes(w.gene)) }, 'Default'),
        h('button', { class: 'btn small', onclick: () => preset((w) => !w.listed) }, 'Control windows'),
        h('button', { class: 'btn small', onclick: () => preset((w) => w.listed) }, 'Listed windows'),
        h('button', { class: 'btn small', onclick: () => preset(() => true) }, 'All')),
      h('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(150px, 1fr))', gap: '4px', maxHeight: '40vh', overflow: 'auto' } }, boxes.map((b) => b.el)),
      h('div', { class: 'grid cols-3', style: { gap: '10px' } }, field('Minimum MAF', minMaf), field('Site spacing (bp)', spacing), field('Components', components)));
    const m = modal('PCA windows and filters', body, [
      h('button', { class: 'btn', onclick: () => m.close() }, 'Cancel'),
      h('button', { class: 'btn primary', onclick: () => {
        const genes = boxes.filter((b) => b.box.checked).map((b) => b.w.gene);
        if (!genes.length) { toast('Choose at least one window', 'error'); return; }
        const same = genes.length === windows.default.length && genes.every((g) => windows.default.includes(g));
        cfg.genes = same ? null : genes;
        cfg.min_maf = Math.min(0.49, Math.max(0.01, Number(minMaf.value) || 0.05));
        cfg.spacing = Math.max(0, Math.round(Number(spacing.value) || 0));
        cfg.components = Math.min(20, Math.max(2, Math.round(Number(components.value) || 10)));
        cfg.x = Math.min(cfg.x, cfg.components - 1); cfg.y = Math.min(cfg.y, cfg.components - 1);
        persist(); m.close(); renderHead(); compute();
      } }, 'Compute')]);
  }

  function pcLabel(i) { return `PC${i + 1} (${fmtPct(pca.explained[i], 1)})`; }

  function renderPca() {
    const k = pca.scores[0] ? pca.scores[0].length : 0;
    const pcOptions = Array.from({ length: k }, (_, i) => ({ value: String(i), label: pcLabel(i) }));
    const xSel = select(pcOptions, String(cfg.x), (v) => { cfg.x = Number(v); persist(); draw(); }, { 'aria-label': 'X axis' });
    const ySel = select(pcOptions, String(cfg.y), (v) => { cfg.y = Number(v); persist(); draw(); }, { 'aria-label': 'Y axis' });
    const colorSel = select([{ value: '', label: 'None' }, ...categoricalFields().map((f) => ({ value: f.name, label: f.label }))], cfg.color, (v) => { cfg.color = v; persist(); draw(); }, { 'aria-label': 'Colour by' });
    const pcFields = derived.find((d) => d.kind === 'pca');
    const countSel = select([2, 4, 6, 10].filter((n) => n <= k).map((n) => ({ value: String(n), label: `PC1–PC${n}` })), String(Math.min(4, k)), () => {});
    const canvasHost = h('div', { class: 'canvas-host', style: { position: 'relative' } });
    const canvas = h('canvas');
    canvasHost.appendChild(canvas);
    const legend = h('div', { class: 'legend', style: { padding: '6px 0 0', borderBottom: '0', background: 'transparent' } });
    const note = h('p', { class: 'muted', style: { margin: '4px 0 0', fontSize: '12px' } });
    clear(pcaBody).append(
      h('div', { style: { display: 'flex', gap: '10px', flexWrap: 'wrap', alignItems: 'end', marginBottom: '8px' }, 'data-export': 'skip' },
        field('X', xSel), field('Y', ySel), field('Colour by', colorSel),
        h('div', { style: { flex: '1' } }),
        h('div', { style: { display: 'flex', gap: '6px', alignItems: 'end' } },
          field('Sample fields', countSel),
          h('button', { class: 'btn small', title: 'Add the PCs as numeric sample fields pc1, pc2, … (Samples columns, CSV export)', onclick: () => setPcFields(Number(countSel.value)) }, pcFields ? 'Update PC fields' : 'Add PC fields'),
          pcFields ? h('button', { class: 'btn small ghost', onclick: () => setPcFields(0) }, 'Remove') : null)),
      canvasHost, legend, note);

    function draw() {
      const width = Math.max(320, canvasHost.clientWidth || pcaBody.clientWidth || 600);
      const cohort = new Set(filterRows(state.filters, state.search).map((r) => r[state.samples.col.sample_id]));
      const colorCol = cfg.color ? state.samples.col[cfg.color] : undefined;
      const values = cfg.color ? (categoricalFields().find((f) => f.name === cfg.color) || { counts: [] }).counts.map((c) => c.value) : [];
      const colorOf = new Map(values.map((v, i) => [v, seriesColor(i)]));
      const grey = getComputedStyle(document.documentElement).getPropertyValue('--ink-3').trim();
      const counts = new Map();
      const points = pca.samples.map((s, i) => {
        const row = state.samples.index.get(s);
        const inCohort = cohort.has(s);
        const value = row && colorCol !== undefined ? String(row[colorCol] ?? '') : '';
        if (inCohort) counts.set(value, (counts.get(value) || 0) + 1);
        const color = !inCohort ? grey : cfg.color ? colorOf.get(value) || grey : null;
        return { x: pca.scores[i][cfg.x], y: pca.scores[i][cfg.y], back: !inCohort, r: inCohort ? 2.6 : 2, color };
      });
      plot = scatterPlot(canvas, { points, width, height: Math.min(560, Math.max(320, width * 0.55)), xLabel: pcLabel(cfg.x), yLabel: pcLabel(cfg.y) });
      clear(legend);
      if (cfg.color) for (const v of values) if (counts.get(v)) legend.append(h('span', { class: 'legend-item' }, h('span', { class: 'swatch', style: { background: colorOf.get(v) } }), `${v} (${fmtInt(counts.get(v))})`));
      const inCohort = pca.samples.filter((s) => cohort.has(s)).length;
      note.textContent = `${fmtInt(pca.samples.length)} of ${fmtInt(pca.sample_count)} samples have PCs (genotyped in all ${pca.genes.length} windows) · ${fmtInt(pca.n_sites)} SNVs · ${fmtInt(inCohort)} in the cohort (coloured; the rest grey) · click a point to pin the sample.`;
    }
    canvas.addEventListener('mousemove', (e) => {
      if (!plot) return;
      const rect = canvas.getBoundingClientRect();
      const i = plot.nearest(e.clientX - rect.left, e.clientY - rect.top);
      if (i < 0) { hideTooltip(); return; }
      const s = pca.samples[i];
      const row = state.samples.index.get(s);
      const field = (name) => (row && state.samples.col[name] !== undefined ? row[state.samples.col[name]] : null);
      const extra = ['population', 'superpopulation', cfg.color].filter((f, j, a) => f && a.indexOf(f) === j && field(f) !== null && field(f) !== undefined);
      showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(s)}</div>${extra.map((f) => `<div class="tt-row"><span class="tt-key">${escapeHtml(f)}</span><span class="tt-val">${escapeHtml(String(field(f)))}</span></div>`).join('')}<div class="tt-sub">PC${cfg.x + 1} ${fmtNum(pca.scores[i][cfg.x])} · PC${cfg.y + 1} ${fmtNum(pca.scores[i][cfg.y])}</div>`);
    });
    canvas.addEventListener('mouseleave', hideTooltip);
    canvas.addEventListener('click', (e) => {
      if (!plot) return;
      const rect = canvas.getBoundingClientRect();
      const i = plot.nearest(e.clientX - rect.left, e.clientY - rect.top);
      if (i >= 0) { togglePinned(pca.samples[i]); toast(`${state.pinned.includes(pca.samples[i]) ? 'Pinned' : 'Unpinned'} ${pca.samples[i]}`); }
    });
    redraw = draw;
    draw();
  }

  let redraw = () => {};

  async function setPcFields(count) {
    try {
      const res = await api(`${ds()}/ancestry/fields`, { method: 'POST', body: { params: query(), count } });
      toast(count ? `Added ${res.fields.join(', ')} to the samples` : 'Removed the PC fields');
      await reloadSamples();
    } catch (err) { toast(err.message, 'error'); }
  }

  // ------------------------------------------------------------------ matching
  function renderMatch() {
    clear(matchBody);
    const cats = categoricalFields();
    if (!pca) { matchBody.appendChild(h('p', { class: 'muted' }, 'Compute the PCA first.')); return; }
    const pick = { field: cats.some((f) => f.name === 'pigmentation') ? 'pigmentation' : (cats.find((f) => f.name === 'superpopulation') || cats[0] || {}).name, a: new Set(), b: new Set() };
    const k = pca.scores[0].length;
    const kSel = select(Array.from({ length: k }, (_, i) => ({ value: String(i + 1), label: `PC1–PC${i + 1}` })).slice(1), String(Math.min(4, k)), () => {});
    const caliper = h('input', { class: 'input', type: 'number', min: '0', step: '0.05', value: '0.2', title: 'Maximum pair distance, in pooled SDs of the PC space; 0 = no limit' });
    const name = h('input', { class: 'input mono', value: 'matched' });
    const valuesHost = h('div', { class: 'grid cols-2', style: { gap: '10px' } });
    const result = h('div');
    const fieldSel = select(cats.map((f) => ({ value: f.name, label: f.label })), pick.field, (v) => { pick.field = v; pick.a.clear(); pick.b.clear(); renderValues(); });
    function renderValues() {
      const f = cats.find((x) => x.name === pick.field) || { counts: [] };
      const list = (set, other) => h('div', { class: 'facet-values', style: { maxHeight: '180px' } }, f.counts.map((c) => checkbox(`${c.value} (${fmtInt(c.count)})`, set.has(c.value), (on) => { if (on) { set.add(c.value); other.delete(c.value); } else set.delete(c.value); renderValues(); }, { disabled: other.has(c.value) })));
      clear(valuesHost).append(field('Group A (each sample gets a partner)', list(pick.a, pick.b)), field('Group B (partners drawn from here)', list(pick.b, pick.a)));
    }
    renderValues();
    const go = h('button', { class: 'btn primary', onclick: async () => {
      if (!pick.a.size || !pick.b.size) { toast('Pick values for both groups', 'error'); return; }
      go.disabled = true;
      try {
        const res = await api(`${ds()}/ancestry/match`, { method: 'POST', body: { params: query(), field: pick.field, group_a: [...pick.a], group_b: [...pick.b], k: Number(kSel.value), caliper: Number(caliper.value), name: name.value.trim(), filters: state.filters, replace: derived.some((d) => d.kind === 'match' && d.name === name.value.trim().toLowerCase()) } });
        showResult(res);
        await reloadSamples();
      } catch (err) { toast(err.message, 'error', 8000); } finally { go.disabled = false; }
    } }, 'Match');
    matchBody.append(
      h('p', { class: 'secondary', style: { margin: '0 0 10px' } }, 'Pairs every sample of group A with its nearest unused sample of group B on the first PCs (greedy 1:1, within the caliper), among samples of the current cohort. The matched samples become a sample field, usable as a cohort filter and for group means on Tracks, so a comparison is not driven by ancestry differences between the groups.'),
      h('div', { class: 'grid cols-3', style: { gap: '10px', marginBottom: '10px' } }, field('Field', fieldSel), field('Match on', kSel), field('Caliper (SD)', caliper)),
      valuesHost,
      h('div', { style: { display: 'flex', gap: '10px', alignItems: 'end', marginTop: '10px' } }, field('Save as field', name), go),
      result);
    renderDerived();

    function showResult(res) {
      const rows = res.balance_before.map((b, i) => h('tr', null, h('td', null, `PC${i + 1}`), h('td', null, fmtNum(b, 2)), h('td', null, fmtNum(res.balance_after[i], 2))));
      clear(result).append(h('div', { class: 'card', style: { marginTop: '12px', padding: '12px', display: 'grid', gap: '8px' } },
        h('div', null, h('b', null, `${fmtInt(res.pairs)} pairs`), ` from ${fmtInt(res.group_a)} ${res.label_a} and ${fmtInt(res.group_b)} ${res.label_b} samples with PCs`, res.caliper_distance ? ` · caliper ${fmtNum(res.caliper_distance, 3)} (median pair distance ${fmtNum(res.median_distance, 3)})` : ''),
        h('table', { class: 'table compact' }, h('thead', null, h('tr', null, h('th', null, ''), h('th', null, 'SMD before'), h('th', null, 'SMD after'))), h('tbody', null, rows)),
        h('p', { class: 'muted', style: { margin: 0, fontSize: '12px' } }, 'Standardised mean difference of each PC between the groups; |SMD| < 0.1 is the usual bar for balance.'),
        h('div', { style: { display: 'flex', gap: '6px' } },
          h('button', { class: 'btn small', onclick: () => { setFilters({ [res.name]: [res.label_a, res.label_b] }); toast('Cohort set to the matched samples'); } }, 'Use as cohort'),
          h('button', { class: 'btn small', onclick: () => { setFilters({ [res.name]: [res.label_a, res.label_b] }); navigate('tracks', { mode: 'groups', groupField: res.name }); } }, 'Group means in Tracks →'))));
    }
  }

  function renderDerived() {
    const sets = derived.filter((d) => d.kind === 'match');
    if (!sets.length) return;
    matchBody.appendChild(h('div', { style: { marginTop: '14px' } }, h('div', { class: 'label', style: { marginBottom: '6px' } }, 'Saved matched sets'),
      h('table', { class: 'table compact' }, h('tbody', null, sets.map((d) => h('tr', null,
        h('td', null, h('b', { class: 'mono' }, d.name)),
        h('td', null, Object.entries(d.counts || {}).map(([k, v]) => `${k} ${fmtInt(v)}`).join(' · ')),
        h('td', { class: 'muted' }, `${d.field}, PC1–PC${d.k}${d.ready ? '' : ' · not applied'}`),
        h('td', null, h('div', { style: { display: 'flex', gap: '4px', justifyContent: 'flex-end' } },
          h('button', { class: 'btn small', onclick: () => setFilters({ [d.name]: [d.label_a, d.label_b] }) }, 'Use as cohort'),
          h('button', { class: 'btn small ghost', onclick: async () => {
            if (!confirm(`Delete the matched set ${d.name}?`)) return;
            try { await api(`${ds()}/scalars/delete`, { method: 'POST', body: { name: d.name } }); await reloadSamples(); } catch (err) { toast(err.message, 'error'); }
          } }, 'Delete')))))))));
  }

  function figureOptions() {
    if (!pca) return null;
    return {
      root: pcaCard,
      redraw: () => redraw(),
      title: `Genotype PCA · ${pcLabel(cfg.x)} vs ${pcLabel(cfg.y)}`,
      caption: [
        `${fmtInt(pca.n_sites)} SNVs (MAF ≥ ${pca.min_maf}, ≥ ${fmtInt(pca.spacing)} bp apart) from ${pca.genes.length} windows: ${pca.genes.join(', ')}. Coloured by ${cfg.color || 'nothing'}; cohort: ${cohortDescription()}.`,
        `Dataset ${state.datasetId} · exported ${new Date().toISOString().slice(0, 10)} from genomics visualize`,
      ],
      filename: slug(`${state.datasetId}_pca_pc${cfg.x + 1}_pc${cfg.y + 1}`),
    };
  }

  disposers.push(onResize(scroll, () => redraw()));
  init();
  return {
    refresh() { redraw(); },
    unmount() { latest.abort(); hideTooltip(); for (const d of disposers) d(); },
  };
}
