// Perturbation Lab: edit an individual's haplotypes in silico, re-predict them with AlphaGenome and
// re-score them with a trained model, all inside the visualizer.
//
// Edits are made on the genomic axis (drag on the sequence lane) and applied to each chosen
// haplotype through its own indels; overwrite / scramble / revert-to-reference / custom sequence
// all keep the haplotype length, so edited predictions line up with the originals.
import { api, apiJob, isAbort, Latest } from '../api.js';
import { navigate, updateRouteParams, watchJobs } from '../app.js';
import { state, setDataset, emit, loadStatus } from '../state.js';
import { h, clear, icon, iconButton, segmented, select, setOptions, field, checkbox, fmtInt, fmtBp, fmtNum, fmtPct, toast, jobOverlay, debounce, errorBox, showTooltip, hideTooltip, escapeHtml } from '../ui.js';
import { Viewport, setupCanvas, theme, ticks, withAlpha, onResize, css } from '../plot.js';
import { loadGeneAnnotations, packGenes, drawGeneLanes, geneLanesHeight, geneAt, geneTooltip } from '../gene_lanes.js';
import { labelGeneOptions, openGeneCard } from '../cards.js';

const GUTTER_L = 92;
const GUTTER_R = 14;
const PANEL_H = 120;
const ROW_H = 18;
const LETTER_SPAN = 4000;
const OPS = [
  { value: 'scramble', label: 'Scramble', title: 'Shuffle the bases of the region (keeps composition)' },
  { value: 'overwrite', label: 'Overwrite', title: 'Replace every base of the region with one base' },
  { value: 'reference', label: 'Revert', title: 'Revert the haplotype to the reference in the region (removes SNVs; indels stay)' },
  { value: 'sequence', label: 'Sequence', title: 'Replace the region with your own sequence (same length)' },
];

export async function mount(root, params) {
  const page = new PerturbPage(root);
  await page.init(params);
  return page;
}

class PerturbPage {
  constructor(root) {
    this.root = root;
    this.latest = new Latest();
    this.seqLatest = new Latest();
    this.disposers = [];
    this.model = null;
    this.info = null;
    this.edits = [];
    this.result = null;
    this.baseline = null;
    this.data = null;
    this.means = null;
    this.seq = null;
    this.selection = null; // [start, end) reference offsets
    this.annotations = null;
    this.geneRows = [];
    this.cfg = { sample: null, gene: null, output: null, tracks: null, haps: 'both', editHaps: ['H1', 'H2'], op: 'scramble', base: 'A', seed: 1, sequence: '', means: false, delta: false, log: false };
  }

  // ------------------------------------------------------------------ setup
  async init(params) {
    this.page = h('div', { class: 'workspace' });
    this.root.appendChild(this.page);
    clear(this.page).appendChild(h('div', { class: 'page-inner' }, h('div', { class: 'muted' }, 'Loading models…')));
    let models;
    try { models = await api('/api/perturb/models'); } catch (err) { clear(this.page).appendChild(h('div', { class: 'page-inner' }, errorBox(err))); return; }
    this.models = models;
    const compatible = models.runs.filter((r) => r.compatible);
    if (!compatible.length) {
      clear(this.page).appendChild(h('div', { class: 'page-inner' }, h('div', { class: 'page-head' }, h('h1', null, 'Perturbation Lab')),
        h('div', { class: 'card empty' }, 'No trained model found under the runs roots. Train one from the Experiments or Jobs page (NN, CNN or CNN2).')));
      return;
    }
    const loaded = models.loaded;
    const runId = params.run || (loaded && loaded.run) || models.default || compatible[0].id;
    this.buildLayout(runId);
    if (loaded && loaded.run === runId) await this.loadModel(runId, loaded.checkpoint, params);
    else this.renderIdle();
  }

  buildLayout(runId) {
    clear(this.page);
    const runs = this.models.runs;
    this.runSelect = select(runs.map((r) => ({ value: r.id, label: `${r.compatible ? '' : '⚠ '}${r.name}`, disabled: !r.compatible })), runId, () => this.onRunPicked(), { style: { maxWidth: '360px' }, 'aria-label': 'Model run' });
    this.ckptSelect = select([], null, () => {}, { 'aria-label': 'Checkpoint' });
    this.loadBtn = h('button', { class: 'btn primary', onclick: () => this.loadModel(this.runSelect.value, this.ckptSelect.value) }, 'Load model');
    this.sampleInput = h('input', { class: 'input mono', list: 'perturb-samples', placeholder: 'sample id', style: { width: '150px' }, 'aria-label': 'Sample' });
    this.sampleList = h('datalist', { id: 'perturb-samples' });
    this.sampleInput.addEventListener('change', () => this.setSample(this.sampleInput.value.trim()));
    this.geneSelect = select([], null, (v) => this.setGene(v), { 'aria-label': 'Gene' });
    this.outputSelect = select([], null, (v) => { this.cfg.output = v; this.cfg.tracks = null; this.renderSide(); this.invalidate(); }, { 'aria-label': 'Output' });
    this.hapSeg = segmented([{ value: 'H1', label: 'H1' }, { value: 'H2', label: 'H2' }, { value: 'both', label: 'H1 & H2' }], this.cfg.haps, (v) => { this.cfg.haps = v; this.invalidate(); this.requestSequence(); });
    this.labelEl = h('span', { class: 'muted', style: { fontSize: '12px' } });
    const toolbar = h('div', { class: 'ws-toolbar' },
      field('Model', h('div', { style: { display: 'flex', gap: '6px' } }, this.runSelect, this.ckptSelect, this.loadBtn)),
      h('div', { class: 'divider' }),
      field('Individual', h('div', { style: { display: 'flex', gap: '6px', alignItems: 'center' } }, this.sampleInput, this.sampleList, this.labelEl)),
      field('Gene', h('div', { style: { display: 'flex', gap: '4px', alignItems: 'center' } }, this.geneSelect, iconButton('info', 'Gene card: HGNC names, database links and Gene Ontology', () => this.geneSelect.value && openGeneCard(this.geneSelect.value), 'icon-btn bordered'))), field('Track output', this.outputSelect), field('Show', this.hapSeg));
    this.locusInput = h('input', { class: 'input locus-input', 'aria-label': 'Locus', spellcheck: 'false' });
    this.locusInput.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.gotoLocus(this.locusInput.value); });
    this.spanEl = h('span', { class: 'locus-span' });
    this.statusEl = h('span', { class: 'muted', style: { fontSize: '12px' } });
    const locus = h('div', { class: 'locus' },
      iconButton('left', 'Pan left', () => this.viewport.pan(-this.viewport.span * 0.5), 'icon-btn bordered'),
      iconButton('right', 'Pan right', () => this.viewport.pan(this.viewport.span * 0.5), 'icon-btn bordered'),
      iconButton('minus', 'Zoom out', () => this.viewport.zoom(2), 'icon-btn bordered'),
      iconButton('plus', 'Zoom in', () => this.viewport.zoom(0.5), 'icon-btn bordered'),
      this.locusInput, this.spanEl,
      h('button', { class: 'btn small', onclick: () => this.gotoModelWindow() }, icon('target', 14), 'Model window'),
      h('button', { class: 'btn small', onclick: () => this.viewport.set(0, this.viewport.length) }, icon('expand', 14), 'Whole window'),
      this.statusEl,
      h('span', { class: 'hint' }, 'Tracks: drag to pan, ⌘/Ctrl+scroll to zoom · Sequence lane: drag to select a region to edit'),
      iconButton('panel', 'Toggle side panel', () => { this.page.classList.toggle('side-collapsed'); this.render(); }, 'icon-btn bordered'));
    this.rulerHost = h('div', { class: 'canvas-host', style: { borderBottom: '1px solid var(--line)' } }, h('canvas'));
    this.panelsHost = h('div', { class: 'canvas-host' });
    this.seqHost = h('div', { class: 'canvas-host seq-lane', style: { borderTop: '1px solid var(--line)' } }, h('canvas'));
    this.legend = h('div', { class: 'legend' });
    this.loadingLine = h('div', { class: 'loading-line', hidden: true });
    this.crosshair = h('div', { class: 'crosshair', hidden: true });
    this.scroll = h('div', { class: 'ws-scroll' }, this.loadingLine, this.rulerHost, this.legend, h('div', { style: { position: 'relative' } }, this.panelsHost, this.seqHost, this.crosshair));
    this.side = h('aside', { class: 'side', 'aria-label': 'Perturbation' });
    this.page.append(h('div', { class: 'ws-main' }, toolbar, locus, this.scroll), this.side);
    this.viewport = new Viewport({ length: 1, start: 0, end: 1, minSpan: 20, onChange: () => this.onViewChange() });
    for (const el of [this.rulerHost, this.panelsHost]) {
      this.viewport.attach(el, () => this.geom(), { onHover: (e, f) => this.onHover(e, f), onLeave: () => this.clearHover() });
    }
    this.attachSequenceLane();
    this.disposers.push(onResize(this.scroll, () => this.render()));
    const onTheme = () => this.render();
    window.addEventListener('themechange', onTheme);
    this.disposers.push(() => window.removeEventListener('themechange', onTheme));
    this.onRunPicked();
  }

  onRunPicked() {
    const run = this.models.runs.find((r) => r.id === this.runSelect.value);
    setOptions(this.ckptSelect, run ? run.checkpoints : [], run ? run.default_checkpoint : null);
    const isLoaded = this.model && this.model.run === this.runSelect.value && this.model.checkpoint === this.ckptSelect.value;
    this.loadBtn.textContent = isLoaded ? 'Loaded' : 'Load model';
    this.loadBtn.disabled = !!isLoaded;
  }

  renderIdle() {
    clear(this.side).appendChild(h('div', { class: 'side-section' },
      h('h3', null, 'Perturbation Lab'),
      h('p', { class: 'help', style: { margin: 0 } }, 'Pick a trained run and load it. Then choose an individual and a gene, select a region on the sequence lane, add edits (scramble, overwrite, revert to reference or your own sequence) and run AlphaGenome + the model to see how the class probabilities and the predicted tracks change.'),
      this.models.training.length ? h('div', { class: 'notice' }, `Training running: ${this.models.training.join(', ')}. Loading a model also uses memory.`) : null,
      h('div', { class: this.models.backend.reasons.length ? 'error-box' : 'notice' }, h('b', null, 'AlphaGenome: '), this.models.backend.label, this.models.backend.reasons.length ? ` — ${this.models.backend.reasons.join('; ')}` : '', ' · ', h('a', { href: '#/alphagenome' }, 'settings'))));
    this.render();
  }

  async loadModel(runId, checkpoint, params = {}) {
    const overlay = jobOverlay(this.scroll, 'Loading the model');
    this.loadBtn.disabled = true;
    try {
      this.model = await apiJob('/api/perturb/model', { params: { run: runId, checkpoint }, onProgress: (job) => overlay.update(job) });
    } catch (err) {
      if (!isAbort(err)) { toast(err.message, 'error', 10000); this.loadBtn.disabled = false; }
      overlay.remove();
      return;
    }
    overlay.remove();
    this.onRunPicked();
    if (state.datasetId !== this.model.dataset_id) {
      await loadStatus();
      emit('datasets');
      await setDataset(this.model.dataset_id);
      emit('datasets');
      toast(`Switched to the model's dataset: ${this.model.dataset_path}`);
    }
    this.byId = new Map(this.model.samples.map((s) => [s.id, s]));
    clear(this.sampleList);
    for (const s of this.model.samples) if (s.label) this.sampleList.appendChild(h('option', { value: s.id }, `${s.label}${s.split ? ` · ${s.split}` : ''}`));
    setOptions(this.geneSelect, this.model.genes, this.model.genes[0]);
    labelGeneOptions(this.geneSelect);
    const sample = [params.sample, ...state.pinned].find((s) => s && this.byId.has(s) && this.byId.get(s).label)
      || (this.model.samples.find((s) => s.split === 'test' && s.label) || this.model.samples.find((s) => s.label) || {}).id;
    this.cfg.sample = sample;
    this.sampleInput.value = sample || '';
    await this.setGene(params.gene && this.model.genes.includes(params.gene) ? params.gene : this.model.genes[0]);
    this.setSample(sample, true);
  }

  async setSample(sample, quiet = false) {
    if (!this.model) return;
    if (!this.byId.has(sample)) { if (!quiet) toast(`${sample} is not in the model's dataset`, 'error'); return; }
    if (sample !== this.cfg.sample) { this.edits = []; this.result = null; }
    this.cfg.sample = sample;
    this.sampleInput.value = sample;
    const s = this.byId.get(sample);
    this.labelEl.textContent = s.label ? `${s.label}${s.split ? ` · ${s.split} split` : ' · not used in training'}` : 'no label';
    this.baseline = null;
    this.renderSide();
    this.invalidate();
    this.requestSequence();
    this.persist();
    try {
      this.baseline = await api('/api/perturb/score', { params: { sample } });
    } catch (err) { toast(err.message, 'error'); }
    this.renderSide();
  }

  async setGene(gene) {
    this.cfg.gene = gene;
    this.geneSelect.value = gene;
    this.edits = []; // edits are positions in one gene's window
    this.result = null;
    this.selection = null;
    this.means = null;
    this.annotations = null;
    try {
      this.info = await api(`/api/d/${encodeURIComponent(this.model.dataset_id)}/genes/${encodeURIComponent(gene)}`);
    } catch (err) { toast(err.message, 'error'); return; }
    loadGeneAnnotations(this.model.dataset_id, gene).then((res) => { if (gene === this.cfg.gene) { this.annotations = res; this.render(); } }).catch(() => {});
    const outputs = Object.keys(this.info.outputs || {});
    if (!outputs.includes(this.cfg.output)) this.cfg.output = outputs.includes(this.model.outputs[0]) ? this.model.outputs[0] : outputs[0];
    setOptions(this.outputSelect, outputs, this.cfg.output);
    this.cfg.tracks = null;
    this.viewport.setDomain(this.info.length || 1);
    this.setModelWindow(false);
    this.renderSide();
    this.invalidate();
    this.requestSequence();
    this.persist();
  }

  persist() {
    if (!this.model) return;
    updateRouteParams({ run: this.model.run, sample: this.cfg.sample, gene: this.cfg.gene });
  }

  // ------------------------------------------------------------------ geometry & locus
  geom() {
    const width = this.scroll.clientWidth;
    return { left: GUTTER_L, width: Math.max(50, width - GUTTER_L - GUTTER_R), total: width };
  }

  xOf(pos) { const g = this.geom(); const v = this.viewport; return g.left + ((pos - v.start) / v.span) * g.width; }

  modelWindow() {
    if (!this.info || !this.model) return null;
    const len = this.info.length || 0;
    const size = Math.min(this.model.window_center_size, len);
    let start = Math.max(0, Math.floor(len / 2) - Math.floor(size / 2));
    const end = Math.min(len, start + size);
    start = Math.max(0, end - size);
    return { start, end };
  }

  setModelWindow(emitChange = true) {
    const mw = this.modelWindow();
    if (mw) this.viewport.set(mw.start, mw.end, emitChange);
  }

  gotoModelWindow() { this.setModelWindow(true); }

  formatPos(pos) { return this.info && this.info.start ? fmtInt(this.info.start + pos) : fmtInt(pos + 1); }

  renderLocus() {
    if (!this.info) return;
    const v = this.viewport;
    if (document.activeElement !== this.locusInput) this.locusInput.value = `${this.info.chromosome}:${fmtInt(this.info.start + v.start)}-${fmtInt(this.info.start + v.end - 1)}`;
    this.spanEl.textContent = `${fmtBp(v.span)} of ${fmtBp(v.length)}`;
  }

  gotoLocus(text) {
    const m = /^(?:[\w.]+:)?\s*(\d+)(?:\s*[-–]\s*(\d+))?$/.exec(text.replace(/,/g, '').trim());
    if (!m || !this.info) { toast('Use chr:start-end or a position', 'error'); this.renderLocus(); return; }
    const conv = (x) => (Number(x) > this.info.start ? Number(x) - this.info.start : Number(x) - 1);
    if (!m[2]) this.viewport.center(conv(m[1]));
    else this.viewport.set(conv(m[1]), conv(m[2]) + 1);
    this.locusInput.blur();
  }

  status(text) { this.statusEl.textContent = text || ''; }

  onViewChange() {
    this.renderLocus();
    this.render();
    this.request();
    this.requestSequence();
  }

  // ------------------------------------------------------------------ data
  haps() { return this.cfg.haps === 'both' ? ['H1', 'H2'] : [this.cfg.haps]; }

  trackMeta() { return ((this.info && this.info.outputs[this.cfg.output]) || { tracks: [] }).tracks; }

  defaultTracks() {
    const meta = this.trackMeta();
    const terms = new Set(this.model ? this.model.ontology_terms : []);
    const used = meta.filter((m) => terms.has((m.metadata || {}).ontology_curie)).map((m) => m.index);
    return (used.length ? used : meta.map((m) => m.index)).slice(0, 6);
  }

  invalidate() {
    this.data = null;
    this.means = null;
    this.request.cancel?.();
    this.requestNow();
  }

  request = debounce(() => this.requestNow(), 140);

  requestKey() { return JSON.stringify([this.cfg.sample, this.cfg.gene, this.cfg.output, this.cfg.tracks, this.haps(), this.result && this.result.key, this.cfg.means]); }

  needsFetch() {
    const d = this.data;
    if (!d || d.key !== this.requestKey()) return true;
    const v = this.viewport;
    if (v.start < d.start || v.end > d.end) return true;
    const binBp = (d.end - d.start) / Math.max(1, d.edges.length - 1);
    return binBp > Math.max(1, v.span / this.geom().width) * 1.6;
  }

  async requestNow() {
    if (!this.model || !this.info || !this.cfg.sample || !this.cfg.output) { this.render(); return; }
    if (!this.cfg.tracks) this.cfg.tracks = this.defaultTracks();
    if (!this.needsFetch()) { this.render(); return; }
    const v = this.viewport;
    const g = this.geom();
    const margin = v.span;
    const start = Math.max(0, Math.floor(v.start - margin));
    const end = Math.min(v.length, Math.ceil(v.end + margin));
    const bins = Math.max(32, Math.min(8192, Math.round(g.width * (end - start) / v.span)));
    const signal = this.latest.next();
    const key = this.requestKey();
    const ds = `/api/d/${encodeURIComponent(this.model.dataset_id)}`;
    const common = { gene: this.cfg.gene, output: this.cfg.output, coords: 'reference', start, end, bins, tracks: this.cfg.tracks.join(',') };
    this.loadingLine.hidden = false;
    try {
      const base = await api(`${ds}/signal`, { params: { ...common, series: this.haps().map((hp) => `${this.cfg.sample}:${hp}`).join(',') }, signal });
      const items = base.series.filter((s) => !s.error).map((s) => ({ hap: s.haplotype, kind: 'baseline', mean: s.mean, min: s.min, max: s.max }));
      const result = this.result;
      if (result && result.sample === this.cfg.sample && result.gene === this.cfg.gene && result.outputs.includes(this.cfg.output)) {
        const edited = await api('/api/perturb/signal', { params: { key: result.key, output: this.cfg.output, tracks: this.cfg.tracks.join(','), start, end, bins }, signal });
        for (const s of edited.series) if (s.kind === 'edited' && this.haps().includes(s.haplotype)) items.push({ hap: s.haplotype, kind: 'edited', mean: s.mean, min: s.min, max: s.max });
      }
      this.data = { key, start: base.start, end: base.end, edges: base.edges, tracks: base.tracks, items };
      if (this.cfg.means) await this.requestMeans(common, signal);
      this.status('');
    } catch (err) {
      if (isAbort(err)) return;
      toast(err.message, 'error', 8000);
      this.status('Request failed');
    } finally {
      if (!signal.aborted) this.loadingLine.hidden = true;
    }
    this.render();
  }

  async requestMeans(common, signal) {
    const ds = `/api/d/${encodeURIComponent(this.model.dataset_id)}`;
    let overlay = null;
    try {
      const res = await apiJob(`${ds}/groups`, {
        params: { ...common, haps: this.haps().join(','), field: this.model.label_field, groups: this.model.classes.join(',') },
        signal,
        onProgress: (job) => { if (!overlay) { overlay = jobOverlay(this.scroll, job.title || 'Class means'); watchJobs(); } overlay.update(job); },
      });
      this.means = { start: res.start, end: res.end, edges: res.edges, groups: res.groups };
    } finally { if (overlay) overlay.remove(); }
  }

  // ------------------------------------------------------------------ sequence lane
  requestSequence = debounce(() => this.requestSequenceNow(), 120);

  async requestSequenceNow() {
    if (!this.model || !this.info || !this.cfg.sample) return;
    const v = this.viewport;
    const span = v.span;
    const start = Math.max(0, Math.floor(v.start - (span <= LETTER_SPAN / 3 ? span : 0)));
    const end = Math.min(v.length, Math.ceil(v.end + (span <= LETTER_SPAN / 3 ? span : 0)));
    const signal = this.seqLatest.next();
    try {
      this.seq = await api('/api/perturb/sequence', { method: 'POST', body: { sample: this.cfg.sample, gene: this.cfg.gene, haplotypes: this.haps(), edits: this.edits, start, end, bins: Math.max(32, Math.round(this.geom().width)) }, signal });
      this.drawSequence();
    } catch (err) {
      if (!isAbort(err)) toast(err.message, 'error', 7000);
    }
  }

  attachSequenceLane() {
    const el = this.seqHost;
    el.style.cursor = 'crosshair';
    let drag = null;
    const posAt = (clientX) => {
      const rect = el.getBoundingClientRect();
      const g = this.geom();
      const frac = (clientX - rect.left - g.left) / g.width;
      return Math.max(0, Math.min(this.viewport.length, this.viewport.start + frac * this.viewport.span));
    };
    el.addEventListener('pointerdown', (e) => {
      if (e.button !== 0 || !this.info) return;
      el.setPointerCapture(e.pointerId);
      drag = { a: posAt(e.clientX) };
    });
    el.addEventListener('pointermove', (e) => {
      if (drag) {
        const b = posAt(e.clientX);
        this.selection = [Math.floor(Math.min(drag.a, b)), Math.ceil(Math.max(drag.a, b))];
        this.render();
        return;
      }
      this.onHover(e, (e.clientX - el.getBoundingClientRect().left - this.geom().left) / this.geom().width);
    });
    const finish = (e) => {
      if (!drag) return;
      const b = posAt(e.clientX);
      const lo = Math.floor(Math.min(drag.a, b)); const hi = Math.ceil(Math.max(drag.a, b));
      drag = null;
      this.selection = hi - lo >= 1 ? [lo, hi] : null;
      this.renderSide();
      this.render();
    };
    el.addEventListener('pointerup', finish);
    el.addEventListener('pointercancel', finish);
    el.addEventListener('pointerleave', () => this.clearHover());
    el.addEventListener('wheel', (e) => {
      if (!(e.ctrlKey || e.metaKey)) return;
      e.preventDefault();
      const g = this.geom();
      const frac = (e.clientX - el.getBoundingClientRect().left - g.left) / g.width;
      this.viewport.zoom(Math.exp(Math.max(-1, Math.min(1, e.deltaY * 0.01))), Math.max(0, Math.min(1, frac)));
    }, { passive: false });
  }

  // ------------------------------------------------------------------ edits & prediction
  addEdit() {
    if (!this.selection) { toast('Drag on the sequence lane to select a region first', 'error'); return; }
    const [start, end] = this.selection;
    const edit = { op: this.cfg.op, start, end, haplotypes: [...this.cfg.editHaps] };
    if (edit.op === 'overwrite') edit.base = this.cfg.base;
    if (edit.op === 'scramble') edit.seed = Number(this.cfg.seed) || 0;
    if (edit.op === 'sequence') {
      edit.sequence = this.cfg.sequence.replace(/\s+/g, '').toUpperCase();
      if (!edit.sequence) { toast('Type the replacement sequence', 'error'); return; }
    }
    this.edits = [...this.edits, edit];
    this.renderSide();
    this.requestSequence();
    this.render();
  }

  removeEdit(i) {
    this.edits = this.edits.filter((_, k) => k !== i);
    this.renderSide();
    this.requestSequence();
    this.render();
  }

  async run() {
    if (!this.edits.length) { toast('Add at least one edit', 'error'); return; }
    const body = { sample: this.cfg.sample, gene: this.cfg.gene, edits: this.edits, outputs: [this.cfg.output] };
    this.running = true;
    this.renderSide();
    const overlay = jobOverlay(this.scroll, 'Re-predicting with AlphaGenome', (jobId) => api(`/api/jobs/${jobId}/cancel`, { method: 'POST', body: {} }).catch(() => {}));
    try {
      let res = await api('/api/perturb/apply', { method: 'POST', body });
      if (res.pending) {
        overlay.update(res.job);
        res = await apiJob('/api/perturb/result', { params: { key: res.key }, onProgress: (job) => overlay.update(job) });
      }
      this.result = res;
      this.resultEdits = JSON.stringify(this.edits);
      toast(`Done in ${fmtNum(res.elapsed, 3)} s`);
    } catch (err) {
      if (!isAbort(err)) toast(err.message, 'error', 10000);
    } finally {
      overlay.remove();
      this.running = false;
    }
    this.renderSide();
    this.invalidate();
  }

  // ------------------------------------------------------------------ side panel
  renderSide() {
    if (!this.model) { this.renderIdle(); return; }
    clear(this.side);
    this.side.append(this.predictionSection(), this.editSection(), this.displaySection());
  }

  predictionSection() {
    const m = this.model;
    const base = this.baseline;
    const r = this.result && this.result.sample === this.cfg.sample && this.result.gene === this.cfg.gene ? this.result : null;
    const stale = r && this.resultEdits !== JSON.stringify(this.edits);
    const rows = m.classes.map((cls, i) => {
      const b = base ? base.probabilities[i] : null;
      const e = r ? r.edited[i] : null;
      const d = r ? r.delta[i] : null;
      return h('div', { class: 'prob-row' },
        h('div', { class: 'prob-name', title: cls }, cls),
        h('div', { class: 'prob-bars' },
          h('div', { class: 'prob-bar base', style: { width: `${Math.round((b ?? 0) * 100)}%` }, title: 'original' }),
          r ? h('div', { class: 'prob-bar edited', style: { width: `${Math.round(e * 100)}%` }, title: 'edited' }) : null),
        h('div', { class: 'prob-val' }, b === null ? '…' : fmtPct(b), r ? h('span', null, ' → ', h('b', null, fmtPct(e))) : null,
          r ? h('div', { class: d >= 0 ? 'delta up' : 'delta down' }, `${d >= 0 ? '+' : ''}${(d * 100).toFixed(1)} pp`) : null));
    });
    const sample = this.byId && this.byId.get(this.cfg.sample);
    return h('div', { class: 'side-section' },
      h('h3', null, 'Model prediction', h('span', { class: 'muted', style: { fontWeight: '500', fontSize: '11.5px' } }, m.target)),
      sample ? h('div', { style: { fontSize: '12px' } }, h('b', null, this.cfg.sample), ` · true class: ${sample.label || '–'}${sample.split ? ` · ${sample.split}` : ''}`) : null,
      h('div', { class: 'prob-list' }, rows),
      h('div', { class: 'legend', style: { padding: '0', border: '0', background: 'none' } },
        h('span', { class: 'legend-item' }, h('span', { class: 'swatch', style: { background: 'var(--ink-3)' } }), 'original'),
        h('span', { class: 'legend-item' }, h('span', { class: 'swatch', style: { background: 'var(--series-2)' } }), 'edited')),
      r ? h('div', { class: 'help' }, `Predicted: ${r.predicted_baseline} → ${r.predicted_edited}. ${r.edits.filter((e) => e.inside_model_window).length ? '' : 'No edit touches the CNN window, so only long-range AlphaGenome effects can change the input. '}`) : null,
      stale ? h('div', { class: 'notice' }, 'Edits changed since this run; run again to update.') : null);
  }

  editSection() {
    const sel = this.selection;
    const opSeg = segmented(OPS, this.cfg.op, (v) => { this.cfg.op = v; this.renderSide(); });
    const opts = [];
    if (this.cfg.op === 'overwrite') opts.push(field('Base', segmented(['A', 'C', 'G', 'T', 'N'].map((b) => ({ value: b, label: b })), this.cfg.base, (v) => { this.cfg.base = v; })));
    if (this.cfg.op === 'scramble') { const i = h('input', { class: 'input', type: 'number', value: this.cfg.seed, style: { width: '90px' } }); i.addEventListener('input', () => { this.cfg.seed = i.value; }); opts.push(field('Seed', i)); }
    if (this.cfg.op === 'sequence') { const t = h('textarea', { class: 'input mono', rows: 3, placeholder: 'ACGT… (same length as the selected haplotype region)' }); t.value = this.cfg.sequence; t.addEventListener('input', () => { this.cfg.sequence = t.value; }); opts.push(t); }
    const hapSeg = segmented([{ value: 'H1', label: 'H1' }, { value: 'H2', label: 'H2' }, { value: 'both', label: 'Both' }], this.cfg.editHaps.length === 2 ? 'both' : this.cfg.editHaps[0], (v) => { this.cfg.editHaps = v === 'both' ? ['H1', 'H2'] : [v]; });
    const list = this.edits.map((e, i) => h('div', { class: 'edit-item' },
      h('div', null, h('b', null, `${i + 1}. ${OPS.find((o) => o.value === e.op).label}`), e.op === 'overwrite' ? ` → ${e.base}` : e.op === 'scramble' ? ` (seed ${e.seed})` : '',
        h('div', { class: 'muted mono', style: { fontSize: '11px' } }, `${this.info ? `${this.info.chromosome}:${this.formatPos(e.start)}-${this.formatPos(e.end - 1)}` : `${e.start}-${e.end}`} · ${fmtBp(e.end - e.start)} · ${e.haplotypes.join('+')}`)),
      h('div', { style: { display: 'flex', gap: '2px' } },
        iconButton('target', 'Show', () => this.viewport.center((e.start + e.end) / 2, Math.max(60, (e.end - e.start) * 3))),
        iconButton('close', 'Remove', () => this.removeEdit(i)))));
    const backend = this.models.backend;
    return h('div', { class: 'side-section' },
      h('h3', null, 'Edits', this.edits.length ? h('button', { class: 'btn small ghost', onclick: () => { this.edits = []; this.renderSide(); this.requestSequence(); this.render(); } }, 'Clear') : null),
      h('div', { class: 'notice', style: { fontSize: '12px' } }, sel ? h('span', null, 'Selection: ', h('span', { class: 'mono' }, `${this.info.chromosome}:${this.formatPos(sel[0])}-${this.formatPos(sel[1] - 1)}`), ` (${fmtBp(sel[1] - sel[0])})`) : 'Drag on the sequence lane (below the tracks) to select a region.'),
      field('Operation', opSeg), ...opts, field('Haplotypes', hapSeg),
      h('div', { style: { display: 'flex', gap: '6px' } },
        h('button', { class: 'btn', disabled: !sel, onclick: () => this.addEdit() }, '+ Add edit'),
        h('button', { class: 'btn ghost', title: 'Select the whole visible region', onclick: () => { this.selection = [this.viewport.start, this.viewport.end]; this.renderSide(); this.render(); } }, 'Select view')),
      list.length ? h('div', { class: 'edit-list' }, list) : null,
      h('button', { class: 'btn primary', disabled: !this.edits.length || this.running || backend.reasons.length > 0, onclick: () => this.run() }, icon('play', 14), this.running ? 'Running…' : 'Run AlphaGenome + model'),
      backend.reasons.length ? h('div', { class: 'error-box', style: { fontSize: '12px' } }, `AlphaGenome: ${backend.reasons.join('; ')} · `, h('a', { href: '#/alphagenome' }, 'settings')) : h('p', { class: 'help' }, `Each edited haplotype is one AlphaGenome call (${backend.label}); results are cached.`));
  }

  displaySection() {
    const meta = this.trackMeta();
    const chosen = new Set(this.cfg.tracks || []);
    const list = h('div', { class: 'checklist', style: { maxHeight: '180px' } }, meta.map((m) => {
      const box = h('input', { type: 'checkbox' });
      box.checked = chosen.has(m.index);
      box.addEventListener('change', () => {
        const set = new Set(this.cfg.tracks);
        if (box.checked) set.add(m.index); else set.delete(m.index);
        this.cfg.tracks = [...set].sort((a, b) => a - b).slice(0, 12);
        this.invalidate();
      });
      return h('label', { title: m.label }, box, h('span', null, m.short || m.label), h('span', { class: 'meta' }, `#${m.index}`));
    }));
    return h('div', { class: 'side-section' },
      h('h3', null, 'Display'),
      list,
      checkbox(`Class means (${this.model.classes.join(' vs ')})`, this.cfg.means, (v) => { this.cfg.means = v; this.invalidate(); }),
      checkbox('Edited − original', this.cfg.delta, (v) => { this.cfg.delta = v; this.render(); }),
      checkbox('log(1+x) scale', this.cfg.log, (v) => { this.cfg.log = v; this.render(); }),
      h('p', { class: 'help' }, 'Blue band: the CNN input window. Class means are computed over every haplotype of each class once and cached.'));
  }

  // ------------------------------------------------------------------ rendering
  render() {
    if (!this.info) { this.drawEmpty(); return; }
    this.renderLocus();
    this.drawRuler();
    this.drawPanels();
    this.drawSequence();
    this.renderLegend();
  }

  drawEmpty() {
    const g = this.geom();
    const ctx = setupCanvas(this.rulerHost.querySelector('canvas'), g.total, 30);
    const t = theme();
    ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.textAlign = 'center';
    ctx.fillText(this.model ? 'Loading…' : 'Load a model to start', g.left + g.width / 2, 20);
    clear(this.panelsHost);
  }

  drawRuler() {
    const g = this.geom();
    const v = this.viewport;
    this.geneRows = packGenes((this.annotations && this.annotations.genes) || [], { start: v.start, end: v.end, pxPerBp: g.width / v.span });
    const H = 34 + geneLanesHeight(this.geneRows);
    const ctx = setupCanvas(this.rulerHost.querySelector('canvas'), g.total, H);
    const t = theme();
    this.drawBands(ctx, t, g, 0, H);
    ctx.font = `11px ${t.font}`; ctx.textAlign = 'center'; ctx.fillStyle = t.ink2;
    for (const tv of ticks(v.start + this.info.start, v.end + this.info.start, Math.max(3, Math.floor(g.width / 120)))) {
      const x = Math.round(this.xOf(tv - this.info.start)) + 0.5;
      if (x < g.left || x > g.left + g.width) continue;
      ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(x, 24); ctx.lineTo(x, 33); ctx.stroke();
      ctx.fillText(fmtInt(tv), x, 18);
    }
    ctx.textAlign = 'right'; ctx.fillStyle = t.ink3; ctx.fillText(this.info.chromosome || '', g.left - 8, 18);
    ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(g.left, 33.5); ctx.lineTo(g.left + g.width, 33.5); ctx.stroke();
    drawGeneLanes(ctx, t, this.geneRows, { top: 34, left: g.left, width: g.width, xOf: (p) => this.xOf(p), highlight: this.cfg.gene });
  }

  /** Model window band, edit regions and the current selection. */
  drawBands(ctx, t, g, top, height) {
    const v = this.viewport;
    const clip = (a, b) => [Math.max(g.left, this.xOf(a)), Math.min(g.left + g.width, this.xOf(b))];
    const mw = this.modelWindow();
    if (mw && mw.end > v.start && mw.start < v.end) { const [a, b] = clip(mw.start, mw.end); ctx.fillStyle = t.modelWindow; ctx.fillRect(a, top, b - a, height); }
    for (const e of this.edits) {
      if (e.end <= v.start || e.start >= v.end) continue;
      const [a, b] = clip(e.start, e.end);
      ctx.fillStyle = withAlpha(css('--series-2'), 0.12); ctx.fillRect(a, top, Math.max(1, b - a), height);
    }
    if (this.selection) {
      const [s0, s1] = this.selection;
      if (s1 > v.start && s0 < v.end) { const [a, b] = clip(s0, s1); ctx.fillStyle = withAlpha(t.accent, 0.16); ctx.fillRect(a, top, Math.max(1, b - a), height); ctx.strokeStyle = t.accent; ctx.lineWidth = 1; ctx.strokeRect(a + 0.5, top + 0.5, Math.max(1, b - a) - 1, height - 1); }
    }
  }

  ensurePanels(n) {
    const host = this.panelsHost;
    if (host.children.length !== n) {
      clear(host);
      for (let i = 0; i < n; i++) host.appendChild(h('div', { class: 'panel' }, h('div', { class: 'lane-label', style: { top: '6px' } }), h('canvas')));
    }
    return [...host.children].map((p) => ({ label: p.querySelector('.lane-label'), canvas: p.querySelector('canvas') }));
  }

  seriesStyle(item) {
    const edited = item.kind === 'edited';
    return { color: edited ? css('--series-2') : css('--ink-2'), dash: item.hap === 'H2', width: edited ? 2 : 1.5 };
  }

  displayItems() {
    const d = this.data;
    if (!d) return [];
    if (!this.cfg.delta) return d.items;
    const out = [];
    for (const it of d.items.filter((x) => x.kind === 'edited')) {
      const base = d.items.find((x) => x.kind === 'baseline' && x.hap === it.hap);
      if (!base) continue;
      out.push({ hap: it.hap, kind: 'delta', mean: it.mean.map((row, ti) => row.map((val, k) => val - base.mean[ti][k])) });
    }
    return out;
  }

  drawPanels() {
    const g = this.geom();
    const tracks = this.cfg.tracks || [];
    const panels = this.ensurePanels(Math.max(1, tracks.length));
    const t = theme();
    const meta = this.trackMeta();
    const d = this.data;
    const items = this.displayItems();
    if (!d || !items.length) {
      panels.forEach((p, i) => {
        const ctx = setupCanvas(p.canvas, g.total, PANEL_H);
        p.label.textContent = tracks.length ? (meta[tracks[i]] || {}).label || '' : 'No tracks';
        ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.textAlign = 'center';
        ctx.fillText(this.cfg.delta && d ? 'Run a perturbation to see the difference' : (this.cfg.sample ? 'Loading…' : 'Choose an individual'), g.left + g.width / 2, PANEL_H / 2);
      });
      return;
    }
    const v = this.viewport;
    const e = d.edges;
    let i0 = 0; let i1 = e.length - 1;
    while (i0 < e.length - 1 && e[i0 + 1] <= v.start) i0++;
    while (i1 > i0 && e[i1 - 1] >= v.end) i1--;
    const tf = (x) => (this.cfg.log ? Math.sign(x) * Math.log1p(Math.abs(x)) : x);
    const means = !this.cfg.delta && this.cfg.means && this.means ? this.means : null;
    d.tracks.forEach((track, ti) => {
      const p = panels[ti];
      if (!p) return;
      p.label.textContent = (meta[track] || {}).label || `Track ${track}`;
      let lo = 0; let hi = 0;
      for (const it of items) for (let k = i0; k < i1; k++) { const val = tf(it.mean[ti][k]); if (val > hi) hi = val; if (val < lo) lo = val; }
      if (means) for (const gr of means.groups) for (const val of gr.mean[ti]) { const x = tf(val); if (x > hi) hi = x; }
      if (hi <= lo) hi = lo + 1;
      hi += (hi - lo) * 0.06;
      const ctx = setupCanvas(p.canvas, g.total, PANEL_H);
      const top = 24; const bottom = PANEL_H - 8;
      const sy = (val) => bottom - ((tf(val) - lo) / (hi - lo)) * (bottom - top);
      this.drawBands(ctx, t, g, top - 4, bottom - top + 4);
      ctx.font = `10.5px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      for (const val of ticks(lo, hi, 3)) {
        const y = Math.round(bottom - ((val - lo) / (hi - lo)) * (bottom - top)) + 0.5;
        if (y < top - 1 || y > bottom + 1) continue;
        ctx.strokeStyle = t.line; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(g.left, y); ctx.lineTo(g.left + g.width, y); ctx.stroke();
        ctx.fillStyle = t.ink3; ctx.fillText(fmtNum(this.cfg.log ? Math.sign(val) * Math.expm1(Math.abs(val)) : val, 2), g.left - 6, y);
      }
      if (lo < 0) { const y0 = Math.round(sy(0)) + 0.5; ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(g.left, y0); ctx.lineTo(g.left + g.width, y0); ctx.stroke(); }
      ctx.save();
      ctx.beginPath(); ctx.rect(g.left, 0, g.width, PANEL_H); ctx.clip();
      const line = (edges, vals, color, dash, width) => {
        ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash ? [5, 4] : []); ctx.lineJoin = 'round';
        ctx.beginPath();
        let pen = false;
        for (let k = 0; k < vals.length; k++) {
          const val = vals[k];
          if (!Number.isFinite(val)) { pen = false; continue; }
          const x = (this.xOf(edges[k]) + this.xOf(edges[k + 1])) / 2;
          if (x < g.left - 50 || x > g.left + g.width + 50) { pen = false; continue; }
          const y = sy(val);
          if (pen) ctx.lineTo(x, y); else { ctx.moveTo(x, y); pen = true; }
        }
        ctx.stroke();
        ctx.setLineDash([]);
      };
      if (means) means.groups.forEach((gr, k) => line(means.edges, gr.mean[ti], withAlpha(css(`--series-${3 + (k % 5)}`), 0.9), true, 1.25));
      for (const it of items) {
        const s = it.kind === 'delta' ? { color: css('--series-2'), dash: it.hap === 'H2', width: 1.75 } : this.seriesStyle(it);
        line(d.edges, it.mean[ti], s.color, s.dash, s.width);
      }
      ctx.restore();
    });
  }

  drawSequence() {
    const g = this.geom();
    const canvas = this.seqHost.querySelector('canvas');
    const t = theme();
    const d = this.seq;
    const rows = d ? d.rows : [];
    const H = 26 + ROW_H * (1 + rows.length * 2) + 10;
    const ctx = setupCanvas(canvas, g.total, H);
    if (!this.info) return;
    const top = 22;
    this.drawBands(ctx, t, g, 0, H);
    ctx.font = `11px ${t.font}`; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillStyle = t.ink3;
    ctx.fillText(d ? (d.mode === 'letters' ? 'Bases (changed bases outlined)' : 'Zoomed out: orange = variants vs reference, accent = edited bases · zoom below 4 kb for letters') : 'Loading sequence…', g.left + 4, 4);
    if (!d) return;
    const v = this.viewport;
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillStyle = t.ink2;
    const labels = ['Reference', ...rows.flatMap((r) => [`${r.haplotype}`, `${r.haplotype} edited`])];
    labels.forEach((label, i) => ctx.fillText(label, g.left - 8, top + i * ROW_H + ROW_H / 2));
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, top, g.width, H - top); ctx.clip();
    if (d.mode === 'letters') {
      const px = g.width / v.span;
      const from = Math.max(v.start, d.start); const to = Math.min(v.end, d.end);
      const showLetters = px >= 7.5;
      ctx.font = `${Math.min(13, Math.max(9, px * 0.75))}px ${t.mono}`; ctx.textAlign = 'center';
      const baseColor = (b) => t.bases[b] || t.bases.N;
      const drawRow = (str, y, ref, changed) => {
        for (let p = from; p < to; p++) {
          const k = p - d.start;
          const b = str[k];
          const x = this.xOf(p);
          if (b === '-') { ctx.fillStyle = t.ink3; ctx.fillRect(x, y + ROW_H / 2 - 1, px, 2); continue; }
          const differs = ref && b !== ref[k];
          if (!ref || differs) { ctx.fillStyle = ref ? baseColor(b) : withAlpha(baseColor(b), 0.22); ctx.fillRect(x + (px > 4 ? 0.5 : 0), y + 1, Math.max(1, px - (px > 4 ? 1 : 0)), ROW_H - 2); }
          if (changed && changed[k]) { ctx.strokeStyle = css('--series-2'); ctx.lineWidth = 2; ctx.strokeRect(x + 1, y + 1.5, Math.max(1, px - 2), ROW_H - 3); }
          if (showLetters) { ctx.fillStyle = differs ? '#fff' : (ref ? t.ink3 : t.ink); ctx.fillText(ref && !differs ? '·' : b, x + px / 2, y + ROW_H / 2 + 0.5); }
        }
      };
      drawRow(d.reference, top, null, null);
      rows.forEach((r, i) => {
        const y = top + ROW_H * (1 + i * 2);
        const changed = Array.from(r.edited, (b, k) => b !== r.original[k]);
        drawRow(r.original, y, d.reference, null);
        drawRow(r.edited, y + ROW_H, d.reference, changed);
      });
    } else {
      const e = d.edges;
      ctx.fillStyle = t.surface2; ctx.fillRect(g.left, top + 3, g.width, ROW_H - 6);
      rows.forEach((r, i) => {
        const y = top + ROW_H * (1 + i * 2);
        for (const [yy, vals, color] of [[y, r.variant_density, css('--series-2')], [y + ROW_H, r.changed_density, t.accent]]) {
          ctx.fillStyle = t.surface2; ctx.fillRect(g.left, yy + 3, g.width, ROW_H - 6);
          for (let k = 0; k < e.length - 1; k++) {
            if (!(vals[k] > 0)) continue;
            const x0 = this.xOf(e[k]); const x1 = this.xOf(e[k + 1]);
            ctx.fillStyle = withAlpha(color, Math.min(1, 0.45 + vals[k] * 20));
            ctx.fillRect(x0, yy + 2, Math.max(1.5, x1 - x0), ROW_H - 4);
          }
        }
      });
    }
    ctx.restore();
  }

  renderLegend() {
    clear(this.legend);
    const add = (label, color, dash) => this.legend.appendChild(h('span', { class: 'legend-item' }, dash ? h('span', { class: 'swatch line dashed', style: { color, borderTopColor: color } }) : h('span', { class: 'swatch line', style: { background: color } }), label));
    if (this.cfg.delta) {
      for (const hp of this.haps()) add(`${hp} edited − original`, css('--series-2'), hp === 'H2');
    } else {
      for (const hp of this.haps()) {
        add(`${this.cfg.sample || ''} ${hp} original`, css('--ink-2'), hp === 'H2');
        if (this.result && this.result.haplotypes.includes(hp)) add(`${hp} edited`, css('--series-2'), hp === 'H2');
      }
      if (this.cfg.means && this.means) this.means.groups.forEach((gr, k) => add(`${gr.group} mean (n=${fmtInt(gr.samples)})`, css(`--series-${3 + (k % 5)}`), true));
    }
    this.legend.appendChild(h('span', { class: 'muted' }, 'H2 dashed · coordinates are genomic (each haplotype remapped through its indels)'));
  }

  // ------------------------------------------------------------------ hover
  clearHover() { this.crosshair.hidden = true; hideTooltip(); }

  onHover(e, frac) {
    if (!this.info || frac < 0 || frac > 1) { this.clearHover(); return; }
    const g = this.geom();
    const pos = Math.floor(this.viewport.start + frac * this.viewport.span);
    this.crosshair.hidden = false;
    this.crosshair.style.left = `${g.left + frac * g.width}px`;
    if (this.rulerHost.contains(e.target) && this.geneRows.length) {
      const rect = this.rulerHost.getBoundingClientRect();
      const gene = geneAt(this.geneRows, e.clientX - rect.left, e.clientY - rect.top, { top: 34, xOf: (p) => this.xOf(p) });
      if (gene) { showTooltip(e.clientX, e.clientY, geneTooltip(gene, escapeHtml)); return; }
    }
    const lines = [`<div class="tt-title">${escapeHtml(this.info.chromosome)}:${fmtInt(this.info.start + pos)}</div>`];
    const d = this.data;
    const panel = e.target.closest ? e.target.closest('.panel') : null;
    if (d && panel) {
      const ti = [...this.panelsHost.children].indexOf(panel);
      let k = -1;
      for (let i = 0; i < d.edges.length - 1; i++) if (d.edges[i] <= pos && pos < d.edges[i + 1]) { k = i; break; }
      if (k >= 0 && ti >= 0) for (const it of this.displayItems()) lines.push(`<div class="tt-row"><span class="tt-key">${it.hap} ${it.kind}</span><span class="tt-val">${fmtNum(it.mean[ti][k])}</span></div>`);
    }
    if (this.seq && this.seq.mode === 'letters' && pos >= this.seq.start && pos < this.seq.end) {
      const k = pos - this.seq.start;
      lines.push(`<div class="tt-row"><span class="tt-key">reference</span><span class="tt-val mono">${this.seq.reference[k]}</span></div>`);
      for (const r of this.seq.rows) lines.push(`<div class="tt-row"><span class="tt-key">${r.haplotype}</span><span class="tt-val mono">${r.original[k]}${r.edited[k] !== r.original[k] ? ` → ${r.edited[k]}` : ''}</span></div>`);
    }
    showTooltip(e.clientX, e.clientY, lines.join(''));
  }

  // ------------------------------------------------------------------ lifecycle
  update(params) {
    if (!this.model) return;
    if (params.sample && params.sample !== this.cfg.sample) this.setSample(params.sample);
    if (params.gene && params.gene !== this.cfg.gene && this.model.genes.includes(params.gene)) this.setGene(params.gene);
  }

  unmount() {
    this.latest.abort();
    this.seqLatest.abort();
    this.request.cancel?.();
    this.requestSequence.cancel?.();
    hideTooltip();
    for (const d of this.disposers) d();
  }
}
