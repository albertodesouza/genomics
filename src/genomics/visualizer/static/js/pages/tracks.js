// Tracks: interactive AlphaGenome signal browser.
//
// Modes:   individuals (pinned samples) · groups (cohort means ± SD by any facet) · population
//          heatmap (every cohort sample as a row, sorted by a facet).
// Coords:  genomic (each haplotype remapped through its own indels) · haplotype (raw prediction
//          index) · training axis (bcftools_chain expanded alignment used for CNN inputs).
import { api, apiJob, isAbort, Latest } from '../api.js';
import { navigate, updateRouteParams, watchJobs } from '../app.js';
import { state, ds, geneInfo, setLocus, setPinned, togglePinned, categoricalFields, filterRows } from '../state.js';
import { h, clear, icon, iconButton, segmented, select, setOptions, field, fmtInt, fmtBp, fmtNum, showTooltip, hideTooltip, toast, jobOverlay, debounce, escapeHtml, errorBox } from '../ui.js';
import { Viewport, ColorSlots, setupCanvas, theme, ticks, withAlpha, seqColor, onResize, css } from '../plot.js';

const GUTTER_L = 64;
const GUTTER_R = 14;
const PANEL_H = 118;
const MAX_INDIVIDUALS = 8;
const COORD_OPTIONS = [
  { value: 'reference', label: 'Genomic', title: 'Reference coordinates: each haplotype is remapped through its own indels (deleted bases are gaps)' },
  { value: 'haplotype', label: 'Haplotype', title: 'Raw prediction index of each haplotype sequence (positions drift after indels)' },
  { value: 'aligned', label: 'Training axis', title: 'bcftools_chain expanded alignment used to build CNN inputs (shared insertion columns)' },
];

export async function mount(root, params) {
  const page = new TracksPage(root);
  await page.init(params);
  return page;
}

class TracksPage {
  constructor(root) {
    this.root = root;
    this.latest = new Latest();
    this.colors = new ColorSlots(8);
    this.hidden = new Set();
    this.data = null;
    this.pop = null;
    this.axis = null;
    this.annotations = null;
    this.hover = null;
    this.disposers = [];
  }

  // ------------------------------------------------------------------ setup
  async init(params) {
    const genes = state.summary.genes.map((g) => g.gene);
    const saved = state.locus || {};
    const gene = (params.gene && genes.includes(params.gene)) ? params.gene : (saved.gene && genes.includes(saved.gene) ? saved.gene : (state.summary.genes.find((g) => g.in_metadata) || state.summary.genes[0] || {}).gene);
    const cats = categoricalFields();
    this.cfg = {
      gene,
      output: params.output || saved.output || null,
      // A URL locus without coords is genomic (coords are omitted from URLs when genomic).
      coords: params.coords || (params.start !== undefined ? 'reference' : saved.coords) || 'reference',
      mode: params.mode || saved.mode || 'individuals',
      hap: saved.hap || 'H1+H2',
      tracks: saved.tracks || null,
      groupField: saved.groupField && cats.some((f) => f.name === saved.groupField) ? saved.groupField : (cats.find((f) => f.name === 'superpopulation') || cats[0] || {}).name,
      groups: saved.groups || null,
      yScale: saved.yScale || 'linear',
      sharedY: saved.sharedY ?? false,
      diff: saved.diff ?? false,
      envelope: saved.envelope ?? true,
      band: saved.band ?? true,
      popTrack: saved.popTrack ?? 0,
      popRows: saved.popRows || 600,
    };
    this.buildLayout();
    this.viewport = new Viewport({ length: 1, start: 0, end: 1, minSpan: 30, onChange: () => this.onViewChange() });
    for (const el of [this.overviewHost, this.rulerHost, this.panelsHost]) {
      const isOverview = el === this.overviewHost;
      if (isOverview) this.attachOverview(el);
      else this.viewport.attach(el, () => this.geom(), {
        onHover: (e, frac) => this.onHover(e, frac),
        onLeave: () => this.clearHover(),
        onSelect: (a, b, done) => this.drawBrush(a, b, done),
      });
    }
    // Clicking (not dragging) a population heatmap row pins that sample.
    let down = null;
    this.panelsHost.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
    this.panelsHost.addEventListener('click', (e) => {
      if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4) return;
      if (this.cfg.mode !== 'population' || !this.hoverRow) return;
      const id = this.hoverRow;
      const wasPinned = state.pinned.includes(id);
      togglePinned(id);
      toast(wasPinned ? `Unpinned ${id}` : `Pinned ${id} — switch to Individuals to compare`);
    });
    this.disposers.push(onResize(this.scroll, () => this.render()));
    const onTheme = () => this.render();
    window.addEventListener('themechange', onTheme);
    this.disposers.push(() => window.removeEventListener('themechange', onTheme));
    await this.setGene(gene, { start: params.start, end: params.end, keepView: !params.gene || params.gene === saved.gene });
  }

  buildLayout() {
    const genes = state.summary.genes.map((g) => ({ value: g.gene, label: g.gene }));
    this.geneSelect = select(genes, this.cfg.gene, (v) => this.setGene(v, {}), { 'aria-label': 'Gene', style: { minWidth: '120px' } });
    this.outputSelect = select([], null, (v) => { this.cfg.output = v; this.cfg.tracks = null; this.applyGeneInfo(); this.invalidate(); }, { 'aria-label': 'Output' });
    this.coordSeg = segmented(COORD_OPTIONS, this.cfg.coords, (v) => this.setCoords(v));
    this.modeSeg = segmented([
      { value: 'individuals', label: 'Individuals', title: 'Pinned individuals as separate lines' },
      { value: 'groups', label: 'Group means', title: 'Mean ± SD over the cohort, split by a sample field' },
      { value: 'population', label: 'Population heatmap', title: 'Every cohort sample as a row (one track)' },
    ], this.cfg.mode, (v) => { this.cfg.mode = v; this.renderSide(); this.invalidate(); });
    this.hapSeg = segmented([
      { value: 'H1', label: 'H1' }, { value: 'H2', label: 'H2' },
      { value: 'both', label: 'H1 & H2', title: 'Both haplotypes as separate lines (H2 dashed)' },
      { value: 'H1+H2', label: 'Mean', title: 'Diploid mean of H1 and H2' },
    ], this.cfg.hap, (v) => { this.cfg.hap = v; this.invalidate(); });

    this.locusInput = h('input', { class: 'input locus-input', 'aria-label': 'Locus', spellcheck: 'false' });
    this.locusInput.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.gotoLocus(this.locusInput.value); });
    this.spanEl = h('span', { class: 'locus-span' });
    this.statusEl = h('span', { class: 'muted', style: { fontSize: '12px' } });

    const toolbar = h('div', { class: 'ws-toolbar' },
      field('Gene', this.geneSelect), field('Output', this.outputSelect),
      field('Coordinates', this.coordSeg),
      h('div', { class: 'divider' }),
      field('Show', this.modeSeg), field('Haplotype', this.hapSeg),
    );
    const locus = h('div', { class: 'locus' },
      iconButton('left', 'Pan left (←)', () => this.viewport.pan(-this.viewport.span * 0.5), 'icon-btn bordered'),
      iconButton('right', 'Pan right (→)', () => this.viewport.pan(this.viewport.span * 0.5), 'icon-btn bordered'),
      iconButton('minus', 'Zoom out (−)', () => this.viewport.zoom(2), 'icon-btn bordered'),
      iconButton('plus', 'Zoom in (+)', () => this.viewport.zoom(0.5), 'icon-btn bordered'),
      this.locusInput, this.spanEl,
      h('button', { class: 'btn small', title: 'Jump to the CNN training window', onclick: () => this.gotoModelWindow() }, icon('target', 14), 'Model window'),
      h('button', { class: 'btn small', title: 'Show the whole prediction window', onclick: () => this.viewport.set(0, this.viewport.length) }, icon('expand', 14), 'Whole window'),
      h('button', { class: 'btn small ghost', title: 'Open this locus on the Sequence page', onclick: () => navigate('sequence', { gene: this.cfg.gene }) }, 'Sequence →'),
      this.statusEl,
      h('span', { class: 'hint' }, 'Drag to pan · ⌘/Ctrl+scroll zoom · Shift+drag region'),
      iconButton('panel', 'Toggle settings panel', () => { this.workspace.classList.toggle('side-collapsed'); this.render(); }, 'icon-btn bordered'));

    this.overviewHost = h('div', { class: 'canvas-host', style: { borderBottom: '1px solid var(--line)' } }, h('canvas'));
    this.rulerHost = h('div', { class: 'canvas-host', style: { borderBottom: '1px solid var(--line)' } }, h('canvas'));
    this.panelsHost = h('div', { class: 'canvas-host' });
    this.crosshair = h('div', { style: { position: 'absolute', top: '0', bottom: '0', width: '1px', background: 'var(--ink-3)', opacity: '.6', pointerEvents: 'none', zIndex: '3' }, hidden: true });
    this.brush = h('div', { style: { position: 'absolute', top: '0', bottom: '0', background: 'var(--accent-wash)', borderLeft: '1px solid var(--accent)', borderRight: '1px solid var(--accent)', pointerEvents: 'none', zIndex: '3' }, hidden: true });
    this.loadingLine = h('div', { class: 'loading-line', hidden: true });
    this.legend = h('div', { class: 'legend' });
    this.footer = h('div', { class: 'muted', style: { padding: '10px 16px', fontSize: '12px' } });
    this.scroll = h('div', { class: 'ws-scroll' }, this.loadingLine, this.overviewHost, this.rulerHost, this.legend, h('div', { style: { position: 'relative' } }, this.panelsHost, this.crosshair, this.brush), this.footer);
    this.side = h('aside', { class: 'side', 'aria-label': 'Track settings' });
    this.workspace = h('div', { class: 'workspace' }, h('div', { class: 'ws-main' }, toolbar, locus, this.scroll), this.side);
    this.root.appendChild(this.workspace);
  }

  geom() {
    const width = this.scroll.clientWidth;
    return { left: GUTTER_L, width: Math.max(50, width - GUTTER_L - GUTTER_R), total: width };
  }

  // ------------------------------------------------------------------ gene / coords
  async setGene(gene, { start, end, keepView } = {}) {
    this.cfg.gene = gene;
    this.geneSelect.value = gene;
    this.data = null; this.pop = null; this.annotations = null; this.axis = null;
    this.status('Loading gene…');
    try {
      this.info = await geneInfo(gene);
    } catch (err) {
      this.showError(err);
      return;
    }
    this.applyGeneInfo();
    if (this.cfg.coords === 'aligned') {
      const ok = await this.loadAxis();
      if (!ok) this.cfg.coords = 'reference';
      this.coordSeg.setValue(this.cfg.coords);
    }
    this.viewport.setDomain(this.domainLength());
    const saved = state.locus || {};
    if (start !== undefined && end !== undefined) this.viewport.set(Number(start), Number(end), false);
    else if (keepView && saved.gene === gene && saved.coords === this.cfg.coords && saved.end > saved.start) this.viewport.set(saved.start, saved.end, false);
    else this.setModelWindow(false);
    this.renderSide();
    this.loadAnnotations();
    this.invalidate();
  }

  applyGeneInfo() {
    const outputs = Object.keys(this.info.outputs || {});
    if (!outputs.length) {
      setOptions(this.outputSelect, [{ value: '', label: 'no predictions' }], '');
      this.cfg.output = null;
      return;
    }
    if (!outputs.includes(this.cfg.output)) this.cfg.output = outputs.includes('rna_seq') ? 'rna_seq' : outputs[0];
    setOptions(this.outputSelect, outputs, this.cfg.output);
    const n = this.trackMeta().length;
    if (!this.cfg.tracks || this.cfg.tracks.some((t) => t >= n)) this.cfg.tracks = Array.from({ length: Math.min(n, 6) }, (_, i) => i);
    if (this.cfg.popTrack >= n) this.cfg.popTrack = 0;
  }

  trackMeta() { return ((this.info && this.info.outputs[this.cfg.output]) || { tracks: [] }).tracks; }

  domainLength() {
    if (this.cfg.coords === 'aligned' && this.axis) return this.axis.expanded_length;
    if (this.cfg.coords === 'haplotype') return (this.info.outputs[this.cfg.output] || {}).length || this.info.length;
    return this.info.length || (this.info.outputs[this.cfg.output] || {}).length || 1;
  }

  async setCoords(coords) {
    const prev = this.cfg.coords;
    const center = this.toGenomicOffset(this.viewport.start + this.viewport.span / 2);
    const span = this.viewport.span;
    this.cfg.coords = coords;
    if (coords === 'aligned' && !this.axis) {
      const ok = await this.loadAxis();
      if (!ok) { this.cfg.coords = prev; this.coordSeg.setValue(prev); return; }
    }
    this.viewport.setDomain(this.domainLength());
    this.data = null; this.pop = null;
    const c = this.fromGenomicOffset(center);
    this.viewport.set(c - span / 2, c + span / 2, false);
    this.invalidate();
  }

  async loadAxis() {
    const overlay = jobOverlay(this.scroll, 'Preparing the training alignment axis');
    try {
      this.axis = await apiJob(`${ds()}/genes/${encodeURIComponent(this.cfg.gene)}/axis`, { onProgress: (job) => overlay.update(job) });
      return true;
    } catch (err) {
      if (!isAbort(err)) toast(`Training axis unavailable: ${err.message}`, 'error', 8000);
      return false;
    } finally { overlay.remove(); }
  }

  /** Reference (window) offset <-> current coordinates. */
  fromGenomicOffset(ref) {
    if (this.cfg.coords !== 'aligned' || !this.axis) return ref;
    const rr = ref - this.axis.ref_start_offset;
    if (rr < 0) return 0;
    const slots = this.axis.insertion_slots;
    let k = 0;
    while (k < slots.length && slots[k] <= rr + k) k++;
    return Math.min(this.axis.expanded_length, rr + k);
  }

  toGenomicOffset(pos) {
    if (this.cfg.coords !== 'aligned' || !this.axis) return pos;
    const slots = this.axis.insertion_slots;
    let lo = 0; let hi = slots.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (slots[m] < pos) lo = m + 1; else hi = m; }
    return this.axis.ref_start_offset + pos - lo;
  }

  isInsertionSlot(pos) {
    if (this.cfg.coords !== 'aligned' || !this.axis) return false;
    const slots = this.axis.insertion_slots;
    let lo = 0; let hi = slots.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (slots[m] < pos) lo = m + 1; else hi = m; }
    return slots[lo] === pos;
  }

  modelWindow() {
    if (this.cfg.coords === 'aligned') return this.axis ? this.axis.model_window : null;
    return this.info.model_window;
  }

  setModelWindow(emit = true) {
    const mw = this.modelWindow();
    if (mw) this.viewport.set(mw.start, mw.end, emit);
    else this.viewport.set(0, this.viewport.length, emit);
  }

  gotoModelWindow() { this.setModelWindow(true); }

  formatPos(pos, withChrom = false) {
    const coords = this.cfg.coords;
    if (coords === 'reference' && this.info.start) {
      const g = this.info.start + pos;
      return withChrom ? `${this.info.chromosome}:${fmtInt(g)}` : fmtInt(g);
    }
    if (coords === 'aligned' && this.info.start) {
      if (this.isInsertionSlot(Math.floor(pos))) return `ins @ ${fmtInt(this.info.start + this.toGenomicOffset(pos))}`;
      const g = this.info.start + this.toGenomicOffset(pos);
      return withChrom ? `${this.info.chromosome}:${fmtInt(g)}` : fmtInt(g);
    }
    return `${fmtInt(pos + 1)}`;
  }

  renderLocus() {
    const v = this.viewport;
    const c = this.cfg.coords;
    if (c === 'reference' && this.info.start) this.locusInput.value = `${this.info.chromosome}:${fmtInt(this.info.start + v.start)}-${fmtInt(this.info.start + v.end - 1)}`;
    else if (c === 'aligned') this.locusInput.value = `axis:${fmtInt(v.start + 1)}-${fmtInt(v.end)}`;
    else this.locusInput.value = `hap:${fmtInt(v.start + 1)}-${fmtInt(v.end)}`;
    this.spanEl.textContent = `${fmtBp(v.span)} of ${fmtBp(v.length)}`;
  }

  gotoLocus(text) {
    const clean = text.replace(/,/g, '').trim();
    const m = /^(?:([\w.]+):)?\s*(\d+)(?:\s*[-–]\s*(\d+))?$/.exec(clean);
    if (!m) { toast('Use chr:start-end, start-end or a single position', 'error'); this.renderLocus(); return; }
    const prefix = m[1];
    let a = Number(m[2]);
    let b = m[3] ? Number(m[3]) : null;
    const genomic = this.info.start && (!prefix || (prefix !== 'hap' && prefix !== 'axis')) && a > (this.info.start - 1);
    const conv = (x) => {
      if (genomic) return this.fromGenomicOffset(x - this.info.start);
      return x - 1;
    };
    if (b === null) { this.viewport.center(conv(a)); return; }
    if (b < a) [a, b] = [b, a];
    this.viewport.set(conv(a), conv(b) + 1);
  }

  status(text) { this.statusEl.textContent = text || ''; }

  showError(err) {
    clear(this.panelsHost).appendChild(h('div', { style: { padding: '16px' } }, errorBox(err)));
  }

  // ------------------------------------------------------------------ data
  configKey() {
    const c = this.cfg;
    const base = [c.gene, c.output, c.coords];
    if (c.mode === 'individuals') return JSON.stringify([...base, 'i', c.tracks, this.seriesSpec()]);
    if (c.mode === 'groups') return JSON.stringify([...base, 'g', c.tracks, c.hap, c.groupField, this.groupValues(), state.filters]);
    return JSON.stringify([...base, 'p', c.popTrack, c.hap, c.groupField, state.filters, c.popRows]);
  }

  individuals() { return state.pinned.slice(0, MAX_INDIVIDUALS); }

  seriesSpec() {
    const haps = this.cfg.hap === 'both' ? ['H1', 'H2'] : [this.cfg.hap];
    return this.individuals().flatMap((s) => haps.map((hp) => `${s}:${hp}`)).join(',');
  }

  cohortCounts() {
    const field = this.cfg.groupField;
    const col = state.samples.col[field];
    const counts = new Map();
    if (col === undefined) return counts;
    for (const row of filterRows(state.filters, '')) {
      const v = row[col];
      if (v !== null && v !== undefined && v !== '') counts.set(String(v), (counts.get(String(v)) || 0) + 1);
    }
    return counts;
  }

  groupValues() {
    const counts = this.cohortCounts();
    const all = [...counts.keys()].sort((a, b) => counts.get(b) - counts.get(a) || a.localeCompare(b));
    const chosen = (this.cfg.groups || []).filter((g) => counts.has(g));
    return (chosen.length ? chosen : all).slice(0, 8);
  }

  invalidate() {
    this.clearHover();
    this.persist();
    this.request.cancel?.();
    this.requestNow();
  }

  persist() {
    const c = this.cfg;
    setLocus({ gene: c.gene, output: c.output, coords: c.coords, mode: c.mode, hap: c.hap, tracks: c.tracks, groupField: c.groupField, groups: c.groups, yScale: c.yScale, sharedY: c.sharedY, diff: c.diff, envelope: c.envelope, band: c.band, popTrack: c.popTrack, popRows: c.popRows, start: this.viewport.start, end: this.viewport.end });
    updateRouteParams({ gene: c.gene, coords: c.coords !== 'reference' ? c.coords : null, mode: c.mode !== 'individuals' ? c.mode : null, start: this.viewport.start, end: this.viewport.end });
  }

  onViewChange() {
    this.renderLocus();
    this.render();
    this.persistDebounced();
    this.request();
  }

  persistDebounced = debounce(() => this.persist(), 300);

  request = debounce(() => this.requestNow(), 120);

  needsFetch(data, key, margin) {
    if (!data || data.key !== key) return true;
    const v = this.viewport;
    if (v.start < data.start || v.end > data.end) return true;
    const bins = Math.max(1, data.edges.length - 1);
    const dataBp = (data.end - data.start) / bins;
    const wantBp = v.span / this.geom().width;
    if (dataBp > Math.max(1, wantBp) * 1.6) return true; // too coarse
    if (margin && dataBp < wantBp / 4 && bins > 3000) return true; // far too fine
    return false;
  }

  async requestNow() {
    if (!this.info || !this.cfg.output) { this.render(); return; }
    const key = this.configKey();
    const mode = this.cfg.mode;
    const current = mode === 'population' ? this.pop : this.data;
    if (!this.needsFetch(current, key, mode !== 'population')) { this.render(); return; }
    const v = this.viewport;
    const g = this.geom();
    const marginFactor = mode === 'population' ? 0.5 : 1.0;
    const margin = v.span * marginFactor;
    const start = Math.max(0, Math.floor(v.start - margin));
    const end = Math.min(v.length, Math.ceil(v.end + margin));
    const bins = Math.max(32, Math.min(mode === 'population' ? 2048 : 8192, Math.round(g.width * (end - start) / v.span)));
    const signal = this.latest.next();
    const common = { gene: this.cfg.gene, output: this.cfg.output, coords: this.cfg.coords, start, end, bins };
    this.loadingLine.hidden = false;
    this.panelsHost.classList.add('stale');
    let overlay = null;
    const onProgress = (job) => {
      if (!overlay) {
        overlay = jobOverlay(this.scroll, job.title || 'Computing', (jobId) => {
          api(`/api/jobs/${jobId}/cancel`, { method: 'POST', body: {} }).catch(() => {});
          this.latest.abort();
          overlay.remove();
          this.loadingLine.hidden = true;
          this.panelsHost.classList.remove('stale');
          toast('Job cancelled');
        });
        watchJobs();
      }
      overlay.update(job);
    };
    try {
      if (mode === 'individuals') {
        if (!this.individuals().length) { this.data = null; this.render(); return; }
        const res = await api(`${ds()}/signal`, { params: { ...common, tracks: this.cfg.tracks.join(','), series: this.seriesSpec() }, signal });
        this.data = this.toItems(res, key, 'series');
        const failed = res.series.filter((s) => s.error);
        if (failed.length) toast(`${failed.length} series failed: ${failed[0].error}`, 'error', 7000);
      } else if (mode === 'groups') {
        const res = await apiJob(`${ds()}/groups`, { params: { ...common, tracks: this.cfg.tracks.join(','), haps: this.cfg.hap === 'both' ? 'H1,H2' : this.cfg.hap, field: this.cfg.groupField, groups: this.groupValues().join(','), filters: state.filters }, signal, onProgress });
        this.data = this.toItems(res, key, 'groups');
      } else {
        const hap = this.cfg.hap === 'both' ? 'H1+H2' : this.cfg.hap;
        const res = await apiJob(`${ds()}/population`, { params: { ...common, track: this.cfg.popTrack, hap, field: this.cfg.groupField, filters: state.filters, max_rows: this.cfg.popRows }, signal, onProgress });
        res.key = key;
        this.pop = res;
      }
      this.status('');
    } catch (err) {
      if (isAbort(err)) return;
      console.error(err);
      toast(err.message, 'error', 8000);
      this.status('Request failed');
    } finally {
      if (overlay) overlay.remove();
      if (!signal.aborted) { this.loadingLine.hidden = true; this.panelsHost.classList.remove('stale'); }
    }
    this.render();
  }

  toItems(res, key, kind) {
    const items = [];
    if (kind === 'series') {
      const keys = [...new Set(res.series.map((s) => s.sample))];
      this.colors.sync(keys);
      for (const s of res.series) {
        if (s.error) continue;
        items.push({ key: `${s.sample}:${s.haplotype}`, entity: s.sample, label: this.cfg.hap === 'both' ? `${s.sample} ${s.haplotype}` : s.sample, sub: this.cfg.hap === 'both' ? '' : (s.haplotype === 'H1+H2' ? 'diploid mean' : s.haplotype), dash: s.haplotype === 'H2' && this.cfg.hap === 'both', mean: s.mean, min: s.min, max: s.max });
      }
    } else {
      this.colors.sync(res.groups.map((g) => g.group));
      for (const g of res.groups) items.push({ key: g.group, entity: g.group, label: g.group, sub: `n=${fmtInt(g.samples)}`, haplotypes: g.haplotypes, mean: g.mean, min: g.min, max: g.max, upper: g.upper, lower: g.lower, failed: g.failed });
    }
    return { key, kind, start: res.start, end: res.end, edges: res.edges, tracks: res.tracks, items };
  }

  // ------------------------------------------------------------------ rendering
  render() {
    if (!this.info) return;
    this.renderLocus();
    this.drawOverview();
    this.drawRuler();
    if (this.cfg.mode === 'population') this.drawPopulation();
    else this.drawPanels();
    this.renderLegend();
    this.renderFooter();
  }

  xOf(pos) {
    const g = this.geom();
    const v = this.viewport;
    return g.left + ((pos - v.start) / v.span) * g.width;
  }

  // Overview strip: whole window, gene models, model window and current view.
  drawOverview() {
    const g = this.geom();
    const canvas = this.overviewHost.querySelector('canvas');
    const H = 38;
    const ctx = setupCanvas(canvas, g.total, H);
    const t = theme();
    const L = this.viewport.length;
    const x = (p) => g.left + (p / L) * g.width;
    ctx.fillStyle = t.surface2;
    ctx.fillRect(g.left, 8, g.width, H - 16);
    const mw = this.modelWindow();
    if (mw) {
      ctx.fillStyle = t.modelWindow; ctx.fillRect(x(mw.start), 8, Math.max(1, x(mw.end) - x(mw.start)), H - 16);
    }
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, 0, g.width, H); ctx.clip();
    for (const gene of this.annotationGenes()) {
      const a = x(this.fromGenomicOffset(gene.start)); const b = x(this.fromGenomicOffset(gene.end));
      ctx.fillStyle = gene.name === this.cfg.gene ? t.ink2 : t.ink3;
      ctx.fillRect(a, H / 2 - 2, Math.max(1, b - a), 4);
    }
    ctx.restore();
    const v = this.viewport;
    const vx = x(v.start); const vw = Math.max(2, x(v.end) - vx);
    ctx.fillStyle = withAlpha(t.accent, 0.12); ctx.fillRect(vx, 4, vw, H - 8);
    ctx.strokeStyle = t.accent; ctx.lineWidth = 1.5; ctx.strokeRect(vx + 0.75, 4.75, vw - 1.5, H - 9.5);
    ctx.font = `10.5px ${t.font}`; ctx.fillStyle = t.ink3; ctx.textBaseline = 'middle';
    ctx.textAlign = 'right'; ctx.fillText('window', g.left - 8, H / 2);
  }

  attachOverview(el) {
    let dragging = false;
    const move = (e) => {
      const g = this.geom();
      const rect = el.getBoundingClientRect();
      const frac = (e.clientX - rect.left - g.left) / g.width;
      this.viewport.center(frac * this.viewport.length);
    };
    el.style.cursor = 'pointer';
    el.addEventListener('pointerdown', (e) => { dragging = true; el.setPointerCapture(e.pointerId); move(e); });
    el.addEventListener('pointermove', (e) => { if (dragging) move(e); });
    el.addEventListener('pointerup', () => { dragging = false; });
    el.title = 'Click or drag to move the view';
  }

  annotationGenes() {
    return (this.annotations && this.annotations.genes) || [];
  }

  drawRuler() {
    const g = this.geom();
    const canvas = this.rulerHost.querySelector('canvas');
    const genes = this.annotationGenes();
    const rows = this.packGenes(genes);
    const laneH = rows.length ? Math.min(rows.length, 4) * 24 + 6 : 0;
    const H = 30 + laneH;
    const ctx = setupCanvas(canvas, g.total, H);
    const t = theme();
    const v = this.viewport;
    // model window band
    const mw = this.modelWindow();
    if (mw && mw.end > v.start && mw.start < v.end) {
      const a = Math.max(g.left, this.xOf(mw.start)); const b = Math.min(g.left + g.width, this.xOf(mw.end));
      ctx.fillStyle = t.modelWindow; ctx.fillRect(a, 0, b - a, H);
      ctx.fillStyle = t.ink3; ctx.font = `10.5px ${t.font}`; ctx.textBaseline = 'top'; ctx.textAlign = 'left';
      if (b - a > 90) ctx.fillText('CNN model window', a + 4, 2);
    }
    // ticks (genomic positions in reference / aligned coords)
    ctx.font = `11px ${t.font}`; ctx.textBaseline = 'alphabetic'; ctx.textAlign = 'center';
    const offset = this.cfg.coords === 'reference' && this.info.start ? this.info.start : (this.cfg.coords === 'aligned' ? 0 : 1);
    const tickVals = ticks(v.start + offset, v.end + offset, Math.max(3, Math.floor(g.width / 120)));
    ctx.strokeStyle = t.lineStrong; ctx.lineWidth = 1;
    for (const tv of tickVals) {
      const pos = tv - offset;
      const x = Math.round(this.xOf(pos)) + 0.5;
      if (x < g.left || x > g.left + g.width) continue;
      ctx.beginPath(); ctx.moveTo(x, 22); ctx.lineTo(x, 29); ctx.stroke();
      ctx.fillStyle = t.ink2;
      ctx.fillText(this.cfg.coords === 'aligned' ? this.formatPos(pos) : fmtInt(tv), x, 18);
    }
    ctx.beginPath(); ctx.moveTo(g.left, 29.5); ctx.lineTo(g.left + g.width, 29.5); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillStyle = t.ink3; ctx.font = `10.5px ${t.font}`;
    ctx.fillText(this.cfg.coords === 'reference' ? (this.info.chromosome || 'pos') : this.cfg.coords === 'aligned' ? 'axis' : 'hap pos', g.left - 8, 18);
    if (!laneH) return;
    ctx.fillText('genes', g.left - 8, 30 + 16);
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, 30, g.width, laneH); ctx.clip();
    rows.slice(0, 4).forEach((row, r) => {
      const y = 30 + 6 + r * 24;
      for (const gene of row) this.drawGeneModel(ctx, t, gene, y);
    });
    ctx.restore();
    this.geneRows = rows;
  }

  packGenes(genes) {
    const v = this.viewport;
    const visible = genes.filter((gn) => this.fromGenomicOffset(gn.end) > v.start && this.fromGenomicOffset(gn.start) < v.end);
    const rows = [];
    const pxPerBp = this.geom().width / v.span;
    for (const gene of visible) {
      const a = this.fromGenomicOffset(gene.start);
      const labelBp = (gene.name.length * 7 + 12) / pxPerBp;
      const b = this.fromGenomicOffset(gene.end) + labelBp;
      let placed = false;
      for (const row of rows) {
        if (row.end < a) { row.items.push(gene); row.end = b; placed = true; break; }
      }
      if (!placed) rows.push({ end: b, items: [gene] });
    }
    return rows.map((r) => r.items);
  }

  drawGeneModel(ctx, t, gene, y) {
    const tx = gene.transcripts && gene.transcripts[0];
    const color = gene.name === this.cfg.gene ? t.accent : t.ink2;
    const x = (p) => this.xOf(this.fromGenomicOffset(p));
    const a = x(gene.start); const b = x(gene.end);
    ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(a, y + 6.5); ctx.lineTo(b, y + 6.5); ctx.stroke();
    // strand chevrons
    const step = 26;
    const dir = gene.strand === '-' ? -1 : 1;
    for (let cx = Math.max(a, this.geom().left) + 10; cx < Math.min(b, this.geom().left + this.geom().width) - 4; cx += step) {
      ctx.beginPath(); ctx.moveTo(cx - 2 * dir, y + 4); ctx.lineTo(cx + 1 * dir, y + 6.5); ctx.lineTo(cx - 2 * dir, y + 9); ctx.stroke();
    }
    if (tx) {
      for (const [s, e] of tx.exons) { const xa = x(s); ctx.fillRect(xa, y + 3, Math.max(1, x(e) - xa), 7); }
      for (const [s, e] of tx.cds) { const xa = x(s); ctx.fillRect(xa, y + 1, Math.max(1, x(e) - xa), 11); }
    } else {
      ctx.fillRect(a, y + 3, Math.max(1, b - a), 7);
    }
    ctx.font = `600 11px ${t.font}`; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
    const right = this.geom().left + this.geom().width;
    const label = `${gene.name} ${gene.strand === '-' ? '←' : '→'}`;
    const lw = ctx.measureText(label).width;
    // Beside the gene when there is room, otherwise on a chip at the left edge of the view.
    let lx = b + 5;
    if (lx + lw > right) lx = Math.max(this.geom().left + 4, a + 4);
    if (lx + lw <= right) {
      ctx.fillStyle = t.surface; ctx.fillRect(lx - 3, y, lw + 6, 13);
      ctx.fillStyle = color; ctx.fillText(label, lx, y + 6.5);
    }
  }

  ensurePanels(n, height) {
    const host = this.panelsHost;
    if (host.dataset.kind !== 'panels' || host.children.length !== n) {
      clear(host);
      host.dataset.kind = 'panels';
      for (let i = 0; i < n; i++) host.appendChild(h('div', { class: 'panel' }, h('div', { class: 'lane-label', style: { top: '6px' } }), h('canvas')));
    }
    return [...host.children].map((p) => ({ el: p, label: p.querySelector('.lane-label'), canvas: p.querySelector('canvas'), height }));
  }

  transform(v) {
    if (this.cfg.yScale === 'log') return Math.sign(v) * Math.log1p(Math.abs(v));
    return v;
  }

  visibleRange(data) {
    const v = this.viewport;
    const e = data.edges;
    let i0 = 0; let i1 = e.length - 1;
    while (i0 < e.length - 1 && e[i0 + 1] <= v.start) i0++;
    while (i1 > i0 && e[i1 - 1] >= v.end) i1--;
    return [i0, i1];
  }

  drawPanels() {
    const g = this.geom();
    const data = this.data;
    const meta = this.trackMeta();
    const tracks = this.cfg.tracks || [];
    const panels = this.ensurePanels(Math.max(1, tracks.length), PANEL_H);
    const t = theme();
    if (!data || !data.items.length) {
      panels.forEach((p, i) => {
        const ctx = setupCanvas(p.canvas, g.total, PANEL_H);
        p.label.textContent = tracks.length ? (meta[tracks[i]] || {}).label || `Track ${tracks[i]}` : 'No tracks selected';
        ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.textAlign = 'center';
        const msg = this.cfg.mode === 'individuals' && !this.individuals().length ? 'Pin individuals on the Samples page (or add them in the panel →)' : (data ? 'No data' : 'Loading…');
        ctx.fillText(msg, g.left + g.width / 2, PANEL_H / 2);
      });
      return;
    }
    const [i0, i1] = this.visibleRange(data);
    const items = this.displayItems(data).filter((it) => !this.hidden.has(it.key));
    const ranges = data.tracks.map((_, ti) => {
      let lo = 0; let hi = 0;
      const band = this.showBand();
      for (const it of items) {
        const top = this.cfg.envelope && it.max ? it.max[ti] : it.mean[ti];
        const upper = band && it.upper ? it.upper[ti] : null;
        const bottom = band && it.lower ? it.lower[ti] : (this.cfg.envelope && it.min ? it.min[ti] : it.mean[ti]);
        for (let i = i0; i < i1; i++) {
          const a = top[i]; if (a > hi) hi = a;
          if (upper) { const u = upper[i]; if (u > hi) hi = u; }
          const b = bottom[i]; if (b < lo) lo = b;
        }
      }
      return [this.transform(lo), this.transform(hi)];
    });
    let shared = null;
    if (this.cfg.sharedY) shared = [Math.min(...ranges.map((r) => r[0])), Math.max(...ranges.map((r) => r[1]))];
    data.tracks.forEach((track, ti) => {
      const p = panels[ti];
      if (!p) return;
      const m = meta[track] || {};
      p.label.textContent = m.label || `Track ${track}`;
      p.label.title = Object.entries(m.metadata || {}).map(([k, val]) => `${k}: ${val}`).join('\n');
      const [ylo, yhiRaw] = shared || ranges[ti];
      const yhi = yhiRaw > ylo ? yhiRaw * 1.06 : ylo + 1;
      this.drawPanel(p.canvas, t, g, data, items, ti, ylo, yhi, i0, i1);
      p.yRange = [ylo, yhi];
    });
    this.panelGeom = { top: 0, height: PANEL_H, count: data.tracks.length };
  }

  /** The ±SD band is within-group spread; it is hidden when showing differences between groups. */
  showBand() { return this.cfg.mode === 'groups' && this.cfg.band && !this.cfg.diff; }

  /** Group means as-is, or (when "difference" is on) minus the haplotype-weighted cohort mean. */
  displayItems(data) {
    if (!(this.cfg.mode === 'groups' && this.cfg.diff && data.kind === 'groups' && data.items.length > 1)) return data.items;
    if (data.diffItems) return data.diffItems;
    const weights = data.items.map((it) => it.haplotypes || 1);
    const baseline = data.tracks.map((_, ti) => {
      const n = data.items[0].mean[ti].length;
      const out = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        let acc = 0; let w = 0;
        data.items.forEach((it, k) => { const v = it.mean[ti][i]; if (Number.isFinite(v)) { acc += v * weights[k]; w += weights[k]; } });
        out[i] = w ? acc / w : NaN;
      }
      return out;
    });
    const sub = (arrs) => arrs && arrs.map((a, ti) => a.map((v, i) => v - baseline[ti][i]));
    data.diffItems = data.items.map((it) => ({ ...it, mean: sub(it.mean), min: sub(it.min), max: sub(it.max), upper: sub(it.upper), lower: sub(it.lower) }));
    return data.diffItems;
  }

  drawPanel(canvas, t, g, data, items, ti, ylo, yhi, i0, i1) {
    const H = PANEL_H;
    const ctx = setupCanvas(canvas, g.total, H);
    const top = 24; const bottom = H - 8;
    const sy = (val) => bottom - ((this.transform(val) - ylo) / (yhi - ylo)) * (bottom - top);
    const v = this.viewport;
    // model window band
    const mw = this.modelWindow();
    if (mw && mw.end > v.start && mw.start < v.end) {
      const a = Math.max(g.left, this.xOf(mw.start)); const b = Math.min(g.left + g.width, this.xOf(mw.end));
      ctx.fillStyle = t.modelWindow; ctx.fillRect(a, top - 4, b - a, bottom - top + 4);
    }
    // grid + y labels
    ctx.font = `10.5px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.lineWidth = 1;
    const yt = ticks(ylo, yhi, 3);
    for (const val of yt) {
      const y = Math.round(bottom - ((val - ylo) / (yhi - ylo)) * (bottom - top)) + 0.5;
      if (y < top - 1 || y > bottom + 1) continue;
      ctx.strokeStyle = t.line; ctx.beginPath(); ctx.moveTo(g.left, y); ctx.lineTo(g.left + g.width, y); ctx.stroke();
      ctx.fillStyle = t.ink3; ctx.fillText(this.cfg.yScale === 'log' ? fmtNum(Math.sign(val) * Math.expm1(Math.abs(val)), 2) : fmtNum(val, 3), g.left - 6, y);
    }
    ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(g.left, bottom + 0.5); ctx.lineTo(g.left + g.width, bottom + 0.5); ctx.stroke();
    if (ylo < 0 && yhi > 0) {
      const y0 = Math.round(sy(0)) + 0.5;
      ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(g.left, y0); ctx.lineTo(g.left + g.width, y0); ctx.stroke();
    }
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, 0, g.width, H); ctx.clip();
    const e = data.edges;
    const binBp = (data.end - data.start) / Math.max(1, e.length - 1);
    const pxPerBin = (binBp / v.span) * g.width;
    const step = pxPerBin > 3;
    const xs = new Float64Array(e.length);
    for (let i = 0; i < e.length; i++) xs[i] = this.xOf(e[i]);
    const from = Math.max(0, i0 - 1); const to = Math.min(e.length - 1, i1 + 1);
    const area = (lo, hi, color) => {
      ctx.fillStyle = color;
      let open = false;
      const pts = [];
      const flush = () => {
        if (pts.length < 2) { pts.length = 0; return; }
        ctx.beginPath();
        ctx.moveTo(pts[0][0], pts[0][1]);
        for (const p of pts) ctx.lineTo(p[0], p[1]);
        for (let k = pts.length - 1; k >= 0; k--) ctx.lineTo(pts[k][0], pts[k][2]);
        ctx.closePath(); ctx.fill(); pts.length = 0;
      };
      for (let i = from; i < to; i++) {
        const a = lo[i]; const b = hi[i];
        if (!Number.isFinite(a) || !Number.isFinite(b)) { flush(); open = false; continue; }
        if (step) { pts.push([xs[i], sy(b), sy(a)]); pts.push([xs[i + 1], sy(b), sy(a)]); }
        else pts.push([(xs[i] + xs[i + 1]) / 2, sy(b), sy(a)]);
        open = true;
      }
      if (open) flush();
    };
    const line = (vals, color, dash, width) => {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash ? [5, 4] : []); ctx.lineJoin = 'round';
      ctx.beginPath();
      let pen = false;
      for (let i = from; i < to; i++) {
        const val = vals[i];
        if (!Number.isFinite(val)) { pen = false; continue; }
        const y = sy(val);
        if (step) {
          if (pen) ctx.lineTo(xs[i], y); else ctx.moveTo(xs[i], y);
          ctx.lineTo(xs[i + 1], y);
        } else {
          const x = (xs[i] + xs[i + 1]) / 2;
          if (pen) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        }
        pen = true;
      }
      ctx.stroke();
      ctx.setLineDash([]);
    };
    const highlight = this.highlight;
    for (const it of items) {
      const color = this.colors.color(it.entity);
      const dim = highlight && highlight !== it.key;
      if (this.showBand() && it.upper) area(it.lower[ti], it.upper[ti], withAlpha(color, dim ? 0.04 : 0.14));
      else if (this.cfg.envelope && binBp > 1.5 && it.min) area(it.min[ti], it.max[ti], withAlpha(color, dim ? 0.04 : 0.16));
    }
    for (const it of items) {
      const color = this.colors.color(it.entity);
      const dim = highlight && highlight !== it.key;
      line(it.mean[ti], dim ? withAlpha(color, 0.25) : color, it.dash, highlight === it.key ? 2.5 : 1.75);
    }
    ctx.restore();
  }

  drawPopulation() {
    const g = this.geom();
    const host = this.panelsHost;
    if (host.dataset.kind !== 'population') {
      clear(host);
      host.dataset.kind = 'population';
      host.appendChild(h('div', { class: 'panel' }, h('div', { class: 'lane-label', style: { top: '6px' } }), h('canvas')));
    }
    const panel = host.firstChild;
    const label = panel.querySelector('.lane-label');
    const canvas = panel.querySelector('canvas');
    const t = theme();
    const pop = this.pop;
    const meta = this.trackMeta()[this.cfg.popTrack] || {};
    label.textContent = `${meta.label || `Track ${this.cfg.popTrack}`} · ${this.cfg.hap === 'both' ? 'H1+H2 mean' : this.cfg.hap}`;
    if (!pop) {
      const ctx = setupCanvas(canvas, g.total, 200);
      ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.textAlign = 'center'; ctx.fillText('Loading…', g.left + g.width / 2, 100);
      return;
    }
    const rows = pop.samples.length;
    const rowH = Math.max(1, Math.min(10, Math.floor(620 / rows)));
    const top = 26;
    const H = top + rows * rowH + 34;
    const ctx = setupCanvas(canvas, g.total, H);
    const cols = pop.matrix[0] ? pop.matrix[0].length : 0;
    // color scale: 99th percentile of finite values
    const vals = [];
    const stride = Math.max(1, Math.floor((rows * cols) / 40000));
    const flat = pop.matrix.flat || null;
    if (flat) for (let i = 0; i < flat.length; i += stride) { const x = flat[i]; if (Number.isFinite(x)) vals.push(this.transform(x)); }
    vals.sort((a, b) => a - b);
    const vmax = vals.length ? Math.max(vals[Math.floor(vals.length * 0.99)] || 0, 1e-9) : 1;
    // render matrix into an offscreen bitmap then scale into the view
    const off = document.createElement('canvas');
    off.width = Math.max(1, cols); off.height = Math.max(1, rows);
    const octx = off.getContext('2d');
    const img = octx.createImageData(off.width, off.height);
    const missing = hexOrRgb(t.surface2);
    for (let r = 0; r < rows; r++) {
      const row = pop.matrix[r];
      for (let c = 0; c < cols; c++) {
        const k = (r * cols + c) * 4;
        const val = row[c];
        const rgb = Number.isFinite(val) ? seqColor(t.seq, this.transform(val) / vmax) : missing;
        img.data[k] = rgb[0]; img.data[k + 1] = rgb[1]; img.data[k + 2] = rgb[2]; img.data[k + 3] = 255;
      }
    }
    octx.putImageData(img, 0, 0);
    const e = pop.edges;
    const x0 = this.xOf(e[0]); const x1 = this.xOf(e[e.length - 1]);
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, top, g.width, rows * rowH); ctx.clip();
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, x0, top, x1 - x0, rows * rowH);
    // model window edges
    const mw = this.modelWindow();
    if (mw) {
      ctx.strokeStyle = t.modelWindowEdge; ctx.lineWidth = 1;
      for (const p of [mw.start, mw.end]) { const x = Math.round(this.xOf(p)) + 0.5; ctx.beginPath(); ctx.moveTo(x, top); ctx.lineTo(x, top + rows * rowH); ctx.stroke(); }
    }
    ctx.restore();
    // group bands
    ctx.font = `11px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    let start = 0;
    const labels = pop.labels || [];
    for (let r = 1; r <= rows; r++) {
      if (r === rows || labels[r] !== labels[start]) {
        const y0 = top + start * rowH; const y1 = top + r * rowH;
        ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(g.left - 4, y1 + 0.5); ctx.lineTo(g.left + g.width, y1 + 0.5); ctx.stroke();
        if (y1 - y0 >= 10) { ctx.fillStyle = t.ink2; ctx.fillText(truncate(labels[start] || '–', 9), g.left - 6, (y0 + y1) / 2); }
        start = r;
      }
    }
    // color legend
    const ly = top + rows * rowH + 12;
    const lw = Math.min(220, g.width / 3);
    for (let i = 0; i < lw; i++) { const c = seqColor(t.seq, i / lw); ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`; ctx.fillRect(g.left + i, ly, 1, 9); }
    ctx.fillStyle = t.ink3; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.font = `11px ${t.font}`;
    const vmaxLabel = this.cfg.yScale === 'log' ? Math.expm1(vmax) : vmax;
    ctx.fillText(`0 → ${fmtNum(vmaxLabel)} (99th pct${this.cfg.yScale === 'log' ? ', log1p' : ''}) · ${fmtInt(rows)} samples${pop.failed_count ? ` · ${pop.failed_count} missing` : ''} · click a row to pin`, g.left + lw + 10, ly + 5);
    this.popGeom = { top, rowH, rows };
  }

  renderLegend() {
    clear(this.legend);
    const data = this.cfg.mode === 'population' ? null : this.data;
    if (this.cfg.mode === 'population') {
      this.legend.append(h('span', null, `Rows sorted by ${this.cfg.groupField || 'sample'}; cohort from the Samples page filters.`));
      return;
    }
    if (!data) { this.legend.append(h('span', { class: 'muted' }, this.cfg.mode === 'groups' ? 'Group means over the cohort' : 'Pinned individuals')); return; }
    for (const it of data.items) {
      const color = this.colors.color(it.entity);
      const sw = it.dash ? h('span', { class: 'swatch line dashed', style: { color, borderTopColor: color } }) : h('span', { class: 'swatch line', style: { background: color } });
      const el = h('span', { class: `legend-item${this.hidden.has(it.key) ? ' off' : ''}`, title: 'Click to hide/show; hover to highlight' }, sw, h('span', null, it.label), h('span', { class: 'muted' }, it.sub || ''));
      el.addEventListener('click', () => { if (this.hidden.has(it.key)) this.hidden.delete(it.key); else this.hidden.add(it.key); this.render(); });
      el.addEventListener('mouseenter', () => { this.highlight = it.key; this.drawPanels(); });
      el.addEventListener('mouseleave', () => { this.highlight = null; this.drawPanels(); });
      this.legend.appendChild(el);
    }
    if (this.cfg.mode === 'groups' && this.cfg.diff) this.legend.append(h('span', { class: 'muted' }, 'Each line: group mean minus the haplotype-weighted mean of the shown groups'));
    else if (this.showBand()) this.legend.append(h('span', { class: 'muted' }, 'Shaded: ±1 SD across haplotypes'));
    else if (this.cfg.envelope) this.legend.append(h('span', { class: 'muted' }, 'Shaded: min–max within each pixel bin'));
  }

  renderFooter() {
    const c = this.cfg;
    const coordText = { reference: 'genomic coordinates (haplotypes remapped through their indels; deletions show as gaps)', haplotype: 'raw haplotype prediction index', aligned: 'bcftools_chain training axis (insertion columns shared across the cohort)' }[c.coords];
    const res = this.data && this.cfg.mode !== 'population' ? ` · resolution ${fmtBp((this.data.end - this.data.start) / Math.max(1, this.data.edges.length - 1))}/bin` : '';
    this.footer.textContent = `${c.gene} · ${c.output || 'no output'} · ${coordText}${res}`;
  }

  drawBrush(a, b, done) {
    if (done) { this.brush.hidden = true; return; }
    const g = this.geom();
    const lo = Math.max(0, Math.min(a, b)); const hi = Math.min(1, Math.max(a, b));
    this.brush.hidden = false;
    this.brush.style.left = `${g.left + lo * g.width}px`;
    this.brush.style.width = `${(hi - lo) * g.width}px`;
  }

  // ------------------------------------------------------------------ hover
  clearHover() { this.crosshair.hidden = true; hideTooltip(); }

  onHover(e, frac) {
    if (frac < 0 || frac > 1) { this.clearHover(); return; }
    const g = this.geom();
    const pos = this.viewport.start + frac * this.viewport.span;
    this.crosshair.hidden = false;
    this.crosshair.style.left = `${g.left + frac * g.width}px`;
    const target = e.target.closest ? e.target.closest('.panel') : null;
    const coordLabel = this.formatPos(Math.floor(pos), true);
    if (this.cfg.mode === 'population') {
      if (!this.pop || !this.popGeom || !target) { showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div>`); return; }
      const rect = target.getBoundingClientRect();
      const r = Math.floor((e.clientY - rect.top - this.popGeom.top) / this.popGeom.rowH);
      if (r < 0 || r >= this.popGeom.rows) { showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div>`); return; }
      const bi = binIndex(this.pop.edges, pos);
      const val = bi >= 0 ? this.pop.matrix[r][bi] : NaN;
      showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(this.pop.samples[r])}</div><div class="tt-row"><span class="tt-key">${escapeHtml(this.cfg.groupField || '')}</span><span class="tt-val">${escapeHtml(this.pop.labels[r] || '–')}</span></div><div class="tt-row"><span class="tt-key">value</span><span class="tt-val">${fmtNum(val)}</span></div><div class="tt-sub">${escapeHtml(coordLabel)}</div>`);
      this.hoverRow = this.pop.samples[r];
      return;
    }
    if (!this.data || !target) { showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div>`); return; }
    const ti = [...this.panelsHost.children].indexOf(target);
    const bi = binIndex(this.data.edges, pos);
    const meta = this.trackMeta()[this.data.tracks[ti]] || {};
    const rows = [];
    const items = this.displayItems(this.data).filter((it) => !this.hidden.has(it.key));
    const sorted = items.map((it) => ({ it, v: bi >= 0 && it.mean[ti] ? it.mean[ti][bi] : NaN })).sort((a, b) => (Number.isFinite(b.v) ? b.v : -Infinity) - (Number.isFinite(a.v) ? a.v : -Infinity));
    for (const { it, v } of sorted.slice(0, 12)) {
      const color = this.colors.color(it.entity);
      const extra = it.upper && bi >= 0 ? ` <span class="muted">±${fmtNum((it.upper[ti][bi] - it.lower[ti][bi]) / 2, 2)}</span>` : '';
      rows.push(`<div class="tt-row"><span class="tt-key"><span class="swatch line" style="background:${color}"></span>${escapeHtml(it.label)}</span><span class="tt-val">${Number.isFinite(v) ? fmtNum(v) : '<span class="muted">gap</span>'}${extra}</span></div>`);
    }
    const binBp = bi >= 0 ? this.data.edges[bi + 1] - this.data.edges[bi] : 0;
    showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div><div class="tt-sub" style="margin:-2px 0 6px">${escapeHtml(meta.short || meta.label || '')}${binBp > 1 ? ` · mean of ${fmtInt(binBp)} bp` : ''}</div>${rows.join('')}`);
  }

  // ------------------------------------------------------------------ side panel
  renderSide() {
    const side = this.side;
    clear(side);
    const meta = this.trackMeta();
    // tracks
    if (this.cfg.mode === 'population') {
      side.appendChild(h('div', { class: 'side-section' },
        h('h3', null, 'Track'),
        select(meta.map((m) => ({ value: m.index, label: m.short || m.label })), this.cfg.popTrack, (v) => { this.cfg.popTrack = Number(v); this.invalidate(); }),
        field('Max rows', select([100, 300, 600, 1200, 2400, 4000].map((n) => ({ value: n, label: `${fmtInt(n)} samples` })), this.cfg.popRows, (v) => { this.cfg.popRows = Number(v); this.invalidate(); }))));
    } else {
      const list = h('div', { class: 'checklist' });
      for (const m of meta) {
        const box = h('input', { type: 'checkbox' });
        box.checked = (this.cfg.tracks || []).includes(m.index);
        box.addEventListener('change', () => {
          const set = new Set(this.cfg.tracks);
          if (box.checked) set.add(m.index); else set.delete(m.index);
          this.cfg.tracks = [...set].sort((a, b) => a - b).slice(0, 16);
          this.invalidate();
        });
        list.appendChild(h('label', { title: m.label }, box, h('span', null, m.short || m.label), h('span', { class: 'meta' }, `#${m.index}`)));
      }
      side.appendChild(h('div', { class: 'side-section' },
        h('h3', null, `Tracks (${(this.cfg.tracks || []).length}/${meta.length})`, h('span', null,
          h('button', { class: 'btn small ghost', onclick: () => { this.cfg.tracks = meta.slice(0, 16).map((m) => m.index); this.renderSide(); this.invalidate(); } }, 'All'),
          h('button', { class: 'btn small ghost', onclick: () => { this.cfg.tracks = meta.length ? [meta[0].index] : []; this.renderSide(); this.invalidate(); } }, 'First'))),
        list));
    }
    // individuals / groups
    if (this.cfg.mode === 'individuals') side.appendChild(this.individualsSection());
    else side.appendChild(this.groupsSection());
    // display
    side.appendChild(h('div', { class: 'side-section' },
      h('h3', null, 'Display'),
      field('Value scale', segmented([{ value: 'linear', label: 'Linear' }, { value: 'log', label: 'log(1+x)' }], this.cfg.yScale, (v) => { this.cfg.yScale = v; this.persist(); this.render(); })),
      this.cfg.mode !== 'population' ? checkbox('Same y-scale for all tracks', this.cfg.sharedY, (v) => { this.cfg.sharedY = v; this.persist(); this.render(); }) : null,
      this.cfg.mode === 'individuals' ? checkbox('Min–max envelope when zoomed out', this.cfg.envelope, (v) => { this.cfg.envelope = v; this.persist(); this.render(); }) : null,
      this.cfg.mode === 'groups' ? checkbox('±1 SD band', this.cfg.band, (v) => { this.cfg.band = v; this.persist(); this.render(); }) : null,
      this.cfg.mode === 'groups' ? checkbox('Show difference from cohort mean', this.cfg.diff, (v) => { this.cfg.diff = v; this.persist(); this.render(); }) : null,
      h('p', { class: 'help' }, 'The blue band marks the CNN model window. Hover a legend entry to highlight it; click to hide it.')));
  }

  individualsSection() {
    const ids = state.pinned;
    this.colors.sync(this.individuals());
    const chips = h('div', { class: 'chips' }, ids.length ? ids.map((id, i) => {
      const row = state.samples.index.get(id);
      const extra = row ? [state.samples.col.superpopulation, state.samples.col.population].filter((c) => c !== undefined).map((c) => row[c]).filter(Boolean).join(' · ') : '';
      return h('span', { class: 'pill removable', title: extra, style: i >= MAX_INDIVIDUALS ? { opacity: 0.5 } : null },
        i < MAX_INDIVIDUALS ? h('span', { class: 'swatch', style: { background: this.colors.color(id) } }) : null, id, h('button', { title: `Remove ${id}`, onclick: () => togglePinned(id) }, '×'));
    }) : h('span', { class: 'muted', style: { fontSize: '12px' } }, 'Nothing pinned yet.'));
    const listId = 'sample-options';
    const input = h('input', { class: 'input', list: listId, placeholder: 'Add sample id…', style: { flex: '1' } });
    const options = h('datalist', { id: listId }, filterRows(state.filters, '').slice(0, 4000).map((r) => h('option', { value: r[state.samples.col.sample_id] })));
    const add = () => {
      const id = input.value.trim();
      if (!state.samples.index.has(id)) { toast(`Unknown sample: ${id}`, 'error'); return; }
      setPinned([...state.pinned, id]);
      input.value = '';
    };
    input.addEventListener('keydown', (e) => { if (e.key === 'Enter') add(); });
    return h('div', { class: 'side-section' },
      h('h3', null, `Individuals (${Math.min(ids.length, MAX_INDIVIDUALS)}/${MAX_INDIVIDUALS})`, h('a', { href: '#/samples', style: { fontSize: '12px', fontWeight: '500' } }, 'Browse…')),
      chips,
      h('div', { style: { display: 'flex', gap: '6px' } }, input, h('button', { class: 'btn', onclick: add }, 'Add'), options),
      ids.length > MAX_INDIVIDUALS ? h('p', { class: 'help' }, `Only the first ${MAX_INDIVIDUALS} pinned samples are drawn (one color each).`) : null);
  }

  groupsSection() {
    const cats = categoricalFields();
    const counts = this.cohortCounts();
    const chosen = new Set(this.groupValues());
    const list = h('div', { class: 'checklist' });
    const values = [...counts.keys()].sort((a, b) => counts.get(b) - counts.get(a) || a.localeCompare(b));
    for (const value of values) {
      const box = h('input', { type: 'checkbox' });
      box.checked = chosen.has(value);
      box.addEventListener('change', () => {
        const next = new Set(this.groupValues());
        if (box.checked) next.add(value); else next.delete(value);
        this.cfg.groups = [...next].slice(0, 8);
        if (!this.cfg.groups.length) this.cfg.groups = null;
        this.renderSide();
        this.invalidate();
      });
      list.appendChild(h('label', null, box, h('span', null, value), h('span', { class: 'meta' }, fmtInt(counts.get(value)))));
    }
    const cohort = filterRows(state.filters, '').length;
    const filterText = Object.entries(state.filters).map(([k, v]) => `${k} ∈ {${v.join(', ')}}`).join('; ');
    return h('div', { class: 'side-section' },
      h('h3', null, this.cfg.mode === 'population' ? 'Rows' : 'Groups'),
      field(this.cfg.mode === 'population' ? 'Sort rows by' : 'Group by', select(cats.map((f) => ({ value: f.name, label: f.label })), this.cfg.groupField, (v) => { this.cfg.groupField = v; this.cfg.groups = null; this.renderSide(); this.invalidate(); })),
      this.cfg.mode === 'groups' ? h('div', null, h('div', { class: 'label', style: { marginBottom: '4px' } }, 'Groups (max 8)'), list) : null,
      h('p', { class: 'help' }, `Cohort: ${fmtInt(cohort)} samples${filterText ? ` (${filterText})` : ' (all)'}. `, h('a', { href: '#/samples' }, 'Edit cohort')),
      this.cfg.mode === 'groups' ? h('p', { class: 'help' }, 'Means are computed once over every haplotype in each group (progress is shown) and cached on disk, so later views are instant.') : null,
      this.cfg.coords === 'aligned' ? h('p', { class: 'help', style: { color: 'var(--ink-2)' } }, 'Training-axis aggregates need a bcftools alignment entry per haplotype; uncached entries are built on first use, which is slow for large cohorts. Genomic coordinates give the same per-base alignment without that cost.') : null);
  }

  // ------------------------------------------------------------------ lifecycle
  onEvent(topic) {
    if (topic === 'pinned') { this.renderSide(); if (this.cfg.mode === 'individuals') this.invalidate(); }
    if (topic === 'cohort' && this.cfg.mode !== 'individuals') { this.renderSide(); this.invalidate(); }
  }

  update(params) {
    if (params.gene && params.gene !== this.cfg.gene) this.setGene(params.gene, {});
  }

  unmount() {
    this.latest.abort();
    this.request.cancel();
    this.persistDebounced.cancel();
    hideTooltip();
    for (const d of this.disposers) d();
  }

  async loadAnnotations() {
    const gene = this.cfg.gene;
    try {
      const res = await apiJob(`${ds()}/genes/${encodeURIComponent(gene)}/annotations`, { onProgress: () => this.status('Loading gene annotations…') });
      if (gene !== this.cfg.gene) return;
      this.annotations = res;
      this.status('');
      this.render();
    } catch (err) {
      if (!isAbort(err)) this.status('');
    }
  }
}

function checkbox(label, value, onChange) {
  const box = h('input', { type: 'checkbox' });
  box.checked = !!value;
  box.addEventListener('change', () => onChange(box.checked));
  return h('label', { class: 'check' }, box, label);
}

function binIndex(edges, pos) {
  let lo = 0; let hi = edges.length - 1;
  if (pos < edges[0] || pos >= edges[hi]) return -1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (edges[m] <= pos) lo = m; else hi = m; }
  return lo;
}

function truncate(text, n) { return text.length > n ? `${text.slice(0, n - 1)}…` : text; }

function hexOrRgb(color) {
  const m = /^#?([0-9a-f]{6})$/i.exec(color || '');
  if (m) { const n = parseInt(m[1], 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
  return [128, 128, 128];
}
