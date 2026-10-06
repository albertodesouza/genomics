// Tracks: interactive AlphaGenome signal browser.
//
// Tracks:  any mix of outputs (RNA-seq, CAGE, DNase, ChIP, …), one panel per track, plus gene models
//          and a sequence lane (per-base letter frequencies of the pinned haplotypes or reference).
// Modes:   individuals (pinned samples) · groups (cohort means ± SD by any facet) · population
//          heatmap (every cohort sample as a row, sorted by a facet).
// Coords:  genomic (each haplotype remapped through its own indels) · haplotype (raw prediction
//          index) · training axis (bcftools_chain expanded alignment used for CNN inputs).
// Compare: AlphaGenome's prediction of the reference genome (a line in every panel) and the
//          observed ENCODE / FANTOM5 signal behind each track (a lane under the panel).
// Click a track name for its card (ontology terms, assay, experiments), a gene model or the gene
// button for the gene card (HGNC, database links, Gene Ontology).
import { api, apiJob, isAbort, Latest } from '../api.js';
import { navigate, updateRouteParams, watchJobs } from '../app.js';
import { state, ds, geneInfo, setLocus, setPinned, togglePinned, categoricalFields, filterRows, cohortDescription } from '../state.js';
import { exportMenu, slug } from '../figure.js';
import { h, clear, icon, iconButton, segmented, select, field, fmtInt, fmtBp, fmtNum, fmtPct, showTooltip, hideTooltip, toast, jobOverlay, debounce, escapeHtml, errorBox } from '../ui.js';
import { Viewport, ColorSlots, setupCanvas, theme, ticks, withAlpha, seqColor, onResize, css } from '../plot.js';
import { loadGeneAnnotations, packGenes, drawGeneLanes, geneLanesHeight, geneAt, geneTooltip } from '../gene_lanes.js';
import { labelGeneOptions, openGeneCard, openTrackCard } from '../cards.js';
import { openPredictForm } from '../forms.js';
import { openScalarForm } from '../scalars.js';

const GUTTER_L = 64;
const GUTTER_R = 14;
const PANEL_H = 118;
const MAX_INDIVIDUALS = 8;
const MAX_TRACKS = 16;
const SEQ_H = 62;
const SEQ_LETTERS_PX = 7; // px per base from which the sequence lane draws letters
const OBS_H = 72; // observed-data lane under each track panel
const REF = '@reference'; // the reference-genome prediction series
const MAX_SEQ_SAMPLES = 24;
const OUTPUT_LABELS = { rna_seq: 'RNA-seq', cage: 'CAGE', procap: 'PRO-cap', dnase: 'DNase', atac: 'ATAC', chip_histone: 'ChIP histone', chip_tf: 'ChIP TF', splice_sites: 'Splice sites', splice_site_usage: 'Splice usage', splice_junctions: 'Splice junctions', contact_maps: 'Contact maps' };
const outputLabel = (name) => OUTPUT_LABELS[String(name).toLowerCase()] || name;
/** Track keys are "output:index" so tracks of several outputs can be shown together. */
const trackKey = (output, index) => `${output}:${index}`;
function parseTrackKey(key) {
  const k = String(key).lastIndexOf(':');
  return k < 0 ? [null, Number(key)] : [key.slice(0, k), Number(key.slice(k + 1))];
}
const isVariantField = (f) => typeof f === 'string' && /^variant:\d+:[A-Z]+:[A-Z]+$/i.test(f);
const groupFieldLabel = (f) => { if (!isVariantField(f)) return f; const [, pos, ref, alt] = f.split(':'); return `genotype at ${Number(pos).toLocaleString('en-US')} ${ref}>${alt}`; };
const defaultGroupField = (cats) => (cats.find((f) => f.name === 'superpopulation') || cats[0] || {}).name;
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
    this.seq = null;
    this.seqLatest = new Latest();
    this.hover = null;
    this.disposers = [];
  }

  // ------------------------------------------------------------------ setup
  async init(params) {
    const genes = state.summary.genes.map((g) => g.gene);
    const saved = state.locus || {};
    const gene = (params.gene && genes.includes(params.gene)) ? params.gene : (saved.gene && genes.includes(saved.gene) ? saved.gene : (state.summary.genes.find((g) => g.in_metadata) || state.summary.genes[0] || {}).gene);
    const cats = categoricalFields();
    // Older saved state kept one output plus numeric track indices.
    const legacyOutput = saved.output || null;
    const migrate = (t) => (typeof t === 'number' ? (legacyOutput ? trackKey(legacyOutput, t) : null) : t);
    this.preferOutput = params.output || legacyOutput;
    this.cfg = {
      gene,
      // A URL locus without coords is genomic (coords are omitted from URLs when genomic).
      coords: params.coords || (params.start !== undefined ? 'reference' : saved.coords) || 'reference',
      mode: params.mode || saved.mode || 'individuals',
      hap: saved.hap || 'H1+H2',
      tracks: Array.isArray(saved.tracks) ? saved.tracks.map(migrate).filter(Boolean) : null,
      // "variant:<pos>:<ref>:<alt>" groups the cohort by genotype at a site (from the Variant page).
      groupField: isVariantField(params.groupField) || cats.some((f) => f.name === params.groupField) ? params.groupField
        : saved.groupField && (cats.some((f) => f.name === saved.groupField) || (isVariantField(saved.groupField) && saved.gene === gene)) ? saved.groupField
          : defaultGroupField(cats),
      groups: saved.groups || null,
      yScale: saved.yScale || 'linear',
      sharedY: saved.sharedY ?? false,
      diff: saved.diff ?? false,
      envelope: saved.envelope ?? true,
      band: saved.band ?? true,
      popTrack: migrate(saved.popTrack ?? null),
      seqSource: saved.seqSource || 'pinned',
      popRows: saved.popRows || 600,
      showRef: saved.showRef ?? true,
      showObserved: saved.showObserved ?? false,
      obsScale: saved.obsScale || 'tpm',
    };
    this.buildLayout();
    this.viewport = new Viewport({ length: 1, start: 0, end: 1, minSpan: 30, onChange: () => this.onViewChange() });
    for (const el of [this.overviewHost, this.rulerHost, this.seqHost, this.panelsHost]) {
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
    // Clicking a gene model on the ruler opens its card.
    let rulerDown = null;
    this.rulerHost.addEventListener('pointerdown', (e) => { rulerDown = [e.clientX, e.clientY]; });
    this.rulerHost.addEventListener('click', (e) => {
      if (!rulerDown || Math.hypot(e.clientX - rulerDown[0], e.clientY - rulerDown[1]) > 4 || !this.geneRows) return;
      const rect = this.rulerHost.getBoundingClientRect();
      const gene = geneAt(this.geneRows, e.clientX - rect.left, e.clientY - rect.top, { top: 30, xOf: this.geneX });
      if (gene) openGeneCard(gene.name);
    });
    this.disposers.push(onResize(this.scroll, () => this.render()));
    const onTheme = () => this.render();
    window.addEventListener('themechange', onTheme);
    this.disposers.push(() => window.removeEventListener('themechange', onTheme));
    await this.setGene(gene, { start: params.start, end: params.end, keepView: !params.gene || params.gene === saved.gene });
  }

  buildLayout() {
    const genes = state.summary.genes.map((g) => ({ value: g.gene, label: g.gene }));
    this.geneSelect = select(genes, this.cfg.gene, (v) => this.setGene(v, {}), { 'aria-label': 'Gene', style: { minWidth: '120px', maxWidth: '300px' } });
    labelGeneOptions(this.geneSelect); // HGNC approved names, when the HGNC table is reachable
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
      field('Gene', h('div', { style: { display: 'flex', gap: '4px', alignItems: 'center' } }, this.geneSelect,
        iconButton('info', 'Gene card: HGNC names, database links (Ensembl, NCBI, UniProt, GTEx, …) and Gene Ontology', () => openGeneCard(this.cfg.gene), 'icon-btn bordered'))),
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
      h('button', { class: 'btn small ghost', title: 'Define one number per sample from the first track over the visible range (a sample field for filters and group means)', onclick: () => this.openScalar() }, 'Region scalar…'),
      exportMenu(() => this.figureOptions()),
      this.statusEl,
      h('span', { class: 'hint' }, 'Drag to pan · ⌘/Ctrl+scroll zoom · Shift+drag region'),
      iconButton('panel', 'Toggle settings panel', () => { this.workspace.classList.toggle('side-collapsed'); this.render(); }, 'icon-btn bordered'));

    this.overviewHost = h('div', { class: 'canvas-host', style: { borderBottom: '1px solid var(--line)' }, 'data-export': 'skip' }, h('canvas'));
    this.rulerHost = h('div', { class: 'canvas-host', style: { borderBottom: '1px solid var(--line)' } }, h('canvas'));
    this.seqHost = h('div', { class: 'canvas-host seq-track', style: { borderBottom: '1px solid var(--line)' } }, h('canvas'));
    this.panelsHost = h('div', { class: 'canvas-host' });
    this.crosshair = h('div', { style: { position: 'absolute', top: '0', bottom: '0', width: '1px', background: 'var(--ink-3)', opacity: '.6', pointerEvents: 'none', zIndex: '3' }, hidden: true });
    this.brush = h('div', { style: { position: 'absolute', top: '0', bottom: '0', background: 'var(--accent-wash)', borderLeft: '1px solid var(--accent)', borderRight: '1px solid var(--accent)', pointerEvents: 'none', zIndex: '3' }, hidden: true });
    this.loadingLine = h('div', { class: 'loading-line', hidden: true });
    this.legend = h('div', { class: 'legend' });
    this.footer = h('div', { class: 'muted', style: { padding: '10px 16px', fontSize: '12px' } });
    this.scroll = h('div', { class: 'ws-scroll' }, this.loadingLine, this.overviewHost, this.legend, h('div', { style: { position: 'relative' } }, this.rulerHost, this.seqHost, this.panelsHost, this.crosshair, this.brush), this.footer);
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
    if (this.cfg.gene !== gene && isVariantField(this.cfg.groupField) && this.info) {
      this.cfg.groupField = defaultGroupField(categoricalFields()); // the site belongs to the previous window
      this.cfg.groups = null;
    }
    this.cfg.gene = gene;
    this.geneSelect.value = gene;
    this.data = null; this.pop = null; this.annotations = null; this.axis = null; this.seq = null;
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

  /** Keep only tracks this gene has; default to the first tracks of the preferred output. */
  applyGeneInfo() {
    const outputs = this.outputs();
    const valid = (key) => { const [o, i] = parseTrackKey(key); return !!(o && this.info.outputs[o] && i >= 0 && i < this.info.outputs[o].tracks.length); };
    const kept = (this.cfg.tracks || []).filter(valid);
    if (kept.length || (this.cfg.tracks && !this.cfg.tracks.length)) this.cfg.tracks = this.sortTracks(kept);
    else if (outputs.length) {
      const o = outputs.includes(this.preferOutput) ? this.preferOutput : outputs.includes('rna_seq') ? 'rna_seq' : outputs[0];
      this.cfg.tracks = this.trackMeta(o).slice(0, 6).map((m) => trackKey(o, m.index));
    } else this.cfg.tracks = [];
    if (!valid(this.cfg.popTrack)) this.cfg.popTrack = this.cfg.tracks[0] || (outputs.length ? trackKey(outputs[0], 0) : null);
  }

  outputs() { return Object.keys((this.info && this.info.outputs) || {}); }

  trackMeta(output) { return ((this.info && this.info.outputs[output]) || { tracks: [] }).tracks; }

  /** {output, index, meta} of a track key. */
  track(key) {
    const [output, index] = parseTrackKey(key);
    return { output, index, meta: this.trackMeta(output)[index] || {} };
  }

  trackLabel(key, short = false) {
    const { output, index, meta } = this.track(key);
    return `${outputLabel(output)} · ${(short && meta.short) || meta.label || `track ${index}`}`;
  }

  /** Panels follow the user's order (drag a panel's grip); new tracks are appended. */
  sortTracks(keys) { return [...new Set(keys)]; }

  /** Move the panel at ``from`` to ``to`` (no refetch: the data key ignores the track order). */
  moveTrack(from, to) {
    const tracks = [...(this.cfg.tracks || [])];
    if (from === to || from < 0 || to < 0 || from >= tracks.length || to >= tracks.length) return;
    tracks.splice(to, 0, ...tracks.splice(from, 1));
    this.cfg.tracks = tracks;
    this.persist();
    this.render();
  }

  /** Data panels in the order of ``cfg.tracks`` (responses come grouped by output). */
  orderPanels() {
    if (!this.data) return;
    const rank = new Map((this.cfg.tracks || []).map((k, i) => [k, i]));
    this.data.panels.sort((a, b) => (rank.get(a.trackKey) ?? 1e9) - (rank.get(b.trackKey) ?? 1e9));
  }

  /** Selected tracks grouped by output: [[output, [indices]], …] in panel order. */
  tracksByOutput() {
    const groups = new Map();
    for (const key of this.cfg.tracks || []) {
      const [o, i] = parseTrackKey(key);
      if (!groups.has(o)) groups.set(o, []);
      groups.get(o).push(i);
    }
    return [...groups.entries()];
  }

  shownOutputs() {
    if (this.cfg.mode === 'population') return this.cfg.popTrack ? [parseTrackKey(this.cfg.popTrack)[0]] : [];
    return this.tracksByOutput().map(([o]) => o);
  }

  domainLength() {
    if (this.cfg.coords === 'aligned' && this.axis) return this.axis.expanded_length;
    const lengths = this.shownOutputs().map((o) => (this.info.outputs[o] || {}).length || 0);
    if (this.cfg.coords === 'haplotype') return Math.max(0, ...lengths) || this.info.length || 1;
    return this.info.length || Math.max(0, ...lengths) || 1;
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
    this.data = null; this.pop = null; this.seq = null;
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
    const base = [c.gene, c.coords];
    const extra = [c.showRef, c.showObserved, c.obsScale];
    const tracks = [...(c.tracks || [])].sort(); // reordering panels must not refetch
    if (c.mode === 'individuals') return JSON.stringify([...base, 'i', tracks, this.seriesSpec(), ...extra]);
    if (c.mode === 'groups') return JSON.stringify([...base, 'g', tracks, c.hap, c.groupField, this.groupValues(), state.filters, ...extra]);
    return JSON.stringify([...base, 'p', c.popTrack, c.hap, c.groupField, state.filters, c.popRows]);
  }

  individuals() { return state.pinned.slice(0, MAX_INDIVIDUALS); }

  seriesSpec() {
    const haps = this.cfg.hap === 'both' ? ['H1', 'H2'] : [this.cfg.hap];
    return this.individuals().flatMap((s) => haps.map((hp) => `${s}:${hp}`)).join(',');
  }

  /** Outputs of this gene with an AlphaGenome prediction of the reference window. */
  hasReference(output) { return ((this.info && this.info.reference_outputs) || []).includes(output); }

  showReference(output) { return this.cfg.showRef && this.cfg.mode !== 'population' && this.hasReference(output); }

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
    if (isVariantField(this.cfg.groupField)) return []; // the server forms the genotype classes
    const counts = this.cohortCounts();
    const all = [...counts.keys()].sort((a, b) => counts.get(b) - counts.get(a) || a.localeCompare(b));
    const chosen = (this.cfg.groups || []).filter((g) => counts.has(g));
    return (chosen.length ? chosen : all).slice(0, 8);
  }

  invalidate() {
    this.clearHover();
    if (this.info) this.viewport.setDomain(this.domainLength());
    this.persist();
    this.request.cancel?.();
    this.requestNow();
    this.requestSequence.cancel?.();
    this.requestSequenceNow();
  }

  persist() {
    const c = this.cfg;
    setLocus({ gene: c.gene, output: null, coords: c.coords, mode: c.mode, hap: c.hap, tracks: c.tracks, groupField: c.groupField, groups: c.groups, yScale: c.yScale, sharedY: c.sharedY, diff: c.diff, envelope: c.envelope, band: c.band, popTrack: c.popTrack, popRows: c.popRows, seqSource: c.seqSource, showRef: c.showRef, showObserved: c.showObserved, obsScale: c.obsScale, start: this.viewport.start, end: this.viewport.end });
    updateRouteParams({ gene: c.gene, coords: c.coords !== 'reference' ? c.coords : null, mode: c.mode !== 'individuals' ? c.mode : null, start: this.viewport.start, end: this.viewport.end });
  }

  onViewChange() {
    this.renderLocus();
    this.render();
    this.persistDebounced();
    this.request();
    this.requestSequence();
  }

  persistDebounced = debounce(() => this.persist(), 300);

  request = debounce(() => this.requestNow(), 120);

  needsFetch(data, key, margin) {
    if (!data || data.key !== key) return true;
    const v = this.viewport;
    if (v.start < data.start || (v.end > data.end && data.end < v.length)) return true;
    const bins = Math.max(1, data.edges.length - 1);
    const dataBp = (data.end - data.start) / bins;
    const wantBp = v.span / this.geom().width;
    if (dataBp > Math.max(1, wantBp) * 1.6) return true; // too coarse
    if (margin && dataBp < wantBp / 4 && bins > 3000) return true; // far too fine
    return false;
  }

  async requestNow() {
    if (!this.info || !this.outputs().length) { this.render(); return; }
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
    const common = { gene: this.cfg.gene, coords: this.cfg.coords, start, end, bins };
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
      const groups = this.tracksByOutput();
      if (mode === 'individuals') {
        // Pinned haplotypes plus the reference prediction (enough on its own when nothing is pinned).
        const specs = groups.map(([output]) => [this.seriesSpec(), this.showReference(output) ? `${REF}:H1` : ''].filter(Boolean).join(','));
        if (!groups.length || specs.some((spec) => !spec)) { this.data = null; this.render(); return; }
        // One request per output, in parallel; each panel keeps its own bin edges.
        const responses = await Promise.all(groups.map(([output, tracks], gi) => api(`${ds()}/signal`, { params: { ...common, output, tracks: tracks.join(','), series: specs[gi] }, signal })));
        this.data = this.toItems(responses, groups, key, 'series');
        const failed = responses.flatMap((res) => res.series.filter((s) => s.error));
        if (failed.length) toast(`${failed.length} series failed: ${failed[0].error}`, 'error', 7000);
      } else if (mode === 'groups') {
        if (!groups.length) { this.data = null; this.render(); return; }
        // Sequential: each output may start a long cohort job with its own progress.
        const responses = [];
        for (const [output, tracks] of groups) {
          responses.push(await apiJob(`${ds()}/groups`, { params: { ...common, output, tracks: tracks.join(','), haps: this.cfg.hap === 'both' ? 'H1,H2' : this.cfg.hap, field: this.cfg.groupField, groups: this.groupValues().join(','), filters: state.filters }, signal, onProgress }));
          if (overlay) { overlay.remove(); overlay = null; }
        }
        // The reference prediction is drawn with the group means (same bins: same start/end/bins).
        const refs = await Promise.all(groups.map(([output, tracks]) => (this.showReference(output)
          ? api(`${ds()}/signal`, { params: { ...common, output, tracks: tracks.join(','), series: `${REF}:H1` }, signal }).catch((err) => { if (isAbort(err)) throw err; return null; })
          : null)));
        this.data = this.toItems(responses, groups, key, 'groups', refs);
      } else {
        if (!this.cfg.popTrack) { this.pop = null; this.render(); return; }
        const hap = this.cfg.hap === 'both' ? 'H1+H2' : this.cfg.hap;
        const [output, track] = parseTrackKey(this.cfg.popTrack);
        const res = await apiJob(`${ds()}/population`, { params: { ...common, output, track, hap, field: this.cfg.groupField, filters: state.filters, max_rows: this.cfg.popRows }, signal, onProgress });
        res.key = key;
        this.pop = res;
      }
      this.status('');
      if (this.data && this.cfg.showObserved && mode !== 'population') {
        this.render();
        if (overlay) { overlay.remove(); overlay = null; }
        await this.loadObserved(groups, common, signal, onProgress);
      }
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

  /** Observed ENCODE / FANTOM5 signal of the shown tracks (same bins as the panels). */
  async loadObserved(groups, common, signal, onProgress) {
    const data = this.data;
    this.status('Loading observed data…');
    const results = await Promise.all(groups.map(([output, tracks]) => apiJob(`${ds()}/observed`, { params: { ...common, output, tracks: tracks.join(','), scale: this.cfg.obsScale }, signal, onProgress })
      .catch((err) => { if (isAbort(err)) throw err; return { error: err.message, tracks: tracks.map((t) => ({ track: t, available: false, reason: err.message })) }; })));
    if (this.data !== data) return;
    data.observed = new Map();
    results.forEach((res, gi) => {
      const [output] = groups[gi];
      for (const item of res.tracks) data.observed.set(trackKey(output, item.track), { ...item, edges: res.edges, start: res.start, end: res.end });
    });
    this.status('');
  }

  /**
   * Responses (one per output) -> one entry per panel, in track order. Each panel item holds the
   * 1-D arrays of that track; ``items`` (legend) is the union of entities across outputs.
   */
  toItems(responses, groups, key, kind, refs = []) {
    const legend = new Map();
    const panels = [];
    const refEntry = (s) => ({ src: s, key: REF, entity: REF, label: 'Reference genome', sub: 'AlphaGenome on hg38', reference: true });
    responses.forEach((res, gi) => {
      const [output, indices] = groups[gi];
      const entries = [];
      if (kind === 'series') {
        for (const s of res.series) {
          if (s.error) continue;
          if (s.sample === REF) { entries.push(refEntry(s)); continue; }
          entries.push({ src: s, key: `${s.sample}:${s.haplotype}`, entity: s.sample, label: this.cfg.hap === 'both' ? `${s.sample} ${s.haplotype}` : s.sample, sub: this.cfg.hap === 'both' ? '' : (s.haplotype === 'H1+H2' ? 'diploid mean' : s.haplotype), dash: s.haplotype === 'H2' && this.cfg.hap === 'both' });
        }
      } else {
        for (const g of res.groups) entries.push({ src: g, key: g.group, entity: g.group, label: g.group, sub: `n=${fmtInt(g.samples)}`, haplotypes: g.haplotypes, failed: g.failed });
        const ref = refs[gi] && refs[gi].series.find((s) => !s.error);
        if (ref) entries.push(refEntry(ref));
      }
      for (const { src, ...it } of entries) if (!legend.has(it.key)) legend.set(it.key, it);
      indices.forEach((index, ti) => {
        const pick = (arr) => (arr ? arr[ti] : null);
        panels.push({
          trackKey: trackKey(output, index), output, index, start: res.start, end: res.end, edges: res.edges,
          items: entries.map(({ src, ...it }) => ({ ...it, mean: pick(src.mean), min: pick(src.min), max: pick(src.max), upper: pick(src.upper), lower: pick(src.lower) })),
        });
      });
    });
    const items = [...legend.values()];
    const entities = items.filter((it) => !it.reference).map((it) => it.entity);
    this.colors.sync(kind === 'series' ? [...new Set(entities)] : entities);
    const first = panels[0] || { start: 0, end: 0, edges: [0] };
    const start = Math.min(...panels.map((pn) => pn.start), first.start);
    const end = Math.max(...panels.map((pn) => pn.end), first.end);
    return { key, kind, start, end, edges: first.edges, panels, items };
  }

  // ------------------------------------------------------------------ rendering
  render() {
    if (!this.info) return;
    this.renderLocus();
    this.drawOverview();
    this.drawRuler();
    this.drawSequence();
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

  /** Gene x position (reference offset -> canvas x in the current coordinates). */
  geneX = (p) => this.xOf(this.fromGenomicOffset(p));

  drawRuler() {
    const g = this.geom();
    const canvas = this.rulerHost.querySelector('canvas');
    const v = this.viewport;
    const rows = packGenes(this.annotationGenes(), { start: v.start, end: v.end, pxPerBp: g.width / v.span, toPos: (p) => this.fromGenomicOffset(p) });
    const H = 30 + geneLanesHeight(rows);
    const ctx = setupCanvas(canvas, g.total, H);
    const t = theme();
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
    drawGeneLanes(ctx, t, rows, { top: 30, left: g.left, width: g.width, xOf: this.geneX, highlight: this.cfg.gene });
    this.geneRows = rows;
  }

  // ------------------------------------------------------------------ sequence lane
  /** Haplotypes counted by the sequence lane ('' = the reference). */
  seqRows() {
    if (this.cfg.seqSource !== 'pinned') return '';
    const hap = this.cfg.hap === 'H1' || this.cfg.hap === 'H2' ? this.cfg.hap : 'H1+H2';
    return state.pinned.slice(0, MAX_SEQ_SAMPLES).map((s) => `${s}:${hap}`).join(',');
  }

  requestSequence = debounce(() => this.requestSequenceNow(), 140);

  async requestSequenceNow() {
    if (!this.info || this.cfg.seqSource === 'off') { this.drawSequence(); return; }
    const v = this.viewport;
    const g = this.geom();
    const key = JSON.stringify([this.cfg.gene, this.cfg.coords, this.seqRows()]);
    const d = this.seq;
    if (d && d.key === key && v.start >= d.start && (v.end <= d.end || d.end >= d.domain)) {
      const wantBp = v.span / g.width;
      const binBp = (d.end - d.start) / Math.max(1, d.edges.length - 1);
      const tooCoarse = binBp > Math.max(1, wantBp) * 1.6;
      const tooFine = binBp < wantBp / 4 && d.edges.length > 3000;
      if (!tooCoarse && !tooFine) { this.drawSequence(); return; }
    }
    const margin = v.span * 0.5;
    const start = Math.max(0, Math.floor(v.start - margin));
    const end = Math.min(v.length, Math.ceil(v.end + margin));
    const bins = Math.max(32, Math.min(8192, Math.round(g.width * (end - start) / v.span)));
    const signal = this.seqLatest.next();
    try {
      const res = await api(`${ds()}/composition`, { params: { gene: this.cfg.gene, coords: this.cfg.coords, rows: this.seqRows(), start, end, bins }, signal });
      res.key = key;
      this.seq = res;
      if (res.errors && res.errors.length) toast(`Sequence lane: ${res.errors[0]}`, 'error', 6000);
    } catch (err) {
      if (isAbort(err)) return;
      this.seq = null;
      this.seqError = err.message;
    }
    this.drawSequence();
  }

  /**
   * Letter frequency per base: a frequency-scaled stack of letters when zoomed in (a logo of the
   * counted haplotypes, or the plain reference), stacked base-composition bars when zoomed out.
   */
  drawSequence() {
    const host = this.seqHost;
    host.hidden = this.cfg.seqSource === 'off';
    if (host.hidden || !this.info) return;
    const g = this.geom();
    const ctx = setupCanvas(host.querySelector('canvas'), g.total, SEQ_H);
    const t = theme();
    const d = this.seq;
    const v = this.viewport;
    const top = 18; const bottom = SEQ_H - 4; const hgt = bottom - top;
    const mw = this.modelWindow();
    if (mw && mw.end > v.start && mw.start < v.end) {
      const a = Math.max(g.left, this.xOf(mw.start)); const b = Math.min(g.left + g.width, this.xOf(mw.end));
      ctx.fillStyle = t.modelWindow; ctx.fillRect(a, 0, b - a, SEQ_H);
    }
    ctx.font = `10.5px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillStyle = t.ink3;
    ctx.fillText('sequence', g.left - 8, top + hgt / 2);
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    if (!d) { ctx.fillText(this.seqError ? `Sequence unavailable: ${this.seqError}` : 'Loading sequence…', g.left + 4, 4); return; }
    const pxPerBp = g.width / v.span;
    const binBp = (d.end - d.start) / Math.max(1, d.edges.length - 1);
    const letters = binBp <= 1 && pxPerBp >= SEQ_LETTERS_PX;
    const n = d.haplotypes.length;
    const what = d.source === 'haplotypes' ? `letter frequency over ${n} pinned haplotype${n === 1 ? '' : 's'}` : 'reference genome';
    const how = letters ? '' : binBp <= 1 ? ' · zoom in for letters' : ` · base composition per ${fmtBp(binBp)} bin`;
    ctx.fillText(`${what}${how}`, g.left + 4, 3);
    const colors = [t.bases.A, t.bases.C, t.bases.G, t.bases.T, t.bases.N, t.surface2];
    const e = d.edges;
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, top, g.width, hgt); ctx.clip();
    let i0 = 0; while (i0 < e.length - 2 && e[i0 + 1] <= v.start) i0++;
    if (letters) {
      ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
      const fontPx = 40;
      ctx.font = `700 ${fontPx}px ${t.mono}`;
      for (let i = i0; i < e.length - 1 && e[i] < v.end; i++) {
        const x = this.xOf(e[i]);
        const cx = x + pxPerBp / 2;
        const sx = Math.min(1, (pxPerBp * 0.92) / (fontPx * 0.62));
        let y = bottom;
        // Largest frequency at the top, as in a sequence logo.
        const order = [0, 1, 2, 3, 4, 5].filter((k) => d.freq[k][i] > 0).sort((a, b) => d.freq[a][i] - d.freq[b][i]);
        for (const k of order) {
          const f = d.freq[k][i];
          const lh = f * hgt;
          if (k === 5) { ctx.fillStyle = t.ink3; ctx.fillRect(x + 1, y - lh / 2 - 1, Math.max(1, pxPerBp - 2), 2); y -= lh; continue; }
          if (lh >= 3) {
            ctx.save();
            ctx.translate(cx, y);
            ctx.scale(sx, lh / (fontPx * 0.72));
            ctx.fillStyle = colors[k];
            ctx.fillText(d.letters[k], 0, 0);
            ctx.restore();
          } else { ctx.fillStyle = colors[k]; ctx.fillRect(x + 1, y - lh, Math.max(1, pxPerBp - 2), lh); }
          y -= lh;
        }
      }
    } else {
      for (let i = i0; i < e.length - 1 && e[i] < v.end; i++) {
        const x0 = this.xOf(e[i]); const x1 = this.xOf(e[i + 1]);
        const w = Math.max(1, x1 - x0 + 0.5);
        let y = bottom;
        for (let k = 0; k < 6; k++) {
          const lh = d.freq[k][i] * hgt;
          if (lh <= 0) continue;
          ctx.fillStyle = colors[k];
          ctx.fillRect(x0, y - lh, w, lh);
          y -= lh;
        }
      }
    }
    ctx.restore();
  }

  seqTooltip(pos) {
    const d = this.seq;
    if (!d) return '';
    const bi = binIndex(d.edges, pos);
    if (bi < 0) return '';
    const binBp = d.edges[bi + 1] - d.edges[bi];
    const names = { '-': 'gap / deleted', N: 'N / other' };
    const rows = [...d.letters].map((l, k) => [l, d.freq[k][bi]]).filter(([, f]) => f > 0).sort((a, b) => b[1] - a[1])
      .map(([l, f]) => `<div class="tt-row"><span class="tt-key">${escapeHtml(names[l] || l)}</span><span class="tt-val">${fmtNum(f * 100, 3)}%</span></div>`);
    const ref = d.reference && binBp === 1 ? d.reference[pos - d.start] : null;
    const sub = d.source === 'haplotypes' ? `${d.haplotypes.length} haplotypes` : 'reference';
    return `<div class="tt-sub" style="margin:-2px 0 6px">${binBp > 1 ? `composition of ${fmtInt(binBp)} bp · ` : ''}${escapeHtml(sub)}${ref ? ` · reference ${escapeHtml(ref)}` : ''}</div>${rows.join('')}`;
  }

  ensurePanels(n, height) {
    const host = this.panelsHost;
    const obs = this.cfg.showObserved ? '1' : '0';
    if (host.dataset.kind !== 'panels' || host.children.length !== n || host.dataset.obs !== obs) {
      clear(host);
      host.dataset.kind = 'panels';
      host.dataset.obs = obs;
      for (let i = 0; i < n; i++) {
        const label = h('div', { class: 'lane-label clickable', style: { top: '6px' }, title: 'Track card: ontology terms, assay, experiments', onpointerdown: (e) => e.stopPropagation(), onclick: (e) => { e.stopPropagation(); this.openTrack(label.dataset.key); } });
        const obsLabel = h('div', { class: 'lane-label clickable obs-label', style: { top: `${height + 4}px` }, title: 'The observed experiments behind this track', onpointerdown: (e) => e.stopPropagation(), onclick: (e) => { e.stopPropagation(); this.openTrack(label.dataset.key); } });
        const obsStats = h('div', { class: 'obs-stats', style: { top: `${height + 4}px` } });
        const panel = h('div', { class: 'panel has-grip' }, label, h('canvas', { class: 'pred' }), this.cfg.showObserved ? [obsLabel, obsStats, h('canvas', { class: 'obs' })] : null);
        panel.prepend(this.panelGrip(panel));
        host.appendChild(panel);
      }
    }
    return [...host.children].map((p) => ({ el: p, label: p.querySelector('.lane-label'), canvas: p.querySelector('canvas.pred'), obsLabel: p.querySelector('.obs-label'), obsStats: p.querySelector('.obs-stats'), obsCanvas: p.querySelector('canvas.obs'), height }));
  }

  /**
   * Drag handle of a track panel: drag it up or down to reorder the panels (the others slide to make
   * room; the page scrolls near its edges; Escape cancels). Arrow keys move a focused handle by one panel.
   */
  panelGrip(panel) {
    const grip = h('button', { class: 'panel-grip', type: 'button', title: 'Drag to reorder (or focus and use ↑ / ↓)', 'aria-label': 'Reorder track' }, icon('grip', 14));
    const indexOf = () => [...this.panelsHost.children].indexOf(panel);
    let drag = null;
    const reset = () => {
      if (!drag) return;
      for (const el of drag.siblings) { el.style.transform = ''; el.classList.remove('reorder-shift'); }
      panel.style.transform = '';
      panel.classList.remove('reordering');
      this.panelsHost.classList.remove('reordering-host');
      window.removeEventListener('keydown', drag.onKey, true);
      cancelAnimationFrame(drag.raf);
      drag = null;
    };
    grip.addEventListener('pointerdown', (e) => {
      e.stopPropagation(); // not a pan of the viewport
      if (e.button !== 0 || (this.cfg.tracks || []).length < 2) return;
      e.preventDefault();
      grip.setPointerCapture(e.pointerId);
      hideTooltip(); this.clearHover();
      const siblings = [...this.panelsHost.children];
      const from = siblings.indexOf(panel);
      const onKey = (k) => { if (k.key === 'Escape') { k.stopPropagation(); reset(); } };
      window.addEventListener('keydown', onKey, true);
      drag = { y: e.clientY, lastY: e.clientY, scroll0: this.scroll.scrollTop, raf: 0, from, to: from, siblings: siblings.filter((el) => el !== panel), mids: siblings.map((el) => { const r = el.getBoundingClientRect(); return r.top + r.height / 2; }), height: panel.getBoundingClientRect().height, onKey };
      for (const el of drag.siblings) el.classList.add('reorder-shift');
      panel.classList.add('reordering');
      this.panelsHost.classList.add('reordering-host');
      drag.raf = requestAnimationFrame(autoScroll);
    });
    // Positions are taken at pointerdown; scrolling since then shifts the pointer by ``ds`` in the page.
    const update = () => {
      const dy = drag.lastY - drag.y + (this.scroll.scrollTop - drag.scroll0);
      panel.style.transform = `translateY(${dy}px)`;
      const center = drag.mids[drag.from] + dy;
      let to = drag.from;
      while (to < drag.mids.length - 1 && center > drag.mids[to + 1]) to++;
      while (to > 0 && center < drag.mids[to - 1]) to--;
      drag.to = to;
      // Panels between the old and new slot slide by one panel height.
      [...this.panelsHost.children].forEach((el, i) => {
        if (el === panel) return;
        const shift = drag.from < to && i > drag.from && i <= to ? -drag.height : drag.from > to && i >= to && i < drag.from ? drag.height : 0;
        el.style.transform = shift ? `translateY(${shift}px)` : '';
      });
    };
    const autoScroll = () => {
      if (!drag) return;
      const r = this.scroll.getBoundingClientRect();
      const edge = 48;
      const step = drag.lastY < r.top + edge ? -Math.ceil((r.top + edge - drag.lastY) / 4) : drag.lastY > r.bottom - edge ? Math.ceil((drag.lastY - r.bottom + edge) / 4) : 0;
      if (step) { const before = this.scroll.scrollTop; this.scroll.scrollTop += step; if (this.scroll.scrollTop !== before) update(); }
      drag.raf = requestAnimationFrame(autoScroll);
    };
    grip.addEventListener('pointermove', (e) => {
      e.stopPropagation();
      if (!drag) return;
      drag.lastY = e.clientY;
      update();
    });
    const finish = (e) => {
      e.stopPropagation();
      if (!drag) return;
      const { from, to } = drag;
      reset();
      this.moveTrack(from, to);
    };
    grip.addEventListener('pointerup', finish);
    grip.addEventListener('pointercancel', (e) => { e.stopPropagation(); reset(); });
    grip.addEventListener('click', (e) => e.stopPropagation());
    grip.addEventListener('keydown', (e) => {
      if (e.key !== 'ArrowUp' && e.key !== 'ArrowDown') return;
      e.preventDefault(); e.stopPropagation();
      const from = indexOf();
      const to = from + (e.key === 'ArrowUp' ? -1 : 1);
      this.moveTrack(from, to);
      const moved = this.panelsHost.children[Math.max(0, Math.min(to, this.panelsHost.children.length - 1))];
      const next = moved && moved.querySelector('.panel-grip');
      if (next) next.focus();
    });
    return grip;
  }

  openTrack(key) {
    if (!key) return;
    const { output, index, meta } = this.track(key);
    const observed = this.data && this.data.observed ? this.data.observed.get(key) : null;
    openTrackCard({ output, index, meta: meta.metadata || {}, label: this.trackLabel(key), gene: this.cfg.gene, observed });
  }

  itemColor(it) { return it.reference ? theme().ink : this.colors.color(it.entity); }

  transform(v) {
    if (this.cfg.yScale === 'log') return Math.sign(v) * Math.log1p(Math.abs(v));
    return v;
  }

  visibleRange(panel) {
    const v = this.viewport;
    const e = panel.edges;
    let i0 = 0; let i1 = e.length - 1;
    while (i0 < e.length - 1 && e[i0 + 1] <= v.start) i0++;
    while (i1 > i0 && e[i1 - 1] >= v.end) i1--;
    return [i0, i1];
  }

  setPanelLabel(p, key) {
    const { meta } = this.track(key);
    p.label.textContent = this.trackLabel(key);
    p.label.dataset.key = key || '';
    p.label.title = `${Object.entries(meta.metadata || {}).map(([k, val]) => `${k}: ${val}`).join('\n')}\n\nClick for the track card (ontology terms, assay, observed experiments)`;
  }

  drawPanels() {
    const g = this.geom();
    const data = this.data;
    const tracks = this.cfg.tracks || [];
    const panels = this.ensurePanels(Math.max(1, tracks.length), PANEL_H);
    const t = theme();
    if (!data || !data.items.length || !data.panels.length) {
      panels.forEach((p, i) => {
        const ctx = setupCanvas(p.canvas, g.total, PANEL_H);
        if (p.obsCanvas) setupCanvas(p.obsCanvas, g.total, OBS_H);
        if (tracks.length) this.setPanelLabel(p, tracks[i]); else p.label.textContent = 'No tracks selected';
        ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.textAlign = 'center';
        const msg = !tracks.length ? 'Choose tracks of any output in the panel →' : this.cfg.mode === 'individuals' && !this.individuals().length ? 'Pin individuals on the Samples page (or add them in the panel →)' : (data ? 'No data' : 'Loading…');
        ctx.fillText(msg, g.left + g.width / 2, PANEL_H / 2);
      });
      return;
    }
    const band = this.showBand();
    this.orderPanels();
    const views = data.panels.map((panel) => {
      const [i0, i1] = this.visibleRange(panel);
      const items = this.displayItems(panel).filter((it) => !this.hidden.has(it.key));
      let lo = 0; let hi = 0;
      for (const it of items) {
        const top = this.cfg.envelope && it.max ? it.max : it.mean;
        const upper = band && it.upper ? it.upper : null;
        const bottom = band && it.lower ? it.lower : (this.cfg.envelope && it.min ? it.min : it.mean);
        for (let i = i0; i < i1; i++) {
          const a = top[i]; if (a > hi) hi = a;
          if (upper) { const u = upper[i]; if (u > hi) hi = u; }
          const b = bottom[i]; if (b < lo) lo = b;
        }
      }
      return { panel, items, i0, i1, range: [this.transform(lo), this.transform(hi)] };
    });
    // A shared y-scale only makes sense within one output (units differ between assays).
    const shared = new Map();
    if (this.cfg.sharedY) {
      for (const { panel, range } of views) {
        const cur = shared.get(panel.output);
        shared.set(panel.output, cur ? [Math.min(cur[0], range[0]), Math.max(cur[1], range[1])] : range);
      }
    }
    views.forEach(({ panel, items, i0, i1, range }, pi) => {
      const p = panels[pi];
      if (!p) return;
      this.setPanelLabel(p, panel.trackKey);
      const [ylo, yhiRaw] = shared.get(panel.output) || range;
      const yhi = yhiRaw > ylo ? yhiRaw * 1.06 : ylo + 1;
      this.drawPanel(p.canvas, t, g, panel, items, ylo, yhi, i0, i1);
      p.yRange = [ylo, yhi];
      if (p.obsCanvas) this.drawObserved(p, t, g, panel);
    });
    this.panelGeom = { top: 0, height: PANEL_H, count: data.panels.length };
  }

  /** Observed signal lane of one track (own y-scale: units are the source's, not AlphaGenome's). */
  drawObserved(p, t, g, panel) {
    const H = OBS_H;
    const ctx = setupCanvas(p.obsCanvas, g.total, H);
    const obs = this.data.observed ? this.data.observed.get(panel.trackKey) : null;
    const color = css('--observed') || t.ink2;
    this.renderAgreement(p, null);
    if (!obs) {
      p.obsLabel.textContent = 'Observed · loading…';
      return;
    }
    if (!obs.available) {
      p.obsLabel.textContent = 'Observed · not available';
      ctx.fillStyle = t.ink3; ctx.font = `11.5px ${t.font}`; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText(obs.reason || 'No observed data', g.left + 8, H / 2 + 6);
      return;
    }
    const n = (obs.items || []).length;
    p.obsLabel.textContent = `Observed · ${obs.provider} · ${obs.units} · mean of ${n}${obs.total > n ? `/${obs.total}` : ''} ${obs.provider === 'FANTOM5' ? 'librar' + (n === 1 ? 'y' : 'ies') : 'experiment' + (n === 1 ? '' : 's')}${obs.match === 'name' ? ' (by name)' : ''}`;
    const e = obs.edges;
    const v = this.viewport;
    let i0 = 0; let i1 = e.length - 1;
    while (i0 < e.length - 1 && e[i0 + 1] <= v.start) i0++;
    while (i1 > i0 && e[i1 - 1] >= v.end) i1--;
    this.renderAgreement(p, this.agreement(panel, obs, i0, i1));
    let hi = 0;
    const top = this.cfg.envelope && obs.max ? obs.max : obs.mean;
    for (let i = i0; i < i1; i++) { const a = top[i]; if (a > hi) hi = a; }
    const ylo = 0; const yhi = hi > 0 ? this.transform(hi) * 1.08 : 1;
    const y0 = 30; const y1 = H - 6;
    const sy = (val) => y1 - ((this.transform(Math.max(0, val)) - ylo) / (yhi - ylo)) * (y1 - y0);
    ctx.fillStyle = t.ink3; ctx.font = `10.5px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillText(fmtNum(this.cfg.yScale === 'log' ? Math.expm1(yhi) : yhi, 2), g.left - 6, y0);
    ctx.fillText('0', g.left - 6, y1);
    ctx.strokeStyle = t.lineStrong; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(g.left, y1 + 0.5); ctx.lineTo(g.left + g.width, y1 + 0.5); ctx.stroke();
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, 0, g.width, H); ctx.clip();
    const binBp = (obs.end - obs.start) / Math.max(1, e.length - 1);
    const pxPerBin = (binBp / v.span) * g.width;
    ctx.fillStyle = withAlpha(color, 0.75);
    for (let i = Math.max(0, i0 - 1); i < Math.min(e.length - 1, i1 + 1); i++) {
      const val = (this.cfg.envelope && obs.max ? obs.max : obs.mean)[i];
      if (!Number.isFinite(val) || val <= 0) continue;
      const xa = this.xOf(e[i]); const xb = this.xOf(e[i + 1]);
      const y = sy(val);
      ctx.fillRect(xa, y, Math.max(pxPerBin > 2 ? xb - xa - 0.5 : 1, 1), y1 - y);
    }
    ctx.restore();
  }

  /**
   * Agreement of a prediction with the observed lane over the bins in view: Pearson r (linear and of
   * log1p values), the share of the observed peak bins (top decile) that are also predicted peak bins,
   * and the ratio of means (predicted / observed). The ratio is only flagged for FANTOM5 CAGE in TPM, the
   * units AlphaGenome's CAGE tracks are trained on (a large ratio there is a magnitude bias, like the
   * exon/promoter CAGE bias on TYRP1); ENCODE signal tracks use other units. Compared against the
   * reference-genome prediction when shown, else the first series.
   */
  agreement(panel, obs, i0, i1) {
    const pred = panel.items.find((it) => it.reference && it.mean) || panel.items.find((it) => it.mean);
    if (!pred || !obs.mean || pred.mean.length !== obs.mean.length) return null;
    const xs = []; const ys = []; const at = [];
    for (let i = i0; i < i1; i++) {
      const x = pred.mean[i]; const y = obs.mean[i];
      if (Number.isFinite(x) && Number.isFinite(y)) { xs.push(x); ys.push(y); at.push(i); }
    }
    if (xs.length < 3) return null;
    const pearson = (a, b) => {
      const n = a.length; let ma = 0; let mb = 0;
      for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
      ma /= n; mb /= n;
      let sab = 0; let saa = 0; let sbb = 0;
      for (let i = 0; i < n; i++) { const da = a[i] - ma; const db = b[i] - mb; sab += da * db; saa += da * da; sbb += db * db; }
      return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : NaN;
    };
    const log = (a) => a.map((v) => Math.log1p(Math.max(0, v)));
    // Top-decile bins (with signal) of each; ties at the threshold are included.
    const top = (a) => {
      const k = Math.max(1, Math.round(a.length * 0.1));
      const threshold = [...a].sort((u, v) => v - u)[k - 1];
      return new Set(a.map((v, i) => (v >= threshold && v > 0 ? i : -1)).filter((i) => i >= 0));
    };
    const tp = top(xs); const to = top(ys);
    const shared = to.size ? [...to].filter((i) => tp.has(i)).length / to.size : NaN;
    const argmax = (a) => a.reduce((best, v, i) => (v > a[best] ? i : best), 0);
    const centre = (i) => (obs.edges[at[i]] + obs.edges[at[i] + 1]) / 2;
    const meanOf = (a) => a.reduce((acc, v) => acc + v, 0) / a.length;
    const mo = meanOf(ys);
    return {
      against: pred.reference ? 'reference genome' : pred.label, n: xs.length,
      r: pearson(xs, ys), rlog: pearson(log(xs), log(ys)), shared,
      peak: Math.abs(centre(argmax(xs)) - centre(argmax(ys))),
      ratio: mo > 0 ? meanOf(xs) / mo : NaN,
      sameUnits: obs.provider === 'FANTOM5' && this.cfg.obsScale === 'tpm',
    };
  }

  renderAgreement(p, a) {
    if (!p.obsStats) return;
    clear(p.obsStats);
    p.obsStats.hidden = !a;
    if (!a) return;
    const weak = !(a.rlog >= 0.3);
    const biased = a.sameUnits && Number.isFinite(a.ratio) && Math.abs(Math.log2(a.ratio)) > 2;
    p.obsStats.classList.toggle('warn', weak || biased);
    p.obsStats.append(
      h('span', { class: weak ? 'bad' : '' }, `r ${fmtNum(a.r, 2)} · log r ${fmtNum(a.rlog, 2)}`),
      h('span', null, `peaks ${fmtPct(a.shared, 0)} shared`),
      h('span', { class: biased ? 'bad' : a.sameUnits ? '' : 'muted' }, `pred/obs ×${fmtNum(a.ratio, 2)}`));
    p.obsStats.title = `Observed vs AlphaGenome (${a.against}) over the ${fmtInt(a.n)} bins in view:\n`
      + `Pearson r ${fmtNum(a.r, 3)} (of log1p values: ${fmtNum(a.rlog, 3)})\n`
      + `${fmtPct(a.shared, 0)} of the observed top-decile bins are also predicted top-decile bins; the two maxima are ${fmtBp(a.peak)} apart\n`
      + `mean predicted / mean observed ${fmtNum(a.ratio, 3)}${a.sameUnits ? '' : ' (different units: not flagged)'}\n`
      + 'Flagged when log r < 0.3, or for FANTOM5 TPM when the ratio is beyond 4× either way. Zoom to a gene or exon to judge it locally.';
  }

  /** The ±SD band is within-group spread; it is hidden when showing differences between groups. */
  showBand() { return this.cfg.mode === 'groups' && this.cfg.band && !this.cfg.diff; }

  /** Group means as-is, or (when "difference" is on) minus the haplotype-weighted cohort mean. */
  displayItems(panel) {
    if (!(this.cfg.mode === 'groups' && this.cfg.diff && this.data.kind === 'groups' && panel.items.length > 1)) return panel.items;
    if (panel.diffItems) return panel.diffItems;
    const weights = panel.items.map((it) => (it.reference ? 0 : it.haplotypes || 1));
    const n = panel.items[0].mean.length;
    const baseline = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      let acc = 0; let w = 0;
      panel.items.forEach((it, k) => { const v = it.mean[i]; if (weights[k] && Number.isFinite(v)) { acc += v * weights[k]; w += weights[k]; } });
      baseline[i] = w ? acc / w : NaN;
    }
    const sub = (a) => a && a.map((v, i) => v - baseline[i]);
    panel.diffItems = panel.items.map((it) => ({ ...it, mean: sub(it.mean), min: sub(it.min), max: sub(it.max), upper: sub(it.upper), lower: sub(it.lower) }));
    return panel.diffItems;
  }

  drawPanel(canvas, t, g, panel, items, ylo, yhi, i0, i1) {
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
    const e = panel.edges;
    const binBp = (panel.end - panel.start) / Math.max(1, e.length - 1);
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
      if (it.reference) continue;
      const color = this.itemColor(it);
      const dim = highlight && highlight !== it.key;
      if (this.showBand() && it.upper) area(it.lower, it.upper, withAlpha(color, dim ? 0.04 : 0.14));
      else if (this.cfg.envelope && binBp > 1.5 && it.min) area(it.min, it.max, withAlpha(color, dim ? 0.04 : 0.16));
    }
    // The reference prediction goes on top, thin and dark, so individuals' deviations stand out.
    for (const it of [...items.filter((x) => !x.reference), ...items.filter((x) => x.reference)]) {
      const color = this.itemColor(it);
      const dim = highlight && highlight !== it.key;
      line(it.mean, dim ? withAlpha(color, 0.25) : color, it.dash || it.reference, highlight === it.key ? 2.5 : it.reference ? 1.25 : 1.75);
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
    label.textContent = `${this.cfg.popTrack ? this.trackLabel(this.cfg.popTrack) : 'No track'} · ${this.cfg.hap === 'both' ? 'H1+H2 mean' : this.cfg.hap}`;
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
      const color = this.itemColor(it);
      const sw = it.dash || it.reference ? h('span', { class: 'swatch line dashed', style: { color, borderTopColor: color } }) : h('span', { class: 'swatch line', style: { background: color } });
      const el = h('span', { class: `legend-item${this.hidden.has(it.key) ? ' off' : ''}`, title: 'Click to hide/show; hover to highlight' }, sw, h('span', null, it.label), h('span', { class: 'muted' }, it.sub || ''));
      el.addEventListener('click', () => { if (this.hidden.has(it.key)) this.hidden.delete(it.key); else this.hidden.add(it.key); this.render(); });
      el.addEventListener('mouseenter', () => { this.highlight = it.key; this.drawPanels(); });
      el.addEventListener('mouseleave', () => { this.highlight = null; this.drawPanels(); });
      this.legend.appendChild(el);
    }
    if (this.cfg.showObserved) this.legend.append(h('span', { class: 'legend-item', title: 'Observed ENCODE / FANTOM5 signal under each track (own scale)' }, h('span', { class: 'swatch', style: { background: css('--observed') } }), h('span', null, 'Observed data'), h('span', { class: 'muted' }, 'lane under each track')));
    if (this.cfg.mode === 'groups' && this.cfg.diff) this.legend.append(h('span', { class: 'muted' }, 'Each line: group mean minus the haplotype-weighted mean of the shown groups'));
    else if (this.showBand()) this.legend.append(h('span', { class: 'muted' }, 'Shaded: ±1 SD across haplotypes'));
    else if (this.cfg.envelope) this.legend.append(h('span', { class: 'muted' }, 'Shaded: min–max within each pixel bin'));
  }

  /** Figure export: the plot area (without the overview strip) with the view described as caption. */
  /** Region scalar form prefilled with the first data track and the visible range (reference offsets). */
  openScalar() {
    const [output, track] = parseTrackKey((this.cfg.tracks || [])[0] || '');
    const prefill = { gene: this.cfg.gene, output: output || undefined, track: Number.isInteger(track) ? track : undefined };
    if (this.cfg.coords !== 'haplotype') {
      prefill.start = Math.round(this.toGenomicOffset(this.viewport.start));
      prefill.end = Math.round(this.toGenomicOffset(this.viewport.end));
    }
    openScalarForm(prefill, { onSaved: (sc) => {
      if (sc.bin_field && this.cfg.mode === 'groups') toast(`Pick ${sc.bin_field} in Group by to split the cohort by it`);
    } });
  }

  figureOptions() {
    if (!this.info) return null;
    const c = this.cfg;
    const outputs = this.shownOutputs().map(outputLabel).join(', ') || 'no output';
    const view = c.mode === 'population'
      ? `Population heatmap of ${this.cfg.popTrack ? this.trackLabel(this.cfg.popTrack) : 'no track'}, rows sorted by ${c.groupField || 'sample'}`
      : c.mode === 'groups' ? `Group means by ${groupFieldLabel(c.groupField)}${c.diff ? ' (difference from the cohort mean)' : ''}` : `Pinned individuals (${c.hap === 'both' ? 'H1 and H2' : c.hap})`;
    const compare = [this.cfg.showRef && c.mode !== 'population' ? 'reference genome prediction (dashed)' : '', this.cfg.showObserved && c.mode !== 'population' ? 'observed ENCODE / FANTOM5 signal' : ''].filter(Boolean).join(' and ');
    const today = new Date().toISOString().slice(0, 10);
    return {
      root: this.scroll,
      redraw: () => this.render(),
      title: `${c.gene} · ${outputs} · ${this.locusInput.value}`,
      caption: [
        `${view}${compare ? `; with ${compare}` : ''}. Cohort: ${cohortDescription()}.`,
        `Dataset ${state.datasetId} · ${this.spanEl.textContent} · exported ${today} from genomics visualize`,
      ],
      filename: slug(`${state.datasetId}_${c.gene}_${this.locusInput.value.replace(/,/g, '')}_${c.mode}`),
    };
  }

  renderFooter() {
    const c = this.cfg;
    const coordText = { reference: 'genomic coordinates (haplotypes remapped through their indels; deletions show as gaps)', haplotype: 'raw haplotype prediction index', aligned: 'bcftools_chain training axis (insertion columns shared across the cohort)' }[c.coords];
    const res = this.data && this.cfg.mode !== 'population' ? ` · resolution ${fmtBp((this.data.end - this.data.start) / Math.max(1, this.data.edges.length - 1))}/bin` : '';
    const outputs = this.shownOutputs().map(outputLabel).join(', ') || 'no output';
    this.footer.textContent = `${c.gene} · ${outputs} · ${coordText}${res}`;
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
    if (this.seqHost.contains(e.target)) {
      showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div>${this.seqTooltip(Math.floor(pos))}`);
      return;
    }
    if (this.rulerHost.contains(e.target) && this.geneRows) {
      const rect = this.rulerHost.getBoundingClientRect();
      const gene = geneAt(this.geneRows, e.clientX - rect.left, e.clientY - rect.top, { top: 30, xOf: this.geneX });
      if (gene) { showTooltip(e.clientX, e.clientY, `${geneTooltip(gene, escapeHtml)}<div class="tt-sub">Click for the gene card (HGNC, Gene Ontology, links)</div>`); return; }
    }
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
    const panel = this.data && target ? this.data.panels[[...this.panelsHost.children].indexOf(target)] : null;
    if (!panel) { showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div>`); return; }
    if (e.target.classList && e.target.classList.contains('obs')) {
      const obs = this.data.observed ? this.data.observed.get(panel.trackKey) : null;
      if (!obs || !obs.available) { showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div><div class="tt-sub">${escapeHtml(obs ? obs.reason || 'No observed data' : 'Loading observed data…')}</div>`); return; }
      const oi = binIndex(obs.edges, pos);
      const ob = oi >= 0 ? obs.edges[oi + 1] - obs.edges[oi] : 0;
      showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div><div class="tt-sub" style="margin:-2px 0 6px">Observed · ${escapeHtml(obs.provider)} · ${escapeHtml(obs.units)}${ob > 1 ? ` · ${fmtInt(ob)} bp bin` : ''}</div>`
        + `<div class="tt-row"><span class="tt-key">mean</span><span class="tt-val">${oi >= 0 ? fmtNum(obs.mean[oi]) : '–'}</span></div>${ob > 1 ? `<div class="tt-row"><span class="tt-key">max in bin</span><span class="tt-val">${fmtNum(obs.max[oi])}</span></div>` : ''}`
        + `<div class="tt-sub">${escapeHtml((obs.items || []).map((x) => x.id).join(', '))} · click the lane name for links</div>`);
      return;
    }
    const bi = binIndex(panel.edges, pos);
    const rows = [];
    const items = this.displayItems(panel).filter((it) => !this.hidden.has(it.key));
    const sorted = items.map((it) => ({ it, v: bi >= 0 && it.mean ? it.mean[bi] : NaN })).sort((a, b) => (Number.isFinite(b.v) ? b.v : -Infinity) - (Number.isFinite(a.v) ? a.v : -Infinity));
    for (const { it, v } of sorted.slice(0, 12)) {
      const color = this.itemColor(it);
      const extra = it.upper && bi >= 0 ? ` <span class="muted">±${fmtNum((it.upper[bi] - it.lower[bi]) / 2, 2)}</span>` : '';
      rows.push(`<div class="tt-row"><span class="tt-key"><span class="swatch line" style="background:${color}"></span>${escapeHtml(it.label)}</span><span class="tt-val">${Number.isFinite(v) ? fmtNum(v) : '<span class="muted">gap</span>'}${extra}</span></div>`);
    }
    const binBp = bi >= 0 ? panel.edges[bi + 1] - panel.edges[bi] : 0;
    showTooltip(e.clientX, e.clientY, `<div class="tt-title">${escapeHtml(coordLabel)}</div><div class="tt-sub" style="margin:-2px 0 6px">${escapeHtml(this.trackLabel(panel.trackKey, true))}${binBp > 1 ? ` · mean of ${fmtInt(binBp)} bp` : ''}</div>${rows.join('')}`);
  }

  // ------------------------------------------------------------------ side panel
  renderSide() {
    const side = this.side;
    clear(side);
    // tracks
    if (this.cfg.mode === 'population') {
      const pick = h('select', { class: 'select', 'aria-label': 'Track', onchange: (e) => { this.cfg.popTrack = e.target.value; this.invalidate(); } },
        this.outputs().map((o) => h('optgroup', { label: outputLabel(o) }, this.trackMeta(o).map((m) => h('option', { value: trackKey(o, m.index) }, m.short || m.label)))));
      if (this.cfg.popTrack) pick.value = this.cfg.popTrack;
      side.appendChild(h('div', { class: 'side-section' },
        h('h3', null, 'Track'), pick,
        field('Max rows', select([100, 300, 600, 1200, 2400, 4000].map((n) => ({ value: n, label: `${fmtInt(n)} samples` })), this.cfg.popRows, (v) => { this.cfg.popRows = Number(v); this.invalidate(); }))));
    } else {
      side.appendChild(this.tracksSection());
    }
    if (this.cfg.mode !== 'population') side.appendChild(this.compareSection());
    side.appendChild(this.sequenceSection());
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

  /** Tracks of every output, one collapsible list per output, at most MAX_TRACKS in total. */
  tracksSection() {
    const chosen = new Set(this.cfg.tracks || []);
    this.openOutputs = this.openOutputs || new Set(this.tracksByOutput().map(([o]) => o));
    if (!this.openOutputs.size && this.outputs().length) this.openOutputs.add(this.outputs()[0]);
    const heading = h('h3', null);
    const counters = new Map();
    const boxes = [];
    const refreshCounts = () => {
      clear(heading).append(`Tracks (${chosen.size}/${MAX_TRACKS})`, chosen.size ? h('button', { class: 'btn small ghost', title: 'Unselect every track', onclick: () => setTracks([]) }, 'Clear') : '');
      for (const [o, el] of counters) {
        const n = [...chosen].filter((k) => parseTrackKey(k)[0] === o).length;
        el.textContent = `${n ? `${n} / ` : ''}${this.trackMeta(o).length}`;
        el.classList.toggle('on', n > 0);
      }
    };
    const setTracks = (keys) => {
      this.cfg.tracks = this.sortTracks(keys);
      chosen.clear(); this.cfg.tracks.forEach((k) => chosen.add(k));
      for (const { box, key } of boxes) box.checked = chosen.has(key);
      refreshCounts();
      this.invalidate();
    };
    const filter = h('input', { class: 'input', type: 'search', placeholder: 'Filter tracks (tissue, strand, mark…)', 'aria-label': 'Filter tracks' });
    const groups = this.outputs().map((o) => {
      const list = h('div', { class: 'checklist', style: { maxHeight: '220px' } });
      for (const m of this.trackMeta(o)) {
        const key = trackKey(o, m.index);
        const box = h('input', { type: 'checkbox' });
        box.checked = chosen.has(key);
        box.addEventListener('change', () => {
          if (box.checked && chosen.size >= MAX_TRACKS) { box.checked = false; toast(`At most ${MAX_TRACKS} tracks at once`, 'error'); return; }
          setTracks(box.checked ? [...chosen, key] : [...chosen].filter((k) => k !== key));
        });
        boxes.push({ box, key });
        const text = `${m.label} ${Object.values(m.metadata || {}).join(' ')}`.toLowerCase();
        list.appendChild(h('label', { title: m.label, 'data-text': text }, box, h('span', null, m.short || m.label), h('span', { class: 'meta' }, `#${m.index}`)));
      }
      const count = h('span', { class: 'track-count' });
      counters.set(o, count);
      const addFirst = () => {
        const shown = [...list.querySelectorAll('label')].filter((l) => !l.hidden).map((l) => l.querySelector('input'));
        const keys = boxes.filter((b) => shown.includes(b.box) && !chosen.has(b.key)).map((b) => b.key);
        const room = MAX_TRACKS - chosen.size;
        if (!room) { toast(`At most ${MAX_TRACKS} tracks at once`, 'error'); return; }
        setTracks([...chosen, ...keys.slice(0, room)]);
      };
      const details = h('details', { class: 'track-group', open: this.openOutputs.has(o) },
        h('summary', null, h('span', null, outputLabel(o)), count),
        h('div', { class: 'track-group-actions' },
          h('button', { class: 'btn small ghost', title: `Add the shown ${outputLabel(o)} tracks (up to ${MAX_TRACKS} in total)`, onclick: addFirst }, 'Add shown'),
          h('button', { class: 'btn small ghost', onclick: () => setTracks([...chosen].filter((k) => parseTrackKey(k)[0] !== o)) }, 'None')),
        list);
      details.addEventListener('toggle', () => { if (details.open) this.openOutputs.add(o); else this.openOutputs.delete(o); });
      return { details, list };
    });
    filter.addEventListener('input', () => {
      const q = filter.value.trim().toLowerCase();
      for (const { details, list } of groups) {
        let hits = 0;
        for (const label of list.querySelectorAll('label')) { label.hidden = !!q && !label.dataset.text.includes(q); if (!label.hidden) hits++; }
        details.hidden = !!q && !hits;
        if (q && hits) details.open = true;
      }
    });
    refreshCounts();
    return h('div', { class: 'side-section' }, heading,
      this.outputs().length > 1 ? h('p', { class: 'help', style: { margin: 0 } }, 'Mix tracks of any outputs; each track gets its own panel.') : null,
      filter, groups.map((g) => g.details));
  }

  /** Reference-genome prediction and observed (ENCODE / FANTOM5) data next to the predictions. */
  compareSection() {
    const shown = this.tracksByOutput().map(([o]) => o);
    const missingRef = shown.filter((o) => !this.hasReference(o));
    const predictRef = () => {
      const terms = [...new Set((this.cfg.tracks || []).map((k) => (this.track(k).meta.metadata || {}).ontology_curie).filter(Boolean))];
      openPredictForm(state.datasetId, { haplotypes: ['ref'], genes: [this.cfg.gene], outputs: (missingRef.length ? missingRef : shown).map((o) => o.toUpperCase()), ontology_terms: terms.length ? terms : undefined });
    };
    const refHelp = !shown.length ? null : missingRef.length
      ? h('p', { class: 'help' }, `No AlphaGenome prediction of the reference window for ${missingRef.map(outputLabel).join(', ')} in ${this.cfg.gene}. `,
        h('a', { href: '#', onclick: (e) => { e.preventDefault(); predictRef(); } }, 'Predict the reference…'))
      : h('p', { class: 'help' }, 'Dashed dark line: AlphaGenome on the reference genome (hg38), the baseline every haplotype deviates from.');
    const scale = shown.includes('cage') && this.cfg.showObserved
      ? field('FANTOM5 CAGE values', segmented([{ value: 'tpm', label: 'TPM', title: 'Tags per million (comparable across libraries)' }, { value: 'counts', label: 'Read counts', title: 'CTSS read counts' }], this.cfg.obsScale, (v) => { this.cfg.obsScale = v; this.invalidate(); }))
      : null;
    return h('div', { class: 'side-section' },
      h('h3', null, 'Compare with'),
      checkbox('Reference genome (AlphaGenome)', this.cfg.showRef, (v) => { this.cfg.showRef = v; this.invalidate(); }),
      refHelp,
      checkbox('Observed data (ENCODE / FANTOM5)', this.cfg.showObserved, (v) => { this.cfg.showObserved = v; this.renderSide(); this.invalidate(); }),
      scale,
      h('p', { class: 'help' }, 'The experiments with the same ontology term, assay and target as each track (CAGE: FANTOM5 libraries; RNA-seq, DNase, ATAC, ChIP: ENCODE GRCh38 bigWigs), read over the window on first use and cached. Units are the source\'s (TPM, signal, fold change), not AlphaGenome\'s: compare shapes and peaks. Click a track name for its ontology terms and experiments.'));
  }

  sequenceSection() {
    return h('div', { class: 'side-section' },
      h('h3', null, 'Sequence lane'),
      segmented([
        { value: 'pinned', label: 'Pinned', title: 'Letter frequency per base over the pinned haplotypes (reference when nothing is pinned)' },
        { value: 'reference', label: 'Reference', title: 'The reference genome of the window' },
        { value: 'off', label: 'Hidden' },
      ], this.cfg.seqSource, (v) => { this.cfg.seqSource = v; this.seq = null; this.persist(); this.render(); this.requestSequenceNow(); }),
      h('p', { class: 'help' }, 'Zoomed in, each base shows its letters scaled by their frequency; zoomed out, stacked bars show base composition per bin (gaps are deleted bases). Haplotypes follow the Haplotype setting.'));
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
      field(this.cfg.mode === 'population' ? 'Sort rows by' : 'Group by', select([...cats.map((f) => ({ value: f.name, label: f.label })), ...(isVariantField(this.cfg.groupField) ? [{ value: this.cfg.groupField, label: groupFieldLabel(this.cfg.groupField) }] : [])], this.cfg.groupField, (v) => { this.cfg.groupField = v; this.cfg.groups = null; this.renderSide(); this.invalidate(); })),
      this.cfg.mode === 'groups' && isVariantField(this.cfg.groupField)
        ? h('p', { class: 'help' }, 'Genotype classes at this site (REF/REF, REF/ALT, ALT/ALT) over the cohort. ', h('a', { href: `#/variant?${new URLSearchParams({ gene: this.cfg.gene, pos: this.cfg.groupField.split(':')[1], ref: this.cfg.groupField.split(':')[2], alt: this.cfg.groupField.split(':')[3] })}` }, 'Variant page'))
        : this.cfg.mode === 'groups' ? h('div', null, h('div', { class: 'label', style: { marginBottom: '4px' } }, 'Groups (max 8)'), list) : null,
      h('p', { class: 'help' }, `Cohort: ${fmtInt(cohort)} samples${filterText ? ` (${filterText})` : ' (all)'}. `, h('a', { href: '#/samples' }, 'Edit cohort')),
      this.cfg.mode === 'groups' ? h('p', { class: 'help' }, 'Means are computed once over every haplotype in each group (progress is shown) and cached on disk, so later views are instant.') : null,
      this.cfg.coords === 'aligned' ? h('p', { class: 'help', style: { color: 'var(--ink-2)' } }, 'Training-axis aggregates need a bcftools alignment entry per haplotype; uncached entries are built on first use, which is slow for large cohorts. Genomic coordinates give the same per-base alignment without that cost.') : null);
  }

  // ------------------------------------------------------------------ lifecycle
  onEvent(topic) {
    if (topic === 'pinned') {
      this.renderSide();
      if (this.cfg.mode === 'individuals') this.invalidate();
      else if (this.cfg.seqSource === 'pinned') this.requestSequenceNow();
    }
    if (topic === 'cohort' && this.cfg.mode !== 'individuals') { this.renderSide(); this.invalidate(); }
    if (topic === 'samples') {
      const cats = categoricalFields();
      if (!isVariantField(this.cfg.groupField) && !cats.some((f) => f.name === this.cfg.groupField)) {
        this.cfg.groupField = defaultGroupField(cats); this.cfg.groups = null;
        if (this.cfg.mode !== 'individuals') this.invalidate();
      }
      this.renderSide();
    }
  }

  update(params) {
    if (params.gene && params.gene !== this.cfg.gene) this.setGene(params.gene, {});
  }

  unmount() {
    this.latest.abort();
    this.request.cancel();
    this.seqLatest.abort();
    this.requestSequence.cancel();
    this.persistDebounced.cancel();
    hideTooltip();
    for (const d of this.disposers) d();
  }

  async loadAnnotations() {
    const gene = this.cfg.gene;
    try {
      const res = await loadGeneAnnotations(state.datasetId, gene, { onProgress: () => this.status('Loading gene annotations…') });
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
