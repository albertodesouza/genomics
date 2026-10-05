// Sequence: pinned haplotypes against the reference (letters when zoomed in, mismatch density
// when zoomed out), with gene models and the variants carried by those samples.
import { api, apiJob, isAbort, Latest } from '../api.js';
import { navigate, updateRouteParams } from '../app.js';
import { state, ds, geneInfo, setLocus, togglePinned } from '../state.js';
import { h, clear, icon, iconButton, segmented, select, field, fmtInt, fmtBp, fmtPct, showTooltip, hideTooltip, toast, jobOverlay, debounce, escapeHtml, errorBox } from '../ui.js';
import { Viewport, setupCanvas, theme, ticks, withAlpha, onResize, css } from '../plot.js';
import { loadGeneAnnotations, packGenes, drawGeneLanes, geneLanesHeight, geneAt, geneTooltip } from '../gene_lanes.js';
import { labelGeneOptions, openGeneCard } from '../cards.js';
import { exportMenu, slug } from '../figure.js';

const GUTTER_L = 120;
const GUTTER_R = 14;
const ROW_H = 18;
const MAX_LETTER_SPAN = 4000;
const MAX_SAMPLES = 24;

export async function mount(root, params) {
  const page = new SequencePage(root);
  await page.init(params);
  return page;
}

class SequencePage {
  constructor(root) {
    this.root = root;
    this.latest = new Latest();
    this.data = null;
    this.axis = null;
    this.annotations = null;
    this.geneRows = [];
    this.genesH = 0;
    this.disposers = [];
  }

  async init(params) {
    const genes = state.summary.genes.map((g) => g.gene);
    const saved = state.locus || {};
    this.cfg = {
      gene: params.gene && genes.includes(params.gene) ? params.gene : (saved.gene && genes.includes(saved.gene) ? saved.gene : genes[0]),
      // A URL locus without coords is genomic (coords are omitted from URLs when genomic).
      coords: params.coords || (params.start !== undefined ? 'reference' : (saved.coords === 'haplotype' ? 'reference' : saved.coords)) || 'reference',
      hap: saved.seqHap || 'H1+H2',
      matches: saved.seqMatches || 'dots',
    };
    this.build();
    this.viewport = new Viewport({ length: 1, minSpan: 20, onChange: () => this.onViewChange() });
    this.viewport.attach(this.host, () => this.geom(), {
      onHover: (e, frac) => this.onHover(e, frac),
      onLeave: () => { hideTooltip(); },
    });
    this.host.addEventListener('click', (e) => this.onClick(e));
    this.disposers.push(onResize(this.scroll, () => this.render()));
    const onTheme = () => this.render();
    window.addEventListener('themechange', onTheme);
    this.disposers.push(() => window.removeEventListener('themechange', onTheme));
    await this.setGene(this.cfg.gene, params);
  }

  build() {
    const genes = state.summary.genes.map((g) => g.gene);
    this.geneSelect = select(genes, this.cfg.gene, (v) => this.setGene(v, {}), { 'aria-label': 'Gene', style: { maxWidth: '300px' } });
    labelGeneOptions(this.geneSelect);
    this.coordSeg = segmented([
      { value: 'reference', label: 'Genomic', title: 'Reference coordinates; insertions drawn as markers between bases' },
      { value: 'aligned', label: 'Training axis', title: 'bcftools_chain expanded alignment used for CNN inputs (insertions become columns)' },
      { value: 'haplotype', label: 'Haplotype', title: 'Raw haplotype sequence positions (not aligned after indels)' },
    ], this.cfg.coords, (v) => this.setCoords(v));
    this.hapSeg = segmented([{ value: 'H1', label: 'H1' }, { value: 'H2', label: 'H2' }, { value: 'H1+H2', label: 'Both' }], this.cfg.hap, (v) => { this.cfg.hap = v; this.persist(); this.fetch(true); });
    this.matchSeg = segmented([{ value: 'dots', label: 'Matches as ·' }, { value: 'bases', label: 'All bases' }], this.cfg.matches, (v) => { this.cfg.matches = v; this.persist(); this.render(); });
    this.locusInput = h('input', { class: 'input locus-input', 'aria-label': 'Locus', spellcheck: 'false' });
    this.locusInput.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.gotoLocus(this.locusInput.value); });
    this.spanEl = h('span', { class: 'locus-span' });
    const toolbar = h('div', { class: 'ws-toolbar' },
      field('Gene', h('div', { style: { display: 'flex', gap: '4px', alignItems: 'center' } }, this.geneSelect, iconButton('info', 'Gene card: HGNC names, database links and Gene Ontology', () => openGeneCard(this.cfg.gene), 'icon-btn bordered'))), field('Coordinates', this.coordSeg), h('div', { class: 'divider' }),
      field('Haplotypes', this.hapSeg), field('Bases', this.matchSeg), h('div', { style: { flex: '1' } }),
      h('a', { class: 'btn small', href: '#/samples' }, icon('pin', 13), `${state.pinned.length} pinned`));
    const locus = h('div', { class: 'locus' },
      iconButton('left', 'Pan left', () => this.viewport.pan(-this.viewport.span * 0.5), 'icon-btn bordered'),
      iconButton('right', 'Pan right', () => this.viewport.pan(this.viewport.span * 0.5), 'icon-btn bordered'),
      iconButton('minus', 'Zoom out', () => this.viewport.zoom(2), 'icon-btn bordered'),
      iconButton('plus', 'Zoom in', () => this.viewport.zoom(0.5), 'icon-btn bordered'),
      this.locusInput, this.spanEl,
      h('button', { class: 'btn small', onclick: () => { const mw = this.modelWindow(); if (mw) this.viewport.set(mw.start, mw.end); } }, icon('target', 14), 'Model window'),
      h('button', { class: 'btn small', onclick: () => this.viewport.set(0, this.viewport.length) }, icon('expand', 14), 'Whole window'),
      h('button', { class: 'btn small ghost', onclick: () => navigate('tracks', { gene: this.cfg.gene }) }, 'Tracks →'),
      exportMenu(() => this.figureOptions()),
      h('span', { class: 'hint' }, 'Zoom below 4 kb for bases · click a variant tick to center it'));
    this.host = h('div', { class: 'seq-host canvas-host' }, h('canvas'));
    this.loadingLine = h('div', { class: 'loading-line', hidden: true });
    this.table = h('div', { style: { padding: '0 16px 24px' } });
    this.scroll = h('div', { class: 'ws-scroll' }, this.loadingLine, this.host, this.table);
    this.root.appendChild(h('div', { class: 'workspace side-collapsed' }, h('div', { class: 'ws-main' }, toolbar, locus, this.scroll), h('aside', { class: 'side' })));
  }

  geom() {
    const w = this.scroll.clientWidth;
    return { left: GUTTER_L, width: Math.max(50, w - GUTTER_L - GUTTER_R), total: w };
  }

  async setGene(gene, params = {}) {
    this.cfg.gene = gene;
    this.geneSelect.value = gene;
    this.data = null;
    this.axis = null;
    this.annotations = null;
    try { this.info = await geneInfo(gene); } catch (err) { clear(this.host).appendChild(errorBox(err)); return; }
    loadGeneAnnotations(state.datasetId, gene).then((res) => { if (gene === this.cfg.gene) { this.annotations = res; this.render(); } }).catch(() => {});
    if (this.cfg.coords === 'aligned' && !(await this.loadAxis())) { this.cfg.coords = 'reference'; this.coordSeg.setValue('reference'); }
    this.viewport.setDomain(this.domain());
    const saved = state.locus || {};
    if (params.start !== undefined && params.end !== undefined) this.viewport.set(Number(params.start), Number(params.end), false);
    else if (saved.gene === gene && saved.coords === this.cfg.coords && saved.end > saved.start) {
      // Arriving from Tracks: keep its locus but zoom to a readable span around its centre.
      const c = (saved.start + saved.end) / 2;
      const span = Math.min(saved.end - saved.start, 600);
      this.viewport.set(c - span / 2, c + span / 2, false);
    } else {
      const mw = this.modelWindow();
      const c = mw ? (mw.start + mw.end) / 2 : this.viewport.length / 2;
      this.viewport.set(c - 300, c + 300, false);
    }
    this.onViewChange();
  }

  domain() {
    if (this.cfg.coords === 'aligned' && this.axis) return this.axis.expanded_length;
    return this.info.length || 1;
  }

  modelWindow() { return this.cfg.coords === 'aligned' ? (this.axis && this.axis.model_window) : this.info.model_window; }

  async loadAxis() {
    const overlay = jobOverlay(this.scroll, 'Preparing the training alignment axis');
    try {
      this.axis = await apiJob(`${ds()}/genes/${encodeURIComponent(this.cfg.gene)}/axis`, { onProgress: (j) => overlay.update(j) });
      return true;
    } catch (err) {
      if (!isAbort(err)) toast(`Training axis unavailable: ${err.message}`, 'error', 8000);
      return false;
    } finally { overlay.remove(); }
  }

  fromGenomic(ref) {
    if (this.cfg.coords !== 'aligned' || !this.axis) return ref;
    const rr = ref - this.axis.ref_start_offset;
    if (rr < 0) return 0;
    const slots = this.axis.insertion_slots;
    let k = 0;
    while (k < slots.length && slots[k] <= rr + k) k++;
    return rr + k;
  }

  toGenomic(pos) {
    if (this.cfg.coords !== 'aligned' || !this.axis) return pos;
    const slots = this.axis.insertion_slots;
    let lo = 0; let hi = slots.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (slots[m] < pos) lo = m + 1; else hi = m; }
    return this.axis.ref_start_offset + pos - lo;
  }

  async setCoords(coords) {
    const center = this.toGenomic(this.viewport.start + this.viewport.span / 2);
    const span = this.viewport.span;
    const prev = this.cfg.coords;
    this.cfg.coords = coords;
    if (coords === 'aligned' && !this.axis && !(await this.loadAxis())) { this.cfg.coords = prev; this.coordSeg.setValue(prev); return; }
    this.viewport.setDomain(this.domain());
    this.data = null;
    const c = this.fromGenomic(center);
    this.viewport.set(c - span / 2, c + span / 2);
  }

  posLabel(pos) {
    if (this.cfg.coords === 'haplotype') return `hap ${fmtInt(pos + 1)}`;
    if (!this.info.start) return fmtInt(pos + 1);
    return `${this.info.chromosome}:${fmtInt(this.info.start + this.toGenomic(pos))}`;
  }

  renderLocus() {
    const v = this.viewport;
    if (this.cfg.coords === 'haplotype') this.locusInput.value = `hap:${fmtInt(v.start + 1)}-${fmtInt(v.end)}`;
    else if (this.cfg.coords === 'aligned') this.locusInput.value = `axis:${fmtInt(v.start + 1)}-${fmtInt(v.end)}`;
    else this.locusInput.value = `${this.info.chromosome}:${fmtInt(this.info.start + v.start)}-${fmtInt(this.info.start + v.end - 1)}`;
    this.spanEl.textContent = `${fmtBp(v.span)} of ${fmtBp(v.length)}`;
  }

  gotoLocus(text) {
    const m = /^(?:([\w.]+):)?\s*(\d+)(?:\s*[-–]\s*(\d+))?$/.exec(text.replace(/,/g, '').trim());
    if (!m) { toast('Use chr:start-end, start-end or a position', 'error'); return; }
    const genomic = this.info.start && !['hap', 'axis'].includes(m[1]) && Number(m[2]) >= this.info.start;
    const conv = (x) => (genomic ? this.fromGenomic(x - this.info.start) : x - 1);
    if (!m[3]) { this.viewport.center(conv(Number(m[2])), Math.min(this.viewport.span, 400)); return; }
    this.viewport.set(conv(Number(m[2])), conv(Number(m[3])) + 1);
  }

  rows() {
    const haps = this.cfg.hap === 'H1+H2' ? 'H1+H2' : this.cfg.hap;
    return state.pinned.slice(0, MAX_SAMPLES).map((s) => `${s}:${haps}`).join(',');
  }

  persist() {
    setLocus({ gene: this.cfg.gene, coords: this.cfg.coords === 'haplotype' ? (state.locus.coords || 'reference') : this.cfg.coords, seqHap: this.cfg.hap, seqMatches: this.cfg.matches, start: this.viewport.start, end: this.viewport.end });
    updateRouteParams({ gene: this.cfg.gene, coords: this.cfg.coords !== 'reference' ? this.cfg.coords : null, start: this.viewport.start, end: this.viewport.end });
  }

  onViewChange() {
    this.renderLocus();
    this.render();
    this.fetchDebounced();
  }

  fetchDebounced = debounce(() => { this.persist(); this.fetch(false); }, 140);

  async fetch(force) {
    const v = this.viewport;
    const g = this.geom();
    const letters = v.span <= MAX_LETTER_SPAN;
    const d = this.data;
    const key = JSON.stringify([this.cfg.gene, this.cfg.coords, this.rows()]);
    if (!force && d && d.key === key && v.start >= d.start && v.end <= d.end) {
      const dataLetters = d.mode === 'letters';
      const bp = d.edges ? (d.end - d.start) / (d.edges.length - 1) : 1;
      if (dataLetters === letters && (letters || bp <= Math.max(1, v.span / g.width) * 1.6)) { this.render(); return; }
    }
    if (!state.pinned.length) { this.data = null; this.render(); return; }
    let start; let end;
    if (letters) {
      const margin = Math.min(v.span * 0.5, (MAX_LETTER_SPAN - v.span) / 2);
      start = Math.max(0, Math.floor(v.start - margin)); end = Math.min(v.length, Math.ceil(v.end + margin));
    } else {
      start = Math.max(0, Math.floor(v.start - v.span * 0.5)); end = Math.min(v.length, Math.ceil(v.end + v.span * 0.5));
    }
    const bins = Math.max(32, Math.min(4096, Math.round(g.width * (end - start) / v.span)));
    const signal = this.latest.next();
    this.loadingLine.hidden = false;
    try {
      const res = await api(`${ds()}/sequence`, { params: { gene: this.cfg.gene, coords: this.cfg.coords, rows: this.rows(), start, end, bins }, signal });
      res.key = key;
      this.data = res;
      const failed = res.rows.filter((r) => r.error);
      if (failed.length) toast(`${failed[0].sample} ${failed[0].haplotype}: ${failed[0].error}`, 'error', 7000);
    } catch (err) {
      if (isAbort(err)) return;
      toast(err.message, 'error');
    } finally {
      if (!signal.aborted) this.loadingLine.hidden = true;
    }
    this.render();
  }

  xOf(pos) { const g = this.geom(); const v = this.viewport; return g.left + ((pos - v.start) / v.span) * g.width; }

  /** Gene x position (reference offset -> canvas x in the current coordinates). */
  geneX = (p) => this.xOf(this.fromGenomic(p));

  /** Figure export: the sequence canvas (not the genotype table), the view described as caption. */
  figureOptions() {
    if (!this.info) return null;
    const shown = state.pinned.slice(0, MAX_SAMPLES);
    const today = new Date().toISOString().slice(0, 10);
    return {
      root: this.host,
      redraw: () => this.render(),
      title: `${this.cfg.gene} · haplotypes against the reference · ${this.locusInput.value}`,
      caption: [
        `${shown.length} pinned individual${shown.length === 1 ? '' : 's'}${shown.length ? `: ${shown.join(', ')}` : ''}.`,
        `Dataset ${state.datasetId} · ${this.spanEl.textContent} · exported ${today} from genomics visualize`,
      ],
      filename: slug(`${state.datasetId}_${this.cfg.gene}_${this.locusInput.value.replace(/,/g, '')}_sequence`),
    };
  }

  render() {
    if (!this.info) return;
    const g = this.geom();
    const canvas = this.host.querySelector('canvas');
    const t = theme();
    const d = this.data;
    const rows = d ? d.rows : [];
    const v = this.viewport;
    this.geneRows = packGenes((this.annotations && this.annotations.genes) || [], { start: v.start, end: v.end, pxPerBp: g.width / v.span, toPos: (p) => this.fromGenomic(p) });
    this.genesH = geneLanesHeight(this.geneRows);
    const top = 52 + this.genesH; // ruler (30) + gene lanes + variant lane (22)
    const H = top + ROW_H + 4 + Math.max(1, rows.length) * ROW_H + 16;
    const ctx = setupCanvas(canvas, g.total, H);
    this.layout = { top, rows: rows.length };
    this.drawRuler(ctx, t, g);
    drawGeneLanes(ctx, t, this.geneRows, { top: 30, left: g.left, width: g.width, xOf: this.geneX, highlight: this.cfg.gene });
    if (!d) {
      ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.textAlign = 'center';
      ctx.fillText(state.pinned.length ? 'Loading…' : 'Pin individuals on the Samples page to compare their haplotypes', g.left + g.width / 2, top + 40);
      this.renderTable();
      return;
    }
    const mw = this.modelWindow();
    if (mw && mw.end > v.start && mw.start < v.end) {
      const a = Math.max(g.left, this.xOf(mw.start)); const b = Math.min(g.left + g.width, this.xOf(mw.end));
      ctx.fillStyle = t.modelWindow; ctx.fillRect(a, top, b - a, H - top);
    }
    this.drawVariants(ctx, t, g, d);
    ctx.font = `11px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillStyle = t.ink2; ctx.fillText('Reference', g.left - 8, top + ROW_H / 2);
    rows.forEach((row, i) => {
      const y = top + ROW_H + 4 + i * ROW_H;
      ctx.fillStyle = row.error ? t.ink3 : t.ink2;
      ctx.fillText(`${row.sample} ${row.haplotype}`, g.left - 8, y + ROW_H / 2);
    });
    ctx.save();
    ctx.beginPath(); ctx.rect(g.left, top, g.width, H - top); ctx.clip();
    if (d.mode === 'letters') this.drawLetters(ctx, t, g, d, top);
    else this.drawDensity(ctx, t, g, d, top);
    ctx.restore();
    this.renderTable();
  }

  drawRuler(ctx, t, g) {
    const v = this.viewport;
    const offset = this.cfg.coords === 'reference' && this.info.start ? this.info.start : 1;
    ctx.font = `11px ${t.font}`; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
    ctx.strokeStyle = t.lineStrong; ctx.lineWidth = 1;
    for (const tv of ticks(v.start + offset, v.end + offset, Math.max(3, Math.floor(g.width / 120)))) {
      const pos = tv - offset;
      const x = Math.round(this.xOf(pos + 0.5)) + 0.5;
      if (x < g.left || x > g.left + g.width) continue;
      ctx.beginPath(); ctx.moveTo(x, 22); ctx.lineTo(x, 29); ctx.stroke();
      ctx.fillStyle = t.ink2;
      ctx.fillText(this.cfg.coords === 'aligned' ? fmtInt(this.info.start + this.toGenomic(pos)) : fmtInt(tv), x, 18);
    }
    ctx.beginPath(); ctx.moveTo(g.left, 29.5); ctx.lineTo(g.left + g.width, 29.5); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillStyle = t.ink3; ctx.font = `10.5px ${t.font}`;
    ctx.fillText(this.cfg.coords === 'haplotype' ? 'hap pos' : this.info.chromosome || '', g.left - 8, 18);
    ctx.fillText('variants', g.left - 8, 44 + this.genesH);
  }

  variantColor(t, type) {
    return { SNV: t.ink2, MNV: t.ink2, INS: css('--series-7'), DEL: css('--series-8') }[type] || css('--series-4');
  }

  drawVariants(ctx, t, g, d) {
    const v = this.viewport;
    const list = (d.variants || []).filter((x) => x.pos >= v.start - 1 && x.pos < v.end + 1);
    this.visibleVariants = list;
    const pxPerBase = g.width / v.span;
    for (const x of list) {
      const cx = this.xOf(x.pos + 0.5);
      const w = Math.max(2, Math.min(pxPerBase - 1, 8));
      ctx.fillStyle = this.variantColor(t, x.type);
      const hgt = 6 + Math.min(10, x.carriers.length * 2);
      ctx.fillRect(cx - w / 2, 48 + this.genesH - hgt, w, hgt);
    }
  }

  drawLetters(ctx, t, g, d, top) {
    const v = this.viewport;
    const px = g.width / v.span;
    const ref = d.reference;
    const from = Math.max(v.start, d.start); const to = Math.min(v.end, d.end);
    const showLetters = px >= 7.5;
    ctx.font = `${Math.min(13, Math.max(9, px * 0.75))}px ${t.mono}`;
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    const baseColor = (b) => t.bases[b] || t.bases.N;
    const cell = (x, y, color, alpha) => { ctx.fillStyle = alpha < 1 ? withAlpha(color, alpha) : color; ctx.fillRect(x + (px > 4 ? 0.5 : 0), y + 1, Math.max(1, px - (px > 4 ? 1 : 0)), ROW_H - 2); };
    // reference row
    for (let p = from; p < to; p++) {
      const b = ref[p - d.start];
      const x = this.xOf(p);
      if (b === '-') { ctx.fillStyle = t.surface2; ctx.fillRect(x, top + 1, px, ROW_H - 2); continue; }
      cell(x, top, baseColor(b), 0.22);
      if (showLetters) { ctx.fillStyle = t.ink; ctx.fillText(b, x + px / 2, top + ROW_H / 2 + 0.5); }
    }
    d.rows.forEach((row, i) => {
      const y = top + ROW_H + 4 + i * ROW_H;
      if (!row.bases) return;
      for (let p = from; p < to; p++) {
        const k = p - d.start;
        const b = row.bases[k]; const r = ref[k];
        const x = this.xOf(p);
        if (b === '-' && r !== '-') {
          ctx.fillStyle = t.ink2; ctx.fillRect(x, y + ROW_H / 2 - 1, px, 2);
        } else if (b !== r && b !== '-' && b !== 'N') {
          cell(x, y, baseColor(b), 1);
          if (showLetters) { ctx.fillStyle = '#fff'; ctx.fillText(b, x + px / 2, y + ROW_H / 2 + 0.5); }
        } else if (b === r && b !== '-') {
          if (showLetters) {
            ctx.fillStyle = this.cfg.matches === 'dots' ? t.ink3 : t.ink2;
            ctx.fillText(this.cfg.matches === 'dots' ? '·' : b, x + px / 2, y + ROW_H / 2 + 0.5);
          } else { ctx.fillStyle = t.surface2; ctx.fillRect(x, y + 3, px, ROW_H - 6); }
        }
      }
      for (const ins of row.insertions || []) {
        if (ins.pos < v.start - 1 || ins.pos >= v.end) continue;
        const x = this.xOf(ins.pos + 1);
        ctx.fillStyle = css('--series-7');
        ctx.fillRect(x - 1, y + 1, 2, ROW_H - 2);
        ctx.fillRect(x - 3, y + 1, 6, 2);
        ctx.fillRect(x - 3, y + ROW_H - 3, 6, 2);
      }
    });
  }

  drawDensity(ctx, t, g, d, top) {
    const e = d.edges;
    const xs = Array.from(e, (p) => this.xOf(p));
    // reference lane: GC content
    for (let i = 0; i < e.length - 1; i++) {
      const gc = d.gc[i];
      ctx.fillStyle = withAlpha(t.accent, 0.1 + 0.6 * Math.max(0, Math.min(1, (gc - 0.3) / 0.4)));
      ctx.fillRect(xs[i], top + 2, Math.max(1, xs[i + 1] - xs[i]), ROW_H - 4);
    }
    const mismatch = css('--series-2');
    const ins = css('--series-7');
    d.rows.forEach((row, r) => {
      const y = top + ROW_H + 4 + r * ROW_H;
      ctx.fillStyle = t.surface2; ctx.fillRect(g.left, y + 3, g.width, ROW_H - 6);
      if (!row.density) return;
      const { mismatch: mm, deletion: del, insertion: inn } = row.density;
      for (let i = 0; i < e.length - 1; i++) {
        const w = Math.max(1, xs[i + 1] - xs[i]);
        if (del[i] > 0) { ctx.fillStyle = withAlpha(t.ink2, Math.min(1, 0.35 + del[i] * 4)); ctx.fillRect(xs[i], y + ROW_H / 2 - 1, w, 2); }
        if (mm[i] > 0) { ctx.fillStyle = withAlpha(mismatch, Math.min(1, 0.45 + mm[i] * 20)); ctx.fillRect(xs[i], y + 2, Math.max(1.5, w), ROW_H - 4); }
        if (inn[i] > 0) { ctx.fillStyle = ins; ctx.fillRect(xs[i], y + 1, Math.max(1.5, w), 3); }
      }
    });
    ctx.fillStyle = t.ink3; ctx.font = `11px ${t.font}`; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText('Zoomed out: orange = mismatches, line = deletions, violet = insertions; reference lane shaded by GC content', g.left + 4, top + ROW_H + 8 + d.rows.length * ROW_H);
  }

  hitTest(e) {
    const rect = this.host.getBoundingClientRect();
    const g = this.geom();
    const x = e.clientX - rect.left; const y = e.clientY - rect.top;
    const pos = Math.floor(this.viewport.start + ((x - g.left) / g.width) * this.viewport.span);
    let row = null;
    if (this.layout && y >= this.layout.top + ROW_H + 4) row = Math.floor((y - this.layout.top - ROW_H - 4) / ROW_H);
    const gh = this.genesH;
    const lane = y >= 30 && y < 30 + gh ? 'genes' : y >= 30 + gh && y < 50 + gh ? 'variants' : y >= this.layout?.top && y < this.layout.top + ROW_H ? 'reference' : 'rows';
    return { x, y, pos, row, lane };
  }

  variantAt(pos, tolerancePx = 4) {
    if (!this.visibleVariants) return null;
    const bpTol = (tolerancePx / this.geom().width) * this.viewport.span;
    let best = null; let bestD = Infinity;
    for (const v of this.visibleVariants) { const dd = Math.abs(v.pos + 0.5 - pos - 0.5); if (dd <= Math.max(0.5, bpTol) && dd < bestD) { best = v; bestD = dd; } }
    return best;
  }

  onHover(e, frac) {
    if (frac < 0 || frac > 1) { hideTooltip(); return; }
    const hit = this.hitTest(e);
    if (hit.lane === 'genes') {
      const gene = geneAt(this.geneRows, hit.x, hit.y, { top: 30, xOf: this.geneX });
      if (gene) showTooltip(e.clientX, e.clientY, geneTooltip(gene, escapeHtml)); else hideTooltip();
      return;
    }
    if (!this.data) { hideTooltip(); return; }
    const d = this.data;
    if (hit.lane === 'variants') {
      const v = this.variantAt(hit.pos);
      if (!v) { hideTooltip(); return; }
      showTooltip(e.clientX, e.clientY, variantHtml(v));
      return;
    }
    const k = hit.pos - d.start;
    const lines = [`<div class="tt-title">${escapeHtml(this.posLabel(hit.pos))}</div>`];
    if (d.mode === 'letters') {
      if (k >= 0 && k < d.reference.length) lines.push(`<div class="tt-row"><span class="tt-key">reference</span><span class="tt-val mono">${escapeHtml(d.reference[k])}</span></div>`);
      const row = hit.row !== null ? d.rows[hit.row] : null;
      if (row && row.bases && k >= 0 && k < row.bases.length) {
        lines.push(`<div class="tt-row"><span class="tt-key">${escapeHtml(row.sample)} ${row.haplotype}</span><span class="tt-val mono">${escapeHtml(row.bases[k] === '-' ? 'deleted' : row.bases[k])}</span></div>`);
        const ins = (row.insertions || []).find((i) => i.pos === hit.pos || i.pos === hit.pos - 1);
        if (ins) lines.push(`<div class="tt-sub">insertion after ${escapeHtml(this.posLabel(ins.pos))}: +${ins.length} bp ${escapeHtml(ins.seq.slice(0, 40))}${ins.seq.length > 40 ? '…' : ''}</div>`);
      }
    } else if (hit.row !== null && d.rows[hit.row] && d.rows[hit.row].density) {
      const row = d.rows[hit.row];
      let i = 0; while (i < d.edges.length - 2 && d.edges[i + 1] <= hit.pos) i++;
      lines.push(`<div class="tt-row"><span class="tt-key">${escapeHtml(row.sample)} ${row.haplotype}</span><span class="tt-val">${fmtBp(d.edges[i + 1] - d.edges[i])} bin</span></div>`);
      lines.push(`<div class="tt-row"><span class="tt-key">mismatches</span><span class="tt-val">${fmtPct(row.density.mismatch[i], 2)}</span></div>`);
      lines.push(`<div class="tt-row"><span class="tt-key">deleted</span><span class="tt-val">${fmtPct(row.density.deletion[i], 2)}</span></div>`);
      lines.push(`<div class="tt-row"><span class="tt-key">insertions</span><span class="tt-val">${fmtInt(row.density.insertion[i])}</span></div>`);
    }
    showTooltip(e.clientX, e.clientY, lines.join(''));
  }

  onClick(e) {
    if (!this.data) return;
    const hit = this.hitTest(e);
    if (hit.lane !== 'variants') return;
    const v = this.variantAt(hit.pos);
    if (v) this.viewport.center(v.pos + 0.5, Math.min(this.viewport.span, 80));
  }

  renderTable() {
    clear(this.table);
    const d = this.data;
    if (!d) return;
    const v = this.viewport;
    const list = (d.variants || []).filter((x) => x.pos >= v.start && x.pos < v.end);
    const rows = d.rows.filter((r) => !r.error);
    this.table.appendChild(h('div', { style: { display: 'flex', gap: '16px', alignItems: 'baseline', margin: '16px 0 8px', flexWrap: 'wrap' } },
      h('h2', null, `Variants in view (${fmtInt(list.length)})`),
      h('span', { class: 'muted' }, rows.map((r) => `${r.sample} ${r.haplotype}: ${fmtInt(r.mismatches)} mismatches, ${fmtInt(r.deletions)} deleted bp, ${(r.insertions || []).length} insertions`).slice(0, 4).join(' · '))));
    if (!list.length) { this.table.appendChild(h('div', { class: 'muted' }, this.cfg.coords === 'haplotype' ? 'Variant positions are not defined in raw haplotype coordinates.' : 'No variants carried by the pinned samples in this view.')); return; }
    const pinned = new Set(state.pinned);
    this.table.appendChild(h('div', { class: 'table-wrap card', style: { maxHeight: '360px' } }, h('table', { class: 'table variants-table' },
      h('thead', null, h('tr', null, h('th', null, 'Position'), h('th', null, 'ID'), h('th', null, 'Ref'), h('th', null, 'Alt'), h('th', null, 'Type'), h('th', null, 'Carriers (GT)'))),
      h('tbody', null, list.slice(0, 400).map((x) => h('tr', { class: 'clickable', onclick: () => this.viewport.center(x.pos + 0.5, Math.min(this.viewport.span, 80)) },
        h('td', null, `${this.info.chromosome}:${fmtInt(x.genomic)}`), h('td', null, x.id || '–'), h('td', null, x.ref), h('td', null, x.alt),
        h('td', null, h('span', { class: 'pill' }, x.type)),
        h('td', null, x.carriers.filter((c) => pinned.has(c.sample)).map((c) => `${c.sample} ${c.gt}`).join(', '))))))));
  }

  onEvent(topic) { if (topic === 'pinned') this.fetch(true); }

  update(params) { if (params.gene && params.gene !== this.cfg.gene) this.setGene(params.gene, {}); }

  unmount() {
    this.latest.abort();
    this.fetchDebounced.cancel();
    hideTooltip();
    for (const d of this.disposers) d();
  }
}

function variantHtml(v) {
  const carriers = v.carriers.map((c) => `${escapeHtml(c.sample)} <span class="muted">${escapeHtml(c.gt)}</span>`).slice(0, 8).join('<br>');
  return `<div class="tt-title">${escapeHtml(v.type)} · ${fmtInt(v.genomic)}</div>`
    + `<div class="tt-row"><span class="tt-key">${escapeHtml(v.id || 'no id')}</span><span class="tt-val mono">${escapeHtml(v.ref)} → ${escapeHtml(v.alt)}</span></div>`
    + `<div class="tt-sub">${carriers}${v.carriers.length > 8 ? '<br>…' : ''}</div>`;
}
