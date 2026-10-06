// Variant: one site across the cohort. Genotype counts and allele frequencies by a sample field,
// AlphaGenome's predicted signal by genotype over a region (an in-silico eQTL: slope per ALT allele,
// ALT vs REF haplotypes) and GTEx's measured eQTL of the same variant (NES and donors' expression by
// genotype per tissue), with whether the two directions agree.
import { api, apiJob, isAbort, Latest } from '../api.js';
import { navigate, updateRouteParams } from '../app.js';
import { state, ds, geneInfo, cohortDescription, categoricalFields } from '../state.js';
import { h, clear, icon, select, field, fmtInt, fmtNum, toast, errorBox, jobOverlay } from '../ui.js';
import { genotypePlot, onResize, css } from '../plot.js';
import { loadGeneAnnotations } from '../gene_lanes.js';
import { openGeneCard } from '../cards.js';
import { dbsnpUrl } from '../links.js';
import { exportMenu, slug } from '../figure.js';

const OUTPUT_LABELS = { rna_seq: 'RNA-seq', cage: 'CAGE', procap: 'PRO-cap', dnase: 'DNase', atac: 'ATAC', chip_histone: 'ChIP histone', chip_tf: 'ChIP TF', splice_sites: 'Splice sites', splice_site_usage: 'Splice usage' };
const REGIONS = [
  { value: 'body', label: 'Gene body' },
  { value: 'promoter', label: 'Promoter (TSS ± 1 kb)' },
  { value: 'variant', label: 'Variant ± 1 kb' },
  { value: 'custom', label: 'Custom…' },
];
const SIGNIFICANT = 0.05;

export async function mount(root, params) {
  const page = new VariantPage(root);
  await page.init(params);
  return page;
}

/** "28120472", "28120472:A:G", "chr15:28,120,472 A>G", "chr15_28120472_A_G_b38" or an rsID. */
export function parseVariantText(text) {
  const s = String(text || '').trim().replace(/,/g, '');
  if (!s) return null;
  let m = /^(rs\d+)$/i.exec(s);
  if (m) return { rsid: m[1].toLowerCase() };
  m = /^(?:chr)?([0-9XYM]+)_(\d+)_([ACGTN]+)_([ACGTN]+)(?:_b38)?$/i.exec(s);
  if (m) return { chrom: m[1], pos: Number(m[2]), ref: m[3].toUpperCase(), alt: m[4].toUpperCase() };
  m = /^(?:(?:chr)?([0-9XYM]+)[:\s])?\s*(\d+)(?:\s*[:\s]\s*([ACGTN]+)\s*(?:>|:|\/|\s)\s*([ACGTN]+))?$/i.exec(s);
  if (m) return { chrom: m[1] || null, pos: Number(m[2]), ref: m[3] ? m[3].toUpperCase() : null, alt: m[4] ? m[4].toUpperCase() : null };
  return null;
}

class VariantPage {
  constructor(root) {
    this.root = root;
    this.latest = new Latest();
    this.disposers = [];
    this.site = null;
    this.effect = null;
    this.gtex = null;
    this.annotations = [];
    this.plots = []; // [{canvas, draw}] redrawn on resize and for figure export
  }

  async init(params) {
    const genes = state.summary.genes.map((g) => g.gene);
    const saved = state.locus || {};
    this.cfg = {
      gene: genes.includes(params.gene) ? params.gene : (genes.includes(saved.gene) ? saved.gene : genes[0]),
      pos: params.pos ? Number(params.pos) : null,
      ref: params.ref || null,
      alt: params.alt || null,
      target: params.target || null,
      region: params.region || 'body',
      custom: params.custom || '',
      output: params.output || null,
      track: params.track !== undefined && params.track !== '' ? Number(params.track) : null,
      trackAuto: params.track === undefined || params.track === '', // pick the strongest track on the gene's strand
      tissues: params.tissues ? params.tissues.split(',') : null,
      field: (categoricalFields().find((f) => f.name === 'superpopulation') || categoricalFields()[0] || {}).name || '',
    };
    this.build(genes);
    if (this.cfg.pos && this.cfg.ref && this.cfg.alt) await this.load();
    else this.showSites();
  }

  build(genes) {
    this.geneSelect = select(genes.map((g) => ({ value: g, label: g })), this.cfg.gene, (v) => { this.cfg.gene = v; this.cfg.pos = null; this.cfg.target = null; this.persist(); this.showSites(); });
    this.input = h('input', { class: 'input mono', placeholder: 'rs12913832 · chr15:28120472 A>G · 28120472', style: { width: '300px' }, spellcheck: 'false' });
    this.input.addEventListener('keydown', (e) => { if (e.key === 'Enter') this.go(); });
    if (this.cfg.pos) this.input.value = `${this.cfg.pos}${this.cfg.ref ? ` ${this.cfg.ref}>${this.cfg.alt}` : ''}`;
    this.content = h('div', { class: 'variant-content' });
    this.scroll = h('div', { class: 'page-inner variant-page' },
      h('div', { class: 'page-head' },
        h('div', null, h('h1', null, 'Variant'), h('p', null, 'Genotypes across the cohort, AlphaGenome’s prediction by genotype and GTEx’s measured eQTL of the same variant.')),
        h('div', { class: 'toolbar' },
          field('Window', this.geneSelect),
          field('Variant', h('div', { style: { display: 'flex', gap: '6px' } }, this.input, h('button', { class: 'btn primary', onclick: () => this.go() }, icon('search', 14), 'Open'))),
          exportMenu(() => this.figureOptions()))),
      this.content);
    this.root.appendChild(this.scroll);
    this.disposers.push(onResize(this.content, () => this.redraw()));
  }

  persist() {
    const c = this.cfg;
    updateRouteParams({ gene: c.gene, pos: c.pos, ref: c.ref, alt: c.alt, target: c.target, region: c.region !== 'body' ? c.region : null, custom: c.region === 'custom' ? c.custom : null, output: c.output, track: c.track, tissues: c.tissues ? c.tissues.join(',') : null });
  }

  redraw() { for (const p of this.plots) p.draw(); }

  // ------------------------------------------------------------------ navigation
  async go() {
    const parsed = parseVariantText(this.input.value);
    if (!parsed) { toast('Type an rsID, chr:pos ref>alt, a GTEx variant id or a position', 'error'); return; }
    let { pos, ref, alt } = parsed;
    if (parsed.rsid) {
      try {
        const r = await api('/api/gtex/resolve', { params: { rsid: parsed.rsid } });
        ({ pos, ref, alt } = r);
        const window = (await geneInfo(this.cfg.gene));
        if (window.chromosome && r.chromosome && window.chromosome.replace(/^chr/, '') !== r.chromosome.replace(/^chr/, '')) {
          toast(`${parsed.rsid} is on ${r.chromosome}; the ${this.cfg.gene} window is on ${window.chromosome}`, 'error', 8000);
          return;
        }
      } catch (err) { toast(err.message, 'error', 8000); return; }
    }
    if (!ref || !alt) { this.showSites(pos); return; }
    Object.assign(this.cfg, { pos, ref, alt });
    this.persist();
    await this.load();
  }

  update(params) {
    if (params.pos && (Number(params.pos) !== this.cfg.pos || params.ref !== this.cfg.ref || params.alt !== this.cfg.alt || (params.gene && params.gene !== this.cfg.gene))) {
      if (params.gene) { this.cfg.gene = params.gene; this.geneSelect.value = params.gene; }
      Object.assign(this.cfg, { pos: Number(params.pos), ref: params.ref, alt: params.alt });
      this.input.value = `${this.cfg.pos} ${this.cfg.ref}>${this.cfg.alt}`;
      this.load();
    }
  }

  /** Sites of the window carried in the cohort, around ``center`` (or the common ones). */
  async showSites(center = null) {
    this.plots = [];
    clear(this.content);
    const body = h('div', { class: 'card-body', style: { position: 'relative', minHeight: '120px' } });
    this.content.appendChild(h('section', { class: 'card' }, h('div', { class: 'card-head' }, h('h2', null, center ? `Sites near ${fmtInt(center)}` : `Common sites in ${this.cfg.gene} (allele frequency ≥ 5%)`)), body));
    const overlay = jobOverlay(body, `Genotypes of the cohort in ${this.cfg.gene}`);
    try {
      const params = center ? { gene: this.cfg.gene, start: center - 50, end: center + 50 } : { gene: this.cfg.gene, min_af: 0.05, limit: 3000 };
      const data = await apiJob(`${ds()}/variant/sites`, { params, signal: this.latest.next(), onProgress: (job) => overlay.update(job) });
      overlay.remove();
      if (!data.sites.length) { body.appendChild(h('div', { class: 'empty' }, center ? 'No sample of the dataset carries a variant there.' : 'No common sites.')); return; }
      body.appendChild(h('p', { class: 'muted', style: { margin: '0 0 8px' } }, `${fmtInt(data.sites.length)}${data.truncated ? '+' : ''} sites · frequencies over all ${fmtInt(data.samples)} samples · click one to open it`));
      body.appendChild(h('div', { class: 'table-wrap', style: { maxHeight: '60vh' } }, h('table', { class: 'table' },
        h('thead', null, h('tr', null, ['Position', 'REF > ALT', 'Type', 'ALT frequency', 'ID'].map((c) => h('th', null, c)))),
        h('tbody', null, data.sites.map((s) => h('tr', { class: 'clickable', onclick: () => { Object.assign(this.cfg, { pos: s.pos, ref: s.ref, alt: s.alt }); this.input.value = `${s.pos} ${s.ref}>${s.alt}`; this.persist(); this.load(); } },
          h('td', { class: 'mono' }, `${data.chromosome}:${fmtInt(s.pos)}`), h('td', { class: 'mono' }, `${trim(s.ref)} > ${trim(s.alt)}`), h('td', null, h('span', { class: 'pill' }, s.type)),
          h('td', null, fmtNum(s.af, 3)), h('td', { class: 'muted mono' }, s.id || '–')))))));
    } catch (err) {
      overlay.remove();
      if (!isAbort(err)) body.appendChild(errorBox(err));
    }
  }

  // ------------------------------------------------------------------ the variant
  async load() {
    this.plots = [];
    this.effect = null; this.gtex = null;
    clear(this.content);
    const signal = this.latest.next();
    this.signal = signal;
    const headBody = h('div', { class: 'card-body', style: { position: 'relative', minHeight: '80px' } });
    this.content.appendChild(h('section', { class: 'card' }, headBody));
    const overlay = jobOverlay(headBody, `Genotypes of the cohort in ${this.cfg.gene}`);
    try {
      [this.site, this.info] = await Promise.all([
        apiJob(`${ds()}/variant/site`, { params: { gene: this.cfg.gene, pos: this.cfg.pos, ref: this.cfg.ref, alt: this.cfg.alt, field: this.cfg.field, filters: state.filters }, signal, onProgress: (job) => overlay.update(job) }),
        geneInfo(this.cfg.gene),
      ]);
    } catch (err) {
      overlay.remove();
      if (!isAbort(err)) headBody.appendChild(errorBox(err));
      return;
    }
    overlay.remove();
    try { this.annotations = (await loadGeneAnnotations(state.datasetId, this.cfg.gene)).genes || []; } catch (err) { this.annotations = []; }
    this.chooseDefaults();
    this.renderHeader(headBody);
    this.effectCard = h('section', { class: 'card', style: { marginTop: '16px' } });
    this.gtexCard = h('section', { class: 'card', style: { marginTop: '16px' } });
    this.verdictCard = h('section', { class: 'card', style: { marginTop: '16px' } });
    this.content.append(this.verdictCard, this.effectCard, this.gtexCard);
    this.renderVerdict();
    await Promise.all([this.loadEffect(), this.loadGtex()]);
  }

  /** Target gene (overlapping or nearest annotated gene), output and a track on its strand. */
  chooseDefaults() {
    const c = this.cfg;
    const offset = c.pos - this.site.window.start;
    const coding = this.annotations.filter((g) => g.type === 'protein_coding');
    const pool = coding.length ? coding : this.annotations;
    if (!c.target || !this.annotations.some((g) => g.name === c.target)) {
      const named = pool.find((g) => g.name === c.gene);
      const containing = pool.filter((g) => g.start <= offset && offset < g.end).sort((a, b) => (a.end - a.start) - (b.end - b.start));
      const nearest = [...pool].sort((a, b) => dist(a, offset) - dist(b, offset))[0];
      c.target = (named || containing[0] || nearest || { name: c.gene }).name;
    }
    const outputs = Object.keys(this.info.outputs).filter((o) => (this.info.outputs[o].tracks || []).length && OUTPUT_LABELS[o]);
    if (!outputs.includes(c.output)) c.output = outputs.includes('rna_seq') ? 'rna_seq' : outputs[0];
    const tracks = (this.info.outputs[c.output] || { tracks: [] }).tracks;
    if (!(Number.isInteger(c.track) && c.track >= 0 && c.track < tracks.length)) {
      const strand = (this.targetGene() || {}).strand;
      const onStrand = tracks.findIndex((t) => (t.metadata || {}).strand === strand);
      c.track = onStrand >= 0 ? onStrand : 0;
    }
  }

  targetGene() { return this.annotations.find((g) => g.name === this.cfg.target) || null; }

  /** Index of the track with the highest cohort mean over the region, on the target gene's strand when known. */
  strongestTrack(means) {
    const tracks = (this.info.outputs[this.cfg.output] || { tracks: [] }).tracks;
    const strand = (this.targetGene() || {}).strand;
    let best = null;
    (means || []).forEach((m, i) => {
      if (m === null || !Number.isFinite(m)) return;
      const s = (tracks[i] && tracks[i].metadata || {}).strand;
      if (strand && s && s !== '.' && s !== strand) return;
      if (best === null || m > means[best]) best = i;
    });
    return best;
  }

  /** [start, end) in window offsets for the chosen region. */
  regionRange() {
    const c = this.cfg;
    const offset = c.pos - this.site.window.start;
    const gene = this.targetGene();
    const length = this.site.window.end - this.site.window.start + 1;
    const clamp = ([a, b]) => [Math.max(0, Math.round(a)), Math.min(length, Math.round(b))];
    if (c.region === 'custom') {
      const m = /(?:chr\w+:)?([\d,]+)\s*-\s*([\d,]+)/.exec(c.custom || '');
      if (m) return clamp([Number(m[1].replace(/,/g, '')) - this.site.window.start, Number(m[2].replace(/,/g, '')) - this.site.window.start + 1]);
    }
    if (c.region === 'promoter' && gene) {
      const tss = gene.strand === '-' ? gene.end : gene.start;
      return clamp([tss - 1000, tss + 1000]);
    }
    if ((c.region === 'body' || c.region === 'promoter') && gene) return clamp([gene.start, gene.end]);
    return clamp([offset - 1000, offset + 1000]);
  }

  regionText() {
    const [a, b] = this.regionRange();
    const w = this.site.window;
    return `${w.chromosome}:${fmtInt(w.start + a)}-${fmtInt(w.start + b - 1)}`;
  }

  trackLabel() {
    const t = (this.info.outputs[this.cfg.output] || { tracks: [] }).tracks[this.cfg.track] || {};
    return `${OUTPUT_LABELS[this.cfg.output] || this.cfg.output} · ${t.label || `track ${this.cfg.track}`}`;
  }

  renderHeader(body) {
    const s = this.site;
    const co = s.cohort;
    const freq = h('table', { class: 'table compact' },
      h('thead', null, h('tr', null, h('th', null, co.field || 'Group'), h('th', null, 'n'), h('th', null, 'ALT freq.'), ...s.labels.map((l) => h('th', null, l)))),
      h('tbody', null,
        h('tr', null, h('td', null, h('b', null, 'Cohort')), h('td', null, fmtInt(co.cohort)), h('td', null, h('b', null, fmtNum(co.af, 3))), ...['0/0', '0/1', '1/1'].map((k) => h('td', null, fmtInt(co.counts[k])))),
        co.by_group.map((g) => h('tr', null, h('td', null, g.group), h('td', null, fmtInt(g.n)), h('td', null, afBar(g.af)), ...['0/0', '0/1', '1/1'].map((k) => h('td', null, fmtInt(g.counts[k])))))));
    this.rsidEl = h('span', { class: 'muted' }, 'rsID: looking up…');
    body.append(h('div', { class: 'variant-head' },
      h('div', null,
        h('div', { class: 'variant-title mono' }, `${s.chromosome}:${fmtInt(s.pos)} ${trim(s.ref)} > ${trim(s.alt)}`),
        h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap', alignItems: 'center', marginTop: '6px' } },
          h('span', { class: 'pill' }, s.type), this.rsidEl,
          h('span', { class: 'muted mono' }, s.variant_id),
          h('a', { href: `https://gnomad.broadinstitute.org/variant/${s.chromosome.replace(/^chr/, '')}-${s.pos}-${s.ref}-${s.alt}?dataset=gnomad_r4`, target: '_blank', rel: 'noopener' }, 'gnomAD ', icon('external', 12)),
          h('a', { href: `https://gtexportal.org/home/snp/${s.variant_id}`, target: '_blank', rel: 'noopener' }, 'GTEx ', icon('external', 12))),
        h('dl', { class: 'kv', style: { marginTop: '12px' } },
          h('dt', null, 'Window'), h('dd', null, h('a', { href: '#', onclick: (e) => { e.preventDefault(); openGeneCard(s.gene); } }, s.gene), ` · ${s.window.chromosome}:${fmtInt(s.window.start)}-${fmtInt(s.window.end)}`),
          h('dt', null, 'Dataset'), h('dd', null, `${state.datasetId} · ALT frequency ${fmtNum(s.af_dataset, 3)} over all samples`),
          h('dt', null, 'Cohort'), h('dd', null, cohortDescription()),
          h('dt', null, 'Genotypes'), h('dd', null, 'from the samples’ phased window VCFs; samples without a call are homozygous REF')),
        h('div', { style: { display: 'flex', gap: '8px', marginTop: '12px', flexWrap: 'wrap' } },
          h('button', { class: 'btn small', title: 'Group means of the cohort by genotype at this site on the Tracks page', onclick: () => navigate('tracks', { gene: s.gene, mode: 'groups', groupField: `variant:${s.pos}:${s.ref}:${s.alt}`, start: Math.max(0, s.pos - s.window.start - 2000), end: s.pos - s.window.start + 2000 }) }, icon('tracks', 14), 'Tracks by genotype →'),
          h('button', { class: 'btn small ghost', onclick: () => navigate('sequence', { gene: s.gene, start: Math.max(0, s.pos - s.window.start - 40), end: s.pos - s.window.start + 40 }) }, 'Sequence →'),
          h('button', { class: 'btn small ghost', onclick: () => this.showSites(s.pos) }, 'Other sites here'))),
      h('div', { class: 'table-wrap' }, freq)));
  }

  // ------------------------------------------------------------------ AlphaGenome by genotype
  effectControls() {
    const c = this.cfg;
    const genes = [...this.annotations].sort((a, b) => a.name.localeCompare(b.name));
    const outputs = Object.keys(this.info.outputs).filter((o) => (this.info.outputs[o].tracks || []).length && OUTPUT_LABELS[o]);
    const tracks = (this.info.outputs[c.output] || { tracks: [] }).tracks;
    const refresh = () => { this.persist(); this.loadEffect(); this.renderVerdict(); };
    const custom = h('input', { class: 'input mono', value: c.custom || this.regionText(), style: { width: '230px' }, hidden: c.region !== 'custom' });
    custom.addEventListener('keydown', (e) => { if (e.key === 'Enter') { c.custom = custom.value; refresh(); } });
    return h('div', { class: 'toolbar', 'data-export': 'skip', style: { padding: '12px 14px', borderBottom: '1px solid var(--line)' } },
      genes.length ? field('Target gene', select(genes.map((g) => ({ value: g.name, label: `${g.name} (${g.strand})` })), c.target, (v) => { c.target = v; c.track = null; c.trackAuto = true; this.chooseDefaults(); this.loadGtex(); this.renderEffectCardShell(); refresh(); })) : null,
      field('Region', h('div', { style: { display: 'flex', gap: '6px' } }, select(REGIONS.filter((r) => genes.length || r.value === 'variant' || r.value === 'custom'), c.region, (v) => { c.region = v; custom.hidden = v !== 'custom'; if (v !== 'custom') refresh(); }), custom)),
      field('Output', select(outputs.map((o) => ({ value: o, label: OUTPUT_LABELS[o] })), c.output, (v) => { c.output = v; c.track = null; c.trackAuto = true; this.chooseDefaults(); this.renderEffectCardShell(); refresh(); })),
      field('Track', select(tracks.map((t, i) => ({ value: String(i), label: t.label || `track ${i}` })), String(c.track), (v) => { c.track = Number(v); c.trackAuto = false; refresh(); }, { style: { maxWidth: '340px' } })));
  }

  renderEffectCardShell() {
    clear(this.effectCard);
    this.effectBody = h('div', { class: 'card-body', style: { position: 'relative', minHeight: '260px' } });
    this.effectCard.append(
      h('div', { class: 'card-head' }, h('h2', null, 'AlphaGenome by genotype'), h('span', { class: 'muted', style: { fontSize: '12px' } }, 'mean predicted signal per sample (H1, H2 averaged) over the region')),
      this.effectControls(), this.effectBody);
  }

  async loadEffect() {
    if (!this.effectBody) this.renderEffectCardShell();
    const body = this.effectBody;
    clear(body);
    this.plots = this.plots.filter((p) => p.kind !== 'effect');
    const [start, end] = this.regionRange();
    if (end <= start) { body.appendChild(errorBox(new Error('Empty region'))); return; }
    const overlay = jobOverlay(body, `${this.cfg.gene}: ${OUTPUT_LABELS[this.cfg.output]} per haplotype over ${this.regionText()}`);
    const c = this.cfg;
    const token = (this.effectToken = (this.effectToken || 0) + 1);
    this.effect = null;
    try {
      const effect = await apiJob(`${ds()}/variant/effect`, { params: { gene: c.gene, pos: c.pos, ref: c.ref, alt: c.alt, output: c.output, track: c.track, start, end, filters: state.filters }, signal: this.signal, onProgress: (job) => overlay.update(job) });
      if (token !== this.effectToken) { overlay.remove(); return; } // superseded by a newer choice
      const best = c.trackAuto ? this.strongestTrack(effect.track_means) : null;
      if (best !== null && best !== c.track) {
        // Region means cover every track, so re-asking for the strongest one is instant.
        overlay.remove();
        c.track = best;
        c.trackAuto = false;
        this.renderEffectCardShell();
        this.persist();
        return this.loadEffect();
      }
      this.effect = effect;
    } catch (err) {
      overlay.remove();
      if (!isAbort(err) && token === this.effectToken) body.appendChild(errorBox(err));
      return;
    }
    overlay.remove();
    const e = this.effect;
    const r = e.regression;
    const base = (e.groups[0].summary || {}).mean;
    const canvas = h('canvas');
    const stats = h('dl', { class: 'kv' },
      h('dt', null, 'Region'), h('dd', null, `${this.regionText()} (${fmtInt(end - start)} bp, ${REGIONS.find((x) => x.value === c.region).label.toLowerCase()}${this.targetGene() && c.region !== 'variant' && c.region !== 'custom' ? ` of ${c.target}` : ''})`),
      h('dt', null, 'Track'), h('dd', null, this.trackLabel(), h('span', { class: 'muted' }, ` · cohort mean ${fmtNum((e.track_means || [])[c.track], 3)}`)),
      h('dt', null, 'Slope per ALT allele'), h('dd', null, h('b', { class: directionClass(r.slope, r.p) }, signed(r.slope)), ` ± ${fmtNum(r.se, 3)} · ${pct(r.slope, base)} of the REF/REF mean · t = ${fmtNum(r.t, 3)} · p = ${fmtP(r.p)} · n = ${fmtInt(r.n)}`),
      h('dt', null, 'ALT vs REF haplotypes'), h('dd', null, `log2 fold change ${signed(e.haplotypes.log2fc)} (${fmtInt(e.haplotypes.alt.n || 0)} ALT, ${fmtInt(e.haplotypes.ref.n || 0)} REF haplotypes)`),
      h('dt', null, 'Reference genome'), h('dd', null, `${fmtNum(e.reference_genome, 4)} (AlphaGenome on hg38, dashed)`),
      h('dt', null, 'Note'), h('dd', { class: 'muted' }, 'An association across the cohort, like an eQTL: variants in linkage disequilibrium with this one move with it.'));
    body.append(h('div', { class: 'variant-grid' }, h('div', { class: 'canvas-host' }, canvas), stats));
    const plot = {
      kind: 'effect', canvas,
      draw: () => genotypePlot(canvas, {
        groups: e.groups.map((g, i) => ({ label: e.labels[i], values: g.values, summary: g.summary })),
        width: Math.max(320, canvas.parentElement.clientWidth || 420), height: 260, yLabel: 'AlphaGenome (mean over region)',
        fit: { intercept: r.intercept, slope: r.slope }, refLine: e.reference_genome, refLabel: 'reference genome', color: css('--accent'),
      }),
    };
    this.plots.push(plot);
    plot.draw();
    this.renderVerdict();
  }

  // ------------------------------------------------------------------ GTEx
  async loadGtex() {
    const card = this.gtexCard;
    clear(card);
    this.plots = this.plots.filter((p) => p.kind !== 'gtex');
    const body = h('div', { class: 'card-body', style: { position: 'relative', minHeight: '120px' } });
    const tissueHost = h('div', { style: { display: 'flex', gap: '6px', flexWrap: 'wrap', alignItems: 'center' } });
    card.append(h('div', { class: 'card-head' }, h('h2', null, `GTEx v8 eQTL · ${this.cfg.target}`), tissueHost), body);
    body.appendChild(h('div', { class: 'muted' }, 'Asking the GTEx portal…'));
    let tissues = [];
    try { tissues = (await api('/api/gtex/tissues')).tissues; } catch (err) { /* offline: defaults only */ }
    const chosen = this.cfg.tissues || ['Skin_Sun_Exposed_Lower_leg', 'Skin_Not_Sun_Exposed_Suprapubic'];
    const name = (id) => (tissues.find((t) => t.id === id) || { name: id.replace(/_/g, ' ') }).name;
    for (const id of chosen) tissueHost.appendChild(h('span', { class: 'pill removable' }, name(id), h('button', { title: 'Remove', onclick: () => { this.cfg.tissues = chosen.filter((x) => x !== id); this.persist(); this.loadGtex(); } }, '×')));
    if (tissues.length) {
      tissueHost.appendChild(select([{ value: '', label: '+ tissue' }, ...tissues.filter((t) => !chosen.includes(t.id)).map((t) => ({ value: t.id, label: t.name }))], '', (v) => { if (v) { this.cfg.tissues = [...chosen, v].slice(0, 8); this.persist(); this.loadGtex(); } }, { style: { maxWidth: '200px' } }));
    }
    const token = (this.gtexToken = (this.gtexToken || 0) + 1);
    this.gtex = null;
    try {
      const gtex = await api('/api/gtex/variant', { params: { variant_id: this.site.variant_id, gene: this.cfg.target, tissues: chosen.join(',') }, signal: this.signal });
      if (token !== this.gtexToken) return;
      this.gtex = gtex;
    } catch (err) {
      if (token !== this.gtexToken) return;
      clear(body);
      if (!isAbort(err)) body.appendChild(errorBox(err));
      return;
    }
    clear(body);
    const g = this.gtex;
    this.rsidEl.replaceChildren(g.rsid ? h('a', { href: dbsnpUrl(g.rsid), target: '_blank', rel: 'noopener' }, g.rsid, ' ', icon('external', 12)) : h('span', { class: 'muted' }, 'no rsID in GTEx'));
    if (!g.in_gtex) {
      body.appendChild(h('div', { class: 'empty' }, 'Not in GTEx v8: only variants with MAF ≥ 1% among GTEx donors were tested.'));
      this.renderVerdict();
      return;
    }
    if (!g.gene || !g.gene.gencode_id) body.appendChild(h('p', { class: 'muted' }, `${this.cfg.target} is not in GTEx's gene annotation (GENCODE v26).`));
    const grid = h('div', { class: 'gtex-grid' });
    for (const d of g.dynamic || []) {
      const box = h('div', { class: 'gtex-tissue' }, h('div', { class: 'gtex-tissue-head' }, h('b', null, name(d.tissue)),
        d.error ? h('span', { class: 'muted' }, d.error) : h('span', null, 'NES ', h('b', { class: directionClass(d.nes, d.p) }, signed(d.nes)), ` · p = ${fmtP(d.p)}`)));
      grid.appendChild(box);
      if (d.error) continue;
      const canvas = h('canvas');
      box.appendChild(h('div', { class: 'canvas-host' }, canvas));
      const plot = {
        kind: 'gtex', canvas,
        draw: () => genotypePlot(canvas, {
          groups: ['0/0', '0/1', '1/1'].map((k, i) => ({ label: this.site.labels[i], values: d.values[k], summary: summarize(d.values[k], d.counts[k]) })),
          width: Math.max(260, canvas.parentElement.clientWidth || 320), height: 220, yLabel: 'normalized expression', color: css('--series-2'),
        }),
      };
      this.plots.push(plot);
      plot.draw();
    }
    if ((g.dynamic || []).length) body.appendChild(grid);
    const sig = g.significant || [];
    body.appendChild(h('h3', { style: { margin: '16px 0 6px', fontSize: '13px' } }, `Significant eQTLs of this variant in GTEx (${sig.length})`));
    if (!sig.length) body.appendChild(h('div', { class: 'muted' }, 'None at GTEx’s per-gene significance threshold, in any tissue.'));
    else {
      body.appendChild(h('div', { class: 'table-wrap', style: { maxHeight: '260px' } }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['Gene', 'Tissue', 'NES', 'p'].map((x) => h('th', null, x)))),
        h('tbody', null, sig.map((r) => h('tr', null, h('td', null, h('a', { href: '#', onclick: (e) => { e.preventDefault(); openGeneCard(r.gene); } }, r.gene)), h('td', null, name(r.tissue)), h('td', { class: directionClass(r.nes, r.p) }, signed(r.nes)), h('td', null, fmtP(r.p))))))));
    }
    this.renderVerdict();
  }

  // ------------------------------------------------------------------ verdict
  renderVerdict() {
    const card = this.verdictCard;
    if (!card) return;
    clear(card);
    const e = this.effect;
    const dyn = ((this.gtex && this.gtex.dynamic) || []).filter((d) => !d.error);
    const rows = dyn.map((d) => {
      const ag = e ? e.regression : null;
      let verdict = 'waiting';
      if (ag && Number.isFinite(ag.slope) && Number.isFinite(d.nes)) {
        const agSig = ag.p < SIGNIFICANT;
        const gtexSig = d.p < SIGNIFICANT;
        if (agSig && gtexSig) verdict = Math.sign(ag.slope) === Math.sign(d.nes) ? 'agree' : 'disagree';
        else verdict = agSig ? 'GTEx n.s.' : gtexSig ? 'AlphaGenome n.s.' : 'both n.s.';
      }
      return { tissue: d.tissue, nes: d.nes, p: d.p, verdict };
    });
    const ag = e && e.regression;
    card.append(h('div', { class: 'card-head' }, h('h2', null, 'Direction: AlphaGenome vs GTEx'), h('span', { class: 'muted', style: { fontSize: '12px' } }, `effect of the ALT allele (${trim(this.cfg.alt)}) on ${this.cfg.target}; significant at p < ${SIGNIFICANT}`)),
      h('div', { class: 'card-body' }, !rows.length
        ? h('div', { class: 'muted' }, this.gtex ? 'No GTEx association to compare (see below).' : 'Waiting for AlphaGenome and GTEx…')
        : h('table', { class: 'table compact' },
          h('thead', null, h('tr', null, ['GTEx tissue', 'GTEx NES', 'AlphaGenome slope', 'Verdict'].map((x) => h('th', null, x)))),
          h('tbody', null, rows.map((r) => h('tr', null,
            h('td', null, r.tissue.replace(/_/g, ' ')),
            h('td', null, h('span', { class: directionClass(r.nes, r.p) }, signed(r.nes)), h('span', { class: 'muted' }, ` p = ${fmtP(r.p)}`)),
            h('td', null, ag ? [h('span', { class: directionClass(ag.slope, ag.p) }, signed(ag.slope)), h('span', { class: 'muted' }, ` p = ${fmtP(ag.p)} · ${this.trackLabel()}`)] : '…'),
            h('td', null, h('span', { class: `verdict ${r.verdict === 'agree' || r.verdict === 'disagree' ? r.verdict : 'ns'}`, title: r.verdict.endsWith('n.s.') ? `not significant (p ≥ ${SIGNIFICANT})` : 'both significant' }, r.verdict))))))));
  }

  figureOptions() {
    if (!this.site) return null;
    const s = this.site;
    const rsid = this.gtex && this.gtex.rsid ? ` (${this.gtex.rsid})` : '';
    return {
      root: this.content,
      redraw: () => this.redraw(),
      title: `${s.chromosome}:${fmtInt(s.pos)} ${trim(s.ref)}>${trim(s.alt)}${rsid} · ${this.cfg.target}`,
      caption: [
        `AlphaGenome ${this.trackLabel()} over ${this.effect ? this.regionText() : 'the region'} by genotype; GTEx v8 eQTL of ${this.cfg.target}. Cohort: ${cohortDescription()}.`,
        `Dataset ${state.datasetId} · exported ${new Date().toISOString().slice(0, 10)} from genomics visualize`,
      ],
      filename: slug(`${state.datasetId}_${s.gene}_${s.pos}_${s.ref}_${s.alt}_${this.cfg.target}`),
    };
  }

  unmount() {
    this.latest.abort();
    for (const d of this.disposers) d();
  }
}

// ------------------------------------------------------------------ helpers
function dist(gene, offset) { return offset < gene.start ? gene.start - offset : offset >= gene.end ? offset - gene.end + 1 : 0; }
function trim(allele) { return allele.length > 12 ? `${allele.slice(0, 10)}…(${allele.length})` : allele; }
function signed(v) { return Number.isFinite(v) ? `${v > 0 ? '+' : ''}${Math.abs(v) >= 1e-2 || v === 0 ? Number(v.toPrecision(3)) : v.toExponential(2)}` : '–'; }
function fmtP(p) { return Number.isFinite(p) ? (p < 1e-3 ? p.toExponential(1) : Number(p.toPrecision(2)).toString()) : '–'; }
function pct(slope, base) { return Number.isFinite(slope) && Number.isFinite(base) && base !== 0 ? `${slope / base > 0 ? '+' : ''}${(100 * slope / base).toFixed(1)}%` : '–'; }
function directionClass(v, p) { return !Number.isFinite(v) || !(p < SIGNIFICANT) ? 'muted' : v > 0 ? 'delta up' : 'delta down'; }
function afBar(af) {
  return h('span', { class: 'af-bar', title: fmtNum(af, 4) }, h('span', { style: { width: `${Math.round(Math.max(0, Math.min(1, af)) * 100)}%` } }), h('span', { class: 'af-val' }, fmtNum(af, 3)));
}
/** Box-plot summary of GTEx donor values (the API sends the values, not their quartiles). */
function summarize(values, n) {
  const v = (values || []).filter(Number.isFinite).sort((a, b) => a - b);
  if (!v.length) return { n: n || 0 };
  const q = (f) => v[Math.min(v.length - 1, Math.max(0, Math.round(f * (v.length - 1))))];
  return { n: n || v.length, min: v[0], q1: q(0.25), median: q(0.5), q3: q(0.75), max: v[v.length - 1], mean: v.reduce((a, b) => a + b, 0) / v.length };
}
