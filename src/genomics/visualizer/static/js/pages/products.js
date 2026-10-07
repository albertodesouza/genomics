// Gene products: one individual's mature transcripts and proteins per haplotype (GENCODE models
// spliced on each haplotype's own sequence and translated), compared with the reference's, plus
// AlphaGenome's splice-junction evidence per tissue track and the candidate isoforms it implies.
import { api, apiJob, isAbort, Latest } from '../api.js';
import { updateRouteParams } from '../app.js';
import { state, ds } from '../state.js';
import { h, clear, icon, select, field, fmtInt, fmtNum, toast, errorBox, jobOverlay, downloadText, checkbox } from '../ui.js';
import { slug } from '../figure.js';
import { transcriptLink } from '../links.js';
import { expressionCard } from '../products_expression.js';
import { StructureCard } from '../products_structure.js';
import { reportCards } from '../products_report.js';
import { tissuePicker } from '../tissue_picker.js';

const DEFAULT_TISSUES = ['CL:1000458', 'CL:2000045', 'EFO:0000572', 'EFO:0002784'];

const HAPS = ['H1', 'H2'];
const ALL = ['ref', 'H1', 'H2'];
const CLASS_LABELS = {
  no_change: 'identical', utr_change: 'UTR only', synonymous: 'synonymous', missense: 'missense', stop_gained: 'stop gained',
  frameshift: 'frameshift', inframe_indel: 'in-frame indel', start_lost: 'start lost', stop_lost: 'stop lost', noncoding_change: 'changed', no_data: '–',
};
const SEVERE = new Set(['stop_gained', 'frameshift', 'start_lost', 'stop_lost']);
const MODERATE = new Set(['missense', 'inframe_indel']);
const QUIET_REGIONS = new Set(['intron', 'upstream', 'downstream']);

const sampleIds = () => (state.samples && state.samples.index ? Array.from(state.samples.index.keys()) : []);

export async function mount(root, params) {
  const page = new ProductsPage(root);
  await page.init(params);
  return page;
}

function changeBadge(change) {
  const cls = (change && change.class) || 'no_data';
  const tone = SEVERE.has(cls) ? 'bad' : MODERATE.has(cls) ? 'warn' : '';
  return h('span', { class: `pill product-class ${tone}`, title: (change && change.hgvs) || '' }, CLASS_LABELS[cls] || cls);
}

function changeCell(product) {
  if (!product) return h('td', { class: 'muted' }, '–');
  const change = product.change || {};
  const nmd = product.nmd && product.nmd.predicted;
  return h('td', null,
    h('div', { class: 'product-cell' }, changeBadge(change), nmd ? h('span', { class: 'pill product-class bad', title: product.nmd.reason }, 'NMD') : null),
    change.hgvs && change.hgvs !== 'p.(=)' ? h('div', { class: 'mono product-hgvs' }, change.hgvs) : null);
}

class ProductsPage {
  constructor(root) {
    this.root = root;
    this.latest = new Latest();
    this.data = null;
    this.selected = null; // {kind: 'tx'|'cand', index}
  }

  async init(params) {
    const genes = state.summary.genes.map((g) => g.gene);
    const saved = state.locus || {};
    const ids = sampleIds();
    this.cfg = {
      gene: genes.includes(params.gene) ? params.gene : (genes.includes(saved.gene) ? saved.gene : genes[0]),
      sample: params.sample || state.pinned[0] || ids[0] || '',
      target: params.target || null,
      track: params.track !== undefined && params.track !== '' ? Number(params.track) : null,
      intronic: params.intronic === '1',
      tissue: /^[A-Za-z]+:\d+$/.test(params.tissue || '') ? params.tissue : 'CL:1000458',
    };
    this.reportState = { allVariants: false, selectedTx: null, selectedHap: null };
    this.build(genes);
    await this.load();
  }

  update(params) {
    if ((params.gene && params.gene !== this.cfg.gene) || (params.sample && params.sample !== this.cfg.sample)) {
      Object.assign(this.cfg, { gene: params.gene || this.cfg.gene, sample: params.sample || this.cfg.sample, target: params.target || null });
      this.geneSelect.value = this.cfg.gene;
      this.sampleInput.value = this.cfg.sample;
      this.load();
    }
  }

  unmount() { this.latest.abort(); }

  build(genes) {
    this.geneSelect = select(genes.map((g) => ({ value: g, label: g })), this.cfg.gene, (v) => { this.cfg.gene = v; this.cfg.target = null; this.load(); });
    const list = h('datalist', { id: 'products-samples' }, sampleIds().slice(0, 5000).map((id) => h('option', { value: id })));
    this.sampleInput = h('input', { class: 'input mono', value: this.cfg.sample, list: 'products-samples', style: { width: '130px' }, spellcheck: 'false' });
    this.sampleInput.addEventListener('change', () => { this.cfg.sample = this.sampleInput.value.trim(); this.load(); });
    this.targetSelect = select([], null, (v) => { this.cfg.target = v; this.load(); });
    this.tissueSelect = tissuePicker(this.cfg.tissue, (v) => { this.cfg.tissue = v; this.persist(); this.exprCard = null; this.loadReport(); if (this.evidenceOpen) this.render(); });
    const exportBtn = h('div', { class: 'toolbar-group' },
      h('button', { class: 'btn', onclick: () => this.download('protein') }, icon('download', 14), 'Proteins FASTA'),
      h('button', { class: 'btn', onclick: () => this.download('mrna') }, icon('download', 14), 'mRNA FASTA'),
      h('button', { class: 'btn', onclick: () => this.downloadJson() }, icon('download', 14), 'JSON'));
    this.content = h('div', { class: 'products-content' });
    this.root.appendChild(h('div', { class: 'page-inner products-page' },
      h('div', { class: 'page-head' },
        h('div', null, h('h1', null, 'Gene products'),
          h('p', null, '(1) What the reference genome makes of a gene in a tissue, transcript by transcript, from GTEx and the Human Protein Atlas. (2) What each of the individual’s two copies changes: how much of each transcript, unstable mRNA, premature stops, different amino acids, protein similarity. (3) The protein’s 3D shape.')),
        h('div', { class: 'toolbar' }, field('Window', this.geneSelect), list, field('Sample', this.sampleInput), field('Gene', this.targetSelect), field('Tissue', this.tissueSelect), exportBtn)),
      this.content));
  }

  persist() {
    const c = this.cfg;
    updateRouteParams({ gene: c.gene, sample: c.sample, target: c.target, track: c.track, intronic: c.intronic ? '1' : null, tissue: c.tissue });
  }

  params(extra = {}) { return { gene: this.cfg.gene, sample: this.cfg.sample, target: this.cfg.target, ...extra }; }

  async load() {
    this.persist();
    clear(this.content);
    this.data = null;
    this.selected = null;
    const host = h('section', { class: 'card', style: { position: 'relative', minHeight: '120px' } });
    this.content.appendChild(host);
    if (!this.cfg.sample) { host.appendChild(h('div', { class: 'empty' }, 'Pick a sample.')); return; }
    const overlay = jobOverlay(host, `Gene products of ${this.cfg.sample} in ${this.cfg.gene}`);
    try {
      this.data = await apiJob(`${ds()}/products`, { params: this.params(), signal: this.latest.next(), onProgress: (job) => overlay.update(job) });
    } catch (err) {
      overlay.remove();
      // 404: no gene annotations, or no window for this sample (data the dataset lacks, not a failure)
      if (!isAbort(err)) host.appendChild(err.status === 404 ? h('div', { class: 'empty' }, err.message) : errorBox(err));
      return;
    }
    overlay.remove();
    host.remove();
    const d = this.data;
    const genes = d.genes || [];
    this.targetSelect.disabled = genes.length < 2;
    clear(this.targetSelect);
    for (const g of genes) this.targetSelect.appendChild(h('option', { value: g }, g));
    if (d.target) this.targetSelect.value = d.target;
    if (!d.transcripts || !d.transcripts.length) {
      this.content.appendChild(h('section', { class: 'card' }, h('div', { class: 'card-body' }, h('div', { class: 'empty' }, d.message || 'No transcripts.'))));
      return;
    }
    const sp = d.splicing || {};
    if (this.cfg.track === null || !sp.tracks || this.cfg.track >= sp.tracks.length) {
      const expr = sp.expression || [];
      this.cfg.track = expr.length ? expr.indexOf(Math.max(...expr)) : 0;
    }
    const base = d.transcripts.findIndex((t) => t.id === d.base);
    this.selected = { kind: 'tx', index: Math.max(base, 0) };
    this.selectedHap = null;
    this.exprCard = null; // built when "More evidence" is first opened
    this.structures = {};
    this.reportHost = h('div', { class: 'report-host' });
    this.loadReport();
    this.render();
  }

  /** The gene report (baseline + haplotypes) for the chosen tissue. */
  loadReport() {
    if (!this.reportHost) return;
    reportCards(this.reportHost, ds(), this.params(), this.cfg.tissue, this.latest.controller && this.latest.controller.signal, (txId, hap) => {
      const index = this.data.transcripts.findIndex((t) => t.id === txId);
      if (index < 0) return;
      this.selected = { kind: 'tx', index };
      this.selectedHap = hap;
      this.render();
      const card = this.structures.protein && this.structures.protein.card.el;
      if (card) card.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }, this.reportState, (report) => {
      if (this.selectedHap) return;
      // open the protein view on the most interesting copy
      const pick = ['H1', 'H2'].find((hap) => { const f = report.haplotypes[hap].flags; return f.premature_stop.length || f.unstable.length || f.missense.length; });
      if (pick) { this.selectedHap = pick; this.render(); }
    });
  }

  render() {
    clear(this.content);
    const evidence = h('details', { class: 'card evidence', open: this.evidenceOpen ? true : null },
      h('summary', { class: 'card-head' }, h('h2', null, 'More evidence'), h('span', { class: 'muted small' }, 'mRNA and protein in every tissue · all transcripts and protein alignment · mRNA folding · AlphaGenome splicing · every variant')));
    const fill = () => {
      if (!this.exprCard) {
        const extra = DEFAULT_TISSUES.includes(this.cfg.tissue) ? {} : { tissues: this.cfg.tissue };
        this.exprCard = expressionCard(ds(), this.params(extra), this.data.target, this.latest.controller && this.latest.controller.signal);
      }
      evidence.append(h('div', { class: 'evidence-body' }, this.summaryCard(), this.exprCard, this.transcriptsCard(), (this.detailCard = h('section', { class: 'card' })),
        this.structureCard('rna'), this.splicingCard(), this.variantsCard()));
      this.renderDetail();
    };
    evidence.addEventListener('toggle', () => {
      this.evidenceOpen = evidence.open;
      if (evidence.open && !evidence.querySelector('.evidence-body')) fill();
    });
    this.content.append(this.reportHost, this.structureCard('protein'), evidence);
    if (this.evidenceOpen) fill();
  }

  // ------------------------------------------------------------------ summary
  summaryCard() {
    const d = this.data;
    const sp = d.splicing || {};
    const base = d.transcripts.find((t) => t.id === d.base);
    const coding = base && base.products.ref.coding;
    const expressed = sp.available ? sp.tracks.map((t, i) => h('span', { class: `pill ${sp.expressed[i] ? 'expressed' : 'muted'}`, title: `strongest annotated junction ${fmtNum(sp.expression[i], 3)}` },
      h('span', { class: `status-dot ${sp.expressed[i] ? 'good' : ''}` }), t.label)) : [h('span', { class: 'muted' }, sp.message || 'No splice predictions')];
    return h('section', { class: 'card' },
      h('div', { class: 'card-head' }, h('h2', null, `${d.target} · ${d.sample}`),
        h('span', { class: 'muted mono' }, `${d.chromosome}:${fmtInt(d.window_start + Math.min(...d.transcripts.map((t) => t.start)))} (${d.strand}) · ${d.gene_type || ''}`)),
      h('div', { class: 'card-body' },
        h('dl', { class: 'kv' },
          h('dt', null, 'Base transcript'), h('dd', null, `${d.base_name} (${d.base})${base && base.mane ? ' · MANE Select' : ''}${coding ? ` · ${fmtInt(base.products.ref.protein.length)} aa` : ' · non-coding'}`),
          ...HAPS.flatMap((hap) => [h('dt', null, hap), h('dd', null, base ? h('div', { class: 'product-cell' }, changeBadge(base.products[hap].change),
            h('span', { class: 'mono' }, base.products[hap].change.hgvs || ''), base.products[hap].nmd && base.products[hap].nmd.predicted ? h('span', { class: 'pill product-class bad' }, 'NMD') : null) : '–')]),
          h('dt', null, 'Spliced in'), h('dd', null, h('div', { class: 'product-cell wrap' }, ...expressed),
            base && base.exons.length < 2 ? h('div', { class: 'muted' }, `${d.base_name} has a single exon, so junction signal does not measure its expression; see RNA-seq/CAGE on Tracks.`) : null),
          d.truncated && d.truncated.length ? h('dt', null, 'Not shown') : null,
          d.truncated && d.truncated.length ? h('dd', { class: 'muted' }, `${d.truncated.length} transcript(s) extend past the window: ${d.truncated.slice(0, 6).join(', ')}`) : null)));
  }

  // ------------------------------------------------------------------ transcripts
  trackSelect() {
    const sp = this.data.splicing || {};
    if (!sp.available) return null;
    return field('Tissue track', select(sp.tracks.map((t) => ({ value: t.index, label: `${t.label}${sp.expressed[t.index] ? '' : ' (not expressed)'}` })), this.cfg.track,
      (v) => { this.cfg.track = Number(v); this.persist(); this.render(); }));
  }

  supportText(txId) {
    const s = ((this.data.splicing || {}).support || {})[txId];
    if (!s) return '–';
    const t = this.cfg.track;
    return ALL.filter((k) => s[k]).map((k) => `${k} ${fmtNum(s[k][t], 2)}`).join(' · ');
  }

  transcriptsCard() {
    const d = this.data;
    const rows = d.transcripts.map((t, i) => {
      const ref = t.products.ref;
      const sel = this.selected && this.selected.kind === 'tx' && this.selected.index === i;
      return h('tr', { class: `clickable${sel ? ' selected' : ''}`, onclick: () => { this.selected = { kind: 'tx', index: i }; this.render(); } },
        h('td', null, transcriptLink(t.id, t.name), t.mane ? h('span', { class: 'pill', style: { marginLeft: '6px' } }, 'MANE') : null,
          t.incomplete && (t.incomplete.start || t.incomplete.end) ? h('span', { class: 'pill muted', style: { marginLeft: '6px' }, title: 'GENCODE marks this model incomplete' }, 'partial') : null),
        h('td', { class: 'muted' }, (t.type || '').replace(/_/g, ' ')),
        h('td', { class: 'num' }, fmtInt(t.exons.length)),
        h('td', { class: 'num' }, fmtInt(ref.mrna_length)),
        h('td', { class: 'num' }, ref.coding ? fmtInt((ref.protein || '').length) : '–'),
        changeCell(t.products.H1), changeCell(t.products.H2),
        h('td', { class: 'mono muted', title: 'Weakest junction usage of the transcript (AlphaGenome), on the selected track' }, this.supportText(t.id)));
    });
    return h('section', { class: 'card' },
      h('div', { class: 'card-head' }, h('h2', null, 'Annotated transcripts'), this.trackSelect()),
      h('div', { class: 'table-wrap', style: { maxHeight: '46vh' } }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['Transcript', 'Type', 'Exons', 'mRNA nt', 'Protein aa', 'H1 vs reference', 'H2 vs reference', 'Junction support'].map((c) => h('th', null, c)))),
        h('tbody', null, rows))));
  }

  // ------------------------------------------------------------------ structure (3D protein, mRNA folds)
  structureCard(part) {
    const item = this.selectedItem();
    if (!item) return h('div');
    const tx = this.selected.kind === 'tx' ? this.data.transcripts[this.selected.index].id : `cand:${this.selected.index}`;
    const changed = ['H1', 'H2'].find((hap) => item.products[hap] && !['no_change', 'utr_change', 'synonymous', 'noncoding_change'].includes((item.products[hap].change || {}).class));
    const hap = this.selectedHap || changed || 'H1';
    const key = `${tx}|${hap}`;
    const hit = this.structures[part];
    if (hit && hit.key === key) return hit.card.el;
    const card = new StructureCard({
      dsPath: ds(), params: this.params(), item: tx, title: item.title, exons: item.exons, payload: this.data, defaultHap: hap, parts: [part],
      signal: this.latest.controller && this.latest.controller.signal,
    });
    this.structures[part] = { key, card };
    return card.el;
  }

  // ------------------------------------------------------------------ detail (protein alignment)
  selectedItem() {
    if (!this.selected) return null;
    if (this.selected.kind === 'tx') {
      const t = this.data.transcripts[this.selected.index];
      return t && { title: t.name, products: t.products, exons: t.exons, reference: t.products.ref };
    }
    const c = (this.data.splicing.candidates || [])[this.selected.index];
    const base = this.data.transcripts.find((t) => t.id === this.data.base);
    return c && { title: `${this.data.base_name} with ${c.kind}`, products: c.products, exons: c.exons, reference: base.products.ref, candidate: c };
  }

  renderDetail() {
    const host = clear(this.detailCard);
    const item = this.selectedItem();
    if (!item) return;
    const lines = [];
    for (const name of ALL) {
      const p = item.products[name];
      if (!p) continue;
      const bits = [`${fmtInt(p.mrna_length)} nt`, `${fmtInt(p.exon_count)} exon${p.exon_count === 1 ? '' : 's'}`];
      if (p.coding) bits.push(`${fmtInt((p.protein || '').length)} aa`, p.stop_found ? `5′UTR ${fmtInt(p.utr5)} · 3′UTR ${fmtInt(p.utr3)}` : 'no stop');
      if (p.nmd) bits.push(p.nmd.predicted ? `NMD: ${p.nmd.reason}` : `no NMD (${p.nmd.reason})`);
      lines.push(h('dt', null, name === 'ref' ? 'Reference' : name), h('dd', null, bits.join(' · '), (p.notes || []).length ? h('div', { class: 'muted' }, p.notes.join('; ')) : null));
    }
    host.append(
      h('div', { class: 'card-head' }, h('h2', null, item.title), h('span', { class: 'muted' }, item.candidate ? 'Candidate isoform: compared with the base transcript’s reference protein' : 'Compared with this transcript’s reference protein')),
      h('div', { class: 'card-body' }, h('dl', { class: 'kv', style: { marginBottom: '12px' } }, ...lines), this.alignment(item)));
  }

  alignment(item) {
    const ref = item.reference.protein || '';
    const seqs = ALL.map((name) => ({ name, seq: (item.products[name] || {}).protein || '' })).filter((s) => item.products[s.name]);
    if (!seqs.some((s) => s.seq)) return h('p', { class: 'muted' }, 'Non-coding transcript: no protein. Download the mRNA FASTA for the sequences.');
    const width = 60;
    const total = Math.max(...seqs.map((s) => s.seq.length), ref.length);
    const blocks = [];
    const showAll = this.showAllBlocks;
    for (let start = 0; start < total; start += width) {
      const differs = seqs.some((s) => s.seq.slice(start, start + width) !== ref.slice(start, start + width));
      if (!differs && !showAll) continue;
      blocks.push(h('div', { class: 'aln-block' },
        ...seqs.map((s) => {
          const row = h('div', { class: 'aln-row' }, h('span', { class: 'aln-name' }, s.name === 'ref' ? (item.candidate ? 'base ref' : 'ref') : s.name), h('span', { class: 'aln-pos' }, fmtInt(start + 1)));
          const chunk = s.seq.slice(start, start + width);
          const seqEl = h('span', { class: 'aln-seq' });
          let run = '';
          let runDiff = null;
          const flush = () => { if (run) seqEl.appendChild(runDiff ? h('mark', null, run) : document.createTextNode(run)); run = ''; };
          for (let i = 0; i < chunk.length; i++) {
            const diff = s.name !== 'ref' && chunk[i] !== ref[start + i];
            if (runDiff !== diff) { flush(); runDiff = diff; }
            run += chunk[i];
          }
          flush();
          if (s.seq.length > start && s.seq.length <= start + width && (item.products[s.name] || {}).stop_found) seqEl.appendChild(h('span', { class: 'aln-stop' }, '*'));
          row.appendChild(seqEl);
          return row;
        })));
    }
    const toggle = checkbox('Show identical blocks', !!showAll, (v) => { this.showAllBlocks = v; this.renderDetail(); });
    if (!blocks.length) return h('div', null, toggle, h('p', { class: 'muted' }, 'The three proteins are identical.'));
    return h('div', null, toggle, h('div', { class: 'aln' }, ...blocks));
  }

  // ------------------------------------------------------------------ splicing
  splicingCard() {
    const sp = this.data.splicing || {};
    const card = h('section', { class: 'card' }, h('div', { class: 'card-head' }, h('h2', null, 'Splicing (AlphaGenome junctions)'), this.trackSelect()));
    const body = h('div', { class: 'card-body' });
    card.appendChild(body);
    if (!sp.available) {
      body.appendChild(h('div', { class: 'empty' }, sp.message || 'No splice-junction predictions.'));
      return card;
    }
    const t = this.cfg.track;
    if (sp.missing && sp.missing.length) body.appendChild(h('p', { class: 'muted' }, `No splice_junctions prediction for ${sp.missing.join(', ')}; run AlphaGenome → Predict with SPLICE_JUNCTIONS for this sample.`));
    if (!sp.expressed[t]) body.appendChild(h('p', { class: 'warn-text' }, `The gene is barely spliced in ${sp.tracks[t].label} (strongest annotated junction ${fmtNum(sp.expression[t], 3)}); usage values there are noise.`));
    body.appendChild(this.sashimi(sp, t));
    const usageCells = (j) => ALL.map((k) => j[k] ? h('td', { class: 'num' }, fmtNum(j[k].usage[t], 2)) : h('td', { class: 'num muted' }, '–'));
    const deltaCell = (j) => {
      const ds_ = HAPS.filter((k) => j[k]).map((k) => j[k].delta[t]);
      const worst = ds_.length ? ds_.reduce((a, b) => (Math.abs(b) > Math.abs(a) ? b : a)) : 0;
      return h('td', { class: `num ${Math.abs(worst) >= 0.1 ? 'delta-strong' : ''}` }, ds_.length ? (worst > 0 ? '+' : '') + fmtNum(worst, 2) : '–');
    };
    const ws = this.data.window_start;
    const coord = (j) => `${fmtInt(ws + j.start)}–${fmtInt(ws + j.end - 1)}`;
    body.append(h('h3', { class: 'section-title' }, `Introns of ${this.data.base_name}`),
      h('div', { class: 'table-wrap', style: { maxHeight: '32vh' } }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['#', 'Intron (genomic)', 'Usage ref', 'H1', 'H2', 'Largest Δ'].map((c) => h('th', null, c)))),
        h('tbody', null, sp.base_introns.map((j, i) => h('tr', null, h('td', null, String(this.data.strand === '-' ? sp.base_introns.length - i : i + 1)), h('td', { class: 'mono' }, coord(j)), ...usageCells(j), deltaCell(j)))))));
    const cands = sp.candidates || [];
    body.append(h('h3', { class: 'section-title' }, 'Candidate isoforms'),
      h('p', { class: 'muted' }, `Junctions in the gene that are not introns of the base transcript, used ≥ ${sp.thresholds.usage} or changed by ≥ ${sp.thresholds.delta} on an expressed track, each spliced into ${this.data.base_name}; and introns whose junction drops to ≤ ${sp.thresholds.retention_ratio}× the reference's. Click one for its proteins.`));
    if (!cands.length) body.appendChild(h('div', { class: 'empty' }, 'No candidate passes the thresholds.'));
    else body.appendChild(h('div', { class: 'table-wrap', style: { maxHeight: '40vh' } }, h('table', { class: 'table compact' },
      h('thead', null, h('tr', null, ['Event', 'Junction (genomic)', 'Annotated in', 'Usage ref', 'H1', 'H2', 'Largest Δ', 'Isoform on ref', 'Isoform on H1', 'Isoform on H2'].map((c) => h('th', null, c)))),
      h('tbody', null, cands.map((c, i) => {
        const sel = this.selected && this.selected.kind === 'cand' && this.selected.index === i;
        return h('tr', { class: `clickable${sel ? ' selected' : ''}`, onclick: () => { this.selected = { kind: 'cand', index: i }; this.render(); this.detailCard.scrollIntoView({ behavior: 'smooth', block: 'start' }); } },
          h('td', null, c.kind, c.note ? h('div', { class: 'muted', style: { fontSize: '11px' } }, c.note) : null),
          h('td', { class: 'mono' }, coord(c.junction)),
          h('td', { class: 'muted' }, (c.junction.annotated || []).join(', ') || 'novel'),
          ...usageCells(c.junction), deltaCell(c.junction),
          changeCell(c.products.ref), changeCell(c.products.H1), changeCell(c.products.H2));
      })))));
    return card;
  }

  /** Base transcript exons on a genomic axis, its introns as arcs above (ref usage) and candidate
   *  junctions below (the haplotype with the largest change), arc weight = usage on the track. */
  sashimi(sp, t) {
    const d = this.data;
    const base = d.transcripts.find((x) => x.id === d.base);
    const juncs = [...sp.base_introns.map((j) => ({ j, cand: false })), ...(sp.candidates || []).map((c) => ({ j: c.junction, cand: true, kind: c.kind }))];
    const lo = Math.min(base.start, ...juncs.map((x) => x.j.start));
    const hi = Math.max(base.end, ...juncs.map((x) => x.j.end));
    const W = 1000;
    const H = 170;
    const mid = 85;
    const x = (p) => 10 + ((p - lo) / Math.max(hi - lo, 1)) * (W - 20);
    const ns = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(ns, 'svg');
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.setAttribute('class', 'sashimi');
    const add = (tag, attrs, text) => { const el = document.createElementNS(ns, tag); for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v); if (text) el.textContent = text; svg.appendChild(el); return el; };
    add('line', { x1: x(base.start), x2: x(base.end), y1: mid, y2: mid, class: 'sashimi-axis' });
    for (const [s, e] of base.exons) add('rect', { x: x(s), y: mid - 7, width: Math.max(x(e) - x(s), 1.5), height: 14, class: 'sashimi-exon' });
    for (const { j, cand, kind } of juncs) {
      let pick = 'ref';
      for (const k of HAPS) if (j[k] && Math.abs(j[k].delta[t]) > (pick === 'ref' ? 0 : Math.abs(j[pick].delta[t]))) pick = k;
      const u = (j[pick] || j.ref).usage[t];
      const x1 = x(j.start);
      const x2 = x(j.end);
      const up = !cand;
      const hgt = Math.min(60, 14 + (x2 - x1) * 0.25);
      const y0 = up ? mid - 8 : mid + 8;
      const yc = up ? y0 - hgt : y0 + hgt;
      const path = add('path', { d: `M${x1},${y0} C${x1},${yc} ${x2},${yc} ${x2},${y0}`, class: `sashimi-arc${cand ? ' cand' : ''}`, 'stroke-width': String(0.6 + 5 * Math.max(0, Math.min(1, u))) });
      const title = document.createElementNS(ns, 'title');
      title.textContent = `${cand ? kind : 'base intron'} ${fmtInt(d.window_start + j.start)}–${fmtInt(d.window_start + j.end - 1)} · usage ref ${fmtNum(j.ref.usage[t], 2)}${HAPS.filter((k) => j[k]).map((k) => ` · ${k} ${fmtNum(j[k].usage[t], 2)} (Δ ${fmtNum(j[k].delta[t], 2)})`).join('')}`;
      path.appendChild(title);
    }
    add('text', { x: 10, y: 14, class: 'sashimi-label' }, `${d.base_name} introns (usage; thicker = more used)`);
    add('text', { x: 10, y: H - 6, class: 'sashimi-label' }, `candidate junctions (haplotype with the largest change) · 5′→3′ is ${d.strand === '-' ? 'right to left' : 'left to right'}`);
    return h('div', { class: 'sashimi-wrap' }, svg);
  }

  // ------------------------------------------------------------------ variants
  variantsCard() {
    const d = this.data;
    const shown = d.variants.filter((v) => this.cfg.intronic || !QUIET_REGIONS.has(v.region));
    const toggle = checkbox(`Intronic and flanking too (${fmtInt(d.variants.length)})`, this.cfg.intronic, (v) => { this.cfg.intronic = v; this.persist(); this.render(); });
    const trim = (s) => (s && s.length > 12 ? `${s.slice(0, 10)}…(${s.length})` : s || '·');
    return h('section', { class: 'card' },
      h('div', { class: 'card-head' }, h('h2', null, `Variants carried in ${d.base_name}`), toggle),
      shown.length ? h('div', { class: 'table-wrap', style: { maxHeight: '36vh' } }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['Position', 'ID', 'REF', 'H1', 'H2', 'Region'].map((c) => h('th', null, c)))),
        h('tbody', null, shown.map((v) => h('tr', null,
          h('td', { class: 'mono' }, `${d.chromosome}:${fmtInt(v.pos)}`), h('td', { class: 'mono muted' }, v.id || '–'), h('td', { class: 'mono' }, trim(v.ref)),
          h('td', { class: 'mono' }, v.H1 ? trim(v.H1) : h('span', { class: 'muted' }, 'ref')), h('td', { class: 'mono' }, v.H2 ? trim(v.H2) : h('span', { class: 'muted' }, 'ref')),
          h('td', null, h('span', { class: `pill${QUIET_REGIONS.has(v.region) ? ' muted' : ''}` }, v.region)))))))
        : h('div', { class: 'card-body' }, h('div', { class: 'empty' }, 'No exonic or splice-region variant on either haplotype.')));
  }

  // ------------------------------------------------------------------ export
  async download(kind) {
    try {
      const res = await fetch(`${ds()}/products/fasta?${new URLSearchParams({ ...this.params(), kind, target: this.cfg.target || '' })}`);
      const text = await res.text();
      if (!res.ok) throw new Error(text);
      downloadText(`${slug(`${this.cfg.sample}_${this.data ? this.data.target : this.cfg.gene}_${kind}`)}.fa`, text);
    } catch (err) { toast(err.message, 'error', 8000); }
  }

  async downloadJson() {
    try {
      const data = await api(`${ds()}/products`, { params: this.params({ sequences: 'mrna' }) });
      downloadText(`${slug(`${this.cfg.sample}_${data.target}_products`)}.json`, JSON.stringify(data, null, 1), 'application/json');
    } catch (err) { toast(err.message, 'error', 8000); }
  }
}
