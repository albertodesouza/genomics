// Structure card of the Gene products page: the protein in 3D (AlphaFold DB model of the reference,
// the haplotype's changes placed on it, or ESMFold when the change needs a new prediction) and the
// mRNA's secondary structure (ViennaRNA) around the start codon, the stop codon or a variant.
import { api, apiJob, isAbort } from './api.js';
import { h, clear, fmtInt, fmtNum, errorBox, jobOverlay, segmented, select, field, toast } from './ui.js';

const MOL_URL = 'https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js';
let molPromise = null;

/** 3Dmol.js, loaded from cdnjs the first time a structure is shown. */
function load3Dmol() {
  if (window.$3Dmol) return Promise.resolve(window.$3Dmol);
  if (!molPromise) {
    molPromise = new Promise((resolve, reject) => {
      const s = document.createElement('script');
      s.src = MOL_URL;
      s.onload = () => resolve(window.$3Dmol);
      s.onerror = () => { molPromise = null; reject(new Error('Could not load the 3D viewer (3Dmol.js from cdnjs.cloudflare.com); is this browser offline?')); };
      document.head.appendChild(s);
    });
  }
  return molPromise;
}

// AlphaFold's pLDDT colours (very high, confident, low, very low)
function plddtColor(v) {
  if (v >= 90) return '#0053D6';
  if (v >= 70) return '#65CBF3';
  if (v >= 50) return '#FFDB13';
  return '#FF7D45';
}

const cssVar = (name, fallback) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;

/**
 * One 3Dmol viewer. ``opts``: {pdb, scale (pLDDT multiplier), offset (model residue = ours + offset),
 * marks: [{pos, label}], lostFrom, changedFrom}.
 */
async function molViewer(host, opts) {
  const $3Dmol = await load3Dmol();
  const el = h('div', { class: 'mol-view' });
  host.appendChild(el);
  const viewer = $3Dmol.createViewer(el, { backgroundColor: cssVar('--surface', '#ffffff'), antialias: true });
  viewer.addModel(opts.pdb, 'pdb');
  const scale = opts.scale || 1;
  const offset = opts.offset || 0;
  viewer.setStyle({}, { cartoon: { colorfunc: (atom) => plddtColor(atom.b * scale) } });
  if (opts.lostFrom) viewer.setStyle({ resi: `${opts.lostFrom + offset}-99999` }, { cartoon: { color: '#bdbdbd', style: 'trace', thickness: 0.25 } });
  if (opts.changedFrom) viewer.addStyle({ resi: `${opts.changedFrom + offset}-99999` }, { stick: { color: '#c026d3', radius: 0.18 } });
  for (const m of opts.marks || []) {
    const resi = m.pos + offset;
    viewer.addStyle({ resi }, { sphere: { color: m.color || '#c026d3', radius: 1.2 } });
    const atom = viewer.selectedAtoms({ resi, atom: 'CA' })[0];
    if (atom && m.label) viewer.addLabel(m.label, { position: { x: atom.x, y: atom.y, z: atom.z }, backgroundColor: '#111111', backgroundOpacity: 0.75, fontColor: '#ffffff', fontSize: 12, inFront: true });
  }
  viewer.zoomTo();
  viewer.render();
  // the canvas is sized once the grid has laid out: fit again then
  requestAnimationFrame(() => { viewer.resize(); viewer.zoomTo(); viewer.zoom(1.25); viewer.render(); });
  return viewer;
}

const legend = () => h('div', { class: 'mol-legend' }, 'pLDDT: ',
  ...[['#0053D6', '> 90'], ['#65CBF3', '70–90'], ['#FFDB13', '50–70'], ['#FF7D45', '< 50']].map(([c, t]) => h('span', { class: 'mol-swatch' }, h('i', { style: { background: c } }), t)),
  h('span', { class: 'mol-swatch' }, h('i', { style: { background: '#c026d3' } }), 'changed'),
  h('span', { class: 'mol-swatch' }, h('i', { style: { background: '#bdbdbd' } }), 'lost'));

export class StructureCard {
  /** ``ctx``: {dsPath, params: {gene, sample, target}, item: tx id | 'cand:i', title, payload, signal}. */
  constructor(ctx) {
    this.ctx = ctx;
    this.hap = ctx.defaultHap || 'H1';
    this.viewers = [];
    this.el = h('section', { class: 'card' });
    this.render();
  }

  render() {
    const c = this.ctx;
    const parts = c.parts || ['protein', 'rna'];
    clear(this.el);
    this.protein = h('div', { class: 'structure-protein', style: { position: 'relative', minHeight: '120px' } });
    this.rna = h('div', { class: 'structure-rna', style: { position: 'relative', minHeight: '120px' } });
    const only = parts.length === 1 ? parts[0] : null;
    const title = only === 'protein' ? h('h2', null, h('span', { class: 'step' }, '3'), ` Protein shape: ${c.title} on ${this.hap}`)
      : only === 'rna' ? h('h2', null, `mRNA folding: ${c.title}`) : h('h2', null, `Structure · ${c.title}`);
    const body = h('div', { class: 'card-body' });
    if (parts.includes('protein')) body.append(only ? '' : h('h3', { class: 'section-title' }, 'Protein (3D)'), this.protein);
    if (parts.includes('rna')) body.append(only ? '' : h('h3', { class: 'section-title' }, 'mRNA (secondary structure)'), this.rna);
    this.el.append(
      h('div', { class: 'card-head' }, title,
        field('Copy', segmented([{ value: 'H1', label: 'H1' }, { value: 'H2', label: 'H2' }], this.hap, (v) => { this.hap = v; this.render(); }))),
      body);
    if (parts.includes('protein')) this.loadProtein();
    if (parts.includes('rna')) this.loadRna();
  }

  // ------------------------------------------------------------------ protein
  async loadProtein() {
    const c = this.ctx;
    const host = this.protein;
    const overlay = jobOverlay(host, 'AlphaFold DB model');
    let data;
    try {
      data = await apiJob(`${c.dsPath}/products/structure`, { params: { ...c.params, tx: c.item }, signal: c.signal, onProgress: (job) => overlay.update(job) });
    } catch (err) {
      overlay.remove();
      if (!isAbort(err)) host.appendChild(errorBox(err));
      return;
    }
    overlay.remove();
    if (!data.available) { host.appendChild(h('div', { class: 'empty' }, data.message)); return; }
    this.structure = data;
    const product = data.products[this.hap];
    const model = data.model || {};
    const usable = model.pdb && model.match !== 'different';
    const info = h('div', { class: 'structure-info' });
    host.appendChild(info);
    const modelText = usable
      ? h('span', null, 'Reference: ', h('a', { href: model.url, target: '_blank', rel: 'noopener' }, `AlphaFold DB ${model.accession}`), model.protein_name ? ` (${model.protein_name})` : '',
        ` · mean pLDDT ${fmtNum(model.plddt_mean, 1)}`, model.match !== 'identical' ? ` · UniProt sequence ${model.match} (numbering offset ${model.offset})` : '')
      : h('span', null, model.error ? `No AlphaFold DB model: ${model.error}.` : `The AlphaFold DB model (${model.accession}) is a different isoform than this transcript's protein.`);
    info.appendChild(h('p', null, modelText));
    const mode = product.mode;
    let describe = '';
    if (mode === 'identical') describe = `The ${this.hap} protein is identical to the reference.`;
    else if (mode === 'substitution') describe = `${this.hap} differs at ${product.substitutions.length} residue${product.substitutions.length === 1 ? '' : 's'}: ${product.substitutions.slice(0, 8).map((s) => `${s.ref}${s.pos}${s.alt}`).join(', ')}. Shown on the reference model: a predicted fold does not resolve a substitution's effect on stability (use AlphaMissense or a ΔΔG predictor for that).`;
    else if (mode === 'truncated') describe = `${this.hap} stops early: ${fmtInt(product.length)} of ${fmtInt(data.reference_length)} residues; the missing part (from residue ${product.lost_from}) is drawn as a grey trace.${product.nmd && product.nmd.predicted ? ' Its mRNA is predicted to undergo NMD, so little of this protein should be made.' : ''}`;
    else describe = `${this.hap}'s protein changes from residue ${product.first_change} on (${(product.change && product.change.hgvs) || 'altered sequence'}); the reference model cannot show it. Fold both with ESMFold to compare.`;
    info.appendChild(h('p', { class: mode === 'identical' ? 'muted' : '' }, describe));
    const sim = product.similarity;
    if (sim && sim.similarity !== null && sim.similarity !== undefined) {
      const tone = sim.similarity >= 0.999 ? 'ok' : sim.similarity >= 0.95 ? 'warn' : 'bad';
      info.appendChild(h('div', { class: `verdict similarity ${tone}` },
        h('div', { class: 'similarity-score' }, `${(sim.similarity * 100).toFixed(sim.similarity >= 0.99 ? 2 : 1)}%`),
        h('div', null, h('div', { class: 'verdict-title' }, 'Sequence similarity to the reference protein'),
          h('div', { class: 'verdict-detail' }, `BLOSUM62 alignment score / reference self-score · identity ${(sim.identity * 100).toFixed(2)}% · ${fmtInt(sim.length)} of ${fmtInt(sim.reference_length)} residues${mode === 'altered' || !usable ? ' · structural similarity (TM-score) after ESMFold below' : ''}`))));
    }
    if (mode === 'substitution') {
      info.appendChild(h('div', { class: 'table-wrap' }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['Residue', 'Model confidence (pLDDT)', 'Position in the fold', 'AlphaMissense'].map((x) => h('th', null, x)))),
        h('tbody', null, product.substitutions.map((x) => h('tr', null,
          h('td', { class: 'mono' }, h('b', null, `${x.ref}${x.pos}${x.alt}`)),
          h('td', null, x.plddt === undefined ? '–' : `${fmtNum(x.plddt, 3)}${x.plddt < 50 ? ' (disordered region)' : x.plddt < 70 ? ' (low confidence)' : ''}`),
          h('td', null, x.burial ? `${x.burial} (${x.neighbours} residues within 10 Å)` : '–'),
          h('td', null, x.alphamissense ? h('span', { class: `pill product-class ${x.alphamissense.class === 'likely pathogenic' ? 'bad' : x.alphamissense.class === 'ambiguous' ? 'warn' : 'good'}` }, `${x.alphamissense.score.toFixed(2)} · ${x.alphamissense.class}`) : '–')))))),
        h('p', { class: 'muted small' }, 'A buried residue in a confident region is more likely to matter for the fold than a surface one in a disordered loop. AlphaMissense estimates pathogenicity, not function (e.g. MC1R R151C, a known loss-of-function red-hair allele, scores likely benign).'));
    }
    const needsFold = mode === 'altered' || !usable;
    if (needsFold) {
      const remote = data.esmfold && data.esmfold.remote;
      const btn = h('button', { class: 'btn primary', disabled: !remote, onclick: () => this.fold(btn) }, `Predict with ESMFold (≤ ${data.esmfold.max} residues)`);
      info.appendChild(h('div', { class: 'structure-fold' }, btn,
        h('span', { class: 'muted small' }, remote ? 'Sends the reference and this haplotype\'s protein sequences to api.esmatlas.com (Meta\'s public ESMFold service). Results are cached.' : 'Unavailable: remote lookups are disabled (--no-remote).')));
    }
    if (!usable) return;
    const grid = h('div', { class: 'mol-grid' });
    const leftHost = h('div', { class: 'mol-cell' }, h('div', { class: 'mol-title' }, 'Reference genome'));
    const rightHost = h('div', { class: 'mol-cell' }, h('div', { class: 'mol-title' }, `${this.hap} copy`));
    grid.append(leftHost, rightHost);
    host.append(grid, legend());
    const marks = mode === 'substitution' ? product.substitutions.map((s) => ({ pos: s.pos, label: `${s.ref}${s.pos}${s.alt}` })) : [];
    try {
      const refMarks = mode === 'substitution' ? product.substitutions.map((s) => ({ pos: s.pos, label: `${s.ref}${s.pos}`, color: '#6b7280' })) : [];
      const left = await molViewer(leftHost, { pdb: model.pdb, scale: 100 / (model.plddt_scale || 100), offset: model.offset || 0, marks: refMarks });
      const right = await molViewer(rightHost, {
        pdb: model.pdb, scale: 100 / (model.plddt_scale || 100), offset: model.offset || 0, marks,
        lostFrom: mode === 'truncated' ? product.lost_from : null,
      });
      if (mode === 'altered') rightHost.appendChild(h('div', { class: 'mol-overlay muted' }, 'Reference model shown; this product needs ESMFold.'));
      if (left.linkViewer) { left.linkViewer(right); right.linkViewer(left); }
      this.viewers = [left, right];
    } catch (err) {
      grid.replaceWith(errorBox(err));
    }
  }

  async fold(btn) {
    const c = this.ctx;
    btn.disabled = true;
    const host = this.protein;
    const overlay = jobOverlay(host, 'ESMFold (api.esmatlas.com)');
    let data;
    try {
      data = await apiJob(`${c.dsPath}/products/fold`, { method: 'POST', params: { ...c.params, tx: c.item, hap: this.hap }, signal: c.signal, onProgress: (job) => overlay.update(job) });
    } catch (err) {
      overlay.remove();
      btn.disabled = false;
      if (!isAbort(err)) toast(err.message, 'error', 9000);
      return;
    }
    overlay.remove();
    const grid = h('div', { class: 'mol-grid' });
    const w = data.window;
    const leftHost = h('div', { class: 'mol-cell' }, h('div', { class: 'mol-title' }, `Reference, ESMFold of residues ${w.reference[0]}–${w.reference[1]} · pLDDT ${Math.round(data.reference.plddt_mean * 100)}`));
    const rightHost = h('div', { class: 'mol-cell' }, h('div', { class: 'mol-title' }, `${this.hap}, ESMFold of residues ${w.product[0]}–${w.product[1]} · pLDDT ${Math.round(data.product.plddt_mean * 100)}`));
    grid.append(leftHost, rightHost);
    const ss = data.structure_similarity || {};
    const tmTone = ss.tm === null || ss.tm === undefined ? '' : ss.tm >= 0.9 ? 'ok' : ss.tm >= 0.5 ? 'warn' : 'bad';
    host.append(h('h4', { class: 'structure-sub' }, 'ESMFold, same window for both'),
      h('div', { class: `verdict similarity ${tmTone}` },
        h('div', { class: 'similarity-score' }, ss.tm === null || ss.tm === undefined ? '–' : ss.tm.toFixed(2)),
        h('div', null, h('div', { class: 'verdict-title' }, 'Structural similarity (TM-score) to the reference fold'),
          h('div', { class: 'verdict-detail' }, ss.tm === null || ss.tm === undefined ? (ss.error || 'not computable') : `1 = same fold, > 0.5 = same fold family, < 0.3 = unrelated · RMSD ${fmtNum(ss.rmsd, 3)} Å over ${fmtInt(ss.aligned)} aligned residues · normalised by the reference window`))),
      grid, legend(),
      h('p', { class: 'muted small' }, `Residues from ${data.first_change} on differ (magenta sticks on the right). ESMFold is single-sequence: its confidence (pLDDT) for a novel frameshift tail is usually low, which itself says the tail is unlikely to fold.`));
    try {
      const left = await molViewer(leftHost, { pdb: data.reference.pdb, scale: 100, offset: -(w.reference[0] - 1) });
      const right = await molViewer(rightHost, { pdb: data.product.pdb, scale: 100, offset: -(w.product[0] - 1), changedFrom: data.first_change });
      this.viewers.push(left, right);
    } catch (err) {
      grid.replaceWith(errorBox(err));
    }
  }

  // ------------------------------------------------------------------ RNA
  /** Start / stop codon, and every variant this haplotype carries inside the product's own exons. */
  windowChoices() {
    const p = this.ctx.payload;
    const exons = this.ctx.exons || [];
    const inExon = (v) => exons.some(([s, e]) => v.offset + Math.max(v.ref.length, 1) > s && v.offset < e);
    const choices = [{ value: 'start', label: 'Start codon (translation initiation)' }, { value: 'stop', label: 'Stop codon' }];
    for (const v of p.variants || []) {
      if (!v[this.hap] || !inExon(v)) continue;
      const short = (x) => (x.length > 6 ? `${x.slice(0, 6)}…` : x);
      choices.push({ value: String(v.offset), label: `${p.chromosome}:${fmtInt(v.pos)} ${short(v.ref)}>${short(v[this.hap])} (exonic in ${this.ctx.title})` });
    }
    return choices;
  }

  async loadRna(center = null, width = null) {
    const c = this.ctx;
    const host = clear(this.rna);
    const choices = this.windowChoices();
    this.rnaCenter = center || this.rnaCenter || (choices[2] ? choices[2].value : 'start');
    if (!choices.some((x) => x.value === this.rnaCenter)) this.rnaCenter = 'start';
    this.rnaWidth = width || this.rnaWidth || 240;
    host.appendChild(h('div', { class: 'toolbar', style: { marginBottom: '8px' } },
      field('Window around', select(choices, this.rnaCenter, (v) => this.loadRna(v))),
      field('Width', select([120, 240, 400, 600].map((n) => ({ value: n, label: `${n} nt` })), this.rnaWidth, (v) => this.loadRna(null, Number(v))))));
    host.appendChild(h('p', { class: 'muted small' }, 'Why not 3D: a full mRNA has no reliable 3D structure prediction (3D RNA methods handle short, isolated RNAs, and mRNA in cells is coated by proteins). The minimum-free-energy secondary structure of a window is the standard proxy; a more negative ΔG is a more stable fold, and a stable fold over the start codon hinders translation.'));
    const panel = h('div', { style: { position: 'relative', minHeight: '80px' } });
    host.appendChild(panel);
    let data;
    try {
      data = await api(`${c.dsPath}/products/rna`, { params: { ...c.params, tx: c.item, hap: this.hap, center: this.rnaCenter, width: this.rnaWidth }, signal: c.signal });
    } catch (err) {
      if (!isAbort(err)) panel.appendChild(errorBox(err));
      return;
    }
    if (!data.available) { panel.appendChild(h('div', { class: 'empty' }, data.message)); return; }
    const ddg = data.ddg;
    panel.appendChild(h('p', null,
      `Transcript positions ${fmtInt(data.reference.start)}–${fmtInt(data.reference.end)} (reference) and ${fmtInt(data.haplotype.start)}–${fmtInt(data.haplotype.end)} (${this.hap}). `,
      data.identical ? h('span', { class: 'muted' }, `Identical sequence on ${this.hap}: same fold.`)
        : h('span', null, `ΔG ${fmtNum(data.reference.mfe, 1)} → ${fmtNum(data.haplotype.mfe, 1)} kcal/mol: `,
          h('b', { class: ddg < -0.5 ? 'dir-up' : ddg > 0.5 ? 'dir-down' : '' }, `ΔΔG ${ddg > 0 ? '+' : ''}${fmtNum(ddg, 2)}`),
          ddg < -0.5 ? ' (more stable on this haplotype)' : ddg > 0.5 ? ' (less stable on this haplotype)' : ' (no meaningful change)')));
    const grid = h('div', { class: 'rna-grid' },
      h('div', null, h('div', { class: 'mol-title' }, `Reference · ${fmtNum(data.reference.mfe, 1)} kcal/mol`), rnaSvg(data.reference, [])),
      h('div', null, h('div', { class: 'mol-title' }, `${this.hap} · ${fmtNum(data.haplotype.mfe, 1)} kcal/mol`), rnaSvg(data.haplotype, data.haplotype.changed)));
    panel.append(grid, h('div', { class: 'mol-legend' }, ...['A', 'C', 'G', 'U'].map((b) => h('span', { class: 'mol-swatch' }, h('i', { class: `rna-base-${b}` }), b)),
      h('span', { class: 'mol-swatch' }, h('i', { class: 'rna-changed-swatch' }), 'differs from the reference')));
  }
}

/** 2D secondary-structure drawing (ViennaRNA layout): backbone, base pairs, bases, changed bases ringed. */
function rnaSvg(fold, changed) {
  const ns = 'http://www.w3.org/2000/svg';
  const xs = fold.xy.map((p) => p[0]);
  const ys = fold.xy.map((p) => p[1]);
  const pad = 12;
  const minX = Math.min(...xs) - pad;
  const minY = Math.min(...ys) - pad;
  const w = Math.max(...xs) - minX + pad;
  const hgt = Math.max(...ys) - minY + pad;
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', `${minX} ${minY} ${w} ${hgt}`);
  svg.setAttribute('class', 'rna-svg');
  const add = (tag, attrs, text) => { const el = document.createElementNS(ns, tag); for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v); if (text !== undefined) el.textContent = text; svg.appendChild(el); return el; };
  const r = Math.max(2.2, Math.min(5, Math.sqrt((w * hgt) / fold.xy.length) / 3.2));
  add('polyline', { points: fold.xy.map((p) => p.join(',')).join(' '), class: 'rna-backbone', 'stroke-width': String(r * 0.25) });
  for (const [i, j] of fold.pairs) add('line', { x1: fold.xy[i][0], y1: fold.xy[i][1], x2: fold.xy[j][0], y2: fold.xy[j][1], class: 'rna-pair', 'stroke-width': String(r * 0.35) });
  const marked = new Set(changed);
  fold.xy.forEach(([x, y], i) => {
    const base = fold.sequence[i];
    const c = add('circle', { cx: x, cy: y, r: marked.has(i) ? r * 1.6 : r, class: `rna-base-${base}${marked.has(i) ? ' rna-changed' : ''}` });
    const t = document.createElementNS(ns, 'title');
    t.textContent = `${base} · transcript position ${fold.start + i}${marked.has(i) ? ' · differs from the reference' : ''}`;
    c.appendChild(t);
  });
  add('text', { x: fold.xy[0][0], y: fold.xy[0][1] - r * 2.2, class: 'rna-end' }, "5′");
  add('text', { x: fold.xy[fold.xy.length - 1][0], y: fold.xy[fold.xy.length - 1][1] - r * 2.2, class: 'rna-end' }, "3′");
  return svg;
}
