// Expression card of the Gene products page: how much mRNA and protein the individual makes in each
// tissue, as a fold change against the reference genome (relative) and as an estimated level in
// TPM-like units (absolute), kept visibly apart because they rest on different evidence.
import { apiJob, isAbort } from './api.js';
import { h, clear, fmtNum, errorBox, jobOverlay } from './ui.js';

const HAPS = ['H1', 'H2'];
const LOG2_RANGE = 1; // bars span 0.5x .. 2x; longer bars are clipped and marked

const fmtFold = (v) => (v === null || v === undefined || !Number.isFinite(v) ? '–' : `${v >= 10 ? v.toFixed(1) : v.toFixed(2)}×`);
const fmtPctChange = (v) => {
  if (v === null || v === undefined || !Number.isFinite(v)) return '';
  const pct = (v - 1) * 100;
  return `${pct >= 0 ? '+' : '−'}${Math.abs(pct) < 10 ? Math.abs(pct).toFixed(1) : Math.round(Math.abs(pct))}%`;
};
const fmtAbs = (v, unit) => (v === null || v === undefined || !Number.isFinite(v) ? '–' : `${v >= 100 ? Math.round(v).toLocaleString('en-US') : v >= 10 ? v.toFixed(1) : v.toFixed(2)}${unit ? ` ${unit}` : ''}`);

function direction(v, similar) {
  if (v === null || v === undefined || !Number.isFinite(v)) return 'unknown';
  if (v >= 1 + similar) return 'up';
  if (v <= 1 / (1 + similar)) return 'down';
  return 'same';
}

const VERB = { up: 'more', down: 'less', same: 'about the same', unknown: '' };

/** Diverging log2 bars around 1× (the reference genome). */
function foldBars(rows, similar, muted) {
  const W = 420;
  const rowH = 22;
  const left = 132;
  const right = 64;
  const top = 16;
  const H = top + rows.length * rowH + 26;
  const plotW = W - left - right;
  const x = (fold) => left + plotW / 2 + (Math.max(-LOG2_RANGE, Math.min(LOG2_RANGE, Math.log2(fold))) / LOG2_RANGE) * (plotW / 2);
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('class', `fold-bars${muted ? ' muted' : ''}`);
  const add = (tag, attrs, text) => { const el = document.createElementNS(ns, tag); for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v); if (text !== undefined) el.textContent = text; svg.appendChild(el); return el; };
  // band of "about the same"
  add('rect', { x: x(1 / (1 + similar)), y: top, width: x(1 + similar) - x(1 / (1 + similar)), height: rows.length * rowH, class: 'fold-band' });
  add('text', { x: x(1), y: top - 5, class: 'fold-tick-label', 'text-anchor': 'middle' }, 'reference genome');
  add('text', { x: x(0.5) + 2, y: top - 5, class: 'fold-tick-label' }, '← less');
  add('text', { x: x(2) - 2, y: top - 5, class: 'fold-tick-label', 'text-anchor': 'end' }, 'more →');
  rows.forEach((r, i) => {
    const y = top + i * rowH;
    add('text', { x: left - 8, y: y + rowH / 2 + 4, class: `fold-label${r.strong ? ' strong' : ''}`, 'text-anchor': 'end' }, r.label);
    if (r.fold === null || r.fold === undefined || !Number.isFinite(r.fold)) {
      add('text', { x: x(1) + 6, y: y + rowH / 2 + 4, class: 'fold-value' }, 'n/a');
      return;
    }
    const x0 = x(1);
    const x1 = x(r.fold);
    const dir = direction(r.fold, similar);
    add('rect', { x: Math.min(x0, x1), y: y + 4, width: Math.max(Math.abs(x1 - x0), 1.5), height: rowH - 8, class: `fold-bar ${dir}${r.strong ? ' strong' : ''}` });
    if (Math.abs(Math.log2(r.fold)) > LOG2_RANGE) add('text', { x: x1 + (r.fold > 1 ? -10 : 2), y: y + rowH / 2 + 4, class: 'fold-clip' }, r.fold > 1 ? '»' : '«');
    add('text', { x: W - right + 6, y: y + rowH / 2 + 4, class: `fold-value${r.strong ? ' strong' : ''}` }, fmtFold(r.fold));
  });
  const axisY = top + rows.length * rowH;
  add('line', { x1: x(1), x2: x(1), y1: top, y2: axisY + 2, class: 'fold-axis' });
  for (const t of [0.5, 0.75, 1, 1.5, 2]) {
    add('line', { x1: x(t), x2: x(t), y1: axisY, y2: axisY + 4, class: 'fold-tick' });
    add('text', { x: x(t), y: axisY + 16, class: 'fold-tick-label', 'text-anchor': 'middle' }, `${t}×`);
  }
  return svg;
}

function tissueCard(t, notes, target) {
  const similar = notes.similar;
  const anchor = t.anchor || {};
  const unit = anchor.unit || '';
  const card = h('div', { class: 'expr-tissue' });
  const sourceChip = t.available
    ? h('span', { class: 'pill muted', title: t.source === 'stored' ? 'AlphaGenome RNA-seq stored in the dataset' : 'AlphaGenome RNA-seq predicted on demand (cached)' }, t.source === 'stored' ? 'stored prediction' : 'predicted on demand')
    : null;
  let expressedChip;
  if (anchor.value === null || anchor.value === undefined) expressedChip = h('span', { class: 'pill muted', title: anchor.error || '' }, 'no observed level');
  else if (t.expressed) expressedChip = h('span', { class: 'pill expressed' }, h('span', { class: 'status-dot good' }), `expressed · ${fmtAbs(anchor.value, unit)}`);
  else expressedChip = h('span', { class: 'pill muted' }, h('span', { class: 'status-dot' }), `not expressed · ${fmtAbs(anchor.value, unit)}`);
  card.appendChild(h('div', { class: 'expr-tissue-head' }, h('h3', null, t.label), h('div', { class: 'product-cell wrap' }, expressedChip, sourceChip)));
  if (!t.available) {
    card.appendChild(h('div', { class: 'empty' }, 'No RNA-seq prediction for this tissue (see the note above).'));
    return card;
  }
  const notExpressed = t.expressed === false;
  if (notExpressed) card.classList.add('not-expressed');
  const mf = t.mrna.fold;
  const pf = t.protein.fold;
  const dM = direction(mf.individual, similar);
  const dP = direction(pf.individual, similar);
  card.appendChild(h('p', { class: `expr-headline ${notExpressed ? 'muted' : ''}` },
    notExpressed
      ? `${target} is not expressed here (${fmtAbs(anchor.value, unit)} observed), so the predicted changes below are noise.`
      : h('span', null,
        'This individual makes ', h('b', { class: `dir-${dM}` }, `${fmtFold(mf.individual)} ${dM === 'same' ? '' : `(${fmtPctChange(mf.individual)})`}`),
        ` the reference genome's ${target} mRNA (${VERB[dM]}) and `, h('b', { class: `dir-${dP}` }, `${fmtFold(pf.individual)}`),
        ` its protein (${VERB[dP]}; relative estimate).`)));
  card.appendChild(foldBars([
    { label: 'mRNA · H1', fold: mf.H1 },
    { label: 'mRNA · H2', fold: mf.H2 },
    { label: 'mRNA · individual', fold: mf.individual, strong: true },
    { label: 'protein · H1', fold: pf.H1 },
    { label: 'protein · H2', fold: pf.H2 },
    { label: 'protein · individual', fold: pf.individual, strong: true },
  ], similar, notExpressed));
  const abs = t.mrna.absolute;
  const kindCell = (kind, what) => h('td', null, h('span', { class: `kind-badge ${kind}` }, kind === 'rel' ? 'RELATIVE' : 'ABSOLUTE'), ' ', what);
  const variantNote = (hap) => {
    const c = t.protein.variant && t.protein.variant[hap];
    return c && !['no_change', 'utr_change', 'synonymous'].includes(c.class) ? h('div', { class: 'mono product-hgvs' }, c.hgvs || c.class) : null;
  };
  card.appendChild(h('div', { class: 'table-wrap' }, h('table', { class: 'table compact expr-table' },
    h('thead', null, h('tr', null, h('th', null, ''), h('th', { class: 'num' }, 'Reference genome'), h('th', { class: 'num' }, 'H1 copy'), h('th', { class: 'num' }, 'H2 copy'), h('th', { class: 'num' }, 'Individual (H1 + H2)'))),
    h('tbody', null,
      h('tr', null, kindCell('rel', 'mRNA'), h('td', { class: 'num' }, '1.00×'), h('td', { class: 'num' }, fmtFold(mf.H1)), h('td', { class: 'num' }, fmtFold(mf.H2)), h('td', { class: `num dir-${dM}` }, h('b', null, fmtFold(mf.individual)))),
      h('tr', null, kindCell('abs', `mRNA (${unit || 'n/a'})`),
        h('td', { class: 'num', title: 'observed level of the gene in this tissue (population median), taken as the reference genome\'s' }, abs ? fmtAbs(abs.ref, unit) : '–', abs ? h('div', { class: 'muted small' }, `2 × ${fmtAbs(abs.reference_haplotype)}`) : null),
        h('td', { class: 'num', title: 'this copy\'s contribution' }, abs ? fmtAbs(abs.H1) : '–'), h('td', { class: 'num', title: 'this copy\'s contribution' }, abs ? fmtAbs(abs.H2) : '–'),
        h('td', { class: `num dir-${dM}` }, h('b', null, abs ? fmtAbs(abs.individual, unit) : '–'))),
      h('tr', null, kindCell('rel', 'protein'), h('td', { class: 'num' }, '1.00×'), h('td', { class: 'num' }, fmtFold(pf.H1), variantNote('H1')), h('td', { class: 'num' }, fmtFold(pf.H2), variantNote('H2')),
        h('td', { class: `num dir-${dP}` }, h('b', null, fmtFold(pf.individual)))),
      h('tr', null, kindCell('abs', 'protein'), h('td', { class: 'num muted', colspan: '4', title: notes.protein_absolute }, 'not available: no tissue-matched quantitative proteomics reference (hover for why)'))))));
  const meta = [
    `Predicted coverage (AlphaGenome units, mean over ${target}'s exons): reference ${fmtNum(t.mrna.predicted.ref, 4)}, H1 ${fmtNum(t.mrna.predicted.H1, 4)}, H2 ${fmtNum(t.mrna.predicted.H2, 4)}.`,
    anchor.value !== null && anchor.value !== undefined ? `Absolute anchor: ${anchor.source}${anchor.proxy ? ' (a proxy for this cell line)' : ''}.` : `No absolute anchor: ${anchor.error || 'no observed level for this tissue'}.`,
    `Protein-making share of transcripts: reference ${fmtNum(t.protein.productive_share.ref, 2)}, H1 ${fmtNum(t.protein.productive_share.H1, 2)}, H2 ${fmtNum(t.protein.productive_share.H2, 2)}.`,
  ];
  card.appendChild(h('ul', { class: 'expr-meta muted' }, meta.map((m) => h('li', null, m))));
  if (t.isoforms && t.isoforms.length) {
    card.appendChild(h('details', { class: 'expr-isoforms' }, h('summary', null, `Isoform shares (${t.isoforms.length} alternative junctions)`),
      h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['Isoform event', 'Share · ref', 'H1', 'H2', 'NMD'].map((c) => h('th', null, c)))),
        h('tbody', null, t.isoforms.map((r) => h('tr', null, h('td', null, r.kind), ...['ref', 'H1', 'H2'].map((k) => h('td', { class: 'num' }, fmtNum(r.share[k], 2))),
          h('td', null, r.nmd.H1 || r.nmd.H2 || r.nmd.ref ? h('span', { class: 'pill product-class bad' }, 'NMD') : '–')))))));
  }
  return card;
}

/** The expression card; ``params`` = {gene, sample, target}. */
export function expressionCard(dsPath, params, target, latestSignal) {
  const body = h('div', { class: 'card-body', style: { position: 'relative', minHeight: '140px' } });
  const card = h('section', { class: 'card expr-card' },
    h('div', { class: 'card-head' }, h('h2', null, 'How much? mRNA and protein against the reference genome')), body);
  const overlay = jobOverlay(body, `Expression of ${target} in ${params.sample}`);
  apiJob(`${dsPath}/products/expression`, { params, signal: latestSignal, onProgress: (job) => overlay.update(job) }).then((data) => {
    overlay.remove();
    clear(body);
    if (!data.available) { body.appendChild(h('div', { class: 'empty' }, data.message || 'No expression estimate.')); return; }
    const notes = data.notes;
    body.appendChild(h('div', { class: 'expr-legend' },
      h('div', { class: 'expr-def' }, h('span', { class: 'kind-badge rel' }, 'RELATIVE'), h('div', null, notes.relative, ' This ratio is what the model can be trusted with.')),
      h('div', { class: 'expr-def' }, h('span', { class: 'kind-badge abs' }, 'ABSOLUTE'), h('div', null, h('b', null, 'Units: TPM or nCPM. '), notes.absolute)),
      h('div', { class: 'expr-def' }, h('span', { class: 'kind-badge prot' }, 'PROTEIN'), h('div', null, notes.protein, ' ', h('span', { class: 'muted' }, notes.protein_absolute)))));
    if (data.problems && data.problems.length) body.appendChild(h('p', { class: 'warn-text' }, data.problems.join(' · ')));
    const groups = {};
    for (const t of data.tissues) (groups[t.group] = groups[t.group] || []).push(t);
    for (const [group, list] of Object.entries(groups)) {
      body.appendChild(h('h3', { class: 'section-title' }, group === 'LCL' ? 'Lymphoblastoid cells' : group === 'melanocyte' ? 'Melanocytes' : `Selected tissue (${group})`));
      body.appendChild(h('div', { class: 'expr-grid' }, list.map((t) => tissueCard(t, notes, data.target))));
    }
    body.appendChild(h('p', { class: 'muted small' }, `“About the same” = within ±${Math.round(notes.similar * 100)}% of the reference (shaded band). Folds are of the whole ${(data.exonic_bases || 0).toLocaleString('en-US')}-bp exonic sequence of ${data.target}, so they include every variant of the haplotype in the predicted window, coding or not.`));
  }).catch((err) => {
    overlay.remove();
    if (!isAbort(err)) { clear(body); body.appendChild(errorBox(err)); }
  });
  return card;
}

export { fmtFold, direction };
