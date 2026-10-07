// Gene report of the Gene products page: (1) what the reference genome makes of the gene in a
// tissue, transcript by transcript (share and absolute level, from GTEx / HPA), and (2) what each
// of the individual's haplotypes changes: transcript levels, unstable mRNA, premature stops,
// different amino acids, protein similarity, and the variants behind them.
import { apiJob, isAbort } from './api.js';
import { h, clear, fmtInt, errorBox, jobOverlay, checkbox } from './ui.js';
import { transcriptLink } from './links.js';
import { knownChip, knownPanel } from './known_variants.js';

const HAPS = ['H1', 'H2'];
const PALETTE = ['--series-1', '--series-2', '--series-3', '--series-4', '--series-5', '--series-6', '--series-7', '--series-8'];
const CATEGORY = {
  frameshift: { label: 'frameshift', tone: 'bad' },
  nonsense: { label: 'nonsense (stop gained)', tone: 'bad' },
  start_lost: { label: 'start lost', tone: 'bad' },
  stop_lost: { label: 'stop lost', tone: 'bad' },
  splice: { label: 'splice site', tone: 'bad' },
  inframe_indel: { label: 'in-frame indel', tone: 'warn' },
  missense: { label: 'missense', tone: 'warn' },
  splice_region: { label: 'splice region', tone: '' },
  synonymous: { label: 'synonymous', tone: '' },
  utr: { label: 'UTR', tone: '' },
  noncoding_exon: { label: 'non-coding exon', tone: '' },
  intron: { label: 'intronic', tone: 'quiet' },
  flank: { label: 'up/downstream', tone: 'quiet' },
  structural: { label: 'structural', tone: 'warn' },
};
const COUNT_ORDER = ['synonymous', 'missense', 'nonsense', 'frameshift', 'inframe_indel', 'start_lost', 'stop_lost', 'splice', 'splice_region', 'utr', 'noncoding_exon', 'intron', 'flank', 'structural'];
const QUIET = new Set(['intron', 'flank']);

const pct = (v, d = 1) => (v === null || v === undefined || !Number.isFinite(v) ? '–' : `${(v * 100).toFixed(d)}%`);
const fold = (v) => (v === null || v === undefined || !Number.isFinite(v) ? '–' : `${v.toFixed(2)}×`);
const level = (v, unit) => (v === null || v === undefined || !Number.isFinite(v) ? '–' : `${v >= 100 ? Math.round(v).toLocaleString('en-US') : v >= 10 ? v.toFixed(1) : v.toFixed(2)}${unit ? ` ${unit}` : ''}`);
const dirClass = (v) => (v === null || v === undefined ? '' : v >= 1.1 ? 'dir-up' : v <= 1 / 1.1 ? 'dir-down' : 'dir-same');
const trim = (s, n = 10) => (s && s.length > n ? `${s.slice(0, n - 1)}…` : s || '');
const AA3 = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter' };
const vkey = (v) => `${v.pos}:${v.ref}:${v.alt}`;
/** The variant of this copy behind a substitution of the base protein (e.g. R151C <- p.Arg151Cys). */
function variantFor(hapData, m) {
  const hgvs = `p.${AA3[m.ref]}${m.pos}${AA3[m.alt]}`;
  return hapData.variants.find((v) => v.base_change && v.base_change.hgvs === hgvs);
}

function amChip(am) {
  if (!am) return null;
  const tone = am.class === 'likely pathogenic' ? 'bad' : am.class === 'ambiguous' ? 'warn' : 'good';
  return h('span', { class: `pill product-class ${tone}`, title: 'AlphaMissense pathogenicity (0 = benign, 1 = pathogenic): a disease-risk proxy, not a measure of function' }, `AlphaMissense ${am.score.toFixed(2)} · ${am.class}`);
}

function categoryChip(cat, count) {
  const c = CATEGORY[cat] || { label: cat, tone: '' };
  return h('span', { class: `pill cat-chip ${c.tone}` }, count !== undefined ? h('b', null, String(count)) : null, ` ${c.label}`);
}

/** The flag rows at the top of a haplotype column. */
function verdicts(hapData, report, hap, known = {}) {
  const f = hapData.flags;
  const unit = report.level.unit || '';
  const rows = [];
  const row = (tone, icon, title, detail) => rows.push(h('div', { class: `verdict ${tone}` }, h('span', { class: 'verdict-icon' }, icon), h('div', null, h('div', { class: 'verdict-title' }, title), detail ? h('div', { class: 'verdict-detail' }, detail) : null)));
  // expression
  const g = hapData.gene_fold;
  if (g === null || g === undefined) row('', '?', 'mRNA level: no prediction for this tissue', (report.expression_problems || []).join(' '));
  else {
    const notExpr = report.level.expressed === false;
    row(notExpr ? 'quiet' : g >= 1.1 ? 'up' : g <= 1 / 1.1 ? 'down' : 'ok', g >= 1.1 ? '▲' : g <= 1 / 1.1 ? '▼' : '=',
      h('span', null, 'mRNA ', h('b', { class: dirClass(g) }, fold(g)), ' the reference copy', notExpr ? ' (gene not expressed here: noise)' : ''),
      report.level.value !== null && report.level.value !== undefined ? `≈ ${level(hapData.level, unit)} from this copy vs ${level(report.level.value / 2, unit)} from a reference copy` : 'no absolute reference level for this tissue');
  }
  // stability
  const share = (x) => (x.reference_share === null || x.reference_share === undefined ? 'expression share unknown' : `${pct(x.reference_share, 0)} of the gene's mRNA`);
  if (f.unstable.length) row('bad', '✕', `Unstable mRNA: nonsense-mediated decay in ${f.unstable.length} transcript${f.unstable.length > 1 ? 's' : ''}`,
    f.unstable.map((x) => `${x.name} (${share(x)}): kept at ~${pct(report.notes.nmd_residual, 0)}`).join('; '));
  else row('ok', '✓', 'mRNA stable: no transcript gains a decay-triggering stop');
  // premature stop
  if (f.premature_stop.length) row('bad', '■', `Premature stop codon in ${f.premature_stop.length} transcript${f.premature_stop.length > 1 ? 's' : ''}`,
    f.premature_stop.map((x) => `${x.name} ${x.change.hgvs || ''} (${share(x)})`).join('; '));
  else row('ok', '✓', 'No premature stop codon');
  if (f.start_lost.length) row('bad', '■', 'Start codon lost', f.start_lost.map((x) => `${x.name} (${share(x)})`).join('; '));
  // amino acids (base transcript)
  const base = f.base_change || {};
  if (f.missense.length) {
    row('warn', '◆', `Different amino acid${f.missense.length > 1 ? 's' : ''} in ${report.base_name}`,
      h('div', { class: 'verdict-list' }, f.missense.map((m) => {
        const v = variantFor(hapData, m);
        const rec = v && known[vkey(v)];
        const traits = rec && (rec.traits || []).slice(0, 3).map((t) => t.trait).join('; ');
        return h('div', null, h('b', { class: 'mono' }, `${m.ref}${m.pos}${m.alt}`), ' ', amChip(m.alphamissense), ' ', rec ? knownChip(rec) : null,
          traits ? h('div', { class: 'muted small' }, `Known for: ${traits}${(rec.trait_count || 0) > 3 ? ` (+${rec.trait_count - 3})` : ''}`) : null);
      })));
  } else if (['frameshift', 'inframe_indel', 'stop_gained', 'stop_lost', 'start_lost'].includes(base.class)) {
    row('bad', '◆', `${report.base_name} protein altered`, base.hgvs);
  } else {
    const syn = hapData.counts.synonymous || 0;
    row('ok', '=', `Same amino acids in ${report.base_name}`, syn ? `${syn} synonymous variant${syn > 1 ? 's' : ''}: DNA changes, protein unchanged` : null);
  }
  // similarity
  const sim = f.base_similarity;
  if (sim && sim.similarity !== null && sim.similarity !== undefined) {
    const tone = sim.similarity >= 0.999 ? 'ok' : sim.similarity >= 0.95 ? 'warn' : 'bad';
    rows.push(h('div', { class: `verdict similarity ${tone}` },
      h('div', { class: 'similarity-score' }, pct(sim.similarity, sim.similarity >= 0.99 ? 2 : 1)),
      h('div', null, h('div', { class: 'verdict-title' }, `Protein similarity to the reference ${report.base_name}`),
        h('div', { class: 'verdict-detail' }, `BLOSUM62 alignment score relative to the reference protein against itself · identity ${pct(sim.identity, 2)} · ${fmtInt(sim.length)} of ${fmtInt(sim.reference_length)} residues`))));
  }
  return h('div', { class: 'verdicts' }, rows);
}

function hapColumn(report, hap, onSelect, state) {
  const d = report.haplotypes[hap];
  const unit = report.level.unit || '';
  const known = state.known || {};
  const col = h('div', { class: 'hap-col' }, h('h3', { class: 'hap-title' }, `${hap} copy`), verdicts(d, report, hap, known));
  // known variants carried on this copy (dbSNP / ClinVar / GWAS / literature), most severe first
  const knownHere = d.variants.filter((v) => known[vkey(v)] && known[vkey(v)].known);
  if (knownHere.length) {
    col.appendChild(h('details', { class: 'known-section', open: true },
      h('summary', null, `Known variants on this copy (${knownHere.length}): clinical records, traits and literature`),
      knownHere.slice(0, 8).map((v) => knownPanel(known[vkey(v)], { chrom: report.chromosome, pos: v.pos, ref: v.ref, alt: v.alt }, { change: v.base_change && v.base_change.hgvs !== 'p.(=)' ? v.base_change.hgvs : null })),
      knownHere.length > 8 ? h('p', { class: 'muted small' }, `${knownHere.length - 8} more in the variant list below.`) : null));
  } else if (state.knownStatus) {
    col.appendChild(h('p', { class: 'muted small' }, state.knownStatus));
  }
  // variant counts
  const counts = COUNT_ORDER.filter((k) => d.counts[k]).map((k) => categoryChip(k, d.counts[k]));
  col.appendChild(h('div', { class: 'cat-row' }, h('span', { class: 'muted small' }, 'Variants on this copy: '), counts.length ? counts : h('span', { class: 'muted' }, 'none in this gene')));
  // transcripts
  col.appendChild(h('div', { class: 'table-wrap' }, h('table', { class: 'table compact hap-table' },
    h('thead', null, h('tr', null, ['Transcript', 'Share', `Level (${unit || '–'})`, 'vs ref', 'mRNA', 'Protein', 'Similarity'].map((c) => h('th', null, c)))),
    h('tbody', null, d.transcripts.map((t) => {
      const base = report.baseline.find((b) => b.id === t.id) || {};
      const cls = (t.change || {}).class;
      const tone = ['stop_gained', 'frameshift', 'start_lost', 'stop_lost'].includes(cls) ? 'bad' : ['missense', 'inframe_indel'].includes(cls) ? 'warn' : '';
      const sel = state.selectedTx === t.id && state.selectedHap === hap;
      return h('tr', { class: `clickable${sel ? ' selected' : ''}`, title: 'Show this protein in 3D below', onclick: () => onSelect(t.id, hap) },
        h('td', null, transcriptLink(t.id, t.name), base.mane ? h('span', { class: 'pill', style: { marginLeft: '4px' } }, 'MANE') : null),
        h('td', { class: 'num' }, t.share === null ? h('span', { class: 'muted', title: 'not in GTEx v8 (newer annotation)' }, '?') : pct(t.share, 0)),
        h('td', { class: 'num' }, level(t.level)),
        h('td', { class: `num ${dirClass(t.fold)}` }, fold(t.fold)),
        h('td', null, t.nmd_gained ? h('span', { class: 'pill product-class bad' }, 'unstable (NMD)') : h('span', { class: 'muted' }, 'stable')),
        h('td', null, base.coding ? h('span', { class: `pill product-class ${tone}`, title: (t.change || {}).hgvs || '' }, cls === 'no_change' || cls === 'utr_change' ? 'identical' : (cls || '').replace('_', ' ')) : h('span', { class: 'muted' }, 'non-coding'),
          t.change && t.change.hgvs && t.change.hgvs !== 'p.(=)' ? h('div', { class: 'mono product-hgvs' }, trim(t.change.hgvs, 34)) : null),
        h('td', { class: 'num' }, t.similarity && t.similarity.similarity !== null ? pct(t.similarity.similarity, t.similarity.similarity >= 0.99 ? 2 : 1) : '–'));
    })))));
  // variants
  const showAll = state.allVariants;
  const shown = d.variants.filter((v) => showAll || !QUIET.has(v.category));
  col.appendChild(h('details', { class: 'hap-variants', open: shown.length <= 12 ? true : null },
    h('summary', null, `Variants (${shown.length}${showAll ? '' : ` coding, splice and UTR; ${d.variants.length - shown.length} intronic / flanking hidden`})`),
    shown.length ? h('div', { class: 'table-wrap', style: { maxHeight: '320px' } }, h('table', { class: 'table compact' },
      h('thead', null, h('tr', null, ['Position', 'Change', 'Known as', 'Effect', `On ${report.base_name}`, 'AlphaMissense', 'Transcripts'].map((c) => h('th', null, c)))),
      h('tbody', null, shown.map((v) => h('tr', null,
        h('td', { class: 'mono' }, fmtInt(v.pos)),
        h('td', { class: 'mono', title: `${v.ref} > ${v.alt}` }, `${trim(v.ref, 8)}>${trim(v.alt, 8)}`),
        h('td', null, state.known ? knownChip(known[vkey(v)]) : h('span', { class: 'muted small' }, '…')),
        h('td', null, categoryChip(v.category)),
        h('td', { class: 'mono' }, v.base_change ? trim(v.base_change.hgvs, 26) : '–'),
        h('td', null, v.alphamissense ? amChip(v.alphamissense) : '–'),
        h('td', { class: 'muted small', title: Object.entries(v.transcripts).map(([k, x]) => `${k}: ${x.class} ${x.hgvs}`).join('\n') },
          Object.keys(v.transcripts).length ? `${Object.keys(v.transcripts).length} coding` : v.region)))))) : h('div', { class: 'empty' }, 'None.')));
  return col;
}

/** Where this tissue's numbers come from: anchor, transcript shares, AlphaGenome RNA-seq and junctions. */
function tissueSourceNote(report) {
  const t = report.tissue || {};
  const bits = [];
  bits.push(t.anchor ? `Gene level: ${t.anchor}${t.anchor_matched_by_name ? ' (matched by name)' : ''}` : 'Gene level: no observed reference for this tissue, so relative changes only');
  bits.push(`Transcript shares: ${(report.transcript_source || {}).label || 'none'}`);
  bits.push(`AlphaGenome RNA-seq: ${report.expression_source === 'stored' ? 'stored in the dataset' : report.expression_source === 'on_demand' ? 'predicted on demand' : 'unavailable'}`);
  bits.push(`Splicing shift: ${t.junctions === 'stored' ? 'stored junction predictions' : t.junctions === 'on_demand' ? 'junctions predicted on demand' : 'no junction prediction (shift not applied)'}`);
  return h('div', { class: 'tissue-sources' }, h('span', { class: 'mono muted small' }, t.ontology), ...bits.map((b) => h('span', { class: 'pill muted' }, b)));
}

function baselineCard(report) {
  const unit = report.level.unit || '';
  const lv = report.level;
  const src = report.transcript_source || {};
  const withShare = report.baseline.filter((b) => b.share !== null && b.share !== undefined);
  const bar = h('div', { class: 'share-bar' });
  withShare.filter((b) => b.share > 0).forEach((b, i) => {
    bar.appendChild(h('div', { class: `share-seg${b.coding ? '' : ' noncoding'}`, style: { width: `${(b.share * 100).toFixed(2)}%`, background: `var(${PALETTE[i % PALETTE.length]})` }, title: `${b.name}: ${pct(b.share)}` }, b.share >= 0.08 ? b.name : ''));
  });
  return h('section', { class: 'card' },
    h('div', { class: 'card-head' }, h('h2', null, h('span', { class: 'step' }, '1'), ` Reference genome: what ${report.target} makes in ${report.tissue.label}`)),
    h('div', { class: 'card-body' },
      tissueSourceNote(report),
      h('p', { class: 'baseline-level' }, lv.value === null || lv.value === undefined
        ? h('span', { class: 'muted' }, `No observed level of ${report.target} for this tissue (${lv.error || 'no reference'}).`)
        : h('span', null, h('b', null, level(lv.value, unit)), lv.expressed ? '' : ' (not expressed: below 1)', ` of ${report.target} mRNA · ${lv.source}`)),
      withShare.length && withShare.some((b) => b.share > 0) ? bar : h('p', { class: 'muted' }, src.error ? `No transcript breakdown: ${src.error}` : 'No transcript of this gene is expressed in the GTEx tissue.'),
      h('div', { class: 'table-wrap' }, h('table', { class: 'table compact' },
        h('thead', null, h('tr', null, ['Transcript', 'Type', 'Protein', 'Share of the gene', `Level (${unit || '–'})`, 'GTEx TPM'].map((c) => h('th', null, c)))),
        h('tbody', null, report.baseline.map((b, i) => h('tr', null,
          h('td', null, h('span', { class: 'swatch', style: { background: b.share > 0 ? `var(${PALETTE[withShare.filter((x) => x.share > 0).indexOf(b) % PALETTE.length]})` : 'transparent' } }), transcriptLink(b.id, b.name), b.mane ? h('span', { class: 'pill', style: { marginLeft: '4px' } }, 'MANE') : null),
          h('td', { class: 'muted' }, (b.type || '').replace(/_/g, ' ')),
          h('td', { class: 'num' }, b.protein_length ? `${fmtInt(b.protein_length)} aa` : '–'),
          h('td', { class: 'num' }, b.in_gtex ? pct(b.share) : h('span', { class: 'muted', title: 'Not in GTEx v8 (GENCODE v26) and no v26 transcript has its intron chain' }, 'not in GTEx'),
            b.gtex_match ? h('div', { class: 'muted small', title: 'GTEx v8 predates this transcript id; the v26 transcript with the identical intron chain stands in' }, `as ${b.gtex_match}`) : null),
          h('td', { class: 'num' }, level(b.level)),
          h('td', { class: 'num muted' }, b.gtex_tpm === null ? '–' : b.gtex_tpm.toFixed(2))))))),
      h('p', { class: 'muted small' }, 'Transcripts: GENCODE v46 models and names (gene symbol + the number Ensembl/HAVANA gave the transcript); each name links to its Ensembl page. ', `Shares: ${src.label || '–'}${src.proxy ? ' (a stand-in: GTEx has no transcript data for this tissue itself)' : ''}. Level = the gene's observed level × share. Observed in a reference population, not predicted.`)));
}

/** Loads the report into ``host``; ``onSelect(txId, hap)`` when a transcript row is clicked. */
export function reportCards(host, dsPath, params, tissue, signal, onSelect, state, onLoaded) {
  clear(host);
  const box = h('section', { class: 'card', style: { position: 'relative', minHeight: '160px' } });
  host.appendChild(box);
  const overlay = jobOverlay(box, 'Gene report');
  apiJob(`${dsPath}/products/report`, { params: { ...params, tissue }, signal, onProgress: (job) => overlay.update(job) }).then((report) => {
    overlay.remove();
    box.remove();
    if (!report.available) { host.appendChild(h('section', { class: 'card' }, h('div', { class: 'card-body' }, h('div', { class: 'empty' }, report.message || 'No report.')))); return; }
    const loadKnown = () => {
      const scope = state.allVariants ? 'all' : 'coding';
      if (state.knownScope === scope || state.knownScope === 'all') return;
      state.knownStatus = 'Looking up known variants (Ensembl VEP, ClinVar, GWAS Catalog)…';
      apiJob(`${dsPath}/products/known`, { params: { ...params, scope }, signal }).then((data) => {
        state.knownScope = scope;
        state.known = { ...(state.known || {}), ...(data.variants || {}) };
        state.knownStatus = data.available === false ? `Known-variant lookup unavailable: ${data.error || data.message}` : null;
        render();
      }).catch((err) => { if (!isAbort(err)) { state.knownStatus = `Known-variant lookup failed: ${err.message}`; render(); } });
    };
    const render = () => {
      clear(host);
      const toggle = checkbox('Show intronic and flanking variants', state.allVariants, (v) => { state.allVariants = v; render(); loadKnown(); });
      host.append(baselineCard(report),
        h('section', { class: 'card' },
          h('div', { class: 'card-head' }, h('h2', null, h('span', { class: 'step' }, '2'), ` ${report.sample}'s two copies of ${report.target}, in ${report.tissue.label}`), toggle),
          h('div', { class: 'card-body' },
            h('div', { class: 'hap-grid' }, HAPS.map((hap) => hapColumn(report, hap, (tx, hp) => { state.selectedTx = tx; state.selectedHap = hp; render(); onSelect(tx, hp); }, state))),
            h('p', { class: 'muted small' }, `Level of a transcript on a copy = reference level / 2 × AlphaGenome's predicted change of the gene's mRNA × the change in the transcript's splice-junction usage × ${report.notes.nmd_residual} if the copy's mRNA gains nonsense-mediated decay. Protein changes combine every variant of the copy (phase-aware); the variant list shows each variant alone on the reference. Click a transcript to see its protein in 3D.`))));
    };
    state.known = null;
    state.knownScope = null;
    render();
    loadKnown();
    if (onLoaded) onLoaded(report);
  }).catch((err) => {
    overlay.remove();
    if (!isAbort(err)) box.appendChild(errorBox(err));
  });
}
