// Info cards (drawers) for technical terms: a gene (HGNC record, database links, Gene Ontology
// annotations, dataset windows) and an AlphaGenome track (biosample ontology term, assay, ChIP
// target, data source and the observed ENCODE / FANTOM5 experiments behind it).
import { api, apiJob } from './api.js';
import { navigate } from './app.js';
import { state } from './state.js';
import { h, clear, drawer, errorBox, fmtInt } from './ui.js';
import {
  GO_EVIDENCE_URL, assayTerm, curieLink, curieLinks, encodeBiosampleUrl, geneLinks, goTermLinks, hgncGroupUrl, isCurie,
  targetLinks, xref, xrefs,
} from './links.js';

const OUTPUT_LABELS = { rna_seq: 'RNA-seq', cage: 'CAGE', procap: 'PRO-cap', dnase: 'DNase-seq', atac: 'ATAC-seq', chip_histone: 'Histone ChIP-seq', chip_tf: 'TF ChIP-seq', splice_sites: 'Splice sites', splice_site_usage: 'Splice site usage', splice_junctions: 'Splice junctions', contact_maps: 'Contact maps' };
const SOURCES = {
  encode: ['ENCODE', 'https://www.encodeproject.org/', 'Encyclopedia of DNA Elements'],
  fantom: ['FANTOM5', 'https://fantom.gsc.riken.jp/5/', 'FANTOM5 CAGE atlas (RIKEN)'],
  gtex: ['GTEx', 'https://gtexportal.org/', 'Genotype-Tissue Expression project (RNA-seq via recount3)'],
  '4dnucleome': ['4D Nucleome', 'https://data.4dnucleome.org/', '4D Nucleome data portal'],
};

const nameCache = new Map();
/** HGNC approved names of gene symbols: {symbol: {symbol, name, hgnc_id, locus_type}} (cached). */
export async function geneNames(symbols) {
  const missing = [...new Set(symbols)].filter((s) => !nameCache.has(s));
  if (missing.length) {
    const res = await api('/api/genes/names', { params: { symbols: missing.join(',') } });
    if (res.hgnc && res.hgnc.loaded) for (const s of missing) nameCache.set(s, res.names[s] || null);
    else return Object.fromEntries(symbols.map((s) => [s, nameCache.get(s) || res.names[s]]).filter(([, v]) => v));
  }
  return Object.fromEntries(symbols.map((s) => [s, nameCache.get(s)]).filter(([, v]) => v));
}

/** Label the options of a gene <select> with HGNC approved names ("TYRP1 — tyrosinase related protein 1"). */
export function labelGeneOptions(selectEl) {
  const symbols = [...selectEl.options].map((o) => o.value).filter(Boolean);
  if (!symbols.length) return;
  geneNames(symbols).then((names) => {
    for (const opt of selectEl.options) { const n = names[opt.value]; if (n) { opt.textContent = `${opt.value} — ${n.name}`; opt.title = `${n.hgnc_id} · ${n.locus_type}`; } }
  }).catch(() => {});
}

const section = (title, ...body) => h('section', { class: 'info-section' }, h('h3', null, title), ...body);
const kv = (rows) => h('dl', { class: 'kv' }, rows.filter((r) => r && r[1] !== null && r[1] !== undefined && r[1] !== '' && !(Array.isArray(r[1]) && !r[1].length)).map(([k, v]) => [h('dt', null, k), h('dd', null, v)]));
const loading = (text = 'Loading…') => h('div', { class: 'muted', style: { fontSize: '12px' } }, h('span', { class: 'spinner', style: { marginRight: '6px' } }), text);

// ------------------------------------------------------------------------------------ gene
/**
 * Gene card. ``symbol`` may be an approved symbol, alias, previous symbol or Ensembl id.
 * ``datasetId`` (default: the active dataset) lists the gene's windows with Tracks / Sequence buttons.
 */
export function openGeneCard(symbol, { datasetId = state.datasetId } = {}) {
  const head = h('div', null, loading());
  const linksHost = h('div');
  const windowsHost = h('div');
  const goHost = h('div', null, loading('Loading Gene Ontology annotations (QuickGO)…'));
  const body = h('div', { style: { display: 'grid', gap: '16px' } }, head, linksHost, windowsHost, goHost);
  const d = drawer(symbol, body);
  d.body.parentElement.classList.add('wide');
  api('/api/genes/info', { params: { symbol, ds: datasetId || '' } }).then((info) => {
    const rec = info.hgnc;
    const gencode = info.gencode;
    const title = d.body.parentElement.querySelector('.drawer-head h2');
    if (rec && title) title.textContent = rec.symbol;
    clear(head);
    if (rec) {
      head.append(
        h('div', { class: 'info-title' }, h('b', null, rec.name), h('span', { class: 'muted' }, ` · ${rec.locus_type}${rec.location ? ` · ${rec.location}` : ''}`)),
        kv([
          ['HGNC', xref({ label: rec.hgnc_id, url: `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${rec.hgnc_id}`, title: 'HGNC symbol report' })],
          ['Aliases', rec.aliases.join(', ')],
          ['Previous symbols', rec.previous.join(', ')],
          ['Gene groups', rec.groups.length ? h('span', { class: 'xrefs' }, rec.groups.map((g) => xref({ label: g.name, url: hgncGroupUrl(g.id), title: 'HGNC gene group' }))) : null],
          ['Ensembl', rec.ensembl_gene_id],
          ['NCBI Gene', rec.entrez_id],
          ['UniProt', rec.uniprot_ids.join(', ')],
          ['MANE Select', rec.mane_select.join(' · ')],
          ['GENCODE locus', gencode ? `${gencode.chrom}:${fmtInt(gencode.start)}-${fmtInt(gencode.end)} (${gencode.strand}) · ${String(gencode.type || '').replace(/_/g, ' ')}` : null],
        ]));
    } else {
      head.append(h('div', { class: 'notice' }, `${symbol} is not an approved HGNC symbol${info.hgnc_status && info.hgnc_status.error ? ` (HGNC unavailable: ${info.hgnc_status.error})` : ''}.`),
        gencode ? kv([['GENCODE locus', `${gencode.chrom}:${fmtInt(gencode.start)}-${fmtInt(gencode.end)} (${gencode.strand})`], ['Ensembl', gencode.id]]) : null);
    }
    const linkGene = { symbol: rec ? rec.symbol : symbol, ...(rec || {}), ...(gencode ? { chrom: gencode.chrom, start: gencode.start, end: gencode.end, id: gencode.id } : {}) };
    linksHost.append(section('Databases', xrefs(geneLinks(linkGene), { grouped: true })));
    if (info.windows && info.windows.length) {
      windowsHost.append(section('In this dataset', h('div', { class: 'chips' }, info.windows.map((w) => h('span', { class: 'chips' },
        h('b', null, w),
        h('button', { class: 'btn small', onclick: () => { d.close(); navigate('tracks', { gene: w }); } }, 'Tracks'),
        h('button', { class: 'btn small', onclick: () => { d.close(); navigate('sequence', { gene: w }); } }, 'Sequence'))))));
    }
    if (!rec) { clear(goHost); return; }
    api('/api/genes/go', { params: { symbol: rec.symbol } }).then((go) => {
      clear(goHost).append(section('Gene Ontology',
        h('p', { class: 'help' }, 'Annotations of ', go.uniprot_ids.map((u, i) => [i ? ', ' : '', xref({ label: `UniProtKB:${u}`, url: `https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=${u}`, title: 'All annotations in QuickGO' })]),
          ' from QuickGO. Evidence codes: ', xref({ label: 'guide', url: GO_EVIDENCE_URL, title: 'GO evidence codes' }), '.'),
        go.aspects.map((a) => goAspect(a))));
    }).catch((err) => { clear(goHost).append(section('Gene Ontology', errorBox(err))); });
  }).catch((err) => { clear(head).append(errorBox(err)); linksHost.append(section('Databases', xrefs(geneLinks({ symbol }), { grouped: true }))); clear(goHost); });
  return d;
}

function goAspect(aspect) {
  const terms = aspect.terms || [];
  const list = h('ul', { class: 'go-list' });
  const render = (n) => {
    clear(list);
    for (const t of terms.slice(0, n)) {
      const [quickgo, amigo] = goTermLinks(t.id);
      list.appendChild(h('li', { class: t.negated ? 'negated' : '', title: t.negated ? 'NOT annotated (negative evidence)' : t.qualifiers.join(', ') },
        h('a', { href: quickgo.url, target: '_blank', rel: 'noopener noreferrer', title: quickgo.title }, t.name || t.id),
        h('span', { class: 'meta mono' }, ' ', t.id),
        h('span', { class: 'meta' }, ` ${t.qualifiers.filter((q) => q !== 'enables' && q !== 'involved_in' && q !== 'located_in').join(', ')}`),
        h('span', { class: 'go-ev' }, t.evidence.join(' ')),
        xref({ ...amigo, label: 'AmiGO' }, { compact: true })));
    }
    if (terms.length > n) list.appendChild(h('li', null, h('button', { class: 'btn small ghost', onclick: () => render(terms.length) }, `Show all ${terms.length}`)));
  };
  render(25);
  return h('details', { class: 'collapsible', open: terms.length > 0 }, h('summary', null, `${aspect.label} (${terms.length})`), terms.length ? list : h('div', { class: 'muted', style: { fontSize: '12px' } }, 'No annotations'));
}

/** A gene symbol that opens its card (plus an optional extra title). */
export function geneChip(symbol, { label = null, title = null } = {}) {
  return h('button', { type: 'button', class: 'gene-link', title: title || `${symbol}: names, database links and Gene Ontology`, onclick: (e) => { e.stopPropagation(); openGeneCard(symbol); }, onpointerdown: (e) => e.stopPropagation() }, label || symbol);
}

// ----------------------------------------------------------------------------------- track
/** Definition of an ontology term from OLS, filled in when it arrives. */
export function termDefinition(curie) {
  const host = h('div', { class: 'muted', style: { fontSize: '12px' } }, isCurie(curie) ? 'Loading definition from OLS…' : '');
  if (!isCurie(curie)) return host;
  api('/api/ontology/term', { params: { curie } }).then((t) => {
    clear(host);
    if (!t.found) { host.textContent = `${curie} was not found in OLS.`; return; }
    host.className = 'term-def';
    host.append(h('div', null, h('b', null, t.label), h('span', { class: 'muted' }, ` (${t.ontology})`)), t.definition ? h('div', null, t.definition) : null,
      t.synonyms && t.synonyms.length ? h('div', { class: 'muted', style: { fontSize: '11.5px' } }, `Synonyms: ${t.synonyms.join('; ')}`) : null);
  }).catch((err) => { host.textContent = `Definition unavailable: ${err.message}`; });
  return host;
}

/**
 * Track card. ``observed`` is the /observed entry of this track when already loaded; otherwise
 * the card offers to look it up (``datasetId`` + ``gene`` are needed for that).
 */
export function openTrackCard({ output, index, meta = {}, label = '', gene = null, datasetId = state.datasetId, observed = null }) {
  const curie = meta.ontology_curie;
  const source = String(meta.data_source || (output === 'cage' ? 'fantom' : '')).toLowerCase();
  const assay = meta['Assay title'] || OUTPUT_LABELS[output] || output;
  const term = assayTerm(meta['Assay title'] || output);
  const target = meta.transcription_factor || meta.histone_mark;
  const biosampleLinks = [...curieLinks(curie)];
  const encodeBio = source === 'encode' ? encodeBiosampleUrl(curie, meta.biosample_type) : null;
  if (encodeBio) biosampleLinks.push({ label: 'ENCODE biosample', url: encodeBio, title: 'ENCODE experiments on this biosample type' });
  const observedHost = h('div');
  const body = h('div', { style: { display: 'grid', gap: '16px' } },
    section('Biosample',
      kv([['Name', meta.biosample_name], ['Ontology term', curie ? curieLink(curie) : null], ['Type', meta.biosample_type ? String(meta.biosample_type).replace(/_/g, ' ') : null], ['Life stage', meta.biosample_life_stage], ['GTEx tissue', meta.gtex_tissue], ['Genetically modified', meta.genetically_modified === undefined ? null : String(meta.genetically_modified)]]),
      termDefinition(curie),
      biosampleLinks.length ? xrefs(biosampleLinks) : null),
    section('Assay',
      kv([['Assay', assay], ['Ontology term', term ? h('span', null, curieLink(term.curie), ` ${term.label}`) : null], ['Strand', meta.strand], ['Endedness', meta.endedness], ['AlphaGenome output', `${OUTPUT_LABELS[output] || output} · track #${index}`]]),
      term ? xrefs(curieLinks(term.curie)) : null),
    target ? section(meta.histone_mark ? 'Histone mark' : 'Transcription factor',
      h('div', { class: 'chips' }, meta.transcription_factor ? h('button', { class: 'btn small', onclick: () => openGeneCard(target) }, `${target}: gene card`) : h('b', null, target), xrefs(targetLinks(target)))) : null,
    section('Data source',
      source && SOURCES[source] ? h('div', null, xref({ label: SOURCES[source][0], url: SOURCES[source][1], title: SOURCES[source][2] }), h('span', { class: 'muted' }, ` · ${SOURCES[source][2]}`)) : h('div', { class: 'muted' }, source || 'unknown'),
      meta.name ? h('div', { class: 'muted mono', style: { fontSize: '11.5px' } }, `AlphaGenome track name: ${meta.name}`) : null),
    observedHost,
    h('details', { class: 'collapsible' }, h('summary', null, 'All track metadata'), kv(Object.entries(meta).map(([k, v]) => [k, String(v)]))));
  const d = drawer(label || `${OUTPUT_LABELS[output] || output} · ${meta.biosample_name || curie || `track ${index}`}`, body);
  d.body.parentElement.classList.add('wide');
  const renderObserved = (item) => {
    clear(observedHost).append(section('Observed data', observedSummary(item)));
  };
  if (observed) renderObserved(observed);
  else if (gene && datasetId) {
    const find = h('button', { class: 'btn small', onclick: async () => {
      clear(observedHost).append(section('Observed data', loading('Looking up the experiments behind this track…')));
      try {
        const res = await apiJob(`/api/d/${encodeURIComponent(datasetId)}/observed`, { params: { gene, output, tracks: index, start: 0, end: 1, bins: 1 } });
        renderObserved(res.tracks[0]);
      } catch (err) { clear(observedHost).append(section('Observed data', errorBox(err))); }
    } }, 'Find the observed experiments');
    observedHost.append(section('Observed data', h('p', { class: 'help' }, 'The ENCODE / FANTOM5 experiments with the same ontology term, assay and target as this track.'), find));
  }
  return d;
}

/** Provider, units and experiments (with links) of an /observed track entry. */
export function observedSummary(item) {
  if (!item) return h('div', { class: 'muted' }, 'Not loaded');
  if (!item.available) return h('div', { class: 'notice' }, item.reason || 'No observed data');
  const rows = (item.items || []).map((it) => h('li', null,
    h('a', { href: it.url, target: '_blank', rel: 'noopener noreferrer', title: 'Experiment / sample page' }, it.id), ' ',
    h('span', null, it.label), h('div', { class: 'muted', style: { fontSize: '11.5px' } }, it.detail,
      it.file_url ? [' · ', h('a', { href: it.file_url, target: '_blank', rel: 'noopener noreferrer' }, 'file')] : null,
      it.loaded === false ? ' · could not be read' : '')));
  return h('div', { style: { display: 'grid', gap: '6px' } },
    h('div', null, h('b', null, item.provider), h('span', { class: 'muted' }, ` · ${item.units}${item.match === 'name' ? ' · matched by cell-line name' : ''}`)),
    h('ul', { class: 'observed-list' }, rows),
    item.total > (item.items || []).length ? h('div', { class: 'muted', style: { fontSize: '12px' } }, `Mean of ${(item.items || []).length} of ${fmtInt(item.total)} matching ${item.provider === 'FANTOM5' ? 'libraries' : 'experiments'}.`) : h('div', { class: 'muted', style: { fontSize: '12px' } }, `Mean of ${(item.items || []).length}.`),
    item.query ? xref({ label: `All matching ${item.provider} experiments`, url: item.query, title: 'Search on the ENCODE portal' }) : null,
    item.errors && item.errors.length ? h('div', { class: 'muted', style: { fontSize: '11.5px' } }, `Skipped: ${item.errors.join('; ')}`) : null);
}
