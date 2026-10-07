// Links to ontology and genome databases for the technical terms the app shows: genes (HGNC,
// Ensembl, NCBI, UniProt, Gene Ontology, ...), ontology terms (CL, UBERON, EFO, OBI, GO via OLS /
// Ontobee), assays, ChIP targets, ENCODE biosamples and experiments, FANTOM5 samples, 1000 Genomes
// samples / populations and dbSNP variants. Pure functions plus a small DOM helper.
import { h } from './ui.js';

const enc = encodeURIComponent;

// IRI namespaces (EFO is not an OBO-PURL ontology).
const IRI_BASE = { EFO: 'http://www.ebi.ac.uk/efo/EFO_', ORPHA: 'http://www.orpha.net/ORDO/Orphanet_' };
export const ONTOLOGY_NAMES = {
  CL: 'Cell Ontology', UBERON: 'Uberon anatomy ontology', EFO: 'Experimental Factor Ontology', CLO: 'Cell Line Ontology',
  OBI: 'Ontology for Biomedical Investigations', GO: 'Gene Ontology', HP: 'Human Phenotype Ontology', MONDO: 'Mondo disease ontology',
  NCIT: 'NCI Thesaurus', BTO: 'BRENDA Tissue Ontology', SO: 'Sequence Ontology', ECO: 'Evidence & Conclusion Ontology', HANCESTRO: 'Human Ancestry Ontology',
};

// Assays of AlphaGenome tracks (ENCODE / FANTOM5 "Assay title") -> ontology terms (checked in OLS).
export const ASSAY_TERMS = {
  CAGE: ['OBI:0001674', 'cap analysis of gene expression assay'],
  hCAGE: ['OBI:0001674', 'cap analysis of gene expression assay'],
  LQhCAGE: ['OBI:0001674', 'cap analysis of gene expression assay'],
  'DNase-seq': ['OBI:0001853', 'DNase I hypersensitive sites sequencing assay'],
  'ATAC-seq': ['OBI:0002039', 'assay for transposase-accessible chromatin using sequencing'],
  'RNA-seq': ['OBI:0001271', 'RNA-seq assay'],
  'polyA plus RNA-seq': ['OBI:0002571', 'polyA-selected RNA sequencing assay'],
  'total RNA-seq': ['EFO:0009653', 'RNA-seq of total RNA'],
  'ChIP-seq': ['OBI:0000716', 'ChIP-seq assay'],
  'TF ChIP-seq': ['OBI:0000716', 'ChIP-seq assay'],
  'Histone ChIP-seq': ['OBI:0002017', 'histone modification identification by ChIP-Seq assay'],
  'PRO-cap': ['OBI:0002753', 'PRO-cap'],
  'in situ Hi-C': ['OBI:0002440', 'Hi-C assay'],
  'Dilution Hi-C': ['OBI:0002440', 'Hi-C assay'],
};
const OUTPUT_ASSAY = { rna_seq: 'RNA-seq', cage: 'CAGE', procap: 'PRO-cap', dnase: 'DNase-seq', atac: 'ATAC-seq', chip_histone: 'Histone ChIP-seq', chip_tf: 'TF ChIP-seq', splice_sites: 'RNA-seq', splice_site_usage: 'RNA-seq', splice_junctions: 'RNA-seq', contact_maps: 'in situ Hi-C' };

export const isCurie = (text) => /^[A-Za-z][A-Za-z0-9_]*:[A-Za-z0-9_.-]+$/.test(String(text || ''));

/** The term's IRI: OBO PURL for most prefixes (``CL:0002567`` -> ``.../obo/CL_0002567``). */
export function curieIri(curie) {
  const [prefix, id] = String(curie).split(':');
  return `${IRI_BASE[prefix] || `http://purl.obolibrary.org/obo/${prefix}_`}${id}`;
}

export function olsUrl(curie, iri = null) {
  const prefix = String(curie).split(':')[0].toLowerCase();
  const ontology = prefix === 'efo' || prefix === 'orpha' ? (prefix === 'orpha' ? 'ordo' : 'efo') : prefix;
  return `https://www.ebi.ac.uk/ols4/ontologies/${ontology}/classes/${enc(enc(iri || curieIri(curie)))}`;
}

/** Links for any ontology CURIE (OLS, Ontobee, and the ontology's own browser when there is one). */
export function curieLinks(curie, { iri = null } = {}) {
  if (!isCurie(curie)) return [];
  const prefix = String(curie).split(':')[0].toUpperCase();
  const links = [{ label: 'OLS', url: olsUrl(curie, iri), title: `EBI Ontology Lookup Service · ${ONTOLOGY_NAMES[prefix] || prefix}` }];
  if (prefix === 'GO') return goTermLinks(curie);
  if (prefix !== 'EFO') links.push({ label: 'Ontobee', url: `https://ontobee.org/ontology/${prefix}?iri=${enc(iri || curieIri(curie))}`, title: 'Ontobee linked-data browser' });
  if (prefix === 'CL') links.push({ label: 'CELLxGENE', url: `https://cellxgene.cziscience.com/cellguide/${curie.replace(':', '_')}`, title: 'CellGuide (CZ CELLxGENE): marker genes and datasets for this cell type' });
  return links;
}

export function goTermLinks(goId) {
  return [
    { label: 'QuickGO', url: `https://www.ebi.ac.uk/QuickGO/term/${goId}`, title: 'QuickGO (EBI): term, ancestors, annotations' },
    { label: 'AmiGO', url: `https://amigo.geneontology.org/amigo/term/${goId}`, title: 'AmiGO 2 (Gene Ontology Consortium)' },
    { label: 'OLS', url: olsUrl(goId), title: 'EBI Ontology Lookup Service' },
  ];
}

export const GO_EVIDENCE_URL = 'https://geneontology.org/docs/guide-go-evidence-codes/';

/** Assay (``Assay title`` or AlphaGenome output name) -> {curie, label} of its OBI / EFO term. */
export function assayTerm(assayOrOutput) {
  const key = ASSAY_TERMS[assayOrOutput] ? assayOrOutput : OUTPUT_ASSAY[String(assayOrOutput || '').toLowerCase()];
  const term = key && ASSAY_TERMS[key];
  return term ? { curie: term[0], label: term[1], assay: key } : null;
}

/**
 * Database links for a gene. ``gene`` may hold ``symbol`` and any of ``hgnc_id``,
 * ``ensembl_gene_id`` (or ``id``), ``entrez_id``, ``uniprot_ids``, ``omim_ids`` and a locus
 * (``chrom``/``start``/``end``); links that need a missing identifier fall back to a search.
 */
export function geneLinks(gene) {
  const symbol = gene.symbol || gene.name || '';
  const ensembl = gene.ensembl_gene_id || (String(gene.id || '').startsWith('ENSG') ? String(gene.id).split('.')[0] : '');
  const uniprot = (gene.uniprot_ids || [])[0];
  const links = [];
  links.push({ group: 'Nomenclature', label: 'HGNC', url: gene.hgnc_id ? `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${gene.hgnc_id}` : `https://www.genenames.org/tools/search/#!/?query=${enc(symbol)}`, title: 'HUGO Gene Nomenclature Committee: approved symbol, aliases, gene groups' });
  links.push({ group: 'Genome', label: 'Ensembl', url: ensembl ? `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${ensembl}` : `https://www.ensembl.org/Homo_sapiens/Search/Results?q=${enc(symbol)}`, title: 'Ensembl gene: transcripts, regulation, variants' });
  links.push({ group: 'Genome', label: 'NCBI Gene', url: gene.entrez_id ? `https://www.ncbi.nlm.nih.gov/gene/${gene.entrez_id}` : `https://www.ncbi.nlm.nih.gov/gene/?term=${enc(`${symbol}[sym] AND human[orgn]`)}`, title: 'NCBI Gene (RefSeq, publications, GeneRIFs)' });
  if (gene.chrom && gene.start && gene.end) links.push({ group: 'Genome', label: 'UCSC', url: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${enc(`${gene.chrom}:${gene.start}-${gene.end}`)}`, title: 'UCSC Genome Browser (hg38) at this gene' });
  links.push({ group: 'Function', label: 'UniProt', url: uniprot ? `https://www.uniprot.org/uniprotkb/${uniprot}/entry` : `https://www.uniprot.org/uniprotkb?query=${enc(`gene_exact:${symbol} AND organism_id:9606`)}`, title: 'UniProtKB protein entry' });
  links.push({ group: 'Gene Ontology', label: 'QuickGO', url: uniprot ? `https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=${uniprot}` : `https://www.ebi.ac.uk/QuickGO/search/${enc(symbol)}`, title: 'Gene Ontology annotations (QuickGO, EBI)' });
  links.push({ group: 'Gene Ontology', label: 'AmiGO', url: uniprot ? `https://amigo.geneontology.org/amigo/gene_product/UniProtKB:${uniprot}` : `https://amigo.geneontology.org/amigo/search/bioentity?q=${enc(symbol)}`, title: 'Gene Ontology annotations (AmiGO 2)' });
  links.push({ group: 'Expression', label: 'GTEx', url: `https://gtexportal.org/home/gene/${enc(symbol)}`, title: 'GTEx Portal: expression across tissues, eQTLs' });
  if (ensembl) links.push({ group: 'Expression', label: 'Protein Atlas', url: `https://www.proteinatlas.org/${ensembl}`, title: 'Human Protein Atlas: tissue / cell type expression' });
  if (ensembl) links.push({ group: 'Variation & disease', label: 'gnomAD', url: `https://gnomad.broadinstitute.org/gene/${ensembl}?dataset=gnomad_r4`, title: 'gnomAD: variants and constraint' });
  for (const omim of (gene.omim_ids || []).slice(0, 2)) links.push({ group: 'Variation & disease', label: `OMIM ${omim}`, url: `https://www.omim.org/entry/${omim}`, title: 'Online Mendelian Inheritance in Man' });
  if (ensembl) links.push({ group: 'Variation & disease', label: 'Open Targets', url: `https://platform.opentargets.org/target/${ensembl}`, title: 'Open Targets: disease associations' });
  links.push({ group: 'Variation & disease', label: 'ClinVar', url: `https://www.ncbi.nlm.nih.gov/clinvar/?term=${enc(`${symbol}[gene]`)}`, title: 'ClinVar variants in this gene' });
  links.push({ group: 'Other', label: 'GeneCards', url: `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${enc(symbol)}`, title: 'GeneCards summary' });
  return links;
}

export const hgncGroupUrl = (id) => `https://www.genenames.org/data/genegroup/#!/group/${id}`;

/** ChIP-seq target (histone mark or TF): ENCODE target page; TFs also link to their gene. */
/**
 * A transcript's GENCODE/Ensembl name (e.g. DDB1-225: gene symbol + the number Ensembl/HAVANA gave
 * the transcript) linked to its Ensembl transcript page. Clicks do not reach the enclosing row.
 */
export function transcriptLink(id, name, { bold = true } = {}) {
  const stable = String(id || '').split('.')[0];
  const label = bold ? h('b', null, name || stable) : (name || stable);
  if (!stable.startsWith('ENST')) return h('span', null, label);
  return h('a', {
    class: 'tx-link', href: `https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=${stable}`, target: '_blank', rel: 'noopener',
    title: `${name ? `${name} = ` : ''}${id} · GENCODE / Ensembl transcript: open it on Ensembl`, onclick: (e) => e.stopPropagation(),
  }, label);
}

export function targetLinks(target) {
  if (!target) return [];
  const links = [{ label: 'ENCODE target', url: `https://www.encodeproject.org/targets/${enc(target)}-human/`, title: 'ENCODE target page (experiments, antibodies)' }];
  if (/^H[1-4]|^H2A|^H2B/i.test(target)) links.push({ label: 'OLS', url: `https://www.ebi.ac.uk/ols4/search?q=${enc(target)}`, title: 'Search the histone mark in the Ontology Lookup Service' });
  return links;
}

/** ENCODE page for a biosample term (``primary_cell_CL_1000458``). */
export function encodeBiosampleUrl(curie, type) {
  if (!isCurie(curie) || !type) return null;
  return `https://www.encodeproject.org/biosample-types/${type}_${curie.replace(':', '_')}/`;
}

export const dbsnpUrl = (rsid) => `https://www.ncbi.nlm.nih.gov/snp/${enc(rsid)}`;
export const igsrSampleUrl = (id) => `https://www.internationalgenome.org/data-portal/sample/${enc(id)}`;
export const igsrPopulationUrl = (code) => `https://www.internationalgenome.org/data-portal/population/${enc(code)}`;
/** 1000 Genomes / HGDP sample ids. */
export const isIgsrSample = (id) => /^(HG|NA)\d{5}$/.test(String(id || ''));
const POPULATIONS = new Set('ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI'.split(' '));
export const isIgsrPopulation = (code) => POPULATIONS.has(String(code || ''));

// ------------------------------------------------------------------------------------- DOM
/** An external link (new tab) with a small arrow. */
export function xref(link, { compact = false } = {}) {
  return h('a', { class: `xref${compact ? ' compact' : ''}`, href: link.url, target: '_blank', rel: 'noopener noreferrer', title: link.title || link.url, onclick: (e) => e.stopPropagation(), onpointerdown: (e) => e.stopPropagation() }, link.label);
}

/** A row of links, optionally grouped (``link.group``). */
export function xrefs(links, { grouped = false } = {}) {
  if (!grouped) return h('span', { class: 'xrefs' }, links.map((l) => xref(l)));
  const groups = new Map();
  for (const l of links) { const g = l.group || ''; if (!groups.has(g)) groups.set(g, []); groups.get(g).push(l); }
  return h('dl', { class: 'kv xref-groups' }, [...groups.entries()].map(([g, ls]) => [h('dt', null, g), h('dd', null, h('span', { class: 'xrefs' }, ls.map((l) => xref(l))))]));
}

/** A CURIE shown as a link to OLS (``CL:0002567 ↗``). */
export function curieLink(curie, label = null) {
  if (!isCurie(curie)) return h('span', { class: 'mono' }, curie || '–');
  return h('a', { class: 'curie mono', href: olsUrl(curie), target: '_blank', rel: 'noopener noreferrer', title: `${curie} in the Ontology Lookup Service`, onclick: (e) => e.stopPropagation() }, label || curie);
}
