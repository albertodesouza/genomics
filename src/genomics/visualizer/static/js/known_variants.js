// Known-variant annotation for the Gene products page: rsID, clinical significance, associated traits
// and literature of a variant (GET /products/known), with links to the reference databases.
import { h } from './ui.js';
import { dbsnpUrl } from './links.js';

const enc = encodeURIComponent;
const MAX_PUBMED_LINK = 150;

/** "chr16" + 89919709 C>T -> links to every database that holds a record of the variant. */
export function knownLinks(rec, { chrom, pos, ref, alt }) {
  const c = String(chrom || '').replace(/^chr/, '');
  const links = [{ label: 'dbSNP', url: dbsnpUrl(rec.rsid), title: `${rec.rsid} on NCBI dbSNP: the reference record of this variant`, primary: true }];
  for (const vcv of rec.clinvar || []) links.push({ label: 'ClinVar', url: `https://www.ncbi.nlm.nih.gov/clinvar/variation/${Number(String(vcv).replace(/^VCV0*/, ''))}/`, title: `ClinVar ${vcv}: clinical interpretations and conditions` });
  for (const omim of rec.omim || []) {
    const [entry, allele] = String(omim).split('.');
    links.push({ label: 'OMIM', url: `https://www.omim.org/entry/${entry}${allele ? `#${allele}` : ''}`, title: `OMIM allelic variant ${omim}: curated description and the studies behind it` });
  }
  if (rec.traits && rec.traits.some((t) => (t.sources || []).some((s) => /GWAS/.test(s)))) links.push({ label: 'GWAS Catalog', url: `https://www.ebi.ac.uk/gwas/variants/${enc(rec.rsid)}`, title: 'NHGRI-EBI GWAS Catalog associations of this variant' });
  links.push({ label: `Literature${rec.pubmed_count ? ` (${rec.pubmed_count})` : ''}`, url: `https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=${enc(`litvar@${rec.rsid}##`)}&query=${enc(rec.rsid)}`, title: 'LitVar: every article that mentions this variant, searchable' });
  if (rec.pubmed && rec.pubmed.length) links.push({ label: 'PubMed', url: `https://pubmed.ncbi.nlm.nih.gov/?term=${enc(rec.pubmed.slice(0, MAX_PUBMED_LINK).join(' OR '))}`, title: `The ${Math.min(rec.pubmed.length, MAX_PUBMED_LINK)} articles Ensembl links to this variant` });
  for (const v of rec.uniprot || []) links.push({ label: 'UniProt', url: `https://web.expasy.org/variant_pages/${enc(v)}.html`, title: `UniProt/Swiss-Prot variant ${v}` });
  for (const p of rec.pharmgkb || []) links.push({ label: 'ClinPGx', url: `https://www.clinpgx.org/variant/${enc(p)}`, title: 'Pharmacogenomics (ClinPGx, formerly PharmGKB)' });
  if (c && pos && ref && alt) links.push({ label: 'gnomAD', url: `https://gnomad.broadinstitute.org/variant/${c}-${pos}-${ref}-${alt}?dataset=gnomad_r4`, title: 'gnomAD: allele frequencies by population' });
  links.push({ label: 'Ensembl', url: `https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${enc(rec.rsid)}`, title: 'Ensembl variant page: consequences, phenotypes, citations' });
  return links;
}

const linkRow = (links) => h('div', { class: 'known-links' }, links.map((l) => h('a', { href: l.url, target: '_blank', rel: 'noopener', title: l.title, class: l.primary ? 'primary' : '' }, l.label)));

/** A compact rsID link (to dbSNP) with a "known" mark when the variant has clinical or literature records. */
export function knownChip(rec) {
  if (!rec) return h('span', { class: 'muted' }, '–');
  const tip = [rec.rsid, rec.clinical_significance.length ? `ClinVar: ${rec.clinical_significance.join(', ')}` : '',
    (rec.traits || []).slice(0, 5).map((t) => t.trait).join('; '), rec.pubmed_count ? `${rec.pubmed_count} publications` : ''].filter(Boolean).join('\n');
  return h('span', { class: 'known-chip' },
    h('a', { class: 'mono', href: dbsnpUrl(rec.rsid), target: '_blank', rel: 'noopener', title: tip, onclick: (e) => e.stopPropagation() }, rec.rsid),
    rec.known ? h('span', { class: 'pill product-class warn', title: tip }, 'known') : null);
}

const fmtP = (p) => (p === null || p === undefined ? '' : p < 1e-300 ? 'p < 1e-300' : `p = ${p.toExponential(0)}`);
const freq = (f) => (f === null || f === undefined ? null : f >= 0.01 ? `${(f * 100).toFixed(1)}%` : `${(f * 100).toFixed(2)}%`);

/** The full record of a known variant: what it is associated with and where that is documented. */
export function knownPanel(rec, variant, { change = null } = {}) {
  const clin = rec.clinical_significance || [];
  const traits = rec.traits || [];
  const clinvarTraits = traits.filter((t) => (t.sources || []).includes('ClinVar'));
  const gwasTraits = traits.filter((t) => !(t.sources || []).includes('ClinVar'));
  const f = rec.frequency || {};
  const freqText = [freq(f.gnomad_genomes) && `gnomAD ${freq(f.gnomad_genomes)}`, freq(f['1000g_eur']) && `1000G EUR ${freq(f['1000g_eur'])}`,
    freq(f['1000g_afr']) && `AFR ${freq(f['1000g_afr'])}`, freq(f['1000g_eas']) && `EAS ${freq(f['1000g_eas'])}`].filter(Boolean).join(' · ');
  return h('div', { class: 'known-panel' },
    h('div', { class: 'known-head' },
      h('a', { class: 'known-rsid mono', href: dbsnpUrl(rec.rsid), target: '_blank', rel: 'noopener', title: 'Open the dbSNP record' }, rec.rsid),
      change ? h('span', { class: 'mono' }, change) : null,
      h('span', { class: 'muted small mono' }, `${variant.chrom}:${variant.pos.toLocaleString('en-US')} ${variant.ref}>${variant.alt}`),
      rec.same_allele ? null : h('span', { class: 'pill muted', title: 'dbSNP records this position with other alleles; check the record' }, 'same position, other allele?'),
      freqText ? h('span', { class: 'muted small' }, freqText) : null),
    clin.length ? h('div', { class: 'known-row' }, h('span', { class: 'known-label' }, 'ClinVar'), h('span', null, ...clin.map((c) => h('span', { class: `pill ${/pathogenic/.test(c) && !/benign/.test(c) ? 'product-class bad' : 'muted'}` }, c)))) : null,
    clinvarTraits.length ? h('div', { class: 'known-row' }, h('span', { class: 'known-label' }, 'Conditions'),
      h('span', null, clinvarTraits.map((t) => t.trait).join(' · '))) : null,
    gwasTraits.length ? h('div', { class: 'known-row' }, h('span', { class: 'known-label' }, 'GWAS traits'),
      h('span', null, ...gwasTraits.slice(0, 8).map((t, i) => h('span', null, i ? ' · ' : '', t.trait, t.best_p !== null ? h('span', { class: 'muted small' }, ` (${fmtP(t.best_p)})`) : null)),
        rec.trait_count > traits.length ? h('span', { class: 'muted small' }, ` · +${rec.trait_count - traits.length} more on the GWAS Catalog`) : null)) : null,
    h('div', { class: 'known-row' }, h('span', { class: 'known-label' }, 'References'), linkRow(knownLinks(rec, variant))));
}
