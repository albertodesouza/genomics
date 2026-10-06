// What the locus box accepts, shared by the Tracks, Sequence and Perturbation Lab pages.
//
// Coordinates are resolved by the page (each has its own coordinate systems and windows); this
// module turns the other two forms into a genomic position:
//   · a gene symbol  -> the gene annotated in the current window, else a window of the dataset
//   · an rsID        -> its GRCh38 position through GTEx (cached server-side)
import { api } from './api.js';

const COORDS_RE = /^(?:([\w.]+):)?\s*([\d,]+)(?:\s*[-–]\s*([\d,]+))?$/;
const RSID_RE = /^rs\d+$/i;
const SYMBOL_RE = /^[A-Za-z][A-Za-z0-9._-]{0,30}$/;
const GENE_PAD = 0.1; // show a gene with a tenth of its length of context on each side

/** What the box understands, for the placeholder and the hint. */
export const LOCUS_HINT = 'chr:start-end, a position, a gene name (OCA2) or an rsID (rs12913832)';

export function isCoordinates(text) {
  return COORDS_RE.test(String(text || '').trim());
}

/**
 * Resolve a gene symbol or rsID typed in the locus box.
 *
 * ``context``: {datasetId, gene, window: {chromosome, start, end}, annotations: [{name, start, end}],
 * genes: [window names]}. Offsets in ``annotations`` are relative to the window, as the pages hold them.
 *
 * Returns one of
 *   {kind: 'range', start, end}        offsets in the current window
 *   {kind: 'window', gene, ...}        another window of the dataset; 'position' when known
 *   {kind: 'position', pos}            a genomic position inside the current window
 * or throws an Error whose message is meant for a toast.
 */
export async function resolveLocus(text, context) {
  const query = String(text || '').trim();
  if (!query) throw new Error('Type a position, a gene name or an rsID');
  if (RSID_RE.test(query)) return resolveRsid(query, context);
  if (!SYMBOL_RE.test(query)) throw new Error(`Use ${LOCUS_HINT}`);
  return resolveSymbol(query, context);
}

function resolveSymbol(symbol, context) {
  const upper = symbol.toUpperCase();
  const annotated = (context.annotations || []).filter((g) => String(g.name || '').toUpperCase() === upper);
  if (annotated.length) {
    // Several transcripts' gene rows can share a name; show them all.
    const start = Math.min(...annotated.map((g) => g.start));
    const end = Math.max(...annotated.map((g) => g.end));
    const pad = Math.max(50, Math.round((end - start) * GENE_PAD));
    return { kind: 'range', start: start - pad, end: end + pad, label: symbol };
  }
  const window = (context.genes || []).find((g) => String(g).toUpperCase() === upper);
  if (window) {
    if (window === context.gene) throw new Error(`${window} is this window, but it is not in the gene annotations`);
    return { kind: 'window', gene: window, label: window };
  }
  const near = (context.genes || []).filter((g) => String(g).toUpperCase().startsWith(upper)).slice(0, 3);
  throw new Error(
    `No gene "${symbol}" in this window's annotations or in the dataset's windows`
    + (near.length ? `. Did you mean ${near.join(', ')}?` : ''),
  );
}

async function resolveRsid(rsid, context) {
  let info;
  try {
    info = await api('/api/gtex/resolve', { params: { rsid: rsid.toLowerCase() } });
  } catch (err) {
    throw new Error(`${rsid}: ${err.message}`);
  }
  const window = context.window || {};
  const sameChromosome = normChrom(info.chromosome) === normChrom(window.chromosome);
  if (sameChromosome && window.start && info.pos >= window.start && info.pos <= window.end) {
    return { kind: 'position', pos: info.pos, label: `${rsid} (${info.chromosome}:${info.pos})` };
  }
  // The dataset's other windows are the only other places this page can show it.
  const elsewhere = (context.windows || []).find((w) => normChrom(w.chromosome) === normChrom(info.chromosome) && info.pos >= w.start && info.pos <= w.end);
  if (elsewhere) return { kind: 'window', gene: elsewhere.gene, position: info.pos, label: `${rsid} (${elsewhere.gene})` };
  throw new Error(`${rsid} is at ${info.chromosome}:${info.pos}, outside every window of this dataset`);
}

function normChrom(value) {
  const text = String(value || '').trim().toLowerCase();
  return text.startsWith('chr') ? text.slice(3) : text;
}
