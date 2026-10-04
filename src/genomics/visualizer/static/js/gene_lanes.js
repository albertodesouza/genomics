// Gene annotation lanes (GENCODE models: exons, CDS, strand chevrons, labels) shared by the
// Tracks, Sequence and Perturbation pages. Gene coordinates are reference (window) offsets; pages
// pass ``toPos`` to map them onto their own axis (e.g. the training alignment axis).
import { apiJob } from './api.js';

export const GENE_ROW_H = 24;
export const MAX_GENE_ROWS = 4;

const cache = new Map();

/** Annotated genes overlapping a dataset window (fetched once per dataset and window). */
export function loadGeneAnnotations(datasetId, gene, { onProgress } = {}) {
  const key = `${datasetId}\u0000${gene}`;
  if (!cache.has(key)) {
    const url = `/api/d/${encodeURIComponent(datasetId)}/genes/${encodeURIComponent(gene)}/annotations`;
    cache.set(key, apiJob(url, { onProgress }).catch((err) => { cache.delete(key); throw err; }));
  }
  return cache.get(key);
}

/** Pack the genes overlapping the view into rows, leaving room for each label. */
export function packGenes(genes, { start, end, pxPerBp, toPos = (p) => p }) {
  const rows = [];
  for (const gene of genes || []) {
    const a = toPos(gene.start);
    const b0 = toPos(gene.end);
    if (b0 <= start || a >= end) continue;
    const b = b0 + (gene.name.length * 7 + 12) / pxPerBp;
    const row = rows.find((r) => r.end < a);
    if (row) { row.items.push(gene); row.end = b; } else rows.push({ end: b, items: [gene] });
  }
  return rows.map((r) => r.items);
}

/** Height of the lanes for packed ``rows`` (0 when there are no genes in view). */
export function geneLanesHeight(rows, maxRows = MAX_GENE_ROWS) {
  return rows.length ? Math.min(rows.length, maxRows) * GENE_ROW_H + 6 : 0;
}

/** One gene model (first transcript) with its label beside it or pinned at the view's left edge. */
export function drawGeneModel(ctx, t, gene, y, { xOf, left, right, highlight }) {
  const tx = gene.transcripts && gene.transcripts[0];
  const color = gene.name === highlight ? t.accent : t.ink2;
  const a = xOf(gene.start); const b = xOf(gene.end);
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(a, y + 6.5); ctx.lineTo(b, y + 6.5); ctx.stroke();
  const dir = gene.strand === '-' ? -1 : 1;
  for (let cx = Math.max(a, left) + 10; cx < Math.min(b, right) - 4; cx += 26) {
    ctx.beginPath(); ctx.moveTo(cx - 2 * dir, y + 4); ctx.lineTo(cx + 1 * dir, y + 6.5); ctx.lineTo(cx - 2 * dir, y + 9); ctx.stroke();
  }
  if (tx) {
    for (const [s, e] of tx.exons) { const xa = xOf(s); ctx.fillRect(xa, y + 3, Math.max(1, xOf(e) - xa), 7); }
    for (const [s, e] of tx.cds) { const xa = xOf(s); ctx.fillRect(xa, y + 1, Math.max(1, xOf(e) - xa), 11); }
  } else {
    ctx.fillRect(a, y + 3, Math.max(1, b - a), 7);
  }
  ctx.font = `600 11px ${t.font}`; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
  const label = `${gene.name} ${gene.strand === '-' ? '←' : '→'}`;
  const lw = ctx.measureText(label).width;
  let lx = b + 5;
  if (lx + lw > right) lx = Math.max(left + 4, a + 4);
  if (lx + lw <= right) {
    ctx.fillStyle = t.surface; ctx.fillRect(lx - 3, y, lw + 6, 13);
    ctx.fillStyle = color; ctx.fillText(label, lx, y + 6.5);
  }
}

/**
 * Draw packed gene rows starting at ``top`` (a "genes" label sits in the left gutter).
 * ``xOf`` maps a reference offset to a canvas x.
 */
export function drawGeneLanes(ctx, t, rows, { top, left, width, xOf, highlight, maxRows = MAX_GENE_ROWS, label = 'genes' }) {
  const height = geneLanesHeight(rows, maxRows);
  if (!height) return 0;
  ctx.font = `10.5px ${t.font}`; ctx.textAlign = 'right'; ctx.textBaseline = 'alphabetic'; ctx.fillStyle = t.ink3;
  ctx.fillText(label, left - 8, top + 16);
  ctx.save();
  ctx.beginPath(); ctx.rect(left, top, width, height); ctx.clip();
  rows.slice(0, maxRows).forEach((row, r) => {
    for (const gene of row) drawGeneModel(ctx, t, gene, top + 6 + r * GENE_ROW_H, { xOf, left, right: left + width, highlight });
  });
  ctx.restore();
  return height;
}

/** The gene under canvas point (x, y) of lanes drawn at ``top`` (for tooltips), or null. */
export function geneAt(rows, x, y, { top, xOf, maxRows = MAX_GENE_ROWS }) {
  const r = Math.floor((y - top - 6 + 4) / GENE_ROW_H);
  if (r < 0 || r >= Math.min(rows.length, maxRows)) return null;
  return rows[r].find((g) => x >= xOf(g.start) - 3 && x <= xOf(g.end) + 3) || null;
}

export function geneTooltip(gene, escapeHtml) {
  const tx = gene.transcripts && gene.transcripts[0];
  const rows = [
    gene.type ? `<div class="tt-row"><span class="tt-key">type</span><span class="tt-val">${escapeHtml(String(gene.type).replace(/_/g, ' '))}</span></div>` : '',
    `<div class="tt-row"><span class="tt-key">strand</span><span class="tt-val">${escapeHtml(gene.strand || '')}</span></div>`,
    tx ? `<div class="tt-row"><span class="tt-key">exons</span><span class="tt-val">${tx.exons.length}</span></div>` : '',
    tx ? `<div class="tt-sub">${escapeHtml(tx.name || tx.id || '')}${tx.id && tx.id !== tx.name ? ` · ${escapeHtml(tx.id)}` : ''}</div>` : '',
  ];
  return `<div class="tt-title">${escapeHtml(gene.name)}</div>${rows.join('')}`;
}
