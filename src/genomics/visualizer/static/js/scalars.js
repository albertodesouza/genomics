// Region scalars: define "one number per sample" from a track over a region; it becomes a sample field.
import { api, apiJob } from './api.js';
import { state, ds, geneInfo, reloadSamples, filterRows } from './state.js';
import { h, clear, select, field, fmtInt, fmtBp, fmtNum, toast, errorBox, jobOverlay, modal, setOptions } from './ui.js';
import { loadGeneAnnotations } from './gene_lanes.js';
import { setupCanvas, theme } from './plot.js';

export const OUTPUT_LABELS = { rna_seq: 'RNA-seq', cage: 'CAGE', procap: 'PRO-cap', dnase: 'DNase', atac: 'ATAC', chip_histone: 'ChIP histone', chip_tf: 'ChIP TF', splice_sites: 'Splice sites', splice_site_usage: 'Splice usage' };
const STAT_TEXT = { mean: 'mean', sum: 'sum' };

export async function loadScalars() {
  return (await api(`${ds()}/scalars`)).scalars || [];
}

/** One line describing a scalar, e.g. "CAGE · melanocyte (+) · sum over chr15:28,000,000-28,001,000 · log2(x+1)". */
export function scalarText(s) {
  const parts = [`${OUTPUT_LABELS[s.output] || s.output} · ${s.track_label}`, `${STAT_TEXT[s.stat] || s.stat} over ${s.region || `${s.gene} ${fmtInt(s.start)}-${fmtInt(s.end)}`}`];
  if (s.transform === 'log2p1') parts.push('log2(x+1)');
  return parts.join(' · ');
}

const slugName = (text) => text.toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '').replace(/^(\d)/, 's$1').slice(0, 48) || 'scalar';

/**
 * Modal to define a region scalar. ``prefill``: { gene, output, track, start, end } (window offsets) from
 * the Tracks view; without it the form starts from the first window. Resolves after saving.
 */
export function openScalarForm(prefill = {}, { onSaved } = {}) {
  const genes = state.summary.genes.map((g) => g.gene);
  const cfg = {
    gene: genes.includes(prefill.gene) ? prefill.gene : genes[0],
    output: prefill.output || null, track: Number.isInteger(prefill.track) ? prefill.track : null,
    mode: prefill.start !== undefined ? 'custom' : 'tss', annot: null, flank: 500, custom: '',
    stat: 'mean', transform: 'none', bins: '3',
  };
  let info = null; let annotations = []; let nameTouched = false;
  const name = h('input', { class: 'input mono', placeholder: 'e.g. cage_tyrp1_tss', oninput: () => { nameTouched = true; } });
  const geneSel = select(genes, cfg.gene, (v) => { cfg.gene = v; cfg.output = null; cfg.track = null; cfg.annot = null; cfg.mode = cfg.mode === 'custom' ? 'tss' : cfg.mode; modeSel.value = cfg.mode; loadGene(); });
  const outputSel = select([], null, (v) => { cfg.output = v; cfg.track = 0; renderTracks(); update(); });
  const trackSel = select([], null, (v) => { cfg.track = Number(v); update(); }, { style: { maxWidth: '420px' } });
  const modeSel = select([
    { value: 'tss', label: 'TSS ± flank of a gene' }, { value: 'body', label: 'Gene body' }, { value: 'custom', label: 'Custom range' },
  ], cfg.mode, (v) => { cfg.mode = v; renderRegion(); update(); });
  const annotSel = select([], null, (v) => { cfg.annot = v; update(); });
  const flank = h('input', { class: 'input', type: 'number', min: '0', step: '100', value: String(cfg.flank), oninput: () => { cfg.flank = Math.max(0, Number(flank.value) || 0); update(); } });
  const custom = h('input', { class: 'input mono', placeholder: 'chr15:28,000,000-28,001,000', oninput: () => { cfg.custom = custom.value; update(); } });
  const statSel = select([{ value: 'mean', label: 'Mean over the region' }, { value: 'sum', label: 'Sum over the region (mean × length)' }], cfg.stat, (v) => { cfg.stat = v; update(); });
  const transformSel = select([{ value: 'none', label: 'None' }, { value: 'log2p1', label: 'log2(x + 1)' }], cfg.transform, (v) => { cfg.transform = v; update(); });
  const binsSel = select([{ value: '0', label: 'No bins' }, { value: '2', label: 'Halves' }, { value: '3', label: 'Tertiles' }, { value: '4', label: 'Quartiles' }, { value: '5', label: 'Quintiles' }], cfg.bins, (v) => { cfg.bins = v; update(); });
  const annotField = field('Gene', annotSel);
  const flankField = field('Flank (bp)', flank);
  const customField = field('Genomic range', custom);
  const preview = h('div', { class: 'muted', style: { fontSize: '12.5px' } });
  const body = h('div', { style: { display: 'grid', gap: '12px', position: 'relative', minWidth: 'min(560px, 80vw)' } },
    h('p', { class: 'secondary', style: { margin: 0 } }, 'One number per sample: the track averaged over the region on each haplotype, then over H1 and H2. It becomes a column on Samples and, with bins, a filter and a Group-by field on Tracks.'),
    h('div', { class: 'grid cols-2', style: { gap: '10px' } }, field('Window', geneSel), field('Output', outputSel)),
    field('Track', trackSel),
    h('div', { class: 'grid cols-2', style: { gap: '10px' } }, field('Region', modeSel), annotField, flankField, customField),
    h('div', { class: 'grid cols-2', style: { gap: '10px' } }, field('Statistic', statSel), field('Transform', transformSel), field('Quantile bins (all samples)', binsSel), field('Field name', name)),
    preview);

  function range() {
    if (!info) return null;
    const clamp = (a, b) => [Math.max(0, Math.round(a)), Math.min(info.length, Math.round(b))];
    if (cfg.mode === 'custom') {
      const m = /(?:(chr\w+):)?([\d,]+)\s*-\s*([\d,]+)/.exec(cfg.custom || '');
      if (!m) return null;
      return clamp(Number(m[2].replace(/,/g, '')) - info.start, Number(m[3].replace(/,/g, '')) - info.start + 1);
    }
    const g = annotations.find((x) => x.name === cfg.annot);
    if (!g) return null;
    if (cfg.mode === 'body') return clamp(g.start, g.end);
    const tss = g.strand === '-' ? g.end : g.start;
    return clamp(tss - cfg.flank, tss + cfg.flank + 1);
  }

  function regionText(r) { return `${info.chromosome}:${fmtInt(info.start + r[0])}-${fmtInt(info.start + r[1] - 1)}`; }

  function suggestedName() {
    const where = cfg.mode === 'custom' ? 'region' : `${cfg.annot || cfg.gene}_${cfg.mode === 'tss' ? 'tss' : 'body'}`;
    return slugName(`${cfg.output || ''}_${where}`);
  }

  function update() {
    if (!nameTouched) name.value = suggestedName();
    const r = range();
    clear(preview);
    if (!info) return;
    if (!r || r[1] <= r[0]) { preview.append(cfg.mode === 'custom' ? `Enter a range inside ${info.chromosome}:${fmtInt(info.start)}-${fmtInt(info.end)}.` : 'Pick a gene.'); return; }
    preview.append(`${regionText(r)} (${fmtBp(r[1] - r[0])}) · first use reads every sample's predictions (minutes on a full cohort); later scalars over the same region and output are instant.`);
  }

  function renderTracks() {
    const out = (info.outputs || {})[cfg.output] || { tracks: [] };
    setOptions(trackSel, out.tracks.map((t, i) => ({ value: String(i), label: t.label || `track ${i}` })), String(cfg.track ?? 0));
    if (cfg.track === null || cfg.track >= out.tracks.length) cfg.track = 0;
  }

  function renderRegion() {
    annotField.hidden = cfg.mode === 'custom';
    flankField.hidden = cfg.mode !== 'tss';
    customField.hidden = cfg.mode !== 'custom';
  }

  async function loadGene() {
    info = await geneInfo(cfg.gene);
    const outputs = Object.keys(info.outputs || {});
    if (!outputs.includes(cfg.output)) cfg.output = outputs.includes('cage') ? 'cage' : outputs[0];
    setOptions(outputSel, outputs.map((o) => ({ value: o, label: OUTPUT_LABELS[o] || o })), cfg.output);
    renderTracks();
    try { annotations = ((await loadGeneAnnotations(state.datasetId, cfg.gene)).genes || []).filter((g) => g.end > 0 && g.start < info.length); } catch (err) { annotations = []; }
    const names = annotations.map((g) => g.name);
    if (!names.includes(cfg.annot)) cfg.annot = names.includes(cfg.gene) ? cfg.gene : names[0] || null;
    setOptions(annotSel, names, cfg.annot);
    if (!annotations.length && cfg.mode !== 'custom') { cfg.mode = 'custom'; modeSel.value = 'custom'; }
    if (prefill.start !== undefined && cfg.gene === prefill.gene && !cfg.custom) {
      cfg.custom = regionText([prefill.start, prefill.end]);
      custom.value = cfg.custom;
    }
    renderRegion();
    update();
  }

  const save = h('button', { class: 'btn primary', onclick: submit }, 'Compute');
  const m = modal('New region scalar', body, [h('button', { class: 'btn', onclick: () => m.close() }, 'Cancel'), save]);

  async function submit() {
    const r = range();
    if (!info || !r || r[1] <= r[0]) { toast('Choose a region inside the window', 'error'); return; }
    const track = ((info.outputs[cfg.output] || {}).tracks || [])[cfg.track] || {};
    const payload = {
      name: name.value.trim(), gene: cfg.gene, output: cfg.output, track: cfg.track, track_label: track.label || `track ${cfg.track}`,
      start: r[0], end: r[1], region: cfg.mode === 'custom' ? regionText(r) : `${cfg.mode === 'tss' ? `${cfg.annot} TSS ± ${fmtInt(cfg.flank)} bp` : `${cfg.annot} gene body`} (${regionText(r)})`,
      stat: cfg.stat, transform: cfg.transform, bins: Number(cfg.bins),
    };
    save.disabled = true;
    const overlay = jobOverlay(body, `${cfg.gene}: ${OUTPUT_LABELS[cfg.output] || cfg.output} per haplotype over the region`);
    try {
      const res = await apiJob(`${ds()}/scalars`, { method: 'POST', body: payload, onProgress: (job) => overlay.update(job) });
      await reloadSamples();
      m.close();
      toast(`Added ${res.scalar.name}${res.scalar.bin_field ? ` and ${res.scalar.bin_field}` : ''} to the samples`);
      if (onSaved) onSaved(res.scalar);
    } catch (err) {
      toast(err.message, 'error', 8000);
    } finally { overlay.remove(); save.disabled = false; }
  }

  loadGene().catch((err) => body.appendChild(errorBox(err)));
  return m;
}

export async function deleteScalar(scalar) {
  await api(`${ds()}/scalars/delete`, { method: 'POST', body: { name: scalar.name } });
  await reloadSamples();
}

/** Histogram of a scalar: all samples (grey) and the current cohort (accent). */
export function drawScalarHistogram(canvas, scalar, { width = 240, height = 64 } = {}) {
  const hist = scalar.histogram || {};
  const edges = hist.edges || [];
  const ctx = setupCanvas(canvas, width, height);
  const t = theme();
  if (edges.length < 2) return;
  const bins = edges.length - 1;
  const lo = edges[0]; const hi = edges[bins];
  const col = state.samples.col[scalar.name];
  const cohort = new Array(bins).fill(0);
  if (col !== undefined) {
    for (const row of filterRows(state.filters, '')) {
      const v = row[col];
      if (typeof v !== 'number') continue;
      cohort[Math.min(bins - 1, Math.max(0, Math.floor(((v - lo) / (hi - lo)) * bins)))] += 1;
    }
  }
  const max = Math.max(1, ...hist.counts);
  const bw = width / bins;
  const bottom = height - 12;
  for (let i = 0; i < bins; i++) {
    const ha = (hist.counts[i] / max) * (bottom - 2);
    const hc = (cohort[i] / max) * (bottom - 2);
    ctx.fillStyle = t.line; ctx.fillRect(i * bw + 0.5, bottom - ha, Math.max(1, bw - 1), ha);
    ctx.fillStyle = t.accent; ctx.fillRect(i * bw + 0.5, bottom - hc, Math.max(1, bw - 1), hc);
  }
  ctx.fillStyle = t.ink3; ctx.font = '10px system-ui, sans-serif'; ctx.textBaseline = 'bottom';
  ctx.textAlign = 'left'; ctx.fillText(fmtNum(lo), 0, height);
  ctx.textAlign = 'right'; ctx.fillText(fmtNum(hi), width, height);
  ctx.strokeStyle = t.ink2; ctx.setLineDash([3, 3]);
  for (const e of scalar.edges || []) { const x = ((e - lo) / (hi - lo)) * width; ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, bottom); ctx.stroke(); }
  ctx.setLineDash([]);
}
