// Overview: dataset summary, cohort composition, genes and outputs; open other datasets.
import { api } from '../api.js';
import { navigate } from '../app.js';
import { state, setDataset, setFilters, emit, loadStatus, geneInfo, categoricalFields } from '../state.js';
import { h, clear, icon, fmtInt, fmtBp, fmtDate, select, toast, errorBox } from '../ui.js';
import { openPredictForm, openTrainForm } from '../forms.js';
import { loadSessions, sessionList } from '../sessions.js';
import { geneNames, openGeneCard, openTrackCard } from '../cards.js';
import { assayTerm, curieLink } from '../links.js';

export async function mount(root) {
  const page = h('div', { class: 'page-inner' });
  root.appendChild(page);
  if (!state.summary) {
    page.append(h('div', { class: 'page-head' }, h('div', null, h('h1', null, 'No dataset loaded'), h('p', null, 'Open a dataset directory, or import one: Import → Quickstart builds a 1000 Genomes or gnomAD HGDP dataset straight from the public servers.'))), datasetsCard());
    return {};
  }
  const s = state.summary;
  const genesWithData = s.genes.filter((g) => g.length);
  page.append(
    h('div', { class: 'page-head' },
      h('div', null,
        h('h1', null, s.name),
        h('p', null, h('span', { class: 'mono' }, s.path), s.last_updated ? h('span', { class: 'muted' }, ` · updated ${fmtDate(s.last_updated)}`) : null)),
      h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap' } },
        h('button', { class: 'btn', onclick: () => navigate('samples') }, 'Browse samples'),
        h('button', { class: 'btn', title: 'Predict tracks for every individual (background job)', onclick: () => openPredictForm() }, icon('alphagenome', 15), 'AlphaGenome predictions'),
        h('button', { class: 'btn', title: 'Train and evaluate a model on this dataset (background job)', onclick: () => openTrainForm() }, icon('play', 14), 'Train a model'),
        h('button', { class: 'btn primary', onclick: () => navigate('tracks') }, 'Open tracks'))),
    runningNotice(),
    tiles(s),
    h('div', { class: 'grid cols-2', style: { marginTop: '16px' } }, compositionCard(), outputsCard(s)),
    h('div', { style: { marginTop: '16px' } }, sessionsCard()),
    h('div', { style: { marginTop: '16px' } }, genesCard(s, genesWithData)),
    h('div', { style: { marginTop: '16px' } }, datasetsCard()),
  );
  return {};
}

function runningNotice() {
  const host = h('div');
  api('/api/tasks').then(({ tasks }) => {
    const mine = tasks.filter((t) => ['starting', 'queued', 'running'].includes(t.status) && [t.params.dataset_path, t.params.output_dir].includes(state.summary.path));
    if (!mine.length) return;
    host.appendChild(h('div', { class: 'notice', style: { marginBottom: '12px', display: 'flex', gap: '10px', alignItems: 'center' } },
      h('span', { class: 'spinner' }), h('span', null, mine.map((t) => `${t.title}: ${Math.round((t.progress || 0) * 100)}%`).join(' · ')), h('a', { href: '#/jobs' }, 'Jobs →')));
  }).catch(() => {});
  return host;
}

function tiles(s) {
  const tile = (label, value, sub) => h('div', { class: 'tile' }, h('div', { class: 'tile-label' }, label), h('div', { class: 'tile-value' }, value), sub ? h('div', { class: 'tile-sub' }, sub) : null);
  const facets = categoricalFields();
  const outputs = s.outputs.length ? s.outputs.join(', ') : 'discovered per gene';
  return h('div', { class: 'tiles' },
    tile('Individuals', fmtInt(s.sample_count), facets.length ? `${facets.length} sample facets: ${facets.map((f) => f.label.toLowerCase()).join(', ')}` : 'no pedigree facets'),
    tile('Gene windows', fmtInt(s.gene_count), `${s.genes.filter((g) => g.in_metadata).length} listed in dataset metadata`),
    tile('Window size', s.window_size ? fmtBp(s.window_size) : '–', 'AlphaGenome input length'),
    tile('Outputs', s.outputs.length ? fmtInt(s.outputs.length) : '–', outputs));
}

/** Saved views of this dataset: click one to restore its locus, tracks, cohort and pinned samples. */
function sessionsCard() {
  const body = h('div', { class: 'card-body' }, h('div', { class: 'muted', style: { fontSize: '12.5px' } }, 'Loading…'));
  const card = h('section', { class: 'card' },
    h('div', { class: 'card-head' }, h('h2', null, 'Saved views'),
      h('span', { class: 'muted', style: { fontSize: '12px' } }, 'locus, tracks, cohort and pinned individuals')),
    body);
  const refresh = () => loadSessions().then((list) => clear(body).appendChild(sessionList(list, { onChanged: refresh })));
  refresh();
  return card;
}

function compositionCard() {
  const fields = categoricalFields();
  const body = h('div', { class: 'card-body' });
  const card = h('section', { class: 'card' });
  if (!fields.length) {
    card.append(h('div', { class: 'card-head' }, h('h2', null, 'Cohort composition')), h('div', { class: 'card-body empty' }, 'No categorical sample fields in this dataset.'));
    return card;
  }
  let fieldName = fields.find((f) => f.name === 'superpopulation') ? 'superpopulation' : fields[0].name;
  const render = () => {
    clear(body);
    const f = fields.find((x) => x.name === fieldName);
    const max = Math.max(...f.counts.map((c) => c.count), 1);
    const shown = f.counts.slice(0, 30);
    body.appendChild(h('div', { class: 'hbars', role: 'list' }, shown.map((c) =>
      h('div', { class: 'hbar', role: 'listitem', title: `Show ${c.value} samples`, onclick: () => { setFilters({ [fieldName]: [c.value] }); navigate('samples'); } },
        h('span', { class: 'hb-name' }, c.value),
        h('span', { class: 'hb-track' }, h('span', { class: 'hb-fill', style: { width: `${(c.count / max) * 100}%` } })),
        h('span', { class: 'hb-val' }, fmtInt(c.count))))));
    if (f.counts.length > shown.length) body.appendChild(h('div', { class: 'muted', style: { marginTop: '8px', fontSize: '12px' } }, `+ ${f.counts.length - shown.length} more values`));
    if (f.missing) body.appendChild(h('div', { class: 'muted', style: { marginTop: '8px', fontSize: '12px' } }, `${fmtInt(f.missing)} samples without a value`));
  };
  card.append(h('div', { class: 'card-head' }, h('h2', null, 'Cohort composition'), select(fields.map((f) => ({ value: f.name, label: f.label })), fieldName, (v) => { fieldName = v; render(); })), body);
  render();
  return card;
}

function outputsCard(s) {
  const body = h('div', { class: 'card-body' }, h('div', { class: 'muted' }, 'Loading track metadata…'));
  const card = h('section', { class: 'card' }, h('div', { class: 'card-head' }, h('h2', null, 'Outputs & tracks')), body);
  const gene = s.genes.find((g) => g.in_metadata) || s.genes[0];
  if (!gene) { clear(body).appendChild(h('div', { class: 'empty' }, 'No genes')); return card; }
  geneInfo(gene.gene).then((info) => {
    clear(body);
    const outputs = [...Object.entries(info.outputs || {}), ...Object.entries(info.other_outputs || {})];
    if (!outputs.length) {
      body.appendChild(h('div', { class: 'empty', style: { display: 'grid', gap: '10px', justifyItems: 'center' } }, `No AlphaGenome predictions found for ${gene.gene}.`,
        h('button', { class: 'btn primary', onclick: () => openPredictForm() }, icon('alphagenome', 15), 'Run AlphaGenome predictions')));
      return;
    }
    const shape = (out) => {
      if (out.kind === 'junctions') return `${fmtInt(out.junctions)} junctions`;
      if (out.kind === 'contact_map') return `${out.bins}×${out.bins} contact map, ${fmtBp(out.resolution)} bins`;
      return `${fmtBp(out.length)}${out.resolution > 1 ? ` in ${out.resolution} bp bins` : ''}`;
    };
    for (const [name, out] of outputs) {
      const plotted = (info.outputs || {})[name];
      body.appendChild(h('details', { class: 'collapsible', style: { marginBottom: '10px' }, open: outputs.length <= 2 && out.tracks.length <= 40 },
        h('summary', { style: { display: 'flex', justifyContent: 'space-between', gap: '12px' } }, h('span', null, name, plotted ? null : h('span', { class: 'pill', style: { marginLeft: '8px' }, title: 'Stored with the predictions; not drawn by the track browser' }, out.tracks.length ? 'not a 1-D track' : 'no tracks')),
          h('span', { class: 'muted', style: { fontWeight: 400 } }, `${fmtInt(out.tracks.length)} tracks · ${shape(out)} · haplotypes ${info.haplotypes.join(', ')}`)),
        h('div', { class: 'table-wrap', style: { maxHeight: '280px' } }, h('table', { class: 'table' },
          h('thead', null, h('tr', null, h('th', null, '#'), h('th', null, 'Biosample'), h('th', null, 'Ontology'), h('th', null, 'Strand'), h('th', null, 'Assay / target'))),
          h('tbody', null, out.tracks.map((t) => {
            const assay = t.metadata['Assay title'] || null;
            const term = assayTerm(assay || name);
            return h('tr', { class: 'clickable', title: 'Track card: ontology terms, assay, observed experiments', onclick: () => openTrackCard({ output: name, index: t.index, meta: t.metadata, gene: gene.gene }) },
              h('td', { class: 'num' }, t.index),
              h('td', null, t.metadata.biosample_name || t.metadata.name || '–'),
              h('td', null, t.metadata.ontology_curie ? curieLink(t.metadata.ontology_curie) : '–'),
              h('td', null, t.metadata.strand || '–'),
              h('td', null, [assay, t.metadata.transcription_factor || t.metadata.histone_mark].filter(Boolean).join(' · ') || '–', term ? [' ', curieLink(term.curie)] : null));
          }))))));
    }
    body.appendChild(h('div', { class: 'muted', style: { fontSize: '12px' } }, `Track metadata read from ${info.representative_sample} / ${gene.gene}.`));
  }).catch((err) => { clear(body).appendChild(errorBox(err)); });
  return card;
}

function genesCard(s) {
  const card = h('section', { class: 'card' });
  const nameCells = new Map();
  const rows = s.genes.map((g) => h('tr', { class: 'clickable', onclick: () => navigate('tracks', { gene: g.gene }) },
    h('td', null, h('b', null, g.gene), g.in_metadata ? null : h('span', { class: 'pill', style: { marginLeft: '8px' }, title: 'Window present on disk but not listed in dataset_metadata.json' }, 'extra')),
    (() => { const cell = h('td', { class: 'muted', style: { fontSize: '12px', maxWidth: '260px', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' } }); nameCells.set(g.gene, cell); return cell; })(),
    h('td', { class: 'mono' }, g.chromosome ? `${g.chromosome}:${fmtInt(g.start)}–${fmtInt(g.end)}` : '–'),
    h('td', { class: 'num' }, g.length ? fmtBp(g.length) : '–'),
    h('td', null, g.strand || '–'),
    h('td', null, h('div', { style: { display: 'flex', gap: '6px', justifyContent: 'flex-end' } },
      h('button', { class: 'btn small', title: 'Gene card: HGNC, database links, Gene Ontology', onclick: (e) => { e.stopPropagation(); openGeneCard(g.gene); } }, icon('info', 14)),
      h('button', { class: 'btn small', onclick: (e) => { e.stopPropagation(); navigate('tracks', { gene: g.gene }); } }, 'Tracks'),
      h('button', { class: 'btn small', onclick: (e) => { e.stopPropagation(); navigate('sequence', { gene: g.gene }); } }, 'Sequence')))));
  card.append(
    h('div', { class: 'card-head' }, h('h2', null, `Gene windows (${s.genes.length})`), h('span', { class: 'muted' }, 'Click a gene to open its tracks')),
    h('div', { class: 'table-wrap', style: { maxHeight: '420px' } }, h('table', { class: 'table' },
      h('thead', null, h('tr', null, h('th', null, 'Gene'), h('th', { title: 'HGNC approved name' }, 'Name (HGNC)'), h('th', null, 'Window'), h('th', { class: 'num' }, 'Size'), h('th', null, 'Strand'), h('th', null, ''))),
      h('tbody', null, rows))));
  geneNames(s.genes.map((g) => g.gene)).then((names) => {
    for (const [gene, cell] of nameCells) { const n = names[gene]; if (n) { cell.textContent = n.name; cell.title = `${n.name} · ${n.hgnc_id} · ${n.locus_type}`; } }
  }).catch(() => {});
  return card;
}

function datasetsCard() {
  const allowed = !state.status || state.status.allow_add_datasets !== false;
  const pathInput = h('input', { class: 'input mono', placeholder: '/path/to/dataset (contains dataset_metadata.json)', style: { flex: '1', minWidth: '260px' } });
  const annInput = h('input', { class: 'input mono', placeholder: 'optional: sample annotations CSV/TSV', style: { flex: '1', minWidth: '220px' } });
  const add = async () => {
    const path = pathInput.value.trim();
    if (!path) return;
    try {
      const d = await api('/api/datasets', { method: 'POST', body: { path, annotations: annInput.value.trim() || null } });
      await loadStatus();
      emit('datasets');
      await setDataset(d.id);
      emit('datasets');
      toast(`Loaded ${d.name}`);
      navigate('overview');
      window.dispatchEvent(new HashChangeEvent('hashchange'));
    } catch (err) { toast(err.message, 'error'); }
  };
  pathInput.addEventListener('keydown', (e) => { if (e.key === 'Enter') add(); });
  const list = h('table', { class: 'table' },
    h('thead', null, h('tr', null, h('th', null, 'Dataset'), h('th', null, 'Path'), h('th', { class: 'num' }, 'Samples'), h('th', { class: 'num' }, 'Genes'), h('th', null, ''))),
    h('tbody', null, state.datasets.map((d) => h('tr', { class: d.id === state.datasetId ? 'selected' : '' },
      h('td', null, h('b', null, d.name)), h('td', { class: 'mono' }, d.path), h('td', { class: 'num' }, fmtInt(d.sample_count)), h('td', { class: 'num' }, fmtInt(d.gene_count)),
      h('td', null, h('div', { style: { display: 'flex', gap: '6px', justifyContent: 'flex-end' } },
        d.id === state.datasetId ? h('span', { class: 'pill' }, 'active') : h('button', { class: 'btn small', onclick: async () => { await setDataset(d.id); emit('datasets'); window.dispatchEvent(new HashChangeEvent('hashchange')); } }, 'Switch'),
        allowed && d.id !== state.datasetId ? h('button', { class: 'btn small ghost', title: 'Remove from this list (files are kept)', onclick: async () => {
          try { await api('/api/datasets/forget', { method: 'POST', body: { id: d.id } }); await loadStatus(); emit('datasets'); window.dispatchEvent(new HashChangeEvent('hashchange')); } catch (err) { toast(err.message, 'error'); }
        } }, icon('trash', 14)) : null)))),
      // Remembered datasets whose directory is gone (moved, deleted, or a cleaned-up temp dir).
      ((state.status && state.status.missing_datasets) || []).map((path) => h('tr', { class: 'muted' },
        h('td', null, h('span', { class: 'pill', title: 'Remembered from an earlier session, but the directory no longer has dataset_metadata.json' }, 'missing')),
        h('td', { class: 'mono' }, path), h('td', { class: 'num' }, '–'), h('td', { class: 'num' }, '–'),
        h('td', null, h('div', { style: { display: 'flex', justifyContent: 'flex-end' } }, allowed ? h('button', { class: 'btn small ghost', title: 'Forget this dataset', onclick: async () => {
          try { await api('/api/datasets/forget', { method: 'POST', body: { path } }); await loadStatus(); emit('datasets'); window.dispatchEvent(new HashChangeEvent('hashchange')); } catch (err) { toast(err.message, 'error'); }
        } }, icon('trash', 14)) : null))))));
  return h('section', { class: 'card' },
    h('div', { class: 'card-head' }, h('h2', null, 'Datasets'), h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center' } },
      h('span', { class: 'muted' }, 'Any directory with the canonical layout works'),
      state.status && state.status.allow_tasks ? h('button', { class: 'btn small', title: 'From your own VCF, or one click for 1000 Genomes / gnomAD HGDP read from the public servers', onclick: () => navigate('import') }, icon('upload', 14), 'Import or quickstart…') : null)),
    h('div', { class: 'table-wrap' }, list),
    allowed ? h('div', { class: 'card-body', style: { display: 'flex', gap: '8px', flexWrap: 'wrap' } }, pathInput, annInput, h('button', { class: 'btn primary', onclick: add }, 'Open dataset')) : null);
}
