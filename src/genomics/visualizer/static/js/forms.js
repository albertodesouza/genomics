// Forms that start background jobs: AlphaGenome predictions, training and evaluation.
// Jobs run as separate processes (see the Jobs page) and survive closing the browser.
import { api } from './api.js';
import { navigate, watchJobs } from './app.js';
import { state, filterRows } from './state.js';
import { h, clear, modal, drawer, toast, field, select, segmented, checkbox, fmtInt, fmtBp, errorBox } from './ui.js';
import { fmtBytes, ontologyPicker, outputPicker, section, storedBytes, windowPicker } from './alphagenome_pickers.js';

export { COMMON_ONTOLOGIES, ontologyPicker } from './alphagenome_pickers.js';

function started(task, what) {
  toast(`${what} started in the background: it keeps running if you close this page.`);
  watchJobs();
  navigate('jobs', { task: task.id });
}

function sampleScope(datasetId) {
  const sameDataset = datasetId === state.datasetId && state.samples;
  const total = sameDataset ? state.samples.rows.length : null;
  const cohort = sameDataset ? filterRows(state.filters, '').map((r) => r[state.samples.col.sample_id]) : [];
  const pinned = sameDataset ? [...state.pinned] : [];
  let scope = 'all';
  const options = [{ value: 'all', label: total ? `All (${fmtInt(total)})` : 'All samples' }];
  if (sameDataset && cohort.length && cohort.length < total) options.push({ value: 'cohort', label: `Cohort (${fmtInt(cohort.length)})`, title: 'The Samples page filters' });
  if (sameDataset && pinned.length) options.push({ value: 'pinned', label: `Pinned (${pinned.length})` });
  const seg = segmented(options, scope, (v) => { scope = v; if (seg.onchange) seg.onchange(); });
  seg.samples = () => (scope === 'cohort' ? cohort : scope === 'pinned' ? pinned : []);
  seg.count = () => (scope === 'cohort' ? cohort.length : scope === 'pinned' ? pinned.length : total);
  return seg;
}

// ------------------------------------------------------------------------- AlphaGenome predictions
export async function openPredictForm(datasetId = state.datasetId, preset = {}) {
  let options;
  try { options = await api(`/api/d/${encodeURIComponent(datasetId)}/predict/options`); } catch (err) { toast(err.message, 'error'); return; }
  const specs = options.output_specs;
  const onDisk = Object.keys(options.existing_outputs).map((o) => o.toUpperCase()).filter((o) => specs[o]);
  let ready = false;
  const update = () => { if (ready) estimate(); };
  const outputs = outputPicker(specs, preset.outputs || (onDisk.length ? onDisk : ['RNA_SEQ']), { existing: options.existing_outputs, onChange: () => { ontologies.refresh(); update(); } });
  const ontologies = ontologyPicker(options.ontologies, preset.ontology_terms || options.ontologies.map((o) => o.curie), { outputs: () => outputs.values(), onChange: update });
  const windows = windowPicker(options, datasetId, { onChange: update, selected: preset.genes });
  const scope = sampleScope(datasetId);
  scope.onchange = update;
  // 'ref' predicts each window's reference sequence (the baseline drawn on the Tracks page).
  const haps = new Set(preset.haplotypes || ['H1', 'H2']);
  let overwrite = false;
  const estimateEl = h('div', { class: 'notice', style: { display: 'grid', gap: '4px' } });
  function estimate() {
    const counts = ontologies.trackCounts();
    outputs.setCounts(counts);
    clear(estimateEl);
    const sampleHaps = [...haps].filter((x) => x !== 'ref').length;
    const n = sampleHaps ? scope.count() : 0;
    const chosen = outputs.values();
    const hapWindows = (n ? n * windows.count() * sampleHaps : 0) + (haps.has('ref') ? windows.count() : 0) || null;
    const lines = [];
    if (hapWindows) lines.push(`${fmtInt(hapWindows)} ${haps.has('ref') ? (sampleHaps ? 'haplotype + reference ' : 'reference ') : 'haplotype '}windows (${fmtBp(options.window_size || 0)} each); windows that already have every chosen output are skipped. Roughly 1–5 s per window on the hosted API.`);
    if (!haps.size) lines.push(h('span', { style: { color: 'var(--warn, #b25e00)' } }, 'Choose H1, H2 and/or the reference genome.'));
    if (counts && hapWindows) {
      const tracks = chosen.reduce((a, o) => a + (counts[o] || 0), 0);
      const bytes = chosen.reduce((a, o) => a + storedBytes(specs[o], options.window_size || 0, counts[o] || 0), 0);
      lines.push(`${fmtInt(tracks)} tracks per window · ≈ ${fmtBytes(bytes)} per haplotype window, ≈ ${fmtBytes(bytes * hapWindows)} in total before compression.`);
      const empty = chosen.filter((o) => counts[o] === 0);
      if (empty.length) lines.push(h('span', { style: { color: 'var(--warn, #b25e00)' } }, `No ${empty.map((o) => specs[o].label).join(', ')} tracks for the chosen tissues / cell types (these outputs would be empty).`));
    }
    if (chosen.some((o) => specs[o].kind !== 'tracks')) lines.push('Contact maps and splice junctions are stored with the predictions; the track browser draws 1-D tracks only.');
    if (windows.newCount()) lines.push(`${windows.newCount()} new window(s) are first built for every sample of the dataset (bcftools; minutes per window for thousands of samples).`);
    estimateEl.append(...lines.map((l) => h('div', null, l)));
  }
  const backend = options.backend;
  const body = h('div', { style: { display: 'grid', gap: '16px' } },
    h('div', { class: backend.reasons.length ? 'error-box' : 'notice' }, h('b', null, 'AlphaGenome: '), backend.label, backend.reasons.length ? ` — ${backend.reasons.join('; ')}. ` : '. ', h('a', { href: '#/alphagenome', onclick: () => d.close() }, 'Change backend')),
    options.running.length ? h('div', { class: 'notice' }, `${options.running.length} prediction job(s) already running for this dataset; a new one waits in the queue.`) : null,
    outputs,
    section('Tissues / cell types', ontologies),
    windows,
    h('div', { class: 'grid cols-2', style: { gap: '12px' } },
      section('Samples', scope),
      section('Haplotypes', h('div', { class: 'chips' },
        ['H1', 'H2'].map((hp) => checkbox(hp, haps.has(hp), (v) => { if (v) haps.add(hp); else haps.delete(hp); update(); })),
        h('span', { title: 'Also predict the reference genome sequence of each window (references/windows/<gene>/predictions_ref): the baseline drawn on the Tracks page' },
          checkbox('Reference genome', haps.has('ref'), (v) => { if (v) haps.add('ref'); else haps.delete('ref'); update(); }))))),
    checkbox('Re-predict windows that already have these outputs (overwrite)', false, (v) => { overwrite = v; }),
    estimateEl);
  const start = h('button', { class: 'btn primary', disabled: backend.reasons.length > 0, onclick: async () => {
    if (!haps.size) { toast('Choose H1, H2 and/or the reference genome', 'error'); return; }
    const payload = {
      outputs: outputs.values(), ontology_terms: ontologies.values(),
      samples: scope.samples(), ...windows.values(), haplotypes: [...haps], overwrite,
    };
    start.disabled = true;
    try {
      const task = await api(`/api/d/${encodeURIComponent(datasetId)}/predict`, { method: 'POST', body: payload });
      d.close();
      started(task, 'AlphaGenome predictions');
    } catch (err) { toast(err.message, 'error', 10000); start.disabled = false; }
  } }, 'Start predictions');
  const d = drawer(`AlphaGenome predictions · ${(state.datasets.find((x) => x.id === datasetId) || {}).name || datasetId}`, body, { actions: [start] });
  d.body.parentElement.classList.add('wide');
  ready = true;
  estimate();
}

// ------------------------------------------------------------------------------------- training
export async function openTrainForm(datasetId = state.datasetId, preset = {}) {
  let options;
  try { options = await api(`/api/d/${encodeURIComponent(datasetId)}/train/options`); } catch (err) { toast(err.message, 'error'); return; }
  const outputNames = Object.keys(options.outputs);
  if (!outputNames.length) {
    const m = modal('Train a model', h('div', { style: { display: 'grid', gap: '10px' } },
      h('p', { style: { margin: 0 } }, 'This dataset has no AlphaGenome predictions yet: training uses them as model input.'),
      h('button', { class: 'btn primary', onclick: () => { m.close(); openPredictForm(datasetId); } }, 'Run AlphaGenome predictions first')), [h('button', { class: 'btn', onclick: () => m.close() }, 'Close')]);
    return;
  }
  const fields = options.fields.filter((f) => f.name !== 'sample_id');
  const cfg = {
    base_config: preset.base_config || options.default_base,
    base_run: preset.base_run || '',
    target_field: preset.target_field || (fields.find((f) => f.name === 'pigmentation') || fields.find((f) => f.name === 'superpopulation') || fields[0] || {}).name,
    class_map: {},
    target_name: '',
    output: outputNames.includes('rna_seq') ? 'rna_seq' : outputNames[0],
    ontology_terms: [],
    genes: [...options.genes],
    model_type: 'CNN2',
    feature_mode: 'signals_only',
    window_center_size: 32768,
    num_epochs: 100,
    batch_size: 64,
    learning_rate: 0.001,
    train_split: 0.7, val_split: 0.15, test_split: 0.15,
    seed: 13,
    run_name: '',
    evaluate_test: true,
  };
  const scope = sampleScope(datasetId);
  const classHost = h('div');
  const tracksHost = h('div');
  const preview = h('div');
  const num = (key, attrs = {}) => {
    const input = h('input', { class: 'input', type: 'number', value: cfg[key], style: { width: '110px' }, ...attrs });
    input.addEventListener('input', () => { cfg[key] = input.value === '' ? '' : Number(input.value); });
    return input;
  };
  const renderClasses = () => {
    clear(classHost);
    const f = fields.find((x) => x.name === cfg.target_field);
    if (!f) { classHost.appendChild(h('div', { class: 'muted' }, 'No categorical sample fields: add metadata to the dataset first.')); return; }
    cfg.class_map = {};
    const rows = f.counts.slice(0, 60).map((c) => {
      cfg.class_map[c.value] = c.value;
      const input = h('input', { class: 'input', value: c.value, placeholder: '(excluded)', style: { width: '100%' } });
      input.addEventListener('input', () => { cfg.class_map[c.value] = input.value.trim(); });
      return h('tr', null, h('td', null, c.value), h('td', { class: 'num' }, fmtInt(c.count)), h('td', null, input));
    });
    classHost.append(
      h('div', { class: 'table-wrap', style: { maxHeight: '200px' } }, h('table', { class: 'table' },
        h('thead', null, h('tr', null, h('th', null, 'Value'), h('th', { class: 'num' }, 'Samples'), h('th', null, 'Class (blank = exclude)'))), h('tbody', null, rows))),
      h('p', { class: 'help', style: { margin: '4px 0 0', fontSize: '12px', color: 'var(--ink-3)' } }, 'Give several values the same class name to group them (e.g. populations into pigmentation classes).'));
  };
  const renderTracks = () => {
    clear(tracksHost);
    const out = options.outputs[cfg.output];
    cfg.ontology_terms = out.ontologies.map((o) => o.curie);
    tracksHost.appendChild(h('div', { class: 'checklist', style: { maxHeight: '140px' } }, out.ontologies.map((o) => {
      const box = h('input', { type: 'checkbox' });
      box.checked = true;
      box.addEventListener('change', () => { cfg.ontology_terms = box.checked ? [...cfg.ontology_terms, o.curie] : cfg.ontology_terms.filter((x) => x !== o.curie); });
      return h('label', null, box, h('span', null, o.name || o.curie), h('span', { class: 'meta mono' }, `${o.curie} ${o.strands.join('/')}`));
    })));
  };
  const geneList = h('div', { class: 'checklist', style: { maxHeight: '140px' } }, options.genes.map((g) => {
    const box = h('input', { type: 'checkbox' });
    box.checked = true;
    box.addEventListener('change', () => { cfg.genes = box.checked ? [...cfg.genes, g] : cfg.genes.filter((x) => x !== g); });
    return h('label', null, box, h('span', null, g));
  }));
  const baseOptions = [
    ...options.base_configs.map((c) => ({ value: `config:${c.path}`, label: c.name })),
    ...options.runs.map((r) => ({ value: `run:${r.id}`, label: `run · ${r.name}` })),
  ];
  const baseSelect = select(baseOptions, cfg.base_run ? `run:${cfg.base_run}` : `config:${cfg.base_config}`, (v) => {
    if (v.startsWith('run:')) { cfg.base_run = v.slice(4); cfg.base_config = ''; } else { cfg.base_config = v.slice(7); cfg.base_run = ''; }
  }, { style: { width: '100%' } });
  const payload = () => ({ ...cfg, sample_ids: scope.samples() });
  const doPreview = async () => {
    clear(preview).appendChild(h('div', { class: 'muted' }, 'Building the config…'));
    try {
      const res = await api(`/api/d/${encodeURIComponent(datasetId)}/train/preview`, { method: 'POST', body: payload() });
      const s = res.summary;
      clear(preview).append(
        h('dl', { class: 'kv' },
          h('dt', null, 'Target'), h('dd', null, `${s.target}: ${s.classes.join(', ')}`),
          h('dt', null, 'Input'), h('dd', null, `${s.genes.length} genes × ${s.output} × ${s.ontology_terms.length} tissues → ${s.input_shape.join(' × ')} per haplotype pair`),
          h('dt', null, 'Samples'), h('dd', null, fmtInt(s.samples)),
          h('dt', null, 'Results'), h('dd', { class: 'mono' }, s.results_dir)),
        h('details', { class: 'collapsible' }, h('summary', null, 'config.yaml'), h('pre', { class: 'code-block' }, res.config)));
      return true;
    } catch (err) { clear(preview).appendChild(errorBox(err)); return false; }
  };
  const body = h('div', { style: { display: 'grid', gap: '14px' } },
    field('Start from', baseSelect),
    h('div', { class: 'grid cols-2', style: { gap: '12px' } },
      field('Predict', select(fields.map((f) => ({ value: f.name, label: `${f.label} (${f.distinct} values)` })), cfg.target_field, (v) => { cfg.target_field = v; renderClasses(); })),
      field('Target name (when grouping)', (() => { const i = h('input', { class: 'input', placeholder: 'e.g. pigmentation' }); i.addEventListener('input', () => { cfg.target_name = i.value; }); return i; })())),
    classHost,
    field('Samples', scope),
    h('div', { class: 'grid cols-2', style: { gap: '12px' } },
      field('AlphaGenome output', select(outputNames, cfg.output, (v) => { cfg.output = v; renderTracks(); })),
      field('Model', select(options.model_types, cfg.model_type, (v) => { cfg.model_type = v; }))),
    field('Tissue tracks', tracksHost),
    field(`Gene windows (${options.genes.length})`, geneList),
    h('div', { class: 'grid cols-2', style: { gap: '12px' } },
      field('Window around each gene', select([8192, 16384, 32768, 65536, 131072].map((n) => ({ value: n, label: fmtBp(n) })), cfg.window_center_size, (v) => { cfg.window_center_size = Number(v); })),
      field('Features', select([{ value: 'signals_only', label: 'Signals' }, { value: 'signals_and_masks', label: 'Signals + variant masks' }, { value: 'masks_only', label: 'Variant masks only' }], cfg.feature_mode, (v) => { cfg.feature_mode = v; }))),
    h('div', { style: { display: 'flex', gap: '12px', flexWrap: 'wrap' } },
      field('Epochs', num('num_epochs', { min: 1 })), field('Batch size', num('batch_size', { min: 1 })), field('Learning rate', num('learning_rate', { step: 'any' })), field('Seed', num('seed'))),
    h('div', { style: { display: 'flex', gap: '12px', flexWrap: 'wrap' } },
      field('Train', num('train_split', { step: 0.05 })), field('Val', num('val_split', { step: 0.05 })), field('Test', num('test_split', { step: 0.05 }))),
    field('Run name (folder under the runs root)', (() => { const i = h('input', { class: 'input', placeholder: 'default: visualizer-<date>' }); i.addEventListener('input', () => { cfg.run_name = i.value; }); return i; })()),
    checkbox('Evaluate the best checkpoint on the test split afterwards', true, (v) => { cfg.evaluate_test = v; }),
    h('div', { class: 'muted', style: { fontSize: '12px' } }, `Runs are written under ${options.runs_root} and appear on the Experiments page. Training jobs share the GPU queue: one runs at a time.`),
    preview);
  renderClasses();
  renderTracks();
  const startBtn = h('button', { class: 'btn primary', onclick: async () => {
    startBtn.disabled = true;
    try {
      const task = await api(`/api/d/${encodeURIComponent(datasetId)}/train`, { method: 'POST', body: payload() });
      d.close();
      started(task, 'Training');
    } catch (err) { toast(err.message, 'error', 10000); clear(preview).appendChild(errorBox(err)); startBtn.disabled = false; }
  } }, 'Start training');
  const d = drawer(`Train a model · ${(state.datasets.find((x) => x.id === datasetId) || {}).name || datasetId}`, body, {
    actions: [h('button', { class: 'btn', onclick: doPreview }, 'Preview config'), startBtn],
  });
}

// ----------------------------------------------------------------------------------- evaluation
export function openEvaluateForm(run) {
  const checkpoints = run.checkpoints && run.checkpoints.length ? run.checkpoints : ['best_accuracy.pt'];
  let checkpoint = checkpoints.includes('best_accuracy.pt') ? 'best_accuracy.pt' : checkpoints[0];
  let split = 'test';
  const body = h('div', { style: { display: 'grid', gap: '12px' } },
    h('div', { class: 'mono muted', style: { fontSize: '12px' } }, run.path),
    field('Checkpoint', select(checkpoints, checkpoint, (v) => { checkpoint = v; })),
    field('Split', segmented([{ value: 'test', label: 'Test' }, { value: 'val', label: 'Validation' }, { value: 'train', label: 'Train' }], split, (v) => { split = v; })),
    h('div', { class: 'muted', style: { fontSize: '12px' } }, 'Writes <split>_<checkpoint>_results.json (metrics, confusion matrix) into the run; it shows up here when done.'));
  const go = h('button', { class: 'btn primary', onclick: async () => {
    go.disabled = true;
    try {
      const task = await api('/api/runs/evaluate', { method: 'POST', body: { run: run.id, checkpoint, split } });
      m.close();
      started(task, 'Evaluation');
    } catch (err) { toast(err.message, 'error', 9000); go.disabled = false; }
  } }, 'Evaluate');
  const m = modal(`Evaluate ${run.name}`, body, [h('button', { class: 'btn', onclick: () => m.close() }, 'Cancel'), go]);
}

