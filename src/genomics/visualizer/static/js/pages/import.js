// Import a dataset from any phased VCF plus sample metadata.
//
// The VCF defines the individuals. Metadata is free-form and optional: a CSV/TSV/JSON/PLINK .fam
// file (uploaded or on the server), pasted text, or typed into the grid; its columns become the
// samples' fields (facets, group means, training targets). Nothing about 1000 Genomes fields is
// assumed. The import itself runs as a background job (see Jobs) and can be followed by AlphaGenome
// predictions. The Quickstart card fills everything in for a public phased cohort (1000 Genomes,
// gnomAD HGDP + 1KG) whose variants are streamed from the source server.
import { api } from '../api.js';
import { navigate, watchJobs } from '../app.js';
import { h, clear, icon, toast, field, select, segmented, checkbox, pathInput, readFileText, fmtInt, fmtBp, errorBox } from '../ui.js';
import { geneSearch, geneSetPicker, ontologyPicker, outputPicker, presetChecklist, section } from '../alphagenome_pickers.js';

const PAGE_SIZE = 100;
// Disk per window base and sample: both haplotypes' raw and fixed FASTA plus the window VCFs
// (measured 4.2-8.7 bytes on 1000 Genomes / gnomAD HGDP imports).
const BYTES_PER_BASE = 6;
import { fmtBytes } from '../alphagenome_pickers.js';
const VCF_FILTER = (name) => /\.(vcf|vcf\.gz|bcf)$/i.test(name);

export async function mount(root) {
  const page = h('div', { class: 'page-inner import-page' });
  root.appendChild(page);
  let defaults;
  try { defaults = await api('/api/import/defaults'); } catch (err) { page.appendChild(errorBox(err)); return {}; }

  const st = {
    vcf: '',
    vcfOverrides: {}, // chromosome -> file outside the {chrom} pattern (quickstart presets)
    reference: defaults.reference_fasta || '',
    gtf: defaults.gtf || '',
    vcfInfo: null,
    columns: [], // metadata fields (excluding the sample id)
    records: new Map(), // sample -> {field: value}
    metaInfo: null,
    onlyWithMetadata: true,
    name: '',
    outputDir: '',
    windowSize: defaults.window_size,
    genes: '',
    regions: '',
    variantFilter: 'all',
    predict: false,
    outputs: new Set(['RNA_SEQ']),
    gridQuery: '',
    gridPage: 0,
    selected: new Set(),
  };

  const sections = {
    quickstart: h('section', { class: 'card quickstart' }),
    vcf: h('section', { class: 'card step' }),
    meta: h('section', { class: 'card step' }),
    windows: h('section', { class: 'card step' }),
    output: h('section', { class: 'card step' }),
    predict: h('section', { class: 'card step' }),
    review: h('section', { class: 'card step' }),
  };
  page.append(
    h('div', { class: 'page-head' },
      h('div', null, h('h1', null, 'Import a dataset'), h('p', null, 'From a phased VCF and any sample metadata. Builds the same layout as the 1000 Genomes dataset, so every page, AlphaGenome predictions and training work on it.')),
      h('button', { class: 'btn ghost', onclick: () => navigate('jobs') }, 'Jobs →')),
    h('div', null, defaults.missing_tools.length ? h('div', { class: 'error-box', style: { marginBottom: '12px' } }, `Missing on the server: ${defaults.missing_tools.join(', ')} (conda install -c bioconda bcftools samtools).`) : null),
    h('div', { class: 'grid', style: { maxWidth: '1100px' } }, ...Object.values(sections)));

  const head = (n, title, sub) => h('div', { class: 'card-head' }, h('div', { style: { display: 'flex', gap: '10px', alignItems: 'baseline' } }, h('span', { class: 'step-num' }, n), h('h2', null, title)), sub ? h('span', { class: 'muted', style: { fontSize: '12px' } }, sub) : null);
  const vcfSamples = () => (st.vcfInfo ? st.vcfInfo.samples : []);
  const importSamples = () => {
    const all = vcfSamples();
    if (st.onlyWithMetadata && st.records.size) return all.filter((s) => st.records.has(s) && Object.keys(st.records.get(s)).length);
    return all;
  };

  // ------------------------------------------------------------------------------- quickstart
  const qs = { presets: null, error: null, busy: null, choice: {} };
  api('/api/import/quickstart').then((res) => { qs.presets = res.presets; renderQuickstart(); }).catch((err) => { qs.error = err.message; renderQuickstart(); });
  const PER_POPULATION = [1, 2, 4, 10, 25, 0];

  function renderQuickstart() {
    const s = sections.quickstart;
    clear(s);
    s.append(h('div', { class: 'card-head' },
      h('div', { style: { display: 'flex', gap: '10px', alignItems: 'baseline' } }, h('h2', null, 'Quickstart'), h('span', { class: 'muted', style: { fontSize: '12px' } }, 'public phased GRCh38 cohorts, nothing to download first')),
      h('span', { class: 'muted', style: { fontSize: '12px' } }, 'or fill in the steps below for your own VCF')));
    const body = h('div', { class: 'card-body' });
    s.appendChild(body);
    if (qs.error) { body.appendChild(errorBox(qs.error)); return; }
    if (!qs.presets) { body.appendChild(h('div', { class: 'muted' }, 'Loading…')); return; }
    body.append(
      h('div', { class: 'quickstart-grid' }, qs.presets.map(presetTile)),
      h('p', { class: 'muted', style: { fontSize: '12px', margin: '10px 0 0' } },
        'Variants are read from the source server region by region (only the chosen windows are fetched; indexes are cached in ~/.cache/genomics/remote_index). ',
        'A balanced subset takes the first unrelated samples of each population; genes default to the pigmentation set of the canonical dataset.'));
  }

  function presetTile(p) {
    const choice = qs.choice[p.id] || (qs.choice[p.id] = { per: p.default_per_population, scope: p.scopes.length ? p.scopes[0].id : null });
    const busy = qs.busy === p.id;
    const perPopulation = select(PER_POPULATION.map((n) => ({ value: n, label: n ? `${n} per population` : 'every sample' })), choice.per, (v) => { choice.per = Number(v); renderQuickstart(); });
    const estimate = () => {
      const scope = p.scopes.find((x) => x.id === choice.scope) || p;
      const samples = choice.per ? Math.min(scope.samples, scope.populations * choice.per) : scope.samples;
      return `≈ ${fmtInt(samples)} samples × ${p.default_genes.length} genes · ≈ ${fmtBytes(samples * p.default_genes.length * defaults.window_size * BYTES_PER_BASE)}`;
    };
    const go = (startNow) => async () => {
      qs.busy = p.id;
      renderQuickstart();
      try {
        const res = await api('/api/import/quickstart', { method: 'POST', body: { preset: p.id, per_population: choice.per, scope: choice.scope } });
        if (startNow) await startQuickstart(res);
        else await fillFromQuickstart(res);
      } catch (err) { toast(err.message, 'error', 10000); }
      qs.busy = null;
      renderQuickstart();
    };
    return h('div', { class: 'tile quickstart-tile' },
      h('div', { style: { display: 'flex', justifyContent: 'space-between', gap: '8px', alignItems: 'baseline' } },
        h('b', null, p.title),
        h('a', { href: p.homepage, target: '_blank', rel: 'noopener', class: 'muted', style: { fontSize: '12px', whiteSpace: 'nowrap' } }, 'source ', icon('external', 12))),
      h('p', { style: { margin: '6px 0', fontSize: '12.5px' } }, p.summary),
      h('div', { class: 'muted', style: { fontSize: '12px' } }, `${fmtInt(p.samples)} samples · ${p.populations} populations · phasing: ${p.phasing} · ${p.citation}`),
      h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap', alignItems: 'center', marginTop: '10px' } },
        perPopulation,
        p.scopes.length ? segmented(p.scopes.map((x) => ({ value: x.id, label: x.label })), choice.scope, (v) => { choice.scope = v; renderQuickstart(); }) : null),
      h('div', { class: 'tile-sub', style: { marginTop: '6px' } }, estimate()),
      h('div', { style: { display: 'flex', gap: '8px', marginTop: 'auto', paddingTop: '10px' } },
        h('button', { class: 'btn primary', disabled: !!qs.busy || missingTools(), onclick: go(true), title: 'Start the import job now with these defaults' }, icon('play', 14), busy ? 'Preparing…' : 'Import now'),
        h('button', { class: 'btn', disabled: !!qs.busy, onclick: go(false), title: 'Fill in the steps below to review or change samples, genes and AlphaGenome options' }, 'Review in the form')));
  }

  const missingTools = () => defaults.missing_tools.length > 0;

  async function startQuickstart(res) {
    const body = {
      name: res.name, vcf: res.vcf, vcf_overrides: res.vcf_overrides, reference_fasta: res.reference_fasta, gtf: defaults.gtf,
      window_size: defaults.window_size, genes: res.genes, regions: [], samples: res.samples, metadata: res.records, variant_filter: 'all', predict: null,
    };
    const task = await api('/api/import/start', { method: 'POST', body });
    toast(`${res.name}: ${fmtInt(res.samples.length)} samples from ${res.populations} populations. Import started; predictions can be run from the Overview when it finishes.`, 'info', 9000);
    watchJobs();
    navigate('jobs', { task: task.id });
  }

  async function fillFromQuickstart(res) {
    Object.assign(st, {
      vcf: res.vcf, vcfOverrides: res.vcf_overrides || {}, reference: res.reference_fasta, vcfInfo: res.inspect,
      name: res.name, outputDir: '', genes: res.genes.join(', '), regions: '', gridQuery: '', gridPage: 0,
    });
    st.selected.clear();
    metaMode = 'path';
    await parseMetadata({ path: res.metadata_path });
    toast(`${fmtInt(res.samples.length)} samples from ${res.populations} populations (${res.vcf_source === 'local' ? 'local copy of the VCFs' : 'streamed from the source'}). Review the steps below, then Start import.`, 'info', 9000);
    sections.vcf.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  // ------------------------------------------------------------------------------- 1. variants
  function renderVcf() {
    const s = sections.vcf;
    clear(s);
    const vcf = pathInput({ value: st.vcf, placeholder: '/data/cohort.vcf.gz   or   /data/chr{chrom}.vcf.gz   or   https://…/chr{chrom}.vcf.gz', filter: VCF_FILTER, onChange: (v) => { if (v !== st.vcf) st.vcfOverrides = {}; st.vcf = v; } });
    const ref = pathInput({ value: st.reference, placeholder: '/refs/GRCh38.fa (with .fai), or its URL', filter: (n) => /\.(fa|fasta|fna)(\.gz)?$/i.test(n), onChange: (v) => { st.reference = v; } });
    vcf.input.addEventListener('keydown', (e) => { if (e.key === 'Enter') inspect(); });
    const btn = h('button', { class: 'btn primary', onclick: () => inspect() }, 'Read samples');
    const info = h('div');
    if (st.vcfInfo) {
      const v = st.vcfInfo;
      info.append(h('div', { class: 'notice', style: { display: 'grid', gap: '4px' } },
        h('div', null, h('b', null, `${fmtInt(v.sample_count)} samples`), ` · ${v.file_count} file${v.file_count > 1 ? 's' : ''} · contigs: ${v.contigs.slice(0, 8).join(', ')}${v.contigs.length > 8 ? '…' : ''}`),
        h('div', null, `Phasing in the first records: ${fmtInt(v.phasing.phased)} phased / ${fmtInt(v.phasing.unphased)} unphased genotypes`),
        v.warnings.map((w) => h('div', { style: { color: 'var(--critical)' } }, w))));
    }
    s.append(head(1, 'Variants', 'VCF (bgzipped, phased) and the reference it was called against'),
      h('div', { class: 'card-body', style: { display: 'grid', gap: '12px' } },
        field('VCF file or per-chromosome pattern (use {chrom})', h('div', { style: { display: 'flex', gap: '8px' } }, h('div', { style: { flex: '1', display: 'grid' } }, vcf), btn)),
        field('Reference FASTA', ref),
        info));
  }

  async function inspect() {
    if (!st.vcf) { toast('Enter the VCF path', 'error'); return; }
    sections.vcf.classList.add('stale');
    try {
      st.vcfInfo = await api('/api/import/inspect', { method: 'POST', body: { vcf: st.vcf } });
      if (!st.name) st.name = (st.vcf.split('/').pop() || 'dataset').replace(/\.(vcf|bcf)(\.gz)?$/i, '').replace(/[._-]*\{chrom\}/g, '').replace(/[._-]+$/, '');
      if (st.metaInfo) await parseMetadata(st.lastMeta);
    } catch (err) { toast(err.message, 'error', 9000); }
    sections.vcf.classList.remove('stale');
    renderAll();
  }

  // ------------------------------------------------------------------------------- 2. metadata
  let metaMode = 'file';
  function renderMeta() {
    const s = sections.meta;
    clear(s);
    s.append(head(2, 'Sample metadata', 'optional · any columns'));
    const body = h('div', { class: 'card-body', style: { display: 'grid', gap: '12px' } });
    s.append(body);
    if (!st.vcfInfo) { body.appendChild(h('div', { class: 'muted' }, 'Read the VCF samples first.')); return; }
    const modeSeg = segmented([
      { value: 'file', label: 'Upload file' }, { value: 'path', label: 'File on server' },
      { value: 'paste', label: 'Paste' }, { value: 'manual', label: 'Enter manually' },
    ], metaMode, (v) => { metaMode = v; renderMeta(); });
    body.appendChild(modeSeg);
    if (metaMode === 'file') {
      const input = h('input', { type: 'file', accept: '.csv,.tsv,.txt,.json,.fam,.ped,.panel', onchange: async () => {
        const file = input.files[0];
        if (!file) return;
        await parseMetadata({ text: await readFileText(file), filename: file.name });
      } });
      body.appendChild(h('label', { class: 'dropzone' }, icon('upload', 18), h('span', null, 'Choose a CSV, TSV, JSON or PLINK .fam file'), input));
    } else if (metaMode === 'path') {
      const p = pathInput({ placeholder: '/data/samples.tsv' });
      body.appendChild(h('div', { style: { display: 'flex', gap: '8px' } }, h('div', { style: { flex: '1', display: 'grid' } }, p), h('button', { class: 'btn', onclick: () => parseMetadata({ path: p.input.value.trim() }) }, 'Load')));
    } else if (metaMode === 'paste') {
      const area = h('textarea', { class: 'input mono', rows: 6, placeholder: 'sample\tgroup\tage\nS001\tcase\t54\nS002\tcontrol\t61\n… (tab, comma or semicolon separated; or JSON)' });
      area.value = st.pasted || '';
      area.addEventListener('input', () => { st.pasted = area.value; });
      body.appendChild(area);
      body.appendChild(h('div', null, h('button', { class: 'btn', onclick: () => parseMetadata({ text: area.value, filename: 'pasted' }) }, 'Use this table')));
    } else {
      body.appendChild(h('p', { class: 'muted', style: { margin: 0 } }, 'Every VCF sample is a row below. Add columns and type values, paste a column from a spreadsheet, or derive a column from the sample ids.'));
    }
    if (st.metaInfo && metaMode !== 'manual') body.appendChild(mappingControls());
    body.appendChild(grid());
  }

  async function parseMetadata(source, roles = {}) {
    if (!source) return;
    try {
      const res = await api('/api/import/metadata', { method: 'POST', body: { ...source, ...roles, samples: vcfSamples() } });
      st.lastMeta = source;
      st.metaInfo = res;
      st.records = new Map(Object.entries(res.records));
      const cols = [];
      for (const rec of st.records.values()) for (const k of Object.keys(rec)) if (!cols.includes(k)) cols.push(k);
      st.columns = cols;
      st.onlyWithMetadata = res.matched > 0 && res.matched < vcfSamples().length;
      st.gridPage = 0;
      if (!res.matched) toast(`No metadata row matched a VCF sample using column "${res.id_column}". Pick the sample-id column below.`, 'error', 9000);
    } catch (err) { toast(err.message, 'error', 9000); }
    renderAll();
  }

  function mappingControls() {
    const m = st.metaInfo;
    const opts = m.columns.map((c) => ({ value: c, label: c }));
    const none = [{ value: '', label: '— none —' }, ...opts];
    const reparse = (patch) => parseMetadata(st.lastMeta, { id_column: m.id_column, family_column: m.family_column || '', sex_column: m.sex_column || '', ...patch });
    return h('div', { style: { display: 'grid', gap: '10px' } },
      h('div', { style: { display: 'flex', gap: '12px', flexWrap: 'wrap' } },
        field('Sample id column', select(opts, m.id_column, (v) => reparse({ id_column: v }))),
        field('Family column (keeps relatives in one split)', select(none, m.family_column || '', (v) => reparse({ family_column: v }))),
        field('Sex column (optional)', select(none, m.sex_column || '', (v) => reparse({ sex_column: v })))),
      h('div', { class: m.matched ? 'notice' : 'error-box' },
        h('b', null, `${fmtInt(m.matched)} of ${fmtInt(m.row_count)} metadata rows match VCF samples`),
        m.samples_without_metadata_count ? ` · ${fmtInt(m.samples_without_metadata_count)} VCF samples have no metadata` : ' · every VCF sample has metadata',
        m.unmatched_metadata.length ? h('div', { class: 'mono', style: { marginTop: '4px', fontSize: '11.5px' } }, `Not in the VCF: ${m.unmatched_metadata.slice(0, 12).join(', ')}${m.unmatched_metadata.length > 12 ? '…' : ''}`) : null));
  }

  // -- editable grid ------------------------------------------------------------------------
  function grid() {
    const wrap = h('div', { style: { display: 'grid', gap: '8px' } });
    const all = vcfSamples();
    const withMeta = all.filter((s) => st.records.has(s));
    const query = st.gridQuery.toLowerCase();
    const rows = all.filter((s) => {
      if (st.onlyWithMetadata && st.records.size && !st.records.has(s)) return false;
      if (!query) return true;
      if (s.toLowerCase().includes(query)) return true;
      const rec = st.records.get(s) || {};
      return Object.values(rec).some((v) => String(v).toLowerCase().includes(query));
    });
    const pages = Math.max(1, Math.ceil(rows.length / PAGE_SIZE));
    st.gridPage = Math.min(st.gridPage, pages - 1);
    const shown = rows.slice(st.gridPage * PAGE_SIZE, (st.gridPage + 1) * PAGE_SIZE);
    const search = h('input', { class: 'input', type: 'search', placeholder: 'Filter rows…', value: st.gridQuery, style: { width: '200px' } });
    search.addEventListener('change', () => { st.gridQuery = search.value; st.gridPage = 0; renderMeta(); });
    const toolbar = h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap', alignItems: 'center' } },
      search,
      st.records.size ? checkbox(`Only samples with metadata (${fmtInt(withMeta.length)})`, st.onlyWithMetadata, (v) => { st.onlyWithMetadata = v; st.gridPage = 0; renderAll(); }) : null,
      h('span', { style: { flex: '1' } }),
      h('button', { class: 'btn small', onclick: addColumn }, '+ Column'),
      h('button', { class: 'btn small', onclick: deriveColumn, title: 'Fill a column from part of the sample id (regular expression)' }, 'Derive from id'),
      h('button', { class: 'btn small', disabled: !st.columns.length, onclick: () => setValues(rows) }, st.selected.size ? `Set value (${st.selected.size} selected)` : `Set value (${fmtInt(rows.length)} shown)`));
    const headRow = h('tr', null,
      h('th', { style: { width: '28px' } }, (() => { const b = h('input', { type: 'checkbox', title: 'Select all matching rows' }); b.checked = rows.length > 0 && rows.every((r) => st.selected.has(r)); b.addEventListener('change', () => { for (const r of rows) { if (b.checked) st.selected.add(r); else st.selected.delete(r); } renderMeta(); }); return b; })()),
      h('th', null, 'Sample'),
      st.columns.map((c) => h('th', null, h('span', { class: 'col-head' }, c, h('button', { class: 'col-del', title: `Remove column ${c}`, onclick: () => { st.columns = st.columns.filter((x) => x !== c); for (const rec of st.records.values()) delete rec[c]; renderAll(); } }, '×')))));
    const body = h('tbody', null, shown.map((sample) => {
      const rec = st.records.get(sample) || {};
      const box = h('input', { type: 'checkbox' });
      box.checked = st.selected.has(sample);
      box.addEventListener('change', () => { if (box.checked) st.selected.add(sample); else st.selected.delete(sample); });
      return h('tr', null, h('td', null, box), h('td', { class: 'mono' }, sample), st.columns.map((c) => {
        const input = h('input', { class: 'cell-input', value: rec[c] ?? '' });
        input.addEventListener('change', () => setCell(sample, c, input.value));
        input.addEventListener('paste', (e) => pasteColumn(e, sample, c, rows));
        return h('td', null, input);
      }));
    }));
    wrap.append(toolbar,
      h('div', { class: 'table-wrap meta-grid' }, h('table', { class: 'table' }, h('thead', null, headRow), body)),
      h('div', { style: { display: 'flex', gap: '8px', alignItems: 'center', fontSize: '12px' }, class: 'muted' },
        `${fmtInt(rows.length)} rows`, pages > 1 ? h('span', null,
          h('button', { class: 'btn small ghost', disabled: st.gridPage === 0, onclick: () => { st.gridPage -= 1; renderMeta(); } }, '‹'),
          ` page ${st.gridPage + 1} / ${pages} `,
          h('button', { class: 'btn small ghost', disabled: st.gridPage >= pages - 1, onclick: () => { st.gridPage += 1; renderMeta(); } }, '›')) : null,
        h('span', null, 'Tip: paste a column copied from a spreadsheet into a cell to fill the rows below it.')),
      fieldSummary());
    return wrap;
  }

  function coerce(value) {
    const text = String(value ?? '').trim();
    if (text === '') return null;
    if (/^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(text)) { const n = Number(text); if (Number.isFinite(n)) return n; }
    return text;
  }

  function setCell(sample, column, value) {
    const rec = st.records.get(sample) || {};
    const v = coerce(value);
    if (v === null) delete rec[column]; else rec[column] = v;
    st.records.set(sample, rec);
    renderSummaries();
    renderReview();
  }

  function pasteColumn(e, sample, column, rows) {
    const text = (e.clipboardData || window.clipboardData).getData('text');
    if (!text.includes('\n')) return;
    e.preventDefault();
    const values = text.replace(/\r/g, '').split('\n');
    if (values[values.length - 1] === '') values.pop();
    const start = rows.indexOf(sample);
    values.forEach((v, i) => { const s = rows[start + i]; if (s) setCell(s, column, v.split('\t')[0]); });
    toast(`Filled ${Math.min(values.length, rows.length - start)} rows of ${column}`);
    renderAll();
  }

  function addColumn() {
    const name = (prompt('New column name (e.g. phenotype, group, batch)') || '').trim().replace(/\s+/g, '_');
    if (!name) return;
    if (name === 'sample_id' || st.columns.includes(name)) { toast('That column already exists', 'error'); return; }
    st.columns.push(name);
    renderAll();
  }

  function deriveColumn() {
    const name = (prompt('Column to fill from the sample ids') || '').trim().replace(/\s+/g, '_');
    if (!name) return;
    const pattern = prompt('Regular expression; the first (group) becomes the value.\nExamples: ^([A-Z]+)  ·  _(case|control)_  ·  ^(..)', '^([A-Za-z]+)');
    if (!pattern) return;
    let re;
    try { re = new RegExp(pattern); } catch (err) { toast(`Invalid expression: ${err.message}`, 'error'); return; }
    if (!st.columns.includes(name)) st.columns.push(name);
    let n = 0;
    for (const sample of vcfSamples()) {
      const m = re.exec(sample);
      if (m) { setCell(sample, name, m[1] ?? m[0]); n += 1; }
    }
    toast(`${name}: ${fmtInt(n)} samples matched`);
    renderAll();
  }

  function setValues(rows) {
    const targets = st.selected.size ? [...st.selected] : rows;
    const column = st.columns.length === 1 ? st.columns[0] : prompt(`Column (${st.columns.join(', ')})`, st.columns[0]);
    if (!column || !st.columns.includes(column)) return;
    const value = prompt(`Value of ${column} for ${targets.length} samples (empty clears it)`);
    if (value === null) return;
    for (const s of targets) setCell(s, column, value);
    st.selected.clear();
    renderAll();
  }

  function fieldStats() {
    const samples = importSamples();
    return st.columns.map((c) => {
      const counts = new Map();
      let missing = 0;
      for (const s of samples) {
        const v = (st.records.get(s) || {})[c];
        if (v === undefined || v === null || v === '') { missing += 1; continue; }
        counts.set(String(v), (counts.get(String(v)) || 0) + 1);
      }
      return { name: c, distinct: counts.size, missing, counts: [...counts.entries()].sort((a, b) => b[1] - a[1]) };
    });
  }

  const summaryHost = h('div');
  function fieldSummary() { renderSummaries(); return summaryHost; }
  function renderSummaries() {
    clear(summaryHost);
    const stats = fieldStats();
    if (!stats.length) return;
    summaryHost.appendChild(h('div', { class: 'field-summary' }, stats.map((f) => h('div', { class: 'tile', style: { padding: '10px 12px' } },
      h('div', { class: 'tile-label' }, f.name),
      h('div', { style: { fontSize: '12px', marginTop: '4px' } }, f.distinct > 1 && f.distinct <= 50 ? `${f.distinct} groups: ${f.counts.slice(0, 4).map(([v, n]) => `${v} (${n})`).join(', ')}${f.counts.length > 4 ? '…' : ''}` : `${f.distinct} distinct values`),
      f.missing ? h('div', { class: 'tile-sub' }, `${fmtInt(f.missing)} missing`) : null))));
  }

  // ------------------------------------------------------------------------------- 3. windows
  function renderWindows() {
    const s = sections.windows;
    clear(s);
    const genes = h('textarea', { class: 'input mono', rows: 3, placeholder: 'OCA2, HERC2, SLC24A5  (HGNC symbols or Ensembl ids; one per line or comma separated)' });
    genes.value = st.genes;
    genes.addEventListener('input', () => { st.genes = genes.value; renderReview(); });
    const regions = h('textarea', { class: 'input mono', rows: 2, placeholder: 'rs12913832=chr15:28120472\nmy_enhancer=chr2:1200000-1203000' });
    regions.value = st.regions;
    regions.addEventListener('input', () => { st.regions = regions.value; renderReview(); });
    const presetLine = (p) => `${p.name}=${p.chrom}:${p.position}`;
    const regionLines = () => st.regions.split('\n').map((x) => x.trim()).filter(Boolean);
    const presets = presetChecklist(defaults.region_presets || [], (p) => regionLines().includes(presetLine(p)), (p, on) => {
      const lines = regionLines().filter((l) => l !== presetLine(p));
      if (on) lines.push(presetLine(p));
      st.regions = lines.join('\n');
      regions.value = st.regions;
      renderReview();
    });
    const addGenes = (list) => {
      const names = parseList(st.genes);
      for (const g of list) if (!names.includes(g.name)) names.push(g.name);
      st.genes = names.join(', ');
      genes.value = st.genes;
      renderReview();
    };
    const search = geneSearch(null, (g) => addGenes([g]), { placeholder: 'Search genes (HGNC symbol, alias, name or Ensembl id) to add them' });
    const gtf = pathInput({ value: st.gtf, placeholder: 'gtf_cache.feather (GENCODE gene table)', onChange: (v) => { st.gtf = v; } });
    s.append(head(3, 'Windows', 'one AlphaGenome window centred on each gene or region'),
      h('div', { class: 'card-body', style: { display: 'grid', gap: '12px' } },
        h('div', { class: 'grid cols-2', style: { gap: '12px' } }, section('Genes', h('div', { style: { display: 'grid', gap: '6px' } }, search, genes)), field('Regions (NAME=chr:start[-end])', regions)),
        h('details', { class: 'collapsible' }, h('summary', null, 'Genes of a Gene Ontology term or HGNC gene group'), geneSetPicker(null, addGenes)),
        h('details', { class: 'collapsible' }, h('summary', null, `Regulatory elements (not genes) · ${(defaults.region_presets || []).length}`),
          h('p', { class: 'muted', style: { fontSize: '12px', margin: '6px 0' } }, 'Windows centred on a causal or lead variant of well-studied enhancers; ticking one adds it to the regions.'),
          presets),
        h('div', { class: 'grid cols-2', style: { gap: '12px' } },
          field('Window size', select(defaults.window_sizes.map((n) => ({ value: n, label: `${fmtBp(n)}${n === 524288 ? ' (like 1000 Genomes)' : ''}` })), st.windowSize, (v) => { st.windowSize = Number(v); renderReview(); })),
          field('Gene table (for gene symbols)', gtf)),
        field('Variants applied to the haplotypes', segmented([{ value: 'all', label: 'SNVs, indels and SVs' }, { value: 'snps', label: 'SNVs only' }], st.variantFilter, (v) => { st.variantFilter = v; }))));
  }

  // ------------------------------------------------------------------------------- 4. output
  function renderOutput() {
    const s = sections.output;
    clear(s);
    const name = h('input', { class: 'input', value: st.name, placeholder: 'my cohort' });
    const slug = () => (st.name || 'dataset').toLowerCase().replace(/[^a-z0-9_.-]+/g, '-').replace(/^-+|-+$/g, '');
    const out = pathInput({ value: st.outputDir, placeholder: `${defaults.output_root}/${slug()}`, onChange: (v) => { st.outputDir = v; renderReview(); } });
    name.addEventListener('input', () => { st.name = name.value; out.input.placeholder = `${defaults.output_root}/${slug()}`; renderReview(); });
    s.append(head(4, 'Dataset'), h('div', { class: 'card-body grid cols-2', style: { gap: '12px' } }, field('Name', name), field('Output directory (new or empty)', out)));
  }

  // ------------------------------------------------------------------------------- 5. AlphaGenome
  let picker = null;
  let outputPick = null;
  function renderPredict() {
    const s = sections.predict;
    clear(s);
    const backend = defaults.backend;
    outputPick = outputPick || outputPicker(defaults.output_specs, [...st.outputs], { onChange: () => { st.outputs = new Set(outputPick.values()); if (picker) { picker.refresh(); outputPick.setCounts(picker.trackCounts()); } renderReview(); } });
    picker = picker || ontologyPicker([], ['CL:1000458'], { outputs: () => outputPick.values(), onChange: () => { if (picker && outputPick) outputPick.setCounts(picker.trackCounts()); } });
    const body = h('div', { class: 'card-body', style: { display: 'grid', gap: '12px' } },
      checkbox('Run AlphaGenome predictions when the import finishes', st.predict, (v) => { st.predict = v; renderPredict(); renderReview(); }));
    if (st.predict) {
      body.append(
        h('div', { class: backend.reasons.length ? 'error-box' : 'notice' }, h('b', null, 'Backend: '), backend.label, backend.reasons.length ? ` — ${backend.reasons.join('; ')}` : '', ' · ', h('a', { href: '#/alphagenome' }, 'change')),
        outputPick,
        section('Tissues / cell types', picker));
    } else {
      body.appendChild(h('div', { class: 'muted', style: { fontSize: '12px' } }, 'You can also run them later from the dataset\'s Overview or the Jobs page.'));
    }
    s.append(head(5, 'AlphaGenome', 'optional'), body);
  }

  // ------------------------------------------------------------------------------- 6. review
  function parseList(text) { return text.split(/[\s,;]+/).map((x) => x.trim()).filter(Boolean); }
  function renderReview() {
    const s = sections.review;
    clear(s);
    const samples = importSamples();
    const genes = parseList(st.genes);
    const regions = st.regions.split('\n').map((x) => x.trim()).filter(Boolean);
    const windows = genes.length + regions.length;
    const bytes = samples.length * windows * st.windowSize * BYTES_PER_BASE;
    const ready = st.vcfInfo && samples.length && windows && st.reference && st.name.trim();
    const start = h('button', { class: 'btn primary', disabled: !ready, onclick: async () => {
      start.disabled = true;
      const metadata = {};
      for (const sample of samples) metadata[sample] = st.records.get(sample) || {};
      const body = {
        name: st.name.trim(), output_dir: st.outputDir, vcf: st.vcf, vcf_overrides: st.vcfOverrides, reference_fasta: st.reference, gtf: st.gtf,
        window_size: st.windowSize, genes, regions, samples, metadata, variant_filter: st.variantFilter,
        predict: st.predict ? { enabled: true, outputs: [...st.outputs], ontology_terms: picker.values() } : null,
      };
      try {
        const task = await api('/api/import/start', { method: 'POST', body });
        toast('Import started in the background: it keeps running if you close this page.');
        watchJobs();
        navigate('jobs', { task: task.id });
      } catch (err) { toast(err.message, 'error', 10000); start.disabled = false; }
    } }, icon('play', 14), 'Start import');
    s.append(head(6, 'Review'), h('div', { class: 'card-body', style: { display: 'grid', gap: '10px' } },
      h('dl', { class: 'kv' },
        h('dt', null, 'Individuals'), h('dd', null, st.vcfInfo ? `${fmtInt(samples.length)} of ${fmtInt(vcfSamples().length)} VCF samples` : '–'),
        h('dt', null, 'Sample fields'), h('dd', null, st.columns.length ? st.columns.join(', ') : 'none (you can add metadata later with --annotations)'),
        h('dt', null, 'Windows'), h('dd', null, windows ? `${windows} × ${fmtBp(st.windowSize)}` : '–'),
        h('dt', null, 'Disk'), h('dd', null, windows && samples.length ? `≈ ${fmtBytes(bytes)} of sequences${st.predict ? ' plus predictions' : ''}` : '–'),
        h('dt', null, 'AlphaGenome'), h('dd', null, st.predict ? `${[...st.outputs].join(', ').toLowerCase()} afterwards (${fmtInt(samples.length * windows * 2)} haplotype windows)` : 'not now')),
      h('div', null, start)));
  }

  function renderAll() {
    renderQuickstart();
    renderVcf();
    renderMeta();
    renderWindows();
    renderOutput();
    renderPredict();
    renderReview();
  }

  renderAll();
  return {};
}
