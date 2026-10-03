// Experiments: run table with generic metrics, training curves, confusion matrices, comparisons.
import { api } from '../api.js';
import { h, clear, icon, select, fmtNum, fmtDate, fmtInt, debounce, errorBox, showTooltip, hideTooltip, escapeHtml } from '../ui.js';
import { lineChart, setupCanvas, theme, seqColor, css, ColorSlots } from '../plot.js';

const PREFERRED = ['best_val_accuracy', 'best_val_loss', 'val_best_accuracy.weighted_f1_score', 'test.weighted_f1_score', 'val.weighted_accuracy'];
const LOWER_IS_BETTER = /loss|mse|mae|error/i;

export async function mount(root) {
  const page = h('div', { class: 'page-inner' });
  root.appendChild(page);
  page.appendChild(h('div', { class: 'muted' }, 'Loading runs…'));
  let data;
  try { data = await api('/api/runs'); } catch (err) { clear(page).appendChild(errorBox(err)); return {}; }
  clear(page);
  if (!data.runs.length) {
    page.append(h('div', { class: 'page-head' }, h('h1', null, 'Experiments')), h('div', { class: 'card empty' }, data.roots.length ? `No runs found under ${data.roots.join(', ')}` : 'No runs root configured. Start the visualizer with --runs-root PATH.'));
    return {};
  }
  const metrics = data.metrics;
  let columns = PREFERRED.filter((m) => metrics.includes(m));
  for (const m of metrics) { if (columns.length >= 4) break; if (!columns.includes(m)) columns.push(m); }
  let sort = { key: columns[0] || 'updated_at', dir: columns[0] && LOWER_IS_BETTER.test(columns[0]) ? 1 : -1 };
  let query = '';
  const compare = new Set();
  const colors = new ColorSlots(8);
  let selected = null;

  const search = h('input', { class: 'input', type: 'search', placeholder: 'Filter runs…', style: { width: '260px' } });
  search.addEventListener('input', debounce(() => { query = search.value.toLowerCase(); renderTable(); }, 120));
  const metricPickers = h('div', { style: { display: 'flex', gap: '6px', flexWrap: 'wrap' } });
  const tableHost = h('div', { class: 'table-wrap', style: { maxHeight: '46vh' } });
  const compareHost = h('div');
  const detailHost = h('div');
  page.append(
    h('div', { class: 'page-head' }, h('div', null, h('h1', null, 'Experiments'), h('p', null, `${fmtInt(data.runs.length)} runs · ${data.roots.join(', ')}`)), search),
    h('section', { class: 'card' }, h('div', { class: 'card-head' }, h('h2', null, 'Runs'), h('div', { style: { display: 'flex', gap: '8px', alignItems: 'center' } }, h('span', { class: 'label' }, 'Metric columns'), metricPickers)), tableHost),
    compareHost, detailHost);

  function renderPickers() {
    clear(metricPickers);
    columns.forEach((m, i) => metricPickers.appendChild(select(metrics.map((x) => ({ value: x, label: x })), m, (v) => { columns[i] = v; renderTable(); }, { style: { maxWidth: '220px' } })));
  }

  function value(run, key) {
    if (key === 'name') return run.name;
    if (key === 'updated_at') return run.updated_at || '';
    return run.metrics[key];
  }

  function renderTable() {
    const runs = data.runs.filter((r) => !query || r.id.toLowerCase().includes(query) || (r.status || '').toLowerCase().includes(query));
    runs.sort((a, b) => {
      const x = value(a, sort.key); const y = value(b, sort.key);
      if (x === undefined || x === null) return 1;
      if (y === undefined || y === null) return -1;
      return (typeof x === 'number' ? x - y : String(x).localeCompare(String(y))) * sort.dir;
    });
    const best = {};
    for (const c of columns) {
      const vals = data.runs.map((r) => r.metrics[c]).filter(Number.isFinite);
      if (vals.length) best[c] = LOWER_IS_BETTER.test(c) ? Math.min(...vals) : Math.max(...vals);
    }
    const th = (label, key, cls = '') => h('th', { class: `sortable ${cls}`, title: key, onclick: () => { sort = { key, dir: sort.key === key ? -sort.dir : (LOWER_IS_BETTER.test(key) || key === 'name' ? 1 : -1) }; renderTable(); } }, label, sort.key === key ? (sort.dir > 0 ? ' ↑' : ' ↓') : '');
    clear(tableHost).appendChild(h('table', { class: 'table' },
      h('thead', null, h('tr', null, h('th', { title: 'Compare training curves' }, 'Cmp'), th('Run', 'name'), h('th', null, 'Status'), th('Updated', 'updated_at'), columns.map((c) => th(shortMetric(c), c, 'num')), h('th', { class: 'num' }, 'Ckpts'))),
      h('tbody', null, runs.map((run) => {
        const box = h('input', { type: 'checkbox', title: 'Compare', onclick: (e) => e.stopPropagation() });
        box.checked = compare.has(run.id);
        box.disabled = !run.has_history;
        box.addEventListener('change', () => { if (box.checked) compare.add(run.id); else compare.delete(run.id); renderCompare(); });
        return h('tr', { class: `clickable${selected === run.id ? ' selected' : ''}`, onclick: () => openRun(run.id) },
          h('td', null, box),
          h('td', { title: run.path, style: { maxWidth: '520px', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' } }, run.group ? h('span', { class: 'muted' }, `${run.group}/`) : null, run.name),
          h('td', null, run.status ? h('span', { style: { display: 'inline-flex', gap: '6px', alignItems: 'center' } }, h('span', { class: `status-dot ${run.status === 'completed' ? 'good' : run.status === 'failed' ? 'bad' : 'warn'}` }), run.status) : h('span', { class: 'muted' }, '–')),
          h('td', { class: 'muted', style: { whiteSpace: 'nowrap' } }, fmtDate(run.updated_at)),
          columns.map((c) => { const v = run.metrics[c]; return h('td', { class: `num${v === best[c] ? ' metric-best' : ''}` }, fmtNum(v, 4)); }),
          h('td', { class: 'num' }, run.checkpoints.length || '–'));
      }))));
  }

  async function renderCompare() {
    clear(compareHost);
    if (compare.size < 1) return;
    const ids = [...compare].slice(0, 8);
    colors.sync(ids);
    const details = await Promise.all(ids.map((id) => api('/api/runs/detail', { params: { id } }).catch(() => null)));
    const card = h('section', { class: 'card', style: { marginTop: '16px' } });
    const keys = ['val_accuracy', 'val_loss', 'train_accuracy', 'train_loss'].filter((k) => details.some((d) => d && d.history && d.history[k]));
    let key = keys[0];
    const chartHost = h('div', { class: 'card-body' });
    const draw = () => {
      clear(chartHost);
      const canvas = h('canvas');
      chartHost.appendChild(canvas);
      const series = details.map((d, i) => d && d.history && d.history[key] ? { y: d.history[key], x: d.history.epoch, color: colors.color(ids[i]) } : null).filter(Boolean);
      lineChart(canvas, { series, width: chartHost.clientWidth - 28 || 800, height: 260, xLabel: 'epoch', yLabel: key });
      chartHost.appendChild(h('div', { class: 'legend', style: { padding: '8px 0 0', border: '0' } }, ids.map((id) => h('span', { class: 'legend-item' }, h('span', { class: 'swatch line', style: { background: colors.color(id) } }), h('span', { title: id }, shortRun(id))))));
    };
    card.append(h('div', { class: 'card-head' }, h('h2', null, `Compare ${ids.length} run${ids.length > 1 ? 's' : ''}`), h('div', { style: { display: 'flex', gap: '8px' } }, select(keys, key, (v) => { key = v; draw(); }), h('button', { class: 'btn small ghost', onclick: () => { compare.clear(); renderTable(); renderCompare(); } }, 'Clear'))), chartHost);
    compareHost.appendChild(card);
    draw();
  }

  async function openRun(id) {
    selected = id;
    renderTable();
    clear(detailHost).appendChild(h('div', { class: 'muted', style: { marginTop: '16px' } }, 'Loading run…'));
    let run;
    try { run = await api('/api/runs/detail', { params: { id } }); } catch (err) { clear(detailHost).appendChild(errorBox(err)); return; }
    clear(detailHost);
    const card = h('section', { class: 'card', style: { marginTop: '16px' } });
    const body = h('div', { class: 'card-body', style: { display: 'grid', gap: '18px' } });
    card.append(h('div', { class: 'card-head' }, h('div', null, h('h2', null, run.name), h('div', { class: 'muted mono', style: { fontSize: '11.5px', marginTop: '2px' } }, run.path)), h('button', { class: 'btn small ghost', onclick: () => { selected = null; clear(detailHost); renderTable(); } }, 'Close')), body);
    detailHost.appendChild(card);
    // headline metrics
    const headline = Object.entries(run.metrics).filter(([k]) => /accuracy|f1|loss|epoch/i.test(k.split('.').pop())).slice(0, 8);
    if (headline.length) body.appendChild(h('div', { class: 'tiles' }, headline.map(([k, v]) => h('div', { class: 'tile' }, h('div', { class: 'tile-label', title: k }, shortMetric(k)), h('div', { class: 'tile-value', style: { fontSize: '22px' } }, fmtNum(v, 4))))));
    // training curves: one chart per measure (never a dual axis)
    if (run.history) {
      const hist = run.history;
      const charts = h('div', { class: 'grid cols-2' });
      const pairs = [['Loss', ['train_loss', 'val_loss']], ['Accuracy', ['train_accuracy', 'val_accuracy']]];
      for (const [title, keys] of pairs) {
        const present = keys.filter((k) => hist[k]);
        if (!present.length) continue;
        const host = h('div', { class: 'chart-box' });
        charts.appendChild(h('div', null, h('h3', { style: { marginBottom: '6px' } }, title), host,
          h('div', { class: 'legend', style: { padding: '6px 0 0', border: '0' } }, present.map((k, i) => h('span', { class: 'legend-item' }, h('span', { class: 'swatch line', style: { background: css(`--series-${i + 1}`) } }), k.replace('_', ' '))))));
        requestAnimationFrame(() => {
          const canvas = h('canvas');
          host.appendChild(canvas);
          const series = present.map((k, i) => ({ y: hist[k], x: hist.epoch, color: css(`--series-${i + 1}`) }));
          const geo = lineChart(canvas, { series, width: host.clientWidth || 500, height: 220, xLabel: 'epoch' });
          attachCurveHover(canvas, geo, present, hist);
        });
      }
      body.appendChild(charts);
    }
    // results
    for (const res of run.results_detail || []) {
      const block = h('div', { style: { display: 'grid', gap: '10px' } }, h('h3', null, res.file));
      if (res.error) { block.appendChild(errorBox(res.error)); body.appendChild(block); continue; }
      const row = h('div', { class: 'grid cols-2' });
      if (Array.isArray(res.confusion_matrix) && res.confusion_matrix.length) row.appendChild(confusionMatrix(res.confusion_matrix, res.class_names || []));
      if (res.per_class_metrics) {
        const names = Object.keys(res.per_class_metrics);
        const cols = [...new Set(names.flatMap((n) => Object.keys(res.per_class_metrics[n])))];
        row.appendChild(h('div', { class: 'table-wrap' }, h('table', { class: 'table' },
          h('thead', null, h('tr', null, h('th', null, 'Class'), cols.map((c) => h('th', { class: 'num' }, c)))),
          h('tbody', null, names.map((n) => h('tr', null, h('td', null, n), cols.map((c) => h('td', { class: 'num' }, fmtNum(res.per_class_metrics[n][c], 4)))))))));
      }
      block.appendChild(row);
      const scalars = Object.entries(res.metrics || {});
      if (scalars.length) block.appendChild(h('div', { class: 'chips' }, scalars.map(([k, v]) => h('span', { class: 'pill', title: k }, `${k}: ${fmtNum(v, 4)}`))));
      body.appendChild(block);
    }
    // plots gallery
    const images = (run.files || []).filter((f) => f.image);
    if (images.length) body.appendChild(h('div', null, h('h3', { style: { marginBottom: '8px' } }, 'Plots'), h('div', { class: 'gallery' }, images.map((f) => {
      const url = `/api/runs/file?id=${encodeURIComponent(run.id)}&path=${encodeURIComponent(f.path)}`;
      return h('figure', null, h('a', { href: url, target: '_blank', rel: 'noopener' }, h('img', { src: url, alt: f.path, loading: 'lazy' })), h('figcaption', null, f.path));
    }))));
    if (run.config) body.appendChild(h('details', { class: 'collapsible' }, h('summary', null, 'config.yaml'), h('pre', { class: 'code-block' }, run.config)));
    if (run.manifest) body.appendChild(h('details', { class: 'collapsible' }, h('summary', null, 'manifest.json'), h('pre', { class: 'code-block' }, JSON.stringify(run.manifest, null, 2))));
    card.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  renderPickers();
  renderTable();
  return {};
}

function attachCurveHover(canvas, geo, keys, hist) {
  if (!geo) return;
  canvas.addEventListener('mousemove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const epoch = geo.xmin + ((x - geo.pad.l) / geo.pw) * (geo.xmax - geo.xmin);
    const xs = hist.epoch || [];
    let i = 0; let best = Infinity;
    xs.forEach((v, k) => { const d = Math.abs(v - epoch); if (d < best) { best = d; i = k; } });
    if (!xs.length || x < geo.pad.l || x > geo.pad.l + geo.pw) { hideTooltip(); return; }
    showTooltip(e.clientX, e.clientY, `<div class="tt-title">epoch ${xs[i]}</div>${keys.map((k, j) => `<div class="tt-row"><span class="tt-key"><span class="swatch line" style="background:${css(`--series-${j + 1}`)}"></span>${escapeHtml(k)}</span><span class="tt-val">${fmtNum(hist[k][i], 4)}</span></div>`).join('')}`);
  });
  canvas.addEventListener('mouseleave', hideTooltip);
}

function confusionMatrix(cm, names) {
  const n = cm.length;
  const cell = Math.max(36, Math.min(72, Math.floor(360 / n)));
  const labelW = 130;
  const canvas = h('canvas');
  const width = labelW + n * cell + 12;
  const height = 26 + n * cell + 30;
  requestAnimationFrame(() => {
    const t = theme();
    const ctx = setupCanvas(canvas, width, height);
    const max = Math.max(1, ...cm.flat());
    ctx.font = `11px ${t.font}`;
    for (let r = 0; r < n; r++) {
      const total = cm[r].reduce((a, b) => a + b, 0) || 1;
      for (let c = 0; c < n; c++) {
        const v = cm[r][c];
        const rgb = seqColor(t.seq, v / max);
        ctx.fillStyle = `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
        const x = labelW + c * cell; const y = 26 + r * cell;
        ctx.fillRect(x + 1, y + 1, cell - 2, cell - 2);
        const lum = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255;
        ctx.fillStyle = lum > 0.55 ? '#0b0b0b' : '#ffffff';
        ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.fillText(String(v), x + cell / 2, y + cell / 2 - 6);
        ctx.font = `10px ${t.font}`;
        ctx.fillText(`${Math.round((v / total) * 100)}%`, x + cell / 2, y + cell / 2 + 8);
        ctx.font = `11px ${t.font}`;
      }
      ctx.fillStyle = t.ink2; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      ctx.fillText(truncate(names[r] || `class ${r}`, 18), labelW - 8, 26 + r * cell + cell / 2);
    }
    ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic'; ctx.fillStyle = t.ink2;
    for (let c = 0; c < n; c++) ctx.fillText(truncate(names[c] || `${c}`, Math.floor(cell / 8)), labelW + c * cell + cell / 2, 18);
    ctx.fillStyle = t.ink3; ctx.textAlign = 'left';
    ctx.fillText('rows: true class · columns: predicted · % of row', labelW, height - 8);
  });
  return h('div', null, canvas);
}

function truncate(s, n) { s = String(s); return s.length > n ? `${s.slice(0, Math.max(1, n - 1))}…` : s; }
function shortMetric(m) { return m.replace(/^val_best_accuracy\./, 'val* ').replace(/weighted_/, 'w-').replace(/_/g, ' '); }
function shortRun(id) { return id.length > 60 ? `${id.slice(0, 28)}…${id.slice(-28)}` : id; }
