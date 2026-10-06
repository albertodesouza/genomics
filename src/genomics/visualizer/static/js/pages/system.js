// What this machine can run (same report as `genomics doctor`): per-feature requirements with how to
// enable what is missing, hardware, and free space where the visualizer writes.
import { api } from '../api.js';
import { navigate } from '../app.js';
import { h, clear, icon, errorBox } from '../ui.js';

const DOT = { ok: 'good', warn: 'warn', missing: 'bad', info: '' };
const LABEL = { ok: 'Ready', warn: 'Partly ready', missing: 'Not available' };
const DOCS = 'https://albertodesouza.github.io/genomics/getting-started/requirements/';

function fmtGb(v) {
  if (v === null || v === undefined) return '?';
  return v >= 1024 ? `${(v / 1024).toFixed(1)} TB` : `${Math.round(v)} GB`;
}

function featureCard(f) {
  const extra = f.key === 'alphagenome'
    ? h('button', { class: 'btn', onclick: () => navigate('alphagenome') }, icon('alphagenome', 14), 'Backend')
    : null;
  return h('section', { class: 'card' },
    h('div', { class: 'card-head' },
      h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center' } },
        h('span', { class: `status-dot ${DOT[f.status]}` }), h('h2', null, f.title),
        h('span', { class: 'muted', style: { fontSize: '12px' } }, LABEL[f.status] || f.status)),
      extra),
    h('div', { class: 'card-body', style: { display: 'grid', gap: '6px' } },
      h('div', { class: 'muted', style: { fontSize: '12px' } }, `For: ${f.purpose}`),
      ...f.checks.map((c) => h('div', { style: { display: 'grid', gridTemplateColumns: '14px minmax(110px, auto) 1fr', gap: '8px', alignItems: 'baseline', fontSize: '13px' } },
        h('span', { class: `status-dot ${DOT[c.status]}` }),
        h('b', null, c.name),
        h('div', { style: { minWidth: 0, overflowWrap: 'anywhere' } }, c.detail,
          c.fix && c.status !== 'ok' ? h('div', { class: 'mono', style: { fontSize: '12px', marginTop: '2px', color: 'var(--ink-2)' } }, c.fix) : null)))));
}

export async function mount(root) {
  const page = h('div', { class: 'page-inner' });
  root.appendChild(page);
  const body = h('div', { class: 'grid', style: { maxWidth: '980px' } }, h('div', { class: 'muted' }, 'Checking…'));
  page.append(
    h('div', { class: 'page-head' }, h('div', null, h('h1', null, 'System'),
      h('p', null, 'What this machine can run, and what to install for the rest. The same report is printed by ',
        h('span', { class: 'mono' }, 'genomics doctor'), '. Sizes and hardware per feature: ',
        h('a', { href: DOCS, target: '_blank', rel: 'noopener' }, 'requirements'), '.'))),
    body);

  let data;
  try { data = await api('/api/system'); } catch (err) { clear(body).appendChild(errorBox(err)); return; }
  clear(body);
  data.features.forEach((f) => body.appendChild(featureCard(f)));

  const hw = data.hardware;
  const gpus = hw.gpus.length ? hw.gpus.map((g) => `${g.name}${g.memory ? ` (${g.memory})` : ''}`).join(', ') : 'no NVIDIA GPU detected';
  const rows = Object.entries(data.storage.locations).filter(([, v]) => v);
  body.appendChild(h('section', { class: 'card' },
    h('div', { class: 'card-head' }, h('h2', null, 'Hardware and storage')),
    h('div', { class: 'card-body', style: { display: 'grid', gap: '10px', fontSize: '13px' } },
      h('div', null, `${hw.platform} · ${hw.cpus} CPUs · ${fmtGb(hw.memory_gb)} RAM · ${gpus}`),
      h('div', { class: 'table-wrap' }, h('table', { class: 'table' },
        h('thead', null, h('tr', null, h('th', null, 'Location'), h('th', null, 'Path'), h('th', { class: 'num' }, 'Free'))),
        h('tbody', null, ...rows.map(([label, v]) => h('tr', null,
          h('td', null, label),
          h('td', { class: 'mono', style: { overflowWrap: 'anywhere' } }, v.path, v.exists ? '' : h('span', { class: 'muted' }, ' (not created yet)')),
          h('td', { class: 'num', style: { whiteSpace: 'nowrap' } }, h('span', { style: v.free_gb < 20 ? { color: 'var(--critical)', fontWeight: '600' } : {} }, fmtGb(v.free_gb)))))))),
      h('div', { class: 'muted', style: { fontSize: '12px' } },
        'A 524 kb window costs ~8–9 MB per sample with one output of ~6 tracks for both haplotypes (sequences, VCFs, predictions); each extra 1 bp-resolution track adds ~0.3 MB per haplotype. ',
        'The training-axis alignment and training caches can grow to hundreds of GB for thousands of samples.'))));
}
