// App shell: navigation, hash router, dataset switcher, theme, cohort chip and job indicator.
import { api } from './api.js';
import { state, loadStatus, setDataset, subscribe, cohortSize, emit } from './state.js';
import { h, clear, icon, toast, progressBar, fmtInt, errorBox } from './ui.js';

const PAGES = [
  { key: 'overview', label: 'Overview', icon: 'overview', load: () => import('./pages/overview.js') },
  { key: 'samples', label: 'Samples', icon: 'samples', load: () => import('./pages/samples.js') },
  { key: 'tracks', label: 'Tracks', icon: 'tracks', load: () => import('./pages/tracks.js') },
  { key: 'sequence', label: 'Sequence', icon: 'sequence', load: () => import('./pages/sequence.js') },
  { key: 'perturb', label: 'Perturbation Lab', icon: 'perturb', load: () => import('./pages/perturb.js') },
  { sep: true },
  { key: 'experiments', label: 'Experiments', icon: 'experiments', load: () => import('./pages/experiments.js') },
  { key: 'jobs', label: 'Jobs', icon: 'jobs', load: () => import('./pages/jobs.js') },
  { key: 'alphagenome', label: 'AlphaGenome', icon: 'alphagenome', load: () => import('./pages/labs.js') },
  { key: 'system', label: 'System', icon: 'info', load: () => import('./pages/system.js') },
  { key: 'import', label: 'Import dataset', icon: 'upload', hidden: true, nav: 'jobs', load: () => import('./pages/import.js') },
];
const ALIASES = { labs: 'alphagenome' };
const NO_DATASET_PAGES = ['experiments', 'alphagenome', 'system', 'overview', 'jobs', 'import', 'perturb'];

const pageEl = document.getElementById('page');
const navEl = document.getElementById('sidenav');
let current = null; // { key, instance }
let routeToken = 0;

export function parseRoute() {
  const raw = location.hash.replace(/^#\/?/, '');
  const [path, query] = raw.split('?');
  const params = Object.fromEntries(new URLSearchParams(query || ''));
  return { page: path || 'overview', params };
}

/** Update the URL query of the current page without re-rendering (shareable links). */
export function updateRouteParams(params) {
  const { page } = parseRoute();
  const qs = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) if (v !== undefined && v !== null && v !== '') qs.set(k, v);
  const next = `#/${page}${qs.toString() ? `?${qs}` : ''}`;
  if (next !== location.hash) history.replaceState(null, '', next);
}

export function navigate(page, params = {}) {
  const qs = new URLSearchParams(params).toString();
  location.hash = `#/${page}${qs ? `?${qs}` : ''}`;
}

function renderNav(active) {
  clear(navEl);
  for (const p of PAGES) {
    if (p.sep) { navEl.appendChild(h('div', { class: 'nav-sep' })); continue; }
    if (p.hidden) continue;
    navEl.appendChild(h('a', { href: `#/${p.key}`, class: p.key === active ? 'active' : '', title: p.label, 'aria-current': p.key === active ? 'page' : null }, icon(p.icon, 18), h('span', null, p.label)));
  }
  navEl.appendChild(h('div', { class: 'nav-foot' }, state.status ? `v${state.status.version}` : ''));
}

async function route() {
  const token = ++routeToken;
  const { page, params } = parseRoute();
  const def = PAGES.find((p) => p.key === (ALIASES[page] || page)) || PAGES[0];
  renderNav(def.nav || def.key);
  document.getElementById('app').classList.remove('nav-open');
  if (current && current.key === def.key && current.instance.update) {
    current.instance.update(params);
    return;
  }
  if (current && current.instance.unmount) current.instance.unmount();
  current = null;
  clear(pageEl);
  if (!state.datasetId && !NO_DATASET_PAGES.includes(def.key)) {
    pageEl.appendChild(h('div', { class: 'page-inner' }, h('div', { class: 'empty' }, 'No dataset loaded. Add one from the Overview page.')));
    return;
  }
  try {
    const mod = await def.load();
    if (token !== routeToken) return;
    const instance = await mod.mount(pageEl, params);
    if (token !== routeToken) { if (instance && instance.unmount) instance.unmount(); return; }
    current = { key: def.key, instance: instance || {} };
    document.title = `${def.label} · Genomics Visualizer`;
  } catch (err) {
    console.error(err);
    pageEl.appendChild(h('div', { class: 'page-inner' }, errorBox(err)));
  }
}

function renderDatasetSelect() {
  const sel = document.getElementById('datasetSelect');
  clear(sel);
  if (!state.datasets.length) sel.appendChild(h('option', { value: '' }, 'No dataset'));
  for (const d of state.datasets) sel.appendChild(h('option', { value: d.id }, `${d.name} · ${fmtInt(d.sample_count)} samples`));
  sel.value = state.datasetId || '';
}

function renderCohortChip() {
  const chip = document.getElementById('cohortChip');
  if (!state.samples) { chip.hidden = true; return; }
  chip.hidden = false;
  clear(chip);
  const total = state.samples.rows.length;
  const size = cohortSize();
  chip.append(icon('samples', 14), h('span', null, h('b', null, fmtInt(size)), size === total ? ' samples' : ` of ${fmtInt(total)} in cohort`), h('span', { class: 'muted' }, '·'), h('span', null, h('b', null, state.pinned.length), ' pinned'));
}

function setupTheme() {
  const btn = document.getElementById('themeToggle');
  const isDark = () => {
    const forced = document.documentElement.getAttribute('data-theme');
    return forced ? forced === 'dark' : window.matchMedia('(prefers-color-scheme: dark)').matches;
  };
  const render = () => { clear(btn).appendChild(icon(isDark() ? 'sun' : 'moon', 17)); };
  btn.addEventListener('click', () => {
    const next = isDark() ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('gv.theme', next); } catch (e) { /* ignore */ }
    render();
    window.dispatchEvent(new Event('themechange'));
  });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { render(); window.dispatchEvent(new Event('themechange')); });
  render();
}

function setupNavToggle() {
  const btn = document.getElementById('navToggle');
  btn.appendChild(icon('menu', 17));
  const app = document.getElementById('app');
  try { if (localStorage.getItem('gv.navCollapsed') === '1') app.classList.add('nav-collapsed'); } catch (e) { /* ignore */ }
  btn.addEventListener('click', () => {
    if (window.innerWidth <= 760) { app.classList.toggle('nav-open'); return; }
    app.classList.toggle('nav-collapsed');
    try { localStorage.setItem('gv.navCollapsed', app.classList.contains('nav-collapsed') ? '1' : '0'); } catch (e) { /* ignore */ }
    window.dispatchEvent(new Event('resize'));
  });
}

// Job indicator: in-server jobs (cohort aggregates; keep running when you leave a page) and
// background jobs (imports, predictions, training; keep running when the browser or server closes).
let jobsTimer = null;
let lastTaskStates = null;
export function watchJobs() {
  if (jobsTimer) return;
  const el = document.getElementById('jobsIndicator');
  const bar = progressBar(0);
  const label = h('a', { href: '#/jobs', class: 'jobs-link' });
  let activeJobs = [];
  const cancel = h('button', { class: 'btn small ghost', title: 'Cancel running in-server jobs', onclick: async () => {
    if (!activeJobs.length || !confirm(`Cancel ${activeJobs.length} running job(s)?\n\n${activeJobs.map((j) => j.title).join('\n')}`)) return;
    await Promise.all(activeJobs.map((j) => api(`/api/jobs/${j.id}/cancel`, { method: 'POST', body: {} }).catch(() => {})));
    toast('Cancelling jobs…');
  } }, '×');
  el.append(h('span', { class: 'spinner' }), label, bar, cancel);
  const tick = async () => {
    try {
      const { jobs, tasks = [] } = await api('/api/jobs');
      activeJobs = jobs;
      notifyFinishedTasks(tasks);
      if (!jobs.length && !tasks.length) { el.hidden = true; clearInterval(jobsTimer); jobsTimer = null; return; }
      el.hidden = false;
      cancel.hidden = !jobs.length;
      const all = [...jobs, ...tasks];
      const first = all[0];
      label.textContent = all.length > 1 ? `${all.length} jobs running` : first.title;
      el.title = [...jobs.map((j) => `${j.title}: ${j.message}`), ...tasks.map((t) => `${t.title}: ${t.message} (background)`)].join('\n');
      bar.set(first.progress || 0);
    } catch (e) { /* server restarting */ }
  };
  jobsTimer = setInterval(tick, 2000);
  tick();
}

// Toast when a background job that was running in this tab finishes.
function notifyFinishedTasks(active) {
  const now = new Map(active.map((t) => [t.id, t]));
  if (lastTaskStates) {
    for (const [id, title] of lastTaskStates) {
      if (!now.has(id)) {
        api(`/api/tasks/${encodeURIComponent(id)}`, { params: { log: 0 } }).then((t) => {
          toast(`${title}: ${t.status === 'done' ? 'finished' : t.status}`, t.status === 'done' ? 'info' : 'error', 8000);
          emit('task', t);
        }).catch(() => {});
      }
    }
  }
  lastTaskStates = new Map(active.map((t) => [t.id, t.title]));
}

// New datasets (imports) and new outputs (predictions) appear without reloading the page.
async function refreshAfterTask(task) {
  try {
    await loadStatus();
    renderDatasetSelect();
    const path = task.params.dataset_path || task.params.output_dir;
    const active = state.datasets.find((d) => d.id === state.datasetId);
    if (active && path && active.path === path) {
      await setDataset(state.datasetId);
      renderCohortChip();
    }
  } catch (e) { /* ignore */ }
}

async function boot() {
  setupTheme();
  setupNavToggle();
  renderNav(parseRoute().page);
  try {
    await loadStatus();
  } catch (err) {
    pageEl.appendChild(h('div', { class: 'page-inner' }, errorBox(err)));
    return;
  }
  let preferred = null;
  try { preferred = localStorage.getItem('gv.dataset'); } catch (e) { /* ignore */ }
  const urlDataset = parseRoute().params.dataset;
  const pick = [urlDataset, preferred].find((id) => id && state.datasets.some((d) => d.id === id)) || (state.datasets[0] && state.datasets[0].id);
  if (pick) {
    try { await setDataset(pick); } catch (err) { toast(`Could not load dataset: ${err.message}`, 'error'); }
  }
  renderDatasetSelect();
  renderCohortChip();
  if ((state.status.jobs && state.status.jobs.length) || (state.status.tasks && state.status.tasks.length)) watchJobs();
  document.getElementById('datasetSelect').addEventListener('change', async (e) => {
    try {
      await setDataset(e.target.value);
      renderCohortChip();
      if (current && current.instance.unmount) current.instance.unmount();
      current = null;
      route();
    } catch (err) { toast(err.message, 'error'); }
  });
  subscribe((topic, detail) => {
    if (topic === 'cohort' || topic === 'pinned' || topic === 'dataset') renderCohortChip();
    if (topic === 'datasets') renderDatasetSelect();
    if (topic === 'job') watchJobs();
    if (topic === 'task' && detail && detail.status === 'done' && ['import', 'predict'].includes(detail.kind)) refreshAfterTask(detail);
    if (current && current.instance.onEvent) current.instance.onEvent(topic, detail);
  });
  window.addEventListener('hashchange', route);
  route();
}

boot();
