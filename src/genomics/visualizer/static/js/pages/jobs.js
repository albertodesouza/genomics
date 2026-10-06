// Jobs: dataset imports, AlphaGenome predictions, training and evaluation running as separate
// processes. They keep running when this tab closes or the visualizer stops, and are listed
// again when it restarts.
import { api } from '../api.js';
import { navigate, updateRouteParams } from '../app.js';
import { state, setDataset, emit, loadStatus } from '../state.js';
import { h, clear, icon, toast, drawer, progressBar, segmented, taskStatus, fmtDate, fmtDuration, errorBox } from '../ui.js';
import { openPredictForm, openTrainForm } from '../forms.js';
import { notificationsBlocked, notificationsEnabled, notificationsSupported, setNotifications } from '../notify.js';

const KIND_LABELS = { import: 'Import', predict: 'AlphaGenome', train: 'Training', evaluate: 'Evaluation' };
const ACTIVE = ['starting', 'queued', 'running'];
const MIN_ETA_PROGRESS = 0.03; // below this the rate estimate is noise

/** Rough time left, assuming the rate so far holds; null when there is nothing to go on. */
function eta(task) {
  if (task.status !== 'running' || !(task.progress > MIN_ETA_PROGRESS) || task.progress >= 1 || !(task.elapsed > 0)) return null;
  return (task.elapsed * (1 - task.progress)) / task.progress;
}

export async function mount(root, params) {
  const page = h('div', { class: 'page-inner' });
  root.appendChild(page);
  let filter = 'all';
  let tasks = [];
  let info = { available: true, allow: true };
  let timer = null;
  let detail = null; // { id, drawer, timer }
  let unmounted = false;
  const listHost = h('div');
  const hasDataset = () => !!state.datasetId;
  page.append(
    h('div', { class: 'page-head' },
      h('div', null, h('h1', null, 'Jobs'), h('p', null, 'Imports, AlphaGenome predictions, training and evaluation run as separate processes: they keep running if you close this tab or stop the visualizer.')),
      h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap' } },
        h('button', { class: 'btn', onclick: () => navigate('import') }, icon('upload', 15), 'Import dataset'),
        h('button', { class: 'btn', disabled: !hasDataset(), title: 'For the dataset selected at the top', onclick: () => openPredictForm() }, icon('alphagenome', 15), 'AlphaGenome predictions'),
        h('button', { class: 'btn primary', disabled: !hasDataset(), title: 'For the dataset selected at the top', onclick: () => openTrainForm() }, icon('play', 14), 'Train a model'))),
    h('div', { style: { display: 'flex', gap: '12px', alignItems: 'center', marginBottom: '12px', flexWrap: 'wrap' } },
      segmented([{ value: 'all', label: 'All' }, { value: 'active', label: 'Running' }, { value: 'finished', label: 'Finished' }], filter, (v) => { filter = v; render(); }),
      notifyToggle()),
    listHost);

  async function refresh() {
    try {
      const data = await api('/api/tasks');
      tasks = data.tasks;
      info = data;
    } catch (err) { clear(listHost).appendChild(errorBox(err)); return; }
    render();
    clearTimeout(timer);
    timer = setTimeout(refresh, tasks.some((t) => ACTIVE.includes(t.status)) ? 2000 : 10000);
  }

  /** Opt-in desktop notification when a job finishes (the browser asks on the first click). */
  function notifyToggle() {
    if (!notificationsSupported()) return null;
    const box = h('input', { type: 'checkbox' });
    box.checked = notificationsEnabled();
    box.disabled = notificationsBlocked() && !box.checked;
    const label = h('label', { class: 'check', title: box.disabled ? 'Notifications are blocked for this site in the browser settings' : 'Only while this tab is in the background; a toast is shown when it is in front' }, box, 'Notify me when a job finishes');
    box.addEventListener('change', async () => {
      const on = await setNotifications(box.checked);
      box.checked = on;
      if (!on && notificationsBlocked()) { box.disabled = true; toast('The browser is blocking notifications for this site; allow them in its site settings', 'error', 9000); }
    });
    return label;
  }

  function render() {
    clear(listHost);
    if (!info.available) { listHost.appendChild(h('div', { class: 'card empty' }, 'Background jobs are not available in this server.')); return; }
    const shown = tasks.filter((t) => filter === 'all' || (filter === 'active' ? ACTIVE.includes(t.status) : !ACTIVE.includes(t.status)));
    if (!shown.length) {
      listHost.appendChild(h('div', { class: 'card empty' }, tasks.length ? 'No jobs match this filter.' : 'No jobs yet. Import a dataset, run AlphaGenome predictions or train a model.'));
      return;
    }
    listHost.append(h('section', { class: 'card' }, h('div', { class: 'table-wrap' }, h('table', { class: 'table jobs-table' },
      h('thead', null, h('tr', null, h('th', null, 'Status'), h('th', null, 'Job'), h('th', { style: { width: '32%' } }, 'Progress'), h('th', null, 'Started'), h('th', { class: 'num' }, 'Time'), h('th', null, ''))),
      h('tbody', null, shown.map(row))))),
    h('p', { class: 'muted', style: { fontSize: '12px' } }, `Job folders (spec, state, full log): ${info.dir}`));
  }

  function row(task) {
    const active = ACTIVE.includes(task.status);
    const bar = progressBar(task.progress);
    const stepText = task.steps > 1 && task.step !== null && task.step !== undefined ? `step ${task.step + 1}/${task.steps}: ${task.step_title || ''} · ` : '';
    return h('tr', { class: `clickable${detail && detail.id === task.id ? ' selected' : ''}`, onclick: () => openDetail(task.id) },
      h('td', null, taskStatus(task)),
      h('td', null, h('div', { style: { display: 'flex', gap: '8px', alignItems: 'center' } }, h('span', { class: 'pill' }, KIND_LABELS[task.kind] || task.kind), h('b', null, task.title))),
      h('td', null, active || task.status === 'done' ? bar : null,
        h('div', { class: 'muted', style: { fontSize: '11.5px', marginTop: '3px', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', maxWidth: '520px' }, title: task.message }, `${stepText}${task.message || ''}${active && task.progress ? ` · ${Math.round(task.progress * 100)}%` : ''}`)),
      h('td', { class: 'muted', style: { whiteSpace: 'nowrap' } }, fmtDate(new Date((task.created || 0) * 1000).toISOString())),
      h('td', { class: 'num', title: eta(task) === null ? '' : 'Estimated from the rate so far' },
        fmtDuration(task.elapsed),
        eta(task) === null ? null : h('div', { class: 'muted', style: { fontSize: '11px' } }, `~${fmtDuration(eta(task))} left`)),
      h('td', null, h('div', { style: { display: 'flex', gap: '6px', justifyContent: 'flex-end' } }, actions(task))));
  }

  function actions(task) {
    const out = [];
    const stop = (e) => e.stopPropagation();
    if (ACTIVE.includes(task.status)) {
      out.push(h('button', { class: 'btn small', onclick: async (e) => { stop(e); if (!confirm(`Cancel "${task.title}"?`)) return; await api(`/api/tasks/${task.id}/cancel`, { method: 'POST', body: {} }).catch((err) => toast(err.message, 'error')); refresh(); } }, 'Cancel'));
    } else {
      const result = resultLink(task);
      if (result) out.push(result);
      if (info.allow && task.status !== 'cancelled') {
        out.push(h('button', {
          class: 'btn small ghost', title: 'Run this job again with the same settings (as a new job)',
          onclick: async (e) => {
            stop(e);
            try {
              const next = await api(`/api/tasks/${task.id}/retry`, { method: 'POST', body: {} });
              toast(`Started again: ${next.title}`);
              refresh();
              openDetail(next.id);
            } catch (err) { toast(err.message, 'error', 8000); }
          },
        }, 'Retry'));
      }
      out.push(h('button', { class: 'btn small ghost', title: 'Remove from the list (deletes the job folder, not its outputs)', onclick: async (e) => { stop(e); await api(`/api/tasks/${task.id}/delete`, { method: 'POST', body: {} }).catch((err) => toast(err.message, 'error')); refresh(); } }, icon('trash', 14)));
    }
    return out;
  }

  function resultLink(task) {
    if (task.status !== 'done') return null;
    if (task.kind === 'import' || task.kind === 'predict') {
      const path = task.params.output_dir || task.params.dataset_path;
      return h('button', { class: 'btn small', onclick: async (e) => { e.stopPropagation(); await openDataset(path); } }, 'Open dataset');
    }
    if (task.kind === 'train' || task.kind === 'evaluate') return h('button', { class: 'btn small', onclick: (e) => { e.stopPropagation(); navigate('experiments'); } }, 'Experiments');
    return null;
  }

  async function openDataset(path) {
    try {
      await loadStatus();
      let d = state.datasets.find((x) => x.path === path);
      if (!d) d = await api('/api/datasets', { method: 'POST', body: { path } });
      await loadStatus();
      await setDataset(d.id);
      emit('datasets');
      navigate('overview');
    } catch (err) { toast(err.message, 'error'); }
  }

  async function openDetail(id) {
    if (detail) { clearTimeout(detail.timer); detail.drawer.close(); }
    updateRouteParams({ task: id });
    const body = h('div', { style: { display: 'grid', gap: '14px' } }, h('div', { class: 'muted' }, 'Loading…'));
    const handle = { id, timer: null };
    detail = handle;
    handle.drawer = drawer('Job', body, { onClose: () => { clearTimeout(handle.timer); if (detail === handle) detail = null; if (!unmounted) { updateRouteParams({}); render(); } } });
    handle.drawer.body.parentElement.classList.add('wide');
    const logEl = h('pre', { class: 'code-block log-view' });
    let stick = true;
    logEl.addEventListener('scroll', () => { stick = logEl.scrollTop + logEl.clientHeight >= logEl.scrollHeight - 30; });
    let built = false;
    const head = h('div');
    const load = async () => {
      let task;
      try { task = await api(`/api/tasks/${encodeURIComponent(id)}`, { params: { log: 400 } }); } catch (err) { clear(body).appendChild(errorBox(err)); return; }
      if (!built) {
        clear(body).append(head, h('div', null, h('h3', { style: { marginBottom: '6px' } }, 'Log'), logEl));
        built = true;
      }
      clear(head).append(
        h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center', justifyContent: 'space-between' } }, h('h2', null, task.title), taskStatus(task)),
        ACTIVE.includes(task.status) ? h('div', { style: { display: 'grid', gap: '4px', marginTop: '8px' } }, progressBar(task.progress), h('div', { class: 'muted', style: { fontSize: '12px' } }, task.message, eta(task) === null ? null : h('span', null, ` · ~${fmtDuration(eta(task))} left`))) : h('div', { class: task.status === 'done' ? 'notice' : 'error-box', style: { marginTop: '8px' } }, task.message || task.status),
        h('dl', { class: 'kv', style: { marginTop: '12px' } },
          h('dt', null, 'Started'), h('dd', null, fmtDate(new Date((task.created || 0) * 1000).toISOString())),
          h('dt', null, 'Duration'), h('dd', null, fmtDuration(task.elapsed)),
          ...Object.entries(task.params || {}).filter(([, v]) => v !== null && v !== '' && typeof v !== 'object').flatMap(([k, v]) => [h('dt', null, k.replace(/_/g, ' ')), h('dd', { class: typeof v === 'string' && v.startsWith('/') ? 'mono' : '' }, String(v))]),
          ...Object.entries(task.params || {}).filter(([, v]) => Array.isArray(v)).flatMap(([k, v]) => [h('dt', null, k.replace(/_/g, ' ')), h('dd', null, v.join(', '))]),
          ...Object.entries(task.result || {}).flatMap(([k, v]) => [h('dt', null, k.replace(/_/g, ' ')), h('dd', null, String(v))]),
          h('dt', null, 'Folder'), h('dd', { class: 'mono' }, task.dir)),
        h('div', { style: { display: 'flex', gap: '8px', marginTop: '10px' } },
          ACTIVE.includes(task.status) ? h('button', { class: 'btn', onclick: async () => { await api(`/api/tasks/${id}/cancel`, { method: 'POST', body: {} }); load(); refresh(); } }, 'Cancel job') : resultLink(task)),
        h('details', { class: 'collapsible', style: { marginTop: '8px' } }, h('summary', null, 'Commands'), h('pre', { class: 'code-block' }, (task.commands || []).map((c) => `# ${c.title}\n${c.command.map((x) => (/\s/.test(x) ? `'${x}'` : x)).join(' ')}`).join('\n\n'))));
      logEl.textContent = (task.log || []).join('\n') || '(no output yet)';
      if (stick) logEl.scrollTop = logEl.scrollHeight;
      if (ACTIVE.includes(task.status) && detail === handle) handle.timer = setTimeout(load, 2000);
    };
    load();
  }

  await refresh();
  if (params.task) openDetail(params.task);
  return {
    update(p) { if (p.task && (!detail || detail.id !== p.task)) openDetail(p.task); },
    onEvent(topic) { if (topic === 'task') refresh(); },
    unmount() { unmounted = true; clearTimeout(timer); if (detail) { clearTimeout(detail.timer); detail.drawer.close(); } },
  };
}
