// Saved sessions: name what a page is showing now, and come back to it later.
//
// A session is the state the pages already keep per dataset in the browser — the locus and track
// settings, the cohort filters, the pinned individuals and the sample search — stored on the server
// so it survives a cleared cache and can be opened on another machine. The payload is written and
// read here only; the server treats it as an opaque blob.
import { api } from './api.js';
import { ds, state, setFilters, setLocus, setPinned } from './state.js';
import { remount } from './app.js';
import { h, clear, icon, modal, toast, field, fmtDate } from './ui.js';

/** The part of the app state a session restores. */
export function currentSession() {
  return { locus: state.locus || {}, filters: state.filters || {}, pinned: state.pinned || [], search: state.search || '' };
}

export async function loadSessions() {
  try {
    return (await api(`${ds()}/sessions`)).sessions || [];
  } catch (err) {
    toast(err.message, 'error');
    return [];
  }
}

/** Apply a saved session to the shared state and re-render the current page. */
export async function openSession(name, { navigateTo = true } = {}) {
  let entry;
  try {
    entry = await api(`${ds()}/sessions`, { params: { name } });
  } catch (err) { toast(err.message, 'error', 8000); return null; }
  const saved = entry.state || {};
  // Pinned ids and filter values can refer to samples or fields the dataset no longer has; the
  // state setters drop what is gone, so a stale session opens as far as it still makes sense.
  if (saved.locus) setLocus(saved.locus);
  if (saved.filters) setFilters(saved.filters);
  if (saved.pinned) setPinned(saved.pinned);
  state.search = saved.search || '';
  const gene = (saved.locus || {}).gene;
  if (navigateTo && gene) {
    const params = new URLSearchParams({ gene });
    location.hash = `#/tracks?${params.toString()}`;
  }
  remount();
  toast(`Opened “${name}”`);
  return entry;
}

/** Save the current state under a name (asks before replacing one). */
export function saveSessionForm({ onSaved } = {}) {
  const name = h('input', { class: 'input', placeholder: 'e.g. TYRP1 promoter, AFR vs EUR', style: { width: '100%' } });
  const description = h('input', { class: 'input', placeholder: 'optional: what this view is for', style: { width: '100%' } });
  const locus = state.locus || {};
  const summary = [
    locus.gene ? `${locus.gene}${locus.start !== undefined ? ` ${Math.round(locus.start)}–${Math.round(locus.end)}` : ''}` : 'no locus yet',
    `${(locus.tracks || []).length} tracks`,
    locus.mode ? `${locus.mode} view` : null,
    `${(state.pinned || []).length} pinned`,
    Object.keys(state.filters || {}).length ? `cohort: ${Object.keys(state.filters).join(', ')}` : 'whole cohort',
  ].filter(Boolean).join(' · ');
  const body = h('div', { style: { display: 'grid', gap: '10px', minWidth: 'min(460px, 80vw)' } },
    field('Name', name), field('Note', description),
    h('p', { class: 'help', style: { margin: 0 } }, `Saves: ${summary}`));
  const go = async (replace) => {
    try {
      await api(`${ds()}/sessions`, { method: 'POST', body: { name: name.value, description: description.value, state: currentSession(), replace } });
    } catch (err) {
      if (!replace && /already exists/.test(err.message)) {
        if (confirm(`“${name.value.trim()}” already exists. Replace it?`)) return go(true);
        return;
      }
      toast(err.message, 'error', 8000);
      return;
    }
    m.close();
    toast(`Saved “${name.value.trim()}”`);
    if (onSaved) onSaved();
  };
  const m = modal('Save this view', body, [
    h('button', { class: 'btn', onclick: () => m.close() }, 'Cancel'),
    h('button', { class: 'btn primary', onclick: () => go(false) }, 'Save'),
  ]);
  name.focus();
  return m;
}

export async function deleteSession(name) {
  if (!confirm(`Delete the saved view “${name}”?`)) return false;
  try {
    await api(`${ds()}/sessions/delete`, { method: 'POST', body: { name } });
    toast(`Deleted “${name}”`);
    return true;
  } catch (err) { toast(err.message, 'error'); return false; }
}

export async function renameSession(name) {
  const next = prompt('New name', name);
  if (next === null || next.trim() === name) return false;
  try {
    await api(`${ds()}/sessions/rename`, { method: 'POST', body: { name, new_name: next } });
    return true;
  } catch (err) { toast(err.message, 'error'); return false; }
}

/** One row per saved session, for the Overview card. */
export function sessionList(sessions, { onChanged } = {}) {
  if (!sessions.length) {
    return h('p', { class: 'muted', style: { margin: 0, fontSize: '12.5px' } },
      'None yet. Set up a view on the Tracks page (locus, tracks, cohort, pinned individuals) and use ', h('b', null, 'Save view…'), ' there.');
  }
  const refresh = () => { if (onChanged) onChanged(); };
  return h('div', { class: 'session-list' }, sessions.map((s) => h('div', { class: 'session-item' },
    h('button', {
      class: 'session-open', title: 'Open this view',
      onclick: () => openSession(s.name),
    }, h('span', { class: 'session-name' }, s.name),
      h('span', { class: 'session-meta' }, s.description ? `${s.description} · ` : '', fmtDate(new Date(s.saved_at * 1000).toISOString()))),
    h('div', { style: { display: 'flex', gap: '2px' } },
      h('button', { class: 'btn small ghost', title: 'Rename', onclick: async () => { if (await renameSession(s.name)) refresh(); } }, 'Rename'),
      h('button', { class: 'btn small ghost', title: 'Delete', onclick: async () => { if (await deleteSession(s.name)) refresh(); } }, icon('close', 13))))));
}
