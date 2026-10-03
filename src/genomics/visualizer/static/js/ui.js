// DOM helpers and small UI components (no framework).

export function h(tag, attrs, ...children) {
  const el = document.createElement(tag);
  if (attrs) {
    for (const [key, value] of Object.entries(attrs)) {
      if (value === undefined || value === null || value === false) continue;
      if (key === 'class') el.className = value;
      else if (key === 'style' && typeof value === 'object') Object.assign(el.style, value);
      else if (key === 'dataset') Object.assign(el.dataset, value);
      else if (key.startsWith('on') && typeof value === 'function') el.addEventListener(key.slice(2), value);
      else if (key === 'html') el.innerHTML = value;
      else if (value === true) el.setAttribute(key, '');
      else el.setAttribute(key, value);
    }
  }
  append(el, children);
  return el;
}

function append(el, children) {
  for (const child of children) {
    if (child === null || child === undefined || child === false) continue;
    if (Array.isArray(child)) append(el, child);
    else el.appendChild(child instanceof Node ? child : document.createTextNode(String(child)));
  }
}

export function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }

const ICONS = {
  menu: '<path d="M3 5h14M3 10h14M3 15h14"/>',
  overview: '<rect x="3" y="3" width="6" height="6" rx="1.5"/><rect x="11" y="3" width="6" height="6" rx="1.5"/><rect x="3" y="11" width="6" height="6" rx="1.5"/><rect x="11" y="11" width="6" height="6" rx="1.5"/>',
  samples: '<circle cx="7" cy="7" r="2.6"/><circle cx="14" cy="8" r="2.1"/><path d="M2.5 16c.6-2.8 2.4-4.2 4.5-4.2s3.9 1.4 4.5 4.2M11.5 12.6c.7-.5 1.5-.8 2.5-.8 1.8 0 3.1 1.2 3.6 3.7"/>',
  tracks: '<path d="M2 15c2 0 2-8 4-8s2 5 4 5 2-9 4-9 2 12 4 12"/>',
  sequence: '<path d="M4 3c0 4 12 4 12 7s-12 3-12 7M16 3c0 4-12 4-12 7s12 3 12 7M6 6.5h8M6 13.5h8"/>',
  experiments: '<path d="M7.5 3v5L3.5 15.5A1.2 1.2 0 0 0 4.6 17h10.8a1.2 1.2 0 0 0 1.1-1.5L12.5 8V3M6.5 3h7M5.5 12.5h9"/>',
  labs: '<path d="M10 2.5l1.8 4.3 4.7.4-3.6 3 1.1 4.6L10 12.4l-4 2.4 1.1-4.6-3.6-3 4.7-.4z"/>',
  sun: '<circle cx="10" cy="10" r="3.5"/><path d="M10 1.8v2M10 16.2v2M1.8 10h2M16.2 10h2M4.2 4.2l1.4 1.4M14.4 14.4l1.4 1.4M4.2 15.8l1.4-1.4M14.4 5.6l1.4-1.4"/>',
  moon: '<path d="M16.5 12.5A7 7 0 0 1 7.5 3.5a7 7 0 1 0 9 9z"/>',
  left: '<path d="M12.5 4.5L7 10l5.5 5.5"/>',
  right: '<path d="M7.5 4.5L13 10l-5.5 5.5"/>',
  plus: '<path d="M10 4v12M4 10h12"/>',
  minus: '<path d="M4 10h12"/>',
  close: '<path d="M5 5l10 10M15 5L5 15"/>',
  target: '<circle cx="10" cy="10" r="6"/><circle cx="10" cy="10" r="1.6"/>',
  expand: '<path d="M3 8V3h5M17 8V3h-5M3 12v5h5M17 12v5h-5"/>',
  panel: '<rect x="3" y="3.5" width="14" height="13" rx="2"/><path d="M12.5 3.5v13"/>',
  download: '<path d="M10 3v10M6 9.5l4 4 4-4M4 16.5h12"/>',
  pin: '<path d="M7 3h6l-1 5 3 3H5l3-3zM10 11v6"/>',
  external: '<path d="M11 3h6v6M17 3l-8 8M14 12v4H4V6h4"/>',
  refresh: '<path d="M16 10a6 6 0 1 1-1.8-4.3M16 3.5v3.5h-3.5"/>',
  search: '<circle cx="9" cy="9" r="5"/><path d="M13 13l4 4"/>',
};

export function icon(name, size = 18) {
  const span = document.createElement('span');
  span.style.display = 'inline-flex';
  span.innerHTML = `<svg width="${size}" height="${size}" viewBox="0 0 20 20" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${ICONS[name] || ''}</svg>`;
  return span.firstChild;
}

export function iconButton(name, title, onClick, cls = 'icon-btn') {
  return h('button', { class: cls, title, 'aria-label': title, type: 'button', onclick: onClick }, icon(name, 16));
}

export function segmented(options, value, onChange) {
  const root = h('div', { class: 'seg', role: 'radiogroup' });
  const buttons = options.map((opt) => {
    const b = h('button', {
      type: 'button', role: 'radio', title: opt.title || '', disabled: opt.disabled || false,
      onclick: () => { if (!opt.disabled) { set(opt.value); onChange(opt.value); } },
    }, opt.label);
    b.dataset.value = opt.value;
    root.appendChild(b);
    return b;
  });
  function set(v) {
    for (const b of buttons) {
      const on = b.dataset.value === String(v);
      b.classList.toggle('on', on);
      b.setAttribute('aria-checked', on ? 'true' : 'false');
    }
  }
  set(value);
  root.setValue = set;
  return root;
}

export function select(options, value, onChange, attrs = {}) {
  const el = h('select', { class: 'select', ...attrs, onchange: (e) => onChange(e.target.value) });
  setOptions(el, options, value);
  return el;
}

export function setOptions(el, options, value) {
  clear(el);
  for (const opt of options) {
    const o = typeof opt === 'object' ? opt : { value: opt, label: opt };
    const node = h('option', { value: o.value }, o.label);
    if (o.disabled) node.disabled = true;
    el.appendChild(node);
  }
  if (value !== undefined && value !== null) el.value = value;
}

export function field(label, control, attrs = {}) {
  return h('label', { class: 'field', ...attrs }, h('span', null, label), control);
}

// ---------- formatting ----------
export const fmtInt = (v) => (v === null || v === undefined || Number.isNaN(v) ? '–' : Math.round(v).toLocaleString('en-US'));

export function fmtNum(v, digits = 3) {
  if (v === null || v === undefined || !Number.isFinite(v)) return '–';
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(2);
  return Number(v.toPrecision(digits)).toLocaleString('en-US', { maximumFractionDigits: 6 });
}

export function fmtBp(n) {
  if (!Number.isFinite(n)) return '–';
  const a = Math.abs(n);
  if (a >= 1e6) return `${(n / 1e6).toFixed(a >= 1e7 ? 1 : 2)} Mb`;
  if (a >= 1e3) return `${(n / 1e3).toFixed(a >= 1e5 ? 0 : 1)} kb`;
  return `${Math.round(n)} bp`;
}

export function fmtPct(v, digits = 1) {
  return Number.isFinite(v) ? `${(v * 100).toFixed(digits)}%` : '–';
}

export function fmtDate(text) {
  if (!text) return '–';
  const d = new Date(text);
  return Number.isNaN(d.getTime()) ? String(text) : d.toLocaleString(undefined, { year: 'numeric', month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' });
}

export function escapeHtml(text) {
  return String(text ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

// ---------- tooltip ----------
const tooltipEl = () => document.getElementById('tooltip');

export function showTooltip(x, y, html) {
  const el = tooltipEl();
  el.innerHTML = html;
  el.hidden = false;
  const pad = 14;
  const rect = el.getBoundingClientRect();
  let left = x + pad;
  let top = y + pad;
  if (left + rect.width > window.innerWidth - 8) left = x - rect.width - pad;
  if (top + rect.height > window.innerHeight - 8) top = Math.max(8, y - rect.height - pad);
  el.style.left = `${Math.max(8, left)}px`;
  el.style.top = `${top}px`;
}

export function hideTooltip() { tooltipEl().hidden = true; }

// ---------- toast ----------
export function toast(message, kind = 'info', timeout = 4500) {
  const el = h('div', { class: `toast ${kind}`, role: kind === 'error' ? 'alert' : 'status' },
    kind === 'error' ? h('span', { class: 'status-dot bad', style: { marginTop: '5px' } }) : null,
    h('div', null, message));
  document.getElementById('toasts').appendChild(el);
  setTimeout(() => el.remove(), timeout);
}

// ---------- drawer & modal ----------
export function drawer(title, body, { onClose, actions } = {}) {
  const backdrop = h('div', { class: 'drawer-backdrop', onclick: () => close() });
  const panel = h('aside', { class: 'drawer', role: 'dialog', 'aria-label': title },
    h('div', { class: 'drawer-head' }, h('h2', null, title), h('div', { style: { display: 'flex', gap: '6px' } }, actions || null, iconButton('close', 'Close', () => close()))),
    h('div', { class: 'drawer-body' }, body));
  const onKey = (e) => { if (e.key === 'Escape') close(); };
  function close() {
    backdrop.remove(); panel.remove(); document.removeEventListener('keydown', onKey);
    if (onClose) onClose();
  }
  document.addEventListener('keydown', onKey);
  document.body.append(backdrop, panel);
  return { close, body: panel.querySelector('.drawer-body') };
}

export function modal(title, body, buttons = []) {
  const backdrop = h('div', { class: 'modal-backdrop', onclick: () => close() });
  const foot = h('div', { class: 'modal-foot' }, buttons);
  const panel = h('div', { class: 'modal', role: 'dialog', 'aria-label': title },
    h('div', { class: 'drawer-head' }, h('h2', null, title), iconButton('close', 'Close', () => close())),
    h('div', { class: 'modal-body' }, body), foot);
  const onKey = (e) => { if (e.key === 'Escape') close(); };
  function close() { backdrop.remove(); panel.remove(); document.removeEventListener('keydown', onKey); }
  document.addEventListener('keydown', onKey);
  document.body.append(backdrop, panel);
  return { close, foot };
}

export function progressBar(fraction) {
  const inner = h('i', { style: { width: `${Math.round((fraction || 0) * 100)}%` } });
  const el = h('div', { class: 'bar', role: 'progressbar', 'aria-valuemin': '0', 'aria-valuemax': '100' }, inner);
  el.set = (f) => { inner.style.width = `${Math.round(Math.max(0, Math.min(1, f)) * 100)}%`; el.setAttribute('aria-valuenow', String(Math.round(f * 100))); };
  return el;
}

/** Overlay with a progress bar for long-running server jobs (optionally cancellable). */
export function jobOverlay(host, title, onCancel = null) {
  const bar = progressBar(0);
  const message = h('div', { class: 'muted', style: { fontSize: '12px' } }, 'Starting…');
  let jobId = null;
  const cancel = onCancel ? h('button', { class: 'btn small', onclick: () => { if (jobId) onCancel(jobId); } }, 'Cancel') : null;
  const card = h('div', { class: 'overlay-card' },
    h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center' } }, h('span', { class: 'spinner' }), h('b', { style: { flex: '1' } }, title), cancel),
    bar, message,
    h('div', { class: 'muted', style: { fontSize: '11.5px' } }, 'Results are cached; the job keeps running if you leave this page.'));
  const el = h('div', { class: 'overlay' }, card);
  host.appendChild(el);
  return {
    update(job) { jobId = job.id; bar.set(job.progress || 0); message.textContent = `${job.message || ''} · ${Math.round((job.progress || 0) * 100)}% · ${Math.round(job.elapsed || 0)}s`; },
    remove() { el.remove(); },
  };
}

export function debounce(fn, ms) {
  let t = null;
  const wrapped = (...args) => { clearTimeout(t); t = setTimeout(() => fn(...args), ms); };
  wrapped.cancel = () => clearTimeout(t);
  return wrapped;
}

export function downloadText(filename, text, type = 'text/plain') {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const a = h('a', { href: url, download: filename });
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function errorBox(err) {
  return h('div', { class: 'error-box' }, err && err.message ? err.message : String(err));
}
