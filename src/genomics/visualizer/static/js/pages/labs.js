// Labs: optional heavyweight tools launched on demand (separate processes).
import { api } from '../api.js';
import { h, clear, icon, toast, errorBox } from '../ui.js';

export async function mount(root) {
  const page = h('div', { class: 'page-inner' });
  root.appendChild(page);
  const list = h('div', { class: 'grid', style: { maxWidth: '900px' } });
  page.append(h('div', { class: 'page-head' }, h('div', null, h('h1', null, 'Labs'), h('p', null, 'Tools with heavy dependencies (trained checkpoints, AlphaGenome API) run as separate processes, started only when you need them.'))), list);

  async function render() {
    let labs;
    try { ({ labs } = await api('/api/labs')); } catch (err) { clear(list).appendChild(errorBox(err)); return; }
    clear(list);
    for (const lab of labs) {
      const start = async () => {
        try {
          await api(`/api/labs/${lab.key}/start`, { method: 'POST', body: {} });
          toast(`${lab.title} starting at ${lab.url} (models load in a few seconds)`);
          setTimeout(() => window.open(lab.url, '_blank', 'noopener'), 2500);
          render();
        } catch (err) { toast(err.message, 'error', 8000); }
      };
      list.appendChild(h('section', { class: 'card' },
        h('div', { class: 'card-head' },
          h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center' } }, h('span', { class: `status-dot ${lab.running ? 'good' : lab.available ? 'warn' : 'bad'}` }), h('h2', null, lab.title)),
          lab.running
            ? h('a', { class: 'btn primary', href: lab.url, target: '_blank', rel: 'noopener' }, 'Open', icon('external', 14))
            : h('button', { class: 'btn primary', disabled: !lab.available, onclick: start }, 'Launch')),
        h('div', { class: 'card-body', style: { display: 'grid', gap: '10px' } },
          h('p', { style: { margin: 0 } }, lab.description),
          lab.running ? h('div', { class: 'notice' }, `Running at ${lab.url}`) : null,
          lab.reasons && lab.reasons.length ? h('div', { class: 'notice' }, h('b', null, 'Unavailable: '), lab.reasons.join('; ')) : null,
          h('div', { class: 'muted', style: { fontSize: '12px' } }, 'Log: ', h('span', { class: 'mono' }, lab.log)))));
    }
  }
  render();
  return {};
}
