// AlphaGenome backend used by the Perturbation Lab and by background prediction jobs: the hosted
// API, a remote server, or a server started on this machine.
import { api } from '../api.js';
import { navigate } from '../app.js';
import { h, clear, icon, toast, errorBox, segmented, field } from '../ui.js';

const LOCAL_STATES = {
  stopped: ['bad', 'Stopped'],
  starting: ['warn', 'Loading model…'],
  ready: ['good', 'Ready'],
  external: ['good', 'Running (started outside the visualizer)'],
  exited: ['bad', 'Exited'],
};

export async function mount(root) {
  const page = h('div', { class: 'page-inner' });
  root.appendChild(page);
  const backendHost = h('div');
  const list = h('div', { class: 'grid' });
  page.append(
    h('div', { class: 'page-head' }, h('div', null, h('h1', null, 'AlphaGenome'), h('p', null, 'Where AlphaGenome predictions are computed: the Perturbation Lab calls it directly and prediction jobs use it in the background.'))),
    h('div', { class: 'grid', style: { maxWidth: '900px' } }, backendHost, list));

  let pollTimer = null;
  let alive = true;
  const schedulePoll = (backend) => {
    clearTimeout(pollTimer);
    if (alive && backend && backend.local && backend.local.state === 'starting') pollTimer = setTimeout(refresh, 2000);
  };

  // ---------------- AlphaGenome backend ----------------
  const form = { mode: 'cloud', address: '', ca_cert: '' };
  let saved = null;
  const address = h('input', { class: 'input mono', placeholder: 'grpc://10.0.0.5:50051  ·  grpcs://host:50051  ·  host:port', oninput: (e) => { form.address = e.target.value; markDirty(); clear(result); } });
  const caCert = h('input', { class: 'input mono', placeholder: '/path/to/alphagenome_research/certs/ca.crt (only for a self-signed TLS server)', oninput: (e) => { form.ca_cert = e.target.value; markDirty(); clear(result); } });
  const modeSeg = segmented([
    { value: 'cloud', label: 'Hosted API' },
    { value: 'remote', label: 'Remote server' },
    { value: 'local', label: 'This machine' },
  ], form.mode, (v) => { form.mode = v; markDirty(); clear(result); renderPanel(); });
  const panel = h('div', { style: { display: 'grid', gap: '10px' } });
  const result = h('div', { class: 'muted', style: { fontSize: '12px' } });
  const statusDot = h('span', { class: 'status-dot' });
  const activeLabel = h('span', { class: 'muted', style: { fontSize: '12px' } });
  const saveBtn = h('button', { class: 'btn', onclick: () => save().catch(() => {}) }, 'Save');
  const testBtn = h('button', { class: 'btn', onclick: () => test(false), title: 'Connect and call the AlphaGenome service (fast)' }, 'Test connection');
  const predictBtn = h('button', { class: 'btn', onclick: () => test(true), title: 'Run a 16 kb prediction (the first one on a fresh local server includes compilation and can take minutes)' }, 'Test prediction');
  let backend = null;

  const dirty = () => saved && (form.mode !== saved.mode || form.address !== saved.address || form.ca_cert !== saved.ca_cert);
  function markDirty() { saveBtn.classList.toggle('primary', !!dirty()); }

  function localPanel(local) {
    const [cls, text] = LOCAL_STATES[local.state] || ['warn', local.state];
    const running = local.state === 'starting' || local.state === 'ready';
    const start = async () => {
      try {
        if (form.mode !== 'local' || dirty()) { form.mode = 'local'; modeSeg.setValue('local'); await save(true); }
        render(await api('/api/alphagenome/local/start', { method: 'POST', body: {} }));
        toast('AlphaGenome server starting: loading the model can take a few minutes.');
      } catch (err) { toast(err.message, 'error', 9000); }
    };
    const stop = async () => {
      try { render(await api('/api/alphagenome/local/stop', { method: 'POST', body: {} })); } catch (err) { toast(err.message, 'error', 8000); }
    };
    const meta = (k, v) => h('div', null, h('span', { class: 'muted' }, `${k}: `), h('span', { class: 'mono' }, v || '–'));
    return h('div', { style: { display: 'grid', gap: '10px' } },
      h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center', justifyContent: 'space-between', flexWrap: 'wrap' } },
        h('div', { style: { display: 'flex', gap: '8px', alignItems: 'center' } }, h('span', { class: `status-dot ${cls}` }), h('b', null, text),
          local.state === 'exited' && local.exit_code !== null ? h('span', { class: 'muted' }, `(code ${local.exit_code})`) : null,
          local.pid ? h('span', { class: 'muted' }, `pid ${local.pid}`) : null),
        running
          ? h('button', { class: 'btn', onclick: stop }, 'Stop server')
          : h('button', { class: 'btn primary', disabled: local.state === 'external' || local.reasons.length > 0, onclick: start }, 'Start server')),
      h('div', { style: { display: 'grid', gap: '2px', fontSize: '12px' } },
        meta('Address', `${local.address}${local.tls ? ' (TLS)' : ' (plaintext)'}`),
        meta('Checkout', local.server_dir),
        meta('Python', local.python),
        local.environment ? meta('Packages', local.environment) : null),
      local.reasons.length ? h('div', { class: 'notice' }, h('b', null, 'Cannot start: '), local.reasons.join('; '),
        h('div', { style: { marginTop: '6px' } }, 'Install everything (checkout, environment with CUDA jax, weights check) from a terminal with ',
          h('span', { class: 'mono' }, local.setup_command || 'genomics alphagenome server setup'), ', then reload this page.')) : null,
      h('div', { class: 'muted', style: { fontSize: '12px' } }, local.host === '0.0.0.0'
        ? 'The server listens on all interfaces (0.0.0.0) and uses this machine\'s GPU; other machines can use it as a remote server (no authentication). It stops when the visualizer exits.'
        : h('span', null, `The server listens on ${local.host} (this machine only; set ALPHAGENOME_SERVER_HOST=0.0.0.0 to share it) and uses this machine's GPU. Loading the model takes 1–3 minutes; the first prediction also compiles it. It stops when the visualizer exits. A server started with `,
          h('span', { class: 'mono' }, 'genomics alphagenome server start'), ' is detected and used as is.')),
      local.log_tail.length ? h('details', { class: 'collapsible', open: local.state !== 'ready' },
        h('summary', null, 'Server log'),
        h('pre', { style: { maxHeight: '220px', overflow: 'auto', padding: '8px', background: 'var(--surface-2)', borderRadius: '6px' } }, local.log_tail.join('\n')),
        h('div', { class: 'muted', style: { fontSize: '12px' } }, 'Full log: ', h('span', { class: 'mono' }, local.log))) : null);
  }

  function renderPanel() {
    clear(panel);
    if (!backend) return;
    if (form.mode === 'cloud') {
      panel.append(h('p', { style: { margin: 0 } }, 'Uses Google\'s hosted AlphaGenome API with ', h('span', { class: 'mono' }, 'ALPHAGENOME_API_KEY'), ' from the environment or ~/.env.'),
        backend.api_key_available ? h('div', { class: 'muted' }, 'API key found.') : h('div', { class: 'notice' }, h('b', null, 'No API key found. '), 'Set ALPHAGENOME_API_KEY before starting the visualizer, or use a self-hosted server.'));
    } else if (form.mode === 'remote') {
      panel.append(
        field('Server address', address),
        field('CA certificate (optional)', caCert),
        h('div', { class: 'muted', style: { fontSize: '12px' } }, h('span', { class: 'mono' }, 'grpc://'), ' = plaintext, ', h('span', { class: 'mono' }, 'grpcs://'), ' = TLS, no scheme = detect. The port defaults to 50051. A self-hosted server ignores the API key.'));
    } else {
      panel.append(localPanel(backend.local));
    }
    const canTest = form.mode !== 'local' || ['ready', 'external'].includes(backend.local.state);
    testBtn.disabled = !canTest;
    predictBtn.disabled = !canTest;
  }

  function render(data) {
    const keepEdits = saved && dirty();
    backend = data;
    saved = { ...data.settings };
    if (!keepEdits) Object.assign(form, saved);
    if (document.activeElement !== address) address.value = form.address;
    if (document.activeElement !== caCert) caCert.value = form.ca_cert;
    modeSeg.setValue(form.mode);
    statusDot.className = `status-dot ${data.reasons.length ? 'bad' : 'good'}`;
    activeLabel.textContent = data.reasons.length ? `${data.label}: ${data.reasons.join('; ')}` : `In use: ${data.label}`;
    markDirty();
    renderPanel();
    schedulePoll(data);
  }

  async function save(quiet = false) {
    try {
      render(await api('/api/alphagenome/settings', { method: 'POST', body: { ...form } }));
      if (!quiet) toast('AlphaGenome backend saved. Jobs already running keep the backend they were started with.');
      renderLabs();
    } catch (err) { toast(err.message, 'error', 8000); throw err; }
  }

  async function test(predict) {
    const buttons = [testBtn, predictBtn];
    buttons.forEach((b) => { b.disabled = true; });
    clear(result).append(h('span', { class: 'muted' }, predict ? 'Running a test prediction…' : 'Connecting…'));
    try {
      const res = await api('/api/alphagenome/test', { method: 'POST', body: { ...form, predict, timeout: 10 } });
      clear(result).append(h('span', { class: `status-dot ${res.ok ? 'good' : 'bad'}`, style: { marginRight: '6px' } }), res.message);
    } catch (err) {
      clear(result).append(h('span', { class: 'status-dot bad', style: { marginRight: '6px' } }), err.message);
    } finally { renderPanel(); }
  }

  backendHost.appendChild(h('section', { class: 'card' },
    h('div', { class: 'card-head' },
      h('div', { style: { display: 'flex', gap: '10px', alignItems: 'center' } }, statusDot, h('h2', null, 'AlphaGenome backend')),
      modeSeg),
    h('div', { class: 'card-body', style: { display: 'grid', gap: '12px' } },
      panel,
      h('div', { style: { display: 'flex', gap: '8px', flexWrap: 'wrap', alignItems: 'center' } }, saveBtn, testBtn, predictBtn),
      result,
      activeLabel)));

  // ---------------- where it is used ----------------
  function renderLabs() {
    clear(list).appendChild(h('section', { class: 'card' },
      h('div', { class: 'card-head' }, h('h2', null, 'Used by')),
      h('div', { class: 'card-body', style: { display: 'grid', gap: '10px' } },
        h('div', { style: { display: 'flex', gap: '12px', alignItems: 'center', justifyContent: 'space-between' } },
          h('div', null, h('b', null, 'Perturbation Lab'), h('div', { class: 'muted', style: { fontSize: '12px' } }, 'Re-predicts edited haplotypes and re-scores them with a trained model, inside this app.')),
          h('button', { class: 'btn', onclick: () => navigate('perturb') }, icon('perturb', 14), 'Open')),
        h('div', { style: { display: 'flex', gap: '12px', alignItems: 'center', justifyContent: 'space-between' } },
          h('div', null, h('b', null, 'Prediction jobs'), h('div', { class: 'muted', style: { fontSize: '12px' } }, 'Predict tracks for every individual of a dataset (Overview → AlphaGenome predictions, or after an import).')),
          h('button', { class: 'btn', onclick: () => navigate('jobs') }, icon('jobs', 14), 'Jobs')))));
  }

  async function refresh() {
    try { render(await api('/api/alphagenome')); } catch (err) { clear(backendHost).appendChild(errorBox(err)); return; }
    renderLabs();
  }

  refresh();
  return { unmount() { alive = false; clearTimeout(pollTimer); } };
}
