// JSON API client: decodes base64 float32 payloads and polls background jobs.

function decodeF32(obj) {
  const bin = atob(obj.$f32);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  const flat = new Float32Array(bytes.buffer);
  const shape = obj.shape || [flat.length];
  if (shape.length <= 1) return flat;
  const rows = shape[0];
  const cols = shape.slice(1).reduce((a, b) => a * b, 1);
  const out = new Array(rows);
  for (let r = 0; r < rows; r++) out[r] = flat.subarray(r * cols, (r + 1) * cols);
  out.flat = flat;
  out.shape = shape;
  return out;
}

function decode(value) {
  if (Array.isArray(value)) return value.map(decode);
  if (value && typeof value === 'object') {
    if (typeof value.$f32 === 'string') return decodeF32(value);
    const out = {};
    for (const [k, v] of Object.entries(value)) out[k] = decode(v);
    return out;
  }
  return value;
}

export class ApiError extends Error {
  constructor(message, status) { super(message); this.status = status; }
}

export async function api(path, { method = 'GET', body, signal, params } = {}) {
  let url = path;
  if (params) {
    const qs = new URLSearchParams();
    for (const [k, v] of Object.entries(params)) {
      if (v === undefined || v === null || v === '') continue;
      qs.set(k, typeof v === 'object' ? JSON.stringify(v) : String(v));
    }
    const s = qs.toString();
    if (s) url += (url.includes('?') ? '&' : '?') + s;
  }
  const res = await fetch(url, {
    method,
    signal,
    headers: body ? { 'Content-Type': 'application/json' } : undefined,
    body: body ? JSON.stringify(body) : undefined,
  });
  let data;
  try { data = await res.json(); } catch (e) { throw new ApiError(`${res.status} ${res.statusText}`, res.status); }
  if (!res.ok || (data && data.error)) throw new ApiError((data && data.error) || res.statusText, res.status);
  return decode(data);
}

const sleep = (ms, signal) => new Promise((resolve, reject) => {
  const t = setTimeout(resolve, ms);
  if (signal) signal.addEventListener('abort', () => { clearTimeout(t); reject(new DOMException('Aborted', 'AbortError')); }, { once: true });
});

/** Request an endpoint that may answer {pending, job}; poll until it returns data. */
export async function apiJob(path, { params, signal, onProgress, method = 'GET', body } = {}) {
  let delay = 350;
  for (;;) {
    const data = await api(path, { params, signal, method, body });
    if (!data || !data.pending) return data;
    if (onProgress) onProgress(data.job);
    await sleep(delay, signal);
    delay = Math.min(1200, delay * 1.25);
  }
}

export const isAbort = (err) => err && (err.name === 'AbortError' || err.code === 20);

/** Keeps only the latest request alive: starting a new one aborts the previous. */
export class Latest {
  constructor() { this.controller = null; }
  next() {
    if (this.controller) this.controller.abort();
    this.controller = new AbortController();
    return this.controller.signal;
  }
  abort() { if (this.controller) this.controller.abort(); this.controller = null; }
}
