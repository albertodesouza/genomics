// Canvas plotting primitives shared by the Tracks, Sequence and Experiments pages.
import { exportState, SvgContext } from './figure.js';

export function css(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

export function theme() {
  return {
    surface: css('--surface'), surface2: css('--surface-2'), ink: css('--ink'), ink2: css('--ink-2'), ink3: css('--ink-3'),
    line: css('--line'), lineStrong: css('--line-strong'), accent: css('--accent'),
    modelWindow: css('--model-window'), modelWindowEdge: css('--model-window-edge'),
    font: css('--font'), mono: css('--mono'),
    bases: { A: css('--base-a'), C: css('--base-c'), G: css('--base-g'), T: css('--base-t'), N: css('--base-n') },
    seq: [0, 1, 2, 3, 4, 5, 6].map((i) => hexToRgb(css(`--seq-${i}`))),
  };
}

export function hexToRgb(hex) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex || '');
  if (!m) return [128, 128, 128];
  const n = parseInt(m[1], 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

export function withAlpha(color, alpha) {
  if (color.startsWith('#')) {
    const [r, g, b] = hexToRgb(color);
    return `rgba(${r},${g},${b},${alpha})`;
  }
  return color;
}

/** Size a canvas for the device pixel ratio; returns a 2D context in CSS pixels. */
export function setupCanvas(canvas, width, height) {
  const w = Math.max(1, Math.floor(width));
  const hgt = Math.max(1, Math.floor(height));
  if (exportState.recorders) {
    // SVG export: record the drawing instead of painting (the canvas keeps its pixels).
    const rec = new SvgContext(w, hgt);
    exportState.recorders.set(canvas, rec);
    return rec;
  }
  const dpr = exportState.dpr || window.devicePixelRatio || 1;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(hgt * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(hgt * dpr);
  }
  canvas.style.width = `${w}px`;
  canvas.style.height = `${hgt}px`;
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, hgt);
  return ctx;
}

export function niceStep(span, target) {
  const raw = span / Math.max(1, target);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10;
  return step * mag;
}

export function ticks(min, max, target = 6) {
  if (!(max > min)) return [min];
  const step = niceStep(max - min, target);
  const out = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) out.push(Number(v.toPrecision(12)));
  return out;
}

/** Sequential color for t in [0,1] along the theme's single-hue ramp. */
export function seqColor(stops, t) {
  if (!Number.isFinite(t)) return null;
  const x = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i;
  const a = stops[i];
  const b = stops[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

/** Stable color slots: an entity keeps its color while it stays on screen. */
export class ColorSlots {
  constructor(size = 8) { this.size = size; this.map = new Map(); }
  sync(keys) {
    const keep = new Set(keys);
    for (const k of [...this.map.keys()]) if (!keep.has(k)) this.map.delete(k);
    const used = new Set(this.map.values());
    for (const k of keys) {
      if (this.map.has(k)) continue;
      let slot = 0;
      while (used.has(slot) && slot < this.size) slot++;
      if (slot >= this.size) slot = this.map.size % this.size;
      this.map.set(k, slot);
      used.add(slot);
    }
  }
  slot(key) { return this.map.get(key) ?? 0; }
  color(key) { return css(`--series-${(this.slot(key) % this.size) + 1}`); }
}

/**
 * Viewport over a 1-D domain [0, length). Handles drag-to-pan, ctrl/⌘+wheel zoom,
 * shift+wheel / horizontal-trackpad pan, double-click zoom, shift+drag region select
 * and keyboard (←/→ pan, +/− zoom) on any element attached via attach().
 */
export class Viewport {
  constructor({ length = 1, start = 0, end = 1, minSpan = 20, onChange }) {
    this.length = Math.max(1, length);
    this.minSpan = minSpan;
    this.onChange = onChange;
    this.start = start;
    this.end = end;
    this.clamp();
    this._raf = null;
  }

  get span() { return this.end - this.start; }

  setDomain(length) { this.length = Math.max(1, length); this.clamp(); }

  clamp() {
    let span = Math.max(this.minSpan, Math.min(this.end - this.start, this.length));
    if (!Number.isFinite(span) || span <= 0) span = this.length;
    let start = Math.round(Math.max(0, Math.min(this.start, this.length - span)));
    this.start = start;
    this.end = Math.round(Math.min(this.length, start + span));
  }

  set(start, end, emit = true) {
    this.start = start;
    this.end = end;
    this.clamp();
    if (emit) this.changed();
  }

  pan(deltaBp) { this.set(this.start + deltaBp, this.end + deltaBp); }

  zoom(factor, anchorFrac = 0.5) {
    const span = this.span;
    const anchor = this.start + span * anchorFrac;
    const next = Math.max(this.minSpan, Math.min(this.length, span * factor));
    this.set(anchor - next * anchorFrac, anchor - next * anchorFrac + next);
  }

  center(pos, span = this.span) { this.set(pos - span / 2, pos + span / 2); }

  changed() {
    if (this._raf) return;
    this._raf = requestAnimationFrame(() => { this._raf = null; if (this.onChange) this.onChange(this); });
  }

  /** Attach interactions. `geom()` returns {left, width} of the plotted x-range within el. */
  attach(el, geom, { onHover, onLeave, onSelect } = {}) {
    let drag = null;
    const fracAt = (clientX) => {
      const rect = el.getBoundingClientRect();
      const g = geom();
      return (clientX - rect.left - g.left) / g.width;
    };
    el.classList.add('interact');
    el.tabIndex = el.tabIndex >= 0 ? el.tabIndex : 0;
    el.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      el.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, start: this.start, end: this.end, select: e.shiftKey || e.altKey, frac0: fracAt(e.clientX), moved: false };
      if (!drag.select) el.classList.add('dragging');
    });
    el.addEventListener('pointermove', (e) => {
      if (drag) {
        const g = geom();
        const dx = e.clientX - drag.x;
        if (Math.abs(dx) > 2) drag.moved = true;
        if (drag.select) {
          if (onSelect) onSelect(drag.frac0, fracAt(e.clientX), false);
        } else {
          const bpPerPx = (drag.end - drag.start) / g.width;
          this.set(drag.start - dx * bpPerPx, drag.end - dx * bpPerPx);
        }
        if (onLeave) onLeave();
        return;
      }
      if (onHover) onHover(e, fracAt(e.clientX));
    });
    const finish = (e) => {
      if (!drag) return;
      el.classList.remove('dragging');
      if (drag.select && drag.moved) {
        const a = Math.min(drag.frac0, fracAt(e.clientX));
        const b = Math.max(drag.frac0, fracAt(e.clientX));
        if (onSelect) onSelect(a, b, true);
        const span = this.span;
        const s = this.start;
        if (b - a > 0.002) this.set(s + a * span, s + b * span);
      }
      drag = null;
    };
    el.addEventListener('pointerup', finish);
    el.addEventListener('pointercancel', finish);
    el.addEventListener('pointerleave', () => { if (!drag && onLeave) onLeave(); });
    el.addEventListener('wheel', (e) => {
      const horizontal = Math.abs(e.deltaX) > Math.abs(e.deltaY) || e.shiftKey;
      if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
        const factor = Math.exp(Math.max(-1, Math.min(1, e.deltaY * 0.01)));
        this.zoom(factor, Math.max(0, Math.min(1, fracAt(e.clientX))));
      } else if (horizontal) {
        e.preventDefault();
        const g = geom();
        const delta = e.shiftKey && !e.deltaX ? e.deltaY : e.deltaX;
        this.pan((delta * this.span) / g.width);
      }
    }, { passive: false });
    el.addEventListener('dblclick', (e) => {
      this.zoom(e.altKey ? 2 : 0.5, Math.max(0, Math.min(1, fracAt(e.clientX))));
    });
    el.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowLeft') { this.pan(-this.span * 0.15); e.preventDefault(); }
      else if (e.key === 'ArrowRight') { this.pan(this.span * 0.15); e.preventDefault(); }
      else if (e.key === '+' || e.key === '=') { this.zoom(0.5); e.preventDefault(); }
      else if (e.key === '-' || e.key === '_') { this.zoom(2); e.preventDefault(); }
    });
  }
}

/** Observe an element's width and call fn on change (debounced to animation frames). */
export function onResize(el, fn) {
  let last = -1;
  let raf = null;
  const ro = new ResizeObserver(() => {
    const w = el.clientWidth;
    if (w === last) return;
    last = w;
    if (raf) cancelAnimationFrame(raf);
    raf = requestAnimationFrame(() => fn(w));
  });
  ro.observe(el);
  return () => ro.disconnect();
}

/** Draw a simple multi-series line chart (experiments: training curves). */
export function lineChart(canvas, { series, width, height, xLabel = '', yLabel = '', yMin = null, yMax = null }) {
  const t = theme();
  const ctx = setupCanvas(canvas, width, height);
  const pad = { l: 48, r: 14, t: 12, b: 30 };
  const pw = width - pad.l - pad.r;
  const ph = height - pad.t - pad.b;
  let xmin = Infinity; let xmax = -Infinity; let ymin = Infinity; let ymax = -Infinity;
  for (const s of series) {
    s.y.forEach((v, i) => {
      if (!Number.isFinite(v)) return;
      const x = s.x ? s.x[i] : i;
      xmin = Math.min(xmin, x); xmax = Math.max(xmax, x); ymin = Math.min(ymin, v); ymax = Math.max(ymax, v);
    });
  }
  if (!Number.isFinite(xmin)) { ctx.fillStyle = t.ink3; ctx.font = `12px ${t.font}`; ctx.fillText('No data', pad.l, pad.t + 20); return null; }
  if (yMin !== null) ymin = yMin;
  if (yMax !== null) ymax = yMax;
  if (ymax === ymin) { ymax += 1; ymin -= 1; }
  const yr = ymax - ymin;
  ymin -= yr * 0.04; ymax += yr * 0.04;
  const sx = (x) => pad.l + ((x - xmin) / Math.max(1e-9, xmax - xmin)) * pw;
  const sy = (y) => pad.t + (1 - (y - ymin) / (ymax - ymin)) * ph;
  ctx.font = `11px ${t.font}`;
  ctx.textBaseline = 'middle';
  ctx.lineWidth = 1;
  for (const v of ticks(ymin, ymax, 4)) {
    const y = Math.round(sy(v)) + 0.5;
    ctx.strokeStyle = t.line; ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(width - pad.r, y); ctx.stroke();
    ctx.fillStyle = t.ink3; ctx.textAlign = 'right'; ctx.fillText(formatTick(v), pad.l - 6, y);
  }
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (const v of ticks(xmin, xmax, 6)) {
    ctx.fillStyle = t.ink3; ctx.fillText(formatTick(v), sx(v), height - pad.b + 7);
  }
  ctx.strokeStyle = t.lineStrong; ctx.beginPath(); ctx.moveTo(pad.l, pad.t + ph + 0.5); ctx.lineTo(width - pad.r, pad.t + ph + 0.5); ctx.stroke();
  if (xLabel) { ctx.fillStyle = t.ink3; ctx.textAlign = 'right'; ctx.fillText(xLabel, width - pad.r, height - 12); }
  if (yLabel) { ctx.save(); ctx.fillStyle = t.ink3; ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillText(yLabel, pad.l + 4, pad.t); ctx.restore(); }
  for (const s of series) {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2; ctx.setLineDash(s.dash || []);
    ctx.beginPath();
    let pen = false;
    s.y.forEach((v, i) => {
      if (!Number.isFinite(v)) { pen = false; return; }
      const x = sx(s.x ? s.x[i] : i); const y = sy(v);
      if (pen) ctx.lineTo(x, y); else { ctx.moveTo(x, y); pen = true; }
    });
    ctx.stroke();
  }
  ctx.setLineDash([]);
  return { sx, sy, pad, xmin, xmax, pw, ph };
}

export function formatTick(v) {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 1e4 || a < 1e-3) return v.toExponential(1);
  return Number(v.toPrecision(4)).toString();
}

/** Axis label with just enough digits to tell ``step``-spaced ticks apart. */
export function tickLabel(v, step) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  const stepExp = Math.floor(Math.log10(Math.abs(step) || 1));
  if (a >= 1e-2 && a < 1e5) return v.toFixed(Math.max(0, -stepExp));
  const digits = Math.max(0, Math.min(6, Math.floor(Math.log10(a)) - stepExp));
  return v.toExponential(digits);
}

/**
 * Values by genotype: box (quartiles, median), whiskers (min–max) and jittered points per group.
 * ``groups``: [{label, values, summary: {n, q1, median, q3, min, max, mean}}].
 * ``fit``: {intercept, slope} drawn across the group centres (dosage 0, 1, 2).
 * ``refLine``: a dashed horizontal line (e.g. AlphaGenome on the reference genome).
 * Returns the geometry (for hover) or null when there is nothing to draw.
 */
export function genotypePlot(canvas, { groups, width, height = 240, yLabel = '', fit = null, refLine = null, refLabel = '', color = null, log = false }) {
  const t = theme();
  const ctx = setupCanvas(canvas, width, height);
  const accent = color || t.accent;
  const L = 64; const R = 16; const T = 14; const B = 40;
  const tf = (v) => (log ? Math.log10(Math.max(v, 1e-12)) : v);
  const all = [];
  for (const g of groups) for (const v of g.values || []) if (Number.isFinite(v)) all.push(tf(v));
  if (refLine !== null && Number.isFinite(refLine)) all.push(tf(refLine));
  ctx.font = `12px ${t.font}`;
  if (!all.length) {
    ctx.fillStyle = t.ink3; ctx.textAlign = 'center'; ctx.fillText('No values', width / 2, height / 2);
    return null;
  }
  let lo = Math.min(...all); let hi = Math.max(...all);
  if (hi === lo) { hi += Math.abs(hi) * 0.1 || 1; lo -= Math.abs(lo) * 0.1 || 1; }
  const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
  const y = (v) => T + (1 - (tf(v) - lo) / (hi - lo)) * (height - T - B);
  const yRaw = (u) => T + (1 - (u - lo) / (hi - lo)) * (height - T - B);
  const plotW = width - L - R;
  const slot = plotW / groups.length;
  const cx = (i) => L + slot * (i + 0.5);
  // grid + y ticks
  ctx.strokeStyle = t.line; ctx.lineWidth = 1; ctx.fillStyle = t.ink3; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.font = `11px ${t.font}`;
  const yTicks = ticks(lo, hi, 5);
  const step = yTicks.length > 1 ? yTicks[1] - yTicks[0] : Math.abs(hi - lo) || 1;
  for (const v of yTicks) {
    const yy = Math.round(yRaw(v)) + 0.5;
    ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(width - R, yy); ctx.stroke();
    ctx.fillText(log ? `1e${Number(v.toFixed(1))}` : tickLabel(v, step), L - 6, yy);
  }
  if (yLabel) {
    ctx.save(); ctx.translate(12, T + (height - T - B) / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center'; ctx.fillStyle = t.ink2; ctx.fillText(yLabel, 0, 0); ctx.restore();
  }
  // points (deterministic jitter so redraws and exports match)
  const boxW = Math.min(70, slot * 0.42);
  groups.forEach((g, i) => {
    const vals = (g.values || []).filter(Number.isFinite);
    ctx.fillStyle = withAlpha(accent, vals.length > 600 ? 0.18 : 0.35);
    for (let k = 0; k < vals.length; k++) {
      const jitter = (((k * 2654435761) >>> 0) / 4294967296 - 0.5) * boxW * 1.5;
      ctx.fillRect(cx(i) + jitter - 1.2, y(vals[k]) - 1.2, 2.4, 2.4);
    }
  });
  // boxes
  groups.forEach((g, i) => {
    const s = g.summary || {};
    if (!s.n) return;
    const x0 = cx(i) - boxW / 2;
    ctx.strokeStyle = t.ink; ctx.lineWidth = 1.25;
    ctx.beginPath(); ctx.moveTo(cx(i), y(s.max)); ctx.lineTo(cx(i), y(s.q3)); ctx.moveTo(cx(i), y(s.q1)); ctx.lineTo(cx(i), y(s.min)); ctx.stroke();
    ctx.fillStyle = withAlpha(t.surface, 0.55);
    ctx.fillRect(x0, y(s.q3), boxW, y(s.q1) - y(s.q3));
    ctx.strokeRect(x0, y(s.q3), boxW, y(s.q1) - y(s.q3));
    ctx.lineWidth = 2.5; ctx.beginPath(); ctx.moveTo(x0, y(s.median)); ctx.lineTo(x0 + boxW, y(s.median)); ctx.stroke();
  });
  // regression line through the dosage positions
  if (fit && Number.isFinite(fit.slope) && Number.isFinite(fit.intercept)) {
    ctx.strokeStyle = accent; ctx.lineWidth = 2; ctx.setLineDash([]);
    ctx.beginPath(); ctx.moveTo(cx(0), y(fit.intercept)); ctx.lineTo(cx(groups.length - 1), y(fit.intercept + fit.slope * (groups.length - 1))); ctx.stroke();
  }
  if (refLine !== null && Number.isFinite(refLine)) {
    ctx.strokeStyle = t.ink; ctx.lineWidth = 1.25; ctx.setLineDash([5, 4]);
    const yy = y(refLine);
    ctx.beginPath(); ctx.moveTo(L, yy); ctx.lineTo(width - R, yy); ctx.stroke(); ctx.setLineDash([]);
    if (refLabel) {
      ctx.font = `11px ${t.font}`;
      const w = ctx.measureText(refLabel).width;
      ctx.fillStyle = withAlpha(t.surface, 0.85); ctx.fillRect(width - R - w - 6, yy - 17, w + 6, 14);
      ctx.fillStyle = t.ink2; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom'; ctx.fillText(refLabel, width - R - 3, yy - 4);
    }
  }
  // labels
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  groups.forEach((g, i) => {
    ctx.fillStyle = t.ink; ctx.font = `600 12px ${t.font}`; ctx.fillText(g.label, cx(i), height - B + 8);
    ctx.fillStyle = t.ink3; ctx.font = `11px ${t.font}`; ctx.fillText(`n=${(g.summary && g.summary.n) || 0}`, cx(i), height - B + 23);
  });
  ctx.strokeStyle = t.lineStrong; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(L + 0.5, T); ctx.lineTo(L + 0.5, height - B); ctx.lineTo(width - R, height - B + 0.5); ctx.stroke();
  return { cx, slot, y, groups };
}
