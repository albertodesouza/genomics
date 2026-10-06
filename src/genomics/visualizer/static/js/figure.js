// Figure export: a page's plot area as a PNG (bitmap at a chosen scale) or an SVG (vector).
//
// Canvases are re-rendered for the export: at a higher device-pixel ratio for PNG, or into SvgContext
// (a recording subset of CanvasRenderingContext2D) for SVG, so plots stay vector in the SVG. The HTML
// around them (track labels, legend, swatches, borders) is drawn from its computed style. A title and
// caption lines (locus, tracks, cohort) are added above the figure. Plot code needs no changes:
// setupCanvas() consults ``exportState``.
import { h, icon, toast, checkbox } from './ui.js';

export const exportState = { dpr: null, recorders: null };

const SKIP = 'button, input, select, textarea, .panel-grip, .hint, .loading-line, .overlay, .tooltip, [data-export="skip"]';
const r2 = (v) => Math.round(v * 100) / 100;
const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

/** "rgba(1,2,3,.5)" / "#rrggbbaa" -> [opaque color, alpha] (SVG 1.1 readers ignore rgba()). */
function splitAlpha(color) {
  const s = String(color || '').trim();
  let m = /^rgba?\(\s*([\d.]+)[,\s]+([\d.]+)[,\s]+([\d.]+)(?:\s*[,/]\s*([\d.]+%?))?\s*\)$/i.exec(s);
  if (m) {
    let a = m[4] === undefined ? 1 : m[4].endsWith('%') ? parseFloat(m[4]) / 100 : parseFloat(m[4]);
    return [`rgb(${Math.round(+m[1])},${Math.round(+m[2])},${Math.round(+m[3])})`, a];
  }
  m = /^#([0-9a-f]{6})([0-9a-f]{2})$/i.exec(s);
  if (m) return [`#${m[1]}`, parseInt(m[2], 16) / 255];
  if (s === 'transparent') return ['#000', 0];
  return [s || '#000', 1];
}

/** CSS font shorthand as written by the plot code ("700 12px Inter, sans-serif") -> SVG attributes. */
function fontAttrs(font) {
  const m = /^\s*(italic\s+)?(?:(normal|bold|\d{3})\s+)?([\d.]+)px\s+(.+)$/.exec(font || '');
  if (!m) return `font-size="10"`;
  const weight = m[2] && m[2] !== 'normal' ? ` font-weight="${m[2]}"` : '';
  return `font-family="${esc(m[4].replace(/"/g, "'"))}" font-size="${m[3]}"${weight}${m[1] ? ' font-style="italic"' : ''}`;
}

const STATE_KEYS = ['fillStyle', 'strokeStyle', 'lineWidth', 'font', 'textAlign', 'textBaseline', 'globalAlpha', 'lineJoin', 'lineCap', 'imageSmoothingEnabled'];

/** Records the Canvas 2D calls the visualizer uses and serialises them as SVG elements. */
export class SvgContext {
  constructor(width, height) {
    this.width = width; this.height = height;
    this.canvas = { width, height };
    this.fillStyle = '#000'; this.strokeStyle = '#000'; this.lineWidth = 1; this.font = '10px sans-serif';
    this.textAlign = 'start'; this.textBaseline = 'alphabetic'; this.globalAlpha = 1; this.lineJoin = 'miter'; this.lineCap = 'butt';
    this.imageSmoothingEnabled = true;
    this._m = [1, 0, 0, 1, 0, 0]; this._dash = []; this._clip = null; this._path = []; this._stack = [];
    this._defs = []; this._out = []; this._ids = 0;
    this._measure = document.createElement('canvas').getContext('2d');
  }

  // state
  save() { this._stack.push({ ...Object.fromEntries(STATE_KEYS.map((k) => [k, this[k]])), m: [...this._m], dash: [...this._dash], clip: this._clip }); }
  restore() {
    const s = this._stack.pop();
    if (!s) return;
    for (const k of STATE_KEYS) this[k] = s[k];
    this._m = s.m; this._dash = s.dash; this._clip = s.clip;
  }
  setLineDash(d) { this._dash = [...(d || [])]; }
  getLineDash() { return [...this._dash]; }
  setTransform(a, b, c, d, e, f) { this._m = typeof a === 'object' ? [a.a, a.b, a.c, a.d, a.e, a.f] : [a, b, c, d, e, f]; }
  resetTransform() { this._m = [1, 0, 0, 1, 0, 0]; }
  transform(a, b, c, d, e, f) {
    const [A, B, C, D, E, F] = this._m;
    this._m = [A * a + C * b, B * a + D * b, A * c + C * d, B * c + D * d, A * e + C * f + E, B * e + D * f + F];
  }
  translate(x, y) { this.transform(1, 0, 0, 1, x, y); }
  scale(x, y) { this.transform(x, 0, 0, y, 0, 0); }
  rotate(a) { const c = Math.cos(a); const s = Math.sin(a); this.transform(c, s, -s, c, 0, 0); }
  measureText(text) { this._measure.font = this.font; return this._measure.measureText(text); }

  // paths (points are stored in device space, like canvas does)
  _pt(x, y) { const [a, b, c, d, e, f] = this._m; return [a * x + c * y + e, b * x + d * y + f]; }
  beginPath() { this._path = []; }
  moveTo(x, y) { const [px, py] = this._pt(x, y); this._path.push(`M${r2(px)} ${r2(py)}`); }
  lineTo(x, y) { const [px, py] = this._pt(x, y); this._path.push(`L${r2(px)} ${r2(py)}`); }
  closePath() { this._path.push('Z'); }
  rect(x, y, w, h) { this.moveTo(x, y); this.lineTo(x + w, y); this.lineTo(x + w, y + h); this.lineTo(x, y + h); this.closePath(); }
  arc(x, y, r, a0, a1, ccw = false) {
    let sweep = ccw ? a0 - a1 : a1 - a0;
    if (sweep < 0) sweep += 2 * Math.PI;
    const n = Math.max(4, Math.ceil(Math.min(sweep, 2 * Math.PI) / (Math.PI / 16)));
    for (let i = 0; i <= n; i++) {
      const a = a0 + (ccw ? -1 : 1) * Math.min(sweep, 2 * Math.PI) * (i / n);
      if (i === 0 && !this._path.length) this.moveTo(x + r * Math.cos(a), y + r * Math.sin(a));
      else this.lineTo(x + r * Math.cos(a), y + r * Math.sin(a));
    }
  }
  clip() {
    if (!this._path.length) return;
    const id = `c${++this._ids}`;
    this._defs.push(`<clipPath id="${id}"${this._clip ? ` clip-path="url(#${this._clip})"` : ''}><path d="${this._path.join('')}"/></clipPath>`);
    this._clip = id;
  }

  // painting
  _paint(kind) {
    const [color, alpha] = splitAlpha(kind === 'fill' ? this.fillStyle : this.strokeStyle);
    const a = alpha * this.globalAlpha;
    if (a <= 0) return '';
    let s = kind === 'fill' ? `fill="${color}"` : `fill="none" stroke="${color}" stroke-width="${r2(this.lineWidth * this._scale())}"`;
    if (a < 1) s += ` ${kind}-opacity="${r2(a)}"`;
    if (kind === 'stroke') {
      if (this._dash.length) s += ` stroke-dasharray="${this._dash.map((v) => r2(v * this._scale())).join(' ')}"`;
      if (this.lineJoin !== 'miter') s += ` stroke-linejoin="${this.lineJoin}"`;
      if (this.lineCap !== 'butt') s += ` stroke-linecap="${this.lineCap}"`;
    }
    if (this._clip) s += ` clip-path="url(#${this._clip})"`;
    return s;
  }
  _scale() { const [a, b, c, d] = this._m; return Math.sqrt(Math.abs(a * d - b * c)) || 1; }
  _emitPath(d, kind) { const p = this._paint(kind); if (p && d) this._out.push(`<path d="${d}" ${p}/>`); }
  fill() { this._emitPath(this._path.join(''), 'fill'); }
  stroke() { this._emitPath(this._path.join(''), 'stroke'); }
  _rectPath(x, y, w, h) {
    const saved = this._path;
    this._path = []; this.rect(x, y, w, h);
    const d = this._path.join('');
    this._path = saved;
    return d;
  }
  fillRect(x, y, w, h) {
    const [a, b, c, d, e, f] = this._m;
    if (b === 0 && c === 0) {
      const p = this._paint('fill');
      if (!p) return;
      let X = a * x + e; let W = a * w; let Y = d * y + f; let H = d * h;
      if (W < 0) { X += W; W = -W; }
      if (H < 0) { Y += H; H = -H; }
      this._out.push(`<rect x="${r2(X)}" y="${r2(Y)}" width="${r2(W)}" height="${r2(H)}" ${p}/>`);
    } else this._emitPath(this._rectPath(x, y, w, h), 'fill');
  }
  strokeRect(x, y, w, h) { this._emitPath(this._rectPath(x, y, w, h), 'stroke'); }
  clearRect() { /* the export starts from an empty group */ }
  fillText(text, x, y) {
    const p = this._paint('fill');
    if (!p || text === '' || text == null) return;
    const anchor = { center: 'middle', right: 'end', end: 'end' }[this.textAlign] || 'start';
    const base = { top: 'text-before-edge', hanging: 'hanging', middle: 'central', bottom: 'text-after-edge', ideographic: 'ideographic' }[this.textBaseline];
    const [a, b, c, d, e, f] = this._m;
    const tf = a === 1 && b === 0 && c === 0 && d === 1 ? '' : ` transform="matrix(${[a, b, c, d, e, f].map(r2).join(' ')})"`;
    const [X, Y] = tf ? [x, y] : [x + e, y + f];
    this._out.push(`<text x="${r2(X)}" y="${r2(Y)}"${tf} ${fontAttrs(this.font)}${anchor !== 'start' ? ` text-anchor="${anchor}"` : ''}${base ? ` dominant-baseline="${base}"` : ''} ${p} xml:space="preserve">${esc(text)}</text>`);
  }
  strokeText() { /* not used by the plots */ }
  drawImage(img, ...args) {
    let sx = 0; let sy = 0; let sw = img.width; let sh = img.height; let dx; let dy; let dw = img.width; let dh = img.height;
    if (args.length === 2) [dx, dy] = args;
    else if (args.length === 4) [dx, dy, dw, dh] = args;
    else [sx, sy, sw, sh, dx, dy, dw, dh] = args;
    let src = img;
    if (sx || sy || sw !== img.width || sh !== img.height) {
      src = document.createElement('canvas');
      src.width = Math.max(1, sw); src.height = Math.max(1, sh);
      src.getContext('2d').drawImage(img, sx, sy, sw, sh, 0, 0, sw, sh);
    }
    const href = src.toDataURL ? src.toDataURL('image/png') : '';
    if (!href) return;
    const [a, b, c, d, e, f] = this._m;
    const pixel = this.imageSmoothingEnabled ? '' : ' style="image-rendering:pixelated" image-rendering="optimizeSpeed"';
    const alpha = this.globalAlpha < 1 ? ` opacity="${r2(this.globalAlpha)}"` : '';
    this._out.push(`<image href="${href}" x="${r2(dx)}" y="${r2(dy)}" width="${r2(dw)}" height="${r2(dh)}" preserveAspectRatio="none" transform="matrix(${[a, b, c, d, e, f].map(r2).join(' ')})"${pixel}${alpha}${this._clip ? ` clip-path="url(#${this._clip})"` : ''}/>`);
  }
  createImageData(w, h) { return new ImageData(Math.max(1, w), Math.max(1, h)); }
  putImageData() { /* only used on real offscreen canvases */ }

  /** Defs and elements, ids prefixed so several recordings can share one SVG document. */
  toSvg(prefix) {
    const fix = (s) => s.replace(/(id="|url\(#)c(\d+)/g, `$1${prefix}c$2`);
    return { defs: this._defs.map(fix).join(''), body: this._out.map(fix).join('') };
  }
}

// ----------------------------------------------------------------------------- DOM walk
function opacityOf(el, root) {
  let o = 1;
  for (let n = el; n && n !== root.parentElement; n = n.parentElement) o *= parseFloat(getComputedStyle(n).opacity || '1');
  return o;
}

function isAbsolute(el, root) {
  for (let n = el; n && n !== root; n = n.parentElement) {
    const p = getComputedStyle(n).position;
    if (p === 'absolute' || p === 'fixed') return true;
  }
  return false;
}

/** Visible elements under ``root`` (document order), with what is needed to draw them. */
function collect(root, ox, oy) {
  const out = [];
  for (const el of root.querySelectorAll('*')) {
    if (el.closest(SKIP) || el.closest('[hidden]')) continue;
    const cs = getComputedStyle(el);
    if (cs.display === 'none' || cs.visibility === 'hidden') continue;
    const r = el.getBoundingClientRect();
    if (r.width <= 0 || r.height <= 0) continue;
    const alpha = opacityOf(el, root);
    if (alpha <= 0.01) continue;
    out.push({ el, cs, alpha, abs: isAbsolute(el, root), x: r.left - ox, y: r.top - oy, w: r.width, h: r.height });
  }
  return out;
}

function boxes(item) {
  const { cs, x, y, w, h, alpha } = item;
  const shapes = [];
  const [bg, ba] = splitAlpha(cs.backgroundColor);
  const radius = Math.min(parseFloat(cs.borderTopLeftRadius) || 0, h / 2, w / 2);
  if (ba > 0 && !item.el.matches('canvas')) shapes.push({ kind: 'rect', x, y, w, h, radius, color: bg, alpha: ba * alpha });
  for (const side of ['Top', 'Right', 'Bottom', 'Left']) {
    const bw = parseFloat(cs[`border${side}Width`]) || 0;
    const style = cs[`border${side}Style`];
    if (!bw || style === 'none' || style === 'hidden') continue;
    const [color, a] = splitAlpha(cs[`border${side}Color`]);
    if (a <= 0) continue;
    const half = bw / 2;
    const seg = { Top: [x, y + half, x + w, y + half], Bottom: [x, y + h - half, x + w, y + h - half], Left: [x + half, y, x + half, y + h], Right: [x + w - half, y, x + w - half, y + h] }[side];
    shapes.push({ kind: 'line', seg, width: bw, color, alpha: a * alpha, dash: style === 'dashed' ? [bw * 3, bw * 2] : style === 'dotted' ? [bw, bw] : [] });
  }
  return shapes;
}

/** Text runs of an element's own text nodes, one per rendered line, positioned and clipped like the page. */
function texts(item, ox, oy, measure) {
  const { el, cs, alpha } = item;
  const runs = [];
  const clipped = cs.overflow === 'hidden' || cs.overflowX === 'hidden';
  const font = `${cs.fontStyle === 'italic' ? 'italic ' : ''}${cs.fontWeight} ${cs.fontSize} ${cs.fontFamily}`;
  const [color, a] = splitAlpha(cs.color);
  const upper = cs.textTransform === 'uppercase';
  measure.font = font;
  const ascent = measure.measureText('Mg').fontBoundingBoxAscent || parseFloat(cs.fontSize) * 0.8;
  const push = (text, rect) => {
    let t = text.replace(/\s+/g, ' ');
    if (upper) t = t.toUpperCase();
    const trimmed = t.trim();
    if (!trimmed) return;
    // A leading space is part of the rect only when the layout rendered it (not collapsed).
    let lead = 0;
    if (/^\s/.test(t)) {
      const withSpace = measure.measureText(` ${trimmed}`).width;
      const without = measure.measureText(trimmed).width;
      if (Math.abs(rect.width - withSpace) < Math.abs(rect.width - without)) lead = withSpace - without;
    }
    runs.push({ text: trimmed, font, color, alpha: a * alpha, x: rect.left - ox + lead, y: rect.top - oy + ascent, clip: clipped ? { x: item.x, y: item.y, w: item.w, h: item.h } : null });
  };
  for (const node of el.childNodes) {
    if (node.nodeType !== 3 || !node.textContent.trim()) continue;
    const range = document.createRange();
    range.selectNodeContents(node);
    const rects = [...range.getClientRects()].filter((r) => r.width > 0);
    if (!rects.length) continue;
    if (rects.length === 1) { push(node.textContent, rects[0]); continue; }
    // Wrapped text: group the words by the line they were laid out on.
    const content = node.textContent;
    // Browsers also break lines after hyphens ("RNA-" / "seq"), so tokens end at a hyphen too.
    const words = [...content.matchAll(/\s*(?:[^\s-]+-?|-)/g)];
    let line = null;
    for (const w of words) {
      range.setStart(node, w.index); range.setEnd(node, w.index + w[0].length);
      const wr = [...range.getClientRects()].filter((r) => r.width > 0);
      const r = wr[wr.length - 1];
      if (!r) continue;
      const top = Math.round(r.top);
      if (line && Math.abs(line.top - top) <= 2) { line.end = w.index + w[0].length; line.right = r.right; continue; }
      if (line) push(content.slice(line.start, line.end), { left: line.left, top: line.top, width: line.right - line.left });
      // The word's first rect may hold only its leading space at the end of the previous line.
      const glyphs = content.slice(w.index).match(/^\s*/)[0].length;
      range.setStart(node, w.index + glyphs);
      const gr = [...range.getClientRects()].filter((x) => x.width > 0);
      const g = gr[0] || r;
      line = { start: w.index + glyphs, end: w.index + w[0].length, top: Math.round(g.top), left: g.left, right: r.right };
    }
    if (line) push(content.slice(line.start, line.end), { left: line.left, top: line.top, width: line.right - line.left });
  }
  return runs;
}

// ----------------------------------------------------------------------------- composition
function header(title, caption, width) {
  const lines = [];
  if (title) lines.push({ text: title, size: 15, weight: 600 });
  for (const c of caption || []) if (c) lines.push({ text: c, size: 12, weight: 400 });
  const pad = 16;
  let y = pad;
  const placed = lines.map((l) => { y += l.size + 6; return { ...l, y: y - 6 }; });
  return { lines: placed, height: lines.length ? y + pad - 6 : 0, width };
}

/** Run ``fn`` with the light theme (when ``light``) and with [data-export="skip"] parts collapsed. */
function exportLayout(root, light, fn) {
  const html = document.documentElement;
  const prev = html.getAttribute('data-theme');
  const skipped = [...root.querySelectorAll('[data-export="skip"]')].map((el) => [el, el.style.display]);
  if (light && prev !== 'light') html.setAttribute('data-theme', 'light');
  for (const [el] of skipped) el.style.display = 'none';
  try { return fn(); } finally {
    for (const [el, display] of skipped) el.style.display = display;
    if (light && prev !== 'light') { if (prev === null) html.removeAttribute('data-theme'); else html.setAttribute('data-theme', prev); }
  }
}

function renderPng({ root, redraw, title, caption, scale }) {
  const ox0 = root.getBoundingClientRect();
  const ox = ox0.left; const oy = ox0.top - root.scrollTop;
  const width = root.clientWidth; const height = root.scrollHeight;
  exportState.dpr = scale;
  try { redraw(); } finally { exportState.dpr = null; }
  const items = collect(root, ox, oy);
  const hd = header(title, caption, width);
  const out = document.createElement('canvas');
  out.width = Math.round(width * scale); out.height = Math.round((height + hd.height) * scale);
  const ctx = out.getContext('2d');
  const t = getComputedStyle(document.documentElement);
  ctx.scale(scale, scale);
  ctx.fillStyle = t.getPropertyValue('--surface').trim() || '#fff';
  ctx.fillRect(0, 0, width, height + hd.height);
  const fontFamily = t.getPropertyValue('--font').trim() || 'sans-serif';
  for (const l of hd.lines) {
    ctx.font = `${l.weight} ${l.size}px ${fontFamily}`;
    ctx.fillStyle = l.weight > 400 ? t.getPropertyValue('--ink').trim() : t.getPropertyValue('--ink-2').trim();
    ctx.fillText(l.text, 16, l.y);
  }
  ctx.translate(0, hd.height);
  const paintBoxes = (item) => {
    for (const s of boxes(item)) {
      ctx.save();
      ctx.globalAlpha = s.alpha;
      if (s.kind === 'rect') {
        ctx.fillStyle = s.color;
        ctx.beginPath();
        if (s.radius && ctx.roundRect) ctx.roundRect(s.x, s.y, s.w, s.h, s.radius); else ctx.rect(s.x, s.y, s.w, s.h);
        ctx.fill();
      } else {
        ctx.strokeStyle = s.color; ctx.lineWidth = s.width; ctx.setLineDash(s.dash);
        ctx.beginPath(); ctx.moveTo(s.seg[0], s.seg[1]); ctx.lineTo(s.seg[2], s.seg[3]); ctx.stroke();
      }
      ctx.restore();
    }
  };
  for (const it of items) if (!it.abs) paintBoxes(it);
  for (const it of items) {
    if (!(it.el instanceof HTMLCanvasElement) || !it.el.width) continue;
    ctx.save(); ctx.globalAlpha = it.alpha; ctx.drawImage(it.el, it.x, it.y, it.w, it.h); ctx.restore();
  }
  for (const it of items) if (it.abs) paintBoxes(it);
  for (const it of items) {
    for (const run of texts(it, ox, oy, ctx)) {
      ctx.save();
      if (run.clip) { ctx.beginPath(); ctx.rect(run.clip.x, run.clip.y, run.clip.w, run.clip.h); ctx.clip(); }
      ctx.globalAlpha = run.alpha; ctx.fillStyle = run.color; ctx.font = run.font; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
      ctx.fillText(run.text, run.x, run.y);
      ctx.restore();
    }
  }
  return out;
}

function renderSvg({ root, redraw, title, caption }) {
  const r0 = root.getBoundingClientRect();
  const ox = r0.left; const oy = r0.top - root.scrollTop;
  const width = root.clientWidth; const height = root.scrollHeight;
  const recorders = new Map();
  exportState.recorders = recorders;
  try { redraw(); } finally { exportState.recorders = null; }
  const items = collect(root, ox, oy);
  const hd = header(title, caption, width);
  const t = getComputedStyle(document.documentElement);
  const fontFamily = (t.getPropertyValue('--font').trim() || 'sans-serif').replace(/"/g, "'");
  const defs = [];
  const body = [];
  let clipN = 0;
  const op = (a) => (a < 1 ? ` opacity="${r2(a)}"` : '');
  const svgBoxes = (item) => {
    for (const s of boxes(item)) {
      if (s.kind === 'rect') body.push(`<rect x="${r2(s.x)}" y="${r2(s.y)}" width="${r2(s.w)}" height="${r2(s.h)}"${s.radius ? ` rx="${r2(s.radius)}"` : ''} fill="${s.color}"${op(s.alpha)}/>`);
      else body.push(`<line x1="${r2(s.seg[0])}" y1="${r2(s.seg[1])}" x2="${r2(s.seg[2])}" y2="${r2(s.seg[3])}" stroke="${s.color}" stroke-width="${r2(s.width)}"${s.dash.length ? ` stroke-dasharray="${s.dash.join(' ')}"` : ''}${op(s.alpha)}/>`);
    }
  };
  for (const it of items) if (!it.abs) svgBoxes(it);
  let k = 0;
  for (const it of items) {
    if (!(it.el instanceof HTMLCanvasElement)) continue;
    const rec = recorders.get(it.el);
    const id = `fig${++k}`;
    if (rec) {
      const { defs: d, body: b } = rec.toSvg(`${id}_`);
      defs.push(d, `<clipPath id="${id}"><rect width="${r2(it.w)}" height="${r2(it.h)}"/></clipPath>`);
      body.push(`<g transform="translate(${r2(it.x)} ${r2(it.y)})" clip-path="url(#${id})"${op(it.alpha)}>${b}</g>`);
    } else if (it.el.width) {
      // A canvas the redraw did not touch: embed its pixels.
      body.push(`<image href="${it.el.toDataURL('image/png')}" x="${r2(it.x)}" y="${r2(it.y)}" width="${r2(it.w)}" height="${r2(it.h)}" preserveAspectRatio="none"${op(it.alpha)}/>`);
    }
  }
  for (const it of items) if (it.abs) svgBoxes(it);
  const measure = document.createElement('canvas').getContext('2d');
  for (const it of items) {
    for (const run of texts(it, ox, oy, measure)) {
      let clip = '';
      if (run.clip) {
        const id = `t${++clipN}`;
        defs.push(`<clipPath id="${id}"><rect x="${r2(run.clip.x)}" y="${r2(run.clip.y)}" width="${r2(run.clip.w)}" height="${r2(run.clip.h)}"/></clipPath>`);
        clip = ` clip-path="url(#${id})"`;
      }
      body.push(`<text x="${r2(run.x)}" y="${r2(run.y)}" ${fontAttrs(run.font)} fill="${run.color}"${op(run.alpha)}${clip} xml:space="preserve">${esc(run.text)}</text>`);
    }
  }
  const surface = t.getPropertyValue('--surface').trim() || '#fff';
  const head = hd.lines.map((l) => `<text x="16" y="${l.y}" font-family="${esc(fontFamily)}" font-size="${l.size}" font-weight="${l.weight}" fill="${l.weight > 400 ? t.getPropertyValue('--ink').trim() : t.getPropertyValue('--ink-2').trim()}">${esc(l.text)}</text>`).join('');
  const H = height + hd.height;
  return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${r2(H)}" viewBox="0 0 ${width} ${r2(H)}">`
    + `<defs>${defs.join('')}</defs><rect width="100%" height="100%" fill="${surface}"/>${head}<g transform="translate(0 ${hd.height})">${body.join('')}</g></svg>\n`;
}

function save(filename, blob) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/**
 * Export ``root`` (the scrollable plot area) as ``format`` ('png' | 'svg').
 * ``redraw`` re-renders the page's canvases synchronously from data already loaded.
 * ``light`` draws with the light theme whatever the page uses.
 */
export async function exportFigure({ root, redraw, title = '', caption = [], filename = 'figure', format = 'png', scale = 2, light = true }) {
  if (format === 'svg') {
    const svg = exportLayout(root, light, () => renderSvg({ root, redraw, title, caption }));
    redraw();
    save(`${filename}.svg`, new Blob([svg], { type: 'image/svg+xml' }));
    return;
  }
  const canvas = exportLayout(root, light, () => renderPng({ root, redraw, title, caption, scale }));
  redraw();
  const blob = await new Promise((resolve) => canvas.toBlob(resolve, 'image/png'));
  if (!blob) throw new Error('The browser could not encode the PNG (figure too large?)');
  save(`${filename}.png`, blob);
}

/** "OCA2_chr15-27719008-28243295" style file name part. */
export function slug(text) { return String(text).replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 120) || 'figure'; }

// ----------------------------------------------------------------------------- toolbar control

const LIGHT_KEY = 'gv.figureLight';

/**
 * An "Export" button with a small menu (PNG 2x / 4x, SVG). ``getOptions()`` returns the
 * exportFigure options except ``format``/``scale``/``light`` (root, redraw, title, caption, filename),
 * or null when there is nothing to export yet.
 */
/** Copy text to the clipboard, falling back to a prompt where the API is unavailable (http://). */
export async function copyText(text, what = 'Copied') {
  try {
    await navigator.clipboard.writeText(text);
    toast(`${what} to the clipboard`);
  } catch (e) {
    window.prompt('Copy this:', text);
  }
}

export function exportMenu(getOptions, { label = 'Export', python = null } = {}) {
  let light = true;
  try { light = localStorage.getItem(LIGHT_KEY) !== '0'; } catch (e) { /* ignore */ }
  const pop = h('div', { class: 'menu-pop', role: 'menu', hidden: true });
  const close = () => { pop.hidden = true; document.removeEventListener('pointerdown', outside, true); };
  const outside = (e) => { if (!wrap.contains(e.target)) close(); };
  const run = async (format, scale) => {
    close();
    const opts = getOptions();
    if (!opts) { toast('Nothing to export yet', 'info'); return; }
    try {
      await exportFigure({ ...opts, format, scale, light });
    } catch (err) {
      toast(`Export failed: ${err.message}`, 'error', 8000);
    }
  };
  const item = (text, sub, fn) => h('button', { class: 'menu-item', type: 'button', role: 'menuitem', onclick: fn }, h('span', null, text), h('span', { class: 'muted' }, sub));
  pop.append(
    item('PNG', '2× pixels', () => run('png', 2)),
    item('PNG, print', '4× pixels', () => run('png', 4)),
    item('SVG', 'vector, editable', () => run('svg', 1)),
    python ? h('div', { class: 'menu-sep' }) : null,
    python ? item('Copy as Python', 'the same arrays in a notebook', () => {
      close();
      const code = python();
      if (!code) { toast('Nothing to copy yet', 'info'); return; }
      copyText(code, 'Python snippet copied');
    }) : null,
    h('div', { class: 'menu-sep' }),
    h('div', { class: 'menu-row' }, checkbox('Light background', light, (v) => {
      light = v;
      try { localStorage.setItem(LIGHT_KEY, v ? '1' : '0'); } catch (e) { /* ignore */ }
    })),
  );
  const btn = h('button', { class: 'btn small', type: 'button', title: 'Download this view as a figure (PNG or SVG) with its locus, tracks and cohort as caption', 'aria-haspopup': 'menu', onclick: () => {
    if (!pop.hidden) { close(); return; }
    pop.hidden = false;
    document.addEventListener('pointerdown', outside, true);
  } }, icon('download', 14), label);
  const wrap = h('span', { class: 'menu-wrap' }, btn, pop);
  pop.addEventListener('keydown', (e) => { if (e.key === 'Escape') { close(); btn.focus(); } });
  return wrap;
}
