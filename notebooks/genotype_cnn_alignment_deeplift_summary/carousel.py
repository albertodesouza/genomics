"""Self-contained, no-kernel-required Prev/Next figure browsers for notebook output."""

import base64
import io
import json
import uuid

from IPython.display import HTML, display


def render_carousel(figures, labels=None, dpi=140):
    """Browse a list of already-built matplotlib Figures one at a time via Prev/Next buttons.

    Renders each Figure to a base64 PNG data URI up front and emits a small self-contained
    HTML/JS widget (plain `<img>` tags toggled by inline JavaScript) instead of an
    `ipywidgets.Image`. ipywidgets' Prev/Next is Python-callback-driven (`Button.on_click` runs a
    handler in the kernel via a comm channel), so it goes dead the moment there's no live kernel
    behind the page -- e.g. after `jupyter nbconvert --to html`, or in a notebook viewer that
    doesn't spin up a kernel. Plain client-side JS has no such dependency: every frame is already
    embedded in the page, so Prev/Next keeps working in the exported static HTML exactly as it
    does live in Jupyter. Figures must have been created with `plt.close(fig)` called right after
    building them.
    """
    if not figures:
        print("(no figures to display)")
        return
    labels = list(labels) if labels is not None else [str(i) for i in range(len(figures))]
    assert len(labels) == len(figures)

    def _to_data_uri(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")

    uid = uuid.uuid4().hex[:8]
    imgs_html = "\n".join(
        f'<img data-idx="{i}" src="{_to_data_uri(fig)}" '
        f'style="display:{"block" if i == 0 else "none"}; max-width:100%;">'
        for i, fig in enumerate(figures)
    )
    labels_json = json.dumps(labels)

    html = f'''
<div id="carousel-{uid}" style="font-family: sans-serif;">
  <script type="application/json" class="carousel-labels">{labels_json}</script>
  <div style="display:flex; align-items:center; gap:8px; margin-bottom:6px;">
    <button onclick="__carouselStep(\'carousel-{uid}\', -1)">&lt; Prev</button>
    <span class="carousel-label"></span>
    <button onclick="__carouselStep(\'carousel-{uid}\', 1)">Next &gt;</button>
  </div>
  <div class="carousel-images">
{imgs_html}
  </div>
</div>
<script>
(function() {{
  if (typeof window.__carouselStep !== "function") {{
    window.__carouselStep = function(containerId, delta) {{
      var c = document.getElementById(containerId);
      var imgs = c.querySelectorAll(".carousel-images img");
      var labels = JSON.parse(c.querySelector(".carousel-labels").textContent);
      var idx = ((parseInt(c.dataset.idx || "0") + delta) % imgs.length + imgs.length) % imgs.length;
      c.dataset.idx = idx;
      imgs.forEach(function(im, i) {{ im.style.display = (i === idx) ? "block" : "none"; }});
      c.querySelector(".carousel-label").textContent = labels[idx] + " (" + (idx + 1) + "/" + imgs.length + ")";
    }};
  }}
  __carouselStep("carousel-{uid}", 0);
}})();
</script>
'''
    display(HTML(html))


def render_carousel_with_toggle(figures_by_view, labels=None, view_labels=None, dpi=140):
    """Browse a list of figures (one per label) via Prev/Next, plus a toggle button that swaps
    between different renderings of the same content per label -- e.g. log-odds vs. probability
    x-axis. `figures_by_view` is a dict {view_key: [fig_per_label, ...]}; every list must have the
    same length and label order. `view_labels` optionally maps view_key -> button text (defaults
    to the key itself). Same no-kernel-required, everything-pre-rendered rationale as
    `render_carousel` -- Prev/Next and the view toggle are both plain client-side JS. Figures must
    have been created with `plt.close(fig)` called right after building them.
    """
    if not figures_by_view:
        print("(no figures to display)")
        return
    view_keys = list(figures_by_view)
    n = len(figures_by_view[view_keys[0]])
    for k in view_keys:
        assert len(figures_by_view[k]) == n, f"view {k!r} has {len(figures_by_view[k])} figures, expected {n}"
    labels = list(labels) if labels is not None else [str(i) for i in range(n)]
    assert len(labels) == n
    view_labels = view_labels or {k: str(k) for k in view_keys}

    def _to_data_uri(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")

    uid = uuid.uuid4().hex[:8]
    imgs_html_parts = []
    for vi, vkey in enumerate(view_keys):
        for i, fig in enumerate(figures_by_view[vkey]):
            show = (i == 0 and vi == 0)
            imgs_html_parts.append(
                f'<img data-idx="{i}" data-view="{vi}" src="{_to_data_uri(fig)}" '
                f'style="display:{"block" if show else "none"}; max-width:100%;">'
            )
    imgs_html = "\n".join(imgs_html_parts)
    labels_json = json.dumps(labels)
    nav_display = "flex" if n > 1 else "none"
    view_buttons_html = "\n".join(
        f'<button data-view="{vi}" onclick="__carouselToggleSetView(\'carousel-toggle-{uid}\', {vi})">'
        f'{view_labels[vkey]}</button>'
        for vi, vkey in enumerate(view_keys)
    )

    html = f'''
<div id="carousel-toggle-{uid}" data-idx="0" data-view="0" style="font-family: sans-serif;">
  <script type="application/json" class="carousel-labels">{labels_json}</script>
  <div style="display:{nav_display}; align-items:center; gap:8px; margin-bottom:6px;">
    <button onclick="__carouselToggleStep(\'carousel-toggle-{uid}\', -1)">&lt; Prev</button>
    <span class="carousel-label"></span>
    <button onclick="__carouselToggleStep(\'carousel-toggle-{uid}\', 1)">Next &gt;</button>
  </div>
  <div style="display:flex; align-items:center; gap:8px; margin-bottom:6px;">
    <span style="font-weight:600;">x-axis:</span>
{view_buttons_html}
  </div>
  <div class="carousel-images">
{imgs_html}
  </div>
</div>
<script>
(function() {{
  if (typeof window.__carouselToggleRender !== "function") {{
    window.__carouselToggleRender = function(containerId) {{
      var c = document.getElementById(containerId);
      var imgs = c.querySelectorAll(".carousel-images img");
      var labels = JSON.parse(c.querySelector(".carousel-labels").textContent);
      var idx = parseInt(c.dataset.idx);
      var view = parseInt(c.dataset.view);
      imgs.forEach(function(im) {{
        var show = (parseInt(im.dataset.idx) === idx && parseInt(im.dataset.view) === view);
        im.style.display = show ? "block" : "none";
      }});
      var labelEl = c.querySelector(".carousel-label");
      if (labelEl) {{ labelEl.textContent = labels[idx] + " (" + (idx + 1) + "/" + labels.length + ")"; }}
      c.querySelectorAll("button[data-view]").forEach(function(btn) {{
        var active = parseInt(btn.dataset.view) === view;
        btn.style.fontWeight = active ? "bold" : "normal";
        btn.style.textDecoration = active ? "underline" : "none";
      }});
    }};
    window.__carouselToggleStep = function(containerId, delta) {{
      var c = document.getElementById(containerId);
      var imgs = c.querySelectorAll(".carousel-images img");
      var n = new Set(Array.from(imgs).map(function(im) {{ return im.dataset.idx; }})).size;
      var idx = ((parseInt(c.dataset.idx) + delta) % n + n) % n;
      c.dataset.idx = idx;
      window.__carouselToggleRender(containerId);
    }};
    window.__carouselToggleSetView = function(containerId, view) {{
      var c = document.getElementById(containerId);
      c.dataset.view = view;
      window.__carouselToggleRender(containerId);
    }};
  }}
  window.__carouselToggleRender("carousel-toggle-{uid}");
}})();
</script>
'''
    display(HTML(html))


def render_carousel_2d(figures, row_labels, col_labels, row_name="gene", col_name="technique", dpi=140):
    """Browse a 2D grid of already-built matplotlib Figures via two independent Prev/Next button
    pairs -- one steps through `row_labels`, the other through `col_labels` -- both updating the
    same displayed image. Pure HTML/JS, same no-kernel-required rationale as `render_carousel`
    above (every (row, col) frame is embedded in the page up front; the buttons just toggle which
    one is visible).

    `figures` is a dict keyed by `(row_label, col_label)` tuples; every combination must be
    present. Figures must have been created with `plt.close(fig)` called right after building
    them.
    """
    if not figures:
        print("(no figures to display)")
        return
    row_labels, col_labels = list(row_labels), list(col_labels)

    def _to_data_uri(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")

    uid = uuid.uuid4().hex[:8]
    imgs_html_parts = []
    for ri, row_label in enumerate(row_labels):
        for ci, col_label in enumerate(col_labels):
            fig = figures[(row_label, col_label)]
            show = (ri == 0 and ci == 0)
            imgs_html_parts.append(
                f'<img data-row="{ri}" data-col="{ci}" src="{_to_data_uri(fig)}" '
                f'style="display:{"block" if show else "none"}; max-width:100%;">'
            )
    imgs_html = "\n".join(imgs_html_parts)
    row_labels_json = json.dumps(row_labels)
    col_labels_json = json.dumps(col_labels)

    html = f'''
<div id="carousel2d-{uid}" data-rowidx="0" data-colidx="0" data-nrows="{len(row_labels)}" data-ncols="{len(col_labels)}" style="font-family: sans-serif;">
  <script type="application/json" class="carousel2d-row-labels">{row_labels_json}</script>
  <script type="application/json" class="carousel2d-col-labels">{col_labels_json}</script>
  <div style="display:flex; align-items:center; gap:8px; margin-bottom:4px;">
    <button onclick="__carousel2dStep(\'carousel2d-{uid}\', \'row\', -1)">&lt; Prev {row_name}</button>
    <span class="carousel2d-row-label"></span>
    <button onclick="__carousel2dStep(\'carousel2d-{uid}\', \'row\', 1)">Next {row_name} &gt;</button>
  </div>
  <div style="display:flex; align-items:center; gap:8px; margin-bottom:6px;">
    <button onclick="__carousel2dStep(\'carousel2d-{uid}\', \'col\', -1)">&lt; Prev {col_name}</button>
    <span class="carousel2d-col-label"></span>
    <button onclick="__carousel2dStep(\'carousel2d-{uid}\', \'col\', 1)">Next {col_name} &gt;</button>
  </div>
  <div class="carousel2d-images">
{imgs_html}
  </div>
</div>
<script>
(function() {{
  if (typeof window.__carousel2dStep !== "function") {{
    window.__carousel2dStep = function(containerId, axis, delta) {{
      var c = document.getElementById(containerId);
      var imgs = c.querySelectorAll(".carousel2d-images img");
      var nRows = parseInt(c.dataset.nrows);
      var nCols = parseInt(c.dataset.ncols);
      var rowIdx = parseInt(c.dataset.rowidx);
      var colIdx = parseInt(c.dataset.colidx);
      if (axis === "row") {{ rowIdx = ((rowIdx + delta) % nRows + nRows) % nRows; }}
      else {{ colIdx = ((colIdx + delta) % nCols + nCols) % nCols; }}
      c.dataset.rowidx = rowIdx;
      c.dataset.colidx = colIdx;
      imgs.forEach(function(im) {{
        var show = (parseInt(im.dataset.row) === rowIdx && parseInt(im.dataset.col) === colIdx);
        im.style.display = show ? "block" : "none";
      }});
      var rowLabels = JSON.parse(c.querySelector(".carousel2d-row-labels").textContent);
      var colLabels = JSON.parse(c.querySelector(".carousel2d-col-labels").textContent);
      c.querySelector(".carousel2d-row-label").textContent = rowLabels[rowIdx] + " (" + (rowIdx + 1) + "/" + nRows + ")";
      c.querySelector(".carousel2d-col-label").textContent = colLabels[colIdx] + " (" + (colIdx + 1) + "/" + nCols + ")";
    }};
  }}
  __carousel2dStep("carousel2d-{uid}", "row", 0);
}})();
</script>
'''
    display(HTML(html))
