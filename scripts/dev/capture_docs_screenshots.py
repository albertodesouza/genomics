#!/usr/bin/env python3
"""Capture the visualizer screenshots used by the documentation (docs/assets/visualizer/).

Drives a *running* visualizer in headless Chromium, sets each view up (pinned samples, cohort,
locus, tracks), optionally draws numbered callouts over the page, and writes one WebP per shot.

    genomics visualize --port 8796 --dataset-id 1kg_high_coverage &
    python3 scripts/dev/capture_docs_screenshots.py --url http://127.0.0.1:8796/
    python3 scripts/dev/capture_docs_screenshots.py --url http://127.0.0.1:8796/ --only tracks-individuals,variant

The shots follow the tutorial's running example on the 1000 Genomes dataset (`1kg_high_coverage`):
HG00096 and the pigmentation genes. Some views need cached cohort aggregates, AlphaGenome calls
(Perturbation Lab, expression in other tissues) or public databases (GTEx, Ensembl, ENCODE); the
first capture can take several minutes, later ones reuse the visualizer's cache.

Needs Playwright with Chromium (`pip install -e ".[test-ui]" && python3 -m playwright install chromium`)
and Pillow.
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import time
from pathlib import Path
from urllib.parse import quote

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "docs" / "assets" / "visualizer"
DATASET = "1kg_high_coverage"
SAMPLE = "HG00096"
PINNED = ["HG00096", "HG00097", "HG00099"]
VIEWPORT = {"width": 1360, "height": 860}

CALLOUT_CSS = """
.doc-callout{position:fixed;z-index:99999;width:26px;height:26px;border-radius:50%;background:#e8590c;color:#fff;
 font:700 14px/26px system-ui,sans-serif;text-align:center;box-shadow:0 0 0 3px #fff,0 2px 8px rgba(0,0,0,.35);pointer-events:none}
.doc-ring{position:fixed;z-index:99998;border:2.5px solid #e8590c;border-radius:8px;pointer-events:none}
"""


class Shooter:
    def __init__(self, page, base: str, out: Path, scale: float, quality: int):
        self.pg = page
        self.base = base.rstrip("/") + "/"
        self.out = out
        self.scale = scale
        self.quality = quality

    # ------------------------------------------------------------------ setup
    def prime(self, pinned=None, filters=None, locus=None, extra=None):
        """Per-dataset state the app reads from localStorage; written by the next go()."""
        items = {"gv.theme": "light", "gv.dataset": DATASET}
        if pinned is not None:
            items[f"gv.{DATASET}.pinned"] = json.dumps(pinned)
        if filters is not None:
            items[f"gv.{DATASET}.filters"] = json.dumps(filters)
        if locus is not None:
            items[f"gv.{DATASET}.locus"] = json.dumps(locus)
        for k, v in (extra or {}).items():
            items[k] = v if isinstance(v, str) else json.dumps(v)
        self.pending = items

    def go(self, path: str, settle: float = 1.5, timeout: float = 240):
        # Leave the current page for one that keeps no locus (pages save their state while they
        # unmount), write the state, then load the target page from scratch so the app reads it.
        self.pg.goto(self.base + "#/system")
        self.pg.reload()
        self.pg.wait_for_load_state("networkidle")
        items = getattr(self, "pending", None) or {"gv.theme": "light", "gv.dataset": DATASET}
        self.pg.evaluate("(items) => { for (const [k, v] of Object.entries(items)) localStorage.setItem(k, v); }", items)
        self.pending = None
        self.pg.goto(self.base + "#/" + path)
        self.pg.reload()
        self.settle(settle, timeout)

    def settle(self, extra: float = 1.0, timeout: float = 240):
        """Wait for the network to go quiet and every spinner / job overlay to disappear."""
        self.pg.wait_for_load_state("networkidle")
        deadline = time.time() + timeout
        busy = "() => [...document.querySelectorAll('.spinner, .overlay, .loading-line')].some((e) => e.offsetParent !== null)"
        while self.pg.evaluate(busy):
            if time.time() > deadline:
                print("  (still busy after timeout; capturing anyway)")
                break
            self.pg.wait_for_timeout(500)
        self.pg.wait_for_timeout(int(extra * 1000))

    def wait_for(self, selector: str, timeout: float = 240):
        self.pg.wait_for_selector(selector, timeout=timeout * 1000)

    # ------------------------------------------------------------------ callouts
    def callouts(self, items):
        """Numbered badges: items are (n, locator, where) with where in {'left','right','top','bottom','center'} or (dx, dy)."""
        self.pg.add_style_tag(content=CALLOUT_CSS)
        for n, loc, where in items:
            try:
                box = loc.bounding_box(timeout=3000)
            except Exception:
                box = None
            if not box:
                print(f"  callout {n}: element not visible")
                continue
            if isinstance(where, tuple):
                x, y = box["x"] + where[0], box["y"] + where[1]
            elif where == "right":
                x, y = box["x"] + box["width"] + 6, box["y"] + box["height"] / 2 - 13
            elif where == "top":
                x, y = box["x"] + box["width"] / 2 - 13, box["y"] - 30
            elif where == "bottom":
                x, y = box["x"] + box["width"] / 2 - 13, box["y"] + box["height"] + 4
            elif where == "center":
                x, y = box["x"] + box["width"] / 2 - 13, box["y"] + box["height"] / 2 - 13
            else:  # left
                x, y = box["x"] - 32, box["y"] + box["height"] / 2 - 13
            self.pg.evaluate(
                "([n, x, y]) => { const d = document.createElement('div'); d.className = 'doc-callout'; d.textContent = n;"
                " d.style.left = x + 'px'; d.style.top = y + 'px'; document.body.appendChild(d); }",
                [n, max(2, x), max(2, y)],
            )

    def ring(self, loc, pad=4):
        box = loc.bounding_box()
        if not box:
            return
        self.pg.add_style_tag(content=CALLOUT_CSS)
        self.pg.evaluate(
            "([x, y, w, h]) => { const d = document.createElement('div'); d.className = 'doc-ring';"
            " Object.assign(d.style, {left: x + 'px', top: y + 'px', width: w + 'px', height: h + 'px'}); document.body.appendChild(d); }",
            [box["x"] - pad, box["y"] - pad, box["width"] + 2 * pad, box["height"] + 2 * pad],
        )

    def clear_callouts(self):
        self.pg.evaluate("() => document.querySelectorAll('.doc-callout, .doc-ring').forEach((e) => e.remove())")

    # ------------------------------------------------------------------ capture
    def snap(self, name: str, locator=None, clip=None, full_main=False):
        """Screenshot the viewport, an element, or a clip rectangle, and save it as WebP."""
        from PIL import Image

        self.pg.mouse.move(2, VIEWPORT["height"] - 2)  # no hover tooltips
        self.pg.evaluate("() => document.querySelectorAll('.toast').forEach((t) => t.remove())")
        self.pg.wait_for_timeout(200)
        if locator is not None:
            # The app scrolls inside its main panel, so an element taller than the viewport would be
            # cut off: grow the viewport to fit it for this one capture.
            locator.scroll_into_view_if_needed()
            box = locator.bounding_box()
            grow = box and box["height"] > VIEWPORT["height"] - 80
            if grow:
                self.pg.set_viewport_size({"width": VIEWPORT["width"], "height": int(box["height"]) + 200})
                self.pg.wait_for_timeout(800)
                locator.scroll_into_view_if_needed()
            png = locator.screenshot()
            if grow:
                self.pg.set_viewport_size(VIEWPORT)
        elif clip is not None:
            png = self.pg.screenshot(clip=clip)
        else:
            png = self.pg.screenshot()
        img = Image.open(io.BytesIO(png)).convert("RGB")
        path = self.out / f"{name}.webp"
        img.save(path, "WEBP", quality=self.quality, method=6)
        print(f"  wrote {path.relative_to(REPO) if path.is_relative_to(REPO) else path} {img.size[0]}x{img.size[1]} {path.stat().st_size // 1024} KB")
        self.clear_callouts()

    def text(self, text, exact=False):
        return self.pg.get_by_text(text, exact=exact).first

    def button(self, name, exact=False):
        return self.pg.get_by_role("button", name=name, exact=exact).first


# ---------------------------------------------------------------------- the shots
TYRP1_TRACKS = ["rna_seq:1", "rna_seq:4"]  # melanocyte of skin, + and - strands
TYRP1_LOCUS = {
    "gene": "TYRP1", "output": None, "coords": "reference", "mode": "individuals", "hap": "H1+H2",
    "tracks": TYRP1_TRACKS, "groupField": "superpopulation", "groups": None, "yScale": "linear",
    "sharedY": False, "lockY": False, "diff": False, "envelope": True, "band": True, "popTrack": "rna_seq:1",
    "popRows": 600, "seqSource": "pinned", "showRef": True, "showObserved": False, "obsScale": "tpm",
    "start": 247235, "end": 277052,
}
# rs12913832 (HERC2 intron 86, the OCA2 enhancer; G = blue-eye allele), GRCh38
PIGMENTATION_GENES = ["DDB1", "EDAR", "HERC2", "MC1R", "MFSD12", "OCA2", "SLC24A5", "SLC45A2", "TCHH", "TYR", "TYRP1"]
RS12913832 = {"gene": "HERC2", "pos": 28120472, "ref": "A", "alt": "G"}


def drawer(s: Shooter):
    s.wait_for("aside.drawer")
    s.settle(2.5, timeout=90)
    return s.pg.locator("aside.drawer").first


# ---- Overview
def shot_overview(s: Shooter):
    s.prime(pinned=PINNED, filters={})
    s.go("overview")
    s.snap("overview")


def shot_app_shell(s: Shooter):
    s.prime(pinned=PINNED, filters={})
    s.go("overview")
    pg = s.pg
    s.callouts([
        (1, pg.locator("#datasetSelect"), "right"),
        (2, pg.locator("#sidenav a").nth(0), (150, 4)),
        (3, pg.locator("#cohortChip"), "left"),
        (4, pg.locator("#themeToggle"), "left"),
        (5, s.button("Browse samples"), "left"),
    ])
    s.snap("app-shell")


def shot_gene_card(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("overview")
    row = s.pg.locator("tr", has_text="TYRP1").first
    row.scroll_into_view_if_needed()
    row.locator("button").first.click()
    drawer(s)
    s.snap("gene-card")


def shot_predict_form(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("overview")
    s.button("AlphaGenome predictions").click()
    drawer(s)
    s.snap("predict-form")


# ---- Samples
def shot_samples(s: Shooter):
    s.prime(pinned=PINNED, filters={"superpopulation": ["AFR", "EUR"]})
    s.go("samples")
    pg = s.pg
    s.callouts([
        (1, pg.get_by_text("Superpopulation", exact=True).first, "left"),
        (2, s.button("Pin first 5"), "left"),
        (3, s.button("Save as training view"), "bottom"),
        (4, s.button("Compare pinned in Tracks"), "right"),
        (5, pg.get_by_text("Region scalars", exact=True).first, "left"),
    ])
    s.snap("samples-cohort")


def shot_samples_pca(s: Shooter):
    s.prime(pinned=PINNED, filters={})
    s.go("samples?view=ancestry", settle=2)
    compute = s.button("Compute PCA")
    if compute.count():
        compute.click()
        s.settle(3, timeout=900)
    s.snap("samples-pca")


def shot_samples_matching(s: Shooter):
    """Match the pigmentation classes on ancestry (no pairs exist), then delete the derived field."""
    s.prime(pinned=PINNED, filters={})
    s.go("samples?view=ancestry", settle=2)
    compute = s.button("Compute PCA")
    if compute.count():
        compute.click()
        s.settle(3, timeout=900)
    card = s.pg.locator("section.card", has_text="Match two groups on ancestry").first
    card.scroll_into_view_if_needed()
    hosts = card.locator(".facet-values")
    hosts.nth(0).get_by_text("weak pigmentation", exact=False).first.click()
    hosts.nth(1).get_by_text("strong pigmentation", exact=False).first.click()
    name = card.locator("input.mono").first
    name.fill("docs_matched")
    card.get_by_role("button", name="Match", exact=True).click()
    s.settle(2, timeout=300)
    print("  weak vs strong pigmentation:", s.pg.locator(".toast").all_inner_texts())  # no pairs exist
    # A comparison that can be matched: AMR vs EUR.
    card.locator("select").first.select_option("superpopulation")
    s.pg.wait_for_timeout(300)
    hosts.nth(0).get_by_text("AMR", exact=False).first.click()
    hosts.nth(1).get_by_text("EUR", exact=False).first.click()
    card.get_by_role("button", name="Match", exact=True).click()
    s.settle(2, timeout=300)
    s.snap("samples-matching", locator=card)
    s.pg.evaluate(f"() => fetch('/api/d/{DATASET}/scalars/delete', {{method: 'POST', headers: {{'Content-Type': 'application/json'}}, body: JSON.stringify({{name: 'docs_matched'}})}})")


def shot_region_scalar(s: Shooter):
    s.prime(pinned=PINNED, filters={})
    s.go("samples")
    s.button("New").click()
    s.pg.wait_for_timeout(1500)
    target = s.pg.locator("aside.drawer, .modal").first
    s.settle(1)
    s.snap("region-scalar-form")


# ---- Tracks
def shot_tracks_individuals(s: Shooter):
    s.prime(pinned=PINNED, filters={}, locus=dict(TYRP1_LOCUS))
    s.go("tracks", settle=3)
    pg = s.pg
    s.callouts([
        (1, pg.get_by_text("Individuals", exact=True).first, "left"),
        (2, pg.locator("input.locus-input").first, "right"),
        (3, pg.locator(".panel .lane-label").first, "right"),
        (4, pg.get_by_text("Compare with", exact=True).first, "left"),
        (5, pg.get_by_text("Tracks (", exact=False).first, "left"),
    ])
    s.snap("tracks-individuals")


def shot_tracks_observed(s: Shooter):
    s.prime(pinned=PINNED, filters={}, locus=dict(TYRP1_LOCUS, showObserved=True))
    s.go("tracks", settle=4)
    s.snap("tracks-observed")


def shot_tracks_groups(s: Shooter):
    s.prime(pinned=PINNED, filters={}, locus=dict(TYRP1_LOCUS, mode="groups", groupField="superpopulation"))
    s.go("tracks", settle=4, timeout=600)
    s.snap("tracks-group-means")


def shot_tracks_heatmap(s: Shooter):
    s.prime(pinned=PINNED, filters={}, locus=dict(TYRP1_LOCUS, mode="population", groupField="superpopulation"))
    s.go("tracks", settle=4, timeout=600)
    s.snap("tracks-heatmap")


def shot_track_card(s: Shooter):
    s.prime(pinned=PINNED, filters={}, locus=dict(TYRP1_LOCUS))
    s.go("tracks", settle=2)
    s.pg.locator(".panel .lane-label.clickable").first.click()
    drawer(s)
    s.snap("track-card")


def shot_tracks_export(s: Shooter):
    """The export menu, and the exported SVG itself (docs/assets/visualizer/figure-tyrp1.svg)."""
    s.prime(pinned=PINNED, filters={}, locus=dict(TYRP1_LOCUS, showObserved=True))
    s.go("tracks", settle=4)
    s.button("Export").click()
    s.pg.wait_for_timeout(500)
    s.snap("tracks-export-menu", clip={"x": 300, "y": 120, "width": 700, "height": 360})
    with s.pg.expect_download(timeout=120000) as dl:
        s.pg.get_by_text("SVG", exact=True).first.click()
    path = s.out / "figure-tyrp1.svg"
    dl.value.save_as(str(path))
    print(f"  wrote {path.name} {path.stat().st_size // 1024} KB")


# ---- Sequence
def shot_sequence(s: Shooter):
    # window offsets around rs12913832 (HERC2 window starts at chr15:27,954,466)
    off = RS12913832["pos"] - 27954466
    s.prime(pinned=PINNED, filters={})
    s.go(f"sequence?gene=HERC2&start={off - 70}&end={off + 70}", settle=3)
    s.snap("sequence-bases")


def shot_sequence_wide(s: Shooter):
    s.prime(pinned=PINNED, filters={})
    s.go("sequence?gene=HERC2", settle=2)
    s.button("Model window").click()
    s.settle(3)
    s.snap("sequence-overview")


# ---- Variant
def shot_variant(s: Shooter):
    v = RS12913832
    s.prime(pinned=PINNED, filters={})
    s.go(f"variant?gene={v['gene']}&pos={v['pos']}&ref={v['ref']}&alt={v['alt']}&target=OCA2", settle=4, timeout=900)
    s.snap("variant-header")
    for name, title in (("variant-direction", "Direction"), ("variant-genotype", "AlphaGenome by genotype"), ("variant-gtex", "GTEx v8 eQTL")):
        card = s.pg.locator(".card", has_text=title).first
        try:
            card.scroll_into_view_if_needed(timeout=5000)
            s.snap(name, locator=card)
        except Exception as exc:
            print(f"  {name}: {exc.__class__.__name__}")


# ---- Gene products
def shot_products(s: Shooter, gene="MC1R"):
    s.prime(pinned=PINNED, filters={})
    s.go(f"products?gene={gene}&sample={SAMPLE}&tissue=CL:1000458", settle=4, timeout=600)
    s.snap(f"products-{gene.lower()}-top")
    cards = s.pg.locator(".report-host > .card")
    for i in range(min(cards.count(), 3)):
        card = cards.nth(i)
        card.scroll_into_view_if_needed()
        s.pg.wait_for_timeout(800)
        s.snap(f"products-{gene.lower()}-step{i + 1}", locator=card)
    shape = s.pg.locator(".card", has_text="Protein shape").last
    try:
        shape.scroll_into_view_if_needed(timeout=10000)
        s.settle(6, timeout=300)  # AlphaFold DB model + 3Dmol.js
        s.snap(f"products-{gene.lower()}-step3", locator=shape)
    except Exception as exc:
        print(f"  step3: {exc.__class__.__name__}")


def shot_products_mc1r(s: Shooter):
    shot_products(s, "MC1R")


def shot_products_tyr(s: Shooter):
    shot_products(s, "TYR")


def shot_products_more(s: Shooter):
    s.prime(pinned=PINNED, filters={})
    s.go(f"products?gene=TYR&sample={SAMPLE}&tissue=CL:1000458", settle=4, timeout=600)
    more = s.pg.get_by_text("More evidence", exact=False).first
    more.scroll_into_view_if_needed()
    more.click()
    s.settle(3, timeout=600)
    for name, title in (("products-how-much", "How much"), ("products-splicing", "Splicing (AlphaGenome junctions)"), ("products-transcripts", "Annotated transcripts")):
        card = s.pg.locator(".card", has_text=title).last
        try:
            card.scroll_into_view_if_needed(timeout=5000)
            s.settle(2, timeout=1200)  # other tissues are predicted on demand
            s.snap(name, locator=card)
        except Exception as exc:
            print(f"  {name}: {exc.__class__.__name__}")


# ---- Experiments
PIGMENTATION_RUN = "cnn2_pigmentation_rna_seq_H1+H2_haplotype_channels_32768_log_s1k6x32f16_s2f32_s3f64_gpavg_fc256_L100-40_relu_0.5_adam"


def shot_experiments(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("experiments", settle=2)
    s.snap("experiments-runs")


def shot_experiment_detail(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("experiments", settle=2)
    s.pg.locator("tr", has_text="cnn2_superpopulation_rna_seq_H1+H2_haplotype_channels_signals_and_masks").first.locator("td").nth(1).click()
    s.settle(3)
    s.snap("experiments-run")
    target = s.pg.locator(".card", has_text="Confusion").first
    try:
        target.scroll_into_view_if_needed(timeout=5000)
        s.pg.wait_for_timeout(600)
        s.snap("experiments-confusion")
    except Exception as exc:
        print(f"  confusion: {exc.__class__.__name__}")


def shot_training_form(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("experiments", settle=2)
    s.button("New training run").click()
    d = drawer(s)
    s.snap("training-form")
    # Train on the 11 pigmentation windows only, then show both negative controls.
    s.pg.evaluate("""(keep) => {
        for (const label of document.querySelectorAll('aside.drawer label')) {
            const box = label.querySelector('input[type=checkbox]');
            const name = label.textContent.trim();
            if (box && /^[A-Z][A-Za-z0-9]+$/.test(name) && box.checked !== keep.includes(name)) box.click();
        }
    }""", PIGMENTATION_GENES)
    d.get_by_text("Shuffle within", exact=False).first.click()
    s.pg.wait_for_timeout(500)
    d.get_by_role("button", name="Matched control windows").click()
    s.settle(2, timeout=600)
    # the innermost block holding both the control switch and the matched pairs
    section = d.locator("div", has=s.pg.get_by_text("Matched control windows in use")).filter(has=s.pg.get_by_text("Negative control", exact=True)).last
    section.scroll_into_view_if_needed()
    s.pg.wait_for_timeout(600)
    s.snap("training-controls", locator=section)


# ---- Perturbation Lab
def _lab(s: Shooter, gene="TYR"):
    s.prime(pinned=PINNED)
    s.go(f"perturb?run={quote(PIGMENTATION_RUN)}&sample={SAMPLE}&gene={gene}", settle=2)
    load = s.pg.get_by_role("button", name="Load model", exact=True)
    if load.count():  # the server keeps the last loaded model ("Loaded")
        load.first.click()
        s.settle(3, timeout=900)


def shot_perturb(s: Shooter):
    _lab(s)
    s.snap("perturb-loaded")
    s.button("Model window").click()
    s.settle(1)
    for _ in range(3):
        s.button("Zoom in").click()
        s.pg.wait_for_timeout(300)
    s.settle(2)
    s.button("Select view").click()
    s.pg.get_by_text("Scramble", exact=True).first.click()
    s.button("+ Add edit").click()
    s.pg.wait_for_timeout(500)
    s.button("Run AlphaGenome + model").click()
    s.settle(4, timeout=900)
    s.snap("perturb-edited")


def shot_perturb_scan(s: Shooter):
    _lab(s)
    s.button("Model window").click()
    s.settle(1)
    s.button("Zoom in").click()
    s.settle(1)
    scan = s.pg.locator(".side-section", has_text="Saturation scan").first
    scan.scroll_into_view_if_needed()
    scan.get_by_text("View", exact=True).first.click()
    size = scan.locator("input[type=number]").nth(0)
    size.fill("1024")
    size.dispatch_event("input")
    scan.get_by_text("H2", exact=True).first.click()
    s.pg.wait_for_timeout(300)
    scan.get_by_role("button", name="Run scan").click()
    s.settle(4, timeout=1800)
    s.snap("perturb-scan")


# ---- Jobs, AlphaGenome, System, Import
def shot_jobs(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("jobs", settle=2)
    s.snap("jobs")


def shot_alphagenome(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("alphagenome", settle=2)
    s.snap("alphagenome-backend")


def shot_system(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("system", settle=3)
    s.snap("system")


def shot_import(s: Shooter):
    s.prime(pinned=PINNED)
    s.go("import", settle=3)
    s.snap("import-quickstart")
    step = s.pg.get_by_text("Sample metadata", exact=False).first
    try:
        step.scroll_into_view_if_needed(timeout=5000)
        s.pg.wait_for_timeout(800)
        s.snap("import-steps")
    except Exception as exc:
        print(f"  steps: {exc.__class__.__name__}")


SHOTS = {
    "overview": shot_overview,
    "app-shell": shot_app_shell,
    "gene-card": shot_gene_card,
    "predict-form": shot_predict_form,
    "samples-cohort": shot_samples,
    "samples-pca": shot_samples_pca,
    "samples-matching": shot_samples_matching,
    "region-scalar-form": shot_region_scalar,
    "tracks-individuals": shot_tracks_individuals,
    "tracks-observed": shot_tracks_observed,
    "tracks-group-means": shot_tracks_groups,
    "tracks-heatmap": shot_tracks_heatmap,
    "track-card": shot_track_card,
    "tracks-export": shot_tracks_export,
    "sequence-bases": shot_sequence,
    "sequence-overview": shot_sequence_wide,
    "variant": shot_variant,
    "products-mc1r": shot_products_mc1r,
    "products-tyr": shot_products_tyr,
    "products-more": shot_products_more,
    "experiments-runs": shot_experiments,
    "experiments-run": shot_experiment_detail,
    "training-form": shot_training_form,
    "perturb": shot_perturb,
    "perturb-scan": shot_perturb_scan,
    "jobs": shot_jobs,
    "alphagenome-backend": shot_alphagenome,
    "system": shot_system,
    "import": shot_import,
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--url", default="http://127.0.0.1:8780/", help="Running visualizer (default: %(default)s)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output directory (default: docs/assets/visualizer)")
    parser.add_argument("--only", default="", help="Comma-separated shot names (default: all)")
    parser.add_argument("--scale", type=float, default=1.5, help="Device scale factor (default: %(default)s)")
    parser.add_argument("--quality", type=int, default=82, help="WebP quality (default: %(default)s)")
    parser.add_argument("--list", action="store_true", help="List the shot names and exit")
    args = parser.parse_args(argv)
    if args.list:
        print("\n".join(SHOTS))
        return 0
    names = [n for n in args.only.split(",") if n] or list(SHOTS)
    unknown = [n for n in names if n not in SHOTS]
    if unknown:
        parser.error(f"unknown shots: {', '.join(unknown)} (see --list)")
    from playwright.sync_api import sync_playwright

    args.out.mkdir(parents=True, exist_ok=True)
    failed = []
    with sync_playwright() as p:
        browser = p.chromium.launch(args=["--use-angle=swiftshader", "--enable-unsafe-swiftshader"])
        ctx = browser.new_context(viewport=VIEWPORT, device_scale_factor=args.scale, accept_downloads=True)
        page = ctx.new_page()
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(args.url)
        shooter = Shooter(page, args.url, args.out, args.scale, args.quality)
        for name in names:
            print(name)
            errors.clear()
            try:
                SHOTS[name](shooter)
            except Exception as exc:  # keep going; report at the end
                failed.append(name)
                print(f"  FAILED: {exc.__class__.__name__}: {str(exc).splitlines()[0]}")
            if errors:
                print(f"  page errors: {errors}")
        browser.close()
    if failed:
        print(f"failed: {', '.join(failed)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
