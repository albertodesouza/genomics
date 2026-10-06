"""Browser smoke tests for the visualizer's single-page app.

Every page is opened in headless Chromium against the synthetic dataset of ``test_visualizer`` and
must render without JavaScript errors; the Tracks page's figure export is round-tripped (PNG and
vector SVG). Skipped unless Playwright and its Chromium are installed:

    python3 -m pip install -e ".[test-ui]" && python3 -m playwright install chromium
"""
import threading

import pytest

from test_visualizer import _free_port, _write_reference_prediction, dataset_dir  # noqa: F401 (fixture)
from test_visualizer_variants import GTEX_ANSWERS, INS, FakeRemote

sync_api = pytest.importorskip("playwright.sync_api")

PAGES = ["overview", "samples", "tracks", "sequence", "variant", "perturb", "experiments", "jobs", "alphagenome", "system", "import"]
# Console noise the app does not control: the browser asks for a favicon, and it logs a line for
# every failed response. A 503 is the server saying this machine cannot provide a feature (no
# AlphaGenome backend on a CI runner); the pages handle it, which
# test_import_page_without_an_alphagenome_backend checks. Any other status still fails a test.
IGNORED_CONSOLE = ("favicon", "status of 503")
APPS = {}


@pytest.fixture(scope="module")
def server_url(dataset_dir, tmp_path_factory):  # noqa: F811 (fixture imported above)
    from genomics.visualizer.alphagenome import AlphaGenomeBackend, LocalAlphaGenomeServer
    from genomics.visualizer.datasets import DatasetCatalog
    from genomics.visualizer.server import Handler, Server, VisualizerApp

    _write_reference_prediction(dataset_dir)
    tmp = tmp_path_factory.mktemp("viz_ui")
    catalog = DatasetCatalog()
    catalog.add(dataset_dir)
    # A backend with no server checkout and settings in the temp dir (never ~/.config).
    backend = AlphaGenomeBackend(LocalAlphaGenomeServer(server_dir=None, python=None, log_dir=tmp), settings_path=tmp / "alphagenome.json")
    app = VisualizerApp(catalog, cache_dir=tmp / "cache", memory_bytes=256 << 20, workers=2, runs_roots=[], remote=False, alphagenome=backend)
    from genomics.visualizer.gtex import GtexClient

    app.gtex = GtexClient(FakeRemote(GTEX_ANSWERS))  # canned GTEx answers (no network in tests)
    handler = type("SmokeHandler", (Handler,), {"app": app})
    port = _free_port()
    server = Server(("127.0.0.1", port), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    APPS[f"http://127.0.0.1:{port}"] = app  # for tests that seed caches
    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        server.shutdown()
        server.server_close()
        app.shutdown()


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as p:
        try:
            b = p.chromium.launch()
        except Exception as exc:  # browser binaries not downloaded
            pytest.skip(f"Chromium for Playwright is not installed ({exc.__class__.__name__}); run `python3 -m playwright install chromium`")
        yield b
        b.close()


@pytest.fixture
def page(browser):
    context = browser.new_context(viewport={"width": 1400, "height": 900}, accept_downloads=True)
    pg = context.new_page()
    pg.errors = []
    pg.on("pageerror", lambda e: pg.errors.append(f"pageerror: {e}"))
    pg.on("console", lambda m: pg.errors.append(f"console: {m.text}") if m.type == "error" and not any(s in m.text for s in IGNORED_CONSOLE) else None)
    yield pg
    context.close()


def _wait_state(page, expression, timeout=15.0):
    """Poll ``expression`` (JS over ``m``, the app's state module) until it is truthy."""
    import time

    deadline = time.time() + timeout
    script = f"async () => {{ const m = await import('/static/js/state.js'); return !!({expression}); }}"
    while not page.evaluate(script):
        if time.time() > deadline:
            raise AssertionError(f"timed out waiting for {expression}")
        page.wait_for_timeout(100)


def _open(page, url):
    page.goto(url)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(300)


# The Import page reports tools the machine lacks; on a runner without bcftools/samtools that is a
# correct report, not a rendering failure.
EXPECTED_BOXES = ("Missing on the server:",)


@pytest.mark.parametrize("name", PAGES)
def test_page_renders_without_errors(page, server_url, name):
    _open(page, f"{server_url}/#/{name}")
    assert page.locator("#page").inner_text().strip(), f"{name} rendered nothing"
    boxes = [t for t in page.locator(".error-box").all_inner_texts() if not t.startswith(EXPECTED_BOXES)]
    assert boxes == [], f"{name}: {boxes}"
    # Error toasts are app failures, except missing data in the fixture (EXTRA has no haplotype FASTA).
    toasts = [t for t in page.locator(".toast.error").all_inner_texts() if "No such file or directory" not in t]
    assert toasts == [], f"{name}: {toasts}"
    assert page.errors == [], f"{name}: {page.errors}"


def test_tracks_draws_canvases_for_pinned_sample(page, server_url):
    _open(page, f"{server_url}/#/samples")
    # Pin a sample through the app's own state so Tracks has a series to draw.
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S1']); }")
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    sizes = page.evaluate("[...document.querySelectorAll('canvas')].map((c) => c.width * c.height)")
    assert sizes and all(s > 1 for s in sizes)
    assert page.errors == []


@pytest.mark.parametrize("fmt,label", [("png", "PNG"), ("svg", "SVG")])
def test_tracks_figure_export(page, server_url, tmp_path, fmt, label):
    _open(page, f"{server_url}/#/samples")
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S1', 'S2']); }")
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    page.get_by_role("button", name="Export").click()
    with page.expect_download() as info:
        page.get_by_role("menuitem").filter(has_text=label).first.click()
    download = info.value
    assert download.suggested_filename.endswith(f".{fmt}") and "GENE1" in download.suggested_filename
    out = tmp_path / download.suggested_filename
    download.save_as(out)
    data = out.read_bytes()
    if fmt == "png":
        assert data[:8] == b"\x89PNG\r\n\x1a\n" and len(data) > 2000
    else:
        svg = data.decode("utf-8")
        assert svg.lstrip().startswith("<?xml") and "<svg" in svg
        assert "GENE1" in svg and "S1" in svg  # title and legend
        assert "<path" in svg and "data:image" not in svg  # plots are vector, not embedded bitmaps
    assert page.errors == []


def test_variant_page_shows_genotypes_effect_and_gtex(page, server_url, tmp_path):
    pos, ref, alt = INS
    _open(page, f"{server_url}/#/variant?gene=GENE1&pos={pos}&ref={ref}&alt={alt}&region=custom&custom=chr1:1001-1300&track=1")
    page.wait_for_selector(".verdict", timeout=20000)
    text = page.locator("#page").inner_text()
    assert f"chr1:{pos:,}" in text and "rs1" in text  # header with the rsID from GTEx
    assert "-0.5" in text  # AlphaGenome slope per ALT allele (track 1 = per-sample constant)
    assert page.locator(".verdict").first.inner_text() in ("agree", "disagree", "AlphaGenome n.s.", "GTEx n.s.", "both n.s.")
    assert page.locator("canvas").count() >= 2  # AlphaGenome by genotype + GTEx tissue
    page.get_by_role("button", name="Export").click()
    with page.expect_download() as info:
        page.get_by_role("menuitem").filter(has_text="SVG").first.click()
    out = tmp_path / info.value.suggested_filename
    info.value.save_as(out)
    svg = out.read_text(encoding="utf-8")
    assert "GENE1" in svg and "<path" in svg
    assert page.errors == []


def test_variant_page_lists_sites_and_opens_one(page, server_url):
    _open(page, f"{server_url}/#/variant?gene=GENE1")
    page.wait_for_selector("table tbody tr")
    page.locator("table tbody tr").first.click()
    page.wait_for_selector(".variant-title")
    assert "pos=" in page.url
    assert page.errors == []


def test_region_scalar_from_tracks_becomes_a_sample_facet(page, server_url):
    _open(page, f"{server_url}/#/samples")
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S1']); }")
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    page.get_by_role("button", name="Region scalar…").click()
    dialog = page.locator(".modal")
    dialog.locator("label.field", has_text="Track").locator("select").select_option("1")  # per-sample constant
    assert dialog.locator("label.field", has_text="Genomic range").locator("input").input_value().startswith("chr1:")
    assert dialog.locator("label.field", has_text="Field name").locator("input").input_value() == "rna_seq_region"
    dialog.get_by_role("button", name="Compute").click()
    page.wait_for_selector(".modal", state="detached", timeout=20000)

    _open(page, f"{server_url}/#/samples")
    facet = page.locator('[data-facet="region-scalars"]')
    facet.locator(".scalar-name", has_text="rna_seq_region").wait_for()
    assert facet.locator("canvas").count() == 1
    bins = page.locator('[data-facet="rna_seq_region_bin"]')
    assert sorted(bins.locator(".fv-name").all_inner_texts()) == ["Q1", "Q2", "Q3"]
    cols = page.evaluate("async () => { const m = await import('/static/js/state.js'); const c = m.state.samples.col; return m.state.samples.rows.map((r) => [r[c.sample_id], r[c.rna_seq_region]]); }")
    assert sorted(cols) == [["S1", 1.0], ["S2", 2.0], ["S3", 3.0]]
    page.once("dialog", lambda d: d.accept())
    facet.get_by_role("button", name="Delete").click()
    page.wait_for_function("!document.querySelector('[data-facet=\"rna_seq_region_bin\"]')")
    assert page.errors == []


def test_ancestry_pca_scatter_pc_fields_and_matching(page, server_url):
    import numpy as np

    app = APPS[server_url]
    dataset = app.catalog.all()[0]
    params = app.ancestry.params(dataset, ["GENE1"], 0.05, 2000, 10)  # the panel's defaults
    app.ancestry.memory.put(app.ancestry.key(dataset, params), {
        "samples": np.array(["S1", "S2", "S3"]), "scores": np.array([[1.0, 0.2], [-1.0, 0.5], [0.8, -0.5]], np.float32),
        "explained": np.array([0.5, 0.2], np.float32), "n_sites": np.int64(12), "sites_per_window": np.array([12]),
    })
    _open(page, f"{server_url}/#/samples?view=ancestry")
    page.locator(".ancestry canvas").wait_for()  # a cached PCA is shown without asking
    assert "PC1 (50.0%)" in page.locator(".ancestry").inner_text()
    assert page.locator(".ancestry .legend-item").count() == 2  # AFR, EUR
    page.get_by_role("button", name="Add PC fields").click()
    _wait_state(page, "'pc1' in m.state.samples.col")
    page.get_by_role("button", name="Update PC fields").wait_for()  # the panel came back after the page remount

    match = page.locator("section.card", has_text="Match two groups")
    match.locator("label.field", has_text="Field").first.locator("select").select_option("superpopulation")
    match.locator("label.field", has_text="Group A").locator("label.check", has_text="EUR").locator("input").check()
    match.locator("label.field", has_text="Group B").locator("label.check", has_text="AFR").locator("input").check()
    match.locator("label.field", has_text="Caliper").locator("input").fill("0")
    match.get_by_role("button", name="Match").click()
    page.locator(".ancestry", has_text="1 pairs").wait_for()
    _wait_state(page, "'matched' in m.state.samples.col")
    values = page.evaluate("async () => { const m = await import('/static/js/state.js'); const c = m.state.samples.col; return m.state.samples.rows.map((r) => [r[c.sample_id], r[c.matched] ?? null]); }")
    assert sorted(values) == [["S1", "AFR"], ["S2", "EUR"], ["S3", None]]  # on PC1-PC2, S1 is nearer S2
    page.get_by_role("button", name="Use as cohort").first.click()
    _wait_state(page, "m.cohortSize() === 2")
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setFilters({}); }")
    assert page.errors == []


def test_perturbation_saturation_scan_lane(page, server_url):
    """The scan section posts the scan, polls it, draws the lane and turns a window into an edit (the
    model and AlphaGenome are stubbed; the scan itself is tested in test_visualizer_workflows)."""
    import json
    import re

    from test_visualizer import L

    dataset = APPS[server_url].catalog.all()[0]
    mw_start = L // 2 - 20  # window_center_size 40, centred like the page's modelWindow()
    model = {"run": "run1", "checkpoint": "best.pt", "model": "CNN2", "target": "superpopulation", "label_field": "superpopulation",
             "classes": ["AFR", "EUR"], "class_counts": {"AFR": 2, "EUR": 1}, "genes": ["GENE1"], "outputs": ["rna_seq"], "ontology_terms": [],
             "window_center_size": 40, "dataset_id": dataset.id, "dataset_path": str(dataset.path), "device": "cpu",
             "samples": [{"id": "S1", "label": "AFR", "split": "test"}, {"id": "S2", "label": "EUR", "split": "train"}]}
    rows = [{"start": mw_start + 10 * k, "end": mw_start + 10 * k + 10, "changed": 0 if k == 3 else 10,
             "edited": [0.6 + d, 0.4 - d], "delta": [d, -d]} for k, d in enumerate([-0.2, 0.05, 0.1, 0.0])]
    posted = []

    def handle(route):
        path = re.sub(r"\?.*$", "", route.request.url).split("/api/perturb/")[1]
        if path == "models":
            body = {"runs": [{"id": "run1", "name": "run1", "compatible": True, "checkpoints": ["best.pt"], "default_checkpoint": "best.pt"}],
                    "loaded": {"run": "run1", "checkpoint": "best.pt"}, "default": "run1", "training": [], "backend": {"label": "stub", "reasons": []}}
        elif path == "model":
            body = model
        elif path == "score":
            body = {"sample": "S1", "label": "AFR", "split": "test", "classes": ["AFR", "EUR"], "probabilities": [0.6, 0.4]}
        elif path == "sequence":
            body = {"gene": "GENE1", "sample": "S1", "start": 0, "end": L, "domain": L, "mode": "density", "edges": [0, L], "rows": []}
        elif path == "scan":
            posted.append(json.loads(route.request.post_data))
            spec = posted[-1]
            body = {"pending": True, "job": {"id": "j1", "progress": 0.5, "message": "Window 2/4", "elapsed": 1}} if len(posted) == 1 else {
                **spec, "key": "k", "classes": ["AFR", "EUR"], "baseline": [0.6, 0.4], "label": "AFR", "model_window": {}, "rows": rows, "calls": 6, "elapsed": 2.5}
        else:
            return route.continue_()
        route.fulfill(status=200, content_type="application/json", body=json.dumps(body))

    page.route(re.compile(r".*/api/perturb/.*"), handle)
    _open(page, f"{server_url}/#/perturb?run=run1&sample=S1&gene=GENE1")
    section = page.locator(".side-section", has_text="Saturation scan")
    section.wait_for()
    assert "0 windows" in section.inner_text()  # 1,024 bp windows do not fit the 40 bp model window
    section.locator("label.field", has_text="Window (bp)").locator("input").fill("10")
    assert "4 windows, up to 8 AlphaGenome calls" in section.inner_text()  # the step followed the window
    section.get_by_role("button", name="Run scan").click()
    page.locator(".side-section .scan-row").first.wait_for()
    assert len(posted) >= 2 and posted[0] == posted[-1]  # polled by re-posting the same body
    assert {k: posted[0][k] for k in ("op", "start", "end", "size", "step", "haplotypes")} == {
        "op": "scramble", "start": mw_start, "end": mw_start + 40, "size": 10, "step": 10, "haplotypes": ["H1", "H2"]}
    section = page.locator(".side-section", has_text="Saturation scan")
    ranked = section.locator(".scan-row").all_inner_texts()
    assert len(ranked) == 3 and ranked[0].startswith("-20.0 pp")  # unchanged window left out; largest |delta| first for the true class
    assert "1 window unchanged" in section.inner_text()
    assert page.evaluate("!document.querySelector('.scan-lane').hidden && document.querySelector('.scan-lane canvas').height > 1")
    section.locator("label.field", has_text="Show class").locator("select").select_option("EUR")
    assert page.locator(".side-section .scan-row").first.inner_text().startswith("+20.0 pp")
    page.locator(".side-section .scan-row").first.get_by_role("button", name="Add as an edit").click()
    edits = page.locator(".side-section", has_text="Edits").locator(".edit-item").all_inner_texts()
    assert len(edits) == 1 and "Scramble" in edits[0] and "H1+H2" in edits[0]
    assert page.errors == []


def test_train_form_negative_controls(page, server_url):
    """The training form's negative-control section: label shuffling and matched control windows."""
    import re

    _open(page, f"{server_url}/#/jobs")
    page.get_by_role("button", name="Train a model").click()
    form = page.locator(".drawer")
    control = form.locator("label.field", has_text="Negative control")
    control.wait_for()
    exactly_within = form.locator("label.field").filter(has=page.locator("span", has_text=re.compile(r"^Within$")))
    assert exactly_within.count() == 0  # only for "shuffle within"

    control.get_by_role("radio", name="Shuffle labels", exact=True).click()
    assert "nothing links genotype to class" in form.inner_text()
    control.get_by_role("radio", name="Shuffle within…", exact=True).click()
    within = exactly_within.locator("select")
    target = form.locator("label.field", has_text="Predict").locator("select").input_value()
    choices = within.evaluate("el => [...el.options].map((o) => o.value)")
    assert target not in choices and len(choices) >= 2  # within the target itself would be a no-op
    assert "no individual keeps their own label" in form.inner_text()
    other = next(c for c in choices if c != within.input_value())
    within.select_option(other)
    assert f"inside each {other}" in form.inner_text()

    # The matched-control button swaps the gene selection; the fixture has one spare window.
    assert "null" not in form.locator("div", has_text="Negative control").last.inner_text()

    # Matching needs windows outside the panel: with every window chosen there are none left.
    page.get_by_role("button", name="Matched control windows…").click()
    page.locator(".toast.error").first.wait_for()
    assert "fewer windows outside its gene panel" in page.locator(".toast.error").first.inner_text()

    # Leaving EXTRA out makes it the only candidate, and it has no reference prediction to measure.
    genes = form.locator("label.field", has_text="Gene windows")
    genes.locator("label", has_text="EXTRA").locator("input").uncheck()
    page.get_by_role("button", name="Matched control windows…").click()
    page.locator(".toast.error").nth(1).wait_for()
    assert "reference" in page.locator(".toast.error").nth(1).inner_text().lower()
    # Both refusals are 400s the form reports: no failed job, and nothing worse in the console.
    assert [e for e in page.errors if "400 (Bad Request)" not in e] == []


def test_locus_box_accepts_gene_names(page, server_url):
    """A gene name in the locus box jumps to it (coordinates and rsIDs are covered elsewhere)."""
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    box = page.locator(".locus-input")
    before = box.input_value()

    # The fixture has no gene annotations, but EXTRA is a window of the dataset: typing it switches.
    box.fill("EXTRA")
    box.press("Enter")
    page.locator(".toast", has_text="EXTRA").wait_for()
    _wait_state(page, "m.state.locus && m.state.locus.gene === 'EXTRA'")
    assert page.locator("select[aria-label='Gene']").input_value() == "EXTRA"

    page.locator("select[aria-label='Gene']").select_option("GENE1")
    _wait_state(page, "m.state.locus && m.state.locus.gene === 'GENE1'")
    box.fill("NOSUCHGENE")
    box.press("Enter")
    page.locator(".toast.error", has_text="No gene").wait_for()
    assert box.input_value() == "NOSUCHGENE" or box.input_value() == before  # the view did not move
    assert page.locator("select[aria-label='Gene']").input_value() == "GENE1"
    assert page.errors == []


def test_tracks_y_scale_lock_keeps_the_range_across_views(page, server_url):
    """With the lock on, moving the view must not rescale a panel; without it, it does."""
    _open(page, f"{server_url}/#/samples")
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S1']); }")
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    panel = page.locator(".panel canvas").first
    shot = lambda: panel.evaluate("c => c.toDataURL()")
    box = page.locator(".locus-input")

    def goto(text):
        box.fill(text)
        box.press("Enter")
        page.wait_for_timeout(400)

    goto("1-40")
    zoomed_auto = shot()
    page.get_by_text("Lock y-scale while panning").click()
    page.wait_for_timeout(300)
    assert shot() == zoomed_auto  # locking alone keeps the range it was switched on with
    _wait_state(page, "m.state.locus && m.state.locus.lockY === true")

    goto("60-100")
    locked_elsewhere = shot()
    page.get_by_text("Lock y-scale while panning").click()  # unlock: the same view, now autoscaled
    page.wait_for_timeout(400)
    assert shot() != locked_elsewhere
    assert page.errors == []


def test_saved_views_round_trip(page, server_url):
    """Save the Tracks view under a name, change things, then restore it from the Overview."""
    _open(page, f"{server_url}/#/samples")
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S1']); m.setFilters({ superpopulation: ['AFR'] }); }")
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    box = page.locator(".locus-input")
    box.fill("20-60")
    box.press("Enter")
    page.wait_for_timeout(300)

    page.get_by_role("button", name="Save view…").click()
    page.locator(".modal input.input").first.fill("my spot")
    page.get_by_role("button", name="Save", exact=True).click()
    page.locator(".toast", has_text="Saved").wait_for()

    # Move away from everything the session captured.
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S2']); m.setFilters({}); }")
    box.fill("1-10")
    box.press("Enter")
    page.wait_for_timeout(300)

    _open(page, f"{server_url}/#/overview")
    card = page.locator("section.card", has_text="Saved views")
    card.get_by_text("my spot").click()
    page.locator(".toast", has_text="Opened").wait_for()
    _wait_state(page, "m.state.pinned.join() === 'S1'")
    restored = page.evaluate("async () => { const m = await import('/static/js/state.js'); return [m.state.filters, m.state.locus.start, m.state.locus.end, m.state.locus.gene]; }")
    assert restored[0] == {"superpopulation": ["AFR"]}
    assert round(restored[1]) == 19 and round(restored[2]) == 60 and restored[3] == "GENE1"

    _open(page, f"{server_url}/#/overview")  # opening a session navigates to its locus on Tracks
    card = page.locator("section.card", has_text="Saved views")
    page.once("dialog", lambda d: d.accept("renamed spot"))
    card.get_by_role("button", name="Rename").click()
    card.get_by_text("renamed spot").wait_for()
    page.once("dialog", lambda d: d.accept())  # the delete confirm takes no text
    card.locator(".session-item button").last.click()
    card.get_by_text("None yet").wait_for()
    assert page.errors == []


def test_jobs_page_eta_retry_and_notification_toggle(page, server_url):
    """The Jobs page: time left for a running job, Retry on a finished one, and the notify opt-in."""
    import json
    import re

    now = 1_700_000_000.0
    tasks = {"tasks": [
        {"id": "t-run", "kind": "train", "title": "Half done", "status": "running", "progress": 0.25,
         "message": "epoch 25/100", "created": now, "started": now, "finished": None, "elapsed": 60.0,
         "params": {}, "result": {}, "steps": 1, "step": 0, "step_title": "", "last_line": "", "dir": "/tmp/t-run"},
        {"id": "t-done", "kind": "predict", "title": "Finished one", "status": "failed", "progress": 1.0,
         "message": "exit 1", "created": now, "started": now, "finished": now + 5, "elapsed": 5.0,
         "params": {}, "result": {}, "steps": 1, "step": 0, "step_title": "", "last_line": "", "dir": "/tmp/t-done"},
    ], "available": True, "allow": True, "dir": "/tmp/tasks"}
    retried = []

    def handle(route):
        url = route.request.url
        if re.search(r"/api/tasks/[^/]+/retry", url):
            retried.append(url)
            route.fulfill(status=200, content_type="application/json",
                          body=json.dumps({**tasks["tasks"][1], "id": "t-new", "title": "Finished one", "status": "running"}))
        elif re.search(r"/api/tasks/[^/]+$", url):
            route.fulfill(status=200, content_type="application/json", body=json.dumps({**tasks["tasks"][0], "log": ["line"], "commands": []}))
        elif url.endswith("/api/tasks"):
            route.fulfill(status=200, content_type="application/json", body=json.dumps(tasks))
        else:
            route.continue_()

    page.route(re.compile(r".*/api/tasks.*"), handle)
    _open(page, f"{server_url}/#/jobs")
    rows = page.locator("table.jobs-table tbody tr")
    rows.first.wait_for()
    # 25% in 60 s implies about 3 minutes left.
    assert "~3m 0s left" in rows.nth(0).inner_text()
    assert "left" not in rows.nth(1).inner_text()  # a finished job has no estimate

    page.locator("label.check", has_text="Notify me when a job finishes").locator("input").wait_for()
    rows.nth(1).get_by_role("button", name="Retry").click()
    page.locator(".toast", has_text="Started again").wait_for()
    assert len(retried) == 1 and "t-done" in retried[0]
    assert rows.nth(0).get_by_role("button", name="Retry").count() == 0  # not while it runs
    assert page.errors == []


def test_tracks_copy_as_python_runs_and_returns_the_view(page, server_url):
    """Export → Copy as Python writes a client call for what is on screen — and it runs."""
    _open(page, f"{server_url}/#/samples")
    page.evaluate("async () => { const m = await import('/static/js/state.js'); m.setPinned(['S1', 'S2']); }")
    _open(page, f"{server_url}/#/tracks?gene=GENE1")
    page.wait_for_function("document.querySelectorAll('.panel canvas').length >= 1")
    box = page.locator(".locus-input")
    box.fill("10-50")
    box.press("Enter")
    page.wait_for_timeout(300)

    page.get_by_role("button", name="Export").click()
    # Headless has no clipboard permission, so copyText falls back to a prompt with the text.
    prompts = []
    page.on("dialog", lambda d: (prompts.append(d.default_value), d.accept()))
    page.get_by_role("menuitem", name="Copy as Python").click()
    page.wait_for_timeout(500)
    assert prompts, "no snippet offered"
    snippet = prompts[0]
    assert "from genomics.visualizer.client import Visualizer" in snippet
    assert "v.signal('GENE1', 'rna_seq'" in snippet and "series=['S1:H1+H2','S2:H1+H2']" in snippet
    assert "start=9, end=50" in snippet and f"Visualizer('{server_url}'" in snippet

    # The point of the snippet is that it runs: execute it and check it returns this view's arrays.
    namespace: dict = {}
    exec(snippet, namespace)  # noqa: S102 - text generated by the app under test
    data = namespace["data"]
    assert [s["label"] for s in data["series"]] == ["S1 H1+H2", "S2 H1+H2"]
    assert data["edges"][0] == 9 and data["edges"][-1] == 50
    assert data["series"][0]["mean"].shape[0] == 2  # both chosen tracks
    assert namespace["reference"]["series"][0]["label"] == "Reference genome"
    assert page.errors == []


def test_import_page_without_an_alphagenome_backend(page, server_url):
    """No AlphaGenome backend (a CI runner, a fresh install): the Import page still works.

    The page prefetches the track catalog to fill the tissue picker; when the server has no backend
    it answers 503, which must stay a missing feature rather than a broken page.
    """
    import json
    import re

    seen = []

    def handle(route):
        seen.append(route.request.url)
        route.fulfill(status=503, content_type="application/json",
                      body=json.dumps({"error": "AlphaGenome backend not ready: no API key"}))

    page.route(re.compile(r".*/api/alphagenome/catalog.*"), handle)
    _open(page, f"{server_url}/#/import")
    assert seen, "the page did not ask for the catalog"
    assert page.locator("#page").inner_text().strip()
    boxes = [t for t in page.locator(".error-box").all_inner_texts() if not t.startswith(EXPECTED_BOXES)]
    assert boxes == [], boxes  # a machine without bcftools may still show that notice
    assert page.locator(".toast.error").count() == 0
    # Only the browser's own line for the refused request, which IGNORED_CONSOLE covers.
    assert [e for e in page.errors if "status of 503" not in e] == []
