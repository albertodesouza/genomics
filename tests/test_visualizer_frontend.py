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
# Resource loads the app does not control (the browser asks for a favicon).
IGNORED_CONSOLE = ("favicon",)
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


@pytest.mark.parametrize("name", PAGES)
def test_page_renders_without_errors(page, server_url, name):
    _open(page, f"{server_url}/#/{name}")
    assert page.locator("#page").inner_text().strip(), f"{name} rendered nothing"
    assert page.locator(".error-box").count() == 0, f"{name}: {page.locator('.error-box').all_inner_texts()}"
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
