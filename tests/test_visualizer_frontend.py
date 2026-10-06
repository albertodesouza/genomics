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


def _open(page, url):
    page.goto(url)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(300)


@pytest.mark.parametrize("name", PAGES)
def test_page_renders_without_errors(page, server_url, name):
    _open(page, f"{server_url}/#/{name}")
    assert page.locator("#page").inner_text().strip(), f"{name} rendered nothing"
    assert page.locator(".error-box").count() == 0, f"{name}: {page.locator('.error-box').all_inner_texts()}"
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
