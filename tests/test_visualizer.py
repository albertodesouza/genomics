"""Tests for the unified genomics visualizer (genomics.visualizer).

Uses a small synthetic canonical-layout dataset whose haplotype sequences are built by applying
each sample's phased VCF to the reference the same way ``bcftools consensus`` does, so the
reference <-> haplotype coordinate mapping can be checked exactly without bcftools installed.
"""
import gzip
import http.client
import json
import socket
import struct
import threading
import time
import zlib
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlencode

import numpy as np
import pytest

from genomics.visualizer.cache import DiskArrayCache, LRUCache
from genomics.visualizer.coords import inflate_gzip, parse_window_vcf, reference_map
from genomics.visualizer.datasets import Dataset, DatasetCatalog
from genomics.visualizer.experiments import ExperimentService
from genomics.visualizer.jobs import JobCancelled, JobManager
from genomics.visualizer.sequences import SequenceService
from genomics.visualizer.signals import SeriesSpec, SignalService, bin_matrix, load_prediction_matrix

WINDOW_START = 1001
REF = ("ACGTTGCAAC" * 30)[:300]
L = len(REF)
TRACKS = 2

# (pos_1based, ref, alts, {sample: gt}) — includes SNVs, an insertion, a deletion, an overlapping
# record bcftools skips, an insertion anchored on a deletion's last base (applied), and a symbolic
# <DEL> carrying INFO/END.
def _r(off: int, n: int = 1) -> str:
    return REF[off:off + n]


VARIANTS = [
    (1011, _r(10), ["T"], {"S1": "1|0", "S2": "0|1", "S3": "0|0"}),
    (1031, _r(30), [_r(30) + "AAA"], {"S1": "1|1", "S2": "0|0", "S3": "1|0"}),
    (1051, _r(50, 4), [_r(50)], {"S1": "0|1", "S2": "1|1", "S3": "0|0"}),
    (1051, _r(50), [_r(50) + "G"], {"S1": "0|1", "S2": "0|0", "S3": "0|0"}),  # overlaps an applied record on S1 H2: skipped
    (1081, _r(80, 4), [_r(80)], {"S1": "0|0", "S2": "0|0", "S3": "1|1"}),
    (1084, _r(83), [_r(83) + "TT"], {"S1": "0|0", "S2": "0|0", "S3": "1|0"}),  # anchored on the deletion's last base: applied
    (1201, _r(200), ["<DEL>"], {"S1": "0|0", "S2": "1|0", "S3": "0|0"}),  # END=1210
]
SAMPLES = {"S1": ("YRI", "AFR"), "S2": ("CEU", "EUR"), "S3": ("YRI", "AFR")}


def _apply(haplotype_index: int, sample: str):
    """Reference application mirroring bcftools consensus' overlap rule; returns (seq, ref->local)."""
    applied = []
    frozen, last_ins = -1, False
    for pos, ref, alts, gts in VARIANTS:
        allele = int(gts[sample].split("|")[haplotype_index])
        if allele == 0:
            continue
        alt = alts[allele - 1]
        off = pos - WINDOW_START
        rlen, alen = len(ref), len(alt)
        if alt == "<DEL>":
            rlen, alen, alt = 1210 - pos + 1, 1, ref[0]
        if off < frozen or (off == frozen and (rlen == alen or last_ins)):
            continue
        # An indel anchored on the last base of the previous record only adds its new bases.
        applied.append((off, rlen, alt, off == frozen))
        frozen, last_ins = off + rlen - 1, alen > rlen
    out, local, i = [], [-1] * L, 0
    for off, rlen, alt, anchored in applied:
        while i < off:
            local[i] = len(out)
            out.append(REF[i])
            i += 1
        if anchored:
            out.extend(alt[rlen:])
            i = max(i, off + rlen)
            continue
        anchor = len(out)
        out.extend(alt)
        for k in range(min(rlen, len(alt))):  # shared prefix keeps its reference bases' positions
            local[off + k] = anchor + k
        i = off + rlen
    while i < L:
        local[i] = len(out)
        out.append(REF[i])
        i += 1
    return "".join(out), local


def _bgzf(data: bytes, block: int = 37) -> bytes:
    """Minimal BGZF writer (multi-member gzip with BC extra field)."""
    out = b""
    for i in range(0, len(data), block):
        chunk = data[i:i + block]
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = comp.compress(chunk) + comp.flush()
        bsize = 12 + 6 + len(payload) + 8 - 1
        header = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, bsize)
        out += header + payload + struct.pack("<II", zlib.crc32(chunk) & 0xFFFFFFFF, len(chunk))
    eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    return out + eof


def _write_vcf(path: Path, sample: str) -> None:
    lines = ["##fileformat=VCFv4.2", f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"]
    for pos, ref, alts, gts in VARIANTS:
        info = "END=1210;SVTYPE=DEL" if alts == ["<DEL>"] else "."
        lines.append(f"chr1\t{pos}\tv{pos}\t{ref}\t{','.join(alts)}\t.\tPASS\t{info}\tGT\t{gts[sample]}")
    path.write_bytes(_bgzf(("\n".join(lines) + "\n").encode()))


@pytest.fixture(scope="module")
def dataset_dir(tmp_path_factory):
    root = tmp_path_factory.mktemp("viz") / "dataset"
    meta = {
        "dataset_name": "synthetic",
        "individuals": list(SAMPLES),
        "individuals_pedigree": {s: {"population": p, "superpopulation": sp, "sex": 1, "sex_label": "Male"} for s, (p, sp) in SAMPLES.items()},
        "genes": ["GENE1"],
        "alphagenome_outputs": ["RNA_SEQ"],
        "window_size": L,
    }
    root.mkdir(parents=True)
    (root / "dataset_metadata.json").write_text(json.dumps(meta))
    ref_dir = root / "references" / "windows" / "GENE1"
    ref_dir.mkdir(parents=True)
    # window_metadata deliberately disagrees with the FASTA header: the header must win.
    (ref_dir / "ref.window.fa").write_text(f">chr1:{WINDOW_START}-{WINDOW_START + L - 1}\n{REF}\n")
    (ref_dir / "window_metadata.json").write_text(json.dumps({"chromosome": "chr1", "start": WINDOW_START - 2, "end": WINDOW_START + L - 1}))
    (root / "references" / "windows" / "EXTRA").mkdir(parents=True)  # on disk only
    (root / "references" / "windows" / "EXTRA" / "ref.window.fa").write_text(">chr2:1-10\nACGTACGTAC\n")
    for sample in SAMPLES:
        gene_dir = root / "individuals" / sample / "windows" / "GENE1"
        gene_dir.mkdir(parents=True)
        _write_vcf(gene_dir / f"{sample}.window.consensus_ready.vcf.gz", sample)
        for h_index, hap in enumerate(("H1", "H2")):
            seq, _local = _apply(h_index, sample)
            fixed = (seq + REF[len(seq):])[:L] if len(seq) < L else seq[:L]
            (gene_dir / f"{sample}.{hap}.window.fixed.fa").write_text(f">{sample}\n{fixed}\n")
            pred = gene_dir / f"predictions_{hap}"
            pred.mkdir()
            # track 0 = haplotype index (so remapping is observable), track 1 = per-sample constant
            values = np.stack([np.arange(L, dtype=np.float32), np.full(L, float(int(sample[1])), np.float32)], axis=1)
            np.savez_compressed(pred / "rna_seq.npz", values=values)
            (pred / "rna_seq_metadata.json").write_text(json.dumps({"metadata": [{"biosample_name": "cell A", "strand": "+"}, {"biosample_name": "cell B", "strand": "-"}]}))
    return root


@pytest.fixture
def dataset(dataset_dir):
    return Dataset(dataset_dir)


def test_inflate_gzip_handles_bgzf_and_plain_gzip():
    data = b"".join(f"line {i}\tvalue\n".encode() for i in range(500))
    assert inflate_gzip(_bgzf(data)) == data
    assert inflate_gzip(gzip.compress(data)) == data


def test_reference_map_matches_bcftools_style_application(dataset_dir):
    for sample in SAMPLES:
        parsed = parse_window_vcf(dataset_dir / "individuals" / sample / "windows" / "GENE1" / f"{sample}.window.consensus_ready.vcf.gz", WINDOW_START)
        for h_index, hap in enumerate(("H1", "H2")):
            seq, expected_local = _apply(h_index, sample)
            ref_map = reference_map(parsed.events(hap), L)
            assert ref_map.local.tolist() == expected_local, (sample, hap)
            kept = [off for off in range(L) if expected_local[off] >= 0]
            assert all(seq[expected_local[off]] in "ACGT" for off in kept)
    s1 = parse_window_vcf(dataset_dir / "individuals" / "S1" / "windows" / "GENE1" / "S1.window.consensus_ready.vcf.gz", WINDOW_START)
    h2 = s1.events("H2")
    assert list(zip(h2.offsets.tolist(), h2.ref_lengths.tolist(), h2.alt_lengths.tolist())) == [(30, 1, 4), (50, 4, 1)]  # SNV is on H1; overlapping insertion skipped
    s3 = parse_window_vcf(dataset_dir / "individuals" / "S3" / "windows" / "GENE1" / "S3.window.consensus_ready.vcf.gz", WINDOW_START)
    assert (83, 1, 3) in list(zip(s3.events("H1").offsets.tolist(), s3.events("H1").ref_lengths.tolist(), s3.events("H1").alt_lengths.tolist()))
    s2 = reference_map(parse_window_vcf(dataset_dir / "individuals" / "S2" / "windows" / "GENE1" / "S2.window.consensus_ready.vcf.gz", WINDOW_START).events("H1"), L)
    assert (s2.local[201:210] < 0).all() and s2.local[200] >= 0 and s2.local[210] >= 0  # symbolic <DEL> applied


def test_load_prediction_matrix_fast_path_and_fallbacks(tmp_path):
    values = np.random.default_rng(0).random((50, 3)).astype(np.float32)
    for saver, name in ((np.savez_compressed, "a.npz"), (np.savez, "b.npz")):
        saver(tmp_path / name, values=values)
        assert np.array_equal(load_prediction_matrix(tmp_path / name), values)
    np.savez(tmp_path / "c.npz", track_0=values[:, 0], track_1=values[:, 1])
    assert np.array_equal(load_prediction_matrix(tmp_path / "c.npz"), values[:, :2])
    np.savez(tmp_path / "d.npz", signal=values[:, 0])
    assert load_prediction_matrix(tmp_path / "d.npz").shape == (50, 1)


def test_bin_matrix_nan_aware_mean_min_max():
    values = np.array([[1, 10], [3, np.nan], [np.nan, np.nan], [5, 20]], dtype=np.float32)
    binned = bin_matrix(values, 2)
    assert binned["edges"].tolist() == [0, 2, 4]
    np.testing.assert_allclose(binned["mean"], [[2, 5], [10, 20]])
    np.testing.assert_allclose(binned["min"], [[1, 5], [10, 20]])
    np.testing.assert_allclose(binned["max"], [[3, 5], [10, 20]])
    small = bin_matrix(values, 10)
    assert small["edges"].tolist() == [0, 1, 2, 3, 4]


def test_dataset_discovery_facets_and_windows(dataset):
    assert dataset.genes == ["EXTRA", "GENE1"]
    window = dataset.window("GENE1")
    assert (window.chromosome, window.start, window.length) == ("chr1", WINDOW_START, L)
    fields = {f["name"]: f for f in dataset.fields}
    assert fields["superpopulation"]["kind"] == "categorical"
    assert "sex" not in fields  # numeric duplicate of sex_label
    assert {row["sample_id"]: row.get("pigmentation") for row in dataset.samples} == {"S1": "strong pigmentation", "S2": "weak pigmentation", "S3": "strong pigmentation"}
    assert dataset.filter_samples({"superpopulation": ["AFR"]}) == ["S1", "S3"]
    info = dataset.gene_info("GENE1")
    assert info["haplotypes"] == ["H1", "H2"]
    assert [t["short"] for t in info["outputs"]["rna_seq"]["tracks"]] == ["cell A (+)", "cell B (-)"]


def test_annotation_table_adds_facets(dataset_dir, tmp_path):
    table = tmp_path / "pheno.tsv"
    table.write_text("sample\tphenotype\nS1\tcase\nS2\tcontrol\nS3\tcase\n")
    ds = Dataset(dataset_dir, annotations=table)
    assert ds.filter_samples({"phenotype": ["case"]}) == ["S1", "S3"]


def _service(tmp_path=None):
    return SignalService(256 << 20, DiskArrayCache(tmp_path), workers=2)


def test_series_reference_coordinates_remap_through_indels(dataset):
    signals = _service()
    payload = signals.series_payload(dataset, "GENE1", "rna_seq", [SeriesSpec("S2", "H1"), SeriesSpec("S2", "H1+H2")], [0, 1], "reference", 0, L, L)
    s2h1 = payload["series"][0]["mean"][0]
    ref_map = signals.ref_map(dataset, "S2", "GENE1", "H1")
    expected = np.where((ref_map.local >= 0) & (ref_map.local < L), ref_map.local, np.nan).astype(np.float32)
    np.testing.assert_array_equal(s2h1, expected)
    assert np.isnan(s2h1[51:54]).all()  # TTGC>T deletes three reference bases on S2 H1
    diploid = payload["series"][1]["mean"][1]
    assert np.isnan(diploid).sum() == 3 and np.nanmin(diploid) == np.nanmax(diploid) == 2.0  # gap only where both haplotypes lack the base
    raw = signals.series_payload(dataset, "GENE1", "rna_seq", [SeriesSpec("S2", "H1")], [0], "haplotype", 0, L, L)
    np.testing.assert_array_equal(raw["series"][0]["mean"][0], np.arange(L, dtype=np.float32))


def test_group_aggregate_matches_numpy_and_is_cached_on_disk(dataset, tmp_path):
    signals = _service(tmp_path)
    progress = lambda *_: None  # noqa: E731
    agg = signals.compute_group_aggregate(dataset, "GENE1", "rna_seq", [1], ["H1", "H2"], "haplotype", ["S1", "S2", "S3"], progress)
    np.testing.assert_allclose(agg["mean"][:, 0], np.full(L, 2.0))
    np.testing.assert_allclose(agg["std"][:, 0], np.full(L, np.std([1, 1, 2, 2, 3, 3])), rtol=1e-5)
    assert int(agg["n_haplotypes"][0]) == 6
    fresh = _service(tmp_path)
    key = fresh.group_aggregate_key(dataset, "GENE1", "rna_seq", [1], ["H1", "H2"], "haplotype", ["S1", "S2", "S3"])
    assert fresh.cached_group_aggregate(key) is not None
    payload = fresh.group_payload({"all": agg}, {"all": 3}, [1], L, 0, L, 10)
    assert len(payload["edges"]) == 11


def test_population_matrix_rows_follow_samples(dataset):
    signals = _service()
    pop = signals.population_matrix(dataset, "GENE1", "rna_seq", 1, "H1+H2", "reference", ["S3", "S1"], 0, L, 30, lambda *_: None)
    assert pop["matrix"].shape == (2, 30)
    np.testing.assert_allclose(np.nanmean(pop["matrix"], axis=1), [3.0, 1.0])


def test_sequence_letters_show_snvs_deletions_and_insertions(dataset):
    seqs = SequenceService(_service())
    res = seqs.window(dataset, "GENE1", [SeriesSpec("S1", "H1+H2")], "reference", 0, 80, 100)
    assert res["mode"] == "letters" and res["reference"] == REF[:80]
    h1, h2 = res["rows"]
    assert h1["bases"][10] == "T" and h2["bases"][10] == REF[10]  # SNV only on H1
    assert h2["bases"][51:54] == "---"  # deletion on H2
    assert [i["pos"] for i in h1["insertions"]] == [30] and h1["insertions"][0]["seq"] == "AAA"
    assert {v["genomic"] for v in res["variants"]} >= {1011, 1031, 1051}
    dense = seqs.window(dataset, "GENE1", [SeriesSpec("S1", "H2")], "reference", 0, L, 10)
    assert dense["mode"] == "letters"  # short windows always come back as letters
    assert dense["rows"][0]["deletions"] == 3


def test_composition_counts_letter_frequencies_per_base(dataset):
    seqs = SequenceService(_service())
    letters = list("ACGTN-")
    ref = seqs.composition(dataset, "GENE1", [], "reference", 0, 80, 100)
    assert ref["source"] == "reference" and ref["letters"] == "".join(letters) and ref["reference"] == REF[:80]
    freq = ref["freq"]
    assert freq.shape == (6, 80) and np.allclose(freq.sum(axis=0), 1.0)
    assert all(freq[letters.index(REF[i]), i] == 1.0 for i in range(80))
    both = seqs.composition(dataset, "GENE1", [SeriesSpec("S1", "H1+H2")], "reference", 0, 80, 100)
    assert both["haplotypes"] == ["S1:H1", "S1:H2"]
    assert both["freq"][letters.index("T"), 10] == 0.5 and both["freq"][letters.index(REF[10]), 10] == 0.5  # SNV on H1 only
    assert both["freq"][letters.index("-"), 52] == 0.5  # deletion on H2
    binned = seqs.composition(dataset, "GENE1", [], "reference", 0, L, 10)
    assert binned["freq"].shape == (6, 10) and "reference" in binned and np.allclose(binned["freq"].sum(axis=0), 1.0)


def test_lru_cache_respects_budget_and_dedupes_loads():
    cache = LRUCache(1000)
    calls = []
    for i in range(5):
        cache.get_or_load(i, lambda i=i: calls.append(i) or np.zeros(50, np.float32))
    assert cache.stats()["bytes"] <= 1000
    cache.get_or_load(4, lambda: calls.append("again"))
    assert "again" not in calls


def test_job_cancellation_stops_work():
    jobs = JobManager(workers=1)
    started = threading.Event()

    def work(progress):
        started.set()
        for i in range(200):
            time.sleep(0.01)
            progress(i / 200, "working")
        return "finished"

    job = jobs.run("k", "slow", work)
    started.wait(2)
    jobs.cancel(job.id)
    for _ in range(100):
        if job.status == "cancelled":
            break
        time.sleep(0.02)
    assert job.status == "cancelled"
    with pytest.raises(JobCancelled):
        from genomics.visualizer.jobs import Progress

        job.cancelled = True
        Progress(job)(0.5, "x")


def test_experiment_service_reads_generic_metrics(tmp_path):
    run = tmp_path / "runs" / "cnn_a"
    (run / "models").mkdir(parents=True)
    (run / "manifest.json").write_text(json.dumps({"status": "completed", "best_val_accuracy": 0.9, "created_at": "2026-01-01T00:00:00"}))
    (run / "val_best_accuracy_results.json").write_text(json.dumps({"weighted_accuracy": 0.91, "confusion_matrix": [[3, 1], [0, 4]], "per_class_metrics": {"a": {"f1": 0.8}, "b": {"f1": 0.9}}}))
    (run / "models" / "training_history.json").write_text(json.dumps({"epoch": [1, 2], "val_loss": [float("inf"), 0.5], "interrupted": False}))
    (run / "config.yaml").write_text("model: cnn\n")
    service = ExperimentService([tmp_path / "runs"])
    listing = service.list()
    assert listing["runs"][0]["metrics"]["val_best_accuracy.weighted_accuracy"] == 0.91
    detail = service.detail("cnn_a")
    assert detail["history"]["val_loss"] == [None, 0.5]
    assert detail["results_detail"][0]["class_names"] == ["a", "b"]
    with pytest.raises(KeyError):
        service.file("cnn_a", "../../etc/passwd")


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_http_api_end_to_end(dataset_dir, tmp_path):
    from genomics.visualizer.server import Handler, Server, VisualizerApp

    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=256 << 20, workers=2, runs_roots=[])
    handler = type("TestHandler", (Handler,), {"app": app})
    port = _free_port()
    server = Server(("127.0.0.1", port), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def get(path, method="GET", body=None):
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
        conn.request(method, path, body=json.dumps(body) if body is not None else None, headers={"Content-Type": "application/json"})
        res = conn.getresponse()
        data = res.read()
        conn.close()
        return res.status, res.getheader("Content-Type"), data

    try:
        status, ctype, body = get("/")
        assert status == 200 and "text/html" in ctype and b"/static/js/app.js" in body
        assert get("/static/js/pages/tracks.js")[0] == 200
        assert get("/samples")[0] == 200  # client-side routes fall back to index.html
        status, _, body = get("/api/status")
        assert json.loads(body)["datasets"][0]["id"] == ds.id
        base = f"/api/d/{ds.id}"
        assert json.loads(get(f"{base}/summary")[2])["sample_count"] == 3
        params = urlencode({"gene": "GENE1", "output": "rna_seq", "series": "S1:H1,S2:H2", "tracks": "0,1", "coords": "reference", "start": 0, "end": L, "bins": 50})
        status, _, body = get(f"{base}/signal?{params}")
        payload = json.loads(body)
        assert status == 200 and payload["series"][0]["mean"]["shape"] == [2, 50]
        decoded = np.frombuffer(__import__("base64").b64decode(payload["series"][1]["mean"]["$f32"]), dtype="<f4").reshape(2, 50)
        np.testing.assert_allclose(decoded[1], 2.0)
        params = urlencode({"gene": "GENE1", "output": "rna_seq", "tracks": "1", "coords": "reference", "field": "superpopulation", "start": 0, "end": L, "bins": 20})
        for _ in range(100):
            data = json.loads(get(f"{base}/groups?{params}")[2])
            if not data.get("pending"):
                break
            time.sleep(0.05)
        assert {g["group"]: g["samples"] for g in data["groups"]} == {"AFR": 2, "EUR": 1}
        assert get(f"{base}/signal?gene=NOPE&output=rna_seq&series=S1:H1")[0] == 404
        status, _, body = get(f"{base}/views/preview", "POST", {"name": "my view", "sample_ids": ["S1"], "genes": ["GENE1"]})
        assert status == 200 and json.loads(body)["view"]["sample_ids"] == ["S1"]
    finally:
        server.shutdown()
        server.server_close()
        app.shutdown()


def test_cli_exposes_visualize_and_workbench_alias(monkeypatch):
    from genomics import cli

    calls = []
    monkeypatch.setattr(cli, "_run_module", lambda module, args: calls.append((module, [str(a) for a in args])) or 0)
    assert cli.main(["visualize", "--dataset", "/tmp/ds", "--port", "9000", "--open"]) == 0
    module, args = calls[-1]
    assert module == "genomics.visualizer" and "--open" in args and args[args.index("--port") + 1] == "9000"
    assert cli.main(["genotype", "workbench", "--dataset-dir", "/tmp/ds"]) == 0
    assert calls[-1][0] == "genomics.visualizer"
    assert cli.main(["genotype", "workbench", "--legacy"]) == 0
    assert calls[-1][0] == "genomics.predictors.genotype_based.apps.genomics_workbench"


# ------------------------------------------------------------------------- reference & links
def _write_reference_prediction(dataset_dir):
    """Reference-window prediction with the dataset's 'cell B (-)' track plus a track it lacks, in another order."""
    pred = dataset_dir / "references" / "windows" / "GENE1" / "predictions_ref"
    pred.mkdir(exist_ok=True)
    values = np.stack([np.full(L, 7.0, np.float32), np.arange(L, dtype=np.float32) * 10], axis=1)
    np.savez_compressed(pred / "rna_seq.npz", values=values)
    (pred / "rna_seq_metadata.json").write_text(json.dumps({"metadata": [{"biosample_name": "cell B", "strand": "-"}, {"biosample_name": "cell X", "strand": "+"}]}))
    return pred


def test_reference_prediction_series_matches_tracks_and_coordinates(dataset_dir, tmp_path):
    import shutil

    from genomics.visualizer.signals import REFERENCE_SAMPLE, match_track_columns

    pred = _write_reference_prediction(dataset_dir)
    try:
        dataset = Dataset(dataset_dir)
        assert dataset.gene_info("GENE1")["reference_outputs"] == ["rna_seq"]
        signals = _service(tmp_path)
        for coords in ("reference", "haplotype"):
            payload = signals.series_payload(dataset, "GENE1", "rna_seq", [SeriesSpec(REFERENCE_SAMPLE, "H1"), SeriesSpec("S1", "H1")], [0, 1], coords, 0, L, L)
            ref = payload["series"][0]
            assert ref["label"] == "Reference genome"
            assert np.isnan(ref["mean"][0]).all()  # 'cell A (+)' was not predicted on the reference
            np.testing.assert_array_equal(ref["mean"][1], np.full(L, 7.0, np.float32))
        # Training-axis coordinates: insertion slots are gaps, the rest follow the reference offsets.
        signals.alignment = SimpleNamespace(axis=lambda ds, gene: {"insertion_slots": np.array([2, 5]), "ref_start_offset": 1, "expanded_length": 10})
        mapped = signals.reference_indexed_window(dataset, "GENE1", np.arange(L, dtype=np.float32).reshape(-1, 1), 1, "aligned", 0, 8)[:, 0]
        np.testing.assert_array_equal(mapped, [1, 2, np.nan, 3, 4, np.nan, 5, 6])
        binned = signals.reference_indexed_window(dataset, "GENE1", np.arange(3, dtype=np.float32).reshape(-1, 1), 128, "reference", 120, 140)[:, 0]
        np.testing.assert_array_equal(binned, [0] * 8 + [1] * 12)  # 128 bp rows expand to bases
    finally:
        shutil.rmtree(pred)
    assert match_track_columns([{"ontology_curie": "CL:1", "strand": "+"}], [{"ontology_curie": "CL:1", "strand": "-"}, {"ontology_curie": "CL:1", "strand": "+"}]) == [-1, 0]


def test_gene_models_use_0_based_gtf_starts():
    pd = pytest.importorskip("pandas")
    from genomics.visualizer.annotations import _features_to_models

    # pyranges tables are 0-based half-open; window offset 0 is the 1-based position 1001.
    df = pd.DataFrame([
        {"Feature": "gene", "Start": 1000, "End": 1010, "Strand": "+", "gene_name": "G", "gene_type": "protein_coding", "transcript_id": None, "transcript_name": None, "transcript_type": None, "tag": None, "exon_number": None},
        {"Feature": "exon", "Start": 1002, "End": 1005, "Strand": "+", "gene_name": "G", "gene_type": "protein_coding", "transcript_id": "T1", "transcript_name": "G-201", "transcript_type": "protein_coding", "tag": "", "exon_number": 1},
    ])
    gene = _features_to_models(df, 1001)[0]
    assert (gene["start"], gene["end"]) == (0, 10)
    assert gene["transcripts"][0]["exons"] == [[2, 5]]


def test_http_api_reference_observed_and_offline_lookups(dataset_dir, tmp_path):
    import shutil

    from genomics.visualizer.server import Handler, Server, VisualizerApp

    pred = _write_reference_prediction(dataset_dir)
    catalog = DatasetCatalog()
    ds = catalog.add(dataset_dir)
    app = VisualizerApp(catalog, cache_dir=tmp_path / "cache", memory_bytes=256 << 20, workers=2, runs_roots=[], remote=False)
    handler = type("TestHandler", (Handler,), {"app": app})
    port = _free_port()
    server = Server(("127.0.0.1", port), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()

    def get(path):
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
        conn.request("GET", path)
        res = conn.getresponse()
        data = res.read()
        conn.close()
        return res.status, json.loads(data)

    try:
        base = f"/api/d/{ds.id}"
        status, info = get(f"{base}/genes/GENE1")
        assert status == 200 and info["reference_outputs"] == ["rna_seq"]
        params = urlencode({"gene": "GENE1", "output": "rna_seq", "series": "@reference:H1", "tracks": "1", "start": 0, "end": L, "bins": 10})
        status, payload = get(f"{base}/signal?{params}")
        mean = np.frombuffer(__import__("base64").b64decode(payload["series"][0]["mean"]["$f32"]), dtype="<f4")
        assert status == 200 and np.allclose(mean, 7.0)
        params = urlencode({"gene": "GENE1", "output": "rna_seq", "tracks": "0,1", "start": 0, "end": L, "bins": 10})
        for _ in range(100):
            status, data = get(f"{base}/observed?{params}")
            if not data.get("pending"):
                break
            time.sleep(0.05)
        # The synthetic tracks have no ontology term / assay: reported, not an error.
        assert status == 200 and [t["available"] for t in data["tracks"]] == [False, False] and "ontology term" in data["tracks"][0]["reason"]
        status, err = get("/api/ontology/term?curie=CL:0000001")
        assert status == 502 and "--no-remote" in err["error"]
        status, names = get("/api/genes/names?symbols=GENE1")
        assert status == 200 and names["names"] == {}
    finally:
        server.shutdown()
        server.server_close()
        app.shutdown()
        shutil.rmtree(pred)
