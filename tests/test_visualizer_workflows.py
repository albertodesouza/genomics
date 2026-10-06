"""Tests for the visualizer's background jobs, VCF import, AlphaGenome dataset predictions,
training launcher and the in-process Perturbation Lab.

Network and GPU are never used: AlphaGenome calls are faked, and the bcftools-based import test is
skipped when bcftools/samtools are not installed.
"""
import gzip
import json
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from genomics.visualizer.cache import DiskArrayCache
from genomics.visualizer.datasets import Dataset, DatasetCatalog, DatasetMemory
from genomics.visualizer.perturb import PerturbError, PerturbService, reorder_columns
from genomics.visualizer.signals import SignalService
from genomics.visualizer.tasks import TaskManager
from genomics.workflows.dataset_builders.vcf_import import metadata as meta
from genomics.workflows.dataset_builders.vcf_import.builder import centred_window

REF = ("ACGTTGCAAC" * 20)[:200]
START = 1001  # genomic position of window offset 0


# ----------------------------------------------------------------------------------- metadata
def test_parse_table_formats_and_id_detection():
    samples = ["HG1", "HG2", "HG3"]
    cols, rows = meta.parse_table("Subject;Eye colour;Age\nHG1;blue;34\nHG2;brown;\n", "x.csv")
    assert cols == ["Subject", "Eye_colour", "Age"] and len(rows) == 2
    assert meta.guess_id_column(cols, rows, samples) == "Subject"
    records = meta.records_by_sample(rows, "Subject")
    assert records == {"HG1": {"Eye_colour": "blue", "Age": 34}, "HG2": {"Eye_colour": "brown"}}

    cols, rows = meta.parse_table(json.dumps({"HG1": {"group": "a"}, "HG3": {"group": "b"}}), "m.json")
    assert meta.guess_id_column(cols, rows, samples) == "sample_id"
    cols, rows = meta.parse_table(json.dumps([{"name": "HG2", "x": 1.5}]))
    assert rows[0]["x"] == 1.5 and meta.guess_id_column(cols, rows, samples) == "name"

    cols, rows = meta.parse_table("F1 HG1 0 0 1 2\nF1 HG2 0 0 2 1\n", "cohort.fam")
    records = meta.records_by_sample(rows, "sample_id", family_column="family_id", sex_column="sex")
    assert records["HG1"]["family_id"] == "F1" and records["HG1"]["sex_label"] == "Male" and records["HG2"]["sex"] == 2

    fields = {f["name"]: f for f in meta.describe_fields({"a": {"g": "x", "n": 1.5}, "b": {"g": "y", "n": 2.5}, "c": {"g": "x"}})}
    assert fields["g"]["kind"] == "categorical" and fields["g"]["counts"][0] == ("x", 2) and fields["n"]["missing"] == 1


def test_centred_window_matches_alphagenome_resize():
    genome = pytest.importorskip("alphagenome.data.genome")
    for start, end, strand in [(100, 201, "+"), (100, 201, "-"), (100, 200, "+"), (100, 200, "-"), (5000, 5001, "+")]:
        interval = genome.Interval("chr1", start, end, strand).resize(524288)
        assert centred_window(start, end, strand, 524288) == (interval.start, interval.end)


# ---------------------------------------------------------------------------------------- tasks
def _wait(tm, task_id, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        task = tm.get(task_id)
        if task["status"] in ("done", "failed", "cancelled", "lost"):
            return task
        time.sleep(0.1)
    raise AssertionError(f"task {task_id} still {tm.get(task_id)['status']}")


def test_task_runner_progress_results_failures_and_hooks(tmp_path):
    tm = TaskManager(tmp_path / "tasks")
    finished = []
    tm.on_finished("demo", finished.append)
    py = sys.executable
    step = "import sys\nprint('@@progress 0.5 halfway', flush=True)\nprint('hello\\rbar 100%')\nprint('@@result {\"n\": 3}')"
    task = tm.create("demo", "Demo", [{"title": "one", "command": [py, "-c", step], "weight": 3},
                                      {"title": "two", "command": [py, "-c", "import sys; open(sys.argv[1], 'w').write('x')", "{task_dir}/out.txt"]}],
                     params={"k": "v"}, files={"in.txt": "data"})
    done = _wait(tm, task["id"])
    assert done["status"] == "done" and done["progress"] == 1.0 and done["result"] == {"n": 3}
    assert (Path(done["dir"]) / "out.txt").read_text() == "x" and (Path(done["dir"]) / "in.txt").read_text() == "data"
    log = tm.log_tail(task["id"])
    assert "bar 100%" in log and not any(line.startswith("@@") for line in log)
    tm.list()
    assert [t["id"] for t in finished] == [task["id"]]
    # A new manager (restarted server) does not repeat the hook.
    again = TaskManager(tmp_path / "tasks")
    repeated = []
    again.on_finished("demo", repeated.append)
    again.list()
    assert repeated == []

    failed = _wait(tm, tm.create("demo", "Fail", [{"title": "f", "command": [py, "-c", "print('boom'); raise SystemExit(3)"]}])["id"])
    assert failed["status"] == "failed" and failed["exit_code"] == 3 and "boom" in failed["message"]
    with pytest.raises(RuntimeError):
        tm.delete(tm.create("demo", "Slow", [{"title": "s", "command": [py, "-c", "import time; time.sleep(30)"]}])["id"])


def test_task_cancel_and_fifo_resource_queue(tmp_path):
    tm = TaskManager(tmp_path / "tasks")
    py = sys.executable
    first = tm.create("demo", "A", [{"title": "a", "command": [py, "-c", "import time; time.sleep(30)"]}], resource="gpu")
    time.sleep(0.2)
    second = tm.create("demo", "B", [{"title": "b", "command": [py, "-c", "print('ran')"]}], resource="gpu")
    deadline = time.time() + 10
    while tm.get(second["id"])["status"] != "queued" and time.time() < deadline:
        time.sleep(0.1)
    assert tm.get(first["id"])["status"] == "running" and tm.get(second["id"])["status"] == "queued"
    tm.cancel(first["id"])
    assert _wait(tm, first["id"])["status"] == "cancelled"
    assert _wait(tm, second["id"])["status"] == "done"
    tm.delete(first["id"])
    assert [t["id"] for t in tm.list()] == [second["id"]]


def test_dataset_memory_roundtrip(tmp_path):
    memory = DatasetMemory(tmp_path / "datasets.json")
    memory.remember(tmp_path / "a", None)
    memory.remember(tmp_path / "b", tmp_path / "ann.tsv")
    memory.remember(tmp_path / "a", None)
    assert [e["path"] for e in memory.entries()] == [str(tmp_path / "b"), str(tmp_path / "a")]
    memory.forget(tmp_path / "b")
    assert [e["path"] for e in memory.entries()] == [str(tmp_path / "a")]
    # Directories that disappeared (e.g. a cleaned-up temp dir) are reported so the UI can forget them.
    assert memory.missing() == [str(tmp_path / "a")]
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "dataset_metadata.json").write_text("{}")
    assert memory.missing() == []


# ------------------------------------------------------------------------------ synthetic dataset
def _haplotypes():
    """S1 H1: SNV at offset 10 (A->T... whatever REF has) + 3 bp insertion after offset 30.
    S1 H2: 3 bp deletion after offset 50."""
    h1 = REF[:10] + ("T" if REF[10] != "T" else "G") + REF[11:31] + "GGG" + REF[31:]
    h2 = REF[:51] + REF[54:]
    return h1, h2


@pytest.fixture
def lab_dataset(tmp_path):
    root = tmp_path / "ds"
    ref_dir = root / "references" / "windows" / "G1"
    ref_dir.mkdir(parents=True)
    (ref_dir / "ref.window.fa").write_text(f">chr1:{START}-{START + len(REF) - 1}\n{REF}\n")
    (ref_dir / "window_metadata.json").write_text(json.dumps({"chromosome": "chr1", "start": START, "end": START + len(REF) - 1}))
    snv_alt = "T" if REF[10] != "T" else "G"
    vcf = "\n".join([
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        f"chr1\t{START + 10}\t.\t{REF[10]}\t{snv_alt}\t.\tPASS\t.\tGT\t1|0",
        f"chr1\t{START + 30}\t.\t{REF[30]}\t{REF[30]}GGG\t.\tPASS\t.\tGT\t1|0",
        f"chr1\t{START + 50}\t.\t{REF[50:54]}\t{REF[50]}\t.\tPASS\t.\tGT\t0|1",
    ]) + "\n"
    gene_dir = root / "individuals" / "S1" / "windows" / "G1"
    gene_dir.mkdir(parents=True)
    (gene_dir / "S1.window.consensus_ready.vcf.gz").write_bytes(gzip.compress(vcf.encode()))
    for hap, seq in zip(("H1", "H2"), _haplotypes()):
        fixed = (seq + REF[len(seq):])[: len(REF)]
        (gene_dir / f"S1.{hap}.window.fixed.fa").write_text(f">S1_{hap}\n{fixed}\n")
        pred = gene_dir / f"predictions_{hap}"
        pred.mkdir()
        values = np.stack([np.arange(len(REF), dtype=np.float32), np.ones(len(REF), np.float32)], axis=1)
        np.savez_compressed(pred / "rna_seq.npz", values=values)
        (pred / "rna_seq_metadata.json").write_text(json.dumps({"metadata": [
            {"ontology_curie": "CL:1", "strand": "+", "biosample_name": "cell"}, {"ontology_curie": "CL:1", "strand": "-", "biosample_name": "cell"}]}))
    (root / "dataset_metadata.json").write_text(json.dumps({
        "dataset_name": "lab", "individuals": ["S1"], "individuals_pedigree": {"S1": {"group": "a"}}, "genes": ["G1"], "window_size": len(REF),
    }))
    return Dataset(root)


def _service(tmp_path):
    signals = SignalService(64 << 20, DiskArrayCache(tmp_path / "cache"))
    return PerturbService(SimpleNamespace(signals=signals))


def test_edits_map_reference_offsets_onto_each_haplotype(lab_dataset, tmp_path):
    service = _service(tmp_path)
    h1, h2 = _haplotypes()
    # Offsets 60-70 lie after H1's 3 bp insertion (+3) and after H2's 3 bp deletion (-3).
    edits = service.normalize_edits([{"op": "overwrite", "start": 60, "end": 70, "base": "a", "haplotypes": ["H1", "H2"]}])
    seq1, done1 = service.edited_haplotype(lab_dataset, "S1", "G1", "H1", edits)
    seq2, done2 = service.edited_haplotype(lab_dataset, "S1", "G1", "H2", edits)
    assert (done1[0]["local_start"], done1[0]["local_end"]) == (63, 73)
    assert (done2[0]["local_start"], done2[0]["local_end"]) == (57, 67)
    assert seq1.tobytes().decode()[63:73] == "A" * 10 and seq1.tobytes().decode()[:63] == h1[:63]
    assert seq2.tobytes().decode()[57:67] == "A" * 10

    # Reverting H1 to the reference over the SNV removes it but keeps the inserted bases.
    revert = service.normalize_edits([{"op": "reference", "start": 0, "end": 40, "haplotypes": ["H1"]}])
    seq, done = service.edited_haplotype(lab_dataset, "S1", "G1", "H1", revert)
    assert done[0]["changed"] == 1 and seq.tobytes().decode()[:34] == REF[:31] + "GGG"

    scrambled, _ = service.edited_haplotype(lab_dataset, "S1", "G1", "H2", service.normalize_edits([{"op": "scramble", "start": 0, "end": 40, "seed": 1}]))
    assert sorted(scrambled.tobytes().decode()[:40]) == sorted(h2[:40]) and len(scrambled) == len(REF)

    with pytest.raises(PerturbError):
        service.edited_haplotype(lab_dataset, "S1", "G1", "H1", service.normalize_edits([{"op": "sequence", "start": 0, "end": 5, "sequence": "ACG"}]))
    with pytest.raises(PerturbError):
        service.normalize_edits([{"op": "delete", "start": 0, "end": 5}])
    with pytest.raises(PerturbError):
        service.normalize_edits([{"op": "overwrite", "start": 5, "end": 5, "base": "A"}])


def test_reorder_columns_follows_stored_track_order():
    stored = [{"ontology_curie": "CL:2", "strand": "+", "name": "x"}, {"ontology_curie": "CL:1", "strand": "-", "name": "y"}]
    predicted = [{"ontology_curie": "CL:1", "strand": "-", "name": "y"}, {"ontology_curie": "CL:3", "strand": "+", "name": "z"}, {"ontology_curie": "CL:2", "strand": "+", "name": "x"}]
    values = np.array([[1.0, 3.0, 2.0]], np.float32)
    np.testing.assert_array_equal(reorder_columns(values, predicted, stored), [[2.0, 1.0]])
    with pytest.raises(PerturbError):
        reorder_columns(values, predicted[:1], stored)


# ------------------------------------------------------------------------- AlphaGenome predictions
class _FakeMetadata:
    def __init__(self, records):
        self.records = records

    def to_dict(self, orient="records"):
        return list(self.records)


def test_predict_dataset_writes_canonical_files_and_skips_done(lab_dataset, monkeypatch):
    pytest.importorskip("alphagenome.models.dna_client")
    from genomics.workflows.alphagenome import predict_dataset
    from genomics.workflows.dataset_builders.non_longevous import build_window_and_predict as bwp

    calls = []

    def fake_predict(client_box, make_client, seq, outputs, terms, timeout_s, max_attempts):
        calls.append((len(seq), [str(o) for o in outputs], terms))
        track = SimpleNamespace(values=np.full((len(seq), 1), len(calls), np.float32), metadata=_FakeMetadata([{"ontology_curie": "CL:9", "strand": ".", "biosample_name": "cell nine", "biosample_type": "primary_cell"}]))
        return SimpleNamespace(cage=track)

    monkeypatch.setattr(bwp, "predict_sequence_resilient", fake_predict)
    monkeypatch.setattr(predict_dataset.DatasetPredictor, "_client", lambda self: object())
    monkeypatch.setattr(predict_dataset, "SUPPORTED_LENGTHS", (len(REF),))
    predictor = predict_dataset.DatasetPredictor(lab_dataset.path, ["cage"], ["CL:9"])
    assert predictor.run() == {"predicted": 2, "skipped": 0, "failed": 0}
    pred = lab_dataset.path / "individuals" / "S1" / "windows" / "G1" / "predictions_H2"
    assert np.load(pred / "cage.npz")["values"].shape == (len(REF), 1)
    assert json.loads((pred / "cage_metadata.json").read_text())["metadata"][0]["ontology_curie"] == "CL:9"
    metadata = json.loads((lab_dataset.path / "dataset_metadata.json").read_text())
    assert metadata["alphagenome_outputs"] == ["CAGE"] and metadata["ontology_details"]["CL:9"]["biosample_name"] == "cell nine"
    assert json.loads((lab_dataset.path / "references" / "windows" / "G1" / "window_metadata.json").read_text())["outputs"] == ["CAGE"]
    # Everything is done: a second run makes no calls.
    assert predict_dataset.DatasetPredictor(lab_dataset.path, ["CAGE"], ["CL:9"]).run()["skipped"] == 2 and len(calls) == 2
    with pytest.raises(ValueError):
        predict_dataset.DatasetPredictor(lab_dataset.path, ["NOT_AN_OUTPUT"], ["CL:9"])


def test_predict_dataset_includes_windows_missing_from_metadata(lab_dataset):
    from genomics.workflows.alphagenome import predict_dataset

    # A window added after the build (not in dataset_metadata.json "genes") is still predicted.
    (lab_dataset.path / "references" / "windows" / "G2").mkdir()
    assert predict_dataset.DatasetPredictor(lab_dataset.path, ["CAGE"], ["CL:9"]).genes == ["G1", "G2"]
    assert predict_dataset.DatasetPredictor(lab_dataset.path, ["CAGE"], ["CL:9"], genes=["G2"]).genes == ["G2"]
    with pytest.raises(ValueError):
        predict_dataset.DatasetPredictor(lab_dataset.path, ["CAGE"], ["CL:9"], genes=["G3"])


def test_predict_dataset_stores_every_output_kind(lab_dataset, monkeypatch, tmp_path):
    """Binned ChIP tracks, contact maps and splice junctions are stored in their own shapes and
    read back: binned tracks expand to bases (through each haplotype's indels), the others are
    listed apart from the plottable tracks."""
    pytest.importorskip("alphagenome.models.dna_client")
    from genomics.workflows.alphagenome import predict_dataset
    from genomics.workflows.dataset_builders.non_longevous import build_window_and_predict as bwp

    res = 8  # stands in for ChIP-seq's 128 bp (the test windows are 200 bp)
    n = len(REF)
    meta = lambda *records: _FakeMetadata(list(records))

    def fake_predict(client_box, make_client, seq, outputs, terms, timeout_s, max_attempts):
        assert terms is None  # --all-tissues
        bins = np.arange(n // res, dtype=np.float32)
        junctions = np.array([SimpleNamespace(start=20, end=90, strand="+"), SimpleNamespace(start=40, end=150, strand="-")], dtype=object)  # as JunctionData holds them
        return SimpleNamespace(
            chip_histone=SimpleNamespace(values=np.stack([bins, -bins], axis=1), resolution=res,
                                         metadata=meta({"ontology_curie": "CL:9", "strand": ".", "histone_mark": "H3K27ac", "biosample_name": "cell nine"},
                                                       {"ontology_curie": "CL:9", "strand": ".", "histone_mark": "H3K4me3", "biosample_name": "cell nine"})),
            contact_maps=SimpleNamespace(values=np.ones((4, 4, 1), np.float32), resolution=50, metadata=meta({"ontology_curie": "EFO:1", "strand": "."})),
            splice_junctions=SimpleNamespace(junctions=junctions, values=np.array([[1.0], [2.0]], np.float32), metadata=meta({"ontology_curie": "CL:9", "strand": "+"})),
            splice_sites=SimpleNamespace(values=np.zeros((n, 0), np.float32), metadata=meta()),
        )

    monkeypatch.setattr(bwp, "predict_sequence_resilient", fake_predict)
    monkeypatch.setattr(predict_dataset.DatasetPredictor, "_client", lambda self: object())
    monkeypatch.setattr(predict_dataset, "SUPPORTED_LENGTHS", (n,))
    outputs = ["chip_histone", "contact_maps", "splice_junctions", "splice_sites"]
    assert predict_dataset.DatasetPredictor(lab_dataset.path, outputs, None).run()["predicted"] == 2
    pred = lab_dataset.path / "individuals" / "S1" / "windows" / "G1" / "predictions_H1"
    with np.load(pred / "chip_histone.npz") as data:
        assert data["values"].shape == (n // res, 2) and int(data["resolution"]) == res
    with np.load(pred / "splice_junctions.npz") as data:
        assert data["starts"].tolist() == [20, 40] and data["strands"].tolist() == ["+", "-"] and data["values"].shape == (2, 1)
    assert json.loads((pred / "contact_maps_metadata.json").read_text())["kind"] == "contact_map"

    dataset = Dataset(lab_dataset.path)
    info = dataset.gene_info("G1")
    assert set(info["outputs"]) == {"rna_seq", "chip_histone"}
    assert info["outputs"]["chip_histone"]["resolution"] == res and info["outputs"]["chip_histone"]["length"] == n
    assert info["outputs"]["chip_histone"]["tracks"][0]["short"] == "H3K27ac · cell nine"
    assert info["outputs"]["rna_seq"]["resolution"] == 1  # legacy file without a stored resolution
    others = info["other_outputs"]
    assert others["contact_maps"]["kind"] == "contact_map" and others["contact_maps"]["bins"] == 4
    assert others["splice_junctions"]["junctions"] == 2 and others["splice_sites"]["tracks"] == []

    signals = SignalService(64 << 20, DiskArrayCache(tmp_path / "cache"))
    hap = signals.haplotype_window(dataset, "S1", "G1", "H1", "chip_histone", "haplotype", 0, n, [0])[:, 0]
    np.testing.assert_array_equal(hap, np.arange(n) // res)
    # Reference offsets after H1's 3 bp insertion (after offset 30) read haplotype bases +3.
    ref = signals.haplotype_window(dataset, "S1", "G1", "H1", "chip_histone", "reference", 60, 70, [1])[:, 0]
    np.testing.assert_array_equal(ref, -(np.arange(63, 73) // res))
    assert signals.domain_length(dataset, "G1", "chip_histone", "haplotype") == n


def test_old_binned_files_infer_resolution_from_the_window(lab_dataset):
    pred = lab_dataset.path / "individuals" / "S1" / "windows" / "G1" / "predictions_H1"
    np.savez_compressed(pred / "chip_tf.npz", values=np.zeros((len(REF) // 8, 3), np.float32))
    described = Dataset(lab_dataset.path).gene_info("G1")["outputs"]["chip_tf"]
    assert described["resolution"] == 8 and described["length"] == len(REF) and len(described["tracks"]) == 3


def test_output_specs_catalog_and_regulatory_presets():
    from genomics.visualizer import launch
    from genomics.workflows.alphagenome import catalog
    from genomics.workflows.alphagenome.outputs import ALL_OUTPUTS, normalize_outputs, stored_bytes
    from genomics.workflows.alphagenome.regulatory_regions import REGULATORY_REGIONS, preset_regions

    assert normalize_outputs(["all"]) == list(ALL_OUTPUTS) and len(ALL_OUTPUTS) == 11
    assert normalize_outputs(["rna-seq", "CHIP_TF", "rna_seq"]) == ["RNA_SEQ", "CHIP_TF"]
    with pytest.raises(ValueError):
        normalize_outputs(["hic"])
    assert stored_bytes("CHIP_TF", 524288, 10) == 4096 * 10 * 4 and stored_bytes("CONTACT_MAPS", 1048576, 1) == 512 * 512 * 4

    class Frame:
        def __init__(self, rows):
            self.rows = rows

        def to_dict(self, orient="records"):
            return self.rows

    metadata = SimpleNamespace(
        rna_seq=Frame([{"name": "a", "strand": "+", "ontology_curie": "CL:1", "biosample_name": "melanocyte", "biosample_type": "primary_cell"},
                       {"name": "b", "strand": "-", "ontology_curie": "CL:1", "biosample_name": "melanocyte", "biosample_type": "primary_cell"}]),
        chip_tf=Frame([{"name": "c", "strand": ".", "ontology_curie": "EFO:2", "biosample_name": "K562", "transcription_factor": "CTCF", "gtex_tissue": float("nan")}]),
        splice_sites=Frame([{"name": "donor", "strand": "+"}]),
    )
    built = catalog.build_catalog(metadata, source="test")
    terms = {t["curie"]: t for t in built["ontologies"]}
    assert terms["CL:1"]["outputs"] == {"RNA_SEQ": 2} and terms["EFO:2"]["marks"] == {"CHIP_TF": ["CTCF"]}
    assert built["outputs"]["SPLICE_SITES"]["tracks"] == 1 and built["outputs"]["SPLICE_SITES"]["tissue_specific"] is False
    assert "gtex_tissue" not in built["tracks"]["CHIP_TF"][0] and "tracks" not in catalog.summary(built)

    names = [r["name"] for r in REGULATORY_REGIONS]
    assert len(names) == len(set(names)) and all(launch.WINDOW_NAME_RE.match(n) for n in names)
    parsed = launch.parse_regions(preset_regions() + ["mine=chr2:1,000-2,000"])
    assert parsed[0]["start"] == parsed[0]["end"] == REGULATORY_REGIONS[0]["position"] and parsed[-1] == {"name": "mine", "chrom": "chr2", "start": 1000, "end": 2000}
    with pytest.raises(launch.LaunchError):
        launch.parse_regions(["bad name=chr1:5"])


def test_predict_task_outputs_tissues_and_new_windows(lab_dataset, tmp_path, monkeypatch):
    from genomics.visualizer import launch

    title, steps, params, files = launch.predict_task(lab_dataset, {"outputs": ["all"], "all_tissues": True})
    command = steps[-1]["command"]
    assert "--all-tissues" in command and command[command.index("--outputs") + 1].count(",") == 10 and params["ontology_terms"] == "all"
    with pytest.raises(launch.LaunchError):
        launch.predict_task(lab_dataset, {"outputs": ["RNA_SEQ"]})  # tissue-specific output without tissues
    assert "--all-tissues" in launch.predict_task(lab_dataset, {"outputs": ["SPLICE_SITES"]})[1][-1]["command"]  # not tissue-specific
    # The synthetic dataset records no source VCF: new windows are refused, existing ones are fine.
    assert not launch.extension_source(lab_dataset)["available"]
    with pytest.raises(launch.LaunchError, match="Cannot build new windows"):
        launch.predict_task(lab_dataset, {"outputs": ["CAGE"], "ontology_terms": ["CL:1"], "new_regions": ["R2=chr1:1100"]})
    assert launch.predict_task(lab_dataset, {"outputs": ["CAGE"], "ontology_terms": ["CL:1"], "new_genes": ["G1"]})[1][0]["title"] == "AlphaGenome predictions"

    # With a source VCF and reference, new regions get a build step before the predictions.
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 5000 + "\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t5000\t6\t5000\t5001\n")
    vcf = tmp_path / "cohort.vcf.gz"
    vcf.write_bytes(b"")
    metadata = json.loads((lab_dataset.path / "dataset_metadata.json").read_text())
    metadata["raw_variant_sources"] = {"vcf_pattern": str(vcf), "reference_fasta": str(fasta)}
    metadata["window_size"] = 16384
    (lab_dataset.path / "dataset_metadata.json").write_text(json.dumps(metadata))
    monkeypatch.setattr("genomics.workflows.dataset_builders.vcf_import.builder.missing_tools", lambda: [])
    dataset = Dataset(lab_dataset.path)
    preset = {"name": "R2", "chrom": "chr1", "start": 1100, "end": 1100}
    title, steps, params, files = launch.predict_task(dataset, {"outputs": ["CAGE"], "ontology_terms": ["CL:1"], "genes": ["G1"], "new_regions": [preset, "R3=chr1:2000-2100"]})
    assert [s["title"] for s in steps] == ["Build 2 new window(s) for every sample", "AlphaGenome predictions"]
    spec = json.loads(files["extend_spec.json"])
    assert spec["extend"] is True and spec["samples"] == ["S1"] and [r["name"] for r in spec["regions"]] == ["R2", "R3"] and spec["window_size"] == 16384
    command = steps[1]["command"]
    assert command[command.index("--genes") + 1] == "G1,R2,R3" and params["new_windows"] == ["R2", "R3"]
    monkeypatch.setattr(launch, "default_gtf", lambda: None)
    with pytest.raises(launch.LaunchError, match="Gene symbols need"):
        launch.predict_task(dataset, {"outputs": ["CAGE"], "ontology_terms": ["CL:1"], "new_genes": ["KITLG"]})


def test_predict_dataset_reference_window(lab_dataset, monkeypatch):
    """Haplotype 'ref' predicts references/windows/<gene>/ref.window.fa into predictions_ref."""
    pytest.importorskip("alphagenome.models.dna_client")
    from genomics.workflows.alphagenome import predict_dataset
    from genomics.workflows.dataset_builders.non_longevous import build_window_and_predict as bwp

    seen = []

    def fake_predict(client_box, make_client, seq, outputs, terms, timeout_s, max_attempts):
        seen.append(seq)
        track = SimpleNamespace(values=np.full((len(seq), 1), 5.0, np.float32), metadata=_FakeMetadata([{"ontology_curie": "CL:1", "strand": "+", "biosample_name": "cell"}]))
        return SimpleNamespace(rna_seq=track)

    monkeypatch.setattr(bwp, "predict_sequence_resilient", fake_predict)
    monkeypatch.setattr(predict_dataset.DatasetPredictor, "_client", lambda self: object())
    monkeypatch.setattr(predict_dataset, "SUPPORTED_LENGTHS", (len(REF),))
    assert predict_dataset.DatasetPredictor(lab_dataset.path, ["RNA_SEQ"], ["CL:1"], haplotypes=["ref"]).run() == {"predicted": 1, "skipped": 0, "failed": 0}
    assert seen == [REF]  # the reference sequence, not a haplotype
    ref_dir = lab_dataset.path / "references" / "windows" / "G1"
    assert np.load(ref_dir / "predictions_ref" / "rna_seq.npz")["values"].shape == (len(REF), 1)
    assert json.loads((ref_dir / "window_metadata.json").read_text())["reference_outputs"] == ["RNA_SEQ"]
    metadata = json.loads((lab_dataset.path / "dataset_metadata.json").read_text())
    assert metadata["window_catalog"]["G1"]["reference_outputs"] == ["RNA_SEQ"] and "alphagenome_outputs" not in metadata
    assert Dataset(lab_dataset.path).gene_info("G1")["reference_outputs"] == ["rna_seq"]
    # Done already: skipped, while H1 + ref predicts only the haplotype.
    assert predict_dataset.DatasetPredictor(lab_dataset.path, ["RNA_SEQ"], ["CL:1"], haplotypes=["ref"]).run()["skipped"] == 1
    from genomics.visualizer import launch

    title, steps, params, files = launch.predict_task(lab_dataset, {"outputs": ["RNA_SEQ"], "ontology_terms": ["CL:1"], "haplotypes": ["ref"], "samples": ["S1"]})
    command = steps[-1]["command"]
    assert command[command.index("--haplotypes") + 1] == "ref" and "samples.txt" not in files and params["windows"] == 1 and "reference genome" in title


# -------------------------------------------------------------------------------------- training
def test_train_config_from_form(lab_dataset, tmp_path):
    from genomics.visualizer import launch

    catalog = DatasetCatalog()
    dataset = catalog.add(lab_dataset.path)
    dataset.set_field("pop", {"S1": "YRI"})
    body = {
        "target_field": "pop", "class_map": {"YRI": "dark", "CEU": "light"}, "target_name": "Pigment class",
        "genes": ["G1"], "output": "rna_seq", "window_center_size": 64, "feature_mode": "signals_only",
        "model_type": "CNN2", "num_epochs": 3, "train_split": 0.6, "val_split": 0.2, "test_split": 0.2, "run_name": "my run",
    }
    config, summary = launch.build_train_config(dataset, body, {}, tmp_path / "runs", tmp_path / "cache")
    di = config["dataset_input"]
    assert di["dataset_dir"] == str(dataset.path) and "dataset_id" not in di and di["ontology_terms"] == ["CL:1"]
    assert config["output"]["prediction_target"] == "pigment_class"
    assert config["output"]["derived_targets"]["pigment_class"]["class_map"] == {"dark": ["YRI"], "light": ["CEU"]}
    assert config["model"]["cnn2"]["kernel_stage1"][0] == 2 == summary["rows_per_gene"]  # 1 tissue x 2 strands
    assert di["results_dir"] == str(tmp_path / "runs" / "my-run") and config["data_split"]["family_split_mode"] == "ignore"
    launch.validate_config(config, tmp_path / "scratch")
    title, steps, params, files = launch.train_task(config, summary, evaluate_test=True)
    assert [s["title"] for s in steps] == ["Training", "Evaluation on the test split"] and "config.yaml" in files
    assert steps[0]["epochs"] == 3 and steps[0]["history_glob"].endswith("models/training_history.json")
    with pytest.raises(launch.LaunchError):
        launch.build_train_config(dataset, {**body, "class_map": {"YRI": "dark"}}, {}, tmp_path, tmp_path)
    with pytest.raises(launch.LaunchError):
        launch.build_train_config(dataset, {**body, "window_center_size": 10_000}, {}, tmp_path, tmp_path)


def test_custom_metadata_fields_reach_training_labels(tmp_path):
    pytest.importorskip("torch")
    from genomics.workflows.dataset_builders.non_longevous.genomic_dataset import GenomicLongevityDataset

    root = tmp_path / "ds"
    (root / "individuals" / "S1").mkdir(parents=True)
    (root / "dataset_metadata.json").write_text(json.dumps({"dataset_name": "x", "individuals": ["S1"]}))
    (root / "individuals" / "S1" / "individual_metadata.json").write_text(json.dumps({"sample_id": "S1", "eye_colour": "blue", "age": 34, "windows": []}))
    _inputs, outputs = GenomicLongevityDataset(root, load_predictions=False, load_sequences=False)[0]
    assert outputs["eye_colour"] == "blue" and outputs["age"] == 34 and outputs["population"] == ""


# ---------------------------------------------------------------------------------------- import
@pytest.mark.skipif(not (shutil.which("bcftools") and shutil.which("samtools") and shutil.which("bgzip")), reason="needs bcftools, samtools and bgzip")
def test_vcf_import_builds_canonical_layout(tmp_path):
    import subprocess

    from genomics.workflows.dataset_builders.vcf_import.builder import DatasetImporter

    chrom = "A" * 100 + REF * 5 + "C" * 100  # 1200 bp
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "\n".join(chrom[i:i + 60] for i in range(0, len(chrom), 60)) + "\n")
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=1,length=1200>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\tA2\tA3"]
    snv = 601
    alt = "T" if chrom[snv - 1] != "T" else "G"
    lines.append(f"1\t{snv}\t.\t{chrom[snv - 1]}\t{alt}\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0")
    lines.append(f"1\t620\t.\t{chrom[619]}\t{chrom[619]}AAAA\t.\tPASS\t.\tGT\t0|1\t0|0\t1|1")
    raw = tmp_path / "cohort.vcf"
    raw.write_text("\n".join(lines) + "\n")
    vcf = tmp_path / "cohort.vcf.gz"
    subprocess.run(["bgzip", "-c", str(raw)], stdout=open(vcf, "wb"), check=True)  # VCF uses "1", FASTA "chr1"
    spec = {
        "name": "tiny", "output_dir": str(tmp_path / "out"), "vcf": str(vcf), "reference_fasta": str(fasta), "window_size": 256,
        "regions": [{"name": "R1", "chrom": "1", "start": 610}], "samples": ["A1", "A2"],
        "metadata": {"A1": {"group": "x", "age": 3}, "A2": {"group": "y"}}, "workers": 2,
    }
    DatasetImporter(spec).run()
    out = tmp_path / "out"
    metadata = json.loads((out / "dataset_metadata.json").read_text())
    assert metadata["individuals"] == ["A1", "A2"] and metadata["genes"] == ["R1"] and metadata["group_distribution"] == {"x": 1, "y": 1}
    window = json.loads((out / "references" / "windows" / "R1" / "window_metadata.json").read_text())
    assert Path(window["raw_variant_source"]["vcf_path"]).exists() and window["chromosome"] == "chr1"
    case = out / "individuals" / "A1" / "windows" / "R1"
    h1 = "".join((case / "A1.H1.window.fixed.fa").read_text().splitlines()[1:])
    h2 = "".join((case / "A1.H2.window.fixed.fa").read_text().splitlines()[1:])
    ref = "".join((out / "references" / "windows" / "R1" / "ref.window.fa").read_text().splitlines()[1:])
    offset = snv - window["start"]
    assert len(h1) == len(h2) == 256 and h1[offset] == alt and h2[offset] == ref[offset] and "AAAA" in h2[offset:offset + 30]
    assert json.loads((out / "individuals" / "A2" / "individual_metadata.json").read_text())["windows"] == ["R1"]
    # The import reads like any canonical dataset, and resuming adds samples without rebuilding.
    dataset = Dataset(out)
    assert {f["name"] for f in dataset.fields} >= {"sample_id", "group", "age"}
    DatasetImporter({**spec, "samples": ["A1", "A2", "A3"], "metadata": {"A3": {"group": "x"}}}).run()
    metadata = json.loads((out / "dataset_metadata.json").read_text())
    assert metadata["individuals"] == ["A1", "A2", "A3"] and metadata["individuals_pedigree"]["A1"]["age"] == 3


    # Extending adds a window for every sample and leaves samples, metadata and provenance alone.
    before = json.loads((out / "dataset_metadata.json").read_text())
    DatasetImporter({"extend": True, "output_dir": str(out), "vcf": str(vcf), "reference_fasta": str(fasta), "window_size": 256,
                     "regions": [{"name": "R2", "chrom": "chr1", "start": 300, "end": 300}], "samples": ["A1", "A2", "A3", "NOT_IN_VCF"]}).run()
    after = json.loads((out / "dataset_metadata.json").read_text())
    assert after["genes"] == ["R1", "R2"] and after["individuals"] == before["individuals"] and after["source"] == before["source"]
    assert after["group_distribution"] == before["group_distribution"] and "R2" in after["window_catalog"]
    assert json.loads((out / "individuals" / "A3" / "individual_metadata.json").read_text())["windows"] == ["R1", "R2"]
    assert (out / "individuals" / "A1" / "windows" / "R2" / "A1.H2.window.fixed.fa").exists()
    assert json.loads((out / "window_extensions.json").read_text())[0]["regions"][0]["name"] == "R2"
    assert [r["name"] for r in json.loads((out / "import_spec.json").read_text())["regions"]] == ["R1", "R2"]
    assert Dataset(out).genes == ["R1", "R2"]


def test_saturation_scan_scores_each_window_and_skips_unchanged(lab_dataset, tmp_path, monkeypatch):
    from genomics.visualizer.jobs import JobCancelled

    service = _service(tmp_path)
    service.app.catalog = SimpleNamespace(add=lambda path: lab_dataset)
    calls = []

    def score(sample, overrides):
        # P(b) grows with the edited haplotypes' signal in column 0 (the fake AlphaGenome counts A bases).
        total = sum(float(o["rna_seq"][0][:, 0].sum()) for o in overrides.values())
        p = min(0.99, 0.5 + total / 1000.0)
        return np.array([1 - p, p])

    service.context = SimpleNamespace(outputs=["rna_seq"], genes=["G1"], class_names=["a", "b"], window_center_size=40, dataset_dir=str(lab_dataset.path),
                                      baseline=lambda sample: np.array([0.5, 0.5]), score=score)
    service.context_id = ("run", "best.pt")
    service._labels = {"S1": "a"}

    def fake_predict(sequence, outputs, terms):
        calls.append(sequence)
        values = np.zeros((len(sequence), 2), np.float32)
        values[:, 0] = np.frombuffer(sequence, np.uint8) == ord("A")
        return {"rna_seq": (values, [{"ontology_curie": "CL:1", "strand": "+"}, {"ontology_curie": "CL:1", "strand": "-"}])}

    service._predict = fake_predict
    spec = service.scan_spec({"sample": "S1", "gene": "G1", "op": "overwrite", "base": "A", "haplotypes": ["H1"], "start": 0, "end": 40, "size": 10})
    assert (spec["start"], spec["end"], spec["step"]) == (0, 40, 10)
    result = service.scan(spec, lambda *a: None)
    assert [r["start"] for r in result["rows"]] == [0, 10, 20, 30] and result["classes"] == ["a", "b"]
    assert all(r["delta"][1] > 0 and r["delta"][0] == pytest.approx(-r["delta"][1]) for r in result["rows"])
    assert result["calls"] == 4 and len(calls) == 4

    # Reverting to the reference only changes windows that hold a variant of H1 (the SNV at offset 10).
    calls.clear()
    revert = service.scan(service.scan_spec({"sample": "S1", "gene": "G1", "op": "reference", "haplotypes": ["H1"], "start": 0, "end": 40, "size": 10}), lambda *a: None)
    assert [r["changed"] for r in revert["rows"]] == [0, 1, 0, 0] and revert["calls"] == 1 and len(calls) == 1
    assert revert["rows"][0]["delta"] == [0.0, 0.0]

    with pytest.raises(PerturbError):
        service.scan_spec({"sample": "S1", "gene": "G1", "op": "delete"})
    with pytest.raises(PerturbError):
        service.scan_spec({"sample": "S2", "gene": "G1"})
    monkeypatch.setattr("genomics.visualizer.perturb.MAX_SCAN_POSITIONS", 10)
    with pytest.raises(PerturbError):  # more windows than allowed
        service.scan_spec({"sample": "S1", "gene": "G1", "start": 0, "end": 40, "size": 1, "step": 1})

    def cancelled(*args):
        return None

    cancelled.cancelled = True
    with pytest.raises(JobCancelled):
        service.scan(spec, cancelled)


def test_negative_control_adds_label_permutation_to_the_config(lab_dataset, tmp_path):
    """The training form's negative controls: a global shuffle, a shuffle within a field, and off."""
    from genomics.visualizer import launch

    dataset = DatasetCatalog().add(lab_dataset.path)
    dataset.set_field("pop", {"S1": "YRI"})
    dataset.set_field("region", {"S1": "AFR"})
    body = {
        "target_field": "pop", "class_map": {"YRI": "dark", "CEU": "light"}, "genes": ["G1"], "output": "rna_seq",
        "window_center_size": 64, "model_type": "CNN2", "num_epochs": 1,
        "train_split": 0.6, "val_split": 0.2, "test_split": 0.2,
    }
    plain, _ = launch.build_train_config(dataset, body, {}, tmp_path / "runs", tmp_path / "cache")
    assert "label_permutation" not in plain

    shuffled, summary = launch.build_train_config(dataset, {**body, "negative_control": {"labels": "permute", "seed": 5}}, {}, tmp_path / "runs", tmp_path / "cache")
    assert shuffled["label_permutation"] == {"enabled": True, "random_seed": 5}
    assert "chance level" in summary["negative_control"]["description"]
    launch.validate_config(shuffled, tmp_path / "scratch")

    within, summary = launch.build_train_config(
        dataset, {**body, "negative_control": {"labels": "permute_within", "stratify_field": "region", "seed": 5}},
        {}, tmp_path / "runs", tmp_path / "cache")
    assert within["label_permutation"] == {"enabled": True, "random_seed": 5, "stratify_field": "region"}
    assert "within region" in summary["negative_control"]["description"]
    assert "Negative control" in within["metadata"]["note"]
    launch.validate_config(within, tmp_path / "scratch")

    # A control panel is just a different gene set, but the run records that it is one.
    _, summary = launch.build_train_config(dataset, {**body, "control_panel": True}, {}, tmp_path / "runs", tmp_path / "cache")
    assert "matched control windows" in summary["negative_control"]["description"]

    for bad in ({"labels": "sideways"}, {"labels": "permute_within", "stratify_field": "nope"},
                {"labels": "permute_within", "stratify_field": "pop"}):  # within the target itself is a no-op
        with pytest.raises(launch.LaunchError):
            launch.build_train_config(dataset, {**body, "negative_control": bad}, {}, tmp_path, tmp_path)


def test_task_retry_runs_the_same_recipe_in_a_new_directory(tmp_path):
    """Retry re-runs a finished task's recipe, with its own task dir and its own copy of the files."""
    tm = TaskManager(tmp_path / "tasks")
    py = sys.executable
    # The step writes into its own task dir, so a retry must not touch the first task's output.
    step = [py, "-c", "import sys; open(sys.argv[1], 'a').write(open(sys.argv[2]).read())", "{task_dir}/out.txt", "{task_dir}/in.txt"]
    first = _wait(tm, tm.create("demo", "Demo", [{"title": "one", "command": step}], params={"k": "v"}, files={"in.txt": "data"})["id"])
    assert first["status"] == "done"

    second = tm.retry(first["id"])
    done = _wait(tm, second["id"])
    assert done["id"] != first["id"] and done["status"] == "done"
    assert done["title"] == first["title"] and done["kind"] == first["kind"] and done["params"] == {"k": "v"}
    assert (Path(done["dir"]) / "out.txt").read_text() == "data"  # its own copy, written once
    assert (Path(first["dir"]) / "out.txt").read_text() == "data"  # the first task's output is untouched

    running = tm.create("demo", "Slow", [{"title": "s", "command": [py, "-c", "import time; time.sleep(30)"]}])
    with pytest.raises(ValueError, match="still running"):
        tm.retry(running["id"])
    tm.cancel(running["id"])

    # A task written before recipes were stored says so instead of running the wrong thing.
    spec_path = Path(first["dir"]) / "task.json"
    spec = json.loads(spec_path.read_text())
    spec.pop("recipe")
    spec_path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="older version"):
        tm.retry(first["id"])
    with pytest.raises(KeyError):
        tm.retry("no-such-task")
