"""Integration test for the pigmentation_sequence_lab app's HTTP surface.

Runs a real ThreadingHTTPServer against a small synthetic dataset directory, with the
alignment (DynamicIndelAligner/BcftoolsChainMapper -- these shell out to the `bcftools` binary,
not available in every environment) and the model/live-API collaborators (torch checkpoint,
AlphaGenome API) replaced by lightweight fakes injected via monkeypatch. This exercises this
app's own request-routing and orchestration logic (individual vs. class dispatch, true RNA-seq
class averaging, CAGE representative-individual fallback, base-frequency LOD selection, and the
400 rejection of edits on a class subject) -- the pieces that are specific to this app rather than
already covered by the analysis-module unit tests.
"""
import gzip
import http.client
import json
import socket
import threading
from pathlib import Path

import numpy as np
import pytest

from genomics.predictors.genotype_based.apps.alphagenome_track_viewer import DEFAULT_POINTS, AlphaGenomeRepository
from genomics.predictors.genotype_based.apps.pigmentation_sequence_lab import (
    PigmentationSequenceLabHandler,
    PigmentationSequenceLabRepository,
)

GENES = ["GENE1"]
WINDOW_LENGTH = 20
HAPLOTYPES = ("H1", "H2")
SAMPLE_VALUES = {"IND1": 1.0, "IND2": 2.0, "IND3": 3.0}  # IND1/IND3 = strong (YRI), IND2 = weak (CEU)
TRACK_METADATA = [{"ontology_curie": "CL:0", "strand": "+"}, {"ontology_curie": "CL:1", "strand": "+"}]
REF_SEQUENCE = "ACGT" * (WINDOW_LENGTH // 4)


def _write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def _write_fasta(path: Path, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">seq\n{sequence}\n")


def _build_synthetic_dataset(root: Path) -> Path:
    dataset_dir = root / "dataset"
    _write_json(
        dataset_dir / "dataset_metadata.json",
        {
            "dataset_name": "synthetic",
            "individuals": list(SAMPLE_VALUES),
            "individuals_pedigree": {
                "IND1": {"population": "YRI", "superpopulation": "AFR", "sex": "male"},
                "IND2": {"population": "CEU", "superpopulation": "EUR", "sex": "female"},
                "IND3": {"population": "YRI", "superpopulation": "AFR", "sex": "female"},
            },
            "genes": GENES,
            "alphagenome_outputs": ["rna_seq"],
            "window_size": WINDOW_LENGTH,
        },
    )
    for gene in GENES:
        _write_fasta(dataset_dir / "references" / "windows" / gene / "ref.window.fa", REF_SEQUENCE)
        _write_json(
            dataset_dir / "references" / "windows" / gene / "window_metadata.json",
            {"chromosome": "chr1", "start": 1000, "end": 1000 + WINDOW_LENGTH - 1},
        )
        for sample_id, value in SAMPLE_VALUES.items():
            base_dir = dataset_dir / "individuals" / sample_id / "windows" / gene
            base_dir.mkdir(parents=True, exist_ok=True)
            # No-variant VCF (both haplotypes match the reference exactly): needed by
            # haplotype_edit_geometry.haplotype_indel_events for the affects_prediction check.
            with gzip.open(base_dir / f"{sample_id}.window.vcf.gz", "wt") as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_id + "\n")
            for haplotype in HAPLOTYPES:
                _write_fasta(base_dir / f"{sample_id}.{haplotype}.window.fixed.fa", REF_SEQUENCE)
                pred_dir = base_dir / f"predictions_{haplotype}"
                pred_dir.mkdir(parents=True, exist_ok=True)
                values = np.full((WINDOW_LENGTH, len(TRACK_METADATA)), value, dtype=np.float32)
                np.savez_compressed(pred_dir / "rna_seq.npz", values=values)
                (pred_dir / "rna_seq_metadata.json").write_text(json.dumps({"metadata": TRACK_METADATA}))
    return dataset_dir


class _FakeAligner:
    def __init__(self, expanded_length: int):
        self.expanded_length = expanded_length

    def build_alignment_axis_for_gene(self, gene, sample_ids):
        return None

    def get_expanded_length(self, gene):
        return self.expanded_length

    def get_reference_centered_expanded_slice(self, gene, window_size):
        size = max(int(window_size), 1)
        center = self.expanded_length // 2
        if self.expanded_length <= size:
            return {
                "center_ref_idx": center, "center_expanded_idx": center,
                "expanded_start": 0, "expanded_end": self.expanded_length,
                "expanded_length": self.expanded_length, "window_size": self.expanded_length,
            }
        half = size // 2
        start = max(0, center - half)
        end = min(self.expanded_length, start + size)
        return {
            "center_ref_idx": center, "center_expanded_idx": center,
            "expanded_start": start, "expanded_end": end,
            "expanded_length": self.expanded_length, "window_size": size,
        }


class _FakeChainMapper:
    def __init__(self, dataset_dir: Path, window_length: int):
        self.dataset_dir = dataset_dir
        self.consensus_dataset_dir = dataset_dir
        self.window_length = window_length

    def get_haplotype_entry(self, gene, sample_id, haplotype):
        fasta_path = self.dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.{haplotype}.window.fixed.fa"
        if not fasta_path.exists():
            return None
        identity = list(range(self.window_length))
        return {
            "copy_from_indices": identity, "expanded_indices": identity,
            "mapping_method": "bcftools_chain", "fasta_path": str(fasta_path), "source_start_idx": 0,
        }


class _FakeModelContext:
    def __init__(self, dataset_dir: Path):
        self.dataset_dir = dataset_dir
        self.genes = GENES
        self.rna_seq_ontology_terms = [rec["ontology_curie"] for rec in TRACK_METADATA]
        self.window_center_size = WINDOW_LENGTH
        self.class_names = ["weak pigmentation", "strong pigmentation"]
        self.strong_pigmentation_idx = 1

    def load_raw_prediction(self, sample_id, gene, haplotype, output="rna_seq"):
        pred_dir = self.dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
        values = np.load(pred_dir / f"{output}.npz")["values"]
        metadata = json.loads((pred_dir / f"{output}_metadata.json").read_text())["metadata"]
        return values, metadata

    def baseline_tensor(self, sample_id):
        return {"sample_id": sample_id, "overrides": {}}

    def build_individual_tensor(self, sample_id, overrides=None):
        return {"sample_id": sample_id, "overrides": dict(overrides or {})}

    def strong_pigmentation_probability(self, tensor):
        return 0.9 if tensor.get("overrides") else 0.3


class _FakeLivePredictor:
    def __init__(self):
        self.calls = []

    def predict_sequence(self, *, cache_key, sequence, interval, output_type, ontology_terms):
        self.calls.append((cache_key, output_type, sequence))
        values = np.full((WINDOW_LENGTH, len(TRACK_METADATA)), 7.0 if output_type == "CAGE" else 5.0, dtype=np.float32)
        return values, [dict(rec) for rec in TRACK_METADATA]


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture
def running_app(tmp_path, monkeypatch):
    dataset_dir = _build_synthetic_dataset(tmp_path)

    repository = object.__new__(PigmentationSequenceLabRepository)
    repository.dataset_dir = dataset_dir
    repository.consensus_dataset_dir = dataset_dir
    repository.ag_repo = AlphaGenomeRepository(dataset_dir)
    repository.model_context = _FakeModelContext(dataset_dir)
    repository._explicit_api_key = "fake-key"
    repository._live_predictor = _FakeLivePredictor()
    repository._full_aligner = None
    repository._full_chain_mapper = None
    repository._representative_sample_cache = {}
    repository._full_axis_built_genes = set()
    ag_options = repository.ag_repo.options_payload()
    repository._individuals_with_class = [
        {**row, "pigmentation_class": repository.ag_repo._class_label_for_row(row, "pigmentation")}
        for row in ag_options["individuals"]
    ]
    repository._defaults = {
        "subject_type": "individual",
        "subject": repository.ag_repo.individuals[0] if repository.ag_repo.individuals else None,
        "gene": repository.model_context.genes[0] if repository.model_context.genes else None,
        "haplotype": "H1",
        "output": "rna_seq",
        "start": ag_options["defaults"]["start"],
        "length": ag_options["defaults"]["length"],
        "points": DEFAULT_POINTS,
    }

    fake_aligner = _FakeAligner(WINDOW_LENGTH)
    fake_chain_mapper = _FakeChainMapper(dataset_dir, WINDOW_LENGTH)
    monkeypatch.setattr(repository.ag_repo, "_get_aligner", lambda: fake_aligner)
    monkeypatch.setattr(repository.ag_repo, "_get_chain_mapper", lambda: fake_chain_mapper)
    monkeypatch.setattr(repository, "_get_full_aligner", lambda: fake_aligner)
    monkeypatch.setattr(repository, "_get_full_chain_mapper", lambda: fake_chain_mapper)
    monkeypatch.setattr(repository, "_get_live_predictor", lambda: repository._live_predictor)

    from http.server import ThreadingHTTPServer

    class Handler(PigmentationSequenceLabHandler):
        pass

    Handler.repository = repository
    port = _free_port()
    server = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield port, repository
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def _get(port, path):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    try:
        conn.request("GET", path)
        resp = conn.getresponse()
        return resp.status, json.loads(resp.read())
    finally:
        conn.close()


def _post(port, path, body):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    try:
        payload = json.dumps(body).encode("utf-8")
        conn.request("POST", path, body=payload, headers={"Content-Type": "application/json"})
        resp = conn.getresponse()
        return resp.status, json.loads(resp.read())
    finally:
        conn.close()


def test_options_lists_genes_and_individuals(running_app):
    port, _ = running_app
    status, data = _get(port, "/api/options")
    assert status == 200
    assert data["genes"] == GENES
    sample_ids = {row["sample_id"] for row in data["individuals"]}
    assert sample_ids == set(SAMPLE_VALUES)
    by_id = {row["sample_id"]: row for row in data["individuals"]}
    assert by_id["IND1"]["pigmentation_class"] == "strong pigmentation"
    assert by_id["IND2"]["pigmentation_class"] == "weak pigmentation"


def test_tracks_individual_rna_seq(running_app):
    port, _ = running_app
    status, data = _get(
        port,
        f"/api/tracks?subject_type=individual&subject=IND1&gene=GENE1&haplotypes=H1&output=rna_seq&track=0&start=1&length={WINDOW_LENGTH}&points=100",
    )
    assert status == 200
    assert len(data["series"]) == 1
    assert all(abs(v - 1.0) < 1e-6 for v in data["series"][0]["y"])


def test_tracks_class_rna_seq_is_true_average(running_app):
    port, _ = running_app
    status, data = _get(
        port,
        f"/api/tracks?subject_type=class&subject=strong&gene=GENE1&haplotypes=H1&output=rna_seq&track=0&start=1&length={WINDOW_LENGTH}&points=100",
    )
    assert status == 200
    assert len(data["series"]) == 1
    # strong pigmentation = IND1 (1.0) + IND3 (3.0) -> true average 2.0
    assert all(abs(v - 2.0) < 1e-6 for v in data["series"][0]["y"])
    assert data["cage_is_true_average"] is None


def test_tracks_class_cage_uses_representative_individual(running_app):
    port, repository = running_app
    status, data = _get(
        port,
        f"/api/tracks?subject_type=class&subject=weak&gene=GENE1&haplotypes=H1&output=cage&track=0&start=1&length={WINDOW_LENGTH}&points=100",
    )
    assert status == 200
    assert data["cage_is_true_average"] is False
    assert data["cage_representative_sample"] == "IND2"
    assert all(abs(v - 7.0) < 1e-6 for v in data["series"][0]["y"])
    # CAGE gets cached into the canonical per-individual location.
    cage_path = repository.dataset_dir / "individuals" / "IND2" / "windows" / "GENE1" / "predictions_H1" / "cage.npz"
    assert cage_path.exists()


def test_base_frequency_letters_lod_for_class(running_app):
    port, _ = running_app
    status, data = _get(
        port,
        f"/api/base-frequency?subject_type=class&subject=strong&gene=GENE1&haplotypes=H1,H2&start=1&length={WINDOW_LENGTH}&points=100",
    )
    assert status == 200
    assert data["lod"] == "letters"
    assert data["mode"] == "frequency"
    assert data["n_used"] == 4  # 2 samples (IND1, IND3) x 2 haplotypes
    for idx in range(WINDOW_LENGTH):
        total = sum(data["frequencies"][sym][idx] for sym in ("A", "C", "G", "T", "gap"))
        assert abs(total - 1.0) < 1e-9


def test_base_frequency_density_lod_when_zoomed_out(running_app):
    port, _ = running_app
    status, data = _get(
        port,
        "/api/base-frequency?subject_type=individual&subject=IND1&gene=GENE1&haplotype=H1&start=1&length=1000&points=50",
    )
    assert status == 200
    assert data["lod"] == "density"
    assert "x" in data and "frequencies" in data


def test_base_frequency_individual_bases_no_gaps(running_app):
    port, _ = running_app
    status, data = _get(
        port,
        f"/api/base-frequency?subject_type=individual&subject=IND1&gene=GENE1&haplotype=H1&start=1&length={WINDOW_LENGTH}&points=100",
    )
    assert status == 200
    assert data["mode"] == "individual"
    assert data["bases"] == REF_SEQUENCE
    assert data["source_indices"] == list(range(WINDOW_LENGTH))


def test_edit_apply_on_individual_changes_prediction(running_app):
    port, repository = running_app
    status, data = _post(
        port,
        "/api/edit/apply",
        {
            "subject_type": "individual", "sample_id": "IND1", "gene": "GENE1", "haplotype": "H1",
            "start": 2, "end": 6, "op": "overwrite", "base": "A", "apply_to_prediction": ["rna_seq"],
        },
    )
    assert status == 200
    assert data["affects_prediction"] is True
    assert data["baseline_strong_pigmentation_probability"] == pytest.approx(0.3)
    assert data["edited_strong_pigmentation_probability"] == pytest.approx(0.9)
    assert data["delta_strong_pigmentation_probability"] == pytest.approx(0.6)
    assert data["edited_segment"] == "AAAA"
    # the fake live predictor was actually invoked with the edited sequence.
    assert repository._live_predictor.calls
    _, output_type, sequence = repository._live_predictor.calls[-1]
    assert output_type == "RNA_SEQ"
    assert sequence[2:6] == "AAAA"


def test_edit_apply_rejected_for_class_subject(running_app):
    port, _ = running_app
    status, data = _post(
        port,
        "/api/edit/apply",
        {
            "subject_type": "class", "sample_id": "strong", "gene": "GENE1", "haplotype": "H1",
            "start": 0, "end": 4, "op": "overwrite", "base": "A",
        },
    )
    assert status == 400
    assert "individual" in data["error"].lower()


def test_edit_preview_does_not_call_live_predictor(running_app):
    port, repository = running_app
    calls_before = len(repository._live_predictor.calls)
    status, data = _post(
        port,
        "/api/edit/preview",
        {"subject_type": "individual", "sample_id": "IND1", "gene": "GENE1", "haplotype": "H1", "start": 0, "end": 4, "op": "scramble", "seed": 1},
    )
    assert status == 200
    assert len(repository._live_predictor.calls) == calls_before
    assert sorted(data["edited_segment"]) == sorted(data["original_segment"])
