import json
from pathlib import Path

import pytest

pytest.importorskip("rich")

from genomics.predictors.genotype_based.data.layout import is_dataset_dir, materialize_dataset


SAMPLES = ("HG00096", "HG00097")
GENE = "TYR"
PER_SAMPLE_FILES = (
    "{s}.H1.window.fixed.fa",
    "{s}.H2.window.fixed.fa",
    "{s}.H1.window.raw.fa",
    "{s}.H2.window.raw.fa",
    "{s}.window.vcf.gz",
    "{s}.window.vcf.gz.tbi",
    "{s}.window.consensus_ready.vcf.gz",
    "{s}.window.consensus_ready.vcf.gz.tbi",
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _make_legacy_dataset(root: Path) -> Path:
    """Old top3 layout: ref.window.fa repeated inside every individual's window."""
    source = root / "top3" / "non_longevous_results_genes_1000_all"
    _write(source / "dataset_metadata.json", json.dumps({"genes": [GENE], "window_size": 524288}))
    _write(source / "selected_samples.csv", "sample\n" + "\n".join(SAMPLES) + "\n")
    _write(source / "gtf_cache.feather", "gtf")
    for sample in SAMPLES:
        ind_dir = source / "individuals" / sample
        _write(
            ind_dir / "individual_metadata.json",
            json.dumps(
                {
                    "sample_id": sample,
                    "population": "GBR",
                    "superpopulation": "EUR",
                    "windows": [GENE],
                    "window_metadata": {GENE: {"chromosome": "chr11", "start": 1, "end": 524288}},
                }
            ),
        )
        window = ind_dir / "windows" / GENE
        _write(window / "ref.window.fa", ">ref\nACGT\n")
        _write(window / "prediction_H1.ok.txt", "ok")
        for pattern in PER_SAMPLE_FILES:
            name = pattern.format(s=sample)
            _write(window / name, name)
        for hap in ("H1", "H2"):
            _write(window / f"predictions_{hap}" / "rna_seq.npz", f"{sample}-{hap}")
            _write(window / f"predictions_{hap}" / "rna_seq_metadata.json", "{}")
    return source


def _assert_canonical(target: Path) -> None:
    assert is_dataset_dir(target)
    assert (target / "references" / "windows" / GENE / "ref.window.fa").read_text(encoding="utf-8") == ">ref\nACGT\n"
    assert json.loads((target / "references" / "windows" / GENE / "window_metadata.json").read_text(encoding="utf-8"))["chromosome"] == "chr11"
    assert (target / "selected_samples.csv").exists()
    assert (target / "gtf_cache.feather").exists()
    for sample in SAMPLES:
        window = target / "individuals" / sample / "windows" / GENE
        meta = json.loads((target / "individuals" / sample / "individual_metadata.json").read_text(encoding="utf-8"))
        assert meta["windows"] == [GENE]
        assert "window_metadata" not in meta
        for pattern in PER_SAMPLE_FILES:
            name = pattern.format(s=sample)
            assert (window / name).read_text(encoding="utf-8") == name
        for hap in ("H1", "H2"):
            assert (window / f"predictions_{hap}" / "rna_seq.npz").read_text(encoding="utf-8") == f"{sample}-{hap}"
        assert not (window / "ref.window.fa").exists()


def _source_data_files(source: Path) -> list:
    leftovers = {"dataset_metadata.json", "individual_metadata.json", "prediction_H1.ok.txt"}
    return sorted(str(p.relative_to(source)) for p in source.rglob("*") if p.is_file() and p.name not in leftovers)


def test_materialize_copy_keeps_source(tmp_path):
    source = _make_legacy_dataset(tmp_path)
    before = _source_data_files(source)
    target = materialize_dataset(source, tmp_path / "v1" / "1kG_high_coverage")
    _assert_canonical(target)
    assert _source_data_files(source) == before


def test_materialize_move_empties_source(tmp_path):
    source = _make_legacy_dataset(tmp_path)
    target = materialize_dataset(source, tmp_path / "v1" / "1kG_high_coverage", move=True)
    _assert_canonical(target)
    # Only metadata and .ok markers remain; the duplicate ref.window.fa copies are gone.
    assert _source_data_files(source) == []
    for sample in SAMPLES:
        assert (source / "individuals" / sample / "individual_metadata.json").exists()


def test_materialize_move_resumes_after_interruption(tmp_path):
    source = _make_legacy_dataset(tmp_path)
    target = tmp_path / "v1" / "1kG_high_coverage"
    # Simulate a run killed after the first individual was moved (no layout marker yet).
    first = SAMPLES[0]
    src_window = source / "individuals" / first / "windows" / GENE
    dst_window = target / "individuals" / first / "windows" / GENE
    for pattern in PER_SAMPLE_FILES:
        name = pattern.format(s=first)
        dst_window.mkdir(parents=True, exist_ok=True)
        (src_window / name).replace(dst_window / name)
    (target / "references" / "windows" / GENE).mkdir(parents=True)
    (src_window / "ref.window.fa").replace(target / "references" / "windows" / GENE / "ref.window.fa")
    assert not is_dataset_dir(target)

    materialize_dataset(source, target, move=True)
    _assert_canonical(target)
    assert _source_data_files(source) == []


def test_materialize_move_keeps_conflicting_source_file(tmp_path):
    source = _make_legacy_dataset(tmp_path)
    odd_ref = source / "individuals" / SAMPLES[1] / "windows" / GENE / "ref.window.fa"
    odd_ref.write_text(">ref\nTTTT\n", encoding="utf-8")
    target = materialize_dataset(source, tmp_path / "v1" / "1kG_high_coverage", move=True)
    _assert_canonical(target)
    assert _source_data_files(source) == [str(odd_ref.relative_to(source))]


def test_materialize_move_across_filesystems(tmp_path, monkeypatch):
    import errno
    import os

    from genomics.predictors.genotype_based.data import layout

    source = _make_legacy_dataset(tmp_path)
    real_replace = os.replace

    def cross_device_replace(src, dst):
        # Renames out of the source tree fail as they would across mounts;
        # the .partial -> final rename inside the target still works.
        if str(src).startswith(str(source)):
            raise OSError(errno.EXDEV, "Invalid cross-device link")
        return real_replace(src, dst)

    monkeypatch.setattr(layout.os, "replace", cross_device_replace)
    target = materialize_dataset(source, tmp_path / "v1" / "1kG_high_coverage", move=True)
    _assert_canonical(target)
    assert _source_data_files(source) == []
    assert not list(target.rglob("*.partial"))
