"""Quickstart presets (1000 Genomes, gnomAD HGDP + 1KG) and remote VCF/reference handling of the
VCF import. The network is never used: downloads, HEAD requests and the VCF inspection are faked."""
import csv
import gzip
import json
from pathlib import Path

import pytest

from genomics.workflows.dataset_builders.vcf_import import builder
from genomics.workflows.dataset_builders.vcf_import import quickstart as qs

PED = """FamilyID SampleID FatherID MotherID Sex Population Superpopulation
F1 P1 0 0 1 YRI AFR
F1 P2 0 0 2 YRI AFR
F1 C1 P1 P2 1 YRI AFR
F2 G1 0 0 2 GBR EUR
F3 G2 0 0 1 GBR EUR
F4 G3 0 0 1 GBR EUR
"""


def _gnomad_tsv(rows):
    header = ["s", "gnomad_sex_imputation", "relatedness_inference", "hgdp_tgp_meta", "high_quality"]
    lines = ["\t".join(header)]
    for sample, project, pop, karyotype, related, hq in rows:
        meta = {"project": project, "population": pop, "genetic_region": "EUR", "study_region": "Europe", "latitude": 44.0, "longitude": 39.0}
        lines.append("\t".join([sample, json.dumps({"sex_karyotype": karyotype}), json.dumps({"related": related}), json.dumps(meta), str(hq).lower()]))
    return gzip.compress(("\n".join(lines) + "\n").encode())


GNOMAD = _gnomad_tsv([
    ("HGDP1", "HGDP", "Adygei", "XX", False, True),
    ("HGDP2", "HGDP", "Adygei", "XY", True, True),
    ("HGDP3", "HGDP", "Basque", "XY", False, False),
    ("HGDP4", "HGDP", "Basque", "XY", False, True),
    ("NA1", "1000 Genomes", "CEU", "XX", False, True),
    ("CHMI", "synthetic_diploid_truth_sample", None, "XX", False, False),
])


def test_parse_1kg_pedigree_marks_trio_children_related():
    records = {r["sample_id"]: r for r in qs._parse_1kg(PED.encode())}
    assert records["C1"]["related"] and not records["P1"]["related"]
    assert records["P2"] == {"sample_id": "P2", "family_id": "F1", "sex": "Female", "population": "YRI", "superpopulation": "AFR", "related": False}


def test_parse_gnomad_flattens_json_columns_and_drops_qc_sample():
    records = {r["sample_id"]: r for r in qs._parse_gnomad(GNOMAD)}
    assert set(records) == {"HGDP1", "HGDP2", "HGDP3", "HGDP4", "NA1"}
    assert records["HGDP1"]["sex"] == "Female" and records["HGDP1"]["population"] == "Adygei" and records["HGDP1"]["latitude"] == 44.0
    assert records["HGDP2"]["related"] and records["HGDP3"]["high_quality"] is False


def test_select_samples_is_balanced_unrelated_and_scoped():
    records = qs._parse_1kg(PED.encode())
    picked = qs.select_samples(records, ["P1", "P2", "C1", "G1", "G2", "G3"], per_population=2)
    assert [r["sample_id"] for r in picked] == ["G1", "G2", "P1", "P2"]  # by population, then id; C1 is a child
    assert [r["sample_id"] for r in qs.select_samples(records, ["P1", "G3"], per_population=5)] == ["G3", "P1"]  # only VCF samples
    assert len(qs.select_samples(records, None, per_population=0)) == 6  # every sample, relatives included
    gnomad = qs._parse_gnomad(GNOMAD)
    assert [r["sample_id"] for r in qs.select_samples(gnomad, None, 4, scope="hgdp")] == ["HGDP1", "HGDP4"]
    assert [r["sample_id"] for r in qs.select_samples(gnomad, None, 4, scope="all")] == ["HGDP1", "HGDP4", "NA1"]


def test_prepare_writes_the_selected_samples_for_the_form(tmp_path, monkeypatch):
    monkeypatch.setattr(qs, "local_vcf", lambda preset: None)
    monkeypatch.setattr(qs, "resolve_vcf_files", lambda vcf: [vcf.replace("{chrom}", "chr1")])
    monkeypatch.setattr(qs, "inspect_vcf", lambda vcf: {"samples": ["HGDP1", "HGDP2", "HGDP3", "HGDP4", "NA1"], "sample_count": 5, "warnings": []})

    def fake_fetch(preset, cache_dir):
        path = Path(cache_dir) / preset.metadata_file
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(GNOMAD)
        return path

    monkeypatch.setattr(qs, "fetch_metadata", fake_fetch)
    res = qs.prepare("hgdp_1kg", tmp_path, per_population=4, scope="hgdp", reference=None)
    assert res["samples"] == ["HGDP1", "HGDP4"] and res["populations"] == 2 and res["vcf_source"] == "remote"
    assert res["reference_fasta"] == qs.REMOTE_REFERENCE and res["vcf_overrides"]["chrX"].endswith("chrX_non_par.full.shapeit5_rare.bcf")
    assert res["name"] == "HGDP quickstart 4 per population"
    # Records match what the form's metadata parser makes of the table (sex -> sex/sex_label).
    assert res["records"]["HGDP1"]["sex_label"] == "Female" and res["records"]["HGDP1"]["sex"] == 2
    with open(res["metadata_path"], encoding="utf-8") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    assert header[0] == "sample_id" and not {"project", "related", "high_quality"} & set(header)  # constant selection flags
    with pytest.raises(builder.ImportSpecError):
        qs.prepare("nope", tmp_path)


def test_remote_vcf_patterns_and_references(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    assert builder.is_remote("https://x/y.vcf.gz") and builder.is_remote("gs://b/o.bcf") and not builder.is_remote("/data/x.vcf.gz")
    served = {"https://h/c.1.vcf.gz", "https://h/c.2.vcf.gz"}
    monkeypatch.setattr(builder, "remote_exists", lambda url, timeout=20.0: url in served)
    assert builder.resolve_vcf_files("https://h/c.{chrom}.vcf.gz") == ["https://h/c.1.vcf.gz", "https://h/c.2.vcf.gz"]  # no "chr"
    assert builder.resolve_vcf_files("https://h/missing.vcf.gz") == []
    assert builder.vcf_indexed("https://h/c.1.vcf.gz")
    assert builder._cwd_for(["bcftools", "view", "https://h/c.1.vcf.gz"]) == str(tmp_path / "cache" / "genomics" / "remote_index")
    assert builder._cwd_for(["bcftools", "view", "/local.vcf.gz"]) is None
    assert builder.reference_fai("/refs/GRCh38.fa") == Path("/refs/GRCh38.fa.fai")
    cached = tmp_path / "cache" / "genomics" / "remote_index" / "GRCh38.fa.fai"
    cached.write_text("chr1\t100\t6\t60\t61\n")
    assert builder.reference_fai("https://h/GRCh38.fa") == cached  # downloaded once, then reused


def test_vcf_overrides_take_precedence_over_the_pattern(monkeypatch, tmp_path):
    for name in ("cohort.chr1.vcf.gz", "cohort.chrX.v2.vcf.gz"):
        (tmp_path / name).write_bytes(b"")
    monkeypatch.setattr(builder, "vcf_indexed", lambda path: True)
    monkeypatch.setattr(builder, "vcf_contigs", lambda path: ["chr1", "chrX"])
    importer = builder.DatasetImporter({
        "output_dir": str(tmp_path / "out"), "vcf": str(tmp_path / "cohort.{chrom}.vcf.gz"), "reference_fasta": "https://h/ref.fa",
        "vcf_overrides": {"chrX": str(tmp_path / "cohort.chrX.v2.vcf.gz")},
    })
    assert importer.reference == "https://h/ref.fa"  # URLs are not resolved as local paths
    assert importer._vcf_for("chrX") == (tmp_path / "cohort.chrX.v2.vcf.gz", "chrX")
    assert importer._vcf_for("chr1") == (tmp_path / "cohort.chr1.vcf.gz", "chr1")
    with pytest.raises(builder.ImportSpecError):
        importer._vcf_for("chr2")


def test_import_task_accepts_remote_sources(monkeypatch, tmp_path):
    from genomics.visualizer import launch

    monkeypatch.setattr(builder, "remote_exists", lambda url, timeout=20.0: True)
    title, steps, params, files = launch.import_task({
        "name": "HGDP quickstart", "output_dir": str(tmp_path / "hgdp"), "vcf": qs.PRESETS["hgdp_1kg"].vcf,
        "vcf_overrides": qs.PRESETS["hgdp_1kg"].vcf_overrides, "reference_fasta": qs.REMOTE_REFERENCE,
        "genes": ["OCA2"], "samples": ["HGDP1"], "metadata": {"HGDP1": {"population": "Adygei"}},
    }, alphagenome_ready=False)
    spec = json.loads(files["import_spec.json"])
    assert spec["reference_fasta"] == qs.REMOTE_REFERENCE and spec["vcf_overrides"]["chrX"].endswith(".bcf")
