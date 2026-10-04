"""Public phased cohorts that can be imported without preparing any file.

Each preset names a phased GRCh38 VCF/BCF set on a public server and its sample metadata. The
variants are read remotely: bcftools fetches only each window's region (see
``builder.remote_index_dir``), so nothing chromosome-sized is downloaded. A local copy is used
instead when one is found (e.g. the 1000 Genomes panel next to the canonical dataset).

:func:`prepare` downloads the metadata once (cached), normalises it to a flat table, picks a
population-balanced subset of unrelated samples and writes it as a TSV for the import form; the
import itself is the ordinary ``vcf-import`` job.
"""
from __future__ import annotations

import csv
import io
import json
import os
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from genomics.workflows.dataset_builders.vcf_import.builder import ImportSpecError, inspect_vcf, resolve_vcf_files
from genomics.workflows.dataset_builders.vcf_import.metadata import parse_table, records_by_sample

EBI_1KG = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
GNOMAD_HGDP = "https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2"
REMOTE_REFERENCE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
DEFAULT_GENES = ("OCA2", "HERC2", "SLC24A5", "SLC45A2", "TYR", "MC1R")
DEFAULT_PER_POPULATION = 4
_SEX = {"1": "Male", "2": "Female", "XY": "Male", "XX": "Female"}
# Selection flags: dropped from the table when every selected sample has the same value.
_FLAG_COLUMNS = ("related", "high_quality", "project")


@dataclass(frozen=True)
class Preset:
    id: str
    title: str
    name: str  # default dataset name
    summary: str
    vcf: str  # remote {chrom} pattern
    vcf_overrides: Dict[str, str]
    metadata_url: str
    metadata_file: str
    homepage: str
    citation: str
    samples: int
    populations: int
    phasing: str
    parse: Callable[[bytes], List[Dict[str, Any]]]
    scopes: Tuple[Tuple[str, str, int, int], ...] = ()  # (id, label, samples, populations); first = default

    def as_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "title": self.title,
            "name": self.name,
            "summary": self.summary,
            "homepage": self.homepage,
            "citation": self.citation,
            "samples": self.samples,
            "populations": self.populations,
            "phasing": self.phasing,
            "scopes": [{"id": s, "label": label, "samples": n, "populations": pops} for s, label, n, pops in self.scopes],
            "default_genes": list(DEFAULT_GENES),
            "default_per_population": DEFAULT_PER_POPULATION,
        }


def _parse_1kg(raw: bytes) -> List[Dict[str, Any]]:
    """``20130606_g1k_3202_samples_ped_population.txt``: whitespace-separated pedigree."""
    lines = [line.split() for line in raw.decode("utf-8").splitlines() if line.strip()]
    header, rows = lines[0], [dict(zip(lines[0], values)) for values in lines[1:]]
    if "SampleID" not in header:
        raise ImportSpecError("Unexpected 1000 Genomes pedigree header: " + " ".join(header))
    ids = {r["SampleID"] for r in rows}
    out = []
    for r in rows:
        parents = [p for p in (r.get("FatherID"), r.get("MotherID")) if p and p != "0"]
        out.append({
            "sample_id": r["SampleID"],
            "family_id": r.get("FamilyID") or r["SampleID"],
            "sex": _SEX.get(r.get("Sex", ""), None),
            "population": r.get("Population"),
            "superpopulation": r.get("Superpopulation"),
            # Children of trios are the related samples (their parents stay in the unrelated set).
            "related": any(p in ids for p in parents),
        })
    return out


def _parse_gnomad(raw: bytes) -> List[Dict[str, Any]]:
    """gnomAD v3.1.2 HGDP + 1KG sample metadata: a TSV whose columns hold JSON structs."""
    import gzip

    text = gzip.decompress(raw).decode("utf-8") if raw[:2] == b"\x1f\x8b" else raw.decode("utf-8")
    csv.field_size_limit(1 << 30)
    out = []
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        try:
            meta = json.loads(row.get("hgdp_tgp_meta") or "null") or {}
            sex = json.loads(row.get("gnomad_sex_imputation") or "null") or {}
            related = json.loads(row.get("relatedness_inference") or "null") or {}
        except ValueError:
            continue
        project = meta.get("project")
        if project not in ("HGDP", "1000 Genomes"):
            continue  # gnomAD's synthetic QC sample
        out.append({
            "sample_id": row["s"],
            "project": project,
            "population": meta.get("population"),
            "genetic_region": meta.get("genetic_region"),
            "study_region": meta.get("study_region"),
            "latitude": meta.get("latitude"),
            "longitude": meta.get("longitude"),
            "sex": _SEX.get(str(sex.get("sex_karyotype") or ""), None),
            "related": bool(related.get("related")),
            "high_quality": str(row.get("high_quality")).lower() == "true",
        })
    return out


PRESETS: Dict[str, Preset] = {
    p.id: p
    for p in (
        Preset(
            id="1kg",
            title="1000 Genomes · high coverage",
            name="1000 Genomes quickstart",
            summary="3,202 samples from 26 populations sequenced at 30x by the New York Genome Center; the source of the canonical dataset.",
            vcf=f"{EBI_1KG}/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{{chrom}}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
            vcf_overrides={"chrX": f"{EBI_1KG}/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"},
            metadata_url=f"{EBI_1KG}/20130606_g1k_3202_samples_ped_population.txt",
            metadata_file="1kg_3202_samples_ped_population.txt",
            homepage="https://www.internationalgenome.org/data-portal/data-collection/30x-grch38",
            citation="Byrska-Bishop et al., Cell 2022",
            samples=3202,
            populations=26,
            phasing="statistical (SHAPEIT2-duohmm), SNVs, indels and SVs",
            parse=_parse_1kg,
        ),
        Preset(
            id="hgdp_1kg",
            title="gnomAD · HGDP + 1000 Genomes",
            name="HGDP quickstart",
            summary="4,091 genomes jointly called by gnomAD v3.1: 925 from the Human Genome Diversity Project's 52 populations plus 1000 Genomes.",
            vcf=f"{GNOMAD_HGDP}/hgdp1kgp_{{chrom}}.filtered.SNV_INDEL.phased.shapeit5.bcf",
            vcf_overrides={"chrX": f"{GNOMAD_HGDP}/hgdp1kgp_chrX_non_par.full.shapeit5_rare.bcf"},
            metadata_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz",
            metadata_file="gnomad_v3.1.2_hgdp_1kg_sample_meta.tsv.bgz",
            homepage="https://gnomad.broadinstitute.org/news/2023-11-gnomad-hgdp-and-1000-genomes-callset/",
            citation="Koenig et al., bioRxiv 2023",
            samples=4091,
            populations=78,
            phasing="statistical (SHAPEIT5), SNVs and indels; chrX outside PAR",
            parse=_parse_gnomad,
            scopes=(("hgdp", "HGDP only", 925, 52), ("all", "HGDP + 1000 Genomes", 4091, 78)),
        ),
    )
}


def get_preset(preset_id: str) -> Preset:
    if preset_id not in PRESETS:
        raise ImportSpecError(f"Unknown quickstart dataset {preset_id!r}; choose {', '.join(PRESETS)}")
    return PRESETS[preset_id]


def presets() -> List[Dict[str, Any]]:
    return [p.as_dict() for p in PRESETS.values()]


# ------------------------------------------------------------------------------------- inputs
def local_vcf(preset: Preset) -> Optional[Tuple[str, Dict[str, str]]]:
    """A local copy of the preset's VCFs: the 1000 Genomes panel recorded by the canonical dataset."""
    if preset.id != "1kg":
        return None
    try:
        from genomics.workspace import DEFAULT_DATASET_DIR

        meta = json.loads((Path(DEFAULT_DATASET_DIR) / "dataset_metadata.json").read_text(encoding="utf-8"))
    except (OSError, ValueError, ImportError):
        return None
    pattern = str((meta.get("raw_variant_sources") or {}).get("vcf_pattern") or os.environ.get("KG1000_VCF_PATTERN") or "")
    if "{chrom}" not in pattern or not Path(pattern.replace("{chrom}", "chr22")).exists():
        return None
    overrides = {}
    chrx = pattern.replace("{chrom}", "chrX").replace(".vcf.gz", ".v2.vcf.gz")
    if not Path(pattern.replace("{chrom}", "chrX")).exists() and Path(chrx).exists():
        overrides["chrX"] = chrx
    return pattern, overrides


def reference_fasta(local: Optional[str]) -> str:
    """The local GRCh38 (same build as both presets) when present, else the 1000 Genomes copy."""
    return local or REMOTE_REFERENCE


def fetch_metadata(preset: Preset, cache_dir: Path) -> Path:
    path = Path(cache_dir) / preset.metadata_file
    if path.exists() and path.stat().st_size:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        with urllib.request.urlopen(preset.metadata_url, timeout=120) as response, open(tmp, "wb") as handle:
            shutil.copyfileobj(response, handle)
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        raise ImportSpecError(f"Could not download the sample metadata from {preset.metadata_url}: {exc}")
    os.replace(tmp, path)
    return path


def select_samples(records: List[Dict[str, Any]], available: Optional[List[str]], per_population: int, scope: Optional[str] = None) -> List[Dict[str, Any]]:
    """Unrelated samples present in the VCF, the first ``per_population`` by id in each population
    (``per_population <= 0`` keeps them all, related ones included)."""
    present = set(available) if available is not None else None
    rows = [r for r in records if present is None or r["sample_id"] in present]
    if scope == "hgdp":
        rows = [r for r in rows if r.get("project") == "HGDP"]
    if per_population <= 0:
        return sorted(rows, key=lambda r: r["sample_id"])
    by_population: Dict[str, List[Dict[str, Any]]] = {}
    for r in sorted(rows, key=lambda r: r["sample_id"]):
        if r.get("related") or r.get("high_quality") is False:
            continue
        by_population.setdefault(str(r.get("population") or "unknown"), []).append(r)
    return [r for pop in sorted(by_population) for r in by_population[pop][:per_population]]


def write_table(rows: List[Dict[str, Any]], path: Path) -> Path:
    columns: List[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    columns = [c for c in columns if c not in _FLAG_COLUMNS or len({str(r.get(c)) for r in rows}) > 1]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(["" if row.get(c) is None else str(row.get(c)).lower() if isinstance(row.get(c), bool) else row.get(c) for c in columns])
    return path


def prepare(preset_id: str, cache_dir: Path, per_population: int = DEFAULT_PER_POPULATION, scope: Optional[str] = None, reference: Optional[str] = None) -> Dict[str, Any]:
    """Everything the import form needs for a preset: VCF, reference, VCF inspection and a TSV
    with the selected samples' metadata."""
    preset = get_preset(preset_id)
    scopes = [s[0] for s in preset.scopes]
    scope = scope if scope in scopes else (scopes[0] if scopes else None)
    cache_dir = Path(cache_dir) / preset.id
    local = local_vcf(preset)
    vcf, overrides = local if local else (preset.vcf, dict(preset.vcf_overrides))
    if not resolve_vcf_files(vcf):
        raise ImportSpecError(f"Could not reach {vcf.replace('{chrom}', 'chr1')}; check the network connection")
    inspection = inspect_vcf(vcf)
    records = preset.parse(fetch_metadata(preset, cache_dir).read_bytes())
    selected = select_samples(records, inspection["samples"], int(per_population), scope)
    if not selected:
        raise ImportSpecError("No sample matched both the metadata and the VCF")
    suffix = f"{scope + '_' if scope else ''}{'all' if int(per_population) <= 0 else f'{int(per_population)}pp'}"
    table = write_table(selected, cache_dir / f"{preset.id}_{suffix}.tsv")
    populations = sorted({str(r.get("population")) for r in selected if r.get("population")})
    columns, rows = parse_table(table.read_text(encoding="utf-8"), table.name)
    records = records_by_sample(rows, "sample_id", family_column="family_id" if "family_id" in columns else None, sex_column="sex")
    name = preset.name + (" all projects" if scope == "all" else "") + (f" {int(per_population)} per population" if int(per_population) > 0 else " all samples")
    return {
        "preset": preset.as_dict(),
        "name": name,
        "vcf": vcf,
        "vcf_overrides": overrides,
        "vcf_source": "local" if local else "remote",
        "reference_fasta": reference_fasta(reference),
        "inspect": inspection,
        "metadata_path": str(table),
        "samples": [r["sample_id"] for r in selected],
        "records": records,
        "populations": len(populations),
        "scope": scope,
        "per_population": int(per_population),
        "genes": list(DEFAULT_GENES),
    }

