from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path
from typing import Iterable, List, Sequence

from genomics.core.data_registry import resolve_dataset


KG1000_HIGH_COVERAGE_BASE_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
)
KG1000_HIGH_COVERAGE_FILENAME = "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
KG1000_HIGH_COVERAGE_FILENAME_V2 = "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"


def normalize_chromosome(chromosome: str) -> str:
    chromosome = str(chromosome)
    return chromosome if chromosome.startswith("chr") else f"chr{chromosome}"


def default_chromosomes() -> List[str]:
    return [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def raw_variant_chromosome_dir(dataset_id: str = "1kg_high_coverage") -> Path:
    return resolve_dataset(dataset_id).path / "raw_variants" / "vcf_chromosomes"


def kg1000_vcf_name(chromosome: str) -> str:
    chrom = normalize_chromosome(chromosome)
    if chrom == "chrX":
        return KG1000_HIGH_COVERAGE_FILENAME_V2.format(chrom=chrom)
    return KG1000_HIGH_COVERAGE_FILENAME.format(chrom=chrom)


def kg1000_vcf_url(chromosome: str, base_url: str = KG1000_HIGH_COVERAGE_BASE_URL) -> str:
    return f"{base_url.rstrip('/')}/{kg1000_vcf_name(chromosome)}"


def _download(url: str, output: Path, *, force: bool = False) -> None:
    if output.exists() and not force:
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".partial")
    with urllib.request.urlopen(url) as response, open(partial, "wb") as f:
        shutil.copyfileobj(response, f)
    partial.replace(output)


def ensure_kg1000_vcfs(
    *,
    chromosomes: Sequence[str],
    dataset_id: str = "1kg_high_coverage",
    output_dir: Path | None = None,
    base_url: str = KG1000_HIGH_COVERAGE_BASE_URL,
    force: bool = False,
    include_index: bool = True,
) -> List[Path]:
    target_dir = output_dir.expanduser().resolve() if output_dir else raw_variant_chromosome_dir(dataset_id)
    outputs: List[Path] = []
    for chromosome in chromosomes:
        name = kg1000_vcf_name(chromosome)
        vcf_path = target_dir / name
        _download(f"{base_url.rstrip('/')}/{name}", vcf_path, force=force)
        outputs.append(vcf_path)
        if include_index:
            _download(f"{base_url.rstrip('/')}/{name}.tbi", Path(f"{vcf_path}.tbi"), force=force)
    return outputs
