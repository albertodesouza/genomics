from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from genomics.core.data_registry import data_root


GRCH38_FULL_ANALYSIS_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
)
GRCH38_FULL_ANALYSIS_NAME = "GRCh38_full_analysis_set_plus_decoy_hla"


@dataclass(frozen=True)
class ReferenceRef:
    reference_id: str
    path: Path
    url: str
    metadata_path: Path
    fai_path: Path


def references_root() -> Path:
    return data_root() / "references"


def grch38_full_analysis_ref() -> ReferenceRef:
    path = references_root() / f"{GRCH38_FULL_ANALYSIS_NAME}.fa"
    return ReferenceRef(
        reference_id=GRCH38_FULL_ANALYSIS_NAME,
        path=path,
        url=GRCH38_FULL_ANALYSIS_URL,
        metadata_path=references_root() / f"{GRCH38_FULL_ANALYSIS_NAME}.metadata.json",
        fai_path=Path(f"{path}.fai"),
    )


def registered_reference_ids() -> tuple[str, ...]:
    return (GRCH38_FULL_ANALYSIS_NAME,)


def resolve_reference(reference_id: str = GRCH38_FULL_ANALYSIS_NAME) -> ReferenceRef:
    if reference_id != GRCH38_FULL_ANALYSIS_NAME:
        known = ", ".join(registered_reference_ids())
        raise KeyError(f"Reference id desconhecido: {reference_id}. Conhecidos: {known}")
    return grch38_full_analysis_ref()


def _download(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output.with_suffix(output.suffix + ".tmp")
    with urllib.request.urlopen(url) as response, open(tmp_path, "wb") as f:
        shutil.copyfileobj(response, f)
    tmp_path.replace(output)


def _decompress_gzip(source: Path, output: Path) -> None:
    tmp_path = output.with_suffix(output.suffix + ".tmp")
    with gzip.open(source, "rb") as src, open(tmp_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    tmp_path.replace(output)


def ensure_grch38_full_analysis(
    *,
    output_path: Optional[Path] = None,
    url: str = GRCH38_FULL_ANALYSIS_URL,
    force: bool = False,
    index: bool = True,
) -> ReferenceRef:
    ref = grch38_full_analysis_ref()
    if output_path is not None:
        output_path = output_path.expanduser().resolve()
        ref = ReferenceRef(
            reference_id=ref.reference_id,
            path=output_path,
            url=url,
            metadata_path=output_path.with_suffix(output_path.suffix + ".metadata.json"),
            fai_path=Path(f"{output_path}.fai"),
        )

    ref.path.parent.mkdir(parents=True, exist_ok=True)
    if force or not ref.path.exists():
        download_target = ref.path
        if url.endswith(".gz") and not ref.path.name.endswith(".gz"):
            download_target = ref.path.with_suffix(ref.path.suffix + ".gz")
        _download(url, download_target)
        if download_target != ref.path:
            _decompress_gzip(download_target, ref.path)

    if index and (force or not ref.fai_path.exists()):
        try:
            subprocess.run(["samtools", "faidx", str(ref.path)], check=True)
        except FileNotFoundError as exc:
            raise RuntimeError("samtools nao encontrado; instale samtools ou use --skip-index") from exc

    payload = {
        "reference_id": ref.reference_id,
        "path": str(ref.path),
        "url": url,
        "fai_path": str(ref.fai_path),
        "fai_exists": ref.fai_path.exists(),
    }
    with open(ref.metadata_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    return ref
