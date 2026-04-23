from __future__ import annotations

from pathlib import Path
import shutil
from typing import Dict, Iterable, Iterator, List, Tuple

from datasets import Dataset, concatenate_datasets, disable_progress_bar, enable_progress_bar, load_from_disk
from tqdm import tqdm

from .features import gene_window_features, sample_features
from .fileio import load_fasta_sequence, load_json, load_npz_matrix, stable_json_digest, write_json
from .models import BuildSummary, ConversionOptions, GeneWindowRecord, SampleRecord, SourceDataset


SUPPORTED_OUTPUTS = (
    "rna_seq",
    "atac",
    "cage",
    "chip_histone",
    "chip_tf",
    "dnase",
    "procap",
)


def build_hf_dataset_store(
    source_dirs: Iterable[Path],
    output_dir: Path,
    options: ConversionOptions | None = None,
) -> BuildSummary:
    options = options or ConversionOptions()
    sources = [SourceDataset(name=path.name, path=path) for path in source_dirs]
    _validate_sources(sources)
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_names = set()
    sample_index: Dict[str, str] = {}
    gene_window_index: Dict[Tuple[str, str], str] = {}
    warnings: List[str] = []
    duplicate_samples_skipped = 0
    duplicate_gene_windows_skipped = 0
    total_samples = 0
    total_gene_windows = 0
    samples_shards_dir = output_dir / ".samples_shards"
    gene_windows_shards_dir = output_dir / ".gene_windows_shards"

    if samples_shards_dir.exists():
        shutil.rmtree(samples_shards_dir)
    if gene_windows_shards_dir.exists():
        shutil.rmtree(gene_windows_shards_dir)
    samples_shards_dir.mkdir(parents=True, exist_ok=True)
    gene_windows_shards_dir.mkdir(parents=True, exist_ok=True)

    sample_buffer: List[Dict] = []
    gene_window_buffer: List[Dict] = []
    samples_shard_count = 0
    gene_windows_shard_count = 0

    total_individuals = 0
    for source in sources:
        dataset_metadata = load_json(source.path / "dataset_metadata.json")
        individuals = list(dataset_metadata.get("individuals", []))
        total_individuals += len(individuals)

    progress = tqdm(total=total_individuals, desc="Converting samples", unit="sample")
    disable_progress_bar()
    try:
        for source in sources:
            print(f"[hf_genomics_dataset] Reading source dataset: {source.path}")
            for sample_row, gene_window_rows in _iter_source_dataset_rows(source, options, progress):
                sample_id = sample_row["sample_id"]
                row_digest = stable_json_digest(sample_row)
                existing_digest = sample_index.get(sample_id)
                if existing_digest is None:
                    sample_index[sample_id] = row_digest
                    sample_buffer.append(sample_row)
                    total_samples += 1
                    if options.chunk_size > 0 and len(sample_buffer) >= options.chunk_size:
                        samples_shard_count = _flush_rows(sample_buffer, sample_features(), samples_shards_dir, samples_shard_count)
                elif existing_digest == row_digest:
                    duplicate_samples_skipped += 1
                else:
                    raise ValueError(
                        f"Conflicting sample metadata for sample_id={sample_id}. "
                        "Resolve the source overlap before building a canonical HF dataset."
                    )

                for row in gene_window_rows:
                    key = (row["sample_id"], row["gene"])
                    row_digest = stable_json_digest(row)
                    existing_digest = gene_window_index.get(key)
                    if existing_digest is None:
                        gene_window_index[key] = row_digest
                        gene_window_buffer.append(row)
                        gene_names.add(row["gene"])
                        total_gene_windows += 1
                        if options.chunk_size > 0 and len(gene_window_buffer) >= options.chunk_size:
                            gene_windows_shard_count = _flush_rows(
                                gene_window_buffer,
                                gene_window_features(),
                                gene_windows_shards_dir,
                                gene_windows_shard_count,
                            )
                        continue
                    if existing_digest == row_digest:
                        duplicate_gene_windows_skipped += 1
                        continue

                    raise ValueError(
                        f"Conflicting gene window payload for sample_id={row['sample_id']} gene={row['gene']}. "
                        "Resolve the source overlap before building a canonical HF dataset."
                    )
    finally:
        progress.close()
        enable_progress_bar()

    if sample_buffer:
        samples_shard_count = _flush_rows(sample_buffer, sample_features(), samples_shards_dir, samples_shard_count)
    if gene_window_buffer:
        gene_windows_shard_count = _flush_rows(
            gene_window_buffer,
            gene_window_features(),
            gene_windows_shards_dir,
            gene_windows_shard_count,
        )

    _merge_shards(samples_shards_dir, output_dir / "samples")
    _merge_shards(gene_windows_shards_dir, output_dir / "gene_windows")
    shutil.rmtree(samples_shards_dir)
    shutil.rmtree(gene_windows_shards_dir)

    summary = BuildSummary(
        source_datasets=[source.name for source in sources],
        total_samples=total_samples,
        total_gene_windows=total_gene_windows,
        genes=sorted(gene_names),
        duplicate_samples_skipped=duplicate_samples_skipped,
        duplicate_gene_windows_skipped=duplicate_gene_windows_skipped,
        warnings=warnings,
    )
    write_json(output_dir / "manifest.json", _build_manifest(summary, sources))
    return summary


def _iter_source_dataset_rows(
    source: SourceDataset,
    options: ConversionOptions,
    progress: tqdm | None = None,
) -> Iterator[Tuple[Dict, List[Dict]]]:
    dataset_metadata = load_json(source.path / "dataset_metadata.json")
    individuals = dataset_metadata.get("individuals", [])

    for sample_id in individuals:
        sample_metadata = load_json(source.path / "individuals" / sample_id / "individual_metadata.json")
        sample_record = _build_sample_record(source, sample_metadata)
        sample_row = _sample_record_to_row(sample_record)
        gene_window_rows: List[Dict] = []

        for gene in sample_record.available_genes:
            gene_record = _build_gene_window_record(source, sample_id, sample_metadata, gene, options)
            gene_window_rows.append(_gene_window_record_to_row(gene_record))

        if progress is not None:
            progress.update(1)

        yield sample_row, gene_window_rows


def _validate_sources(sources: List[SourceDataset]) -> None:
    if not sources:
        raise ValueError("At least one source dataset is required")

    seen_names = set()
    for source in sources:
        if source.name in seen_names:
            raise ValueError(f"Duplicate source dataset name: {source.name}")
        seen_names.add(source.name)

        if not source.path.exists():
            raise FileNotFoundError(f"Source dataset not found: {source.path}")

        metadata_path = source.path / "dataset_metadata.json"
        individuals_dir = source.path / "individuals"
        if not metadata_path.exists():
            raise FileNotFoundError(f"dataset_metadata.json not found in {source.path}")
        if not individuals_dir.exists():
            raise FileNotFoundError(f"individuals directory not found in {source.path}")


def _build_sample_record(source: SourceDataset, sample_metadata: Dict) -> SampleRecord:
    sample_id = sample_metadata["sample_id"]
    return SampleRecord(
        sample_id=sample_id,
        family_id=_as_optional_string(sample_metadata.get("family_id")),
        sex=_as_optional_int(sample_metadata.get("sex")),
        population=_as_optional_string(sample_metadata.get("population")),
        superpopulation=_as_optional_string(sample_metadata.get("superpopulation")),
        longevity=sample_metadata.get("longevity"),
        frog_likelihood=[float(value) for value in sample_metadata.get("frog_likelihood", [])],
        frog_population_names=[str(value) for value in sample_metadata.get("frog_population_names", [])],
        available_genes=[str(gene) for gene in sample_metadata.get("windows", [])],
        created_at=_as_optional_string(sample_metadata.get("created_at")),
        last_updated=_as_optional_string(sample_metadata.get("last_updated")),
        source_dataset=source.name,
        source_path=str(source.path),
        source_sample_path=str(source.path / "individuals" / sample_id),
    )


def _build_gene_window_record(
    source: SourceDataset,
    sample_id: str,
    sample_metadata: Dict,
    gene: str,
    options: ConversionOptions,
) -> GeneWindowRecord:
    window_metadata = sample_metadata.get("window_metadata", {}).get(gene, {})
    window_dir = source.path / "individuals" / sample_id / "windows" / gene
    if not window_dir.exists():
        raise FileNotFoundError(f"Window directory not found for sample_id={sample_id} gene={gene}: {window_dir}")

    ref_sequence = ""
    h1_sequence = ""
    h2_sequence = ""
    if options.include_sequences:
        ref_sequence = _load_optional_sequence(window_dir / "ref.window.fa")
        h1_sequence = _load_optional_sequence(window_dir / f"{sample_id}.H1.window.fixed.fa")
        h2_sequence = _load_optional_sequence(window_dir / f"{sample_id}.H2.window.fixed.fa")

    predictions_h1 = _load_predictions(window_dir / "predictions_H1")
    predictions_h2 = _load_predictions(window_dir / "predictions_H2")

    return GeneWindowRecord(
        sample_id=sample_id,
        gene=gene,
        source_dataset=source.name,
        chromosome=_as_optional_string(window_metadata.get("chromosome")),
        start=_as_optional_int(window_metadata.get("start")),
        end=_as_optional_int(window_metadata.get("end")),
        window_size=_as_optional_int(window_metadata.get("window_size")),
        window_type=_as_optional_string(window_metadata.get("type")),
        ref_sequence=ref_sequence,
        h1_sequence=h1_sequence,
        h2_sequence=h2_sequence,
        outputs=[str(value) for value in window_metadata.get("outputs", [])],
        ontologies=[str(value) for value in window_metadata.get("ontologies", [])],
        predictions_h1=predictions_h1,
        predictions_h2=predictions_h2,
        source_path=str(source.path),
        source_sample_path=str(source.path / "individuals" / sample_id),
        source_window_path=str(window_dir),
        source_metadata_path=str(source.path / "individuals" / sample_id / "individual_metadata.json"),
    )


def _load_predictions(predictions_dir: Path) -> Dict[str, List[List[float]]]:
    payload = {name: [] for name in SUPPORTED_OUTPUTS}
    if not predictions_dir.exists():
        return payload

    for npz_path in sorted(predictions_dir.glob("*.npz")):
        output_name = npz_path.stem.lower()
        if output_name not in payload:
            continue
        payload[output_name] = load_npz_matrix(npz_path)

    return payload


def _sample_record_to_row(record: SampleRecord) -> Dict:
    return {
        "sample_id": record.sample_id,
        "family_id": record.family_id or "",
        "sex": record.sex if record.sex is not None else -1,
        "population": record.population or "",
        "superpopulation": record.superpopulation or "",
        "longevity": bool(record.longevity) if record.longevity is not None else False,
        "frog_likelihood": record.frog_likelihood,
        "frog_population_names": record.frog_population_names,
        "available_genes": record.available_genes,
        "created_at": record.created_at or "",
        "last_updated": record.last_updated or "",
        "source_dataset": record.source_dataset,
        "source_path": record.source_path,
        "source_sample_path": record.source_sample_path,
    }


def _gene_window_record_to_row(record: GeneWindowRecord) -> Dict:
    return {
        "sample_id": record.sample_id,
        "gene": record.gene,
        "source_dataset": record.source_dataset,
        "chromosome": record.chromosome or "",
        "start": record.start if record.start is not None else -1,
        "end": record.end if record.end is not None else -1,
        "window_size": record.window_size if record.window_size is not None else -1,
        "window_type": record.window_type or "",
        "ref_sequence": record.ref_sequence or "",
        "h1_sequence": record.h1_sequence or "",
        "h2_sequence": record.h2_sequence or "",
        "outputs": record.outputs,
        "ontologies": record.ontologies,
        "source_path": record.source_path,
        "source_sample_path": record.source_sample_path,
        "source_window_path": record.source_window_path,
        "source_metadata_path": record.source_metadata_path,
        "predictions_h1": record.predictions_h1,
        "predictions_h2": record.predictions_h2,
    }


def _build_manifest(summary: BuildSummary, sources: List[SourceDataset]) -> Dict:
    return {
        "schema_version": 1,
        "layout": "normalized",
        "tables": ["samples", "gene_windows"],
        "source_datasets": [
            {
                "name": source.name,
                "path": str(source.path),
            }
            for source in sources
        ],
        "summary": {
            "total_samples": summary.total_samples,
            "total_gene_windows": summary.total_gene_windows,
            "genes": summary.genes,
            "duplicate_samples_skipped": summary.duplicate_samples_skipped,
            "duplicate_gene_windows_skipped": summary.duplicate_gene_windows_skipped,
            "warnings": summary.warnings,
        },
    }


def _load_optional_sequence(path: Path) -> str:
    if not path.exists():
        return ""
    return load_fasta_sequence(path)


def _flush_rows(rows: List[Dict], features, shards_dir: Path, shard_count: int) -> int:
    if not rows:
        return shard_count

    shard_dir = shards_dir / f"shard_{shard_count:05d}"
    Dataset.from_list(rows, features=features).save_to_disk(str(shard_dir))
    rows.clear()
    return shard_count + 1


def _merge_shards(shards_dir: Path, output_dir: Path) -> None:
    shard_paths = sorted(path for path in shards_dir.iterdir() if path.is_dir())
    if not shard_paths:
        Dataset.from_list([]).save_to_disk(str(output_dir))
        return

    datasets = [load_from_disk(str(path)) for path in shard_paths]
    merged = datasets[0] if len(datasets) == 1 else concatenate_datasets(datasets)
    merged.save_to_disk(str(output_dir))


def _as_optional_string(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _as_optional_int(value: object) -> int | None:
    if value is None:
        return None
    return int(value)
