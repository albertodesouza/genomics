from __future__ import annotations

from pathlib import Path
import shutil
from typing import Dict, Iterable, Iterator, List, Tuple

from datasets import Dataset, concatenate_datasets, disable_progress_bar, enable_progress_bar, load_from_disk
from tqdm import tqdm

from .features import gene_reference_features, gene_window_features, sample_features
from .fileio import (
    detect_vcf_chr_prefix,
    load_fasta_sequence,
    load_json,
    load_npz_matrix,
    load_vcf_variants,
    stable_json_digest,
    write_json,
)
from .models import BuildSummary, ConversionOptions, GeneReferenceRecord, GeneWindowRecord, SampleRecord, SourceDataset


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
    gene_reference_index: Dict[str, str] = {}
    warnings: List[str] = []
    duplicate_samples_skipped = 0
    duplicate_gene_windows_skipped = 0
    total_samples = 0
    total_gene_windows = 0
    total_gene_references = 0
    samples_shards_dir = output_dir / ".samples_shards"
    gene_windows_shards_dir = output_dir / ".gene_windows_shards"
    gene_references_shards_dir = output_dir / ".gene_references_shards"

    if samples_shards_dir.exists():
        shutil.rmtree(samples_shards_dir)
    if gene_windows_shards_dir.exists():
        shutil.rmtree(gene_windows_shards_dir)
    if gene_references_shards_dir.exists():
        shutil.rmtree(gene_references_shards_dir)
        
    samples_shards_dir.mkdir(parents=True, exist_ok=True)
    gene_windows_shards_dir.mkdir(parents=True, exist_ok=True)
    gene_references_shards_dir.mkdir(parents=True, exist_ok=True)

    sample_buffer: List[Dict] = []
    gene_window_buffer: List[Dict] = []
    gene_reference_buffer: List[Dict] = []
    samples_shard_count = 0
    gene_windows_shard_count = 0
    gene_references_shard_count = 0

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
            for sample_row, gene_window_rows, gene_reference_rows in _iter_source_dataset_rows(source, options, progress):
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

                for row in gene_reference_rows:
                    gene_symbol = row["gene_symbol"]
                    row_digest = stable_json_digest(row)
                    existing_digest = gene_reference_index.get(gene_symbol)
                    if existing_digest is None:
                        gene_reference_index[gene_symbol] = row_digest
                        gene_reference_buffer.append(row)
                        total_gene_references += 1
                        if options.chunk_size > 0 and len(gene_reference_buffer) >= options.chunk_size:
                            gene_references_shard_count = _flush_rows(
                                gene_reference_buffer,
                                gene_reference_features(),
                                gene_references_shards_dir,
                                gene_references_shard_count,
                            )
                    elif existing_digest != row_digest:
                         warnings.append(f"Conflicting reference data for gene {gene_symbol}. Using the first one encountered.")

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
    if gene_reference_buffer:
        gene_references_shard_count = _flush_rows(
            gene_reference_buffer,
            gene_reference_features(),
            gene_references_shards_dir,
            gene_references_shard_count,
        )
    if gene_window_buffer:
        gene_windows_shard_count = _flush_rows(
            gene_window_buffer,
            gene_window_features(),
            gene_windows_shards_dir,
            gene_windows_shard_count,
        )

    _merge_shards(samples_shards_dir, output_dir / "samples")
    _merge_shards(gene_references_shards_dir, output_dir / "gene_reference")
    _merge_shards(gene_windows_shards_dir, output_dir / "gene_windows")
    shutil.rmtree(samples_shards_dir)
    shutil.rmtree(gene_references_shards_dir)
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
) -> Iterator[Tuple[Dict, List[Dict], List[Dict]]]:
    dataset_metadata = load_json(source.path / "dataset_metadata.json")
    individuals = dataset_metadata.get("individuals", [])
    source_genes_seen = set()

    for sample_id in individuals:
        sample_metadata = load_json(source.path / "individuals" / sample_id / "individual_metadata.json")
        sample_record = _build_sample_record(source, sample_metadata)
        sample_row = _sample_record_to_row(sample_record)
        gene_window_rows: List[Dict] = []
        gene_reference_rows: List[Dict] = []

        for gene in sample_record.available_genes:
            gene_record = _build_gene_window_record(source, sample_id, sample_metadata, gene, options)
            gene_window_rows.append(_gene_window_record_to_row(gene_record))

            if gene not in source_genes_seen:
                ref_record = _build_gene_reference_record(source, sample_id, sample_metadata, gene, options)
                if ref_record:
                    gene_reference_rows.append(_gene_reference_record_to_row(ref_record))
                    source_genes_seen.add(gene)

        if progress is not None:
            progress.update(1)

        yield sample_row, gene_window_rows, gene_reference_rows


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

    h1_sequence = ""
    h2_sequence = ""
    if options.include_sequences:
        h1_sequence = _load_optional_sequence(window_dir / f"{sample_id}.H1.window.fixed.fa")
        h2_sequence = _load_optional_sequence(window_dir / f"{sample_id}.H2.window.fixed.fa")

    predictions_h1 = _load_predictions(window_dir / "predictions_H1")
    predictions_h2 = _load_predictions(window_dir / "predictions_H2")

    vcf_text = ""
    if options.include_vcf:
        chromosome = _as_optional_string(window_metadata.get("chromosome"))
        start = _as_optional_int(window_metadata.get("start"))
        end = _as_optional_int(window_metadata.get("end"))
        
        # Use cohort VCF if available for better filtering (genotypes)
        cohort_vcf_pattern = options.vcf_pattern
        if cohort_vcf_pattern and chromosome and start is not None and end is not None:
            vcf_path = None
            if '{chrom}' in cohort_vcf_pattern:
                path_naked = Path(cohort_vcf_pattern.replace('{chrom}', chromosome)).resolve()
                path_chr = Path(cohort_vcf_pattern.replace('{chrom}', 'chr' + chromosome.replace('chr',''))).resolve()
                if path_chr.exists():
                    vcf_path = path_chr
                elif path_naked.exists():
                    vcf_path = path_naked
            else:
                vcf_path = Path(cohort_vcf_pattern).resolve()
            
            if vcf_path and vcf_path.exists():
                from .fileio import detect_vcf_chr_prefix
                vcf_prefix = detect_vcf_chr_prefix(vcf_path)
                vcf_chrom = chromosome
                if vcf_prefix == 'chr' and not vcf_chrom.startswith('chr'):
                    vcf_chrom = 'chr' + vcf_chrom
                elif vcf_prefix == '' and vcf_chrom.startswith('chr'):
                    vcf_chrom = vcf_chrom.replace('chr', '', 1)
                
                region = f"{vcf_chrom}:{start+1}-{end}"
                vcf_text = load_vcf_variants(vcf_path, sample_id=sample_id, region=region)
        
        # Fallback to local window VCF if cohort extraction failed or wasn't requested
        if not vcf_text:
            vcf_path = window_dir / f"{sample_id}.window.consensus_ready.vcf.gz"
            if not vcf_path.exists():
                vcf_path = window_dir / f"{sample_id}.window.vcf.gz"
            if vcf_path.exists():
                vcf_text = load_vcf_variants(vcf_path, sample_id=sample_id)

    return GeneWindowRecord(
        sample_id=sample_id,
        gene=gene,
        source_dataset=source.name,
        chromosome=_as_optional_string(window_metadata.get("chromosome")),
        start=_as_optional_int(window_metadata.get("start")),
        end=_as_optional_int(window_metadata.get("end")),
        window_size=_as_optional_int(window_metadata.get("window_size")),
        window_type=_as_optional_string(window_metadata.get("type")),
        h1_sequence=h1_sequence,
        h2_sequence=h2_sequence,
        vcf_text=vcf_text,
        outputs=[str(value) for value in window_metadata.get("outputs", [])],
        ontologies=[str(value) for value in window_metadata.get("ontologies", [])],
        predictions_h1=predictions_h1,
        predictions_h2=predictions_h2,
        source_path=str(source.path),
        source_sample_path=str(source.path / "individuals" / sample_id),
        source_window_path=str(window_dir),
        source_metadata_path=str(source.path / "individuals" / sample_id / "individual_metadata.json"),
    )


def _build_gene_reference_record(
    source: SourceDataset,
    sample_id: str,
    sample_metadata: Dict,
    gene: str,
    options: ConversionOptions,
) -> GeneReferenceRecord | None:
    window_metadata = sample_metadata.get("window_metadata", {}).get(gene, {})
    window_dir = source.path / "individuals" / sample_id / "windows" / gene
    if not window_dir.exists():
        return None

    ref_sequence = _load_optional_sequence(window_dir / "ref.window.fa")
    
    chromosome = _as_optional_string(window_metadata.get("chromosome"))
    start = _as_optional_int(window_metadata.get("start"))
    end = _as_optional_int(window_metadata.get("end"))
    
    return GeneReferenceRecord(
        gene_symbol=gene,
        chromosome=chromosome or "",
        start=start if start is not None else -1,
        end=end if end is not None else -1,
        ref_sequence=ref_sequence,
    )


def _load_predictions(predictions_dir: Path) -> Dict[str, bytes]:
    payload = {name: b"" for name in SUPPORTED_OUTPUTS}
    if not predictions_dir.exists():
        return payload

    for npz_path in sorted(predictions_dir.glob("*.npz")):
        output_name = npz_path.stem.lower()
        if output_name in payload:
            try:
                with open(npz_path, "rb") as f:
                    payload[output_name] = f.read()
            except Exception as e:
                print(f"[WARN] Error reading {npz_path}: {e}")
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
        "h1_sequence": record.h1_sequence or "",
        "h2_sequence": record.h2_sequence or "",
        "vcf_text": record.vcf_text or "",
        "outputs": record.outputs,
        "ontologies": record.ontologies,
        "source_path": record.source_path,
        "source_sample_path": record.source_sample_path,
        "source_window_path": record.source_window_path,
        "source_metadata_path": record.source_metadata_path,
        "predictions_h1": record.predictions_h1,
        "predictions_h2": record.predictions_h2,
    }


def _gene_reference_record_to_row(record: GeneReferenceRecord) -> Dict:
    return {
        "gene_symbol": record.gene_symbol,
        "chromosome": record.chromosome,
        "start": record.start,
        "end": record.end,
        "ref_sequence": record.ref_sequence,
    }


def _build_manifest(summary: BuildSummary, sources: List[SourceDataset]) -> Dict:
    return {
        "schema_version": 1,
        "layout": "normalized",
        "tables": ["samples", "gene_windows", "gene_reference"],
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
