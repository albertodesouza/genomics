from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Mapping, Optional, Sequence, Tuple

import yaml


DEFAULT_WINDOW_BP = 1_048_576
DEFAULT_STRIDE_BP = 524_288
DEFAULT_ONTOLOGY_TERMS = ["CL:1000458", "CL:0000346", "CL:2000092"]
DEFAULT_OUTPUTS = ["RNA_SEQ", "CAGE", "ATAC", "DNASE", "CHIP_HISTONE", "CHIP_TF"]


@dataclass(frozen=True)
class WindowSpec:
    index: int
    start: int
    end: int
    keep_start: int
    keep_end: int


@dataclass(frozen=True)
class VariantAlleles:
    position0: int
    h1: str
    h2: str


@dataclass(frozen=True)
class PredictionTask:
    sample_id: str
    haplotype: str
    strand: str
    window: WindowSpec
    sequence: str


def _load_yaml(path: Path) -> Mapping[str, object]:
    with open(path, "r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Config YAML invalido: {path}")
    return payload


def _as_mapping(payload: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = payload.get(key, {})
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"Config section '{key}' deve ser um mapa")
    return value


def _as_list(value: object, default: Sequence[str]) -> List[str]:
    if value is None:
        return list(default)
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    if isinstance(value, list):
        return [str(item) for item in value]
    raise ValueError("Valor de lista invalido")


def _resolve_config_path(value: object, *, config_path: Path, label: str) -> Path:
    if value is None:
        raise ValueError(f"{label} nao configurado")
    path_text = str(value)
    if path_text.startswith("/path/to/"):
        raise FileNotFoundError(
            f"{label} ainda usa placeholder no YAML: {path_text}. "
            "Passe o caminho real pela CLI ou edite o config."
        )
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = (config_path.parent / path).resolve()
    if not path.exists():
        raise FileNotFoundError(f"{label} nao encontrado: {path}")
    return path


def _resolve_reference_fasta(value: object, *, config_path: Path) -> Path:
    from genomics.core.reference_registry import grch38_full_analysis_ref

    if value is None or str(value).startswith("/path/to/"):
        ref = grch38_full_analysis_ref()
        if ref.path.exists():
            return ref.path
        raise FileNotFoundError(
            "reference.fasta nao configurado e a referencia canonica nao existe. "
            f"Esperado: {ref.path}. Rode: genomics references ensure-grch38"
        )
    return _resolve_config_path(value, config_path=config_path, label="reference.fasta")


def _resolve_dataset_chromosome_vcf(dataset_id: str, chromosome: str) -> Path:
    from genomics.core.data_registry import resolve_dataset

    dataset_dir = resolve_dataset(dataset_id).path
    chrom = chromosome if chromosome.startswith("chr") else f"chr{chromosome}"
    candidates = sorted((dataset_dir / "raw_variants" / "vcf_chromosomes").glob(f"*{chrom}*.vcf.gz"))
    if not candidates:
        raise FileNotFoundError(
            f"VCF do {chromosome} nao encontrado para dataset_id={dataset_id} em "
            f"{dataset_dir / 'raw_variants' / 'vcf_chromosomes'}"
        )
    exact = [path for path in candidates if f".{chrom}." in path.name]
    return exact[0] if exact else candidates[0]


def _read_fasta_chromosome(path: Path, chromosome: str) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    wanted = chromosome[3:] if chromosome.startswith("chr") else f"chr{chromosome}"
    aliases = {chromosome, wanted}
    chunks: List[str] = []
    active = False
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                active = name in aliases
                if chunks and not active:
                    break
                continue
            if active:
                chunks.append(line.strip().upper())
    if not chunks:
        raise ValueError(f"Cromossomo {chromosome} nao encontrado em {path}")
    return "".join(chunks)


def _window_plan(chromosome_length: int, window_bp: int, stride_bp: int) -> List[WindowSpec]:
    if window_bp <= 0 or stride_bp <= 0:
        raise ValueError("window_bp e stride_bp devem ser positivos")
    if stride_bp > window_bp:
        raise ValueError("stride_bp nao deve ser maior que window_bp")
    if chromosome_length <= 0:
        raise ValueError("chromosome_length deve ser positivo")

    starts = list(range(0, max(chromosome_length - window_bp + 1, 1), stride_bp))
    last_start = max(chromosome_length - window_bp, 0)
    if starts[-1] != last_start:
        starts.append(last_start)

    centers = [start + min(window_bp, chromosome_length - start) / 2.0 for start in starts]
    windows: List[WindowSpec] = []
    for i, start in enumerate(starts):
        end = min(start + window_bp, chromosome_length)
        keep_start = 0 if i == 0 else int(round((centers[i - 1] + centers[i]) / 2.0))
        keep_end = chromosome_length if i == len(starts) - 1 else int(round((centers[i] + centers[i + 1]) / 2.0))
        windows.append(WindowSpec(i, start, end, keep_start, keep_end))
    return windows


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path, "r", encoding="utf-8")


def _read_vcf_samples(vcf_path: Path) -> List[str]:
    with _open_text(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return fields[9:]
    raise ValueError(f"Cabecalho #CHROM nao encontrado em {vcf_path}")


def _select_samples(vcf_samples: Sequence[str], requested: object) -> List[str]:
    if requested is None or requested == "all":
        return list(vcf_samples)
    if isinstance(requested, str):
        path = Path(requested)
        if path.exists():
            return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        return [item.strip() for item in requested.split(",") if item.strip()]
    if isinstance(requested, list):
        return [str(item) for item in requested]
    raise ValueError("variants.samples deve ser 'all', lista, CSV ou arquivo")


def _parse_gt(gt_field: str, require_phased: bool) -> Optional[Tuple[int, int]]:
    gt = gt_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    if "|" in gt:
        parts = gt.split("|")
    elif "/" in gt:
        if require_phased:
            return None
        parts = gt.split("/")
    else:
        return None
    if len(parts) != 2 or "." in parts:
        return None
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        return None


def _load_sample_snvs(
    vcf_path: Path,
    chromosome: str,
    sample_id: str,
    *,
    require_phased: bool,
    include_indels: bool,
) -> List[VariantAlleles]:
    if include_indels:
        raise NotImplementedError(
            "include_indels=true esta reservado no config, mas a primeira versao "
            "implementa apenas SNVs coordenada-preservada."
        )

    sample_names = _read_vcf_samples(vcf_path)
    try:
        sample_col = sample_names.index(sample_id) + 9
    except ValueError as exc:
        raise ValueError(f"Sample {sample_id} nao encontrado no VCF") from exc

    aliases = {chromosome, chromosome[3:] if chromosome.startswith("chr") else f"chr{chromosome}"}
    variants: List[VariantAlleles] = []
    with _open_text(vcf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= sample_col or fields[0] not in aliases:
                continue
            ref = fields[3].upper()
            alts = [alt.upper() for alt in fields[4].split(",")]
            if len(ref) != 1 or any(len(alt) != 1 for alt in alts):
                continue
            gt = _parse_gt(fields[sample_col], require_phased=require_phased)
            if gt is None or gt == (0, 0):
                continue

            def allele_base(index: int) -> str:
                if index == 0:
                    return ref
                if 1 <= index <= len(alts):
                    return alts[index - 1]
                return ref

            h1 = allele_base(gt[0])
            h2 = allele_base(gt[1])
            if h1 == ref and h2 == ref:
                continue
            variants.append(VariantAlleles(position0=int(fields[1]) - 1, h1=h1, h2=h2))
    return variants


def _build_haplotype_window(reference: str, window: WindowSpec, variants: Sequence[VariantAlleles], haplotype: str) -> str:
    seq = bytearray(reference[window.start : window.end].encode("ascii"))
    attr = "h1" if haplotype == "H1" else "h2"
    for variant in variants:
        if window.start <= variant.position0 < window.end:
            seq[variant.position0 - window.start] = ord(getattr(variant, attr))
    return seq.decode("ascii")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def _import_alphagenome(alphagenome_research_dir: Optional[Path]):
    if alphagenome_research_dir is not None:
        src = alphagenome_research_dir / "src"
        if src.exists() and str(src) not in sys.path:
            sys.path.insert(0, str(src))
    try:
        from alphagenome_research.model import dna_model  # type: ignore
    except ModuleNotFoundError as exc:
        if exc.name != "alphagenome_research":
            raise
        configured = str(alphagenome_research_dir) if alphagenome_research_dir else "nao configurado"
        raise ModuleNotFoundError(
            "alphagenome_research nao esta importavel neste ambiente. "
            f"Diretorio configurado: {configured}. "
            "Na DGX, rode: cd ~/I2CA && git clone https://github.com/google-deepmind/alphagenome_research.git "
            "&& python -m pip install -e ./alphagenome_research"
        ) from exc

    return dna_model


def _output_types(dna_model, names: Sequence[str]):
    output_types = []
    for name in names:
        normalized = name.upper().replace("-", "_").replace(" ", "_")
        try:
            output_types.append(getattr(dna_model.OutputType, normalized))
        except AttributeError as exc:
            raise ValueError(f"Output AlphaGenome desconhecido: {name}") from exc
    return output_types


def _create_model(dna_model, alphagenome_cfg: Mapping[str, object]):
    checkpoint_path = alphagenome_cfg.get("checkpoint_path")
    checkpoint_source = str(alphagenome_cfg.get("checkpoint_source", "huggingface"))
    model_version = str(alphagenome_cfg.get("model_version", "all_folds"))
    if checkpoint_path:
        return dna_model.create(str(checkpoint_path))
    if checkpoint_source == "kaggle":
        return dna_model.create_from_kaggle(model_version)
    if checkpoint_source == "huggingface":
        return dna_model.create_from_huggingface(model_version)
    raise ValueError("alphagenome.checkpoint_source deve ser 'huggingface', 'kaggle' ou use checkpoint_path")


def _save_output_npz(output_dir: Path, outputs, output_names: Sequence[str]) -> List[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    saved: List[str] = []
    for output_name in output_names:
        attr_name = output_name.lower()
        track_data = getattr(outputs, attr_name, None)
        if track_data is None or not hasattr(track_data, "values"):
            continue
        npz_path = output_dir / f"{attr_name}.npz"
        import numpy as np

        np.savez_compressed(npz_path, values=track_data.values)
        metadata = getattr(track_data, "metadata", None)
        if metadata is not None:
            meta_path = output_dir / f"{attr_name}_metadata.json"
            if hasattr(metadata, "to_dict"):
                payload = metadata.to_dict(orient="records")
            else:
                payload = str(metadata)
            with open(meta_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
        saved.append(attr_name)
    return saved


def _save_prediction_arrays(output_dir: Path, predictions: Mapping[str, object], metadata: Mapping[str, object]) -> List[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    saved: List[str] = []
    import numpy as np

    for output_name, values in predictions.items():
        attr_name = output_name.lower()
        np.savez_compressed(output_dir / f"{attr_name}.npz", values=values)
        track_metadata = metadata.get(output_name)
        if track_metadata is not None:
            if hasattr(track_metadata, "to_dict"):
                payload = track_metadata.to_dict(orient="records")
            else:
                payload = str(track_metadata)
            with open(output_dir / f"{attr_name}_metadata.json", "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
        saved.append(attr_name)
    return saved


def _predict_sequences_batched(dna_model, model, sequences: Sequence[str], requested_outputs, ontology_terms):
    import numpy as np
    import jax
    from alphagenome import tensor_utils
    from alphagenome.data import ontology
    from alphagenome.models import dna_model as ag_dna_model
    from alphagenome_research.model import dna_model as ag_research_dna_model
    from alphagenome_research.model.metadata import metadata as metadata_lib

    organism = ag_dna_model.Organism.HOMO_SAPIENS
    requested_outputs = tuple(set(requested_outputs))
    if ontology_terms is not None:
        ontology_terms = set(ontology.from_curie(o) if isinstance(o, str) else o for o in ontology_terms)
    ag_metadata = model._metadata[organism]
    track_masks = metadata_lib.create_track_masks(
        ag_metadata,
        requested_outputs=requested_outputs,
        requested_ontologies=ontology_terms,
    )

    with model._device_context as device, jax.transfer_guard("disallow"):
        sequence_batch = np.stack([np.asarray(model._one_hot_encoder.encode(sequence)) for sequence in sequences], axis=0)
        organism_index = np.full((len(sequences),), ag_research_dna_model.convert_to_organism_index(organism), dtype=np.int32)
        predictions = model._predict(
            model._params,
            model._state,
            jax.device_put(sequence_batch, device),
            jax.device_put(organism_index, device),
            requested_outputs=requested_outputs,
            negative_strand_mask=jax.device_put(np.zeros((len(sequences),), dtype=bool), device),
            strand_reindexing=jax.device_put(ag_metadata.strand_reindexing, device),
        )
        predictions = ag_research_dna_model._filter_predictions(
            predictions,
            track_masks=jax.device_put(track_masks, device),
        )
        predictions = jax.tree.map(tensor_utils.upcast_floating, predictions)
        predictions = jax.device_get(predictions)

    metadata_by_name = {}
    for output_type in requested_outputs:
        track_metadata = ag_metadata.get(output_type)
        if track_metadata is not None:
            metadata_by_name[output_type.name] = track_metadata[track_masks[output_type]]

    rows = []
    for row_index in range(len(sequences)):
        row = {}
        for output_type in requested_outputs:
            values = predictions.get(output_type)
            if values is not None and not isinstance(values, dict):
                row[output_type.name] = values[row_index]
        rows.append(row)
    return rows, metadata_by_name


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def _append_jsonl(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(payload, sort_keys=True) + "\n")


def _marker_path(root: Path, sample_id: str, haplotype: str, strand: str, window_index: int) -> Path:
    return root / "markers" / sample_id / haplotype / strand / f"window_{window_index:04d}.ok"


def _prediction_path(root: Path, sample_id: str, haplotype: str, strand: str, window_index: int) -> Path:
    return root / "predictions" / "by_sample" / sample_id / haplotype / strand / f"window_{window_index:04d}"


def _shard(items: Sequence[str], shard_index: Optional[int], num_shards: Optional[int]) -> List[str]:
    if shard_index is None and num_shards is None:
        return list(items)
    if shard_index is None or num_shards is None:
        raise ValueError("Use shard_index e num_shards juntos")
    if num_shards <= 0 or not (0 <= shard_index < num_shards):
        raise ValueError("Shard invalido")
    return [item for i, item in enumerate(items) if i % num_shards == shard_index]


def _chunks(items: Sequence[PredictionTask], size: int) -> Sequence[Sequence[PredictionTask]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def run(
    config_path: Path,
    *,
    sample_override: Optional[str] = None,
    ref_fasta_override: Optional[Path] = None,
    vcf_override: Optional[Path] = None,
    output_dir_override: Optional[Path] = None,
    outputs_override: Optional[str] = None,
    shard_index: Optional[int] = None,
    num_shards: Optional[int] = None,
    window_shard_index: Optional[int] = None,
    num_window_shards: Optional[int] = None,
    max_windows: Optional[int] = None,
    haplotype_filter: Optional[str] = None,
    strand_filter: Optional[str] = None,
    batch_size_override: Optional[int] = None,
) -> int:
    config = _load_yaml(config_path)
    alphagenome_cfg = _as_mapping(config, "alphagenome")
    reference_cfg = _as_mapping(config, "reference")
    variants_cfg = _as_mapping(config, "variants")
    windows_cfg = _as_mapping(config, "windows")
    outputs_cfg = _as_mapping(config, "outputs")
    runtime_cfg = _as_mapping(config, "runtime")
    output_cfg = _as_mapping(config, "output")

    chromosome = str(reference_cfg.get("chromosome", "chr15"))
    reference_fasta = (
        ref_fasta_override.expanduser().resolve()
        if ref_fasta_override
        else _resolve_reference_fasta(reference_cfg.get("fasta"), config_path=config_path)
    )
    if vcf_override:
        vcf_path = vcf_override.expanduser().resolve()
    elif variants_cfg.get("vcf"):
        vcf_path = _resolve_config_path(variants_cfg.get("vcf"), config_path=config_path, label="variants.vcf")
    else:
        vcf_path = _resolve_dataset_chromosome_vcf(str(variants_cfg.get("dataset_id", "1kg_high_coverage")), chromosome)
    output_root = output_dir_override.expanduser().resolve() if output_dir_override else Path(str(output_cfg.get("dataset_dir", "results/alphagenome_chr15_1kg"))).expanduser()
    window_bp = int(windows_cfg.get("window_bp", DEFAULT_WINDOW_BP))
    stride_bp = int(windows_cfg.get("stride_bp", DEFAULT_STRIDE_BP))
    include_indels = bool(variants_cfg.get("include_indels", False))
    require_phased = bool(variants_cfg.get("phased", True))
    resume = bool(runtime_cfg.get("resume", True))
    batch_size = int(batch_size_override or runtime_cfg.get("batch_size", 1))
    strands = _as_list(outputs_cfg.get("strands", ["plus", "minus"]), ["plus", "minus"])
    output_names = _as_list(outputs_override or outputs_cfg.get("requested_outputs"), DEFAULT_OUTPUTS)
    ontology_terms = _as_list(outputs_cfg.get("ontology_terms"), DEFAULT_ONTOLOGY_TERMS)

    if window_bp != DEFAULT_WINDOW_BP or stride_bp != DEFAULT_STRIDE_BP:
        raise ValueError("Esta configuracao inicial exige window_bp=1048576 e stride_bp=524288")
    if batch_size <= 0:
        raise ValueError("runtime.batch_size deve ser positivo")

    reference = _read_fasta_chromosome(reference_fasta, chromosome)
    chromosome_length = int(reference_cfg.get("chromosome_length") or len(reference))
    if chromosome_length > len(reference):
        raise ValueError("reference.chromosome_length excede o tamanho carregado do FASTA")
    reference = reference[:chromosome_length]
    windows = _window_plan(chromosome_length, window_bp, stride_bp)
    if max_windows is not None:
        if max_windows <= 0:
            raise ValueError("max_windows deve ser positivo")
        windows = windows[:max_windows]
    windows = _shard(windows, window_shard_index, num_window_shards)

    haplotypes = [haplotype_filter] if haplotype_filter else ["H1", "H2"]
    if any(haplotype not in {"H1", "H2"} for haplotype in haplotypes):
        raise ValueError("haplotype deve ser H1 ou H2")
    if strand_filter:
        strands = [strand_filter]

    vcf_samples = _read_vcf_samples(vcf_path)
    samples = [sample_override] if sample_override else _select_samples(vcf_samples, variants_cfg.get("samples", "all"))
    missing_samples = sorted(set(samples) - set(vcf_samples))
    if missing_samples:
        raise ValueError(f"Samples ausentes no VCF: {', '.join(missing_samples[:10])}")
    samples = _shard(samples, shard_index, num_shards)

    _write_json(
        output_root / "dataset_metadata.json",
        {
            "workflow": "genomics.workflows.alphagenome.chr15_local",
            "chromosome": chromosome,
            "chromosome_length": chromosome_length,
            "window_bp": window_bp,
            "stride_bp": stride_bp,
            "variant_policy": {"snvs": True, "include_indels": include_indels, "phased": require_phased},
            "requested_outputs": output_names,
            "ontology_terms": ontology_terms,
            "strands": strands,
            "haplotypes": haplotypes,
            "max_windows": max_windows,
            "window_shard_index": window_shard_index,
            "num_window_shards": num_window_shards,
            "samples": samples,
            "reference_fasta": str(reference_fasta),
            "vcf": str(vcf_path),
            "batching_note": (
                "runtime.batch_size > 1 usa a chamada JAX interna model._predict com tensores "
                "[B, 1048576, 4]. Aumente gradualmente ate saturar VRAM sem OOM."
            ),
        },
    )
    _write_json(output_root / "window_plan.json", [window.__dict__ for window in windows])

    dna_model = _import_alphagenome(Path(str(alphagenome_cfg.get("alphagenome_research_dir"))).expanduser() if alphagenome_cfg.get("alphagenome_research_dir") else None)
    model = _create_model(dna_model, alphagenome_cfg)
    requested_outputs = _output_types(dna_model, output_names)

    manifest_path = output_root / "manifest.jsonl"
    for sample_id in samples:
        print(f"[INFO] Carregando SNVs phased para {sample_id}...", flush=True)
        variants = _load_sample_snvs(
            vcf_path,
            chromosome,
            sample_id,
            require_phased=require_phased,
            include_indels=include_indels,
        )
        print(f"[INFO] {sample_id}: {len(variants):,} SNVs nao-referencia", flush=True)
        tasks: List[PredictionTask] = []
        for haplotype in haplotypes:
            for window in windows:
                plus_sequence = _build_haplotype_window(reference, window, variants, haplotype)
                for strand in strands:
                    if strand not in {"plus", "minus"}:
                        raise ValueError("outputs.strands deve conter apenas 'plus' e/ou 'minus'")
                    marker = _marker_path(output_root, sample_id, haplotype, strand, window.index)
                    if resume and marker.exists():
                        continue
                    sequence = plus_sequence if strand == "plus" else _reverse_complement(plus_sequence)
                    tasks.append(PredictionTask(sample_id, haplotype, strand, window, sequence))

        for task_batch in _chunks(tasks, batch_size):
            start = time.time()
            if batch_size == 1:
                task = task_batch[0]
                outputs = model.predict_sequence(
                    task.sequence,
                    requested_outputs=requested_outputs,
                    ontology_terms=ontology_terms,
                )
                elapsed = time.time() - start
                pred_dir = _prediction_path(output_root, task.sample_id, task.haplotype, task.strand, task.window.index)
                saved_outputs = _save_output_npz(pred_dir, outputs, output_names)
                task_results = [(task, saved_outputs, elapsed, elapsed)]
            else:
                prediction_rows, metadata_by_name = _predict_sequences_batched(
                    dna_model,
                    model,
                    [task.sequence for task in task_batch],
                    requested_outputs,
                    ontology_terms,
                )
                elapsed = time.time() - start
                task_results = []
                per_task_elapsed = elapsed / max(len(task_batch), 1)
                for task, prediction_row in zip(task_batch, prediction_rows):
                    pred_dir = _prediction_path(output_root, task.sample_id, task.haplotype, task.strand, task.window.index)
                    saved_outputs = _save_prediction_arrays(pred_dir, prediction_row, metadata_by_name)
                    task_results.append((task, saved_outputs, per_task_elapsed, elapsed))

            for task, saved_outputs, task_elapsed, batch_elapsed in task_results:
                pred_dir = _prediction_path(output_root, task.sample_id, task.haplotype, task.strand, task.window.index)
                _write_json(
                    pred_dir / "window_metadata.json",
                    {
                        "sample_id": task.sample_id,
                        "haplotype": task.haplotype,
                        "strand": task.strand,
                        "window": task.window.__dict__,
                        "saved_outputs": saved_outputs,
                        "elapsed_seconds": task_elapsed,
                        "batch_elapsed_seconds": batch_elapsed,
                        "batch_size": len(task_batch),
                        "minus_strand_orientation": "raw_reverse_complement_sequence" if task.strand == "minus" else None,
                    },
                )
                marker = _marker_path(output_root, task.sample_id, task.haplotype, task.strand, task.window.index)
                marker.parent.mkdir(parents=True, exist_ok=True)
                marker.write_text("ok\n", encoding="utf-8")
                _append_jsonl(
                    manifest_path,
                    {
                        "sample_id": task.sample_id,
                        "haplotype": task.haplotype,
                        "strand": task.strand,
                        "window_index": task.window.index,
                        "prediction_dir": str(pred_dir),
                        "elapsed_seconds": round(task_elapsed, 3),
                        "batch_elapsed_seconds": round(batch_elapsed, 3),
                        "batch_size": len(task_batch),
                    },
                )
                print(
                    f"[INFO] {task.sample_id} {task.haplotype} {task.strand} "
                    f"window {task.window.index + 1}/{len(windows)} em {task_elapsed:.1f}s/janela "
                    f"(batch={len(task_batch)}, batch_total={batch_elapsed:.1f}s)",
                    flush=True,
                )
    print(f"[DONE] Dataset inicial escrito em {output_root}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Predicoes AlphaGenome locais para o chr15 do 1000 Genomes")
    parser.add_argument("--config", "-c", type=Path, required=True)
    parser.add_argument("--sample", default=None, help="Executa somente um sample ID")
    parser.add_argument("--ref-fasta", type=Path, default=None, help="Override de reference.fasta")
    parser.add_argument("--vcf", type=Path, default=None, help="Override de variants.vcf")
    parser.add_argument("--output-dir", type=Path, default=None, help="Override de output.dataset_dir")
    parser.add_argument("--outputs", default=None, help="CSV de outputs AlphaGenome, ex: RNA_SEQ")
    parser.add_argument("--shard-index", type=int, default=None)
    parser.add_argument("--num-shards", type=int, default=None)
    parser.add_argument("--window-shard-index", type=int, default=None)
    parser.add_argument("--num-window-shards", type=int, default=None)
    parser.add_argument("--max-windows", type=int, default=None, help="Limita janelas para smoke tests")
    parser.add_argument("--haplotype", choices=["H1", "H2"], default=None, help="Limita a um haplotipo")
    parser.add_argument("--strand", choices=["plus", "minus"], default=None, help="Limita a uma strand")
    parser.add_argument("--batch-size", type=int, default=None, help="Batch JAX real de sequencias 1 MiB")
    args = parser.parse_args(argv)
    return run(
        args.config,
        sample_override=args.sample,
        ref_fasta_override=args.ref_fasta,
        vcf_override=args.vcf,
        output_dir_override=args.output_dir,
        outputs_override=args.outputs,
        shard_index=args.shard_index,
        num_shards=args.num_shards,
        window_shard_index=args.window_shard_index,
        num_window_shards=args.num_window_shards,
        max_windows=args.max_windows,
        haplotype_filter=args.haplotype,
        strand_filter=args.strand,
        batch_size_override=args.batch_size,
    )


if __name__ == "__main__":
    raise SystemExit(main())
