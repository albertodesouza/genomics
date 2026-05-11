from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

from genotype_based_predictor.config import load_config
from genotype_based_predictor.dynamic_indel_alignment import (
    DEFAULT_ALIGNMENT_CENTER_WINDOW_SIZE,
    DynamicIndelAligner,
    _allele_sequence,
    _parse_gt,
    _resolve_ref_index,
    _stream_gene_variants,
)


DEFAULT_BATCH_SIZE = 16


def _read_fasta_sequence(path: Path) -> str:
    lines: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines).upper()


def _read_sample_ids(path: Path) -> List[str]:
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def _load_sample_ids(config, limit: Optional[int], all_samples: bool) -> List[str]:
    sample_ids: List[str] = list(config.dataset_input.sample_ids or [])
    sample_ids_path = config.dataset_input.sample_ids_path
    if sample_ids_path:
        sample_ids.extend(_read_sample_ids(Path(sample_ids_path)))

    if not sample_ids:
        metadata_path = Path(config.dataset_input.dataset_dir) / "dataset_metadata.json"
        with open(metadata_path) as f:
            metadata = json.load(f)
        sample_ids.extend(str(sample_id) for sample_id in metadata.get("individuals", []))

    seen = set()
    unique = []
    for sample_id in sample_ids:
        if sample_id in seen:
            continue
        seen.add(sample_id)
        unique.append(sample_id)
    if all_samples or limit is None or limit <= 0:
        return unique
    return unique[:limit]


def _build_reference_alignment(ref_sequence: str, axis: Dict[str, object]) -> str:
    expanded_length = int(axis["expanded_length"])
    aligned = ["X"] * expanded_length
    expanded_index_map = axis["expanded_index_map"]
    for ref_idx_raw, expanded_idx in expanded_index_map.items():
        ref_idx = int(ref_idx_raw)
        if 0 <= ref_idx < len(ref_sequence):
            aligned[int(expanded_idx)] = ref_sequence[ref_idx]
    return "".join(aligned)


def _load_window_metadata(dataset_dir: Path, gene_name: str) -> Dict[str, object]:
    with open(dataset_dir / "references" / "windows" / gene_name / "window_metadata.json") as f:
        return json.load(f)


def _build_batch_haplotype_alignments(
    *,
    ref_aligned: str,
    ref_sequence: str,
    axis: Dict[str, object],
    vcf_path: str,
    sample_ids: List[str],
) -> Dict[str, Dict[str, str]]:
    expanded_index_map = {int(k): int(v) for k, v in axis["expanded_index_map"].items()}
    insertion_slots_by_ref = {int(k): [int(x) for x in v] for k, v in axis.get("insertion_slots_by_ref", {}).items()}
    alignment_start_1based = int(axis["alignment_start_1based"])
    alignment_end_1based = int(axis["alignment_end_1based"])
    region = str(axis.get("region") or f"{axis.get('chromosome', '')}:{alignment_start_1based}-{alignment_end_1based}")
    if not region or region.startswith(":"):
        raise ValueError("Regiao de alinhamento ausente no axis")

    result = {
        sample_id: {"H1": list(ref_aligned), "H2": list(ref_aligned)}
        for sample_id in sample_ids
    }

    query_region = region
    # The axis region can be wider when it was built to compute haplotype offsets.
    chrom = query_region.split(":", 1)[0]
    query_region = f"{chrom}:{alignment_start_1based}-{alignment_end_1based}"
    for row in _stream_gene_variants(vcf_path, sample_ids, query_region):
        ref = str(row["ref"])
        alt = str(row["alt"])
        ref_idx = _resolve_ref_index(ref_sequence, int(row["pos_1based"]) - alignment_start_1based, ref)
        if ref_idx < 0 or ref_idx >= len(ref_sequence):
            continue
        gts = row.get("gts", [])
        for sample_idx, sample_id in enumerate(sample_ids):
            if sample_idx >= len(gts):
                continue
            gt_left, gt_right = _parse_gt(str(gts[sample_idx]))
            for haplotype, token in (("H1", gt_left), ("H2", gt_right)):
                if token in {"0", ".", "nan", "None"}:
                    continue
                allele = _allele_sequence(ref, alt, token)
                aligned = result[sample_id][haplotype]
                for offset in range(len(ref)):
                    target_ref_idx = ref_idx + offset
                    target = expanded_index_map.get(target_ref_idx)
                    if target is None or not (0 <= target < len(aligned)):
                        continue
                    aligned[target] = allele[offset] if offset < len(allele) else "X"
                if len(allele) > len(ref):
                    inserted = allele[len(ref):]
                    for slot, base in zip(insertion_slots_by_ref.get(ref_idx, []), inserted):
                        if 0 <= slot < len(aligned):
                            aligned[slot] = base

    return {
        sample_id: {hap: "".join(chars) for hap, chars in haps.items()}
        for sample_id, haps in result.items()
    }


def export_aligned_dna(
    config_path: Path,
    output_path: Path,
    gene_name: Optional[str],
    sample_limit: Optional[int],
    all_samples: bool = False,
    batch_size: int = DEFAULT_BATCH_SIZE,
    center_window_size: Optional[int] = DEFAULT_ALIGNMENT_CENTER_WINDOW_SIZE,
) -> Path:
    config = load_config(config_path)
    dataset_dir = Path(config.dataset_input.dataset_dir)
    genes = list(config.dataset_input.genes_to_use or [])
    if gene_name is None:
        if not genes:
            raise ValueError("Informe --gene ou configure dataset_input.genes_to_use")
        gene_name = genes[0]

    sample_ids = _load_sample_ids(config, sample_limit, all_samples)
    if not sample_ids:
        raise ValueError("Nenhum sample_id encontrado na config/view")

    aligner = DynamicIndelAligner(
        dataset_dir,
        selected_sample_ids=sample_ids,
        sample_cache_limit=max(batch_size, 8),
        sample_prefetch_size=max(batch_size, 8),
        center_window_size=center_window_size,
    )
    aligner.build_alignment_axis_for_gene(gene_name, sample_ids)

    axis = aligner.get_alignment_axis(gene_name)
    expanded_length = int(axis["expanded_length"])
    ref_start_offset = int(axis.get("ref_start_offset", 0))
    alignment_start_1based = int(axis.get("alignment_start_1based", 0))
    alignment_end_1based = int(axis.get("alignment_end_1based", 0))
    ref_path = dataset_dir / "references" / "windows" / gene_name / "ref.window.fa"
    ref_sequence = _read_fasta_sequence(ref_path)
    if ref_start_offset or center_window_size:
        ref_sequence = ref_sequence[ref_start_offset:ref_start_offset + int(axis.get("ref_length", len(ref_sequence)))]
    ref_aligned = _build_reference_alignment(ref_sequence, axis)
    window_meta = _load_window_metadata(dataset_dir, gene_name)
    raw_variant_source = window_meta.get("raw_variant_source") or {}
    vcf_path = raw_variant_source.get("vcf_path")
    if not vcf_path:
        raise KeyError(f"raw_variant_source.vcf_path ausente para {gene_name}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        out.write(f"# gene={gene_name}\n")
        out.write(f"# expanded_length={expanded_length}\n")
        out.write(f"# alignment_ref_start_offset={ref_start_offset}\n")
        out.write(f"# alignment_start_1based={alignment_start_1based}\n")
        out.write(f"# alignment_end_1based={alignment_end_1based}\n")
        out.write(f"# center_window_size={'full' if center_window_size is None else center_window_size}\n")
        out.write("sample_id\tH1_aligned\tH2_aligned\n")
        out.write(f"REF\t{ref_aligned}\t{ref_aligned}\n")

        for batch_start in range(0, len(sample_ids), max(batch_size, 1)):
            batch = sample_ids[batch_start:batch_start + max(batch_size, 1)]
            aligner._build_sample_entries_for_gene(gene_name, batch)
            batch_alignments = _build_batch_haplotype_alignments(
                ref_aligned=ref_aligned,
                ref_sequence=ref_sequence,
                axis=axis,
                vcf_path=str(vcf_path),
                sample_ids=batch,
            )
            for sample_id in batch:
                h1_aligned = batch_alignments[sample_id]["H1"]
                h2_aligned = batch_alignments[sample_id]["H2"]
                out.write(f"{sample_id}\t{h1_aligned}\t{h2_aligned}\n")
            aligner._sample_entry_cache[gene_name].clear()

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Export aligned DNA sequences using DynamicIndelAligner")
    parser.add_argument("config_path", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--gene", type=str, default=None)
    parser.add_argument("--sample-limit", type=int, default=5)
    parser.add_argument("--all-samples", action="store_true", help="export all samples from the view or dataset metadata")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="number of individuals processed at a time")
    parser.add_argument("--center-window-size", type=int, default=DEFAULT_ALIGNMENT_CENTER_WINDOW_SIZE, help="reference-centered region to align; default 32768")
    parser.add_argument("--full-window", action="store_true", help="align the full reference window instead of the centered region")
    args = parser.parse_args()
    output_path = export_aligned_dna(
        args.config_path.resolve(),
        args.output_path.resolve(),
        args.gene,
        args.sample_limit,
        all_samples=args.all_samples,
        batch_size=args.batch_size,
        center_window_size=None if args.full_window else args.center_window_size,
    )
    print(output_path)


if __name__ == "__main__":
    main()
