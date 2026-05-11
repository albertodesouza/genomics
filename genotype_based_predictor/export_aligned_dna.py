from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

from genotype_based_predictor.config import load_config
from genotype_based_predictor.dynamic_indel_alignment import DynamicIndelAligner


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


def _build_haplotype_alignment(sequence: str, entry: Dict[str, object], expanded_length: int) -> str:
    aligned = ["X"] * expanded_length
    copy_from = entry.get("copy_from_indices", [])
    copy_to = entry.get("expanded_indices", [])
    for source_idx_raw, target_idx_raw in zip(copy_from, copy_to):
        source_idx = int(source_idx_raw)
        target_idx = int(target_idx_raw)
        if 0 <= source_idx < len(sequence) and 0 <= target_idx < expanded_length:
            aligned[target_idx] = sequence[source_idx]
    return "".join(aligned)


def export_aligned_dna(
    config_path: Path,
    output_path: Path,
    gene_name: Optional[str],
    sample_limit: Optional[int],
    all_samples: bool = False,
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

    aligner = DynamicIndelAligner(dataset_dir, selected_sample_ids=sample_ids)
    for sample_id in sample_ids:
        aligner.get_haplotype_entry(gene_name, sample_id, "H1")
        aligner.get_haplotype_entry(gene_name, sample_id, "H2")

    axis = aligner.get_alignment_axis(gene_name)
    expanded_length = int(axis["expanded_length"])
    ref_path = dataset_dir / "references" / "windows" / gene_name / "ref.window.fa"
    ref_sequence = _read_fasta_sequence(ref_path)
    ref_aligned = _build_reference_alignment(ref_sequence, axis)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        out.write(f"# gene={gene_name}\n")
        out.write(f"# expanded_length={expanded_length}\n")
        out.write("sample_id\tH1_aligned\tH2_aligned\n")
        out.write(f"REF\t{ref_aligned}\t{ref_aligned}\n")

        for sample_id in sample_ids:
            sample_dir = dataset_dir / "individuals" / sample_id / "windows" / gene_name
            h1_sequence = _read_fasta_sequence(sample_dir / f"{sample_id}.H1.window.fixed.fa")
            h2_sequence = _read_fasta_sequence(sample_dir / f"{sample_id}.H2.window.fixed.fa")
            h1_entry = aligner.get_haplotype_entry(gene_name, sample_id, "H1")
            h2_entry = aligner.get_haplotype_entry(gene_name, sample_id, "H2")
            if h1_entry is None or h2_entry is None:
                raise RuntimeError(f"Entrada de alinhamento ausente para {sample_id}")
            h1_aligned = _build_haplotype_alignment(h1_sequence, h1_entry, expanded_length)
            h2_aligned = _build_haplotype_alignment(h2_sequence, h2_entry, expanded_length)
            out.write(f"{sample_id}\t{h1_aligned}\t{h2_aligned}\n")

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Export aligned DNA sequences using DynamicIndelAligner")
    parser.add_argument("config_path", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--gene", type=str, default=None)
    parser.add_argument("--sample-limit", type=int, default=5)
    parser.add_argument("--all-samples", action="store_true", help="export all samples from the view or dataset metadata")
    args = parser.parse_args()
    output_path = export_aligned_dna(
        args.config_path.resolve(),
        args.output_path.resolve(),
        args.gene,
        args.sample_limit,
        all_samples=args.all_samples,
    )
    print(output_path)


if __name__ == "__main__":
    main()
