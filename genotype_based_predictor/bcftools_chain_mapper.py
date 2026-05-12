from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from genotype_based_predictor.dynamic_indel_alignment import DynamicIndelAligner


DEFAULT_DATASET_DIR = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
DEFAULT_CONSENSUS_ROOT = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_all")
DEFAULT_TSV_ROOT = Path("genotype_based_predictor/aligned_dna_genes_1000_all")
BCFTOOLS_CHAIN_MAPPER_VERSION = "bcftools_chain_mapper_v3"


@dataclass(frozen=True)
class ChainBlock:
    size: int
    dt: int
    dq: int


@dataclass(frozen=True)
class ChainHeader:
    t_start: int
    t_end: int
    q_start: int
    q_end: int


def read_fasta_sequence(path: Path) -> str:
    lines: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines).upper()


def adjust_to_target_size(consensus_seq: str, ref_window_seq: str, target_size: int) -> str:
    length = len(consensus_seq)
    if length == target_size:
        return consensus_seq
    if length > target_size:
        return consensus_seq[:target_size]
    needed = target_size - length
    pad_from_ref = ref_window_seq[length:length + needed]
    if len(pad_from_ref) < needed:
        pad_from_ref += "N" * (needed - len(pad_from_ref))
    return consensus_seq + pad_from_ref


def run_bcftools_consensus_with_chain(
    *,
    case_dir: Path,
    sample_id: str,
    haplotype: str,
    cache_dir: Path,
    force: bool = False,
) -> Tuple[Path, Path]:
    if shutil.which("bcftools") is None:
        raise RuntimeError("bcftools nao encontrado no PATH. Rode source scripts/start_genomics_universal.sh")
    hap_arg = "1" if haplotype == "H1" else "2"
    ref_path = case_dir / "ref.window.fa"
    vcf_path = case_dir / f"{sample_id}.window.consensus_ready.vcf.gz"
    if not ref_path.exists():
        raise FileNotFoundError(f"ref.window.fa nao encontrado: {ref_path}")
    if not vcf_path.exists():
        raise FileNotFoundError(f"consensus_ready VCF nao encontrado: {vcf_path}")
    cache_dir.mkdir(parents=True, exist_ok=True)
    raw_out = cache_dir / f"{sample_id}.{haplotype}.window.raw.rebuilt.fa"
    chain_out = cache_dir / f"{sample_id}.{haplotype}.window.raw.chain"
    if raw_out.exists() and chain_out.exists() and not force:
        return raw_out, chain_out
    if force:
        raw_out.unlink(missing_ok=True)
        chain_out.unlink(missing_ok=True)
    with open(raw_out, "w") as out:
        proc = subprocess.run(
            ["bcftools", "consensus", "-H", hap_arg, "-f", str(ref_path), "-c", str(chain_out), str(vcf_path)],
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
        )
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, proc.args, stderr=proc.stderr)
    return raw_out, chain_out


def parse_chain(path: Path) -> Tuple[ChainHeader, List[ChainBlock]]:
    with open(path) as f:
        header_line = ""
        for line in f:
            line = line.strip()
            if line:
                header_line = line
                break
        if not header_line.startswith("chain "):
            raise ValueError(f"Arquivo chain invalido: {path}")
        parts = header_line.split()
        header = ChainHeader(
            t_start=int(parts[5]),
            t_end=int(parts[6]),
            q_start=int(parts[10]),
            q_end=int(parts[11]),
        )
        blocks: List[ChainBlock] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = [int(x) for x in line.split()]
            if len(fields) == 1:
                blocks.append(ChainBlock(size=fields[0], dt=0, dq=0))
            elif len(fields) == 3:
                blocks.append(ChainBlock(size=fields[0], dt=fields[1], dq=fields[2]))
            else:
                raise ValueError(f"Linha chain invalida em {path}: {line}")
    return header, blocks


def chain_to_fixed_fasta_map(
    chain_path: Path,
    raw_length: int,
    fixed_length: int,
) -> Tuple[Dict[int, int], Dict[int, Tuple[int, int]], List[int]]:
    """Return maps from fixed FASTA index to full-window reference/insertion/deletion coordinates.

    The chain maps target=reference to query=raw consensus. Coordinates are converted
    to 0-based offsets relative to the original reference window and raw consensus.
    """
    header, blocks = parse_chain(chain_path)
    consensus_to_ref: Dict[int, int] = {}
    consensus_insertions: Dict[int, Tuple[int, int]] = {}
    deletion_ref_indices: List[int] = []

    t_pos = header.t_start
    q_pos = header.q_start
    for block in blocks:
        for offset in range(block.size):
            consensus_idx = q_pos + offset - header.q_start
            ref_idx = t_pos + offset - header.t_start
            consensus_to_ref[consensus_idx] = ref_idx
        t_pos += block.size
        q_pos += block.size

        if block.dt:
            for offset in range(block.dt):
                deletion_ref_indices.append(t_pos + offset - header.t_start)
        if block.dq:
            insertion_after_ref_idx = t_pos - header.t_start - 1
            for order in range(block.dq):
                consensus_insertions[q_pos + order - header.q_start] = (insertion_after_ref_idx, order)
        t_pos += block.dt
        q_pos += block.dq

    # adjust_to_target_size pads short consensus with reference tail indexed by consensus length.
    if fixed_length > raw_length:
        for fasta_idx in range(raw_length, fixed_length):
            consensus_to_ref[fasta_idx] = fasta_idx
    elif fixed_length < raw_length:
        consensus_to_ref = {k: v for k, v in consensus_to_ref.items() if k < fixed_length}
        consensus_insertions = {k: v for k, v in consensus_insertions.items() if k < fixed_length}

    return consensus_to_ref, consensus_insertions, deletion_ref_indices


def build_chain_entry(
    *,
    axis: Dict[str, object],
    chain_path: Path,
    raw_length: int,
    fixed_length: int,
) -> Dict[str, object]:
    consensus_to_ref, consensus_insertions, deletion_ref_indices = chain_to_fixed_fasta_map(chain_path, raw_length, fixed_length)
    ref_start_offset = int(axis.get("ref_start_offset", 0))
    ref_length = int(axis.get("ref_length", 0))
    expanded_index_map = {int(k): int(v) for k, v in axis.get("expanded_index_map", {}).items()}
    insertion_slots_by_ref = {int(k): [int(x) for x in v] for k, v in axis.get("insertion_slots_by_ref", {}).items()}

    copy_from_indices: List[int] = []
    expanded_indices: List[int] = []
    insertion_indices: List[int] = []
    deletion_indices: List[int] = []
    skipped_insertions = 0
    skipped_ref = 0

    for fasta_idx, ref_idx_full in sorted(consensus_to_ref.items()):
        local_ref_idx = int(ref_idx_full) - ref_start_offset
        expanded_idx = expanded_index_map.get(local_ref_idx)
        if expanded_idx is None or not (0 <= local_ref_idx < ref_length):
            skipped_ref += 1
            continue
        copy_from_indices.append(int(fasta_idx))
        expanded_indices.append(expanded_idx)

    for fasta_idx, (after_ref_full, order) in sorted(consensus_insertions.items()):
        local_ref_idx = int(after_ref_full) - ref_start_offset
        slots = insertion_slots_by_ref.get(local_ref_idx, [])
        if 0 <= order < len(slots):
            copy_from_indices.append(int(fasta_idx))
            expanded_indices.append(int(slots[order]))
            insertion_indices.append(int(slots[order]))
        else:
            skipped_insertions += 1

    for ref_idx_full in deletion_ref_indices:
        local_ref_idx = int(ref_idx_full) - ref_start_offset
        expanded_idx = expanded_index_map.get(local_ref_idx)
        if expanded_idx is not None and 0 <= local_ref_idx < ref_length:
            deletion_indices.append(int(expanded_idx))

    ordered = sorted(zip(copy_from_indices, expanded_indices), key=lambda item: item[1])
    return {
        "source_start_idx": 0,
        "copy_from_indices": [int(a) for a, _b in ordered],
        "expanded_indices": [int(b) for _a, b in ordered],
        "insertion_indices": sorted(set(insertion_indices)),
        "deletion_indices": sorted(set(deletion_indices)),
        "mapping_method": "bcftools_chain",
        "skipped_ref_positions": skipped_ref,
        "skipped_insertions": skipped_insertions,
    }


class BcftoolsChainMapper:
    def __init__(
        self,
        *,
        dataset_dir: Path,
        consensus_dataset_dir: Path,
        aligner: DynamicIndelAligner,
        cache_dir: Optional[Path] = None,
    ):
        self.dataset_dir = Path(dataset_dir)
        self.consensus_dataset_dir = Path(consensus_dataset_dir)
        self.aligner = aligner
        self.cache_dir = Path(cache_dir) if cache_dir else self.dataset_dir / "alignment_cache" / BCFTOOLS_CHAIN_MAPPER_VERSION
        self._entry_cache: Dict[Tuple[str, str, str, str], Dict[str, object]] = {}

    def _case_dir(self, gene: str, sample_id: str) -> Path:
        return self.consensus_dataset_dir / "individuals" / sample_id / "windows" / gene

    def _axis_cache_key(self, gene: str) -> str:
        axis = self.aligner.get_alignment_axis(gene)
        return (
            f"expanded_{int(axis.get('expanded_length', 0))}"
            f"_ref_{int(axis.get('ref_length', 0))}"
            f"_offset_{int(axis.get('ref_start_offset', 0))}"
            f"_{axis.get('sample_set_key', 'samples') }"
        )

    def _entry_cache_path(self, gene: str, sample_id: str, haplotype: str, axis_key: str) -> Path:
        return self.cache_dir / gene / axis_key / sample_id / f"{haplotype}.entry.json"

    def get_haplotype_entry(self, gene: str, sample_id: str, haplotype: str) -> Optional[Dict[str, object]]:
        axis_key = self._axis_cache_key(gene)
        key = (gene, sample_id, haplotype, axis_key)
        if key in self._entry_cache:
            return self._entry_cache[key]
        path = self._entry_cache_path(gene, sample_id, haplotype, axis_key)
        if path.exists():
            try:
                with open(path) as f:
                    payload = json.load(f)
                if payload.get("version") == BCFTOOLS_CHAIN_MAPPER_VERSION:
                    entry = payload.get("entry")
                    if isinstance(entry, dict):
                        self._entry_cache[key] = entry
                        return entry
            except Exception:
                pass

        case_dir = self._case_dir(gene, sample_id)
        fixed_path = case_dir / f"{sample_id}.{haplotype}.window.fixed.fa"
        if not fixed_path.exists():
            return None
        raw_path = case_dir / f"{sample_id}.{haplotype}.window.raw.fa"
        raw_rebuilt, chain_path = self._build_validated_consensus(case_dir, sample_id, haplotype)
        raw_rebuilt_seq = read_fasta_sequence(raw_rebuilt)
        fixed_seq = read_fasta_sequence(fixed_path)
        ref_seq = read_fasta_sequence(case_dir / "ref.window.fa")

        axis = self.aligner.get_alignment_axis(gene)
        entry = build_chain_entry(
            axis=axis,
            chain_path=chain_path,
            raw_length=len(raw_rebuilt_seq),
            fixed_length=len(fixed_seq),
        )
        entry["fasta_path"] = str(fixed_path)
        entry["chain_path"] = str(chain_path)
        entry["raw_rebuilt_path"] = str(raw_rebuilt)
        entry["version"] = BCFTOOLS_CHAIN_MAPPER_VERSION
        entry["axis_cache_key"] = axis_key
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump({"version": BCFTOOLS_CHAIN_MAPPER_VERSION, "entry": entry}, f, indent=2)
        self._entry_cache[key] = entry
        return entry

    def _build_validated_consensus(self, case_dir: Path, sample_id: str, haplotype: str) -> Tuple[Path, Path]:
        cache_dir = self.cache_dir / case_dir.name / sample_id / haplotype
        raw_path = case_dir / f"{sample_id}.{haplotype}.window.raw.fa"
        fixed_path = case_dir / f"{sample_id}.{haplotype}.window.fixed.fa"
        ref_path = case_dir / "ref.window.fa"
        raw_rebuilt, chain_path = run_bcftools_consensus_with_chain(
            case_dir=case_dir,
            sample_id=sample_id,
            haplotype=haplotype,
            cache_dir=cache_dir,
        )

        def is_valid(raw_file: Path) -> Tuple[bool, str]:
            raw_rebuilt_seq = read_fasta_sequence(raw_file)
            if raw_path.exists() and raw_rebuilt_seq != read_fasta_sequence(raw_path):
                return False, f"Consenso reconstruido difere do raw existente: {raw_path}"
            fixed_seq = read_fasta_sequence(fixed_path)
            ref_seq = read_fasta_sequence(ref_path)
            rebuilt_fixed = adjust_to_target_size(raw_rebuilt_seq, ref_seq, len(fixed_seq))
            if rebuilt_fixed != fixed_seq:
                return False, f"Consenso reconstruido difere do fixed existente: {fixed_path}"
            return True, ""

        valid, reason = is_valid(raw_rebuilt)
        if valid:
            return raw_rebuilt, chain_path

        raw_rebuilt, chain_path = run_bcftools_consensus_with_chain(
            case_dir=case_dir,
            sample_id=sample_id,
            haplotype=haplotype,
            cache_dir=cache_dir,
            force=True,
        )
        valid, reason = is_valid(raw_rebuilt)
        if not valid:
            raise ValueError(reason)
        return raw_rebuilt, chain_path


def load_tsv_sample_ids(tsv_root: Path, gene: str) -> List[str]:
    idx_path = tsv_root / f"{gene}.tsv.idx.json"
    with open(idx_path) as f:
        payload = json.load(f)
    return [str(x) for x in payload.get("samples", [])]


def main() -> None:
    parser = argparse.ArgumentParser(description="Prototype bcftools chain -> expanded-axis mapper")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--haplotype", choices=["H1", "H2"], required=True)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument("--consensus-root", type=Path, default=DEFAULT_CONSENSUS_ROOT)
    parser.add_argument("--tsv-root", type=Path, default=DEFAULT_TSV_ROOT)
    parser.add_argument("--cache-dir", type=Path, default=None)
    args = parser.parse_args()

    case_dir = args.consensus_root / "individuals" / args.sample / "windows" / args.gene
    if args.cache_dir is None:
        cache_dir = Path(tempfile.gettempdir()) / "bcftools_chain_mapper_probe" / args.gene / args.sample / args.haplotype
    else:
        cache_dir = args.cache_dir / args.gene / args.sample / args.haplotype
    raw_rebuilt, chain_path = run_bcftools_consensus_with_chain(
        case_dir=case_dir,
        sample_id=args.sample,
        haplotype=args.haplotype,
        cache_dir=cache_dir,
    )

    raw_rebuilt_seq = read_fasta_sequence(raw_rebuilt)
    raw_existing_path = case_dir / f"{args.sample}.{args.haplotype}.window.raw.fa"
    fixed_existing_path = case_dir / f"{args.sample}.{args.haplotype}.window.fixed.fa"
    ref_seq = read_fasta_sequence(case_dir / "ref.window.fa")
    fixed_existing_seq = read_fasta_sequence(fixed_existing_path)
    fixed_rebuilt_seq = adjust_to_target_size(raw_rebuilt_seq, ref_seq, len(fixed_existing_seq))
    raw_match = raw_existing_path.exists() and raw_rebuilt_seq == read_fasta_sequence(raw_existing_path)
    fixed_match = fixed_rebuilt_seq == fixed_existing_seq

    sample_ids = load_tsv_sample_ids(args.tsv_root, args.gene)
    aligner = DynamicIndelAligner(args.dataset_dir, selected_sample_ids=sample_ids, center_window_size=32768)
    aligner.build_alignment_axis_for_gene(args.gene, sample_ids)
    axis = aligner.get_alignment_axis(args.gene)
    entry = build_chain_entry(axis=axis, chain_path=chain_path, raw_length=len(raw_rebuilt_seq), fixed_length=len(fixed_existing_seq))

    print(json.dumps({
        "gene": args.gene,
        "sample": args.sample,
        "haplotype": args.haplotype,
        "case_dir": str(case_dir),
        "raw_rebuilt": str(raw_rebuilt),
        "chain": str(chain_path),
        "raw_length": len(raw_rebuilt_seq),
        "fixed_length": len(fixed_existing_seq),
        "raw_match": raw_match,
        "fixed_match": fixed_match,
        "entry_summary": {
            "copy_positions": len(entry["copy_from_indices"]),
            "insertions": len(entry["insertion_indices"]),
            "deletions": len(entry["deletion_indices"]),
            "skipped_ref_positions": entry["skipped_ref_positions"],
            "skipped_insertions": entry["skipped_insertions"],
            "mapping_method": entry["mapping_method"],
        },
    }, indent=2))


if __name__ == "__main__":
    main()
