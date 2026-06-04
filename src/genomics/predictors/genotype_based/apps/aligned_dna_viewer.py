from __future__ import annotations

import argparse
import gzip
import html
import json
import shutil
import subprocess
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse

from genomics.predictors.genotype_based.alignment.bcftools_chain_mapper import BcftoolsChainMapper
from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner
from genomics.workspace import DEFAULT_CONSENSUS_DATASET_DIR, DEFAULT_DATASET_DIR

DEFAULT_MAX_CELLS = 200_000
DEFAULT_MAX_VARIANTS = 1_000
DEFAULT_MODEL_WINDOW_SIZE = 32_768
DEFAULT_VCF_TIMEOUT_SECONDS = 20


def _parse_region(region: str) -> Tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start_raw, end_raw = span.split("-", 1)
    return chrom, int(start_raw), int(end_raw)


def _variant_type(ref: str, alt: str) -> str:
    alts = [item for item in str(alt).split(",") if item and item != "."]
    if not alts:
        return "missing"
    if any(item.startswith("<") or "[" in item or "]" in item for item in alts):
        return "SV"
    if len(ref) == 1 and all(len(item) == 1 for item in alts):
        return "SNV"
    if all(len(item) == len(ref) for item in alts):
        return "MNV"
    return "INDEL"


def _stream_vcf_variants(vcf_path: str, sample_ids: List[str], region: str):
    if shutil.which("bcftools") is None:
        yield from _stream_vcf_variants_from_gzip(vcf_path, sample_ids, region)
        return

    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n"
    cmd = ["bcftools", "query"]
    if sample_ids:
        cmd.extend(["-s", ",".join(sample_ids)])
    cmd.extend(["-r", region, "-f", fmt, vcf_path])

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            chrom, pos, variant_id, ref, alt, qual, filt = fields[:7]
            yield {
                "chrom": chrom,
                "pos_1based": int(pos),
                "id": "" if variant_id == "." else variant_id,
                "ref": ref,
                "alt": alt,
                "qual": "" if qual == "." else qual,
                "filter": filt,
                "gts": fields[7:],
            }
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
        try:
            returncode = proc.wait(timeout=DEFAULT_VCF_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=2)
            raise TimeoutError(f"Consulta VCF excedeu {DEFAULT_VCF_TIMEOUT_SECONDS}s para {region}")
        if returncode != 0:
            stderr = proc.stderr.read() if proc.stderr is not None else ""
            raise subprocess.CalledProcessError(returncode, cmd, stderr=stderr)
        if proc.stderr is not None:
            proc.stderr.close()


def _stream_vcf_variants_from_gzip(vcf_path: str, sample_ids: List[str], region: str):
    chrom, start_1based, end_1based = _parse_region(region)
    sample_column_indices: List[int] = []

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_to_col = {sample_id: idx for idx, sample_id in enumerate(header[9:])}
                missing = [sample_id for sample_id in sample_ids if sample_id not in sample_to_col]
                if missing:
                    raise KeyError(f"Samples ausentes no VCF: {missing}")
                sample_column_indices = [9 + sample_to_col[sample_id] for sample_id in sample_ids]
                break

        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            row_chrom = fields[0]
            if row_chrom != chrom and row_chrom.removeprefix("chr") != chrom.removeprefix("chr"):
                continue
            pos = int(fields[1])
            if pos < start_1based:
                continue
            if pos > end_1based:
                break

            fmt_keys = fields[8].split(":") if len(fields) > 8 else []
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                gt_idx = 0

            gts = []
            for col_idx in sample_column_indices:
                sample_field = fields[col_idx] if col_idx < len(fields) else "./."
                parts = sample_field.split(":")
                gts.append(parts[gt_idx] if gt_idx < len(parts) else parts[0])

            yield {
                "chrom": fields[0],
                "pos_1based": pos,
                "id": "" if fields[2] == "." else fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "qual": "" if fields[5] == "." else fields[5],
                "filter": fields[6],
                "gts": gts,
            }


def _decode(value: bytes) -> str:
    return value.decode("utf-8", errors="replace")


def _read_fasta_sequence(path: Path) -> str:
    lines: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines).upper()


def _gene_from_path(path: Path) -> str:
    name = path.name
    for suffix in (".color.columns.txt", ".columns.txt", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


class TsvAlignmentIndex:
    def __init__(self, tsv_path: Path):
        self.tsv_path = Path(tsv_path)
        self.index_path = self.tsv_path.with_suffix(self.tsv_path.suffix + ".idx.json")
        self.data = self._load_or_build()

    @property
    def gene(self) -> str:
        return str(self.data["gene"])

    @property
    def expanded_length(self) -> int:
        return int(self.data["expanded_length"])

    @property
    def samples(self) -> List[str]:
        return list(self.data["samples"])

    def _load_or_build(self) -> Dict:
        stat = self.tsv_path.stat()
        if self.index_path.exists():
            try:
                with open(self.index_path) as f:
                    cached = json.load(f)
                if cached.get("mtime_ns") == stat.st_mtime_ns and cached.get("size") == stat.st_size:
                    return cached
            except Exception:
                pass

        data = self._build_index(stat.st_mtime_ns, stat.st_size)
        with open(self.index_path, "w") as f:
            json.dump(data, f, indent=2)
        return data

    def _build_index(self, mtime_ns: int, size: int) -> Dict:
        comments: List[str] = []
        offsets: Dict[str, int] = {}
        samples: List[str] = []
        gene: Optional[str] = None
        expanded_length: Optional[int] = None
        alignment_start_1based: Optional[int] = None

        with open(self.tsv_path, "rb") as f:
            while True:
                offset = f.tell()
                line = f.readline()
                if not line:
                    break
                line = line.rstrip(b"\n")
                if not line:
                    continue
                if line.startswith(b"#"):
                    text = _decode(line)
                    comments.append(text)
                    if text.startswith("# gene="):
                        gene = text.split("=", 1)[1].strip()
                    elif text.startswith("# expanded_length="):
                        expanded_length = int(text.split("=", 1)[1].strip())
                    elif text.startswith("# alignment_start_1based="):
                        alignment_start_1based = int(text.split("=", 1)[1].strip())
                    continue

                first_tab = line.find(b"\t")
                if first_tab < 0:
                    continue
                sample_id = _decode(line[:first_tab])
                if sample_id == "sample_id":
                    continue
                offsets[sample_id] = offset
                if sample_id != "REF":
                    samples.append(sample_id)

        if "REF" not in offsets:
            raise ValueError(f"Arquivo TSV sem linha REF: {self.tsv_path}")

        if gene is None:
            gene = _gene_from_path(self.tsv_path)
        if expanded_length is None:
            ref_h1, _ref_h2 = self.read_sample_sequences("REF")
            expanded_length = len(ref_h1)

        return {
            "path": str(self.tsv_path.resolve()),
            "mtime_ns": mtime_ns,
            "size": size,
            "gene": gene,
            "expanded_length": expanded_length,
            "comments": comments,
            "alignment_start_1based": alignment_start_1based,
            "offsets": offsets,
            "samples": samples,
        }

    def read_sample_sequences(self, sample_id: str) -> Tuple[str, str]:
        offset = self.data["offsets"].get(sample_id)
        if offset is None:
            raise KeyError(f"Sample ausente no TSV de {self.gene}: {sample_id}")
        with open(self.tsv_path, "rb") as f:
            f.seek(int(offset))
            line = f.readline().rstrip(b"\n")
        parts = line.split(b"\t", 2)
        if len(parts) != 3:
            raise ValueError(f"Linha invalida para {sample_id} em {self.tsv_path}")
        return _decode(parts[1]), _decode(parts[2])

    def read_sample_window(self, sample_id: str, start: int, end: int) -> Tuple[str, str]:
        h1, h2 = self.read_sample_sequences(sample_id)
        return h1[start - 1:end], h2[start - 1:end]

    def expanded_window_ref_offsets(self, start: int, end: int) -> List[int]:
        ref_h1, _ref_h2 = self.read_sample_sequences("REF")
        offsets = []
        ref_offset = 0
        for expanded_pos, base in enumerate(ref_h1[:end], start=1):
            current_offset = ref_offset if base != "X" else max(ref_offset - 1, 0)
            if expanded_pos >= start:
                offsets.append(current_offset)
            if base != "X":
                ref_offset += 1
        return offsets

    def expanded_window_genomic_positions(self, start: int, end: int, reference_start_1based: int) -> List[Optional[int]]:
        ref_h1, _ref_h2 = self.read_sample_sequences("REF")
        positions: List[Optional[int]] = []
        ref_offset = 0
        for expanded_pos, base in enumerate(ref_h1[:end], start=1):
            current_pos = reference_start_1based + ref_offset if base != "X" else None
            if expanded_pos >= start:
                positions.append(current_pos)
            if base != "X":
                ref_offset += 1
        return positions

    def reference_centered_model_window(self, window_size: int = DEFAULT_MODEL_WINDOW_SIZE) -> Dict[str, int]:
        ref_h1, _ref_h2 = self.read_sample_sequences("REF")
        ref_length = sum(1 for base in ref_h1 if base != "X")
        center_ref_idx = ref_length // 2
        ref_offset = 0
        center_expanded_pos = 1
        for expanded_pos, base in enumerate(ref_h1, start=1):
            if base != "X":
                if ref_offset == center_ref_idx:
                    center_expanded_pos = expanded_pos
                    break
                ref_offset += 1
        size = max(int(window_size), 1)
        half = size // 2
        start = center_expanded_pos - half
        end = start + size - 1
        if start < 1:
            end += 1 - start
            start = 1
        if end > self.expanded_length:
            start = max(1, start - (end - self.expanded_length))
            end = self.expanded_length
        return {
            "policy": "reference_center_to_expanded_axis",
            "window_size": end - start + 1,
            "requested_window_size": size,
            "center_ref_idx_0based": center_ref_idx,
            "center_expanded_pos": center_expanded_pos,
            "start": start,
            "end": end,
        }


class AlignmentRepository:
    def __init__(
        self,
        tsv_root: Path,
        dataset_dir: Optional[Path] = None,
        alignment_mapping: str = "bcftools_chain",
        consensus_dataset_dir: Optional[Path] = None,
    ):
        self.tsv_root = Path(tsv_root)
        self.dataset_dir = Path(dataset_dir) if dataset_dir else DEFAULT_DATASET_DIR
        self.alignment_mapping = alignment_mapping
        self.consensus_dataset_dir = Path(consensus_dataset_dir) if consensus_dataset_dir else DEFAULT_CONSENSUS_DATASET_DIR
        self.sample_metadata = self._load_sample_metadata()
        self.indexes = self._discover_indexes()

    def _load_sample_metadata(self) -> Dict[str, Dict[str, str]]:
        metadata_path = self.dataset_dir / "dataset_metadata.json"
        if not metadata_path.exists():
            return {}
        try:
            with open(metadata_path) as f:
                payload = json.load(f)
        except Exception:
            return {}
        pedigree = payload.get("individuals_pedigree", {}) or {}
        if not isinstance(pedigree, dict):
            return {}
        result: Dict[str, Dict[str, str]] = {}
        for sample_id, meta in pedigree.items():
            if not isinstance(meta, dict):
                meta = {}
            result[str(sample_id)] = {
                "sample_id": str(sample_id),
                "population": str(meta.get("population", "")),
                "superpopulation": str(meta.get("superpopulation", "")),
                "sex": str(meta.get("sex_label", meta.get("sex", ""))),
            }
        return result

    def _discover_tsvs(self) -> List[Path]:
        if self.tsv_root.is_file():
            return [self.tsv_root]
        return sorted(self.tsv_root.glob("*.tsv"))

    def _discover_indexes(self) -> Dict[str, TsvAlignmentIndex]:
        indexes: Dict[str, TsvAlignmentIndex] = {}
        for path in self._discover_tsvs():
            idx = TsvAlignmentIndex(path)
            key = idx.gene
            if key in indexes:
                key = path.stem
            indexes[key] = idx
        if not indexes:
            raise FileNotFoundError(f"Nenhum .tsv encontrado em {self.tsv_root}")
        return indexes

    def genes_payload(self) -> List[Dict[str, object]]:
        return [
            {
                "id": gene,
                "gene": idx.gene,
                "expanded_length": idx.expanded_length,
                "samples": len(idx.samples),
                "path": str(idx.tsv_path),
            }
            for gene, idx in sorted(self.indexes.items())
        ]

    def get(self, gene: str) -> TsvAlignmentIndex:
        if gene not in self.indexes:
            raise KeyError(f"Gene nao encontrado: {gene}")
        return self.indexes[gene]

    def samples_payload(self, gene: str) -> Dict[str, object]:
        idx = self.get(gene)
        rows = []
        pop_counts: Dict[str, int] = {}
        super_counts: Dict[str, int] = {}
        for sample_id in idx.samples:
            meta = dict(self.sample_metadata.get(sample_id, {"sample_id": sample_id}))
            meta.setdefault("sample_id", sample_id)
            rows.append(meta)
            pop = meta.get("population") or ""
            sup = meta.get("superpopulation") or ""
            if pop:
                pop_counts[pop] = pop_counts.get(pop, 0) + 1
            if sup:
                super_counts[sup] = super_counts.get(sup, 0) + 1
        return {
            "gene": gene,
            "samples": idx.samples,
            "sample_rows": rows,
            "populations": dict(sorted(pop_counts.items())),
            "superpopulations": dict(sorted(super_counts.items())),
            "expanded_length": idx.expanded_length,
        }

    def _load_window_metadata(self, gene: str) -> Dict[str, object]:
        metadata_path = self.dataset_dir / "references" / "windows" / gene / "window_metadata.json"
        if not metadata_path.exists():
            raise FileNotFoundError(f"window_metadata.json nao encontrado para {gene}: {metadata_path}")
        with open(metadata_path) as f:
            return json.load(f)

    def variants_payload(self, gene: str, samples: List[str], start: int, length: int, max_variants: int) -> Dict[str, object]:
        idx = self.get(gene)
        start = max(start, 1)
        end = min(start + max(length, 1) - 1, idx.expanded_length)
        offsets = idx.expanded_window_ref_offsets(start, end)
        if not offsets:
            return {"gene": gene, "start": start, "end": end, "variants": [], "columns": samples, "truncated": False}

        ref_start_offset = min(offsets)
        ref_end_offset = max(offsets)
        window_meta = self._load_window_metadata(gene)
        raw_variant_source = window_meta.get("raw_variant_source") or {}
        vcf_path = raw_variant_source.get("vcf_path")
        if not vcf_path:
            raise KeyError(f"raw_variant_source.vcf_path ausente para {gene}")

        chrom = str(window_meta["chromosome"])
        reference_start_1based = idx.data.get("alignment_start_1based") or int(window_meta["start"])
        genomic_start = reference_start_1based + ref_start_offset
        genomic_end = reference_start_1based + ref_end_offset
        region = f"{chrom}:{genomic_start}-{genomic_end}"
        offset_to_expanded_positions: Dict[int, List[int]] = {}
        for expanded_pos, ref_offset in zip(range(start, end + 1), offsets):
            offset_to_expanded_positions.setdefault(ref_offset, []).append(expanded_pos)

        selected = [sample for sample in samples if sample in idx.samples]
        variants = []
        truncated = False
        for row in _stream_vcf_variants(str(vcf_path), selected, region):
            ref_offset = int(row["pos_1based"]) - reference_start_1based
            expanded_positions = offset_to_expanded_positions.get(ref_offset, [])
            if not expanded_positions:
                continue
            gts = [str(gt) for gt in row.get("gts", [])]
            called_samples = []
            alt_samples = []
            for sample_id, gt in zip(selected, gts):
                if gt in {".", "./.", ".|."}:
                    continue
                called_samples.append(sample_id)
                alleles = gt.replace("/", "|").split("|")
                if any(allele not in {"0", "."} for allele in alleles):
                    alt_samples.append(sample_id)
            variants.append({
                "expanded_pos": expanded_positions[0],
                "expanded_span": f"{expanded_positions[0]}" if len(expanded_positions) == 1 else f"{expanded_positions[0]}-{expanded_positions[-1]}",
                "genomic_pos": int(row["pos_1based"]),
                "id": row.get("id", ""),
                "ref": row["ref"],
                "alt": row["alt"],
                "type": _variant_type(str(row["ref"]), str(row["alt"])),
                "filter": row.get("filter", ""),
                "qual": row.get("qual", ""),
                "alt_count": len(alt_samples),
                "called_count": len(called_samples),
                "alt_samples": alt_samples[:12],
                "gts": {sample_id: gt for sample_id, gt in zip(selected, gts)},
            })
            if len(variants) >= max_variants:
                truncated = True
                break

        return {
            "gene": gene,
            "start": start,
            "end": end,
            "reference_region": region,
            "vcf_path": str(vcf_path),
            "columns": selected,
            "variants": variants,
            "truncated": truncated,
            "max_variants": max_variants,
        }

    def reference_genomic_positions(self, gene: str, start: int, end: int) -> List[Optional[int]]:
        idx = self.get(gene)
        reference_start_1based = idx.data.get("alignment_start_1based")
        if not reference_start_1based:
            window_meta = self._load_window_metadata(gene)
            reference_start_1based = int(window_meta["start"])
        return idx.expanded_window_genomic_positions(start, end, reference_start_1based)

    def model_window_payload(self, gene: str, window_size: int = DEFAULT_MODEL_WINDOW_SIZE) -> Dict[str, int]:
        return self.get(gene).reference_centered_model_window(window_size)

    def _build_chain_context(self, gene: str, sample_ids: List[str]) -> Tuple[DynamicIndelAligner, BcftoolsChainMapper]:
        idx = self.get(gene)
        aligner = DynamicIndelAligner(self.dataset_dir, selected_sample_ids=idx.samples, center_window_size=DEFAULT_MODEL_WINDOW_SIZE)
        aligner.build_alignment_axis_for_gene(gene, idx.samples)
        mapper = BcftoolsChainMapper(
            dataset_dir=self.dataset_dir,
            consensus_dataset_dir=self.consensus_dataset_dir,
            aligner=aligner,
        )
        return aligner, mapper

    def _read_chain_aligned_window(
        self,
        mapper: BcftoolsChainMapper,
        gene: str,
        sample_id: str,
        haplotype: str,
        start: int,
        end: int,
    ) -> str:
        entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
        if entry is None:
            return "X" * max(end - start + 1, 0)
        fasta_path = entry.get("fasta_path") or self.consensus_dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.{haplotype}.window.fixed.fa"
        fasta = _read_fasta_sequence(Path(fasta_path))
        expanded_to_source = {int(target) + 1: int(source) for source, target in zip(entry.get("copy_from_indices", []), entry.get("expanded_indices", []))}
        bases = []
        for pos in range(start, end + 1):
            source_idx = expanded_to_source.get(pos)
            bases.append("X" if source_idx is None or source_idx < 0 or source_idx >= len(fasta) else fasta[source_idx])
        return "".join(bases)

    def _read_chain_aligned_window_with_masks(
        self,
        mapper: BcftoolsChainMapper,
        gene: str,
        sample_id: str,
        haplotype: str,
        start: int,
        end: int,
    ) -> Tuple[str, List[Dict[str, object]]]:
        entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
        if entry is None:
            length = max(end - start + 1, 0)
            return "X" * length, [{"valid": 0, "ins": 0, "del": 0, "source_idx": None} for _ in range(length)]
        fasta_path = entry.get("fasta_path") or self.consensus_dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.{haplotype}.window.fixed.fa"
        fasta = _read_fasta_sequence(Path(fasta_path))
        expanded_to_source = {int(target) + 1: int(source) for source, target in zip(entry.get("copy_from_indices", []), entry.get("expanded_indices", []))}
        insertion_positions = {int(v) + 1 for v in entry.get("insertion_indices", [])}
        deletion_positions = {int(v) + 1 for v in entry.get("deletion_indices", [])}
        bases = []
        masks = []
        for pos in range(start, end + 1):
            source_idx = expanded_to_source.get(pos)
            valid = source_idx is not None and 0 <= source_idx < len(fasta)
            bases.append(fasta[source_idx] if valid else "X")
            masks.append({
                "valid": 1 if valid else 0,
                "ins": 1 if pos in insertion_positions else 0,
                "del": 1 if pos in deletion_positions else 0,
                "source_idx": source_idx if valid else None,
            })
        return "".join(bases), masks

class AlignmentViewerHandler(BaseHTTPRequestHandler):
    repository: AlignmentRepository
    max_cells: int

    def log_message(self, fmt: str, *args):
        return

    def _send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, text: str, content_type: str = "text/html; charset=utf-8", status: int = 200) -> None:
        body = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/":
                self._send_text(INDEX_HTML)
            elif parsed.path == "/api/genes":
                self._send_json({"genes": self.repository.genes_payload(), "max_cells": self.max_cells})
            elif parsed.path == "/api/samples":
                self._handle_samples(parsed.query)
            elif parsed.path == "/api/view":
                self._handle_view(parsed.query)
            elif parsed.path == "/api/variants":
                self._handle_variants(parsed.query)
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _handle_samples(self, query: str) -> None:
        qs = parse_qs(query)
        gene = qs.get("gene", [""])[0]
        self._send_json(self.repository.samples_payload(gene))

    def _handle_variants(self, query: str) -> None:
        qs = parse_qs(query)
        gene = qs.get("gene", [""])[0]
        selected = [s for raw in qs.get("samples", []) for s in raw.split(",") if s]
        start = int(qs.get("start", ["1"])[0])
        length = int(qs.get("length", ["200"])[0])
        max_variants = int(qs.get("max_variants", [str(DEFAULT_MAX_VARIANTS)])[0])
        self._send_json(self.repository.variants_payload(gene, selected, start, length, max_variants))

    def _handle_view(self, query: str) -> None:
        qs = parse_qs(query)
        gene = qs.get("gene", [""])[0]
        idx = self.repository.get(gene)
        selected = [s for raw in qs.get("samples", []) for s in raw.split(",") if s]
        if not selected:
            selected = idx.samples[:5]

        haplotype = qs.get("haplotype", ["both"])[0]
        if haplotype not in {"both", "H1", "H2"}:
            raise ValueError("haplotype deve ser both, H1 ou H2")

        start = max(int(qs.get("start", ["1"])[0]), 1)
        length = max(int(qs.get("length", ["200"])[0]), 1)
        end = min(start + length - 1, idx.expanded_length)
        variant_only = qs.get("variant_only", ["false"])[0].lower() in {"1", "true", "yes", "on"}
        model_window_size = int(qs.get("model_window_size", [str(DEFAULT_MODEL_WINDOW_SIZE)])[0])
        model_window = self.repository.model_window_payload(gene, model_window_size)

        columns = ["REF"]
        selected_column_names: List[Tuple[str, str]] = []
        for sample_id in selected:
            if haplotype in {"both", "H1"}:
                columns.append(f"{sample_id}_H1")
                selected_column_names.append((sample_id, "H1"))
            if haplotype in {"both", "H2"}:
                columns.append(f"{sample_id}_H2")
                selected_column_names.append((sample_id, "H2"))

        visible_positions = end - start + 1
        if visible_positions * len(columns) > self.max_cells:
            raise ValueError(
                f"Janela grande demais: {visible_positions} posicoes x {len(columns)} colunas. "
                f"Limite atual: {self.max_cells} celulas. Reduza length ou individuos."
            )

        ref_h1, _ref_h2 = idx.read_sample_window("REF", start, end)
        genomic_positions = self.repository.reference_genomic_positions(gene, start, end)
        window_sequences: List[Tuple[str, str]] = []
        window_masks: List[Tuple[str, List[Dict[str, object]]]] = []
        _aligner, chain_mapper = self.repository._build_chain_context(gene, selected)
        for sample_id in selected:
            h1_window, h1_masks = self.repository._read_chain_aligned_window_with_masks(chain_mapper, gene, sample_id, "H1", start, end)
            h2_window, h2_masks = self.repository._read_chain_aligned_window_with_masks(chain_mapper, gene, sample_id, "H2", start, end)
            if haplotype in {"both", "H1"}:
                window_sequences.append((f"{sample_id}_H1", h1_window))
                window_masks.append((f"{sample_id}_H1", h1_masks))
            if haplotype in {"both", "H2"}:
                window_sequences.append((f"{sample_id}_H2", h2_window))
                window_masks.append((f"{sample_id}_H2", h2_masks))

        rows = []
        for offset, pos in enumerate(range(start, end + 1)):
            ref_base = ref_h1[offset]
            bases = [ref_base]
            masks = [{"valid": 1 if ref_base != "X" else 0, "ins": 0, "del": 0, "source_idx": None}]
            has_variant = False
            for (_name, seq), (_mask_name, mask_seq) in zip(window_sequences, window_masks):
                base = seq[offset]
                bases.append(base)
                masks.append(mask_seq[offset])
                if base != ref_base:
                    has_variant = True
            if variant_only and not has_variant:
                continue
            rows.append({"pos": pos, "genomic_pos": genomic_positions[offset], "ref": ref_base, "bases": bases, "masks": masks})

        self._send_json({
            "gene": gene,
            "start": start,
            "end": end,
            "expanded_length": idx.expanded_length,
            "model_window": model_window,
            "columns": columns,
            "rows": rows,
        })


INDEX_HTML = r"""
<!doctype html>
<html lang="pt-BR">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Aligned DNA Viewer</title>
  <style>
    :root { color-scheme: dark; --bg:#11151c; --panel:#171d27; --line:#2c3544; --text:#e8edf5; --mut:#ff5d5d; --gap:#ffd166; --accent:#7cc7ff; }
    * { box-sizing: border-box; }
    body { margin: 0; background: var(--bg); color: var(--text); font-family: system-ui, sans-serif; }
    header { padding: 16px 20px; border-bottom: 1px solid var(--line); background: #0d1117; }
    h1 { margin: 0; font-size: 20px; }
    main { display: grid; grid-template-columns: 330px 1fr; min-height: calc(100vh - 58px); }
    body.sidebar-hidden main { grid-template-columns: 1fr; }
    body.sidebar-hidden aside { display: none; }
    aside { padding: 16px; border-right: 1px solid var(--line); background: var(--panel); overflow: auto; }
    section { padding: 16px; overflow: auto; }
    label { display: block; margin-top: 12px; font-size: 13px; color: #b8c2d2; }
    select, input, button { width: 100%; margin-top: 5px; padding: 8px; border: 1px solid var(--line); border-radius: 8px; background: #0d1117; color: var(--text); }
    select[multiple] { height: 220px; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
    button { cursor: pointer; background: #1f6feb; border-color: #388bfd; font-weight: 700; }
    button.secondary { background: #21262d; border-color: var(--line); }
    .sidebar-toggle { position: fixed; top: 12px; right: 16px; z-index: 20; width: auto; margin: 0; padding: 7px 11px; font-size: 12px; }
    .row { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
    .filter-panel { border: 1px solid var(--line); border-radius: 10px; padding: 9px; background: #0d1117; margin-top: 10px; }
    .filter-panel h3 { margin: 0 0 8px; font-size: 13px; color: #dce8fb; display:flex; justify-content:space-between; }
    .check-list { max-height: 150px; overflow: auto; }
    .check-list label { display:flex; gap:7px; align-items:center; padding:3px 0; color:#dce8fb; font-size:12px; }
    .check-list input { width:auto; margin:0; }
    .bottom-actions { padding: 12px; border-top: 1px solid var(--line); background: #0d1117; position: sticky; bottom: 0; }
    .hint { color: #9aa8bb; font-size: 12px; line-height: 1.4; }
    .status { margin: 0 0 12px; color: #b8c2d2; }
    .table-wrap { border: 1px solid var(--line); border-radius: 10px; overflow: auto; max-height: calc(100vh - 150px); background: #0d1117; }
    .variants-wrap { margin-top: 14px; border: 1px solid var(--line); border-radius: 10px; overflow: auto; background: #0d1117; max-height: 360px; }
    .sticky-variants { position: sticky; top: 0; z-index: 10; padding-bottom: 12px; background: var(--bg); }
    .sticky-variants .variants-head { margin-top: 0; }
    .variants-head { display:flex; justify-content:space-between; gap:12px; align-items:center; margin-top: 18px; }
    .variants-head h2 { font-size:16px; margin:0; }
    .pill { display:inline-block; padding:2px 7px; border-radius:999px; background:#1f2937; color:#dce8fb; font-size:12px; }
    .snv { color:#7ee787; font-weight:800; }
    .indel { color:#ffd166; font-weight:900; }
    .sv { color:#ff9ad5; font-weight:900; }
    table { border-collapse: collapse; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 13px; min-width: 100%; }
    th, td { border-bottom: 1px solid #202938; border-right: 1px solid #202938; padding: 3px 8px; text-align: center; white-space: nowrap; }
    th { position: sticky; top: 0; z-index: 2; background: #162033; color: var(--accent); }
    th:first-child, td:first-child { position: sticky; left: 0; z-index: 1; background: #121927; color: #b8c2d2; }
    th:first-child { z-index: 3; }
    .mut { color: var(--mut); font-weight: 800; }
    .gap { color: var(--gap); font-weight: 900; }
    .ref { color: #dce5f2; }
    tr.model-window td { background: rgba(46, 160, 67, 0.18); }
    tr.model-center td { background: rgba(124, 199, 255, 0.24); box-shadow: inset 0 1px 0 #7cc7ff, inset 0 -1px 0 #7cc7ff; }
    tr.variant-row { cursor: pointer; }
    tr.variant-row:hover td { background: rgba(124, 199, 255, 0.14); }
    tr.jump-target td { background: rgba(255, 209, 102, 0.26) !important; box-shadow: inset 0 1px 0 #ffd166, inset 0 -1px 0 #ffd166; }
    .model-badge { margin: 0 0 12px; padding: 10px 12px; border: 1px solid #2ea043; border-radius: 10px; background: rgba(46, 160, 67, 0.12); color: #c7f7d4; font-size: 13px; }
    @media (max-width: 860px) { main { grid-template-columns: 1fr; } aside { border-right: 0; border-bottom: 1px solid var(--line); } }
  </style>
</head>
<body>
  <button class="secondary sidebar-toggle" id="sidebarToggle">Ocultar filtros</button>
  <header><h1>Aligned DNA Viewer</h1></header>
  <main>
    <aside>
      <label>Gene</label>
      <select id="gene"></select>

      <div class="filter-panel">
        <h3><span>Superpopulações</span><span id="superpopCount"></span></h3>
        <input id="superpopSearch" placeholder="filtrar superpopulação" />
        <div class="row"><button class="secondary" id="superpopAll">Todas</button><button class="secondary" id="superpopNone">Nenhuma</button></div>
        <div class="check-list" id="superpopChecks"></div>
      </div>

      <div class="filter-panel">
        <h3><span>Populações</span><span id="populationCount"></span></h3>
        <input id="populationSearch" placeholder="filtrar população" />
        <div class="row"><button class="secondary" id="populationAll">Todas</button><button class="secondary" id="populationNone">Nenhuma</button></div>
        <div class="check-list" id="populationChecks"></div>
      </div>

      <div class="filter-panel">
        <h3><span>Indivíduos</span><span id="sampleCount"></span></h3>
        <input id="sampleFilter" placeholder="ex: HG00096" />
        <div class="check-list" id="sampleChecks"></div>
      </div>
      <div class="row">
        <button class="secondary" id="first5">Primeiros 5</button>
        <button class="secondary" id="selectVisible">Visíveis</button>
      </div>
      <div class="row">
        <button class="secondary" id="clear">Limpar</button>
        <button class="secondary" id="invertVisible">Inverter</button>
      </div>

      <label>Haplótipo</label>
      <select id="haplotype">
        <option value="both">H1 + H2</option>
        <option value="H1">H1</option>
        <option value="H2">H2</option>
      </select>

      <div class="row">
        <div><label>Posição inicial</label><input id="start" type="number" value="1" min="1" /></div>
        <div><label>Tamanho</label><input id="length" type="number" value="300" min="1" /></div>
      </div>

      <label><input id="variantOnly" type="checkbox" style="width:auto;margin-right:6px" /> mostrar apenas posições diferentes da referência</label>
      <label><input id="highlightModelWindow" type="checkbox" checked style="width:auto;margin-right:6px" /> destacar faixa usada no treinamento CNN</label>
      <label>Tamanho da janela do modelo</label>
      <input id="modelWindowSize" type="number" value="32768" min="1" />

      <button id="load">Carregar janela</button>
      <div class="row">
        <button class="secondary" id="prevWindow">Janela anterior</button>
        <button class="secondary" id="nextWindow">Próxima janela</button>
      </div>
      <button class="secondary" id="goModelWindow">Ir para janela da CNN</button>
      <p class="hint">Cores: vermelho = base diferente da referência; amarelo = X/gap. A tabela carrega apenas a janela selecionada.</p>
    </aside>
    <section>
      <p class="status" id="status">Carregando...</p>
      <div class="model-badge" id="modelWindowInfo">Faixa CNN: -</div>
      <div class="sticky-variants">
        <div class="variants-head"><h2>Mutações no VCF da janela selecionada</h2><span class="pill" id="variantStatus">-</span></div>
        <div class="variants-wrap"><table id="variantTable"></table></div>
      </div>
      <div class="table-wrap" id="tableWrap"><table id="table"></table><div class="bottom-actions"><button id="nextWindowBottom">Próxima janela</button></div></div>
    </section>
  </main>
  <script>
    const geneSelect = document.getElementById('gene');
    const sampleFilter = document.getElementById('sampleFilter');
    const statusEl = document.getElementById('status');
    const table = document.getElementById('table');
    const tableWrap = document.getElementById('tableWrap');
    const modelWindowInfo = document.getElementById('modelWindowInfo');
    const variantTable = document.getElementById('variantTable');
    const variantStatus = document.getElementById('variantStatus');
    const sidebarToggle = document.getElementById('sidebarToggle');
    let allSamples = [];
    let sampleRows = [];
    let populationRows = [];
    let superpopulationRows = [];
    let maxCells = 0;

    async function fetchJson(url) {
      const res = await fetch(url);
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || res.statusText);
      return data;
    }

    function setStatus(text) { statusEl.textContent = text; }

    function setSidebarHidden(hidden) {
      document.body.classList.toggle('sidebar-hidden', hidden);
      sidebarToggle.textContent = hidden ? 'Mostrar filtros' : 'Ocultar filtros';
    }

    function checkedValues(id) { return Array.from(document.querySelectorAll(`#${id} input:checked`)).map(el => el.value); }
    function updateCounts() {
      document.getElementById('superpopCount').textContent = `${checkedValues('superpopChecks').length}/${document.querySelectorAll('#superpopChecks input').length}`;
      document.getElementById('populationCount').textContent = `${checkedValues('populationChecks').length}/${document.querySelectorAll('#populationChecks input').length}`;
      document.getElementById('sampleCount').textContent = `${checkedValues('sampleChecks').length}/${document.querySelectorAll('#sampleChecks input').length}`;
    }
    function setChecks(id, checked) {
      document.querySelectorAll(`#${id} input`).forEach(el => { el.checked = checked; });
      updateCounts();
      if (id === 'superpopChecks') renderPopulationOptions();
      else if (id !== 'sampleChecks') renderSampleOptions();
    }
    function renderCheckList(id, rows, searchId, countId, formatter) {
      const previous = new Set(checkedValues(id));
      const hadPrevious = document.querySelectorAll(`#${id} input`).length > 0;
      const needle = document.getElementById(searchId).value.trim().toLowerCase();
      const filtered = rows.filter(row => !needle || String(row.id).toLowerCase().includes(needle));
      document.getElementById(id).innerHTML = filtered.map(row => {
        const checked = (!hadPrevious || previous.has(row.id)) ? 'checked' : '';
        return `<label><input type="checkbox" value="${row.id}" ${checked}>${formatter(row)}</label>`;
      }).join('');
      document.getElementById(countId).textContent = `${checkedValues(id).length}/${filtered.length}`;
      document.querySelectorAll(`#${id} input`).forEach(el => el.addEventListener('change', () => {
        updateCounts();
        if (id === 'superpopChecks') renderPopulationOptions();
        else renderSampleOptions();
      }));
    }

    function renderPopulationOptions() {
      const selected = new Set(checkedValues('populationChecks'));
      const hadPrevious = document.querySelectorAll('#populationChecks input').length > 0;
      const superpops = new Set(checkedValues('superpopChecks'));
      const counts = new Map();
      for (const row of sampleRows) {
        if (superpops.size && !superpops.has(row.superpopulation || '')) continue;
        const pop = row.population || '';
        if (!pop) continue;
        counts.set(pop, (counts.get(pop) || 0) + 1);
      }
      const rows = Array.from(counts.entries()).sort((a, b) => a[0].localeCompare(b[0])).map(([id, count]) => ({id, count}));
      const needle = document.getElementById('populationSearch').value.trim().toLowerCase();
      const filtered = rows.filter(row => !needle || row.id.toLowerCase().includes(needle));
      document.getElementById('populationChecks').innerHTML = filtered.map(row => {
        const checked = (!hadPrevious || selected.has(row.id)) ? 'checked' : '';
        return `<label><input type="checkbox" value="${row.id}" ${checked}>${row.id} (${row.count})</label>`;
      }).join('');
      document.getElementById('populationCount').textContent = `${checkedValues('populationChecks').length}/${filtered.length}`;
      document.querySelectorAll('#populationChecks input').forEach(el => el.addEventListener('change', () => { updateCounts(); renderSampleOptions(); }));
      renderSampleOptions();
    }
    function renderSampleOptions() {
      const previous = new Set(checkedValues('sampleChecks'));
      const hadPrevious = document.querySelectorAll('#sampleChecks input').length > 0;
      const needle = sampleFilter.value.trim().toLowerCase();
      const superpops = new Set(checkedValues('superpopChecks'));
      const pops = new Set(checkedValues('populationChecks'));
      let rows = sampleRows;
      if (superpops.size) rows = rows.filter(row => superpops.has(row.superpopulation || ''));
      if (pops.size) rows = rows.filter(row => pops.has(row.population || ''));
      if (needle) rows = rows.filter(row => row.sample_id.toLowerCase().includes(needle));
      const visible = rows.slice(0, 500);
      document.getElementById('sampleChecks').innerHTML = visible.map((row, idx) => {
        const checked = (!hadPrevious ? idx < 5 : previous.has(row.sample_id)) ? 'checked' : '';
        return `<label><input type="checkbox" value="${row.sample_id}" ${checked}>${row.sample_id} | ${row.superpopulation || '-'} / ${row.population || '-'}</label>`;
      }).join('');
      document.querySelectorAll('#sampleChecks input').forEach(el => el.addEventListener('change', updateCounts));
      updateCounts();
    }

    async function loadGenes() {
      const data = await fetchJson('/api/genes');
      maxCells = data.max_cells;
      geneSelect.innerHTML = '';
      for (const g of data.genes) {
        const opt = document.createElement('option');
        opt.value = g.id;
        opt.textContent = `${g.gene} (${g.samples} indivíduos, L=${g.expanded_length})`;
        geneSelect.appendChild(opt);
      }
      await loadSamples();
    }

    async function loadSamples() {
      const gene = geneSelect.value;
      const data = await fetchJson(`/api/samples?gene=${encodeURIComponent(gene)}`);
      allSamples = data.samples;
      sampleRows = data.sample_rows || data.samples.map(sample => ({sample_id: sample, population: '', superpopulation: ''}));
      populationRows = Object.entries(data.populations || {}).map(([id, count]) => ({id, count}));
      superpopulationRows = Object.entries(data.superpopulations || {}).map(([id, count]) => ({id, count}));
      renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`);
      renderPopulationOptions();
      setStatus(`${gene}: ${allSamples.length} indivíduos disponíveis. Limite de renderização: ${maxCells} células.`);
    }

    function selectFirstN(n) {
      document.querySelectorAll('#sampleChecks input').forEach((el, i) => { el.checked = i < n; });
      updateCounts();
    }

    function selectVisible(checked) {
      setChecks('sampleChecks', checked);
    }

    function invertVisible() {
      document.querySelectorAll('#sampleChecks input').forEach(el => { el.checked = !el.checked; });
      updateCounts();
    }

    function moveWindow(direction) {
      const startEl = document.getElementById('start');
      const length = Math.max(1, Number(document.getElementById('length').value || 1));
      const current = Math.max(1, Number(startEl.value || 1));
      startEl.value = Math.max(1, current + direction * length);
      loadView().catch(e => setStatus(e.message));
    }

    function selectedSamples() {
      return checkedValues('sampleChecks');
    }

    function cellClass(base, ref, colIndex) {
      if (base === 'X') return 'gap';
      if (colIndex > 0 && base !== ref) return 'mut';
      return colIndex === 0 ? 'ref' : '';
    }

    function variantTypeClass(type) {
      if (type === 'SNV') return 'snv';
      if (type === 'INDEL') return 'indel';
      if (type === 'SV') return 'sv';
      return '';
    }

    function firstExpandedPosition(span) {
      const text = String(span || '').split('-', 1)[0];
      const pos = Number(text);
      return Number.isFinite(pos) ? pos : null;
    }

    function jumpToAlignedPosition(pos) {
      const target = table.querySelector(`tr[data-pos="${pos}"]`);
      if (!target) {
        setStatus(`Posição ${pos} não está renderizada na janela atual.`);
        return;
      }
      table.querySelectorAll('tr.jump-target').forEach(row => row.classList.remove('jump-target'));
      target.classList.add('jump-target');
      const top = target.offsetTop - tableWrap.clientHeight / 2 + target.clientHeight / 2;
      tableWrap.scrollTo({top: Math.max(0, top), behavior: 'smooth'});
    }

    function maskTooltip(columnName, row, idx) {
      const mask = (row.masks || [])[idx] || {};
      const source = mask.source_idx == null ? '-' : mask.source_idx;
      return [
        `${columnName}`,
        `index: ${row.pos}`,
        `ref genome pos: ${row.genomic_pos || '-'}`,
        `base: ${(row.bases || [])[idx] || '-'}`,
        `valid_mask: ${mask.valid ?? '-'}`,
        `insertion_mask: ${mask.ins ?? '-'}`,
        `deletion_mask: ${mask.del ?? '-'}`,
        `fasta/prediction index: ${source}`,
      ].join('\n');
    }

    function renderTable(data) {
      table.innerHTML = '';
      const modelWindow = data.model_window || {};
      const highlightModel = document.getElementById('highlightModelWindow').checked;
      modelWindowInfo.textContent = `Faixa CNN (${modelWindow.policy || '-'}) | index ${modelWindow.start || '-'}-${modelWindow.end || '-'} | centro ${modelWindow.center_expanded_pos || '-'} | tamanho ${modelWindow.window_size || '-'}`;
      const thead = document.createElement('thead');
      const hrow = document.createElement('tr');
      ['index', 'ref genome pos', ...data.columns].forEach(name => {
        const th = document.createElement('th');
        th.textContent = name;
        hrow.appendChild(th);
      });
      thead.appendChild(hrow);
      table.appendChild(thead);

      const tbody = document.createElement('tbody');
      for (const row of data.rows) {
        const tr = document.createElement('tr');
        tr.dataset.pos = row.pos;
        if (highlightModel && modelWindow.start && row.pos >= modelWindow.start && row.pos <= modelWindow.end) {
          tr.classList.add('model-window');
        }
        if (highlightModel && modelWindow.center_expanded_pos && row.pos === modelWindow.center_expanded_pos) {
          tr.classList.add('model-center');
        }
        const pos = document.createElement('td');
        pos.textContent = row.pos;
        tr.appendChild(pos);
        const genomicPos = document.createElement('td');
        genomicPos.textContent = row.genomic_pos || '';
        tr.appendChild(genomicPos);
        row.bases.forEach((base, idx) => {
          const td = document.createElement('td');
          td.textContent = base;
          td.className = cellClass(base, row.ref, idx);
          td.title = maskTooltip(data.columns[idx], row, idx);
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      }
      table.appendChild(tbody);
    }

    function renderVariants(data) {
      variantTable.innerHTML = '';
      variantStatus.textContent = `${data.variants.length}${data.truncated ? '+' : ''} variantes | ${data.reference_region || 'sem região'}`;
      const thead = document.createElement('thead');
      const hrow = document.createElement('tr');
      ['pos alinhada', 'posição VCF', 'tipo', 'REF', 'ALT', 'GT alternativo', 'genótipos selecionados', 'ID', 'FILTER'].forEach(name => {
        const th = document.createElement('th');
        th.textContent = name;
        hrow.appendChild(th);
      });
      thead.appendChild(hrow);
      variantTable.appendChild(thead);

      const tbody = document.createElement('tbody');
      if (!data.variants.length) {
        const tr = document.createElement('tr');
        const td = document.createElement('td');
        td.colSpan = 9;
        td.textContent = 'Nenhuma variante do VCF encontrada nesta janela.';
        tr.appendChild(td);
        tbody.appendChild(tr);
      }
      for (const row of data.variants) {
        const tr = document.createElement('tr');
        const alignedPos = firstExpandedPosition(row.expanded_span);
        if (alignedPos !== null) {
          tr.classList.add('variant-row');
          tr.title = `Ir para posição alinhada ${alignedPos}`;
          tr.addEventListener('click', () => jumpToAlignedPosition(alignedPos));
        }
        const gtText = Object.entries(row.gts || {}).map(([sample, gt]) => `${sample}:${gt}`).join('  ');
        const altText = (row.alt_samples || []).join(', ');
        const values = [row.expanded_span, row.genomic_pos, row.type, row.ref, row.alt, `${row.alt_count}/${row.called_count}${altText ? ' | ' + altText : ''}`, gtText, row.id || '-', row.filter || '-'];
        values.forEach((value, idx) => {
          const td = document.createElement('td');
          td.textContent = value;
          if (idx === 2) td.className = variantTypeClass(row.type);
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      }
      variantTable.appendChild(tbody);
    }

    async function loadView() {
      const samples = selectedSamples();
      if (!samples.length) { setStatus('Selecione pelo menos um indivíduo.'); return; }
      const params = new URLSearchParams();
      params.set('gene', geneSelect.value);
      params.set('samples', samples.join(','));
      params.set('haplotype', document.getElementById('haplotype').value);
      params.set('start', document.getElementById('start').value);
      params.set('length', document.getElementById('length').value);
      params.set('variant_only', document.getElementById('variantOnly').checked ? 'true' : 'false');
      params.set('model_window_size', document.getElementById('modelWindowSize').value);
      setStatus('Carregando janela...');
      const data = await fetchJson(`/api/view?${params.toString()}`);
      renderTable(data);
      variantStatus.textContent = 'Carregando VCF...';
      const variants = await fetchJson(`/api/variants?${params.toString()}`);
      renderVariants(variants);
      setStatus(`${data.gene}: posições ${data.start}-${data.end}; ${data.rows.length} linhas renderizadas.`);
    }

    async function goToModelWindow() {
      const params = new URLSearchParams();
      params.set('gene', geneSelect.value);
      params.set('samples', selectedSamples().slice(0, 1).join(',') || (allSamples[0] || ''));
      params.set('start', '1');
      params.set('length', '1');
      params.set('model_window_size', document.getElementById('modelWindowSize').value);
      const data = await fetchJson(`/api/view?${params.toString()}`);
      const mw = data.model_window;
      document.getElementById('start').value = mw.start;
      document.getElementById('length').value = mw.window_size;
      await loadView();
    }

    geneSelect.addEventListener('change', loadSamples);
    document.getElementById('superpopSearch').addEventListener('input', () => {
      renderCheckList('superpopChecks', superpopulationRows, 'superpopSearch', 'superpopCount', row => `${row.id} (${row.count})`);
      renderPopulationOptions();
    });
    document.getElementById('populationSearch').addEventListener('input', renderPopulationOptions);
    sampleFilter.addEventListener('input', renderSampleOptions);
    document.getElementById('superpopAll').addEventListener('click', () => setChecks('superpopChecks', true));
    document.getElementById('superpopNone').addEventListener('click', () => setChecks('superpopChecks', false));
    document.getElementById('populationAll').addEventListener('click', () => setChecks('populationChecks', true));
    document.getElementById('populationNone').addEventListener('click', () => setChecks('populationChecks', false));
    document.getElementById('first5').addEventListener('click', () => selectFirstN(5));
    document.getElementById('selectVisible').addEventListener('click', () => selectVisible(true));
    document.getElementById('clear').addEventListener('click', () => selectFirstN(0));
    document.getElementById('invertVisible').addEventListener('click', invertVisible);
    document.getElementById('load').addEventListener('click', () => loadView().catch(e => setStatus(e.message)));
    document.getElementById('prevWindow').addEventListener('click', () => moveWindow(-1));
    document.getElementById('nextWindow').addEventListener('click', () => moveWindow(1));
    document.getElementById('nextWindowBottom').addEventListener('click', () => moveWindow(1));
    document.getElementById('highlightModelWindow').addEventListener('change', () => loadView().catch(e => setStatus(e.message)));
    document.getElementById('goModelWindow').addEventListener('click', () => goToModelWindow().catch(e => setStatus(e.message)));
    sidebarToggle.addEventListener('click', () => {
      const hidden = !document.body.classList.contains('sidebar-hidden');
      setSidebarHidden(hidden);
      localStorage.setItem('alignedDnaSidebarHidden', hidden ? '1' : '0');
    });
    setSidebarHidden(localStorage.getItem('alignedDnaSidebarHidden') === '1');

    loadGenes().then(loadView).catch(e => setStatus(e.message));
  </script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Web viewer for aligned DNA TSV files")
    parser.add_argument(
        "tsv_root",
        type=Path,
        help="TSV file or directory containing one .tsv per gene",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--max-cells", type=int, default=DEFAULT_MAX_CELLS)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR, help="Dataset dir used to load population/superpopulation metadata")
    parser.set_defaults(alignment_mapping="bcftools_chain")
    parser.add_argument("--consensus-dataset-dir", type=Path, default=DEFAULT_CONSENSUS_DATASET_DIR)
    args = parser.parse_args()

    repository = AlignmentRepository(
        args.tsv_root.resolve(),
        args.dataset_dir.resolve() if args.dataset_dir else None,
        alignment_mapping=args.alignment_mapping,
        consensus_dataset_dir=args.consensus_dataset_dir.resolve() if args.consensus_dataset_dir else None,
    )

    class Handler(AlignmentViewerHandler):
        pass

    Handler.repository = repository
    Handler.max_cells = args.max_cells

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Aligned DNA viewer: http://{args.host}:{args.port}")
    print(f"TSV root: {args.tsv_root.resolve()}")
    print(f"Alignment mapping: {args.alignment_mapping}")
    print(f"Genes: {', '.join(repository.indexes.keys())}")
    server.serve_forever()


if __name__ == "__main__":
    main()
