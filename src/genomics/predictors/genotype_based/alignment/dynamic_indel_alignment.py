"""Dynamic INDEL-aware expanded-axis construction from window VCFs."""

from __future__ import annotations

import json
import gzip
import hashlib
import shutil
import subprocess
import tempfile
import time
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ALIGNMENT_ALGORITHM_VERSION = "dynamic_indel_ref_window_v4"
CENTRAL_ALIGNMENT_CACHE_VERSION = "bcftools_chain_alignment_v1"
DEFAULT_ALIGNMENT_CENTER_WINDOW_SIZE = 32_768


def _normalize_alt_alleles(alt_raw: str) -> List[str]:
    return [a for a in str(alt_raw).split(",") if a and a != "."]


def _classify_length_behavior(ref: str, alt: str) -> str:
    alt_alleles = _normalize_alt_alleles(alt)
    if not alt_alleles or len(alt_alleles) != 1:
        return "special_unhandled"
    alt0 = alt_alleles[0]
    if str(ref).startswith("<"):
        return "special_unhandled"
    if alt0 in {"<DEL>", "<del>"}:
        return "length_change"
    if str(alt0).startswith("<"):
        return "special_unhandled"
    if "[" in alt0 or "]" in alt0:
        return "special_unhandled"
    if len(ref) == len(alt0):
        return "same_length"
    return "length_change"


def _parse_gt(gt: str) -> Tuple[str, str]:
    text = str(gt).replace("/", "|")
    if "|" not in text:
        return "0", "0"
    left, right = text.split("|", 1)
    return left, right


def _allele_sequence(ref: str, alt: str, token: str) -> str:
    if token in {".", "nan", "None"}:
        return "N"
    if token == "0":
        return str(ref)
    try:
        idx = int(token) - 1
    except ValueError:
        return str(ref)
    alt_alleles = _normalize_alt_alleles(alt)
    if 0 <= idx < len(alt_alleles):
        allele = str(alt_alleles[idx])
        if allele in {"<DEL>", "<del>"}:
            return str(ref[:1])
        return allele
    return str(ref)


def _read_fasta_sequence(path: Path) -> str:
    lines: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines).upper()


def _path_fingerprint(path: Path) -> Dict[str, object]:
    path = Path(path)
    try:
        stat = path.stat()
    except FileNotFoundError:
        return {"path": str(path), "exists": False}
    return {
        "path": str(path.resolve()),
        "exists": True,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _json_key_dict_to_int(payload: Dict) -> Dict[int, object]:
    return {int(key): value for key, value in payload.items()}


def _resolve_ref_index(ref_sequence: str, pos0: int, ref: str) -> int:
    """Use the local FASTA as the source of truth when VCF POS is off by a small amount."""
    ref = str(ref).upper()
    if not ref_sequence or not ref or ref.startswith("<"):
        return pos0
    if 0 <= pos0 <= len(ref_sequence) - len(ref) and ref_sequence[pos0:pos0 + len(ref)] == ref:
        return pos0
    for delta in (1, -1, 2, -2, 3, -3):
        candidate = pos0 + delta
        if 0 <= candidate <= len(ref_sequence) - len(ref) and ref_sequence[candidate:candidate + len(ref)] == ref:
            return candidate
    return pos0


def _stream_gene_variants(vcf_path: str, sample_ids: List[str], region: str):
    if not sample_ids:
        return

    if shutil.which("bcftools") is None:
        yield from _stream_gene_variants_from_gzip(vcf_path, sample_ids, region)
        return

    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"
    sample_file = tempfile.NamedTemporaryFile("w", prefix="genomics_bcftools_samples_", suffix=".txt", delete=False)
    sample_file_path = Path(sample_file.name)
    try:
        sample_file.write("\n".join(sample_ids))
        sample_file.write("\n")
        sample_file.close()

        cmd = [
            "bcftools",
            "query",
            "-S",
            str(sample_file_path),
            "-r",
            region,
            "-f",
            fmt,
            vcf_path,
        ]

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        try:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                chrom, pos, variant_id, ref, alt = fields[:5]
                yield {
                    "chrom": chrom,
                    "pos_1based": int(pos),
                    "id": variant_id,
                    "ref": ref,
                    "alt": alt,
                    "gts": fields[5:],
                }
        finally:
            if proc.stdout is not None:
                proc.stdout.close()
            returncode = proc.wait()
            stderr = proc.stderr.read() if proc.stderr is not None else ""
            if proc.stderr is not None:
                proc.stderr.close()
            if returncode != 0:
                raise RuntimeError(
                    "bcftools query falhou ao construir eixo INDEL "
                    f"(region={region}, samples={len(sample_ids)}, vcf={vcf_path}): {stderr.strip()}"
                )
    finally:
        sample_file_path.unlink(missing_ok=True)


def _parse_region(region: str) -> Tuple[str, int, int]:
    chrom, span = region.split(":", 1)
    start_raw, end_raw = span.split("-", 1)
    return chrom, int(start_raw), int(end_raw)


def _stream_gene_variants_from_gzip(vcf_path: str, sample_ids: List[str], region: str):
    chrom, start_1based, end_1based = _parse_region(region)
    wanted_samples = list(sample_ids)
    sample_column_indices: List[int] = []

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_to_col = {sample_id: idx for idx, sample_id in enumerate(header[9:])}
                missing = [sample_id for sample_id in wanted_samples if sample_id not in sample_to_col]
                if missing:
                    raise KeyError(f"Samples ausentes no VCF: {missing}")
                sample_column_indices = [9 + sample_to_col[sample_id] for sample_id in wanted_samples]
                break

        if not sample_column_indices:
            raise RuntimeError(f"Header #CHROM nao encontrado em {vcf_path}")

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
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "gts": gts,
            }


class DynamicIndelAligner:
    """Builds expanded-axis mappings dynamically for the active experiment view."""

    def __init__(
        self,
        dataset_dir: Path,
        selected_sample_ids: Optional[Iterable[str]] = None,
        sample_cache_limit: int = 8,
        sample_prefetch_size: int = 8,
        persistent_cache_dir: Optional[Path] = None,
        center_window_size: Optional[int] = DEFAULT_ALIGNMENT_CENTER_WINDOW_SIZE,
    ):
        self.dataset_dir = Path(dataset_dir)
        self.selected_sample_ids = set(selected_sample_ids or [])
        self.center_window_size = center_window_size
        self.persistent_cache_root = (
            Path(persistent_cache_dir)
            if persistent_cache_dir
            else self.dataset_dir / "alignment_cache" / CENTRAL_ALIGNMENT_CACHE_VERSION / "sample_sets"
        )
        self._gene_cache: Dict[str, Dict[str, object]] = {}
        self._sample_entry_cache: Dict[str, OrderedDict[Tuple[str, str], Dict[str, object]]] = {}
        self._sample_cache_limit = sample_cache_limit
        self._sample_prefetch_size = sample_prefetch_size
        self.profile_stats = {
            "bcftools_query_s": 0.0,
            "entry_build_s": 0.0,
            "query_calls": 0,
            "entry_build_calls": 0,
        }

    def _sample_set_key(self) -> str:
        if not self.selected_sample_ids:
            return "samples_unspecified"
        payload = "\n".join(sorted(str(sample_id) for sample_id in self.selected_sample_ids))
        digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:16]
        return f"samples_{len(self.selected_sample_ids)}_{digest}"

    def _gene_cache_dir(self, gene_name: str) -> Path:
        return self.persistent_cache_root / self._sample_set_key() / "genes" / gene_name

    def _sample_set_manifest_path(self) -> Path:
        return self.persistent_cache_root / self._sample_set_key() / "sample_set.json"

    def _write_sample_set_manifest(self) -> None:
        path = self._sample_set_manifest_path()
        sample_ids = sorted(str(sample_id) for sample_id in self.selected_sample_ids)
        payload = {
            "cache_version": CENTRAL_ALIGNMENT_CACHE_VERSION,
            "axis_algorithm_version": ALIGNMENT_ALGORITHM_VERSION,
            "sample_set_key": self._sample_set_key(),
            "sample_count": len(sample_ids),
            "sample_ids_sha1": hashlib.sha1("\n".join(sample_ids).encode("utf-8")).hexdigest(),
            "center_window_size": self.center_window_size,
            "dataset_dir": str(self.dataset_dir.resolve()),
        }
        self._atomic_write_json(path, payload)

    def _axis_cache_path(self, gene_name: str) -> Path:
        return self._gene_cache_dir(gene_name) / "axis.json"

    def _sample_cache_path(self, gene_name: str, sample_id: str) -> Path:
        return self._gene_cache_dir(gene_name) / "samples" / f"{sample_id}.json"

    def _atomic_write_json(self, path: Path, payload: Dict[str, object]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = path.with_suffix(path.suffix + f".tmp.{time.time_ns()}")
        with open(tmp_path, "w") as f:
            json.dump(payload, f, separators=(",", ":"))
        tmp_path.replace(path)

    def _fingerprints_for_gene(self, gene_name: str, gene_payload: Optional[Dict[str, object]] = None) -> Dict[str, object]:
        if gene_payload is None:
            window_meta = self._load_window_metadata(gene_name)
            raw_variant_source = window_meta.get("raw_variant_source") or {}
            vcf_path = raw_variant_source.get("vcf_path", "")
        else:
            vcf_path = str(gene_payload.get("vcf_path", ""))
        ref_sequence_path = self.dataset_dir / "references" / "windows" / gene_name / "ref.window.fa"
        window_metadata_path = self.dataset_dir / "references" / "windows" / gene_name / "window_metadata.json"
        return {
            "algorithm_version": ALIGNMENT_ALGORITHM_VERSION,
            "dataset_dir": str(self.dataset_dir.resolve()),
            "gene": gene_name,
            "sample_set_key": self._sample_set_key(),
            "sample_count": len(self.selected_sample_ids),
            "center_window_size": self.center_window_size,
            "window_metadata": _path_fingerprint(window_metadata_path),
            "ref_window_fasta": _path_fingerprint(ref_sequence_path),
            "vcf": _path_fingerprint(Path(vcf_path)) if vcf_path else {"path": "", "exists": False},
        }

    def _cache_meta_matches(self, meta: Dict[str, object], expected: Dict[str, object]) -> bool:
        return meta == expected

    def _load_window_metadata(self, gene_name: str) -> Dict[str, object]:
        meta_path = self.dataset_dir / "references" / "windows" / gene_name / "window_metadata.json"
        with open(meta_path) as f:
            return json.load(f)

    def _build_region(self, window_meta: Dict[str, object]) -> str:
        start_1based = int(window_meta["start"])
        end_1based = int(window_meta["end"])
        return f"{window_meta['chromosome']}:{start_1based}-{end_1based}"

    def _alignment_ref_span(self, window_meta: Dict[str, object]) -> Tuple[int, int, int, int]:
        full_start_1based = int(window_meta["start"])
        full_end_1based = int(window_meta["end"])
        full_ref_length = full_end_1based - full_start_1based + 1
        if self.center_window_size is None or int(self.center_window_size) <= 0:
            return full_start_1based, full_end_1based, 0, full_ref_length
        size = min(int(self.center_window_size), full_ref_length)
        center_offset = full_ref_length // 2
        ref_start_offset = max(0, center_offset - (size // 2))
        ref_end_offset = min(full_ref_length, ref_start_offset + size)
        if ref_end_offset - ref_start_offset < size:
            ref_start_offset = max(0, ref_end_offset - size)
        region_start_1based = full_start_1based + ref_start_offset
        region_end_1based = full_start_1based + ref_end_offset - 1
        return region_start_1based, region_end_1based, ref_start_offset, ref_end_offset - ref_start_offset

    def _build_alignment_region(self, window_meta: Dict[str, object]) -> str:
        start_1based, end_1based, _offset, _length = self._alignment_ref_span(window_meta)
        return f"{window_meta['chromosome']}:{start_1based}-{end_1based}"

    def _load_gene_mapping(self, gene_name: str) -> Dict[str, object]:
        cached = self._gene_cache.get(gene_name)
        if cached is not None:
            return cached

        window_meta = self._load_window_metadata(gene_name)
        raw_variant_source = window_meta.get("raw_variant_source") or {}
        vcf_path = raw_variant_source.get("vcf_path")
        if not vcf_path:
            raise KeyError(f"raw_variant_source.vcf_path ausente para {gene_name}")
        vcf_path = str(vcf_path)

        start_1based = int(window_meta["start"])
        end_1based = int(window_meta["end"])
        full_ref_length = end_1based - start_1based + 1
        alignment_start_1based, alignment_end_1based, ref_start_offset, ref_length = self._alignment_ref_span(window_meta)
        ref_sequence_path = self.dataset_dir / "references" / "windows" / gene_name / "ref.window.fa"
        ref_sequence = _read_fasta_sequence(ref_sequence_path) if ref_sequence_path.exists() else ""
        if ref_sequence and ref_length != full_ref_length:
            ref_sequence = ref_sequence[ref_start_offset:ref_start_offset + ref_length]

        payload = {
            "vcf_path": vcf_path,
            "region": self._build_region(window_meta) if ref_start_offset else self._build_alignment_region(window_meta),
            "full_start_1based": start_1based,
            "full_ref_length": full_ref_length,
            "alignment_start_1based": alignment_start_1based,
            "alignment_end_1based": alignment_end_1based,
            "ref_start_offset": ref_start_offset,
            "start_1based": alignment_start_1based,
            "ref_length": ref_length,
            "ref_sequence": ref_sequence,
            "expanded_length": ref_length,
            "expanded_index_map": {i: i for i in range(ref_length)},
            "insertion_slots_by_ref": {},
            "fingerprints": self._fingerprints_for_gene(gene_name, {
                "vcf_path": vcf_path,
            }),
        }
        self._gene_cache[gene_name] = payload
        self._sample_entry_cache[gene_name] = OrderedDict()
        return payload

    def _load_axis_from_persistent_cache(self, gene_name: str, gene_payload: Dict[str, object]) -> bool:
        path = self._axis_cache_path(gene_name)
        if not path.exists():
            return False
        try:
            with open(path) as f:
                payload = json.load(f)
        except Exception:
            return False
        if not self._cache_meta_matches(payload.get("fingerprints", {}), gene_payload.get("fingerprints", {})):
            return False
        axis = payload.get("axis", {})
        try:
            gene_payload["expanded_length"] = int(axis["expanded_length"])
            gene_payload["expanded_index_map"] = _json_key_dict_to_int(axis.get("expanded_index_map", {}))
            gene_payload["insertion_slots_by_ref"] = _json_key_dict_to_int(axis.get("insertion_slots_by_ref", {}))
            gene_payload["axis_finalized"] = True
        except Exception:
            return False
        return True

    def _save_axis_to_persistent_cache(self, gene_name: str, gene_payload: Dict[str, object]) -> None:
        payload = {
            "fingerprints": gene_payload.get("fingerprints", {}),
            "axis": {
                "expanded_length": int(gene_payload["expanded_length"]),
                "ref_length": int(gene_payload["ref_length"]),
                "expanded_index_map": gene_payload.get("expanded_index_map", {}),
                "insertion_slots_by_ref": gene_payload.get("insertion_slots_by_ref", {}),
            },
        }
        self._atomic_write_json(self._axis_cache_path(gene_name), payload)

    def _load_sample_from_persistent_cache(self, gene_name: str, sample_id: str, gene_payload: Dict[str, object]) -> Optional[Dict[str, object]]:
        if not gene_payload.get("axis_finalized"):
            self._load_axis_from_persistent_cache(gene_name, gene_payload)
        path = self._sample_cache_path(gene_name, sample_id)
        if not path.exists():
            return None
        try:
            with open(path) as f:
                payload = json.load(f)
        except Exception:
            return None
        if not self._cache_meta_matches(payload.get("fingerprints", {}), gene_payload.get("fingerprints", {})):
            return None
        if int(payload.get("expanded_length", -1)) != int(gene_payload.get("expanded_length", -2)):
            return None
        entries = payload.get("entries")
        return entries if isinstance(entries, dict) else None

    def _save_sample_to_persistent_cache(self, gene_name: str, sample_id: str, gene_payload: Dict[str, object], entries: Dict[str, object]) -> None:
        payload = {
            "fingerprints": gene_payload.get("fingerprints", {}),
            "expanded_length": int(gene_payload["expanded_length"]),
            "sample_id": sample_id,
            "entries": entries,
        }
        self._atomic_write_json(self._sample_cache_path(gene_name, sample_id), payload)

    def _build_sample_entries_for_gene(self, gene_name: str, sample_ids: List[str]) -> None:
        gene_payload = self._load_gene_mapping(gene_name)
        if not gene_payload.get("axis_finalized"):
            self._load_axis_from_persistent_cache(gene_name, gene_payload)
        cache = self._sample_entry_cache[gene_name]
        missing_sample_ids = [sample_id for sample_id in sample_ids if sample_id not in cache]
        for sample_id in list(missing_sample_ids):
            cached_entry = self._load_sample_from_persistent_cache(gene_name, sample_id, gene_payload)
            if cached_entry is None:
                continue
            cache[sample_id] = cached_entry
            cache.move_to_end(sample_id)
            while len(cache) > self._sample_cache_limit:
                cache.popitem(last=False)
        missing_sample_ids = [sample_id for sample_id in sample_ids if sample_id not in cache]
        if not missing_sample_ids:
            return

        selected_ids = sorted(self.selected_sample_ids) if self.selected_sample_ids else missing_sample_ids
        first_missing = missing_sample_ids[0]
        if first_missing in selected_ids:
            start_idx = selected_ids.index(first_missing)
            prefetch_ids = [sid for sid in selected_ids[start_idx:start_idx + self._sample_prefetch_size] if sid not in cache]
        else:
            prefetch_ids = missing_sample_ids[:self._sample_prefetch_size]
        missing_sample_ids = prefetch_ids

        t0 = time.perf_counter()
        insertion_after_ref_index: Dict[int, int] = {}
        start_1based = int(gene_payload["start_1based"])
        full_start_1based = int(gene_payload.get("full_start_1based", start_1based))
        ref_length = int(gene_payload["ref_length"])
        ref_sequence = str(gene_payload.get("ref_sequence", ""))
        axis_finalized = bool(gene_payload.get("axis_finalized"))
        sample_entries: Dict[str, Dict[str, Dict[str, object]]] = {}

        for sample_idx, sample_id in enumerate(missing_sample_ids):
            hap_entries: Dict[str, Dict[str, object]] = {}
            for haplotype, gt_index in (("H1", 0), ("H2", 1)):
                deletion_positions = set()
                insertion_specs: List[Tuple[int, int]] = []
                hap_entries[haplotype] = {
                    "deletion_positions": deletion_positions,
                    "insertion_specs": insertion_specs,
                    "source_start_idx": int(gene_payload.get("ref_start_offset", 0)),
                }
            sample_entries[sample_id] = hap_entries

        query_t0 = time.perf_counter()
        for row in _stream_gene_variants(gene_payload["vcf_path"], missing_sample_ids, gene_payload["region"]):
            if _classify_length_behavior(str(row["ref"]), str(row["alt"])) != "length_change":
                continue
            row_gts = row.get("gts", [])
            ref = str(row["ref"])
            pos0_full = int(row["pos_1based"]) - full_start_1based
            ref_start_offset = int(gene_payload.get("ref_start_offset", 0))
            pos0 = _resolve_ref_index(ref_sequence, pos0_full - ref_start_offset, ref)
            before_alignment = pos0 < 0
            if not before_alignment and pos0 >= ref_length:
                continue

            for sample_idx, sample_id in enumerate(missing_sample_ids):
                if sample_idx >= len(row_gts):
                    continue
                gt_left, gt_right = _parse_gt(str(row_gts[sample_idx]))
                for haplotype, token in (("H1", gt_left), ("H2", gt_right)):
                    allele = _allele_sequence(ref, str(row["alt"]), token)
                    delta_len = len(allele) - len(ref)
                    if before_alignment:
                        hap_state = sample_entries[sample_id][haplotype]
                        span_end = pos0 + len(ref)
                        if span_end <= 0:
                            hap_state["source_start_idx"] += delta_len
                        elif delta_len < 0:
                            hap_state["source_start_idx"] += max(delta_len, -max(0, -pos0))
                        continue
                    if len(allele) == len(ref):
                        continue
                    hap_state = sample_entries[sample_id][haplotype]
                    if len(allele) < len(ref):
                        for offset in range(1, len(ref)):
                            ref_idx = pos0 + offset
                            if 0 <= ref_idx < ref_length:
                                hap_state["deletion_positions"].add(ref_idx)
                    else:
                        ins_len = delta_len
                        if ins_len > 0:
                            hap_state["insertion_specs"].append((pos0, ins_len))
                            if not axis_finalized:
                                insertion_after_ref_index[pos0] = max(insertion_after_ref_index.get(pos0, 0), ins_len)

        self.profile_stats["bcftools_query_s"] += time.perf_counter() - query_t0
        self.profile_stats["query_calls"] += 1

        if axis_finalized:
            expanded_index_map = dict(gene_payload["expanded_index_map"])
            insertion_slots_by_ref = dict(gene_payload["insertion_slots_by_ref"])
        else:
            expanded_index_map: Dict[int, int] = {}
            insertion_slots_by_ref: Dict[int, List[int]] = {}
            expanded_cursor = 0
            for ref_idx in range(ref_length):
                expanded_index_map[ref_idx] = expanded_cursor
                expanded_cursor += 1
                ins_slots = insertion_after_ref_index.get(ref_idx, 0)
                if ins_slots > 0:
                    insertion_slots_by_ref[ref_idx] = list(range(expanded_cursor, expanded_cursor + ins_slots))
                    expanded_cursor += ins_slots

            gene_payload["expanded_length"] = expanded_cursor
            gene_payload["expanded_index_map"] = expanded_index_map
            gene_payload["insertion_slots_by_ref"] = insertion_slots_by_ref

        for sample_id, hap_entries in sample_entries.items():
            per_hap: Dict[str, object] = {}
            for haplotype, hap_state in hap_entries.items():
                deletion_positions = hap_state["deletion_positions"]
                insertion_specs = hap_state["insertion_specs"]
                insertion_slots: List[int] = []
                for ref_idx, ins_len in insertion_specs:
                    insertion_slots.extend(insertion_slots_by_ref.get(ref_idx, [])[:ins_len])

                copy_from_indices: List[int] = []
                expanded_indices: List[int] = []
                source_cursor = 0
                insertion_len_by_ref = {ref_idx: ins_len for ref_idx, ins_len in insertion_specs}
                for ref_idx in range(ref_length):
                    if ref_idx not in deletion_positions:
                        copy_from_indices.append(source_cursor)
                        expanded_indices.append(expanded_index_map[ref_idx])
                        source_cursor += 1
                    ins_len = insertion_len_by_ref.get(ref_idx, 0)
                    if ins_len > 0:
                        slots = insertion_slots_by_ref.get(ref_idx, [])[:ins_len]
                        for slot in slots:
                            copy_from_indices.append(source_cursor)
                            expanded_indices.append(slot)
                            source_cursor += 1
                        source_cursor += ins_len - len(slots)

                per_hap[haplotype] = {
                    "source_start_idx": int(hap_state.get("source_start_idx", gene_payload.get("ref_start_offset", 0))),
                    "copy_from_indices": copy_from_indices,
                    "expanded_indices": expanded_indices,
                    "insertion_indices": insertion_slots,
                    "deletion_indices": sorted(expanded_index_map[idx] for idx in deletion_positions),
                }

            cache[sample_id] = per_hap
            if gene_payload.get("axis_finalized"):
                self._save_sample_to_persistent_cache(gene_name, sample_id, gene_payload, per_hap)
            cache.move_to_end(sample_id)
            while len(cache) > self._sample_cache_limit:
                cache.popitem(last=False)

        self.profile_stats["entry_build_s"] += time.perf_counter() - t0
        self.profile_stats["entry_build_calls"] += 1

    def build_alignment_axis_for_gene(self, gene_name: str, sample_ids: List[str]) -> None:
        self._write_sample_set_manifest()
        gene_payload = self._load_gene_mapping(gene_name)
        if self._load_axis_from_persistent_cache(gene_name, gene_payload):
            return
        if gene_payload.get("axis_finalized"):
            return

        start_1based = int(gene_payload["start_1based"])
        full_start_1based = int(gene_payload.get("full_start_1based", start_1based))
        ref_length = int(gene_payload["ref_length"])
        ref_sequence = str(gene_payload.get("ref_sequence", ""))
        insertion_after_ref_index: Dict[int, int] = {}

        query_t0 = time.perf_counter()
        for row in _stream_gene_variants(gene_payload["vcf_path"], sample_ids, gene_payload["region"]):
            if _classify_length_behavior(str(row["ref"]), str(row["alt"])) != "length_change":
                continue
            ref = str(row["ref"])
            pos0_full = int(row["pos_1based"]) - full_start_1based
            pos0 = _resolve_ref_index(ref_sequence, pos0_full - int(gene_payload.get("ref_start_offset", 0)), ref)
            if pos0 < 0 or pos0 >= ref_length:
                continue
            for gt in row.get("gts", []):
                gt_left, gt_right = _parse_gt(str(gt))
                for token in (gt_left, gt_right):
                    allele = _allele_sequence(ref, str(row["alt"]), token)
                    ins_len = len(allele) - len(ref)
                    if ins_len > 0:
                        insertion_after_ref_index[pos0] = max(insertion_after_ref_index.get(pos0, 0), ins_len)

        self.profile_stats["bcftools_query_s"] += time.perf_counter() - query_t0
        self.profile_stats["query_calls"] += 1

        expanded_index_map: Dict[int, int] = {}
        insertion_slots_by_ref: Dict[int, List[int]] = {}
        expanded_cursor = 0
        for ref_idx in range(ref_length):
            expanded_index_map[ref_idx] = expanded_cursor
            expanded_cursor += 1
            ins_slots = insertion_after_ref_index.get(ref_idx, 0)
            if ins_slots > 0:
                insertion_slots_by_ref[ref_idx] = list(range(expanded_cursor, expanded_cursor + ins_slots))
                expanded_cursor += ins_slots

        gene_payload["expanded_length"] = expanded_cursor
        gene_payload["expanded_index_map"] = expanded_index_map
        gene_payload["insertion_slots_by_ref"] = insertion_slots_by_ref
        gene_payload["axis_finalized"] = True
        self._save_axis_to_persistent_cache(gene_name, gene_payload)

    def get_expanded_length(self, gene_name: str) -> int:
        return int(self._load_gene_mapping(gene_name)["expanded_length"])

    def get_alignment_axis(self, gene_name: str) -> Dict[str, object]:
        """Return the current expanded reference axis metadata for a gene."""
        payload = self._load_gene_mapping(gene_name)
        if not payload.get("axis_finalized"):
            self._load_axis_from_persistent_cache(gene_name, payload)
        return {
            "expanded_length": int(payload["expanded_length"]),
            "ref_length": int(payload["ref_length"]),
            "full_ref_length": int(payload.get("full_ref_length", payload["ref_length"])),
            "ref_start_offset": int(payload.get("ref_start_offset", 0)),
            "alignment_start_1based": int(payload.get("alignment_start_1based", payload.get("start_1based", 0))),
            "alignment_end_1based": int(payload.get("alignment_end_1based", payload.get("start_1based", 0) + int(payload["ref_length"]) - 1)),
            "region": str(payload.get("region", "")),
            "expanded_index_map": dict(payload.get("expanded_index_map", {})),
            "insertion_slots_by_ref": dict(payload.get("insertion_slots_by_ref", {})),
            "algorithm_version": ALIGNMENT_ALGORITHM_VERSION,
            "sample_set_key": self._sample_set_key(),
            "cache_dir": str(self._gene_cache_dir(gene_name)),
            "fingerprints": payload.get("fingerprints", {}),
        }

    def get_reference_centered_expanded_slice(self, gene_name: str, window_size: int) -> Dict[str, int]:
        """Return a fixed-size expanded-axis slice centered on the reference gene center."""
        payload = self._load_gene_mapping(gene_name)
        if not payload.get("axis_finalized"):
            self._load_axis_from_persistent_cache(gene_name, payload)
        ref_length = int(payload["ref_length"])
        expanded_length = int(payload["expanded_length"])
        expanded_index_map = payload.get("expanded_index_map", {})
        center_ref_idx = ref_length // 2
        center_expanded_idx = int(expanded_index_map.get(center_ref_idx, center_ref_idx))
        size = max(int(window_size), 1)
        if expanded_length <= size:
            return {
                "center_ref_idx": center_ref_idx,
                "center_expanded_idx": center_expanded_idx,
                "expanded_start": 0,
                "expanded_end": expanded_length,
                "expanded_length": expanded_length,
                "window_size": expanded_length,
            }
        half = size // 2
        start = center_expanded_idx - half
        end = start + size
        if start < 0:
            end -= start
            start = 0
        if end > expanded_length:
            start = max(0, start - (end - expanded_length))
            end = expanded_length
        return {
            "center_ref_idx": center_ref_idx,
            "center_expanded_idx": center_expanded_idx,
            "expanded_start": start,
            "expanded_end": end,
            "expanded_length": expanded_length,
            "window_size": end - start,
        }

    def get_haplotype_entry(self, gene_name: str, sample_id: str, haplotype: str) -> Optional[Dict[str, object]]:
        self._build_sample_entries_for_gene(gene_name, [sample_id])
        sample_entry = self._sample_entry_cache[gene_name].get(sample_id)
        if sample_entry is None:
            return None
        return sample_entry.get(haplotype)
