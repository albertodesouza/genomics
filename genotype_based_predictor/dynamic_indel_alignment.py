"""Dynamic INDEL-aware expanded-axis construction from window VCFs."""

from __future__ import annotations

import json
import subprocess
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


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


def _query_gene_variants(vcf_path: str, sample_ids: List[str], region: str) -> List[Dict[str, object]]:
    if not sample_ids:
        return []

    fmt = "%POS\t%REF\t%ALT[\t%GT]\n"
    cmd = [
        "bcftools",
        "query",
        "-s",
        ",".join(sample_ids),
        "-r",
        region,
        "-f",
        fmt,
        vcf_path,
    ]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    rows: List[Dict[str, object]] = []
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        pos, ref, alt = fields[:3]
        rows.append({
            "pos_1based": int(pos),
            "ref": ref,
            "alt": alt,
            "gts": fields[3:],
        })
    return rows


class DynamicIndelAligner:
    """Builds expanded-axis mappings dynamically for the active experiment view."""

    def __init__(self, dataset_dir: Path, selected_sample_ids: Optional[Iterable[str]] = None):
        self.dataset_dir = Path(dataset_dir)
        self.selected_sample_ids = set(selected_sample_ids or [])
        self._gene_cache: Dict[str, Dict[str, object]] = {}
        self._sample_entry_cache: Dict[str, OrderedDict[Tuple[str, str], Dict[str, object]]] = {}
        self._sample_cache_limit = 32

    def _load_window_metadata(self, gene_name: str) -> Dict[str, object]:
        meta_path = self.dataset_dir / "references" / "windows" / gene_name / "window_metadata.json"
        with open(meta_path) as f:
            return json.load(f)

    def _build_region(self, window_meta: Dict[str, object]) -> str:
        start_1based = int(window_meta["start"]) + 1
        end_1based = int(window_meta["end"]) + 1
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

        start_1based = int(window_meta["start"]) + 1
        end_1based = int(window_meta["end"]) + 1
        ref_length = end_1based - start_1based + 1

        payload = {
            "vcf_path": vcf_path,
            "region": self._build_region(window_meta),
            "start_1based": start_1based,
            "ref_length": ref_length,
            "expanded_length": ref_length,
        }
        self._gene_cache[gene_name] = payload
        self._sample_entry_cache[gene_name] = OrderedDict()
        return payload

    def _build_sample_entries_for_gene(self, gene_name: str, sample_ids: List[str]) -> None:
        gene_payload = self._load_gene_mapping(gene_name)
        cache = self._sample_entry_cache[gene_name]
        missing_sample_ids = [sample_id for sample_id in sample_ids if sample_id not in cache]
        if not missing_sample_ids:
            return

        variants = _query_gene_variants(gene_payload["vcf_path"], missing_sample_ids, gene_payload["region"])
        insertion_after_ref_index: Dict[int, int] = {}
        start_1based = int(gene_payload["start_1based"])
        ref_length = int(gene_payload["ref_length"])
        sample_entries: Dict[str, Dict[str, Dict[str, object]]] = {}

        for sample_idx, sample_id in enumerate(missing_sample_ids):
            hap_entries: Dict[str, Dict[str, object]] = {}
            for haplotype, gt_index in (("H1", 0), ("H2", 1)):
                deletion_positions = set()
                insertion_specs: List[Tuple[int, int]] = []
                for row in variants:
                    if _classify_length_behavior(str(row["ref"]), str(row["alt"])) != "length_change":
                        continue
                    row_gts = row.get("gts", [])
                    if sample_idx >= len(row_gts):
                        continue
                    gt_left, gt_right = _parse_gt(str(row_gts[sample_idx]))
                    token = gt_left if gt_index == 0 else gt_right
                    ref = str(row["ref"])
                    allele = _allele_sequence(ref, str(row["alt"]), token)
                    pos0 = int(row["pos_1based"]) - start_1based
                    if pos0 < 0 or pos0 >= ref_length:
                        continue
                    if len(allele) == len(ref):
                        continue
                    if len(allele) < len(ref):
                        for offset in range(1, len(ref)):
                            ref_idx = pos0 + offset
                            if 0 <= ref_idx < ref_length:
                                deletion_positions.add(ref_idx)
                    else:
                        ins_len = len(allele) - len(ref)
                        if ins_len > 0:
                            insertion_specs.append((pos0, ins_len))
                            insertion_after_ref_index[pos0] = max(insertion_after_ref_index.get(pos0, 0), ins_len)

                hap_entries[haplotype] = {
                    "deletion_positions": deletion_positions,
                    "insertion_specs": insertion_specs,
                }
            sample_entries[sample_id] = hap_entries

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

                per_hap[haplotype] = {
                    "copy_from_indices": copy_from_indices,
                    "expanded_indices": expanded_indices,
                    "insertion_indices": insertion_slots,
                    "deletion_indices": sorted(expanded_index_map[idx] for idx in deletion_positions),
                }

            cache[sample_id] = per_hap
            cache.move_to_end(sample_id)
            while len(cache) > self._sample_cache_limit:
                cache.popitem(last=False)

    def get_expanded_length(self, gene_name: str) -> int:
        return int(self._load_gene_mapping(gene_name)["expanded_length"])

    def get_haplotype_entry(self, gene_name: str, sample_id: str, haplotype: str) -> Optional[Dict[str, object]]:
        self._build_sample_entries_for_gene(gene_name, [sample_id])
        sample_entry = self._sample_entry_cache[gene_name].get(sample_id)
        if sample_entry is None:
            return None
        return sample_entry.get(haplotype)
