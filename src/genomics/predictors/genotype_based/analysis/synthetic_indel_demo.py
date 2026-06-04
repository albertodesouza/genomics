from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SyntheticVariant:
    pos_1based: int
    ref: str
    alt: str
    genotypes: Dict[str, str]


def classify_length_behavior(ref: str, alt: str) -> str:
    # Symbolic structural alleles like <DEL> are intentionally excluded from
    # this didactic demo because their exact effect depends on extra VCF
    # fields such as END/SVLEN and on an explicit policy about whether an
    # anchor base is preserved. To keep the demo human-verifiable, we only use
    # unambiguous explicit alleles here.
    if alt in {"<DEL>", "<del>"}:
        return "special_unhandled"
    if alt.startswith("<") and alt.endswith(">"):
        return "special_unhandled"
    if len(ref) == len(alt):
        return "same_length"
    return "length_change"


def parse_gt(gt: str) -> Tuple[str, str]:
    text = str(gt).replace("/", "|")
    if "|" not in text:
        return "0", "0"
    left, right = text.split("|", 1)
    return left, right


def allele_sequence(ref: str, alt: str, token: str) -> str:
    if token == "0":
        return ref
    if token == "1":
        return alt
    return ref


def build_demo_definition() -> Dict[str, object]:
    reference = "ACGTACGTAA"
    sample_ids = ["IND1", "IND2"]
    haplotypes = ["H1", "H2"]
    variants = [
        SyntheticVariant(pos_1based=2, ref="C", alt="T", genotypes={"IND1": "0|1", "IND2": "1|1"}),
        SyntheticVariant(pos_1based=4, ref="TA", alt="T", genotypes={"IND1": "1|1", "IND2": "0|0"}),
        SyntheticVariant(pos_1based=6, ref="C", alt="CGG", genotypes={"IND1": "0|1", "IND2": "1|0"}),
        SyntheticVariant(pos_1based=8, ref="T", alt="G", genotypes={"IND1": "1|0", "IND2": "0|1"}),
    ]
    track = np.arange(1, len(reference) + 1, dtype=np.float32)
    return {
        "reference": reference,
        "sample_ids": sample_ids,
        "haplotypes": haplotypes,
        "variants": variants,
        "track": track,
    }


def build_haplotype_sequence(reference: str, variants: Sequence[SyntheticVariant], sample_id: str, haplotype: str) -> str:
    chars: List[str] = []
    ref_cursor = 0
    hap_index = 0 if haplotype == "H1" else 1
    sorted_variants = sorted(variants, key=lambda item: item.pos_1based)
    for variant in sorted_variants:
        pos0 = variant.pos_1based - 1
        chars.append(reference[ref_cursor:pos0])
        gt = variant.genotypes[sample_id]
        token = parse_gt(gt)[hap_index]
        chars.append(allele_sequence(variant.ref, variant.alt, token))
        ref_cursor = pos0 + len(variant.ref)
    chars.append(reference[ref_cursor:])
    return "".join(chars)


def build_haplotype_track(reference_track: np.ndarray, variants: Sequence[SyntheticVariant], sample_id: str, haplotype: str) -> np.ndarray:
    values: List[float] = []
    ref_cursor = 0
    hap_index = 0 if haplotype == "H1" else 1
    sorted_variants = sorted(variants, key=lambda item: item.pos_1based)
    for variant in sorted_variants:
        pos0 = variant.pos_1based - 1
        values.extend(reference_track[ref_cursor:pos0].tolist())
        gt = variant.genotypes[sample_id]
        token = parse_gt(gt)[hap_index]
        allele = allele_sequence(variant.ref, variant.alt, token)
        if len(allele) == len(variant.ref):
            values.extend(reference_track[pos0:pos0 + len(variant.ref)].tolist())
        elif len(allele) < len(variant.ref):
            values.append(float(reference_track[pos0]))
        else:
            values.append(float(reference_track[pos0]))
            inserted_count = len(allele) - len(variant.ref)
            values.extend([float(reference_track[pos0])] * inserted_count)
        ref_cursor = pos0 + len(variant.ref)
    values.extend(reference_track[ref_cursor:].tolist())
    return np.asarray(values, dtype=np.float32)


def build_demo_alignment(reference: str, variants: Sequence[SyntheticVariant], sample_ids: Sequence[str], haplotypes: Sequence[str]) -> Dict[str, object]:
    ref_length = len(reference)
    insertion_after_ref_index: Dict[int, int] = {}
    sample_entries: Dict[str, Dict[str, object]] = {}

    for sample_id in sample_ids:
        hap_entries: Dict[str, object] = {}
        for haplotype, hap_idx in (("H1", 0), ("H2", 1)):
            deletion_positions = set()
            insertion_specs: List[Tuple[int, int]] = []
            for variant in variants:
                if classify_length_behavior(variant.ref, variant.alt) != "length_change":
                    continue
                token = parse_gt(variant.genotypes[sample_id])[hap_idx]
                allele = allele_sequence(variant.ref, variant.alt, token)
                pos0 = variant.pos_1based - 1
                if len(allele) < len(variant.ref):
                    for offset in range(1, len(variant.ref)):
                        deletion_positions.add(pos0 + offset)
                elif len(allele) > len(variant.ref):
                    ins_len = len(allele) - len(variant.ref)
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

    expanded_length = expanded_cursor
    payload: Dict[str, object] = {"expanded_length": expanded_length, "samples": {}}
    for sample_id, hap_entries in sample_entries.items():
        sample_payload: Dict[str, object] = {}
        for haplotype in haplotypes:
            hap_state = hap_entries[haplotype]
            deletion_positions = hap_state["deletion_positions"]
            insertion_specs = hap_state["insertion_specs"]
            insertion_slots: List[int] = []
            insertion_len_by_ref = {ref_idx: ins_len for ref_idx, ins_len in insertion_specs}
            for ref_idx, ins_len in insertion_specs:
                insertion_slots.extend(insertion_slots_by_ref.get(ref_idx, [])[:ins_len])

            copy_from_indices: List[int] = []
            expanded_indices: List[int] = []
            source_cursor = 0
            for ref_idx in range(ref_length):
                if ref_idx not in deletion_positions:
                    copy_from_indices.append(source_cursor)
                    expanded_indices.append(expanded_index_map[ref_idx])
                    source_cursor += 1
                ins_len = insertion_len_by_ref.get(ref_idx, 0)
                if ins_len > 0:
                    for slot in insertion_slots_by_ref.get(ref_idx, [])[:ins_len]:
                        copy_from_indices.append(source_cursor)
                        expanded_indices.append(slot)
                        source_cursor += 1

            sample_payload[haplotype] = {
                "copy_from_indices": copy_from_indices,
                "expanded_indices": expanded_indices,
                "insertion_indices": insertion_slots,
                "deletion_indices": sorted(expanded_index_map[idx] for idx in deletion_positions),
            }
        payload["samples"][sample_id] = sample_payload
    return payload


def build_demo_tensor(track: np.ndarray, entry: Dict[str, object], expanded_length: int) -> np.ndarray:
    values = np.zeros(expanded_length, dtype=np.float32)
    ins_mask = np.zeros(expanded_length, dtype=np.float32)
    del_mask = np.zeros(expanded_length, dtype=np.float32)
    valid_mask = np.zeros(expanded_length, dtype=np.float32)

    for source_idx, target_idx in zip(entry["copy_from_indices"], entry["expanded_indices"]):
        values[int(target_idx)] = float(track[int(source_idx)])
        valid_mask[int(target_idx)] = 1.0
    for idx in entry["insertion_indices"]:
        ins_mask[int(idx)] = 1.0
    for idx in entry["deletion_indices"]:
        del_mask[int(idx)] = 1.0
    return np.stack([values, ins_mask, del_mask, valid_mask], axis=0)


def build_demo_global_tensor() -> Dict[str, object]:
    demo = build_demo_definition()
    reference = demo["reference"]
    sample_ids = demo["sample_ids"]
    haplotypes = demo["haplotypes"]
    variants = demo["variants"]
    track = demo["track"]
    alignment = build_demo_alignment(reference, variants, sample_ids, haplotypes)
    expanded_length = alignment["expanded_length"]

    tensor = np.zeros((len(sample_ids), len(haplotypes), 4, expanded_length), dtype=np.float32)
    hap_sequences: Dict[Tuple[str, str], str] = {}
    for sample_idx, sample_id in enumerate(sample_ids):
        for hap_idx, haplotype in enumerate(haplotypes):
            entry = alignment["samples"][sample_id][haplotype]
            hap_track = build_haplotype_track(track, variants, sample_id, haplotype)
            tensor[sample_idx, hap_idx] = build_demo_tensor(hap_track, entry, expanded_length)
            hap_sequences[(sample_id, haplotype)] = build_haplotype_sequence(reference, variants, sample_id, haplotype)

    variant_rows = []
    for variant in variants:
        variant_rows.append({
            "pos_1based": variant.pos_1based,
            "ref": variant.ref,
            "alt": variant.alt,
            "vcf_type": classify_length_behavior(variant.ref, variant.alt),
            **{sample_id: variant.genotypes[sample_id] for sample_id in sample_ids},
        })

    return {
        "reference": reference,
        "track": track,
        "sample_ids": sample_ids,
        "haplotypes": haplotypes,
        "variants": variants,
        "variants_df": pd.DataFrame(variant_rows),
        "alignment": alignment,
        "tensor": tensor,
        "hap_sequences": hap_sequences,
    }
