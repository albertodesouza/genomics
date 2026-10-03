"""Haplotype <-> reference coordinate mapping from per-sample window VCFs.

AlphaGenome predictions are made on each haplotype's own consensus sequence (``*.window.fixed.fa``,
``bcftools consensus`` output truncated/padded to the window size), so position ``i`` of a
prediction array is position ``i`` of that haplotype, not of the reference. This module rebuilds
the haplotype's indel drift from its phased window VCF -- mirroring ``bcftools consensus``'s rule
of skipping a record that overlaps an already-applied one -- so any haplotype-indexed array can be
remapped onto reference (genomic) coordinates, shared by every sample, with deleted reference
bases reported as missing.
"""
from __future__ import annotations

import gzip
import struct
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HAPLOTYPE_INDEX = {"H1": 0, "H2": 1}


@dataclass
class HaplotypeEvents:
    """Applied (non-reference) variants for one haplotype, as window-relative 0-based offsets."""

    offsets: np.ndarray  # int64, 0-based offset of the record's POS within the window
    ref_lengths: np.ndarray  # int32
    alt_lengths: np.ndarray  # int32

    @classmethod
    def empty(cls) -> "HaplotypeEvents":
        return cls(np.zeros(0, np.int64), np.zeros(0, np.int32), np.zeros(0, np.int32))

    @property
    def indel_count(self) -> int:
        return int(np.count_nonzero(self.ref_lengths != self.alt_lengths))


@dataclass
class SampleVariants:
    """All non-reference records of a sample's window VCF plus per-haplotype applied events."""

    positions: np.ndarray  # int64 genomic 1-based POS
    ids: List[str]
    refs: List[str]
    alts: List[str]  # the specific ALT allele carried (per record, first non-ref allele carried)
    genotypes: List[str]
    carried: np.ndarray  # int8 (n, 2): allele index carried by H1/H2 (0 = ref, -1 = missing)
    alt_alleles: List[List[str]]
    haplotypes: Dict[str, HaplotypeEvents] = field(default_factory=dict)

    def events(self, haplotype: str) -> HaplotypeEvents:
        return self.haplotypes.get(haplotype) or HaplotypeEvents.empty()


def _parse_allele(token: bytes) -> int:
    if token in (b".", b""):
        return -1
    try:
        return int(token)
    except ValueError:
        return -1


def _info_end(info: bytes) -> Optional[int]:
    for item in info.split(b";"):
        if item.startswith(b"END="):
            try:
                return int(item[4:])
            except ValueError:
                return None
    return None


def inflate_gzip(data: bytes) -> bytes:
    """Decompress (BGZF or plain) gzip data.

    ``gzip.decompress`` copies the remaining buffer for every member, which is quadratic for BGZF
    files (one member per 64 KiB block) and holds the GIL while doing so. BGZF blocks carry their
    own size, so each is inflated directly with a GIL-free ``zlib.decompress`` call.
    """
    view = memoryview(data)
    total = len(data)
    pos = 0
    parts: List[bytes] = []
    while pos + 18 <= total and data[pos] == 0x1F and data[pos + 1] == 0x8B:
        block_size = None
        if data[pos + 3] & 0x04:  # FEXTRA
            xlen = struct.unpack_from("<H", data, pos + 10)[0]
            i = pos + 12
            while i + 4 <= pos + 12 + xlen:
                sub_len = struct.unpack_from("<H", data, i + 2)[0]
                if data[i] == 66 and data[i + 1] == 67 and sub_len == 2:
                    block_size = struct.unpack_from("<H", data, i + 4)[0] + 1
                i += 4 + sub_len
            payload_start = pos + 12 + xlen
        if block_size is None:
            parts.append(gzip.decompress(bytes(view[pos:])))
            return b"".join(parts)
        payload_end = pos + block_size - 8
        size = struct.unpack_from("<I", data, payload_end + 4)[0]
        if payload_end > payload_start:
            parts.append(zlib.decompress(view[payload_start:payload_end], -15, max(size, 1)))
        pos += block_size
    return b"".join(parts)


def parse_window_vcf(path: Path, window_start_1based: int, sample_column: int = 0) -> SampleVariants:
    """Parse a single-sample (or ``sample_column``-th sample) phased window VCF."""
    with open(path, "rb") as f:
        raw = f.read()
    text = inflate_gzip(raw) if str(path).endswith(".gz") else raw

    positions: List[int] = []
    ids: List[str] = []
    refs: List[str] = []
    alts: List[str] = []
    genotypes: List[str] = []
    carried: List[Tuple[int, int]] = []
    alt_alleles: List[List[str]] = []
    applied: Dict[str, List[Tuple[int, int, int]]] = {"H1": [], "H2": []}
    frozen_end = {"H1": -1, "H2": -1}
    last_insertion = {"H1": False, "H2": False}
    column = 9 + sample_column

    header_at = text.find(b"\n#CHROM")
    header_end = text.find(b"\n", header_at + 1)
    single_sample = header_at >= 0 and text[header_at + 1:header_end].count(b"\t") == 9
    homref_suffix = (b"\t0|0", b"\t0/0") if single_sample else ()

    for line in text.split(b"\n"):
        # Fast path: in single-sample VCFs most records end with a homozygous-reference GT.
        if not line or line[:1] == b"#" or (homref_suffix and line.endswith(homref_suffix)):
            continue
        cols = line.split(b"\t", column + 1)
        if len(cols) <= column:
            continue
        gt_field = cols[column].split(b":", 1)[0]
        if gt_field in (b"0|0", b"0/0", b"0", b".|.", b"./.", b"."):
            continue
        sep = b"|" if b"|" in gt_field else b"/"
        tokens = gt_field.split(sep)
        a1 = _parse_allele(tokens[0])
        a2 = _parse_allele(tokens[1]) if len(tokens) > 1 else a1
        if a1 <= 0 and a2 <= 0:
            continue
        pos = int(cols[1])
        ref = cols[3].decode("ascii", "replace").upper()
        alt_list = cols[4].decode("ascii", "replace").upper().split(",")
        carried_idx = a1 if a1 > 0 else a2
        positions.append(pos)
        ids.append("" if cols[2] == b"." else cols[2].decode("ascii", "replace"))
        refs.append(ref)
        alts.append(alt_list[carried_idx - 1] if 0 < carried_idx <= len(alt_list) else ".")
        genotypes.append(gt_field.decode("ascii", "replace"))
        carried.append((a1, a2))
        alt_alleles.append(alt_list)

        offset = pos - window_start_1based
        for hap, allele in (("H1", a1), ("H2", a2)):
            if allele <= 0 or allele > len(alt_list):
                continue
            alt = alt_list[allele - 1]
            ref_len, alt_len = len(ref), len(alt)
            if alt == "<DEL>":
                # bcftools consensus applies symbolic deletions over POS..INFO/END.
                end = _info_end(cols[7])
                if end is None or end <= pos:
                    continue
                ref_len, alt_len = end - pos + 1, 1
            elif not alt or alt[0] in "<*" or "[" in alt or "]" in alt:
                continue
            # bcftools consensus skips records overlapping a previously applied variant, except an
            # indel anchored on the last base of a preceding record that was not an insertion.
            if offset < frozen_end[hap]:
                continue
            is_insertion = alt_len > ref_len
            if offset == frozen_end[hap] and (ref_len == alt_len or last_insertion[hap]):
                continue
            applied[hap].append((offset, ref_len, alt_len))
            frozen_end[hap] = offset + ref_len - 1
            last_insertion[hap] = is_insertion

    haplotypes = {}
    for hap, rows in applied.items():
        if rows:
            arr = np.asarray(rows, dtype=np.int64)
            haplotypes[hap] = HaplotypeEvents(arr[:, 0], arr[:, 1].astype(np.int32), arr[:, 2].astype(np.int32))
        else:
            haplotypes[hap] = HaplotypeEvents.empty()
    return SampleVariants(
        positions=np.asarray(positions, dtype=np.int64),
        ids=ids,
        refs=refs,
        alts=alts,
        genotypes=genotypes,
        carried=np.asarray(carried, dtype=np.int8).reshape(-1, 2),
        alt_alleles=alt_alleles,
        haplotypes=haplotypes,
    )


def parse_events_batch(jobs: List[Tuple[str, str, int]]) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    """Worker for process pools: ``[(sample, vcf_path, window_start)] -> [(sample, events, hap_ids)]``.

    ``events`` is ``(n, 3)`` int64 (offset, ref_len, alt_len); ``hap_ids`` is ``(n,)`` int8 with
    0 = H1 and 1 = H2. Failed files are returned with ``events = None``-like empty arrays of
    shape ``(0, 3)`` and ``hap_ids`` of ``-1`` length 1 so the caller can tell them apart.
    """
    out = []
    for sample, path, window_start in jobs:
        try:
            parsed = parse_window_vcf(Path(path), window_start)
        except Exception:
            out.append((sample, np.zeros((0, 3), np.int64), np.full(1, -1, np.int8)))
            continue
        events, haps = pack_events(parsed.haplotypes)
        out.append((sample, events, haps))
    return out


def pack_events(haplotypes: Dict[str, HaplotypeEvents]) -> Tuple[np.ndarray, np.ndarray]:
    """Concatenate the length-changing events of both haplotypes (SNVs do not move coordinates)."""
    blocks, ids = [], []
    for hap, idx in HAPLOTYPE_INDEX.items():
        ev = haplotypes.get(hap)
        if ev is None or not ev.offsets.size:
            continue
        keep = ev.ref_lengths != ev.alt_lengths
        if keep.any():
            blocks.append(np.stack([ev.offsets[keep], ev.ref_lengths[keep].astype(np.int64), ev.alt_lengths[keep].astype(np.int64)], axis=1))
            ids.append(np.full(int(keep.sum()), idx, dtype=np.int8))
    if not blocks:
        return np.zeros((0, 3), np.int64), np.zeros(0, np.int8)
    return np.concatenate(blocks), np.concatenate(ids)


def unpack_events(events: np.ndarray, haps: np.ndarray) -> Dict[str, HaplotypeEvents]:
    out = {}
    for hap, idx in HAPLOTYPE_INDEX.items():
        sel = events[haps == idx]
        out[hap] = HaplotypeEvents(sel[:, 0].astype(np.int64), sel[:, 1].astype(np.int32), sel[:, 2].astype(np.int32))
    return out


@dataclass
class ReferenceMap:
    """For each reference offset in ``[0, length)``: haplotype-local index, or -1 if deleted."""

    local: np.ndarray  # int64, -1 where the reference base is deleted on this haplotype
    insertion_offsets: np.ndarray  # reference offsets after which bases are inserted
    insertion_lengths: np.ndarray

    def valid(self, haplotype_length: int) -> np.ndarray:
        return (self.local >= 0) & (self.local < haplotype_length)


def reference_map(events: HaplotypeEvents, length: int) -> ReferenceMap:
    """Vectorised reference-offset -> haplotype-index map for a window of ``length`` bases."""
    change = np.zeros(length + 1, dtype=np.int64)
    deleted = np.zeros(length + 1, dtype=np.int32)
    ins_off: List[int] = []
    ins_len: List[int] = []
    for offset, ref_len, alt_len in zip(events.offsets.tolist(), events.ref_lengths.tolist(), events.alt_lengths.tolist()):
        delta = alt_len - ref_len
        if delta == 0:
            continue
        first_shifted = offset + min(ref_len, alt_len)
        if first_shifted < 0:
            change[0] += delta
        elif first_shifted <= length:
            change[first_shifted] += delta
        if delta < 0:
            lo = max(offset + alt_len, 0)
            hi = min(offset + ref_len, length)
            if lo < hi:
                deleted[lo] += 1
                deleted[hi] -= 1
        elif 0 <= offset + ref_len - 1 < length:
            ins_off.append(offset + ref_len - 1)
            ins_len.append(delta)
    local = np.arange(length, dtype=np.int64) + np.cumsum(change)[:length]
    local[np.cumsum(deleted)[:length] > 0] = -1
    return ReferenceMap(local, np.asarray(ins_off, dtype=np.int64), np.asarray(ins_len, dtype=np.int64))


def remap_to_reference(
    values: np.ndarray, ref_map: ReferenceMap, start: int, end: int
) -> np.ndarray:
    """Gather ``values`` (haplotype-indexed, first axis) at reference offsets ``[start, end)``.

    Missing positions (deleted on the haplotype, or beyond the truncated haplotype) become NaN.
    """
    local = ref_map.local[start:end]
    valid = (local >= 0) & (local < values.shape[0])
    out_shape = (end - start,) + tuple(values.shape[1:])
    out = np.full(out_shape, np.nan, dtype=np.float32)
    out[valid] = values[local[valid]]
    return out


def bases_on_reference(sequence: bytes, ref_map: ReferenceMap, start: int, end: int) -> bytes:
    """Haplotype bases at reference offsets ``[start, end)``; ``-`` where deleted/missing."""
    arr = np.frombuffer(sequence, dtype=np.uint8)
    local = ref_map.local[start:end]
    valid = (local >= 0) & (local < arr.size)
    out = np.full(end - start, ord("-"), dtype=np.uint8)
    out[valid] = arr[local[valid]]
    return out.tobytes()


def merge_insertions(ref_maps: Sequence[ReferenceMap], start: int, end: int) -> Dict[int, int]:
    """Max inserted length after each reference offset in ``[start, end)`` across haplotypes."""
    merged: Dict[int, int] = {}
    for ref_map in ref_maps:
        mask = (ref_map.insertion_offsets >= start) & (ref_map.insertion_offsets < end)
        for off, ln in zip(ref_map.insertion_offsets[mask].tolist(), ref_map.insertion_lengths[mask].tolist()):
            merged[off] = max(merged.get(off, 0), ln)
    return merged


def window_start_from_metadata(window_meta: Dict[str, object]) -> Optional[int]:
    try:
        return int(window_meta["start"])  # type: ignore[index]
    except (KeyError, TypeError, ValueError):
        return None
