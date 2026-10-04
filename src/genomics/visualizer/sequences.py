"""Haplotype sequences against the reference, in any coordinate system."""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from genomics.visualizer.datasets import Dataset
from genomics.visualizer.signals import SeriesSpec, SignalService, _clamp_range

MAX_LETTER_SPAN = 4000
MAX_ROWS = 64
GAP = ord("-")
LETTERS = ("A", "C", "G", "T", "N", "-")


def _variant_type(ref: str, alt: str) -> str:
    if alt.startswith("<"):
        return alt.strip("<>").split(":")[0] or "SV"
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    return "INS" if len(alt) > len(ref) else "DEL"


class SequenceService:
    def __init__(self, signals: SignalService):
        self.signals = signals

    def _row_reference(self, dataset: Dataset, gene: str, spec: SeriesSpec, start: int, end: int) -> Tuple[np.ndarray, List[Dict[str, object]]]:
        hap = np.frombuffer(self.signals.haplotype_sequence(dataset, spec.sample, gene, spec.haplotype), dtype=np.uint8)
        ref_map = self.signals.ref_map(dataset, spec.sample, gene, spec.haplotype)
        local = ref_map.local[start:end]
        valid = (local >= 0) & (local < hap.size)
        row = np.full(end - start, GAP, dtype=np.uint8)
        row[valid] = hap[local[valid]]
        insertions = []
        mask = (ref_map.insertion_offsets >= start) & (ref_map.insertion_offsets < end)
        for off, length in zip(ref_map.insertion_offsets[mask].tolist(), ref_map.insertion_lengths[mask].tolist()):
            anchor = ref_map.local[off]
            seq = hap[anchor + 1: anchor + 1 + length].tobytes().decode("ascii", "replace") if anchor >= 0 else ""
            insertions.append({"pos": off, "length": length, "seq": seq[:200]})
        return row, insertions

    def _row_aligned(self, dataset: Dataset, gene: str, spec: SeriesSpec, start: int, end: int) -> np.ndarray:
        hap = np.frombuffer(self.signals.haplotype_sequence(dataset, spec.sample, gene, spec.haplotype), dtype=np.uint8)
        expanded, source = self.signals.alignment.entry_arrays(dataset, gene, spec.sample, spec.haplotype)
        row = np.full(end - start, GAP, dtype=np.uint8)
        mask = (expanded >= start) & (expanded < end) & (source >= 0) & (source < hap.size)
        row[expanded[mask] - start] = hap[source[mask]]
        return row

    def _reference_row(self, dataset: Dataset, gene: str, coords: str, start: int, end: int) -> np.ndarray:
        ref = np.frombuffer(self.signals.reference_sequence(dataset, gene), dtype=np.uint8)
        if coords == "aligned":
            axis = self.signals.alignment.axis(dataset, gene)
            ref_idx = axis["ref_of_expanded"][start:end]
            row = np.full(end - start, GAP, dtype=np.uint8)
            ok = ref_idx >= 0
            row[ok] = ref[ref_idx[ok] + axis["ref_start_offset"]]
            return row
        row = np.full(end - start, ord("N"), dtype=np.uint8)
        hi = min(end, ref.size)
        if hi > start:
            row[: hi - start] = ref[start:hi]
        return row

    def _row(self, dataset: Dataset, gene: str, spec: SeriesSpec, coords: str, start: int, end: int) -> Tuple[np.ndarray, List[Dict[str, object]]]:
        """One haplotype's bases over ``[start, end)`` of ``coords`` (``-`` where it has none)."""
        if coords == "reference":
            return self._row_reference(dataset, gene, spec, start, end)
        if coords == "aligned":
            return self._row_aligned(dataset, gene, spec, start, end), []
        hap = np.frombuffer(self.signals.haplotype_sequence(dataset, spec.sample, gene, spec.haplotype), dtype=np.uint8)
        row = np.full(end - start, GAP, dtype=np.uint8)
        hi = min(end, hap.size)
        if hi > start:
            row[: hi - start] = hap[start:hi]
        return row, []

    def _expand_rows(self, rows: List[SeriesSpec]) -> List[SeriesSpec]:
        expanded: List[SeriesSpec] = []
        for spec in rows:
            haps = ["H1", "H2"] if spec.haplotype == "H1+H2" else [spec.haplotype]
            expanded.extend(SeriesSpec(spec.sample, h) for h in haps)
        return expanded[:MAX_ROWS]

    def _domain(self, dataset: Dataset, gene: str, coords: str, rows: List[SeriesSpec]) -> int:
        if coords == "aligned":
            return int(self.signals.alignment.axis(dataset, gene)["expanded_length"])
        if coords == "reference":
            return self.signals.reference_length(dataset, gene)
        return max((len(self.signals.haplotype_sequence(dataset, s.sample, gene, s.haplotype)) for s in rows), default=0) or self.signals.reference_length(dataset, gene)

    def composition(
        self,
        dataset: Dataset,
        gene: str,
        rows: List[SeriesSpec],
        coords: str,
        start: int,
        end: int,
        bins: int,
    ) -> Dict[str, object]:
        """Per-base letter frequencies (A, C, G, T, N, gap) over haplotypes, or of the reference.

        Without ``rows`` the reference sequence is counted (one letter per base). Each bin holds the
        fraction of (haplotype, base) cells with each letter, so zoomed-in bins are per-base
        frequencies across the chosen haplotypes and zoomed-out bins are base composition.
        """
        expanded = self._expand_rows(rows)
        domain = self._domain(dataset, gene, coords, expanded)
        start, end = _clamp_range(start, end, domain)
        span = end - start
        ref_row = self._reference_row(dataset, gene, coords, start, end)
        counts = np.zeros((len(LETTERS), span), dtype=np.uint16)
        used: List[str] = []
        errors: List[str] = []
        sources = []
        for spec in expanded:
            try:
                sources.append(self._row(dataset, gene, spec, coords, start, end)[0])
                used.append(f"{spec.sample}:{spec.haplotype}")
            except Exception as exc:
                errors.append(f"{spec.sample} {spec.haplotype}: {exc}")
        if not sources:
            sources = [ref_row]
        for row in sources:
            upper = row & 0xDF  # ASCII upper case (soft-masked bases count as their letter)
            known = np.zeros(span, dtype=bool)
            for i, letter in enumerate(LETTERS[:4]):
                hit = upper == ord(letter)
                counts[i] += hit
                known |= hit
            gap = row == GAP
            counts[5] += gap
            counts[4] += ~(known | gap)
        edges = np.unique(np.linspace(0, span, min(max(1, int(bins)), span) + 1).astype(np.int64))
        totals = np.add.reduceat(counts.astype(np.float32), edges[:-1], axis=1)
        freq = totals / (np.diff(edges).astype(np.float32) * len(sources))
        payload: Dict[str, object] = {
            "gene": gene,
            "coords": coords,
            "start": start,
            "end": end,
            "domain": domain,
            "edges": (edges + start).astype(np.int64),
            "letters": "".join(LETTERS),
            "freq": freq.astype(np.float32),
            "source": "haplotypes" if used else "reference",
            "haplotypes": used,
            "errors": errors[:5],
        }
        if span <= MAX_LETTER_SPAN:
            payload["reference"] = ref_row.tobytes().decode("ascii", "replace")
        return payload

    def window(
        self,
        dataset: Dataset,
        gene: str,
        rows: List[SeriesSpec],
        coords: str,
        start: int,
        end: int,
        bins: int,
    ) -> Dict[str, object]:
        expanded_rows = self._expand_rows(rows)
        domain = self._domain(dataset, gene, coords, expanded_rows)
        start, end = _clamp_range(start, end, domain)
        span = end - start
        ref_row = self._reference_row(dataset, gene, coords, start, end)
        letters = span <= MAX_LETTER_SPAN
        out_rows = []
        for spec in expanded_rows:
            item: Dict[str, object] = {"sample": spec.sample, "haplotype": spec.haplotype}
            try:
                row, insertions = self._row(dataset, gene, spec, coords, start, end)
                diff = (row != ref_row) & (row != GAP) & (ref_row != GAP) & (row != ord("N"))
                deleted = (row == GAP) & (ref_row != GAP)
                item["mismatches"] = int(diff.sum())
                item["deletions"] = int(deleted.sum())
                item["insertions"] = insertions
                if letters:
                    item["bases"] = row.tobytes().decode("ascii", "replace")
                else:
                    item["density"] = _bin_flags(diff, deleted, insertions, start, bins)
            except Exception as exc:
                item["error"] = str(exc)
            out_rows.append(item)
        payload: Dict[str, object] = {
            "gene": gene,
            "coords": coords,
            "start": start,
            "end": end,
            "domain": domain,
            "mode": "letters" if letters else "density",
            "rows": out_rows,
        }
        if letters:
            payload["reference"] = ref_row.tobytes().decode("ascii", "replace")
        else:
            gc = np.isin(ref_row, np.frombuffer(b"GC", dtype=np.uint8)).astype(np.float32)
            edges = np.unique(np.linspace(0, span, min(bins, span) + 1).astype(np.int64))
            payload["edges"] = (edges + start).astype(np.int64)
            payload["gc"] = (np.add.reduceat(gc, edges[:-1]) / np.diff(edges)).astype(np.float32)
        if coords != "haplotype":
            payload["variants"] = self.variants(dataset, gene, [s.sample for s in rows], coords, start, end)
        return payload

    def variants(self, dataset: Dataset, gene: str, samples: List[str], coords: str, start: int, end: int, limit: int = 5000) -> List[Dict[str, object]]:
        window = dataset.window(gene)
        if window.start is None:
            return []
        to_expanded: Optional[np.ndarray] = None
        if coords == "aligned":
            to_expanded = self.signals.alignment.reference_offset_to_expanded(dataset, gene)
        records: Dict[Tuple[int, str, str], Dict[str, object]] = {}
        seen_samples = []
        for sample in dict.fromkeys(samples):
            try:
                parsed = self.signals.variants(dataset, sample, gene)
            except Exception:
                continue
            seen_samples.append(sample)
            offsets = parsed.positions - window.start
            if to_expanded is not None:
                in_axis = (offsets >= 0) & (offsets < to_expanded.size)
                coord = np.full(offsets.size, -1, dtype=np.int64)
                coord[in_axis] = to_expanded[offsets[in_axis]]
            else:
                coord = offsets
            for i in np.nonzero((coord >= start) & (coord < end))[0].tolist():
                ref, alt = parsed.refs[i], parsed.alts[i]
                key = (int(parsed.positions[i]), ref, alt)
                record = records.get(key)
                if record is None:
                    record = records[key] = {
                        "pos": int(coord[i]),
                        "genomic": int(parsed.positions[i]),
                        "id": parsed.ids[i],
                        "ref": ref[:50],
                        "alt": alt[:50],
                        "type": _variant_type(ref, alt),
                        "carriers": [],
                    }
                record["carriers"].append({"sample": sample, "gt": parsed.genotypes[i]})
            if len(records) > limit:
                break
        return sorted(records.values(), key=lambda r: (r["pos"], r["genomic"]))[:limit]


def _bin_flags(diff: np.ndarray, deleted: np.ndarray, insertions: List[Dict[str, object]], start: int, bins: int) -> Dict[str, np.ndarray]:
    span = diff.size
    edges = np.unique(np.linspace(0, span, min(bins, span) + 1).astype(np.int64))
    widths = np.diff(edges).astype(np.float32)
    mism = np.add.reduceat(diff.astype(np.float32), edges[:-1]) / widths
    dels = np.add.reduceat(deleted.astype(np.float32), edges[:-1]) / widths
    ins = np.zeros(edges.size - 1, dtype=np.float32)
    if insertions:
        idx = np.searchsorted(edges, np.asarray([int(i["pos"]) - start for i in insertions]), side="right") - 1
        np.add.at(ins, np.clip(idx, 0, ins.size - 1), 1.0)
    return {"mismatch": mism.astype(np.float32), "deletion": dels.astype(np.float32), "insertion": ins}
