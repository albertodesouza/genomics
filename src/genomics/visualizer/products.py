"""Gene products of one individual: each haplotype's mature transcripts and proteins.

For every annotated transcript of a gene (GENCODE, from the dataset's GTF cache) the exons are
mapped onto each haplotype through its indel map, spliced, and, for coding transcripts, translated
from the annotated start codon to the first in-frame stop. Phasing is therefore respected: two
variants in one codon on the same haplotype give one combined amino-acid change. Each product is
compared with the reference's (an HGVS-style ``p.`` description) and checked for nonsense-mediated
decay with the 50-nt rule.

AlphaGenome's splice-junction predictions (``splice_junctions.npz`` of the haplotype and of the
reference window) add the splicing evidence: per tissue track, the usage of every junction, PSI5
(share of its donor's reads) and PSI3 (share of its acceptor's reads), on the reference and on the
haplotype; a junction's usage is the smaller of the two. Junctions that are well used, or whose usage changes, and that are not introns of the
base transcript (MANE Select, else Ensembl canonical) become candidate isoforms: the base
transcript with that intron spliced in (exon skipping, cryptic or alternative donor/acceptor), and
an annotated intron whose junction collapses becomes an intron-retention candidate. Candidates
are translated the same way.

All coordinates are reference-window offsets, 0-based and half-open (as in
:mod:`genomics.visualizer.annotations`); junction offsets on a haplotype are mapped back through
its indel map.
"""
from __future__ import annotations

import json
import re
import threading
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from genomics.visualizer.coords import HaplotypeEvents, ReferenceMap, reference_map
from genomics.visualizer.datasets import Dataset

HAPLOTYPES = ("H1", "H2")
REFERENCE = "ref"

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {a + b + c: _AMINO[16 * i + 4 * j + k] for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa",
}
_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

NMD_DISTANCE = 50  # nt between the stop codon and the last exon-exon junction (the 50-nt rule)
NMD_START_PROXIMAL = 150  # PTCs this close to the start codon tend to escape NMD (re-initiation)
NMD_LONG_EXON = 407  # PTCs in exons longer than this tend to escape NMD (Lindeboom et al. 2016)
SPLICE_REGION_EXON = 3
SPLICE_REGION_INTRON = 8
FLANK = 1000  # variants this far up/downstream of a transcript are listed
EXPRESSED_JUNCTION = 0.1  # a track "expresses" the gene when its strongest annotated junction reaches this
SITE_FLOOR = 0.1  # PSI denominators are at least this share of the gene's strongest junction (damps noise at weak sites)
JUNCTION_FLOOR = 0.05  # candidates need this share of the gene's strongest junction somewhere
CANDIDATE_USAGE = 0.2  # an alternative junction this well used is listed even without a change
CANDIDATE_DELTA = 0.1
RETENTION_RATIO = 0.5  # an intron whose junction keeps this little of its signal (relative to the gene's other introns)
RETENTION_MIN_SIGNAL = 0.25  # ... and carries at least this share of the gene's strongest junction on the reference
MAX_CANDIDATES = 12


# ---------------------------------------------------------------------------------- sequences
def reverse_complement(seq: str) -> str:
    return seq.encode("ascii").translate(_COMPLEMENT)[::-1].decode("ascii")


def translate(cds: str) -> Tuple[str, bool]:
    """Amino acids up to (not including) the first stop codon, and whether a stop was reached."""
    out = []
    for i in range(0, len(cds) - 2, 3):
        aa = CODONS.get(cds[i:i + 3].upper(), "X")
        if aa == "*":
            return "".join(out), True
        out.append(aa)
    return "".join(out), False


def aa3(residue: str) -> str:
    return THREE_LETTER.get(residue, "Xaa")


def introns_of(exons: Sequence[Sequence[int]]) -> List[Tuple[int, int]]:
    ordered = sorted((int(s), int(e)) for s, e in exons)
    return [(a[1], b[0]) for a, b in zip(ordered[:-1], ordered[1:]) if b[0] > a[1]]


def exons_from_introns(start: int, end: int, introns: Iterable[Tuple[int, int]]) -> List[List[int]]:
    """The exons of a transcript spanning ``[start, end)`` with the given introns spliced out."""
    exons: List[List[int]] = []
    cursor = start
    for s, e in sorted(introns):
        if s > cursor:
            exons.append([cursor, s])
        cursor = max(cursor, e)
    if end > cursor:
        exons.append([cursor, end])
    return exons


def splice_in(exons: Sequence[Sequence[int]], junction: Tuple[int, int]) -> List[List[int]]:
    """``exons`` with ``junction`` (an intron) spliced: introns it overlaps are replaced by it.

    Overlapped introns are dropped, so an end that falls inside one extends the neighbouring exon
    into it (an alternative site in the intron) and an end inside an exon shortens that exon.
    Exons lying entirely inside the junction are skipped.
    """
    js, je = junction
    start = min(int(s) for s, _ in exons)
    end = max(int(e) for _, e in exons)
    kept = [(s, e) for s, e in introns_of(exons) if e <= js or s >= je]
    if js <= start or je >= end:
        return [list(x) for x in sorted(exons)]
    return exons_from_introns(start, end, kept + [(js, je)])


def retain_intron(exons: Sequence[Sequence[int]], intron: Tuple[int, int]) -> List[List[int]]:
    start = min(int(s) for s, _ in exons)
    end = max(int(e) for _, e in exons)
    return exons_from_introns(start, end, [i for i in introns_of(exons) if i != tuple(intron)])


# ---------------------------------------------------------------------------------- haplotypes
@dataclass
class HaplotypeView:
    """A haplotype's bases plus its reference-offset <-> haplotype-index maps."""

    name: str
    seq: np.ndarray  # uint8 bases
    local: np.ndarray  # reference offset -> haplotype index (-1: deleted)
    inverse: np.ndarray  # haplotype index -> reference offset (-1: inserted)

    @classmethod
    def reference(cls, seq: bytes) -> "HaplotypeView":
        arr = np.frombuffer(seq, dtype=np.uint8)
        idx = np.arange(arr.size, dtype=np.int64)
        return cls(REFERENCE, arr, idx, idx)

    @classmethod
    def from_map(cls, name: str, seq: bytes, ref_map: ReferenceMap) -> "HaplotypeView":
        arr = np.frombuffer(seq, dtype=np.uint8)
        local = np.where((ref_map.local >= 0) & (ref_map.local < arr.size), ref_map.local, -1)
        inverse = np.full(arr.size, -1, dtype=np.int64)
        ok = local >= 0
        inverse[local[ok]] = np.nonzero(ok)[0]
        return cls(name, arr, local, inverse)

    def to_local(self, offset: int) -> int:
        return int(self.local[offset]) if 0 <= offset < self.local.size else -1

    def spliced(self, exons: Sequence[Sequence[int]], strand: str) -> Tuple[str, np.ndarray, List[int]]:
        """Mature transcript (5'->3'), the haplotype index of each of its bases, and exon lengths.

        Bases inserted inside an exon are included; insertions at an exon edge are taken as intronic.
        """
        pieces: List[np.ndarray] = []
        lengths: List[int] = []
        for s, e in sorted((int(a), int(b)) for a, b in exons):
            lo = self.local[max(s, 0):min(e, self.local.size)]
            lo = lo[lo >= 0]
            if lo.size == 0:
                continue
            pieces.append(np.arange(lo[0], lo[-1] + 1, dtype=np.int64))
            lengths.append(int(lo[-1] + 1 - lo[0]))
        idx = np.concatenate(pieces) if pieces else np.zeros(0, dtype=np.int64)
        seq = self.seq[idx].tobytes().decode("ascii", "replace")
        if strand == "-":
            return reverse_complement(seq), idx[::-1].copy(), lengths[::-1]
        return seq, idx, lengths


# ---------------------------------------------------------------------------------- products
def _first_atg(mrna: str, frm: int) -> int:
    return mrna.find("ATG", max(frm, 0))


def build_product(view: HaplotypeView, exons: Sequence[Sequence[int]], strand: str, cds_start: Optional[int], ref_stop: Optional[int] = None,
                  incomplete_start: bool = False, with_sequence: bool = False, phase: int = 0) -> Dict[str, Any]:
    """Splice ``exons`` on ``view`` and, if ``cds_start`` (reference offset of the A of the start
    codon) is given, translate from it. ``ref_stop`` (reference offset of the stop codon's first
    base on the reference) lets the reading frame at the original stop be compared (frameshift)."""
    mrna, idx, lengths = view.spliced(exons, strand)
    out: Dict[str, Any] = {"mrna_length": len(mrna), "exon_count": len(lengths)}
    if with_sequence:
        out["mrna"] = mrna
    if cds_start is None:
        return out
    notes: List[str] = []
    hap_start = view.to_local(cds_start)
    hits = np.nonzero(idx == hap_start)[0] if hap_start >= 0 else np.zeros(0, dtype=np.int64)
    start_lost = False
    if hits.size == 0:
        start_lost = True
        k = _first_atg(mrna, 0)
        notes.append("start codon not in this transcript" if hap_start >= 0 else "start codon deleted")
    else:
        k = int(hits[0])
        if incomplete_start:
            k += phase  # the annotated CDS starts mid-codon
        elif mrna[k:k + 3] != "ATG":
            start_lost = True
            notes.append(f"start codon {mrna[k:k + 3]}")
            k = _first_atg(mrna, k)
    if k < 0:
        out.update(coding=True, start_lost=True, protein="", stop_found=False, notes=notes + ["no ATG downstream"], cds_start=None)
        return out
    if start_lost:
        notes.append(f"translation assumed from the next AUG at transcript position {k + 1}")
    protein, stop_found = translate(mrna[k:])
    stop_at = k + 3 * len(protein)  # first base of the stop codon (or end of the ORF)
    out.update(
        coding=True,
        start_lost=start_lost,
        cds_start=k,
        cds_end=stop_at + 3 if stop_found else len(mrna),
        protein=protein,
        stop_found=stop_found,
        utr5=k,
        utr3=max(len(mrna) - (stop_at + 3), 0) if stop_found else 0,
        notes=notes,
    )
    if not stop_found:
        notes.append("no in-frame stop before the transcript ends (non-stop)")
    out["nmd"] = nmd_status(lengths, k, stop_at, stop_found)
    if ref_stop is not None:
        hap_stop = view.to_local(ref_stop)
        pos = np.nonzero(idx == hap_stop)[0] if hap_stop >= 0 else np.zeros(0, dtype=np.int64)
        out["frame_at_ref_stop"] = int((int(pos[0]) - k) % 3) if pos.size else None
    return out


def nmd_status(exon_lengths: Sequence[int], cds_start: int, stop_at: int, stop_found: bool) -> Dict[str, Any]:
    """The 50-nt rule plus the two common escape rules (start-proximal, long exon)."""
    if not stop_found:
        return {"predicted": False, "reason": "no stop codon"}
    if len(exon_lengths) < 2:
        return {"predicted": False, "reason": "single exon"}
    bounds = np.cumsum(exon_lengths)
    last_junction = int(bounds[-2])
    stop_end = stop_at + 3
    distance = last_junction - stop_end
    if distance <= NMD_DISTANCE:
        return {"predicted": False, "reason": "stop in the last exon or within 50 nt of the last junction", "distance": distance}
    exon = int(np.searchsorted(bounds, stop_at, side="right"))
    if stop_at - cds_start < NMD_START_PROXIMAL:
        return {"predicted": False, "reason": f"start-proximal stop (<{NMD_START_PROXIMAL} nt), may escape NMD", "distance": distance, "escape": True}
    if exon < len(exon_lengths) and exon_lengths[exon] > NMD_LONG_EXON:
        return {"predicted": False, "reason": f"stop in a long exon (>{NMD_LONG_EXON} nt), may escape NMD", "distance": distance, "escape": True}
    return {"predicted": True, "reason": f"stop {distance} nt upstream of the last exon junction", "distance": distance}


def _cds(product: Dict[str, Any]) -> Optional[str]:
    mrna, start, end = product.get("mrna"), product.get("cds_start"), product.get("cds_end")
    return mrna[start:end] if mrna is not None and start is not None else None


def describe_change(ref: Dict[str, Any], alt: Dict[str, Any]) -> Dict[str, Any]:
    """HGVS-style protein change of ``alt`` against ``ref`` (both from :func:`build_product`)."""
    if not ref.get("coding") or not alt.get("coding"):
        if ref.get("mrna_length") == alt.get("mrna_length") and ref.get("mrna") == alt.get("mrna"):
            return {"class": "no_change", "hgvs": ""}
        return {"class": "noncoding_change", "hgvs": "", "length_change": (alt.get("mrna_length") or 0) - (ref.get("mrna_length") or 0)}
    r, a = ref.get("protein") or "", alt.get("protein") or ""
    r_stop, a_stop = bool(ref.get("stop_found")), bool(alt.get("stop_found"))
    if alt.get("start_lost") and not ref.get("start_lost"):
        return {"class": "start_lost", "hgvs": "p.Met1?", "length": len(a)}
    if r == a and r_stop == a_stop:
        if _cds(ref) != _cds(alt):
            return {"class": "synonymous", "hgvs": "p.(=)"}
        if ref.get("mrna") != alt.get("mrna"):
            return {"class": "utr_change", "hgvs": "p.(=)", "length_change": (alt.get("mrna_length") or 0) - (ref.get("mrna_length") or 0)}
        return {"class": "no_change", "hgvs": "p.(=)"}
    i = 0
    n = min(len(r), len(a))
    while i < n and r[i] == a[i]:
        i += 1
    frame = alt.get("frame_at_ref_stop")
    if frame not in (None, 0) and i < len(r):
        new = aa3(a[i]) if i < len(a) else "Ter"
        if i == len(a) and a_stop:
            return {"class": "stop_gained", "hgvs": f"p.{aa3(r[i])}{i + 1}Ter", "position": i + 1, "frameshift": True}
        tail = f"Ter{len(a) - i + 1}" if a_stop else "Ter?"
        return {"class": "frameshift", "hgvs": f"p.{aa3(r[i])}{i + 1}{new}fs{tail}", "position": i + 1}
    if i == len(a) and a_stop and i < len(r):
        return {"class": "stop_gained", "hgvs": f"p.{aa3(r[i])}{i + 1}Ter", "position": i + 1}
    if i == len(r) and len(a) > len(r):
        tail = f"Ter{len(a) - len(r) + 1}" if a_stop else "Ter?"
        return {"class": "stop_lost", "hgvs": f"p.Ter{i + 1}{aa3(a[i])}ext{tail}", "position": i + 1}
    if len(r) == len(a) and r_stop == a_stop:
        subs = [f"{aa3(x)}{p + 1}{aa3(y)}" for p, (x, y) in enumerate(zip(r, a)) if x != y]
        shown = ";".join(subs[:6]) + (f";(+{len(subs) - 6})" if len(subs) > 6 else "")
        return {"class": "missense", "hgvs": f"p.[{shown}]" if len(subs) > 1 else f"p.{shown}", "count": len(subs), "position": i + 1}
    j = 0
    while j < min(len(r), len(a)) - i and r[len(r) - 1 - j] == a[len(a) - 1 - j]:
        j += 1
    rs, as_ = r[i:len(r) - j], a[i:len(a) - j]
    first = f"{aa3(r[i])}{i + 1}" if i < len(r) else f"Ter{i + 1}"
    last = f"_{aa3(r[i + len(rs) - 1])}{i + len(rs)}" if len(rs) > 1 else ""
    if not rs:
        flank = f"{aa3(r[i - 1])}{i}_{aa3(r[i])}{i + 1}" if 0 < i < len(r) else first
        hgvs = f"p.{flank}ins{''.join(aa3(x) for x in as_[:10])}{'…' if len(as_) > 10 else ''}"
    elif not as_:
        hgvs = f"p.{first}{last}del"
    else:
        hgvs = f"p.{first}{last}delins{''.join(aa3(x) for x in as_[:10])}{'…' if len(as_) > 10 else ''}"
    return {"class": "inframe_indel", "hgvs": hgvs, "position": i + 1, "length_change": len(a) - len(r)}


# ---------------------------------------------------------------------------------- variants
def variant_region(tx: Dict[str, Any], start: int, end: int) -> Optional[str]:
    """Where ``[start, end)`` falls in ``tx`` (reference offsets), or None if beyond the flanks."""
    strand = tx["strand"]
    exons = sorted((int(s), int(e)) for s, e in tx["exons"])
    t0, t1 = exons[0][0], exons[-1][1]
    if end <= t0 - FLANK or start >= t1 + FLANK:
        return None
    if end <= t0:
        return "upstream" if strand == "+" else "downstream"
    if start >= t1:
        return "downstream" if strand == "+" else "upstream"
    cds = sorted((int(s), int(e)) for s, e in tx.get("cds") or [])
    for s, e in exons:
        if start < e and end > s:
            if cds and any(start < ce and end > cs for cs, ce in cds):
                label = "coding"
            elif cds:
                before = end <= cds[0][0]
                label = "5' UTR" if before == (strand == "+") else "3' UTR"
            else:
                label = "non-coding exon"
            near = min(abs(start - s), abs(e - end))
            internal_edge = (s != t0 and abs(start - s) < SPLICE_REGION_EXON) or (e != t1 and abs(e - end) < SPLICE_REGION_EXON)
            return f"{label}, splice region" if near < SPLICE_REGION_EXON and internal_edge else label
    for a, b in introns_of(exons):
        if start < b and end > a:
            d_left = start - a  # into the intron from its 5' (genomic left) end
            d_right = b - end
            left_kind, right_kind = ("donor", "acceptor") if strand == "+" else ("acceptor", "donor")
            if d_left < 2:
                return f"splice {left_kind}"
            if d_right < 2:
                return f"splice {right_kind}"
            if d_left < SPLICE_REGION_INTRON or d_right < SPLICE_REGION_INTRON:
                return "splice region (intron)"
            return "intron"
    return "intron"


# ---------------------------------------------------------------------------------- junctions
@dataclass
class Junctions:
    """Predicted junctions on reference offsets: ``(start, end) -> values per track``."""

    tracks: List[Dict[str, Any]]
    values: Dict[Tuple[int, int], np.ndarray]

    @classmethod
    def load(cls, npz: Path, view: HaplotypeView, strand: str) -> Optional["Junctions"]:
        if not npz.exists():
            return None
        with np.load(npz) as z:
            arrays = {k: z[k] for k in ("starts", "ends", "strands", "values")}
        meta_path = npz.with_name(f"{npz.stem}_metadata.json")
        tracks: List[Dict[str, Any]] = []
        if meta_path.exists():
            with open(meta_path, "r", encoding="utf-8") as f:
                tracks = list(json.load(f).get("metadata") or [])
        return cls.from_arrays(arrays, tracks, view, strand)

    @classmethod
    def from_arrays(cls, arrays: Dict[str, np.ndarray], tracks: List[Dict[str, Any]], view: HaplotypeView, strand: str) -> "Junctions":
        """Junctions predicted on a haplotype (offsets in its own coordinates), mapped to reference offsets."""
        starts, ends, strands, values = arrays["starts"], arrays["ends"], arrays["strands"], arrays["values"]
        out: Dict[Tuple[int, int], np.ndarray] = {}
        n = view.inverse.size
        for s, e, st, v in zip(starts.tolist(), ends.tolist(), strands.tolist(), values):
            if st != strand or s < 0 or e < 1 or e > n:
                continue
            rs, re_last = int(view.inverse[s]), int(view.inverse[e - 1])
            if rs < 0 or re_last < 0:
                continue
            key = (rs, re_last + 1)
            out[key] = out[key] + v if key in out else np.asarray(v, dtype=np.float64)
        return cls(tracks, out)

    def psi(self, strand: str, floor: np.ndarray) -> Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]]:
        """PSI5 (share of the donor's signal) and PSI3 (share of the acceptor's) per junction.

        Denominators are at least ``floor`` (per track), so a junction at a site with almost no
        signal does not get a large share by chance."""
        donor: Dict[int, np.ndarray] = {}
        acceptor: Dict[int, np.ndarray] = {}
        for (s, e), v in self.values.items():
            d, a = (s, e) if strand == "+" else (e, s)
            donor[d] = donor.get(d, 0) + v
            acceptor[a] = acceptor.get(a, 0) + v
        out = {}
        for (s, e), v in self.values.items():
            d, a = (s, e) if strand == "+" else (e, s)
            dt, at = donor[d], acceptor[a]
            out[(s, e)] = (v / np.maximum(dt, floor), v / np.maximum(at, floor))
        return out


def transcript_support(transcripts: List[Dict[str, Any]], junctions: Dict[str, "Junctions"], strand: str) -> Dict[str, Dict[str, np.ndarray]]:
    """Per transcript and haplotype: the usage of its weakest junction per track (as :func:`splicing`)."""
    ref = junctions[REFERENCE]
    n_tracks = len(next(iter(ref.values.values()))) if ref.values else max(len(ref.tracks), 1)
    zeros = np.zeros(n_tracks)
    expression = np.zeros(n_tracks)
    for t in transcripts:
        for key in introns_of(t["exons"]):
            expression = np.maximum(expression, ref.values.get(key, zeros))
    floor = np.maximum(SITE_FLOOR * expression, 1e-6)
    psi = {name: j.psi(strand, floor) for name, j in junctions.items()}
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for t in transcripts:
        ints = introns_of(t["exons"])
        if not ints:
            continue
        out[t["id"]] = {name: np.min(np.stack([np.minimum(*psi[name][k]) if k in psi[name] else zeros for k in ints]), axis=0) for name in junctions}
    return out


def track_labels(tracks: List[Dict[str, Any]], count: int) -> List[Dict[str, Any]]:
    out = []
    for i in range(count):
        t = tracks[i] if i < len(tracks) else {}
        name = t.get("biosample_name") or t.get("name") or f"track {i}"
        assay = t.get("Assay title") or ""
        out.append({"index": i, "ontology": t.get("ontology_curie"), "biosample": name, "assay": assay, "source": t.get("data_source"),
                    "label": f"{name}{' · ' + assay if assay else ''}"})
    return out


def _r(values: np.ndarray, digits: int = 3) -> List[float]:
    return [round(float(x), digits) for x in np.asarray(values, dtype=np.float64)]


def junction_kind(junction: Tuple[int, int], base_introns: List[Tuple[int, int]], strand: str) -> str:
    s, e = junction
    starts = {a for a, _ in base_introns}
    ends = {b for _, b in base_introns}
    if s in starts and e in ends:
        skipped = sum(1 for a, _ in base_introns if s < a < e)
        return f"exon skipping ({skipped} exon{'s' if skipped != 1 else ''})" if skipped else "annotated intron"
    left_known, right_known = s in starts, e in ends
    five, three = ("donor", "acceptor") if strand == "+" else ("acceptor", "donor")
    if left_known and not right_known:
        return f"alternative {three}"
    if right_known and not left_known:
        return f"alternative {five}"
    return "novel junction (both sites)"


# ---------------------------------------------------------------------------------- service
class ProductService:
    """Builds :func:`products` payloads for the server, with a small in-memory cache."""

    def __init__(self, signals, annotations, max_items: int = 32):
        self.signals = signals
        self.annotations = annotations
        self._cache: "OrderedDict[Tuple, Dict[str, Any]]" = OrderedDict()
        self._lock = threading.Lock()
        self.max_items = max_items

    def views(self, dataset: Dataset, sample: str, gene: str) -> Dict[str, HaplotypeView]:
        views = {REFERENCE: HaplotypeView.reference(self.signals.reference_sequence(dataset, gene))}
        for hap in HAPLOTYPES:
            path = dataset.haplotype_fasta(sample, gene, hap)
            if not path.exists():
                raise FileNotFoundError(f"No {hap} sequence for {sample}/{gene}")
            views[hap] = HaplotypeView.from_map(hap, self.signals.haplotype_sequence(dataset, sample, gene, hap), self.signals.ref_map(dataset, sample, gene, hap))
        return views

    def junction_paths(self, dataset: Dataset, sample: str, gene: str) -> Dict[str, Path]:
        paths = {REFERENCE: dataset.reference_prediction_path(gene, "splice_junctions")}
        for hap in HAPLOTYPES:
            paths[hap] = dataset.prediction_path(sample, gene, hap, "splice_junctions")
        return paths

    def products(self, dataset: Dataset, sample: str, gene: str, models: List[Dict[str, Any]], target: Optional[str] = None,
                 with_sequence: bool = False) -> Dict[str, Any]:
        paths = self.junction_paths(dataset, sample, gene)
        stamp = tuple(p.stat().st_mtime_ns if p.exists() else 0 for p in paths.values())
        key = (dataset.fingerprint, sample, gene, target, with_sequence, stamp, id(models))
        with self._lock:
            hit = self._cache.get(key)
            if hit is not None:
                self._cache.move_to_end(key)
                return hit
        payload = products(dataset, sample, gene, models, self.views(dataset, sample, gene), paths, self.signals.variants(dataset, sample, gene),
                           target=target, with_sequence=with_sequence)
        with self._lock:
            self._cache[key] = payload
            while len(self._cache) > self.max_items:
                self._cache.popitem(last=False)
        return payload


def coding_start_stop(tx: Dict[str, Any]) -> Tuple[Optional[int], Optional[int]]:
    """Reference offsets of the start codon's A and the stop codon's first base (GTF CDS excludes the stop)."""
    cds = sorted((int(s), int(e)) for s, e in tx.get("cds") or [])
    if not cds:
        return None, None
    if tx["strand"] == "+":
        return cds[0][0], cds[-1][1]
    return cds[-1][1] - 1, cds[0][0] - 1


def _incomplete(tx: Dict[str, Any]) -> Dict[str, bool]:
    tags = set(tx.get("tags") or [])
    return {"start": bool(tags & {"cds_start_NF", "mRNA_start_NF"}), "end": bool(tags & {"cds_end_NF", "mRNA_end_NF"})}


def _variant_rows(variants, window_start: int, base: Dict[str, Any]) -> List[Dict[str, Any]]:
    rows = []
    for i, pos in enumerate(variants.positions.tolist()):
        carried = variants.carried[i].tolist()
        if not any(c > 0 for c in carried):
            continue
        ref = variants.refs[i]
        start = int(pos) - window_start
        region = variant_region(base, start, start + max(len(ref), 1))
        if region is None:
            continue
        alts = variants.alt_alleles[i]
        rows.append({
            "pos": int(pos), "offset": start, "id": variants.ids[i], "ref": ref,
            "H1": alts[carried[0] - 1] if carried[0] > 0 else None,
            "H2": alts[carried[1] - 1] if len(carried) > 1 and carried[1] > 0 else None,
            "genotype": variants.genotypes[i], "region": region,
        })
    return rows


def products(dataset: Dataset, sample: str, gene: str, models: List[Dict[str, Any]], views: Dict[str, HaplotypeView],
             junction_paths: Dict[str, Path], variants, target: Optional[str] = None, with_sequence: bool = False) -> Dict[str, Any]:
    window = dataset.window(gene)
    length = int(window.length or views[REFERENCE].seq.size)
    names = [g["name"] for g in models]
    chosen = next((g for g in models if g["name"] == (target or gene)), None) or next((g for g in models if g["name"] == gene), None)
    if chosen is None:
        coding = [g for g in models if any(t.get("cds") for t in g["transcripts"])]
        chosen = (coding or models or [None])[0]
    if chosen is None:
        return {"gene": gene, "sample": sample, "genes": names, "target": None, "transcripts": [], "message": "No gene model in this window"}
    strand = chosen["strand"]
    complete = [t for t in chosen["transcripts"] if t["exons"] and min(s for s, _ in t["exons"]) >= 0 and max(e for _, e in t["exons"]) <= length]
    truncated = [t["name"] for t in chosen["transcripts"] if t not in complete]
    if not complete:
        return {"gene": gene, "sample": sample, "genes": names, "target": chosen["name"], "transcripts": [], "truncated": truncated,
                "message": "Every transcript of this gene extends past the window"}
    base = next((t for t in complete if t.get("mane")), None) or next((t for t in complete if t.get("canonical")), None) or complete[0]

    # -- annotated transcripts on every haplotype
    rows = []
    reference_products: Dict[str, Dict[str, Any]] = {}
    for tx in complete:
        start, stop = coding_start_stop(tx)
        inc = _incomplete(tx)
        ref = build_product(views[REFERENCE], tx["exons"], strand, start, stop, inc["start"], with_sequence=True, phase=tx.get("cds_phase", 0))
        reference_products[tx["id"]] = ref
        row = {k: tx.get(k) for k in ("id", "name", "type", "mane", "canonical", "start", "end")}
        row.update(exons=tx["exons"], cds=tx.get("cds") or [], incomplete=inc, products={REFERENCE: _public(ref, with_sequence)})
        for hap in HAPLOTYPES:
            alt = build_product(views[hap], tx["exons"], strand, start, stop, inc["start"], with_sequence=True, phase=tx.get("cds_phase", 0))
            public = _public(alt, with_sequence)
            public["change"] = describe_change(ref, alt)
            row["products"][hap] = public
        rows.append(row)

    # -- variants carried, relative to the base transcript
    variant_rows = _variant_rows(variants, int(window.start or 1), base) if variants is not None else []

    payload: Dict[str, Any] = {
        "gene": gene, "sample": sample, "genes": names, "target": chosen["name"], "strand": strand, "chromosome": window.chromosome,
        "window_start": window.start, "gene_type": chosen.get("type"), "base": base["id"], "base_name": base["name"],
        "transcripts": rows, "truncated": truncated, "variants": variant_rows,
    }
    payload["splicing"] = splicing(chosen, base, complete, views, junction_paths, strand, reference_products, with_sequence)
    return payload


def _public(product: Dict[str, Any], with_sequence: bool) -> Dict[str, Any]:
    return {k: v for k, v in product.items() if k != "mrna" or with_sequence}


def splicing(gene_model: Dict[str, Any], base: Dict[str, Any], transcripts: List[Dict[str, Any]], views: Dict[str, HaplotypeView],
             junction_paths: Dict[str, Path], strand: str, reference_products: Dict[str, Dict[str, Any]], with_sequence: bool) -> Dict[str, Any]:
    loaded = {name: Junctions.load(path, views[name], strand) for name, path in junction_paths.items()}
    missing = [name for name, j in loaded.items() if j is None]
    if loaded.get(REFERENCE) is None:
        return {"available": False, "missing": missing,
                "message": "No splice_junctions prediction for the reference window; run predict-dataset with SPLICE_JUNCTIONS for --haplotypes ref"}
    gs = min(min(s for s, _ in t["exons"]) for t in transcripts)
    ge = max(max(e for _, e in t["exons"]) for t in transcripts)
    ref_j = loaded[REFERENCE]
    n_tracks = len(next(iter(ref_j.values.values()))) if ref_j.values else len(ref_j.tracks)
    zeros = np.zeros(n_tracks)
    annotated: Dict[Tuple[int, int], List[str]] = {}
    for t in transcripts:
        for intron in introns_of(t["exons"]):
            annotated.setdefault(intron, []).append(t["name"])
    base_introns = introns_of(base["exons"])

    def vec(name: str, key) -> np.ndarray:
        j = loaded.get(name)
        return j.values.get(key, zeros) if j is not None else zeros

    expression = np.zeros(n_tracks)
    for key in annotated:
        expression = np.maximum(expression, vec(REFERENCE, key))
    expressed = expression >= EXPRESSED_JUNCTION
    floor = np.maximum(SITE_FLOOR * expression, 1e-6)
    psi = {name: j.psi(strand, floor) for name, j in loaded.items() if j is not None}

    def usage(name: str, key) -> np.ndarray:
        p = psi.get(name, {}).get(key)
        # the share at the junction's more contested end: exon skipping lowers PSI5 of the intron
        # before the skipped exon but leaves its PSI3 at 1 (nothing else uses that acceptor)
        return np.minimum(p[0], p[1]) if p is not None else zeros

    def junction_row(key, haps: Sequence[str]) -> Dict[str, Any]:
        row = {"start": key[0], "end": key[1], "annotated": annotated.get(key, []), "ref": {"signal": _r(vec(REFERENCE, key)), "usage": _r(usage(REFERENCE, key))}}
        for hap in haps:
            u = usage(hap, key)
            row[hap] = {"signal": _r(vec(hap, key)), "usage": _r(u), "delta": _r(u - usage(REFERENCE, key))}
        return row

    haps = [h for h in HAPLOTYPES if loaded.get(h) is not None]
    base_rows = [junction_row(key, haps) for key in base_introns]

    # -- annotated transcript support: weakest junction usage per track
    support = {}
    for t in transcripts:
        ints = introns_of(t["exons"])
        if not ints:
            continue
        support[t["id"]] = {name: _r(np.min(np.stack([usage(name, k) for k in ints]), axis=0)) for name in [REFERENCE] + haps}

    # -- candidate isoforms (junctions in the gene span, on its strand, not introns of the base)
    keys = set()
    for name in [REFERENCE] + haps:
        keys.update(k for k in loaded[name].values if gs <= k[0] and k[1] <= ge)
    scored = []
    base_set = set(base_introns)
    for key in keys:
        if key in base_set:
            continue
        strong = expressed & (np.max(np.stack([vec(name, key) for name in [REFERENCE] + haps]), axis=0) >= JUNCTION_FLOOR * expression)
        if not strong.any():
            continue
        best_usage, best_delta = 0.0, 0.0
        for hap in haps:
            u = usage(hap, key)
            d = u - usage(REFERENCE, key)
            best_usage = max(best_usage, float(np.max(np.where(strong, u, 0))))
            best_delta = max(best_delta, float(np.max(np.where(strong, np.abs(d), 0))))
        ref_usage = float(np.max(np.where(strong, usage(REFERENCE, key), 0)))
        if max(best_usage, ref_usage) >= CANDIDATE_USAGE or best_delta >= CANDIDATE_DELTA:
            scored.append((best_delta, max(best_usage, ref_usage), key))
    scored.sort(key=lambda x: (-x[0], -x[1]))
    start, stop = coding_start_stop(base)
    inc = _incomplete(base)
    base_ref = reference_products[base["id"]]
    candidates = []
    for delta, top, key in scored[:MAX_CANDIDATES]:
        exons = splice_in(base["exons"], key)
        candidates.append(_candidate(views, exons, strand, start, stop, inc, base_ref, haps, with_sequence, base.get("cds_phase", 0),
                                     kind=junction_kind(key, base_introns, strand), junction=junction_row(key, haps), max_delta=round(delta, 3)))
    for hap in haps:
        if not base_introns:
            break
        # per-intron hap/ref signal, relative to the gene's median intron: a whole-gene expression
        # change moves every intron alike and is not retention
        ratios = np.stack([vec(hap, k) / np.maximum(vec(REFERENCE, k), 1e-12) for k in base_introns])
        relative = ratios / np.maximum(np.median(ratios, axis=0), 1e-12)
        for i, key in enumerate(base_introns):
            ref_signal = vec(REFERENCE, key)
            ratio = np.where(expressed & (ref_signal >= RETENTION_MIN_SIGNAL * expression), relative[i], 1.0)
            if np.min(ratio) <= RETENTION_RATIO and not any(c["junction"]["start"] == key[0] and c["junction"]["end"] == key[1] for c in candidates):
                candidates.append(_candidate(views, retain_intron(base["exons"], key), strand, start, stop, inc, base_ref, haps, with_sequence, base.get("cds_phase", 0),
                                             kind="intron retention (junction lost)", junction=junction_row(key, haps), max_delta=round(float(1 - np.min(ratio)), 3),
                                             note="inferred from the loss of the junction; AlphaGenome does not predict retention directly"))
    return {
        "available": True, "missing": missing, "tracks": track_labels(ref_j.tracks, n_tracks),
        "expression": _r(expression), "expressed": [bool(x) for x in expressed],
        "base_introns": base_rows, "support": support, "candidates": candidates,
        "thresholds": {"expressed": EXPRESSED_JUNCTION, "usage": CANDIDATE_USAGE, "delta": CANDIDATE_DELTA, "retention_ratio": RETENTION_RATIO,
                       "retention_min_signal": RETENTION_MIN_SIGNAL, "site_floor": SITE_FLOOR, "junction_floor": JUNCTION_FLOOR},
    }


def _candidate(views, exons, strand, start, stop, inc, base_ref, haps, with_sequence, phase, **info) -> Dict[str, Any]:
    out = dict(info)
    out["exons"] = exons
    out["products"] = {}
    for name in [REFERENCE] + list(haps):
        p = build_product(views[name], exons, strand, start, stop, inc["start"], with_sequence=True, phase=phase)
        public = _public(p, with_sequence)
        public["change"] = describe_change(base_ref, p)
        out["products"][name] = public
    return out


def fasta(payload: Dict[str, Any], kind: str = "protein") -> str:
    """FASTA of every product in a :func:`products` payload (``kind``: protein or mrna)."""
    lines: List[str] = []
    sample, gene = payload.get("sample"), payload.get("target")

    def emit(name: str, hap: str, product: Dict[str, Any], extra: str = "") -> None:
        seq = product.get("protein") if kind == "protein" else product.get("mrna")
        if not seq:
            return
        change = (product.get("change") or {}).get("hgvs") or ""
        label = REFERENCE if hap == REFERENCE else f"{sample}.{hap}"
        lines.append(f">{gene}|{name}|{label}{'|' + extra if extra else ''}{'|' + change if change else ''}")
        lines.extend(seq[i:i + 60] for i in range(0, len(seq), 60))

    for tx in payload.get("transcripts", []):
        for hap, product in tx["products"].items():
            emit(tx["name"], hap, product)
    for i, cand in enumerate((payload.get("splicing") or {}).get("candidates") or []):
        for hap, product in cand["products"].items():
            if hap != REFERENCE:
                emit(f"{payload.get('base_name')}~cand{i + 1}", hap, product, cand.get("kind", ""))
    return "\n".join(lines) + ("\n" if lines else "")


# ---------------------------------------------------------------------------------- per-variant effects
# Most to least severe; a variant's consequence is the worst over the gene's coding transcripts.
SEVERITY = ["frameshift", "nonsense", "start_lost", "stop_lost", "splice", "inframe_indel", "missense", "splice_region",
            "synonymous", "utr", "noncoding_exon", "intron", "flank", "structural"]
_CLASS_TO_CATEGORY = {"frameshift": "frameshift", "stop_gained": "nonsense", "start_lost": "start_lost", "stop_lost": "stop_lost",
                      "inframe_indel": "inframe_indel", "missense": "missense", "synonymous": "synonymous", "utr_change": "utr",
                      "no_change": None, "noncoding_change": "noncoding_exon"}


def _region_category(region: str) -> str:
    if region.startswith("splice donor") or region.startswith("splice acceptor"):
        return "splice"
    if "splice region" in region:
        return "splice_region"
    if region == "coding":
        return "synonymous"
    if "UTR" in region:
        return "utr"
    if region == "non-coding exon":
        return "noncoding_exon"
    if region == "intron":
        return "intron"
    return "flank"


def single_variant_view(reference: HaplotypeView, offset: int, ref: str, alt: str) -> HaplotypeView:
    """The reference window with one variant applied (VCF alleles, ``offset`` = POS - window start)."""
    seq = reference.seq.tobytes().decode("ascii")
    edited = seq[:offset] + alt + seq[offset + len(ref):]
    events = HaplotypeEvents(np.array([offset], np.int64), np.array([len(ref)], np.int32), np.array([len(alt)], np.int32))
    return HaplotypeView.from_map("variant", edited.encode("ascii"), reference_map(events, reference.local.size))


def variant_effects(reference: HaplotypeView, transcripts: List[Dict[str, Any]], base: Dict[str, Any], strand: str,
                    reference_products: Dict[str, Dict[str, Any]], variants, window_start: int) -> List[Dict[str, Any]]:
    """Each variant the sample carries in the gene, applied alone to the reference genome.

    Per variant: its consequence on every coding transcript (synonymous, missense, nonsense, ...),
    the worst of them (``category``), the amino-acid change on the base transcript, and which
    haplotypes carry it. The haplotype's combined effect (all its variants together) is in the
    products themselves; this is the per-variant view of the same changes.
    """
    if variants is None:
        return []
    coding = [t for t in transcripts if t.get("cds")]
    exons = [(int(s), int(e)) for t in transcripts for s, e in t["exons"]]
    rows = []
    for i, pos in enumerate(variants.positions.tolist()):
        carried = variants.carried[i].tolist()
        if not any(c > 0 for c in carried):
            continue
        ref = variants.refs[i]
        offset = int(pos) - window_start
        region = variant_region(base, offset, offset + max(len(ref), 1))
        if region is None:
            continue
        alleles = sorted({c for c in carried if c > 0})
        for allele in alleles:
            alt = variants.alt_alleles[i][allele - 1]
            haps = [h for h, c in zip(HAPLOTYPES, carried) if c == allele]
            row = {"pos": int(pos), "offset": offset, "id": variants.ids[i], "ref": ref, "alt": alt, "haplotypes": haps,
                   "genotype": variants.genotypes[i], "region": region, "transcripts": {}}
            category = _region_category(region)
            near_exon = any(s - 10 <= offset < e + 10 for s, e in exons)
            if alt.startswith("<") or not re.fullmatch(r"[ACGTN]+", alt) or not re.fullmatch(r"[ACGTN]+", ref):
                category = "structural"
            elif near_exon and coding:
                view = single_variant_view(reference, offset, ref, alt)
                worst = None
                for t in coding:
                    start, stop = coding_start_stop({**t, "strand": strand})
                    inc = _incomplete(t)
                    product = build_product(view, t["exons"], strand, start, stop, inc["start"], with_sequence=True, phase=t.get("cds_phase", 0))
                    change = describe_change(reference_products[t["id"]], product)
                    cat = _CLASS_TO_CATEGORY.get(change["class"])
                    if cat is None:
                        continue
                    if (product.get("nmd") or {}).get("predicted") and not (reference_products[t["id"]].get("nmd") or {}).get("predicted"):
                        change = {**change, "nmd": True}
                    row["transcripts"][t["name"]] = {"class": change["class"], "hgvs": change.get("hgvs", ""), "nmd": bool(change.get("nmd"))}
                    if worst is None or SEVERITY.index(cat) < SEVERITY.index(worst):
                        worst = cat
                if worst is not None and SEVERITY.index(worst) < SEVERITY.index(category):
                    category = worst
                elif worst is None and category == "synonymous":
                    category = "utr"  # coding in the base model but changes no transcript's CDS
                base_change = row["transcripts"].get(base["name"])
                row["base_change"] = base_change
            row["category"] = category
            rows.append(row)
    rows.sort(key=lambda r: (SEVERITY.index(r["category"]), r["pos"]))
    return rows
