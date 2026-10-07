"""3D protein structures and mRNA secondary structures for the gene-products view.

**Proteins.** The reference protein's structure comes from the AlphaFold Protein Structure Database:
the transcript's UniProt entry is found through UniProt's Ensembl cross-reference, and its model
is used when its sequence is the transcript's protein (else the offset between the two is used, or
the model is not shown). A haplotype's product is then placed on that model:

``identical``      same protein: the reference model
``substitution``   same length, some residues differ: the reference model with them marked (a
                   predicted fold does not resolve a point mutation's effect on stability)
``truncated``      the product is a prefix of the reference (stop gained, NMD-escaping): the model
                   with the residues it lacks greyed out
``altered``        anything else (frameshift tail, in-frame indel, another isoform): only a new
                   prediction can show it. ESMFold folds a window of at most 400 residues around the
                   first change, for the product and the reference alike so the two are comparable.

ESMFold runs on ``api.esmatlas.com``: the product's sequence is sent there, so it runs only when
the user asks and never with ``--no-remote``. Results are cached by sequence under
``<cache>/structures``.

**mRNA.** Full-length mRNA has no reliable 3D structure prediction (3D RNA predictors handle short,
isolated RNAs; mRNA in cells is bound by proteins). The secondary structure (minimum free energy,
ViennaRNA) of a window of the mature transcript is the established proxy: around the start codon
(translation initiation) or around a variant, for the reference and the haplotype, with ΔG of
both. Needs the ``ViennaRNA`` package (``RNA``).
"""
from __future__ import annotations

import hashlib
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.visualizer.products import HaplotypeView
from genomics.visualizer.remote import DAY, USER_AGENT, RemoteCache, RemoteError

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAX = 400
ESMFOLD_CONTEXT = 150  # residues of context before the first change in an ESMFold window
RNA_WINDOW = 240
RNA_WINDOW_MAX = 600


class StructureError(RuntimeError):
    pass


def compare_proteins(ref: str, alt: str) -> Dict[str, Any]:
    """How ``alt`` relates to ``ref``, residue numbers 1-based."""
    if alt == ref:
        return {"mode": "identical"}
    if len(alt) == len(ref):
        subs = [{"pos": i + 1, "ref": a, "alt": b} for i, (a, b) in enumerate(zip(ref, alt)) if a != b]
        return {"mode": "substitution", "substitutions": subs}
    if alt and ref.startswith(alt):
        return {"mode": "truncated", "length": len(alt), "lost_from": len(alt) + 1}
    first = next((i for i, (a, b) in enumerate(zip(ref, alt)) if a != b), min(len(ref), len(alt)))
    return {"mode": "altered", "first_change": first + 1}


def fold_window(length: int, first_change: int, limit: int = ESMFOLD_MAX) -> List[int]:
    """1-based inclusive residue window of at most ``limit`` with context before ``first_change``."""
    if length <= limit:
        return [1, length]
    start = max(1, first_change - ESMFOLD_CONTEXT)
    end = min(length, start + limit - 1)
    start = max(1, end - limit + 1)
    return [start, end]


def mean_b_factor(pdb: str) -> Optional[float]:
    values = [float(line[60:66]) for line in pdb.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[60:66].strip()]
    return round(float(np.mean(values)), 2) if values else None


class StructureService:
    def __init__(self, remote: RemoteCache, cache_dir: Optional[Path]):
        self.remote = remote
        self.cache_dir = Path(cache_dir) / "structures" if cache_dir else None

    # -- reference model ---------------------------------------------------------------------
    def uniprot(self, transcript_id: str, protein: str) -> Optional[Dict[str, Any]]:
        tid = transcript_id.split(".")[0]
        data = self.remote.get_json(UNIPROT_SEARCH, ttl=90 * DAY, params={
            "query": f"xref:ensembl-{tid}", "fields": "accession,sequence,reviewed,protein_name", "format": "json", "size": 10})
        results = (data or {}).get("results") or []
        if not results:
            return None

        def seq(r):
            return ((r.get("sequence") or {}).get("value")) or ""

        exact = [r for r in results if seq(r) == protein]
        reviewed = [r for r in results if "Swiss-Prot" in str(r.get("entryType", ""))]
        best = (exact or reviewed or results)[0]
        name = (((best.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}).get("value")
        return {"accession": best["primaryAccession"], "sequence": seq(best), "name": name, "reviewed": best in reviewed}

    def reference_model(self, transcript_id: str, protein: str) -> Dict[str, Any]:
        """AlphaFold DB model of the transcript's protein, with how its numbering maps to ours."""
        try:
            entry = self.uniprot(transcript_id, protein)
        except RemoteError as exc:
            raise StructureError(f"UniProt: {exc}")
        if entry is None:
            raise StructureError(f"No UniProt entry cross-references {transcript_id}")
        try:
            models = self.remote.get_json(AFDB_API.format(accession=entry["accession"]), ttl=90 * DAY)
        except RemoteError as exc:
            raise StructureError(f"AlphaFold DB has no model for {entry['accession']} ({exc})")
        if not models:
            raise StructureError(f"AlphaFold DB has no model for {entry['accession']}")
        model = models[0]
        afdb_seq = model.get("uniprotSequence") or model.get("sequence") or ""
        start = int(model.get("uniprotStart") or 1)
        # offset: our residue i is model residue i + offset
        if afdb_seq == protein:
            offset, match = 0, "identical"
        elif protein and protein in afdb_seq:
            offset, match = afdb_seq.index(protein), "contained"
        elif afdb_seq and afdb_seq in protein:
            offset, match = -protein.index(afdb_seq), "partial"
        else:
            offset, match = None, "different"
        pdb_url = model.get("pdbUrl")
        pdb = self.remote.get_text(pdb_url, ttl=180 * DAY) if pdb_url else ""
        return {
            "source": "AlphaFold DB", "accession": entry["accession"], "protein_name": entry.get("name"), "model": model.get("entryId"),
            "version": model.get("latestVersion"), "url": f"https://alphafold.ebi.ac.uk/entry/{entry['accession']}",
            "match": match, "offset": offset, "model_start": start, "pdb": pdb, "plddt_mean": mean_b_factor(pdb), "plddt_scale": 100,
        }

    # -- ESMFold ---------------------------------------------------------------------------
    def esmfold(self, sequence: str, allow_remote: bool) -> Dict[str, Any]:
        if not sequence:
            raise StructureError("Empty sequence")
        if len(sequence) > ESMFOLD_MAX:
            raise StructureError(f"ESMFold's public API folds at most {ESMFOLD_MAX} residues")
        key = hashlib.sha1(sequence.encode()).hexdigest()
        path = self.cache_dir / "esmfold" / f"{key}.pdb" if self.cache_dir else None
        if path is not None and path.exists():
            pdb = path.read_text(encoding="utf-8")
        else:
            if not allow_remote:
                raise StructureError("ESMFold needs api.esmatlas.com, and remote lookups are disabled (--no-remote)")
            request = urllib.request.Request(ESMFOLD_API, data=sequence.encode("ascii"), headers={"User-Agent": USER_AGENT, "Content-Type": "text/plain"})
            try:
                with urllib.request.urlopen(request, timeout=300) as response:
                    pdb = response.read().decode("utf-8")
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                raise StructureError(f"ESMFold (api.esmatlas.com) failed: {exc}")
            if "ATOM" not in pdb:
                raise StructureError(f"ESMFold returned no structure: {pdb[:200]}")
            if path is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(pdb, encoding="utf-8")
        return {"source": "ESMFold", "pdb": pdb, "plddt_mean": mean_b_factor(pdb), "plddt_scale": 1, "length": len(sequence)}

    def fold_pair(self, ref: str, alt: str, allow_remote: bool) -> Dict[str, Any]:
        """ESMFold of the same window (≤ 400 residues around the first change) of both proteins."""
        info = compare_proteins(ref, alt)
        first = info.get("first_change") or info.get("lost_from") or 1
        win_alt = fold_window(len(alt), first)
        win_ref = fold_window(len(ref), first)
        ref_seq, alt_seq = ref[win_ref[0] - 1:win_ref[1]], alt[win_alt[0] - 1:win_alt[1]]
        ref_model, alt_model = self.esmfold(ref_seq, allow_remote), self.esmfold(alt_seq, allow_remote)
        return {
            "window": {"reference": win_ref, "product": win_alt},
            "reference": ref_model,
            "product": alt_model,
            "first_change": first,
            "structure_similarity": structure_similarity(ref_model["pdb"], alt_model["pdb"], ref_seq, alt_seq),
            "sequence_similarity": protein_similarity(ref, alt),
        }


# ---------------------------------------------------------------------------------- RNA
def rna_available() -> bool:
    import importlib.util

    return importlib.util.find_spec("RNA") is not None


def fold_rna(sequence: str) -> Dict[str, Any]:
    """Minimum-free-energy secondary structure (ViennaRNA) and a 2D layout of ``sequence``."""
    try:
        import RNA
    except ImportError as exc:
        raise StructureError("RNA folding needs ViennaRNA: pip install ViennaRNA (or pip install -e '.[visualizer]')") from exc
    rna = sequence.upper().replace("T", "U")
    structure, mfe = RNA.fold(rna)
    coords = RNA.simple_xy_coordinates(structure)
    xy = [[round(float(c.X), 2), round(float(c.Y), 2)] for c in list(coords)[:len(rna)]]
    pairs = []
    stack: List[int] = []
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            pairs.append([stack.pop(), i])
    return {"sequence": rna, "structure": structure, "mfe": round(float(mfe), 2), "xy": xy, "pairs": pairs,
            "paired_fraction": round(2 * len(pairs) / max(len(rna), 1), 3)}


def mrna_window(view: HaplotypeView, exons: Sequence[Sequence[int]], strand: str, center_offset: int, width: int) -> Dict[str, Any]:
    """A window of the spliced transcript centred on reference offset ``center_offset``."""
    mrna, idx, _lengths = view.spliced(exons, strand)
    ref_offsets = view.inverse[idx] if idx.size else np.zeros(0, dtype=np.int64)
    if not len(mrna):
        raise StructureError("Empty transcript")
    hits = np.nonzero(ref_offsets == center_offset)[0]
    if hits.size:
        center = int(hits[0])
    else:  # the centre base is deleted on this haplotype: nearest base still present
        valid = np.nonzero(ref_offsets >= 0)[0]
        center = int(valid[np.argmin(np.abs(ref_offsets[valid] - center_offset))]) if valid.size else 0
    half = width // 2
    start = max(0, min(center - half, len(mrna) - width))
    end = min(len(mrna), start + width)
    return {"sequence": mrna[start:end], "start": start, "end": end, "ref_offsets": ref_offsets[start:end].tolist(), "center": center}


def rna_comparison(ref_view: HaplotypeView, hap_view: HaplotypeView, exons, strand: str, center_offset: int, width: int = RNA_WINDOW) -> Dict[str, Any]:
    width = int(max(40, min(width, RNA_WINDOW_MAX)))
    ref = mrna_window(ref_view, exons, strand, center_offset, width)
    hap = mrna_window(hap_view, exons, strand, center_offset, width)
    ref_by_offset = {o: b for o, b in zip(ref["ref_offsets"], ref["sequence"])}
    changed = [i for i, (o, b) in enumerate(zip(hap["ref_offsets"], hap["sequence"])) if o < 0 or ref_by_offset.get(o, b) != b]
    ref_fold, hap_fold = fold_rna(ref["sequence"]), fold_rna(hap["sequence"])
    return {
        "width": width,
        "reference": {**ref_fold, "start": ref["start"] + 1, "end": ref["end"], "ref_offsets": ref["ref_offsets"]},
        "haplotype": {**hap_fold, "start": hap["start"] + 1, "end": hap["end"], "changed": changed, "ref_offsets": hap["ref_offsets"]},
        "ddg": round(hap_fold["mfe"] - ref_fold["mfe"], 2),
        "identical": ref["sequence"] == hap["sequence"],
    }


# ---------------------------------------------------------------------------------- similarity
_BLOSUM62_ROWS = """\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def _blosum62() -> Tuple[Dict[str, int], np.ndarray]:
    lines = _BLOSUM62_ROWS.strip("\n").splitlines()
    letters = lines[0].split()
    matrix = np.array([[int(v) for v in line.split()[1:]] for line in lines[1:]], dtype=np.int16)
    return {a: i for i, a in enumerate(letters)}, matrix


BLOSUM_INDEX, BLOSUM62 = _blosum62()
GAP = -6
MAX_ALIGN_CELLS = 6_000_000


def _codes(seq: str) -> np.ndarray:
    x = BLOSUM_INDEX["X"]
    return np.array([BLOSUM_INDEX.get(c, x) for c in seq], dtype=np.int16)


def global_alignment(a: str, b: str, gap: int = GAP) -> List[Tuple[int, int]]:
    """Needleman-Wunsch (BLOSUM62, linear gaps): aligned residue pairs ``(i, j)`` (0-based).

    The common prefix and suffix are matched directly and only the middle is aligned; a middle
    larger than ``MAX_ALIGN_CELLS`` is left unaligned (gaps), which only underestimates similarity.
    """
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    s = 0
    while s < min(len(a), len(b)) - p and a[len(a) - 1 - s] == b[len(b) - 1 - s]:
        s += 1
    pairs = [(i, i) for i in range(p)]
    ma, mb = a[p:len(a) - s], b[p:len(b) - s]
    if ma and mb and len(ma) * len(mb) <= MAX_ALIGN_CELLS:
        ca, cb = _codes(ma), _codes(mb)
        n, m = len(ma), len(mb)
        H = np.zeros((n + 1, m + 1), dtype=np.float32)
        H[0] = gap * np.arange(m + 1)
        H[:, 0] = gap * np.arange(n + 1)
        js = np.arange(m + 1, dtype=np.float32)
        for i in range(1, n + 1):
            diag = H[i - 1, :-1] + BLOSUM62[ca[i - 1], cb]
            up = H[i - 1, 1:] + gap
            d = np.empty(m + 1, dtype=np.float32)
            d[0] = H[i, 0]
            d[1:] = np.maximum(diag, up)
            # H[j] = max(d[j], H[j-1] + gap) = gap*j + cummax(d[k] - gap*k)
            H[i] = gap * js + np.maximum.accumulate(d - gap * js)
        i, j = n, m
        middle = []
        while i > 0 and j > 0:
            if H[i, j] == H[i - 1, j - 1] + BLOSUM62[ca[i - 1], cb[j - 1]]:
                middle.append((p + i - 1, p + j - 1))
                i, j = i - 1, j - 1
            elif H[i, j] == H[i - 1, j] + gap:
                i -= 1
            else:
                j -= 1
        pairs.extend(reversed(middle))
    pairs.extend((len(a) - s + k, len(b) - s + k) for k in range(s))
    return pairs


def protein_similarity(ref: str, alt: str) -> Dict[str, Any]:
    """How similar ``alt`` is to ``ref``: identity and BLOSUM62 score over the global alignment.

    ``similarity`` = alignment score of ref vs alt (gaps and unaligned reference residues score
    nothing) / score of ref vs itself, in [0, 1]; ``identity`` = identical aligned residues / length
    of the reference; ``coverage`` = aligned reference residues / reference length.
    """
    if not ref:
        return {"similarity": None, "identity": None, "coverage": None}
    if ref == alt:
        return {"similarity": 1.0, "identity": 1.0, "coverage": 1.0, "length": len(alt), "reference_length": len(ref), "substitutions": 0}
    cr, ca = _codes(ref), _codes(alt) if alt else np.zeros(0, dtype=np.int16)
    self_score = float(BLOSUM62[cr, cr].sum())
    if len(ref) == len(alt):
        pairs = [(i, i) for i in range(len(ref))]
    else:
        pairs = global_alignment(ref, alt) if alt else []
    if pairs:
        ia = np.array([p[0] for p in pairs])
        ib = np.array([p[1] for p in pairs])
        score = float(BLOSUM62[cr[ia], ca[ib]].sum())
        identical = int(np.count_nonzero(cr[ia] == ca[ib]))
    else:
        score, identical = 0.0, 0
    return {
        "similarity": round(max(0.0, min(1.0, score / self_score)), 4) if self_score > 0 else None,
        "identity": round(identical / len(ref), 4),
        "coverage": round(len(pairs) / len(ref), 4),
        "length": len(alt), "reference_length": len(ref),
        "substitutions": int(len(pairs) - identical),
    }


def ca_coordinates(pdb: str) -> Tuple[List[int], np.ndarray, np.ndarray]:
    """Residue numbers, CA coordinates and B-factors (pLDDT) of the first chain of a PDB file."""
    resi, xyz, b = [], [], []
    chain = None
    for line in pdb.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if chain is None:
            chain = line[21]
        if line[21] != chain:
            continue
        resi.append(int(line[22:26]))
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        b.append(float(line[60:66]) if line[60:66].strip() else 0.0)
    return resi, np.array(xyz, dtype=np.float64).reshape(-1, 3), np.array(b, dtype=np.float64)


def _kabsch(p: np.ndarray, q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotation and translation that best map ``p`` onto ``q``."""
    pc, qc = p.mean(axis=0), q.mean(axis=0)
    u, _s, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return rot, qc - pc @ rot


def tm_score(x: np.ndarray, y: np.ndarray, length: int) -> Dict[str, Any]:
    """TM-score of the aligned CA pairs ``x[k] <-> y[k]``, normalised by ``length`` (the reference).

    Superpositions are searched from fragments of the alignment and refined on the pairs within
    the distance cut-off, as in TM-score's fixed-alignment search. 1 = same fold; > 0.5 = same fold
    family; < 0.3 = unrelated. Also returns the RMSD of all pairs under the best superposition.
    """
    n = len(x)
    if n < 3 or length < 3:
        return {"tm": None, "rmsd": None, "aligned": n}
    d0 = max(0.5, 1.24 * (length - 15) ** (1.0 / 3.0) - 1.8) if length > 21 else 0.5
    best, best_rot = -1.0, None
    sizes = sorted({n, max(4, n // 2), max(4, n // 4), max(4, min(n, 12))}, reverse=True)
    for size in sizes:
        step = max(1, size // 2)
        for start in range(0, n - size + 1, step):
            idx = np.arange(start, start + size)
            for _ in range(20):
                rot, shift = _kabsch(x[idx], y[idx])
                dist = np.linalg.norm(x @ rot + shift - y, axis=1)
                score = float(np.sum(1.0 / (1.0 + (dist / d0) ** 2)) / length)
                if score > best:
                    best, best_rot = score, (rot, shift)
                cutoff = d0
                new = np.nonzero(dist < cutoff)[0]
                while len(new) < 3 and cutoff < 20:
                    cutoff += 0.5
                    new = np.nonzero(dist < cutoff)[0]
                if len(new) < 3 or np.array_equal(new, idx):
                    break
                idx = new
    rot, shift = best_rot
    rmsd = float(np.sqrt(np.mean(np.sum((x @ rot + shift - y) ** 2, axis=1))))
    return {"tm": round(best, 4), "rmsd": round(rmsd, 2), "aligned": n, "d0": round(d0, 2)}


def structure_similarity(ref_pdb: str, alt_pdb: str, ref_seq: str, alt_seq: str) -> Dict[str, Any]:
    """TM-score between two models of (windows of) the reference and the altered protein."""
    _ri, rx, _rb = ca_coordinates(ref_pdb)
    _ai, ax, _ab = ca_coordinates(alt_pdb)
    if len(rx) != len(ref_seq) or len(ax) != len(alt_seq):
        return {"tm": None, "rmsd": None, "aligned": 0, "error": "model and sequence lengths differ"}
    pairs = global_alignment(ref_seq, alt_seq)
    if len(pairs) < 3:
        return {"tm": None, "rmsd": None, "aligned": len(pairs)}
    ia = np.array([p[0] for p in pairs])
    ib = np.array([p[1] for p in pairs])
    return tm_score(rx[ia], ax[ib], len(ref_seq))


def residue_context(pdb: str, positions: Sequence[int], offset: int = 0) -> Dict[int, Dict[str, Any]]:
    """Per residue (our numbering): model confidence (pLDDT) and burial (CA neighbours within 10 Å)."""
    resi, xyz, b = ca_coordinates(pdb)
    index = {r: k for k, r in enumerate(resi)}
    out = {}
    for pos in positions:
        k = index.get(pos + offset)
        if k is None:
            continue
        neighbours = int(np.count_nonzero(np.linalg.norm(xyz - xyz[k], axis=1) < 10.0)) - 1
        out[pos] = {"plddt": round(float(b[k]), 1), "neighbours": neighbours,
                    "burial": "buried" if neighbours >= 16 else "intermediate" if neighbours >= 10 else "surface"}
    return out
