"""Observed (experimental) signal for AlphaGenome tracks: the public data each track was trained on.

AlphaGenome's tracks come from ENCODE (RNA-seq, DNase, ATAC, histone / TF ChIP-seq, PRO-cap),
FANTOM5 (CAGE), GTEx (RNA-seq) and 4D Nucleome (contact maps). A track's metadata carries its
ontology term (CL / UBERON / EFO), assay, strand and data source, which is enough to find the
matching experiments again:

``fantom``  FANTOM5 CAGE libraries whose sample derives from the same ontology term, read from the
            FANTOM5 ontology (``ff-phase2-170801.obo``) and the hg38 data hub (``trackDb.txt``;
            CTSS read counts or TPM bigWigs per strand). Cell lines fall back to their name.
``encode``  ENCODE portal experiments with the same biosample term, assay and target; one GRCh38
            bigWig per experiment (the strand's signal of unique reads for RNA-seq, read-depth
            normalized signal for DNase, fold change over control for ChIP / ATAC, ...).

Values are read per base over the window with :mod:`genomics.visualizer.bigwig` (HTTP range
requests, no download of whole files), averaged over at most ``MAX_SOURCES`` experiments, and
cached on disk (``<cache>/observed``). Observed data is in reference coordinates and in the units
of the source (TPM, RPM-like signal, fold change), not in AlphaGenome's units: compare shapes.
"""
from __future__ import annotations

import re
import threading
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from genomics.visualizer.bigwig import BigWigError, open_bigwig
from genomics.visualizer.cache import DiskArrayCache, LRUCache, stable_key
from genomics.visualizer.remote import DAY, RemoteCache, RemoteError

FANTOM_OBO = "https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Ontology/ff-phase2-170801.obo.txt"
FANTOM_TRACKDB = "https://fantom.gsc.riken.jp/5/datahub/hg38/trackDb.txt"
FANTOM_SSTAR = "https://fantom.gsc.riken.jp/5/sstar/FF:{ff}"
ENCODE = "https://www.encodeproject.org"
MAX_SOURCES = 4
SCALES = ("tpm", "counts")
# ENCODE bigWig output types, best first (strand words are filtered separately).
ENCODE_OUTPUT_RANK = (
    "read-depth normalized signal",
    "signal of unique reads",
    "fold change over control",
    "signal of all reads",
    "signal p-value",
    "signal",
    "raw signal",
)
# Library names of treated / fractionated / time-course samples are ranked after plain ones.
_PERTURBED = re.compile(r"response|treat|infect|fraction|stimul|differentiat|\b\d+\s*(?:hr|h|min|day|d)\b|revived|SLAM", re.I)
ProgressFn = Callable[[float, str], None]


class ObservedUnavailable(Exception):
    pass


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (text or "").lower()).strip()


# ----------------------------------------------------------------------------------- FANTOM5
class Fantom5Index:
    """FANTOM5 hg38 CAGE libraries (bigWig URLs per strand / scale) and the ontology terms they derive from."""

    def __init__(self, obo_text: str, trackdb_text: str):
        classes: Dict[str, Dict[str, Any]] = {}
        library_parents: Dict[str, List[str]] = {}
        for block in obo_text.split("\n[Term]"):
            m = re.search(r"^id: (FF:\S+)", block, re.M)
            if not m:
                continue
            ff = m.group(1)
            parents = re.findall(r"^is_a: (FF:\S+)", block, re.M)
            derives = set(re.findall(r"^(?:relationship|intersection_of): derives_from (\S+)", block, re.M))
            if re.match(r"FF:\d+-\d+[A-Z]\d+$", ff):  # a library sample (e.g. FF:11274-116H5)
                library_parents[ff[3:]] = parents
            else:
                classes[ff] = {"derives": derives}
        self.libraries: List[Dict[str, Any]] = []
        by_ff: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for block in re.split(r"\n\s*\n", trackdb_text):
            url = re.search(r"^\s*bigDataUrl (\S+)", block, re.M)
            meta = re.search(r"^\s*metadata ontology_id=(\S+) sequence_tech=(\S+)", block, re.M)
            if not url or not meta:
                continue
            href = url.group(1).replace("http://", "https://", 1)
            scale = "counts" if "/ctss/" in href else "tpm" if "/tpm/" in href else None
            strand = "+" if href.endswith(".fwd.bw") else "-" if href.endswith(".rev.bw") else None
            if scale is None or strand is None:
                continue
            ff, tech = meta.group(1), meta.group(2)
            label = re.search(r"^\s*longLabel (.*)$", block, re.M)
            category = re.search(r"category=(\S+)", block)
            name = urllib.parse.unquote(href.rsplit("/", 1)[-1]).split(".CNhs")[0]
            cnhs = re.search(r"(CNhs\d+)", href)
            lib = by_ff.get((ff, tech))
            if lib is None:
                curies = set()
                for parent in library_parents.get(ff, []):
                    curies |= classes.get(parent, {}).get("derives", set())
                lib = {"ff": ff, "library": cnhs.group(1) if cnhs else "", "name": name or (label.group(1) if label else ff), "tech": tech,
                       "category": category.group(1) if category else "", "urls": {}, "curies": sorted(curies)}
                by_ff[(ff, tech)] = lib
                self.libraries.append(lib)
            lib["urls"].setdefault(scale, {})[strand] = href
        self.by_curie: Dict[str, List[Dict[str, Any]]] = {}
        for lib in self.libraries:
            for curie in lib["curies"]:
                self.by_curie.setdefault(curie, []).append(lib)

    def match(self, curie: str, name: str = "", cell_line: bool = False) -> Tuple[List[Dict[str, Any]], str]:
        """Libraries for an ontology term (``ontology``), or by cell-line name (``name``)."""
        libs = self.by_curie.get(curie or "", [])
        if libs:
            return self._rank(libs), "ontology"
        if cell_line and name and len(name) >= 3:
            pattern = re.compile(r"(?<![A-Za-z0-9])" + re.escape(name).replace(r"\-", r"[- ]?") + r"(?![A-Za-z0-9])", re.I)
            libs = [lib for lib in self.libraries if pattern.search(lib["name"])]
            if libs:
                return self._rank(libs), "name"
        return [], ""

    @staticmethod
    def _rank(libs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return sorted(libs, key=lambda lib: (bool(_PERTURBED.search(lib["name"])), lib["tech"] != "hCAGE", lib["name"]))


# ------------------------------------------------------------------------------------ ENCODE
def encode_search_params(meta: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    curie = meta.get("ontology_curie")
    assay = meta.get("Assay title")
    if not curie or not assay:
        return None
    params: Dict[str, Any] = {
        "type": "Experiment",
        "status": "released",
        "biosample_ontology.term_id": curie,
        "assay_title": assay,
        "replicates.library.biosample.donor.organism.scientific_name": "Homo sapiens",
    }
    target = meta.get("transcription_factor") or meta.get("histone_mark")
    if target:
        params["target.label"] = target
    return params


def encode_search_url(meta: Dict[str, Any]) -> Optional[str]:
    """The ENCODE portal search page listing every experiment behind a track (for people)."""
    params = encode_search_params(meta)
    return f"{ENCODE}/search/?{urllib.parse.urlencode(params)}" if params else None


def pick_encode_file(files: List[Dict[str, Any]], strand: Optional[str]) -> Optional[Dict[str, Any]]:
    """The GRCh38 bigWig of an experiment that matches the track's strand, best output type first."""
    def strand_ok(output_type: str) -> bool:
        has_plus, has_minus = "plus strand" in output_type, "minus strand" in output_type
        if strand == "+":
            return has_plus
        if strand == "-":
            return has_minus
        return not has_plus and not has_minus

    def rank(output_type: str) -> int:
        bare = output_type.replace("plus strand ", "").replace("minus strand ", "")
        return ENCODE_OUTPUT_RANK.index(bare) if bare in ENCODE_OUTPUT_RANK else len(ENCODE_OUTPUT_RANK)

    candidates = [f for f in files if f.get("file_format") == "bigWig" and f.get("assembly") == "GRCh38" and f.get("status") == "released"
                  and strand_ok(str(f.get("output_type") or "")) and f.get("href")]
    if not candidates and strand in ("+", "-"):  # unstranded library for a stranded track
        return pick_encode_file(files, None)
    candidates.sort(key=lambda f: (rank(str(f.get("output_type") or "")), not f.get("preferred_default"), -len(f.get("biological_replicates") or [])))
    return candidates[0] if candidates else None


# ----------------------------------------------------------------------------------- service
class ObservedService:
    def __init__(self, remote: RemoteCache, disk: DiskArrayCache, cache_bytes: int = 512 << 20):
        self.remote = remote
        self.disk = disk
        self.arrays = LRUCache(cache_bytes)
        self._fantom: Optional[Fantom5Index] = None
        self._lock = threading.Lock()

    # -- sources --------------------------------------------------------------------------------
    def fantom(self) -> Fantom5Index:
        with self._lock:
            if self._fantom is None:
                obo = self.remote.get_text(FANTOM_OBO, ttl=365 * DAY)
                trackdb = self.remote.get_text(FANTOM_TRACKDB, ttl=180 * DAY)
                self._fantom = Fantom5Index(obo, trackdb)
            return self._fantom

    def sources(self, output: str, meta: Dict[str, Any], scale: str = "tpm") -> Dict[str, Any]:
        """Where the observed signal of one track comes from: provider, experiments, files, links."""
        output = output.lower()
        source = str(meta.get("data_source") or ("fantom" if output == "cage" else "encode" if meta.get("Assay title") else "")).lower()
        strand = meta.get("strand") if meta.get("strand") in ("+", "-") else None
        curie = meta.get("ontology_curie") or ""
        if output in ("splice_sites", "splice_site_usage", "splice_junctions", "contact_maps"):
            raise ObservedUnavailable("No per-base observed track for this output (splicing / 3D outputs are derived from RNA-seq or Hi-C reads)")
        if source == "fantom" or output == "cage":
            if not strand:
                raise ObservedUnavailable("CAGE track without a strand")
            index = self.fantom()
            cell_line = meta.get("biosample_type") == "cell_line" or curie.startswith("EFO:")
            libs, how = index.match(curie, str(meta.get("biosample_name") or ""), cell_line)
            if not libs:
                raise ObservedUnavailable(f"No FANTOM5 CAGE library derives from {curie or 'this sample'}")
            scale = scale if scale in SCALES else "tpm"
            chosen = [lib for lib in libs if strand in lib["urls"].get(scale, {})][:MAX_SOURCES]
            return {
                "provider": "FANTOM5",
                "match": how,
                "units": "TPM (CTSS)" if scale == "tpm" else "CTSS read counts",
                "total": len(libs),
                "query": None,
                "items": [{
                    "id": lib["library"] or lib["ff"],
                    "label": lib["name"],
                    "detail": f"{lib['tech']} · FF:{lib['ff']}{' · ' + lib['category'] if lib['category'] else ''}",
                    "url": FANTOM_SSTAR.format(ff=lib["ff"]),
                    "file": lib["urls"][scale][strand],
                    "strand": strand,
                } for lib in chosen],
            }
        if source == "gtex":
            raise ObservedUnavailable("GTEx RNA-seq has no public per-base coverage track (see the GTEx portal for expression)")
        if source in ("encode", ""):
            params = encode_search_params(meta)
            if params is None:
                raise ObservedUnavailable("Track metadata has no ontology term / assay to look up in ENCODE")
            fields = ["accession", "description", "biosample_summary", "assay_title", "assay_term_id", "target.label", "lab.title", "date_released",
                      "files.accession", "files.href", "files.output_type", "files.file_format", "files.assembly", "files.status",
                      "files.preferred_default", "files.biological_replicates"]
            data = self.remote.get_json(f"{ENCODE}/search/", ttl=30 * DAY, params={**params, "format": "json", "limit": 25, "field": fields})
            experiments = data.get("@graph") or []
            items = []
            for exp in sorted(experiments, key=lambda e: str(e.get("date_released") or ""), reverse=True):
                f = pick_encode_file(exp.get("files") or [], strand)
                if f is None:
                    continue
                items.append({
                    "id": exp.get("accession"),
                    "label": exp.get("biosample_summary") or exp.get("description") or exp.get("accession"),
                    "detail": f"{f.get('output_type')} · {f.get('accession')}{' · ' + (exp.get('lab') or {}).get('title', '') if exp.get('lab') else ''}",
                    "url": f"{ENCODE}/experiments/{exp.get('accession')}/",
                    "file": f"{ENCODE}{f['href']}",
                    "file_url": f"{ENCODE}/files/{f.get('accession')}/",
                    "output_type": f.get("output_type"),
                    "assay_term_id": exp.get("assay_term_id"),
                })
            if not items:
                raise ObservedUnavailable(f"No released ENCODE GRCh38 bigWig for {params['assay_title']} in {curie}" + (f" ({params['target.label']})" if params.get("target.label") else ""))
            units = sorted({i["output_type"] for i in items[:MAX_SOURCES]})
            return {
                "provider": "ENCODE",
                "match": "ontology",
                "units": units[0] if len(units) == 1 else "mixed: " + ", ".join(units),
                "total": len(items),
                "query": encode_search_url(meta),
                "items": items[:MAX_SOURCES],
            }
        raise ObservedUnavailable(f"No observed-data lookup for data source {source!r}")

    # -- values ---------------------------------------------------------------------------------
    @staticmethod
    def _region(window) -> Tuple[str, int, int]:
        if window.chromosome is None or window.start is None or not window.length:
            raise ObservedUnavailable("The window has no genomic coordinates")
        start0 = int(window.start) - 1  # window.start is the 1-based position of offset 0
        return str(window.chromosome), start0, start0 + int(window.length)

    def _file_values(self, url: str, chrom: str, start: int, end: int) -> np.ndarray:
        key = stable_key(["observed", 1, url, chrom, start, end])
        return self.arrays.get_or_load(("file", key), lambda: self._load_file(key, url, chrom, start, end))

    def _load_file(self, key: str, url: str, chrom: str, start: int, end: int) -> np.ndarray:
        hit = self.disk.load("observed", key)
        if hit is not None and "values" in hit:
            return hit["values"].astype(np.float32, copy=False)
        values = open_bigwig(url).values(chrom, start, end)
        self.disk.save("observed", key, {"values": values}, compress=True)
        return values

    def _key(self, dataset, gene: str, output: str, track: int, scale: str) -> Tuple:
        return ("track", dataset.id, gene, output, int(track), scale if output.lower() == "cage" else "")

    def cached(self, dataset, gene: str, output: str, track: int, scale: str) -> Optional[Dict[str, Any]]:
        return self.arrays.get(self._key(dataset, gene, output, track, scale))

    def load(self, dataset, gene: str, output: str, track: int, meta: Dict[str, Any], scale: str, progress: Optional[ProgressFn] = None) -> Dict[str, Any]:
        """Observed signal of one track over the whole window: ``{"info": {...}, "values": (L,) or None}``."""
        key = self._key(dataset, gene, output, track, scale)
        hit = self.arrays.get(key)
        if hit is not None:
            return hit
        try:
            info = self.sources(output, meta, scale)
        except ObservedUnavailable as exc:
            result = {"info": {"available": False, "reason": str(exc)}, "values": None}
            self.arrays.put(key, result)
            return result
        except RemoteError as exc:  # not cached: retried next time
            return {"info": {"available": False, "reason": f"Lookup failed: {exc}"}, "values": None}
        chrom, start, end = self._region(dataset.window(gene))
        items = info["items"]
        columns: List[Optional[np.ndarray]] = [None] * len(items)
        errors: List[str] = []
        done = [0]

        def fetch(i: int) -> None:
            try:
                columns[i] = self._file_values(items[i]["file"], chrom, start, end)
            except (BigWigError, OSError, ValueError) as exc:
                errors.append(f"{items[i]['id']}: {exc}")
            done[0] += 1
            if progress is not None:
                progress(done[0] / max(1, len(items)), f"{info['provider']} {items[i]['id']} ({done[0]}/{len(items)})")

        with ThreadPoolExecutor(max_workers=min(4, max(1, len(items)))) as pool:
            list(pool.map(fetch, range(len(items))))
        ok = [c for c in columns if c is not None]
        for item, column in zip(items, columns):
            item["loaded"] = column is not None
        if not ok:
            return {"info": {"available": False, "reason": "; ".join(errors) or "No observed data could be read"}, "values": None}
        stacked = np.stack(ok, axis=1)
        values = np.abs(stacked).mean(axis=1).astype(np.float32)  # FANTOM5 reverse-strand files are negative
        info = {**info, "available": True, "errors": errors, "sources_loaded": len(ok), "region": f"{chrom}:{start + 1}-{end}"}
        result = {"info": info, "values": values}
        self.arrays.put(key, result)
        return result
