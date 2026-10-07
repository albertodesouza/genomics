"""A small client for a running visualizer, so notebooks stop re-implementing the loaders.

The visualizer already resolves windows, matches track metadata, remaps haplotype coordinates,
aggregates groups and caches the results. A notebook that opens the dataset directly repeats all of
that, usually a little differently, which is how a figure and the app end up disagreeing. This asks
the running app for the same arrays it is drawing:

    from genomics.visualizer.client import Visualizer

    v = Visualizer("http://127.0.0.1:8791", dataset="1kg_high_coverage")
    sig = v.signal("TYRP1", "rna_seq", series=["HG00096:H1+H2"], tracks=[0])
    sig["edges"], sig["series"][0]["mean"]      # numpy arrays, ready to plot

    means = v.group_means("TYRP1", "rna_seq", field="pigmentation", groups=["weak", "strong"])
    v.samples()                                 # a DataFrame when pandas is installed, else dicts

*Copy as Python* on the Tracks page writes the call for the view on screen.

Only the standard library and numpy are needed. Endpoints that compute in the background return
``{"pending": true}`` until done; every method here polls until the result arrives, printing
progress when ``progress=True`` (the default in a terminal is quiet).
"""
from __future__ import annotations

import base64
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np

REFERENCE_SAMPLE = "@reference"  # signals.REFERENCE_SAMPLE; spelled out so this module stays importable alone
DEFAULT_URL = "http://127.0.0.1:8791"
POLL_SECONDS = 0.5
MAX_POLL_SECONDS = 3.0
ARRAY_KEYS = ("edges", "mean", "min", "max", "upper", "lower", "values", "scores", "positions", "af")


class VisualizerError(RuntimeError):
    """The visualizer refused a request, or is not reachable."""


def _arrays(value: Any, key: str = "") -> Any:
    """Decode the payload's arrays into numpy.

    Float arrays travel as ``{"$f32": <base64 little-endian float32>, "shape": [...]}`` (see the
    server's JSONEncoder); integer arrays travel as plain lists, so the ones worth having as arrays
    are named in ``ARRAY_KEYS``.
    """
    if isinstance(value, dict):
        if isinstance(value.get("$f32"), str):
            flat = np.frombuffer(base64.b64decode(value["$f32"]), dtype="<f4")
            shape = value.get("shape") or [flat.size]
            return flat.reshape(tuple(int(n) for n in shape))
        return {k: _arrays(v, k) for k, v in value.items()}
    if isinstance(value, list):
        if key in ARRAY_KEYS and value and all(isinstance(x, (int, float, type(None))) for x in value):
            return np.asarray([np.nan if x is None else x for x in value], dtype=np.float64)
        return [_arrays(v) for v in value]
    return value


class Visualizer:
    """Read-only access to a running ``genomics visualize`` server.

    ``url`` is where the app is serving; ``dataset`` is the dataset id (the first one if omitted).
    Nothing here writes to the dataset: the client only reads what the app can draw.
    """

    def __init__(self, url: str = DEFAULT_URL, dataset: Optional[str] = None, timeout: float = 120.0, progress: bool = False):
        self.url = url.rstrip("/")
        self.timeout = timeout
        self.progress = progress
        self._status = self._request("/api/status")
        ids = [d["id"] for d in self._status.get("datasets", [])]
        if dataset is None:
            if not ids:
                raise VisualizerError(f"{self.url} has no dataset open")
            dataset = ids[0]
        elif dataset not in ids:
            raise VisualizerError(f"Unknown dataset {dataset!r}; this server has {', '.join(ids) or 'none'}")
        self.dataset = dataset

    def __repr__(self) -> str:
        return f"Visualizer({self.url!r}, dataset={self.dataset!r})"

    # -- plumbing ---------------------------------------------------------------------------
    def _request(self, path: str, params: Optional[Dict[str, Any]] = None) -> Any:
        query = urllib.parse.urlencode({k: v for k, v in (params or {}).items() if v is not None})
        url = f"{self.url}{path}{'?' + query if query else ''}"
        try:
            with urllib.request.urlopen(url, timeout=self.timeout) as response:  # noqa: S310 - a local app
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", "replace")
            try:
                detail = json.loads(detail).get("error", detail)
            except ValueError:
                pass
            raise VisualizerError(f"{path}: {detail}") from None
        except urllib.error.URLError as exc:
            raise VisualizerError(f"Cannot reach {self.url} ({exc.reason}). Is `genomics visualize` running?") from None

    def _get(self, path: str, params: Optional[Dict[str, Any]] = None) -> Any:
        """A request that waits out the background job the endpoint may start."""
        delay = POLL_SECONDS
        said = False
        while True:
            payload = self._request(path, params)
            if not (isinstance(payload, dict) and payload.get("pending")):
                if said and self.progress:
                    print("done")
                return payload
            job = payload.get("job") or {}
            if self.progress:
                print(f"\r{job.get('title', 'working')}: {job.get('message', '')} {round(100 * (job.get('progress') or 0))}%  ", end="")
                said = True
            time.sleep(delay)
            delay = min(MAX_POLL_SECONDS, delay * 1.3)

    def _ds(self, path: str) -> str:
        return f"/api/d/{urllib.parse.quote(self.dataset)}{path}"

    @staticmethod
    def _csv(values: Optional[Iterable[Any]]) -> Optional[str]:
        return None if values is None else ",".join(str(v) for v in values)

    # -- the dataset ------------------------------------------------------------------------
    def datasets(self) -> List[Dict[str, Any]]:
        """Every dataset this server has open."""
        return self._request("/api/status").get("datasets", [])

    def summary(self) -> Dict[str, Any]:
        """Windows, outputs, ontologies and sample fields of the dataset."""
        return self._get(self._ds("/summary"))

    def genes(self) -> List[str]:
        return [g["gene"] for g in self.summary().get("genes", [])]

    def samples(self, as_frame: bool = True) -> Any:
        """The sample table, including any region scalars and PC fields.

        A pandas DataFrame when pandas is installed and ``as_frame``, else ``{"columns", "rows"}``.
        """
        payload = self._get(self._ds("/samples"))
        if not as_frame:
            return payload
        try:
            import pandas as pd
        except ImportError:
            return payload
        return pd.DataFrame(payload["rows"], columns=payload["columns"])

    def cohort(self, filters: Optional[Dict[str, Sequence[str]]] = None) -> List[str]:
        """Sample ids matching ``filters`` ({field: [values]}) — the cohort the Samples page builds."""
        payload = self._get(self._ds("/samples"))
        index = {name: i for i, name in enumerate(payload["columns"])}
        wanted = {k: {str(x) for x in v} for k, v in (filters or {}).items() if k in index}
        out = []
        for row in payload["rows"]:
            if all(str(row[index[field]]) in values for field, values in wanted.items()):
                out.append(row[index["sample_id"]])
        return out

    # -- what the Tracks page draws ---------------------------------------------------------
    def signal(self, gene: str, output: str, series: Sequence[str], tracks: Sequence[int] = (0,),
               start: int = 0, end: Optional[int] = None, bins: int = 1000, coords: str = "reference") -> Dict[str, Any]:
        """Per-haplotype predicted signal, binned exactly as the Tracks page bins it.

        ``series`` are ``"<sample>:<haplotype>"`` (``H1``, ``H2`` or ``H1+H2``), ``tracks`` are track
        indices of ``output``. Returns ``{"edges", "tracks", "series": [{"label", "mean", "min", ...}]}``
        with numpy arrays; ``mean`` is (tracks, bins).
        """
        params = {"gene": gene, "output": output, "series": self._csv(series), "tracks": self._csv(tracks),
                  "coords": coords, "start": start, "bins": bins}
        params["end"] = self._window_end(gene, end)
        return _arrays(self._get(self._ds("/signal"), params))

    def group_means(self, gene: str, output: str, field: str, groups: Optional[Sequence[str]] = None,
                    tracks: Sequence[int] = (0,), haps: Sequence[str] = ("H1", "H2"), start: int = 0,
                    end: Optional[int] = None, bins: int = 1000, coords: str = "reference",
                    filters: Optional[Dict[str, Sequence[str]]] = None) -> Dict[str, Any]:
        """Mean ± SD per group of a sample field, over the cohort — the Tracks "group means" view.

        Computing these reads every haplotype of each group once, so the first call for a window can
        take minutes; the server caches it, and this polls until it is done.
        """
        params = {"gene": gene, "output": output, "field": field, "groups": self._csv(groups),
                  "tracks": self._csv(tracks), "haps": self._csv(haps), "coords": coords,
                  "start": start, "bins": bins, "filters": json.dumps(filters or {})}
        params["end"] = self._window_end(gene, end)
        return _arrays(self._get(self._ds("/groups"), params))

    def reference_signal(self, gene: str, output: str, tracks: Sequence[int] = (0,), start: int = 0,
                         end: Optional[int] = None, bins: int = 1000) -> Dict[str, Any]:
        """AlphaGenome's prediction of the reference window — the dashed baseline on Tracks."""
        return self.signal(gene, output, series=[f"{REFERENCE_SAMPLE}:H1"], tracks=tracks, start=start, end=end, bins=bins)

    def observed(self, gene: str, output: str, tracks: Sequence[int] = (0,), start: int = 0, end: Optional[int] = None,
                 bins: int = 1000, scale: str = "tpm") -> Dict[str, Any]:
        """The measured ENCODE / FANTOM5 signal behind tracks, in the source's own units."""
        params = {"gene": gene, "output": output, "tracks": self._csv(tracks), "start": start, "bins": bins, "scale": scale}
        params["end"] = self._window_end(gene, end)
        return _arrays(self._get(self._ds("/observed"), params))

    def _window_end(self, gene: str, end: Optional[int]) -> int:
        if end is not None:
            return end
        for window in self.summary().get("genes", []):
            if window["gene"] == gene and window.get("length"):
                return int(window["length"])
        raise VisualizerError(f"Give an end offset: the length of {gene} is unknown")

    # -- what the Variant page computes ------------------------------------------------------
    def variant_sites(self, gene: str, min_af: float = 0.05, start: Optional[int] = None,
                      end: Optional[int] = None, limit: int = 2000) -> List[Dict[str, Any]]:
        """Sites carried in the cohort in a window, with ALT frequency (genomic positions)."""
        params = {"gene": gene, "min_af": min_af, "limit": limit, "start": start, "end": end}
        return self._get(self._ds("/variant/sites"), params).get("sites", [])

    def variant_effect(self, gene: str, pos: int, ref: str, alt: str, output: str, track: int,
                       start: int, end: int, field: Optional[str] = None,
                       filters: Optional[Dict[str, Sequence[str]]] = None) -> Dict[str, Any]:
        """Predicted signal over a region by genotype at one site, with the slope per ALT allele.

        The in-silico eQTL of the Variant page: each sample's mean over ``[start, end)`` grouped by
        0/0, 0/1 and 1/1, plus the OLS slope, its standard error and p.
        """
        params = {"gene": gene, "pos": pos, "ref": ref, "alt": alt, "output": output, "track": track,
                  "start": start, "end": end, "field": field, "filters": json.dumps(filters or {})}
        return _arrays(self._get(self._ds("/variant/effect"), params))

    def genotypes(self, gene: str, pos: int, ref: str, alt: str, field: Optional[str] = None) -> Dict[str, Any]:
        """Genotype counts and ALT frequency at one site, overall and by a sample field."""
        params = {"gene": gene, "pos": pos, "ref": ref, "alt": alt, "field": field}
        return self._get(self._ds("/variant/site"), params)

    # -- gene products ---------------------------------------------------------------------
    def products(self, sample: str, gene: str, target: Optional[str] = None, sequences: bool = False) -> Dict[str, Any]:
        """One sample's mature transcripts and proteins per haplotype (the Gene products page).

        ``transcripts[i]["products"][hap]`` holds ``protein``, ``mrna_length``, ``nmd`` and
        ``change`` (class and HGVS ``p.`` against the reference) for ``ref``, ``H1`` and ``H2``;
        ``splicing`` holds AlphaGenome's junction usage per tissue track and the candidate isoforms.
        ``target`` picks another gene in the window; ``sequences`` adds each mRNA.
        """
        params = {"gene": gene, "sample": sample, "target": target, "sequences": "mrna" if sequences else None}
        return self._get(self._ds("/products"), params)

    def products_fasta(self, sample: str, gene: str, kind: str = "protein", target: Optional[str] = None) -> str:
        """FASTA of every product (``kind``: protein or mrna), headers ``gene|transcript|sample.H|change``."""
        from genomics.visualizer.products import fasta

        return fasta(self.products(sample, gene, target=target, sequences=kind == "mrna"), kind)

    def report(self, sample: str, gene: str, tissue: str = "CL:1000458", target: Optional[str] = None) -> Dict[str, Any]:
        """The gene report of the Gene products page for one tissue (any AlphaGenome ontology term with RNA-seq).

        ``baseline``: each transcript's share and absolute level for the reference genome (GTEx /
        HPA). ``haplotypes[H]``: the copy's level per transcript (``share``, ``level``, ``fold``),
        ``flags`` (``unstable``, ``premature_stop``, ``missense`` with AlphaMissense,
        ``base_similarity``), and every variant with its consequence (``variants``, ``counts``).
        """
        return self._get(self._ds("/products/report"), {"gene": gene, "sample": sample, "tissue": tissue, "target": target})

    def expression(self, sample: str, gene: str, target: Optional[str] = None) -> Dict[str, Any]:
        """mRNA and protein of a sample against the reference genome, per tissue (the Gene products page).

        ``tissues[i]["mrna"]["fold"]`` holds the relative change (``H1``, ``H2``, ``individual``),
        ``["mrna"]["absolute"]`` the estimated level in the anchor's unit (``anchor``), and
        ``["protein"]["fold"]`` the relative protein change; ``notes`` says what each one rests on.
        May run AlphaGenome for tissues the dataset has no stored RNA-seq of.
        """
        return self._get(self._ds("/products/expression"), {"gene": gene, "sample": sample, "target": target})

    # -- derived sample values ---------------------------------------------------------------
    def scalars(self) -> List[Dict[str, Any]]:
        """Region scalars defined on this dataset (their values are columns of ``samples()``)."""
        return self._get(self._ds("/scalars")).get("scalars", [])

    def pca(self, genes: Optional[Sequence[str]] = None, min_maf: float = 0.05, spacing: int = 2000,
            components: int = 10) -> Dict[str, Any]:
        """Genotype PCA of the cohort — the Samples page's Ancestry PCA.

        Returns ``{"samples", "scores" (n, k), "explained", "n_sites", ...}``. The first call for a
        window set reads every sample's VCF, which this polls through.
        """
        params = {"genes": self._csv(genes), "min_maf": min_maf, "spacing": spacing, "components": components}
        return _arrays(self._get(self._ds("/ancestry/pca"), params))

    def sessions(self) -> List[Dict[str, Any]]:
        """Saved views of this dataset (name, note, when); ``session(name)`` returns one's state."""
        return self._get(self._ds("/sessions")).get("sessions", [])

    def session(self, name: str) -> Dict[str, Any]:
        return self._get(self._ds("/sessions"), {"name": name})
