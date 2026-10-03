"""Single-process HTTP server for the genomics visualizer (stdlib only).

All pages are served by one static single-page app (``static/``) talking to a JSON API under
``/api``. Numeric arrays travel as base64 little-endian float32 (``{"$f32": ..., "shape": ...}``)
instead of JSON number lists, and responses are gzip-compressed when the client accepts it.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import math
import mimetypes
import re
import signal
import sys
import threading
import time
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, unquote, urlparse

import numpy as np

from genomics.visualizer.alignment import AlignmentService
from genomics.visualizer.alphagenome import create_backend
from genomics.visualizer.annotations import AnnotationService
from genomics.visualizer.cli_args import build_arg_parser
from genomics.visualizer.cache import DiskArrayCache, stable_key
from genomics.visualizer.datasets import DatasetCatalog
from genomics.visualizer.experiments import ExperimentService
from genomics.visualizer.jobs import JobManager
from genomics.visualizer.labs import LabsService
from genomics.visualizer.sequences import SequenceService
from genomics.visualizer.signals import COORDINATE_SYSTEMS, DIPLOID, MAX_GROUP_TRACKS, SignalService, parse_series
from genomics.visualizer.views import ViewService

STATIC_DIR = Path(__file__).resolve().parent / "static"
MAX_GROUPS = 12
MAX_POPULATION_ROWS = 4000
VERSION = "1.0"


class JSONEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Any:  # noqa: D401 - json API
        if isinstance(o, np.ndarray):
            if o.dtype.kind == "f":
                arr = np.ascontiguousarray(o, dtype="<f4")
                return {"$f32": base64.b64encode(arr.tobytes()).decode("ascii"), "shape": list(arr.shape)}
            if o.dtype.kind in "iub":
                return o.astype(np.int64).tolist()
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            value = float(o)
            return value if math.isfinite(value) else None
        if isinstance(o, Path):
            return str(o)
        if isinstance(o, set):
            return sorted(o)
        return super().default(o)


def encode_json(payload: Any) -> bytes:
    return json.dumps(to_jsonable_shallow(payload), cls=JSONEncoder, separators=(",", ":"), allow_nan=False).encode("utf-8")


def to_jsonable_shallow(value: Any) -> Any:
    """Replace non-finite Python floats (numpy arrays are handled by the encoder)."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): to_jsonable_shallow(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable_shallow(v) for v in value]
    return value


class HttpError(Exception):
    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status


class Query:
    def __init__(self, raw: str):
        self.params = parse_qs(raw, keep_blank_values=True)

    def str(self, name: str, default: Optional[str] = None, required: bool = False) -> str:
        values = self.params.get(name)
        if not values or values[0] == "":
            if required:
                raise HttpError(HTTPStatus.BAD_REQUEST, f"Missing parameter: {name}")
            return default or ""
        return values[0]

    def int(self, name: str, default: int, lo: Optional[int] = None, hi: Optional[int] = None) -> int:
        raw = self.str(name)
        try:
            value = int(float(raw)) if raw else default
        except ValueError:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Invalid integer for {name}: {raw}")
        if lo is not None:
            value = max(lo, value)
        if hi is not None:
            value = min(hi, value)
        return value

    def list(self, name: str) -> List[str]:
        return [item.strip() for raw in self.params.get(name, []) for item in raw.split(",") if item.strip()]

    def ints(self, name: str) -> List[int]:
        try:
            return [int(v) for v in self.list(name)]
        except ValueError:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Invalid integer list for {name}")

    def json(self, name: str) -> Any:
        raw = self.str(name)
        if not raw:
            return None
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Invalid JSON for {name}")


class VisualizerApp:
    def __init__(
        self,
        catalog: DatasetCatalog,
        cache_dir: Optional[Path],
        memory_bytes: int,
        workers: int,
        runs_roots: List[Path],
        model_window: int = 32768,
        gtf: Optional[Path] = None,
        consensus_dirs: Optional[Dict[str, str]] = None,
        labs: Optional[LabsService] = None,
        allow_add_datasets: bool = True,
    ):
        self.catalog = catalog
        self.cache_dir = cache_dir
        self.alignment = AlignmentService(model_window=model_window, cache_bytes=max(memory_bytes // 8, 64 << 20), consensus_dirs=consensus_dirs)
        self.signals = SignalService(memory_bytes, DiskArrayCache(cache_dir), alignment=self.alignment, workers=workers)
        self.sequences = SequenceService(self.signals)
        self.annotations = AnnotationService(cache_dir, gtf=gtf)
        self.experiments = ExperimentService(runs_roots)
        self.views = ViewService()
        self.labs = labs or LabsService(None)
        self.jobs = JobManager(workers=2)
        self.allow_add_datasets = allow_add_datasets
        self.started = time.time()
        self.routes: List[Tuple[str, re.Pattern, Callable[..., Any]]] = []
        self._register_routes()

    # -- routing ---------------------------------------------------------------------
    def route(self, method: str, pattern: str, handler: Callable[..., Any]) -> None:
        self.routes.append((method, re.compile("^" + pattern + "$"), handler))

    def _register_routes(self) -> None:
        r = self.route
        r("GET", r"/api/status", self.api_status)
        r("GET", r"/api/datasets", self.api_datasets)
        r("POST", r"/api/datasets", self.api_add_dataset)
        r("GET", r"/api/jobs", self.api_jobs)
        r("POST", r"/api/jobs/(?P<job_id>[^/]+)/cancel", self.api_cancel_job)
        r("GET", r"/api/d/(?P<ds>[^/]+)/summary", self.api_summary)
        r("GET", r"/api/d/(?P<ds>[^/]+)/samples", self.api_samples)
        r("GET", r"/api/d/(?P<ds>[^/]+)/samples/(?P<sample>[^/]+)", self.api_sample)
        r("GET", r"/api/d/(?P<ds>[^/]+)/genes/(?P<gene>[^/]+)", self.api_gene)
        r("GET", r"/api/d/(?P<ds>[^/]+)/genes/(?P<gene>[^/]+)/annotations", self.api_annotations)
        r("GET", r"/api/d/(?P<ds>[^/]+)/genes/(?P<gene>[^/]+)/axis", self.api_axis)
        r("GET", r"/api/d/(?P<ds>[^/]+)/signal", self.api_signal)
        r("GET", r"/api/d/(?P<ds>[^/]+)/groups", self.api_groups)
        r("GET", r"/api/d/(?P<ds>[^/]+)/population", self.api_population)
        r("GET", r"/api/d/(?P<ds>[^/]+)/sequence", self.api_sequence)
        r("POST", r"/api/d/(?P<ds>[^/]+)/views/preview", self.api_view_preview)
        r("POST", r"/api/d/(?P<ds>[^/]+)/views/save", self.api_view_save)
        r("GET", r"/api/runs", self.api_runs)
        r("GET", r"/api/runs/detail", self.api_run_detail)
        r("GET", r"/api/runs/file", self.api_run_file)
        r("GET", r"/api/labs", self.api_labs)
        r("POST", r"/api/labs/(?P<key>[^/]+)/start", self.api_lab_start)
        r("POST", r"/api/labs/(?P<key>[^/]+)/stop", self.api_lab_stop)
        r("GET", r"/api/alphagenome", self.api_alphagenome)
        r("POST", r"/api/alphagenome/settings", self.api_alphagenome_settings)
        r("POST", r"/api/alphagenome/test", self.api_alphagenome_test)
        r("POST", r"/api/alphagenome/local/start", self.api_alphagenome_local_start)
        r("POST", r"/api/alphagenome/local/stop", self.api_alphagenome_local_stop)

    def dispatch(self, method: str, path: str, query: Query, body: Any) -> Any:
        for route_method, pattern, handler in self.routes:
            if route_method != method:
                continue
            match = pattern.match(path)
            if match:
                kwargs = {k: unquote(v) for k, v in match.groupdict().items()}
                return handler(query=query, body=body, **kwargs)
        raise HttpError(HTTPStatus.NOT_FOUND, f"No route for {method} {path}")

    def dataset(self, ds: str):
        try:
            return self.catalog.get(ds)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    def _gene(self, dataset, gene: str) -> str:
        if gene not in dataset.genes:
            raise HttpError(HTTPStatus.NOT_FOUND, f"Unknown gene: {gene}")
        return gene

    @staticmethod
    def _coords(query: Query) -> str:
        coords = query.str("coords", "reference")
        if coords not in COORDINATE_SYSTEMS:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"coords must be one of {', '.join(COORDINATE_SYSTEMS)}")
        return coords

    def _job_response(self, key: str, title: str, fn: Callable[[Callable[[float, str], None]], Any]) -> Any:
        job = self.jobs.run(key, title, fn)
        if job.status == "done":
            return job.result
        if job.status == "error":
            raise HttpError(HTTPStatus.INTERNAL_SERVER_ERROR, job.error or "Job failed")
        if job.status == "cancelled":
            raise HttpError(HTTPStatus.CONFLICT, "Cancelled")
        return {"pending": True, "job": job.as_dict()}

    # -- general --------------------------------------------------------------------
    def api_status(self, query: Query, body: Any) -> Any:
        return {
            "version": VERSION,
            "uptime": round(time.time() - self.started, 1),
            "datasets": [d.listing() for d in self.catalog.all()],
            "alignment": self.alignment.available(),
            "caches": {
                "arrays": self.signals.arrays.stats(),
                "small": self.signals.small.stats(),
                "results": self.signals.results.stats(),
                "alignment_entries": self.alignment.entries.stats(),
            },
            "cache_dir": str(self.cache_dir) if self.cache_dir else None,
            "runs_roots": [str(r) for r in self.experiments.roots],
            "allow_add_datasets": self.allow_add_datasets,
            "jobs": self.jobs.active(),
        }

    def api_datasets(self, query: Query, body: Any) -> Any:
        return {"datasets": [d.listing() for d in self.catalog.all()]}

    def api_add_dataset(self, query: Query, body: Any) -> Any:
        if not self.allow_add_datasets:
            raise HttpError(HTTPStatus.FORBIDDEN, "Adding datasets is disabled")
        path = str((body or {}).get("path") or "").strip()
        if not path:
            raise HttpError(HTTPStatus.BAD_REQUEST, "Provide a dataset directory path")
        annotations = (body or {}).get("annotations") or None
        try:
            dataset = self.catalog.add(Path(path), annotations=Path(annotations) if annotations else None)
        except FileNotFoundError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        return dataset.listing()

    def api_jobs(self, query: Query, body: Any) -> Any:
        return {"jobs": self.jobs.active()}

    def api_cancel_job(self, query: Query, body: Any, job_id: str) -> Any:
        job = self.jobs.cancel(job_id)
        if job is None:
            raise HttpError(HTTPStatus.NOT_FOUND, f"Unknown job: {job_id}")
        return job.as_dict()

    # -- dataset ----------------------------------------------------------------------
    def api_summary(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        summary = dataset.summary()
        summary["annotations"] = self.annotations.status(dataset)
        return summary

    def api_samples(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        columns = [f["name"] for f in dataset.fields]
        return {
            "fields": dataset.fields,
            "columns": columns,
            "rows": [[row.get(c) for c in columns] for row in dataset.samples],
        }

    def api_sample(self, query: Query, body: Any, ds: str, sample: str) -> Any:
        dataset = self.dataset(ds)
        try:
            return dataset.sample_detail(sample)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    def api_gene(self, query: Query, body: Any, ds: str, gene: str) -> Any:
        dataset = self.dataset(ds)
        info = dict(dataset.gene_info(self._gene(dataset, gene)))
        info["annotations"] = self.annotations.status(dataset)
        info["alignment"] = self.alignment.available()
        return info

    def api_annotations(self, query: Query, body: Any, ds: str, gene: str) -> Any:
        dataset = self.dataset(ds)
        self._gene(dataset, gene)
        cached = self.annotations.cached(dataset)
        if cached is None:
            result = self._job_response(f"annotations:{dataset.id}", "Loading gene annotations", lambda progress: self.annotations.build(dataset, progress))
            if isinstance(result, dict) and result.get("pending"):
                return result
            cached = result
        return {"source": cached.get("source"), "genes": cached.get("genes", {}).get(gene, [])}

    def api_axis(self, query: Query, body: Any, ds: str, gene: str) -> Any:
        dataset = self.dataset(ds)
        self._gene(dataset, gene)

        def build(progress):
            progress(0.05, "Building training alignment axis (first use can take minutes)")
            axis = self.alignment.axis(dataset, gene)
            window = dataset.window(gene)
            return {
                "expanded_length": axis["expanded_length"],
                "ref_start_offset": axis["ref_start_offset"],
                "ref_length": axis["ref_length"],
                "genomic_start": (window.start or 0) + axis["ref_start_offset"],
                "insertion_slots": axis["insertion_slots"].tolist(),
                "model_window": axis["model_window"],
                "sample_set_key": axis["sample_set_key"],
            }

        if self.alignment.has_axis(dataset, gene):
            return build(lambda *_: None)
        return self._job_response(f"axis:{dataset.id}:{gene}", f"Alignment axis for {gene}", build)

    # -- signals ------------------------------------------------------------------------
    def _common_signal_params(self, dataset, query: Query) -> Tuple[str, str, str, int, int, int]:
        gene = self._gene(dataset, query.str("gene", required=True))
        output = query.str("output", required=True)
        coords = self._coords(query)
        start = query.int("start", 0, lo=0)
        end = query.int("end", start + 1, lo=start + 1)
        bins = query.int("bins", 1000, lo=1, hi=8192)
        return gene, output, coords, start, end, bins

    def api_signal(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        gene, output, coords, start, end, bins = self._common_signal_params(dataset, query)
        series = parse_series(query.str("series", required=True))
        tracks = query.ints("tracks") or [0]
        try:
            return self.signals.series_payload(dataset, gene, output, series, tracks, coords, start, end, bins)
        except (FileNotFoundError, KeyError) as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        except (ValueError, IndexError) as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def _cohort(self, dataset, query: Query) -> List[str]:
        filters = query.json("filters") or {}
        if not isinstance(filters, dict):
            raise HttpError(HTTPStatus.BAD_REQUEST, "filters must be a JSON object")
        return dataset.filter_samples({str(k): [str(x) for x in (v or [])] for k, v in filters.items()})

    def api_groups(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        gene, output, coords, start, end, bins = self._common_signal_params(dataset, query)
        tracks = query.ints("tracks") or [0]
        if len(tracks) > MAX_GROUP_TRACKS:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"At most {MAX_GROUP_TRACKS} tracks for group aggregates")
        haps = query.list("haps") or ["H1", "H2"]
        haps = sorted({h for hap in haps for h in (("H1", "H2") if hap == DIPLOID else (hap,))})
        field = query.str("field", required=True)
        cohort = self._cohort(dataset, query)
        values = query.list("groups")
        if not values:
            counts: Dict[str, int] = {}
            for sample in cohort:
                key = str(dataset.samples[dataset.sample_index[sample]].get(field, ""))
                if key:
                    counts[key] = counts.get(key, 0) + 1
            values = sorted(counts, key=lambda k: (-counts[k], k))[:MAX_GROUPS]
        values = values[:MAX_GROUPS]
        groups = dataset.group_samples(field, values, cohort)
        groups = {k: v for k, v in groups.items() if v}
        if not groups:
            raise HttpError(HTTPStatus.BAD_REQUEST, "No samples in the selected groups")
        keys = {name: self.signals.group_aggregate_key(dataset, gene, output, tracks, haps, coords, members) for name, members in groups.items()}
        cached = {name: self.signals.cached_group_aggregate(key) for name, key in keys.items()}
        if any(v is None for v in cached.values()):
            missing = [name for name, v in cached.items() if v is None]
            total = sum(len(groups[n]) for n in missing)

            def compute(progress):
                done = 0
                for name in missing:
                    share = len(groups[name]) / max(total, 1)
                    scaled = progress.sub(done / max(total, 1), share, f"{name}: ")
                    # compute_group_aggregate returns immediately for groups cached meanwhile.
                    self.signals.compute_group_aggregate(dataset, gene, output, tracks, haps, coords, groups[name], scaled)
                    done += len(groups[name])
                return True

            # Keyed by the whole request (not just the missing groups) so polls keep hitting the
            # same job while it fills the cache group by group.
            job_key = "groups:" + stable_key(sorted(keys.values()))
            result = self._job_response(job_key, f"Group means for {gene} by {field}", compute)
            if isinstance(result, dict) and result.get("pending"):
                return result
            cached = {name: self.signals.cached_group_aggregate(key) for name, key in keys.items()}
        domain = self.signals.domain_length(dataset, gene, output, coords)
        payload = self.signals.group_payload({k: v for k, v in cached.items() if v is not None}, {k: len(v) for k, v in groups.items()}, tracks, domain, start, end, bins)
        payload.update({"gene": gene, "output": output, "coords": coords, "field": field, "haplotypes": haps})
        return payload

    def api_population(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        gene, output, coords, start, end, bins = self._common_signal_params(dataset, query)
        track = query.int("track", 0, lo=0)
        hap = query.str("hap", DIPLOID)
        field = query.str("field", "")
        max_rows = query.int("max_rows", 600, lo=1, hi=MAX_POPULATION_ROWS)
        cohort = self._cohort(dataset, query)
        if field:
            cohort.sort(key=lambda s: (str(dataset.samples[dataset.sample_index[s]].get(field, "")), s))
        if len(cohort) > max_rows:
            picks = np.linspace(0, len(cohort) - 1, max_rows).round().astype(int)
            cohort = [cohort[i] for i in sorted(set(picks.tolist()))]
        if not cohort:
            raise HttpError(HTTPStatus.BAD_REQUEST, "The cohort is empty")
        labels = [str(dataset.samples[dataset.sample_index[s]].get(field, "")) if field else "" for s in cohort]
        key = "population:" + stable_key([dataset.fingerprint, gene, output, track, hap, coords, cohort, start, end, bins])

        def compute(progress):
            payload = dict(self.signals.population_matrix(dataset, gene, output, track, hap, coords, cohort, start, end, bins, progress))
            payload.update({"labels": labels, "field": field, "track": track, "haplotype": hap, "coords": coords, "gene": gene, "output": output})
            return payload

        return self._job_response(key, f"Population heatmap for {gene}", compute)

    def api_sequence(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        gene = self._gene(dataset, query.str("gene", required=True))
        coords = self._coords(query)
        start = query.int("start", 0, lo=0)
        end = query.int("end", start + 1, lo=start + 1)
        bins = query.int("bins", 1000, lo=1, hi=4096)
        rows = parse_series(query.str("rows", ""))
        try:
            return self.sequences.window(dataset, gene, rows, coords, start, end, bins)
        except FileNotFoundError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    # -- views ---------------------------------------------------------------------------
    def api_view_preview(self, query: Query, body: Any, ds: str) -> Any:
        try:
            return self.views.preview(self.dataset(ds), body or {})
        except ValueError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def api_view_save(self, query: Query, body: Any, ds: str) -> Any:
        try:
            return self.views.save(self.dataset(ds), body or {})
        except FileExistsError as exc:
            raise HttpError(HTTPStatus.CONFLICT, str(exc))
        except ValueError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    # -- experiments ------------------------------------------------------------------------
    def api_runs(self, query: Query, body: Any) -> Any:
        return self.experiments.list()

    def api_run_detail(self, query: Query, body: Any) -> Any:
        try:
            return self.experiments.detail(query.str("id", required=True))
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    def api_run_file(self, query: Query, body: Any) -> Any:
        try:
            data, content_type = self.experiments.file(query.str("id", required=True), query.str("path", required=True))
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        return RawResponse(data, content_type)

    # -- labs --------------------------------------------------------------------------------
    def api_labs(self, query: Query, body: Any) -> Any:
        return {"labs": self.labs.list()}

    def api_lab_start(self, query: Query, body: Any, key: str) -> Any:
        try:
            return self.labs.start(key)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        except RuntimeError as exc:
            raise HttpError(HTTPStatus.CONFLICT, str(exc))

    def api_lab_stop(self, query: Query, body: Any, key: str) -> Any:
        try:
            return self.labs.stop(key)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    # -- AlphaGenome backend -------------------------------------------------------------------
    def _alphagenome(self):
        if self.labs.alphagenome is None:
            raise HttpError(HTTPStatus.NOT_FOUND, "AlphaGenome backend settings are not available")
        return self.labs.alphagenome

    def api_alphagenome(self, query: Query, body: Any) -> Any:
        return self._alphagenome().describe()

    def api_alphagenome_settings(self, query: Query, body: Any) -> Any:
        try:
            return self._alphagenome().update(body or {})
        except ValueError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def api_alphagenome_test(self, query: Query, body: Any) -> Any:
        return self._alphagenome().test(body or {})

    def api_alphagenome_local_start(self, query: Query, body: Any) -> Any:
        try:
            self._alphagenome().local.start()
        except RuntimeError as exc:
            raise HttpError(HTTPStatus.CONFLICT, str(exc))
        return self._alphagenome().describe()

    def api_alphagenome_local_stop(self, query: Query, body: Any) -> Any:
        self._alphagenome().local.stop()
        return self._alphagenome().describe()

    def shutdown(self) -> None:
        self.jobs.shutdown()
        self.labs.stop_all()


class RawResponse:
    def __init__(self, data: bytes, content_type: str):
        self.data = data
        self.content_type = content_type


class Handler(BaseHTTPRequestHandler):
    app: VisualizerApp
    protocol_version = "HTTP/1.1"
    verbose = False

    def log_message(self, fmt: str, *args: Any) -> None:
        if self.verbose:
            sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    def _send(self, status: int, body: bytes, content_type: str, extra: Optional[Dict[str, str]] = None) -> None:
        accepts_gzip = "gzip" in (self.headers.get("Accept-Encoding") or "")
        compressible = content_type.startswith(("application/json", "text/", "application/javascript"))
        headers = dict(extra or {})
        if accepts_gzip and compressible and len(body) > 1400:
            body = gzip.compress(body, compresslevel=5)
            headers["Content-Encoding"] = "gzip"
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        for key, value in headers.items():
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _send_json(self, payload: Any, status: int = 200) -> None:
        self._send(status, encode_json(payload), "application/json; charset=utf-8", {"Cache-Control": "no-store"})

    def _serve_static(self, rel: str) -> None:
        rel = rel.lstrip("/") or "index.html"
        target = (STATIC_DIR / rel).resolve()
        if STATIC_DIR not in target.parents or not target.is_file():
            target = STATIC_DIR / "index.html"
        stat = target.stat()
        etag = f'"{stat.st_mtime_ns:x}-{stat.st_size:x}"'
        if self.headers.get("If-None-Match") == etag:
            self.send_response(HTTPStatus.NOT_MODIFIED)
            self.send_header("ETag", etag)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        content_type = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        if target.suffix == ".js":
            content_type = "application/javascript"
        if content_type.startswith("text/") or content_type == "application/javascript":
            content_type += "; charset=utf-8"
        self._send(200, target.read_bytes(), content_type, {"ETag": etag, "Cache-Control": "no-cache"})

    def _handle(self, method: str) -> None:
        parsed = urlparse(self.path)
        if not parsed.path.startswith("/api/"):
            if method != "GET" and method != "HEAD":
                self._send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)
                return
            # Unknown paths fall back to index.html so client-side routes survive a reload.
            rel = parsed.path[len("/static/"):] if parsed.path.startswith("/static/") else ""
            self._serve_static(rel)
            return
        body = None
        if method == "POST":
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b""
            try:
                body = json.loads(raw.decode("utf-8")) if raw else {}
            except json.JSONDecodeError:
                self._send_json({"error": "Invalid JSON body"}, HTTPStatus.BAD_REQUEST)
                return
        try:
            result = self.app.dispatch(method, parsed.path, Query(parsed.query), body)
            if isinstance(result, RawResponse):
                self._send(200, result.data, result.content_type, {"Cache-Control": "max-age=60"})
            else:
                self._send_json(result)
        except HttpError as exc:
            self._send_json({"error": str(exc)}, exc.status)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as exc:  # unexpected: report without killing the server
            import traceback

            traceback.print_exc()
            try:
                self._send_json({"error": f"{type(exc).__name__}: {exc}"}, HTTPStatus.INTERNAL_SERVER_ERROR)
            except (BrokenPipeError, ConnectionResetError):
                pass

    def do_GET(self) -> None:
        self._handle("GET")

    def do_HEAD(self) -> None:
        self._handle("GET")

    def do_POST(self) -> None:
        self._handle("POST")


class Server(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True


def create_app(args: argparse.Namespace) -> VisualizerApp:
    from genomics.workspace import DEFAULT_DATASET_DIR, DEFAULT_GENOTYPE_RUNS_ROOT, cache_path

    catalog = DatasetCatalog()
    for dataset_id in args.dataset_id:
        from genomics.core.data_registry import resolve_dataset

        ref = resolve_dataset(dataset_id)
        catalog.add(ref.path, dataset_id=dataset_id, annotations=args.annotations)
    for path in args.dataset:
        catalog.add(path, annotations=args.annotations)
    if not catalog.all() and Path(DEFAULT_DATASET_DIR, "dataset_metadata.json").exists():
        catalog.add(DEFAULT_DATASET_DIR, annotations=args.annotations)
    if not catalog.all():
        print("No dataset found; add one with --dataset PATH or from the UI.", file=sys.stderr)

    runs_roots = list(args.runs_root) or ([DEFAULT_GENOTYPE_RUNS_ROOT] if Path(DEFAULT_GENOTYPE_RUNS_ROOT).exists() else [])
    cache_dir = None if args.no_disk_cache else Path(args.cache_dir or cache_path("visualizer")).resolve()
    consensus_dirs = {}
    if args.consensus_dataset_dir:
        consensus_dirs = {d.id: str(Path(args.consensus_dataset_dir).resolve()) for d in catalog.all()}
    primary = catalog.all()[0] if catalog.all() else None
    log_dir = (cache_dir / "logs") if cache_dir else None
    labs = LabsService(
        primary.path if primary else None,
        consensus_dir=Path(args.consensus_dataset_dir) if args.consensus_dataset_dir else None,
        pigmentation_config=args.pigmentation_config,
        port=args.lab_port,
        log_dir=log_dir,
        alphagenome=create_backend(args, log_dir),
    )
    return VisualizerApp(
        catalog,
        cache_dir=cache_dir,
        memory_bytes=max(args.memory_mb, 256) << 20,
        workers=max(1, args.workers),
        runs_roots=runs_roots,
        model_window=args.model_window,
        gtf=args.gtf,
        consensus_dirs=consensus_dirs,
        labs=labs,
        allow_add_datasets=not args.no_add_datasets,
    )


def serve(app: VisualizerApp, host: str, port: int, open_browser: bool = False, verbose: bool = False) -> None:
    handler = type("VisualizerHandler", (Handler,), {"app": app, "verbose": verbose})
    server = Server((host, port), handler)
    url = f"http://{host}:{port}/"
    print(f"Genomics Visualizer: {url}")
    for dataset in app.catalog.all():
        print(f"  dataset {dataset.id}: {dataset.path} ({len(dataset.samples)} samples, {len(dataset.genes)} genes)")
    if app.cache_dir:
        print(f"  cache: {app.cache_dir}")

    def stop(_signum=None, _frame=None):
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGTERM, stop)
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:
        pass
    finally:
        app.shutdown()
        server.server_close()


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    app = create_app(args)
    serve(app, args.host, args.port, open_browser=args.open, verbose=args.verbose)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
