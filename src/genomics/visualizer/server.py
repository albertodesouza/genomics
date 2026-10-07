"""Single-process HTTP server for the genomics visualizer (stdlib only).

All pages are served by one static single-page app (``static/``) talking to a JSON API under
``/api``. Numeric arrays travel as base64 little-endian float32 (``{"$f32": ..., "shape": ...}``)
instead of JSON number lists, and responses are gzip-compressed when the client accepts it.
"""
from __future__ import annotations

import argparse
import base64
import errno
import gzip
import json
import math
import mimetypes
import os
import re
import shutil
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

from genomics.visualizer import launch
from genomics.visualizer.alignment import AlignmentService
from genomics.visualizer.alphagenome import AlphaGenomeBackend, create_backend
from genomics.visualizer.gene_index import GeneIndex, find_table
from genomics.visualizer.genotypes import GENOTYPE_LABELS, GenotypeService
from genomics.visualizer.controls import ControlError, ControlService
from genomics.visualizer.scalars import ScalarError, ScalarService, ScalarStore
from genomics.visualizer.sessions import SessionError, SessionStore
from genomics.visualizer.ancestry import AncestryService, default_windows, match_groups, window_coverage
from genomics.visualizer.gtex import DEFAULT_TISSUES, GtexClient, GtexError
from genomics.visualizer.annotations import AnnotationService
from genomics.visualizer.cli_args import build_arg_parser
from genomics.visualizer.cache import DiskArrayCache, stable_key
from genomics.visualizer.datasets import DatasetCatalog, DatasetMemory
from genomics.visualizer.experiments import ExperimentService
from genomics.visualizer.jobs import JobManager
from genomics.visualizer.knowledge import KnowledgeBase, public_record
from genomics.visualizer.observed import SCALES as OBSERVED_SCALES, ObservedService
from genomics.visualizer.perturb import PerturbError, PerturbService
from genomics.visualizer.expression import DEFAULT_ONTOLOGIES, TISSUES as EXPRESSION_TISSUES, ExpressionService, tissue_spec
from genomics.visualizer.report import ReportService
from genomics.visualizer.variant_knowledge import VariantKnowledge
from genomics.visualizer.products import ProductService, coding_start_stop, fasta as products_fasta
from genomics.visualizer.structures import (ESMFOLD_MAX, StructureError, StructureService, compare_proteins, protein_similarity, residue_context,
                                            rna_available, rna_comparison)
from genomics.visualizer.sequences import SequenceService
from genomics.visualizer import startup
from genomics.visualizer.startup import APP_NAME, DEFAULT_PORT, PORT_ATTEMPTS, StartupError
from genomics.visualizer.remote import RemoteCache, RemoteError
from genomics.visualizer.signals import COORDINATE_SYSTEMS, DIPLOID, MAX_GROUP_TRACKS, SignalService, _clamp_range, bin_matrix, parse_series
from genomics.visualizer.tasks import TaskManager
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
        alphagenome: Optional[AlphaGenomeBackend] = None,
        allow_add_datasets: bool = True,
        tasks_dir: Optional[Path] = None,
        allow_tasks: bool = True,
        dataset_memory: Optional[DatasetMemory] = None,
        perturb_default_config: Optional[Path] = None,
        remote: bool = True,
        scalar_store: Optional[ScalarStore] = None,
        session_store: Optional[SessionStore] = None,
    ):
        self.catalog = catalog
        self.cache_dir = cache_dir
        self.alignment = AlignmentService(model_window=model_window, cache_bytes=max(memory_bytes // 8, 64 << 20), consensus_dirs=consensus_dirs)
        self.signals = SignalService(memory_bytes, DiskArrayCache(cache_dir), alignment=self.alignment, workers=workers)
        self.sequences = SequenceService(self.signals)
        self.annotations = AnnotationService(cache_dir, gtf=gtf)
        self.products = ProductService(self.signals, self.annotations)
        self.gtf = gtf
        self.gene_index = GeneIndex()
        self.remote = RemoteCache(cache_dir, offline=not remote)
        self.knowledge = KnowledgeBase(self.remote)
        self.observed = ObservedService(self.remote, DiskArrayCache(cache_dir), cache_bytes=max(memory_bytes // 8, 64 << 20))
        self.genotypes = GenotypeService(self.signals, DiskArrayCache(cache_dir), cache_bytes=max(memory_bytes // 8, 64 << 20))
        self.gtex = GtexClient(self.remote)
        self.ancestry = AncestryService(self.genotypes, DiskArrayCache(cache_dir))
        self.scalars = ScalarService(self.genotypes, scalar_store, ancestry=self.ancestry)
        self.controls = ControlService(self.signals)
        self.sessions = session_store if session_store is not None else SessionStore()
        self._ag_catalog: Optional[Dict[str, Any]] = None
        self._ag_catalog_lock = threading.Lock()
        self.experiments = ExperimentService(runs_roots)
        self.views = ViewService()
        self.alphagenome = alphagenome
        self.jobs = JobManager(workers=2)
        self.allow_add_datasets = allow_add_datasets
        self.allow_tasks = allow_tasks and tasks_dir is not None
        self.tasks = TaskManager(tasks_dir) if tasks_dir is not None else None
        self.dataset_memory = dataset_memory
        self.perturb = PerturbService(self, default_config=perturb_default_config)
        self.structures = StructureService(self.remote, cache_dir)
        self.expression = ExpressionService(self)
        self.report = ReportService(self)
        self.variant_knowledge = VariantKnowledge(self.remote, cache_dir)
        if self.tasks is not None:
            self.tasks.on_finished("import", self._on_import_finished)
            self.tasks.on_finished("predict", self._on_predict_finished)
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
        r("GET", r"/api/d/(?P<ds>[^/]+)/observed", self.api_observed)
        r("GET", r"/api/d/(?P<ds>[^/]+)/sequence", self.api_sequence)
        r("GET", r"/api/d/(?P<ds>[^/]+)/variant/sites", self.api_variant_sites)
        r("GET", r"/api/d/(?P<ds>[^/]+)/variant/site", self.api_variant_site)
        r("GET", r"/api/d/(?P<ds>[^/]+)/variant/effect", self.api_variant_effect)
        r("GET", r"/api/d/(?P<ds>[^/]+)/scalars", self.api_scalars)
        r("POST", r"/api/d/(?P<ds>[^/]+)/scalars", self.api_scalar_save)
        r("POST", r"/api/d/(?P<ds>[^/]+)/scalars/delete", self.api_scalar_delete)
        r("GET", r"/api/d/(?P<ds>[^/]+)/ancestry/windows", self.api_ancestry_windows)
        r("GET", r"/api/d/(?P<ds>[^/]+)/ancestry/pca", self.api_ancestry_pca)
        r("POST", r"/api/d/(?P<ds>[^/]+)/ancestry/fields", self.api_ancestry_fields)
        r("POST", r"/api/d/(?P<ds>[^/]+)/ancestry/match", self.api_ancestry_match)
        r("GET", r"/api/gtex/variant", self.api_gtex_variant)
        r("GET", r"/api/gtex/resolve", self.api_gtex_resolve)
        r("GET", r"/api/gtex/tissues", self.api_gtex_tissues)
        r("GET", r"/api/d/(?P<ds>[^/]+)/composition", self.api_composition)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products", self.api_products)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products/fasta", self.api_products_fasta)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products/expression", self.api_products_expression)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products/report", self.api_products_report)
        r("GET", r"/api/products/tissues", self.api_products_tissues)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products/known", self.api_products_known)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products/structure", self.api_products_structure)
        r("POST", r"/api/d/(?P<ds>[^/]+)/products/fold", self.api_products_fold)
        r("GET", r"/api/d/(?P<ds>[^/]+)/products/rna", self.api_products_rna)
        r("POST", r"/api/d/(?P<ds>[^/]+)/views/preview", self.api_view_preview)
        r("POST", r"/api/d/(?P<ds>[^/]+)/views/save", self.api_view_save)
        r("GET", r"/api/runs", self.api_runs)
        r("GET", r"/api/runs/detail", self.api_run_detail)
        r("GET", r"/api/runs/file", self.api_run_file)
        r("POST", r"/api/runs/evaluate", self.api_run_evaluate)
        r("POST", r"/api/datasets/forget", self.api_forget_dataset)
        r("POST", r"/api/d/(?P<ds>[^/]+)/refresh", self.api_refresh_dataset)
        r("GET", r"/api/fs/complete", self.api_fs_complete)
        r("GET", r"/api/tasks", self.api_tasks)
        r("GET", r"/api/tasks/(?P<task_id>[^/]+)", self.api_task)
        r("POST", r"/api/tasks/(?P<task_id>[^/]+)/cancel", self.api_task_cancel)
        r("POST", r"/api/tasks/(?P<task_id>[^/]+)/retry", self.api_task_retry)
        r("POST", r"/api/tasks/(?P<task_id>[^/]+)/delete", self.api_task_delete)
        r("GET", r"/api/import/defaults", self.api_import_defaults)
        r("POST", r"/api/import/inspect", self.api_import_inspect)
        r("POST", r"/api/import/metadata", self.api_import_metadata)
        r("POST", r"/api/import/start", self.api_import_start)
        r("GET", r"/api/import/quickstart", self.api_quickstart_presets)
        r("POST", r"/api/import/quickstart", self.api_quickstart_prepare)
        r("GET", r"/api/d/(?P<ds>[^/]+)/predict/options", self.api_predict_options)
        r("POST", r"/api/d/(?P<ds>[^/]+)/predict", self.api_predict_start)
        r("GET", r"/api/d/(?P<ds>[^/]+)/train/options", self.api_train_options)
        r("POST", r"/api/d/(?P<ds>[^/]+)/train/preview", self.api_train_preview)
        r("POST", r"/api/d/(?P<ds>[^/]+)/train/controls", self.api_train_controls)
        r("GET", r"/api/d/(?P<ds>[^/]+)/sessions", self.api_sessions)
        r("POST", r"/api/d/(?P<ds>[^/]+)/sessions", self.api_session_save)
        r("POST", r"/api/d/(?P<ds>[^/]+)/sessions/rename", self.api_session_rename)
        r("POST", r"/api/d/(?P<ds>[^/]+)/sessions/delete", self.api_session_delete)
        r("POST", r"/api/d/(?P<ds>[^/]+)/train", self.api_train_start)
        r("GET", r"/api/perturb/models", self.api_perturb_models)
        r("GET", r"/api/perturb/model", self.api_perturb_model)
        r("GET", r"/api/perturb/score", self.api_perturb_score)
        r("POST", r"/api/perturb/sequence", self.api_perturb_sequence)
        r("POST", r"/api/perturb/apply", self.api_perturb_apply)
        r("GET", r"/api/perturb/result", self.api_perturb_result)
        r("POST", r"/api/perturb/scan", self.api_perturb_scan)
        r("GET", r"/api/perturb/signal", self.api_perturb_signal)
        r("GET", r"/api/system", self.api_system)
        r("GET", r"/api/alphagenome", self.api_alphagenome)
        r("GET", r"/api/alphagenome/catalog", self.api_alphagenome_catalog)
        r("GET", r"/api/genes/search", self.api_gene_search)
        r("GET", r"/api/genes/info", self.api_gene_info)
        r("GET", r"/api/genes/go", self.api_gene_go)
        r("GET", r"/api/genes/names", self.api_gene_names)
        r("GET", r"/api/genesets/search", self.api_geneset_search)
        r("GET", r"/api/genesets/genes", self.api_geneset_genes)
        r("GET", r"/api/ontology/term", self.api_ontology_term)
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
            "app": APP_NAME,
            "version": VERSION,
            "uptime": round(time.time() - self.started, 1),
            "datasets": [d.listing() for d in self.catalog.all()],
            "missing_datasets": self.dataset_memory.missing() if self.dataset_memory is not None else [],
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
            "allow_tasks": self.allow_tasks,
            "tasks_dir": str(self.tasks.root) if self.tasks else None,
            "jobs": self.jobs.active(),
            "tasks": self.tasks.active() if self.tasks else [],
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
        if self.dataset_memory is not None:
            self.dataset_memory.remember(dataset.path, dataset.annotations_path)
        return dataset.listing()

    def api_forget_dataset(self, query: Query, body: Any) -> Any:
        path = str((body or {}).get("path") or "")
        if path and not (body or {}).get("id"):  # a remembered dataset whose directory is gone
            if self.dataset_memory is not None:
                self.dataset_memory.forget(Path(path))
            return {"datasets": [d.listing() for d in self.catalog.all()]}
        dataset = self.dataset(str((body or {}).get("id") or ""))
        if self.dataset_memory is not None:
            self.dataset_memory.forget(dataset.path)
        self.catalog.remove(dataset.id)
        return {"datasets": [d.listing() for d in self.catalog.all()]}

    def api_refresh_dataset(self, query: Query, body: Any, ds: str) -> Any:
        try:
            return self.catalog.reload(ds).listing()
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        except (FileNotFoundError, ValueError) as exc:
            raise HttpError(HTTPStatus.CONFLICT, str(exc))

    def api_fs_complete(self, query: Query, body: Any) -> Any:
        """Directory entries for path inputs (local paths are how datasets and VCFs are chosen)."""
        if not (self.allow_add_datasets or self.allow_tasks):
            raise HttpError(HTTPStatus.FORBIDDEN, "Browsing paths is disabled")
        raw = query.str("path", "") or str(Path.home()) + "/"
        path = Path(raw).expanduser()
        directory, prefix = (path, "") if raw.endswith("/") else (path.parent, path.name)
        if not directory.is_dir():
            return {"directory": str(directory), "entries": []}
        entries = []
        try:
            for child in sorted(directory.iterdir(), key=lambda p: p.name.lower()):
                if child.name.startswith(".") and not prefix.startswith("."):
                    continue
                if prefix and not child.name.startswith(prefix):
                    continue
                is_dir = child.is_dir()
                entries.append({"name": child.name, "path": str(child) + ("/" if is_dir else ""), "dir": is_dir, "dataset": is_dir and (child / "dataset_metadata.json").exists()})
                if len(entries) >= 200:
                    break
        except PermissionError:
            pass
        return {"directory": str(directory), "entries": entries}

    def api_jobs(self, query: Query, body: Any) -> Any:
        return {"jobs": self.jobs.active(), "tasks": self.tasks.active() if self.tasks else []}

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
        self.scalars.sync(dataset)
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

    def _term_names(self, dataset) -> Dict[str, Dict[str, str]]:
        """{CURIE: {biosample_name, biosample_type}} from the dataset and the cached AlphaGenome catalog."""
        names: Dict[str, Dict[str, str]] = {}
        catalog = self._ag_catalog
        if catalog is None and self.cache_dir:
            from genomics.workflows.alphagenome import catalog as ag_catalog

            with self._ag_catalog_lock:
                if self._ag_catalog is None:
                    self._ag_catalog = ag_catalog.load_catalog(Path(self.cache_dir) / "alphagenome_catalog.json")
                catalog = self._ag_catalog
        for term in (catalog or {}).get("ontologies") or []:
            if term.get("name"):
                names[term["curie"]] = {"biosample_name": term["name"], "biosample_type": term.get("type") or ""}
        for curie, detail in (dataset.metadata.get("ontology_details") or {}).items():
            if (detail or {}).get("biosample_name"):
                names[curie] = {"biosample_name": detail["biosample_name"], "biosample_type": detail.get("biosample_type") or ""}
        return names

    def _gene_info(self, dataset, gene: str) -> Dict[str, Any]:
        """``dataset.gene_info`` with tissue names filled in for track metadata that only has a CURIE
        (older predictions stored ``ontology_curie`` and ``strand`` only). Idempotent, in place."""
        from genomics.visualizer.datasets import track_label, track_short_label

        info = dataset.gene_info(gene)
        if info.get("_named"):
            return info
        names = self._term_names(dataset)
        for out in (info.get("outputs") or {}).values():
            for track in out.get("tracks") or []:
                meta = track.get("metadata") or {}
                extra = names.get(meta.get("ontology_curie") or "")
                if extra and not meta.get("biosample_name"):
                    meta.update({k: v for k, v in extra.items() if v and not meta.get(k)})
                    track["label"] = track_label(track["index"], meta)
                    track["short"] = track_short_label(track["index"], meta)
        info["_named"] = True
        return info

    def api_gene(self, query: Query, body: Any, ds: str, gene: str) -> Any:
        dataset = self.dataset(ds)
        info = dict(self._gene_info(dataset, self._gene(dataset, gene)))
        info.pop("_named", None)
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
        if field.startswith("variant:"):
            # Genotype classes at one site: "variant:<pos>:<ref>:<alt>" (from the Variant page).
            geno, pending = self._genotypes_or_pending(dataset, gene)
            if pending is not None:
                return pending
            site = self._site(geno, field)
            labels = self._genotype_labels(geno, site)
            dose = geno.dosage(site)
            gindex = {sample: i for i, sample in enumerate(geno.samples)}
            groups = {labels[k]: [smp for smp in cohort if smp in gindex and dose[gindex[smp]] == k] for k in range(3)}
            values = list(groups)
        elif not values:
            counts: Dict[str, int] = {}
            for sample in cohort:
                key = str(dataset.samples[dataset.sample_index[sample]].get(field, ""))
                if key:
                    counts[key] = counts.get(key, 0) + 1
            values = sorted(counts, key=lambda k: (-counts[k], k))[:MAX_GROUPS]
        values = values[:MAX_GROUPS]
        if not field.startswith("variant:"):
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

    # -- region scalars ---------------------------------------------------------------
    # -- saved sessions ------------------------------------------------------------------------
    def _session(self, fn: Callable[[], Any]) -> Any:
        try:
            return fn()
        except SessionError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def api_sessions(self, query: Query, body: Any, ds: str) -> Any:
        """Saved sessions of a dataset; ``name`` returns one with its stored state."""
        dataset = self.dataset(ds)
        name = query.str("name", "")
        if name:
            entry = self.sessions.get(dataset, name)
            if entry is None:
                raise HttpError(HTTPStatus.NOT_FOUND, f"No saved session called {name!r}")
            return entry
        # The listing leaves out the state blobs: the pages only need them when one is opened.
        return {"sessions": [{k: v for k, v in s.items() if k != "state"} for s in self.sessions.list(dataset)]}

    def api_session_save(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        body = body or {}
        return self._session(lambda: self.sessions.save(
            dataset, body.get("name"), body.get("state") or {}, body.get("description") or "", replace=bool(body.get("replace"))))

    def api_session_rename(self, query: Query, body: Any, ds: str) -> Any:
        body = body or {}
        return self._session(lambda: self.sessions.rename(self.dataset(ds), str(body.get("name") or ""), body.get("new_name")))

    def api_session_delete(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        name = str((body or {}).get("name") or "")
        if not self.sessions.delete(dataset, name):
            raise HttpError(HTTPStatus.NOT_FOUND, f"No saved session called {name!r}")
        return {"deleted": name}

    def api_scalars(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        self.scalars.sync(dataset)
        return {"scalars": [self.scalars.describe(dataset, spec) for spec in self.scalars.specs(dataset)]}

    def api_scalar_save(self, query: Query, body: Any, ds: str) -> Any:
        """Define (or with ``replace``, redefine) a region scalar and compute it; pending while the
        per-haplotype region means are computed. The sample fields appear on the next samples listing."""
        dataset = self.dataset(ds)
        body = body or {}
        gene = self._gene(dataset, str(body.get("gene") or ""))
        try:
            spec = self.scalars.validate(dataset, body, self._gene_info(dataset, gene), replace=bool(body.get("replace")))
        except ScalarError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        if not self.scalars.is_cached(dataset, spec):
            key = "region:" + self.scalars.region_key(dataset, spec)
            title = f"{gene} {spec['output']} per haplotype over the region"
            result = self._job_response(key, title, lambda progress: self.genotypes.region_means(dataset, gene, spec["output"], spec["start"], spec["end"], progress))
            if isinstance(result, dict) and result.get("pending"):
                return result
        try:
            self.scalars.save(dataset, spec)
        except ScalarError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        return {"scalar": self.scalars.describe(dataset, spec)}

    def api_scalar_delete(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        name = str((body or {}).get("name") or "")
        if not self.scalars.store.delete(dataset, name):
            raise HttpError(HTTPStatus.NOT_FOUND, f"No region scalar named {name!r}")
        self.scalars.unapply(dataset, name)
        return {"deleted": name}

    # -- ancestry PCA and matching ----------------------------------------------------
    def api_ancestry_windows(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        coverage = window_coverage(dataset)
        listed = {str(g) for g in dataset.metadata.get("genes") or []}
        return {
            "samples": len(dataset.samples),
            "windows": [{"gene": g, "coverage": round(coverage.get(g, 0.0), 3), "listed": g in listed,
                         "cached": self.genotypes.cached(dataset, g) is not None} for g in dataset.genes],
            "default": default_windows(dataset, coverage),
            "derived": [self.scalars.describe(dataset, sp) for sp in self.scalars.specs(dataset) if sp.get("kind") in ("pca", "match")],
        }

    def _pca_params(self, dataset, source: Dict[str, Any]) -> Dict[str, Any]:
        genes = source.get("genes")
        if isinstance(genes, str):
            genes = [g for g in genes.split(",") if g]
        try:
            return self.ancestry.params(dataset, genes or None, float(source.get("min_maf") or 0.05),
                                        int(source.get("spacing") if source.get("spacing") not in (None, "") else 2000),
                                        int(source.get("components") or 10))
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc.args[0]))
        except (TypeError, ValueError) as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def _pca_or_pending(self, dataset, params: Dict[str, Any]):
        result = self.ancestry.cached(dataset, params)
        if result is not None:
            return result, None
        key = "pca:" + self.ancestry.key(dataset, params)
        title = f"Genotype PCA over {len(params['genes'])} window{'s' if len(params['genes']) != 1 else ''}"
        result = self._job_response(key, title, lambda progress: self.ancestry.compute(dataset, params, progress))
        if isinstance(result, dict) and result.get("pending"):
            return None, result
        return result, None

    def api_ancestry_pca(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        params = self._pca_params(dataset, {k: query.str(k, "") for k in ("genes", "min_maf", "spacing", "components")})
        if query.str("cached_only", "") in ("1", "true") and self.ancestry.cached(dataset, params) is None:
            return {"cached": False}
        result, pending = self._pca_or_pending(dataset, params)
        if pending is not None:
            return pending
        payload = self.ancestry.payload(result, params)
        payload["sample_count"] = len(dataset.samples)
        return payload

    def api_ancestry_fields(self, query: Query, body: Any, ds: str) -> Any:
        """Store the PCs of a computed PCA as numeric sample fields pc1..pcK (``count``); count 0 removes them."""
        dataset = self.dataset(ds)
        body = body or {}
        params = self._pca_params(dataset, body.get("params") or {})
        count = int(body.get("count") or 0)
        if count <= 0:
            self.scalars.store.delete(dataset, "pcs")
            self.scalars.unapply(dataset, "pcs")
            return {"fields": []}
        if self.ancestry.cached(dataset, params) is None:
            raise HttpError(HTTPStatus.CONFLICT, "Compute the PCA first")
        spec = {"kind": "pca", "name": "pcs", "params": params, "count": max(1, min(count, params["components"]))}
        try:
            self.scalars.check_name(dataset, "pcs", ScalarService.spec_fields(spec), replace=True, what="PC set")
            self.scalars.save(dataset, spec)
        except ScalarError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        return {"fields": ScalarService.spec_fields(spec), "derived": self.scalars.describe(dataset, spec)}

    def api_ancestry_match(self, query: Query, body: Any, ds: str) -> Any:
        """Pair group A with group B (values of one sample field, within the cohort) on the first k PCs and
        store the matched samples as a categorical field (``name``) with the two group labels."""
        dataset = self.dataset(ds)
        body = body or {}
        params = self._pca_params(dataset, body.get("params") or {})
        result = self.ancestry.cached(dataset, params)
        if result is None:
            raise HttpError(HTTPStatus.CONFLICT, "Compute the PCA first")
        field_name = str(body.get("field") or "")
        a_values = [str(v) for v in body.get("group_a") or []]
        b_values = [str(v) for v in body.get("group_b") or []]
        if not field_name or not a_values or not b_values:
            raise HttpError(HTTPStatus.BAD_REQUEST, "field, group_a and group_b are required")
        if set(a_values) & set(b_values):
            raise HttpError(HTTPStatus.BAD_REQUEST, "The two groups overlap")
        name = str(body.get("name") or "").strip().lower()
        label_a = str(body.get("label_a") or "+".join(a_values))[:40]
        label_b = str(body.get("label_b") or "+".join(b_values))[:40]
        if label_a == label_b:
            raise HttpError(HTTPStatus.BAD_REQUEST, "The two groups need different labels")
        filters = body.get("filters") or {}
        if not isinstance(filters, dict):
            raise HttpError(HTTPStatus.BAD_REQUEST, "filters must be a JSON object")
        cohort = set(dataset.filter_samples({str(k): [str(x) for x in (v or [])] for k, v in filters.items()}))
        samples = result["samples"].tolist()

        def members(values):
            wanted = set(values)
            return [i for i, s in enumerate(samples) if s in cohort and str(dataset.samples[dataset.sample_index[s]].get(field_name, "")) in wanted]

        group_a, group_b = members(a_values), members(b_values)
        if not group_a or not group_b:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"No samples with PCs in {'group A' if not group_a else 'group B'} (samples need genotypes in every PCA window)")
        try:
            k = int(body.get("k") or 4)
            caliper = float(body.get("caliper") if body.get("caliper") not in (None, "") else 0.2)
            matched = match_groups(result["scores"], group_a, group_b, k, caliper)
        except (TypeError, ValueError) as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        values = {}
        for ia, ib, _ in matched["pairs"]:
            values[samples[ia]] = label_a
            values[samples[ib]] = label_b
        spec = {"kind": "match", "name": name, "params": params, "field": field_name, "group_a": a_values, "group_b": b_values,
                "label_a": label_a, "label_b": label_b, "k": matched["k"], "caliper": caliper, "filters": filters, "values": values}
        try:
            self.scalars.check_name(dataset, name, [name], replace=bool(body.get("replace")), what="matched set")
            if not values:
                raise ScalarError("No pairs within the caliper; widen it or use fewer PCs")
            self.scalars.save(dataset, spec)
        except ScalarError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        distances = [d for _, _, d in matched["pairs"]]
        return {
            "name": name, "label_a": label_a, "label_b": label_b, "pairs": len(matched["pairs"]),
            "group_a": len(group_a), "group_b": len(group_b), "k": matched["k"], "caliper_distance": matched["caliper_distance"],
            "median_distance": float(np.median(distances)) if distances else None,
            "balance_before": matched["balance_before"], "balance_after": matched["balance_after"],
        }

    # -- variants (Variant page) -------------------------------------------------------
    def _genotypes_or_pending(self, dataset, gene: str):
        """(CohortGenotypes, None) once built, else (None, pending job response)."""
        geno = self.genotypes.cached(dataset, gene)
        if geno is not None:
            return geno, None
        key = "genotypes:" + self.genotypes.key(dataset, gene)
        result = self._job_response(key, f"Genotypes of the cohort in {gene}", lambda progress: self.genotypes.genotypes(dataset, gene, progress))
        if isinstance(result, dict) and result.get("pending"):
            return None, result
        return result, None

    @staticmethod
    def _site(geno, spec: str) -> int:
        """Site index of "variant:<pos>:<ref>:<alt>" or "<pos>:<ref>:<alt>"."""
        parts = spec.split(":")
        if parts[0] == "variant":
            parts = parts[1:]
        if len(parts) != 3 or not parts[0].isdigit():
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Variant must be <pos>:<ref>:<alt>, got {spec!r}")
        try:
            return geno.find(int(parts[0]), parts[1].upper(), parts[2].upper())
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc.args[0]))

    @staticmethod
    def _genotype_labels(geno, site: int) -> List[str]:
        ref, alt = geno.refs[site], geno.alts[site]
        if len(ref) <= 3 and len(alt) <= 3:
            return [f"{ref}/{ref}", f"{ref}/{alt}", f"{alt}/{alt}"]
        return list(GENOTYPE_LABELS)

    def _variant_query(self, dataset, query: Query):
        gene = self._gene(dataset, query.str("gene", required=True))
        geno, pending = self._genotypes_or_pending(dataset, gene)
        if pending is not None:
            return gene, None, None, pending
        spec = f"{query.int('pos', 0, lo=1)}:{query.str('ref', required=True)}:{query.str('alt', required=True)}"
        return gene, geno, self._site(geno, spec), None

    def api_variant_sites(self, query: Query, body: Any, ds: str) -> Any:
        """Sites carried in the cohort over genomic positions [start, end] (default: the window), by position."""
        dataset = self.dataset(ds)
        gene = self._gene(dataset, query.str("gene", required=True))
        geno, pending = self._genotypes_or_pending(dataset, gene)
        if pending is not None:
            return pending
        window = dataset.window(gene)
        lo = query.int("start", window.start or 0)
        hi = query.int("end", window.end or 2 ** 62)
        min_af = float(query.str("min_af", "0") or 0)
        limit = query.int("limit", 2000, lo=1, hi=20000)
        n_haps = 2 * max(len(geno.samples), 1)
        i0 = int(np.searchsorted(geno.positions, lo, side="left"))
        i1 = int(np.searchsorted(geno.positions, hi, side="right"))
        counts = geno.allele_counts()
        sites = []
        for i in range(i0, i1):
            af = float(counts[i]) / n_haps
            if af < min_af:
                continue
            rec = self.genotypes.site_record(geno, i, window.chromosome)
            rec["af"] = af
            sites.append(rec)
            if len(sites) >= limit:
                break
        return {"gene": gene, "chromosome": window.chromosome, "samples": len(geno.samples), "sites": sites, "truncated": len(sites) >= limit}

    def api_variant_site(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        gene, geno, site, pending = self._variant_query(dataset, query)
        if pending is not None:
            return pending
        window = dataset.window(gene)
        field = query.str("field", "") or None
        if field and not any(f["name"] == field for f in dataset.fields):
            field = None
        record = self.genotypes.site_record(geno, site, window.chromosome)
        record.update(gene=gene, window={"chromosome": window.chromosome, "start": window.start, "end": window.end},
                      labels=self._genotype_labels(geno, site), af_dataset=float(geno.allele_counts()[site]) / (2 * len(geno.samples)))
        record["cohort"] = self.genotypes.site_summary(dataset, geno, site, self._cohort(dataset, query), field)
        return record

    def api_variant_effect(self, query: Query, body: Any, ds: str) -> Any:
        """AlphaGenome track summarised over [start, end) (reference offsets) by genotype at the site."""
        dataset = self.dataset(ds)
        gene, geno, site, pending = self._variant_query(dataset, query)
        if pending is not None:
            return pending
        output = query.str("output", required=True)
        track = query.int("track", 0, lo=0)
        start = query.int("start", 0, lo=0)
        end = query.int("end", start + 1, lo=start + 1)
        if end - start > 2_000_000:
            raise HttpError(HTTPStatus.BAD_REQUEST, "Region too long")
        region = self.genotypes.cached_region_means(dataset, gene, output, start, end)
        if region is None:
            key = "region:" + self.genotypes.region_key(dataset, gene, output, start, end)
            result = self._job_response(key, f"{gene} {output} per haplotype over the region", lambda progress: self.genotypes.region_means(dataset, gene, output, start, end, progress))
            if isinstance(result, dict) and result.get("pending"):
                return result
            region = result
        if track >= region["means"].shape[2]:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"{output} has {region['means'].shape[2]} tracks")
        payload = self.genotypes.effect(geno, site, region, track, self._cohort(dataset, query))
        payload.update(gene=gene, output=output, track=track, start=start, end=end, labels=self._genotype_labels(geno, site))
        return payload

    def api_gtex_variant(self, query: Query, body: Any) -> Any:
        variant_id = query.str("variant_id", required=True)
        if not re.match(r"^chr[0-9XYM]+_\d+_[ACGTN]+_[ACGTN]+_b38$", variant_id):
            raise HttpError(HTTPStatus.BAD_REQUEST, "variant_id must look like chr15_28120472_A_G_b38")
        tissues = query.list("tissues") or list(DEFAULT_TISSUES)
        try:
            return self.gtex.report(variant_id, query.str("gene", "") or None, tissues[:8])
        except GtexError as exc:
            raise HttpError(HTTPStatus.BAD_GATEWAY, f"GTEx: {exc}")

    def api_gtex_resolve(self, query: Query, body: Any) -> Any:
        rsid = query.str("rsid", required=True).strip()
        if not re.match(r"^rs\d+$", rsid, re.I):
            raise HttpError(HTTPStatus.BAD_REQUEST, "Expected an rsID such as rs12913832")
        try:
            info = self.gtex.variant(snp_id=rsid.lower())
        except GtexError as exc:
            raise HttpError(HTTPStatus.BAD_GATEWAY, f"GTEx: {exc}")
        if info is None:
            raise HttpError(HTTPStatus.NOT_FOUND, f"{rsid} is not in GTEx v8 (only variants with MAF >= 1% in GTEx donors are)")
        return {"rsid": info.get("snpId"), "variant_id": info.get("variantId"), "chromosome": info.get("chromosome"), "pos": info.get("pos"), "ref": info.get("ref"), "alt": info.get("alt")}

    def api_gtex_tissues(self, query: Query, body: Any) -> Any:
        try:
            return {"tissues": self.gtex.tissues(), "default": list(DEFAULT_TISSUES)}
        except GtexError as exc:
            raise HttpError(HTTPStatus.BAD_GATEWAY, f"GTEx: {exc}")

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

    def api_observed(self, query: Query, body: Any, ds: str) -> Any:
        """Observed (ENCODE / FANTOM5) signal behind AlphaGenome tracks, in the requested coordinates."""
        dataset = self.dataset(ds)
        gene, output, coords, start, end, bins = self._common_signal_params(dataset, query)
        tracks = query.ints("tracks") or [0]
        if len(tracks) > MAX_GROUP_TRACKS:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"At most {MAX_GROUP_TRACKS} tracks")
        scale = query.str("scale", "tpm")
        if scale not in OBSERVED_SCALES:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"scale must be one of {', '.join(OBSERVED_SCALES)}")
        described = self._gene_info(dataset, gene)["outputs"].get(output)
        if described is None:
            raise HttpError(HTTPStatus.NOT_FOUND, f"No {output} predictions for {gene}")
        metas = {t: (described["tracks"][t].get("metadata") or {}) for t in tracks if 0 <= t < len(described["tracks"])}
        missing = [t for t in metas if self.observed.cached(dataset, gene, output, t, scale) is None]
        if missing:
            def load(progress):
                for i, t in enumerate(missing):
                    self.observed.load(dataset, gene, output, t, metas[t], scale, progress.sub(i / len(missing), 1 / len(missing), f"track {t}: "))
                return True

            result = self._job_response("observed:" + stable_key([dataset.id, gene, output, sorted(missing), scale]), f"Observed data for {gene} ({output})", load)
            if isinstance(result, dict) and result.get("pending"):
                return result
        domain = self.signals.domain_length(dataset, gene, output, coords)
        start, end = _clamp_range(start, end, domain)
        bins = max(1, min(int(bins), 8192))
        items = []
        edges = None
        for t in tracks:
            if t not in metas:
                items.append({"track": t, "available": False, "reason": "Unknown track"})
                continue
            hit = self.observed.cached(dataset, gene, output, t, scale) or self.observed.load(dataset, gene, output, t, metas[t], scale)
            item = {"track": t, **hit["info"]}
            if hit["values"] is not None:
                try:
                    window = self.signals.reference_indexed_window(dataset, gene, hit["values"].reshape(-1, 1), 1, coords, start, end)
                except RuntimeError as exc:
                    item.update(available=False, reason=str(exc))
                else:
                    binned = bin_matrix(window, bins)
                    edges = binned["edges"]
                    item.update(mean=binned["mean"][0], min=binned["min"][0], max=binned["max"][0])
            items.append(item)
        if edges is None:
            edges = bin_matrix(np.zeros((end - start, 1), np.float32), bins)["edges"]
        return {"gene": gene, "output": output, "coords": coords, "start": start, "end": end, "domain": domain,
                "edges": (edges + start).astype(np.int64), "scale": scale, "tracks": items}

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

    def _products_context(self, query: Query, ds: str) -> Any:
        """(dataset, gene, sample, gene models, target) of a products request, or a response to return as is."""
        dataset = self.dataset(ds)
        gene = self._gene(dataset, query.str("gene", required=True))
        sample = query.str("sample", required=True)
        if not re.match(r"^[A-Za-z0-9._-]+$", sample):
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Invalid sample id: {sample}")

        def empty(message: str) -> Dict[str, Any]:  # data the dataset lacks: an answer, not an error
            return {"gene": gene, "sample": sample, "genes": [], "target": None, "transcripts": [], "message": message}

        if not dataset.gene_dir(sample, gene).is_dir():
            return empty(f"{sample} has no {gene} window in this dataset")
        cached = self.annotations.cached(dataset)
        if cached is None:
            result = self._job_response(f"annotations:{dataset.id}", "Loading gene annotations", lambda progress: self.annotations.build(dataset, progress))
            if isinstance(result, dict) and result.get("pending"):
                return result
            cached = result
        if not cached.get("source"):
            return empty("No gene annotations for this dataset: add gtf_cache.feather to it or start with --gtf")
        return dataset, gene, sample, cached.get("genes", {}).get(gene, []), query.str("target", "") or None

    def _products(self, query: Query, ds: str, with_sequence: bool) -> Any:
        ctx = self._products_context(query, ds)
        if isinstance(ctx, dict):
            return ctx
        dataset, gene, sample, models, target = ctx
        try:
            return self.products.products(dataset, sample, gene, models, target=target, with_sequence=with_sequence)
        except FileNotFoundError as exc:
            return {"gene": gene, "sample": sample, "genes": [], "target": None, "transcripts": [], "message": str(exc)}

    @staticmethod
    def _product_item(payload: Dict[str, Any], tx: str) -> Dict[str, Any]:
        """An annotated transcript (``tx`` = its id) or a candidate isoform (``cand:<i>``) of a products payload."""
        base = next((t for t in payload["transcripts"] if t["id"] == payload.get("base")), None)
        if tx.startswith("cand:"):
            cands = (payload.get("splicing") or {}).get("candidates") or []
            try:
                cand = cands[int(tx[5:])]
            except (ValueError, IndexError):
                raise HttpError(HTTPStatus.NOT_FOUND, f"No candidate {tx}")
            return {"exons": cand["exons"], "products": cand["products"], "transcript": base, "reference_protein": base["products"]["ref"].get("protein") or "",
                    "title": f"{payload.get('base_name')} with {cand.get('kind')}"}
        t = next((t for t in payload["transcripts"] if t["id"] == tx), None)
        if t is None:
            raise HttpError(HTTPStatus.NOT_FOUND, f"No transcript {tx}")
        return {"exons": t["exons"], "products": t["products"], "transcript": t, "reference_protein": t["products"]["ref"].get("protein") or "", "title": t["name"]}

    def api_products_expression(self, query: Query, body: Any, ds: str) -> Any:
        """mRNA and protein of one sample relative to the reference genome, and absolute estimates, per tissue."""
        ctx = self._products_context(query, ds)
        if isinstance(ctx, dict):
            return {"available": False, "message": ctx.get("message")} if not ctx.get("pending") else ctx
        dataset, gene, sample, models, target = ctx
        extra = [c for c in query.list("tissues") if c and c not in DEFAULT_ONTOLOGIES] if query.str("tissues", "") else []
        stub = [{"ontology": c} for c in DEFAULT_ONTOLOGIES + extra]
        hit = self.expression.cached(dataset, sample, gene, target, stub)
        if hit is not None:
            return hit

        def build(progress):
            specs = [self._tissue_spec(c) for c in DEFAULT_ONTOLOGIES + extra]
            return self.expression.compute(dataset, sample, gene, models, target, progress, tissues=specs)

        return self._job_response(self.expression.key(dataset, sample, gene, target, stub), f"Expression of {target or gene} in {sample}", build)

    def api_products_report(self, query: Query, body: Any, ds: str) -> Any:
        """The gene report: reference baseline per transcript, and each haplotype's levels, flags, variants and similarity."""
        ctx = self._products_context(query, ds)
        if isinstance(ctx, dict):
            return {"available": False, "message": ctx.get("message")} if not ctx.get("pending") else ctx
        dataset, gene, sample, models, target = ctx
        tissue = query.str("tissue", DEFAULT_ONTOLOGIES[0])
        if not re.match(r"^[A-Za-z]+:\d+$", tissue):
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Not an ontology term: {tissue}")
        hit = self.report.cached(dataset, sample, gene, target, tissue)
        if hit is not None:
            return hit

        def build(progress):
            progress(0.01, "Resolving the tissue (AlphaGenome catalog, GTEx, Human Protein Atlas)")
            return self.report.compute(dataset, sample, gene, models, target, self._tissue_spec(tissue), progress)

        return self._job_response(self.report.key(dataset, sample, gene, target, tissue), f"Gene report of {target or gene} in {sample} ({tissue})", build)

    def api_products_known(self, query: Query, body: Any, ds: str) -> Any:
        """What is known about the sample's variants in the gene (rsID, ClinVar, traits, literature)."""
        payload = self._products(query, ds, with_sequence=False)
        if not payload.get("transcripts"):
            return {"available": False, "message": payload.get("message") or "No transcripts"} if not payload.get("pending") else payload
        scope = query.str("scope", "coding")
        quiet = {"intron", "upstream", "downstream"}
        alleles = sorted({(v["pos"], v["ref"], v[h]) for v in payload["variants"] for h in ("H1", "H2")
                          if v.get(h) and (scope == "all" or v["region"] not in quiet)})
        key = f"known:{ds}:{query.str('sample')}:{query.str('gene')}:{payload.get('target')}:{scope}"
        return self._job_response(key, f"Known variants in {payload.get('target')}",
                                  lambda progress: {**self.variant_knowledge.annotate(payload["chromosome"], alleles), "scope": scope, "queried": len(alleles)})

    # -- tissues for the products pipeline -------------------------------------------------------
    def _full_catalog(self) -> Optional[Dict[str, Any]]:
        """The AlphaGenome track catalog with per-track rows (cached file or memory; fetched when a backend is ready)."""
        from genomics.workflows.alphagenome import catalog as ag_catalog

        path = Path(self.cache_dir) / "alphagenome_catalog.json" if self.cache_dir else None
        with self._ag_catalog_lock:
            if self._ag_catalog is None and path is not None:
                self._ag_catalog = ag_catalog.load_catalog(path)
            if self._ag_catalog is not None:
                return self._ag_catalog
            if self._backend_summary()["reasons"]:
                return None
            try:
                client = self.alphagenome.create_client(timeout=30.0) if self.alphagenome is not None else None
                if client is None:
                    return None
                self._ag_catalog = ag_catalog.fetch_catalog(client, source=self._backend_summary()["label"])
            except Exception:
                return None
            if path is not None:
                ag_catalog.save_catalog(self._ag_catalog, path)
            return self._ag_catalog

    def _tissue_terms(self) -> Dict[str, Dict[str, Any]]:
        """{ontology: {name, type, rna_seq, junctions, gtex_tissue}} for every term with AlphaGenome RNA-seq."""
        catalog = self._full_catalog()
        if catalog is None:
            return {}
        if getattr(self, "_tissue_terms_cache", None) is not None and self._tissue_terms_cache[0] is catalog:
            return self._tissue_terms_cache[1]
        terms: Dict[str, Dict[str, Any]] = {}
        for row in (catalog.get("tracks") or {}).get("RNA_SEQ", []):
            curie = row.get("ontology_curie")
            if not curie:
                continue
            t = terms.setdefault(curie, {"curie": curie, "name": row.get("biosample_name") or curie, "type": row.get("biosample_type") or "",
                                         "rna_seq": 0, "junctions": 0, "gtex_tissue": None})
            t["rna_seq"] += 1
            if row.get("gtex_tissue"):
                t["gtex_tissue"] = row["gtex_tissue"]
        for row in (catalog.get("tracks") or {}).get("SPLICE_JUNCTIONS", []):
            if row.get("ontology_curie") in terms:
                terms[row["ontology_curie"]]["junctions"] += 1
        self._tissue_terms_cache = (catalog, terms)
        return terms

    def _gtex_ontologies(self) -> Dict[str, str]:
        try:
            return self.gtex.ontology_tissues()
        except Exception:  # offline: GTEx tissues only through AlphaGenome's own GTEx tracks
            return {}

    def _gtex_names(self) -> Dict[str, str]:
        try:
            return self.gtex.name_tissues()
        except Exception:
            return {}

    def _tissue_spec(self, curie: str) -> Dict[str, Any]:
        return tissue_spec(curie, self._tissue_terms().get(curie), self.expression.references, self._gtex_ontologies(), self._gtex_names())

    def api_products_tissues(self, query: Query, body: Any) -> Any:
        """Every ontology term with AlphaGenome RNA-seq, with where its observed baseline would come from."""
        from genomics.visualizer.gtex import tissue_name_key

        terms = self._tissue_terms()
        gtex = self._gtex_ontologies()
        names = self._gtex_names()
        rows = []
        for curie, t in terms.items():
            tissue = t["gtex_tissue"] or gtex.get(curie) or (names.get(tissue_name_key(t["name"])) if t["type"] == "tissue" else None)
            rows.append({**t, "gtex_tissue": tissue, "default": curie in DEFAULT_ONTOLOGIES,
                         "baseline": "gtex" if tissue else ("hpa_cell_line" if t["type"] == "cell_line" else "hpa_single_cell"
                                                            if t["type"] in ("primary_cell", "in_vitro_differentiated_cells") else "none")})
        for spec in EXPRESSION_TISSUES:
            if spec["ontology"] not in terms:
                rows.append({"curie": spec["ontology"], "name": spec["label"], "type": "", "rna_seq": None, "junctions": None,
                             "gtex_tissue": None, "default": True, "baseline": "curated"})
        rows.sort(key=lambda r: (not r["default"], (r["name"] or "").lower()))
        return {"catalog": bool(terms), "tissues": rows}

    def api_products_structure(self, query: Query, body: Any, ds: str) -> Any:
        """The reference protein's AlphaFold DB model and how each haplotype's product maps onto it."""
        payload = self._products(query, ds, with_sequence=False)
        if not payload.get("transcripts"):
            return {"available": False, "message": payload.get("message") or "No transcripts"} if not payload.get("pending") else payload
        item = self._product_item(payload, query.str("tx", required=True))
        ref_protein = item["reference_protein"]
        if not ref_protein:
            return {"available": False, "message": "Non-coding transcript: no protein"}
        tid = item["transcript"]["id"]

        def build(progress):
            progress(0.1, "UniProt and AlphaFold DB")
            try:
                model = self.structures.reference_model(tid, ref_protein)
            except StructureError as exc:
                model = {"error": str(exc)}
            am = self.report.alphamissense(tid, ref_protein)
            products = {}
            for hap, product in item["products"].items():
                info = compare_proteins(ref_protein, product.get("protein") or "")
                info.update(length=len(product.get("protein") or ""), change=product.get("change"), nmd=product.get("nmd"),
                            similarity=protein_similarity(ref_protein, product.get("protein") or ""))
                subs = info.get("substitutions") or []
                if subs and model.get("pdb") and model.get("match") != "different":
                    context = residue_context(model["pdb"], [x["pos"] for x in subs], model.get("offset") or 0)
                    for x in subs:
                        x.update(context.get(x["pos"], {}))
                        score = am["scores"].get(f"{x['ref']}{x['pos']}{x['alt']}")
                        x["alphamissense"] = {"score": score[0], "class": score[1]} if score else None
                products[hap] = info
            return {"available": True, "title": item["title"], "transcript": tid, "reference_length": len(ref_protein), "model": model,
                    "products": products, "esmfold": {"max": ESMFOLD_MAX, "remote": not self.remote.offline}}

        return self._job_response(f"structure:{ds}:{query.str('sample')}:{query.str('gene')}:{payload.get('target')}:{query.str('tx')}", "Protein structure", build)

    def api_products_fold(self, query: Query, body: Any, ds: str) -> Any:
        """ESMFold of a haplotype's product and of the reference protein (same window; sends both to api.esmatlas.com)."""
        payload = self._products(query, ds, with_sequence=False)
        if not payload.get("transcripts"):
            return payload
        item = self._product_item(payload, query.str("tx", required=True))
        hap = query.str("hap", "H1")
        product = item["products"].get(hap)
        if product is None or not product.get("protein") or not item["reference_protein"]:
            raise HttpError(HTTPStatus.BAD_REQUEST, "Nothing to fold: no protein")
        ref_protein, alt = item["reference_protein"], product["protein"]

        def build(progress):
            progress(0.1, "ESMFold (api.esmatlas.com), about 10-30 s per protein")
            try:
                return self.structures.fold_pair(ref_protein, alt, allow_remote=not self.remote.offline)
            except StructureError as exc:
                raise RuntimeError(str(exc))

        key = f"fold:{ds}:{query.str('sample')}:{query.str('gene')}:{payload.get('target')}:{query.str('tx')}:{hap}"
        return self._job_response(key, f"ESMFold of {item['title']} ({hap})", build)

    def api_products_rna(self, query: Query, body: Any, ds: str) -> Any:
        """Secondary structure (ViennaRNA) of a window of the reference and a haplotype's mature transcript."""
        ctx = self._products_context(query, ds)
        if isinstance(ctx, dict):
            return {"available": False, "message": ctx.get("message")} if not ctx.get("pending") else ctx
        dataset, gene, sample, models, target = ctx
        if not rna_available():
            return {"available": False, "message": "RNA folding needs ViennaRNA: pip install ViennaRNA (part of pip install -e '.[visualizer]')"}
        payload = self.products.products(dataset, sample, gene, models, target=target)
        item = self._product_item(payload, query.str("tx", required=True))
        hap = query.str("hap", "H1")
        if hap not in ("H1", "H2"):
            raise HttpError(HTTPStatus.BAD_REQUEST, "hap must be H1 or H2")
        views = self.products.views(dataset, sample, gene)
        center = query.str("center", "start")
        tx = item["transcript"]
        if center in ("start", "stop"):
            start, stop = coding_start_stop({**tx, "strand": payload["strand"]})
            if start is None:  # non-coding: the transcript's 5' end
                exons = sorted(item["exons"])
                start = exons[0][0] if payload["strand"] == "+" else exons[-1][1] - 1
                stop = start
            offset = start if center == "start" else stop
        else:
            try:
                offset = int(center)
            except ValueError:
                raise HttpError(HTTPStatus.BAD_REQUEST, "center must be start, stop or a window offset")
        try:
            result = rna_comparison(views["ref"], views[hap], item["exons"], payload["strand"], offset, query.int("width", 240, lo=40, hi=600))
        except StructureError as exc:
            return {"available": False, "message": str(exc)}
        result.update(available=True, title=item["title"], hap=hap, center=center, offset=offset, window_start=payload.get("window_start"))
        return result

    def api_products(self, query: Query, body: Any, ds: str) -> Any:
        """Mature transcripts and proteins of one sample's haplotypes, with AlphaGenome splicing evidence."""
        return self._products(query, ds, with_sequence=query.str("sequences", "") in ("1", "true", "mrna"))

    def api_products_fasta(self, query: Query, body: Any, ds: str) -> Any:
        kind = query.str("kind", "protein")
        if kind not in ("protein", "mrna"):
            raise HttpError(HTTPStatus.BAD_REQUEST, "kind must be protein or mrna")
        payload = self._products(query, ds, with_sequence=kind == "mrna")
        if isinstance(payload, dict) and payload.get("pending"):
            return payload
        return RawResponse(products_fasta(payload, kind).encode("utf-8"), "text/plain; charset=utf-8")

    def api_composition(self, query: Query, body: Any, ds: str) -> Any:
        """Per-base letter frequencies over the given haplotypes (the reference without ``rows``)."""
        dataset = self.dataset(ds)
        gene = self._gene(dataset, query.str("gene", required=True))
        coords = self._coords(query)
        start = query.int("start", 0, lo=0)
        end = query.int("end", start + 1, lo=start + 1)
        bins = query.int("bins", 1000, lo=1, hi=8192)
        rows = parse_series(query.str("rows", ""))
        try:
            return self.sequences.composition(dataset, gene, rows, coords, start, end, bins)
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

    def api_run_evaluate(self, query: Query, body: Any) -> Any:
        body = body or {}
        run_id = str(body.get("run") or "")
        run_dir = self.experiments.run_dirs().get(run_id)
        if run_dir is None:
            raise HttpError(HTTPStatus.NOT_FOUND, f"Unknown run: {run_id}")
        try:
            title, steps, params, files = launch.evaluate_task(run_id, run_dir, body)
        except launch.LaunchError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        return self._start_task("evaluate", title, steps, params, files, resource="gpu")

    # -- detached tasks ----------------------------------------------------------------------------
    def _task_manager(self) -> TaskManager:
        if self.tasks is None:
            raise HttpError(HTTPStatus.NOT_FOUND, "Background tasks are not available (no tasks directory)")
        return self.tasks

    def _start_task(self, kind: str, title: str, steps: List[Dict[str, Any]], params: Dict[str, Any], files: Dict[str, str], resource: Optional[str] = None, env: Optional[Dict[str, str]] = None) -> Any:
        if not self.allow_tasks:
            raise HttpError(HTTPStatus.FORBIDDEN, "Starting background tasks is disabled (--no-jobs)")
        from genomics.workspace import repo_root

        return self._task_manager().create(kind, title, steps, params=params, cwd=repo_root(), env=env, resource=resource, files=files)

    def api_tasks(self, query: Query, body: Any) -> Any:
        if self.tasks is None:
            return {"tasks": [], "available": False}
        return {"tasks": self.tasks.list(kind=query.str("kind") or None), "available": True, "allow": self.allow_tasks, "dir": str(self.tasks.root)}

    def api_task(self, query: Query, body: Any, task_id: str) -> Any:
        try:
            return self._task_manager().get(task_id, log_lines=query.int("log", 300, lo=0, hi=2000))
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    def api_task_cancel(self, query: Query, body: Any, task_id: str) -> Any:
        try:
            return self._task_manager().cancel(task_id)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    def api_task_retry(self, query: Query, body: Any, task_id: str) -> Any:
        """Run a finished task again, with the same recipe, as a new task."""
        if not self.allow_tasks:
            raise HttpError(HTTPStatus.FORBIDDEN, "Background tasks are disabled (--no-jobs)")
        try:
            return self._task_manager().retry(task_id)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        except ValueError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def api_task_delete(self, query: Query, body: Any, task_id: str) -> Any:
        try:
            self._task_manager().delete(task_id)
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))
        except RuntimeError as exc:
            raise HttpError(HTTPStatus.CONFLICT, str(exc))
        return {"deleted": task_id}

    def _on_import_finished(self, task: Dict[str, Any]) -> None:
        path = Path(task["params"].get("output_dir") or "")
        if task["status"] != "done" or not (path / "dataset_metadata.json").exists():
            return
        existing = next((d for d in self.catalog.all() if d.path == path.resolve()), None)
        dataset = self.catalog.reload(existing.id) if existing else self.catalog.add(path)
        if self.dataset_memory is not None:
            self.dataset_memory.remember(dataset.path, None)

    def _on_predict_finished(self, task: Dict[str, Any]) -> None:
        path = Path(task["params"].get("dataset_path") or "").resolve()
        for dataset in self.catalog.all():
            if dataset.path == path:
                self.catalog.reload(dataset.id)

    # -- import ------------------------------------------------------------------------------------
    def api_import_defaults(self, query: Query, body: Any) -> Any:
        defaults = launch.import_defaults()
        from genomics.workflows.dataset_builders.vcf_import.builder import missing_tools

        defaults["missing_tools"] = missing_tools()
        defaults["backend"] = self._backend_summary()
        return defaults

    def api_import_inspect(self, query: Query, body: Any) -> Any:
        from genomics.workflows.dataset_builders.vcf_import.builder import ImportSpecError, inspect_vcf

        vcf = str((body or {}).get("vcf") or "").strip()
        if not vcf:
            raise HttpError(HTTPStatus.BAD_REQUEST, "Give a VCF path (or a pattern with {chrom})")
        try:
            return inspect_vcf(vcf)
        except ImportSpecError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        except RuntimeError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Could not read the VCF: {exc}")

    def api_import_metadata(self, query: Query, body: Any) -> Any:
        """Parse a metadata table and match it to the VCF samples (nothing is stored)."""
        from genomics.workflows.dataset_builders.vcf_import import metadata as meta

        body = body or {}
        text = body.get("text")
        if not text and body.get("path"):
            path = Path(str(body["path"])).expanduser()
            if not path.is_file():
                raise HttpError(HTTPStatus.NOT_FOUND, f"File not found: {path}")
            text = path.read_text(encoding="utf-8", errors="replace")
            body.setdefault("filename", path.name)
        if not text:
            raise HttpError(HTTPStatus.BAD_REQUEST, "Provide the metadata file contents")
        samples = [str(s) for s in body.get("samples") or []]
        try:
            columns, rows = meta.parse_table(str(text), str(body.get("filename") or ""))
        except (ValueError, json.JSONDecodeError) as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, f"Could not parse the metadata: {exc}")
        if not columns:
            raise HttpError(HTTPStatus.BAD_REQUEST, "The metadata file is empty")
        id_column = body.get("id_column") if body.get("id_column") in columns else meta.guess_id_column(columns, rows, samples)
        family = body.get("family_column")
        family = family if family in columns else (None if "family_column" in body else meta.guess_role(columns, meta.FAMILY_COLUMN_NAMES))
        sex = body.get("sex_column")
        sex = sex if sex in columns else (None if "sex_column" in body else meta.guess_role(columns, meta.SEX_COLUMN_NAMES))
        records = meta.records_by_sample(rows, id_column, family_column=family, sex_column=sex)
        sample_set = set(samples)
        matched = [s for s in records if s in sample_set] if samples else list(records)
        return {
            "columns": columns,
            "row_count": len(rows),
            "id_column": id_column,
            "family_column": family,
            "sex_column": sex,
            "records": records,
            "matched": len(matched),
            "unmatched_metadata": [s for s in records if samples and s not in sample_set][:200],
            "samples_without_metadata": [s for s in samples if s not in records][:200],
            "samples_without_metadata_count": sum(1 for s in samples if s not in records),
            "fields": meta.describe_fields({s: records[s] for s in matched}),
        }

    def api_import_start(self, query: Query, body: Any) -> Any:
        try:
            title, steps, params, files = launch.import_task(body or {}, alphagenome_ready=not self._backend_summary()["reasons"])
        except launch.LaunchError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        env = self.alphagenome.child_env() if self.alphagenome is not None else None
        return self._start_task("import", title, steps, params, files, env=env)

    def api_quickstart_presets(self, query: Query, body: Any) -> Any:
        from genomics.workflows.dataset_builders.vcf_import import quickstart

        return {"presets": quickstart.presets(), "reference_fasta": launch.default_reference() or quickstart.REMOTE_REFERENCE}

    def api_quickstart_prepare(self, query: Query, body: Any) -> Any:
        """Metadata (downloaded once), a balanced sample subset and the VCF inspection of a public
        cohort, for the import form; nothing is imported here."""
        from genomics.workflows.dataset_builders.vcf_import import quickstart
        from genomics.workflows.dataset_builders.vcf_import.builder import ImportSpecError

        body = body or {}
        try:
            per_population = int(body.get("per_population", quickstart.DEFAULT_PER_POPULATION))
        except (TypeError, ValueError):
            raise HttpError(HTTPStatus.BAD_REQUEST, "per_population must be a number (0 = every sample)")
        cache = Path(self.cache_dir or self._task_manager().root) / "quickstart"
        try:
            return quickstart.prepare(str(body.get("preset") or ""), cache, per_population, body.get("scope"), launch.default_reference())
        except ImportSpecError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        except RuntimeError as exc:
            raise HttpError(HTTPStatus.BAD_GATEWAY, f"Could not read the remote VCF: {exc}")

    # -- AlphaGenome predictions ------------------------------------------------------------------
    def _backend_summary(self) -> Dict[str, Any]:
        if self.alphagenome is None:
            from genomics.core.alphagenome_connection import backend_configured

            return {"label": "environment", "reasons": [] if backend_configured() else ["ALPHAGENOME_API_KEY not set"]}
        return {"label": self.alphagenome.label, "reasons": self.alphagenome.reasons()}

    def api_predict_options(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        ontologies: Dict[str, Dict[str, Any]] = {}
        for curie, detail in (dataset.metadata.get("ontology_details") or {}).items():
            ontologies[curie] = {"curie": curie, "name": (detail or {}).get("biosample_name") or ""}
        for curie in dataset.metadata.get("ontologies") or []:
            ontologies.setdefault(str(curie), {"curie": str(curie), "name": ""})
        existing: Dict[str, int] = {}
        for gene in dataset.genes[:3]:
            info = dataset.gene_info(gene)
            for name, out in {**(info.get("other_outputs") or {}), **(info.get("outputs") or {})}.items():
                existing[name] = len(out.get("tracks") or [])
        windows = []
        for gene in dataset.genes:
            window = dataset.window(gene)
            windows.append({"name": gene, "type": (window.metadata or {}).get("type") or "gene", "chromosome": window.chromosome, "start": window.start, "end": window.end})
        return {
            "outputs": list(launch.PREDICT_OUTPUTS),
            "output_specs": launch.import_defaults_outputs(),
            "existing_outputs": existing,
            "ontologies": sorted(ontologies.values(), key=lambda o: o["curie"]),
            "genes": dataset.genes,
            "windows": windows,
            "window_size": dataset.metadata.get("window_size"),
            "extension": launch.extension_source(dataset),
            "region_presets": launch.region_presets(),
            "gene_table": str(find_table(dataset.path, self.gtf) or "") or None,
            "backend": self._backend_summary(),
            "running": [t for t in (self.tasks.active() if self.tasks else []) if t["kind"] == "predict" and t["params"].get("dataset_path") == str(dataset.path)],
        }

    def api_predict_start(self, query: Query, body: Any, ds: str) -> Any:
        dataset = self.dataset(ds)
        backend = self._backend_summary()
        if backend["reasons"]:
            raise HttpError(HTTPStatus.CONFLICT, f"AlphaGenome backend not ready: {'; '.join(backend['reasons'])}")
        try:
            title, steps, params, files = launch.predict_task(dataset, body or {})
        except launch.LaunchError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        params["backend"] = backend["label"]
        env = self.alphagenome.child_env() if self.alphagenome is not None else None
        return self._start_task("predict", title, steps, params, files, resource="alphagenome", env=env)

    # -- training ----------------------------------------------------------------------------------
    def _runs_root(self) -> Path:
        if self.experiments.roots:
            return self.experiments.roots[0]
        from genomics.workspace import DEFAULT_GENOTYPE_RUNS_ROOT

        return Path(DEFAULT_GENOTYPE_RUNS_ROOT)

    def _train_config(self, dataset, body: Any):
        from genomics.workspace import cache_path

        try:
            config, summary = launch.build_train_config(dataset, body or {}, self.experiments.run_dirs(), self._runs_root(), Path(cache_path("genotype_based_predictor", "visualizer")))
            launch.validate_config(config, Path(self.cache_dir or self._task_manager().root) / "tmp")
        except launch.LaunchError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        return config, summary

    def api_train_options(self, query: Query, body: Any, ds: str) -> Any:
        options = launch.train_options(self.dataset(ds))
        options["runs_root"] = str(self._runs_root())
        options["runs"] = [{"id": run_id, "name": path.name} for run_id, path in self.experiments.run_dirs().items() if (path / "config.yaml").exists()]
        return options

    def api_train_preview(self, query: Query, body: Any, ds: str) -> Any:
        import yaml

        config, summary = self._train_config(self.dataset(ds), body)
        return {"summary": summary, "config": yaml.safe_dump(config, sort_keys=False, allow_unicode=True)}

    def api_train_controls(self, query: Query, body: Any, ds: str) -> Any:
        """A control window per chosen window, matched on the signal the model reads (a background job)."""
        body = body or {}
        dataset = self.dataset(ds)
        genes = [str(g) for g in (body.get("genes") or [])]
        output = str(body.get("output") or "rna_seq").lower()
        terms = [str(t) for t in (body.get("ontology_terms") or [])]
        size = max(1, int(body.get("window_center_size") or 32768))
        if not terms:
            info = dataset.gene_info(genes[0]) if genes else {}
            terms = sorted({str((t.get("metadata") or {}).get("ontology_curie"))
                            for t in ((info.get("outputs") or {}).get(output) or {}).get("tracks", [])
                            if (t.get("metadata") or {}).get("ontology_curie")})
        key = stable_key({"v": 1, "dataset": dataset.fingerprint, "genes": sorted(genes), "output": output, "terms": sorted(terms), "size": size})
        # Cheap checks first, so a bad request is a 400 rather than a failed background job.
        try:
            unknown = [g for g in genes if g not in dataset.genes]
            if not genes:
                raise ControlError("Choose at least one window for the panel")
            if unknown:
                raise ControlError(f"Unknown window(s): {', '.join(unknown)}")
            candidates = self.controls.candidates(dataset, genes)
            if len(candidates) < len(set(genes)):
                raise ControlError("This dataset has fewer windows outside its gene panel than the panel needs")
            # Which windows have a reference prediction is a file check, so it need not cost a job.
            missing = [g for g in dict.fromkeys(genes) if not dataset.reference_prediction_path(g, output).exists()]
            if missing:
                raise ControlError(
                    f"No reference {output} prediction for {', '.join(sorted(missing))}. "
                    "Predict the reference window of every chosen gene first (haplotype 'ref').")
            usable = [g for g in candidates if dataset.reference_prediction_path(g, output).exists()]
            if len(usable) < len(set(genes)):
                raise ControlError(
                    f"Only {len(usable)} control window(s) have a reference {output} prediction, for a panel of "
                    f"{len(set(genes))}. Predict the reference window of the control windows first.")
        except ControlError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        return self._job_response(
            f"controls:{key}", f"Matching control windows ({len(genes)} windows)",
            lambda progress: self.controls.matched_panel(dataset, genes, output, terms, size, progress))

    def api_train_start(self, query: Query, body: Any, ds: str) -> Any:
        config, summary = self._train_config(self.dataset(ds), body)
        title, steps, params, files = launch.train_task(config, summary, evaluate_test=bool((body or {}).get("evaluate_test", True)))
        return self._start_task("train", title, steps, params, files, resource="gpu")

    # -- perturbation lab --------------------------------------------------------------------------
    def _perturb(self, fn: Callable[[], Any]) -> Any:
        try:
            return fn()
        except PerturbError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))
        except FileNotFoundError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc))

    def api_perturb_models(self, query: Query, body: Any) -> Any:
        data = self.perturb.models()
        data["backend"] = self._backend_summary()
        data["training"] = [t["title"] for t in (self.tasks.active() if self.tasks else []) if t["kind"] in ("train", "evaluate")]
        return data

    def api_perturb_model(self, query: Query, body: Any) -> Any:
        run_id = query.str("run", required=True)
        checkpoint = query.str("checkpoint") or None
        if self.perturb.context is not None and self.perturb.context_id == (run_id, checkpoint or "best_accuracy.pt"):
            return self._perturb(lambda: self.perturb.describe(self.perturb._label_field(self.perturb.dataset())))
        return self._job_response(f"perturb-load:{run_id}:{checkpoint}", "Loading model for the Perturbation Lab", lambda progress: self.perturb.load(run_id, checkpoint, progress))

    def api_perturb_score(self, query: Query, body: Any) -> Any:
        sample = query.str("sample", required=True)
        return self._perturb(lambda: self.perturb.score(sample))

    def api_perturb_sequence(self, query: Query, body: Any) -> Any:
        body = body or {}
        return self._perturb(lambda: self.perturb.sequence_view(
            str(body.get("sample") or ""), str(body.get("gene") or ""), body.get("haplotypes") or ["H1", "H2"], body.get("edits") or [],
            int(body.get("start") or 0), int(body.get("end") or 1), int(body.get("bins") or 1000)))

    def api_perturb_apply(self, query: Query, body: Any) -> Any:
        body = body or {}
        sample, gene = str(body.get("sample") or ""), str(body.get("gene") or "")
        outputs = [str(o) for o in body.get("outputs") or []]
        edits = self._perturb(lambda: self.perturb.normalize_edits(body.get("edits") or []))
        key = self._perturb(lambda: self.perturb.apply_key(sample, gene, edits, outputs))
        cached = self.perturb.result(key)
        if cached is not None:
            return cached
        job = self.jobs.run(f"perturb:{key}", f"Perturbation: {sample} {gene}", lambda progress: self.perturb.apply(sample, gene, edits, outputs, progress))
        return {"pending": True, "key": key, "job": job.as_dict()}

    def api_perturb_scan(self, query: Query, body: Any) -> Any:
        """Saturation scan of the loaded model (a background job; poll by posting the same body)."""
        spec = self._perturb(lambda: self.perturb.scan_spec(body or {}))
        key = self.perturb.scan_key(spec)
        positions = (spec["end"] - spec["start"] - spec["size"]) // spec["step"] + 1
        return self._perturb(lambda: self._job_response(f"scan:{key}", f"Saturation scan: {spec['sample']} {spec['gene']} ({positions} windows)", lambda progress: self.perturb.scan(spec, progress)))

    def api_perturb_result(self, query: Query, body: Any) -> Any:
        key = query.str("key", required=True)
        cached = self.perturb.result(key)
        if cached is not None:
            return cached
        job = self.jobs.peek(f"perturb:{key}")
        if job is None:
            raise HttpError(HTTPStatus.NOT_FOUND, "This perturbation is no longer running; apply it again")
        if job.status == "error":
            raise HttpError(HTTPStatus.INTERNAL_SERVER_ERROR, job.error or "Perturbation failed")
        if job.status == "cancelled":
            raise HttpError(HTTPStatus.CONFLICT, "Cancelled")
        if job.status == "done" and isinstance(job.result, dict):
            return job.result
        return {"pending": True, "key": key, "job": job.as_dict()}

    def api_perturb_signal(self, query: Query, body: Any) -> Any:
        return self._perturb(lambda: self.perturb.signal(
            query.str("key", required=True), query.str("output", required=True), query.ints("tracks"),
            query.int("start", 0, lo=0), query.int("end", 1, lo=1), query.int("bins", 1000, lo=1, hi=8192)))

    # -- AlphaGenome backend -------------------------------------------------------------------
    def _alphagenome(self):
        if self.alphagenome is None:
            raise HttpError(HTTPStatus.NOT_FOUND, "AlphaGenome backend settings are not available")
        return self.alphagenome

    def api_system(self, query: Query, body: Any) -> Any:
        """``genomics doctor`` for this visualizer: features, hardware, and free space where it writes."""
        from genomics import doctor

        data = doctor.report(include_local=False, backend=self._backend_summary())
        locations = data["storage"]["locations"]
        if self.cache_dir:
            locations["visualizer cache"] = doctor.free_space(Path(self.cache_dir))
        else:
            locations.pop("visualizer cache", None)
        if self.tasks is not None:
            locations["background jobs"] = doctor.free_space(self.tasks.root)
        for dataset in self.catalog.all():
            locations[f"dataset {dataset.id}"] = doctor.free_space(Path(dataset.path))
        return data

    def api_alphagenome(self, query: Query, body: Any) -> Any:
        return self._alphagenome().describe()

    def api_alphagenome_catalog(self, query: Query, body: Any) -> Any:
        """Every AlphaGenome output and ontology term with track counts (fetched once, then cached)."""
        from genomics.workflows.alphagenome import catalog as ag_catalog

        path = Path(self.cache_dir) / "alphagenome_catalog.json" if self.cache_dir else None
        refresh = query.str("refresh", "") in ("1", "true")
        with self._ag_catalog_lock:
            if self._ag_catalog is None and path is not None and not refresh:
                self._ag_catalog = ag_catalog.load_catalog(path)
            if self._ag_catalog is not None and not refresh:
                return ag_catalog.summary(self._ag_catalog)
            backend = self._backend_summary()
            if backend["reasons"]:
                raise HttpError(HTTPStatus.SERVICE_UNAVAILABLE, f"AlphaGenome backend not ready: {'; '.join(backend['reasons'])}")
            try:
                if self.alphagenome is not None:
                    client = self.alphagenome.create_client(timeout=30.0)
                else:
                    from genomics.core.alphagenome_connection import create_dna_client

                    client = create_dna_client(timeout=30.0)
                fetched = ag_catalog.fetch_catalog(client, source=backend["label"])
            except ImportError as exc:
                raise HttpError(HTTPStatus.SERVICE_UNAVAILABLE, f"{exc}; install the client with: pip install -e '.[alphagenome]'")
            except Exception as exc:
                raise HttpError(HTTPStatus.BAD_GATEWAY, f"Could not read the AlphaGenome track catalog: {type(exc).__name__}: {exc}")
            self._ag_catalog = fetched
            if path is not None:
                ag_catalog.save_catalog(fetched, path)
            return ag_catalog.summary(fetched)

    def _gene_table(self, ds: str) -> Optional[Path]:
        dataset_path = self.dataset(ds).path if ds else None
        return find_table(dataset_path, self.gtf)

    def _hgnc(self):
        """The HGNC table, or None when it cannot be fetched (searches then use GENCODE names only)."""
        try:
            return self.knowledge.hgnc(required=False)
        except RemoteError:
            return None

    def api_gene_search(self, query: Query, body: Any) -> Any:
        """Genes by HGNC approved symbol, previous symbol, alias, name or id, with GENCODE coordinates."""
        table = self._gene_table(query.str("ds", ""))
        if table is None:
            raise HttpError(HTTPStatus.NOT_FOUND, "No gene table (gtf_cache.feather) for this dataset; use regions instead")
        text = query.str("q", "")
        limit = query.int("limit", 25, lo=1, hi=200)
        try:
            gencode = self.gene_index.search(table, text, limit=limit)
        except ImportError as exc:
            raise HttpError(HTTPStatus.SERVICE_UNAVAILABLE, f"Reading the gene table needs pandas/pyarrow: {exc}")
        hgnc = self._hgnc()
        hits: List[Dict[str, Any]] = []
        if hgnc is not None:
            for rec, why in hgnc.search(text, limit):
                row = self.gene_index.locate(table, rec["symbol"], rec["ensembl_gene_id"])
                if row is None:
                    continue
                hits.append({**row, "symbol": rec["symbol"], "full_name": rec["name"], "hgnc_id": rec["hgnc_id"],
                             "locus_type": rec["locus_type"], "location": rec["location"], "matched": why})
            seen = {h["id"] for h in hits}
            hits += [g for g in gencode if g["id"] not in seen]
        else:
            hits = gencode
        return {"table": str(table), "genes": hits[:limit], "source": "HGNC + GENCODE" if hgnc is not None else "GENCODE", "hgnc": self.knowledge.hgnc_status()}

    def api_gene_info(self, query: Query, body: Any) -> Any:
        """HGNC record (names, aliases, groups, cross-references) and GENCODE coordinates of a gene."""
        text = query.str("symbol", required=True)
        ds = query.str("ds", "")
        hgnc = self._hgnc()
        rec = hgnc.get(text) if hgnc is not None else None
        table = self._gene_table(ds)
        row = None
        if table is not None:
            try:
                row = self.gene_index.locate(table, rec["symbol"] if rec else text, rec["ensembl_gene_id"] if rec else "")
            except ImportError:
                row = None
        in_dataset = []
        if ds:
            dataset = self.dataset(ds)
            in_dataset = [g for g in dataset.genes if g.upper() in {text.upper(), (rec or {}).get("symbol", "").upper(), (row or {}).get("name", "").upper()}]
        return {"query": text, "hgnc": public_record(rec) if rec else None, "gencode": row, "windows": in_dataset, "hgnc_status": self.knowledge.hgnc_status()}

    def _remote(self, fn: Callable[[], Any]) -> Any:
        try:
            return fn()
        except RemoteError as exc:
            raise HttpError(HTTPStatus.BAD_GATEWAY, str(exc))
        except KeyError as exc:
            raise HttpError(HTTPStatus.NOT_FOUND, str(exc).strip("'\""))
        except ValueError as exc:
            raise HttpError(HTTPStatus.BAD_REQUEST, str(exc))

    def api_gene_go(self, query: Query, body: Any) -> Any:
        symbol = query.str("symbol", required=True)
        return self._remote(lambda: self.knowledge.gene_go(symbol))

    def api_gene_names(self, query: Query, body: Any) -> Any:
        symbols = query.list("symbols")[:2000]
        hgnc = self._hgnc()
        return {"names": self.knowledge.gene_names(symbols) if hgnc is not None else {}, "hgnc": self.knowledge.hgnc_status()}

    def api_geneset_search(self, query: Query, body: Any) -> Any:
        """Gene sets to pick genes from: Gene Ontology terms (QuickGO) and HGNC gene groups."""
        text = query.str("q", "")
        sources = set(query.list("sources") or ["go", "hgnc"])
        out: Dict[str, Any] = {"go": [], "hgnc_groups": [], "errors": []}
        if "go" in sources:
            try:
                out["go"] = self.knowledge.go_search(text, limit=query.int("limit", 20, lo=1, hi=100))
            except RemoteError as exc:
                out["errors"].append(f"Gene Ontology: {exc}")
        if "hgnc" in sources:
            hgnc = self._hgnc()
            if hgnc is not None:
                out["hgnc_groups"] = hgnc.search_groups(text, limit=query.int("limit", 20, lo=1, hi=100))
            else:
                out["errors"].append(f"HGNC: {self.knowledge.hgnc_status().get('error')}")
        return out

    def api_geneset_genes(self, query: Query, body: Any) -> Any:
        """Genes of a GO term (and its descendants) or of an HGNC gene group, with coordinates."""
        source = query.str("source", "go")
        ident = query.str("id", required=True)
        ds = query.str("ds", "")
        if source == "go":
            data = self._remote(lambda: self.knowledge.go_genes(ident, include_regulation=query.str("regulation", "") in ("1", "true")))
            title = (data.get("term") or {}).get("name") or ident
            symbols = [(g["symbol"], g.get("terms") or []) for g in data["genes"]]
            extra = {"term": data.get("term"), "relations": data.get("relations")}
        elif source == "hgnc":
            symbols = [(sym, []) for sym in self._remote(lambda: self.knowledge.group_genes(ident))]
            title = ident
            extra = {}
        else:
            raise HttpError(HTTPStatus.BAD_REQUEST, "source must be go or hgnc")
        hgnc = self._hgnc()
        table = self._gene_table(ds)
        windows = set(g.upper() for g in self.dataset(ds).genes) if ds else set()
        genes = []
        for symbol, terms in symbols:
            rec = hgnc.get(symbol) if hgnc is not None else None
            row = self.gene_index.locate(table, symbol, (rec or {}).get("ensembl_gene_id", "")) if table is not None else None
            genes.append({
                "symbol": symbol,
                "full_name": (rec or {}).get("name", ""),
                "hgnc_id": (rec or {}).get("hgnc_id", ""),
                "locus_type": (rec or {}).get("locus_type", ""),
                "gencode": row,
                "in_dataset": symbol.upper() in windows or bool(row and row["name"].upper() in windows),
                "terms": terms,
            })
        return {"source": source, "id": ident, "title": title, "genes": genes, **extra}

    def api_ontology_term(self, query: Query, body: Any) -> Any:
        curie = query.str("curie", required=True)
        return self._remote(lambda: self.knowledge.ontology_term(curie))

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
        if self.alphagenome is None:
            return
        # Background predictions outlive the visualizer; keep a local AlphaGenome server they use.
        users = [t for t in (self.tasks.active() if self.tasks else []) if t["kind"] in ("predict", "import")]
        if users and self.alphagenome.settings.get("mode") == "local" and self.alphagenome.local.owned_running:
            print(f"Leaving the local AlphaGenome server running for {len(users)} background task(s) (pid {self.alphagenome.local._process.pid}).", file=sys.stderr)
            return
        self.alphagenome.local.stop()


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
    from genomics.workspace import DEFAULT_DATASET_DIR, DEFAULT_GENOTYPE_RUNS_ROOT, cache_path, results_path

    catalog = DatasetCatalog()
    for dataset_id in args.dataset_id:
        from genomics.core.data_registry import resolve_dataset

        ref = resolve_dataset(dataset_id)
        catalog.add(ref.path, dataset_id=dataset_id, annotations=args.annotations)
    for path in args.dataset:
        catalog.add(path, annotations=args.annotations)
    if not catalog.all() and Path(DEFAULT_DATASET_DIR, "dataset_metadata.json").exists():
        catalog.add(DEFAULT_DATASET_DIR, annotations=args.annotations)
    memory = None if args.no_add_datasets else DatasetMemory()
    for entry in memory.entries() if memory else []:
        path = Path(str(entry["path"]))
        if not (path / "dataset_metadata.json").exists():
            print(f"  skipping remembered dataset (missing): {path}", file=sys.stderr)
            continue
        try:
            catalog.add(path, annotations=Path(entry["annotations"]) if entry.get("annotations") else None)
        except Exception as exc:
            print(f"  skipping remembered dataset {path}: {exc}", file=sys.stderr)
    if not catalog.all():
        print("No dataset found; add one with --dataset PATH or from the UI.", file=sys.stderr)

    runs_roots = list(args.runs_root) or ([DEFAULT_GENOTYPE_RUNS_ROOT] if Path(DEFAULT_GENOTYPE_RUNS_ROOT).exists() else [])
    cache_dir = None if args.no_disk_cache else Path(args.cache_dir or cache_path("visualizer")).resolve()
    consensus_dirs = {}
    if args.consensus_dataset_dir:
        consensus_dirs = {d.id: str(Path(args.consensus_dataset_dir).resolve()) for d in catalog.all()}
    log_dir = (cache_dir / "logs") if cache_dir else None
    tasks_dir = Path(args.jobs_dir or results_path("visualizer", "jobs")).expanduser().resolve()
    default_config = args.pigmentation_config
    if default_config is None:
        from genomics.workspace import repo_root

        default_config = repo_root() / launch.DEFAULT_BASE_CONFIG
    return VisualizerApp(
        catalog,
        cache_dir=cache_dir,
        memory_bytes=max(args.memory_mb, 256) << 20,
        workers=max(1, args.workers),
        runs_roots=runs_roots,
        model_window=args.model_window,
        gtf=args.gtf,
        consensus_dirs=consensus_dirs,
        alphagenome=create_backend(args, log_dir),
        allow_add_datasets=not args.no_add_datasets,
        tasks_dir=tasks_dir,
        allow_tasks=not args.no_jobs,
        dataset_memory=memory,
        perturb_default_config=Path(default_config),
        scalar_store=ScalarStore.default(),
        session_store=SessionStore.default(),
        remote=not getattr(args, "no_remote", False),
    )


def ensure_tools_on_path() -> None:
    """bcftools/samtools from this interpreter's environment even when it was not activated."""
    env_bin = str(Path(sys.executable).parent)
    if shutil.which("bcftools") is None and (Path(env_bin) / "bcftools").exists():
        os.environ["PATH"] = os.pathsep.join([env_bin, os.environ.get("PATH", "")])


def bind_server(host: str, port: Optional[int]) -> Server:
    """Listen before the datasets load, so a busy port is reported at once. Without ``--port`` the
    next free port is used when the default one is taken."""
    factory = lambda candidate: Server((host, candidate), Handler)  # noqa: E731
    try:
        if port is not None:
            return factory(port)
        server, chosen = startup.bind_first_free(factory, DEFAULT_PORT, PORT_ATTEMPTS)
    except OSError as exc:
        if port is None:
            raise StartupError(f"Ports {DEFAULT_PORT}-{DEFAULT_PORT + PORT_ATTEMPTS - 1} are in use; choose one with --port") from exc
        if exc.errno == errno.EADDRINUSE:
            raise StartupError(f"Port {port} is in use by another program; choose another with --port") from exc
        raise StartupError(f"Cannot listen on {host}:{port}: {exc.strerror or exc}") from exc
    if chosen != DEFAULT_PORT:
        print(f"Port {DEFAULT_PORT} is in use; using {chosen}", file=sys.stderr)
    return server


def open_in_browser(url: str, wait: bool = False) -> None:
    if not startup.display_available():
        print(f"  no graphical display here: open {url} in a browser", file=sys.stderr)
        return
    if wait:
        webbrowser.open(url)
        return
    # webbrowser can block until the browser exits; the server is already listening.
    threading.Thread(target=webbrowser.open, args=(url,), name="open-browser", daemon=True).start()


def print_remote_hint(host: str, port: int) -> None:
    hint = startup.ssh_tunnel_hint(host, port)
    if hint:
        print(f"  over SSH: run `{hint}` on your computer, then open http://localhost:{port}/")


def reuse_running(args: argparse.Namespace) -> bool:
    """Point at a visualizer already running on the requested (or default) port when it has
    everything this invocation asks for; ``--open`` then opens it instead of failing."""
    port = args.port if args.port is not None else DEFAULT_PORT
    status = startup.probe_visualizer(args.host, port)
    if status is None or not startup.can_reuse(args, status):
        if status is not None and args.port is not None:
            raise StartupError(
                f"A visualizer at {startup.server_url(args.host, port)} is already using port {port} without the requested datasets; "
                "open them from its Overview page, stop it, or choose another --port"
            )
        return False
    url = startup.server_url(args.host, port)
    datasets = ", ".join(str(d.get("id")) for d in status.get("datasets") or []) or "no datasets"
    print(f"Genomics Visualizer already running: {url} ({datasets})")
    print_remote_hint(args.host, port)
    if args.open:
        open_in_browser(url, wait=True)
    return True


def serve(app: VisualizerApp, server: Server, open_browser: bool = False, verbose: bool = False) -> None:
    server.RequestHandlerClass = type("VisualizerHandler", (Handler,), {"app": app, "verbose": verbose})
    host, port = server.server_address[:2]
    url = startup.server_url(str(host), int(port))
    print(f"Genomics Visualizer: {url}")
    for dataset in app.catalog.all():
        print(f"  dataset {dataset.id}: {dataset.path} ({len(dataset.samples)} samples, {len(dataset.genes)} genes)")
    if app.cache_dir:
        print(f"  cache: {app.cache_dir}")
    if not app.remote.offline:
        # Warm the HGNC table (one ~17 MB download, then cached) so the first gene search is quick.
        threading.Thread(target=lambda: app._hgnc(), name="hgnc-prefetch", daemon=True).start()
    if app.tasks is not None:
        active = app.tasks.active()
        print(f"  background jobs: {app.tasks.root}" + (f" ({len(active)} running)" if active else ""))
    print_remote_hint(str(host), int(port))

    def stop(_signum=None, _frame=None):
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGTERM, stop)
    if open_browser:
        open_in_browser(url)
    try:
        server.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:
        pass
    finally:
        app.shutdown()
        server.server_close()


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        if reuse_running(args):
            return 0
        server = bind_server(args.host, args.port)
    except StartupError as exc:
        print(f"genomics visualize: {exc}", file=sys.stderr)
        return 1
    try:
        ensure_tools_on_path()
        app = create_app(args)
    except BaseException:
        server.server_close()
        raise
    serve(app, server, open_browser=args.open, verbose=args.verbose)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
