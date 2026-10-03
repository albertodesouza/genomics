"""Argument definitions for ``genomics visualize`` (stdlib only, so ``--help`` stays light)."""
from __future__ import annotations

import argparse
import os
from pathlib import Path


def add_visualizer_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    data = parser.add_argument_group("data")
    data.add_argument("--dataset", action="append", default=[], type=Path, metavar="DIR", help="Dataset directory with dataset_metadata.json (repeatable)")
    data.add_argument("--dataset-id", action="append", default=[], metavar="ID", help="Registered dataset id (repeatable), e.g. 1kg_high_coverage")
    data.add_argument("--annotations", type=Path, default=None, metavar="TABLE", help="CSV/TSV of extra sample columns (first column = sample id) used as facets")
    data.add_argument("--consensus-dataset-dir", type=Path, default=None, metavar="DIR", help="Consensus dataset for the training alignment axis (default: the dataset itself)")
    data.add_argument("--runs-root", action="append", default=[], type=Path, metavar="DIR", help="Experiment runs directory (repeatable; default: results/genotype_based_predictor/runs)")
    data.add_argument("--gtf", type=Path, default=None, metavar="TABLE", help="GTF feature table (feather/parquet) for gene models (default: <dataset>/gtf_cache.feather)")

    perf = parser.add_argument_group("performance")
    perf.add_argument("--cache-dir", type=Path, default=None, metavar="DIR", help="On-disk cache for cohort aggregates and indel indexes (default: results/cache/visualizer)")
    perf.add_argument("--no-disk-cache", action="store_true", help="Keep aggregates in memory only")
    perf.add_argument("--memory-mb", type=int, default=4096, help="In-memory cache budget in MB (default: 4096)")
    perf.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 4), help="Threads for interactive requests (bulk jobs use up to 16)")
    perf.add_argument("--model-window", type=int, default=32768, help="CNN training window size for the training-axis coordinates (default: 32768)")

    labs = parser.add_argument_group("labs")
    labs.add_argument("--pigmentation-config", type=Path, default=None, help="Config for the optional Pigmentation Sequence Lab")
    labs.add_argument("--lab-port", type=int, default=8781, help="Port used when the Pigmentation Sequence Lab is launched")
    labs.add_argument("--alphagenome-address", default=None, metavar="URL", help="Self-hosted AlphaGenome server for the Labs, e.g. grpc://host:50051 or grpcs://host:50051 (overrides the setting saved from the UI)")
    labs.add_argument("--alphagenome-ca-cert", type=Path, default=None, metavar="PEM", help="CA certificate for a self-hosted AlphaGenome server using TLS")
    labs.add_argument("--alphagenome-server-dir", type=Path, default=None, metavar="DIR", help="alphagenome_research checkout with server.py, for 'Start server on this machine' (default: $ALPHAGENOME_SERVER_DIR or ../alphagenome_research)")
    labs.add_argument("--alphagenome-server-python", type=Path, default=None, metavar="PYTHON", help="Interpreter with alphagenome_research + JAX for the local server (default: $ALPHAGENOME_SERVER_PYTHON, else the first conda env with alphagenome_research and CUDA jax)")

    server = parser.add_argument_group("server")
    server.add_argument("--host", default="127.0.0.1", help="Bind address (default: 127.0.0.1)")
    server.add_argument("--port", type=int, default=8780, help="Port (default: 8780)")
    server.add_argument("--open", action="store_true", help="Open the browser after starting")
    server.add_argument("--no-add-datasets", action="store_true", help="Disallow opening other dataset paths from the UI")
    server.add_argument("--verbose", action="store_true", help="Log every request")
    return parser


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="genomics visualize", description="Interactive genomics dataset and AlphaGenome prediction visualizer")
    return add_visualizer_arguments(parser)
