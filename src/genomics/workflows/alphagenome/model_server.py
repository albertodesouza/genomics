"""Serve the AlphaGenome model from an ``alphagenome_research`` checkout over the public gRPC API.

This file runs **inside the AlphaGenome server environment** (JAX + ``alphagenome_research``), not
inside the ``genomics`` environment, so it imports nothing from ``genomics`` and only the standard
library at module level. Start it with ``genomics alphagenome server start`` (or *Start server* on the
visualizer's AlphaGenome page), which runs::

    <server python> .../model_server.py --server-dir <alphagenome_research checkout> --port 50051

It reuses the ``AlphaGenomeServer`` servicer of the checkout's ``server.py``
(https://github.com/FeLiPeOLi7/alphagenome_research) and changes what running ``server.py`` directly
cannot do:

* ``GetMetadata`` returns the model's track metadata (``server.py`` answers with an empty message, so
  ``DnaClient.output_metadata()`` -- the visualizer's track catalog and tissue picker -- came back empty);
* bind address and port are options (``server.py`` always listens on ``0.0.0.0:50051``; the default
  here is ``127.0.0.1``, reachable only from this machine);
* cached Hugging Face weights load without network access or a login prompt (``create_from_huggingface``
  calls ``whoami()`` first and opens an interactive ``login()`` without a token), and ``--checkpoint``
  loads a local checkpoint (e.g. downloaded from Kaggle);
* it refuses to start on CPU-only JAX unless ``--allow-cpu`` is given, instead of loading for minutes
  and serving predictions that take hours.

Clients are the official SDK: ``dna_client.create(api_key="local", address="127.0.0.1:50051")`` for a
plaintext server, or :func:`genomics.core.alphagenome_connection.create_dna_client`.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 50051
DEFAULT_MODEL_VERSION = "all_folds"
READY_MARKER = "ALPHAGENOME_SERVER_READY"


def log(message: str) -> None:
    print(f"[model-server {time.strftime('%H:%M:%S')}] {message}", flush=True)


def hf_repo_id(model_version: str) -> str:
    return f"google/alphagenome-{model_version.replace('_', '-').lower()}"


def resolve_checkpoint(model_version: str, checkpoint: Optional[str]) -> str:
    """A local checkpoint directory: ``checkpoint`` as given, the Hugging Face cache, or a download."""
    if checkpoint:
        path = Path(checkpoint).expanduser()
        if not path.is_dir():
            raise SystemExit(f"checkpoint directory not found: {path}")
        return str(path)
    import huggingface_hub

    repo_id = hf_repo_id(model_version)
    try:
        path = huggingface_hub.snapshot_download(repo_id=repo_id, local_files_only=True)
        log(f"weights: {repo_id} from the Hugging Face cache ({path})")
        return path
    except Exception:  # not cached (LocalEntryNotFoundError and friends)
        pass
    log(f"weights: downloading {repo_id} from Hugging Face (~700 MB, once)")
    try:
        return huggingface_hub.snapshot_download(repo_id=repo_id)
    except Exception as exc:
        raise SystemExit(
            f"could not download {repo_id}: {type(exc).__name__}: {exc}\n"
            f"The weights are gated: accept the terms at https://huggingface.co/{repo_id}, then log in once "
            "with `hf auth login` (or set HF_TOKEN) in this environment. A checkpoint downloaded elsewhere "
            "(e.g. from Kaggle) can be passed with --checkpoint DIR."
        ) from None


def check_device(allow_cpu: bool) -> Any:
    import jax

    devices = jax.devices()
    gpus = [d for d in devices if d.platform in ("gpu", "cuda", "rocm", "tpu")]
    if gpus:
        log(f"device: {gpus[0].device_kind} ({gpus[0].platform}); jax {jax.__version__}")
        return gpus[0]
    message = (
        f"jax {jax.__version__} sees no GPU (devices: {', '.join(str(d) for d in devices)}). Install a CUDA-enabled jax "
        "in this environment (`genomics alphagenome server setup`, or pip install -U \"jax[cuda12]\")"
    )
    if not allow_cpu:
        raise SystemExit(message + "; --allow-cpu runs on CPU anyway (very slow).")
    log("WARNING: " + message + "; running on CPU because of --allow-cpu")
    return devices[0]


def metadata_response(model: Any, organism: Any) -> Any:
    """``MetadataResponse`` with the model's (padding-free) metadata for every output, as the hosted API sends."""
    from alphagenome.models import dna_output, junction_data_utils, track_data_utils
    from alphagenome.protos import dna_model_pb2, dna_model_service_pb2

    output_metadata = model.output_metadata(organism)
    response = dna_model_service_pb2.MetadataResponse()
    for output_type in dna_output.OutputType:
        frame = output_metadata.get(output_type)
        if frame is None:
            continue
        proto = dna_model_pb2.OutputMetadata(output_type=output_type.to_proto())
        if output_type == dna_output.OutputType.SPLICE_JUNCTIONS:
            proto.junctions.CopyFrom(junction_data_utils.metadata_to_proto(frame))
        else:
            proto.tracks.CopyFrom(track_data_utils.metadata_to_proto(frame))
        response.output_metadata.append(proto)
    return response


def build_servicer(server_module: Any, model: Any) -> Any:
    base = server_module.AlphaGenomeServer

    class Servicer(base):  # type: ignore[misc, valid-type]
        def __init__(self, model: Any) -> None:
            super().__init__(model=model)
            self._metadata_cache: Dict[Any, Any] = {}

        def GetMetadata(self, request: Any, context: Any) -> Any:  # noqa: N802 (gRPC method name)
            organism = self._parse_organism(request.organism)
            if organism not in self._metadata_cache:
                self._metadata_cache[organism] = metadata_response(self.model, organism)
            yield self._metadata_cache[organism]

    return Servicer(model)


def tls_credentials(server_dir: Path, cert: Optional[str], key: Optional[str], plaintext: bool) -> Optional[tuple]:
    """(key bytes, cert bytes) when TLS is configured: explicit files, else ``certs/server.{crt,key}`` like server.py."""
    if plaintext:
        return None
    cert_path = Path(cert).expanduser() if cert else server_dir / os.environ.get("ALPHAGENOME_TLS_CERT", "certs/server.crt")
    key_path = Path(key).expanduser() if key else server_dir / os.environ.get("ALPHAGENOME_TLS_KEY", "certs/server.key")
    if cert or key:
        for path in (cert_path, key_path):
            if not path.is_file():
                raise SystemExit(f"TLS file not found: {path}")
    if cert_path.is_file() and key_path.is_file():
        return key_path.read_bytes(), cert_path.read_bytes()
    return None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Serve AlphaGenome (alphagenome_research) over the public gRPC API")
    parser.add_argument("--server-dir", type=Path, default=Path.cwd(), help="alphagenome_research checkout with server.py (default: current directory)")
    parser.add_argument("--host", default=os.environ.get("ALPHAGENOME_SERVER_HOST", DEFAULT_HOST), help=f"Bind address (default: {DEFAULT_HOST}; 0.0.0.0 serves the whole network, without authentication)")
    parser.add_argument("--port", type=int, default=int(os.environ.get("ALPHAGENOME_SERVER_PORT", DEFAULT_PORT)), help=f"Port (default: {DEFAULT_PORT})")
    parser.add_argument("--model-version", default=DEFAULT_MODEL_VERSION, help="Hugging Face model version: all_folds (default) or fold_0 ... fold_3")
    parser.add_argument("--checkpoint", default=None, metavar="DIR", help="Local checkpoint directory instead of Hugging Face")
    parser.add_argument("--tls-cert", default=None, metavar="PEM", help="Server certificate (default: certs/server.crt in the checkout, when present)")
    parser.add_argument("--tls-key", default=None, metavar="PEM", help="Server private key (default: certs/server.key in the checkout, when present)")
    parser.add_argument("--plaintext", action="store_true", help="Do not use TLS even when certs/server.crt and certs/server.key exist")
    parser.add_argument("--max-workers", type=int, default=10, help="gRPC worker threads (default: 10)")
    parser.add_argument("--allow-cpu", action="store_true", help="Start even when jax has no GPU (predictions take minutes each)")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    server_dir = args.server_dir.expanduser().resolve()
    if not (server_dir / "server.py").is_file():
        raise SystemExit(f"server.py not found in {server_dir} (pass the alphagenome_research checkout with --server-dir)")
    os.chdir(server_dir)
    sys.path.insert(0, str(server_dir))
    # Same settings server.py applies when imported, but they must be in place before jax is first
    # imported (check_device below): otherwise XLA preallocates 75% of GPU memory -- ~90 GB of the
    # unified memory of a DGX Spark, which the visualizer and training also need.
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
    os.environ.setdefault("HF_HUB_DISABLE_XET", "1")
    credentials = tls_credentials(server_dir, args.tls_cert, args.tls_key, args.plaintext)

    started = time.perf_counter()
    device = check_device(args.allow_cpu)
    import server as server_module  # the checkout's server.py: sets XLA env vars, defines AlphaGenomeServer
    from alphagenome_research.model import dna_model

    import grpc
    from alphagenome.protos import dna_model_service_pb2_grpc

    checkpoint = resolve_checkpoint(args.model_version, args.checkpoint)
    log("loading the model (1-3 minutes) ...")
    model = dna_model.create(checkpoint, device=device)
    servicer = build_servicer(server_module, model)

    max_message = getattr(server_module, "MAX_MESSAGE_LENGTH", 100 * 1024 * 1024)
    from concurrent import futures

    server = grpc.server(
        futures.ThreadPoolExecutor(max_workers=args.max_workers),
        options=[("grpc.max_receive_message_length", max_message), ("grpc.max_send_message_length", max_message)],
    )
    dna_model_service_pb2_grpc.add_DnaModelServiceServicer_to_server(servicer, server)
    target = f"{args.host}:{args.port}"
    if credentials is not None:
        bound = server.add_secure_port(target, grpc.ssl_server_credentials((credentials,)))
        scheme = "grpcs"
    else:
        bound = server.add_insecure_port(target)
        scheme = "grpc"
    if not bound:
        raise SystemExit(f"could not listen on {target} (port in use?)")
    server.start()
    shown = "127.0.0.1" if args.host in ("0.0.0.0", "::", "") else args.host
    log(f"model loaded in {time.perf_counter() - started:.0f} s")
    log(f"{READY_MARKER} listening on {target} ({'TLS' if credentials else 'plaintext'}); clients use {scheme}://{shown}:{args.port}")
    try:
        server.wait_for_termination()
    except KeyboardInterrupt:
        log("stopping")
        server.stop(grace=5)
    return 0


if __name__ == "__main__":
    sys.exit(main())
