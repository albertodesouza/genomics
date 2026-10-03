"""Connect to AlphaGenome: the hosted API or a self-hosted gRPC server.

``alphagenome.models.dna_client.create`` only opens TLS channels trusted by the system CA store,
so it cannot talk to a self-hosted ``alphagenome_research`` ``server.py`` that runs in plaintext
or with a self-signed certificate. :func:`create_dna_client` builds the channel itself when a
server address is configured and otherwise defers to ``dna_client.create`` for the hosted API.

Configuration (explicit arguments win over the environment):

* ``ALPHAGENOME_ADDRESS`` -- self-hosted server, e.g. ``grpc://10.0.0.5:50051`` (plaintext),
  ``grpcs://10.0.0.5:50051`` (TLS) or ``10.0.0.5:50051`` (TLS detected automatically).
  Unset means the hosted API, which needs ``ALPHAGENOME_API_KEY`` (env or ``~/.env``).
* ``ALPHAGENOME_TLS_CA_CERT`` -- CA certificate (PEM) that signed a self-hosted server's TLS
  certificate (``certs/ca.crt`` from ``alphagenome_research/scripts/generate_certs.sh``).

grpc and alphagenome are imported lazily so this module stays cheap to import.
"""
from __future__ import annotations

import os
import socket
import ssl
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

ADDRESS_ENV = "ALPHAGENOME_ADDRESS"
CA_CERT_ENV = "ALPHAGENOME_TLS_CA_CERT"
API_KEY_ENV = "ALPHAGENOME_API_KEY"
DEFAULT_SERVER_PORT = 50051
# The self-hosted server ignores the key, but the client always sends one.
PLACEHOLDER_API_KEY = "self-hosted"

_TLS_SCHEMES = {"grpcs": "tls", "https": "tls", "grpc": "insecure", "http": "insecure"}


def resolve_api_key(explicit: Optional[str] = None) -> str:
    """Resolve the AlphaGenome API key from an explicit argument, the environment, or ``~/.env``,
    raising a clear error if none is found."""
    if explicit:
        return explicit
    api_key = os.environ.get(API_KEY_ENV)
    if api_key:
        return api_key
    try:
        from dotenv import dotenv_values

        api_key = dotenv_values(Path.home() / ".env").get(API_KEY_ENV)
    except ImportError:
        api_key = None
    if not api_key:
        raise RuntimeError(f"{API_KEY_ENV} not found in the environment, an explicit argument, or ~/.env.")
    return api_key


def api_key_available(explicit: Optional[str] = None) -> bool:
    try:
        resolve_api_key(explicit)
        return True
    except RuntimeError:
        return False


@dataclass(frozen=True)
class Endpoint:
    """A self-hosted server: ``target`` is ``host:port``; ``transport`` is ``"tls"``,
    ``"insecure"`` or ``"auto"`` (probe the server for TLS)."""

    target: str
    transport: str = "auto"
    ca_cert: Optional[str] = None

    @property
    def host_port(self) -> Tuple[str, int]:
        host, _, port = self.target.rpartition(":")
        return host.strip("[]"), int(port)

    @property
    def url(self) -> str:
        scheme = {"tls": "grpcs://", "insecure": "grpc://"}.get(self.transport, "")
        return scheme + self.target


def parse_address(address: str, ca_cert: Optional[str] = None) -> Endpoint:
    """Parse ``[grpc|grpcs|http|https://]host[:port][/]`` into an :class:`Endpoint`.

    A missing port defaults to 50051 (443 for ``https``). Without a scheme the transport is
    probed on first use.
    """
    text = (address or "").strip()
    if not text:
        raise ValueError("empty AlphaGenome server address")
    transport = "auto"
    if "://" in text:
        scheme, text = text.split("://", 1)
        scheme = scheme.lower()
        if scheme not in _TLS_SCHEMES:
            raise ValueError(f"unsupported scheme {scheme!r} (use grpc://, grpcs://, http:// or https://)")
        transport = _TLS_SCHEMES[scheme]
    text = text.split("/", 1)[0]
    if not text:
        raise ValueError(f"no host in AlphaGenome server address {address!r}")
    if text.startswith("["):  # IPv6 literal
        host, _, rest = text[1:].partition("]")
        port_text = rest[1:] if rest.startswith(":") else ""
        host = f"[{host}]"
    elif text.count(":") == 1:
        host, port_text = text.split(":")
    else:
        host, port_text = text, ""
    if port_text:
        if not port_text.isdigit() or not 0 < int(port_text) < 65536:
            raise ValueError(f"invalid port in AlphaGenome server address {address!r}")
        port = int(port_text)
    else:
        port = 443 if address.strip().lower().startswith("https://") else DEFAULT_SERVER_PORT
    ca = str(Path(ca_cert).expanduser()) if ca_cert else None
    return Endpoint(f"{host}:{port}", transport, ca)


def endpoint_from_env(address: Optional[str] = None, ca_cert: Optional[str] = None) -> Optional[Endpoint]:
    """The configured self-hosted endpoint, or ``None`` for the hosted API."""
    address = address or os.environ.get(ADDRESS_ENV)
    if not address:
        return None
    return parse_address(address, ca_cert or os.environ.get(CA_CERT_ENV) or None)


def backend_configured(api_key: Optional[str] = None) -> bool:
    """True when predictions can be requested: a self-hosted server is set, or a key exists."""
    return bool(os.environ.get(ADDRESS_ENV)) or api_key_available(api_key)


def speaks_tls(host: str, port: int, timeout: float = 3.0) -> bool:
    """Whether ``host:port`` completes a TLS handshake (certificate not verified)."""
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    context.set_alpn_protocols(["h2"])
    try:
        with socket.create_connection((host, port), timeout=timeout) as sock:
            with context.wrap_socket(sock, server_hostname=host):
                return True
    except (ssl.SSLError, OSError):
        return False


def resolve_transport(endpoint: Endpoint, timeout: float = 3.0) -> Endpoint:
    if endpoint.transport != "auto":
        return endpoint
    host, port = endpoint.host_port
    transport = "tls" if speaks_tls(host, port, timeout) else "insecure"
    return Endpoint(endpoint.target, transport, endpoint.ca_cert)


def open_channel(endpoint: Endpoint, timeout: Optional[float] = 10.0) -> Tuple[Any, Endpoint]:
    """Open a ready gRPC channel to ``endpoint`` (resolving ``auto`` transport first)."""
    import grpc

    endpoint = resolve_transport(endpoint, timeout=min(timeout or 3.0, 3.0))
    options = [("grpc.max_send_message_length", -1), ("grpc.max_receive_message_length", -1)]
    if endpoint.transport == "tls":
        root = Path(endpoint.ca_cert).read_bytes() if endpoint.ca_cert else None
        channel = grpc.secure_channel(endpoint.target, grpc.ssl_channel_credentials(root_certificates=root), options=options)
    else:
        channel = grpc.insecure_channel(endpoint.target, options=options)
    try:
        grpc.channel_ready_future(channel).result(timeout=timeout)
    except grpc.FutureTimeoutError:
        channel.close()
        raise ConnectionError(_unreachable_message(endpoint, timeout)) from None
    return channel, endpoint


def _unreachable_message(endpoint: Endpoint, timeout: Optional[float]) -> str:
    message = f"could not connect to {endpoint.url} within {timeout:g}s"
    if endpoint.transport == "tls":
        hint = "check that the CA certificate signed the server certificate and that it lists this host"
        message += f" over TLS ({hint})" if endpoint.ca_cert else " over TLS (self-signed server? set the CA certificate)"
    return message


def create_dna_client(
    api_key: Optional[str] = None,
    address: Optional[str] = None,
    ca_cert: Optional[str] = None,
    timeout: Optional[float] = None,
) -> Any:
    """``DnaClient`` for the configured self-hosted server, else for the hosted API."""
    from alphagenome.models import dna_client

    endpoint = endpoint_from_env(address, ca_cert)
    if endpoint is None:
        if timeout is None:
            return dna_client.create(api_key=resolve_api_key(api_key))
        return dna_client.create(api_key=resolve_api_key(api_key), timeout=timeout)
    channel, _ = open_channel(endpoint, timeout=timeout if timeout is not None else 30.0)
    key = api_key or os.environ.get(API_KEY_ENV) or PLACEHOLDER_API_KEY
    return dna_client.DnaClient(channel=channel, metadata=[("x-goog-api-key", key)])


def check_connection(
    address: Optional[str] = None,
    ca_cert: Optional[str] = None,
    api_key: Optional[str] = None,
    timeout: float = 10.0,
    predict: bool = False,
) -> Dict[str, Any]:
    """Reach the server and call the AlphaGenome service.

    Calls ``GetMetadata`` (cheap); with ``predict`` also runs a 16 kb ``predict_sequence``, which
    proves the model itself runs (the first call on a fresh server includes JIT compilation).
    Never raises: returns ``{"ok": bool, "message": str, ...}``.
    """
    import grpc

    started = time.perf_counter()
    result: Dict[str, Any] = {"ok": False, "target": "AlphaGenome API (hosted)", "transport": "tls"}
    try:
        endpoint = parse_address(address, ca_cert) if address else None
    except ValueError as exc:
        return {**result, "message": str(exc)}
    if endpoint is not None:
        result["target"] = endpoint.url
        if endpoint.ca_cert and not Path(endpoint.ca_cert).is_file():
            return {**result, "message": f"CA certificate not found: {endpoint.ca_cert}"}
    try:
        if endpoint is None:
            key = resolve_api_key(api_key)
            channel, _ = open_channel(Endpoint("gdmscience.googleapis.com:443", "tls"), timeout=timeout)
        else:
            key = api_key or os.environ.get(API_KEY_ENV) or PLACEHOLDER_API_KEY
            channel, endpoint = open_channel(endpoint, timeout=timeout)
            result.update(target=endpoint.url, transport=endpoint.transport)
    except (RuntimeError, ConnectionError, OSError) as exc:
        return {**result, "message": str(exc)}
    result["connect_ms"] = round((time.perf_counter() - started) * 1000)
    try:
        from alphagenome.models import dna_client
        from alphagenome.protos import dna_model_service_pb2, dna_model_service_pb2_grpc

        metadata = [("x-goog-api-key", key)]
        stub = dna_model_service_pb2_grpc.DnaModelServiceStub(channel)
        request = dna_model_service_pb2.MetadataRequest(organism=dna_client.Organism.HOMO_SAPIENS.to_proto())
        responses = stub.GetMetadata(request, metadata=metadata, timeout=timeout)
        try:
            next(iter(responses), None)
        finally:
            responses.cancel()
        result["rpc_ms"] = round((time.perf_counter() - started) * 1000)
        if predict:
            t0 = time.perf_counter()
            client = dna_client.DnaClient(channel=channel, metadata=metadata)
            output = client.predict_sequence(
                "ACGT" * (2 ** 14 // 4),
                requested_outputs=[dna_client.OutputType.RNA_SEQ],
                ontology_terms=["UBERON:0002107"],
            )
            shape = list(getattr(getattr(output, "rna_seq", None), "values", []).shape)
            result.update(predict_ms=round((time.perf_counter() - t0) * 1000), predict_shape=shape)
    except grpc.RpcError as exc:
        code = exc.code() if hasattr(exc, "code") else None
        detail = exc.details() if hasattr(exc, "details") else str(exc)
        hints = {
            grpc.StatusCode.UNAUTHENTICATED: "the server rejected the API key",
            grpc.StatusCode.PERMISSION_DENIED: "the server rejected the API key",
            grpc.StatusCode.UNIMPLEMENTED: "the server is reachable but does not serve AlphaGenome",
            grpc.StatusCode.DEADLINE_EXCEEDED: "the server did not answer in time",
        }
        name = code.name if code is not None else "error"
        return {**result, "message": f"{hints.get(code, 'request failed')} ({name}: {detail})"}
    except Exception as exc:  # predict_sequence on a misbehaving server
        return {**result, "message": f"request failed: {exc}"}
    finally:
        channel.close()
    result["ok"] = True
    parts = [f"connected to {result['target']} in {result['connect_ms']} ms"]
    if "predict_ms" in result:
        parts.append(f"test prediction {result['predict_shape']} in {result['predict_ms'] / 1000:.1f} s")
    result["message"] = "; ".join(parts)
    return result
