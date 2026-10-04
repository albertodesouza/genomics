"""HTTP GETs to public bioinformatics services with an on-disk response cache (stdlib only).

Used for HGNC, QuickGO, OLS, ENCODE and FANTOM5 lookups. Responses are cached under
``<cache>/remote/`` for ``ttl`` seconds (a stale copy is still served when the service is
unreachable), so the visualizer keeps working offline once something has been looked up.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional

USER_AGENT = "genomics-visualizer/1.0 (python urllib)"
DAY = 24 * 3600


class RemoteError(RuntimeError):
    pass


class RemoteCache:
    def __init__(self, root: Optional[Path], timeout: float = 60.0, offline: bool = False):
        self.root = Path(root) / "remote" if root else None
        self.timeout = timeout
        self.offline = offline  # serve cached responses only (``genomics visualize --no-remote``)
        self._locks: Dict[str, threading.Lock] = {}
        self._guard = threading.Lock()

    def _lock(self, key: str) -> threading.Lock:
        with self._guard:
            return self._locks.setdefault(key, threading.Lock())

    def path(self, url: str, suffix: str = "") -> Optional[Path]:
        if self.root is None:
            return None
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()[:32]
        return self.root / f"{digest}{suffix}"

    def get_bytes(self, url: str, ttl: float = 7 * DAY, headers: Optional[Dict[str, str]] = None, params: Optional[Dict[str, Any]] = None) -> bytes:
        if params:
            url = f"{url}{'&' if '?' in url else '?'}{urllib.parse.urlencode({k: v for k, v in params.items() if v is not None}, doseq=True)}"
        path = self.path(url, ".bin.gz")
        with self._lock(url):
            if path is not None and path.exists() and (self.offline or time.time() - path.stat().st_mtime < ttl):
                return gzip.decompress(path.read_bytes())
            if self.offline:
                raise RemoteError(f"{urllib.parse.urlsplit(url).netloc}: not cached and remote lookups are disabled (--no-remote)")
            try:
                data = self._download(url, headers or {})
            except RemoteError:
                if path is not None and path.exists():  # stale but usable offline
                    return gzip.decompress(path.read_bytes())
                raise
            if path is not None:
                try:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
                    tmp.write_bytes(gzip.compress(data, compresslevel=5))
                    os.replace(tmp, path)
                except OSError:
                    pass
            return data

    def get_json(self, url: str, ttl: float = 7 * DAY, params: Optional[Dict[str, Any]] = None) -> Any:
        data = self.get_bytes(url, ttl, headers={"Accept": "application/json"}, params=params)
        try:
            return json.loads(data.decode("utf-8"))
        except ValueError as exc:
            raise RemoteError(f"{url}: not JSON ({exc})") from exc

    def get_text(self, url: str, ttl: float = 7 * DAY, params: Optional[Dict[str, Any]] = None, headers: Optional[Dict[str, str]] = None) -> str:
        return self.get_bytes(url, ttl, headers=headers, params=params).decode("utf-8", "replace")

    def _download(self, url: str, headers: Dict[str, str]) -> bytes:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **headers})
        last: Optional[Exception] = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    data = response.read()
                    if response.headers.get("Content-Encoding") == "gzip":
                        data = gzip.decompress(data)
                    return data
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    # ENCODE answers an empty search with 404 and a JSON body.
                    body = exc.read()
                    if body.startswith(b"{"):
                        return body
                    raise RemoteError(f"{url}: HTTP 404") from exc
                last = exc
                if exc.code < 500 and exc.code != 429:
                    break
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last = exc
            time.sleep(0.5 * (attempt + 1))
        raise RemoteError(f"{urllib.parse.urlsplit(url).netloc} unreachable: {last}")
