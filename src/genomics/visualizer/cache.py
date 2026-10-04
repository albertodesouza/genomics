"""Thread-safe, byte-bounded in-memory LRU cache and a small on-disk array cache."""
from __future__ import annotations

import hashlib
import json
import os
import threading
from collections import OrderedDict
from pathlib import Path
from typing import Any, Callable, Dict, Hashable, Optional, Tuple

import numpy as np


def _sizeof(value: Any) -> int:
    """Approximate retained size, good enough to bound the cache."""
    if isinstance(value, np.ndarray):
        return int(value.nbytes) + 112
    if isinstance(value, (bytes, bytearray, str)):
        return len(value) + 49
    if isinstance(value, (tuple, list)):
        return sum(_sizeof(item) for item in value) + 8 * len(value) + 56
    if isinstance(value, dict):
        return sum(_sizeof(item) for item in value.values()) + 100 * len(value) + 64
    if hasattr(value, "__dict__"):
        return sum(_sizeof(item) for item in vars(value).values()) + 64
    return 64


class LRUCache:
    """LRU cache bounded by an approximate byte budget.

    ``get_or_load`` deduplicates concurrent loads of the same key: the first caller loads while
    the others wait for its result, so a burst of identical requests does a single disk read.
    """

    def __init__(self, max_bytes: int):
        self.max_bytes = max(int(max_bytes), 0)
        self._items: "OrderedDict[Hashable, Tuple[Any, int]]" = OrderedDict()
        self._bytes = 0
        self._lock = threading.Lock()
        self._loading: Dict[Hashable, threading.Event] = {}
        self.hits = 0
        self.misses = 0

    def get(self, key: Hashable) -> Optional[Any]:
        with self._lock:
            item = self._items.get(key)
            if item is None:
                return None
            self._items.move_to_end(key)
            self.hits += 1
            return item[0]

    def put(self, key: Hashable, value: Any) -> None:
        size = _sizeof(value)
        if size > self.max_bytes:
            return
        with self._lock:
            old = self._items.pop(key, None)
            if old is not None:
                self._bytes -= old[1]
            self._items[key] = (value, size)
            self._bytes += size
            while self._bytes > self.max_bytes and self._items:
                _key, (_value, evicted) = self._items.popitem(last=False)
                self._bytes -= evicted

    def get_or_load(self, key: Hashable, loader: Callable[[], Any]) -> Any:
        while True:
            with self._lock:
                item = self._items.get(key)
                if item is not None:
                    self._items.move_to_end(key)
                    self.hits += 1
                    return item[0]
                pending = self._loading.get(key)
                if pending is None:
                    pending = threading.Event()
                    self._loading[key] = pending
                    self.misses += 1
                    owner = True
                else:
                    owner = False
            if not owner:
                pending.wait()
                continue
            try:
                value = loader()
                self.put(key, value)
                return value
            finally:
                with self._lock:
                    self._loading.pop(key, None)
                pending.set()

    def stats(self) -> Dict[str, int]:
        with self._lock:
            return {
                "items": len(self._items),
                "bytes": self._bytes,
                "max_bytes": self.max_bytes,
                "hits": self.hits,
                "misses": self.misses,
            }


def stable_key(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:32]


class DiskArrayCache:
    """Stores dicts of numpy arrays as ``.npz`` files keyed by a stable hash of their inputs."""

    def __init__(self, root: Optional[Path]):
        self.root = Path(root) if root else None

    def _path(self, namespace: str, key: str) -> Optional[Path]:
        if self.root is None:
            return None
        return self.root / namespace / f"{key}.npz"

    def load(self, namespace: str, key: str) -> Optional[Dict[str, np.ndarray]]:
        path = self._path(namespace, key)
        if path is None or not path.exists():
            return None
        try:
            with np.load(path, allow_pickle=False) as data:
                return {name: np.asarray(data[name]) for name in data.files}
        except Exception:
            return None

    def save(self, namespace: str, key: str, arrays: Dict[str, np.ndarray], compress: bool = False) -> None:
        path = self._path(namespace, key)
        if path is None:
            return
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(f"{path.stem}.{os.getpid()}.{threading.get_ident()}.tmp.npz")
            (np.savez_compressed if compress else np.savez)(tmp, **arrays)
            os.replace(tmp, path)
        except OSError:
            pass
