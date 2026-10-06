"""Saved sessions: a named snapshot of what a page is showing, to return to or hand to someone else.

A session holds the locus and track settings (gene, coordinates, mode, tracks and their order, group
field, y-scale, compare-with lanes), the cohort filters and the pinned individuals — everything the
pages already keep per dataset in the browser. Saving it moves that state out of one browser's
local storage into a file beside the region scalars, so a view survives a cleared cache, can be
opened on another machine and can be listed for a collaborator.

The payload is stored as the pages produce it (an opaque blob to this module, bounded in size); only
the name, the description and the bookkeeping are interpreted here, so adding a page setting needs no
change on the server. Stored per dataset path in ``~/.config/genomics/visualizer_sessions.json``.
"""
from __future__ import annotations

import json
import re
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from genomics.visualizer.datasets import Dataset, load_json

MAX_SESSIONS = 100
MAX_PAYLOAD = 256 * 1024  # a session is settings, not data
NAME_RE = re.compile(r"^[^\x00-\x1f]{1,80}$")


class SessionError(ValueError):
    """Bad name, oversized payload, or too many saved sessions."""


class SessionStore:
    """Saved sessions per dataset path, persisted to a JSON file (or in memory without a path)."""

    def __init__(self, path: Optional[Path] = None):
        self.path = Path(path) if path is not None else None
        self._memory: Dict[str, List[Dict[str, Any]]] = {}
        self._lock = threading.Lock()

    @classmethod
    def default(cls) -> "SessionStore":
        import os

        base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
        return cls(Path(base) / "genomics" / "visualizer_sessions.json")

    def _read(self) -> Dict[str, List[Dict[str, Any]]]:
        if self.path is None:
            return self._memory
        try:
            data = load_json(self.path)
        except (OSError, ValueError):
            return {}
        items = data.get("datasets") if isinstance(data, dict) else None
        return {str(k): [s for s in v if isinstance(s, dict)] for k, v in (items or {}).items() if isinstance(v, list)}

    def _write(self, data: Dict[str, List[Dict[str, Any]]]) -> None:
        if self.path is None:
            self._memory = data
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"datasets": data}, indent=2) + "\n", encoding="utf-8")
        tmp.replace(self.path)

    def list(self, dataset: Dataset) -> List[Dict[str, Any]]:
        """Saved sessions of a dataset, most recently saved first."""
        with self._lock:
            items = [dict(s) for s in self._read().get(str(dataset.path), [])]
        return sorted(items, key=lambda s: s.get("saved_at") or 0, reverse=True)

    def get(self, dataset: Dataset, name: str) -> Optional[Dict[str, Any]]:
        for item in self.list(dataset):
            if item.get("name") == name:
                return item
        return None

    def save(self, dataset: Dataset, name: str, payload: Dict[str, Any], description: str = "", replace: bool = False) -> Dict[str, Any]:
        name = str(name or "").strip()
        if not NAME_RE.match(name):
            raise SessionError("A session needs a name of at most 80 characters")
        if not isinstance(payload, dict):
            raise SessionError("The session state must be an object")
        encoded = json.dumps(payload)
        if len(encoded) > MAX_PAYLOAD:
            raise SessionError(f"This session is {len(encoded) // 1024} kB; at most {MAX_PAYLOAD // 1024} kB is stored")
        entry = {"name": name, "description": str(description or "").strip()[:400],
                 "saved_at": time.time(), "state": json.loads(encoded)}
        with self._lock:
            data = self._read()
            items = data.get(str(dataset.path), [])
            if any(s.get("name") == name for s in items) and not replace:
                raise SessionError(f"A session called {name!r} already exists; save it under another name or replace it")
            items = [s for s in items if s.get("name") != name]
            if len(items) >= MAX_SESSIONS:
                raise SessionError(f"At most {MAX_SESSIONS} saved sessions per dataset; delete one first")
            data[str(dataset.path)] = items + [entry]
            self._write(data)
        return entry

    def rename(self, dataset: Dataset, name: str, new_name: str) -> Dict[str, Any]:
        new_name = str(new_name or "").strip()
        if not NAME_RE.match(new_name):
            raise SessionError("A session needs a name of at most 80 characters")
        with self._lock:
            data = self._read()
            items = data.get(str(dataset.path), [])
            entry = next((s for s in items if s.get("name") == name), None)
            if entry is None:
                raise SessionError(f"No saved session called {name!r}")
            if new_name != name and any(s.get("name") == new_name for s in items):
                raise SessionError(f"A session called {new_name!r} already exists")
            entry["name"] = new_name
            data[str(dataset.path)] = items
            self._write(data)
        return dict(entry)

    def delete(self, dataset: Dataset, name: str) -> bool:
        with self._lock:
            data = self._read()
            items = data.get(str(dataset.path), [])
            kept = [s for s in items if s.get("name") != name]
            if len(kept) == len(items):
                return False
            data[str(dataset.path)] = kept
            self._write(data)
        return True
