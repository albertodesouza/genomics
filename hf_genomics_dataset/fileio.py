from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def load_npz_matrix(path: Path) -> List[List[float]]:
    with np.load(path) as data:
        if "values" in data:
            array = data["values"]
        else:
            keys = list(data.keys())
            if not keys:
                raise ValueError(f"No arrays found in {path}")
            array = data[keys[0]]

    array = np.asarray(array, dtype=np.float32)
    if array.ndim == 1:
        array = array.reshape(-1, 1)
    elif array.ndim != 2:
        raise ValueError(f"Unsupported array rank {array.ndim} in {path}")

    return array.tolist()


def load_fasta_sequence(path: Path) -> str:
    sequence_lines: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_lines.append(line)
    return "".join(sequence_lines)


def stable_json_digest(payload: Dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
