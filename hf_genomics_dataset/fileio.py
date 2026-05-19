from __future__ import annotations

import gzip
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


def load_vcf_variants(path: Path, sample_id: Optional[str] = None, region: Optional[str] = None) -> str:
    """Read a .vcf or .vcf.gz file and return only the variant lines (skipping headers).
    If sample_id is provided, filters to keep only sites where the sample has a variant allele.
    If region is provided, uses bcftools to extract only that region (e.g. 'chr1:100-200').
    """
    import subprocess

    if not path.exists():
        return ""

    # Try using bcftools for efficient filtering if sample_id is provided
    if sample_id:
        try:
            # -H: no header
            # -s: specific sample
            # -c 1: minimum allele count 1 (keep only variants for this sample)
            cmd = ["bcftools", "view", "-H"]
            if sample_id:
                cmd.extend(["-s", sample_id, "-c", "1"])
            if region:
                cmd.append("-r" if ":" in region else "-R")
                cmd.append(region)
            
            cmd.append(str(path))
            
            proc = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
            )
            return proc.stdout.strip()
        except Exception as e:
            # Fallback to manual reading if bcftools fails
            pass

    lines: List[str] = []
    # Handle gzipped files
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"

    try:
        with opener(path, mode, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                lines.append(line.strip())
    except Exception:
        return ""

    return "\n".join(lines)


def detect_vcf_chr_prefix(vcf_path: Path) -> str:
    """Detect whether VCF uses 'chr' prefix for contigs in its header."""
    import subprocess

    try:
        proc = subprocess.run(
            ["bcftools", "view", "-h", str(vcf_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        for line in proc.stdout.splitlines():
            if line.startswith("##contig=<ID=chr"):
                return "chr"
            if line.startswith("##contig=<ID=1,") or line.startswith("##contig=<ID=1>"):
                return ""
    except Exception:
        pass

    if ".chr" in vcf_path.name:
        return "chr"

    return ""
