"""Minimal bigWig reader (stdlib + numpy) for local files and remote URLs (HTTP Range requests).

Only what the visualizer needs: chromosome sizes and per-base values over a region, read from the
full-resolution data through the R-tree index. Remote files are read in a few ranged requests
(header, index nodes, then the overlapping data blocks merged into contiguous ranges), so a 500 kb
window of an ENCODE or FANTOM5 track costs a handful of requests instead of downloading the file.

Format: https://genome.ucsc.edu/goldenPath/help/bigWig.html (Kent et al. 2010, Bioinformatics).
"""
from __future__ import annotations

import struct
import threading
import urllib.error
import urllib.request
import zlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

BIGWIG_MAGIC = 0x888FFC26
BPT_MAGIC = 0x78CA8C91
RTREE_MAGIC = 0x2468ACE0
USER_AGENT = "genomics-visualizer/1.0 (+https://github.com/; bigWig range reader)"
HEADER_PREFETCH = 1 << 16
MERGE_GAP = 1 << 16  # data blocks closer than this are fetched in one request
TIMEOUT = 60.0

_BEDGRAPH = np.dtype([("start", "<u4"), ("end", "<u4"), ("value", "<f4")])
_VARSTEP = np.dtype([("start", "<u4"), ("value", "<f4")])


class BigWigError(RuntimeError):
    pass


def is_url(source: str) -> bool:
    return str(source).startswith(("http://", "https://"))


class _Source:
    """Random access to a local file or a URL (follows redirects once, re-resolves on 403/expiry)."""

    def __init__(self, source: Union[str, Path], timeout: float = TIMEOUT):
        self.source = str(source)
        self.remote = is_url(self.source)
        self.timeout = timeout
        self.resolved: Optional[str] = None
        self.requests = 0

    def read(self, offset: int, size: int) -> bytes:
        if size <= 0:
            return b""
        if not self.remote:
            with open(self.source, "rb") as handle:
                handle.seek(offset)
                return handle.read(size)
        for attempt in range(3):
            url = self.resolved or self.source
            request = urllib.request.Request(url, headers={"Range": f"bytes={offset}-{offset + size - 1}", "User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    self.requests += 1
                    status = getattr(response, "status", 200)
                    data = response.read()
                    if self.resolved is None:
                        self.resolved = response.geturl()
                    if status == 200 and len(data) > size:  # server ignored Range
                        data = data[offset:offset + size]
                    return data
            except urllib.error.HTTPError as exc:
                # Pre-signed redirect targets (ENCODE -> S3) expire: resolve again from the source.
                if exc.code in (401, 403) and self.resolved and attempt < 2:
                    self.resolved = None
                    continue
                if exc.code == 416:
                    return b""
                raise BigWigError(f"{self.source}: HTTP {exc.code} {exc.reason}") from exc
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                if attempt < 2:
                    continue
                raise BigWigError(f"{self.source}: {exc}") from exc
        raise BigWigError(f"{self.source}: could not read bytes {offset}-{offset + size - 1}")


class BigWig:
    """A bigWig file; thread-safe (index nodes are cached, data is read on demand)."""

    def __init__(self, source: Union[str, Path], timeout: float = TIMEOUT):
        self._src = _Source(source, timeout)
        self._lock = threading.Lock()
        self._nodes: Dict[int, Tuple[bool, list]] = {}
        head = self._src.read(0, HEADER_PREFETCH)
        if len(head) < 64:
            raise BigWigError(f"{source}: not a bigWig file (too short)")
        magic = struct.unpack_from("<I", head, 0)[0]
        if magic == BIGWIG_MAGIC:
            self.e = "<"
        elif struct.unpack_from(">I", head, 0)[0] == BIGWIG_MAGIC:
            self.e = ">"
        else:
            raise BigWigError(f"{source}: not a bigWig file (magic {magic:#x})")
        e = self.e
        (_, self.version, self.zoom_levels, self.chrom_tree_offset, self.data_offset, self.index_offset,
         _, _, _, _, self.uncompress_buf_size, _) = struct.unpack_from(e + "IHHQQQHHQQIQ", head, 0)
        self._head = head
        self.chroms: Dict[str, Tuple[int, int]] = self._read_chrom_tree()
        self._rtree = self._read_rtree_header()

    @property
    def source(self) -> str:
        return self._src.source

    @property
    def requests(self) -> int:
        return self._src.requests

    def _bytes(self, offset: int, size: int) -> bytes:
        if offset + size <= len(self._head):
            return self._head[offset:offset + size]
        return self._src.read(offset, size)

    # -- chromosome B+ tree -------------------------------------------------------------
    def _read_chrom_tree(self) -> Dict[str, Tuple[int, int]]:
        e = self.e
        header = self._bytes(self.chrom_tree_offset, 32)
        magic, block_size, key_size, val_size, item_count, _ = struct.unpack_from(e + "IIIIQQ", header, 0)
        if magic != BPT_MAGIC:
            raise BigWigError(f"{self.source}: bad chromosome tree")
        chroms: Dict[str, Tuple[int, int]] = {}

        def node(offset: int) -> None:
            head = self._bytes(offset, 4)
            is_leaf, _, count = struct.unpack_from(e + "BBH", head, 0)
            item = key_size + (8 if is_leaf else 8)
            body = self._bytes(offset + 4, count * item)
            for i in range(count):
                base = i * item
                key = body[base:base + key_size].split(b"\0", 1)[0].decode("ascii", "replace")
                if is_leaf:
                    chrom_id, chrom_size = struct.unpack_from(e + "II", body, base + key_size)
                    chroms[key] = (chrom_id, chrom_size)
                else:
                    node(struct.unpack_from(e + "Q", body, base + key_size)[0])

        node(self.chrom_tree_offset + 32)
        return chroms

    # -- R-tree index -------------------------------------------------------------------
    def _read_rtree_header(self) -> Dict[str, int]:
        head = self._bytes(self.index_offset, 48)
        magic, block_size, item_count = struct.unpack_from(self.e + "IIQ", head, 0)
        if magic != RTREE_MAGIC:
            raise BigWigError(f"{self.source}: bad R-tree index")
        return {"block_size": block_size, "root": self.index_offset + 48}

    def _node(self, offset: int) -> Tuple[bool, list]:
        with self._lock:
            hit = self._nodes.get(offset)
        if hit is not None:
            return hit
        e = self.e
        head = self._bytes(offset, 4)
        is_leaf, _, count = struct.unpack_from(e + "BBH", head, 0)
        item = 32 if is_leaf else 24
        body = self._bytes(offset + 4, count * item)
        items = []
        for i in range(count):
            if is_leaf:
                items.append(struct.unpack_from(e + "IIIIQQ", body, i * item))
            else:
                items.append(struct.unpack_from(e + "IIIIQ", body, i * item))
        result = (bool(is_leaf), items)
        with self._lock:
            self._nodes[offset] = result
        return result

    def _blocks(self, chrom_id: int, start: int, end: int) -> List[Tuple[int, int]]:
        """(offset, size) of every data block overlapping [start, end) on chrom_id."""
        out: List[Tuple[int, int]] = []

        def overlaps(sc: int, sb: int, ec: int, eb: int) -> bool:
            return (sc, sb) < (chrom_id, end) and (ec, eb) > (chrom_id, start)

        def walk(offset: int) -> None:
            is_leaf, items = self._node(offset)
            for item in items:
                if not overlaps(*item[:4]):
                    continue
                if is_leaf:
                    out.append((item[4], item[5]))
                else:
                    walk(item[4])

        walk(self._rtree["root"])
        return sorted(set(out))

    # -- data ---------------------------------------------------------------------------
    def values(self, chrom: str, start: int, end: int, missing: float = 0.0) -> np.ndarray:
        """Per-base float32 values over 0-based half-open [start, end); bases without data get ``missing``."""
        start, end = max(0, int(start)), int(end)
        out = np.full(max(0, end - start), missing, dtype=np.float32)
        name = self.resolve_chrom(chrom)
        if name is None or end <= start:
            return out
        chrom_id, size = self.chroms[name]
        end = min(end, size)
        blocks = self._blocks(chrom_id, start, end)
        for offset, raw in self._fetch_blocks(blocks):
            self._decode_into(out, start, end, chrom_id, raw)
        return out

    def resolve_chrom(self, chrom: str) -> Optional[str]:
        """``chr15`` / ``15`` naming differences are bridged."""
        if chrom in self.chroms:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        if alt in self.chroms:
            return alt
        if chrom in ("chrM", "MT", "M"):
            for candidate in ("chrM", "MT", "M"):
                if candidate in self.chroms:
                    return candidate
        return None

    def _fetch_blocks(self, blocks: List[Tuple[int, int]]):
        groups: List[List[Tuple[int, int]]] = []
        for block in blocks:
            if groups and block[0] - (groups[-1][-1][0] + groups[-1][-1][1]) <= MERGE_GAP:
                groups[-1].append(block)
            else:
                groups.append([block])
        for group in groups:
            first = group[0][0]
            last = group[-1][0] + group[-1][1]
            data = self._src.read(first, last - first)
            for offset, size in group:
                chunk = data[offset - first:offset - first + size]
                yield offset, (zlib.decompress(chunk) if self.uncompress_buf_size else chunk)

    def _decode_into(self, out: np.ndarray, start: int, end: int, chrom_id: int, raw: bytes) -> None:
        e = self.e
        sec_chrom, sec_start, sec_end, step, span, kind, _, count = struct.unpack_from(e + "IIIIIBBH", raw, 0)
        if sec_chrom != chrom_id or sec_end <= start or sec_start >= end:
            return
        body = raw[24:]
        if kind == 1:  # bedGraph
            rec = np.frombuffer(body, dtype=_BEDGRAPH.newbyteorder(e) if e == ">" else _BEDGRAPH, count=count)
            starts, ends, vals = rec["start"].astype(np.int64), rec["end"].astype(np.int64), rec["value"]
        elif kind == 2:  # variableStep
            rec = np.frombuffer(body, dtype=_VARSTEP.newbyteorder(e) if e == ">" else _VARSTEP, count=count)
            starts = rec["start"].astype(np.int64)
            ends, vals = starts + span, rec["value"]
        elif kind == 3:  # fixedStep
            vals = np.frombuffer(body, dtype=np.dtype(e + "f4"), count=count)
            starts = sec_start + step * np.arange(count, dtype=np.int64)
            ends = starts + span
        else:
            raise BigWigError(f"{self.source}: unknown section type {kind}")
        a = np.clip(starts, start, end) - start
        b = np.clip(ends, start, end) - start
        keep = b > a
        a, b, vals = a[keep], b[keep], vals[keep]
        if not a.size:
            return
        if np.all(b - a == 1):
            out[a] = vals
            return
        lengths = b - a
        idx = np.repeat(a - np.cumsum(np.concatenate(([0], lengths[:-1]))), lengths) + np.arange(lengths.sum())
        out[idx] = np.repeat(vals, lengths)


_open_lock = threading.Lock()
_open: Dict[str, BigWig] = {}


def open_bigwig(source: Union[str, Path], timeout: float = TIMEOUT) -> BigWig:
    """Opened files are kept (their header and index nodes are reused across reads)."""
    key = str(source)
    with _open_lock:
        hit = _open.get(key)
    if hit is not None:
        return hit
    bw = BigWig(source, timeout)
    with _open_lock:
        if len(_open) > 256:
            _open.clear()
        _open[key] = bw
    return bw
