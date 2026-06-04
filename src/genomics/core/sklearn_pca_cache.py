#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cache em disco de StandardScaler + IncrementalPCA para baselines sklearn.

O cache é compartilhado entre experimentos com o mesmo dataset processado
(mesmo diretório em datasets/) e o mesmo pca_components pedido.

CLI (após o cache do dataset existir):
    genomics neural pca-cache legacy/neural_ancestry_predictor_deprecated/configs/genes_1000.yaml
    genomics neural pca-cache legacy/neural_ancestry_predictor_deprecated/configs/genes_1000.yaml --force
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import torch
from torch.utils.data import DataLoader

if TYPE_CHECKING:
    from sklearn.decomposition import IncrementalPCA
    from sklearn.preprocessing import StandardScaler

try:
    from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn

    _RICH_PROGRESS_AVAILABLE = True
except ImportError:
    _RICH_PROGRESS_AVAILABLE = False

SKLEARN_PCA_CACHE_VERSION = 1
SKLEARN_PCA_CACHE_SUBDIR = "pca_cache"
SCALER_PCA_FILENAME = "scaler_pca.joblib"
METADATA_FILENAME = "pca_metadata.json"
COMPLETE_FLAG = ".pca_cache_complete"

LAPACK_INT32_INDEX_MAX = 2147483647


class StreamingRandomizedPCA:
    """PCA randomizada ajustada em streaming, com componentes em memmap."""

    def __init__(
        self,
        *,
        components_path: str,
        n_components: int,
        n_features: int,
        dtype: str,
        feature_chunk_size: int,
        explained_variance: np.ndarray,
        explained_variance_ratio: np.ndarray,
        singular_values: np.ndarray,
    ) -> None:
        self.components_path = str(components_path)
        self.n_components = int(n_components)
        self.n_components_ = int(n_components)
        self.n_features_in_ = int(n_features)
        self.dtype = str(dtype)
        self.feature_chunk_size = int(feature_chunk_size)
        self.explained_variance_ = np.asarray(explained_variance, dtype=np.float64)
        self.explained_variance_ratio_ = np.asarray(explained_variance_ratio, dtype=np.float64)
        self.singular_values_ = np.asarray(singular_values, dtype=np.float64)

    def _components_memmap(self, mode: str = "r") -> np.memmap:
        return np.memmap(
            self.components_path,
            mode=mode,
            dtype=np.dtype(self.dtype),
            shape=(self.n_components_, self.n_features_in_),
        )

    @property
    def components_(self) -> np.memmap:
        return self._components_memmap("r")

    def transform(self, X: np.ndarray) -> np.ndarray:
        if X.shape[1] != self.n_features_in_:
            raise ValueError(f"X tem {X.shape[1]} features; esperado {self.n_features_in_}")
        dtype = np.dtype(self.dtype)
        components = self._components_memmap("r")
        out = np.zeros((X.shape[0], self.n_components_), dtype=np.float32)
        for start in range(0, self.n_features_in_, self.feature_chunk_size):
            end = min(start + self.feature_chunk_size, self.n_features_in_)
            out += X[:, start:end].astype(dtype, copy=False) @ components[:, start:end].T
        return out


def max_incremental_pca_components_lapack_safe(n_features: int) -> int:
    nf = max(int(n_features), 1)
    max_rows = LAPACK_INT32_INDEX_MAX // nf
    return max(1, (max_rows - 1) // 2)


def largest_k_dividing_n_train(n_train: int, k_cap: int) -> int:
    """
    Maior k em [1, min(k_cap, n_train)] tal que n_train % k == 0.
    Garante que, em IncrementalPCA por blocos de k linhas, não sobra resto < k
    (evita padding com linhas repetidas no último partial_fit).
    """
    n_train = int(n_train)
    k_cap = int(k_cap)
    if n_train <= 0:
        return 1
    hi = min(k_cap, n_train)
    if hi < 1:
        return 1
    for k in range(hi, 0, -1):
        if n_train % k == 0:
            return k
    return 1


def sklearn_flatten_batch(features: torch.Tensor) -> np.ndarray:
    x = features.detach().cpu().numpy()
    return np.reshape(x, (x.shape[0], -1)).astype(np.float64, copy=False)


class _NoRichProgress:
    """Context manager stand-in when Rich is unavailable or console is None."""

    def __enter__(self) -> "_NoRichProgress":
        return self

    def __exit__(self, *args: Any) -> bool:
        return False

    def add_task(self, description: str, total: Optional[int] = None) -> int:
        return 0

    def update(self, task_id: int, advance: int = 1, **kwargs: Any) -> None:
        pass


def _rich_progress_cm(rich_console: Optional[Any]) -> Any:
    if rich_console is not None and _RICH_PROGRESS_AVAILABLE:
        return Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=rich_console,
        )
    return _NoRichProgress()


def _dataloader_len(loader: DataLoader) -> Optional[int]:
    try:
        return len(loader)
    except TypeError:
        return None


def fit_standard_scaler_incremental(
    train_loader: DataLoader,
    *,
    rich_console: Optional[Any] = None,
    progress_desc: str = "PCA: StandardScaler (treino)",
) -> StandardScaler:
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    total = _dataloader_len(train_loader)
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task(progress_desc, total=total)
        for features, _targets, _idx in train_loader:
            scaler.partial_fit(sklearn_flatten_batch(features))
            progress.update(tid, advance=1)
    return scaler


def fit_incremental_pca_on_train(
    train_loader: DataLoader,
    scaler: StandardScaler,
    n_components: int,
    log: Optional[Callable[[str], None]] = None,
    forbid_tail_padding: bool = False,
    *,
    rich_console: Optional[Any] = None,
    progress_desc: str = "PCA: IncrementalPCA (treino)",
) -> IncrementalPCA:
    from sklearn.decomposition import IncrementalPCA

    log = log or print
    pca = IncrementalPCA(n_components=n_components)
    buf_list: List[np.ndarray] = []

    def _nbuf() -> int:
        return sum(x.shape[0] for x in buf_list)

    total = _dataloader_len(train_loader)
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task(progress_desc, total=total)
        for features, _targets, _idx in train_loader:
            X = scaler.transform(sklearn_flatten_batch(features))
            buf_list.append(X)
            while _nbuf() >= n_components:
                big = np.vstack(buf_list) if len(buf_list) > 1 else buf_list[0]
                pca.partial_fit(big[:n_components])
                rest = big[n_components:]
                buf_list = [rest] if rest.shape[0] > 0 else []
            progress.update(tid, advance=1)

        if _nbuf() > 0:
            big = np.vstack(buf_list) if len(buf_list) > 1 else buf_list[0]
            if big.shape[0] >= n_components:
                pca.partial_fit(big)
            else:
                r = big.shape[0]
                if forbid_tail_padding:
                    raise RuntimeError(
                        f"IncrementalPCA: resto final r={r} < k={n_components} com forbid_tail_padding=True. "
                        f"Ative model.sklearn.pca_align_n_train ou verifique n_train % k == 0."
                    )
                pad = big[np.random.RandomState(0).choice(r, n_components - r, replace=True)]
                log(
                    "[dim]IncrementalPCA: último partial_fit precisa de ≥k linhas (regra do sklearn). "
                    f"Sobraram {r} amostras após blocos de k={n_components}; "
                    f"completamos com {n_components - r} linhas repetidas (amostragem com reposição, seed fixa). "
                    "Isto é esperado quando n_train não é múltiplo de k; todas as amostras reais entram nos blocos "
                    "anteriores e neste resto; na etapa seguinte (transform) cada amostra é projetada uma vez. "
                    "Para evitar repetições, use pca_align_n_train: true no YAML.[/dim]"
                )
                pca.partial_fit(np.vstack([big, pad]))

    return pca


def stack_scaled_pca_batches(
    loader: DataLoader,
    scaler: StandardScaler,
    pca: Any,
    *,
    rich_console: Optional[Any] = None,
    progress_desc: str = "PCA: transform",
) -> Tuple[np.ndarray, np.ndarray]:
    xs: List[np.ndarray] = []
    ys: List[np.ndarray] = []
    total = _dataloader_len(loader)
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task(progress_desc, total=total)
        for features, targets, _idx in loader:
            X = scaler.transform(sklearn_flatten_batch(features))
            Xr = pca.transform(X)
            xs.append(Xr)
            t = targets.detach().cpu().numpy()
            ys.append(t.reshape(-1))
            progress.update(tid, advance=1)
    if not xs:
        k = getattr(pca, "n_components_", pca.n_components)
        return np.empty((0, int(k))), np.empty((0,), dtype=np.int64)
    return np.vstack(xs), np.concatenate(ys)


def _random_omega_chunk(
    start: int,
    end: int,
    ell: int,
    *,
    dtype: np.dtype,
    random_state: int,
) -> np.ndarray:
    # Seed por chunk: reproduzível sem materializar a matriz Omega D x ell.
    rng = np.random.default_rng(np.random.SeedSequence([int(random_state), int(start)]))
    return rng.standard_normal((end - start, ell)).astype(dtype, copy=False)


def _batch_indices_array(indices: Any) -> np.ndarray:
    if hasattr(indices, "detach"):
        return indices.detach().cpu().numpy().reshape(-1)
    return np.asarray(indices).reshape(-1)


def fit_streaming_randomized_pca_on_train(
    train_loader: DataLoader,
    scaler: StandardScaler,
    n_components: int,
    output_dir: Path,
    *,
    oversampling: int = 32,
    n_iter: int = 2,
    feature_chunk_size: int = 16384,
    dtype: str = "float32",
    random_state: int = 42,
    log: Optional[Callable[[str], None]] = None,
    rich_console: Optional[Any] = None,
) -> StreamingRandomizedPCA:
    """
    Ajusta PCA randomizada sem chamar SVD/LAPACK sobre matrizes n x D grandes.

    A única decomposição densa é em l x l, onde l = k + oversampling. Os vetores
    de componentes são gravados como memmap para não inflar o scaler_pca.joblib.
    """
    log = log or print
    output_dir = Path(output_dir)
    np_dtype = np.dtype(dtype)
    if np_dtype not in (np.dtype("float32"), np.dtype("float64")):
        raise ValueError("dtype deve ser float32 ou float64")

    first = next(iter(train_loader))
    n_features = int(np.prod(first[0].shape[1:]))
    n_train = len(train_loader.dataset)
    k = int(n_components)
    ell = min(n_features, n_train, k + int(oversampling))
    if ell < k:
        raise ValueError(f"randomized PCA: ell={ell} menor que k={k}")

    log(
        f"[cyan]PCA randomizada streaming: k={k}, oversampling={oversampling}, "
        f"n_iter={n_iter}, ell={ell}, n_train={n_train}, D={n_features}, dtype={np_dtype.name}[/cyan]"
    )

    Y = np.zeros((n_train, ell), dtype=np_dtype)
    index_to_row: Dict[int, int] = {}
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task("PCA randomizada: X @ Omega", total=_dataloader_len(train_loader))
        row_start = 0
        for features, _targets, idx in train_loader:
            X = scaler.transform(sklearn_flatten_batch(features)).astype(np_dtype, copy=False)
            row_end = row_start + X.shape[0]
            batch_indices = _batch_indices_array(idx)
            if len(batch_indices) != X.shape[0]:
                raise RuntimeError("Batch PCA randomizada tem número de índices diferente do número de amostras")
            for offset, sample_idx in enumerate(batch_indices):
                sample_idx = int(sample_idx)
                if sample_idx in index_to_row:
                    raise RuntimeError(f"Índice de treino duplicado no PCA randomizado: {sample_idx}")
                index_to_row[sample_idx] = row_start + offset
            Y_batch = np.zeros((X.shape[0], ell), dtype=np_dtype)
            for start in range(0, n_features, feature_chunk_size):
                end = min(start + feature_chunk_size, n_features)
                omega = _random_omega_chunk(start, end, ell, dtype=np_dtype, random_state=random_state)
                Y_batch += X[:, start:end] @ omega
            Y[row_start:row_end] = Y_batch
            row_start = row_end
            progress.update(tid, advance=1)

    Q, _ = np.linalg.qr(Y.astype(np.float64, copy=False), mode="reduced")
    Q = Q[:, :ell]
    del Y
    if len(index_to_row) != n_train:
        raise RuntimeError(f"PCA randomizada viu {len(index_to_row)} índices únicos; esperado {n_train}")

    for iteration in range(int(n_iter)):
        Z_path = output_dir / f"randomized_pca_power_Z_iter{iteration}.{np_dtype.name}.memmap"
        Z = np.memmap(Z_path, mode="w+", dtype=np_dtype, shape=(ell, n_features))
        Z[:] = 0
        with _rich_progress_cm(rich_console) as progress:
            tid = progress.add_task(
                f"PCA randomizada: power {iteration + 1}/{n_iter} X.T @ Q",
                total=_dataloader_len(train_loader),
            )
            for features, _targets, idx in train_loader:
                X = scaler.transform(sklearn_flatten_batch(features)).astype(np_dtype, copy=False)
                batch_rows = [index_to_row[int(i)] for i in _batch_indices_array(idx)]
                Q_batch = Q[batch_rows]
                for start in range(0, n_features, feature_chunk_size):
                    end = min(start + feature_chunk_size, n_features)
                    Z[:, start:end] += (Q_batch.T @ X[:, start:end]).astype(np_dtype, copy=False)
                progress.update(tid, advance=1)

        Y_power = np.zeros((n_train, ell), dtype=np_dtype)
        with _rich_progress_cm(rich_console) as progress:
            tid = progress.add_task(
                f"PCA randomizada: power {iteration + 1}/{n_iter} X @ Z.T",
                total=_dataloader_len(train_loader),
            )
            for features, _targets, idx in train_loader:
                X = scaler.transform(sklearn_flatten_batch(features)).astype(np_dtype, copy=False)
                batch_rows = [index_to_row[int(i)] for i in _batch_indices_array(idx)]
                Y_batch = np.zeros((X.shape[0], ell), dtype=np_dtype)
                for start in range(0, n_features, feature_chunk_size):
                    end = min(start + feature_chunk_size, n_features)
                    Y_batch += X[:, start:end] @ np.asarray(Z[:, start:end]).T
                Y_power[batch_rows] = Y_batch
                progress.update(tid, advance=1)

        try:
            Z._mmap.close()
        except Exception:
            pass
        try:
            Z_path.unlink()
        except OSError:
            pass

        Q, _ = np.linalg.qr(Y_power.astype(np.float64, copy=False), mode="reduced")
        Q = Q[:, :ell]
        del Y_power

    B_path = output_dir / f"randomized_pca_B.{np_dtype.name}.memmap"
    B = np.memmap(B_path, mode="w+", dtype=np_dtype, shape=(ell, n_features))
    B[:] = 0
    total_sum_sq = 0.0
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task("PCA randomizada: Q.T @ X", total=_dataloader_len(train_loader))
        for features, _targets, idx in train_loader:
            X = scaler.transform(sklearn_flatten_batch(features)).astype(np_dtype, copy=False)
            batch_rows = [index_to_row[int(i)] for i in _batch_indices_array(idx)]
            Q_batch = Q[batch_rows]
            total_sum_sq += float(np.square(X.astype(np.float64, copy=False)).sum())
            for start in range(0, n_features, feature_chunk_size):
                end = min(start + feature_chunk_size, n_features)
                B[:, start:end] += (Q_batch.T @ X[:, start:end]).astype(np_dtype, copy=False)
            progress.update(tid, advance=1)
    del Q

    C = np.zeros((ell, ell), dtype=np.float64)
    total_chunks = (n_features + feature_chunk_size - 1) // feature_chunk_size
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task("PCA randomizada: B @ B.T", total=total_chunks)
        for start in range(0, n_features, feature_chunk_size):
            end = min(start + feature_chunk_size, n_features)
            block = np.asarray(B[:, start:end], dtype=np.float64)
            C += block @ block.T
            progress.update(tid, advance=1)

    evals, U = np.linalg.eigh(C)
    order = np.argsort(evals)[::-1]
    evals = np.maximum(evals[order], 0.0)
    U = U[:, order]
    U_k = U[:, :k]

    components_path = output_dir / f"randomized_pca_components.{np_dtype.name}.memmap"
    components = np.memmap(components_path, mode="w+", dtype=np_dtype, shape=(k, n_features))
    singular_values = np.sqrt(evals[:k])
    inv_singular_values = np.zeros_like(singular_values)
    nonzero = singular_values > np.finfo(np.float64).eps
    inv_singular_values[nonzero] = 1.0 / singular_values[nonzero]
    with _rich_progress_cm(rich_console) as progress:
        tid = progress.add_task("PCA randomizada: componentes", total=total_chunks)
        for start in range(0, n_features, feature_chunk_size):
            end = min(start + feature_chunk_size, n_features)
            block_components = U_k.T @ np.asarray(B[:, start:end], dtype=np.float64)
            block_components *= inv_singular_values[:, None]
            components[:, start:end] = block_components.astype(np_dtype)
            progress.update(tid, advance=1)
    components.flush()

    denom = max(n_train - 1, 1)
    explained_variance = (singular_values ** 2) / denom
    total_variance = total_sum_sq / denom if n_train > 1 else 0.0
    if total_variance > 0:
        explained_variance_ratio = explained_variance / total_variance
    else:
        explained_variance_ratio = np.zeros_like(explained_variance)

    try:
        B._mmap.close()
    except Exception:
        pass
    try:
        B_path.unlink()
    except OSError:
        pass
    return StreamingRandomizedPCA(
        components_path=str(components_path.resolve()),
        n_components=k,
        n_features=n_features,
        dtype=np_dtype.name,
        feature_chunk_size=feature_chunk_size,
        explained_variance=explained_variance,
        explained_variance_ratio=explained_variance_ratio,
        singular_values=singular_values,
    )


def get_sklearn_pca_cache_dir_path(
    dataset_cache_dir: Path,
    pca_components_requested: int,
    pca_backend: str = "incremental",
) -> Path:
    """
    dataset_cache_dir = .../processed_cache/datasets/<dataset_name>
    Retorno: .../processed_cache/pca_cache/<dataset_name>_pca<req>
    """
    dataset_cache_dir = Path(dataset_cache_dir)
    base = dataset_cache_dir.parent.parent
    tag = dataset_cache_dir.name
    suffix = f"_pca{int(pca_components_requested)}"
    if pca_backend != "incremental":
        suffix += f"_{pca_backend}"
    return base / SKLEARN_PCA_CACHE_SUBDIR / f"{tag}{suffix}"


def _effective_pca_k(
    pca_req: int,
    n_train: int,
    n_features: int,
    log: Callable[[str], None],
) -> int:
    effective_k = max(1, min(pca_req, n_train, n_features))
    k_lapack = max_incremental_pca_components_lapack_safe(n_features)
    if effective_k > k_lapack:
        log(
            f"[yellow]⚠ IncrementalPCA: n_components {effective_k} → {k_lapack} "
            f"(evita overflow int32 no LAPACK: matrizes ~2k×D com D={n_features} "
            f"excedem {LAPACK_INT32_INDEX_MAX} elementos).[/yellow]"
        )
        effective_k = k_lapack
    return effective_k


def compute_sklearn_pca_effective_k(
    config: Dict[str, Any],
    *,
    n_train: int,
    n_features: int,
    log: Optional[Callable[[str], None]] = None,
) -> Tuple[int, int]:
    """
    Resolve k efetivo e o pca_components pedido no YAML.
    Com ``model.sklearn.pca_align_n_train: true``, reduz k (se preciso) para o maior
    divisor de n_train que não exceda o k já limitado por n_train/n_features/LAPACK.
    """
    log = log or print
    sk = config.get("model", {}).get("sklearn", {})
    pca_req = sk.get("pca_components")
    if pca_req is None:
        pca_req = min(500, n_train, n_features)
    pca_req = int(pca_req)
    backend = str(sk.get("pca_backend", "incremental"))
    align = bool(sk.get("pca_align_n_train", False))
    if backend == "incremental":
        k = _effective_pca_k(pca_req, n_train, n_features, log)
    elif backend == "randomized_streaming":
        k = max(1, min(pca_req, n_train, n_features))
    else:
        raise ValueError(f"pca_backend não suportado: {backend}")
    if align:
        k_aligned = largest_k_dividing_n_train(n_train, k)
        if k_aligned != k:
            log(
                f"[cyan]pca_align_n_train: k {k} → {k_aligned} "
                f"(n_train={n_train} múltiplo de k; evita padding no último partial_fit)[/cyan]"
            )
        k = k_aligned
    return k, pca_req


def _read_metadata(cache_dir: Path) -> Optional[Dict[str, Any]]:
    p = cache_dir / METADATA_FILENAME
    if not p.exists():
        return None
    with open(p, "r") as f:
        return json.load(f)


def pca_cache_is_valid(
    cache_dir: Path,
    dataset_cache_dir: Path,
    pca_components_requested: int,
    n_features: int,
    n_train: int,
    n_val: int,
    n_test: int,
    prediction_target: str,
    pca_align_n_train: bool,
    pca_backend: str,
    randomized_pca_oversampling: int,
    randomized_pca_n_iter: int,
    randomized_pca_feature_chunk_size: int,
    randomized_pca_dtype: str,
) -> bool:
    cache_dir = Path(cache_dir)
    if not (cache_dir / COMPLETE_FLAG).exists():
        return False
    if not (cache_dir / SCALER_PCA_FILENAME).exists():
        return False
    for name in ("X_train.npy", "y_train.npy", "X_val.npy", "y_val.npy", "X_test.npy", "y_test.npy"):
        if not (cache_dir / name).exists():
            return False
    meta = _read_metadata(cache_dir)
    if not meta or meta.get("cache_version") != SKLEARN_PCA_CACHE_VERSION:
        return False
    if meta.get("dataset_cache_stem") != Path(dataset_cache_dir).name:
        return False
    if int(meta.get("pca_components_requested", -1)) != int(pca_components_requested):
        return False
    if int(meta.get("n_features_original", -1)) != int(n_features):
        return False
    if int(meta.get("n_train", -1)) != int(n_train):
        return False
    if int(meta.get("n_val", -1)) != int(n_val):
        return False
    if int(meta.get("n_test", -1)) != int(n_test):
        return False
    if meta.get("prediction_target") != prediction_target:
        return False
    if bool(meta.get("pca_align_n_train", False)) != bool(pca_align_n_train):
        return False
    if str(meta.get("pca_backend", "incremental")) != str(pca_backend):
        return False
    if str(pca_backend) == "randomized_streaming":
        if int(meta.get("randomized_pca_oversampling", -1)) != int(randomized_pca_oversampling):
            return False
        if int(meta.get("randomized_pca_n_iter", -1)) != int(randomized_pca_n_iter):
            return False
        if int(meta.get("randomized_pca_feature_chunk_size", -1)) != int(randomized_pca_feature_chunk_size):
            return False
        if str(meta.get("randomized_pca_dtype", "")) != str(randomized_pca_dtype):
            return False
        components_path = cache_dir / f"randomized_pca_components.{randomized_pca_dtype}.memmap"
        if not components_path.exists():
            return False
    return True


def build_sklearn_pca_cache(
    config: Dict[str, Any],
    dataset_cache_dir: Path,
    train_loader: DataLoader,
    val_loader: DataLoader,
    test_loader: DataLoader,
    *,
    force: bool = False,
    log: Callable[[str], None] = print,
    rich_console: Optional[Any] = None,
) -> Path:
    """
    Ajusta StandardScaler + IncrementalPCA só no treino; transforma train/val/test;
    grava arrays .npy + scaler_pca.joblib + metadata.
    """
    dataset_cache_dir = Path(dataset_cache_dir)
    sk = config.get("model", {}).get("sklearn", {})
    align = bool(sk.get("pca_align_n_train", False))
    pca_backend = str(sk.get("pca_backend", "incremental"))
    randomized_oversampling = int(sk.get("randomized_pca_oversampling", 32))
    randomized_n_iter = int(sk.get("randomized_pca_n_iter", 2))
    randomized_feature_chunk_size = int(sk.get("randomized_pca_feature_chunk_size", 16384))
    randomized_dtype = str(sk.get("randomized_pca_dtype", "float32"))

    n_train = len(train_loader.dataset)
    n_val = len(val_loader.dataset)
    n_test = len(test_loader.dataset)
    first = next(iter(train_loader))
    n_features = int(np.prod(first[0].shape[1:]))

    effective_k, pca_req = compute_sklearn_pca_effective_k(
        config, n_train=n_train, n_features=n_features, log=log
    )

    out_dir = get_sklearn_pca_cache_dir_path(dataset_cache_dir, pca_req, pca_backend)
    pred_target = config.get("output", {}).get("prediction_target", "")

    if (
        not force
        and pca_cache_is_valid(
            out_dir,
            dataset_cache_dir,
            pca_req,
            n_features,
            n_train,
            n_val,
            n_test,
            pred_target,
            align,
            pca_backend,
            randomized_oversampling,
            randomized_n_iter,
            randomized_feature_chunk_size,
            randomized_dtype,
        )
    ):
        log(f"[green]✓ Cache PCA já existe e é válido: {out_dir}[/green]")
        return out_dir

    if force and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log(
        f"[cyan]Construindo cache PCA → {out_dir} | backend={pca_backend}, "
        f"k={effective_k} (pedido={pca_req}, n_train={n_train}, D={n_features})[/cyan]"
    )

    if rich_console is None:
        log("[cyan]PCA cache: StandardScaler (treino)...[/cyan]")
    scaler = fit_standard_scaler_incremental(
        train_loader,
        rich_console=rich_console,
        progress_desc="PCA cache: StandardScaler (treino)",
    )
    if pca_backend == "incremental":
        if rich_console is None:
            log("[cyan]PCA cache: IncrementalPCA (treino)...[/cyan]")
        pca = fit_incremental_pca_on_train(
            train_loader,
            scaler,
            effective_k,
            log=log,
            forbid_tail_padding=align,
            rich_console=rich_console,
            progress_desc="PCA cache: IncrementalPCA (treino)",
        )
    elif pca_backend == "randomized_streaming":
        if rich_console is None:
            log("[cyan]PCA cache: PCA randomizada streaming (treino)...[/cyan]")
        pca = fit_streaming_randomized_pca_on_train(
            train_loader,
            scaler,
            effective_k,
            out_dir,
            oversampling=randomized_oversampling,
            n_iter=randomized_n_iter,
            feature_chunk_size=randomized_feature_chunk_size,
            dtype=randomized_dtype,
            random_state=42,
            log=log,
            rich_console=rich_console,
        )
    else:
        raise ValueError(f"pca_backend não suportado: {pca_backend}")

    if rich_console is None:
        log("[cyan]PCA cache: transformando splits...[/cyan]")
    X_train, y_train = stack_scaled_pca_batches(
        train_loader,
        scaler,
        pca,
        rich_console=rich_console,
        progress_desc="PCA cache: transform treino",
    )
    X_val, y_val = stack_scaled_pca_batches(
        val_loader,
        scaler,
        pca,
        rich_console=rich_console,
        progress_desc="PCA cache: transform val",
    )
    X_test, y_test = stack_scaled_pca_batches(
        test_loader,
        scaler,
        pca,
        rich_console=rich_console,
        progress_desc="PCA cache: transform test",
    )

    np.save(out_dir / "X_train.npy", X_train.astype(np.float32))
    np.save(out_dir / "y_train.npy", y_train.astype(np.int64))
    np.save(out_dir / "X_val.npy", X_val.astype(np.float32))
    np.save(out_dir / "y_val.npy", y_val.astype(np.int64))
    np.save(out_dir / "X_test.npy", X_test.astype(np.float32))
    np.save(out_dir / "y_test.npy", y_test.astype(np.int64))

    import joblib

    joblib.dump({"scaler": scaler, "pca": pca}, out_dir / SCALER_PCA_FILENAME)

    meta = {
        "cache_version": SKLEARN_PCA_CACHE_VERSION,
        "dataset_cache_dir": str(dataset_cache_dir.resolve()),
        "dataset_cache_stem": dataset_cache_dir.name,
        "pca_components_requested": pca_req,
        "pca_n_components_effective": int(effective_k),
        "n_features_original": int(n_features),
        "n_train": int(n_train),
        "n_val": int(n_val),
        "n_test": int(n_test),
        "prediction_target": pred_target,
        "pca_align_n_train": align,
        "pca_backend": pca_backend,
        "randomized_pca_oversampling": randomized_oversampling,
        "randomized_pca_n_iter": randomized_n_iter,
        "randomized_pca_feature_chunk_size": randomized_feature_chunk_size,
        "randomized_pca_dtype": randomized_dtype,
    }
    with open(out_dir / METADATA_FILENAME, "w") as f:
        json.dump(meta, f, indent=2)

    (out_dir / COMPLETE_FLAG).touch()
    log(f"[green]✓ Cache PCA concluído: {out_dir}[/green]")
    return out_dir


def ensure_sklearn_pca_cache(
    config: Dict[str, Any],
    dataset_cache_dir: Path,
    train_loader: DataLoader,
    val_loader: DataLoader,
    test_loader: DataLoader,
    *,
    force: bool = False,
    log: Callable[[str], None] = print,
    rich_console: Optional[Any] = None,
) -> Path:
    """Garante cache em disco; reconstrói se inválido ou force=True."""
    return build_sklearn_pca_cache(
        config,
        dataset_cache_dir,
        train_loader,
        val_loader,
        test_loader,
        force=force,
        log=log,
        rich_console=rich_console,
    )


def _load_neural_ancestry_predictor_deprecated_module():
    """Carrega neural_ancestry_predictor_deprecated.py pelo caminho (evita falha de import pelo cwd)."""
    import importlib.util

    repo_root = Path(__file__).resolve().parents[3]
    nap_path = repo_root / "legacy" / "neural_ancestry_predictor_deprecated" / "neural_ancestry_predictor_deprecated.py"
    if not nap_path.is_file():
        raise FileNotFoundError(f"neural_ancestry_predictor_deprecated.py não encontrado em {nap_path.parent}")
    spec = importlib.util.spec_from_file_location("neural_ancestry_predictor_deprecated", nap_path)
    if spec is None or spec.loader is None:
        raise ImportError("Não foi possível carregar neural_ancestry_predictor_deprecated")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["neural_ancestry_predictor_deprecated"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_neural_ancestry_predictor_deprecated_for_cli():
    """Carrega ``neural_ancestry_predictor_deprecated`` pelo caminho do pacote (scripts CLI no mesmo diretório)."""
    return _load_neural_ancestry_predictor_deprecated_module()


def _load_neural_ancestry_predictor_module():
    """Compatibility alias for the deprecated package rename."""
    return _load_neural_ancestry_predictor_deprecated_module()


def load_neural_ancestry_predictor_for_cli():
    """Compatibility alias for the deprecated package rename."""
    return load_neural_ancestry_predictor_deprecated_for_cli()


def run_cli() -> None:
    parser = argparse.ArgumentParser(
        description="Gera cache em disco de StandardScaler + IncrementalPCA para baselines sklearn."
    )
    parser.add_argument("--config", type=str, required=True, help="YAML de configuração")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reconstrói o cache mesmo se existir e for válido",
    )
    args = parser.parse_args()

    nap = _load_neural_ancestry_predictor_deprecated_module()
    console = nap.console
    get_dataset_cache_dir = nap.get_dataset_cache_dir
    load_config = nap.load_config
    prepare_data = nap.prepare_data
    validate_cache = nap.validate_cache

    config_path = Path(args.config)
    config = load_config(config_path)
    dataset_cache_dir = Path(get_dataset_cache_dir(config))
    if not dataset_cache_dir.exists() or not validate_cache(dataset_cache_dir, config):
        console.print(
            f"[red]Cache do dataset inválido ou ausente: {dataset_cache_dir}[/red]\n"
            "[yellow]Execute primeiro o neural_ancestry_predictor_deprecated em modo train para gerar o cache do dataset.[/yellow]"
        )
        sys.exit(1)

    exp_dir = Path(config["dataset_input"]["processed_cache_dir"]) / "_pca_cache_cli_build"
    exp_dir.mkdir(parents=True, exist_ok=True)
    _full, train_loader, val_loader, test_loader = prepare_data(config, exp_dir)

    ensure_sklearn_pca_cache(
        config,
        dataset_cache_dir,
        train_loader,
        val_loader,
        test_loader,
        force=args.force,
        log=console.print,
        rich_console=console,
    )


if __name__ == "__main__":
    run_cli()
