#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cache em disco de StandardScaler + IncrementalPCA para baselines sklearn.

O cache é compartilhado entre experimentos com o mesmo dataset processado
(mesmo diretório em datasets/) e o mesmo pca_components pedido.

CLI (após o cache do dataset existir) — este mesmo ficheiro:
    python3 sklearn_pca_cache.py --config configs/genes_1000.yaml
    python3 sklearn_pca_cache.py --config configs/genes_1000.yaml --force
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import joblib
import numpy as np
import torch
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader

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
    pca: IncrementalPCA,
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


def get_sklearn_pca_cache_dir_path(dataset_cache_dir: Path, pca_components_requested: int) -> Path:
    """
    dataset_cache_dir = .../processed_cache/datasets/<dataset_name>
    Retorno: .../processed_cache/pca_cache/<dataset_name>_pca<req>
    """
    dataset_cache_dir = Path(dataset_cache_dir)
    base = dataset_cache_dir.parent.parent
    tag = dataset_cache_dir.name
    return base / SKLEARN_PCA_CACHE_SUBDIR / f"{tag}_pca{int(pca_components_requested)}"


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
    align = bool(sk.get("pca_align_n_train", False))
    k = _effective_pca_k(pca_req, n_train, n_features, log)
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

    n_train = len(train_loader.dataset)
    n_val = len(val_loader.dataset)
    n_test = len(test_loader.dataset)
    first = next(iter(train_loader))
    n_features = int(np.prod(first[0].shape[1:]))

    effective_k, pca_req = compute_sklearn_pca_effective_k(
        config, n_train=n_train, n_features=n_features, log=log
    )

    out_dir = get_sklearn_pca_cache_dir_path(dataset_cache_dir, pca_req)
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
        )
    ):
        log(f"[green]✓ Cache PCA já existe e é válido: {out_dir}[/green]")
        return out_dir

    if force and out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log(
        f"[cyan]Construindo cache PCA → {out_dir} | k={effective_k} (pedido={pca_req}, n_train={n_train}, D={n_features})[/cyan]"
    )

    if rich_console is None:
        log("[cyan]PCA cache: StandardScaler (treino)...[/cyan]")
    scaler = fit_standard_scaler_incremental(
        train_loader,
        rich_console=rich_console,
        progress_desc="PCA cache: StandardScaler (treino)",
    )
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


def _load_neural_ancestry_predictor_module():
    """Carrega neural_ancestry_predictor.py pelo caminho (evita falha de import pelo cwd)."""
    import importlib.util

    here = Path(__file__).resolve().parent
    nap_path = here / "neural_ancestry_predictor.py"
    if not nap_path.is_file():
        raise FileNotFoundError(f"neural_ancestry_predictor.py não encontrado em {here}")
    spec = importlib.util.spec_from_file_location("neural_ancestry_predictor", nap_path)
    if spec is None or spec.loader is None:
        raise ImportError("Não foi possível carregar neural_ancestry_predictor")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["neural_ancestry_predictor"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_neural_ancestry_predictor_for_cli():
    """Carrega ``neural_ancestry_predictor`` pelo caminho do pacote (scripts CLI no mesmo diretório)."""
    return _load_neural_ancestry_predictor_module()


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

    nap = _load_neural_ancestry_predictor_module()
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
            "[yellow]Execute primeiro o neural_ancestry_predictor em modo train para gerar o cache do dataset.[/yellow]"
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
