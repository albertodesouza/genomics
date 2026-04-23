from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Tuple

import joblib
import numpy as np
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler

from .config import PredictionConfig
from .data import create_dataloaders, generate_dataset_name


METADATA_FILENAME = "pca_metadata.json"
SCALER_PCA_FILENAME = "scaler_pca.joblib"
COMPLETE_FLAG = ".pca_cache_complete"


def sklearn_flatten_batch(features) -> np.ndarray:
    array = features.detach().cpu().numpy()
    return np.reshape(array, (array.shape[0], -1)).astype(np.float64, copy=False)


def get_pca_cache_dir(config: PredictionConfig) -> Path:
    cache_root = Path(config.dataset_input.processed_cache_dir or ".") / "pca_cache"
    requested = int(config.model.sklearn.get("pca_components", 100))
    return cache_root / f"{generate_dataset_name(config)}_pca{requested}"


def pca_cache_is_valid(cache_dir: Path, requested_components: int) -> bool:
    if not (cache_dir / COMPLETE_FLAG).exists():
        return False
    if not (cache_dir / SCALER_PCA_FILENAME).exists():
        return False
    if not (cache_dir / METADATA_FILENAME).exists():
        return False
    metadata = json.loads((cache_dir / METADATA_FILENAME).read_text(encoding="utf-8"))
    return int(metadata.get("pca_components_requested", -1)) == int(requested_components)


def ensure_sklearn_pca_cache(config: PredictionConfig, force: bool = False) -> Path:
    requested_components = int(config.model.sklearn.get("pca_components", 100))
    cache_dir = get_pca_cache_dir(config)
    if cache_dir.exists() and pca_cache_is_valid(cache_dir, requested_components) and not force:
        return cache_dir

    if cache_dir.exists():
        for child in cache_dir.iterdir():
            if child.is_file():
                child.unlink()

    cache_dir.mkdir(parents=True, exist_ok=True)
    data_bundle = create_dataloaders(config)
    train_loader = data_bundle.train_loader
    val_loader = data_bundle.val_loader
    test_loader = data_bundle.test_loader

    scaler = StandardScaler()
    train_batches = []
    train_targets = []
    for features, targets in train_loader:
        batch = sklearn_flatten_batch(features)
        scaler.partial_fit(batch)
        train_batches.append(batch)
        train_targets.append(targets.detach().cpu().numpy())

    X_train = np.vstack([scaler.transform(batch) for batch in train_batches])
    y_train = np.concatenate(train_targets)
    effective_k = max(1, min(requested_components, X_train.shape[0], X_train.shape[1]))
    pca = IncrementalPCA(n_components=effective_k)
    pca.fit(X_train)

    X_train_pca = pca.transform(X_train)
    X_val_pca, y_val = _transform_loader(val_loader, scaler, pca)
    X_test_pca, y_test = _transform_loader(test_loader, scaler, pca)

    np.save(cache_dir / "X_train.npy", X_train_pca)
    np.save(cache_dir / "y_train.npy", y_train)
    np.save(cache_dir / "X_val.npy", X_val_pca)
    np.save(cache_dir / "y_val.npy", y_val)
    np.save(cache_dir / "X_test.npy", X_test_pca)
    np.save(cache_dir / "y_test.npy", y_test)
    joblib.dump({"scaler": scaler, "pca": pca}, cache_dir / SCALER_PCA_FILENAME)
    metadata = {
        "pca_components_requested": requested_components,
        "pca_components_effective": effective_k,
        "n_train": int(X_train_pca.shape[0]),
        "n_val": int(X_val_pca.shape[0]),
        "n_test": int(X_test_pca.shape[0]),
        "n_features_original": int(X_train.shape[1]),
    }
    (cache_dir / METADATA_FILENAME).write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (cache_dir / COMPLETE_FLAG).write_text("ok\n", encoding="utf-8")
    return cache_dir


def _transform_loader(loader, scaler: StandardScaler, pca: IncrementalPCA) -> Tuple[np.ndarray, np.ndarray]:
    xs = []
    ys = []
    for features, targets in loader:
        batch = scaler.transform(sklearn_flatten_batch(features))
        xs.append(pca.transform(batch))
        ys.append(targets.detach().cpu().numpy())
    if not xs:
        return np.empty((0, pca.n_components_)), np.empty((0,))
    return np.vstack(xs), np.concatenate(ys)
