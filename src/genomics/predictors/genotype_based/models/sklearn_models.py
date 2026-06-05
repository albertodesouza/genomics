# -*- coding: utf-8 -*-
"""
sklearn_models.py
=================

Baselines sklearn: StandardScaler + IncrementalPCA + SVM / Logistic Regression / RF / XGBoost.

Funções principais
------------------
build_sklearn_classifier     : Instancia classificador baseado no config
train_sklearn_baseline       : Treina pipeline completo (PCA + classificador)
sklearn_predict_labels       : Predições em DataLoader
sklearn_metrics_dict         : Métricas de avaliação
run_sklearn_eval_and_save    : Salva resultados JSON + log W&B
run_sklearn_test_mode        : Avaliação em modo test usando artefato joblib
load_sklearn_baseline_artifact: Carrega artefato salvo
"""
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import joblib
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    precision_recall_fscore_support,
)
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from torch.utils.data import DataLoader
from rich.console import Console

from genomics.predictors.genotype_based.config import PipelineConfig
from genomics.core import update_manifest
from genomics.core.metrics import print_classification_metrics, save_classification_plots, save_results_json

console = Console()

SKLEARN_ARTIFACT_FILENAME = "sklearn_baseline.joblib"
SCALER_PCA_FILENAME = "scaler_pca.joblib"


def _effective_seed(config: PipelineConfig) -> int:
    seed = config.data_split.random_seed
    return 42 if seed is None or seed == -1 else int(seed)


def sklearn_flatten_batch(features) -> np.ndarray:
    """Converte batch de tensores para matriz numpy 2D."""
    return features.detach().cpu().numpy().reshape(features.shape[0], -1)


def build_sklearn_classifier(config: PipelineConfig, model_type: str, random_seed: int, n_train: int = None) -> Any:
    """Instancia classificador (SVM, LOGREG, RF ou XGBoost) a partir do config."""
    model_type = model_type.upper()
    sk = config.model.sklearn

    if model_type == "SVM":
        svm_cfg = sk.svm
        cw = svm_cfg.class_weight
        if cw not in (None, "balanced"):
            cw = None
        base = LinearSVC(C=svm_cfg.C, max_iter=svm_cfg.max_iter,
                         dual=False, class_weight=cw, random_state=random_seed)
        if svm_cfg.calibrate_probabilities:
            cv = max(2, min(svm_cfg.calibration_cv, n_train or 100))
            return CalibratedClassifierCV(base, cv=cv, method="sigmoid")
        return base

    if model_type == "LOGREG":
        lr = sk.logistic_regression
        cw = lr.class_weight
        if cw not in (None, "balanced"):
            cw = None
        return LogisticRegression(C=lr.C,
                                  penalty=lr.penalty,
                                  solver=lr.solver,
                                  max_iter=lr.max_iter,
                                  class_weight=cw,
                                  random_state=random_seed,
                                  n_jobs=lr.n_jobs,
                                  multi_class="auto")

    if model_type == "RF":
        rf = sk.random_forest
        cw = rf.class_weight
        if cw not in (None, "balanced"):
            cw = None
        return RandomForestClassifier(n_estimators=rf.n_estimators,
                                      max_depth=rf.max_depth,
                                      class_weight=cw, random_state=rf.random_state if rf.random_state else random_seed,
                                      n_jobs=rf.n_jobs)

    if model_type == "XGBOOST":
        try:
            import xgboost as xgb
        except ImportError as e:
            raise ImportError("XGBOOST requer 'xgboost'. Instale com: pip install xgboost") from e
        xgb_cfg = sk.xgboost
        return xgb.XGBClassifier(n_estimators=xgb_cfg.n_estimators,
                                   max_depth=xgb_cfg.max_depth,
                                   learning_rate=xgb_cfg.learning_rate,
                                   subsample=xgb_cfg.subsample,
                                   colsample_bytree=xgb_cfg.colsample_bytree,
                                   tree_method=xgb_cfg.tree_method,
                                   eval_metric=xgb_cfg.eval_metric,
                                   random_state=xgb_cfg.random_state if xgb_cfg.random_state else random_seed,
                                   n_jobs=xgb_cfg.n_jobs)

    raise ValueError(f"Tipo sklearn não suportado: {model_type}. Use SVM, LOGREG, RF ou XGBOOST.")


def sklearn_predict_labels(loader: DataLoader, scaler: Any, pca: Any, clf: Any) -> Tuple[np.ndarray, np.ndarray]:
    """Gera predições e targets para todo o DataLoader."""
    preds, ys = [], []
    for features, targets, _idx in loader:
        X = scaler.transform(sklearn_flatten_batch(features))
        Xr = pca.transform(X)
        preds.append(np.asarray(clf.predict(Xr)).reshape(-1))
        ys.append(targets.detach().cpu().numpy().reshape(-1))
    if not preds:
        return np.empty((0,), dtype=np.int64), np.empty((0,), dtype=np.int64)
    return np.concatenate(preds), np.concatenate(ys)


def sklearn_metrics_dict(y_true: np.ndarray, y_pred: np.ndarray, full_dataset: Any) -> Dict[str, Any]:
    """Calcula métricas padronizadas (mesmas chaves que Tester.test())."""
    labels = list(range(full_dataset.get_num_classes()))
    target_names = [full_dataset.idx_to_target[i] for i in labels]
    p, r, f1, _ = precision_recall_fscore_support(y_true, y_pred, average="weighted", zero_division=0)
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(p), "recall": float(r), "f1": float(f1),
        "confusion_matrix": confusion_matrix(y_true, y_pred, labels=labels).tolist(),
        "classification_report": classification_report(y_true, y_pred, labels=labels, target_names=target_names, zero_division=0),
    }


def run_sklearn_eval_and_save(results: Dict, experiment_dir: Path, dataset_name: str,
                              wandb_run: Optional[Any], split_name: str) -> None:
    """Salva resultados JSON e loga métricas no W&B."""
    serializable = {}
    for k, v in results.items():
        serializable[k] = v.tolist() if isinstance(v, np.ndarray) else v
    json_file = experiment_dir / f"{dataset_name}_results.json"
    save_results_json(serializable, json_file, console)
    save_classification_plots(serializable, experiment_dir / "plots", dataset_name, console)
    print_classification_metrics(serializable, f"📊 {dataset_name}", console)
    console.print(
        f"[green]✓ {dataset_name}: "
        f"acc={results['accuracy']:.4f} "
        f"precision={results['precision']:.4f} "
        f"recall={results['recall']:.4f} "
        f"f1={results['f1']:.4f}[/green]"
    )
    console.print(f"[green]✓ Resultados de {dataset_name} salvos em: {json_file}[/green]")
    if wandb_run:
        try:
            wandb_run.log({f"{split_name}_{k}": results[k] for k in ("accuracy", "precision", "recall", "f1")})
        except Exception:
            pass


def load_sklearn_baseline_artifact(experiment_dir: Path) -> Dict[str, Any]:
    """Carrega artefato sklearn treinado (.joblib)."""
    import sys

    predictor_dir = Path(__file__).resolve().parents[5] / "legacy" / "neural_ancestry_predictor_deprecated"
    if str(predictor_dir) not in sys.path:
        sys.path.insert(0, str(predictor_dir))

    path = experiment_dir / "models" / SKLEARN_ARTIFACT_FILENAME
    if not path.exists():
        raise FileNotFoundError(f"Artefato sklearn não encontrado: {path}")
    data = joblib.load(path)
    if not isinstance(data, dict) or "classifier" not in data:
        raise ValueError(f"Artefato inválido em {path}")
    if "pca_cache_dir" in data:
        bundle = joblib.load(Path(data["pca_cache_dir"]) / SCALER_PCA_FILENAME)
        data = {**data, "scaler": bundle["scaler"], "pca": bundle["pca"]}
    elif "scaler" not in data or "pca" not in data:
        raise ValueError(f"Artefato inválido: faltam scaler/pca")
    return data


def train_sklearn_baseline(config: PipelineConfig, model_type: str, train_loader: DataLoader,
                           val_loader: DataLoader, test_loader: DataLoader,
                           full_dataset: Any, experiment_dir: Path,
                           wandb_run: Optional[Any] = None) -> Dict:
    """
    Treina pipeline sklearn completo: StandardScaler → PCA → classificador.

    Usa cache de PCA em disco (se configurado) para acelerar experimentos
    subsequentes. Salva artefato em experiment_dir/models/sklearn_baseline.joblib.
    """
    from genomics.core.sklearn_pca_cache import (
        ensure_sklearn_pca_cache, METADATA_FILENAME as SKLEARN_PCA_METADATA_FILENAME,
        fit_standard_scaler_incremental, fit_incremental_pca_on_train,
        fit_streaming_randomized_pca_on_train, stack_scaled_pca_batches, compute_sklearn_pca_effective_k,
    )

    from genomics.predictors.genotype_based.config import get_dataset_cache_dir

    if config.output.prediction_target == "frog_likelihood":
        raise ValueError("Baselines sklearn suportam apenas classificação.")

    model_type = model_type.upper()
    sk = config.model.sklearn
    random_seed = _effective_seed(config)
    use_pca_cache = sk.use_pca_cache
    force_pca = config.debug.force_pca_cache_rebuild

    first = next(iter(train_loader))
    n_features = int(np.prod(first[0].shape[1:]))
    n_train = len(train_loader.dataset)

    models_dir = experiment_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = models_dir / SKLEARN_ARTIFACT_FILENAME

    if use_pca_cache:
        # A função de cache legada recebe config serializado e o diretório do dataset processado.
        dataset_cache_dir = get_dataset_cache_dir(config)
        pca_dir = ensure_sklearn_pca_cache(config.model_dump(), dataset_cache_dir, train_loader, val_loader,
                                            test_loader, force=force_pca, log=console.print, rich_console=console)
        with open(pca_dir / SKLEARN_PCA_METADATA_FILENAME) as f:
            pca_meta = json.load(f)
        k = int(pca_meta["pca_n_components_effective"])

        X_train = np.load(pca_dir / "X_train.npy")
        y_train = np.load(pca_dir / "y_train.npy")
        valid = y_train >= 0
        X_train, y_train = X_train[valid], y_train[valid]

        clf = build_sklearn_classifier(config, model_type, random_seed, n_train=len(y_train))
        console.print(
            f"[cyan]Treinando {model_type}: "
            f"{len(y_train)} amostras x {X_train.shape[1]} componentes PCA...[/cyan]"
        )
        clf.fit(X_train, y_train)
        console.print(f"[green]✓ {model_type} treinado[/green]")
        joblib.dump({"classifier": clf, "model_type": model_type, "pca_n_components": k,
                      "pca_cache_dir": str(pca_dir.resolve()),
                      "dataset_cache_dir": str(dataset_cache_dir.resolve())}, artifact_path)
        console.print(f"[green]✓ Artefato sklearn salvo em: {artifact_path}[/green]")

        for name, xk, yk in [("val", "X_val", "y_val")]:
            Xs = np.load(pca_dir / f"{xk}.npy")
            yt = np.load(pca_dir / f"{yk}.npy")
            results = sklearn_metrics_dict(yt, clf.predict(Xs), full_dataset)
            run_sklearn_eval_and_save(results, experiment_dir, name, wandb_run, name)

        update_manifest(
            experiment_dir,
            status="completed",
            model_type=model_type,
            sklearn_artifact=str(artifact_path),
            pca_cache_dir=str(pca_dir.resolve()),
        )
        return {"model_type": model_type, "pca_n_components": k, "artifact_path": str(artifact_path)}

    # Sem cache: PCA in-memory
    k, _ = compute_sklearn_pca_effective_k(config.model_dump(), n_train=n_train, n_features=n_features, log=console.print)
    scaler = fit_standard_scaler_incremental(train_loader, rich_console=console)
    if sk.pca_backend == "randomized_streaming":
        pca = fit_streaming_randomized_pca_on_train(
            train_loader,
            scaler,
            k,
            models_dir,
            oversampling=sk.randomized_pca_oversampling,
            n_iter=sk.randomized_pca_n_iter,
            feature_chunk_size=sk.randomized_pca_feature_chunk_size,
            dtype=sk.randomized_pca_dtype,
            random_state=random_seed,
            log=console.print,
            rich_console=console,
        )
    else:
        pca = fit_incremental_pca_on_train(train_loader, scaler, k, log=console.print, rich_console=console)
    X_train, y_train = stack_scaled_pca_batches(train_loader, scaler, pca, rich_console=console)

    valid = y_train >= 0
    X_train, y_train = X_train[valid], y_train[valid]
    clf = build_sklearn_classifier(config, model_type, random_seed, n_train=len(y_train))
    console.print(
        f"[cyan]Treinando {model_type}: "
        f"{len(y_train)} amostras x {X_train.shape[1]} componentes PCA...[/cyan]"
    )
    clf.fit(X_train, y_train)
    console.print(f"[green]✓ {model_type} treinado[/green]")
    joblib.dump({"scaler": scaler, "pca": pca, "classifier": clf, "model_type": model_type, "pca_n_components": k}, artifact_path)
    console.print(f"[green]✓ Artefato sklearn salvo em: {artifact_path}[/green]")

    for name, loader in [("val", val_loader)]:
        y_pred, y_true = sklearn_predict_labels(loader, scaler, pca, clf)
        run_sklearn_eval_and_save(sklearn_metrics_dict(y_true, y_pred, full_dataset), experiment_dir, name, wandb_run, name)

    update_manifest(
        experiment_dir,
        status="completed",
        model_type=model_type,
        sklearn_artifact=str(artifact_path),
        pca_n_components=k,
    )
    return {"model_type": model_type, "pca_n_components": k, "artifact_path": str(artifact_path)}


def run_sklearn_test_mode(config: PipelineConfig, train_loader: DataLoader, val_loader: DataLoader,
                          test_loader: DataLoader, full_dataset: Any,
                          experiment_dir: Path, wandb_run: Optional[Any] = None,
                          output_name: Optional[str] = None) -> None:
    """Avaliação em modo test usando artefato joblib salvo."""
    if config.output.prediction_target == "frog_likelihood":
        raise ValueError("Baselines sklearn suportam apenas classificação.")
    art = load_sklearn_baseline_artifact(experiment_dir)
    scaler, pca, clf = art["scaler"], art["pca"], art["classifier"]
    choice = config.test_dataset.lower()
    loader_map = {"train": train_loader, "val": val_loader, "test": test_loader}
    loader = loader_map.get(choice, test_loader)
    y_pred, y_true = sklearn_predict_labels(loader, scaler, pca, clf)
    results = sklearn_metrics_dict(y_true, y_pred, full_dataset)
    run_sklearn_eval_and_save(results, experiment_dir, output_name or choice, wandb_run, choice)
