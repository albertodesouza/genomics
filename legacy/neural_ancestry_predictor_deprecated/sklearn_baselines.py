#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sklearn baselines for neural ancestry experiments."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import json

import joblib
import numpy as np
from rich.console import Console
from sklearn.calibration import CalibratedClassifierCV
from sklearn.decomposition import IncrementalPCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, precision_recall_fscore_support
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from torch.utils.data import DataLoader

from genomics.core.sklearn_pca_cache import (
    METADATA_FILENAME as SKLEARN_PCA_METADATA_FILENAME,
    SCALER_PCA_FILENAME,
    compute_sklearn_pca_effective_k,
    ensure_sklearn_pca_cache,
    fit_incremental_pca_on_train,
    fit_standard_scaler_incremental,
    sklearn_flatten_batch,
    stack_scaled_pca_batches,
)
from neural_ancestry_predictor_deprecated.config import get_dataset_cache_dir


console = Console()

SKLEARN_BASELINE_TYPES = frozenset({'SVM', 'RF', 'XGBOOST'})
SKLEARN_ARTIFACT_FILENAME = 'sklearn_baseline.joblib'


def _sklearn_effective_random_seed(config: Dict) -> int:
    seed = config['data_split'].get('random_seed')
    if seed is None or seed == -1:
        return 42
    return int(seed)


def _ensure_sklearn_classification_target(config: Dict) -> None:
    if config['output'].get('prediction_target') == 'frog_likelihood':
        raise ValueError(
            "Baselines sklearn (SVM/RF/XGBOOST) suportam apenas classificação. "
            "Altere output.prediction_target ou use um modelo NN/CNN/CNN2."
        )


def sklearn_predict_labels(
    loader: DataLoader,
    scaler: StandardScaler,
    pca: IncrementalPCA,
    clf: Any
) -> Tuple[np.ndarray, np.ndarray]:
    """Predições em todo o loader (rótulos int)."""
    preds: List[np.ndarray] = []
    ys: List[np.ndarray] = []
    for features, targets, _idx in loader:
        X = scaler.transform(sklearn_flatten_batch(features))
        Xr = pca.transform(X)
        pr = clf.predict(Xr)
        preds.append(np.asarray(pr).reshape(-1))
        t = targets.detach().cpu().numpy()
        ys.append(t.reshape(-1))
    if not preds:
        return np.empty((0,), dtype=np.int64), np.empty((0,), dtype=np.int64)
    return np.concatenate(preds), np.concatenate(ys)


def build_sklearn_classifier(
    config: Dict,
    model_type: str,
    random_seed: int,
    n_train_samples: Optional[int] = None
) -> Any:
    """Instancia o classificador final (após PCA)."""
    model_type = model_type.upper()
    sk = config['model'].get('sklearn', {})

    if model_type == 'SVM':
        svm_cfg = sk.get('svm', {})
        cw = svm_cfg.get('class_weight')
        if cw is not None and cw != 'balanced':
            cw = None
        base = LinearSVC(
            C=float(svm_cfg.get('C', 1.0)),
            max_iter=int(svm_cfg.get('max_iter', 20000)),
            dual=False,
            class_weight=cw,
            random_state=random_seed,
        )
        if svm_cfg.get('calibrate_probabilities', False):
            if n_train_samples is not None and n_train_samples < 3:
                console.print(
                    "[yellow]⚠ calibrate_probabilities ignorado: n_train < 3 (use LinearSVC sem calibração)[/yellow]"
                )
                return base
            cv = int(svm_cfg.get('calibration_cv', 3))
            if n_train_samples is not None:
                cv = max(2, min(cv, n_train_samples))
            return CalibratedClassifierCV(base, cv=cv, method='sigmoid')
        return base

    if model_type == 'RF':
        rf = sk.get('random_forest', {})
        cw = rf.get('class_weight')
        if cw is not None and cw != 'balanced':
            cw = None
        md = rf.get('max_depth')
        return RandomForestClassifier(
            n_estimators=int(rf.get('n_estimators', 200)),
            max_depth=None if md is None else int(md),
            class_weight=cw,
            random_state=int(rf.get('random_state', random_seed)),
            n_jobs=int(rf.get('n_jobs', -1)),
        )

    if model_type == 'XGBOOST':
        try:
            import xgboost as xgb
        except ImportError as e:
            raise ImportError(
                "XGBOOST requer o pacote 'xgboost'. Instale com: pip install xgboost"
            ) from e
        xgb_cfg = sk.get('xgboost', {})
        return xgb.XGBClassifier(
            n_estimators=int(xgb_cfg.get('n_estimators', 200)),
            max_depth=int(xgb_cfg.get('max_depth', 6)),
            learning_rate=float(xgb_cfg.get('learning_rate', 0.1)),
            subsample=float(xgb_cfg.get('subsample', 0.8)),
            colsample_bytree=float(xgb_cfg.get('colsample_bytree', 0.8)),
            tree_method=str(xgb_cfg.get('tree_method', 'hist')),
            random_state=int(xgb_cfg.get('random_state', random_seed)),
            n_jobs=int(xgb_cfg.get('n_jobs', -1)),
            eval_metric=str(xgb_cfg.get('eval_metric', 'mlogloss')),
        )

    raise ValueError(f"Tipo sklearn baseline não suportado: {model_type}")


def _print_svm_convergence(clf: Any, max_iter: int) -> None:
    """Imprime se o LinearSVC convergiu ou atingiu o limite de iterações."""
    base = clf
    if hasattr(clf, 'calibrated_classifiers_'):
        try:
            base = clf.calibrated_classifiers_[0].estimator
        except (AttributeError, IndexError):
            pass
    elif hasattr(clf, 'estimator'):
        base = clf.estimator

    if not hasattr(base, 'n_iter_'):
        console.print("[yellow]  ⚠ SVM: informação de convergência não disponível[/yellow]")
        return

    n_iter = int(np.max(base.n_iter_))
    if n_iter < max_iter:
        console.print(
            f"[green]  ✓ SVM convergiu em {n_iter} iterações "
            f"(limite: {max_iter})[/green]"
        )
    else:
        console.print(
            f"[bold red]  ✗ SVM NÃO convergiu! Atingiu o limite de {max_iter} iterações. "
            f"Aumente max_iter ou reduza C.[/bold red]"
        )


def sklearn_metrics_dict(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    full_dataset: Any
) -> Dict[str, Any]:
    """Mesmas chaves que Tester.test() para classificação."""
    labels = list(range(full_dataset.get_num_classes()))
    target_names = [full_dataset.idx_to_target[i] for i in labels]
    results: Dict[str, Any] = {}
    results['accuracy'] = float(accuracy_score(y_true, y_pred))
    results['precision'], results['recall'], results['f1'], _ = precision_recall_fscore_support(
        y_true, y_pred, average='weighted', zero_division=0
    )
    results['precision'] = float(results['precision'])
    results['recall'] = float(results['recall'])
    results['f1'] = float(results['f1'])
    results['confusion_matrix'] = confusion_matrix(y_true, y_pred, labels=labels)
    results['classification_report'] = classification_report(
        y_true, y_pred, labels=labels, target_names=target_names, zero_division=0
    )
    return results


def train_sklearn_baseline(
    config: Dict,
    model_type: str,
    train_loader: DataLoader,
    val_loader: DataLoader,
    test_loader: DataLoader,
    full_dataset: Any,
    experiment_dir: Path,
    wandb_run: Optional[Any] = None
) -> Dict[str, Any]:
    """
    Treina classificador sklearn após redução PCA.

    Com ``model.sklearn.use_pca_cache`` (default True), StandardScaler + IncrementalPCA
    e matrizes reduzidas são gravados em disco sob ``processed_cache_dir/pca_cache/``,
    reutilizáveis entre experimentos. Caso contrário, PCA é ajustado só em memória
    nesta execução (comportamento antigo).

    Salva artefato em experiment_dir/models/sklearn_baseline.joblib
    """
    _ensure_sklearn_classification_target(config)
    model_type = model_type.upper()
    sk = config['model'].get('sklearn', {})
    random_seed = _sklearn_effective_random_seed(config)
    use_pca_cache = sk.get('use_pca_cache', True)
    force_pca = config.get('debug', {}).get('force_pca_cache_rebuild', False)

    first = next(iter(train_loader))
    features0 = first[0]
    n_features = int(np.prod(features0.shape[1:]))
    n_train = len(train_loader.dataset)
    align_n_train = bool(sk.get('pca_align_n_train', False))

    models_dir = experiment_dir / 'models'
    models_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = models_dir / SKLEARN_ARTIFACT_FILENAME

    if use_pca_cache:
        dataset_cache_dir = Path(get_dataset_cache_dir(config))
        pca_dir = ensure_sklearn_pca_cache(
            config,
            dataset_cache_dir,
            train_loader,
            val_loader,
            test_loader,
            force=force_pca,
            log=console.print,
            rich_console=console,
        )
        with open(pca_dir / SKLEARN_PCA_METADATA_FILENAME, 'r') as f:
            pca_meta = json.load(f)
        effective_k = int(pca_meta['pca_n_components_effective'])
        pca_req = int(pca_meta.get('pca_components_requested', effective_k))

        console.print(
            f"[cyan]Sklearn baseline: {model_type} | PCA k={effective_k} "
            f"(pedido={pca_req}, cache={pca_dir.name})[/cyan]"
        )

        X_train = np.load(pca_dir / 'X_train.npy')
        y_train = np.load(pca_dir / 'y_train.npy')

        console.print("[cyan]Passo classificador: fit em matrizes do cache PCA...[/cyan]")
        valid_mask = y_train >= 0
        if not np.all(valid_mask):
            console.print(f"[yellow]⚠ Removendo {int((~valid_mask).sum())} labels inválidos do treino PCA[/yellow]")
            X_train = X_train[valid_mask]
            y_train = y_train[valid_mask]
        if len(np.unique(y_train)) < 2:
            raise ValueError(
                f"Treino inválido após filtrar labels: classes presentes={np.unique(y_train).tolist()}. "
                "Isso indica cache/split incorreto ou apenas uma classe no treino."
            )
        clf = build_sklearn_classifier(
            config, model_type, random_seed, n_train_samples=len(y_train)
        )
        clf.fit(X_train, y_train)
        if model_type == 'SVM':
            _print_svm_convergence(clf, int(sk.get('svm', {}).get('max_iter', 20000)))

        artifact = {
            'classifier': clf,
            'model_type': model_type,
            'pca_n_components': effective_k,
            'pca_cache_dir': str(pca_dir.resolve()),
        }
        joblib.dump(artifact, artifact_path)
        console.print(f"[green]✓ Artefato sklearn salvo em {artifact_path}[/green]")

        history = {
            'model_type': model_type,
            'pca_n_components': effective_k,
            'pca_cache_dir': str(pca_dir.resolve()),
            'artifact_path': str(artifact_path),
        }

        console.print("[cyan]Avaliação train/val/test (vetores do cache PCA)...[/cyan]")
        for name, x_key, y_key in [
            ('train', 'X_train', 'y_train'),
            ('val', 'X_val', 'y_val'),
            ('test', 'X_test', 'y_test'),
        ]:
            Xs = np.load(pca_dir / f'{x_key}.npy')
            y_true = np.load(pca_dir / f'{y_key}.npy')
            y_pred = clf.predict(Xs)
            results = sklearn_metrics_dict(y_true, y_pred, full_dataset)
            run_sklearn_eval_and_save(results, experiment_dir, name, wandb_run, split_name=name)

        return history

    effective_k, pca_req = compute_sklearn_pca_effective_k(
        config, n_train=n_train, n_features=n_features, log=console.print
    )

    console.print(
        f"[cyan]Sklearn baseline: {model_type} | PCA k={effective_k} "
        f"(pedido={pca_req}, n_train={n_train}, D={n_features}, sem pca_cache)[/cyan]"
    )

    scaler = fit_standard_scaler_incremental(
        train_loader,
        rich_console=console,
        progress_desc="Sklearn baseline: StandardScaler (1/4)",
    )

    pca = fit_incremental_pca_on_train(
        train_loader,
        scaler,
        effective_k,
        log=console.print,
        forbid_tail_padding=align_n_train,
        rich_console=console,
        progress_desc="Sklearn baseline: IncrementalPCA (2/4)",
    )

    console.print("[cyan]Passo 3/4: Montando matriz de treino reduzida e fit do classificador...[/cyan]")
    X_train, y_train = stack_scaled_pca_batches(
        train_loader,
        scaler,
        pca,
        rich_console=console,
        progress_desc="Sklearn baseline: PCA transform treino (3/4)",
    )
    valid_mask = y_train >= 0
    if not np.all(valid_mask):
        console.print(f"[yellow]⚠ Removendo {int((~valid_mask).sum())} labels inválidos do treino[/yellow]")
        X_train = X_train[valid_mask]
        y_train = y_train[valid_mask]
    if len(np.unique(y_train)) < 2:
        raise ValueError(
            f"Treino inválido após filtrar labels: classes presentes={np.unique(y_train).tolist()}. "
            "Isso indica split incorreto ou apenas uma classe no treino."
        )
    clf = build_sklearn_classifier(
        config, model_type, random_seed, n_train_samples=len(y_train)
    )
    clf.fit(X_train, y_train)
    if model_type == 'SVM':
        _print_svm_convergence(clf, int(sk.get('svm', {}).get('max_iter', 20000)))

    artifact = {
        'scaler': scaler,
        'pca': pca,
        'classifier': clf,
        'model_type': model_type,
        'pca_n_components': effective_k,
    }
    joblib.dump(artifact, artifact_path)
    console.print(f"[green]✓ Artefato sklearn salvo em {artifact_path}[/green]")

    history = {
        'model_type': model_type,
        'pca_n_components': effective_k,
        'artifact_path': str(artifact_path),
    }

    console.print("[cyan]Passo 4/4: Avaliação train/val/test...[/cyan]")
    for name, loader in [('train', train_loader), ('val', val_loader), ('test', test_loader)]:
        y_pred, y_true = sklearn_predict_labels(loader, scaler, pca, clf)
        results = sklearn_metrics_dict(y_true, y_pred, full_dataset)
        run_sklearn_eval_and_save(results, experiment_dir, name, wandb_run, split_name=name)

    return history


def run_sklearn_eval_and_save(
    results: Dict[str, Any],
    experiment_dir: Path,
    dataset_name: str,
    wandb_run: Optional[Any],
    split_name: str
) -> None:
    """Serializa JSON como run_test_and_save; log opcional no W&B."""
    results_serializable = {}
    for key, value in results.items():
        if isinstance(value, np.ndarray):
            results_serializable[key] = value.tolist()
        elif isinstance(value, (list, tuple)):
            results_serializable[key] = [
                item.tolist() if isinstance(item, np.ndarray) else item
                for item in value
            ]
        else:
            results_serializable[key] = value

    json_file = experiment_dir / f'{dataset_name}_results.json'
    with open(json_file, 'w') as f:
        json.dump(results_serializable, f, indent=2)
    console.print(f"[green]✓ Resultados de {dataset_name} salvos em: {json_file}[/green]\n")

    if wandb_run:
        try:
            wandb_run.log({
                f'{split_name}_accuracy': results['accuracy'],
                f'{split_name}_precision': results['precision'],
                f'{split_name}_recall': results['recall'],
                f'{split_name}_f1': results['f1'],
            })
        except Exception:
            pass


def load_sklearn_baseline_artifact(experiment_dir: Path) -> Dict[str, Any]:
    path = experiment_dir / 'models' / SKLEARN_ARTIFACT_FILENAME
    if not path.exists():
        raise FileNotFoundError(f"Artefato sklearn não encontrado: {path}")
    data = joblib.load(path)
    if not isinstance(data, dict) or 'classifier' not in data:
        raise ValueError(f"Artefato sklearn inválido em {path}")
    if 'pca_cache_dir' in data:
        bundle_path = Path(data['pca_cache_dir']) / SCALER_PCA_FILENAME
        if not bundle_path.exists():
            raise FileNotFoundError(f"Cache PCA ausente (esperado {bundle_path})")
        bundle = joblib.load(bundle_path)
        data = {**data, 'scaler': bundle['scaler'], 'pca': bundle['pca']}
    elif 'scaler' not in data or 'pca' not in data:
        raise ValueError(f"Artefato sklearn inválido em {path} (faltam scaler/pca)")
    return data


def run_sklearn_test_mode(
    config: Dict,
    train_loader: DataLoader,
    val_loader: DataLoader,
    test_loader: DataLoader,
    full_dataset: Any,
    experiment_dir: Path,
    wandb_run: Optional[Any] = None
) -> None:
    """Avaliação em modo --mode test usando artefato joblib."""
    _ensure_sklearn_classification_target(config)
    art = load_sklearn_baseline_artifact(experiment_dir)
    scaler = art['scaler']
    pca = art['pca']
    clf = art['classifier']

    test_dataset_choice = config.get('test_dataset', 'test').lower()
    if test_dataset_choice == 'train':
        selected_loader = train_loader
        dataset_name = 'train'
        label = 'Train'
    elif test_dataset_choice == 'val':
        selected_loader = val_loader
        dataset_name = 'val'
        label = 'Validation'
    else:
        selected_loader = test_loader
        dataset_name = 'test'
        label = 'Test'

    console.print(f"[cyan]Teste sklearn no conjunto: {label}[/cyan]")
    y_pred, y_true = sklearn_predict_labels(selected_loader, scaler, pca, clf)
    results = sklearn_metrics_dict(y_true, y_pred, full_dataset)
    run_sklearn_eval_and_save(results, experiment_dir, dataset_name, wandb_run, split_name=dataset_name)
