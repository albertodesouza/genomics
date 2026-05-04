# -*- coding: utf-8 -*-
"""
normalization.py
================

Funções de normalização aplicadas às features genômicas antes de entrar
nos modelos. Todas as funções operam sobre tensores PyTorch e preservam
o device original (CPU ou GPU).

Funções disponíveis
-------------------
zscore_normalize        : Normalização Z-score (µ=0, σ=1)
minmax_keep_zero        : Min-Max preservando zeros exatos
log_normalize           : Normalização logarítmica via log1p

Notas de design
---------------
- Parâmetros de normalização (mean, std, max, log_max) são sempre
  pre-computados APENAS nos dados de treino+validação para evitar
  data leakage com o conjunto de teste.
- As funções são stateless: recebem o tensor e os parâmetros e retornam
  o tensor normalizado, sem efeitos colaterais.
"""

from typing import Optional

import torch


def zscore_normalize(x: torch.Tensor, mean: float, std: float) -> torch.Tensor:
    """
    Normalização Z-score padrão: (x - µ) / σ.

    Se std == 0 (feature constante), retorna x sem modificação para evitar
    divisão por zero.

    Parameters
    ----------
    x : torch.Tensor
        Tensor de entrada (qualquer shape, qualquer device).
    mean : float
        Média global pré-computada no conjunto de treino+validação.
    std : float
        Desvio padrão global pré-computado no conjunto de treino+validação.

    Returns
    -------
    torch.Tensor
        Tensor normalizado com mesma shape e device de `x`.
    """
    if std == 0:
        return x
    return (x - mean) / std


def minmax_keep_zero(x: torch.Tensor, xmax: float) -> torch.Tensor:
    """
    Normalização Min-Max onde zeros são preservados como zeros.

    Divide todos os valores por `xmax` (máximo dos valores não-zero).
    Zeros verdadeiros permanecem como 0.0 após a normalização, o que é
    útil para dados esparsos (ex: ATAC-seq) onde ausência de sinal tem
    significado biológico.

    Parameters
    ----------
    x : torch.Tensor
        Tensor de entrada.
    xmax : float
        Máximo dos valores não-zero pré-computado. Se 0, retorna x intacto.

    Returns
    -------
    torch.Tensor
        Tensor normalizado no intervalo [0, 1] (exceto se houver valores
        acima do máximo usado no fit).
    """
    if xmax == 0:
        return x
    return x / xmax


def log_normalize(x: torch.Tensor, log_max: float) -> torch.Tensor:
    """
    Normalização logarítmica via log1p, dividida pelo log1p do máximo.

    Transforma x → log1p(x) / log1p(max), resultando em valores em [0, 1].
    Zeros permanecem zero automaticamente porque log1p(0) = 0.

    Útil para dados com distribuição altamente assimétrica (ex: RNA-seq
    com alguns genes muito expressos).

    Parameters
    ----------
    x : torch.Tensor
        Tensor de entrada com valores ≥ 0.
    log_max : float
        log1p do valor máximo pré-computado. Se 0, retorna x intacto.

    Returns
    -------
    torch.Tensor
        Tensor normalizado.
    """
    if log_max == 0:
        return x
    return torch.log1p(x) / log_max


def apply_normalization(
    x: torch.Tensor,
    normalization_params: dict,
    track_idx: Optional[int] = None,
) -> torch.Tensor:
    """
    Aplica o método de normalização correto com base em `normalization_params`.

    Suporta tanto normalização global (um conjunto de parâmetros para todas
    as tracks) quanto normalização por-track (parâmetros independentes por track).

    Parameters
    ----------
    x : torch.Tensor
        Tensor a normalizar.
    normalization_params : dict
        Dicionário com 'method' e os parâmetros correspondentes.
        Formato global:   {'method': 'zscore', 'mean': ..., 'std': ...}
        Formato per-track: {'method': 'zscore', 'per_track': True,
                            'track_params': [{'mean': ..., 'std': ...}, ...]}
    track_idx : int | None
        Índice da track quando usando normalização per-track.
        Ignorado em normalização global.

    Returns
    -------
    torch.Tensor
        Tensor normalizado.
    """
    method = normalization_params.get("method", "zscore")
    per_track = normalization_params.get("per_track", False)

    if per_track and track_idx is not None:
        params = normalization_params["track_params"][track_idx]
    else:
        params = normalization_params

    if method == "zscore":
        return zscore_normalize(x, params.get("mean", 0.0), params.get("std", 1.0))
    elif method == "minmax_keep_zero":
        return minmax_keep_zero(x, params.get("max", 1.0))
    elif method == "log":
        return log_normalize(x, params.get("log_max", 1.0))
    else:
        # Fallback para zscore
        return zscore_normalize(x, params.get("mean", 0.0), params.get("std", 1.0))
