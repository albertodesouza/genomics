# -*- coding: utf-8 -*-
"""
utils.py
========

Funções utilitárias usadas em todo o módulo `genotype_based_predictor`.

Conteúdo
--------
set_random_seeds    : Configura todas as seeds para reprodutibilidade
worker_init_fn      : Inicializa seeds de workers do DataLoader
block_reduce_2d     : Redimensionamento 2D com block-reduce/interpolação
taint_sample        : Debug – insere sinal sentinela em amostras
"""

import os
import random
from typing import Tuple

import numpy as np
import torch

from rich.console import Console

console = Console()


# ===========================================================================
# Reprodutibilidade
# ===========================================================================

def set_random_seeds(seed: int, strict_determinism: bool = True) -> None:
    """
    Configura todas as seeds randômicas para garantir reprodutibilidade.

    Afeta: Python random, NumPy, PyTorch CPU, PyTorch GPU (se disponível),
    variáveis de ambiente (PYTHONHASHSEED, CUBLAS_WORKSPACE_CONFIG).

    Parameters
    ----------
    seed : int
        Valor da seed randômica.
    strict_determinism : bool
        Se True, força determinismo total via `torch.use_deterministic_algorithms`.
        Pode deixar o treinamento 10-30% mais lento. Se False, permite pequenas
        variações em troca de melhor performance (~99% reprodutível).
    """
    os.environ["PYTHONHASHSEED"] = str(seed)

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.cuda.empty_cache()

    deterministic_algorithms_enabled = False

    if strict_determinism:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

        try:
            torch.use_deterministic_algorithms(True)
            deterministic_algorithms_enabled = True
        except AttributeError:
            torch.set_deterministic(True)
            deterministic_algorithms_enabled = True
        except Exception as e:
            console.print(
                "[yellow]Nao foi possivel ativar torch.use_deterministic_algorithms(True). "
                f"Continuando com determinismo parcial. Detalhe: {e}[/yellow]"
            )

        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"

        if deterministic_algorithms_enabled:
            console.print(
                f"[green]🎲 Seed configurada: {seed} "
                f"(determinismo ESTRITO — 100% reprodutível)[/green]"
            )
            console.print(
                "[yellow]   ⚠ Treinamento pode ser 10-30% mais lento "
                "devido ao determinismo estrito[/yellow]"
            )
        else:
            console.print(
                f"[green]🎲 Seed configurada: {seed} "
                f"(determinismo PARCIAL — ~99% reprodutível)[/green]"
            )
    else:
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True

        try:
            torch.use_deterministic_algorithms(False)
        except AttributeError:
            try:
                torch.set_deterministic(False)
            except AttributeError:
                pass
        except Exception:
            pass

        console.print(
            f"[green]🎲 Seed configurada: {seed} "
            f"(determinismo PARCIAL — ~99% reprodutível)[/green]"
        )


def worker_init_fn(worker_id: int) -> None:
    """
    Inicializa a seed de cada worker do DataLoader de forma determinística.

    Deve ser passada como `worker_init_fn` ao `DataLoader` quando
    `num_workers > 0` e `persistent_workers=True` para garantir
    reprodutibilidade entre execuções.

    Parameters
    ----------
    worker_id : int
        ID do worker (0, 1, 2, ...) — fornecido automaticamente pelo PyTorch.
    """
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


# ===========================================================================
# Processamento de arrays 2D
# ===========================================================================

def block_reduce_2d(
    arr: np.ndarray,
    target_shape: Tuple[int, int],
    func: str = "max",
) -> np.ndarray:
    """
    Redimensiona um array 2D para `target_shape` usando estratégias distintas
    para upscale (interpolação) e downscale (block-reduce), aplicadas
    independentemente em cada dimensão.

    Isso permite preservar picos de sinal (ex: picos de ATAC-seq) quando
    há downscale em uma dimensão e upscale em outra.

    Parameters
    ----------
    arr : np.ndarray
        Array 2D de entrada com shape (h, w).
    target_shape : Tuple[int, int]
        Shape desejada (target_h, target_w).
    func : str
        Função de agregação para downscale: 'max', 'min' ou 'mean'.

    Returns
    -------
    np.ndarray
        Array redimensionado com shape `target_shape`.

    Raises
    ------
    ValueError
        Se `func` não for 'max', 'min' ou 'mean'.
    """
    h, w = arr.shape
    th, tw = target_shape

    agg_map = {"max": np.max, "min": np.min, "mean": np.mean}
    if func not in agg_map:
        raise ValueError(
            f"Função de agregação desconhecida: '{func}'. Use 'max', 'min' ou 'mean'."
        )
    agg_func = agg_map[func]

    # Passo 1: LARGURA (colunas)
    if tw >= w:
        from scipy import ndimage
        temp = ndimage.zoom(arr, (1, tw / w), order=1)
    else:
        bw = w / tw
        temp = np.zeros((h, tw), dtype=arr.dtype)
        for j in range(tw):
            x_start = int(j * bw)
            x_end = min(int((j + 1) * bw), w)
            if x_end > x_start:
                temp[:, j] = agg_func(arr[:, x_start:x_end], axis=1)

    # Passo 2: ALTURA (linhas)
    h_temp = temp.shape[0]
    if th >= h_temp:
        from scipy import ndimage
        result = ndimage.zoom(temp, (th / h_temp, 1), order=1)
    else:
        bh = h_temp / th
        result = np.zeros((th, tw), dtype=arr.dtype)
        for i in range(th):
            y_start = int(i * bh)
            y_end = min(int((i + 1) * bh), h_temp)
            if y_end > y_start:
                result[i, :] = agg_func(temp[y_start:y_end, :], axis=0)

    return result


# ===========================================================================
# Debug / Taint
# ===========================================================================

def taint_sample(
    features: torch.Tensor,
    target_class: int,
    num_classes: int,
    taint_type: str = "additive",
    taint_value: float = 1.0,
    taint_horizontal_size: int = 6,
    taint_vertical_size: int = 6,
    taint_horizontal_step: int = 100,
    taint_vertical_step: int = 6,
) -> torch.Tensor:
    """
    Insere um sinal sentinela numa amostra para fins de debug.

    Permite verificar se a rede "vê" um sinal de classe antes de treinar
    com dados reais. Útil para confirmar que o pipeline de dados e a
    arquitetura do modelo estão corretos.

    Modos
    -----
    additive (padrão)
        Mantém os dados originais e adiciona o bloco sentinela sobre eles.
        Útil para verificar se a rede consegue detectar o sinal mesmo com
        ruído de fundo.
    override
        Zera toda a entrada e coloca apenas o bloco sentinela.
        Teste mais radical: a rede SÓ vê o taint. Se não aprender com
        isso, há problema estrutural na arquitetura ou no pipeline.

    Posicionamento do bloco
    -----------------------
    O bloco é posicionado em função da classe:
        start_col = target_class * taint_horizontal_step
        start_row = target_class * taint_vertical_step

    Parameters
    ----------
    features : torch.Tensor
        Tensor de entrada (1D, 2D ou mais dimensões).
    target_class : int
        Classe-alvo (0 a num_classes-1).
    num_classes : int
        Número total de classes.
    taint_type : str
        'additive' ou 'override'.
    taint_value : float
        Valor inserido no bloco sentinela.
    taint_horizontal_size : int
        Largura do bloco sentinela em colunas.
    taint_vertical_size : int
        Altura do bloco sentinela em linhas (apenas para dados 2D).
    taint_horizontal_step : int
        Deslocamento horizontal por classe.
    taint_vertical_step : int
        Deslocamento vertical por classe (apenas para dados 2D).

    Returns
    -------
    torch.Tensor
        Tensor com o sinal sentinela aplicado.
    """
    if taint_type == "override":
        features = torch.zeros_like(features)
    else:
        features = features.clone()

    if features.ndim == 1:
        input_size = features.shape[0]
        position = int(target_class * taint_horizontal_step)
        end_position = min(position + taint_horizontal_size, input_size)
        if position < input_size:
            features[position:end_position] = taint_value

    elif features.ndim == 2:
        num_rows, num_cols = features.shape
        start_col = int(target_class * taint_horizontal_step)
        start_row = int(target_class * taint_vertical_step)
        end_row = min(start_row + taint_vertical_size, num_rows)
        end_col = min(start_col + taint_horizontal_size, num_cols)
        if start_row < num_rows and start_col < num_cols:
            features[start_row:end_row, start_col:end_col] = taint_value

    else:
        flat = features.view(-1)
        input_size = flat.numel()
        position = int(target_class * taint_horizontal_step)
        end_position = min(position + taint_horizontal_size, input_size)
        if position < input_size:
            flat[position:end_position] = taint_value

    return features
