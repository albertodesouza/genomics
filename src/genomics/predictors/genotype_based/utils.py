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

import torch

from genomics.core.arrays import block_reduce_2d
from genomics.core.reproducibility import set_random_seeds, worker_init_fn


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
