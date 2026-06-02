# -*- coding: utf-8 -*-
"""
data_splitting.py
=================

Funções para divisão de dados em conjuntos de treino/validação/teste com
consciência de grupos familiares (family-aware splitting).

O objetivo é garantir que indivíduos da mesma família permaneçam no mesmo
split, evitando data leakage por similaridade genética entre parentes.

Funções
-------
build_family_aware_sample_groups : Agrupa amostras por família
build_valid_sample_index_map     : Filtra amostras sem target válido
"""

from pathlib import Path
from typing import Any, Dict, List, Tuple

from genomics_pipeline.splitting import build_family_groups
from genotype_based_predictor.config import PipelineConfig


def _extract_family_links(pedigree: Dict) -> List[str]:
    from genomics_pipeline.splitting import extract_family_links

    return extract_family_links(pedigree)


def build_family_aware_sample_groups(
    base_dataset: Any,
    config: PipelineConfig,
) -> Tuple[List[List[int]], Dict[str, Any]]:
    """
    Constrói grupos de amostras que devem permanecer no mesmo split.

    Tenta identificar grupos familiares através de três estratégias, em
    ordem de prioridade:
    1. `family_id` presente em `individual_metadata.json` (mais direto)
    2. Links de pedigree (father/mother/etc.) via union-find
    3. Fallback: cada indivíduo é seu próprio grupo (singleton)

    O parâmetro `family_split_mode` no config controla o comportamento:
    - ``'family_aware'`` (padrão): agrupa por família
    - ``'ignore'``: trata cada indivíduo como grupo independente

    Parameters
    ----------
    base_dataset : GenomicLongevityDataset
        Dataset base com atributo `dataset_metadata` e `dataset_dir`.
    config : PipelineConfig
        Configuração completa do experimento.

    Returns
    -------
    Tuple[List[List[int]], Dict[str, Any]]
        - Lista de grupos, onde cada grupo é uma lista de índices de amostras
          que devem permanecer juntos no mesmo split.
        - Dicionário com informações sobre o agrupamento (para logging).
    """
    dataset_metadata = getattr(base_dataset, "dataset_metadata", {}) or {}
    individuals = dataset_metadata.get("individuals", [])
    pedigree_map = dataset_metadata.get("individuals_pedigree", {})
    if not individuals:
        individuals = [str(idx) for idx in range(len(base_dataset))]
    return build_family_groups(
        sample_ids=individuals,
        pedigree_map=pedigree_map,
        individuals_dir=Path(base_dataset.dataset_dir) / "individuals",
        family_split_mode=config.data_split.family_split_mode,
    )


def build_valid_sample_index_map(base_dataset: Any, config: PipelineConfig) -> List[int]:
    """
    Retorna os índices do dataset base que têm um target válido para a
    tarefa de predição configurada.

    Utiliza `ProcessedGenomicDataset` em modo "probe" (sem normalização)
    apenas para descobrir quais amostras têm target válido, sem processar
    os dados de predição.

    Parameters
    ----------
    base_dataset : GenomicLongevityDataset
        Dataset base.
    config : PipelineConfig
        Configuração completa do experimento.

    Returns
    -------
    List[int]
        Lista de índices no dataset base com target válido.
    """
    from genotype_based_predictor.dataset import ProcessedGenomicDataset

    probe_dataset = ProcessedGenomicDataset(
        base_dataset=base_dataset,
        config=config,
        normalization_params={"mean": 0.0, "std": 1.0},
        compute_normalization=False,
    )
    return list(probe_dataset.valid_sample_indices)
