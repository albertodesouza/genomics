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

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

from rich.console import Console

from genotype_based_predictor.config import PipelineConfig

console = Console()


def _extract_family_links(pedigree: Dict) -> List[str]:
    """
    Extrai IDs de parentes a partir de um dicionário de pedigree.

    Procura por chaves que contenham termos relacionados a parentesco
    (father, mother, parent, child, sibling, family).

    Parameters
    ----------
    pedigree : Dict
        Dicionário com metadados de pedigree de um indivíduo.

    Returns
    -------
    List[str]
        Lista de IDs de indivíduos relacionados.
    """
    related = []
    for key, value in pedigree.items():
        key_lower = str(key).lower()
        if not any(
            token in key_lower
            for token in ["father", "mother", "parent", "child", "sibling", "family"]
        ):
            continue
        if value in [None, "", "0"]:
            continue
        if isinstance(value, (list, tuple, set)):
            related.extend(str(v) for v in value if v not in [None, "", "0"])
        else:
            related.append(str(value))
    return related


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
    family_mode = config.data_split.family_split_mode

    dataset_metadata = getattr(base_dataset, "dataset_metadata", {}) or {}
    individuals = dataset_metadata.get("individuals", [])
    pedigree_map = dataset_metadata.get("individuals_pedigree", {})

    # Sem lista de indivíduos: tratar como singletons
    if not individuals:
        groups = [[idx] for idx in range(len(base_dataset))]
        return groups, {
            "family_split_mode": family_mode,
            "num_groups": len(groups),
            "num_family_groups": 0,
            "num_singletons": len(groups),
            "grouping_source": "individual",
        }

    sample_to_idx = {sample_id: idx for idx, sample_id in enumerate(individuals)}

    # Tentativa 1: usar family_id de individual_metadata.json
    family_ids: Dict[str, str] = {}
    individuals_dir = Path(base_dataset.dataset_dir) / "individuals"

    for sample_id in individuals:
        metadata_file = individuals_dir / sample_id / "individual_metadata.json"
        if not metadata_file.exists():
            continue
        try:
            with open(metadata_file, "r") as f:
                individual_metadata = json.load(f)
            family_id = individual_metadata.get("family_id")
            if family_id not in [None, "", "0"]:
                family_ids[sample_id] = str(family_id)
        except Exception:
            continue

    if family_ids:
        groups_by_family_id: Dict[str, List[int]] = {}
        for sample_id, idx in sample_to_idx.items():
            fid = family_ids.get(sample_id, sample_id)
            groups_by_family_id.setdefault(fid, []).append(idx)

        if family_mode != "ignore":
            groups = list(groups_by_family_id.values())
            return groups, {
                "family_split_mode": family_mode,
                "grouping_source": "family_id",
                "num_groups": len(groups),
                "num_family_groups": sum(1 for g in groups if len(g) > 1),
                "num_singletons": sum(1 for g in groups if len(g) == 1),
            }

    # Modo ignore: singletons
    if family_mode == "ignore":
        groups = [[idx] for idx in range(len(individuals))]
        return groups, {
            "family_split_mode": family_mode,
            "grouping_source": "individual",
            "num_groups": len(groups),
            "num_family_groups": 0,
            "num_singletons": len(groups),
        }

    # Tentativa 2: union-find sobre links de pedigree
    parent = {sample_id: sample_id for sample_id in individuals}

    def find(sample_id: str) -> str:
        while parent[sample_id] != sample_id:
            parent[sample_id] = parent[parent[sample_id]]
            sample_id = parent[sample_id]
        return sample_id

    def union(left: str, right: str) -> None:
        if left not in parent or right not in parent:
            return
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for sample_id, pedigree in pedigree_map.items():
        if sample_id not in parent:
            continue
        for related_id in _extract_family_links(pedigree):
            if related_id in parent:
                union(sample_id, related_id)

    groups_by_root: Dict[str, List[int]] = {}
    for sample_id, idx in sample_to_idx.items():
        root = find(sample_id)
        groups_by_root.setdefault(root, []).append(idx)

    groups = list(groups_by_root.values())
    return groups, {
        "family_split_mode": family_mode,
        "grouping_source": "pedigree",
        "num_groups": len(groups),
        "num_family_groups": sum(1 for g in groups if len(g) > 1),
        "num_singletons": sum(1 for g in groups if len(g) == 1),
    }


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
