"""Family-aware train/test split for this notebook's individuals.

Same grouping approach `genomics.core.splitting.build_family_groups`/`split_groups` use
elsewhere in this repo (group related individuals by `family_id`, shuffle whole groups with
a seeded RNG, then slice by fraction) -- reimplemented here, rather than imported, since this
notebook is deliberately kept independent of the rest of the `genomics` package (see the
notebook's own setup cell).

1000 Genomes trios/families share a `family_id` (see `notebooks/data/individuals.json`); a
plain per-individual random split would let relatives straddle the train/test boundary,
leaking genotype (and therefore RNA-seq signal) between the two.
"""

from __future__ import annotations

import random
from typing import Dict, List, Tuple

from .data import population_to_class


def family_aware_train_test_split(
    individuals: Dict[str, dict],
    class_map: Dict[str, List[str]],
    test_size: float = 0.2,
    random_seed: int = 13,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Splits individuals into train/test, keeping every member of the same family
    (`individuals[sample_id]["family_id"]`) on the same side, and grouping/splitting each
    *population* separately (not just each class) so the population mix within a class is
    preserved across train and test. A split done at the class level alone can, purely by
    chance, shift (e.g.) YRI toward train and GWD toward test within "strong pigmentation" --
    family clustering already shrinks the effective sample size per split, and population is
    itself a source of RNA-seq variation this notebook cares about, so it shouldn't be left to
    chance. 1000 Genomes families never span more than one population (verified against
    `notebooks/data/individuals.json`), so splitting per population is safe.

    The split is applied at the group level (like `split_groups`), so `test_size` is a
    fraction of *families* within each population, not of individuals -- exact counts will
    vary slightly depending on family sizes.
    """
    pop_to_class = population_to_class(class_map)

    by_population: Dict[str, List[str]] = {}
    for sample_id, info in individuals.items():
        population = info["population"]
        if population in pop_to_class:
            by_population.setdefault(population, []).append(sample_id)

    train_by_class: Dict[str, List[str]] = {class_name: [] for class_name in class_map}
    test_by_class: Dict[str, List[str]] = {class_name: [] for class_name in class_map}

    for population, sample_ids in by_population.items():
        class_name = pop_to_class[population]

        groups: Dict[str, List[str]] = {}
        for sample_id in sample_ids:
            family_id = individuals[sample_id].get("family_id") or sample_id
            groups.setdefault(str(family_id), []).append(sample_id)

        group_list = list(groups.values())
        random.Random(random_seed).shuffle(group_list)

        n_test_groups = int(round(test_size * len(group_list)))
        test_groups, train_groups = group_list[:n_test_groups], group_list[n_test_groups:]

        train_by_class[class_name].extend(sid for group in train_groups for sid in group)
        test_by_class[class_name].extend(sid for group in test_groups for sid in group)

    for class_name in class_map:
        train_by_class[class_name].sort()
        test_by_class[class_name].sort()

    return train_by_class, test_by_class
