"""Loaders for this notebook's bundled, self-contained data (`notebooks/data/`).

Replaces what previously came from `genomics.predictors.genotype_based.config` +
`genomics.core.data_registry` (the pigmentation YAML config, the shared dataset
registry) so the notebook has no import-time or read-time dependency on the rest of
the `genomics` package or on the external dataset mount.
"""

from __future__ import annotations

import json
import random
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

NOTEBOOK_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = NOTEBOOK_DIR / "data"


def load_experiment() -> dict:
    """Returns {"genes": [...], "ontology_terms": [...], "class_map": {...}}."""
    with open(DATA_DIR / "experiment.json") as f:
        return json.load(f)


def load_genes() -> Dict[str, dict]:
    """Returns {gene: {"chrom": ..., "start": ..., "end": ...}} (0-based half-open)."""
    with open(DATA_DIR / "genes.json") as f:
        return json.load(f)


def load_individuals() -> Dict[str, dict]:
    """Returns {sample_id: {"population", "superpopulation", "sex", "family_id"}},
    restricted to individuals in the class_map. Sourced from the 1000 Genomes pedigree
    CSV -- see `notebooks/utils/prepare_data.py`."""
    with open(DATA_DIR / "individuals.json") as f:
        return json.load(f)


def population_to_class(class_map: Dict[str, List[str]]) -> Dict[str, str]:
    return {population: class_name for class_name, populations in class_map.items() for population in populations}


def samples_by_class(individuals: Dict[str, dict], class_map: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """Returns {class_name: [sample_id, ...]} (sorted, deterministic order)."""
    pop_to_class = population_to_class(class_map)
    grouped: Dict[str, List[str]] = {class_name: [] for class_name in class_map}
    for sample_id, info in individuals.items():
        class_name = pop_to_class.get(info["population"])
        if class_name is not None:
            grouped[class_name].append(sample_id)
    for sample_ids in grouped.values():
        sample_ids.sort()
    return grouped


def gene_paths(gene: str) -> Tuple[Path, Path]:
    """Returns (ref_window_fa, gene_vcf_path) for a gene."""
    return DATA_DIR / "references" / gene / "ref.window.fa", DATA_DIR / "variants" / f"{gene}.vcf.gz"


def sample_background_variants_from_vcf(
    vcf_path: Path,
    exclude_positions: set,
    min_allele_frequency: float,
    n: int,
    random_seed: int = 13,
) -> list:
    """Real 1000 Genomes background variants -- not synthetic substitutions: every biallelic
    SNV in this gene's sliced 1000 Genomes VCF (`notebooks/data/variants/<gene>.vcf.gz`, the
    same file this notebook series already uses to build individual consensus sequences) with
    at least `min_allele_frequency` in the panel, excluding the curated literature positions,
    then a reproducible random sample of `n` of them. Shared across notebooks so they compare
    against the exact same background set (same seed/threshold)."""
    from alphagenome.data import genome

    if shutil.which("bcftools") is None:
        raise RuntimeError(
            "bcftools not found on PATH. Run scripts/env/start_genomics_universal.sh "
            "or otherwise install bcftools."
        )

    view = subprocess.Popen(
        ["bcftools", "view", "-v", "snps", "-m2", "-M2", "-i", f"INFO/AF>={min_allele_frequency}", str(vcf_path)],
        stdout=subprocess.PIPE,
    )
    query = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%AF\n"],
        stdin=view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    view.stdout.close()
    view.wait()
    if query.returncode != 0:
        raise RuntimeError(f"bcftools query failed: {query.stderr}")

    candidates = []
    for line in query.stdout.splitlines():
        chrom, pos, ref, alt, af = line.split("\t")
        pos = int(pos)
        if pos in exclude_positions:
            continue
        candidates.append((chrom, pos, ref, alt))

    rng = random.Random(random_seed)
    sampled = rng.sample(candidates, min(n, len(candidates)))
    return [
        genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
        for chrom, pos, ref, alt in sampled
    ]
