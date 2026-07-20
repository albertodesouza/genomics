"""Programmatic literature-known variants with a signed direction of effect on RNA-seq
expression (enhance/diminish), sourced from GTEx skin-tissue eQTLs via the repo's
`gtex_database` skill (`.claude/skills/gtex_database/scripts/gtex_cli.py`) -- replaces
hand-curated, hardcoded per-gene variant lists.

GTEx's normalized effect size (NES) is a signed, quantitative measurement of a variant's
effect on a gene's expression (positive = higher expression, negative = lower) -- the
closest available match to what this notebook series already measures (AlphaGenome's
melanocyte-of-skin RNA-seq track). "Significant eQTL" filtering is already GTEx's own
threshold; no extra p-value filtering is applied here.

Follows `annotations.py`'s cache-to-`notebooks/.cache/annotations/` convention: fetching
once from any notebook warms the cache for every other notebook (`ANNOTATIONS_CACHE_DIR`
resolves to the same physical directory everywhere).
"""

from __future__ import annotations

import json
import re
import subprocess
import warnings
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

GTEX_CLI = Path(__file__).resolve().parent.parent.parent / ".claude" / "skills" / "gtex_database" / "scripts" / "gtex_cli.py"

# Skin tissues -- the closest GTEx match to AlphaGenome's melanocyte-of-skin RNA-seq track
# already used throughout this notebook series (ONTOLOGY_TERMS[0]).
SKIN_TISSUES = "Skin - Sun Exposed (Lower leg),Skin - Not Sun Exposed (Suprapubic)"

# score_variants scoring is per-variant, live-API work -- cap candidates per gene so a
# high-eQTL-count gene (e.g. TYR, 67 unique skin eQTLs) doesn't dominate the API budget.
DEFAULT_TOP_N = 10

_VARIANT_ID_RE = re.compile(r"^(chr[0-9XYM]+)_(\d+)_([ACGT]+)_([ACGT]+)_b\d+$")


def _run_gtex_cli(args: List[str], output_path: Path) -> dict:
    result = subprocess.run(
        ["uv", "run", str(GTEX_CLI), *args, "--output", str(output_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"gtex_cli.py {args} failed:\n{result.stderr}")
    return json.loads(output_path.read_text())


def get_gencode_id(gene: str, cache_dir: Path) -> str:
    """Gene symbol -> Versioned GENCODE ID, cached in a single dict file."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / "gtex_gencode_ids.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    if gene in cache:
        return cache[gene]

    tmp_path = cache_dir / f"_tmp_gencode_{gene}.json"
    data = _run_gtex_cli(["resolve-gencode-id", gene], tmp_path)
    tmp_path.unlink(missing_ok=True)
    gencode_id = data["gencode_id"]
    cache[gene] = gencode_id
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))
    return gencode_id


def fetch_gene_eqtls_raw(
    gene: str, cache_dir: Path, tissues: str = SKIN_TISSUES, force_refresh: bool = False,
) -> List[dict]:
    """Raw GTEx `singleTissueEqtl` records for a gene, cached verbatim per gene."""
    eqtl_cache_dir = cache_dir / "gtex_eqtls"
    eqtl_cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = eqtl_cache_dir / f"{gene}.json"

    if not force_refresh and cache_path.exists():
        return json.loads(cache_path.read_text())

    gencode_id = get_gencode_id(gene, cache_dir)
    data = _run_gtex_cli(
        ["get-gene-eqtls", gencode_id, "--tissues", tissues], cache_path,
    )
    return data


def _parse_variant_id(variant_id: str) -> tuple:
    m = _VARIANT_ID_RE.match(variant_id)
    if not m:
        raise ValueError(f"Unrecognized GTEx variantId format: {variant_id!r}")
    chrom, pos, ref, alt = m.groups()
    return chrom, int(pos), ref, alt


def get_literature_variants_for_window(
    gene: str,
    window_chrom: str,
    window_start: int,
    window_end: int,
    cache_dir: Path,
    tissues: str = SKIN_TISSUES,
    top_n: Optional[int] = DEFAULT_TOP_N,
    force_refresh: bool = False,
) -> pd.DataFrame:
    """Literature-known variants for `gene` with a signed direction of effect on RNA-seq
    expression, from GTEx skin-tissue eQTLs, restricted to an explicit window.

    Returns a DataFrame with columns: gene, variant ("chr19:3565601:T>C"), direction
    ("enhance" if nes>0 else "diminish"), nes, p_value, tissue, n_tissues -- sorted by
    p_value ascending (most statistically confident first), capped to `top_n` rows (pass
    `top_n=None` for the full count, e.g. when using this as a gene-ranking metric).

    Variants seen in multiple tissues are deduplicated (keeping the most-significant
    tissue's nes/p_value); a variant with a sign *conflict* across tissues is dropped
    with a warning rather than silently picking one direction.

    Candidates are filtered to `(window_chrom, window_start, window_end)`, since
    `score_variants` scores each variant against an interval centered on *the variant*,
    not the gene -- an eQTL far outside the gene's own window would silently vanish from
    downstream `tidy_scores` filtering by gene name. Callers with an already-built
    `windows_df` (the staged, active gene panel) should use `get_literature_variants`
    instead; this window-explicit form also works for candidate genes that haven't been
    staged into the pipeline yet (e.g. during gene-selection ranking), since it only needs
    a GTF-derived window, not a fully materialized dataset.
    """
    raw = fetch_gene_eqtls_raw(gene, cache_dir, tissues=tissues, force_refresh=force_refresh)

    by_variant: Dict[str, dict] = {}
    conflicts = set()
    for rec in raw:
        variant_id = rec["variantId"]
        chrom, pos, ref, alt = _parse_variant_id(variant_id)
        if chrom != window_chrom or not (window_start <= pos < window_end):
            continue

        direction = "enhance" if rec["nes"] > 0 else "diminish"
        existing = by_variant.get(variant_id)
        if existing is None:
            by_variant[variant_id] = {
                "gene": gene,
                "variant": f"{chrom}:{pos}:{ref}>{alt}",
                "direction": direction,
                "nes": rec["nes"],
                "p_value": rec["pValue"],
                "tissue": rec["tissueSiteDetailId"],
                "n_tissues": 1,
            }
        else:
            existing["n_tissues"] += 1
            if existing["direction"] != direction:
                conflicts.add(variant_id)
            elif rec["pValue"] < existing["p_value"]:
                existing["nes"] = rec["nes"]
                existing["p_value"] = rec["pValue"]
                existing["tissue"] = rec["tissueSiteDetailId"]

    if conflicts:
        warnings.warn(
            f"{gene}: dropping {len(conflicts)} variant(s) with conflicting effect "
            f"direction across tissues: {sorted(conflicts)}"
        )
    for variant_id in conflicts:
        del by_variant[variant_id]

    columns = ["gene", "variant", "direction", "nes", "p_value", "tissue", "n_tissues"]
    if not by_variant:
        return pd.DataFrame(columns=columns)

    df = pd.DataFrame(by_variant.values(), columns=columns).sort_values("p_value").reset_index(drop=True)
    if top_n is not None:
        df = df.head(top_n).reset_index(drop=True)
    return df


def get_literature_variants(
    gene: str,
    windows_df: pd.DataFrame,
    cache_dir: Path,
    tissues: str = SKIN_TISSUES,
    top_n: Optional[int] = DEFAULT_TOP_N,
    force_refresh: bool = False,
) -> pd.DataFrame:
    """`get_literature_variants_for_window` for a gene already staged in `windows_df` (the
    active gene panel) -- resolves the window from there instead of taking it explicitly."""
    gene_row = windows_df.loc[windows_df["gene"] == gene]
    if gene_row.empty:
        raise ValueError(f"{gene!r} not found in windows_df")
    gene_row = gene_row.iloc[0]
    return get_literature_variants_for_window(
        gene, gene_row["chrom"], gene_row["start"], gene_row["end"], cache_dir,
        tissues=tissues, top_n=top_n, force_refresh=force_refresh,
    )


def get_literature_variants_all_genes(
    genes: List[str],
    windows_df: pd.DataFrame,
    cache_dir: Path,
    tissues: str = SKIN_TISSUES,
    top_n: Optional[int] = DEFAULT_TOP_N,
    force_refresh: bool = False,
) -> Dict[str, pd.DataFrame]:
    """`get_literature_variants` for every gene in `genes` -- one independent per-gene
    cache file, so re-running after the first successful fetch costs zero network calls."""
    return {
        gene: get_literature_variants(
            gene, windows_df, cache_dir, tissues=tissues, top_n=top_n, force_refresh=force_refresh,
        )
        for gene in genes
    }


def literature_variant_strings(
    gene: str,
    windows_df: pd.DataFrame,
    cache_dir: Path,
    tissues: str = SKIN_TISSUES,
    top_n: Optional[int] = DEFAULT_TOP_N,
    force_refresh: bool = False,
) -> List[str]:
    """Drop-in replacement for a hardcoded `variants = [...]` list: just the variant
    strings, sorted by p_value ascending."""
    df = get_literature_variants(
        gene, windows_df, cache_dir, tissues=tissues, top_n=top_n, force_refresh=force_refresh,
    )
    return df["variant"].tolist()


def add_phenotype_direction(df: pd.DataFrame, gene_coefficient_sign: float) -> pd.DataFrame:
    """Adds a `phenotype_direction` column ("pro-pigmentation"/"anti-pigmentation") derived from
    each variant's expression-direction (the `direction` column: "enhance"/"diminish", from GTEx
    NES) combined with `gene_coefficient_sign` -- the sign of this gene's own coefficient in a
    classifier trained to predict "strong pigmentation" (the African-ancestry-population-proxy
    class used throughout this notebook series) from this gene's expression. A positive
    coefficient means higher expression of this gene is, empirically, in this study's own data,
    associated with the "strong pigmentation" class; negative means the opposite.

    "pro-pigmentation" = a variant whose effect on expression (enhance/diminish) points the same
    way this gene's expression already points toward "strong pigmentation" in the trained model
    (i.e. `direction == "enhance"` and `gene_coefficient_sign > 0`, or `direction == "diminish"`
    and `gene_coefficient_sign < 0`). Every other combination is "anti-pigmentation".

    Note on interpretation: this labels a variant's *predicted phenotype direction according to
    this study's own trained classifier*, not an independent ground truth -- there is no external
    per-gene "does more expression mean more pigmentation" table involved (deliberately: that
    would need error-prone manual special-casing, e.g. ASIP is an MC1R antagonist and would flip
    sign relative to most other genes here). Cross-checking a variant's predicted RNA-seq shift
    against this label is by construction non-circular for `score_variants`' raw_score and for
    `gene_span_model`/`nn_model` (different features/model than the one the sign convention was
    drawn from), but partially circular for `mane_exon_model`'s own logit, since that is the exact
    model `gene_coefficient_sign` came from.
    """
    df = df.copy()
    increases_pigmentation_class = (df["direction"] == "enhance") == (gene_coefficient_sign > 0)
    df["phenotype_direction"] = np.where(increases_pigmentation_class, "pro-pigmentation", "anti-pigmentation")
    return df
