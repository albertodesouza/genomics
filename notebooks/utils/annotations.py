"""Caches the small network lookups `simplified_notebook.ipynb`'s setup cells make
(gene descriptions from mygene.info, ontology labels from EBI OLS, and the GENCODE
GTF itself) under `notebooks/.cache/annotations/`, so re-running the notebook
top-to-bottom -- routine during development -- doesn't re-fetch the same responses
(or re-download the ~tens-of-MB GTF feather) every time.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import requests

GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)


def _cached_lookup(cache_dir: Path, cache_name: str, key: str, fetch) -> Optional[str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / cache_name
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    if key in cache:
        return cache[key]

    value = fetch()
    cache[key] = value
    cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))
    return value


def get_gene_description(gene: str, cache_dir: Path) -> Optional[str]:
    def fetch():
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": gene, "species": "human", "fields": "name,summary"},
            timeout=10,
        )
        r.raise_for_status()
        hits = r.json().get("hits", [])
        return (hits[0].get("summary") or hits[0].get("name")) if hits else None

    return _cached_lookup(cache_dir, "gene_descriptions.json", gene, fetch)


def get_ontology_label(term_id: str, cache_dir: Path) -> Optional[str]:
    def fetch():
        r = requests.get(
            "https://www.ebi.ac.uk/ols4/api/search",
            params={"q": term_id, "ontology": "cl", "exact": "true"},
            timeout=10,
        )
        r.raise_for_status()
        docs = r.json()["response"]["docs"]
        return docs[0]["label"] if docs else None

    return _cached_lookup(cache_dir, "ontology_labels.json", term_id, fetch)


def load_gtf(cache_dir: Path):
    """Downloads the GENCODE v46 GTF feather once to `cache_dir`, reusing it on
    every subsequent call -- this reference file is static and not worth re-fetching
    on every notebook run."""
    import pandas as pd

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / "gencode.v46.annotation.gtf.gz.feather"
    if not cache_path.exists():
        r = requests.get(GTF_URL, timeout=60)
        r.raise_for_status()
        cache_path.write_bytes(r.content)
    return pd.read_feather(cache_path)
