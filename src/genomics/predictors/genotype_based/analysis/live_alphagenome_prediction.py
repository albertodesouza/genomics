"""Thin, cached wrapper around the live AlphaGenome API for re-predicting a (possibly edited)
haplotype sequence.

Promotes notebooks/genotype_cnn_alignment_deeplift.ipynb's ``predict_modified_sequence`` (Section
9) from a notebook-local, RNA-seq-only helper into a reusable, output-type-agnostic one, so the
same code path serves both RNA-seq re-prediction after an edit and CAGE track fetching for the
interactive viewer (neither of which the notebook needed together).

Deliberately does **not** reuse ``genomics.workflows.alphagenome.neural_module.AlphaGenomeAnalyzer``:
despite its ``predict_sequence`` method name, that class ignores the sequence argument's content
and calls ``predict_interval`` on a synthetic interval built from the sequence's length -- it
cannot actually score an edited sequence. This module talks to
``alphagenome.models.dna_client`` directly, exactly as the notebook does.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_OUTPUT_ATTR_BY_NAME = {
    "RNA_SEQ": "rna_seq",
    "CAGE": "cage",
    "ATAC": "atac",
    "DNASE": "dnase",
    "CHIP_HISTONE": "chip_histone",
    "CHIP_TF": "chip_tf",
}


def resolve_api_key(explicit: Optional[str] = None) -> str:
    """Resolve the AlphaGenome API key from an explicit argument, the environment, or ``~/.env``
    (matching the notebook's own resolution order), raising a clear error if none is found."""
    if explicit:
        return explicit
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if api_key:
        return api_key
    try:
        from dotenv import dotenv_values

        values = dotenv_values(Path.home() / ".env")
        api_key = values.get("ALPHAGENOME_API_KEY")
    except ImportError:
        api_key = None
    if not api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY not found in the environment, an explicit argument, or ~/.env.")
    return api_key


def reorder_to_canonical(
    values: np.ndarray,
    metadata_records: Sequence[Dict[str, Any]],
    canonical_order: Sequence[Tuple[str, str]],
) -> np.ndarray:
    """Reindex AlphaGenome output columns (``values``, with one column per row of
    ``metadata_records``) into ``canonical_order`` = [(ontology_curie, strand), ...], the fixed
    column order the trained model's cached on-disk predictions use.

    Promoted from the notebook's ``reorder_to_canonical``, generalized to take a plain list of
    dicts (as returned by :func:`LiveAlphaGenomePredictor.predict_sequence`) instead of requiring
    a pandas DataFrame.
    """
    lookup = {(rec["ontology_curie"], rec["strand"]): i for i, rec in enumerate(metadata_records)}
    try:
        col_idx = [lookup[key] for key in canonical_order]
    except KeyError as exc:
        raise KeyError(f"canonical track {exc} not present in predicted metadata") from exc
    return values[:, col_idx]


class LiveAlphaGenomePredictor:
    """Fetches (and disk-caches) AlphaGenome predictions for an explicit DNA sequence, generalizing
    the notebook's RNA-seq-only ``predict_modified_sequence`` to any requested output type."""

    def __init__(self, api_key: Optional[str] = None, organism: Optional[Any] = None, cache_dir: Optional[Path] = None):
        from alphagenome.models import dna_client

        self._dna_client_module = dna_client
        self.client = dna_client.create(api_key=resolve_api_key(api_key))
        self.organism = organism if organism is not None else dna_client.Organism.HOMO_SAPIENS
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        if self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_paths(self, cache_key: str) -> Tuple[Optional[Path], Optional[Path]]:
        if self.cache_dir is None:
            return None, None
        values_path = self.cache_dir / f"{cache_key}.npz"
        meta_path = self.cache_dir / f"{cache_key}_meta.json"
        return values_path, meta_path

    def predict_sequence(
        self,
        *,
        cache_key: str,
        sequence: str,
        interval: Any,
        output_type: str,
        ontology_terms: Optional[List[str]],
    ) -> Tuple[np.ndarray, List[Dict[str, Any]]]:
        """Run (or fetch from cache) an AlphaGenome ``predict_sequence`` call for ``sequence`` over
        ``interval``, restricted to ``output_type`` (e.g. "RNA_SEQ", "CAGE") and ``ontology_terms``.

        Returns ``(values, metadata_records)``, where ``metadata_records`` is a list of
        ``{"ontology_curie": ..., "strand": ...}`` dicts aligned to ``values``'s columns -- the
        same shape :func:`reorder_to_canonical` expects.
        """
        values_path, meta_path = self._cache_paths(cache_key)
        if values_path is not None and values_path.exists() and meta_path is not None and meta_path.exists():
            values = np.load(values_path)["values"]
            metadata_records = json.loads(meta_path.read_text())
            return values, metadata_records

        output_type_upper = output_type.upper()
        if output_type_upper not in _OUTPUT_ATTR_BY_NAME:
            raise ValueError(f"output_type nao suportado: {output_type!r}")
        attr_name = _OUTPUT_ATTR_BY_NAME[output_type_upper]

        output = self.client.predict_sequence(
            sequence,
            organism=self.organism,
            requested_outputs=[getattr(self._dna_client_module.OutputType, output_type_upper)],
            ontology_terms=ontology_terms,
            interval=interval,
        )
        track_output = getattr(output, attr_name)
        values = np.asarray(track_output.values)
        metadata_records = track_output.metadata[["ontology_curie", "strand"]].to_dict("records")

        if values_path is not None and meta_path is not None:
            tmp_values_path = values_path.with_name(values_path.name + ".tmp.npz")
            tmp_meta_path = meta_path.with_suffix(meta_path.suffix + ".tmp")
            np.savez_compressed(tmp_values_path, values=values)
            tmp_meta_path.write_text(json.dumps(metadata_records))
            os.replace(tmp_values_path, values_path)
            os.replace(tmp_meta_path, meta_path)

        return values, metadata_records
