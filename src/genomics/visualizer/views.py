"""Cohort -> ``.view.json`` export (the former standalone View Builder app)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from genomics.predictors.genotype_based.apps.view_builder import (
    DEFAULT_DOWNSAMPLE_FACTOR,
    DEFAULT_NORMALIZATION_METHOD,
    DEFAULT_WINDOW_CENTER_SIZE,
    NORMALIZATION_METHODS,
    _dataset_input_yaml_snippet,
    _default_output_path,
    _safe_view_name,
)
from genomics.visualizer.datasets import Dataset


def _csv(value: Any) -> List[str]:
    if value is None:
        return []
    parts = value if isinstance(value, list) else str(value).split(",")
    return [str(p).strip() for p in parts if str(p).strip()]


class ViewService:
    def options(self) -> Dict[str, Any]:
        return {
            "normalization_methods": NORMALIZATION_METHODS,
            "defaults": {
                "window_center_size": DEFAULT_WINDOW_CENTER_SIZE,
                "downsample_factor": DEFAULT_DOWNSAMPLE_FACTOR,
                "normalization_method": DEFAULT_NORMALIZATION_METHOD,
            },
        }

    def preview(self, dataset: Dataset, payload: Dict[str, Any]) -> Dict[str, Any]:
        name = _safe_view_name(str(payload.get("name") or "new_view"))
        genes = _csv(payload.get("genes")) or list(dataset.genes)
        unknown = sorted(set(genes) - set(dataset.genes))
        if unknown:
            raise ValueError(f"Unknown genes: {', '.join(unknown[:10])}")
        samples = _csv(payload.get("sample_ids"))
        unknown_samples = [s for s in samples if s not in dataset.sample_index]
        if unknown_samples:
            raise ValueError(f"Unknown samples: {', '.join(unknown_samples[:10])}")
        outputs = _csv(payload.get("alphagenome_outputs")) or ["rna_seq"]
        window = int(payload.get("window_center_size") or DEFAULT_WINDOW_CENTER_SIZE)
        downsample = int(payload.get("downsample_factor") or DEFAULT_DOWNSAMPLE_FACTOR)
        normalization = str(payload.get("normalization_method") or DEFAULT_NORMALIZATION_METHOD)
        if window < 1 or downsample < 1:
            raise ValueError("window_center_size and downsample_factor must be >= 1")
        if normalization not in NORMALIZATION_METHODS:
            raise ValueError(f"normalization_method must be one of {', '.join(NORMALIZATION_METHODS)}")
        view = {
            "name": name,
            "description": str(payload.get("description") or ""),
            "dataset_dir": str(dataset.path),
            "alphagenome_outputs": outputs,
            "haplotype_mode": "H1+H2",
            "tensor_layout": "haplotype_channels",
            "window_center_size": window,
            "downsample_factor": downsample,
            "genes_to_use": genes,
            "sample_ids": samples,
            "sample_ids_path": None,
            "superpopulations_to_use": None,
            "populations_to_use": None,
            "normalization_method": normalization,
            "ontology_terms": _csv(payload.get("ontology_terms")) or None,
        }
        output_path = str(payload.get("output_path") or _default_output_path(name))
        return {
            "view": view,
            "output_path": output_path,
            "yaml_snippet": _dataset_input_yaml_snippet(view),
            "counts": {"genes": len(genes), "samples": len(samples), "outputs": len(outputs)},
        }

    def save(self, dataset: Dataset, payload: Dict[str, Any]) -> Dict[str, Any]:
        preview = self.preview(dataset, payload)
        path = Path(preview["output_path"]).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
        path = path.resolve()
        if path.exists() and not payload.get("overwrite"):
            raise FileExistsError(f"File already exists: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(preview["view"], f, indent=2)
            f.write("\n")
        preview["saved_path"] = str(path)
        return preview
