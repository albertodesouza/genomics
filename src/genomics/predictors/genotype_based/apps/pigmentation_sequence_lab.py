"""Interactive web app turning notebooks/genotype_cnn_alignment_deeplift.ipynb's analysis into a
browser tool: pick an individual or a pigmentation-class average, see RNA-seq/CAGE tracks over the
full 524,288 bp AlphaGenome window with zoom/pan, see a bcftools_chain-aligned base-frequency
"sequence logo" (or, for a single individual, their actual bases with ``X`` for gaps), and
interactively overwrite/scramble a selected region of an individual's sequence and re-run
AlphaGenome + the trained CNN2 pigmentation classifier to see the effect.

Follows this repo's established convention (see ``alphagenome_track_viewer.py``,
``aligned_dna_viewer.py``): a hand-rolled ``http.server`` app serving inline HTML/vanilla JS, no
Flask/FastAPI. Kept as a *sibling* of ``alphagenome_track_viewer.py`` rather than a modification of
it, since this app additionally needs a loaded CNN2 checkpoint (torch) and a live AlphaGenome API
key for editing -- neither of which the read-only, dependency-light track viewer should require.

Reuses, rather than reimplements: ``AlphaGenomeRepository`` (RNA-seq track/heatmap serving,
downsampling, pigmentation class helpers), ``DynamicIndelAligner``/``BcftoolsChainMapper``
(bcftools_chain alignment), and the promoted notebook logic in
``genomics.predictors.genotype_based.analysis`` (edit geometry, live AlphaGenome re-prediction,
alignment frequency, CNN2 inference context).
"""
from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse

import numpy as np

from genomics.predictors.genotype_based.alignment.bcftools_chain_mapper import BcftoolsChainMapper
from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner
from genomics.predictors.genotype_based.analysis.alignment_frequency import (
    class_base_frequency,
    individual_aligned_bases_with_source,
)
from genomics.predictors.genotype_based.analysis.haplotype_edit_geometry import (
    EditResult,
    apply_overwrite,
    apply_scramble_region,
    haplotype_indel_events,
    haplotype_local_idx,
)
from genomics.predictors.genotype_based.analysis.live_alphagenome_prediction import (
    LiveAlphaGenomePredictor,
    reorder_to_canonical,
    resolve_api_key,
)
from genomics.predictors.genotype_based.analysis.pigmentation_model_context import PigmentationModelContext
from genomics.predictors.genotype_based.apps.alphagenome_track_viewer import (
    DEFAULT_POINTS,
    MAX_POINTS,
    AlphaGenomeRepository,
    _align_signal_window,
    _downsample,
    _extract_full_signal,
    _json_scalar,
    _load_prediction_array,
    _safe_int,
    _sanitize_json,
)
from genomics.workspace import DEFAULT_CONSENSUS_DATASET_DIR

DEFAULT_PORT = 8781
# Building an individual's bcftools_chain entry (first time only, then disk-cached) shells out to
# `bcftools consensus` -- an I/O-bound subprocess call that releases the GIL, so a class-average
# over hundreds of never-before-seen individuals is worth parallelizing across a thread pool
# rather than fetching each one serially.
CLASS_AVERAGE_MAX_WORKERS = 12
LETTER_LOD_THRESHOLD = 400  # max span (bp) rendered as literal letters instead of density bars
CAGE_ONTOLOGY_TERMS = ["CL:0002566", "CL:0002567"]  # dark / light melanocyte, matches the notebook
CLASS_ALIASES = {
    "strong": "strong pigmentation",
    "weak": "weak pigmentation",
    "strong pigmentation": "strong pigmentation",
    "weak pigmentation": "weak pigmentation",
}


def _resolve_class_value(subject: str) -> str:
    resolved = CLASS_ALIASES.get(subject.strip().lower())
    if resolved is None:
        raise ValueError(f"subject de classe desconhecido: {subject!r} (use 'strong' ou 'weak')")
    return resolved


def _read_fasta_sequence(path: Path) -> str:
    lines: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines).upper()


def _gene_interval(dataset_dir: Path, gene: str):
    from alphagenome.data import genome

    window_meta = json.loads((dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text())
    chrom = window_meta["chromosome"]
    start_1based = int(window_meta["start"])
    full_length = int(window_meta["end"]) - start_1based + 1
    return genome.Interval(chromosome=chrom, start=start_1based - 1, end=start_1based - 1 + full_length)


def _haplotype_fixed_fasta_path(dataset_dir: Path, sample_id: str, gene: str, haplotype: str) -> Path:
    return dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.{haplotype}.window.fixed.fa"


def _cage_prediction_paths(dataset_dir: Path, sample_id: str, gene: str, haplotype: str) -> Tuple[Path, Path]:
    pred_dir = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
    return pred_dir / "cage.npz", pred_dir / "cage_metadata.json"


def _downsample_frequency(freq: Dict[str, Any], start: int, points: int) -> Tuple[List[int], Dict[str, List[float]]]:
    """Bucket-mean downsample per symbol, mirroring alphagenome_track_viewer's ``_downsample`` but
    applied independently to each of the A/C/G/T/gap arrays in a frequency payload."""
    symbols = ("A", "C", "G", "T", "gap")
    length = len(freq["A"]) if len(freq["A"]) else 0
    if length == 0:
        return [], {sym: [] for sym in symbols}
    if length <= points:
        xs = list(range(start, start + length))
        return xs, {sym: [float(v) for v in freq[sym]] for sym in symbols}

    bucket_count = max(1, points)
    edges = np.linspace(0, length, num=bucket_count + 1, dtype=np.int64)
    xs: List[int] = []
    out: Dict[str, List[float]] = {sym: [] for sym in symbols}
    for i in range(bucket_count):
        left, right = int(edges[i]), int(edges[i + 1])
        if right <= left:
            continue
        xs.append(start + (left + right - 1) // 2)
        for sym in symbols:
            out[sym].append(float(np.mean(freq[sym][left:right])))
    return xs, out


def _stats_for(signal: np.ndarray) -> Dict[str, Any]:
    finite = signal[np.isfinite(signal)]
    return {
        "min": None if finite.size == 0 else _json_scalar(np.min(finite)),
        "max": None if finite.size == 0 else _json_scalar(np.max(finite)),
        "mean": None if finite.size == 0 else _json_scalar(np.mean(finite)),
    }


def _aligned_track_series(
    mapper: BcftoolsChainMapper,
    expanded_length: int,
    dataset_dir: Path,
    sample_id: str,
    gene: str,
    haplotype: str,
    output: str,
    track: int,
    start: int,
    length: int,
    points: int,
    label: Optional[str] = None,
) -> Dict[str, Any]:
    """One sample/haplotype's track, scattered onto the gene's *full* (524,288bp) bcftools_chain
    expanded axis at the caller's requested [start, start+length-1] window -- unlike
    ``AlphaGenomeRepository.compare_payload``, which (via its own hardcoded
    ``DEFAULT_VIEW_LENGTH``-centered aligner) always returns the same fixed 32,768bp CNN-input
    window regardless of the requested start/length, this supports genuinely panning/zooming
    across the whole AlphaGenome window while staying bcftools_chain-aligned."""
    pred_dir = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
    array, _array_key, _npz_keys = _load_prediction_array(pred_dir / f"{output}.npz")
    full_signal = _extract_full_signal(array, track=track)
    entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
    if entry is None:
        raise RuntimeError(f"Alignment entry ausente para {sample_id} {haplotype}")
    signal = _align_signal_window(full_signal, entry, expanded_length, start, length)
    x, y, method = _downsample(signal, start=start, points=points)
    return {
        "label": label or f"{sample_id}_{haplotype}",
        "sample": sample_id,
        "haplotype": haplotype,
        "x": x,
        "y": y,
        "returned_points": len(y),
        "downsample": method,
        "stats": _stats_for(signal),
    }


def _aligned_signal_for_sample(
    mapper: BcftoolsChainMapper,
    expanded_length: int,
    dataset_dir: Path,
    sample_id: str,
    gene: str,
    haplotype: str,
    output: str,
    track: int,
    start: int,
    length: int,
) -> Optional[np.ndarray]:
    pred_dir = dataset_dir / "individuals" / sample_id / "windows" / gene / f"predictions_{haplotype}"
    array, _array_key, _npz_keys = _load_prediction_array(pred_dir / f"{output}.npz")
    full_signal = _extract_full_signal(array, track=track)
    entry = mapper.get_haplotype_entry(gene, sample_id, haplotype)
    if entry is None:
        return None
    return _align_signal_window(full_signal, entry, expanded_length, start, length)


def _aligned_class_mean_track_series(
    mapper: BcftoolsChainMapper,
    expanded_length: int,
    dataset_dir: Path,
    sample_ids: List[str],
    gene: str,
    haplotype: str,
    output: str,
    track: int,
    start: int,
    length: int,
    points: int,
    label: str,
) -> Dict[str, Any]:
    """True average of every sample's aligned track over [start, start+length-1] -- the
    bcftools_chain-aligned, arbitrary-window analog of
    ``AlphaGenomeRepository.heatmap_payload(mode="class_mean", track_count=1)``.

    Fetches each sample's aligned signal concurrently (see ``CLASS_AVERAGE_MAX_WORKERS``): most
    individuals' bcftools_chain entries are already disk-cached from prior use, but any that
    aren't require a real, uncached ``bcftools consensus`` subprocess call, and a class can span
    hundreds of individuals.
    """
    aligned_signals = []
    missing = 0
    with ThreadPoolExecutor(max_workers=min(CLASS_AVERAGE_MAX_WORKERS, len(sample_ids)) or 1) as pool:
        futures = {
            pool.submit(
                _aligned_signal_for_sample, mapper, expanded_length, dataset_dir, sample_id, gene, haplotype, output, track, start, length
            ): sample_id
            for sample_id in sample_ids
        }
        for future in as_completed(futures):
            try:
                signal = future.result()
            except Exception:
                missing += 1
                continue
            if signal is None:
                missing += 1
                continue
            aligned_signals.append(signal)
    if not aligned_signals:
        raise RuntimeError(f"Nenhum sinal alinhado disponivel para a media da classe (gene={gene}, haplotype={haplotype})")

    stack = np.stack(aligned_signals, axis=0)
    with np.errstate(invalid="ignore"):
        finite_count = np.sum(np.isfinite(stack), axis=0)
        sums = np.nansum(stack, axis=0)
        mean_signal = np.divide(
            sums, finite_count, out=np.full(stack.shape[1], np.nan, dtype=np.float32), where=finite_count > 0
        )
    x, y, method = _downsample(mean_signal, start=start, points=points)
    return {
        "label": label,
        "sample": label,
        "haplotype": haplotype,
        "x": x,
        "y": y,
        "returned_points": len(y),
        "downsample": method,
        "stats": _stats_for(mean_signal),
        "samples_used": len(aligned_signals),
        "samples_requested": len(sample_ids),
        "missing_count": missing,
    }


class PigmentationSequenceLabRepository:
    def __init__(
        self,
        dataset_dir: Path,
        config_path: Path,
        checkpoint: str = "best_accuracy",
        consensus_dataset_dir: Optional[Path] = None,
        api_key: Optional[str] = None,
    ):
        self.dataset_dir = Path(dataset_dir).resolve()
        self.consensus_dataset_dir = (
            Path(consensus_dataset_dir).resolve() if consensus_dataset_dir else DEFAULT_CONSENSUS_DATASET_DIR
        )
        self.ag_repo = AlphaGenomeRepository(self.dataset_dir, consensus_dataset_dir=self.consensus_dataset_dir)
        self.model_context = PigmentationModelContext(config_path, checkpoint=checkpoint)
        self._explicit_api_key = api_key
        self._live_predictor: Optional[LiveAlphaGenomePredictor] = None
        self._full_aligner: Optional[DynamicIndelAligner] = None
        self._full_chain_mapper: Optional[BcftoolsChainMapper] = None
        self._representative_sample_cache: Dict[str, str] = {}
        # DynamicIndelAligner.build_alignment_axis_for_gene() unconditionally rewrites a manifest
        # file and re-reads/re-parses the gene's persistent axis.json (multi-MB across a large
        # cohort) even when the axis is already finalized in this aligner's own in-memory cache --
        # so call it at most once per gene per process, here, instead of before every request.
        self._full_axis_built_genes: set = set()

        ag_options = self.ag_repo.options_payload()
        self._individuals_with_class = [
            {**row, "pigmentation_class": self.ag_repo._class_label_for_row(row, "pigmentation")}
            for row in ag_options["individuals"]
        ]
        self._defaults = {
            "subject_type": "individual",
            "subject": self.ag_repo.individuals[0] if self.ag_repo.individuals else None,
            "gene": self.model_context.genes[0] if self.model_context.genes else None,
            "haplotype": "H1",
            "output": "rna_seq",
            "start": ag_options["defaults"]["start"],
            "length": ag_options["defaults"]["length"],
            "points": DEFAULT_POINTS,
        }

    # -- shared collaborators -------------------------------------------------
    def api_key_available(self) -> bool:
        try:
            resolve_api_key(self._explicit_api_key)
            return True
        except RuntimeError:
            return False

    @property
    def defaults(self) -> Dict[str, Any]:
        return self._defaults

    def _get_live_predictor(self) -> LiveAlphaGenomePredictor:
        if self._live_predictor is None:
            cache_dir = self.dataset_dir / ".pigmentation_sequence_lab_cache" / "live_predictions"
            self._live_predictor = LiveAlphaGenomePredictor(api_key=self._explicit_api_key, cache_dir=cache_dir)
        return self._live_predictor

    def _get_full_aligner(self) -> DynamicIndelAligner:
        if self._full_aligner is None:
            self._full_aligner = DynamicIndelAligner(
                self.dataset_dir, selected_sample_ids=self.ag_repo.individuals, center_window_size=None
            )
        return self._full_aligner

    def _ensure_full_axis_built(self, gene: str) -> DynamicIndelAligner:
        aligner = self._get_full_aligner()
        if gene not in self._full_axis_built_genes:
            aligner.build_alignment_axis_for_gene(gene, self.ag_repo.individuals)
            self._full_axis_built_genes.add(gene)
        return aligner

    def _get_full_chain_mapper(self) -> BcftoolsChainMapper:
        if self._full_chain_mapper is None:
            self._full_chain_mapper = BcftoolsChainMapper(
                dataset_dir=self.dataset_dir,
                consensus_dataset_dir=self.consensus_dataset_dir,
                aligner=self._get_full_aligner(),
            )
        return self._full_chain_mapper

    def representative_sample_for_class(self, class_value: str) -> str:
        if class_value not in self._representative_sample_cache:
            samples = self.ag_repo._samples_for_class("pigmentation", class_value)
            if not samples:
                raise ValueError(f"Nenhum individuo encontrado para a classe {class_value!r}")
            self._representative_sample_cache[class_value] = samples[0]
        return self._representative_sample_cache[class_value]

    # -- options -----------------------------------------------------------
    def options_payload(self) -> Dict[str, Any]:
        return {
            "dataset_dir": str(self.dataset_dir),
            "genes": self.model_context.genes,
            "individuals": self._individuals_with_class,
            "pigmentation_classes": ["strong pigmentation", "weak pigmentation"],
            "rna_seq_ontology_terms": self.model_context.rna_seq_ontology_terms,
            "cage_ontology_terms": CAGE_ONTOLOGY_TERMS,
            "window_center_size": self.model_context.window_center_size,
            "haplotypes": ["H1", "H2"],
            "letter_lod_threshold": LETTER_LOD_THRESHOLD,
            "api_key_available": self.api_key_available(),
            "defaults": self._defaults,
        }

    def model_window_payload(self, gene: str) -> Dict[str, Any]:
        """The CNN2's centered window_center_size-bp region, expressed in the *full* 524,288bp
        bcftools_chain expanded-axis coordinate system used everywhere else in this app.

        Deliberately does not reuse ``AlphaGenomeRepository.model_window_payload``: that method's
        aligner is hardcoded to a 32,768bp-centered axis (``alphagenome_track_viewer.DEFAULT_VIEW_LENGTH``),
        so its "expanded_start"/"expanded_end" are relative to an axis that already starts near the
        middle of the full window -- a different, incompatible coordinate system from the one this
        app's tracks/base-frequency endpoints use for full-window pan/zoom.
        """
        aligner = self._ensure_full_axis_built(gene)
        center_slice = aligner.get_reference_centered_expanded_slice(gene, self.model_context.window_center_size)
        return {
            "gene": gene,
            "start": int(center_slice["expanded_start"]) + 1,
            "length": int(center_slice["expanded_end"]) - int(center_slice["expanded_start"]),
            "x_axis": "bcftools_chain_expanded_alignment",
            "model_window": center_slice,
        }

    def track_options_payload(self, subject_type: str, subject: str, gene: str, haplotype: str, output: str) -> Dict[str, Any]:
        if subject_type == "class":
            sample = self.representative_sample_for_class(_resolve_class_value(subject))
        else:
            sample = subject
        if output == "cage":
            self.ensure_cage_cached(sample, gene, haplotype)
        return self.ag_repo.track_options_payload(sample, gene, haplotype, output)

    # -- CAGE fetch-or-cache-hit into the canonical per-individual layout ------
    def ensure_cage_cached(self, sample_id: str, gene: str, haplotype: str) -> None:
        values_path, meta_path = _cage_prediction_paths(self.dataset_dir, sample_id, gene, haplotype)
        if values_path.exists() and meta_path.exists():
            return

        predictor = self._get_live_predictor()
        fasta_path = _haplotype_fixed_fasta_path(self.dataset_dir, sample_id, gene, haplotype)
        sequence = _read_fasta_sequence(fasta_path)
        interval = _gene_interval(self.dataset_dir, gene)
        values, metadata_records = predictor.predict_sequence(
            cache_key=f"{sample_id}_{gene}_{haplotype}_cage_canonical",
            sequence=sequence,
            interval=interval,
            output_type="CAGE",
            ontology_terms=CAGE_ONTOLOGY_TERMS,
        )
        values_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_values_path = values_path.with_name(values_path.name + ".tmp.npz")
        tmp_meta_path = meta_path.with_suffix(meta_path.suffix + ".tmp")
        np.savez_compressed(tmp_values_path, values=values)
        tmp_meta_path.write_text(json.dumps({"metadata": metadata_records}))
        os.replace(tmp_values_path, values_path)
        os.replace(tmp_meta_path, meta_path)

    # -- tracks (RNA-seq true average / CAGE representative-individual) -------
    def tracks_payload(
        self,
        subject_type: str,
        subject: str,
        gene: str,
        haplotypes: List[str],
        output: str,
        track: int,
        start: int,
        length: int,
        points: int,
    ) -> Dict[str, Any]:
        if output not in {"rna_seq", "cage"}:
            raise ValueError(f"output deve ser 'rna_seq' ou 'cage', recebeu {output!r}")
        haplotypes = [h for h in haplotypes if h in {"H1", "H2"}] or ["H1"]

        aligner = self._ensure_full_axis_built(gene)
        expanded_length = aligner.get_expanded_length(gene)
        mapper = self._get_full_chain_mapper()

        cage_representative_sample = None
        cage_is_true_average: Optional[bool] = None
        series: List[Dict[str, Any]] = []

        if subject_type == "individual":
            if output == "cage":
                for haplotype in haplotypes:
                    self.ensure_cage_cached(subject, gene, haplotype)
            for haplotype in haplotypes:
                series.append(_aligned_track_series(
                    mapper, expanded_length, self.dataset_dir, subject, gene, haplotype, output, track, start, length, points,
                ))
        elif subject_type == "class":
            class_value = _resolve_class_value(subject)
            if output == "cage":
                cage_is_true_average = False
                cage_representative_sample = self.representative_sample_for_class(class_value)
                for haplotype in haplotypes:
                    self.ensure_cage_cached(cage_representative_sample, gene, haplotype)
                    series.append(_aligned_track_series(
                        mapper, expanded_length, self.dataset_dir, cage_representative_sample, gene, haplotype,
                        output, track, start, length, points, label=f"{class_value}_{haplotype}_representative",
                    ))
            else:
                samples = self.ag_repo._samples_for_class("pigmentation", class_value)
                if not samples:
                    raise ValueError(f"Nenhum individuo encontrado para a classe {class_value!r}")
                for haplotype in haplotypes:
                    series.append(_aligned_class_mean_track_series(
                        mapper, expanded_length, self.dataset_dir, samples, gene, haplotype,
                        output, track, start, length, points, label=f"{class_value}_{haplotype}_mean",
                    ))
        else:
            raise ValueError(f"subject_type deve ser 'individual' ou 'class', recebeu {subject_type!r}")

        return {
            "gene": gene,
            "output": output,
            "track": track,
            "model_window": self.model_window_payload(gene)["model_window"],
            "x_axis": "bcftools_chain_expanded_alignment",
            "series": series,
            "subject_type": subject_type,
            "subject": subject,
            "cage_representative_sample": cage_representative_sample,
            "cage_is_true_average": cage_is_true_average,
        }

    # -- base-frequency / letters plot -----------------------------------------
    def base_frequency_payload(
        self,
        subject_type: str,
        subject: str,
        gene: str,
        haplotypes: List[str],
        start: int,
        length: int,
        points: int,
    ) -> Dict[str, Any]:
        haplotypes = [h for h in haplotypes if h in {"H1", "H2"}] or ["H1"]
        end = start + max(length, 1) - 1
        self._ensure_full_axis_built(gene)
        mapper = self._get_full_chain_mapper()
        lod = "letters" if length <= LETTER_LOD_THRESHOLD else "density"

        if subject_type == "individual":
            haplotype = haplotypes[0]
            bases, source_indices = individual_aligned_bases_with_source(mapper, gene, subject, haplotype, start, end)
            if lod == "letters":
                payload: Dict[str, Any] = {
                    "mode": "individual",
                    "haplotype": haplotype,
                    "bases": bases,
                    "source_indices": source_indices,
                    "start": start,
                    "end": end,
                }
            else:
                freq = {sym: np.zeros(len(bases)) for sym in ("A", "C", "G", "T", "gap")}
                for idx, base in enumerate(bases):
                    freq[base if base in "ACGT" else "gap"][idx] = 1.0
                xs, downsampled = _downsample_frequency(freq, start, points)
                payload = {
                    "mode": "individual",
                    "haplotype": haplotype,
                    "x": xs,
                    "frequencies": downsampled,
                    "n_used": 1,
                    "n_requested": 1,
                    "start": start,
                    "end": end,
                }
        elif subject_type == "class":
            class_value = _resolve_class_value(subject)
            samples = self.ag_repo._samples_for_class("pigmentation", class_value)
            if not samples:
                raise ValueError(f"Nenhum individuo encontrado para a classe {class_value!r}")
            freq = class_base_frequency(mapper, gene, samples, haplotypes, start, end)
            n_used, n_requested = int(freq.pop("n_used")), int(freq.pop("n_requested"))
            if lod == "letters":
                payload = {
                    "mode": "frequency",
                    "frequencies": {sym: [float(v) for v in arr] for sym, arr in freq.items()},
                    "n_used": n_used,
                    "n_requested": n_requested,
                    "start": start,
                    "end": end,
                }
            else:
                xs, downsampled = _downsample_frequency(freq, start, points)
                payload = {
                    "mode": "frequency",
                    "x": xs,
                    "frequencies": downsampled,
                    "n_used": n_used,
                    "n_requested": n_requested,
                    "start": start,
                    "end": end,
                }
        else:
            raise ValueError(f"subject_type deve ser 'individual' ou 'class', recebeu {subject_type!r}")

        payload["subject_type"] = subject_type
        payload["subject"] = subject
        payload["gene"] = gene
        payload["lod"] = lod
        payload["letter_lod_threshold"] = LETTER_LOD_THRESHOLD
        return payload

    # -- editing ----------------------------------------------------------------
    def _load_haplotype_sequence(self, sample_id: str, gene: str, haplotype: str) -> str:
        return _read_fasta_sequence(_haplotype_fixed_fasta_path(self.dataset_dir, sample_id, gene, haplotype))

    def edit_preview_payload(
        self, sample_id: str, gene: str, haplotype: str, start: int, end: int, op: str, base: Optional[str], seed: int
    ) -> Dict[str, Any]:
        sequence = self._load_haplotype_sequence(sample_id, gene, haplotype)
        edit = self._apply_edit(sequence, start, end, op, base, seed)
        return {
            "sample_id": sample_id, "gene": gene, "haplotype": haplotype,
            "start": edit.start, "end": edit.end, "op": op,
            "original_segment": edit.original_segment, "edited_segment": edit.edited_segment,
        }

    @staticmethod
    def _apply_edit(sequence: str, start: int, end: int, op: str, base: Optional[str], seed: int) -> EditResult:
        if op == "overwrite":
            if not base:
                raise ValueError("op='overwrite' requer 'base'")
            return apply_overwrite(sequence, start, end, base.upper())
        if op == "scramble":
            return apply_scramble_region(sequence, start, end, seed=seed)
        raise ValueError(f"op deve ser 'overwrite' ou 'scramble', recebeu {op!r}")

    def _model_window_local_range(self, sample_id: str, gene: str, haplotype: str) -> Optional[Tuple[int, int]]:
        """The CNN2's centered window_center_size-bp region, translated into this haplotype's own
        raw FASTA-local [start, end) range -- used to tell whether an edit (given in that same
        local coordinate space) can possibly change the CNN's prediction at all.

        Computed directly from the reference window's own center (mirroring
        ``DynamicIndelAligner.get_reference_centered_expanded_slice``'s ``ref_length // 2``
        centering, but expressed in plain reference coordinates) and this haplotype's own indel
        events (``haplotype_edit_geometry.haplotype_indel_events``/``haplotype_local_idx``) --
        deliberately avoids the bcftools_chain aligner/mapper entirely, since their "expanded axis"
        coordinate system depends on which ``center_window_size`` built the axis (see
        ``model_window_payload``'s docstring) and is not needed here: this only ever needs a single
        haplotype's own indel drift, which ``haplotype_local_idx`` already gives directly.
        """
        window_meta = json.loads((self.dataset_dir / "references" / "windows" / gene / "window_metadata.json").read_text())
        start_1based = int(window_meta["start"])
        ref_length = int(window_meta["end"]) - start_1based + 1
        window_size = min(self.model_context.window_center_size, ref_length)
        center_idx = ref_length // 2
        half = window_size // 2
        lo = max(0, center_idx - half)
        hi = min(ref_length, lo + window_size)
        if hi <= lo:
            return None

        vcf_path = self.dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.window.vcf.gz"
        if not vcf_path.exists():
            return None
        events = haplotype_indel_events(vcf_path, haplotype)
        local_start = haplotype_local_idx(events, start_1based, start_1based + lo)
        local_end = haplotype_local_idx(events, start_1based, start_1based + hi - 1) + 1
        return local_start, local_end

    def edit_apply_payload(
        self,
        sample_id: str,
        gene: str,
        haplotype: str,
        start: int,
        end: int,
        op: str,
        base: Optional[str],
        seed: int,
        apply_to_prediction: List[str],
    ) -> Dict[str, Any]:
        apply_to_prediction = [o for o in apply_to_prediction if o in {"rna_seq", "cage"}] or ["rna_seq"]

        sequence = self._load_haplotype_sequence(sample_id, gene, haplotype)
        edit = self._apply_edit(sequence, start, end, op, base, seed)

        local_range = self._model_window_local_range(sample_id, gene, haplotype)
        affects_prediction = local_range is not None and not (edit.end <= local_range[0] or edit.start >= local_range[1])

        baseline_tensor = self.model_context.baseline_tensor(sample_id)
        baseline_prob = self.model_context.strong_pigmentation_probability(baseline_tensor)

        # Only RNA-seq feeds the CNN2's input tensor (config.dataset_input.alphagenome_outputs ==
        # ["rna_seq"]); CAGE is re-predicted here purely so the CAGE track panel can reflect the
        # edit too, it never enters `overrides`/the model's prediction.
        overrides: Dict[Tuple[str, str], Tuple[np.ndarray, list]] = {}
        edited_tracks: Dict[str, Any] = {}
        if affects_prediction:
            predictor = self._get_live_predictor()
            interval = _gene_interval(self.dataset_dir, gene)
            _baseline_array, baseline_meta = self.model_context.load_raw_prediction(sample_id, gene, haplotype, "rna_seq")
            canonical_order = [(m["ontology_curie"], m["strand"]) for m in baseline_meta]

            for output_name in apply_to_prediction:
                ontology_terms = (
                    self.model_context.rna_seq_ontology_terms if output_name == "rna_seq" else CAGE_ONTOLOGY_TERMS
                )
                cache_key = f"{sample_id}_{gene}_{haplotype}_{op}_{start}_{end}_{output_name}"
                values, metadata_records = predictor.predict_sequence(
                    cache_key=cache_key, sequence=edit.modified_sequence, interval=interval,
                    output_type=output_name.upper(), ontology_terms=ontology_terms,
                )
                if output_name == "rna_seq":
                    reordered = reorder_to_canonical(values, metadata_records, canonical_order)
                    overrides[(gene, haplotype)] = (reordered, baseline_meta)
                edited_tracks[output_name] = {"values_shape": list(values.shape), "metadata": metadata_records}

        edited_tensor = (
            self.model_context.build_individual_tensor(sample_id, overrides=overrides) if overrides else baseline_tensor
        )
        edited_prob = self.model_context.strong_pigmentation_probability(edited_tensor)

        return {
            "sample_id": sample_id, "gene": gene, "haplotype": haplotype,
            "start": edit.start, "end": edit.end, "op": op,
            "original_segment": edit.original_segment, "edited_segment": edit.edited_segment,
            "affects_prediction": affects_prediction,
            "baseline_strong_pigmentation_probability": baseline_prob,
            "edited_strong_pigmentation_probability": edited_prob,
            "delta_strong_pigmentation_probability": edited_prob - baseline_prob,
            "edited_tracks": edited_tracks,
        }


class PigmentationSequenceLabHandler(BaseHTTPRequestHandler):
    repository: PigmentationSequenceLabRepository

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _send_json(self, payload: object, status: int = 200) -> None:
        body = json.dumps(_sanitize_json(payload), allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, text: str, content_type: str = "text/html; charset=utf-8", status: int = 200) -> None:
        body = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        try:
            if parsed.path == "/":
                self._send_text(INDEX_HTML)
            elif parsed.path == "/api/options":
                self._send_json(self.repository.options_payload())
            elif parsed.path == "/api/model-window":
                self._handle_model_window(parsed.query)
            elif parsed.path == "/api/tracks":
                self._handle_tracks(parsed.query)
            elif parsed.path == "/api/track-options":
                self._handle_track_options(parsed.query)
            elif parsed.path == "/api/base-frequency":
                self._handle_base_frequency(parsed.query)
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        try:
            length = int(self.headers.get("Content-Length", "0") or "0")
            raw_body = self.rfile.read(length) if length else b"{}"
            body = json.loads(raw_body.decode("utf-8")) if raw_body else {}
            if parsed.path == "/api/edit/preview":
                self._handle_edit_preview(body)
            elif parsed.path == "/api/edit/apply":
                self._handle_edit_apply(body)
            elif parsed.path == "/api/edit/reset":
                self._send_json({"ok": True})
            else:
                self._send_json({"error": "not found"}, status=HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _handle_model_window(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.defaults
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        self._send_json(self.repository.model_window_payload(gene))

    def _handle_tracks(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.defaults
        subject_type = qs.get("subject_type", [defaults["subject_type"]])[0]
        subject = qs.get("subject", [defaults["subject"] or ""])[0]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        haplotype = qs.get("haplotype", [defaults["haplotype"]])[0]
        haplotypes = [h for raw in qs.get("haplotypes", []) for h in raw.split(",") if h] or [haplotype]
        output = qs.get("output", [defaults["output"]])[0]
        track = _safe_int(qs.get("track", ["0"])[0], 0, 0, 1_000_000)
        start = _safe_int(qs.get("start", [str(defaults["start"])])[0], defaults["start"], 1, 1_000_000_000)
        length = _safe_int(qs.get("length", [str(defaults["length"])])[0], defaults["length"], 1, 50_000_000)
        points = _safe_int(qs.get("points", [str(DEFAULT_POINTS)])[0], DEFAULT_POINTS, 1, MAX_POINTS)
        self._send_json(
            self.repository.tracks_payload(subject_type, subject, gene, haplotypes, output, track, start, length, points)
        )

    def _handle_track_options(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.defaults
        subject_type = qs.get("subject_type", [defaults["subject_type"]])[0]
        subject = qs.get("subject", [defaults["subject"] or ""])[0]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        haplotype = qs.get("haplotype", [defaults["haplotype"]])[0]
        output = qs.get("output", [defaults["output"]])[0]
        self._send_json(self.repository.track_options_payload(subject_type, subject, gene, haplotype, output))

    def _handle_base_frequency(self, query: str) -> None:
        qs = parse_qs(query)
        defaults = self.repository.defaults
        subject_type = qs.get("subject_type", [defaults["subject_type"]])[0]
        subject = qs.get("subject", [defaults["subject"] or ""])[0]
        gene = qs.get("gene", [defaults["gene"] or ""])[0]
        haplotype = qs.get("haplotype", [defaults["haplotype"]])[0]
        haplotypes = [h for raw in qs.get("haplotypes", []) for h in raw.split(",") if h] or [haplotype]
        start = _safe_int(qs.get("start", [str(defaults["start"])])[0], defaults["start"], 1, 1_000_000_000)
        length = _safe_int(qs.get("length", [str(defaults["length"])])[0], defaults["length"], 1, 50_000_000)
        points = _safe_int(qs.get("points", [str(DEFAULT_POINTS)])[0], DEFAULT_POINTS, 1, MAX_POINTS)
        self._send_json(self.repository.base_frequency_payload(subject_type, subject, gene, haplotypes, start, length, points))

    def _handle_edit_preview(self, body: Dict[str, Any]) -> None:
        if body.get("subject_type", "individual") != "individual":
            self._send_json({"error": "editing is only available for a single individual"}, status=HTTPStatus.BAD_REQUEST)
            return
        self._send_json(self.repository.edit_preview_payload(
            sample_id=body["sample_id"], gene=body["gene"], haplotype=body["haplotype"],
            start=int(body["start"]), end=int(body["end"]), op=body["op"],
            base=body.get("base"), seed=int(body.get("seed", 0)),
        ))

    def _handle_edit_apply(self, body: Dict[str, Any]) -> None:
        if body.get("subject_type", "individual") != "individual":
            self._send_json({"error": "editing is only available for a single individual"}, status=HTTPStatus.BAD_REQUEST)
            return
        self._send_json(self.repository.edit_apply_payload(
            sample_id=body["sample_id"], gene=body["gene"], haplotype=body["haplotype"],
            start=int(body["start"]), end=int(body["end"]), op=body["op"],
            base=body.get("base"), seed=int(body.get("seed", 0)),
            apply_to_prediction=body.get("apply_to_prediction", ["rna_seq"]),
        ))


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Pigmentation Sequence Lab</title>
  <style>
    :root { color-scheme: dark; }
    * { box-sizing: border-box; }
    body { margin: 0; font-family: -apple-system, Segoe UI, Roboto, sans-serif; background: #0f1420; color: #dbe7f3; }
    header { padding: 10px 16px; background: #141b2b; border-bottom: 1px solid #25344e; }
    h1 { font-size: 16px; margin: 0 0 8px 0; }
    .controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; }
    .controls label { font-size: 12px; color: #9ca9bf; display: flex; flex-direction: column; gap: 2px; }
    select, input[type=text], input[type=number] { background: #0f1420; color: #dbe7f3; border: 1px solid #25344e; border-radius: 4px; padding: 4px 6px; }
    button { background: #1c2740; color: #dbe7f3; border: 1px solid #334467; border-radius: 4px; padding: 5px 10px; cursor: pointer; }
    button:hover { background: #24314f; }
    button:disabled { opacity: 0.5; cursor: not-allowed; }
    main { padding: 12px 16px; display: flex; flex-direction: column; gap: 16px; }
    .panel { background: #141b2b; border: 1px solid #25344e; border-radius: 6px; padding: 10px; }
    .panel h2 { font-size: 13px; margin: 0 0 8px 0; color: #9ca9bf; text-transform: uppercase; letter-spacing: .04em; }
    svg { width: 100%; height: 220px; display: block; background: #0b101c; border-radius: 4px; }
    #logoSvg, #basesSvg { height: 130px; }
    .row { display: flex; gap: 8px; align-items: center; flex-wrap: wrap; }
    .hint { font-size: 11px; color: #7684a0; }
    .badge { display: inline-block; padding: 2px 8px; border-radius: 10px; font-size: 11px; background: #25344e; }
    .badge.warn { background: #4a3620; color: #ffd166; }
    .toolbar { display: none; gap: 8px; align-items: center; margin-top: 8px; padding: 8px; background: #0b101c; border-radius: 4px; }
    .toolbar.visible { display: flex; }
    #predictionCard { margin-top: 8px; font-size: 12px; }
    #predictionCard b { color: #ffd166; }
    #error { color: #ff6b6b; font-size: 12px; min-height: 16px; }
    .hidden { display: none !important; }
    svg.plot { cursor: grab; user-select: none; }
    svg.plot.no-pan { cursor: crosshair; }
    .panel-head { display: flex; align-items: baseline; justify-content: space-between; gap: 8px; }
  </style>
</head>
<body>
<header>
  <h1>Pigmentation Sequence Lab</h1>
  <div class="controls">
    <label>Subject
      <select id="subjectType">
        <option value="individual">Individual</option>
        <option value="class">Pigmentation class average</option>
      </select>
    </label>
    <label id="individualPickerLabel">Individual
      <select id="individualPicker"></select>
    </label>
    <label id="classPickerLabel" class="hidden">Class
      <select id="classPicker">
        <option value="strong">strong pigmentation</option>
        <option value="weak">weak pigmentation</option>
      </select>
    </label>
    <label>Gene
      <select id="gene"></select>
    </label>
    <label>Haplotypes
      <span class="row"><label><input type="checkbox" id="hapH1" checked> H1</label><label><input type="checkbox" id="hapH2" checked> H2</label></span>
    </label>
    <label>RNA-seq track
      <select id="trackRnaSeq"></select>
    </label>
    <label>CAGE track
      <select id="trackCage"></select>
    </label>
  </div>
  <div class="row" style="margin-top:8px;">
    <button id="centerModelWindow">Center CNN window</button>
    <button id="fullWindow">Full 512 kbp</button>
    <span class="hint" id="windowLabel"></span>
    <span class="hint">&middot; drag to pan &middot; shift+drag or scroll to zoom</span>
  </div>
  <div id="error"></div>
</header>
<main>
  <div class="panel">
    <h2>RNA-seq track</h2>
    <svg id="tracksRnaSeqSvg" class="plot"></svg>
  </div>
  <div class="panel">
    <div class="panel-head">
      <h2>CAGE track</h2>
      <span class="badge" id="cageBadge" style="display:none"></span>
    </div>
    <svg id="tracksCageSvg" class="plot"></svg>
  </div>
  <!-- Base-frequency (sequence logo) panel disabled for now: computing it for a whole
       pigmentation class (~500-700 individuals) via real bcftools_chain alignment was taking
       several minutes on the first view of a gene, and competing for CPU/GIL with concurrent
       CAGE fetches, slowing those down too. Re-enable once class-frequency caching/sampling is
       improved (see loadLogo() below, also commented out).
  <div class="panel">
    <h2>Base frequency (sequence logo) &mdash; <span id="logoClassLabel"></span></h2>
    <svg id="logoSvg" class="plot"></svg>
    <div class="hint" id="logoHint"></div>
  </div>
  -->
  <div class="panel" id="individualBasesPanel">
    <h2>Individual actual bases (bcftools_chain aligned, X = no aligned base)</h2>
    <svg id="basesSvg" class="plot no-pan"></svg>
    <div class="hint" id="basesHint"></div>
    <div class="toolbar" id="editToolbar">
      <span id="selectionLabel" class="hint"></span>
      <label>Overwrite with
        <select id="overwriteBase"><option>A</option><option>C</option><option>G</option><option>T</option></select>
      </label>
      <button id="overwriteBtn">Overwrite selection</button>
      <button id="scrambleBtn">Scramble selection</button>
      <button id="resetEditBtn">Reset edits</button>
    </div>
    <div id="predictionCard"></div>
  </div>
</main>
<script>
const $ = (id) => document.getElementById(id);
const COLORS = ['#5aa9ff', '#ff9f5a', '#7ee787', '#f778ba', '#ffd166'];
const BASE_COLORS = {A: '#4caf50', C: '#2196f3', G: '#ffb300', T: '#e53935', gap: '#555b6e'};
const FULL_WINDOW_LENGTH = 524288;
const MIN_VIEW_LENGTH = 20;
const NS = 'http://www.w3.org/2000/svg';

let options = null;
let view = {start: 1, length: 32768};
let editSelection = null; // {fastaStart, fastaEnd}
let lastIndividualBases = null; // {bases, source_indices, start, end}
let viewRefreshTimer = null;

function showError(err) { $('error').textContent = err && err.message ? err.message : String(err); console.error(err); }
function clearError() { $('error').textContent = ''; }

async function getJSON(url) {
  const resp = await fetch(url);
  const data = await resp.json();
  if (!resp.ok) throw new Error(data.error || `HTTP ${resp.status}`);
  return data;
}
async function postJSON(url, body) {
  const resp = await fetch(url, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body)});
  const data = await resp.json();
  if (!resp.ok) throw new Error(data.error || `HTTP ${resp.status}`);
  return data;
}

function currentSubjectType() { return $('subjectType').value; }
function currentSubject() {
  return currentSubjectType() === 'individual' ? $('individualPicker').value : $('classPicker').value;
}
function currentHaplotypes() {
  const hs = [];
  if ($('hapH1').checked) hs.push('H1');
  if ($('hapH2').checked) hs.push('H2');
  return hs.length ? hs : ['H1'];
}
function pigmentationClassForIndividual(sampleId) {
  const row = (options.individuals || []).find(r => r.sample_id === sampleId);
  return row ? row.pigmentation_class : null;
}

function populateSelectors() {
  $('individualPicker').innerHTML = (options.individuals || []).map(r =>
    `<option value="${r.sample_id}">${r.sample_id} (${r.population || '?'}, ${r.pigmentation_class || 'unknown'})</option>`
  ).join('');
  $('gene').innerHTML = (options.genes || []).map(g => `<option value="${g}">${g}</option>`).join('');
  if (options.defaults.subject) $('individualPicker').value = options.defaults.subject;
  if (options.defaults.gene) $('gene').value = options.defaults.gene;
  view = {start: options.defaults.start, length: options.defaults.length};
}

function toggleSubjectControls() {
  const isIndividual = currentSubjectType() === 'individual';
  $('individualPickerLabel').classList.toggle('hidden', !isIndividual);
  $('classPickerLabel').classList.toggle('hidden', isIndividual);
  $('individualBasesPanel').classList.toggle('hidden', !isIndividual);
}

function formatBp(n) { return Math.abs(n) >= 1000 ? `${(n / 1000).toFixed(1)}kb` : `${Math.round(n)}bp`; }

function updateWindowLabel() {
  $('windowLabel').textContent = `${formatBp(view.start)} - ${formatBp(view.start + view.length - 1)} (${formatBp(view.length)})`;
}

// -- view (pan/zoom) state ----------------------------------------------------
function setView(newStart, newLength) {
  newLength = Math.max(MIN_VIEW_LENGTH, Math.min(FULL_WINDOW_LENGTH, Math.round(newLength)));
  newStart = Math.max(1, Math.min(FULL_WINDOW_LENGTH - newLength + 1, Math.round(newStart)));
  view = {start: newStart, length: newLength};
  updateWindowLabel();
}

function scheduleViewRefresh(delay) {
  clearTimeout(viewRefreshTimer);
  viewRefreshTimer = setTimeout(() => { refreshView().catch(showError); }, delay);
}

function panByFraction(deltaFrac) {
  setView(view.start + deltaFrac * view.length, view.length);
  scheduleViewRefresh(30);
}

function zoomToFracRange(fracLo, fracHi) {
  setView(view.start + fracLo * view.length, (fracHi - fracLo) * view.length);
  scheduleViewRefresh(30);
}

function zoomByFactorAtFrac(frac, factor) {
  const anchor = view.start + frac * view.length;
  const newLength = view.length * factor;
  setView(anchor - frac * newLength, newLength);
  scheduleViewRefresh(150); // wheel events can fire rapidly; coalesce into the last one
}

// -- generic drag-to-pan / shift-drag-to-zoom / scroll-to-zoom interactions ---
function drawSelectionOverlay(svg, x1, x2) {
  let rect = svg.querySelector('.selection-overlay');
  if (!rect) {
    rect = document.createElementNS(NS, 'rect');
    rect.setAttribute('class', 'selection-overlay');
    rect.setAttribute('fill', 'rgba(255,209,102,0.25)');
    rect.setAttribute('stroke', '#ffd166');
    rect.setAttribute('y', 0);
    rect.setAttribute('height', '100%');
    svg.appendChild(rect);
  }
  rect.setAttribute('x', Math.min(x1, x2));
  rect.setAttribute('width', Math.max(1, Math.abs(x2 - x1)));
}
function clearSelectionOverlay(svg) {
  const rect = svg.querySelector('.selection-overlay');
  if (rect) rect.remove();
}

function attachInteractions(svg, {allowPlainDrag, onPan, onZoomRect, onWheelZoom, onPlainDragSelect}) {
  let drag = null; // {clientX, mode}

  function fracAt(clientX) {
    const rect = svg.getBoundingClientRect();
    return Math.min(1, Math.max(0, (clientX - rect.left) / rect.width));
  }

  svg.addEventListener('mousedown', (ev) => {
    const mode = ev.shiftKey ? 'zoom' : (allowPlainDrag ? 'pan' : (onPlainDragSelect ? 'select' : null));
    if (!mode) return;
    drag = {clientX: ev.clientX, mode};
    if (mode === 'pan') svg.style.cursor = 'grabbing';
    ev.preventDefault();
  });

  window.addEventListener('mousemove', (ev) => {
    if (!drag) return;
    if (drag.mode === 'pan') {
      const content = svg.querySelector('.plot-content');
      if (content) content.setAttribute('transform', `translate(${ev.clientX - drag.clientX},0)`);
    } else {
      const rect = svg.getBoundingClientRect();
      drawSelectionOverlay(svg, drag.clientX - rect.left, ev.clientX - rect.left);
    }
  });

  window.addEventListener('mouseup', (ev) => {
    if (!drag) return;
    svg.style.cursor = '';
    if (drag.mode === 'pan') {
      const content = svg.querySelector('.plot-content');
      if (content) content.removeAttribute('transform');
      const deltaFrac = (ev.clientX - drag.clientX) / svg.getBoundingClientRect().width;
      if (Math.abs(deltaFrac) > 0.002) onPan(-deltaFrac);
    } else {
      clearSelectionOverlay(svg);
      const fracLo = Math.min(fracAt(drag.clientX), fracAt(ev.clientX));
      const fracHi = Math.max(fracAt(drag.clientX), fracAt(ev.clientX));
      if (fracHi - fracLo > 0.01) {
        if (drag.mode === 'zoom') onZoomRect(fracLo, fracHi);
        else if (drag.mode === 'select' && onPlainDragSelect) onPlainDragSelect(fracLo, fracHi);
      }
    }
    drag = null;
  });

  svg.addEventListener('wheel', (ev) => {
    ev.preventDefault();
    onWheelZoom(fracAt(ev.clientX), ev.deltaY > 0 ? 1.3 : 1 / 1.3);
  }, {passive: false});
}

// -- line plots (RNA-seq / CAGE tracks) --------------------------------------
function drawLinePlot(svgId, series, modelWindow) {
  const svg = $(svgId);
  const width = svg.clientWidth || 900, height = svg.clientHeight || 220;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.innerHTML = '';
  const content = document.createElementNS(NS, 'g');
  content.setAttribute('class', 'plot-content');
  svg.appendChild(content);

  const pad = {l: 50, r: 16, t: 14, b: 24};
  const pts = [];
  for (const s of series) for (let i = 0; i < s.y.length; i++) if (Number.isFinite(s.y[i])) pts.push({x: s.x[i], y: s.y[i]});
  if (!pts.length) {
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', '50%'); t.setAttribute('y', '50%'); t.setAttribute('text-anchor', 'middle'); t.setAttribute('fill', '#7684a0');
    t.textContent = 'No finite values';
    content.appendChild(t);
    return;
  }
  const xmin = Math.min(...pts.map(p => p.x)), xmax = Math.max(...pts.map(p => p.x));
  let ymin = Math.min(...pts.map(p => p.y)), ymax = Math.max(...pts.map(p => p.y));
  if (ymin === ymax) { ymin -= 1; ymax += 1; }
  const sx = v => pad.l + ((v - xmin) / Math.max(1, xmax - xmin)) * (width - pad.l - pad.r);
  const sy = v => pad.t + (1 - (v - ymin) / (ymax - ymin)) * (height - pad.t - pad.b);
  if (modelWindow) {
    const mStart = (modelWindow.expanded_start || 0) + 1, mEnd = modelWindow.expanded_end || 0;
    if (mEnd >= xmin && mStart <= xmax) {
      const rect = document.createElementNS(NS, 'rect');
      rect.setAttribute('x', sx(Math.max(mStart, xmin))); rect.setAttribute('y', pad.t);
      rect.setAttribute('width', Math.max(0, sx(Math.min(mEnd, xmax)) - sx(Math.max(mStart, xmin))));
      rect.setAttribute('height', height - pad.t - pad.b);
      rect.setAttribute('fill', '#ffd16622');
      content.appendChild(rect);
    }
  }
  series.forEach((s, i) => {
    const p = s.y.map((v, idx) => ({x: s.x[idx], y: v})).filter(q => Number.isFinite(q.y));
    const d = p.map((q, idx) => `${idx ? 'L' : 'M'}${sx(q.x).toFixed(2)},${sy(q.y).toFixed(2)}`).join(' ');
    const path = document.createElementNS(NS, 'path');
    path.setAttribute('d', d); path.setAttribute('fill', 'none');
    path.setAttribute('stroke', COLORS[i % COLORS.length]); path.setAttribute('stroke-width', '1.8');
    content.appendChild(path);
    const label = document.createElementNS(NS, 'text');
    label.setAttribute('x', width - 160); label.setAttribute('y', 14 + i * 14);
    label.setAttribute('fill', COLORS[i % COLORS.length]); label.setAttribute('font-size', '11');
    label.textContent = s.label || `series ${i}`;
    content.appendChild(label);
  });
  [[pad.l, height - 6, xmin], [width - pad.r, height - 6, xmax]].forEach(([x, y, text]) => {
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', x); t.setAttribute('y', y); t.setAttribute('fill', '#7684a0'); t.setAttribute('font-size', '10');
    if (x > width / 2) t.setAttribute('text-anchor', 'end');
    t.textContent = formatBp(text);
    content.appendChild(t);
  });
}

const TRACK_PANELS = {
  rna_seq: {svg: 'tracksRnaSeqSvg', trackSelect: 'trackRnaSeq', badge: null},
  cage: {svg: 'tracksCageSvg', trackSelect: 'trackCage', badge: 'cageBadge'},
};

async function refreshTrackOptions(output) {
  const panel = TRACK_PANELS[output];
  const params = new URLSearchParams({
    subject_type: currentSubjectType(), subject: currentSubject(), gene: $('gene').value,
    haplotype: currentHaplotypes()[0], output,
  });
  const data = await getJSON(`/api/track-options?${params}`);
  const select = $(panel.trackSelect);
  const previous = select.value;
  select.innerHTML = data.tracks.map(t => `<option value="${t.index}">${t.label}</option>`).join('');
  if (previous && data.tracks.some(t => String(t.index) === previous)) select.value = previous;
}

async function loadTracks(output) {
  const panel = TRACK_PANELS[output];
  const params = new URLSearchParams({
    subject_type: currentSubjectType(), subject: currentSubject(), gene: $('gene').value,
    output, track: $(panel.trackSelect).value || '0',
    start: String(view.start), length: String(view.length), points: '600',
  });
  currentHaplotypes().forEach(h => params.append('haplotypes', h));
  const data = await getJSON(`/api/tracks?${params}`);
  drawLinePlot(panel.svg, data.series || [], data.model_window);
  if (panel.badge) {
    const badge = $(panel.badge);
    if (data.cage_is_true_average === false) {
      badge.style.display = 'inline-block';
      badge.textContent = `representative individual ${data.cage_representative_sample} (not a true class average yet)`;
    } else {
      badge.style.display = 'none';
    }
  }
}

// -- sequence logo / individual bases plots ----------------------------------
function drawFrequencyLogo(svgId, payload) {
  const svg = $(svgId);
  const width = svg.clientWidth || 900, height = svg.clientHeight || 130;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.innerHTML = '';
  const content = document.createElementNS(NS, 'g');
  content.setAttribute('class', 'plot-content');
  svg.appendChild(content);
  const symbols = ['A', 'C', 'G', 'T', 'gap'];

  if (payload.lod === 'letters') {
    const length = payload.mode === 'individual' ? payload.bases.length : payload.frequencies.A.length;
    if (!length) return;
    const colWidth = width / length;
    for (let i = 0; i < length; i++) {
      let order;
      if (payload.mode === 'individual') {
        const base = payload.bases[i];
        order = [{sym: base, freq: 1.0}];
      } else {
        order = symbols.map(sym => ({sym, freq: payload.frequencies[sym][i]})).filter(e => e.freq > 0).sort((a, b) => a.freq - b.freq);
      }
      let yCursor = height;
      for (const {sym, freq} of order) {
        const h = Math.max(1, freq * (height - 4));
        yCursor -= h;
        const text = document.createElementNS(NS, 'text');
        text.setAttribute('x', i * colWidth + colWidth / 2);
        text.setAttribute('y', yCursor + h * 0.85);
        text.setAttribute('text-anchor', 'middle');
        text.setAttribute('fill', BASE_COLORS[sym] || '#888');
        text.setAttribute('font-family', 'monospace');
        text.setAttribute('font-weight', 'bold');
        text.style.fontSize = `${h}px`;
        text.textContent = sym === 'gap' ? 'X' : sym;
        content.appendChild(text);
      }
    }
    return;
  }

  // density LOD: stacked frequency bars
  const xs = payload.x || [];
  if (!xs.length) return;
  const colWidth = width / xs.length;
  for (let i = 0; i < xs.length; i++) {
    let yCursor = height;
    for (const sym of symbols) {
      const freq = payload.frequencies[sym][i] || 0;
      const h = freq * height;
      yCursor -= h;
      if (h <= 0) continue;
      const rect = document.createElementNS(NS, 'rect');
      rect.setAttribute('x', i * colWidth); rect.setAttribute('y', yCursor);
      rect.setAttribute('width', Math.max(1, colWidth - 0.5)); rect.setAttribute('height', h);
      rect.setAttribute('fill', BASE_COLORS[sym] || '#888');
      content.appendChild(rect);
    }
  }
}

// Disabled for now (see the commented-out panel markup above): computing the group frequency
// logo for a whole pigmentation class (~500-700 individuals) was taking several minutes on the
// first view of a gene, and the concurrent bcftools_chain work was competing for CPU/GIL with
// the CAGE fetch running in parallel, slowing that down too.
// async function loadLogo() {
//   const classForLogo = currentSubjectType() === 'class' ? currentSubject() : pigmentationClassForIndividual(currentSubject());
//   $('logoClassLabel').textContent = classForLogo || '(unknown class)';
//   if (!classForLogo) { $('logoSvg').innerHTML = ''; return; }
//   const params = new URLSearchParams({
//     subject_type: 'class', subject: classForLogo, gene: $('gene').value,
//     start: String(view.start), length: String(view.length), points: '600',
//   });
//   currentHaplotypes().forEach(h => params.append('haplotypes', h));
//   const data = await getJSON(`/api/base-frequency?${params}`);
//   drawFrequencyLogo('logoSvg', data);
//   $('logoHint').textContent = data.lod === 'letters'
//     ? `letters (n=${data.n_used}/${data.n_requested} individuals)`
//     : `zoom in to <= ${options.letter_lod_threshold}bp to see individual letters (n=${data.n_used}/${data.n_requested} individuals)`;
// }

async function loadIndividualBases() {
  if (currentSubjectType() !== 'individual') return;
  const params = new URLSearchParams({
    subject_type: 'individual', subject: currentSubject(), gene: $('gene').value,
    haplotype: currentHaplotypes()[0], start: String(view.start), length: String(view.length), points: '600',
  });
  const data = await getJSON(`/api/base-frequency?${params}`);
  drawFrequencyLogo('basesSvg', data);
  if (data.lod === 'letters') {
    lastIndividualBases = {bases: data.bases, source_indices: data.source_indices, start: data.start, end: data.end};
    $('basesHint').textContent = `${data.bases.length}bp shown — drag to select a region to edit, shift+drag or scroll to zoom`;
  } else {
    lastIndividualBases = null;
    $('editToolbar').classList.remove('visible');
    $('basesHint').textContent = `zoom in to <= ${options.letter_lod_threshold}bp to select/edit bases (shift+drag or scroll to zoom)`;
  }
}

function onIndividualBasesPlainDrag(fracLo, fracHi) {
  if (!lastIndividualBases) return;
  const n = lastIndividualBases.bases.length;
  const lo = Math.max(0, Math.min(n - 1, Math.floor(fracLo * n)));
  const hi = Math.max(0, Math.min(n - 1, Math.floor(fracHi * n) - 1));
  if (hi < lo) return;
  const src = lastIndividualBases.source_indices.slice(lo, hi + 1).filter(v => v !== null);
  if (!src.length) { showError(new Error('Selection has no editable bases (entirely a gap/deletion).')); return; }
  editSelection = {fastaStart: Math.min(...src), fastaEnd: Math.max(...src) + 1};
  $('editToolbar').classList.add('visible');
  $('selectionLabel').textContent = `selected ${hi - lo + 1} rendered bp -> FASTA[${editSelection.fastaStart}:${editSelection.fastaEnd}]`;
}

async function applyEdit(op) {
  if (!editSelection) return;
  clearError();
  const body = {
    subject_type: 'individual', sample_id: currentSubject(), gene: $('gene').value,
    haplotype: currentHaplotypes()[0], start: editSelection.fastaStart, end: editSelection.fastaEnd,
    op, base: $('overwriteBase').value, apply_to_prediction: ['rna_seq', 'cage'],
  };
  $('overwriteBtn').disabled = true; $('scrambleBtn').disabled = true;
  try {
    const result = await postJSON('/api/edit/apply', body);
    const delta = result.delta_strong_pigmentation_probability;
    $('predictionCard').innerHTML = `P(strong pigmentation): <b>${result.baseline_strong_pigmentation_probability.toFixed(3)}</b> &rarr; <b>${result.edited_strong_pigmentation_probability.toFixed(3)}</b> (&Delta; ${delta >= 0 ? '+' : ''}${delta.toFixed(3)})` +
      (result.affects_prediction ? '' : ' <span class="badge warn">edit is outside the CNN window: prediction unchanged</span>');
    await refreshView();
  } catch (err) { showError(err); } finally {
    $('overwriteBtn').disabled = false; $('scrambleBtn').disabled = false;
  }
}

async function resetEdits() {
  await postJSON('/api/edit/reset', {});
  editSelection = null;
  $('editToolbar').classList.remove('visible');
  $('predictionCard').innerHTML = '';
  await refreshView();
}

// -- top-level refresh: "view" (pan/zoom only) vs. "all" (subject/gene/etc.) -
async function refreshView() {
  clearError();
  updateWindowLabel();
  await Promise.all([loadTracks('rna_seq'), loadTracks('cage'), /* loadLogo(), */ loadIndividualBases()]);
}

async function refreshAll() {
  clearError();
  toggleSubjectControls();
  updateWindowLabel();
  try {
    await Promise.all([refreshTrackOptions('rna_seq'), refreshTrackOptions('cage')]);
    await refreshView();
  } catch (err) { showError(err); }
}

async function centerModelWindow() {
  const data = await getJSON(`/api/model-window?gene=${encodeURIComponent($('gene').value)}`);
  setView(data.start, data.length);
  await refreshView();
}
async function fullWindowView() {
  setView(1, FULL_WINDOW_LENGTH);
  await refreshView();
}

async function init() {
  options = await getJSON('/api/options');
  populateSelectors();
  toggleSubjectControls();
  ['subjectType', 'individualPicker', 'classPicker', 'gene', 'hapH1', 'hapH2'].forEach(id => {
    $(id).addEventListener('change', () => refreshAll().catch(showError));
  });
  $('trackRnaSeq').addEventListener('change', () => loadTracks('rna_seq').catch(showError));
  $('trackCage').addEventListener('change', () => loadTracks('cage').catch(showError));
  $('centerModelWindow').addEventListener('click', () => centerModelWindow().catch(showError));
  $('fullWindow').addEventListener('click', () => fullWindowView().catch(showError));
  $('overwriteBtn').addEventListener('click', () => applyEdit('overwrite'));
  $('scrambleBtn').addEventListener('click', () => applyEdit('scramble'));
  $('resetEditBtn').addEventListener('click', () => resetEdits().catch(showError));

  const zoomWheelPan = {onZoomRect: zoomToFracRange, onWheelZoom: zoomByFactorAtFrac};
  attachInteractions($('tracksRnaSeqSvg'), {allowPlainDrag: true, onPan: panByFraction, ...zoomWheelPan});
  attachInteractions($('tracksCageSvg'), {allowPlainDrag: true, onPan: panByFraction, ...zoomWheelPan});
  // attachInteractions($('logoSvg'), ...) removed: the logo panel is commented out above.
  attachInteractions($('basesSvg'), {allowPlainDrag: false, onPlainDragSelect: onIndividualBasesPlainDrag, ...zoomWheelPan});

  await refreshAll();
}

init().catch(showError);
</script>
</body>
</html>
"""


def serve(
    dataset_dir: Path,
    config_path: Path,
    host: str,
    port: int,
    checkpoint: str = "best_accuracy",
    consensus_dataset_dir: Optional[Path] = None,
    api_key: Optional[str] = None,
) -> None:
    repository = PigmentationSequenceLabRepository(
        dataset_dir, config_path, checkpoint=checkpoint, consensus_dataset_dir=consensus_dataset_dir, api_key=api_key,
    )
    PigmentationSequenceLabHandler.repository = repository
    server = ThreadingHTTPServer((host, port), PigmentationSequenceLabHandler)
    print(f"Pigmentation Sequence Lab: http://{host}:{port}")
    print(f"Dataset: {repository.dataset_dir}")
    print(f"Config: {config_path}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Interactive pigmentation CNN sequence/prediction lab")
    parser.add_argument("dataset_dir", type=Path, help="Dataset directory containing dataset_metadata.json")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml"),
        help="genotype_based pigmentation config (bcftools_chain alignment_mapping, CNN2 model)",
    )
    parser.add_argument("--checkpoint", default="best_accuracy", help="Checkpoint alias/path: best_accuracy, best_loss, final or .pt path")
    parser.add_argument("--host", default="127.0.0.1", help="Bind host")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help="Bind port")
    parser.add_argument("--consensus-dataset-dir", type=Path, default=DEFAULT_CONSENSUS_DATASET_DIR)
    parser.add_argument("--api-key", default=None, help="AlphaGenome API key (defaults to ALPHAGENOME_API_KEY env or ~/.env)")
    args = parser.parse_args(argv)
    serve(
        args.dataset_dir, args.config, args.host, args.port,
        checkpoint=args.checkpoint, consensus_dataset_dir=args.consensus_dataset_dir, api_key=args.api_key,
    )


if __name__ == "__main__":
    main()
