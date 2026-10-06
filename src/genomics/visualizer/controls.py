"""Matched control windows: a negative control for "is this gene panel special?".

A model trained on a panel of phenotype genes is only evidence about those genes if a panel of
phenotype-irrelevant windows, reading comparable signal, does worse. Picking irrelevant windows is
easy; picking windows that deliver a *comparable* amount of signal to the network is the part that
decides whether the control can be read at all, since a window the model barely sees is flat for
reasons that have nothing to do with biology.

So candidates are matched to the panel on the statistic this project already uses for it (see
``scripts/experiments/specificity_control_preflight.py``): the total predicted signal inside the
crop the model actually reads, on the *reference* window, summed over the chosen ontology terms and
strands. It needs no AlphaGenome calls — the reference predictions are already on disk — and over
the published panel it tracks the realised perturbation at Spearman rho = 0.855.

Matching is on log10 of that total, so pairs are compared by ratio rather than by difference, and
is an exact minimum-total-cost 1:1 assignment (``match_windows``).
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from genomics.visualizer.datasets import Dataset
from genomics.visualizer.jobs import JobCancelled, is_cancelled


class ControlError(ValueError):
    """No reference predictions, or no candidate windows to match against."""


def crop_signal(signals, dataset: Dataset, gene: str, output: str, terms: Sequence[str], window_center_size: int) -> float:
    """Total predicted signal inside the model's crop of the reference window of ``gene``.

    Sums every track whose ontology term is in ``terms`` (both strands), over the centred crop of
    ``window_center_size`` bases — the same rows and columns the trained model reads.
    """
    matrix = signals.reference_prediction(dataset, gene, output)
    tracks = (dataset.gene_info(gene)["outputs"].get(output) or {}).get("tracks") or []
    wanted = {str(t) for t in terms}
    columns = [t["index"] for t in tracks if str((t.get("metadata") or {}).get("ontology_curie")) in wanted]
    if not columns:
        raise ControlError(f"{gene} has no {output} tracks for the chosen ontology terms")
    window = dataset.model_window(gene, window_center_size) or {}
    start, end = int(window.get("start", 0)), int(window.get("end", matrix.shape[0]))
    block = matrix[start:end, [c for c in columns if c < matrix.shape[1]]]
    return float(np.nansum(block))


def match_windows(panel: Dict[str, float], candidates: Dict[str, float]) -> List[Dict[str, Any]]:
    """Assign each panel gene one candidate window, minimising the total |log10(ratio)|.

    Costs are absolute differences on a line (log10 of the signal), where an optimal assignment
    never crosses, so sorting both sides and choosing which candidates to skip (a O(n*m) dynamic
    program) gives the exact minimum rather than a greedy approximation.
    """
    if not panel:
        raise ControlError("No panel windows to match")
    if len(candidates) < len(panel):
        raise ControlError(f"{len(candidates)} candidate window(s) for {len(panel)} panel window(s); need at least as many")
    logs = lambda value: math.log10(max(float(value), 1e-9))  # noqa: E731 - a 0 total is a floor, not an error
    p = sorted(panel.items(), key=lambda kv: (logs(kv[1]), kv[0]))
    c = sorted(candidates.items(), key=lambda kv: (logs(kv[1]), kv[0]))
    n, m = len(p), len(c)
    inf = float("inf")
    # dp[i][j]: best cost matching the first i panel genes within the first j candidates.
    dp = [[inf] * (m + 1) for _ in range(n + 1)]
    for j in range(m + 1):
        dp[0][j] = 0.0
    for i in range(1, n + 1):
        for j in range(i, m + 1):
            take = dp[i - 1][j - 1] + abs(logs(p[i - 1][1]) - logs(c[j - 1][1]))
            dp[i][j] = min(dp[i][j - 1], take)
    pairs: List[Dict[str, Any]] = []
    i, j = n, m
    while i > 0:
        if dp[i][j] == dp[i][j - 1]:  # candidate j-1 was skipped
            j -= 1
            continue
        gene, value = p[i - 1]
        control, control_value = c[j - 1]
        pairs.append({
            "gene": gene, "signal": value, "control": control, "control_signal": control_value,
            "ratio": (float(control_value) / float(value)) if value else None,
        })
        i -= 1
        j -= 1
    pairs.reverse()
    return pairs


class ControlService:
    """Builds matched control panels for the training form (reference predictions only)."""

    def __init__(self, signals):
        self.signals = signals

    @staticmethod
    def candidates(dataset: Dataset, panel: Sequence[str]) -> List[str]:
        """Windows that can stand in for ``panel``: every other window not named in the metadata.

        The metadata's ``genes`` are the dataset's phenotype panel, so windows outside it were added
        as controls; a window already in ``panel`` cannot also be its own control.
        """
        listed = {str(g) for g in (dataset.metadata.get("genes") or [])}
        chosen = {str(g) for g in panel}
        return [g for g in dataset.genes if g not in chosen and g not in listed]

    def signals_for(self, dataset: Dataset, genes: Sequence[str], output: str, terms: Sequence[str],
                    window_center_size: int, progress=None) -> Tuple[Dict[str, float], Dict[str, str]]:
        """Crop signal per gene, plus the reason each skipped gene has none."""
        values: Dict[str, float] = {}
        skipped: Dict[str, str] = {}
        for i, gene in enumerate(genes):
            if is_cancelled(progress):
                raise JobCancelled()
            if progress is not None:
                progress(i / max(1, len(genes)), f"{gene} ({i + 1}/{len(genes)})")
            try:
                values[gene] = crop_signal(self.signals, dataset, gene, output, terms, window_center_size)
            except (FileNotFoundError, ControlError, KeyError) as exc:
                skipped[gene] = str(exc)
        return values, skipped

    def matched_panel(self, dataset: Dataset, panel: Sequence[str], output: str, terms: Sequence[str],
                      window_center_size: int, progress=None) -> Dict[str, Any]:
        """A control window per panel window, matched on crop signal, with the pairing table."""
        panel = [g for g in dict.fromkeys(str(g) for g in panel)]
        if not panel:
            raise ControlError("Choose at least one window for the panel")
        unknown = [g for g in panel if g not in dataset.genes]
        if unknown:
            raise ControlError(f"Unknown window(s): {', '.join(unknown)}")
        candidates = self.candidates(dataset, panel)
        if not candidates:
            raise ControlError("This dataset has no windows outside its gene panel to use as controls")
        panel_values, panel_skipped = self.signals_for(dataset, panel, output, terms, window_center_size, progress)
        if panel_skipped:
            raise ControlError(
                f"No reference {output} prediction for {', '.join(sorted(panel_skipped))}. "
                "Predict the reference window of every chosen gene first (haplotype 'ref')."
            )
        candidate_values, candidate_skipped = self.signals_for(dataset, candidates, output, terms, window_center_size, progress)
        if len(candidate_values) < len(panel_values):
            raise ControlError(
                f"Only {len(candidate_values)} control window(s) have a reference {output} prediction, "
                f"for a panel of {len(panel_values)}. Predict the reference window of the control windows first."
            )
        pairs = match_windows(panel_values, candidate_values)
        total_panel = sum(panel_values.values())
        total_control = sum(p["control_signal"] for p in pairs)
        ratios = [abs(math.log10(p["ratio"])) for p in pairs if p.get("ratio")]
        return {
            "output": output,
            "ontology_terms": list(terms),
            "window_center_size": int(window_center_size),
            "pairs": pairs,
            "controls": [p["control"] for p in pairs],
            "panel_total": total_panel,
            "control_total": total_control,
            "total_ratio": (total_control / total_panel) if total_panel else None,
            "worst_ratio": max((p["ratio"] for p in pairs if p.get("ratio")), key=lambda r: abs(math.log10(r)), default=None),
            "mean_abs_log_ratio": (sum(ratios) / len(ratios)) if ratios else None,
            "candidates": len(candidate_values),
            "skipped": candidate_skipped,
        }
