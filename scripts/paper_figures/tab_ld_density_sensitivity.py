#!/usr/bin/env python3
"""Reviewer ask (Q7): are the Sidak/GATES/mean-chi2 gene scores sensitive to per-window
variant density or LD heterogeneity, and does either differ systematically between panel
and control (a confound that would undermine the whole GWAS-arm comparison)?

Reads the already-computed gene_scores JSON (gwas_ranking/gene_ranking_uncorrected_44.json)
-- no new experiment, purely a reanalysis of numbers the paper already reports. m_eff/m
(Li & Ji effective-test count over raw variant count) is the LD-heterogeneity proxy: 1 means
every variant in the window is independent, lower means heavier redundancy.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import mannwhitneyu, pearsonr

from _common import CONTROL, PANEL, RES, TABDIR

GWAS_JSON = RES / "gwas_ranking" / "gene_ranking_uncorrected_44.json"


def main():
    gs = json.loads(GWAS_JSON.read_text())["gene_scores"]["window"]
    genes = PANEL + CONTROL
    m = np.array([gs[g]["m"] for g in genes], float)
    ld_ratio = np.array([gs[g]["m_eff_li_ji"] for g in genes], float) / m
    scores = {"S_sidak": np.array([gs[g]["S_sidak"] for g in genes], float),
              "S_meanchi2": np.array([gs[g]["S_meanchi2"] for g in genes], float),
              "S_gates": np.array([gs[g]["S_gates"] for g in genes], float)}
    is_panel = np.array([g in PANEL for g in genes])

    rows = []
    for name, s in scores.items():
        r_m, p_m = pearsonr(m, s)
        r_ld, p_ld = pearsonr(ld_ratio, s)
        rows.append((name, r_m, p_m, r_ld, p_ld))

    lines = [r"\begin{tabular}{lrrrr}", r"\toprule",
             r"estimator & $r$ vs.\ $m$ & $p$ & $r$ vs.\ $m_{\mathrm{eff}}/m$ & $p$ \\",
             r"\midrule"]
    label = {"S_sidak": r"{\v S}id{\'a}k", "S_meanchi2": r"mean $\chi^2$", "S_gates": "GATES"}
    for name, r_m, p_m, r_ld, p_ld in rows:
        lines.append(f"{label[name]} & ${r_m:+.3f}$ & {p_m:.3f} & ${r_ld:+.3f}$ & {p_ld:.3f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    out = TABDIR / "tab-ld-density-sensitivity.tex"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")
    print()

    for name in ("m", "ld_ratio"):
        arr = m if name == "m" else ld_ratio
        u, p = mannwhitneyu(arr[is_panel], arr[~is_panel])
        print(f"{name}: panel_mean={arr[is_panel].mean():.4f} control_mean={arr[~is_panel].mean():.4f} "
              f"MWU p={p:.3f}")

    gap = scores["S_sidak"] - scores["S_meanchi2"]
    for name, arr in (("m", m), ("m_eff/m", ld_ratio)):
        r, p = pearsonr(arr, gap)
        print(f"gap(sidak - meanchi2) vs {name}: r={r:+.3f} p={p:.3f}")


if __name__ == "__main__":
    main()
