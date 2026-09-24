#!/usr/bin/env python3
"""Tables 1 and 2: per-arm accuracy, and per-arm knockdown readouts.

Both are written as LaTeX `tabular` bodies under paper/tables/ and pulled in with
\\input, so a rerun of the sweep updates the paper without an edit.
"""
from __future__ import annotations

import numpy as np

import pandas as pd

from _common import (CONTROL, DRAW1, DRAW2, DRAW3, KD_DIR, PANEL, TABDIR,
                     load_final_table, load_gwas)

BASELINE = 0.7099


def fmt(x, n=3, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    return f"{x:+.{n}f}" if sign else f"{x:.{n}f}"


def tex_escape(g):
    return g.replace("_", r"\_")


def thousands(x, n=0):
    return f"{x:,.{n}f}".replace(",", r"\,")


def accuracy_table(df):
    """Balanced accuracy per arm on both splits, panel and controls, ordered within class.

    Ordered by VALIDATION accuracy. No selection is marked: the multi-gene probe this
    table's rows once fed belongs to a version of the paper that no longer exists.
    """
    try:
        gw = load_gwas("uncorrected", "44")["gene_scores"]["window"]
    except FileNotFoundError:
        gw = {}
    lines = [r"\begin{tabular}{llrrrrrr}", r"\toprule",
             r"& gene & \multicolumn{2}{c}{bal.\ acc.} & acc. & rec.\ strong "
             r"& rec.\ weak & $m$ \\",
             r"\cmidrule(lr){3-4}",
             r"& & val. & test & \multicolumn{3}{c}{\textit{test}} & \\",
             r"\midrule"]
    for label, genes in (("Panel", PANEL), ("Control, draw 1", DRAW1),
                         ("Control, draw 2", DRAW2), ("Control, draw 3", DRAW3)):
        sub = df[df["gene"].isin(genes)].sort_values("bal_acc_val", ascending=False)
        if sub.empty:
            continue
        lines.append(rf"\multicolumn{{8}}{{l}}{{\textit{{{label}}} ({len(sub)} arms)}} \\")
        for _, r in sub.iterrows():
            m = gw.get(r["gene"], {}).get("m")
            lines.append(
                f"& {tex_escape(r['gene'])} & {fmt(r['bal_acc_val'])} & "
                f"{fmt(r['bal_acc'])} & {fmt(r['acc'])} & "
                f"{fmt(r['rec_strong'])} & {fmt(r['rec_weak'])} & "
                f"{'--' if m is None else thousands(m)} \\\\")
        lines.append(r"\addlinespace")
    lines += [r"\midrule",
              rf"& \textit{{majority-class baseline}} & {BASELINE:.4f} & {BASELINE:.4f} "
              rf"& {BASELINE:.4f} & 1.000 & 0.000 & -- \\",
              r"\bottomrule", r"\end{tabular}"]
    out = TABDIR / "tab-balacc.tex"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out} ({len(df)} arms)")


def log2fc_per_arm():
    """Mean exonic log2 fold change per single-gene arm, from that arm's own replay CSV.

    Not in poolmax_final_table.csv, which predates the column. Read here rather than
    added there because several finished figures read that table and would silently
    change shape.
    """
    out = {}
    for f in sorted(KD_DIR.glob("single_*.csv")):
        d = pd.read_csv(f)
        d = d[d["method"] == "biology_tss"]
        if d.empty or "expr_log2fc" not in d.columns:
            continue
        g = str(d["gene"].iloc[0])
        v = pd.to_numeric(d["expr_log2fc"], errors="coerce").dropna()
        if len(v):
            out[g] = float(v.mean())
    return out


def perturbation_table(df):
    """The knockdown readouts per arm, plus the flip rate and the spread of Delta."""
    fc = log2fc_per_arm()
    lines = [r"\begin{tabular}{llrrrrrrr}", r"\toprule",
             # Delta_in is raw AlphaGenome units and is NOT comparable across genes,
             # because normalisation divides each gene by its own log-max; the tilde
             # column is the same perturbation measured on the normalised tensor the
             # network consumes, which is the comparable one. The old header called it
             # Delta_in/log(1+m), which it is not: the log is applied per position
             # before the division, so it is not a rescaling of the raw sum.
             r"& gene & $\Delta_{\mathrm{in}}$ & $\tilde{\Delta}_{\mathrm{in}}$ & "
             r"$\Delta_{\mathrm{expr}}$ & $\Delta_{\mathrm{expr}}^{\mathrm{rel}}$ & "
             r"$\log_2$FC & $\Delta$ (log-odds) & flip \\",
             r"\midrule"]
    for label, genes in (("Panel", PANEL), ("Control, draw 1", DRAW1),
                         ("Control, draw 2", DRAW2), ("Control, draw 3", DRAW3)):
        sub = df[df["gene"].isin(genes)].copy()
        if sub.empty:
            continue
        sub = sub.reindex(sub["delta"].abs().sort_values(ascending=False).index)
        lines.append(rf"\multicolumn{{9}}{{l}}{{\textit{{{label}}} ({len(sub)} arms)}} \\")
        for _, r in sub.iterrows():
            v = fc.get(r["gene"])
            # Red for a positive fold change: the scramble RAISED predicted expression,
            # so the intervention failed on its own terms. Same convention as the
            # top-10 table.
            fctxt = "--" if v is None else (
                rf"\textcolor{{red}}{{{fmt(v, 3, True)}}}" if v > 0 else fmt(v, 3, True))
            lines.append(
                f"& {tex_escape(r['gene'])} & {thousands(r['delta_in'])} & "
                f"{thousands(r['delta_in_norm'])} & {thousands(r['delta_expr'])} & "
                f"{fmt(r['delta_expr_rel'], 3, True)} & {fctxt} & "
                f"{fmt(r['delta'], 3, True)} $\\pm$ {fmt(r['delta_sd'], 2)} & "
                f"{fmt(r['flip_rate'], 3)} \\\\")
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}"]
    out = TABDIR / "tab-perturbation.tex"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out} ({len(df)} arms)")


if __name__ == "__main__":
    df = load_final_table()
    missing = sorted(set(PANEL + CONTROL) - set(df["gene"]))
    if missing:
        print(f"NOTE: {len(missing)} arms not yet in the final table: {missing}")
    accuracy_table(df)
    perturbation_table(df)
