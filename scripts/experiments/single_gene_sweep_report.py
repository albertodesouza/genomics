#!/usr/bin/env python3
"""Aggregate the one-classifier-per-gene sweep into the panel-vs-control comparison.

The panel probe measures RELIANCE: a property of the trained classifier that depends on
what else is in its input, so a control gene can look flat merely because the panel genes
crowded it out. Training one classifier per gene pins reliance at its maximum -- with
nothing else to route through, the classifier must use that locus or fail -- so a flat
response there cannot be explained by competition. That is what makes this the specificity
control Section 5 says the paper does not have.

Two readouts per arm, on the same footing for every gene:

  acc   the accuracy the window alone supports, against the majority-class baseline. This
        is NOT a measure of pigmentation information: a control gene tops the table, so
        the axis cannot separate panel from control and is reported for context only.
  |D|   the same 100 bp promoter scramble, replayed through that arm's own checkpoint.
        This is the axis the design is for.

|D_in| is reported on two scales because they disagree about everything. It is summed
BEFORE normalisation, but each arm divides its tracks by a log_max fitted on that gene
alone, and across these arms that divisor spans ~760x (LRRC36 2.3e-03 to DDB1 1.8e+00).
So the raw column is not the perturbation the network receives. On the raw scale the
response looks free of the delivered-magnitude confound (rho ~ +0.19, n.s.); divided by
log_max it is not (rho ~ +0.57, p ~ 0.006). The consequence is that |D| is NOT comparable
across genes without this column: a gene whose window is nearly unexpressed has its
scramble amplified into the same input range as a gene delivering 400x more raw signal.
log_max varies per track within a gene too, so the median used here is an approximation --
the exact figure needs per-track deltas, which the knockdown CSVs do not carry.

Delivery matching is the crux. Section 5 prescribes that control genes be matched to the
responding genes on delivered perturbation rather than merely chosen for biological
irrelevance, so each control is paired here with the panel gene nearest it in |D_in| and
the response ratio at matched delivery is reported.

Reads only CSVs already on disk. No AlphaGenome calls, no training.
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
RUNS = REPO_ROOT / "results/genotype_based_predictor"
FOURTEEN = KO_DIR / "pigmentation_test_split_knockout_14gene.csv"
# |Delta| as published in Table II, for the reliance-vs-response contrast.
PANEL_DELTA = {"SLC24A5": .503, "TYR": .337, "SLC45A2": .259, "TYRP1": .251, "DDB1": .240,
               "MFSD12": .240, "MC1R": .055, "HERC2": .018, "OCA2": .002, "EDAR": .002,
               "TCHH": .000}
OUT = RUNS / "single_gene_sweep_report.json"
# Class is read from the control gene lists, not from a knockdown CSV: the second batch of
# eight was built after the last knockdown run and appears in none of them.
CONTROL_LISTS = [REPO_ROOT / "docs/gene_lists/specificity_control_genes.txt",
                 REPO_ROOT / "docs/gene_lists/specificity_control_genes_8.txt",
                 REPO_ROOT / "docs/gene_lists/specificity_control_genes_11b.txt"]

# EDAR and TCHH are not pigmentation genes. EDAR (rs3827760, V370A) is a hair-thickness,
# tooth-shape and sweat-gland locus and the cause of hypohidrotic ectodermal dysplasia;
# TCHH (trichohyalin) is a hair-shape locus. Neither has a pigmentation phenotype of any
# evidence class -- Mendelian, GWAS or functional -- and neither carries a pigmentation GO
# annotation. They entered the panel from the predecessor pipeline, whose three retained
# ontologies include two hair-follicle cell types, and are excluded here because this
# study is about pigmentation.
#
# The criterion is human phenotype, NOT GO annotation. Screening the panel for the
# PRESENCE of the same GO terms used to screen the controls for their ABSENCE would also
# drop SLC24A5, HERC2 and DDB1: GO records curated molecular process, and SLC24A5 is
# annotated as a cation exchanger even though its variant rs1426654 is the largest-effect
# light-skin allele known and its loss causes oculocutaneous albinism type 6. That
# criterion is unfit on this side of the comparison, so it is not used.
#
# The controls' irrelevance is established by three screens that agree -- no pigmentation
# GO annotation, no Mendelian pigmentary disorder, no membership of KEGG melanogenesis
# (hsa04916) -- plus a locus-level check that no control WINDOW overlaps any of 58
# established pigmentation loci, which matters because the classifier reads 524,288 bp and
# not a gene.
NOT_PIGMENTATION = {"EDAR", "TCHH"}


def p_weak(strong: np.ndarray, weak: np.ndarray) -> np.ndarray:
    z = np.stack([strong, weak], 1)
    z = z - z.max(1, keepdims=True)
    e = np.exp(z)
    return e[:, 1] / e.sum(1)


def main() -> int:
    if not FOURTEEN.exists():
        raise SystemExit(f"ABORT: {FOURTEEN} missing; the 14-gene knockdown must run first.")
    ref = pd.read_csv(FOURTEEN)
    ref = ref[ref.method == "biology_tss"]
    delivery = ref.groupby("gene")["delta_in"].mean().to_dict()
    controls = set()
    for lp in CONTROL_LISTS:
        if not lp.exists():
            raise SystemExit(f"ABORT: control gene list {lp} missing; class would be wrong.")
        controls |= {g.strip().upper() for g in lp.read_text().split() if g.strip()}
    lab = ref[ref.gene == ref.gene.iloc[0]].true_label.value_counts()
    baseline = float(lab.max() / lab.sum())

    rows = []
    for f in sorted(glob.glob(str(KO_DIR / "pigmentation_test_split_knockout_single_*.csv"))):
        arm = re.search(r"single_(\w+)\.csv", f).group(1)
        gene = arm.upper()
        d = pd.read_csv(f)
        d = d[d.method == "biology_tss"]
        if d.empty:
            print(f"  WARN: {gene} has no biology_tss rows; skipped", file=sys.stderr)
            continue
        dd = (p_weak(d.perturbed_strong_logit.values, d.perturbed_weak_logit.values)
              - p_weak(d.baseline_strong_logit.values, d.baseline_weak_logit.values))
        acc = float("nan")
        for j in glob.glob(str(RUNS / f"runs_single_gene_{arm}" / "*" / "test_best_accuracy_results.json")):
            acc = float(json.load(open(j)).get("weighted_accuracy", float("nan")))
        # Delivery from this arm's own rows, falling back to the 14-gene run. |D_in| is
        # summed before normalisation, so the two agree wherever both exist (checked:
        # TYR 55330, PSMC4 12605 either way) -- but only the arm's own CSV has the eight.
        din = float(d["delta_in"].mean()) if "delta_in" in d else float("nan")
        if din != din:
            din = float(delivery.get(gene, float("nan")))
        # Median per-track log_max of this arm: the scale its input was divided by.
        logmax = float("nan")
        for npj in glob.glob(str(RUNS / f"runs_single_gene_{arm}" / "*" / "models"
                                 / "normalization_params.json")):
            logmax = float(np.median([t["log_max"]
                                      for t in json.load(open(npj))["track_params"]]))
        # Per-arm null band, where the arm replayed the random-window draw as well.
        nullrows = pd.read_csv(f)
        nullrows = nullrows[nullrows.method.astype(str).str.startswith("random_window")]
        if nullrows.empty:
            null99 = float("nan")
        else:
            nd = (p_weak(nullrows.perturbed_strong_logit.values, nullrows.perturbed_weak_logit.values)
                  - p_weak(nullrows.baseline_strong_logit.values, nullrows.baseline_weak_logit.values))
            null99 = float(np.percentile(np.abs(nd), 99))
        rows.append({
            "gene": gene, "class": "control" if gene in controls else "panel", "n": int(len(d)),
            "acc": acc, "acc_over_baseline": acc - baseline,
            "abs_delta_single": float(abs(dd.mean())),
            "abs_delta_panel": PANEL_DELTA.get(gene, float("nan")),
            "delta_in": din, "log_max_median": logmax,
            "delta_in_normalised": din / logmax if logmax == logmax and logmax else float("nan"),
            "null_p99": null99,
        })
    df = pd.DataFrame(rows).sort_values("abs_delta_single", ascending=False)

    missing = set(PANEL_DELTA) | {"PSMC4", "SMCR8", "TPM2"}
    absent = sorted(missing - set(df.gene))
    if absent:
        print(f"  WARN: arms not yet present: {', '.join(absent)}", file=sys.stderr)

    print(f"majority-class baseline on the test split = {baseline:.4f} "
          f"({int(lab.max())}/{int(lab.sum())})\n")
    print(f"{'gene':9}{'class':9}{'acc':>8}{'vs base':>9}{'|D| 1g':>9}{'|D| panel':>11}"
          f"{'|D_in| raw':>11}{'log_max':>10}{'|D_in|/lm':>11}{'null p99':>10}{'clears':>8}")
    print("-" * 105)
    for _, r in df.iterrows():
        cls = r["class"] + ("*" if r.gene in NOT_PIGMENTATION else "")
        cl = "-" if r.null_p99 != r.null_p99 else ("yes" if r.abs_delta_single > r.null_p99 else "no")
        n99 = "        -" if r.null_p99 != r.null_p99 else f"{r.null_p99:9.4f}"
        print(f"{r.gene:9}{cls:9}{r.acc:8.4f}{r.acc_over_baseline:+9.4f}"
              f"{r.abs_delta_single:9.3f}{r.abs_delta_panel:11.3f}{r.delta_in:11.0f}"
              f"{r.log_max_median:10.2e}{r.delta_in_normalised:11.0f}{n99}{cl:>8}")

    pan = df[df["class"] == "panel"]
    ctl = df[df["class"] == "control"]
    pan_pig = pan[~pan.gene.isin(NOT_PIGMENTATION)]
    # Both are reported, always, and the pre-registered one first. Reporting only the
    # focused set would hide that the exclusion was decided after the measurement.
    for label, P in [(f"pre-registered  panel (n={len(pan)}) vs control (n={len(ctl)})", pan),
                     (f"pigmentation    panel (n={len(pan_pig)}) vs control (n={len(ctl)})"
                      f"  [-{', '.join(sorted(NOT_PIGMENTATION))}]", pan_pig)]:
        print(f"\n--- {label} ---")
        for col in ["acc", "abs_delta_single"]:
            u = mannwhitneyu(P[col], ctl[col], alternative="greater")
            print(f"  {col:18} panel median={P[col].median():.3f}  control median={ctl[col].median():.3f}"
                  f"   Mann-Whitney U={u.statistic:.0f} p={u.pvalue:.4f} (one-sided, panel>control)")

    print("\n--- does the single-gene response still track delivered magnitude? ---")
    r1 = spearmanr(df.abs_delta_single, df.delta_in)
    r1n = spearmanr(df.abs_delta_single, df.delta_in_normalised)
    r2 = spearmanr(df.abs_delta_panel.dropna(),
                   df.loc[df.abs_delta_panel.notna(), "delta_in"])
    print(f"  single-gene |D| vs |D_in| RAW        : rho={r1.statistic:+.3f} p={r1.pvalue:.4g}  (n={len(df)})")
    print(f"  single-gene |D| vs |D_in|/log_max    : rho={r1n.statistic:+.3f} p={r1n.pvalue:.4g}"
          f"  <- the scale the network sees")
    print(f"  panel     |D| vs |D_in| : rho={r2.statistic:+.3f} p={r2.pvalue:.3f}  (published limitation)")

    print("\n--- delivery-matched pairs, matched on |D_in|/log_max (the network's scale) ---")
    pairs = []
    for _, c in ctl.iterrows():
        j = (pan.delta_in_normalised - c.delta_in_normalised).abs().idxmin()
        p = pan.loc[j]
        ratio_del = c.delta_in_normalised / p.delta_in_normalised
        ratio_resp = p.abs_delta_single / c.abs_delta_single if c.abs_delta_single else float("inf")
        pairs.append({"control": c.gene, "panel": p.gene, "delivery_ratio": ratio_del,
                      "response_ratio": ratio_resp})
        print(f"  {c.gene:9}{c.delta_in_normalised:9.0f} -> {c.abs_delta_single:.3f}   vs   "
              f"{p.gene:9}{p.delta_in_normalised:9.0f} -> {p.abs_delta_single:.3f}   "
              f"delivery {ratio_del:.2f}x, panel/control response {ratio_resp:.2f}x")

    OUT.write_text(json.dumps({
        "baseline": baseline, "arms": df.to_dict("records"), "delivery_matched_pairs": pairs,
        "not_pigmentation_excluded": sorted(NOT_PIGMENTATION),
        "mannwhitney": {c: {"U": float(mannwhitneyu(pan[c], ctl[c], alternative="greater").statistic),
                            "p": float(mannwhitneyu(pan[c], ctl[c], alternative="greater").pvalue)}
                        for c in ["acc", "abs_delta_single"]},
        "mannwhitney_pigmentation_only": {
            c: {"U": float(mannwhitneyu(pan_pig[c], ctl[c], alternative="greater").statistic),
                "p": float(mannwhitneyu(pan_pig[c], ctl[c], alternative="greater").pvalue)}
            for c in ["acc", "abs_delta_single"]},
        "spearman_single_delta_vs_delivery_raw": {"rho": float(r1.statistic), "p": float(r1.pvalue)},
        "spearman_single_delta_vs_delivery_normalised": {"rho": float(r1n.statistic),
                                                         "p": float(r1n.pvalue)},
        "arms_absent": absent,
    }, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
