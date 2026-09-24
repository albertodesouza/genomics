#!/usr/bin/env python3
"""Final table for the max-pool arms: classifier quality and knockdown response per gene.

Joins two measurements that answer different questions and must not be confused:
  bal_acc / acc / rec      -- how much the gene's window says about the label (test split)
  bal_acc_val              -- the same on validation, which is the split the ten-gene
                              selector must rank on: choosing the ten on test would pick
                              them on the split their accuracy is then reported on
  delta / delta_in / delta_expr -- what happens when its promoter is scrambled

`delta` follows the convention already in readout_normalisation.json: the MEAN SIGNED change in
(strong - weak) log-odds over the 162 test individuals, not the mean absolute change. The sign is
the direction the classifier moves, and averaging absolute values would turn per-individual noise
into an apparent response, so the sd is reported next to it.

Class 0 is 'strong pigmentation' (support 115) and class 1 is 'weak pigmentation' (support 47),
verified from per_class_metrics rather than assumed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/home/breno/I2CA/genomics")
RUNS = REPO / "results/genotype_based_predictor"
KD = RUNS / "knockout_bulk/poolmax"

PANEL = ["SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"]
CONTROL = ["SPRED2", "FRA10AC1", "PSMC4", "LRRC36", "PPP1R3E", "PRSS55",
           "EIF1B", "SMCR8", "LACTB2", "ECHDC3", "TPM2"]
# random_11_2, the second pre-existing random draw, taken whole. Both draws predate the
# specificity question and are independent of it, which is the only thing making the
# comparison interpretable; a draw is never reordered, substituted or trimmed.
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]
# The third pre-existing draw, built 2026-09-12 by control_expansion_11c.sh. Behind a
# flag and OFF by default: the paper reports 31 arms, and a plain run of this script
# regenerates the table the paper reads, so switching the default would substitute a
# 42-arm table into the paper as a side effect of running the script.
DRAW3 = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6", "KIAA0319", "RIDA",
         "SEM1", "SFMBT2", "SUMF2"]
CONTROL = CONTROL + DRAW2



def split_metrics(gene: str, split: str) -> dict:
    """Balanced accuracy and recalls for one arm on one split.

    Recomputed from per_class_metrics rather than read off a summary field, because the
    checkpoint's own reported accuracy is unbalanced and the classes are 115/47 on test
    and 108/41 on validation.
    """
    hits = sorted(RUNS.glob(
        f"runs_poolmax*/{gene.lower()}_poolmax*/*/{split}_best_accuracy_results.json"))
    if len(hits) != 1:
        raise SystemExit(f"ABORT: {gene} has {len(hits)} {split} result files")
    d = json.loads(hits[0].read_text())
    pcm = d["per_class_metrics"]
    if "strong pigmentation" not in pcm or "weak pigmentation" not in pcm:
        raise SystemExit(f"ABORT: {gene} per_class_metrics keys are {list(pcm)}")
    cm = np.asarray(d["confusion_matrix"], dtype=float)
    rec_strong = float(pcm["strong pigmentation"]["recall"])
    rec_weak = float(pcm["weak pigmentation"]["recall"])
    return {
        "bal_acc": (rec_strong + rec_weak) / 2.0,
        "acc": float(cm.trace() / cm.sum()),
        "rec_weak": rec_weak,
        "rec_strong": rec_strong,
    }


def clf_metrics(gene: str) -> dict:
    """Test-split metrics under the historical names, plus the validation split.

    The two splits answer different questions and the paper uses them for different
    things: the selector that picks the ten genes for the multi-gene probe must rank on
    VALIDATION, or the ten genes are chosen on the same split their accuracy is later
    reported on. The test columns keep their original names so every figure reading this
    table is unaffected; `bal_acc_val` is additive.
    """
    r = split_metrics(gene, "test")
    v = split_metrics(gene, "val")
    r["bal_acc_val"] = v["bal_acc"]
    r["acc_val"] = v["acc"]
    return r


def kd_metrics(gene: str) -> dict:
    f = KD / f"single_{gene.lower()}.csv"
    if not f.exists():
        return {"n_kd": 0}
    d = pd.read_csv(f)
    if d.empty:
        return {"n_kd": 0}
    return {
        "n_kd": len(d),
        "delta": float(d["delta_log_odds"].mean()),
        "delta_sd": float(d["delta_log_odds"].std(ddof=1)),
        "delta_in": float(d["delta_in"].mean()),
        "delta_in_norm": float(d["delta_in_norm"].mean()),
        "delta_expr": float(d["delta_expr"].mean()),
        "delta_expr_rel": float(d["delta_expr_rel"].mean()),
        "flip_rate": float(d["flipped"].mean()),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-out", default=str(RUNS / "poolmax_final_table.json"))
    ap.add_argument("--csv-out", default=str(RUNS / "poolmax_final_table.csv"))
    ap.add_argument("--include-draw3", action="store_true",
                    help="add the third control draw (11 arms); write to a PREVIEW path, "
                         "not over the table the paper reads")
    a = ap.parse_args()

    controls = CONTROL + (DRAW3 if a.include_draw3 else [])
    if a.include_draw3 and Path(a.csv_out) == RUNS / "poolmax_final_table.csv":
        raise SystemExit("ABORT: --include-draw3 with the default --csv-out would "
                         "overwrite the 31-arm table the paper reads. Pass --csv-out and "
                         "--json-out explicitly.")
    rows = []
    for gene in PANEL + controls:
        r = {"gene": gene, "classe": "painel" if gene in PANEL else "controle"}
        r.update(clf_metrics(gene))
        r.update(kd_metrics(gene))
        rows.append(r)

    df = pd.DataFrame(rows).sort_values("bal_acc", ascending=False).reset_index(drop=True)
    incomplete = df[df["n_kd"] != 162]
    hdr = (f"{'gene':9s} {'classe':8s} {'bal_val':>7s} {'bal_acc':>7s} {'acc':>6s} "
           f"{'rec_wk':>6s} {'rec_st':>6s} "
           f"{'delta':>8s} {'sd':>7s} {'d/sd':>5s} {'delta_in':>10s} {'delta_expr':>11s} {'n':>4s}")
    print(hdr); print("-" * len(hdr))
    for _, r in df.iterrows():
        d = f"{r['delta']:+8.4f}" if r["n_kd"] else f"{'--':>8s}"
        sd = f"{r['delta_sd']:7.3f}" if r["n_kd"] else f"{'--':>7s}"
        rat = f"{abs(r['delta'])/r['delta_sd']:5.2f}" if r["n_kd"] else f"{'--':>5s}"
        di = f"{r['delta_in_norm']:10.1f}" if r["n_kd"] else f"{'--':>10s}"
        de = f"{r['delta_expr']:11.1f}" if r["n_kd"] else f"{'--':>11s}"
        print(f"{r['gene']:9s} {r['classe']:8s} {r['bal_acc_val']:7.4f} "
              f"{r['bal_acc']:7.4f} {r['acc']:6.4f} "
              f"{r['rec_weak']:6.3f} {r['rec_strong']:6.3f} {d} {sd} {rat} {di} {de} {int(r['n_kd']):4d}")
    print("-" * len(hdr))
    print("classe 0 = strong (n=115), classe 1 = weak (n=47); preditor constante = 0.5000 bal / 0.7099 acc")
    print("d/sd = |delta| / sd entre individuos; abaixo de 1 a media nao se distingue de zero.")
    print("delta_in = perturbacao no tensor normalizado (comparavel entre genes);")
    print("delta_in_raw = mesma perturbacao em unidades AlphaGenome cruas (convencao publicada,")
    print("  NAO comparavel entre genes: normalizacao divide pelo log_max do proprio gene).")
    if len(incomplete):
        print(f"INCOMPLETO: {list(incomplete['gene'])} nao tem 162 linhas de knockdown -- "
              f"nao usar essas celulas")

    # The selector diagnostic. Printed every run because the two orderings disagree and
    # the disagreement is the reason the arm was retrained: any claim about "the ten
    # highest-accuracy arms" has to say which split ranked them.
    k = 10
    by_val = list(df.sort_values("bal_acc_val", ascending=False)["gene"].head(k))
    by_test = list(df.sort_values("bal_acc", ascending=False)["gene"].head(k))
    both = [g for g in by_val if g in by_test]
    print(f"\ntop-{k} on VAL  (the selector): {', '.join(by_val)}")
    print(f"top-{k} on TEST (not usable)  : {', '.join(by_test)}")
    print(f"overlap {len(both)}/{k}; val-only {sorted(set(by_val) - set(by_test))}; "
          f"test-only {sorted(set(by_test) - set(by_val))}")
    rv = df["bal_acc_val"].rank(ascending=False)
    rt = df["bal_acc"].rank(ascending=False)
    print(f"Spearman(val, test) over {len(df)} arms = {rv.corr(rt, method='spearman'):+.3f}")
    cut = df.sort_values("bal_acc_val", ascending=False)["bal_acc_val"].to_numpy()
    print(f"val cut at rank {k}: {cut[k - 1]:.4f} vs rank {k + 1}: {cut[k]:.4f} "
          f"(margin {cut[k - 1] - cut[k]:.4f})")

    Path(a.json_out).write_text(json.dumps(rows, indent=1))
    df.to_csv(a.csv_out, index=False)
    print(f"\njson -> {a.json_out}\ncsv  -> {a.csv_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
