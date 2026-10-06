#!/usr/bin/env python3
"""Is the panel-vs-control comparison robust to an expression floor, or to any one gene?

An expression floor is the right inclusion policy for the pipeline -- it is what should
run before spending 2,144 AlphaGenome calls building a gene's windows, and it turns the
paper's post-hoc withdrawal of EDAR and TCHH into a criterion declared in advance. This
script asks the separate question of whether applying it changes what the comparison
reports, and the answer is that it does not, for two reasons this script exists to make
checkable.

FIRST, the p-value is non-monotonic in the threshold. A threshold chosen after seeing
the outcome could therefore be made to produce almost anything, so the sweep is printed
in full rather than at one chosen point.

SECOND, the separation is carried by a single gene. Dropping LRRC36 alone moves the
comparison further than any floor does; no other gene comes close. A result that
survives only the removal of its most inconvenient observation is not a result, and
printing the leave-one-out column is the honest way to say so.

The floor itself is read from expression_floor_calibration.json, which anchors it to the
distribution of crop_total_signal over random genomic windows rather than to anything in
this comparison. Two anchors are reported. The p95 anchor is REJECTED: the calibration
shows random 524 kbp windows frequently contain expressed genes (18% exceed 5,000), so
crop_total_signal measures the neighbourhood rather than the target gene, and a p95 floor
excludes five of eleven panel genes including SLC24A5. Only the weak form survives -- a
window below the median random window is not distinguishable from arbitrary genome -- and
it is reported as a necessary, not sufficient, condition.

Reads only JSON already on disk. No AlphaGenome calls, no training.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REPO_ROOT = Path("/home/breno/I2CA/genomics")
os.chdir(REPO_ROOT)
RES = REPO_ROOT / "results/genotype_based_predictor"
SWEEP = RES / "single_gene_sweep_report.json"
PREFLIGHT = RES / "specificity_control_preflight.json"
CALIB = RES / "expression_floor_calibration.json"
OUT = RES / "single_gene_floor_sensitivity.json"
MIN_PER_GROUP = 4      # below this the test reports noise; the sweep says so rather than a p


def mw(d):
    p_, c_ = d[d["class"] == "panel"], d[d["class"] == "control"]
    if len(p_) < MIN_PER_GROUP or len(c_) < MIN_PER_GROUP:
        return len(p_), len(c_), float("nan")
    return len(p_), len(c_), float(
        mannwhitneyu(p_.abs_delta_single, c_.abs_delta_single, alternative="greater").pvalue)


def main():
    for p in (SWEEP, PREFLIGHT):
        if not p.exists():
            raise SystemExit(f"ABORT: {p} missing.")
    rep = json.loads(SWEEP.read_text())
    d = pd.DataFrame(rep["arms"])
    d = d[d["class"].isin(("panel", "control"))].copy()
    pre = json.loads(PREFLIGHT.read_text())["per_gene"]
    d["ref_signal"] = [pre.get(g, {}).get("crop_total_signal", np.nan) for g in d.gene]

    out = {"generated_at": datetime.now(timezone.utc).isoformat()}
    n_p, n_c, p0 = mw(d)
    print(f"baseline, no floor: {n_p} panel vs {n_c} control, p = {p0:.4f}\n")
    out["baseline"] = {"n_panel": n_p, "n_control": n_c, "p": p0}

    # ---- the two calibrated anchors ------------------------------------------------
    if CALIB.exists():
        cal = json.loads(CALIB.read_text())
        v = np.array([x["crop_total_signal"] for x in cal["draws"].values()])
        anchors = {"median random window (weak form, RETAINED)": float(np.median(v)),
                   "p95 random window (REJECTED, see docstring)": float(np.percentile(v, 95))}
        print(f"floor anchored to {len(v)} random genomic windows:")
        out["anchors"] = {}
        for lab, thr in anchors.items():
            k = d[d.ref_signal >= thr]
            a, b, p = mw(k)
            ex = sorted(set(d.gene) - set(k.gene))
            ex_pan = [g for g in ex if d.loc[d.gene == g, "class"].iloc[0] == "panel"]
            print(f"  {lab}")
            print(f"    floor = {thr:.1f}   {a} panel vs {b} control   "
                  f"p = {'  n/a' if p != p else f'{p:.4f}'}")
            print(f"    excludes {len(ex)}: {', '.join(ex) if ex else '(none)'}")
            print(f"    of which panel genes: {', '.join(ex_pan) if ex_pan else '(none)'}")
            out["anchors"][lab] = {"floor": thr, "n_panel": a, "n_control": b, "p": p,
                                   "excluded": ex, "excluded_panel": ex_pan}
        print()
    else:
        print(f"(no {CALIB.name}; anchors skipped)\n")

    # ---- threshold sweep: printed whole, because one point could be shopped ----------
    print("threshold sweep on the reference signal -- note the non-monotonicity:")
    print(f"{'floor':>9}{'panel':>7}{'ctrl':>6}{'med pan':>9}{'med ctl':>9}{'p':>9}")
    sweep = []
    for thr in [0, 50, 100, 145, 250, 500, 1000, 2000, 3000, 5000, 10000]:
        k = d[d.ref_signal >= thr]
        a, b, p = mw(k)
        mp = k[k["class"] == "panel"].abs_delta_single.median()
        mc = k[k["class"] == "control"].abs_delta_single.median()
        sweep.append({"floor": thr, "n_panel": a, "n_control": b, "p": p})
        ps = "     n/a" if p != p else f"{p:9.4f}"
        print(f"{thr:>9}{a:>7}{b:>6}{mp:>9.3f}{mc:>9.3f}{ps}")
    out["sweep"] = sweep
    fin = [s["p"] for s in sweep if s["p"] == s["p"]]
    if len(fin) > 2:
        turns = sum(1 for i in range(1, len(fin) - 1)
                    if (fin[i] - fin[i - 1]) * (fin[i + 1] - fin[i]) < 0)
        print(f"\n  the p-value changes direction {turns} times across the sweep; "
              f"range {min(fin):.4f}-{max(fin):.4f}")
        out["sweep_direction_changes"] = turns

    # ---- leave-one-out: is any single gene carrying the comparison? -----------------
    print("\nleave-one-out, every gene, sorted by how much its removal moves p:")
    loo = []
    for g in d.gene:
        a, b, p = mw(d[d.gene != g])
        loo.append({"gene": g, "class": d.loc[d.gene == g, "class"].iloc[0],
                    "p_without": p, "shift": p - p0})
    loo.sort(key=lambda r: r["p_without"])
    print(f"{'gene':10}{'class':9}{'p without it':>14}{'shift':>9}")
    for r in loo[:5] + [{"gene": "...", "class": "", "p_without": float("nan"),
                         "shift": float("nan")}] + loo[-2:]:
        pw = "   n/a" if r["p_without"] != r["p_without"] else f"{r['p_without']:14.4f}"
        sh = "" if r["shift"] != r["shift"] else f"{r['shift']:+9.4f}"
        print(f"{r['gene']:10}{r['class']:9}{pw}{sh}")
    out["leave_one_out"] = loo

    top = loo[0]
    print(f"\n  removing {top['gene']} alone gives p = {top['p_without']:.4f}, "
          f"against {p0:.4f} with every arm kept.")
    below = [r for r in loo if r["p_without"] < 0.05]
    print(f"  {len(below)} of {len(loo)} single removals put p below 0.05"
          f"{': ' + ', '.join(r['gene'] for r in below) if below else ''}.")
    print("  A separation that requires dropping its most extreme observation is not one.")

    OUT.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
