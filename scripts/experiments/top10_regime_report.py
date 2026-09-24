#!/usr/bin/env python3
"""Did the perturbation return to the local regime in the 10-gene model?

The single-gene arms broke a property the published panel had: every individual moved the same
way. Measured cause was that scrambling the only channel does not perturb the input, it erases
it, so the output collapses to a constant L* and delta_i = L* - L_i -- whose sign is just the
negative of each individual's own margin, and therefore inverts between AFR and EUR.

Three signatures separate the two regimes, and this reports all three per gene:

  slope of (delta ~ baseline margin)   ~0 local          ~-1 erasure
  sd(L_pert) / sd(L_base)              ~1 local          ~0  erasure
  fraction moving the same direction   ~1 local          ~0.7/0.3 split by class under erasure

The last one is reported alongside the AFR and EUR means, because a pooled fraction near 0.7 is
what class imbalance (115 strong / 47 weak) produces from a pure erasure, and must not be read
as partial uniformity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress

REPO = Path("/home/breno/I2CA/genomics")
RUNS = REPO / "results/genotype_based_predictor"
PANEL_GENES = {"SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12", "TYR", "OCA2"}


def regime_rows(d: pd.DataFrame) -> list[dict]:
    d = d.copy()
    d["lb"] = d.baseline_strong_logit - d.baseline_weak_logit
    d["lp"] = d.perturbed_strong_logit - d.perturbed_weak_logit
    out = []
    for g, s in d.groupby("gene"):
        lr = linregress(s.lb, s.delta_log_odds)
        a = s[s.superpopulation == "AFR"]
        e = s[s.superpopulation == "EUR"]
        frac_pos = float((s.delta_log_odds > 0).mean())
        out.append({
            "gene": g,
            "classe": "painel" if g in PANEL_GENES else "controle",
            "n": len(s),
            "delta": float(s.delta_log_odds.mean()),
            "delta_sd": float(s.delta_log_odds.std(ddof=1)),
            "delta_afr": float(a.delta_log_odds.mean()) if len(a) else float("nan"),
            "delta_eur": float(e.delta_log_odds.mean()) if len(e) else float("nan"),
            "mesmo_sinal": bool(len(a) and len(e)
                                and np.sign(a.delta_log_odds.mean()) == np.sign(e.delta_log_odds.mean())),
            "frac_pos": frac_pos,
            "uniformidade": max(frac_pos, 1 - frac_pos),
            "slope": float(lr.slope),
            "r2": float(lr.rvalue ** 2),
            "sd_ratio": float(s.lp.std() / s.lb.std()) if s.lb.std() > 0 else float("nan"),
            "delta_in": float(s.delta_in_norm.mean()),
            "delta_expr": float(s.delta_expr.mean()),
            "flip_rate": float(s.flipped.mean()),
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", type=Path,
                    default=RUNS / "knockout_bulk/top10/top10_knockdown.csv")
    ap.add_argument("--json-out", default=str(RUNS / "top10_regime.json"))
    a = ap.parse_args()

    d = pd.read_csv(a.in_path)
    rows = sorted(regime_rows(d), key=lambda r: -abs(r["delta"]))

    hdr = (f"{'gene':9s} {'classe':8s} {'delta':>9s} {'AFR':>9s} {'EUR':>9s} {'unif':>5s} "
           f"{'slope':>7s} {'r2':>5s} {'sd_rat':>6s} {'regime':>9s}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        # erasure = the output stopped depending on the individual
        if r["sd_ratio"] < 0.25 or r["slope"] < -0.75:
            regime = "APAGAM."
        elif r["sd_ratio"] > 0.7 and abs(r["slope"]) < 0.35:
            regime = "local"
        else:
            regime = "misto"
        print(f"{r['gene']:9s} {r['classe']:8s} {r['delta']:+9.4f} {r['delta_afr']:+9.4f} "
              f"{r['delta_eur']:+9.4f} {r['uniformidade']:5.3f} {r['slope']:+7.3f} "
              f"{r['r2']:5.3f} {r['sd_ratio']:6.3f} {regime:>9s}")
    print("-" * len(hdr))
    n_same = sum(r["mesmo_sinal"] for r in rows)
    n_unif = sum(r["uniformidade"] > 0.95 for r in rows)
    print(f"mesmo sinal em AFR e EUR: {n_same}/{len(rows)}   |   "
          f"uniformidade > 0.95: {n_unif}/{len(rows)}")
    print("unif = fracao movida na direcao majoritaria (1.000 = todos os individuos juntos)")
    print("regime: 'local' = perturbacao (painel publicado), 'APAGAM.' = entrada apagada")

    Path(a.json_out).write_text(json.dumps(rows, indent=1))
    print(f"\njson -> {a.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
