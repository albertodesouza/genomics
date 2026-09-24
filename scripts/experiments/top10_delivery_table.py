#!/usr/bin/env python3
"""Top-10 knockdown: the response next to how much the intervention actually delivered.

Three different questions, one row each:

  RESPOSTA   delta (log-odds), sd, d/sd, unif -- how the classifier moved
  ENTREGA    delta_in          -- perturbation that reached the CNN input, measured on the
                                  NORMALIZED tensor. The raw-unit version is not comparable
                                  between genes: normalization divides by each gene's own
                                  log_max, so a near-silent window amplifies a tiny raw change.
  EXPRESSAO  log2fc, delta_expr_rel -- whether AlphaGenome's predicted expression of the gene
                                  actually fell. log2fc is signed and is the only column that
                                  can say the intervention was a knockdown rather than merely
                                  a change; delta_expr on its own is a magnitude.

`ratio` = |delta| / delta_in is the efficiency the paper reports: response per unit of
delivered perturbation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

RUNS = Path("/home/breno/I2CA/genomics/results/genotype_based_predictor")
PANEL = {"SLC24A5", "SLC45A2", "HERC2", "TYRP1", "MC1R", "DDB1", "MFSD12"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", type=Path,
                    default=RUNS / "knockout_bulk/top10/top10_knockdown.csv")
    ap.add_argument("--json-out", default=str(RUNS / "top10_delivery_table.json"))
    a = ap.parse_args()

    d = pd.read_csv(a.in_path)
    need = {"expr_log2fc", "delta_expr_signed"}
    if not need <= set(d.columns):
        raise SystemExit(f"ABORT: {a.in_path} nao tem {sorted(need - set(d.columns))}; "
                         f"rode o knockdown novamente com a versao assinada.")

    rows = []
    for g, s in d.groupby("gene"):
        lo = s.delta_log_odds
        sd = float(lo.std(ddof=1))
        fp = float((lo > 0).mean())
        din = float(s.delta_in_norm.mean())
        rows.append({
            "gene": g,
            "classe": "painel" if g in PANEL else "controle",
            "n": len(s),
            "delta": float(lo.mean()),
            "delta_sd": sd,
            "d_sd": abs(float(lo.mean())) / sd if sd > 0 else float("nan"),
            "unif": max(fp, 1 - fp),
            "delta_in": din,
            "delta_in_raw": float(s.delta_in.mean()),
            "ratio": abs(float(lo.mean())) / din if din > 0 else float("nan"),
            "delta_expr": float(s.delta_expr.mean()),
            "delta_expr_rel": float(s.delta_expr_rel.mean()),
            "expr_log2fc": float(s.expr_log2fc.mean()),
            "frac_expr_down": float((s.delta_expr_signed < 0).mean()),
            "flip": float(s.flipped.mean()),
        })
    rows.sort(key=lambda r: -abs(r["delta"]))

    hdr = (f"{'gene':9s} {'classe':8s} {'delta':>9s} {'d/sd':>5s} {'unif':>5s} | "
           f"{'delta_in':>9s} {'ratio':>9s} | {'log2fc':>7s} {'%down':>6s} {'expr_rel':>8s}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['gene']:9s} {r['classe']:8s} {r['delta']:+9.4f} {r['d_sd']:5.2f} "
              f"{r['unif']:5.3f} | {r['delta_in']:9.1f} {r['ratio']:9.2e} | "
              f"{r['expr_log2fc']:+7.3f} {r['frac_expr_down']:6.1%} {r['delta_expr_rel']:8.4f}")
    print("-" * len(hdr))
    print("delta_in = perturbacao no tensor normalizado (comparavel entre genes)")
    print("ratio    = |delta| / delta_in : resposta por unidade entregue")
    print("log2fc   = log2(expressao perturbada / basal) sobre exons, ASSINADO")
    print("%down    = fracao de individuos em que a expressao predita CAIU")
    print("expr_rel = soma de |mudanca| sobre exons, relativa a basal (magnitude)")

    Path(a.json_out).write_text(json.dumps(rows, indent=1))
    print(f"\njson -> {a.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
