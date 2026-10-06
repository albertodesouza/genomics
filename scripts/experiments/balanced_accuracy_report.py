#!/usr/bin/env python3
"""Relatorio por braco em acuracia balanceada.

Por que balanceada e nao bruta: o split de teste e 115/47 (71%/29%), entao um
preditor constante da 0.7099 de acuracia bruta -- um numero que parece um
resultado. Em acuracia balanceada o mesmo preditor da exatamente 0.5, e o
colapso fica visivel sem inspecionar a matriz de confusao.

  acuracia balanceada = media dos recalls por classe = (TN/(TN+FP) + TP/(TP+FN)) / 2

Le a matriz de confusao de test_best_accuracy_results.json e a curva de
train_accuracy de models/training_history.json (para a epoca de escape do plato).
"""
import argparse, json
from pathlib import Path

import numpy as np

REPO = Path("/home/breno/I2CA/genomics")
RUNS = REPO / "results/genotype_based_predictor"
# fracao da classe majoritaria no split de TREINO (480/761): valor em que a
# train_accuracy fica presa enquanto a rede preve so a majoritaria
TRAIN_MAJORITY = 480 / 761


def balanced_accuracy(cm) -> float:
    cm = np.asarray(cm, dtype=float)
    recalls = [cm[i, i] / cm[i].sum() for i in range(cm.shape[0]) if cm[i].sum() > 0]
    return float(np.mean(recalls)) if recalls else float("nan")


def scan(run_dir: Path):
    res = list(run_dir.glob("*/test_best_accuracy_results.json"))
    if not res:
        return None
    d = json.loads(res[0].read_text())
    cm = d.get("confusion_matrix")
    if cm is None:
        return None
    out = {
        "run": run_dir.name,
        "acc_raw": float(d.get("weighted_accuracy", float("nan"))),
        "bal_acc": balanced_accuracy(cm),
        "confusion": cm,
        "n": int(d.get("num_samples", 0)),
    }
    cmn = np.asarray(cm, dtype=float)
    out["recall_per_class"] = [float(cmn[i, i] / cmn[i].sum()) if cmn[i].sum() else float("nan")
                              for i in range(cmn.shape[0])]
    hist = list(run_dir.glob("*/models/training_history.json"))
    if hist:
        h = json.loads(hist[0].read_text())
        tr = np.asarray(h.get("train_accuracy", []), dtype=float)
        va = np.asarray(h.get("val_accuracy", []), dtype=float)
        esc = np.flatnonzero(tr > TRAIN_MAJORITY + 0.005)
        out["epochs_run"] = int(len(tr))
        out["escape_epoch"] = int(esc[0]) + 1 if esc.size else None
        out["train_acc_final"] = float(tr[-1]) if tr.size else float("nan")
        out["val_acc_best"] = float(va.max()) if va.size else float("nan")
        out["val_best_epoch"] = int(np.argmax(va)) + 1 if va.size else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="runs_single_gene_*_melstranddita",
                    help="padrao de diretorios de run sob results/genotype_based_predictor")
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args()

    rows = [r for p in sorted(RUNS.glob(a.glob)) if p.is_dir() and (r := scan(p))]
    rows.sort(key=lambda r: (-r["bal_acc"], r["run"]))

    print(f"\n{'braco':<26} {'bal_acc':>8} {'acc_bruta':>10} {'rec_c0':>7} {'rec_c1':>7} "
          f"{'escape':>7} {'ep':>5}  estado")
    print("-" * 92)
    for r in rows:
        esc = r.get("escape_epoch")
        collapsed = abs(r["bal_acc"] - 0.5) < 1e-9 or (r["recall_per_class"][1] == 0.0)
        state = "COLAPSOU (preditor constante)" if collapsed else "ok"
        print(f"{r['run']:<26} {r['bal_acc']:>8.4f} {r['acc_raw']:>10.4f} "
              f"{r['recall_per_class'][0]:>7.3f} {r['recall_per_class'][1]:>7.3f} "
              f"{str(esc) if esc else 'NUNCA':>7} {r.get('epochs_run', 0):>5}  {state}")

    ok = [r for r in rows if abs(r["bal_acc"] - 0.5) >= 1e-9 and r["recall_per_class"][1] > 0]
    print("-" * 92)
    print(f"{len(ok)}/{len(rows)} bracos convergiram  |  "
          f"bal_acc mediana dos convergidos = "
          f"{np.median([r['bal_acc'] for r in ok]):.4f}" if ok else "nenhum braco convergiu")
    print("referencia: preditor constante = 0.5000 balanceada / 0.7099 bruta\n")

    if a.json_out:
        Path(a.json_out).write_text(json.dumps(rows, indent=1))
        print(f"json -> {a.json_out}")


if __name__ == "__main__":
    main()
