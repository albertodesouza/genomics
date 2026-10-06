#!/usr/bin/env python3
"""Assemble the 2x2 of the 3-vs-3 magnitude-matched specificity control.

The test asks whether the *tracks* pipeline is specific to the phenotype panel,
which is the question the predicted-transcriptome ("the representation filters
out variants with no transcriptional effect") argument turns on. Audit test 2 of
the paper answered it only for the dosage logistic regression.

Two 3-gene panels, matched on reference-window crop signal (the proxy that tracks
realised delivered perturbation over the published panel at Spearman rho = 0.855):

    TPM2   48,071  <->  TYR      48,068   (1.00x)
    SMCR8  23,951  <->  MFSD12   26,326   (0.91x)
    PSMC4  18,620  <->  MC1R     14,995   (1.24x)

Reading the table this script prints:

  * The **genotype row** is the control on the control. If both triples saturate,
    each carries essentially all the label information at genotype level, so a
    gap in the tracks row cannot be dismissed as "the control windows are just
    less variable". If the genotype row itself has a gap, the tracks row is
    confounded by it and the comparison is not clean -- say so rather than
    reporting the tracks row alone.
  * The **tracks row** is the measurement. Control approaching panel refutes the
    filter argument; control falling materially below is the first evidence of
    phenotype specificity in this work.

Reads only files already on disk; runs nothing. Cells that have not been produced
yet print as `--`, so this is safe to run while the arms are still training.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
RESULTS = REPO_ROOT / "results" / "genotype_based_predictor"
# Each arm writes to its own tree: the run name is derived from architecture and
# layout only, so both arms -- and the published 11-gene no_alignment run -- would
# otherwise resolve to one directory and clobber each other.
RUNS = RESULTS / "runs_three_vs_three"

TRIPLES = {
    "panel": ["TYR", "MFSD12", "MC1R"],
    "control": ["TPM2", "SMCR8", "PSMC4"],
}
# Signal proxy from specificity_control_preflight.py, for the header line.
CROP_SIGNAL = {
    "TYR": 48068, "MFSD12": 26326, "MC1R": 14995,
    "TPM2": 48071, "SMCR8": 23951, "PSMC4": 18620,
}
LR_JSON = {arm: RESULTS / f"three_vs_three_lr_{arm}_pigmentation.json" for arm in TRIPLES}
TEST_OUTPUT_NAME = "three_vs_three_{arm}_test_results"   # keep in sync with the .sh


def _load(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def lr_cell(arm: str, variant: str = "pca"):
    """(accuracy, balanced accuracy, n_variants) for the dosage-LR arm, or None.

    `random_gene_genotype_control.py` writes `results: {raw: {...}, pca: {...}}`.
    PCA is the published variant of the genotype baseline, so it is the default.
    """
    d = _load(LR_JSON[arm])
    if not d:
        return None
    fit = (d.get("results") or {}).get(variant)
    if not isinstance(fit, dict):
        return None
    return fit.get("accuracy"), fit.get("balanced_accuracy"), d.get("n_variants")


def tracks_cell(arm: str):
    """(accuracy, balanced accuracy, run dir) for the tracks CNN arm, or None."""
    name = TEST_OUTPUT_NAME.format(arm=arm)
    # `genotype test --output-name X` writes X_results.json, so match both spellings.
    hits = (sorted(RUNS.glob(f"*/*/{name}.json"))
            + sorted(RUNS.glob(f"*/*/{name}_results.json"))
            + sorted(RUNS.glob(f"*/{name}.json"))
            + sorted(RUNS.glob(f"*/{name}_results.json")))
    for h in hits:
        d = _load(h)
        if not d:
            continue
        acc = d.get("weighted_accuracy") or d.get("accuracy")
        if acc is None:
            continue
        bal = d.get("balanced_accuracy")
        if bal is None:
            pcm = d.get("per_class_metrics") or {}
            recalls = [v.get("recall") for v in pcm.values() if isinstance(v, dict)]
            recalls = [r for r in recalls if r is not None]
            bal = sum(recalls) / len(recalls) if recalls else None
        return acc, bal, h.parent.parent.name
    return None


def fmt(cell, width=22):
    if cell is None:
        return "--".rjust(width)
    acc, bal, _extra = cell
    a = f"{acc:.4f}" if isinstance(acc, (int, float)) else "?"
    b = f"{bal:.4f}" if isinstance(bal, (int, float)) else "?"
    return f"{a} / {b}".rjust(width)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", type=Path, default=None,
                    help="also write the assembled table as JSON to this path")
    args = ap.parse_args()

    cells = {}
    for arm in TRIPLES:
        cells[("genotype", arm)] = lr_cell(arm)
        cells[("tracks", arm)] = tracks_cell(arm)

    print("3-vs-3 magnitude-matched specificity control -- pigmentation, test split (n=162)")
    print("acc / balanced acc\n")
    for arm, genes in TRIPLES.items():
        total = sum(CROP_SIGNAL[g] for g in genes)
        print(f"  {arm:<8} {', '.join(genes):<24} crop signal {total:>8,}")
    print()
    print(f"  {'representation':<26}{'panel triple':>22}{'control triple':>22}")
    print("  " + "-" * 70)
    for row, label in (("genotype", "dosage LR (PCA)"), ("tracks", "tracks-only CNN")):
        print(f"  {label:<26}{fmt(cells[(row, 'panel')])}{fmt(cells[(row, 'control')])}")
    print()

    geno = [cells[("genotype", a)] for a in TRIPLES]
    tr = [cells[("tracks", a)] for a in TRIPLES]
    if all(g is not None for g in geno):
        gap = abs(geno[0][0] - geno[1][0])
        print(f"  genotype-row gap: {gap:.4f}"
              f"  ({'saturated in both -- tracks row is interpretable' if gap < 0.05 else 'NOT flat: the tracks comparison inherits this gap'})")
    if any(t is None for t in tr):
        missing = [a for a in TRIPLES if cells[("tracks", a)] is None]
        print(f"  tracks row incomplete: {', '.join(missing)} not yet evaluated.")
        print("  Run scripts/experiments/three_vs_three_specificity.sh.")
    else:
        d = tr[0][0] - tr[1][0]
        print(f"  tracks-row gap (panel - control): {d:+.4f}")
        print("  Pre-registered: control approaching panel refutes the filter argument;")
        print("  control materially below is the first evidence of phenotype specificity.")

    if args.json:
        args.json.write_text(json.dumps(
            {f"{row}_{arm}": cells[(row, arm)] for row, arm in cells},
            indent=2), encoding="utf-8")
        print(f"\n  wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
