#!/usr/bin/env python3
"""E4 -- genotype-only baseline over the same eleven windows the CNNs see.

Reads each individual's per-gene window VCF directly, keeps only variants whose reference
position falls inside the SAME central 32,768 bp crop the CNN consumes, and encodes them as
allele dosage (0/1/2). No AlphaGenome value is involved at any point. Trains logistic regression
on the family-aware split of the corresponding predictor config, with and without a PCA
projection, and reports the same metrics as Table I.

This is the baseline the paper needs to close its central argument: the channels-only ablation
shows the frozen representation is unnecessary, but without this we cannot say the whole pipeline
beats counting variants. The published maximum-likelihood SNP classifier is restricted to one
haplotype and to genotyping-array sites, so it is a lower bound, not this.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/genotype_baseline.py --task pigmentation
  ... --task superpopulation
  ... --limit-individuals 40      # smoke test
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _path in (REPO_ROOT / "src", REPO_ROOT / "notebooks"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
os.chdir(REPO_ROOT)

import numpy as np

CONFIGS = {
    "pigmentation": "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml",
    "superpopulation": "configs/predictors/genotype_based/superpopulation/superpopulation.yaml",
}
CROP = 32768
FULL = 524288


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


def crop_bounds(window_start_1based, full_len):
    """The reference span of the central CROP bp of the window -- the region the CNN sees."""
    c = full_len // 2
    lo = c - CROP // 2
    hi = lo + CROP
    return window_start_1based + lo, window_start_1based + hi


def read_dosages(vcf_path, lo_pos, hi_pos):
    """{(pos, ref, alt): dosage} for variants inside [lo_pos, hi_pos)."""
    out = {}
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            pos = int(cols[1])
            if pos < lo_pos or pos >= hi_pos:
                continue
            ref, alts = cols[3], cols[4].split(",")
            gt = cols[9].split(":")[0]
            sep = "|" if "|" in gt else "/"
            alleles = gt.split(sep)
            if len(alleles) != 2:
                continue
            for a in alleles:
                if a in ("0", "."):
                    continue
                n = int(a)
                if not (1 <= n <= len(alts)):
                    continue
                key = (pos, ref, alts[n - 1])
                out[key] = out.get(key, 0) + 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", choices=sorted(CONFIGS), required=True)
    ap.add_argument("--limit-individuals", type=int, default=None)
    ap.add_argument("--min-count", type=int, default=5,
                    help="drop variants carried by fewer than this many training individuals")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix
    from sklearn.preprocessing import StandardScaler

    from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
    from genomics.predictors.genotype_based.data.pipeline import _load_split_index

    cfg = load_config(REPO_ROOT / CONFIGS[args.task])
    genes = list(cfg.dataset_input.genes_to_use)
    dataset_dir = Path(cfg.dataset_input.dataset_dir)
    split_index = _load_split_index(get_dataset_cache_dir(cfg))
    meta = json.loads((dataset_dir / "dataset_metadata.json").read_text())
    ped = meta.get("individuals_pedigree", {})

    splits = {k: list(split_index[k]) for k in ("train", "val", "test")}
    if args.limit_individuals:
        for k in splits:
            splits[k] = splits[k][: max(4, args.limit_individuals // 3)]
    all_ids = [s for k in ("train", "val", "test") for s in splits[k]]
    _log(f"task={args.task}  train/val/test = {len(splits['train'])}/{len(splits['val'])}/{len(splits['test'])}")

    win = {}
    for g in genes:
        wm = json.loads((dataset_dir / "references" / "windows" / g / "window_metadata.json").read_text())
        st = int(wm["start"]); full = int(wm["end"]) - st + 1
        win[g] = crop_bounds(st, full)
    _log(f"Reading window VCFs for {len(all_ids)} individuals x {len(genes)} genes "
         f"(central {CROP} bp only)...")

    per_ind = {}
    t0 = time.monotonic()
    for i, sid in enumerate(all_ids, 1):
        d = {}
        for g in genes:
            vcf = dataset_dir / "individuals" / sid / "windows" / g / f"{sid}.window.vcf.gz"
            if not vcf.exists():
                continue
            lo, hi = win[g]
            for k, v in read_dosages(vcf, lo, hi).items():
                d[(g,) + k] = v
        per_ind[sid] = d
        if i % 200 == 0:
            _log(f"  {i}/{len(all_ids)} individuals ({i/(time.monotonic()-t0):.1f}/s)")

    train_set = set(splits["train"])
    counts = {}
    for sid in splits["train"]:
        for k in per_ind.get(sid, {}):
            counts[k] = counts.get(k, 0) + 1
    variants = sorted(k for k, c in counts.items() if c >= args.min_count)
    vindex = {k: j for j, k in enumerate(variants)}
    _log(f"{len(counts)} variants seen in train; {len(variants)} kept at min_count={args.min_count}")

    def matrix(ids):
        X = np.zeros((len(ids), len(variants)), dtype=np.float32)
        for r, sid in enumerate(ids):
            for k, v in per_ind.get(sid, {}).items():
                j = vindex.get(k)
                if j is not None:
                    X[r, j] = v
        return X

    def labels(ids):
        if args.task == "pigmentation":
            strong = {"YRI", "ESN", "LWK", "MSL", "GWD"}
            return np.array([0 if ped[s]["population"] in strong else 1 for s in ids])
        names = sorted({ped[s]["superpopulation"] for s in all_ids})
        idx = {n: i for i, n in enumerate(names)}
        return np.array([idx[ped[s]["superpopulation"]] for s in ids])

    Xtr, Xva, Xte = (matrix(splits[k]) for k in ("train", "val", "test"))
    ytr, yva, yte = (labels(splits[k]) for k in ("train", "val", "test"))
    _log(f"design matrix: train {Xtr.shape}, val {Xva.shape}, test {Xte.shape}")

    results = {}
    for mode in ("raw", "pca"):
        if mode == "pca":
            n_comp = int(min(300, Xtr.shape[0] - 1, Xtr.shape[1]))
            sc = StandardScaler().fit(Xtr)
            p = PCA(n_components=n_comp, random_state=13).fit(sc.transform(Xtr))
            A, B, C = (p.transform(sc.transform(M)) for M in (Xtr, Xva, Xte))
            _log(f"PCA: {n_comp} components, {p.explained_variance_ratio_.sum():.3f} variance")
        else:
            A, B, C = Xtr, Xva, Xte
        best = None
        for Creg in (0.01, 0.1, 1.0, 10.0):
            clf = LogisticRegression(C=Creg, max_iter=5000, n_jobs=-1)
            clf.fit(A, ytr)
            va = accuracy_score(yva, clf.predict(B))
            if best is None or va > best[0]:
                best = (va, Creg, clf)
        va, Creg, clf = best
        pred = clf.predict(C)
        results[mode] = {
            "C": Creg, "val_accuracy": float(va),
            "n_features": int(A.shape[1]),
            "accuracy": float(accuracy_score(yte, pred)),
            "balanced_accuracy": float(balanced_accuracy_score(yte, pred)),
            "macro_f1": float(f1_score(yte, pred, average="macro")),
            "confusion_matrix": confusion_matrix(yte, pred).tolist(),
        }
        r = results[mode]
        _log(f"=== genotype LR ({mode}) === C={Creg} val={va:.4f} | TEST acc={r['accuracy']:.4f} "
             f"bal={r['balanced_accuracy']:.4f} macroF1={r['macro_f1']:.4f}")
        _log(f"    CM: {r['confusion_matrix']}")

    out = args.out or REPO_ROOT / f"results/genotype_based_predictor/genotype_baseline_{args.task}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "task": args.task, "crop_bp": CROP, "min_count": args.min_count,
        "n_variants": len(variants),
        "n_train": len(splits["train"]), "n_val": len(splits["val"]), "n_test": len(splits["test"]),
        "results": results,
    }, indent=2))
    _log(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
