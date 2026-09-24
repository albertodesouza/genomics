#!/usr/bin/env python3
"""E19 -- is the unaligned pipeline's coordinate frame itself an ancestry channel?

`analysis.individual_consensus._adjust_to_target_size` forces every personal sequence to exactly
524 288 bp -- truncating on the right when the haplotype is net-inserted, padding from the
reference when it is net-deleted -- so the sequence is LEFT-ANCHORED at the window start and the
crop indices `ProcessedDataset._extract_center_raw` takes are the same for every individual. The
crop bounds therefore do NOT move. What moves is the content: tensor index `u` holds reference
offset `r` where `r + cumdelta(r) = u`, so an individual carrying a net insertion upstream has
everything from that point on displaced leftward by the cumulative offset. The fixed crop covers a
DIFFERENT reference span in every individual, displaced by their upstream indel load. Indel load
is not ancestry-neutral, so an unaligned tracks-only tensor may carry ancestry information in its
geometry alone, with no AlphaGenome value involved.

This script measures that channel directly. For every individual, gene and haplotype it reads the
window VCF and derives ONLY length/geometry quantities -- never an allele identity, never a
dosage, never a track value:

  net_delta        sum of (len(ALT) - len(REF)) over the whole 512 kbp window on that haplotype;
                   also the number of bases truncated (>0) or reference-padded (<0) at the right
                   end by the target-size adjustment
  disp_crop_start  reference displacement of the content at the crop's left edge, i.e. minus the
                   cumulative indel offset upstream of it
  disp_crop_centre the same at the crop's centre index -- the locus the classifier's middle column
                   actually holds
  crop_ref_span    how many reference bp the fixed 32 768-index crop actually covers

and fits the same logistic regression, on the same family-aware split, as the E4 genotype
baseline. A model well above chance means the unaligned tensor has a genotype pathway that DITA
removes by construction -- and that the DITA-vs-no-alignment accuracy gap is at least partly the
paper's own confound argument rather than a defect of the alignment.

Scope, stated explicitly: this establishes that the channel EXISTS and is ancestry-informative.
It does not establish that the trained CNN uses it. That would need the permutation control (E8).

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/indel_phase_ancestry.py --task superpopulation
  ... --task pigmentation
  ... --limit-individuals 60      # smoke test
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
FEATURES = ("net_delta", "disp_crop_start", "disp_crop_centre", "crop_ref_span")


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


def read_haplotype_deltas(vcf_path, window_start):
    """[(ref_offset, delta_h1, delta_h2)] sorted by offset, for indels only.

    delta = len(ALT) - len(REF) for the allele that haplotype carries; 0 for reference or SNV.
    Offsets are 0-based from the window start. Only rows where at least one haplotype has a
    non-zero length change are retained -- SNVs cannot displace anything.
    """
    rows = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            ref, alts = cols[3], cols[4].split(",")
            gt = cols[9].split(":")[0]
            sep = "|" if "|" in gt else "/"
            alleles = gt.split(sep)
            if len(alleles) != 2:
                continue
            d = []
            for a in alleles:
                if a in ("0", "."):
                    d.append(0)
                    continue
                n = int(a)
                if not (1 <= n <= len(alts)):
                    d.append(0)
                    continue
                d.append(len(alts[n - 1]) - len(ref))
            if d[0] == 0 and d[1] == 0:
                continue
            rows.append((int(cols[1]) - window_start, d[0], d[1]))
    rows.sort(key=lambda r: r[0])
    return rows


def ref_offset_at(deltas, u):
    """Reference offset held by personal-sequence index `u`.

    The personal index of reference offset r is r + cumdelta(r), where cumdelta accumulates the
    length changes of every variant strictly upstream of r. Inverting that map at `u` gives the
    locus the tensor actually holds at index `u`.
    """
    cum, prev_r = 0, 0
    for r, d in deltas:
        if prev_r + cum <= u < r + cum:
            return u - cum
        cum += d
        prev_r = r
    return u - cum


def geometry(deltas, full_len):
    """The four geometry scalars for one haplotype, under fixed-length left-anchored sequences."""
    mid = full_len // 2
    crop_lo, crop_hi = mid - CROP // 2, mid + CROP // 2
    lo_ref, hi_ref = ref_offset_at(deltas, crop_lo), ref_offset_at(deltas, crop_hi)
    return {
        "net_delta": float(sum(d for _, d in deltas)),
        "disp_crop_start": float(lo_ref - crop_lo),
        "disp_crop_centre": float(ref_offset_at(deltas, mid) - mid),
        "crop_ref_span": float(hi_ref - lo_ref),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", choices=sorted(CONFIGS), required=True)
    ap.add_argument("--limit-individuals", type=int, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

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
    _log(f"task={args.task}  train/val/test = "
         f"{len(splits['train'])}/{len(splits['val'])}/{len(splits['test'])}")

    win = {}
    for g in genes:
        wm = json.loads((dataset_dir / "references" / "windows" / g / "window_metadata.json").read_text())
        st = int(wm["start"])
        win[g] = (st, int(wm["end"]) - st + 1)
    _log(f"Deriving indel geometry for {len(all_ids)} individuals x {len(genes)} genes "
         f"(length arithmetic only, no allele identities)...")

    cols = [f"{g}_h{h}_{f}" for g in genes for h in (1, 2) for f in FEATURES]
    per_ind, missing = {}, 0
    t0 = time.monotonic()
    for i, sid in enumerate(all_ids, 1):
        vals = []
        for g in genes:
            vcf = dataset_dir / "individuals" / sid / "windows" / g / f"{sid}.window.vcf.gz"
            st, full = win[g]
            if not vcf.exists():
                missing += 1
                vals.extend([0.0] * (2 * len(FEATURES)))
                continue
            rows = read_haplotype_deltas(vcf, st)
            for h in (0, 1):
                gm = geometry([(r, d[h]) for r, *d in rows], full)
                vals.extend(gm[f] for f in FEATURES)
        per_ind[sid] = np.asarray(vals, dtype=np.float32)
        if i % 200 == 0:
            _log(f"  {i}/{len(all_ids)} individuals ({i/(time.monotonic()-t0):.1f}/s)")
    _log(f"done in {time.monotonic()-t0:.0f}s; {missing} missing window VCFs; "
         f"{len(cols)} geometry features")

    def matrix(ids):
        return np.vstack([per_ind[s] for s in ids])

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

    # Descriptive: does indel load itself separate the classes at all?
    net_idx = [j for j, c in enumerate(cols) if c.endswith("_net_delta")]
    load = np.abs(Xte[:, net_idx]).sum(axis=1)
    by_class = {}
    class_names = (["strong", "weak"] if args.task == "pigmentation"
                   else sorted({ped[s]["superpopulation"] for s in all_ids}))
    for c, name in enumerate(class_names):
        m = yte == c
        if m.any():
            by_class[name] = {"n": int(m.sum()), "mean_abs_net_indel_bp": float(load[m].mean())}
    _log(f"mean |net indel| per individual over 11 windows, by class: "
         + ", ".join(f"{k}={v['mean_abs_net_indel_bp']:.0f}bp" for k, v in by_class.items()))

    results = {}
    for mode, keep in (("geometry_all", None), ("displacement_only", "disp_crop_centre")):
        sel = ([j for j, c in enumerate(cols) if c.endswith(f"_{keep}")] if keep
               else list(range(len(cols))))
        sc = StandardScaler().fit(Xtr[:, sel])
        A, B, C = (sc.transform(M[:, sel]) for M in (Xtr, Xva, Xte))
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
            "C": Creg, "val_accuracy": float(va), "n_features": len(sel),
            "accuracy": float(accuracy_score(yte, pred)),
            "balanced_accuracy": float(balanced_accuracy_score(yte, pred)),
            "macro_f1": float(f1_score(yte, pred, average="macro")),
            "confusion_matrix": confusion_matrix(yte, pred).tolist(),
        }
        r = results[mode]
        _log(f"=== indel geometry LR ({mode}, {len(sel)} feats) === C={Creg} val={va:.4f} | "
             f"TEST acc={r['accuracy']:.4f} bal={r['balanced_accuracy']:.4f} "
             f"macroF1={r['macro_f1']:.4f}")

    majority = float(max(np.bincount(yte)) / len(yte))
    _log(f"majority-class baseline on test: {majority:.4f}")

    out = args.out or (REPO_ROOT / "results" / "genotype_based_predictor" /
                       f"indel_phase_ancestry_{args.task}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "task": args.task,
        "crop_bp": CROP,
        "genes": genes,
        "features_per_gene_haplotype": list(FEATURES),
        "n_features": len(cols),
        "n_train": len(splits["train"]), "n_val": len(splits["val"]), "n_test": len(splits["test"]),
        "missing_window_vcfs": missing,
        "majority_class_baseline": majority,
        "indel_load_by_class": by_class,
        "class_names": list(class_names),
        "results": results,
        "generated_at": _now(),
    }, indent=2))
    _log(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
