#!/usr/bin/env python3
"""E20 -- does DITA's own coordinate frame leak ancestry too?

E19 (`indel_phase_ancestry.py`) showed the UNALIGNED frame is ancestry-informative: its fixed
crop holds a different reference span in every individual, displaced by their cumulative upstream
indel load, and 88 length scalars derived from that reach 0.861 on superpopulation.

It is tempting to read DITA's 14.5-point accuracy deficit (0.646 vs 0.791, tracks-only
superpopulation) as the price of removing that channel. That inference is only valid if DITA does
not substitute a comparably informative channel of its own -- and it plainly substitutes
SOMETHING. In the tracks-only configuration the event channels M_i are dropped, but X_i is still
zero wherever pi_i(u) is undefined, so the pattern of unoccupied units on the expanded axis
remains an implicit trace of the individual's deletions and of how much shorter their insertions
are than the cohort maximum at each anchor. That pattern is exactly the `valid_mask` that
`indel_tensor_builder.build_aligned_haplotype_tensor` computes.

This script measures how ancestry-informative that substitute is, using the SAME design as E19 so
the two numbers are directly comparable: the same family-aware split, the same logistic
regression, the same metrics, and the same feature budget (4 per gene per haplotype = 88, plus a
restricted 22-feature variant). Features are derived from the aligner's own per-sample entries --
no AlphaGenome track is ever loaded, and no allele identity or dosage is read:

  n_unocc    unoccupied units in the crop the classifier reads
  n_runs     number of contiguous unoccupied runs (roughly, distinct deletion/slack events)
  max_run    longest unoccupied run
  centroid   mean normalised position of the unoccupied units (spatial asymmetry)

Reading the outcome:
  ~0.86  DITA swapped one genotype channel for another of similar strength; the accuracy deficit
         is NOT explained by leak removal and that reading must be dropped from the paper.
  <<0.86 DITA genuinely reduces the leak, and the accuracy reading becomes defensible.

Either way the knockdown is unaffected: a scramble leaves the occupancy pattern bit-identical, so
it cancels in Delta.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/dita_occupancy_ancestry.py --task superpopulation
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
FEATURES = ("n_unocc", "n_runs", "max_run", "centroid")


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


AXIS_REL = "alignment_cache/dynamic_indel_ref_window_v4"


def load_axis(dataset_dir, sample_set_key, gene):
    """Read the cached cohort axis directly.

    We bypass DynamicIndelAligner's fingerprint validation deliberately: this is a read-only
    analysis, the axis content is what the pipeline uses, and the only mismatching field is
    `center_window_size` (cached as null), which does not affect the axis itself.
    """
    p = Path(dataset_dir) / AXIS_REL / sample_set_key / gene / "axis.json"
    a = json.loads(p.read_text())["axis"]
    return {
        "expanded_length": int(a["expanded_length"]),
        "ref_length": int(a["ref_length"]),
        # JSON keys are strings; convert as _json_key_dict_to_int does.
        "map": {int(k): int(v) for k, v in a["expanded_index_map"].items()},
        # anchor ref index -> explicit list of expanded-axis units allocated as insertion slots
        "slots": {int(k): [int(x) for x in v] for k, v in a.get("insertion_slots_by_ref", {}).items()},
    }


def haplotype_occupancy(vcf_path, window_start, ref_start_offset, axis, hap_idx, slice_bounds):
    """Unoccupied-unit statistics for one haplotype over the cropped expanded axis.

    A unit is occupied iff the haplotype supplies a base for it: every reference position it
    retains occupies its mapped unit, and an insertion of length L at an anchor occupies L of
    the slots allocated there. Unoccupied units are therefore this haplotype's deletions plus
    the insertion slack it does not use -- exactly the zeros `valid_mask` carries.
    """
    ref_len = axis["ref_length"]
    amap, slots = axis["map"], axis["slots"]
    lo, hi = slice_bounds

    deleted_ref = set()
    used_slots = {}
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            r0 = int(c[1]) - window_start - ref_start_offset   # 0-based within the axis ref span
            if r0 < -1000 or r0 >= ref_len + 1000:
                continue
            ref, alts = c[3], c[4].split(",")
            gt = c[9].split(":")[0]
            sep = "|" if "|" in gt else "/"
            al = gt.split(sep)
            if len(al) != 2:
                continue
            a = al[hap_idx]
            if a in ("0", "."):
                continue
            n = int(a)
            if not (1 <= n <= len(alts)):
                continue
            alt = alts[n - 1]
            d = len(alt) - len(ref)
            if d < 0:                                  # deletion: bases after the kept prefix
                for k in range(len(alt), len(ref)):
                    rp = r0 + k
                    if 0 <= rp < ref_len:
                        deleted_ref.add(rp)
            elif d > 0:                                # insertion: consumes slots at the anchor
                if 0 <= r0 < ref_len:
                    used_slots[r0] = max(used_slots.get(r0, 0), min(d, len(slots.get(r0, ()))))

    unocc = []
    for rp in deleted_ref:
        u = amap.get(rp)
        if u is not None and lo <= u < hi:
            unocc.append(u - lo)
    for rp, slot_units in slots.items():
        for uu in slot_units[used_slots.get(rp, 0):]:
            if lo <= uu < hi:
                unocc.append(uu - lo)

    local_length = hi - lo
    n = len(unocc)
    if n == 0:
        return {"n_unocc": 0.0, "n_runs": 0.0, "max_run": 0.0, "centroid": 0.0}
    idx = np.sort(np.asarray(unocc, dtype=np.int64))
    brk = np.flatnonzero(np.diff(idx) > 1)
    runs = np.diff(np.r_[0, brk + 1, len(idx)])
    return {
        "n_unocc": float(n),
        "n_runs": float(len(runs)),
        "max_run": float(runs.max()),
        "centroid": float(idx.mean() / max(local_length, 1)),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", choices=sorted(CONFIGS), required=True)
    ap.add_argument("--limit-individuals", type=int, default=None)
    ap.add_argument("--allow-axis-build", action="store_true",
                    help="permit building a missing cohort axis (SLOW -- hours). Off by default.")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix
    from sklearn.preprocessing import StandardScaler

    from genomics.predictors.genotype_based.alignment.dynamic_indel_alignment import DynamicIndelAligner
    from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
    from genomics.predictors.genotype_based.data.pipeline import _load_split_index

    cfg = load_config(REPO_ROOT / CONFIGS[args.task])
    di = cfg.dataset_input
    genes = list(di.genes_to_use)
    dataset_dir = Path(di.dataset_dir)
    window_center_size = int(di.window_center_size)
    split_index = _load_split_index(get_dataset_cache_dir(cfg))
    meta = json.loads((dataset_dir / "dataset_metadata.json").read_text())
    ped = meta.get("individuals_pedigree", {})

    splits = {k: list(split_index[k]) for k in ("train", "val", "test")}
    all_ids_full = [s for k in ("train", "val", "test") for s in splits[k]]

    # The aligner's cache key is a function of the FULL selected sample set, so it must be built
    # from the complete cohort even when we only evaluate a subset.
    aligner = DynamicIndelAligner(dataset_dir, selected_sample_ids=set(all_ids_full),
                                  center_window_size=window_center_size)
    key = aligner._sample_set_key()
    axis_root = Path(dataset_dir) / "alignment_cache" / "dynamic_indel_ref_window_v4" / key
    _log(f"task={args.task}  sample-set key = {key}")
    if not axis_root.exists() and not args.allow_axis_build:
        _log(f"ABORT: no cached cohort axis at {axis_root}")
        _log("Building it is a multi-hour job. Re-run with --allow-axis-build to permit it.")
        return 2
    _log(f"cohort axis cache present: {axis_root}")

    if args.limit_individuals:
        for k in splits:
            splits[k] = splits[k][: max(4, args.limit_individuals // 3)]
    all_ids = [s for k in ("train", "val", "test") for s in splits[k]]
    _log(f"train/val/test = {len(splits['train'])}/{len(splits['val'])}/{len(splits['test'])}")

    # Per gene: the cached cohort axis and the crop the classifier actually reads.
    #
    # NOTE: get_reference_centered_expanded_slice looks up expanded_index_map with an int key.
    # When the axis cache fails its fingerprint check (it does here -- the cache records
    # center_window_size=null against the config's 32768) the map retains string keys, the
    # lookup silently falls back to its default, and the crop starts at 0 instead of being
    # reference-centred. We reproduce that behaviour rather than correct it, because it is what
    # the trained model consumed.
    axes, slices, ref_off, win_start = {}, {}, {}, {}
    for g in genes:
        axes[g] = load_axis(dataset_dir, key, g)
        info = aligner.get_reference_centered_expanded_slice(g, window_center_size)
        slices[g] = (int(info["expanded_start"]), int(info["expanded_end"]))
        wm = json.loads((dataset_dir / "references" / "windows" / g / "window_metadata.json").read_text())
        win_start[g] = int(wm["start"])
        full = int(wm["end"]) - int(wm["start"]) + 1
        ref_off[g] = full // 2 - axes[g]["ref_length"] // 2
        _log(f"  {g}: axis ref={axes[g]['ref_length']} expanded={axes[g]['expanded_length']} "
             f"(+{axes[g]['expanded_length']-axes[g]['ref_length']}), crop {slices[g]}, "
             f"ref_start_offset={ref_off[g]}")

    cols = [f"{g}_h{h}_{f}" for g in genes for h in (1, 2) for f in FEATURES]
    per_ind, missing = {}, 0
    t0 = time.monotonic()
    for i, sid in enumerate(all_ids, 1):
        vals = []
        for g in genes:
            vcf = dataset_dir / "individuals" / sid / "windows" / g / f"{sid}.window.vcf.gz"
            if not vcf.exists():
                missing += 1
                vals.extend([0.0] * (2 * len(FEATURES)))
                continue
            for hap in (0, 1):
                f = haplotype_occupancy(vcf, win_start[g], ref_off[g], axes[g], hap, slices[g])
                vals.extend(f[k] for k in FEATURES)
        per_ind[sid] = np.asarray(vals, dtype=np.float32)
        if i % 200 == 0:
            _log(f"  {i}/{len(all_ids)} individuals ({i/(time.monotonic()-t0):.1f}/s)")
    _log(f"done in {time.monotonic()-t0:.0f}s; {missing} missing window VCFs; "
         f"{len(cols)} occupancy features")

    def matrix(ids): return np.vstack([per_ind[s] for s in ids])

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

    unocc_idx = [j for j, c in enumerate(cols) if c.endswith("_n_unocc")]
    load = Xte[:, unocc_idx].sum(axis=1)
    class_names = (["strong", "weak"] if args.task == "pigmentation"
                   else sorted({ped[s]["superpopulation"] for s in all_ids}))
    by_class = {}
    for c, name in enumerate(class_names):
        m = yte == c
        if m.any():
            by_class[name] = {"n": int(m.sum()), "mean_unoccupied_units": float(load[m].mean())}
    _log("mean unoccupied units per individual over 11 windows, by class: "
         + ", ".join(f"{k}={v['mean_unoccupied_units']:.0f}" for k, v in by_class.items()))

    results = {}
    for mode, keep in (("occupancy_all", None), ("n_unocc_only", "n_unocc")):
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
        _log(f"=== DITA occupancy LR ({mode}, {len(sel)} feats) === C={Creg} val={va:.4f} | "
             f"TEST acc={r['accuracy']:.4f} bal={r['balanced_accuracy']:.4f} "
             f"macroF1={r['macro_f1']:.4f}")

    majority = float(max(np.bincount(yte)) / len(yte))
    _log(f"majority-class baseline on test: {majority:.4f}")

    out = args.out or (REPO_ROOT / "results" / "genotype_based_predictor" /
                       f"dita_occupancy_ancestry_{args.task}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "task": args.task,
        "sample_set_key": key,
        "window_center_size": window_center_size,
        "genes": genes,
        "features_per_gene_haplotype": list(FEATURES),
        "n_features": len(cols),
        "n_train": len(splits["train"]), "n_val": len(splits["val"]), "n_test": len(splits["test"]),
        "missing_window_vcfs": missing,
        "majority_class_baseline": majority,
        "unoccupied_by_class": by_class,
        "class_names": list(class_names),
        "results": results,
        "generated_at": _now(),
    }, indent=2))
    _log(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
