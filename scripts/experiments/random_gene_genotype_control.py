#!/usr/bin/env python3
"""E2 (genotype tier) -- random-gene negative control for the E4 genotype baseline.

Refits the E4 logistic regression (results/genotype_based_predictor/genotype_baseline_*.json,
scripts/experiments/genotype_baseline.py) on eleven gene windows with no documented pigmentation
role, instead of the eleven pigmentation-panel genes (MC1R, TYRP1, TYR, SLC45A2, DDB1, EDAR,
MFSD12, OCA2, HERC2, SLC24A5, TCHH). If accuracy stays near saturation on this control panel, the
task is being solved by ancestry-informative variants anywhere in the genome rather than by
pigmentation biology specifically, and Table I's genotype-LR numbers must be read as a statement
about ancestry, not about pigmentation (see paper Section V-C / V-B).

Control panel: ECHDC3, EIF1B, FRA10AC1, LACTB2, LRRC36, PPP1R3E, PRSS55, PSMC4, SMCR8, SPRED2,
TPM2 -- the "random_11_1" panel already drawn for the unrelated non-longevous-dataset project
(/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_random_11_1/), reused here rather than
redrawn. Verified against QuickGO (EBI) immediately before use: none of the 11 has a GO annotation
under "pigmentation" (GO:0043473), "melanin biosynthetic process" (GO:0042438), "melanocyte
differentiation" (GO:0030318), or any GO term whose name contains pigment/melanin/melanocyte.
Caveat kept in view rather than hidden: SPRED2's paralog SPRED1 causes Legius syndrome, which
includes hyperpigmented macules (a MAPK-pathway effect) -- SPRED2 itself carries no such GO
annotation, so it is kept, but this is exactly the kind of case a hostile reviewer would raise, and
this note exists so nobody has to rediscover it.

No AlphaGenome call is involved anywhere in this script. Genotype dosages are read directly from
the raw per-chromosome 1000 Genomes phased VCFs (the same files the main dataset's window VCFs were
themselves cut from, per each real gene's window_metadata.json -> raw_variant_source.vcf_pattern),
via a single bcftools region query per gene covering every individual in the split at once -- there
is no per-individual window VCF for these genes on disk, so genotype_baseline.py's approach (read
one pre-cut mini-VCF per individual per gene) does not apply here and is not reused.

Window definition, matched to how the real gene windows were built
(src/genomics/workflows/dataset_builders/non_longevous/build_window_and_predict.py, which resizes
the GENCODE gene-body interval to a fixed 524,288 bp centered on the gene body's own midpoint): the
524,288 bp outer window and the network's 32,768 bp central crop share the same center by
construction, so this script skips the outer window and goes straight to the crop, computed as
gene_body_midpoint +/- 16,384 bp from the same cached GENCODE v46 GTF
(gtf_cache.feather, already downloaded for the non-longevous-dataset project, reused here rather
than re-fetched). This is arithmetically equivalent to what the real pipeline does, not an
approximation of it.

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/random_gene_genotype_control.py --task pigmentation
  ... --task superpopulation
  ... --limit-individuals 40      # smoke test
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
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
CONTROL_GENES = [
    "ECHDC3", "EIF1B", "FRA10AC1", "LACTB2", "LRRC36",
    "PPP1R3E", "PRSS55", "PSMC4", "SMCR8", "SPRED2", "TPM2",
]
GTF_CACHE = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_random_11_1/gtf_cache.feather")
VCF_PATTERN = ("/dados/GENOMICS_DATA/top3/longevity_dataset/vcf_chromosomes/"
               "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
BCFTOOLS = "/home/breno/miniforge3/envs/genomics/bin/bcftools"


def _now(): return datetime.now(timezone.utc).isoformat()
def _log(m): print(f"[{_now()}] {m}", flush=True)


def control_gene_crops(genes):
    """{gene: (chrom, lo_pos_1based, hi_pos_1based_exclusive)} centered on each gene body's own
    midpoint, from the cached GENCODE v46 GTF -- see module docstring for why this is equivalent
    to resizing the outer 524,288 bp window and then cropping its center."""
    import pandas as pd
    gtf = pd.read_feather(GTF_CACHE)
    genes_df = gtf[gtf["Feature"] == "gene"].drop_duplicates("gene_name").set_index("gene_name")
    out = {}
    for g in genes:
        if g not in genes_df.index:
            raise ValueError(f"gene {g!r} not found in cached GTF ({GTF_CACHE})")
        row = genes_df.loc[g]
        start, end = int(row["Start"]), int(row["End"])
        center = (start + end) // 2
        out[g] = (str(row["Chromosome"]), center - CROP // 2, center + CROP // 2)
    return out


def read_dosages_bcftools(vcf_path, chrom, lo_pos, hi_pos, sample_ids):
    """{sample_id: {(pos, ref, alt): dosage}} for variants in [lo_pos, hi_pos), read with one
    bcftools query call for every sample at once (bcftools preserves the -s order exactly, verified
    empirically against a differing genotype before writing this)."""
    region = f"{chrom}:{lo_pos}-{hi_pos - 1}"
    cmd = [BCFTOOLS, "query", "-r", region, "-s", ",".join(sample_ids),
           "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    per_sample = {sid: {} for sid in sample_ids}
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        pos, ref, alt_field = int(cols[0]), cols[1], cols[2]
        alts = alt_field.split(",")
        gts = cols[3:]
        for sid, gt in zip(sample_ids, gts):
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
                d = per_sample[sid]
                d[key] = d.get(key, 0) + 1
    return per_sample


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", choices=sorted(CONFIGS), required=True)
    ap.add_argument("--limit-individuals", type=int, default=None)
    ap.add_argument("--min-count", type=int, default=5,
                     help="drop variants carried by fewer than this many training individuals")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--genes", default=None,
                    help="Comma-separated gene symbols to fit on, overriding the eleven-gene "
                         "control panel. Any gene present in the cached GENCODE feather works, "
                         "panel genes included -- which is what lets the same code supply both "
                         "arms of the 3-vs-3 magnitude-matched control "
                         "(--genes TPM2,SMCR8,PSMC4 against --genes TYR,MFSD12,MC1R). Default "
                         "keeps the published eleven-gene panel, so existing invocations are "
                         "unchanged.")
    args = ap.parse_args()

    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix
    from sklearn.preprocessing import StandardScaler

    from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
    from genomics.predictors.genotype_based.data.pipeline import _load_split_index

    cfg = load_config(REPO_ROOT / CONFIGS[args.task])
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
    genes = ([g.strip() for g in args.genes.split(",") if g.strip()]
             if args.genes else list(CONTROL_GENES))
    if not genes:
        raise SystemExit("--genes was given but parsed to an empty list")
    _log(f"control genes: {genes}")

    crops = control_gene_crops(genes)
    for g, (chrom, lo, hi) in crops.items():
        _log(f"  {g}: {chrom}:{lo}-{hi} (crop {hi - lo} bp)")

    per_ind = {sid: {} for sid in all_ids}
    t0 = time.monotonic()
    for g, (chrom, lo, hi) in crops.items():
        vcf = Path(VCF_PATTERN.format(chrom=chrom))
        if not vcf.exists():
            raise FileNotFoundError(vcf)
        by_sample = read_dosages_bcftools(vcf, chrom, lo, hi, all_ids)
        n_variants = len({k for d in by_sample.values() for k in d})
        _log(f"  {g}: {n_variants} distinct variant alleles across {len(all_ids)} individuals "
             f"({time.monotonic() - t0:.1f}s elapsed)")
        for sid, d in by_sample.items():
            for k, v in d.items():
                per_ind[sid][(g,) + k] = v

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
        _log(f"=== control genotype LR ({mode}) === C={Creg} val={va:.4f} | TEST acc={r['accuracy']:.4f} "
             f"bal={r['balanced_accuracy']:.4f} macroF1={r['macro_f1']:.4f}")
        _log(f"    CM: {r['confusion_matrix']}")

    out = args.out or REPO_ROOT / f"results/genotype_based_predictor/random_gene_control_{args.task}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "task": args.task, "crop_bp": CROP, "min_count": args.min_count,
        "control_genes": genes,
        "control_gene_crops": {g: list(v) for g, v in crops.items()},
        "n_variants": len(variants),
        "n_train": len(splits["train"]), "n_val": len(splits["val"]), "n_test": len(splits["test"]),
        "results": results,
    }, indent=2))
    _log(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
