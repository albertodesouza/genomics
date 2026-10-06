#!/usr/bin/env python3
"""Variant-level gene ranking for panel AND control windows, against the CNN ranking.

This is the arm Section "The same genes ranked by variant-level association" of the
draft specifies and the published GWAS run does not cover: the association scan was
run on the eleven panel windows only, so the panel-against-control comparison it
promises -- the one statistic that puts the two estimators on identical footing -- was
never computable. Here the same plink2 machinery, the same cohort and the same label
are applied to all twenty-two windows, the three gene scores the literature offers are
computed over each, and the resulting ordering is tested head to head against the
ordering the per-gene CNNs induce.

Isolated from the published arm: separate work dir, separate results dir, separate
output JSON. Nothing here writes into results/genotype_based_predictor/gwas/.

Gene scores, all over the same MAF-filtered variants of one 524,288 bp window:

  min         the raw minimum p-value. Reported because the two corrected scores are
              adjustments of it, not as a gene score: the smallest of m uniforms has
              expectation 1/(m+1), so it ranks variant density in part.
  sidak       1 - (1 - p_min)^m_eff, with m_eff the effective number of independent
              tests from the eigenvalues of the variant correlation matrix by Li & Ji
              (2005). A p-value under a single-signal null, comparable across windows
              of different size. Declared primary in advance: pigmentation architecture
              is concentrated in single large-effect alleles, which is the regime where
              the best-SNP family is the more powerful, so this choice is the harder one
              for the CNN to beat.
  meanchi2    the aggregating estimator MAGMA applies by default. The window's chi2 sum
              is a quadratic form z'z with z ~ N(0, R), hence distributed as the
              eigenvalue-weighted sum of independent chi2(1); its tail is taken by
              Satterthwaite moment matching to a scaled chi2.
  gates       the extended Simes procedure of Li et al. (2011): min over k of
              m_eff * p_(k) / m_eff(k). The nested effective counts m_eff(k) use the
              Cheverud/Nyholt variance-of-eigenvalues estimator, which reduces to a
              cumulative sum of squared off-diagonal correlations and so is computable
              for every k in O(m^2); the eigenvalue estimator of Li & Ji would need one
              decomposition per k. The total m_eff in the numerator uses the same
              estimator as the nested ones, so the ratio is internally consistent.

Every score is returned as S_g = -log10 p_gene, computed in log space, so a window
whose p_gene underflows a double is still ranked rather than tied at zero.

Usage (the conda env carries scipy; the repo venv does not):
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/gwas_gene_ranking.py --stage all
"""
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO_ROOT = Path("/home/breno/I2CA/genomics")
DATASET = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
RESULTS = REPO_ROOT / "results" / "genotype_based_predictor" / "gwas_ranking"
CNN_TABLE = REPO_ROOT / "results" / "genotype_based_predictor" / "poolmax_final_table.csv"
CONTROL_JSON = (REPO_ROOT / "results" / "genotype_based_predictor"
                / "random_gene_control_pigmentation.json")

PANEL = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR",
         "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"]
# EDAR and TCHH are in the predecessor's panel but are hair loci with no pigmentation
# phenotype of any evidence class, and the draft excludes them from the analysis. They
# are scanned here so their scores are on record, and dropped from the head-to-head.
NOT_PIGMENTATION = {"EDAR", "TCHH"}

STRONG = {"YRI", "ESN", "LWK", "MSL", "GWD"}
WEAK = {"FIN", "CEU", "GBR"}
CROP = 32768
MAF = 0.01
GW = 5e-8
ACGT = set("ACGT")

_SITES: dict = {}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


def plink(*args, cwd: Path) -> None:
    cmd = ["plink2", *[str(a) for a in args], "--threads", "16", "--memory", "32000"]
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stdout[-4000:] + "\n" + p.stderr[-4000:] + "\n")
        raise RuntimeError(f"plink2 failed: {' '.join(cmd)}")


# random_11_2: the second pre-existing random draw, taken whole. Both draws predate
# and are independent of the specificity question, which is what makes the comparison
# interpretable; a draw is never reordered, substituted or trimmed.
DRAW2 = ["CD47", "COA1", "EGF", "FAM234B", "FYB1", "HSH2D",
         "LYNX1", "OR11H12", "OR51S1", "TRHR", "TSPAN11"]

# random_11_3 (in spirit; source directory carries no numeric suffix -- see
# scripts/experiments/link_draw3_preview.py). The third and earliest of the three
# pre-existing draws, made before and independently of the specificity question, for the
# unrelated non-longevous-dataset project. Taken whole, same as DRAW2.
DRAW3 = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6",
         "KIAA0319", "RIDA", "SEM1", "SFMBT2", "SUMF2"]


def controls(include_draw2: bool = True, include_draw3: bool = False) -> list[str]:
    c = list(json.loads(CONTROL_JSON.read_text())["control_genes"])
    if include_draw2:
        c = c + DRAW2
    if include_draw3:
        c = c + DRAW3
    return c


# ---------------------------------------------------------------- cohort


def cohort(work: Path) -> dict:
    import csv
    rows = list(csv.DictReader(open(DATASET / "selected_samples.csv")))
    for r in rows:
        assert r["Population"] in STRONG or r["Population"] in WEAK, r["Population"]
    ids = [r["SampleID"] for r in rows]
    pheno = {r["SampleID"]: (2 if r["Population"] in STRONG else 1) for r in rows}
    founders = [r["SampleID"] for r in rows if r["FatherID"] == "0" and r["MotherID"] == "0"]
    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(13)
    shuffled = list(pheno.values())
    rng.shuffle(shuffled)
    random_pheno = dict(zip(ids, shuffled))
    with open(work / "pheno.tsv", "w") as f:
        f.write("#IID\tPIGM\tRANDOM\n")
        for s in ids:
            f.write(f"{s}\t{pheno[s]}\t{random_pheno[s]}\n")
    for name, subset in (("keep_all", ids), ("keep_founders", founders)):
        with open(work / f"{name}.txt", "w") as f:
            f.write("#IID\n")
            for s in subset:
                f.write(s + "\n")
    _log(f"cohort: {len(ids)} individuals, {sum(1 for s in ids if pheno[s]==2)} strong / "
         f"{sum(1 for s in ids if pheno[s]==1)} weak, {len(founders)} founders")
    return {"ids": ids, "pheno": pheno, "founders": founders,
            "populations": sorted({r["Population"] for r in rows})}


# ---------------------------------------------------------------- build


def window_bounds(gene: str) -> tuple[str, int, int]:
    wm = json.loads((DATASET / "references" / "windows" / gene / "window_metadata.json").read_text())
    return wm["chromosome"], int(wm["start"]), int(wm["end"])


def site_list(gene: str, probe_id: str) -> list[tuple[str, int, str, str]]:
    vcf = DATASET / "individuals" / probe_id / "windows" / gene / f"{probe_id}.window.vcf.gz"
    out = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 6)
            ref, alt = c[3], c[4]
            if len(ref) != 1 or len(alt) != 1 or ref not in ACGT or alt not in ACGT:
                continue
            out.append((c[0], int(c[1]), ref, alt))
    return out


def _read_dosages(sid: str) -> tuple[np.ndarray, np.ndarray]:
    """Sparse (union column, dosage) for one individual over every site of every window.

    A site inside the OCA2/HERC2 overlap is reached twice, through both gene
    directories, and maps to the same union column with the same dosage, so the
    second write is idempotent.
    """
    idx, dos = [], []
    for gene in _SITES["genes"]:
        vcf = DATASET / "individuals" / sid / "windows" / gene / f"{sid}.window.vcf.gz"
        index = _SITES["index"]
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                c = line.rstrip("\n").split("\t")
                col = index.get((c[0], int(c[1]), c[3], c[4]))
                if col is None:
                    continue
                gt = c[9].split(":")[0]
                d = gt.count("1")
                if d and "." not in gt:
                    idx.append(col)
                    dos.append(d)
    return np.asarray(idx, dtype=np.int32), np.asarray(dos, dtype=np.int8)


def _init_worker(sites):
    global _SITES
    _SITES = sites


def write_bed(prefix: Path, geno: np.ndarray, bim_rows: list, fam_ids: list, pheno: dict) -> None:
    n_var, n_samp = geno.shape
    codes = np.empty_like(geno, dtype=np.uint8)
    codes[geno == 2] = 0
    codes[geno == 1] = 2
    codes[geno == 0] = 3
    pad = (-n_samp) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros((n_var, pad), np.uint8)], axis=1)
    q = codes.reshape(n_var, -1, 4)
    packed = (q[:, :, 0] | (q[:, :, 1] << 2) | (q[:, :, 2] << 4) | (q[:, :, 3] << 6)).astype(np.uint8)
    with open(prefix.with_suffix(".bed"), "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        f.write(packed.tobytes())
    with open(prefix.with_suffix(".bim"), "w") as f:
        for chrom, pos, ref, alt in bim_rows:
            f.write(f"{chrom}\t{chrom}:{pos}:{ref}:{alt}\t0\t{pos}\t{alt}\t{ref}\n")
    with open(prefix.with_suffix(".fam"), "w") as f:
        for s in fam_ids:
            f.write(f"0\t{s}\t0\t0\t0\t{pheno[s]}\n")


def build(work: Path, coh: dict, genes: list[str], jobs: int) -> dict:
    ids = coh["ids"]
    bounds = {g: window_bounds(g) for g in genes}

    # The union site table is de-duplicated, so every variant is tested exactly once;
    # window membership is kept separately and is NOT de-duplicated, so a variant in
    # the OCA2/HERC2 overlap belongs to both windows and each window is complete. That
    # is what the gene score is defined over.
    union: dict[tuple, int] = {}
    order: list[tuple] = []
    members: dict[str, list[tuple]] = {}
    for g in genes:
        sites = site_list(g, ids[0])
        members[g] = sites
        for k in sites:
            if k not in union:
                union[k] = -1
                order.append(k)
    order.sort(key=lambda k: (int(k[0][3:].replace("X", "23").replace("Y", "24")), k[1], k[3]))
    for i, k in enumerate(order):
        union[k] = i
    _log(f"union: {len(order)} biallelic SNV sites over {len(genes)} windows "
         f"({sum(len(v) for v in members.values())} window memberships)")

    geno = np.zeros((len(ids), len(order)), dtype=np.int8)
    t0 = time.monotonic()
    payload = {"genes": genes, "index": union}
    with mp.Pool(jobs, initializer=_init_worker, initargs=(payload,)) as pool:
        for r, (idx, dos) in enumerate(pool.imap(_read_dosages, ids, chunksize=8)):
            geno[r, idx] = dos
            if (r + 1) % 200 == 0:
                _log(f"  {r+1}/{len(ids)} individuals ({(r+1)/(time.monotonic()-t0):.1f}/s)")
    _log(f"genotype matrix {geno.shape}, {int((geno>0).sum())} non-reference calls")

    np.save(work / "geno.npy", geno)
    crop_of = {}
    for g in genes:
        _, start, end = bounds[g]
        c = (start + end) // 2
        crop_of[g] = (c - CROP // 2, c - CROP // 2 + CROP)
    with open(work / "union_sites.tsv", "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tCOL\n")
        for k in order:
            f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{union[k]}\n")

    bim_rows = [(k[0][3:], k[1], k[2], k[3]) for k in order]
    write_bed(work / "union_cohort", geno.T, bim_rows, ids, coh["pheno"])
    plink("--bfile", "union_cohort", "--maf", MAF, "--make-pgen", "--out", "union_maf", cwd=work)
    kept = [l.split("\t")[2] for l in open(work / "union_maf.pvar") if not l.startswith("#")]
    _log(f"{len(kept)} sites at MAF >= {MAF}")

    info = {
        "genes": genes,
        "n_sites_union_raw": len(order),
        "n_sites_union_maf": len(kept),
        "windows": {g: {"chrom": bounds[g][0], "start": bounds[g][1], "end": bounds[g][2],
                        "crop_start": crop_of[g][0], "crop_end": crop_of[g][1],
                        "n_sites_raw": len(members[g])} for g in genes},
    }
    return info


# ---------------------------------------------------------------- association


def run_arm(work: Path) -> None:
    plink("--pfile", "union_maf", "--indep-pairwise", "200kb", 0.1,
          "--out", "union_prune", cwd=work)
    plink("--pfile", "union_maf", "--extract", "union_prune.prune.in", "--pca", 10, "approx",
          "--out", "union_pcs", cwd=work)
    plink("--pfile", "union_maf", "--pheno", "pheno.tsv", "--pheno-name", "PIGM",
          "--glm", "firth-fallback", "allow-no-covars", "omit-ref", "log10",
          "--out", "union_uncorrected", cwd=work)
    plink("--pfile", "union_maf", "--pheno", "pheno.tsv", "--pheno-name", "PIGM",
          "--covar", "union_pcs.eigenvec", "--covar-variance-standardize",
          "--glm", "firth-fallback", "hide-covar", "omit-ref", "log10",
          "--out", "union_pc10", cwd=work)
    plink("--pfile", "union_maf", "--keep", "keep_founders.txt", "--pheno", "pheno.tsv",
          "--pheno-name", "PIGM", "--maf", MAF,
          "--glm", "firth-fallback", "allow-no-covars", "omit-ref", "log10",
          "--out", "union_founders", cwd=work)
    plink("--pfile", "union_maf", "--pheno", "pheno.tsv", "--pheno-name", "RANDOM",
          "--glm", "firth-fallback", "allow-no-covars", "omit-ref", "log10",
          "--out", "union_randomlabel", cwd=work)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("build", "run", "all"), default="all")
    ap.add_argument("--work", type=Path,
                    default=Path(os.environ.get("GWAS_RANK_WORK", "/tmp/gwas_gene_ranking")))
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--draw1-only", action="store_true",
                    help="eleven controls instead of twenty-two (reproduces the earlier run)")
    ap.add_argument("--include-draw3", action="store_true",
                    help="thirty-three controls instead of twenty-two (adds the third "
                         "pre-existing random draw)")
    args = ap.parse_args()

    os.chdir(REPO_ROOT)
    if RESULTS.resolve() == (REPO_ROOT / "results" / "genotype_based_predictor" / "gwas").resolve():
        raise SystemExit("ABORT: results dir collides with the published GWAS arm")
    RESULTS.mkdir(parents=True, exist_ok=True)
    work = args.work
    work.mkdir(parents=True, exist_ok=True)
    if args.include_draw3 and args.draw1_only:
        raise SystemExit("ABORT: --draw1-only and --include-draw3 are mutually exclusive")
    genes = PANEL + controls(not args.draw1_only, args.include_draw3)
    n_expected = 22 if args.draw1_only else (44 if args.include_draw3 else 33)
    assert len(set(genes)) == len(genes) == n_expected, genes
    coh = cohort(work)
    info_path = work / "build_info.json"

    if args.stage in ("build", "all"):
        info_path.write_text(json.dumps(build(work, coh, genes, args.jobs), indent=2))
    if args.stage in ("run", "all"):
        _log("association: union_maf")
        run_arm(work)
        _log("association done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
