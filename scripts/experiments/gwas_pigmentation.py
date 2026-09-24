#!/usr/bin/env python3
"""GWAS on the pigmentation proxy label -- the estimator the probe is meant to replace.

Two arms, both on the 1072-individual eight-population cohort of the pigmentation task
(703 strong = AFR populations, 369 weak = EUR populations):

  panel  the eleven 524,288 bp gene windows the frozen model annotates -- i.e. exactly
         the sequence the downstream classifier is shown. Genotypes come from the
         per-individual window VCFs already in the dataset; the site list is identical
         across individuals, so the union is one individual's site list.
  chr15  all of chromosome 15, from the 1kGP high-coverage phased panel. Three of the
         eleven genes (OCA2, HERC2, SLC24A5) sit on it, so it supplies the
         chromosome-scale background against which their windows can be read.

Each arm is tested three ways with identical plink2 machinery: uncorrected, corrected
with ten genotype principal components (the standard remedy for stratification), and
uncorrected on the 840 pedigree founders (to separate relatedness from ancestry).

Usage:
  /home/breno/miniforge3/envs/genomics/bin/python3 scripts/experiments/gwas_pigmentation.py --stage all
  ... --stage build     # genotype matrices only
  ... --stage run       # plink2 association only
  ... --stage summarise # JSON summary only
"""
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO_ROOT = Path("/home/breno/I2CA/genomics")
DATASET = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
CHR15_VCF = (DATASET / "raw_variants" / "vcf_chromosomes"
             / "1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
RESULTS = REPO_ROOT / "results" / "genotype_based_predictor" / "gwas"

GENES = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR",
         "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"]
STRONG = {"YRI", "ESN", "LWK", "MSL", "GWD"}
WEAK = {"FIN", "CEU", "GBR"}
CROP = 32768
MAF = 0.01
GW = 5e-8
ACGT = set("ACGT")

# Fine-mapped pigmentation variants on chr15, for the rank lookup. Coordinates are the
# GRCh38 positions already resolved from Ensembl in rsid_coordinates_grch38.json.
KNOWN_CHR15 = {"rs1426654": 48134287, "rs12913832": 28120472}

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
    # RANDOM is a calibration control: the same 703/369 case-control split assigned
    # independently of ancestry. It shares every genotype, filter and test with PIGM,
    # so any inflation it does not show is attributable to the label, not the machinery.
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
    with open(work / "sample_info.tsv", "w") as f:
        f.write("IID\tPOP\tSUPERPOP\tSEX\tFOUNDER\tPIGM\n")
        for r in rows:
            s = r["SampleID"]
            fo = int(r["FatherID"] == "0" and r["MotherID"] == "0")
            f.write(f"{s}\t{r['Population']}\t{r['Superpopulation']}\t{r['Sex']}\t{fo}\t{pheno[s]}\n")
    _log(f"cohort: {len(ids)} individuals, {sum(1 for s in ids if pheno[s]==2)} strong / "
         f"{sum(1 for s in ids if pheno[s]==1)} weak, {len(founders)} founders")
    return {"ids": ids, "pheno": pheno, "founders": founders,
            "populations": sorted({r["Population"] for r in rows})}


# ---------------------------------------------------------------- panel arm


def window_bounds(gene: str) -> tuple[str, int, int]:
    wm = json.loads((DATASET / "references" / "windows" / gene / "window_metadata.json").read_text())
    return wm["chromosome"], int(wm["start"]), int(wm["end"])


def site_list(gene: str, probe_id: str) -> list[tuple[int, str, str]]:
    """Biallelic SNV sites in a gene window. Identical across individuals, so one file suffices."""
    vcf = DATASET / "individuals" / probe_id / "windows" / gene / f"{probe_id}.window.vcf.gz"
    out = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 6)
            pos, ref, alt = int(c[1]), c[3], c[4]
            if len(ref) != 1 or len(alt) != 1 or ref not in ACGT or alt not in ACGT:
                continue
            out.append((pos, ref, alt))
    return out


def _read_dosages(sid: str) -> tuple[np.ndarray, np.ndarray]:
    """Sparse (column index, dosage) for one individual over every kept panel site."""
    idx, dos = [], []
    for gene, index in _SITES.items():
        vcf = DATASET / "individuals" / sid / "windows" / gene / f"{sid}.window.vcf.gz"
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                c = line.rstrip("\n").split("\t")
                col = index.get((int(c[1]), c[3], c[4]))
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
    """PLINK1 binary, SNP-major. A1 = ALT, A2 = REF, no missing calls."""
    n_var, n_samp = geno.shape
    codes = np.empty_like(geno, dtype=np.uint8)
    codes[geno == 2] = 0    # hom A1
    codes[geno == 1] = 2    # het
    codes[geno == 0] = 3    # hom A2
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
    # FID is left at 0 so plink2 drops the column, matching the VCF-derived chr15 arm
    # and letting one --pheno/--keep file with a bare #IID header serve both.
    with open(prefix.with_suffix(".fam"), "w") as f:
        for s in fam_ids:
            f.write(f"0\t{s}\t0\t0\t0\t{pheno[s]}\n")


def build_panel(work: Path, coh: dict, jobs: int) -> dict:
    ids = coh["ids"]
    bounds = {g: window_bounds(g) for g in GENES}
    # One flat, de-duplicated site table: the OCA2 and HERC2 windows overlap, so a site
    # can appear in two windows. Each site is attributed to the window whose centre is
    # nearest, which is also the window the CNN reads it in at highest weight.
    seen: dict[tuple, tuple[str, int]] = {}
    for g in GENES:
        chrom, start, end = bounds[g]
        centre = (start + end) / 2.0
        for pos, ref, alt in site_list(g, ids[0]):
            key = (chrom, pos, ref, alt)
            prev = seen.get(key)
            if prev is None or abs(pos - centre) < prev[1]:
                seen[key] = (g, abs(pos - centre))
    keys = sorted(seen, key=lambda k: (int(k[0][3:]), k[1], k[3]))
    col_of = {k: i for i, k in enumerate(keys)}
    per_gene_index: dict[str, dict] = {g: {} for g in GENES}
    for k in keys:
        chrom, pos, ref, alt = k
        per_gene_index[seen[k][0]][(pos, ref, alt)] = col_of[k]
    _log(f"panel: {len(keys)} biallelic SNV sites over {len(GENES)} windows "
         f"({sum(len(v) for v in per_gene_index.values())} after de-duplication)")

    geno = np.zeros((len(ids), len(keys)), dtype=np.int8)
    t0 = time.monotonic()
    with mp.Pool(jobs, initializer=_init_worker, initargs=(per_gene_index,)) as pool:
        for r, (idx, dos) in enumerate(pool.imap(_read_dosages, ids, chunksize=8)):
            geno[r, idx] = dos
            if (r + 1) % 200 == 0:
                _log(f"  {r+1}/{len(ids)} individuals ({(r+1)/(time.monotonic()-t0):.1f}/s)")
    _log(f"panel genotype matrix {geno.shape}, {int((geno>0).sum())} non-reference calls")

    crop_of = {}
    for g in GENES:
        chrom, start, end = bounds[g]
        c = (start + end) // 2
        crop_of[g] = (c - CROP // 2, c - CROP // 2 + CROP)
    with open(work / "panel_sites.tsv", "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tGENE\tIN_CROP\n")
        for k in keys:
            chrom, pos, ref, alt = k
            g = seen[k][0]
            lo, hi = crop_of[g]
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{g}\t{int(lo <= pos < hi)}\n")

    bim_rows = [(k[0][3:], k[1], k[2], k[3]) for k in keys]
    write_bed(work / "panel_cohort", geno.T, bim_rows, ids, coh["pheno"])
    plink("--bfile", "panel_cohort", "--maf", MAF, "--make-pgen", "--out", "panel_maf", cwd=work)
    n_kept = sum(1 for _ in open(work / "panel_maf.pvar")) - 1
    _log(f"panel: {n_kept} sites at MAF >= {MAF}")
    return {"n_sites_raw": len(keys), "n_sites_maf": n_kept,
            "windows": {g: {"chrom": bounds[g][0], "start": bounds[g][1], "end": bounds[g][2],
                            "crop_start": crop_of[g][0], "crop_end": crop_of[g][1]} for g in GENES}}


def build_chr15(work: Path) -> dict:
    if not (work / "chr15_maf.pgen").exists():
        plink("--vcf", CHR15_VCF, "--keep", "keep_all.txt", "--snps-only", "just-acgt",
              "--max-alleles", 2, "--min-alleles", 2,
              "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", 60, "missing",
              "--maf", MAF, "--geno", 0.05, "--make-pgen", "--out", "chr15_maf", cwd=work)
    n = sum(1 for _ in open(work / "chr15_maf.pvar")) - 1
    _log(f"chr15: {n} biallelic SNVs at MAF >= {MAF}")
    return {"n_sites_maf": n}


# ---------------------------------------------------------------- association


def run_arm(work: Path, prefix: str) -> None:
    plink("--pfile", prefix, "--indep-pairwise", "200kb", 0.1, "--out", f"{prefix}_prune", cwd=work)
    plink("--pfile", prefix, "--extract", f"{prefix}_prune.prune.in", "--pca", 10, "approx",
          "--out", f"{prefix}_pcs", cwd=work)
    plink("--pfile", prefix, "--pheno", "pheno.tsv", "--pheno-name", "PIGM",
          "--glm", "firth-fallback", "allow-no-covars", "omit-ref",
          "--out", f"{prefix}_uncorrected", cwd=work)
    plink("--pfile", prefix, "--pheno", "pheno.tsv", "--pheno-name", "PIGM",
          "--covar", f"{prefix}_pcs.eigenvec", "--covar-variance-standardize",
          "--glm", "firth-fallback", "hide-covar", "omit-ref",
          "--out", f"{prefix}_pc10", cwd=work)
    plink("--pfile", prefix, "--keep", "keep_founders.txt", "--pheno", "pheno.tsv",
          "--pheno-name", "PIGM", "--maf", MAF,
          "--glm", "firth-fallback", "allow-no-covars", "omit-ref",
          "--out", f"{prefix}_founders", cwd=work)
    plink("--pfile", prefix, "--pheno", "pheno.tsv", "--pheno-name", "RANDOM",
          "--glm", "firth-fallback", "allow-no-covars", "omit-ref",
          "--out", f"{prefix}_randomlabel", cwd=work)


# ---------------------------------------------------------------- summary


def read_sumstats(path: Path):
    import pandas as pd
    d = pd.read_csv(path, sep="\t", dtype={"#CHROM": str})
    d = d[d.TEST == "ADD"].copy()
    d["ERRCODE"] = d["ERRCODE"].astype(str)
    return d


def lambda_gc(p) -> float:
    from scipy import stats
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return float("nan")
    chisq = stats.chi2.isf(np.clip(p, 1e-300, 1.0), 1)
    return float(np.median(chisq) / stats.chi2.ppf(0.5, 1))


def arm_summary(work: Path, prefix: str, tag: str) -> dict:
    name = "RANDOM" if tag == "randomlabel" else "PIGM"
    f = work / f"{prefix}_{tag}.{name}.glm.logistic.hybrid"
    d = read_sumstats(f)
    ok = d[(d.ERRCODE == ".") & d.P.notna()]
    out = {
        "n_tests": int(len(d)),
        "n_converged": int(len(ok)),
        "n_unfinished": int((d.ERRCODE != ".").sum()),
        "lambda_gc": lambda_gc(ok.P) if len(ok) else None,
        "lambda_gc_all_reported": lambda_gc(d.P.dropna()) if d.P.notna().any() else None,
        "n_genomewide": int((ok.P < GW).sum()) if len(ok) else 0,
        "frac_genomewide": float((ok.P < GW).mean()) if len(ok) else 0.0,
        "min_p": float(ok.P.min()) if len(ok) else None,
    }
    return out


def summarise(work: Path, coh: dict, build_info: dict) -> dict:
    import pandas as pd
    from sklearn.metrics import roc_auc_score

    summary = {"generated": _now(), "cohort": {
        "n": len(coh["ids"]), "n_strong": sum(1 for s in coh["ids"] if coh["pheno"][s] == 2),
        "n_weak": sum(1 for s in coh["ids"] if coh["pheno"][s] == 1),
        "n_founders": len(coh["founders"]), "populations": coh["populations"],
        "strong_populations": sorted(STRONG), "weak_populations": sorted(WEAK)},
        "maf_threshold": MAF, "genomewide_threshold": GW, "build": build_info, "arms": {}}

    info = pd.read_csv(work / "sample_info.tsv", sep="\t")
    for prefix in ("panel_maf", "chr15_maf"):
        arm = {t: arm_summary(work, prefix, t)
               for t in ("uncorrected", "pc10", "founders", "randomlabel")}
        ev = [float(x.split()[0]) for x in open(work / f"{prefix}_pcs.eigenval")]
        pcs = pd.read_csv(work / f"{prefix}_pcs.eigenvec", sep="\t").rename(columns={"#IID": "IID"})
        m = info.merge(pcs, on="IID")
        y = (m.PIGM == 2).astype(int).values
        s, w = m.PC1[y == 1], m.PC1[y == 0]
        arm["structure"] = {
            "eigenvalues": [round(e, 4) for e in ev],
            # PC sign is arbitrary, so report the orientation-invariant AUC.
            "pc1_auc_vs_label": float(max(roc_auc_score(y, m.PC1), 1.0 - roc_auc_score(y, m.PC1))),
            "pc1_point_biserial_r": float(abs(np.corrcoef(m.PC1, y)[0, 1])),
            "pc1_ranges_disjoint": bool(s.min() > w.max() or w.min() > s.max()),
        }
        summary["arms"][prefix.replace("_maf", "")] = arm

    # Where the panel arm's significant sites fall, gene by gene.
    sites = pd.read_csv(work / "panel_sites.tsv", sep="\t")
    d = read_sumstats(work / "panel_maf_uncorrected.PIGM.glm.logistic.hybrid")
    d = d[(d.ERRCODE == ".") & d.P.notna()]
    d["POS"] = d.POS.astype(int)
    j = d.merge(sites, on="POS", how="left")
    per_gene = {}
    for g, sub in j.groupby("GENE"):
        per_gene[g] = {
            "n_tested": int(len(sub)),
            "n_genomewide": int((sub.P < GW).sum()),
            "frac_genomewide": float((sub.P < GW).mean()),
            "min_p": float(sub.P.min()),
            "min_p_pos": int(sub.loc[sub.P.idxmin(), "POS"]),
            "n_tested_in_crop": int((sub.IN_CROP == 1).sum()),
            "n_genomewide_in_crop": int(((sub.IN_CROP == 1) & (sub.P < GW)).sum()),
        }
    summary["panel_per_gene"] = dict(sorted(per_gene.items(), key=lambda kv: kv[1]["min_p"]))

    # The panel ranking the GWAS induces, against the one the knockdown probe induces.
    # Neither is a claim about pigmentation biology on this cohort; the comparison is
    # between the spread the two estimators produce over the same eleven windows.
    probe_path = REPO_ROOT / "results" / "genotype_based_predictor" / "readout_normalisation.json"
    if probe_path.exists():
        from scipy import stats as _st
        probe = {r["gene"]: r for r in json.loads(probe_path.read_text())["per_gene"]}
        genes = [g for g in per_gene if g in probe]
        gw_score = np.array([-np.log10(per_gene[g]["min_p"]) for g in genes])
        gw_frac = np.array([per_gene[g]["frac_genomewide"] for g in genes])
        pr = np.array([probe[g]["abs_delta"] for g in genes])
        summary["probe_comparison"] = {
            "genes": genes,
            "gwas_neglog10_min_p": [round(float(x), 2) for x in gw_score],
            "gwas_frac_genomewide": [round(float(x), 4) for x in gw_frac],
            "probe_abs_delta": [float(probe[g]["abs_delta"]) for g in genes],
            "probe_null_mean_abs": [float(probe[g]["null_mean_abs"]) for g in genes],
            "spearman_minp_vs_probe": {
                "rho": float(_st.spearmanr(gw_score, pr).statistic),
                "p": float(_st.spearmanr(gw_score, pr).pvalue)},
            "spearman_frac_vs_probe": {
                "rho": float(_st.spearmanr(gw_frac, pr).statistic),
                "p": float(_st.spearmanr(gw_frac, pr).pvalue)},
            "gwas_spread": {
                "neglog10_min_p_range": [float(gw_score.min()), float(gw_score.max())],
                "frac_genomewide_range": [float(gw_frac.min()), float(gw_frac.max())],
                "all_below_threshold": bool(all(per_gene[g]["min_p"] < GW for g in genes))},
            "probe_spread": {"abs_delta_range": [float(pr.min()), float(pr.max())],
                             "null_mean_abs": float(np.mean(
                                 [probe[g]["null_mean_abs"] for g in genes]))},
        }

    # Rank of the two fine-mapped chr15 variants inside the chr15 arm.
    c = read_sumstats(work / "chr15_maf_uncorrected.PIGM.glm.logistic.hybrid")
    c = c[(c.ERRCODE == ".") & c.P.notna()].sort_values("P").reset_index(drop=True)
    c["rank"] = np.arange(1, len(c) + 1)
    known = {}
    for rsid, pos in KNOWN_CHR15.items():
        hit = c[c.POS == pos]
        known[rsid] = ({"pos": pos, "present": False} if hit.empty else {
            "pos": pos, "present": True, "rank": int(hit["rank"].iloc[0]),
            "n_tests": int(len(c)), "p": float(hit.P.iloc[0]),
            "odds_ratio": float(hit.OR.iloc[0]), "a1": str(hit.A1.iloc[0]),
            "percentile": float(hit["rank"].iloc[0] / len(c))})
    summary["chr15_known_variants"] = known

    # chr15 significance inside vs outside the three panel windows on that chromosome.
    wins = [(v["start"], v["end"]) for g, v in build_info["panel"]["windows"].items()
            if v["chrom"] == "chr15"]
    inside = np.zeros(len(c), dtype=bool)
    for lo, hi in wins:
        inside |= (c.POS.values >= lo) & (c.POS.values <= hi)
    summary["chr15_window_partition"] = {
        "n_windows": len(wins),
        "bp_in_windows": int(sum(hi - lo + 1 for lo, hi in wins)),
        "inside": {"n_tested": int(inside.sum()),
                   "n_genomewide": int(((c.P.values < GW) & inside).sum())},
        "outside": {"n_tested": int((~inside).sum()),
                    "n_genomewide": int(((c.P.values < GW) & ~inside).sum())},
    }
    return summary


def export_sumstats(work: Path) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    for prefix in ("panel_maf", "chr15_maf"):
        for tag in ("uncorrected", "pc10", "founders", "randomlabel"):
            name = "RANDOM" if tag == "randomlabel" else "PIGM"
            d = read_sumstats(work / f"{prefix}_{tag}.{name}.glm.logistic.hybrid")
            keep = ["#CHROM", "POS", "ID", "A1", "A1_FREQ", "OR", "P", "ERRCODE"]
            out = RESULTS / f"{prefix.replace('_maf','')}_{tag}.sumstats.tsv.gz"
            d[keep].to_csv(out, sep="\t", index=False, compression="gzip")
    shutil.copy(work / "panel_sites.tsv", RESULTS / "panel_sites.tsv")
    shutil.copy(work / "sample_info.tsv", RESULTS / "sample_info.tsv")
    for prefix in ("panel_maf", "chr15_maf"):
        shutil.copy(work / f"{prefix}_pcs.eigenvec", RESULTS / f"{prefix.replace('_maf','')}_pcs.eigenvec")
    _log(f"sumstats exported to {RESULTS}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("build", "run", "summarise", "all"), default="all")
    ap.add_argument("--work", type=Path,
                    default=Path(os.environ.get("GWAS_WORK", "/tmp/gwas_pigmentation")))
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--out", type=Path, default=RESULTS / "gwas_pigmentation.json")
    args = ap.parse_args()

    os.chdir(REPO_ROOT)
    work = args.work
    work.mkdir(parents=True, exist_ok=True)
    coh = cohort(work)
    info_path = work / "build_info.json"

    if args.stage in ("build", "all"):
        build = {"panel": build_panel(work, coh, args.jobs), "chr15": build_chr15(work)}
        info_path.write_text(json.dumps(build, indent=2))
    build = json.loads(info_path.read_text())

    if args.stage in ("run", "all"):
        for prefix in ("panel_maf", "chr15_maf"):
            _log(f"association: {prefix}")
            run_arm(work, prefix)

    if args.stage in ("summarise", "all"):
        summary = summarise(work, coh, build)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(summary, indent=2))
        export_sumstats(work)
        _log(f"wrote {args.out}")
        for k in ("arms", "panel_per_gene", "probe_comparison",
                  "chr15_known_variants", "chr15_window_partition"):
            if k in summary:
                print(f"--- {k}\n" + json.dumps(summary[k], indent=1)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
