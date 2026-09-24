#!/usr/bin/env python3
"""Does carrying a variant INSIDE the scrambled window change the knockdown response?

knockdown_site_variant_landscape.py establishes where the 100 bp scramble windows sit relative to
the cohort's variation. For most genes the window is variant-free, which is what the probe's
identifiability argument wants: the edit is not a variant, so the contrast it induces cannot be
confounded with the ancestry structure of real genotypes.

Two windows are not variant-free. MC1R's TSS+/-50 contains a common, ancestry-differentiated SNV
(chr16:89918814 G>A, AF 0.453, AF_AFR 0.429 vs AF_EUR 0.265) carried by 70% of the cohort, and the
control gene TPM2's contains a three-base A>T run at AF 0.226. For those genes the scramble does
two things at once -- it destroys the promoter motif AND it overwrites an allele the individual
actually has -- and the second part is exactly the confound the design is meant to exclude.

This tests whether that matters, per gene, on the individuals the probe was run on: split the test
split by carrier status at each variant inside the window and compare the knockdown Delta.

    no difference  -> the scramble's effect does not depend on which allele it overwrote; the
                      identifiability argument survives even at the polymorphic windows.
    difference     -> at that gene the readout mixes the promoter edit with an allele contrast,
                      and the gene's Delta cannot be read as a pure non-variant intervention.

Reads only the knockdown CSV, the landscape JSON and the on-disk per-window VCFs (tabix).
No AlphaGenome calls, no GPU.

Usage:
  python3 scripts/experiments/knockdown_site_carrier_effect.py
  ... --knockout <csv>   # default: the published 11-gene raw_center_crop run
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path("/home/breno/I2CA/genomics")
DATASET = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
LANDSCAPE = REPO_ROOT / "results/genotype_based_predictor/knockdown_site_variant_landscape.json"
OUT = REPO_ROOT / "results/genotype_based_predictor/knockdown_site_carrier_effect.json"
METHOD = "biology_tss"


def p_weak(strong, weak):
    s = np.vstack([np.asarray(strong, float), np.asarray(weak, float)])
    s = s - s.max(axis=0, keepdims=True)
    e = np.exp(s)
    return e[1] / e.sum(axis=0)


def carried_alleles(sample_id, gene, chrom, lo, hi):
    """{pos: n_alt_alleles} for variants whose POS is inside [lo, hi] on this individual."""
    vcf = DATASET / "individuals" / sample_id / "windows" / gene / f"{sample_id}.window.vcf.gz"
    out = subprocess.run(["tabix", str(vcf), f"{chrom}:{lo}-{hi}"], capture_output=True, text=True)
    res = {}
    for line in out.stdout.splitlines():
        c = line.split("\t")
        pos = int(c[1])
        if not (lo <= pos <= hi):
            continue  # spanning structural records whose POS lies outside the window
        gt = c[9].split(":")[0]
        sep = "|" if "|" in gt else "/"
        alleles = [a for a in gt.split(sep)]
        n = sum(1 for a in alleles if a not in ("0", "."))
        res[pos] = n
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--knockout", type=Path,
                    default=KO_DIR / "pigmentation_test_split_knockout_no_alignment.csv")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    land = {g["gene"]: g for g in json.loads(LANDSCAPE.read_text())["per_gene"]}
    ko = pd.read_csv(args.knockout)
    ko = ko[ko["method"] == METHOD].copy() if "method" in ko.columns else ko
    ko["delta"] = (p_weak(ko["perturbed_strong_logit"], ko["perturbed_weak_logit"])
                   - p_weak(ko["baseline_strong_logit"], ko["baseline_weak_logit"]))

    results = []
    print("=" * 96)
    print("DOES CARRYING A VARIANT INSIDE THE SCRAMBLED WINDOW CHANGE THE RESPONSE?")
    print(f"knockout = {args.knockout.name}, method = {METHOD}")
    print("=" * 96)
    for gene, sub in ko.groupby("gene"):
        info = land.get(gene)
        if info is None or info["panel_sites_in_scramble_window"] == 0:
            continue
        chrom, (lo, hi) = info["chrom"], info["scramble_window"]
        geno = {sid: carried_alleles(sid, gene, chrom, lo, hi) for sid in sub["sample_id"].unique()}
        positions = sorted({p for d in geno.values() for p in d})
        printed = False
        for pos in positions:
            n_alt = np.array([geno[s].get(pos, 0) for s in sub["sample_id"]])
            carrier = n_alt > 0
            if carrier.sum() < 5 or (~carrier).sum() < 5:
                continue  # too few on one side to compare
            d_car, d_non = sub["delta"].to_numpy()[carrier], sub["delta"].to_numpy()[~carrier]
            u, pu = stats.mannwhitneyu(d_car, d_non, alternative="two-sided")
            rho, prho = stats.spearmanr(n_alt, sub["delta"].to_numpy())
            if not printed:
                print(f"\n{gene}  {chrom}:{lo}-{hi}   n = {len(sub)}   overall Delta = {sub['delta'].mean():+.4f}")
                printed = True
            print(f"   pos {pos} (offset {pos - info['tss_1based']:+d}): "
                  f"carriers {carrier.sum():>3} Delta = {d_car.mean():+.4f} | "
                  f"non {(~carrier).sum():>3} Delta = {d_non.mean():+.4f} | "
                  f"MWU p = {pu:.4f}  dosage rho = {rho:+.3f} (p = {prho:.4f})")
            # Both of the moderate-frequency sites here are AFR-specific, and every AFR
            # individual carries the "strong pigmentation" label by construction, so a pooled
            # carrier/non-carrier difference is partly a label difference. The clean contrast is
            # WITHIN a superpopulation, where the label is constant.
            within = {}
            for sup, idx in sub.groupby("superpopulation").groups.items():
                mask = sub.index.isin(idx)
                c, n = carrier & mask, (~carrier) & mask
                if c.sum() < 5 or n.sum() < 5:
                    continue
                uu, pp = stats.mannwhitneyu(sub["delta"].to_numpy()[c],
                                            sub["delta"].to_numpy()[n], alternative="two-sided")
                within[sup] = {"n_carriers": int(c.sum()), "n_non_carriers": int(n.sum()),
                               "delta_carriers": float(sub["delta"].to_numpy()[c].mean()),
                               "delta_non_carriers": float(sub["delta"].to_numpy()[n].mean()),
                               "p": float(pp)}
                print(f"        within {sup}: carriers {c.sum():>3} Delta = {within[sup]['delta_carriers']:+.4f} | "
                      f"non {n.sum():>3} Delta = {within[sup]['delta_non_carriers']:+.4f} | p = {pp:.4f}")
            if not within:
                print("        (no superpopulation has >= 5 on both sides -- carrier status is "
                      "collinear with ancestry here)")
            results.append({
                "within_superpopulation": within,
                "gene": gene, "chrom": chrom, "pos": pos,
                "offset_from_tss": pos - info["tss_1based"],
                "n_carriers": int(carrier.sum()), "n_non_carriers": int((~carrier).sum()),
                "delta_carriers": float(d_car.mean()), "delta_non_carriers": float(d_non.mean()),
                "mannwhitney_p": float(pu), "dosage_spearman_rho": float(rho),
                "dosage_spearman_p": float(prho),
                "overall_delta": float(sub["delta"].mean()),
            })

    if results:
        ps = [r["mannwhitney_p"] for r in results]
        order = np.argsort(ps)
        m = len(ps)
        holm = np.empty(m)
        run_max = 0.0
        for rank, idx in enumerate(order):
            v = min(1.0, (m - rank) * ps[idx])
            run_max = max(run_max, v)
            holm[idx] = run_max
        for r, h in zip(results, holm):
            r["mannwhitney_p_holm"] = float(h)
        print()
        sig = [f"{r['gene']}@{r['pos']}" for r in results if r["mannwhitney_p_holm"] < 0.05]
        print(f"{m} comparisons, Holm-corrected; significant at 0.05: {sig or 'none'}")

    args.out.write_text(json.dumps({"knockout": args.knockout.name, "method": METHOD,
                                    "comparisons": results}, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
