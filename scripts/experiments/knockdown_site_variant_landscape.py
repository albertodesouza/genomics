#!/usr/bin/env python3
"""Are the knockdown sites places where the cohort actually varies?

The probe's identifiability argument rests on the intervention NOT being a variant: a 100 bp
scramble at a promoter location is an edit no individual carries, so the contrast it induces
cannot be confounded with the ancestry structure of real genotypes. That argument is only
airtight if the scrambled window is also not a place where the cohort happens to differ --
otherwise part of the delivered perturbation is "erasing an allele the individual has" rather
than "erasing a promoter".

This measures it directly, per gene and per individual, against a matched null:

  scramble window   the 100 bp actually scrambled, in REFERENCE coordinates
                    - biology_tss:  the MANE TSS, individual-invariant.
                    - cage_*:       the haplotype-local summit from the knockdown CSV, mapped
                                    back to reference coordinates by inverting the same indel
                                    drift model used to place the edit (haplotype.py).
  matched null      every other 100 bp window inside the same 32,768 bp CNN crop for the same
                    individual -- so the comparison holds gene, individual and crop fixed and
                    only moves the position.

Reported per gene: how many of the cohort's individuals carry >= 1 variant inside the scrambled
window, the SNV/INDEL split, and where the scramble window's carried-variant count falls in the
crop-wide 100 bp distribution (empirical percentile).

Also checks whether any of the fine-mapped causal variants in snps.csv falls inside a scramble
window, which is the link to the causal-variant-footprint experiment.

Reads only on-disk per-window VCFs (tabix) plus the GENCODE cache. No AlphaGenome calls, no GPU.

Usage:
  python3 scripts/experiments/knockdown_site_variant_landscape.py
  ... --individuals test          # 162 test-split individuals only (default: all)
  ... --genes SLC24A5,TYR
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _p in (REPO_ROOT / "src", REPO_ROOT / "notebooks"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
import pandas as pd

DATASET = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
KO_CSV = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout_no_alignment.csv"
OUT = REPO_ROOT / "results/genotype_based_predictor/knockdown_site_variant_landscape.json"
PANEL = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH"]
CONTROLS = ["TPM2", "SMCR8", "PSMC4"]
CROP = 32768
NULL_HALF = 16384   # matched-null half-width around the scramble centre
SCRAMBLE = 100


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def window_meta(gene):
    m = json.loads((DATASET / "references" / "windows" / gene / "window_metadata.json").read_text())
    return m["chromosome"], int(m["start"]), int(m["end"])


def crop_bounds(full_len, size=CROP):
    """Identical centering formula to knockout.cnn_crop_bounds -- the slice the CNN sees."""
    size = min(size, full_len)
    center = full_len // 2
    lo = max(0, center - size // 2)
    hi = min(full_len, lo + size)
    if hi - lo < size:
        lo = max(0, hi - size)
    return lo, hi


def parse_records(text, hap_index=None):
    """(pos_1based, is_indel, carried_h1, carried_h2) for each VCF line in `text`."""
    out = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        pos, ref, alts, gt = int(c[1]), c[3], c[4].split(","), c[9].split(":")[0]
        sep = "|" if "|" in gt else "/"
        alleles = gt.split(sep)
        if len(alleles) != 2:
            continue
        carried = []
        is_indel = False
        for a in alleles:
            if a in ("0", "."):
                carried.append(False)
                continue
            n = int(a)
            if not (1 <= n <= len(alts)):
                carried.append(False)
                continue
            carried.append(True)
            alt = alts[n - 1]
            if alt.startswith("<"):
                is_indel = True
            elif len(alt) != len(ref):
                is_indel = True
        out.append((pos, is_indel, carried[0], carried[1]))
    return out


def tabix(vcf, region):
    r = subprocess.run(["tabix", str(vcf), region], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"tabix failed on {vcf} {region}: {r.stderr[:200]}")
    return r.stdout


def indel_events(vcf_text, hap_index):
    """Replays haplotype.py's drift model from an already-read VCF region: (pos, delta) list."""
    events = []
    for line in vcf_text.splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        pos, ref, alts, info, gt = int(c[1]), c[3], c[4].split(","), c[7], c[9].split(":")[0]
        sep = "|" if "|" in gt else "/"
        alleles = gt.split(sep)
        if len(alleles) != 2:
            continue
        a = alleles[hap_index]
        if a in ("0", "."):
            continue
        n = int(a)
        if not (1 <= n <= len(alts)):
            continue
        alt = alts[n - 1]
        if alt.startswith("<"):
            if alt != "<DEL>":
                continue
            delta = None
            for f in info.split(";"):
                if f.startswith("SVLEN="):
                    delta = int(f.split("=", 1)[1].split(",")[0])
                    break
            if delta is None:
                end = next((int(f.split("=", 1)[1]) for f in info.split(";") if f.startswith("END=")), None)
                if end is None:
                    continue
                delta = -(end - pos)
        else:
            delta = len(alt) - len(ref)
        if delta:
            events.append((pos, delta))
    events.sort()
    return events


def ref_pos_from_local(events, start_1based, local_idx):
    """Inverse of haplotype.haplotype_local_idx: local 0-based index -> reference 1-based pos.

    haplotype_local_idx(p) = (p - start) + sum(delta for pos<p), monotone non-decreasing in p, so
    a scan over the (few) drift events brackets the answer exactly.
    """
    drift = 0
    prev_pos = start_1based
    for pos, delta in events + [(10**12, 0)]:
        # on [prev_pos, pos) the drift is constant
        lo_local = (prev_pos - start_1based) + drift
        hi_local = (pos - start_1based) + drift
        if lo_local <= local_idx < hi_local:
            return local_idx - drift + start_1based
        drift += delta
        prev_pos = pos
    return local_idx + start_1based  # past the last event


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--individuals", choices=["test", "all"], default="all")
    ap.add_argument("--genes", default=None)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    genes = args.genes.split(",") if args.genes else PANEL + CONTROLS

    from genotype_cnn_alignment_deeplift_summary.annotations import build_tss_tables, build_transcript_extractors, load_gtf
    _log("loading GENCODE...")
    gtf = load_gtf(REPO_ROOT / "notebooks" / ".cache" / "annotations")
    gtf_mane, gtf_pc, _, _, _ = build_transcript_extractors(gtf)
    tss_mane, tss_pc = build_tss_tables(gtf_mane, gtf_pc)

    def gene_tss(g):
        rows = tss_mane[tss_mane["gene_name"] == g]
        if rows.empty:
            rows = tss_pc[tss_pc["gene_name"] == g]
        return int(rows.iloc[0]["Start"])  # 0-based

    ko = pd.read_csv(KO_CSV)
    test_ids = sorted(ko["sample_id"].unique())
    if args.individuals == "test":
        sample_ids = test_ids
    else:
        sample_ids = sorted(p.name for p in (DATASET / "individuals").iterdir() if p.is_dir())
    _log(f"{len(sample_ids)} individuals, {len(genes)} genes")

    # CAGE summit reference positions, per (sample, gene, haplotype, method), from the knockdown CSV.
    cage_local = defaultdict(dict)
    for _, r in ko.iterrows():
        if r["method"] == "biology_tss":
            continue
        cage_local[(r["sample_id"], r["gene"], r["method"])] = (
            int(r["h1_target_local_idx"]), int(r["h2_target_local_idx"]))

    causal = pd.read_csv(REPO_ROOT / "snps.csv") if (REPO_ROOT / "snps.csv").exists() else None
    rsid_coords = {}
    p = REPO_ROOT / "results/genotype_based_predictor/rsid_coordinates_grch38.json"
    if p.exists():
        rsid_coords = json.loads(p.read_text())

    results = {"generated_at": datetime.now(timezone.utc).isoformat(),
               "n_individuals": len(sample_ids), "scramble_window_bp": SCRAMBLE,
               "crop_bp": CROP, "per_gene": [], "cage_methods": {}}

    for gene in genes:
        chrom, start_1b, end_1b = window_meta(gene)
        full_len = end_1b - start_1b + 1
        lo, hi = crop_bounds(full_len)
        crop_start_1b, crop_end_1b = start_1b + lo, start_1b + hi - 1
        tss0 = gene_tss(gene)
        tss_1b = tss0 + 1
        sc_lo, sc_hi = tss_1b - SCRAMBLE // 2, tss_1b - SCRAMBLE // 2 + SCRAMBLE - 1
        in_crop = crop_start_1b <= tss_1b <= crop_end_1b

        # Matched null region: +/- NULL_HALF bp centred on the scramble window, so the null is
        # defined whether or not the scramble site falls inside the CNN crop (for genes whose
        # window is centred on a long gene body, the TSS does not).
        null_lo, null_hi = tss_1b - NULL_HALF, tss_1b + NULL_HALF
        scramble_in_crop = crop_start_1b <= sc_lo and sc_hi <= crop_end_1b

        n_carrier = 0            # individuals with >=1 carried variant in the scramble window
        n_carrier_indel = 0
        carried_counts = []      # per individual, carried variants inside the scramble window
        null_pct = []            # percentile of that count in the individual's own local 100bp dist
        cohort_sites = None
        local_dens = []
        crop_dens = []
        cage_ref_offsets = defaultdict(list)
        cage_carrier = defaultdict(int)
        cage_n = defaultdict(int)

        for sid in sample_ids:
            vcf = DATASET / "individuals" / sid / "windows" / gene / f"{sid}.window.vcf.gz"
            if not vcf.exists():
                continue
            text = tabix(vcf, f"{chrom}:{null_lo}-{null_hi}")
            recs = parse_records(text)
            if cohort_sites is None:
                cohort_sites = np.array([r[0] for r in recs])  # panel sites are identical across individuals
            carried = np.array([(r[0], r[1]) for r in recs if (r[2] or r[3])], dtype=int).reshape(-1, 2)
            local_dens.append(len(carried))
            pos = carried[:, 0]
            k = int(((pos >= sc_lo) & (pos <= sc_hi)).sum())
            k_ind = int(((pos >= sc_lo) & (pos <= sc_hi) & (carried[:, 1] == 1)).sum())
            carried_counts.append(k)
            n_carrier += (k > 0)
            n_carrier_indel += (k_ind > 0)
            # matched null: disjoint 100 bp windows tiling the local region, same individual
            edges = np.arange(null_lo, null_hi + 1, SCRAMBLE)
            counts = np.histogram(pos, bins=np.append(edges, null_hi + 1))[0]
            null_pct.append(float((counts <= k).mean()))

            crop_txt = tabix(vcf, f"{chrom}:{crop_start_1b}-{crop_end_1b}")
            crop_recs = parse_records(crop_txt)
            crop_dens.append(sum(1 for r in crop_recs if (r[2] or r[3])))

            # CAGE summit locations, test individuals only (that is where the CSV has rows)
            full_text = None
            for method in ("cage_melanocyte", "cage_gene_start"):
                key = (sid, gene, method)
                if key not in cage_local:
                    continue
                if full_text is None:
                    full_text = tabix(vcf, f"{chrom}:{start_1b}-{end_1b}")
                    full_recs = parse_records(full_text)
                    full_carried = np.array([(r[0], r[1]) for r in full_recs if (r[2] or r[3])],
                                            dtype=int).reshape(-1, 2)
                for hap_i, loc in enumerate(cage_local[key]):
                    ev = indel_events(full_text, hap_i)
                    ref_pos = ref_pos_from_local(ev, start_1b, loc)
                    cage_ref_offsets[method].append(ref_pos - tss_1b)
                    a, b = ref_pos - SCRAMBLE // 2, ref_pos - SCRAMBLE // 2 + SCRAMBLE - 1
                    fp = full_carried[:, 0]
                    kk = int(((fp >= a) & (fp <= b)).sum())
                    cage_carrier[method] += (kk > 0)
                    cage_n[method] += 1

        carried_counts = np.array(carried_counts)
        local_dens = np.array(local_dens)
        crop_dens = np.array(crop_dens)
        n = len(carried_counts)
        # cohort-level: how many panel sites are polymorphic inside the scramble window at all
        n_sites_window = int(((cohort_sites >= sc_lo) & (cohort_sites <= sc_hi)).sum()) if cohort_sites is not None else 0

        rec = {
            "gene": gene, "chrom": chrom, "tss_1based": tss_1b,
            "scramble_window": [int(sc_lo), int(sc_hi)], "tss_inside_cnn_crop": bool(in_crop),
            "crop_1based": [int(crop_start_1b), int(crop_end_1b)],
            "n_individuals_measured": int(n),
            "panel_sites_in_scramble_window": n_sites_window,
            "individuals_carrying_any": int(n_carrier),
            "fraction_carrying_any": float(n_carrier / n) if n else None,
            "individuals_carrying_indel": int(n_carrier_indel),
            "mean_carried_in_scramble_window": float(carried_counts.mean()) if n else None,
            "max_carried_in_scramble_window": int(carried_counts.max()) if n else None,
            "expected_under_local_density": float(local_dens.mean() * SCRAMBLE / (2 * NULL_HALF + 1)) if n else None,
            "mean_carried_per_local_region": float(local_dens.mean()) if n else None,
            "local_region_bp": 2 * NULL_HALF + 1,
            "mean_carried_per_crop": float(crop_dens.mean()) if n else None,
            "carried_per_100bp_in_crop": float(crop_dens.mean() * SCRAMBLE / CROP) if n else None,
            "mean_percentile_vs_local_null": float(np.mean(null_pct)) if n else None,
            "scramble_window_inside_crop": bool(scramble_in_crop),
            "tss_distance_to_crop_bp": (0 if scramble_in_crop else
                                        int(min(abs(tss_1b - crop_start_1b), abs(tss_1b - crop_end_1b)))),
        }
        if rec["expected_under_local_density"]:
            rec["enrichment_vs_local"] = rec["mean_carried_in_scramble_window"] / rec["expected_under_local_density"]
        for method in ("cage_melanocyte", "cage_gene_start"):
            if cage_n.get(method):
                offs = np.array(cage_ref_offsets[method])
                rec[method] = {
                    "n_haplotypes": int(cage_n[method]),
                    "median_offset_from_tss_bp": float(np.median(offs)),
                    "iqr_offset_bp": [float(np.percentile(offs, 25)), float(np.percentile(offs, 75))],
                    "fraction_carrying_any": float(cage_carrier[method] / cage_n[method]),
                }
        # causal variants inside the scramble window
        if causal is not None:
            hits = []
            for _, row in causal.iterrows():
                rsid = str(row.get("rsid", ""))
                c = rsid_coords.get(rsid)
                if not c:
                    continue
                cp = c.get("position") or c.get("pos")
                cc = c.get("chromosome") or c.get("chrom")
                if cc and str(cc).replace("chr", "") == chrom.replace("chr", "") and cp and sc_lo <= int(cp) <= sc_hi:
                    hits.append(rsid)
            rec["causal_variants_in_scramble_window"] = hits
        results["per_gene"].append(rec)
        _log(f"{gene:8s} carriers={rec['fraction_carrying_any']:.3f} "
             f"mean={rec['mean_carried_in_scramble_window']:.3f} "
             f"exp={rec['expected_under_local_density']:.3f} "
             f"sites={n_sites_window} crop={in_crop}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2))
    _log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
