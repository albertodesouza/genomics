#!/usr/bin/env python3
"""
SNP Ancestry Predictor
======================

Predicts genetic ancestry from SNP data using allele frequency-based methods.

Pipeline:
  Step 1: Generate 23andMe-format files from multi-sample VCFs
  Step 2: Compute per-population allele frequency statistics from reference subset
  Step 3: Predict ancestry using Maximum Likelihood or Admixture models

Usage:
    source ../scripts/start_genomics_universal.sh
    python3 snp_ancestry_predictor.py --config configs/default.yaml

Author: Alberto F. De Souza
"""

import argparse
import glob as glob_module
import json
import math
import os
import resource
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import yaml
from scipy.optimize import minimize

from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeRemainingColumn,
    SpinnerColumn,
    MofNCompleteColumn,
)
from rich.table import Table

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))
from vcf_to_23andme.vcf_to_23andme import normalize_chrom, generate_header


DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

GT_DOSE = {
    "0|0": 0, "0/0": 0,
    "0|1": 1, "0/1": 1, "1|0": 1, "1/0": 1,
    "1|1": 2, "1/1": 2,
}


# ═══════════════════════════════════════════════════════════════
# Configuration & Helpers
# ═══════════════════════════════════════════════════════════════


def load_config(path: Path) -> dict:
    """Load YAML configuration file."""
    with open(path, "r") as f:
        return yaml.safe_load(f)


def load_splits(path: str) -> dict:
    """Load splits_metadata.json."""
    with open(path, "r") as f:
        return json.load(f)


def sample_metadata(splits: dict) -> Dict[str, dict]:
    """Build {sample_id: {superpopulation, population, sex, split}} from splits."""
    meta: Dict[str, dict] = {}
    for split_name in ("train", "val", "test"):
        for entry in splits.get(split_name, []):
            meta[entry["sample_id"]] = {
                "superpopulation": entry["superpopulation"],
                "population": entry["population"],
                "sex": entry.get("sex", 0),
                "split": split_name,
            }
    return meta


def find_vcf(pattern: str, chrom: str) -> Optional[str]:
    """Find VCF file for *chrom*, handling naming variations (e.g. chrX .v2)."""
    path = pattern.format(chrom=chrom)
    if os.path.exists(path):
        return path
    base = path.replace(".vcf.gz", "")
    hits = sorted(glob_module.glob(f"{base}*.vcf.gz"))
    return hits[0] if hits else None


def require_bcftools(console: Console) -> None:
    """Abort with a helpful message when bcftools is not on PATH."""
    try:
        subprocess.run(
            ["bcftools", "--version"], capture_output=True, check=True
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        console.print(
            "[red]Error: bcftools not found. Activate the genomics environment:[/red]"
        )
        console.print("  source ../scripts/start_genomics_universal.sh")
        sys.exit(1)


def load_snp_panel(path: str) -> Set[str]:
    """Load a SNP panel file (one rsID per line)."""
    panel: Set[str] = set()
    with open(path) as f:
        for line in f:
            tok = line.strip()
            if tok and not tok.startswith("#"):
                panel.add(tok)
    return panel


def _ensure_fd_limit(needed: int, console: Console) -> int:
    """Raise the soft file-descriptor limit; return usable batch size."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        if soft < needed:
            new_soft = min(needed, hard)
            resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft, hard))
            soft = new_soft
        if soft >= needed:
            return needed - 256
    except Exception:
        pass
    usable = max(100, soft - 256) if soft > 256 else 100
    console.print(
        f"  [yellow]File descriptor limit ({soft}) below {needed}; "
        f"using batched mode (batch={usable})[/yellow]"
    )
    return usable


def _gt_to_23andme(gt_str: str, ref: str, alt: str) -> str:
    """Fast genotype conversion for biallelic SNPs."""
    dose = GT_DOSE.get(gt_str)
    if dose is None:
        return "--"
    if dose == 0:
        return ref + ref
    if dose == 2:
        return alt + alt
    sep = "|" if "|" in gt_str else "/"
    indices = gt_str.split(sep)
    alleles = (ref, alt)
    try:
        return alleles[int(indices[0])] + alleles[int(indices[1])]
    except (ValueError, IndexError):
        return "--"


# ═══════════════════════════════════════════════════════════════
# Step 1 — Generate 23andMe files
# ═══════════════════════════════════════════════════════════════


def _process_chromosome_batch(
    vcf_path: str,
    batch_ids: List[str],
    ind_dir: Path,
    fname_tpl: str,
    inc_expr: str,
    panel: Optional[Set[str]],
) -> int:
    """Run bcftools query for *batch_ids* on one chromosome and write results."""
    n = len(batch_ids)
    sample_csv = ",".join(batch_ids)

    cmd = [
        "bcftools", "query",
        "-s", sample_csv,
        "-i", inc_expr,
        "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n",
        vcf_path,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20
    )

    fhs: Dict[str, "IO"] = {}
    for sid in batch_ids:
        fhs[sid] = open(
            ind_dir / sid / fname_tpl.format(sample_id=sid), "a"
        )

    written = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5 + n:
            continue

        rsid = parts[2]
        if panel is not None and rsid not in panel:
            continue

        ref, alt = parts[3], parts[4]
        nc = normalize_chrom(parts[0])
        prefix = f"{rsid}\t{nc}\t{parts[1]}\t"

        for i in range(n):
            geno = _gt_to_23andme(parts[5 + i], ref, alt)
            fhs[batch_ids[i]].write(prefix + geno + "\n")

        written += 1

    proc.wait()
    for fh in fhs.values():
        fh.close()
    return written


def step1_generate_23andme(config: dict, console: Console) -> None:
    """Generate 23andMe-format files from multi-sample 1000 Genomes VCFs."""
    inp = config["input"]
    conv = config.get("conversion", {})

    splits = load_splits(inp["splits_metadata"])
    meta = sample_metadata(splits)
    all_ids = list(meta.keys())
    ind_dir = Path(inp["individuals_dir"])
    fname_tpl = conv.get("output_filename", "{sample_id}_23andme.txt")

    pending = [
        s for s in all_ids
        if not (ind_dir / s / fname_tpl.format(sample_id=s)).exists()
    ]
    if not pending:
        console.print("[green]All 23andMe files already exist. Skipping.[/green]")
        return

    console.print(
        f"  Pending: {len(pending):,} / {len(all_ids):,} individuals"
    )

    panel: Optional[Set[str]] = None
    if conv.get("snp_panel"):
        panel = load_snp_panel(conv["snp_panel"])
        console.print(f"  SNP panel: {len(panel):,} SNPs")

    chroms = inp.get("chromosomes", DEFAULT_CHROMS)
    vcf_pat = inp["vcf_pattern"]

    first_vcf = next(
        (find_vcf(vcf_pat, c) for c in chroms if find_vcf(vcf_pat, c)), None
    )
    if not first_vcf:
        console.print("[red]No VCF files found.[/red]")
        return

    proc = subprocess.run(
        ["bcftools", "query", "-l", first_vcf],
        capture_output=True, text=True,
    )
    vcf_ids = set(proc.stdout.strip().split("\n"))
    valid = [s for s in pending if s in vcf_ids]
    if len(valid) < len(pending):
        console.print(
            f"  [yellow]{len(pending) - len(valid)} samples not found in VCF[/yellow]"
        )
    if not valid:
        console.print("[yellow]No valid pending samples.[/yellow]")
        return

    fmt_version = conv.get("format_version", "V5")
    header = generate_header("GRCh38", fmt_version)
    for sid in valid:
        out = ind_dir / sid / fname_tpl.format(sample_id=sid)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as f:
            f.write(header)

    skip_rsid = conv.get("skip_no_rsid", True)
    inc_parts = ['TYPE="snp"', "N_ALT=1"]
    if skip_rsid:
        inc_parts.append('ID!="."')
    inc_expr = " && ".join(inc_parts)

    batch_size = _ensure_fd_limit(len(valid) + 256, console)

    batches = [
        valid[i : i + batch_size]
        for i in range(0, len(valid), batch_size)
    ]
    n_batches = len(batches)
    if n_batches > 1:
        console.print(
            f"  Processing in {n_batches} batches of up to {batch_size}"
        )

    total_var = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as prog:
        task = prog.add_task("Step 1", total=len(chroms) * n_batches)

        for bi, batch in enumerate(batches):
            for chrom in chroms:
                vcf = find_vcf(vcf_pat, chrom)
                if not vcf:
                    prog.advance(task)
                    continue

                label = (
                    f"Step 1 — {chrom}"
                    if n_batches == 1
                    else f"Step 1 — {chrom} (batch {bi + 1}/{n_batches})"
                )
                prog.update(task, description=label)

                written = _process_chromosome_batch(
                    vcf, batch, ind_dir, fname_tpl, inc_expr, panel
                )
                if bi == 0:
                    total_var += written
                prog.advance(task)

    console.print(
        f"[green]Step 1 done — ~{total_var:,} variants/individual, "
        f"{len(valid)} files written[/green]"
    )


# ═══════════════════════════════════════════════════════════════
# Step 2 — Compute allele frequency statistics
# ═══════════════════════════════════════════════════════════════


def step2_compute_statistics(config: dict, console: Console) -> None:
    """Compute per-population allele frequencies from VCFs."""
    inp = config["input"]
    stat_cfg = config.get("statistics", {})
    pred_cfg = config.get("prediction", {})

    output_file = stat_cfg["output_file"]
    if os.path.exists(output_file):
        console.print(f"[green]Statistics file already exists: {output_file}[/green]")
        console.print("  Delete the file to recompute.")
        return

    splits = load_splits(inp["splits_metadata"])
    meta = sample_metadata(splits)

    ref_subsets = stat_cfg.get("reference_subsets", ["train"])
    ref_ids = sorted(
        sid for sid, m in meta.items() if m["split"] in ref_subsets
    )

    level = pred_cfg.get("level", "superpopulation")
    pop_map = {sid: meta[sid][level] for sid in ref_ids}
    populations = sorted(set(pop_map.values()))
    pop_idx = {p: i for i, p in enumerate(populations)}
    K = len(populations)

    pop_sizes: Dict[str, int] = defaultdict(int)
    for sid in ref_ids:
        pop_sizes[pop_map[sid]] += 1

    console.print(f"  Reference subsets: {ref_subsets} ({len(ref_ids)} individuals)")
    console.print(f"  Level: {level} ({K} classes)")
    for p in populations:
        console.print(f"    {p}: {pop_sizes[p]}")

    ind_pop_idx = {sid: pop_idx[pop_map[sid]] for sid in ref_ids}

    stat_panel: Optional[Set[str]] = None
    if stat_cfg.get("snp_panel"):
        stat_panel = load_snp_panel(stat_cfg["snp_panel"])
        console.print(f"  SNP panel filter: {len(stat_panel):,} SNPs")

    chroms = inp.get("chromosomes", DEFAULT_CHROMS)
    vcf_pat = inp["vcf_pattern"]
    min_maf = stat_cfg.get("min_maf", 0.01)
    max_snps = stat_cfg.get("max_snps", None)

    inc_expr = 'TYPE="snp" && N_ALT=1 && ID!="."'
    sample_csv = ",".join(ref_ids)
    n_ref = len(ref_ids)

    snp_ref_allele: Dict[str, str] = {}
    snp_chrom_pos: Dict[str, Tuple[str, str]] = {}
    snp_ref_count: Dict[str, List[int]] = {}
    snp_total_count: Dict[str, List[int]] = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as prog:
        task = prog.add_task("Step 2 — chromosomes", total=len(chroms))

        for chrom in chroms:
            vcf = find_vcf(vcf_pat, chrom)
            if not vcf:
                prog.advance(task)
                continue

            prog.update(task, description=f"Step 2 — {chrom}")

            cmd = [
                "bcftools", "query",
                "-s", sample_csv,
                "-i", inc_expr,
                "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n",
                vcf,
            ]
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20
            )

            for line in proc.stdout:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5 + n_ref:
                    continue

                rsid = parts[2]
                if stat_panel is not None and rsid not in stat_panel:
                    continue

                ref = parts[3]
                nc = normalize_chrom(parts[0])
                pos = parts[1]

                rc = [0] * K
                tc = [0] * K

                for i in range(n_ref):
                    dose = GT_DOSE.get(parts[5 + i])
                    if dose is None:
                        continue
                    pidx = ind_pop_idx[ref_ids[i]]
                    rc[pidx] += 2 - dose
                    tc[pidx] += 2

                snp_ref_allele[rsid] = ref
                snp_chrom_pos[rsid] = (nc, pos)
                snp_ref_count[rsid] = rc
                snp_total_count[rsid] = tc

            proc.wait()
            prog.advance(task)

    console.print(f"  Total SNPs collected: {len(snp_ref_count):,}")

    freq_data: Dict[str, List[float]] = {}
    fst_values: Dict[str, float] = {}

    for rsid in snp_ref_count:
        rc = snp_ref_count[rsid]
        tc = snp_total_count[rsid]

        freqs: List[float] = []
        tot_ref = 0
        tot_all = 0
        for k in range(K):
            f = rc[k] / tc[k] if tc[k] > 0 else 0.5
            freqs.append(round(f, 6))
            tot_ref += rc[k]
            tot_all += tc[k]

        if tot_all == 0:
            continue

        overall = tot_ref / tot_all
        maf = min(overall, 1.0 - overall)
        if maf < min_maf:
            continue

        p_bar = sum(freqs) / K
        if 0.0 < p_bar < 1.0:
            var_b = sum((f - p_bar) ** 2 for f in freqs) / K
            fst = var_b / (p_bar * (1.0 - p_bar))
        else:
            fst = 0.0

        freq_data[rsid] = freqs
        fst_values[rsid] = fst

    console.print(f"  After MAF >= {min_maf}: {len(freq_data):,} SNPs")

    if max_snps and len(freq_data) > max_snps:
        top_ids = sorted(fst_values, key=fst_values.get, reverse=True)[:max_snps]
        top_set = set(top_ids)
        freq_data = {r: freq_data[r] for r in top_ids}
        console.print(
            f"  After Fst top-{max_snps:,} selection: {len(freq_data):,} SNPs"
        )

    output = {
        "metadata": {
            "reference_subsets": ref_subsets,
            "level": level,
            "populations": populations,
            "pop_sizes": dict(pop_sizes),
            "n_individuals": len(ref_ids),
            "n_snps": len(freq_data),
            "min_maf": min_maf,
            "max_snps": max_snps,
            "created_at": datetime.now().isoformat(),
        },
        "ref_alleles": {r: snp_ref_allele[r] for r in freq_data},
        "snp_info": {r: list(snp_chrom_pos[r]) for r in freq_data},
        "allele_frequencies": freq_data,
    }

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(output, f)

    size_mb = os.path.getsize(output_file) / (1024 * 1024)
    console.print(
        f"[green]Step 2 done — {len(freq_data):,} SNPs saved to "
        f"{output_file} ({size_mb:.1f} MB)[/green]"
    )


# ═══════════════════════════════════════════════════════════════
# Step 3 — Predict ancestry
# ═══════════════════════════════════════════════════════════════


def _load_23andme_genotypes(filepath: str, needed: Set[str]) -> Dict[str, str]:
    """Load genotypes from a 23andMe file, keeping only *needed* rsIDs."""
    genos: Dict[str, str] = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[0] in needed:
                genos[parts[0]] = parts[3]
    return genos


def _genotype_dose(genotype: str, ref_allele: str) -> Optional[int]:
    """Count of *ref_allele* copies in a two-letter genotype string."""
    if genotype == "--" or len(genotype) != 2:
        return None
    return sum(1 for c in genotype if c == ref_allele)


def classify_mle(
    genotypes: Dict[str, str],
    ref_alleles: Dict[str, str],
    allele_freqs: Dict[str, List[float]],
    populations: List[str],
) -> Tuple[str, Dict[str, float], int]:
    """Maximum Likelihood Classification.

    Returns (predicted_pop, {pop: log_likelihood}, n_snps_used).
    """
    K = len(populations)
    eps = 1e-6

    valid_rsids = []
    doses: List[int] = []
    freq_rows: List[List[float]] = []

    for rsid, freqs in allele_freqs.items():
        gt = genotypes.get(rsid)
        if gt is None:
            continue
        d = _genotype_dose(gt, ref_alleles.get(rsid, ""))
        if d is None:
            continue
        valid_rsids.append(rsid)
        doses.append(d)
        freq_rows.append(freqs)

    if not valid_rsids:
        return populations[0], {p: 0.0 for p in populations}, 0

    dose_arr = np.array(doses, dtype=np.int32)
    freq_arr = np.clip(np.array(freq_rows, dtype=np.float64), eps, 1.0 - eps)

    log_p = np.log(freq_arr)
    log_1p = np.log(1.0 - freq_arr)

    log_prob = np.where(
        dose_arr[:, None] == 2,
        2.0 * log_p,
        np.where(dose_arr[:, None] == 1, np.log(2.0) + log_p + log_1p, 2.0 * log_1p),
    )

    log_liks = log_prob.sum(axis=0)
    best = int(np.argmax(log_liks))

    return (
        populations[best],
        {populations[k]: float(log_liks[k]) for k in range(K)},
        len(valid_rsids),
    )


def _admixture_nll(
    alphas_raw: np.ndarray,
    dose_arr: np.ndarray,
    freq_matrix: np.ndarray,
) -> float:
    """Negative log-likelihood for the admixture model (vectorised)."""
    alphas = alphas_raw / alphas_raw.sum()
    eps = 1e-10

    p_mixed = np.clip(freq_matrix @ alphas, eps, 1.0 - eps)

    log_prob = np.where(
        dose_arr == 2,
        2.0 * np.log(p_mixed),
        np.where(
            dose_arr == 1,
            np.log(2.0) + np.log(p_mixed) + np.log(1.0 - p_mixed),
            2.0 * np.log(1.0 - p_mixed),
        ),
    )
    return -float(np.sum(log_prob))


def estimate_admixture(
    genotypes: Dict[str, str],
    ref_alleles: Dict[str, str],
    allele_freqs: Dict[str, List[float]],
    populations: List[str],
    n_restarts: int = 20,
) -> Tuple[Dict[str, float], int]:
    """Admixture proportion estimation via constrained MLE.

    Returns ({pop: proportion}, n_snps_used).
    """
    K = len(populations)

    doses_list: List[int] = []
    freq_list: List[List[float]] = []

    for rsid, freqs in allele_freqs.items():
        gt = genotypes.get(rsid)
        if gt is None:
            continue
        d = _genotype_dose(gt, ref_alleles.get(rsid, ""))
        if d is None:
            continue
        doses_list.append(d)
        freq_list.append(freqs)

    if not doses_list:
        return {p: 1.0 / K for p in populations}, 0

    dose_arr = np.array(doses_list, dtype=np.int32)
    freq_matrix = np.array(freq_list, dtype=np.float64)

    best_nll = float("inf")
    best_alpha = np.ones(K) / K
    rng = np.random.RandomState(42)

    constraints = {"type": "eq", "fun": lambda x: x.sum() - 1.0}
    bounds = [(1e-6, 1.0)] * K

    for i in range(n_restarts):
        x0 = np.ones(K) / K if i == 0 else rng.dirichlet(np.ones(K))

        result = minimize(
            _admixture_nll,
            x0,
            args=(dose_arr, freq_matrix),
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 500, "ftol": 1e-12},
        )

        if result.fun < best_nll:
            best_nll = result.fun
            best_alpha = result.x / result.x.sum()

    return (
        {populations[k]: float(best_alpha[k]) for k in range(K)},
        len(doses_list),
    )


def step3_predict_ancestry(config: dict, console: Console) -> None:
    """Predict ancestry for evaluation individuals."""
    inp = config["input"]
    stat_cfg = config.get("statistics", {})
    pred_cfg = config.get("prediction", {})

    output_file = stat_cfg["output_file"]
    if not os.path.exists(output_file):
        console.print(
            "[red]Statistics file not found. Run Step 2 first.[/red]"
        )
        return

    with open(output_file) as f:
        stats = json.load(f)

    populations = stats["metadata"]["populations"]
    ref_alleles: Dict[str, str] = stats["ref_alleles"]
    allele_freqs: Dict[str, List[float]] = stats["allele_frequencies"]
    needed_snps = set(allele_freqs.keys())

    console.print(
        f"  Loaded {len(needed_snps):,} SNPs, "
        f"{len(populations)} {stats['metadata']['level']} classes"
    )

    splits = load_splits(inp["splits_metadata"])
    meta = sample_metadata(splits)
    eval_subsets = pred_cfg.get("evaluation_subsets", ["test"])
    eval_ids = [sid for sid, m in meta.items() if m["split"] in eval_subsets]

    ind_dir = Path(inp["individuals_dir"])
    fname_tpl = config.get("conversion", {}).get(
        "output_filename", "{sample_id}_23andme.txt"
    )

    method = pred_cfg.get("method", "mle")
    level = pred_cfg.get("level", "superpopulation")
    results_dir = Path(
        pred_cfg.get("results_dir", str(SCRIPT_DIR / "results"))
    )
    results_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"  Evaluation: {len(eval_ids)} individuals from {eval_subsets}")
    console.print(f"  Method: {method}")

    predictions: List[dict] = []

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as prog:
        task = prog.add_task("Predicting", total=len(eval_ids))

        for sid in eval_ids:
            prog.update(task, description=f"Step 3 — {sid}")

            fpath = ind_dir / sid / fname_tpl.format(sample_id=sid)
            if not fpath.exists():
                console.print(
                    f"  [yellow]File not found for {sid}, skipping[/yellow]"
                )
                prog.advance(task)
                continue

            genos = _load_23andme_genotypes(str(fpath), needed_snps)
            true_label = meta[sid][level]

            if method == "mle":
                pred, log_liks, n_used = classify_mle(
                    genos, ref_alleles, allele_freqs, populations
                )
                predictions.append(
                    {
                        "sample_id": sid,
                        "true": true_label,
                        "predicted": pred,
                        "correct": pred == true_label,
                        "snps_used": n_used,
                    }
                )

            elif method == "admixture_mle":
                n_restarts = (
                    pred_cfg.get("admixture_mle", {}).get("n_restarts", 20)
                )
                props, n_used = estimate_admixture(
                    genos, ref_alleles, allele_freqs, populations, n_restarts
                )
                pred = max(props, key=props.get)
                predictions.append(
                    {
                        "sample_id": sid,
                        "true": true_label,
                        "predicted": pred,
                        "correct": pred == true_label,
                        "snps_used": n_used,
                        "proportions": props,
                    }
                )
            else:
                console.print(f"[red]Unknown method: {method}[/red]")
                return

            prog.advance(task)

    if not predictions:
        console.print("[yellow]No predictions produced.[/yellow]")
        return

    true_labels = [p["true"] for p in predictions]
    pred_labels = [p["predicted"] for p in predictions]
    total = len(predictions)
    correct = sum(1 for p in predictions if p["correct"])
    accuracy = correct / total

    class_metrics: Dict[str, dict] = {}
    for pop in populations:
        tp = sum(1 for t, p in zip(true_labels, pred_labels) if t == pop and p == pop)
        fp = sum(1 for t, p in zip(true_labels, pred_labels) if t != pop and p == pop)
        fn = sum(1 for t, p in zip(true_labels, pred_labels) if t == pop and p != pop)
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        support = sum(1 for t in true_labels if t == pop)
        class_metrics[pop] = {
            "precision": round(prec, 4),
            "recall": round(rec, 4),
            "f1": round(f1, 4),
            "support": support,
        }

    confusion: Dict[str, Dict[str, int]] = {
        t: {p: 0 for p in populations} for t in populations
    }
    for t, p in zip(true_labels, pred_labels):
        confusion[t][p] += 1

    w_prec = sum(
        class_metrics[p]["precision"] * class_metrics[p]["support"]
        for p in populations
    ) / total
    w_rec = sum(
        class_metrics[p]["recall"] * class_metrics[p]["support"]
        for p in populations
    ) / total
    w_f1 = sum(
        class_metrics[p]["f1"] * class_metrics[p]["support"]
        for p in populations
    ) / total

    # ── Display ──────────────────────────────────────
    console.print(f"\n[bold]Results — {method}, {level}[/bold]")

    mt = Table(title="Overall Metrics")
    mt.add_column("Metric")
    mt.add_column("Value", justify="right")
    mt.add_row("Accuracy", f"{accuracy:.4f}")
    mt.add_row("Precision (weighted)", f"{w_prec:.4f}")
    mt.add_row("Recall (weighted)", f"{w_rec:.4f}")
    mt.add_row("F1 (weighted)", f"{w_f1:.4f}")
    mt.add_row("Samples", str(total))
    console.print(mt)

    ct = Table(title="Per-Class Metrics")
    ct.add_column("Class")
    ct.add_column("Precision", justify="right")
    ct.add_column("Recall", justify="right")
    ct.add_column("F1", justify="right")
    ct.add_column("Support", justify="right")
    for pop in populations:
        m = class_metrics[pop]
        ct.add_row(
            pop,
            f"{m['precision']:.4f}",
            f"{m['recall']:.4f}",
            f"{m['f1']:.4f}",
            str(m["support"]),
        )
    console.print(ct)

    cmt = Table(title="Confusion Matrix  (rows = true, columns = predicted)")
    cmt.add_column("True \\ Pred")
    for p in populations:
        cmt.add_column(p, justify="right")
    for t in populations:
        cmt.add_row(t, *(str(confusion[t][p]) for p in populations))
    console.print(cmt)

    # ── Save ─────────────────────────────────────────
    results_path = results_dir / f"predictions_{method}_{level}.json"

    results_output = {
        "metadata": {
            "method": method,
            "level": level,
            "evaluation_subsets": eval_subsets,
            "n_snps": len(needed_snps),
            "n_individuals": total,
            "created_at": datetime.now().isoformat(),
        },
        "metrics": {
            "accuracy": round(accuracy, 4),
            "weighted_precision": round(w_prec, 4),
            "weighted_recall": round(w_rec, 4),
            "weighted_f1": round(w_f1, 4),
        },
        "per_class_metrics": class_metrics,
        "confusion_matrix": confusion,
        "predictions": predictions,
    }

    with open(results_path, "w") as f:
        json.dump(results_output, f, indent=2)

    console.print(
        f"\n[green]Step 3 done — results saved to {results_path}[/green]"
    )


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════


def main() -> None:
    parser = argparse.ArgumentParser(
        description="SNP Ancestry Predictor — predict genetic ancestry from SNP data"
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to YAML configuration file",
    )
    args = parser.parse_args()

    config = load_config(Path(args.config))
    console = Console()

    console.print(
        Panel(
            "[bold]SNP Ancestry Predictor[/bold]\n"
            f"Config: {args.config}",
            title="Genomics",
            border_style="blue",
        )
    )

    require_bcftools(console)

    steps = config.get("pipeline", {}).get("steps", {})

    if steps.get("generate_23andme", True):
        console.print(
            "\n[bold blue]══════ Step 1: Generate 23andMe Files ══════[/bold blue]"
        )
        step1_generate_23andme(config, console)

    if steps.get("compute_statistics", True):
        console.print(
            "\n[bold blue]══════ Step 2: Compute Statistics ══════[/bold blue]"
        )
        step2_compute_statistics(config, console)

    if steps.get("predict_ancestry", True):
        console.print(
            "\n[bold blue]══════ Step 3: Predict Ancestry ══════[/bold blue]"
        )
        step3_predict_ancestry(config, console)

    console.print("\n[bold green]All steps completed.[/bold green]")


if __name__ == "__main__":
    main()
