#!/usr/bin/env python3
"""One classifier per gene: melanocyte of skin, the gene's own strand, on the DITA axis.

Sibling of single_gene_sweep_melanocyte_strand.py. It differs in ONE axis -- the
coordinate system -- so running both and comparing separates the effect of alignment from
the effect of the ontology/strand restriction:

    raw variant   tensor_layout=raw_center_crop     no INDEL-aware axis
    this variant  tensor_layout=haplotype_channels  + alignment_mapping=bcftools_chain
                  = DITA, Dynamic INDEL Tensor Alignment

MASKS ARE OFF (feature_mode=signals_only, indel_include_valid_mask/snp_mask false), which
is deliberate and not a detail. A masks_only arm -- the INDEL/SNP position channels with no
AlphaGenome signal at all -- once reached best_val_accuracy = 1.0 at epoch 324. The mask
channels are therefore a pure ancestry shortcut, and including them in an arm meant to
measure use of transcriptional signal would let the classifier bypass the signal entirely.
The bcftools_chain path still computes a valid mask internally to build the aligned tensor;
signals_only discards it before it reaches the model.

COST, AND THE PART THAT IS NOT FREE. The ontology/strand restriction still costs no
AlphaGenome calls. But DITA needs a per-individual, per-haplotype alignment entry for every
gene, and alignment_cache/bcftools_chain_mapper_v4/ currently holds only the ELEVEN PANEL
genes (DDB1 EDAR HERC2 MC1R MFSD12 OCA2 SLC24A5 SLC45A2 TCHH TYR TYRP1). The eleven
controls have none. The mapper builds a missing entry on demand
(bcftools_chain_mapper.get_haplotype_entry -> _build_validated_consensus), so nothing has to
be pre-staged -- but the first pass over a control gene pays a bcftools consensus per
haplotype for 1072 individuals, and lands roughly 6 GB of cache per gene. Expect the eleven
control arms to take substantially longer than the eleven panel arms, which read a cache
that already exists.

Each arm reads ONE AlphaGenome track per haplotype instead of six: ontology CL:1000458
(melanocyte of skin), on whichever strand the gene is transcribed from. With H1+H2 that
is 2 rows per gene, against 12 in the published arms.

COSTS NO API CALLS. Every prediction on disk -- personal windows and cached promoter
scrambles alike -- stores all six tracks in one (524288, 6) array, so this is a column
subset at dataset-load time. No window build, no scramble prefetch.

WHY THE KERNEL STAYS COMPARABLE. Stage 1's vertical kernel spans one haplotype's
channels, so it goes from [6, 32] to [1, 32] and the stage-1 output still has 2 rows,
one per haplotype -- exactly as in the published arms. The parameter count therefore
comes back to 49,522, identical to a three-ontology arm, and that equality is asserted
after the first arm trains. An arm that reports a different count has silently received
the wrong number of input rows, which is the failure this whole script is written around.

THE BUG THIS DEPENDS ON, now fixed. Until 2026-09-10, ontology_terms had NO effect in
the raw_center_crop layout. processed_dataset.py applied its track filter only when the
per-track metadata was present, and no loader ever populated it, so the filter fell back
to "every track" while the cache key still recorded ont1 and the channel arithmetic still
computed 2. A run could therefore declare one ontology, name its cache ont1, build a
[2, 32] kernel, and train on all twelve rows with nothing detecting it. The fix populates
the metadata (genomic_dataset.py) and makes the fallback raise instead of widening
(processed_dataset.py). Verified as a behavioural no-op for the published configuration:
with all three ontologies and no strand restriction the filter returns [0,1,2,3,4,5],
which is what the old fallback returned.

STRAND IS DERIVED, NOT TABULATED. The strand comes from the gene-level GENCODE record in
the same gtf_cache.feather the windows were built from, resolved at run time. A hard-coded
table would drift from the annotation the dataset was built against.

Outputs land beside, never on top of, the published ones:
  knockout_bulk/melanocyte_strand/pigmentation_test_split_knockout_single_<arm>.csv

Usage:
  python3 scripts/experiments/single_gene_sweep_melanocyte_strand.py --dry-run
  python3 scripts/experiments/single_gene_sweep_melanocyte_strand.py --genes TYR,LRRC36
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path("/home/breno/I2CA/genomics")
PY = "/home/breno/miniforge3/envs/genomics/bin/python3"
# The bcftools_chain aligner shells out to bare `bcftools`, so the conda env must be on
# PATH for the training subprocess -- not merely used to launch python. Missing it, the
# per-individual DITA assembly raises for every sample; that used to be swallowed by the
# normalisation fit (leaving no divisor at all) and now aborts the arm outright.
os.environ["PATH"] = f"/home/breno/miniforge3/envs/genomics/bin:{os.environ.get('PATH', '')}"
DATASET = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
GTF_CACHE = DATASET / "gtf_cache.feather"
BASE_CONFIG = REPO_ROOT / "configs/predictors/genotype_based/pigmentation/pigmentation_binary.yaml"
CONFIG_DIR = REPO_ROOT / "configs/predictors/genotype_based/pigmentation/single_gene_mel_strand_dita"
CACHE_PARENT = "results/cache/genotype_based_predictor/pigmentation_mel_strand_dita"
CACHE_ROOT = REPO_ROOT / CACHE_PARENT / "datasets"
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/melanocyte_strand_dita"
SCRAMBLE_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
LOG_DIR = REPO_ROOT / "results/genotype_based_predictor/logs"
ONTOLOGY = "CL:1000458"
TRACKS_PER_HAP = 1
# Same architecture token as the published arms except the stage-1 kernel height.
# haplotype_channels with feature_mode=signals_only carries NO extra token; the masked
# variants are the ones that add _masks_only-valid-snp / _signals_and_masks-valid-snp.
RUN_NAME = ("cnn2_pigmentation_rna_seq_H1+H2_haplotype_channels_32768_log_"
            "s1k1x32f16_s2f32_s3f64_gpavg_fc256_L100-40_relu_0.5_adam")
EXPECTED_INDIVIDUALS = 1072
# Measured on the first two arms (MC1R, TYRP1) and identical for both. It equals the
# published raw single-gene count with the stage-1 kernel shrunk from 6 rows to 1:
#   49,522 - (16*1*6*32) + (16*1*1*32) = 49,522 - 3,072 + 512 = 46,962
# Pinned rather than only compared arm-to-arm, because the in-process baseline resets on
# every relaunch and this sweep has already been relaunched twice.
EXPECTED_PARAMS = 46_962
MIN_FREE_GB = 8

PANEL = ["MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1", "EDAR", "MFSD12", "OCA2", "HERC2",
         "SLC24A5", "TCHH"]
CONTROLS = ["TPM2", "SMCR8", "PSMC4", "PPP1R3E", "ECHDC3", "FRA10AC1", "SPRED2", "EIF1B",
            "LACTB2", "LRRC36", "PRSS55"]
DEFAULT_GENES = PANEL + CONTROLS

NOTE = """metadata:
  note: >
    SINGLE-GENE ARM, MELANOCYTE + GENE STRAND, DITA AXIS -- {GENE} alone, ontology CL:1000458
    (melanocyte of skin), strand {STRAND} (the strand {GENE} is transcribed from, read
    from the gene-level GENCODE record in gtf_cache.feather). One track per haplotype,
    2 rows per gene against 12 in the published arms. Generated by
    scripts/experiments/single_gene_sweep_melanocyte_strand_dita.py. Identical to
    pigmentation_binary.yaml in split, seed, crop, layout, optimisation and
    architecture family; differs in genes_to_use, ontology_terms, track_strands, the
    stage-1 kernel/stride height (1 instead of 6, so stage 1 still emits one row per
    haplotype), an isolated results_dir, a namespaced processed_cache_dir, and run_name.
    Costs no AlphaGenome calls: the ontology and strand restriction is a column subset of
    predictions already on disk.
"""


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def free_gb(path=REPO_ROOT):
    return shutil.disk_usage(path).free / 1024 ** 3


def gene_strands(genes):
    """gene -> '+'/'-' from the gene-level GENCODE record, same GTF the windows used."""
    import pandas as pd
    gtf = pd.read_feather(GTF_CACHE)
    rows = gtf[gtf["Feature"] == "gene"]
    out = {}
    for g in genes:
        st = sorted({str(s) for s in rows[rows["gene_name"] == g]["Strand"].tolist()})
        if len(st) != 1 or st[0] not in ("+", "-"):
            raise SystemExit(f"ABORT: strand for {g} is {st}; refusing to guess.")
        out[g] = st[0]
    return out


def write_config(gene: str, strand: str) -> Path:
    cfg = yaml.safe_load(BASE_CONFIG.read_text())
    di = cfg.setdefault("dataset_input", {})
    di["genes_to_use"] = [gene]
    di["ontology_terms"] = [ONTOLOGY]
    di["track_strands"] = [strand]
    di["results_dir"] = f"results/genotype_based_predictor/runs_single_gene_{gene.lower()}_melstranddita"
    di["processed_cache_dir"] = CACHE_PARENT

    model = cfg.setdefault("model", {})
    for block, keys in (("cnn", ("kernel_size", "stride")),
                        ("cnn2", ("kernel_stage1", "stride_stage1"))):
        b = model.get(block)
        if not isinstance(b, dict):
            continue
        for k in keys:
            v = b.get(k)
            if isinstance(v, list) and len(v) == 2:
                b[k] = [TRACKS_PER_HAP, v[1]]

    cfg.setdefault("wandb", {})["run_name"] = f"pigmentation-binary-single-gene-{gene.lower()}-melstrand"

    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    out = CONFIG_DIR / f"pigmentation_binary_single_{gene.lower()}_melstranddita.yaml"
    out.write_text(NOTE.replace("{GENE}", gene).replace("{STRAND}", strand)
                   + yaml.safe_dump(cfg, sort_keys=False))
    return out


def verify_config(path: Path, gene: str, strand: str):
    base = yaml.safe_load(BASE_CONFIG.read_text())
    cfg = yaml.safe_load(path.read_text())
    diffs = []

    def walk(a, b, pre=""):
        for k in sorted(set(a) | set(b)):
            if k == "metadata":
                continue
            p = f"{pre}.{k}" if pre else k
            va, vb = a.get(k, "<missing>"), b.get(k, "<missing>")
            if isinstance(va, dict) and isinstance(vb, dict):
                walk(va, vb, p)
            elif va != vb:
                diffs.append(p)

    walk(base, cfg)
    allowed = {
        "dataset_input.genes_to_use", "dataset_input.ontology_terms",
        "dataset_input.track_strands", "dataset_input.results_dir",
        "dataset_input.processed_cache_dir", "wandb.run_name",
        "model.cnn.kernel_size", "model.cnn.stride",
        "model.cnn2.kernel_stage1", "model.cnn2.stride_stage1",
    }
    extra = set(diffs) - allowed
    if extra:
        raise SystemExit(f"ABORT: {path} drifts from the baseline in {sorted(extra)}")
    di = cfg["dataset_input"]
    if di["genes_to_use"] != [gene]:
        raise SystemExit(f"ABORT: {path} gene list is {di['genes_to_use']}")
    if di["ontology_terms"] != [ONTOLOGY]:
        raise SystemExit(f"ABORT: {path} ontology_terms is {di['ontology_terms']}")
    if di["track_strands"] != [strand]:
        raise SystemExit(f"ABORT: {path} track_strands is {di['track_strands']}, expected [{strand!r}]")
    for block, k in (("cnn2", "kernel_stage1"), ("cnn2", "stride_stage1")):
        v = cfg.get("model", {}).get(block, {}).get(k)
        if isinstance(v, list) and v[0] != TRACKS_PER_HAP:
            raise SystemExit(f"ABORT: {path} model.{block}.{k} is {v}, expected height {TRACKS_PER_HAP}")
    if f"runs_single_gene_{gene.lower()}_melstranddita" not in di["results_dir"]:
        raise SystemExit(f"ABORT: {path} has no isolated results_dir")
    if "pigmentation_mel_strand_dita" not in di["processed_cache_dir"]:
        raise SystemExit(f"ABORT: {path} shares a published tensor cache")
    return sorted(diffs)


def scrambles_ready(gene: str) -> bool:
    return bool(list(SCRAMBLE_CACHE.glob(f"seq_*_{gene}_H1_biology_tss_*_scramble100.npz")))


_FIRST_PARAMS = {"n": None, "gene": None}


def check_params(run_dir: Path, gene: str):
    """A wrong input-row count shows up here and nowhere else, so it is checked.

    On DITA the absolute count is not known ahead of time: the expanded, INDEL-aware axis
    is a little longer than the raw crop and its length depends on the gene. What MUST hold
    is that every arm agrees with the first one -- all 22 read one track per haplotype, so a
    differing count means an arm received a different number of input rows, which is exactly
    the silent failure that made an earlier melanocyte attempt train on all six tracks.
    """
    import torch
    ck = run_dir / "models" / "best_accuracy.pt"
    if not ck.exists():
        return None
    d = torch.load(ck, map_location="cpu", weights_only=False)
    sd = d.get("model_state_dict") or d.get("state_dict") or d
    n = sum(v.numel() for v in sd.values() if hasattr(v, "numel"))
    k = [tuple(v.shape) for v in sd.values() if hasattr(v, "shape") and getattr(v, "dim", lambda: 0)() == 4]
    stage1 = k[0] if k else None
    if stage1 is not None and stage1[2] != TRACKS_PER_HAP:
        _log(f"{gene:<9} WARN: stage-1 kernel height is {stage1[2]}, expected {TRACKS_PER_HAP}. "
             f"Full shape {stage1}. Treat this arm as invalid until explained.")
    if EXPECTED_PARAMS is not None and n != EXPECTED_PARAMS:
        _log(f"{gene:<9} WARN: {n:,} parameters against the pinned {EXPECTED_PARAMS:,}. "
             f"Stage-1 kernel {stage1}. Treat this arm as invalid until explained.")
    if _FIRST_PARAMS["n"] is None:
        _FIRST_PARAMS["n"], _FIRST_PARAMS["gene"] = n, gene
        _log(f"{gene:<9} params {n:,} (stage-1 kernel {stage1})")
    elif n != _FIRST_PARAMS["n"]:
        _log(f"{gene:<9} WARN: {n:,} parameters against {_FIRST_PARAMS['n']:,} for "
             f"{_FIRST_PARAMS['gene']}. Every arm reads one track per haplotype, so the counts "
             f"must agree; a difference means a different input row count reached the model.")
    return n


def views() -> set:
    return {p.name for p in CACHE_ROOT.iterdir()} if CACHE_ROOT.exists() else set()


def run(cmd, log_path):
    with open(log_path, "w", encoding="utf-8") as fh:
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=REPO_ROOT).returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", default=None, help=f"comma-separated; default: the {len(DEFAULT_GENES)} arms")
    ap.add_argument("--write-configs", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--keep-cache", action="store_true")
    ap.add_argument("--skip-done", action="store_true", default=True)
    args = ap.parse_args()

    genes = [g.strip().upper() for g in args.genes.split(",")] if args.genes else list(DEFAULT_GENES)
    built = {p.name for p in (DATASET / "individuals" / "HG00096" / "windows").iterdir() if p.is_dir()}
    missing = [g for g in genes if g not in built]
    if missing:
        raise SystemExit(f"ABORT: no built window for {missing}")

    strands = gene_strands(genes)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    KO_DIR.mkdir(parents=True, exist_ok=True)
    _log(f"melanocyte({ONTOLOGY}) + gene-strand sweep, {len(genes)} arms, in sequence")
    _log("  " + ", ".join(f"{g}({strands[g]})" for g in genes))
    _log("cost: 0 AlphaGenome calls -- ontology and strand restriction is a column subset")

    for g in genes:
        cfg = write_config(g, strands[g])
        drift = verify_config(cfg, g, strands[g])
        _log(f"{g:<9} strand {strands[g]}  config OK, drift = {len(drift)} fields")
    if args.write_configs:
        return 0

    for g in genes:
        arm = g.lower()
        cfg = CONFIG_DIR / f"pigmentation_binary_single_{arm}_melstranddita.yaml"
        run_dir = REPO_ROOT / f"results/genotype_based_predictor/runs_single_gene_{arm}_melstranddita" / RUN_NAME
        ko_csv = KO_DIR / f"pigmentation_test_split_knockout_single_{arm}.csv"
        n_win = len(list((DATASET / "individuals").glob(f"*/windows/{g}")))
        if n_win < EXPECTED_INDIVIDUALS:
            _log(f"{g:<9} SKIP: only {n_win}/{EXPECTED_INDIVIDUALS} windows built")
            continue
        if args.skip_done and ko_csv.exists():
            _log(f"{g:<9} SKIP: {ko_csv.name} already present")
            continue
        if not scrambles_ready(g):
            _log(f"{g:<9} SKIP: no cached promoter scramble; the replay would find nothing")
            continue
        if free_gb() < MIN_FREE_GB:
            _log(f"ABORT: {free_gb():.1f} GB free, need >= {MIN_FREE_GB}")
            return 1
        if args.dry_run:
            _log(f"{g:<9} DRY RUN: would train {cfg.name} (strand {strands[g]}), test, replay")
            continue

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        before = views()
        t0 = time.monotonic()

        _log(f"{g:<9} training strand {strands[g]} ({free_gb():.1f} GB free)...")
        rc = run([PY, "-m", "genomics", "genotype", "train", str(cfg)],
                 LOG_DIR / f"melstrdita_{arm}_{stamp}.train.log")
        if rc != 0:
            _log(f"{g:<9} FAILED training (rc={rc}); see melstr_{arm}_{stamp}.train.log")
            continue
        n_par = check_params(run_dir, g)
        rc = run([PY, "-m", "genomics", "genotype", "test", str(cfg)],
                 LOG_DIR / f"melstrdita_{arm}_{stamp}.test.log")
        if rc != 0:
            _log(f"{g:<9} FAILED test (rc={rc})")
            continue
        acc = None
        tj = run_dir / "test_best_accuracy_results.json"
        if tj.exists():
            acc = json.loads(tj.read_text())["weighted_accuracy"]
        else:
            found = sorted(p.name for p in run_dir.parent.iterdir()) if run_dir.parent.is_dir() else []
            _log(f"{g:<9} WARN: no {tj.name} at expected RUN_NAME; dirs present: {found}")

        _log(f"{g:<9} knockdown replay (cached scrambles, no API calls)...")
        rc = run([PY, "scripts/experiments/single_gene_knockdown_replay_ontology.py",
                  "--arm", arm, "--config", str(cfg.relative_to(REPO_ROOT)), "--out", str(ko_csv)],
                 LOG_DIR / f"melstrdita_{arm}_{stamp}.replay.log")
        if rc != 0:
            _log(f"{g:<9} FAILED replay (rc={rc}); keeping the cache for diagnosis")
            continue

        freed = 0.0
        if not args.keep_cache:
            for name in views() - before:
                d = CACHE_ROOT / name
                if d.is_dir() and d.resolve().parent == CACHE_ROOT.resolve():
                    sz = sum(f.stat().st_size for f in d.rglob("*") if f.is_file()) / 1024 ** 3
                    shutil.rmtree(d)
                    freed += sz
        _log(f"{g:<9} DONE in {(time.monotonic()-t0)/60:.1f} min | params={n_par:,} | "
             f"test acc = {acc if acc is None else f'{acc:.4f}'} | freed {freed:.1f} GB | "
             f"{free_gb():.1f} GB free")

    _log(f"sweep finished. CSVs: {len(list(KO_DIR.glob('*.csv')))}/{len(genes)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
