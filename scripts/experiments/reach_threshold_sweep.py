#!/usr/bin/env python3
"""Is `reach` a real measurement or a tautology?

The published `reach` counts an input unit as changed when the perturbation moves it by
more than 1e-6, i.e. by more than float noise. That threshold is deliberate -- reach was
built to answer "did the edit arrive at all", not "did it arrive big" -- but it makes the
statistic nearly saturated (93.8-99.9% for every gene), and a statistic that is 0.99
everywhere separates nothing. This script asks what survives if the threshold is raised.

Raising it is not one choice but three, because "changes by at least 10%" needs a
denominator, and the three defensible denominators do not agree:

  rel     |d_u| >= q * b_u          10% of THIS unit's own baseline. Degenerate where the
                                    baseline is ~0, which in log-normalised space is most
                                    of the window, so it flatters low-signal genes.
  dyn     |d_u| >= q * max_u(b_u)   10% of the gene's own dynamic range in this haplotype.
                                    Scale-aware and well-defined at b_u = 0.
  glob    |d_u| >= q * max over all genes of that individual's per-gene max
                                    Same yardstick for every gene, so it does not let a
                                    weak gene grade itself on its own curve.

`dyn` and `glob` differ exactly where it matters: `dyn` normalises away the very
magnitude difference that |Delta_in| exists to expose, so a gene delivering 59 units can
still score high on `dyn` by delivering them relative to its own small maximum. `glob` is
the one that cannot do that.

Also reported, and threshold-free: the concentration of the delivered change -- the
smallest fraction of units carrying 50% / 90% of the total |Delta|. If the movement is
spread evenly, reach and magnitude say the same thing; if it is concentrated in a handful
of units, a high reach at 1e-6 is actively misleading.

Read-only w.r.t. the AlphaGenome API: every perturbed prediction comes from the on-disk
cache, exactly as knockdown_expression_delta.py does. Rows without a cache entry are
counted and skipped, never re-billed.

Usage:
  python3 scripts/experiments/reach_threshold_sweep.py --individuals 40
  python3 scripts/experiments/reach_threshold_sweep.py            # all 162
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _path in (REPO_ROOT / "src", REPO_ROOT / "notebooks", REPO_ROOT / "scripts" / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd

DEFAULT_IN = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/pigmentation_test_split_knockout.csv"
DEFAULT_OUT = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk/reach_threshold_sweep.csv"
AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
METHOD = "biology_tss"

# The published reach threshold, then the ladder we are testing it against.
QS = [0.001, 0.01, 0.05, 0.10, 0.25, 0.50]

FIELDS = (
    ["sample_id", "gene", "haplotype", "crop_len", "abs_delta", "base_max", "base_sum",
     "reach_1e6", "frac_base_zero", "reach_rel10", "conc50", "conc90"]
    + [f"reach_dyn{q:g}" for q in QS]
    + [f"reach_glob{q:g}" for q in QS]
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _log(msg: str) -> None:
    print(f"[{_now()}] {msg}", flush=True)


def concentration(d: np.ndarray, frac: float) -> float:
    """Smallest fraction of units carrying `frac` of the total absolute change.

    Threshold-free companion to reach: reach says how many units moved, this says how few
    of them account for the movement. 1.0 means perfectly even, ~0 means a spike.
    """
    tot = d.sum()
    if tot <= 0:
        return float("nan")
    s = np.sort(d)[::-1]
    k = int(np.searchsorted(np.cumsum(s), frac * tot) + 1)
    return k / d.size


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_path", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--individuals", type=int, default=None,
                    help="subsample the first N test individuals (deterministic order)")
    args = ap.parse_args()

    import torch
    from bulk_knockout_pigmentation import _build_context
    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction, reorder_to_canonical

    device = torch.device("cpu")
    _log("Building context (cache only, no AlphaGenome calls)...")
    ctx, _ko_genes, _test_ids, _ped, _cn, _si, _wi = _build_context(device)

    df = pd.read_csv(args.in_path)
    df = df[df["method"] == METHOD].copy()
    if args.individuals:
        keep = list(dict.fromkeys(df["sample_id"]))[: args.individuals]
        df = df[df["sample_id"].isin(keep)]
    _log(f"{len(df)} rows x 2 haplotypes = {2 * len(df)} measurements "
         f"({df['sample_id'].nunique()} individuals, {df['gene'].nunique()} genes)")

    done = set()
    if args.out.exists():
        prev = pd.read_csv(args.out)
        done = set(zip(prev["sample_id"], prev["gene"], prev["haplotype"]))
        _log(f"Resuming: {len(done)} measurements already present.")
    write_header = not args.out.exists()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # Pass 1 collects per-(sample,haplotype) the max baseline over genes, which `glob` needs
    # as its yardstick; it is filled in as genes are seen and applied on the fly, so a
    # partial run still writes correct dyn/rel columns and a glob column that only becomes
    # final once every gene of that individual has been processed.
    rows_by_key: dict[tuple[str, str], list[dict]] = {}

    n_ok = n_missing = n_failed = 0
    t0 = time.monotonic()
    with open(args.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if write_header:
            w.writeheader(); fh.flush()

        for _, row in df.iterrows():
            sample_id, gene = row["sample_id"], row["gene"]
            for hap, idx_col in (("H1", "h1_target_local_idx"), ("H2", "h2_target_local_idx")):
                if (sample_id, gene, hap) in done:
                    continue
                tgt = int(row[idx_col])
                key = f"{sample_id}_{gene}_{hap}_{METHOD}_{tgt}_scramble{int(row['scramble_window'])}"
                npz, meta_json = AG_CACHE / f"seq_{key}.npz", AG_CACHE / f"seq_{key}_meta.json"
                if not (npz.exists() and meta_json.exists()):
                    n_missing += 1
                    continue
                try:
                    orig_array, orig_meta = load_raw_prediction(ctx.dataset_dir, sample_id, gene, hap)
                    canonical = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
                    mod_array = reorder_to_canonical(
                        np.load(npz)["values"], pd.DataFrame(json.loads(meta_json.read_text())), canonical)

                    # Same code path the published reach uses, so the 1e-6 column here must
                    # reproduce crop_frac_changed in knockdown_expression_delta.csv.
                    base = ctx.full_ds._process_window_haplotype_channels(
                        sample_id, gene, hap, {"rna_seq": orig_array}, {"rna_seq": orig_meta})
                    pert = ctx.full_ds._process_window_haplotype_channels(
                        sample_id, gene, hap, {"rna_seq": mod_array}, {"rna_seq": orig_meta})
                    if base is None or pert is None:
                        raise RuntimeError("processed window returned None")

                    sb = np.asarray(base[0], dtype=np.float64)
                    sp = np.asarray(pert[0], dtype=np.float64)
                    d = np.abs(sp - sb)
                    b = sb
                    if d.ndim > 1:
                        d = d.sum(axis=0); b = b.sum(axis=0)
                    n = d.size
                    bmax = float(b.max())
                    nz = b > 0

                    rec = {
                        "sample_id": sample_id, "gene": gene, "haplotype": hap,
                        "crop_len": n, "abs_delta": float(d.sum()),
                        "base_max": bmax, "base_sum": float(b.sum()),
                        "reach_1e6": float((d > 1e-6).mean()),
                        "frac_base_zero": float((~nz).mean()),
                        # per-unit relative, evaluated only where a relative change is defined
                        "reach_rel10": float((d[nz] >= 0.10 * b[nz]).mean()) if nz.any() else float("nan"),
                        "conc50": concentration(d, 0.50),
                        "conc90": concentration(d, 0.90),
                    }
                    for q in QS:
                        rec[f"reach_dyn{q:g}"] = float((d >= q * bmax).mean()) if bmax > 0 else float("nan")
                    rows_by_key.setdefault((sample_id, hap), []).append({"rec": rec, "d": d, "bmax": bmax})
                    n_ok += 1
                except Exception as exc:  # noqa: BLE001
                    n_failed += 1
                    _log(f"FAIL {sample_id}/{gene}/{hap}: {type(exc).__name__}: {exc}")

            if n_ok and n_ok % 50 == 0:
                el = time.monotonic() - t0
                _log(f"{n_ok} measurements ({n_ok/el:.2f}/s), missing={n_missing}, failed={n_failed}")

        # `glob` needs every gene of an individual+haplotype in hand before it can be
        # written, so the whole pass is held and flushed here rather than streamed.
        _log("computing the cross-gene yardstick and writing rows...")
        for (sample_id, hap), items in rows_by_key.items():
            gmax = max(it["bmax"] for it in items)
            for it in items:
                rec, d = it["rec"], it["d"]
                for q in QS:
                    rec[f"reach_glob{q:g}"] = float((d >= q * gmax).mean()) if gmax > 0 else float("nan")
                w.writerow(rec)
        fh.flush()

    _log(f"done: ok={n_ok} missing={n_missing} failed={n_failed} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
