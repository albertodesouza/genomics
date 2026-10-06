#!/usr/bin/env python3
"""Promoter knockdown on the top-10 multi-gene CNN, one gene at a time.

Why this is not ten runs of single_gene_knockdown_replay.py: in a 10-gene model the baseline
tensor is the SAME object for all ten arms, and building it means reading ten 524k x 6 raw
prediction files per individual. Ten separate processes would rebuild it ten times. Here each
individual's ten gene signals are built once, then each gene's rows are swapped for their
scrambled version in turn -- one load instead of ten, and the baseline forward pass is shared.

What this measures is RELIANCE, the panel-style quantity: how much the model's decision moves
when one of its ten inputs is scrambled while the other nine stay at baseline. That is a
different question from the single-gene arms, which measured the information one window carries.

Read-only w.r.t. the AlphaGenome API: every scrambled prediction comes from the on-disk cache,
and a missing entry is reported, never re-billed.

Correctness guard: rebuilds one individual's unperturbed tensor by hand and asserts it matches
what the dataset itself produces, before anything is written.
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
for _p in (REPO_ROOT / "src", REPO_ROOT / "notebooks"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
os.chdir(REPO_ROOT)

import numpy as np
import pandas as pd
import torch

AG_CACHE = REPO_ROOT / "notebooks" / ".cache" / "promoter_knockout_predictions"
KO_DIR = REPO_ROOT / "results/genotype_based_predictor/knockout_bulk"
LOCI_SOURCES = [
    KO_DIR / "pigmentation_test_split_knockout_no_alignment.csv",
    KO_DIR / "pigmentation_test_split_knockout_14gene.csv",
    # The null (random-window) loci. Panel genes carry method "random_window" with a
    # draw column; control genes carry "random_windowd1" with no draw column. Both
    # resolve to the same cache-key spelling through key_method() below. The two panel
    # null CSVs agree on every overlapping row (verified: 1782 rows, identical targets),
    # so their order here cannot change a result.
    KO_DIR / "pigmentation_test_split_null_scramble.csv",
    KO_DIR / "pigmentation_test_split_null_14gene.csv",
    KO_DIR / "scramble_manifest.csv",
]
DEFAULT_CFG = ("configs/predictors/genotype_based/pigmentation/"
               "pigmentation_binary_top10_multigene.yaml")

FIELDS = [
    "sample_id", "population", "superpopulation", "true_label", "gene", "n_genes_in_model",
    "method", "scramble_window", "h1_target_local_idx", "h2_target_local_idx",
    "baseline_strong_logit", "baseline_weak_logit", "perturbed_strong_logit", "perturbed_weak_logit",
    "delta_log_odds", "baseline_pred", "perturbed_pred", "flipped",
    "correct_baseline", "correct_perturbed",
    "delta_in", "delta_in_norm", "delta_expr", "delta_expr_rel", "expr_baseline",
    "expr_mod", "delta_expr_signed", "expr_log2fc",
]


def _log(m):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {m}", flush=True)


def key_method(sub: pd.DataFrame, method: str) -> str:
    """The method name as it appears in the AlphaGenome cache key.

    null_knockout_pigmentation.py writes its key as f"{METHOD}d{draw}", so a loci table
    carrying a draw column names the method WITHOUT the draw suffix and the key needs it
    appended; a table that already carries the suffixed name (the control manifest) has no
    draw column and is used verbatim. Getting this wrong does not produce a wrong number,
    it produces a cache miss, which the replay reports rather than re-billing.
    """
    if "draw" not in sub.columns:
        return method
    draws = sorted(sub["draw"].unique())
    if len(draws) != 1:
        raise SystemExit(f"ABORT: loci carry {len(draws)} draws {draws}; expected exactly one")
    return f"{method}d{int(draws[0])}"


def load_loci(genes, methods):
    """{gene: (rows, key_method)} -- where the scramble was applied, and under which
    cache-key spelling. `methods` is the ordered list of acceptable method values; the
    first one with rows for a gene wins, so one invocation can cover a gene set whose
    panel and control halves were written under different method names.
    """
    frames = {}
    for g in genes:
        hit = None
        for src in LOCI_SOURCES:
            if not src.exists():
                continue
            d = pd.read_csv(src)
            for m in methods:
                sub = d[(d["gene"] == g) & (d["method"] == m)].copy()
                if not sub.empty:
                    sub["__src"] = src.name
                    hit = (sub, key_method(sub, m), m)
                    break
            if hit:
                break
        if not hit:
            raise SystemExit(f"ABORT: no scramble loci for {g} with any of {methods}")
        frames[g] = hit
    return frames


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=DEFAULT_CFG)
    ap.add_argument("--method", default="biology_tss",
                    help="comma-separated candidate method names; the first with rows for a "
                         "gene wins. Use 'random_window,random_windowd1' for the null band, "
                         "whose panel and control halves were written under different names.")
    ap.add_argument("--out", type=Path,
                    default=KO_DIR / "top10" / "top10_knockdown.csv")
    ap.add_argument("--limit", type=int, default=None, help="first N individuals (smoke test)")
    ap.add_argument("--no-expr", action="store_true")
    a = ap.parse_args()

    from genotype_cnn_alignment_deeplift_summary.knockout import load_raw_prediction, reorder_to_canonical
    from genotype_cnn_alignment_deeplift_summary.logging_utils import quiet_pipeline_logs
    from genotype_cnn_alignment_deeplift_summary.model_loading import build_model, load_checkpoint
    from genomics.predictors.genotype_based.config import (
        generate_experiment_name, get_dataset_cache_dir, get_experiment_runs_dir, load_config)
    from genomics.predictors.genotype_based.data.pipeline import (
        _make_runtime_processed_datasets, _resolve_runtime_dataset_dir)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    cfg = load_config(REPO_ROOT / a.config)
    genes_all = list(cfg.dataset_input.genes_to_use)
    if len(genes_all) < 2:
        raise SystemExit(f"ABORT: this script is for a multi-gene model; config has {genes_all}")
    _log(f"device={device}  {len(genes_all)} genes: {genes_all}")

    exp = get_experiment_runs_dir(cfg) / generate_experiment_name(cfg)
    ck = exp / "models" / "best_accuracy.pt"
    if not ck.exists():
        raise SystemExit(f"ABORT: no checkpoint at {ck}. Train the top-10 model first.")

    methods = [m.strip() for m in a.method.split(",") if m.strip()]
    loci = load_loci(genes_all, methods)
    for g in genes_all:
        _sub, _key, _m = loci[g]
        print(f"  loci {g:9s} method={_m:16s} key={_key:16s} src={_sub['__src'].iloc[0]} "
              f"n={len(_sub)}", flush=True)

    cache_dir = get_dataset_cache_dir(cfg)
    ds_dir = _resolve_runtime_dataset_dir(cfg)
    with quiet_pipeline_logs():
        full_ds, _tr, _va, _te = _make_runtime_processed_datasets(ds_dir, cache_dir, cfg)
    model = build_model(cfg, full_ds, device)
    model = load_checkpoint(model, ck, device)
    model.eval()

    class_names = full_ds.get_class_names()
    strong_idx = [i for i, n in full_ds.idx_to_target.items() if n == "strong pigmentation"][0]
    weak_idx = next(i for i in range(2) if i != strong_idx)
    ped = full_ds.dataset_metadata.get("individuals_pedigree", {})
    dataset_dir = Path(cfg.dataset_input.dataset_dir)

    def dita_signal(sid, g, hap, arr, meta):
        r = full_ds._process_window_haplotype_channels(sid, g, hap, {"rna_seq": arr}, {"rna_seq": meta})
        if r is None:
            raise RuntimeError(f"haplotype_channels returned None for {sid}/{g}/{hap}")
        signals, masks = r
        if full_ds.feature_mode == "signals_only":
            return signals
        if full_ds.feature_mode == "masks_only":
            return masks
        return np.concatenate([signals, masks], axis=0)

    def tensor_from(rows):
        """rows[(gene, hap)] -> the model's normalized input, gene order fixed by the config."""
        stacked = np.stack([
            np.concatenate([rows[(g, "H1")] for g in genes_all], axis=0),
            np.concatenate([rows[(g, "H2")] for g in genes_all], axis=0),
        ])
        return full_ds._normalize_features_tensor(torch.FloatTensor(stacked.astype(np.float32)))

    # exon machinery for delta_expr
    exon_ctx = None
    if not a.no_expr:
        from genotype_cnn_alignment_deeplift_summary.annotations import (
            build_transcript_extractors, gene_exon_haplotype_local_boxes, load_gtf)
        _log("loading GENCODE annotations for delta_expr ...")
        _gtf = load_gtf(REPO_ROOT / "notebooks" / ".cache" / "annotations")
        _m, _pc, _exm, _exc, _gid = build_transcript_extractors(_gtf)
        exon_ctx = (_gid, _exm, _exc, gene_exon_haplotype_local_boxes)
        _log("annotations ready")

    def exon_mask(sid, g, hap, n):
        gid, exm, exc, fn = exon_ctx
        boxes, _s = fn(dataset_dir, sid, g, hap, gid, exm, exc)
        m = np.zeros(n, dtype=bool)
        for lo, hi in boxes:
            lo, hi = max(0, int(lo)), min(n, int(hi))
            if hi > lo:
                m[lo:hi] = True
        if not m.any():
            raise RuntimeError(f"no exonic positions for {g}/{sid}/{hap}")
        return m

    sid_list = [full_ds._sample_id_for_base_index(b) for b in full_ds.valid_sample_indices]
    individuals = sorted(set.intersection(*[set(t[0]["sample_id"]) for t in loci.values()]))
    individuals = [s for s in individuals if s in sid_list]
    if a.limit:
        individuals = individuals[:a.limit]
    _log(f"{len(individuals)} individuals present in both the loci and this dataset index")

    # ---- correctness guard ----
    probe = individuals[0]
    base_rows = {(g, h): dita_signal(probe, g, h, *load_raw_prediction(dataset_dir, probe, g, h))
                 for g in genes_all for h in ("H1", "H2")}
    mine = tensor_from(base_rows).squeeze()
    theirs = full_ds[sid_list.index(probe)][0].squeeze()
    if mine.shape != theirs.shape:
        raise SystemExit(f"ABORT: shape {tuple(mine.shape)} vs {tuple(theirs.shape)}")
    md = float((mine - theirs).abs().max())
    if md > 1e-4:
        raise SystemExit(f"ABORT: reconstructed baseline differs (max|diff|={md:.3e})")
    _log(f"guard OK: max|diff|={md:.2e}, tensor {tuple(mine.shape)}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if a.out.exists():
        prev = pd.read_csv(a.out)
        done = set(zip(prev["sample_id"], prev["gene"]))
    n_ok = n_missing = n_fail = 0
    t0 = time.monotonic()
    with open(a.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if not done:
            w.writeheader(); fh.flush()
        for i, sid in enumerate(individuals):
            if all((sid, g) in done for g in genes_all):
                continue
            raw = {(g, h): load_raw_prediction(dataset_dir, sid, g, h)
                   for g in genes_all for h in ("H1", "H2")}
            base_rows = {k: dita_signal(sid, k[0], k[1], *v) for k, v in raw.items()}
            base_t = tensor_from(base_rows)
            with torch.no_grad():
                bl = model(base_t.unsqueeze(0).float().to(device))[0].cpu().numpy()
            pr = ped.get(sid, {})
            true_label = full_ds._get_target_value(pr)
            bpred = class_names[int(bl.argmax())]

            for g in genes_all:
                if (sid, g) in done:
                    continue
                g_rows, g_key, _ = loci[g]
                row = g_rows[g_rows["sample_id"] == sid]
                if row.empty:
                    n_missing += 1
                    continue
                row = row.iloc[0]
                try:
                    pert_rows = dict(base_rows)
                    din = dexpr = dbase = dmod = 0.0
                    miss = False
                    for hap, col in (("H1", "h1_target_local_idx"), ("H2", "h2_target_local_idx")):
                        tgt = int(row[col])
                        key = f"{sid}_{g}_{hap}_{g_key}_{tgt}_scramble{int(row['scramble_window'])}"
                        npz, mj = AG_CACHE / f"seq_{key}.npz", AG_CACHE / f"seq_{key}_meta.json"
                        if not (npz.exists() and mj.exists()):
                            miss = True
                            break
                        orig_arr, orig_meta = raw[(g, hap)]
                        canonical = [(m["ontology_curie"], m["strand"]) for m in orig_meta]
                        mod = reorder_to_canonical(np.load(npz)["values"],
                                                   pd.DataFrame(json.loads(mj.read_text())), canonical)
                        pr_rows = dita_signal(sid, g, hap, mod, orig_meta)
                        din += float(np.abs(np.asarray(pr_rows, float)
                                            - np.asarray(base_rows[(g, hap)], float)).sum())
                        pert_rows[(g, hap)] = pr_rows
                        if exon_ctx is not None:
                            o64 = np.asarray(orig_arr, dtype=np.float64)
                            m64 = np.asarray(mod, dtype=np.float64)
                            em = exon_mask(sid, g, hap, o64.shape[0])
                            # dexpr is a MAGNITUDE (sum of |change|) and cannot say whether
                            # predicted expression fell, which is the only direction that makes
                            # this a knockdown rather than merely a perturbation. dmod keeps the
                            # perturbed exon sum so the signed change and log2 fold change are
                            # available downstream.
                            dexpr += float(np.abs(m64[em] - o64[em]).sum())
                            dbase += float(o64[em].sum())
                            dmod += float(m64[em].sum())
                    if miss:
                        n_missing += 1
                        continue
                    pert_t = tensor_from(pert_rows)
                    with torch.no_grad():
                        pl = model(pert_t.unsqueeze(0).float().to(device))[0].cpu().numpy()
                    ppred = class_names[int(pl.argmax())]
                    w.writerow({
                        "sample_id": sid, "population": pr.get("population"),
                        "superpopulation": pr.get("superpopulation"), "true_label": true_label,
                        "gene": g, "n_genes_in_model": len(genes_all), "method": g_key,
                        "scramble_window": int(row["scramble_window"]),
                        "h1_target_local_idx": int(row["h1_target_local_idx"]),
                        "h2_target_local_idx": int(row["h2_target_local_idx"]),
                        "baseline_strong_logit": float(bl[strong_idx]),
                        "baseline_weak_logit": float(bl[weak_idx]),
                        "perturbed_strong_logit": float(pl[strong_idx]),
                        "perturbed_weak_logit": float(pl[weak_idx]),
                        "delta_log_odds": float((pl[strong_idx] - pl[weak_idx])
                                                - (bl[strong_idx] - bl[weak_idx])),
                        "baseline_pred": bpred, "perturbed_pred": ppred,
                        "flipped": int(bpred != ppred),
                        "correct_baseline": int(bpred == true_label),
                        "correct_perturbed": int(ppred == true_label),
                        "delta_in": din,
                        "delta_in_norm": float((pert_t - base_t).abs().sum()),
                        "delta_expr": dexpr if exon_ctx is not None else float("nan"),
                        "expr_baseline": dbase if exon_ctx is not None else float("nan"),
                        "delta_expr_rel": (dexpr / dbase) if (exon_ctx is not None and dbase > 0)
                                          else float("nan"),
                        "expr_mod": dmod if exon_ctx is not None else float("nan"),
                        "delta_expr_signed": (dmod - dbase) if exon_ctx is not None else float("nan"),
                        "expr_log2fc": (float(np.log2(dmod / dbase))
                                        if (exon_ctx is not None and dbase > 0 and dmod > 0)
                                        else float("nan")),
                    })
                    n_ok += 1
                except Exception as exc:  # noqa: BLE001
                    n_fail += 1
                    _log(f"FAILED {sid}/{g}: {exc!r}")
            fh.flush()
            if (i + 1) % 10 == 0:
                el = time.monotonic() - t0
                _log(f"{i+1}/{len(individuals)} individuals, {n_ok} rows, {el/60:.1f} min "
                     f"(ETA {(el/(i+1)*(len(individuals)-i-1))/60:.1f} min), "
                     f"missing={n_missing} failed={n_fail}")
    _log(f"DONE: ok={n_ok} missing={n_missing} failed={n_fail} -> {a.out}")
    if n_missing or n_fail:
        _log("WARNING: incomplete -- do not table these numbers as-is.")
    return 0 if (n_fail == 0 and n_missing == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
