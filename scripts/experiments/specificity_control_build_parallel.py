#!/usr/bin/env python3
"""Parallel cohort driver for the specificity-control window build.

`build_non_longevous_dataset.py` honours `pipeline.parallel.n_workers` for some
steps but processes samples in step 4 (`run_predictions`) strictly serially, at
roughly 2.5 min/sample here. Over the 1072 individuals of the pigmentation
cohort that is ~45 hours, almost all of it spent waiting on AlphaGenome rather
than on CPU.

This driver keeps the same per-sample work -- it calls the builder's own
`run_build_window_predict`, `IndividualDatasetBuilder` and window-organising
code, so the on-disk result is byte-for-byte what the serial path produces --
and runs samples across a process pool. Per-sample outputs are written to
disjoint directories, so the only shared state is the checkpoint, which the
parent owns: workers return results and never write it.

Shared reference windows (references/windows/<gene>/) are created once up front,
before any worker starts, so no two workers race to write them.

Resumable: samples already in the checkpoint are skipped, so this can take over
from an interrupted serial run and vice versa.

Usage:
  python3 scripts/experiments/specificity_control_build_parallel.py \
      --config configs/workflows/non_longevous_dataset/specificity_control_genes.yaml \
      --workers 8
  ... --limit 4 --workers 2          # smoke test
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path("/home/breno/I2CA/genomics")
for _p in (REPO_ROOT / "src",):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
os.chdir(REPO_ROOT)

from genomics.workflows.dataset_builders.non_longevous.build_non_longevous_dataset import (  # noqa: E402
    load_checkpoint,
    load_config,
    load_metadata_csv,
    parse_fasta_header,
    run_build_window_predict,
    save_checkpoint,
    select_samples,
)
from genomics.workflows.dataset_builders.non_longevous.dataset_builder import (  # noqa: E402
    IndividualDatasetBuilder,
)


def _log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc).isoformat()}] {msg}", flush=True)


def _process_one(sample_row: dict, config: dict, output_dir_str: str):
    """Build all configured windows for one sample. Runs in a worker process."""
    output_dir = Path(output_dir_str)
    sample_id = sample_row["SampleID"]
    params = config["build_window_params"]

    try:
        sample_info = {
            "FamilyID": sample_row.get("FamilyID", "0"),
            "SampleID": sample_id,
            "Sex": int(sample_row["Sex"]),
            "Population": sample_row["Population"],
            "Superpopulation": sample_row["Superpopulation"],
        }
        # FROG likelihoods are omitted here: they are already recorded in each
        # individual's metadata from the original eleven-window build, and
        # IndividualDatasetBuilder preserves existing fields it is not given.
        ind_builder = IndividualDatasetBuilder(
            base_dir=output_dir, sample_id=sample_id, sample_info=sample_info,
        )
        ind_builder.create_structure()

        success, _target = run_build_window_predict(sample_id, config, output_dir)
        if not success:
            return sample_id, False, "run_build_window_predict returned failure"

        windows_dir = output_dir / "individuals" / sample_id / "windows"
        if windows_dir.exists():
            mode = params.get("mode", "gene")
            outputs_str = params.get("outputs", "") or ""
            outputs = [o.strip() for o in outputs_str.split(",") if o.strip()]
            ontology_str = params.get("ontology", "") or ""
            ontologies = [o.strip() for o in ontology_str.split(",") if o.strip()]

            for window_dir in sorted(d for d in windows_dir.iterdir() if d.is_dir()):
                window_name = window_dir.name
                ref_fasta = output_dir / "references" / "windows" / window_name / "ref.window.fa"
                chromosome, start, end = "unknown", 0, params.get("window_size", 1000000)
                if ref_fasta.exists():
                    c, s, e = parse_fasta_header(ref_fasta)
                    if c is not None:
                        chromosome, start, end = c, s, e
                ind_builder.add_window(
                    target_name=window_name, window_type=mode, chromosome=chromosome,
                    start=start, end=end, outputs=outputs, ontologies=ontologies,
                )
        ind_builder.save_metadata()
        return sample_id, True, None
    except Exception as exc:                                          # noqa: BLE001
        return sample_id, False, f"{type(exc).__name__}: {exc}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--workers", type=int, default=8,
                    help="Parallel samples in flight. The binding resource is the "
                         "AlphaGenome API, not CPU (default: 8)")
    ap.add_argument("--limit", type=int, help="Process at most this many pending samples (smoke test)")
    args = ap.parse_args()

    config = load_config(Path(args.config).resolve())
    output_dir = Path(config["project"]["output_dir"]).resolve()
    checkpoint_file = output_dir / config["pipeline"]["checkpoint_file"]

    df = load_metadata_csv(Path(config["data_sources"]["metadata_csv"]).resolve())
    selected = select_samples(df, config)
    _log(f"Cohort: {len(selected)} individuals")

    checkpoint = load_checkpoint(checkpoint_file)
    done = set(checkpoint["completed_samples"])
    pending = [r for _, r in selected.iterrows() if r["SampleID"] not in done]
    if args.limit:
        pending = pending[: args.limit]
    _log(f"Already complete: {len(done)}  |  pending: {len(pending)}  |  workers: {args.workers}")
    if not pending:
        _log("Nothing to do.")
        return 0

    # Create the shared reference windows once, serially, so workers never race
    # on references/windows/<gene>/. Processing a single sample is enough.
    genes_needed = []
    gene_list_file = config["build_window_params"].get("gene", {}).get("gene_list_file")
    if gene_list_file:
        genes_needed = [g.strip() for g in Path(gene_list_file).read_text(encoding="utf-8").split() if g.strip()]
    missing = [g for g in genes_needed
               if not (output_dir / "references" / "windows" / g / "ref.window.fa").exists()]
    if missing:
        _log(f"Reference windows missing for {missing}; building them via one serial sample first")
        seed = pending.pop(0)
        sid, ok, err = _process_one(dict(seed), config, str(output_dir))
        if ok:
            checkpoint["completed_samples"].append(sid)
        else:
            checkpoint["failed_samples"].append(sid)
            _log(f"Seed sample {sid} FAILED: {err}")
        save_checkpoint(checkpoint_file, checkpoint)

    n_ok = n_fail = 0
    total = len(pending)
    started = datetime.now(timezone.utc)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_process_one, dict(r), config, str(output_dir)): r["SampleID"]
                   for r in pending}
        for i, fut in enumerate(as_completed(futures), 1):
            sample_id, ok, err = fut.result()
            if ok:
                checkpoint["completed_samples"].append(sample_id)
                n_ok += 1
            else:
                checkpoint["failed_samples"].append(sample_id)
                n_fail += 1
                _log(f"  {sample_id} FAILED: {err}")
            save_checkpoint(checkpoint_file, checkpoint)

            if i % 10 == 0 or i == total:
                elapsed = (datetime.now(timezone.utc) - started).total_seconds()
                rate = i / elapsed if elapsed else 0.0
                eta_h = (total - i) / rate / 3600 if rate else float("nan")
                _log(f"{i}/{total} done ({n_ok} ok, {n_fail} failed) "
                     f"| {rate*60:.1f} samples/min | ETA {eta_h:.1f} h")

    _log(f"Finished: {n_ok} succeeded, {n_fail} failed. Checkpoint: {checkpoint_file}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
