#!/usr/bin/env python3
"""Expose the third random draw's windows inside the working dataset, by symlink.

WHAT THIS IS FOR, AND WHAT IT IS NOT. The eleven genes of the third pre-existing random
draw (`random_11_3` in spirit, though its directory carries no suffix) already have
AlphaGenome predictions on disk from the non-longevous-dataset project. Their input
sequences are byte-identical to what this project's pipeline would send -- verified by
md5 on `*.window.fixed.fa` -- but the predictions are from a DIFFERENT CALL: AlphaGenome
is not reproducible call to call, and on a spot check the two differ at 944,619 of
3,145,728 positions, up to 0.125 in absolute value, about 0.1% in per-track sums.

So the arms built on these links are a PREVIEW. They answer "what does the ranking do at
33 controls" for no API spend. They are NOT the definitive arms: those need
control_expansion-style rebuild with fresh calls, so that all 42 arms come from one
prediction vintage. Nothing here should be substituted into the paper.

WHY SYMLINKS AND NOT COPIES. Reading is all training does with a window directory: the
bcftools alignment cache lives at the dataset root (`alignment_cache/`), not inside the
per-gene window, so nothing is written back through the link. Symlinks cost no disk where
copies would cost ~50 GB, and `--unlink` removes every one of them, leaving the working
dataset exactly as it was. The source tree belongs to another user and is only ever read.

Idempotent: an existing correct link is left alone, and a path that exists and is NOT a
link to the expected target is reported and skipped rather than replaced.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

WORK = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")
SRC = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes_1000_random")
GENES = ["ATP11B", "BCL3", "C6orf52", "FBXO5", "FOXN2", "HERC6", "KIAA0319", "RIDA",
         "SEM1", "SFMBT2", "SUMF2"]
# Any gene already measured: the individuals that have it are exactly the cohort the
# paper's arms are trained on, so the link set is defined by what is already there
# rather than by a hardcoded count that could drift.
REFERENCE_GENE = "CD47"
MANIFEST = Path("/home/breno/I2CA/genomics/results/genotype_based_predictor/"
                "draw3_preview_links.json")


def cohort() -> list[str]:
    ids = sorted(d.name for d in (WORK / "individuals").iterdir()
                 if (d / "windows" / REFERENCE_GENE).exists())
    if not ids:
        raise SystemExit(f"ABORT: no individual carries {REFERENCE_GENE}; wrong dataset?")
    return ids


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unlink", action="store_true", help="remove the links this made")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    ids = cohort()
    print(f"cohort: {len(ids)} individuals carrying {REFERENCE_GENE}")

    if a.unlink:
        removed = kept = 0
        for s in ids:
            for g in GENES:
                p = WORK / "individuals" / s / "windows" / g
                if p.is_symlink():
                    if not a.dry_run:
                        p.unlink()
                    removed += 1
                elif p.exists():
                    kept += 1
        print(f"removed {removed} symlinks; left {kept} real directories untouched")
        if not a.dry_run and MANIFEST.exists():
            MANIFEST.unlink()
        return 0

    made = ok = missing = conflict = 0
    for s in ids:
        for g in GENES:
            src = SRC / "individuals" / s / "windows" / g
            dst = WORK / "individuals" / s / "windows" / g
            if not src.is_dir():
                missing += 1
                continue
            if dst.is_symlink():
                ok += 1 if dst.resolve() == src.resolve() else 0
                conflict += 0 if dst.resolve() == src.resolve() else 1
                continue
            if dst.exists():
                # A real directory here means the definitive build already ran for this
                # gene. Never shadow it with a preview link.
                conflict += 1
                continue
            if not a.dry_run:
                dst.symlink_to(src)
            made += 1
    print(f"links made {made}, already correct {ok}, source missing {missing}, "
          f"conflicts left alone {conflict}")
    if conflict:
        print("  WARNING: conflicts are paths that exist but are not the expected link; "
              "inspect before trusting any arm built over them")

    # Completeness per gene, judged from the working tree the trainer will read, not from
    # the source: a link that failed to appear must not be discovered at training time.
    print("\nper gene, visible in the working dataset:")
    bad = []
    for g in GENES:
        n = sum(1 for s in ids if (WORK / "individuals" / s / "windows" / g).exists())
        flag = "" if n == len(ids) else "  <-- INCOMPLETE"
        print(f"  {g:9s} {n}/{len(ids)}{flag}")
        if n != len(ids):
            bad.append(g)
    if bad and not a.dry_run:
        print(f"\nABORT: incomplete genes {bad}; do not train over these")
        return 1

    if not a.dry_run:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_text(json.dumps({
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "purpose": "preview arms for the third random draw; predictions are a "
                       "different AlphaGenome call vintage than the other 31 arms",
            "source": str(SRC), "work": str(WORK), "genes": GENES,
            "n_individuals": len(ids), "links": made + ok,
            "remove_with": "python3 scripts/experiments/link_draw3_preview.py --unlink",
        }, indent=1))
        print(f"\nwrote {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
