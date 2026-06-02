from __future__ import annotations

import argparse
import json
from pathlib import Path

from genomics_workspace import DEFAULT_DATASET_DIR, LEGACY_TOP3_DATASET_DIR
from genotype_based_predictor.data.layout import sync_bcftools_chain_artifacts


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Sincroniza artefatos legacy necessarios para alignment_mapping=bcftools_chain no dataset canonico."
    )
    parser.add_argument("--source-dir", type=Path, default=LEGACY_TOP3_DATASET_DIR)
    parser.add_argument("--target-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument("--link-mode", choices=["hardlink", "symlink", "copy"], default="hardlink")
    parser.add_argument("--apply", action="store_true", help="Efetiva links/copias. Sem esta flag roda em dry-run.")
    parser.add_argument("--sample-limit", type=int, default=None)
    parser.add_argument("--genes", nargs="*", default=None)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    stats = sync_bcftools_chain_artifacts(
        args.source_dir,
        args.target_dir,
        link_mode=args.link_mode,
        dry_run=not args.apply,
        sample_limit=args.sample_limit,
        genes=args.genes,
    )
    if args.json:
        print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
