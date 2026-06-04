from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List, Optional

import pandas as pd


OUTPUT_METADATA_ATTRS = [
    "rna_seq",
    "cage",
    "dnase",
    "atac",
    "chip_histone",
    "chip_tf",
    "splice_sites",
    "splice_junctions",
    "splice_site_usage",
    "procap",
    "contact_maps",
]


def export_tracks(api_key: str, output: Path) -> pd.DataFrame:
    from alphagenome.models import dna_client as alphagenome_client

    alphagenome_model = alphagenome_client.create(api_key=api_key)
    metadata = alphagenome_model.output_metadata(organism=alphagenome_client.Organism.HOMO_SAPIENS)

    frames = []
    for attr in OUTPUT_METADATA_ATTRS:
        if hasattr(metadata, attr):
            df = getattr(metadata, attr).copy()
            df["output_type"] = attr.upper()
            frames.append(df)

    if not frames:
        raise RuntimeError("No AlphaGenome output metadata frames were returned")

    all_tracks = pd.concat(frames, ignore_index=True)
    all_tracks.to_csv(output, index=False)
    return all_tracks


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Export AlphaGenome output metadata tracks to CSV")
    parser.add_argument("--api-key", default=os.environ.get("ALPHAGENOME_API_KEY"))
    parser.add_argument("--output", type=Path, default=Path("alpha_genome_all_tracks.csv"))
    args = parser.parse_args(argv)

    if not args.api_key:
        parser.error("--api-key is required or set ALPHAGENOME_API_KEY")

    tracks = export_tracks(args.api_key, args.output)
    print(f"Saved {len(tracks)} AlphaGenome tracks to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
