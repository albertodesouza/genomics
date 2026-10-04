"""The AlphaGenome output types, their resolution and how they are stored in a dataset.

Importable without the ``alphagenome`` client (the visualizer's forms and the CLI help use it).

Storage under ``individuals/<sample>/windows/<window>/predictions_<H>/``:

``tracks``       ``<output>.npz``: ``values`` (positions / resolution, tracks) float32 and
                 ``resolution`` (bp per row). 1 bp for most outputs, 128 bp for ChIP-seq.
``contact_map``  ``contact_maps.npz``: ``values`` (bins, bins, tracks) and ``resolution`` (2048).
``junctions``    ``splice_junctions.npz``: ``starts``/``ends`` (0-based window offsets, half-open),
                 ``strands`` and ``values`` (junctions, tracks).

Every file has a ``<output>_metadata.json`` with ``{"metadata": [one record per track],
"output_type", "kind", "resolution"}``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional


@dataclass(frozen=True)
class OutputSpec:
    name: str  # dna_client.OutputType name
    label: str
    group: str
    kind: str  # tracks | contact_map | junctions
    resolution: Optional[int]  # bp per stored row (None for junctions)
    description: str
    tissue_specific: bool = True

    @property
    def attr(self) -> str:
        """Attribute on ``dna_output.Output`` and the file stem on disk."""
        return self.name.lower()

    def as_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "attr": self.attr,
            "label": self.label,
            "group": self.group,
            "kind": self.kind,
            "resolution": self.resolution,
            "description": self.description,
            "tissue_specific": self.tissue_specific,
        }


OUTPUT_SPECS: Dict[str, OutputSpec] = {
    spec.name: spec
    for spec in (
        OutputSpec("RNA_SEQ", "RNA-seq", "Expression", "tracks", 1, "Gene expression (polyA+ and total RNA-seq coverage)"),
        OutputSpec("CAGE", "CAGE", "Expression", "tracks", 1, "Transcription start sites (cap analysis of gene expression)"),
        OutputSpec("PROCAP", "PRO-cap", "Expression", "tracks", 1, "Nascent transcription initiation (few cell lines)"),
        OutputSpec("DNASE", "DNase-seq", "Accessibility", "tracks", 1, "Chromatin accessibility (DNase I hypersensitivity)"),
        OutputSpec("ATAC", "ATAC-seq", "Accessibility", "tracks", 1, "Chromatin accessibility (transposase)"),
        OutputSpec("CHIP_HISTONE", "Histone ChIP-seq", "ChIP-seq", "tracks", 128, "Histone marks (H3K27ac, H3K4me3, H3K27me3, …) in 128 bp bins"),
        OutputSpec("CHIP_TF", "TF ChIP-seq", "ChIP-seq", "tracks", 128, "Transcription factor binding (CTCF, POLR2A, …) in 128 bp bins"),
        OutputSpec("SPLICE_SITES", "Splice sites", "Splicing", "tracks", 1, "Donor/acceptor probability per strand (not tissue-specific)", tissue_specific=False),
        OutputSpec("SPLICE_SITE_USAGE", "Splice site usage", "Splicing", "tracks", 1, "Fraction of transcripts using each splice site"),
        OutputSpec("SPLICE_JUNCTIONS", "Splice junctions", "Splicing", "junctions", None, "Split-read counts for predicted intron junctions"),
        OutputSpec("CONTACT_MAPS", "Contact maps", "3D genome", "contact_map", 2048, "Hi-C/Micro-C contact frequencies in 2048 bp bins (few cell lines)"),
    )
}
ALL_OUTPUTS = tuple(OUTPUT_SPECS)
TRACK_OUTPUTS = tuple(name for name, spec in OUTPUT_SPECS.items() if spec.kind == "tracks")
BASE_RESOLUTION_OUTPUTS = tuple(name for name in TRACK_OUTPUTS if OUTPUT_SPECS[name].resolution == 1)
SUPPORTED_LENGTHS = (16384, 131072, 524288, 1048576)


def normalize_outputs(values: Iterable[str]) -> List[str]:
    """Upper-case output names (``all`` expands to every output); raises ``ValueError`` on unknown names."""
    names: List[str] = []
    for value in values:
        name = str(value).strip().upper().replace("-", "_")
        if not name:
            continue
        if name == "ALL":
            names.extend(ALL_OUTPUTS)
            continue
        if name not in OUTPUT_SPECS:
            raise ValueError(f"Unknown AlphaGenome output {value!r}; choose among {', '.join(ALL_OUTPUTS)} (or 'all')")
        names.append(name)
    return list(dict.fromkeys(names))


def spec_for(name: str) -> Optional[OutputSpec]:
    return OUTPUT_SPECS.get(str(name).strip().upper().replace("-", "_"))


def stored_bytes(name: str, length: int, tracks: int) -> int:
    """Uncompressed size of one haplotype window's stored output (junctions: a rough guess)."""
    spec = OUTPUT_SPECS[name]
    if spec.kind == "junctions":
        return (length // 64) * (tracks * 4 + 17)  # ~9k candidate junctions per 524 kb window
    bins = max(1, length // (spec.resolution or 1))
    if spec.kind == "contact_map":
        return bins * bins * tracks * 4
    return bins * tracks * 4
