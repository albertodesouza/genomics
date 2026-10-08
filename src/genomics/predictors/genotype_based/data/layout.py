"""Utilities to validate and materialize the on-disk dataset layout."""

from __future__ import annotations

import errno
import filecmp
import json
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

console = Console()


LAYOUT_VERSION = 1


def is_dataset_dir(dataset_dir: Path) -> bool:
    marker = Path(dataset_dir) / "layout_metadata.json"
    if not marker.exists():
        return False
    try:
        with open(marker) as f:
            meta = json.load(f)
        return meta.get("layout_version") == LAYOUT_VERSION
    except Exception:
        return False


def _load_json(path: Path) -> Dict:
    with open(path) as f:
        return json.load(f)


def _write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def _same_content(a: Path, b: Path) -> bool:
    return filecmp.cmp(a, b, shallow=False)


def _move_file(src: Path, dst: Path) -> None:
    """Moves ``src`` to ``dst`` so that ``dst`` only ever appears complete.

    On the same filesystem this is a rename (no extra space). Across
    filesystems the file is copied to a temporary name, renamed into place and
    only then removed from the source, so the extra space is one file at a time
    and an interrupted run never leaves a truncated ``dst``.
    """
    try:
        os.replace(src, dst)
        return
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
    tmp = dst.with_name(dst.name + ".partial")
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)
    src.unlink()


def _transfer_file(src: Path, dst: Path, *, move: bool, stats: Counter) -> None:
    """Copies (default) or moves ``src`` to ``dst``; an existing ``dst`` is kept.

    When moving and ``dst`` already exists, an identical ``src`` (e.g. the
    per-individual copies of a window's ``ref.window.fa``) is deleted so the
    source shrinks as the target grows; a differing ``src`` is left in place
    and counted as a conflict.
    """
    if dst.exists():
        if move:
            if _same_content(src, dst):
                src.unlink()
                stats["duplicates_removed"] += 1
            else:
                stats["conflicts"] += 1
                console.print(f"[yellow]Conflito, origem mantida:[/yellow] {src} difere de {dst}")
        else:
            stats["already_present"] += 1
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if move:
        _move_file(src, dst)
        stats["moved"] += 1
    else:
        shutil.copy2(src, dst)
        stats["copied"] += 1


def _materialize_file(src: Path, dst: Path, *, link_mode: str, dry_run: bool) -> str:
    if dst.exists():
        return "exists"
    if dry_run:
        return "missing"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if link_mode == "hardlink":
        try:
            os.link(src, dst)
            return "linked"
        except OSError:
            shutil.copy2(src, dst)
            return "copied"
    if link_mode == "symlink":
        dst.symlink_to(src)
        return "linked"
    shutil.copy2(src, dst)
    return "copied"


def _format_vcf_pattern(vcf_pattern: str, chromosome: str) -> str:
    try:
        return vcf_pattern.format(chrom=chromosome)
    except (KeyError, IndexError, ValueError):
        return vcf_pattern.replace("{chrom}", chromosome)


def _candidate_chromosomes(chromosome: str) -> Iterable[str]:
    chromosome = str(chromosome)
    yield chromosome
    if chromosome.startswith("chr"):
        yield chromosome[3:]
    else:
        yield f"chr{chromosome}"


def _resolve_vcf_path(raw_variant_sources: Dict, chromosome: str) -> Optional[str]:
    vcf_pattern = raw_variant_sources.get("vcf_pattern") or os.environ.get("KG1000_VCF_PATTERN")
    if vcf_pattern:
        candidates = [_format_vcf_pattern(str(vcf_pattern), chrom) for chrom in _candidate_chromosomes(chromosome)]
        for candidate in candidates:
            if Path(candidate).exists():
                return candidate
        return candidates[0]

    vcf_root_dir = raw_variant_sources.get("vcf_root_dir") or os.environ.get("KG1000_VCF_ROOT_DIR")
    if not vcf_root_dir:
        return None

    root = Path(vcf_root_dir)
    for chrom in _candidate_chromosomes(chromosome):
        candidates = [
            root / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
            root / f"{chrom}.vcf.gz",
        ]
        for candidate in candidates:
            if candidate.exists():
                return str(candidate)
    return None


def update_window_variant_sources(dataset_dir: Path, vcf_pattern: Optional[str] = None, vcf_root_dir: Optional[str] = None) -> None:
    """Adds per-window VCF paths required by the indel alignment code."""
    dataset_dir = Path(dataset_dir)
    dataset_meta_path = dataset_dir / "dataset_metadata.json"
    dataset_meta = _load_json(dataset_meta_path)
    raw_variant_sources = dict(dataset_meta.get("raw_variant_sources", {}) or {})
    if vcf_pattern:
        raw_variant_sources["vcf_pattern"] = vcf_pattern
    else:
        raw_variant_sources.setdefault("vcf_pattern", os.environ.get("KG1000_VCF_PATTERN"))
    if vcf_root_dir:
        raw_variant_sources["vcf_root_dir"] = vcf_root_dir
    else:
        raw_variant_sources.setdefault("vcf_root_dir", os.environ.get("KG1000_VCF_ROOT_DIR"))
    raw_variant_sources.setdefault(
        "notes",
        "Set KG1000_VCF_PATTERN or KG1000_VCF_ROOT_DIR to preserve the chromosome-level 1000G VCF source used to generate consensus FASTAs.",
    )
    dataset_meta["raw_variant_sources"] = raw_variant_sources

    window_catalog = dict(dataset_meta.get("window_catalog", {}) or {})
    ref_windows_dir = dataset_dir / "references" / "windows"
    if not window_catalog and ref_windows_dir.exists():
        for meta_path in sorted(ref_windows_dir.glob("*/window_metadata.json")):
            window_catalog[meta_path.parent.name] = _load_json(meta_path)

    for window_name, meta in window_catalog.items():
        chromosome = str(meta.get("chromosome", ""))
        raw_variant_source = dict(meta.get("raw_variant_source", {}) or {})
        raw_variant_source["chromosome"] = chromosome
        raw_variant_source["vcf_pattern"] = raw_variant_sources.get("vcf_pattern")
        raw_variant_source["vcf_root_dir"] = raw_variant_sources.get("vcf_root_dir")
        vcf_path = _resolve_vcf_path(raw_variant_sources, chromosome) if chromosome else None
        if vcf_path:
            raw_variant_source["vcf_path"] = vcf_path
        meta["raw_variant_source"] = raw_variant_source
        _write_json(ref_windows_dir / window_name / "window_metadata.json", meta)

    dataset_meta["window_catalog"] = window_catalog
    _write_json(dataset_meta_path, dataset_meta)


def _iter_individual_dirs(dataset_dir: Path) -> Iterable[Path]:
    individuals_dir = dataset_dir / "individuals"
    if not individuals_dir.exists():
        return []
    return sorted([p for p in individuals_dir.iterdir() if p.is_dir()])


# Dataset-level files carried over when present in the source.
_DATASET_LEVEL_FILES = (
    "gtf_cache.feather",
    "metadata_statistics.json",
    "selected_samples.csv",
    "vcf_validation_report.json",
)


def materialize_dataset(source_dir: Path, target_dir: Path, *, move: bool = False) -> Path:
    """Creates the normalized dataset view consumed by this package.

    With ``move=True`` the data files are moved out of ``source_dir`` instead
    of copied, so the conversion needs almost no extra disk space when source
    and target share a filesystem. The source's metadata JSONs are only read,
    never moved, so an interrupted run can simply be started again; what is
    left in the source afterwards is metadata, ``prediction_H*.ok.txt``
    markers and any conflicting files reported at the end.
    """
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)

    layout_marker = target_dir / "layout_metadata.json"
    if layout_marker.exists() and is_dataset_dir(target_dir):
        console.print(f"[green]Dataset já existe:[/green] {target_dir}")
        return target_dir

    dataset_meta = _load_json(source_dir / "dataset_metadata.json")
    target_dir.mkdir(parents=True, exist_ok=True)
    stats: Counter = Counter()

    for file_name in _DATASET_LEVEL_FILES:
        src = source_dir / file_name
        if src.exists():
            _transfer_file(src, target_dir / file_name, move=move, stats=stats)

    ref_windows_dir = target_dir / "references" / "windows"
    individuals_out_dir = target_dir / "individuals"

    global_window_metadata: Dict[str, Dict] = {}
    gene_order = list(dataset_meta.get("genes", []))

    individual_dirs = list(_iter_individual_dirs(source_dir))
    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Materializando individuos...", total=len(individual_dirs))
        for ind_dir in individual_dirs:
            sample_id = ind_dir.name
            ind_meta_path = ind_dir / "individual_metadata.json"
            if not ind_meta_path.exists():
                progress.update(task, advance=1)
                continue

            ind_meta = _load_json(ind_meta_path)
            windows = list(ind_meta.get("windows", []))
            window_metadata = ind_meta.get("window_metadata", {})

            new_ind_meta = {
                "sample_id": ind_meta.get("sample_id", sample_id),
                "family_id": ind_meta.get("family_id", "0"),
                "sex": ind_meta.get("sex", 0),
                "population": ind_meta.get("population", ""),
                "superpopulation": ind_meta.get("superpopulation", ""),
                "longevity": ind_meta.get("longevity", False),
                "windows": windows,
                "created_at": ind_meta.get("created_at"),
                "last_updated": ind_meta.get("last_updated"),
            }
            if "frog_likelihood" in ind_meta:
                new_ind_meta["frog_likelihood"] = ind_meta["frog_likelihood"]
            if "frog_population_names" in ind_meta:
                new_ind_meta["frog_population_names"] = ind_meta["frog_population_names"]

            _write_json(individuals_out_dir / sample_id / "individual_metadata.json", new_ind_meta)

            for window_name in windows:
                src_window_dir = ind_dir / "windows" / window_name
                if not src_window_dir.exists():
                    continue

                if window_name in window_metadata and window_name not in global_window_metadata:
                    global_window_metadata[window_name] = dict(window_metadata[window_name])
                    _write_json(ref_windows_dir / window_name / "window_metadata.json", global_window_metadata[window_name])

                ref_src = src_window_dir / "ref.window.fa"
                if ref_src.exists():
                    _transfer_file(ref_src, ref_windows_dir / window_name / "ref.window.fa", move=move, stats=stats)

                dst_window_dir = individuals_out_dir / sample_id / "windows" / window_name
                # Haplotype FASTAs plus the raw per-window VCF and bcftools
                # chain artefacts (raw.fa, consensus_ready VCF) when available.
                for file_name in (
                    f"{sample_id}.H1.window.fixed.fa",
                    f"{sample_id}.H2.window.fixed.fa",
                    f"{sample_id}.H1.window.raw.fa",
                    f"{sample_id}.H2.window.raw.fa",
                    f"{sample_id}.window.vcf.gz",
                    f"{sample_id}.window.vcf.gz.tbi",
                    f"{sample_id}.window.consensus_ready.vcf.gz",
                    f"{sample_id}.window.consensus_ready.vcf.gz.tbi",
                ):
                    src = src_window_dir / file_name
                    if src.exists():
                        _transfer_file(src, dst_window_dir / file_name, move=move, stats=stats)

                for pred_dir_name in ("predictions_H1", "predictions_H2"):
                    src_pred_dir = src_window_dir / pred_dir_name
                    dst_pred_dir = dst_window_dir / pred_dir_name
                    if not src_pred_dir.exists():
                        continue
                    dst_pred_dir.mkdir(parents=True, exist_ok=True)
                    for pred_file in sorted(src_pred_dir.iterdir()):
                        if pred_file.is_file():
                            _transfer_file(pred_file, dst_pred_dir / pred_file.name, move=move, stats=stats)

            progress.update(task, advance=1)

    out_meta = dict(dataset_meta)
    out_meta["layout_version"] = LAYOUT_VERSION
    out_meta["source_dataset_dir"] = str(source_dir.resolve())
    out_meta["window_catalog"] = global_window_metadata
    out_meta.setdefault("raw_variant_sources", {
        "vcf_pattern": os.environ.get("KG1000_VCF_PATTERN"),
        "vcf_root_dir": os.environ.get("KG1000_VCF_ROOT_DIR"),
        "notes": "Set KG1000_VCF_PATTERN or KG1000_VCF_ROOT_DIR to preserve the chromosome-level 1000G VCF source used to generate consensus FASTAs.",
    })
    if gene_order:
        out_meta["genes"] = gene_order

    for window_name, meta in global_window_metadata.items():
        chromosome = str(meta.get("chromosome", ""))
        chrom_token = chromosome.replace("chr", "") if chromosome else None
        if chrom_token:
            raw_variant_source = dict(meta.get("raw_variant_source", {}) or {})
            raw_variant_source.setdefault("chromosome", chromosome)
            raw_variant_source.setdefault("vcf_pattern", out_meta.get("raw_variant_sources", {}).get("vcf_pattern"))
            raw_variant_source.setdefault("vcf_root_dir", out_meta.get("raw_variant_sources", {}).get("vcf_root_dir"))
            raw_variant_source.setdefault("vcf_pattern_example", None if not out_meta.get("raw_variant_sources", {}).get("vcf_pattern") else _format_vcf_pattern(out_meta["raw_variant_sources"]["vcf_pattern"], chrom_token))
            vcf_path = _resolve_vcf_path(out_meta.get("raw_variant_sources", {}), chromosome)
            if vcf_path:
                raw_variant_source.setdefault("vcf_path", vcf_path)
            meta["raw_variant_source"] = raw_variant_source
            _write_json(ref_windows_dir / window_name / "window_metadata.json", meta)
    _write_json(target_dir / "dataset_metadata.json", out_meta)
    _write_json(layout_marker, {
        "layout_version": LAYOUT_VERSION,
        "source_dataset_dir": str(source_dir.resolve()),
        "target_dataset_dir": str(target_dir.resolve()),
    })

    summary = ", ".join(f"{key}={stats[key]}" for key in sorted(stats)) or "nenhum arquivo transferido"
    console.print(f"[green]Dataset materializado em[/green] {target_dir} ({summary})")
    if stats["conflicts"]:
        console.print(f"[yellow]{stats['conflicts']} arquivo(s) em conflito ficaram na origem:[/yellow] {source_dir}")
    return target_dir


def sync_bcftools_chain_artifacts(
    source_dir: Path,
    target_dir: Path,
    *,
    link_mode: str = "hardlink",
    dry_run: bool = True,
    sample_limit: Optional[int] = None,
    genes: Optional[Sequence[str]] = None,
) -> Dict[str, int]:
    """Populate target dataset with legacy artifacts required by bcftools_chain.

    The canonical v1 layout stores reference windows under ``references/`` and
    omits per-sample raw FASTA/VCF intermediates. ``bcftools_chain`` still needs
    the per-sample raw FASTAs and consensus-ready VCFs. This helper links or
    copies only those missing artifacts from a legacy source dataset.
    """
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)
    if link_mode not in {"hardlink", "symlink", "copy"}:
        raise ValueError("link_mode deve ser 'hardlink', 'symlink' ou 'copy'")
    if not (source_dir / "individuals").exists():
        raise FileNotFoundError(f"Diretorio de individuos ausente no source: {source_dir}")
    if not (target_dir / "individuals").exists():
        raise FileNotFoundError(f"Diretorio de individuos ausente no target: {target_dir}")

    selected_genes = set(genes or [])
    stats = {
        "samples": 0,
        "windows": 0,
        "exists": 0,
        "missing": 0,
        "linked": 0,
        "copied": 0,
        "source_missing": 0,
    }
    suffix_templates = (
        "{sample}.H1.window.raw.fa",
        "{sample}.H2.window.raw.fa",
        "{sample}.window.consensus_ready.vcf.gz",
        "{sample}.window.consensus_ready.vcf.gz.tbi",
        "{sample}.window.vcf.gz",
        "{sample}.window.vcf.gz.tbi",
    )

    sample_dirs = list(_iter_individual_dirs(source_dir))
    if sample_limit is not None:
        sample_dirs = sample_dirs[:sample_limit]

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Sincronizando artefatos bcftools_chain...", total=len(sample_dirs))
        for source_sample_dir in sample_dirs:
            sample_id = source_sample_dir.name
            source_windows_dir = source_sample_dir / "windows"
            target_windows_dir = target_dir / "individuals" / sample_id / "windows"
            if not source_windows_dir.exists() or not target_windows_dir.exists():
                progress.update(task, advance=1)
                continue
            stats["samples"] += 1
            for source_window_dir in sorted(p for p in source_windows_dir.iterdir() if p.is_dir()):
                gene = source_window_dir.name
                if selected_genes and gene not in selected_genes:
                    continue
                target_window_dir = target_windows_dir / gene
                if not target_window_dir.exists():
                    continue
                stats["windows"] += 1
                for template in suffix_templates:
                    filename = template.format(sample=sample_id)
                    src = source_window_dir / filename
                    dst = target_window_dir / filename
                    if not src.exists():
                        stats["source_missing"] += 1
                        continue
                    result = _materialize_file(src, dst, link_mode=link_mode, dry_run=dry_run)
                    stats[result] += 1
            progress.update(task, advance=1)

    console.print(
        "[green]bcftools_chain artifact sync[/green] "
        f"dry_run={dry_run} link_mode={link_mode} stats={stats}"
    )
    return stats
