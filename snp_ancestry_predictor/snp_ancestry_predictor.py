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
import concurrent.futures
import glob as glob_module
import hashlib
import json
import math
import multiprocessing
import os
import resource
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

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
from rich.rule import Rule
from rich.table import Table

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))
from vcf_to_23andme.vcf_to_23andme import (
    normalize_chrom,
    generate_header,
    ensure_chip_panel_available,
    ensure_dbsnp_available,
    detect_vcf_chrom_style,
)


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
        cfg = yaml.safe_load(f)
    cfg["__config_dir"] = str(path.resolve().parent)
    return cfg


def load_splits(path: str) -> dict:
    """Load splits_metadata.json."""
    with open(path, "r") as f:
        return json.load(f)


def load_split_index(path: str, individuals_dir: str) -> dict:
    """Build split metadata from genotype_based_predictor split_index.json."""
    with open(path, "r") as f:
        split_index = json.load(f)
    ind_root = Path(individuals_dir)
    out = {"train": [], "val": [], "test": []}
    for split_name in ("train", "val", "test"):
        for sample_id in split_index.get(split_name, []):
            meta_path = ind_root / sample_id / "individual_metadata.json"
            if not meta_path.exists():
                continue
            with open(meta_path) as f:
                ind_meta = json.load(f)
            out[split_name].append({
                "sample_id": sample_id,
                "superpopulation": ind_meta.get("superpopulation"),
                "population": ind_meta.get("population"),
                "sex": ind_meta.get("sex", 0),
                "family_id": ind_meta.get("family_id", sample_id),
            })
    return out


def load_configured_splits(config: dict) -> dict:
    """Load either legacy splits_metadata or genotype cache split_index."""
    inp = config["input"]
    if inp.get("split_index"):
        return load_split_index(inp["split_index"], inp["individuals_dir"])
    if inp.get("genotype_cache_dir"):
        split_index = Path(inp["genotype_cache_dir"]) / "split_index.json"
        if split_index.exists():
            return load_split_index(str(split_index), inp["individuals_dir"])
    return load_splits(inp["splits_metadata"])


def _resolve_derived_target(source_value: Optional[str], dt_cfg: dict) -> Optional[str]:
    """Resolve one derived target value from its source field."""
    class_map = dt_cfg.get("class_map", {})
    exclude_unmapped = dt_cfg.get("exclude_unmapped", False)
    for class_name, source_values in class_map.items():
        if source_value in source_values:
            return class_name
    if exclude_unmapped:
        return None
    return source_value


def sample_metadata(splits: dict, config: Optional[dict] = None) -> Dict[str, dict]:
    """Build sample metadata, optionally including derived targets."""
    pred_cfg = (config or {}).get("prediction", {})
    derived_targets = pred_cfg.get("derived_targets", {})

    meta: Dict[str, dict] = {}
    for split_name in ("train", "val", "test"):
        for entry in splits.get(split_name, []):
            sample_id = entry["sample_id"]
            item = {
                "superpopulation": entry["superpopulation"],
                "population": entry["population"],
                "sex": entry.get("sex", 0),
                "split": split_name,
            }
            for level_name, dt_cfg in derived_targets.items():
                source_field = dt_cfg.get("source_field")
                source_value = entry.get(source_field) if source_field else None
                item[level_name] = _resolve_derived_target(source_value, dt_cfg)
            meta[sample_id] = item
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


def _resolve_path(path: Optional[str], config: dict) -> Optional[Path]:
    """Resolve config paths relative to the YAML file directory."""
    if not path:
        return None
    p = Path(path).expanduser()
    if p.is_absolute():
        return p
    base = Path(config.get("__config_dir", "."))
    return (base / p).resolve()


def _as_list(value, default: Optional[List[str]] = None) -> List[str]:
    if value is None:
        return list(default or [])
    if isinstance(value, str):
        return [value]
    return list(value)


def _region_source_cfg(config: dict, section: dict) -> dict:
    """Merge global input region settings with a step-specific section."""
    inp = config.get("input", {})
    merged = {
        k: inp.get(k)
        for k in (
            "region_bed",
            "genotype_dataset_dir",
            "genotype_cache_dir",
            "genotype_config",
            "genotype_view",
            "window_center_size",
            "genes_to_use",
        )
        if inp.get(k) is not None
    }
    for key, value in section.items():
        if key in (
            "region_bed",
            "genotype_dataset_dir",
            "genotype_cache_dir",
            "genotype_config",
            "genotype_view",
            "window_center_size",
            "genes_to_use",
        ):
            merged[key] = value
    return merged


def _load_genotype_dataset_input(genotype_config: Path) -> dict:
    """Load dataset_input from a genotype_based_predictor YAML/view pair."""
    with open(genotype_config) as f:
        raw = yaml.safe_load(f) or {}
    dataset_input = dict(raw.get("dataset_input", {}) or {})
    view_path = dataset_input.get("view_path")
    if view_path:
        view_file = Path(view_path).expanduser()
        if not view_file.is_absolute():
            view_file = (genotype_config.parent / view_file).resolve()
        with open(view_file) as f:
            view_payload = json.load(f)
        merged = dict(view_payload)
        for key, value in dataset_input.items():
            if key != "view_path":
                merged[key] = value
        dataset_input = merged
    return dataset_input


def _regions_from_genotype_metadata(
    dataset_dir: Path,
    genes_to_use: Optional[Iterable[str]],
    window_center_size: int,
) -> List[Tuple[str, int, int, str]]:
    """Replicate genotype_based_predictor's centered 32k window calculation."""
    meta_path = dataset_dir / "dataset_metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"dataset_metadata.json not found: {meta_path}")
    with open(meta_path) as f:
        dataset_meta = json.load(f)

    catalog = dataset_meta.get("window_catalog") or {}
    gene_order = list(genes_to_use or dataset_meta.get("genes", []) or catalog.keys())
    half = int(window_center_size) // 2
    regions: List[Tuple[str, int, int, str]] = []

    for gene in gene_order:
        gm = catalog.get(gene)
        if not gm:
            continue
        chrom = str(gm.get("chromosome", ""))
        orig_start = int(gm.get("start", 0))
        orig_size = int(gm.get("window_size", 0))
        if not chrom or orig_size <= 0:
            continue
        center = orig_start + (orig_size + 1) // 2
        start = max(0, center - half)
        end_inclusive = center + half - 1
        regions.append((chrom, start, end_inclusive + 1, str(gene)))

    if not regions:
        raise ValueError(f"No genotype windows found in {meta_path}")
    return regions


def _json_key_dict_to_int(payload: Dict) -> Dict[int, object]:
    return {int(key): value for key, value in payload.items()}


def _expanded_center_slice(axis: Dict[str, object], window_size: int) -> Tuple[int, int]:
    """Mirror DynamicIndelAligner.get_reference_centered_expanded_slice."""
    ref_length = int(axis.get("ref_length", 0))
    expanded_length = int(axis.get("expanded_length", ref_length))
    expanded_index_map = _json_key_dict_to_int(axis.get("expanded_index_map", {}))
    center_ref_idx = ref_length // 2
    center_expanded_idx = int(expanded_index_map.get(center_ref_idx, center_ref_idx))
    size = max(int(window_size), 1)
    if expanded_length <= size:
        return 0, expanded_length
    half = size // 2
    start = center_expanded_idx - half
    end = start + size
    if start < 0:
        end -= start
        start = 0
    if end > expanded_length:
        start = max(0, start - (end - expanded_length))
        end = expanded_length
    return start, end


def _axis_effective_reference_intervals(
    axis: Dict[str, object],
    genomic_start: int,
    window_size: int,
) -> List[Tuple[int, int]]:
    """Map the expanded-axis center slice back to BED-like reference intervals."""
    expanded_start, expanded_end = _expanded_center_slice(axis, window_size)
    expanded_index_map = _json_key_dict_to_int(axis.get("expanded_index_map", {}))
    selected = sorted(
        int(ref_idx)
        for ref_idx, expanded_idx in expanded_index_map.items()
        if expanded_start <= int(expanded_idx) < expanded_end
    )
    if not selected:
        return []

    intervals: List[Tuple[int, int]] = []
    run_start = selected[0]
    prev = selected[0]
    for ref_idx in selected[1:]:
        if ref_idx == prev + 1:
            prev = ref_idx
            continue
        intervals.append((genomic_start + run_start, genomic_start + prev + 1))
        run_start = prev = ref_idx
    intervals.append((genomic_start + run_start, genomic_start + prev + 1))
    return intervals


def _regions_from_genotype_cache(cache_dir: Path) -> List[Tuple[str, int, int, str]]:
    """Load effective expanded-axis windows from a genotype_based_predictor cache."""
    meta_path = cache_dir / "metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"metadata.json not found in genotype cache: {meta_path}")
    with open(meta_path) as f:
        meta = json.load(f)
    window_meta = meta.get("gene_window_metadata") or {}
    gene_order = list(meta.get("gene_order") or window_meta.keys())
    processing = meta.get("processing_params") or {}
    window_size = int(processing.get("window_center_size", 32768))
    alignment_signature = meta.get("alignment_cache_signature") or {}
    regions: List[Tuple[str, int, int, str]] = []
    for gene in gene_order:
        gm = window_meta.get(gene)
        if not gm:
            continue
        chrom = str(gm.get("chromosome", ""))
        start = int(gm.get("start", 0))
        end_value = int(gm.get("end", 0))
        if not chrom or end_value <= start:
            continue

        axis_path = None
        axis_cache_dir = (alignment_signature.get(gene) or {}).get("cache_dir")
        if axis_cache_dir:
            axis_path = Path(axis_cache_dir) / "axis.json"

        if axis_path and axis_path.exists():
            with open(axis_path) as f:
                axis_payload = json.load(f)
            axis = axis_payload.get("axis", {})
            effective_intervals = _axis_effective_reference_intervals(axis, max(0, start), window_size)
            for idx, (interval_start, interval_end) in enumerate(effective_intervals):
                suffix = f"#{idx + 1}" if len(effective_intervals) > 1 else ""
                regions.append((chrom, interval_start, interval_end, f"{gene}{suffix}"))
            continue

        regions.append((chrom, max(0, start), end_value + 1, str(gene)))
    if not regions:
        raise ValueError(f"No gene_window_metadata found in {meta_path}")
    return regions


def _write_regions_bed(regions: List[Tuple[str, int, int, str]], tag_payload: dict) -> Path:
    digest = hashlib.sha1(
        json.dumps(tag_payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:12]
    out_dir = SCRIPT_DIR / "refs" / "generated_regions"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(tag_payload.get("prefix", "genotype_windows"))
    bed_path = out_dir / f"{prefix}_{digest}.bed"
    if not bed_path.exists():
        tmp = bed_path.with_suffix(bed_path.suffix + f".tmp.{os.getpid()}")
        with open(tmp, "w") as f:
            for chrom, start, end, name in regions:
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
        os.replace(tmp, bed_path)
    return bed_path


def prepare_region_bed(config: dict, section: dict, console: Optional[Console] = None) -> Optional[Path]:
    """Return a BED with analysis windows, if configured."""
    region_cfg = _region_source_cfg(config, section)
    genotype_config = _resolve_path(region_cfg.get("genotype_config"), config)
    genotype_view = _resolve_path(region_cfg.get("genotype_view"), config)
    genotype_cache_dir = _resolve_path(region_cfg.get("genotype_cache_dir"), config)
    if genotype_cache_dir:
        regions = _regions_from_genotype_cache(genotype_cache_dir)
        bed_path = _write_regions_bed(
            regions,
            {
                "prefix": "genotype_effective_windows",
                "genotype_cache_dir": str(genotype_cache_dir),
            },
        )
        if console:
            console.print(
                f"  Genotype cache windows: {len(regions):,} regions from {genotype_cache_dir}"
            )
        return bed_path

    bed = _resolve_path(region_cfg.get("region_bed"), config)
    if bed:
        if not bed.exists():
            raise FileNotFoundError(f"Region BED not found: {bed}")
        return bed

    dataset_input: Dict = {}
    if genotype_config:
        dataset_input = _load_genotype_dataset_input(genotype_config)
    elif genotype_view:
        with open(genotype_view) as f:
            dataset_input = json.load(f)

    dataset_dir_value = region_cfg.get("genotype_dataset_dir") or dataset_input.get("dataset_dir")
    if not dataset_dir_value:
        return None

    dataset_dir = Path(str(dataset_dir_value)).expanduser()
    if not dataset_dir.is_absolute():
        base = genotype_config.parent if genotype_config else Path(config.get("__config_dir", "."))
        dataset_dir = (base / dataset_dir).resolve()

    genes_to_use = region_cfg.get("genes_to_use", dataset_input.get("genes_to_use"))
    window_center_size = int(region_cfg.get(
        "window_center_size",
        dataset_input.get("window_center_size", 32768),
    ))
    regions = _regions_from_genotype_metadata(dataset_dir, genes_to_use, window_center_size)
    bed_path = _write_regions_bed(
        regions,
        {
            "dataset_dir": str(dataset_dir),
            "genes_to_use": list(genes_to_use or []),
            "window_center_size": window_center_size,
        },
    )
    if console:
        console.print(
            f"  Genotype windows: {len(regions):,} regions from {dataset_dir} "
            f"(w={window_center_size:,})"
        )
    return bed_path


def _load_region_intervals(region_bed: Optional[Path]) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    if not region_bed:
        return intervals
    with open(region_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = normalize_chrom(parts[0])
            intervals[chrom].append((int(parts[1]), int(parts[2])))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def _position_in_regions(chrom: str, pos_1based: int, intervals: Dict[str, List[Tuple[int, int]]]) -> bool:
    if not intervals:
        return True
    pos0 = int(pos_1based) - 1
    for start, end in intervals.get(normalize_chrom(chrom), []):
        if start <= pos0 < end:
            return True
        if start > pos0:
            break
    return False


def _region_chromosomes(region_bed: Optional[Path]) -> Optional[Set[str]]:
    if not region_bed:
        return None
    return set(_load_region_intervals(region_bed).keys())


def _build_variant_include_expr(section: dict) -> str:
    """Build bcftools include expression for biallelic variant classes."""
    variant_types = [str(v).lower() for v in _as_list(section.get("variant_types"), ["snp"])]
    type_exprs = []
    if "snp" in variant_types or "snv" in variant_types:
        type_exprs.append('TYPE="snp"')
    if "indel" in variant_types or "indels" in variant_types:
        type_exprs.append('TYPE="indel"')
    if not type_exprs:
        raise ValueError(f"Unsupported variant_types: {variant_types}")
    inc_parts = ["(" + " || ".join(type_exprs) + ")", "N_ALT=1"]
    if section.get("skip_no_rsid", True):
        inc_parts.append('ID~"^rs"')
    return " && ".join(inc_parts)


def _default_genotype_encoding(section: dict) -> str:
    variant_types = [str(v).lower() for v in _as_list(section.get("variant_types"), ["snp"])]
    return "gt" if any(v not in ("snp", "snv") for v in variant_types) else "alleles"


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


def _gt_to_numeric(gt_str: str) -> str:
    """Return the biallelic genotype code (0/1) without phasing payload."""
    gt = gt_str.split(":", 1)[0]
    if gt in GT_DOSE:
        return gt
    return "--"


# ═══════════════════════════════════════════════════════════════
# Step 1 — Generate 23andMe files
# ═══════════════════════════════════════════════════════════════


def _start_bcftools_query(
    vcf_path: str,
    sample_csv: str,
    inc_expr: str,
    dbsnp_vcf: Optional[str] = None,
    region_bed: Optional[str] = None,
) -> Tuple[subprocess.Popen, Optional[subprocess.Popen]]:
    """Start a bcftools query pipeline.

    When *dbsnp_vcf* is given, pipes ``bcftools annotate -c ID`` into the
    query so that rsIDs are added on-the-fly without creating intermediate
    files on disk.
    """
    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"

    if dbsnp_vcf:
        ann_cmd = [
            "bcftools", "annotate",
            "-a", dbsnp_vcf, "-c", "ID",
        ]
        if region_bed:
            ann_cmd.extend(["-R", region_bed])
        ann_cmd.append(vcf_path)
        q_cmd = [
            "bcftools", "query",
            "-s", sample_csv,
            "-i", inc_expr,
            "-f", fmt,
        ]
        q_cmd.extend([
            "-",
        ])
        p_ann = subprocess.Popen(ann_cmd, stdout=subprocess.PIPE)
        proc = subprocess.Popen(
            q_cmd, stdin=p_ann.stdout,
            stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
        )
        p_ann.stdout.close()
        return proc, p_ann

    cmd = [
        "bcftools", "query",
        "-s", sample_csv,
        "-i", inc_expr,
        "-f", fmt,
    ]
    if region_bed:
        cmd.extend(["-R", region_bed])
    cmd.append(vcf_path)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
    )
    return proc, None


def _partial_path(
    sid: str, ind_dir: Path, base_tpl: str, chrom: str
) -> Path:
    """Path of the per-chromosome ``.partial`` checkpoint file for *sid*.

    Example: ``{sid}_v5_23andme.txt`` -> ``{sid}_v5_23andme.chr5.partial``.
    """
    base = base_tpl.format(sample_id=sid)
    stem, _ = os.path.splitext(base)
    return ind_dir / sid / f"{stem}.{chrom}.partial"


def _tmp_path(
    sid: str, ind_dir: Path, base_tpl: str, chrom: str
) -> Path:
    """Path of the in-flight ``.tmp`` file for the chrom checkpoint."""
    base = base_tpl.format(sample_id=sid)
    stem, _ = os.path.splitext(base)
    return ind_dir / sid / f"{stem}.{chrom}.tmp"


def _final_path(sid: str, ind_dir: Path, base_tpl: str) -> Path:
    return ind_dir / sid / base_tpl.format(sample_id=sid)


def _concatenate_partials(
    sid: str,
    ind_dir: Path,
    base_tpl: str,
    chroms: List[str],
    header: str,
) -> bool:
    """Assemble final 23andMe file from all per-chrom ``.partial`` files.

    Writes ``header + chr1.partial + ... + chrX.partial`` atomically to
    ``{sid}_*.txt`` via a ``.assembling`` staging file and ``os.replace``.
    Removes the ``.partial`` files only after the rename succeeds.

    Returns True iff the final file ended up in place.
    """
    final = _final_path(sid, ind_dir, base_tpl)
    asm = final.with_suffix(final.suffix + ".assembling")

    missing = [
        c for c in chroms
        if not _partial_path(sid, ind_dir, base_tpl, c).exists()
    ]
    if missing:
        return False

    with open(asm, "w") as out:
        out.write(header)
        for c in chroms:
            p = _partial_path(sid, ind_dir, base_tpl, c)
            with open(p) as src:
                shutil.copyfileobj(src, out)

    os.replace(asm, final)

    for c in chroms:
        _partial_path(sid, ind_dir, base_tpl, c).unlink(missing_ok=True)
    return True


def _cleanup_orphan_tmps(
    sids: List[str], ind_dir: Path, base_tpl: str
) -> int:
    """Remove ``.tmp`` files left by interrupted previous runs.

    Returns the number of files removed.
    """
    n = 0
    for sid in sids:
        sample_dir = ind_dir / sid
        if not sample_dir.is_dir():
            continue
        base = base_tpl.format(sample_id=sid)
        stem, _ = os.path.splitext(base)
        for tmp in sample_dir.glob(f"{stem}.*.tmp"):
            try:
                tmp.unlink()
                n += 1
            except OSError:
                pass
    return n


def _process_chromosome_batch(
    vcf_path: str,
    batch_ids: List[str],
    ind_dir: Path,
    fname_tpl: str,
    inc_expr: str,
    panel: Optional[Set[str]],
    chrom_label: str,
    dbsnp_vcf: Optional[str] = None,
    region_bed: Optional[str] = None,
    genotype_encoding: str = "alleles",
) -> Dict[str, int]:
    """Run bcftools query for *batch_ids* on one chromosome.

    Output for each sample is written first to ``{sid}_*.{chrom}.tmp`` and,
    upon successful exit of bcftools, atomically renamed to
    ``{sid}_*.{chrom}.partial`` (the per-chrom checkpoint).

    If interrupted, the ``.tmp`` files are left behind and cleaned by the
    next run; ``.partial`` files from earlier successful chromosomes stay
    intact, providing per-chromosome resume.

    Returns per-individual SNP counts for this chromosome.
    """
    n = len(batch_ids)
    sample_csv = ",".join(batch_ids)

    proc, p_ann = _start_bcftools_query(
        vcf_path, sample_csv, inc_expr, dbsnp_vcf, region_bed,
    )

    fhs: Dict[str, "IO"] = {}
    tmp_paths: Dict[str, Path] = {}
    for sid in batch_ids:
        tmp = _tmp_path(sid, ind_dir, fname_tpl, chrom_label)
        tmp.parent.mkdir(parents=True, exist_ok=True)
        tmp_paths[sid] = tmp
        fhs[sid] = open(tmp, "w")

    counts: Dict[str, int] = {sid: 0 for sid in batch_ids}

    try:
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 + n:
                continue

            var_id = parts[2]
            ref, alt = parts[3], parts[4]
            nc = normalize_chrom(parts[0])
            if not var_id or var_id == ".":
                var_id = f"{nc}:{parts[1]}:{ref}:{alt}"
            if panel is not None and var_id not in panel:
                continue

            prefix = f"{var_id}\t{nc}\t{parts[1]}\t"

            for i in range(n):
                if genotype_encoding == "gt":
                    geno = _gt_to_numeric(parts[5 + i])
                else:
                    geno = _gt_to_23andme(parts[5 + i], ref, alt)
                fhs[batch_ids[i]].write(prefix + geno + "\n")
                counts[batch_ids[i]] += 1

        proc.wait()
        if p_ann is not None:
            p_ann.wait()
    finally:
        for fh in fhs.values():
            fh.close()

    if proc.returncode != 0 or (p_ann is not None and p_ann.returncode != 0):
        for tmp in tmp_paths.values():
            tmp.unlink(missing_ok=True)
        rc = proc.returncode
        raise RuntimeError(
            f"bcftools failed for {chrom_label} (exit code {rc})"
        )

    for sid, tmp in tmp_paths.items():
        partial = _partial_path(sid, ind_dir, fname_tpl, chrom_label)
        os.replace(tmp, partial)

    return counts


# ── Step 1 multiprocessing workers (must be module-level for pickle) ──

_WORKER_PANEL: Optional[Set[str]] = None
_WORKER_INC_EXPR: str = ""
_WORKER_REGION_BED: Optional[str] = None
_WORKER_GENOTYPE_ENCODING: str = "alleles"


def _worker_init(
    panel_path: Optional[str],
    inc_expr: str,
    region_bed: Optional[str],
    genotype_encoding: str,
) -> None:
    """Load SNP panel once per worker process for parallel Step 1."""
    global _WORKER_PANEL, _WORKER_INC_EXPR, _WORKER_REGION_BED, _WORKER_GENOTYPE_ENCODING
    _WORKER_PANEL = load_snp_panel(panel_path) if panel_path else None
    _WORKER_INC_EXPR = inc_expr
    _WORKER_REGION_BED = region_bed
    _WORKER_GENOTYPE_ENCODING = genotype_encoding


def _worker_run_chrom(
    args: Tuple[str, List[str], str, str, str, Optional[str]],
) -> Tuple[str, Dict[str, int]]:
    """Run ``_process_chromosome_batch`` in a worker process.

    *args*: ``(vcf_path, batch_ids, ind_dir_str, fname_tpl, chrom_label,
    dbsnp_vcf)``
    """
    vcf, batch, ind_dir_str, fname_tpl, chrom_label, dbsnp_vcf = args
    counts = _process_chromosome_batch(
        vcf,
        batch,
        Path(ind_dir_str),
        fname_tpl,
        _WORKER_INC_EXPR,
        _WORKER_PANEL,
        chrom_label=chrom_label,
        dbsnp_vcf=dbsnp_vcf,
        region_bed=_WORKER_REGION_BED,
        genotype_encoding=_WORKER_GENOTYPE_ENCODING,
    )
    return chrom_label, counts


def step1_generate_23andme(config: dict, console: Console) -> None:
    """Generate 23andMe-format files from multi-sample 1000 Genomes VCFs."""
    inp = config["input"]
    conv = config.get("conversion", {})

    splits = load_configured_splits(config)
    meta = sample_metadata(splits, config)
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

    fmt_version = conv.get("format_version", "V5")

    panel: Optional[Set[str]] = None
    panel_path_for_workers: Optional[str] = None
    if conv.get("snp_panel"):
        panel_path_for_workers = str(conv["snp_panel"])
        panel = load_snp_panel(panel_path_for_workers)
        console.print(f"  SNP panel: {len(panel):,} SNPs")
    elif conv.get("filter_by_chip_panel", False):
        ref_dir = conv.get("ref_dir", str(SCRIPT_DIR / "refs"))
        panel_file = ensure_chip_panel_available(fmt_version, ref_dir)
        panel_path_for_workers = str(panel_file)
        panel = load_snp_panel(panel_path_for_workers)
        console.print(f"  Chip panel ({fmt_version}): {len(panel):,} SNPs")

    chroms = inp.get("chromosomes", DEFAULT_CHROMS)
    vcf_pat = inp["vcf_pattern"]
    region_bed = prepare_region_bed(config, conv, console)
    if region_bed:
        region_chroms = _region_chromosomes(region_bed)
        chroms = [c for c in chroms if normalize_chrom(c) in region_chroms]
        console.print(f"  Region filter: {region_bed} ({len(chroms)} chromosome(s))")
    genotype_encoding = conv.get("genotype_encoding", _default_genotype_encoding(conv))
    if genotype_encoding not in ("alleles", "gt"):
        raise ValueError("conversion.genotype_encoding must be 'alleles' or 'gt'")

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

    dbsnp_vcf: Optional[str] = None
    if conv.get("annotate_dbsnp", False):
        build = conv.get("genome_build", "GRCh38")
        ref_dir = conv.get("ref_dir", str(SCRIPT_DIR / "refs"))
        vcf_style = detect_vcf_chrom_style(first_vcf)
        console.print(f"  Preparing dbSNP ({build}) for rsID annotation...")
        dbsnp_path = ensure_dbsnp_available(build, Path(ref_dir), vcf_style)
        dbsnp_vcf = str(dbsnp_path)
        console.print(f"  dbSNP ready: {dbsnp_path.name}")

    header = generate_header("GRCh38", fmt_version)

    for sid in valid:
        (ind_dir / sid).mkdir(parents=True, exist_ok=True)

    n_tmp_removed = _cleanup_orphan_tmps(valid, ind_dir, fname_tpl)
    if n_tmp_removed:
        console.print(
            f"  [yellow]Cleaned {n_tmp_removed:,} orphan .tmp file(s) "
            f"from a previous run[/yellow]"
        )

    n_partials = sum(
        1
        for s in valid
        for c in chroms
        if _partial_path(s, ind_dir, fname_tpl, c).exists()
    )
    total_chr_jobs = len(valid) * len(chroms)
    if n_partials > 0:
        console.print(
            f"  [yellow]Resuming: {n_partials:,}/{total_chr_jobs:,} "
            f"chromosome jobs already done from previous run "
            f"({100.0 * n_partials / total_chr_jobs:.1f}%)[/yellow]"
        )

    inc_expr = _build_variant_include_expr(conv)

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

    # Parallel chromosome jobs: resolve worker count (default auto = all CPUs).
    raw_pc = conv.get("parallel_chroms", "auto")
    if raw_pc == "auto" or raw_pc is None:
        n_workers_parallel = os.cpu_count() or 1
    else:
        try:
            n_workers_parallel = max(1, int(raw_pc))
        except (TypeError, ValueError):
            n_workers_parallel = os.cpu_count() or 1
            console.print(
                "  [yellow]parallel_chroms invalid; using auto[/yellow]"
            )

    total_jobs = len(chroms) * n_batches
    tasks: List[Tuple[str, List[str], str, str, str, Optional[str]]] = []
    for batch in batches:
        for chrom in chroms:
            vcf = find_vcf(vcf_pat, chrom)
            if not vcf:
                continue
            samples_needing = [
                s for s in batch
                if not _partial_path(s, ind_dir, fname_tpl, chrom).exists()
            ]
            if not samples_needing:
                continue
            tasks.append(
                (
                    vcf,
                    samples_needing,
                    str(ind_dir.resolve()),
                    fname_tpl,
                    chrom,
                    dbsnp_vcf,
                )
            )

    skipped_jobs = total_jobs - len(tasks)
    use_parallel = n_workers_parallel > 1 and len(tasks) > 1

    if use_parallel:
        console.print(
            f"  Parallel chromosome jobs: {n_workers_parallel} workers, "
            f"{len(tasks)} pending partition jobs ({skipped_jobs} slots skipped)"
        )
    elif n_workers_parallel > 1 and len(tasks) <= 1:
        console.print(
            "  [dim]parallel_chroms>1 but only one pending job; "
            "running sequentially[/dim]"
        )

    per_ind_snps: Dict[str, int] = {sid: 0 for sid in valid}

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as prog:
        task_prog = prog.add_task("Step 1", total=total_jobs)

        if skipped_jobs:
            prog.advance(task_prog, skipped_jobs)

        if use_parallel:
            mp_ctx = multiprocessing.get_context(
                "spawn" if sys.platform in ("win32", "darwin") else "fork"
            )
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_workers_parallel,
                mp_context=mp_ctx,
                initializer=_worker_init,
                initargs=(
                    panel_path_for_workers,
                    inc_expr,
                    str(region_bed) if region_bed else None,
                    genotype_encoding,
                ),
            ) as ex:
                futures = {
                    ex.submit(_worker_run_chrom, t): t[4] for t in tasks
                }
                try:
                    for fut in concurrent.futures.as_completed(futures):
                        chrom_label, batch_counts = fut.result()
                        for sid, cnt in batch_counts.items():
                            per_ind_snps[sid] += cnt
                        prog.update(
                            task_prog,
                            description=f"Step 1 — {chrom_label} done",
                        )
                        prog.advance(task_prog)
                except KeyboardInterrupt:
                    try:
                        ex.shutdown(wait=False, cancel_futures=True)
                    except TypeError:
                        ex.shutdown(wait=False)
                    raise
        else:
            for t in tasks:
                chrom_label, batch_counts = _process_chromosome_batch(
                    t[0],
                    t[1],
                    Path(t[2]),
                    t[3],
                    inc_expr,
                    panel,
                    chrom_label=t[4],
                    dbsnp_vcf=t[5],
                    region_bed=str(region_bed) if region_bed else None,
                    genotype_encoding=genotype_encoding,
                )
                for sid, cnt in batch_counts.items():
                    per_ind_snps[sid] += cnt
                prog.update(
                    task_prog,
                    description=(
                        f"Step 1 — {chrom_label}"
                        if n_batches == 1
                        else f"Step 1 — {chrom_label} (batch)"
                    ),
                )
                prog.advance(task_prog)

        prog.update(task_prog, description="Step 1 — assembling final files")
        for sid in valid:
            _concatenate_partials(
                sid, ind_dir, fname_tpl, chroms, header,
            )

    for sid in valid:
        final = _final_path(sid, ind_dir, fname_tpl)
        if not final.exists():
            per_ind_snps[sid] = 0
            continue
        try:
            with open(final) as f:
                per_ind_snps[sid] = sum(
                    1 for ln in f if ln and not ln.startswith("#")
                )
        except OSError:
            pass

    counts = list(per_ind_snps.values())
    min_snps, max_snps, mean_snps = min(counts), max(counts), sum(counts) / len(counts)

    panel_info = ""
    if panel is not None:
        hit_rate = mean_snps / len(panel) * 100
        panel_info = f" (panel coverage: {hit_rate:.1f}% of {len(panel):,})"

    if min_snps == max_snps:
        console.print(f"  SNPs per individual: {min_snps:,}{panel_info}")
    else:
        console.print(
            f"  SNPs per individual: min={min_snps:,}, max={max_snps:,}, "
            f"mean={mean_snps:,.0f}{panel_info}"
        )

    zero_ids = [sid for sid, c in per_ind_snps.items() if c == 0]
    if zero_ids:
        console.print(
            f"[red]Error: {len(zero_ids)} individual(s) have 0 panel SNPs:[/red]"
        )
        for sid in zero_ids[:20]:
            console.print(f"  [red]{sid}[/red]")
        if len(zero_ids) > 20:
            console.print(f"  [red]... and {len(zero_ids) - 20} more[/red]")
        sys.exit(1)

    console.print(
        f"[green]Step 1 done — {len(valid)} files written[/green]"
    )


# ═══════════════════════════════════════════════════════════════
# Step 2 — Compute allele frequency statistics
# ═══════════════════════════════════════════════════════════════


def _build_statistics_path(config: dict) -> str:
    """Build the statistics filename from the parameters that affect it.

    The generated name encodes reference_subsets, level, min_maf, max_snps,
    snp_panel, and haplotype_mode so that different configurations produce
    distinct files.
    """
    stat_cfg = config.get("statistics", {})
    pred_cfg = config.get("prediction", {})

    output_dir = stat_cfg.get(
        "output_dir",
        stat_cfg.get("output_file", str(SCRIPT_DIR / "results")),
    )
    if output_dir.endswith(".json"):
        output_dir = os.path.dirname(output_dir)

    ref_subsets = sorted(stat_cfg.get("reference_subsets", ["train"]))
    level = pred_cfg.get("level", "superpopulation")
    min_maf = stat_cfg.get("min_maf", 0.01)
    max_snps = stat_cfg.get("max_snps", None)
    haplotype_mode = pred_cfg.get("haplotype_mode", "H1+H2")
    split_source = config.get("input", {}).get("split_index") or config.get("input", {}).get("splits_metadata")
    split_tag = "splits" + hashlib.sha1(str(split_source).encode("utf-8")).hexdigest()[:8]

    panel_tag = "nopanel"
    snp_panel = stat_cfg.get("snp_panel")
    if snp_panel:
        panel_tag = Path(snp_panel).stem

    region_tag = "allregions"
    region_cfg = _region_source_cfg(config, stat_cfg)
    region_source = (
        region_cfg.get("region_bed")
        or region_cfg.get("genotype_cache_dir")
        or region_cfg.get("genotype_dataset_dir")
        or region_cfg.get("genotype_config")
        or region_cfg.get("genotype_view")
    )
    if region_source:
        payload = {
            "region_source": str(region_source),
            "window_center_size": region_cfg.get("window_center_size"),
            "genes_to_use": region_cfg.get("genes_to_use"),
        }
        region_tag = "regions" + hashlib.sha1(
            json.dumps(payload, sort_keys=True).encode("utf-8")
        ).hexdigest()[:8]

    variant_tag = "+".join(sorted(
        str(v).lower() for v in _as_list(stat_cfg.get("variant_types"), ["snp"])
    ))
    genotype_encoding = stat_cfg.get(
        "genotype_encoding",
        config.get("conversion", {}).get("genotype_encoding", _default_genotype_encoding(stat_cfg)),
    )

    parts = [
        "snp_ancestry_statistics",
        "+".join(ref_subsets),
        level,
        split_tag,
        f"maf{min_maf}",
        f"max{max_snps}" if max_snps else "maxall",
        panel_tag,
        region_tag,
        variant_tag,
        genotype_encoding,
    ]
    if haplotype_mode != "H1+H2":
        parts.append(haplotype_mode)
    filename = "_".join(parts) + ".json"
    return os.path.join(output_dir, filename)


def step2_compute_statistics(config: dict, console: Console) -> None:
    """Compute per-population allele frequencies from 23andMe files."""
    inp = config["input"]
    stat_cfg = config.get("statistics", {})
    pred_cfg = config.get("prediction", {})

    output_file = _build_statistics_path(config)
    if os.path.exists(output_file):
        console.print(f"[green]Statistics file already exists: {output_file}[/green]")
        console.print("  Reusing. Delete the file to force recomputation.")
        return

    splits = load_configured_splits(config)
    meta = sample_metadata(splits, config)

    level = pred_cfg.get("level", "superpopulation")
    ref_subsets = stat_cfg.get("reference_subsets", ["train"])
    ref_ids = sorted(
        sid
        for sid, m in meta.items()
        if m["split"] in ref_subsets and m.get(level) is not None
    )

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

    min_maf = stat_cfg.get("min_maf", 0.01)
    max_snps = stat_cfg.get("max_snps", None)
    haplotype_mode = pred_cfg.get("haplotype_mode", "H1+H2")
    genotype_encoding = stat_cfg.get(
        "genotype_encoding",
        config.get("conversion", {}).get("genotype_encoding", _default_genotype_encoding(stat_cfg)),
    )
    if genotype_encoding not in ("alleles", "gt"):
        raise ValueError("statistics.genotype_encoding must be 'alleles' or 'gt'")
    region_bed = prepare_region_bed(config, stat_cfg, console)
    region_intervals = _load_region_intervals(region_bed)

    console.print(f"  Haplotype mode: {haplotype_mode}")
    console.print(f"  Genotype encoding: {genotype_encoding}")
    if region_bed:
        console.print(f"  Region filter: {region_bed}")

    ind_dir = Path(inp["individuals_dir"])
    fname_tpl = config.get("conversion", {}).get(
        "output_filename", "{sample_id}_23andme.txt"
    )

    tracked_allele: Dict[str, str] = {}
    snp_chrom_pos: Dict[str, Tuple[str, str]] = {}
    snp_count: Dict[str, List[int]] = {}
    snp_total: Dict[str, List[int]] = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as prog:
        task = prog.add_task("Step 2 — reading 23andMe files", total=len(ref_ids))

        for sid in ref_ids:
            prog.update(task, description=f"Step 2 — {sid}")
            pidx = ind_pop_idx[sid]

            fpath = ind_dir / sid / fname_tpl.format(sample_id=sid)
            if not fpath.exists():
                console.print(
                    f"  [yellow]23andMe file missing for {sid}, skipping[/yellow]"
                )
                prog.advance(task)
                continue

            with open(fpath) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue

                    rsid = parts[0]
                    geno = parts[3]

                    if geno == "--":
                        continue
                    if stat_panel is not None and rsid not in stat_panel:
                        continue
                    if not _position_in_regions(parts[1], int(parts[2]), region_intervals):
                        continue

                    if rsid not in tracked_allele:
                        tracked_allele[rsid] = "0" if genotype_encoding == "gt" else geno[0]
                        snp_chrom_pos[rsid] = (parts[1], parts[2])
                        snp_count[rsid] = [0] * K
                        snp_total[rsid] = [0] * K

                    ref = tracked_allele[rsid]
                    dose = _genotype_dose(geno, ref, haplotype_mode)
                    if dose is None:
                        continue
                    snp_count[rsid][pidx] += dose
                    snp_total[rsid][pidx] += 1 if haplotype_mode in ("H1", "H2") else 2

            prog.advance(task)

    console.print(f"  Total SNPs collected: {len(snp_count):,}")

    if len(snp_count) == 0:
        console.print(
            "[red]No SNPs collected from any 23andMe file.[/red]\n"
            "  Common causes:\n"
            "    • conversion.output_filename points at files generated with a\n"
            "      different SNP panel (no rsIDs overlap with statistics.snp_panel).\n"
            "    • statistics.snp_panel is empty or unreadable.\n"
            "    • The 23andMe files exist but are empty (e.g. failed Step 1).\n"
            "  Aborting before saving an empty statistics file."
        )
        sys.exit(1)

    snp_ref_allele = tracked_allele
    snp_ref_count = snp_count
    snp_total_count = snp_total

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
            "region_bed": str(region_bed) if region_bed else None,
            "variant_types": _as_list(stat_cfg.get("variant_types"), ["snp"]),
            "genotype_encoding": genotype_encoding,
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


def _genotype_dose(
    genotype: str, ref_allele: str, haplotype_mode: str = "H1+H2",
) -> Optional[int]:
    """Count of *ref_allele* copies in a genotype string.

    *haplotype_mode* selects which alleles to consider:
      "H1"    — only genotype[0]  (dose 0 or 1)
      "H2"    — only genotype[1]  (dose 0 or 1)
      "H1+H2" — both alleles      (dose 0, 1, or 2)
    """
    if genotype == "--":
        return None

    if "/" in genotype or "|" in genotype:
        sep = "|" if "|" in genotype else "/"
        alleles = genotype.split(sep)
        if len(alleles) != 2 or "." in alleles:
            return None
        if haplotype_mode == "H1":
            return 1 if alleles[0] == ref_allele else 0
        if haplotype_mode == "H2":
            return 1 if alleles[1] == ref_allele else 0
        return sum(1 for a in alleles if a == ref_allele)

    if len(genotype) != 2:
        return None
    if haplotype_mode == "H1":
        return 1 if genotype[0] == ref_allele else 0
    if haplotype_mode == "H2":
        return 1 if genotype[1] == ref_allele else 0
    return sum(1 for c in genotype if c == ref_allele)


def classify_mle(
    genotypes: Dict[str, str],
    ref_alleles: Dict[str, str],
    allele_freqs: Dict[str, List[float]],
    populations: List[str],
    haplotype_mode: str = "H1+H2",
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
        d = _genotype_dose(gt, ref_alleles.get(rsid, ""), haplotype_mode)
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

    if haplotype_mode in ("H1", "H2"):
        log_prob = np.where(
            dose_arr[:, None] == 1, log_p, log_1p,
        )
    else:
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


def _admixture_nll_softmax(
    theta: np.ndarray,
    dose_arr: np.ndarray,
    freq_matrix: np.ndarray,
    freq_matrix_T: np.ndarray,
    ploidy: int = 2,
) -> Tuple[float, np.ndarray]:
    """NLL and gradient for the admixture model with softmax reparametrisation.

    Uses unconstrained parameters *theta* mapped to the simplex via softmax,
    allowing the use of L-BFGS-B (no equality constraints needed).
    *ploidy* is 1 for haploid (single-haplotype) or 2 for diploid.
    Returns (nll, grad_theta) so ``jac=True`` can be used.
    """
    eps = 1e-10

    theta_s = theta - theta.max()
    exp_t = np.exp(theta_s)
    alpha = exp_t / exp_t.sum()

    p_mixed = np.clip(freq_matrix @ alpha, eps, 1.0 - eps)

    if ploidy == 1:
        log_prob = dose_arr * np.log(p_mixed) + (1.0 - dose_arr) * np.log(1.0 - p_mixed)
        dlog_dp = dose_arr / p_mixed - (1.0 - dose_arr) / (1.0 - p_mixed)
    else:
        log_prob = np.where(
            dose_arr == 2,
            2.0 * np.log(p_mixed),
            np.where(
                dose_arr == 1,
                np.log(2.0) + np.log(p_mixed) + np.log(1.0 - p_mixed),
                2.0 * np.log(1.0 - p_mixed),
            ),
        )
        dlog_dp = np.where(
            dose_arr == 2,
            2.0 / p_mixed,
            np.where(
                dose_arr == 1,
                1.0 / p_mixed - 1.0 / (1.0 - p_mixed),
                -2.0 / (1.0 - p_mixed),
            ),
        )

    nll = -float(np.sum(log_prob))
    grad_alpha = -(freq_matrix_T @ dlog_dp)
    grad_theta = alpha * (grad_alpha - np.dot(grad_alpha, alpha))

    return nll, grad_theta


def estimate_admixture(
    genotypes: Dict[str, str],
    ref_alleles: Dict[str, str],
    allele_freqs: Dict[str, List[float]],
    populations: List[str],
    n_restarts: int = 20,
    haplotype_mode: str = "H1+H2",
) -> Tuple[Dict[str, float], int]:
    """Admixture proportion estimation via L-BFGS-B with softmax reparametrisation.

    Parameters are mapped to the probability simplex via softmax, eliminating
    the need for equality constraints and enabling the much faster L-BFGS-B
    optimiser with analytical gradient.

    Returns ({pop: proportion}, n_snps_used).
    """
    K = len(populations)
    ploidy = 1 if haplotype_mode in ("H1", "H2") else 2

    doses_list: List[int] = []
    freq_list: List[List[float]] = []

    for rsid, freqs in allele_freqs.items():
        gt = genotypes.get(rsid)
        if gt is None:
            continue
        d = _genotype_dose(gt, ref_alleles.get(rsid, ""), haplotype_mode)
        if d is None:
            continue
        doses_list.append(d)
        freq_list.append(freqs)

    if not doses_list:
        return {p: 1.0 / K for p in populations}, 0

    dose_arr = np.array(doses_list, dtype=np.int32)
    freq_matrix = np.array(freq_list, dtype=np.float64)
    freq_matrix_T = np.ascontiguousarray(freq_matrix.T)

    best_nll = float("inf")
    best_alpha = np.ones(K) / K
    rng = np.random.RandomState(42)
    no_improve = 0

    for i in range(n_restarts):
        if i == 0:
            theta0 = np.zeros(K)
        else:
            d = rng.dirichlet(np.ones(K))
            theta0 = np.log(np.clip(d, 1e-10, None))

        result = minimize(
            _admixture_nll_softmax,
            theta0,
            args=(dose_arr, freq_matrix, freq_matrix_T, ploidy),
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": 200, "ftol": 1e-10},
        )

        if result.fun < best_nll - 1e-6:
            best_nll = result.fun
            t = result.x
            t_s = t - t.max()
            e = np.exp(t_s)
            best_alpha = e / e.sum()
            no_improve = 0
        else:
            no_improve += 1
            if no_improve >= 3:
                break

    return (
        {populations[k]: float(best_alpha[k]) for k in range(K)},
        len(doses_list),
    )


def estimate_admixture_em(
    genotypes: Dict[str, str],
    ref_alleles: Dict[str, str],
    allele_freqs: Dict[str, List[float]],
    populations: List[str],
    max_iter: int = 1000,
    tol: float = 1e-7,
    haplotype_mode: str = "H1+H2",
) -> Tuple[Dict[str, float], int]:
    """Admixture proportion estimation via the EM algorithm.

    Uses the classical Expectation-Maximization algorithm from FRAPPE
    (Tang et al., 2005) to estimate admixture proportions.  Each allele
    copy is treated as belonging latently to one of K populations; the
    E-step computes posterior responsibilities and the M-step updates
    the proportions.

    This is much faster than numerical optimisers (L-BFGS-B, SLSQP)
    because each iteration is a few pure-numpy matrix–vector multiplies
    with no scipy overhead, and convergence is monotonic (the likelihood
    never decreases).

    Returns ({pop: proportion}, n_snps_used).
    """
    K = len(populations)
    ploidy = 1 if haplotype_mode in ("H1", "H2") else 2

    doses_list: List[int] = []
    freq_list: List[List[float]] = []

    for rsid, freqs in allele_freqs.items():
        gt = genotypes.get(rsid)
        if gt is None:
            continue
        d = _genotype_dose(gt, ref_alleles.get(rsid, ""), haplotype_mode)
        if d is None:
            continue
        doses_list.append(d)
        freq_list.append(freqs)

    if not doses_list:
        return {p: 1.0 / K for p in populations}, 0

    dose_arr = np.array(doses_list, dtype=np.float64)
    freq_matrix = np.array(freq_list, dtype=np.float64)
    freq_matrix_T = np.ascontiguousarray(freq_matrix.T)

    eps = 1e-10
    dose_alt = float(ploidy) - dose_arr

    alpha = np.ones(K, dtype=np.float64) / K

    for _ in range(max_iter):
        p_mixed = np.clip(freq_matrix @ alpha, eps, 1.0 - eps)
        q_mixed = 1.0 - p_mixed

        ratio_ref = dose_arr / p_mixed
        ratio_alt = dose_alt / q_mixed

        c_ref = freq_matrix_T @ ratio_ref
        c_alt = ratio_alt.sum() - freq_matrix_T @ ratio_alt

        alpha_new = alpha * (c_ref + c_alt)
        alpha_new /= alpha_new.sum()

        if np.max(np.abs(alpha_new - alpha)) < tol:
            alpha = alpha_new
            break
        alpha = alpha_new

    return (
        {populations[k]: float(alpha[k]) for k in range(K)},
        len(doses_list),
    )


def predict_single_individual(
    config: dict, individual_path: str, console: Console,
) -> None:
    """Predict ancestry for a single individual from a 23andMe file.

    Loads pre-computed statistics from the YAML-configured path and runs
    the configured method on the provided file.  No splits_metadata or
    individuals_dir is required.
    """
    stat_cfg = config.get("statistics", {})
    pred_cfg = config.get("prediction", {})

    output_file = _build_statistics_path(config)
    if not os.path.exists(output_file):
        console.print(
            f"[red]Statistics file not found: {output_file}[/red]\n"
            "  Run the full pipeline (Steps 1-2) first to compute statistics."
        )
        return

    ind_path = Path(individual_path)
    if not ind_path.exists():
        console.print(f"[red]File not found: {individual_path}[/red]")
        return

    with open(output_file) as f:
        stats = json.load(f)

    populations = stats["metadata"]["populations"]
    ref_alleles: Dict[str, str] = stats["ref_alleles"]
    allele_freqs: Dict[str, List[float]] = stats["allele_frequencies"]
    needed_snps = set(allele_freqs.keys())

    method = pred_cfg.get("method", "admixture_em")
    level = pred_cfg.get("level", "superpopulation")

    sample_id = ind_path.stem
    console.print(
        f"  Statistics: {len(needed_snps):,} SNPs, "
        f"{len(populations)} {stats['metadata']['level']} classes"
    )
    haplotype_mode = pred_cfg.get("haplotype_mode", "H1+H2")
    console.print(f"  File: {individual_path}")
    console.print(f"  Method: {method}")
    console.print(f"  Haplotype mode: {haplotype_mode}")

    genos = _load_23andme_genotypes(str(ind_path), needed_snps)
    if not genos:
        console.print("[red]No matching SNPs found in the file.[/red]")
        return

    if method == "mle":
        pred, log_liks, n_used = classify_mle(
            genos, ref_alleles, allele_freqs, populations, haplotype_mode,
        )
        ll_arr = np.array([log_liks[p] for p in populations])
        ll_arr -= ll_arr.max()
        exp_ll = np.exp(ll_arr)
        probs = exp_ll / exp_ll.sum()
        props = {
            populations[k]: round(float(probs[k]), 6)
            for k in range(len(populations))
        }
    elif method == "admixture_em":
        em_cfg = pred_cfg.get("admixture_em", {})
        props, n_used = estimate_admixture_em(
            genos, ref_alleles, allele_freqs, populations,
            max_iter=em_cfg.get("max_iter", 1000),
            tol=em_cfg.get("tol", 1e-7),
            haplotype_mode=haplotype_mode,
        )
        pred = max(props, key=props.get)
    elif method == "admixture_mle":
        n_restarts = pred_cfg.get("admixture_mle", {}).get("n_restarts", 20)
        props, n_used = estimate_admixture(
            genos, ref_alleles, allele_freqs, populations, n_restarts,
            haplotype_mode=haplotype_mode,
        )
        pred = max(props, key=props.get)
    else:
        console.print(f"[red]Unknown method: {method}[/red]")
        return

    entry = {
        "sample_id": sample_id,
        "predicted": pred,
        "snps_used": n_used,
        "proportions": props,
        "method": method,
        "level": level,
    }

    ind_file = Path(f"{sample_id}_{method}_{level}.json")
    with open(ind_file, "w") as fout:
        json.dump(entry, fout, indent=2)

    console.print(Rule(
        f"Prediction — {sample_id}", style="bold cyan",
    ))

    pt = Table(title="Ancestry Proportions")
    pt.add_column("Population")
    pt.add_column("Proportion", justify="right")
    sorted_props = sorted(props.items(), key=lambda x: x[1], reverse=True)
    for pop, frac in sorted_props:
        pt.add_row(pop, f"{frac:.6f}")
    console.print(pt)

    console.print(f"  Predicted ancestry: [bold]{pred}[/bold]")
    console.print(f"  SNPs used: {n_used:,}")
    console.print(f"\n[green]Result saved to {ind_file}[/green]")
    console.print(
        f"  [dim]Generate a pie chart with:[/dim]  "
        f"python3 plot_ancestry_pie.py {ind_file}"
    )


def step3_predict_ancestry(config: dict, console: Console) -> None:
    """Predict ancestry for evaluation individuals."""
    inp = config["input"]
    stat_cfg = config.get("statistics", {})
    pred_cfg = config.get("prediction", {})

    output_file = _build_statistics_path(config)
    if not os.path.exists(output_file):
        console.print(
            f"[red]Statistics file not found: {output_file}[/red]\n"
            "  Run Step 2 first."
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

    splits = load_configured_splits(config)
    meta = sample_metadata(splits, config)
    level = pred_cfg.get("level", "superpopulation")
    eval_subsets = pred_cfg.get("evaluation_subsets", ["test"])
    eval_ids = [
        sid for sid, m in meta.items()
        if m["split"] in eval_subsets and m.get(level) is not None
    ]

    ind_dir = Path(inp["individuals_dir"])
    fname_tpl = config.get("conversion", {}).get(
        "output_filename", "{sample_id}_23andme.txt"
    )

    method = pred_cfg.get("method", "mle")
    results_dir = Path(
        pred_cfg.get("results_dir", str(SCRIPT_DIR / "results"))
    )
    results_dir.mkdir(parents=True, exist_ok=True)
    ind_results_dir = results_dir / "individuals"
    ind_results_dir.mkdir(parents=True, exist_ok=True)

    haplotype_mode = pred_cfg.get("haplotype_mode", "H1+H2")
    console.print(f"  Evaluation: {len(eval_ids)} individuals from {eval_subsets}")
    console.print(f"  Method: {method}")
    console.print(f"  Haplotype mode: {haplotype_mode}")

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
            sid_split = meta[sid]["split"]

            if method == "mle":
                pred, log_liks, n_used = classify_mle(
                    genos, ref_alleles, allele_freqs, populations,
                    haplotype_mode,
                )
                ll_arr = np.array([log_liks[p] for p in populations])
                ll_arr -= ll_arr.max()
                exp_ll = np.exp(ll_arr)
                probs = exp_ll / exp_ll.sum()
                props = {
                    populations[k]: round(float(probs[k]), 6)
                    for k in range(len(populations))
                }
                entry = {
                    "sample_id": sid,
                    "split": sid_split,
                    "true": true_label,
                    "predicted": pred,
                    "correct": pred == true_label,
                    "snps_used": n_used,
                    "proportions": props,
                }
                predictions.append(entry)

            elif method == "admixture_mle":
                n_restarts = (
                    pred_cfg.get("admixture_mle", {}).get("n_restarts", 20)
                )
                props, n_used = estimate_admixture(
                    genos, ref_alleles, allele_freqs, populations, n_restarts,
                    haplotype_mode=haplotype_mode,
                )
                pred = max(props, key=props.get)
                entry = {
                    "sample_id": sid,
                    "split": sid_split,
                    "true": true_label,
                    "predicted": pred,
                    "correct": pred == true_label,
                    "snps_used": n_used,
                    "proportions": props,
                }
                predictions.append(entry)

            elif method == "admixture_em":
                em_cfg = pred_cfg.get("admixture_em", {})
                em_max_iter = em_cfg.get("max_iter", 1000)
                em_tol = em_cfg.get("tol", 1e-7)
                props, n_used = estimate_admixture_em(
                    genos, ref_alleles, allele_freqs, populations,
                    max_iter=em_max_iter, tol=em_tol,
                    haplotype_mode=haplotype_mode,
                )
                pred = max(props, key=props.get)
                entry = {
                    "sample_id": sid,
                    "split": sid_split,
                    "true": true_label,
                    "predicted": pred,
                    "correct": pred == true_label,
                    "snps_used": n_used,
                    "proportions": props,
                }
                predictions.append(entry)
            else:
                console.print(f"[red]Unknown method: {method}[/red]")
                return

            ind_file = ind_results_dir / f"{sid}_{method}_{level}.json"
            ind_data = {**entry, "method": method, "level": level}
            with open(ind_file, "w") as fout:
                json.dump(ind_data, fout, indent=2)

            prog.advance(task)

    if not predictions:
        console.print("[yellow]No predictions produced.[/yellow]")
        return

    def _compute_metrics(
        preds: List[dict], populations: List[str],
    ) -> Tuple[float, float, float, float, Dict[str, dict], Dict[str, Dict[str, int]]]:
        true_labels = [p["true"] for p in preds]
        pred_labels = [p["predicted"] for p in preds]
        total = len(preds)
        correct = sum(1 for p in preds if p["correct"])
        accuracy = correct / total

        class_met: Dict[str, dict] = {}
        for pop in populations:
            tp = sum(1 for t, p in zip(true_labels, pred_labels) if t == pop and p == pop)
            fp = sum(1 for t, p in zip(true_labels, pred_labels) if t != pop and p == pop)
            fn = sum(1 for t, p in zip(true_labels, pred_labels) if t == pop and p != pop)
            prec = tp / (tp + fp) if (tp + fp) else 0.0
            rec = tp / (tp + fn) if (tp + fn) else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            support = sum(1 for t in true_labels if t == pop)
            class_met[pop] = {
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
            class_met[p]["precision"] * class_met[p]["support"]
            for p in populations
        ) / total
        w_rec = sum(
            class_met[p]["recall"] * class_met[p]["support"]
            for p in populations
        ) / total
        w_f1 = sum(
            class_met[p]["f1"] * class_met[p]["support"]
            for p in populations
        ) / total

        return accuracy, w_prec, w_rec, w_f1, class_met, confusion

    def _display_metrics(
        label: str,
        preds: List[dict],
        populations: List[str],
        console: Console,
    ) -> Tuple[float, float, float, float, Dict[str, dict], Dict[str, Dict[str, int]]]:
        accuracy, w_prec, w_rec, w_f1, class_met, confusion = _compute_metrics(
            preds, populations
        )
        total = len(preds)

        console.print()
        console.print(Rule(label, style="bold cyan"))

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
            m = class_met[pop]
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

        return accuracy, w_prec, w_rec, w_f1, class_met, confusion

    # ── Per-subset results ───────────────────────────
    subset_results: Dict[str, dict] = {}
    for subset in eval_subsets:
        subset_preds = [p for p in predictions if p["split"] == subset]
        if not subset_preds:
            console.print(f"\n[yellow]No predictions for subset '{subset}'.[/yellow]")
            continue
        acc, wp, wr, wf, cm, conf = _display_metrics(
            f"Results — {method}, {level} — subset: {subset} — {len(subset_preds)} samples",
            subset_preds,
            populations,
            console,
        )
        subset_results[subset] = {
            "accuracy": round(acc, 4),
            "weighted_precision": round(wp, 4),
            "weighted_recall": round(wr, 4),
            "weighted_f1": round(wf, 4),
            "per_class_metrics": cm,
            "confusion_matrix": conf,
            "n_individuals": len(subset_preds),
        }

    # ── Aggregate results (when more than one subset) ─
    if len(eval_subsets) > 1:
        acc, wp, wr, wf, class_metrics, confusion = _display_metrics(
            f"Results — {method}, {level} — ALL ({', '.join(eval_subsets)}) — {len(predictions)} samples",
            predictions,
            populations,
            console,
        )
    else:
        acc, wp, wr, wf, class_metrics, confusion = _compute_metrics(
            predictions, populations
        )

    # ── Save ─────────────────────────────────────────
    results_path = results_dir / f"predictions_{method}_{level}.json"

    results_output = {
        "metadata": {
            "method": method,
            "level": level,
            "evaluation_subsets": eval_subsets,
            "n_snps": len(needed_snps),
            "n_individuals": len(predictions),
            "created_at": datetime.now().isoformat(),
        },
        "metrics": {
            "accuracy": round(acc, 4),
            "weighted_precision": round(wp, 4),
            "weighted_recall": round(wr, 4),
            "weighted_f1": round(wf, 4),
        },
        "per_class_metrics": class_metrics,
        "confusion_matrix": confusion,
        "per_subset_metrics": subset_results,
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
    parser.add_argument(
        "--individual",
        type=str,
        default=None,
        help="Path to a 23andMe file for single-individual prediction. "
             "Skips Steps 1-2 and uses pre-computed statistics.",
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

    if args.individual:
        console.print(
            "\n[bold blue]══════ Single-Individual Prediction ══════[/bold blue]"
        )
        predict_single_individual(config, args.individual, console)
        return

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
