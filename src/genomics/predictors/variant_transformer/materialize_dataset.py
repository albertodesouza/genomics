from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import torch
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

from genomics.core.config_io import load_json, write_json
from genomics.core.dataset_metadata import (
    iter_sample_metadata,
    load_dataset_metadata,
    load_vcf_sources_from_dataset as load_common_vcf_sources_from_dataset,
    load_window_catalog,
)
from genomics.core.splitting import SplitSpec, records_from_objects, split_sample_records
from genomics.core.targets import target_value
from .allele_codec import classify_variant, encode_allele, normalize_length
from .constants import HAPLOTYPE_TO_ID, VARIANT_TYPE_TO_ID
from .variant_schema import Region, SampleMetadata, VariantToken, sort_tokens

console = Console()
_WARNED_VCF_FALLBACK_ROOTS: set[str] = set()

DEFAULT_VCF_ROOT_CANDIDATES: List[Path] = []


def _load_json(path: Path) -> Dict:
    return load_json(path)


def _write_json(path: Path, payload: Dict) -> None:
    write_json(path, payload)


def load_regions(path: Path) -> List[Region]:
    regions: List[Region] = []
    with open(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED invalido em {path}:{line_no}: esperado chrom start end gene_id")
            chrom = fields[0]
            start0 = int(fields[1])
            end0 = int(fields[2])
            gene_id = fields[3]
            # BED is 0-based half-open; bcftools regions are 1-based inclusive.
            regions.append(Region(chrom=chrom, start=start0 + 1, end=end0, gene_id=gene_id))
    if not regions:
        raise ValueError(f"Nenhuma regiao encontrada em {path}")
    return regions


def load_regions_from_dataset(dataset_dir: Path, genes: Optional[set[str]] = None) -> List[Region]:
    catalog = load_window_catalog(dataset_dir, load_dataset_metadata(dataset_dir))
    regions: List[Region] = []
    for gene_id, item in sorted(catalog.items()):
        if genes and gene_id not in genes:
            continue
        chrom = item.get("chromosome") or item.get("chrom")
        start = item.get("start") or item.get("start_1based") or item.get("window_start")
        end = item.get("end") or item.get("end_1based") or item.get("window_end")
        if chrom is None or start is None or end is None:
            continue
        regions.append(Region(chrom=str(chrom), start=int(start), end=int(end), gene_id=str(gene_id)))
    if not regions:
        raise ValueError(f"Nenhuma regiao valida encontrada em {dataset_dir}")
    return regions


def load_sample_metadata(path: Optional[Path], dataset_dir: Optional[Path], target: str, classes: Sequence[str]) -> List[SampleMetadata]:
    class_set = set(classes)
    samples: List[SampleMetadata] = []
    if path is not None:
        if path.suffix.lower() == ".json":
            payload = _load_json(path)
            rows = payload.get("samples", payload if isinstance(payload, list) else [])
            for row in rows:
                sid = str(row.get("sample_id") or row.get("id"))
                value = str(row.get(target, row.get("target", "")))
                if value in class_set:
                    samples.append(SampleMetadata(sid, value, row.get("family_id"), row.get("population"), row.get("superpopulation"), row))
        else:
            with open(path) as f:
                header = f.readline().rstrip("\n").split("\t")
                col = {name: idx for idx, name in enumerate(header)}
                sid_idx = col.get("sample_id", col.get("id", 0))
                target_idx = col.get(target, col.get("target"))
                if target_idx is None:
                    raise ValueError(f"Coluna target '{target}' ausente em {path}")
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    sid = fields[sid_idx]
                    value = fields[target_idx]
                    if value not in class_set:
                        continue
                    family_id = fields[col["family_id"]] if "family_id" in col and col["family_id"] < len(fields) else None
                    pop = fields[col["population"]] if "population" in col and col["population"] < len(fields) else None
                    superpop = fields[col["superpopulation"]] if "superpopulation" in col and col["superpopulation"] < len(fields) else None
                    samples.append(SampleMetadata(sid, value, family_id, pop, superpop))
    elif dataset_dir is not None:
        for row in iter_sample_metadata(dataset_dir):
            sid = str(row.get("sample_id"))
            value = target_value(row, target)
            if value is None:
                value = row.get("target")
            if value in class_set:
                samples.append(SampleMetadata(str(sid), str(value), row.get("family_id"), row.get("population"), row.get("superpopulation"), row))
    else:
        raise ValueError("Informe --samples-metadata ou --dataset-dir")
    if not samples:
        raise ValueError("Nenhuma amostra com target valido encontrada")
    return samples


def _format_vcf_pattern(vcf_pattern: str, chromosome: str) -> str:
    try:
        return vcf_pattern.format(chrom=chromosome)
    except Exception:
        return vcf_pattern.replace("{chrom}", chromosome)


def _candidate_chromosomes(chromosome: str) -> Iterable[str]:
    chromosome = str(chromosome)
    yield chromosome
    if chromosome.startswith("chr"):
        yield chromosome[3:]
    else:
        yield f"chr{chromosome}"


def _strip_chr(chromosome: str) -> str:
    text = str(chromosome)
    return text[3:] if text.startswith("chr") else text


def resolve_vcf_path(chromosome: str, vcf_pattern: Optional[str], vcf_root_dir: Optional[str]) -> str:
    vcf_pattern = vcf_pattern or os.environ.get("KG1000_VCF_PATTERN")
    if vcf_pattern:
        candidates = [_format_vcf_pattern(vcf_pattern, chrom) for chrom in _candidate_chromosomes(chromosome)]
        for candidate in candidates:
            if Path(candidate).exists():
                return candidate
        return candidates[0]
    vcf_root_dir = vcf_root_dir or os.environ.get("KG1000_VCF_ROOT_DIR")
    roots: List[Path] = []
    if vcf_root_dir:
        roots.append(Path(vcf_root_dir))
    roots.extend(root for root in DEFAULT_VCF_ROOT_CANDIDATES if root not in roots)

    attempted: List[str] = []
    for root in roots:
        if not root.exists():
            attempted.append(str(root))
            continue
        match = _find_vcf_in_root(root, chromosome)
        if match is not None:
            if vcf_root_dir and root != Path(vcf_root_dir):
                warn_key = f"{vcf_root_dir}->{root}"
                if warn_key not in _WARNED_VCF_FALLBACK_ROOTS:
                    console.print(f"[yellow]VCF root informado nao existe/nao contem os cromossomos requisitados; usando {root}[/yellow]")
                    _WARNED_VCF_FALLBACK_ROOTS.add(warn_key)
            return str(match)
        attempted.append(str(root))
    raise FileNotFoundError(f"VCF nao encontrado para cromossomo {chromosome}. Roots testados: {attempted}")


def resolve_region_vcf_path(region: Region, vcf_pattern: Optional[str], vcf_root_dir: Optional[str], dataset_vcf_sources: Optional[Dict[str, str]] = None) -> str:
    if dataset_vcf_sources:
        exact = dataset_vcf_sources.get(region.gene_id)
        if exact:
            return exact
        chrom = _strip_chr(region.chrom)
        chrom_match = dataset_vcf_sources.get(chrom) or dataset_vcf_sources.get(f"chr{chrom}")
        if chrom_match:
            return chrom_match
    return resolve_vcf_path(region.chrom, vcf_pattern, vcf_root_dir)


def load_vcf_sources_from_dataset(dataset_dir: Optional[Path]) -> Dict[str, str]:
    return load_common_vcf_sources_from_dataset(dataset_dir)


def _find_vcf_in_root(root: Path, chromosome: str) -> Optional[Path]:
    for chrom in _candidate_chromosomes(chromosome):
        for candidate in (
            root / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
            root / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz",
            root / f"{chrom}.vcf.gz",
        ):
            if candidate.exists():
                return candidate

    candidates: List[Path] = []
    for chrom in _candidate_chromosomes(chromosome):
        candidates.extend(root.glob(f"*{chrom}*.vcf.gz"))
    if not candidates:
        return None

    def score(path: Path) -> tuple[int, str]:
        name = path.name
        priority = 0
        if "phased_panel" in name:
            priority -= 20
        if "filtered" in name:
            priority -= 10
        if name.startswith("1kGP_high_coverage_Illumina"):
            priority -= 5
        return priority, name

    return sorted(candidates, key=score)[0]


def _stream_region_variants_bcftools(vcf_path: str, sample_ids: Sequence[str], region: Region) -> Iterator[Dict]:
    fmt = "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n"
    cmd = ["bcftools", "query", "-s", ",".join(sample_ids), "-r", region.bcftools_region, "-f", fmt, vcf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            yield {"chrom": fields[0], "pos": int(fields[1]), "ref": fields[2], "alt": fields[3], "gts": fields[4:]}
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
        returncode = proc.wait()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        if proc.stderr is not None:
            proc.stderr.close()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, cmd, stderr=stderr)


def _stream_region_variants_gzip(vcf_path: str, sample_ids: Sequence[str], region: Region) -> Iterator[Dict]:
    sample_to_col: Dict[str, int] = {}
    sample_cols: List[int] = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_to_col = {sid: idx for idx, sid in enumerate(header[9:])}
                missing = [sid for sid in sample_ids if sid not in sample_to_col]
                if missing:
                    raise KeyError(f"Samples ausentes no VCF: {missing[:10]}")
                sample_cols = [9 + sample_to_col[sid] for sid in sample_ids]
                break
        if not sample_cols:
            raise RuntimeError(f"Header #CHROM nao encontrado em {vcf_path}")
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if _strip_chr(fields[0]) != _strip_chr(region.chrom):
                continue
            pos = int(fields[1])
            if pos < region.start:
                continue
            if pos > region.end:
                break
            fmt_keys = fields[8].split(":") if len(fields) > 8 else []
            gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else 0
            gts = []
            for col_idx in sample_cols:
                parts = fields[col_idx].split(":") if col_idx < len(fields) else ["./."]
                gts.append(parts[gt_idx] if gt_idx < len(parts) else parts[0])
            yield {"chrom": fields[0], "pos": pos, "ref": fields[3], "alt": fields[4], "gts": gts}


def stream_region_variants(vcf_path: str, sample_ids: Sequence[str], region: Region) -> Iterator[Dict]:
    if shutil.which("bcftools") is None:
        console.print("[yellow]bcftools nao encontrado no PATH; usando fallback gzip lento[/yellow]")
        yield from _stream_region_variants_gzip(vcf_path, sample_ids, region)
        return
    try:
        yield from _stream_region_variants_bcftools(vcf_path, sample_ids, region)
    except subprocess.CalledProcessError as exc:
        stderr = str(getattr(exc, "stderr", "") or "")
        if "BGZF EOF" not in stderr and "may be truncated" not in stderr:
            raise
        console.print(f"[yellow]bcftools falhou por BGZF EOF; fallback gzip para {Path(vcf_path).name}[/yellow]")
        yield from _stream_region_variants_gzip(vcf_path, sample_ids, region)


def _normalize_alt_alleles(alt_raw: str) -> List[str]:
    return [item for item in str(alt_raw).split(",") if item and item != "." and not item.startswith("<")]


def _parse_gt(gt: str, unphased_policy: str) -> Optional[Tuple[str, str]]:
    text = str(gt).split(":", 1)[0]
    if "|" in text:
        left, right = text.split("|", 1)
        return left, right
    if "/" in text:
        if unphased_policy == "error":
            raise ValueError(f"Genotipo nao faseado: {gt}")
        return None
    return None


def _allele_for_token(ref: str, alt_raw: str, token: str) -> Optional[str]:
    if token in {"0", ".", ""}:
        return None
    try:
        alt_idx = int(token) - 1
    except ValueError:
        return None
    alts = _normalize_alt_alleles(alt_raw)
    if 0 <= alt_idx < len(alts):
        return alts[alt_idx].upper()
    return None


def _split_sample_groups(samples: List[SampleMetadata], train: float, val: float, test: float, seed: Optional[int], family_mode: str) -> Tuple[Dict[str, List[str]], Dict[str, object]]:
    spec = SplitSpec(
        train_split=train,
        val_split=val,
        test_split=test,
        random_seed=seed,
        family_split_mode=family_mode,
    )
    return split_sample_records(records_from_objects(samples), spec)


def _tokens_to_tensors(tokens: List[VariantToken], gene_to_idx: Dict[str, int], l_max: int, max_indel_size: int) -> Dict[str, torch.Tensor]:
    length_norm = torch.tensor([[normalize_length(t.length, max_indel_size)] for t in tokens], dtype=torch.float32)
    if not tokens:
        length_norm = torch.empty((0, 1), dtype=torch.float32)
    return {
        "variant_type": torch.tensor([VARIANT_TYPE_TO_ID[t.variant_type] for t in tokens], dtype=torch.long),
        "haplotype": torch.tensor([HAPLOTYPE_TO_ID[t.haplotype] for t in tokens], dtype=torch.long),
        "gene": torch.tensor([gene_to_idx[t.gene_id] for t in tokens], dtype=torch.long),
        "length_norm": length_norm,
        "position_relative": torch.tensor([t.position_relative for t in tokens], dtype=torch.long),
        "position": torch.tensor([t.position for t in tokens], dtype=torch.long),
        "ref_allele": torch.stack([encode_allele(t.reference_allele, l_max) for t in tokens]) if tokens else torch.empty((0, l_max), dtype=torch.long),
        "alt_allele": torch.stack([encode_allele(t.alternate_allele, l_max) for t in tokens]) if tokens else torch.empty((0, l_max), dtype=torch.long),
    }


def materialize_variant_dataset(
    *,
    output_dir: Path,
    regions: List[Region],
    samples: List[SampleMetadata],
    classes: Sequence[str],
    target: str,
    vcf_pattern: Optional[str],
    vcf_root_dir: Optional[str],
    l_max: int,
    max_indel_size: int,
    unphased_policy: str,
    sample_batch_size: int,
    train_split: float,
    val_split: float,
    test_split: float,
    random_seed: Optional[int],
    family_split_mode: str,
    max_samples: Optional[int] = None,
    central_window_size: Optional[int] = None,
    dataset_vcf_sources: Optional[Dict[str, str]] = None,
) -> Path:
    output_dir = Path(output_dir)
    samples_dir = output_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)
    class_to_idx = {name: idx for idx, name in enumerate(classes)}
    gene_to_idx = {gene: idx for idx, gene in enumerate(sorted({r.gene_id for r in regions}))}
    if central_window_size is not None:
        centered_regions: List[Region] = []
        for region in regions:
            region_len = int(region.end) - int(region.start) + 1
            width = min(int(central_window_size), region_len)
            start = int(region.start) + max((region_len - width) // 2, 0)
            end = start + width - 1
            centered_regions.append(Region(chrom=region.chrom, start=start, end=end, gene_id=region.gene_id))
        regions = centered_regions
    if max_samples is not None:
        samples = samples[:max_samples]
    sample_by_id = {sample.sample_id: sample for sample in samples}
    sample_ids = [sample.sample_id for sample in samples]
    batches = [sample_ids[i:i + sample_batch_size] for i in range(0, len(sample_ids), sample_batch_size)]
    region_vcfs = {
        (region.chrom, region.start, region.end, region.gene_id): resolve_region_vcf_path(region, vcf_pattern, vcf_root_dir, dataset_vcf_sources)
        for region in regions
    }
    sample_index = []
    with Progress(SpinnerColumn(), TextColumn("{task.description}"), BarColumn(), TextColumn("{task.completed}/{task.total}"), TimeElapsedColumn(), console=console) as progress:
        task = progress.add_task("Extraindo variantes VCF...", total=len(regions) * len(batches))
        for batch_ids in batches:
            token_store: Dict[str, List[VariantToken]] = {sid: [] for sid in batch_ids}
            for region in regions:
                vcf_path = region_vcfs[(region.chrom, region.start, region.end, region.gene_id)]
                for row in stream_region_variants(vcf_path, batch_ids, region):
                    ref = str(row["ref"]).upper()
                    alt_raw = str(row["alt"]).upper()
                    if not ref or ref.startswith("<"):
                        continue
                    for sample_idx, gt in enumerate(row.get("gts", [])):
                        if sample_idx >= len(batch_ids):
                            break
                        parsed = _parse_gt(gt, unphased_policy)
                        if parsed is None:
                            continue
                        for haplotype, token in (("H1", parsed[0]), ("H2", parsed[1])):
                            alt = _allele_for_token(ref, alt_raw, token)
                            if alt is None or alt == ref:
                                continue
                            variant_type = classify_variant(ref, alt)
                            length = len(alt) - len(ref)
                            token_store[batch_ids[sample_idx]].append(VariantToken(
                                chrom=str(row["chrom"]),
                                position=int(row["pos"]),
                                position_relative=int(row["pos"]) - int(region.start),
                                gene_id=region.gene_id,
                                haplotype=haplotype,
                                reference_allele=ref,
                                alternate_allele=alt,
                                variant_type=variant_type,
                                length=length,
                            ))
                progress.advance(task)

            for sid, tokens in token_store.items():
                tokens = sort_tokens(tokens)
                sample_meta = sample_by_id[sid]
                tensors = _tokens_to_tensors(tokens, gene_to_idx, l_max, max_indel_size)
                tensors.update({
                    "sample_id": sid,
                    "target": torch.tensor(class_to_idx[sample_meta.target], dtype=torch.long),
                    "num_tokens": len(tokens),
                })
                torch.save(tensors, samples_dir / f"{sid}.pt")
                sample_index.append({
                    "sample_id": sid,
                    "target": sample_meta.target,
                    "target_idx": class_to_idx[sample_meta.target],
                    "num_tokens": len(tokens),
                    "family_id": sample_meta.family_id,
                    "population": sample_meta.population,
                    "superpopulation": sample_meta.superpopulation,
                    "path": f"samples/{sid}.pt",
                })

    splits, split_info = _split_sample_groups(samples, train_split, val_split, test_split, random_seed, family_split_mode)
    metadata = {
        "format_version": 1,
        "target": target,
        "classes": list(classes),
        "class_to_idx": class_to_idx,
        "num_samples": len(samples),
        "num_genes": len(gene_to_idx),
        "l_max": l_max,
        "max_indel_size": max_indel_size,
        "central_window_size": central_window_size,
        "regions": [region.__dict__ for region in regions],
        "vcf_sources": {"|".join(map(str, key)): value for key, value in region_vcfs.items()},
        "splits": {key: len(value) for key, value in splits.items()},
        "split_info": split_info,
    }
    _write_json(output_dir / "metadata.json", metadata)
    _write_json(output_dir / "gene_vocab.json", gene_to_idx)
    _write_json(output_dir / "sample_index.json", {"samples": sample_index})
    _write_json(output_dir / "splits.json", splits)
    console.print(f"[bold green]Dataset materializado:[/bold green] {output_dir}")
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description="Materializa tokens esparsos de variantes a partir de VCFs faseados")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target", default="superpopulation")
    parser.add_argument("--classes", nargs="+", default=["AFR", "AMR", "EAS", "EUR", "SAS"])
    parser.add_argument("--regions-bed", type=Path)
    parser.add_argument("--dataset-dir", type=Path, help="Dataset canonico para ler metadados, regioes e fontes VCF")
    parser.add_argument("--genes", nargs="*", default=None)
    parser.add_argument("--samples-metadata", type=Path)
    parser.add_argument("--vcf-pattern")
    parser.add_argument("--vcf-root-dir")
    parser.add_argument("--l-max", type=int, default=16)
    parser.add_argument("--max-indel-size", type=int, default=50)
    parser.add_argument("--unphased-policy", choices=["skip", "error"], default="skip")
    parser.add_argument("--sample-batch-size", type=int, default=128)
    parser.add_argument("--train-split", type=float, default=0.7)
    parser.add_argument("--val-split", type=float, default=0.15)
    parser.add_argument("--test-split", type=float, default=0.15)
    parser.add_argument("--random-seed", type=int, default=13)
    parser.add_argument("--family-split-mode", choices=["family_aware", "ignore"], default="family_aware")
    parser.add_argument("--max-samples", type=int, default=None, help="Limita amostras para smoke tests/debug")
    parser.add_argument("--central-window-size", type=int, default=None, help="Se definido, materializa apenas a janela central N bp de cada regiao")
    args = parser.parse_args()

    genes = set(args.genes or []) or None
    if args.regions_bed:
        regions = load_regions(args.regions_bed)
        if genes:
            regions = [region for region in regions if region.gene_id in genes]
    elif args.dataset_dir:
        regions = load_regions_from_dataset(args.dataset_dir, genes=genes)
    else:
        raise ValueError("Informe --regions-bed ou --dataset-dir")
    samples = load_sample_metadata(args.samples_metadata, args.dataset_dir, args.target, args.classes)
    dataset_vcf_sources = load_vcf_sources_from_dataset(args.dataset_dir)
    materialize_variant_dataset(
        output_dir=Path(args.output_dir),
        regions=regions,
        samples=samples,
        classes=args.classes,
        target=args.target,
        vcf_pattern=args.vcf_pattern,
        vcf_root_dir=args.vcf_root_dir,
        l_max=args.l_max,
        max_indel_size=args.max_indel_size,
        unphased_policy=args.unphased_policy,
        sample_batch_size=max(1, args.sample_batch_size),
        train_split=args.train_split,
        val_split=args.val_split,
        test_split=args.test_split,
        random_seed=args.random_seed,
        family_split_mode=args.family_split_mode,
        max_samples=args.max_samples,
        central_window_size=args.central_window_size,
        dataset_vcf_sources=dataset_vcf_sources,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
