"""Build a canonical-layout dataset from a phased VCF and free-form sample metadata.

Produces exactly the layout the 1000 Genomes builder writes (so the visualizer, AlphaGenome
prediction, training and the bcftools_chain aligner work unchanged)::

    dataset_metadata.json, layout_metadata.json, import_spec.json, import_report.json
    references/windows/<gene>/ref.window.fa + window_metadata.json
    individuals/<sample>/individual_metadata.json
    individuals/<sample>/windows/<gene>/<sample>.window.vcf.gz(.tbi)
    individuals/<sample>/windows/<gene>/<sample>.window.consensus_ready.vcf.gz(.tbi)
    individuals/<sample>/windows/<gene>/<sample>.H{1,2}.window.raw.fa / .fixed.fa

Windows are centred on each gene exactly like ``build_window_and_predict`` (AlphaGenome's
``Interval.resize``) and per-sample VCFs/consensus sequences are made with the same bcftools
commands. Each gene's region is extracted from the input VCF once for all samples, then split per
sample in parallel. Re-running resumes: finished sample windows are skipped, and new genes or
samples are merged into the existing metadata.

``"extend": true`` adds windows to an existing dataset (an import or the 1000 Genomes dataset):
only the window lists change (``genes``, ``window_catalog``, ``gene_strands`` and each sample's
``windows``); samples, their metadata and the dataset's provenance are left as they are, and the
request is appended to ``window_extensions.json``.

The VCF and the reference may be URLs (``https://``, ``ftp://``, ``s3://``, ``gs://``): bcftools and
samtools then fetch only each window's region through the remote ``.tbi``/``.csi``/``.fai``
index, which htslib caches in :func:`remote_index_dir`. ``vcf_overrides`` maps chromosomes to
files that do not follow the ``{chrom}`` pattern (e.g. a differently named chrX).

Predictions are not made here; run ``genomics alphagenome predict-dataset`` (or the visualizer's
Jobs page) afterwards.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import re
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from genomics.workflows.dataset_builders.vcf_import.metadata import (
    MAX_CATEGORIES,
    describe_fields,
    guess_id_column,
    guess_role,
    FAMILY_COLUMN_NAMES,
    SEX_COLUMN_NAMES,
    parse_table,
    records_by_sample,
)

ALPHAGENOME_WINDOW_SIZES = (16384, 131072, 524288, 1048576)
REQUIRED_TOOLS = ("bcftools", "samtools")
CONSENSUS_EXCLUDE = 'ALT~"<" && ALT!="<DEL>" && ALT!="<NON_REF>"'
REMOTE_RE = re.compile(r"^(https?|ftp|s3|gs)://", re.IGNORECASE)
CHROMOSOMES = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"])


class ImportSpecError(RuntimeError):
    """A problem with the import inputs, reported to the user as-is."""


def emit_progress(fraction: float, message: str) -> None:
    print(f"@@progress {max(0.0, min(1.0, fraction)):.4f} {message}", flush=True)


def now_iso() -> str:
    return _dt.datetime.now().isoformat()


def tool_path(name: str) -> Optional[str]:
    """A tool on PATH, else next to this interpreter (an un-activated conda env)."""
    found = shutil.which(name)
    if found:
        return found
    sibling = Path(sys.executable).parent / name
    return str(sibling) if sibling.exists() else None


def is_remote(path: Any) -> bool:
    return bool(REMOTE_RE.match(str(path or "")))


def remote_index_dir() -> Path:
    """Where htslib keeps the indexes of remote files (it saves them in the working directory)."""
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    path = Path(base) / "genomics" / "remote_index"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cwd_for(cmd: Sequence[str]) -> Optional[str]:
    return str(remote_index_dir()) if any(is_remote(c) for c in cmd) else None


def remote_exists(url: str, timeout: float = 20.0) -> bool:
    """HEAD request for http(s) URLs; other schemes are assumed to exist."""
    if not url.lower().startswith(("http://", "https://")):
        return True
    import urllib.error
    import urllib.request

    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=timeout) as response:
            return 200 <= response.status < 400
    except (urllib.error.URLError, OSError, ValueError):
        return False


def path_exists(path: str) -> bool:
    return remote_exists(path) if is_remote(path) else Path(path).expanduser().exists()


def reference_fai(reference: str) -> Path:
    """Local ``.fai`` of a reference FASTA; for a URL it is downloaded once into the index cache."""
    if not is_remote(reference):
        return Path(f"{Path(reference).expanduser()}.fai")
    local = remote_index_dir() / f"{reference.rstrip('/').rsplit('/', 1)[-1]}.fai"
    if not local.exists():
        import urllib.request

        tmp = local.with_suffix(".fai.tmp")
        try:
            with urllib.request.urlopen(f"{reference}.fai", timeout=60) as response, open(tmp, "wb") as handle:
                shutil.copyfileobj(response, handle)
        except OSError as exc:
            tmp.unlink(missing_ok=True)
            raise ImportSpecError(f"Could not download the reference index {reference}.fai: {exc}")
        os.replace(tmp, local)
    return local


def reference_available(reference: str) -> bool:
    return path_exists(reference)


def run(cmd: Sequence[str], stdout=None) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    cmd[0] = tool_path(cmd[0]) or cmd[0]
    proc = subprocess.run(cmd, stdout=stdout or subprocess.PIPE, stderr=subprocess.PIPE, cwd=_cwd_for(cmd))
    if proc.returncode != 0:
        detail = (proc.stderr or b"").decode("utf-8", "replace").strip().splitlines()[-3:]
        raise RuntimeError(f"{' '.join(str(c) for c in cmd[:3])} failed: {' | '.join(detail)}")
    return proc


def missing_tools() -> List[str]:
    return [tool for tool in REQUIRED_TOOLS if tool_path(tool) is None]


# ------------------------------------------------------------------------------------------- VCF
def resolve_vcf_files(vcf: str) -> List[Any]:
    """Existing files for a VCF path or a ``{chrom}`` pattern (URLs are kept as strings)."""
    if is_remote(vcf):
        if "{chrom}" not in vcf:
            return [vcf] if remote_exists(vcf) else []
        candidates = [vcf.replace("{chrom}", c) for c in CHROMOSOMES]
        with ThreadPoolExecutor(max_workers=8) as pool:
            found = list(pool.map(remote_exists, candidates))
        if not any(found):  # contigs named without "chr"
            candidates = [vcf.replace("{chrom}", c[3:]) for c in CHROMOSOMES]
            with ThreadPoolExecutor(max_workers=8) as pool:
                found = list(pool.map(remote_exists, candidates))
        return [c for c, ok in zip(candidates, found) if ok]
    if "{chrom}" in vcf:
        import glob

        return sorted(Path(p) for p in glob.glob(vcf.replace("{chrom}", "*")) if not p.endswith((".tbi", ".csi")))
    path = Path(vcf).expanduser()
    return [path] if path.exists() else []


def vcf_samples(path: Path) -> List[str]:
    out = run(["bcftools", "query", "-l", path]).stdout.decode("utf-8", "replace")
    return [line.strip() for line in out.splitlines() if line.strip()]


def vcf_contigs(path: Path) -> List[str]:
    try:
        out = run(["bcftools", "index", "-s", path]).stdout.decode()
        names = [line.split("\t")[0] for line in out.splitlines() if line.strip()]
        if names:
            return names
    except RuntimeError:
        pass
    header = run(["bcftools", "view", "-h", path]).stdout.decode("utf-8", "replace")
    return [line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">") for line in header.splitlines() if line.startswith("##contig=<ID=")]


def vcf_indexed(path: Any) -> bool:
    if is_remote(path):  # read through the remote index; bcftools reports a missing one
        return True
    return Path(f"{path}.tbi").exists() or Path(f"{path}.csi").exists()


def phasing_sample(path: Path, region: Optional[str] = None, records: int = 300) -> Dict[str, int]:
    """Count phased/unphased heterozygous-capable genotypes in the first records."""
    cmd = ["bcftools", "query", "-f", "[%GT\\t]\\n"]
    if region:
        cmd += ["-r", region]
    cmd[0] = tool_path("bcftools") or "bcftools"
    proc = subprocess.Popen([*cmd, str(path)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=_cwd_for([str(path)]))
    counts = {"phased": 0, "unphased": 0, "records": 0}
    try:
        for line in proc.stdout:
            for gt in line.split("\t"):
                if "|" in gt:
                    counts["phased"] += 1
                elif "/" in gt and gt not in ("./.",):
                    counts["unphased"] += 1
            counts["records"] += 1
            if counts["records"] >= records:
                break
    finally:
        proc.kill()
        proc.wait()
    return counts


def inspect_vcf(vcf: str) -> Dict[str, Any]:
    """Samples, contigs, index and phasing of a VCF (or ``{chrom}`` pattern) for the import form."""
    missing = missing_tools()
    if missing:
        raise ImportSpecError(f"Missing tools on PATH: {', '.join(missing)} (install bcftools/samtools, e.g. conda install -c bioconda bcftools samtools)")
    files = resolve_vcf_files(vcf)
    if not files:
        raise ImportSpecError(f"No VCF found at {vcf}")
    first = files[0]
    samples = vcf_samples(first)
    warnings = []
    unindexed = [str(f) for f in files if not vcf_indexed(f)]
    if unindexed:
        warnings.append(f"{len(unindexed)} file(s) have no .tbi/.csi index; they will be indexed during the import (needs write access)")
    contigs = vcf_contigs(first) if vcf_indexed(first) else []
    phasing = phasing_sample(first)
    if phasing["unphased"] and not phasing["phased"]:
        warnings.append("Genotypes look unphased: H1/H2 haplotypes will follow allele order, not true phase")
    elif phasing["unphased"]:
        warnings.append("Some genotypes are unphased; those use allele order for H1/H2")
    return {
        "files": [str(f) for f in files[:50]],
        "file_count": len(files),
        "samples": samples,
        "sample_count": len(samples),
        "contigs": contigs[:60],
        "chr_prefix": any(c.startswith("chr") for c in contigs) if contigs else None,
        "phasing": phasing,
        "warnings": warnings,
    }


# --------------------------------------------------------------------------------------- windows
@dataclass
class Target:
    name: str
    chromosome: str  # FASTA naming
    start: int  # window, 0-based inclusive
    end: int  # window, 0-based exclusive
    strand: str = "+"
    feature_start: Optional[int] = None  # gene/region, 0-based
    feature_end: Optional[int] = None
    kind: str = "gene"


def centred_window(start: int, end: int, strand: str, size: int) -> Tuple[int, int]:
    """AlphaGenome ``Interval.resize`` on a 0-based half-open interval."""
    width = end - start
    negative = strand == "-"
    centre = (start + end) // 2 + (0 if negative else width % 2)
    if negative:
        return centre - size // 2, centre + (size + 1) // 2
    return centre - (size + 1) // 2, centre + size // 2


def fasta_chrom_names(fai: Path) -> List[str]:
    with open(fai, "r", encoding="utf-8") as handle:
        return [line.split("\t", 1)[0] for line in handle if line.strip()]


def match_chrom(name: str, available: Iterable[str]) -> Optional[str]:
    available = list(available)
    candidates = [name, name[3:] if name.startswith("chr") else f"chr{name}"]
    if name in ("chrM", "MT", "M"):
        candidates += ["chrM", "MT", "M"]
    for candidate in candidates:
        if candidate in available:
            return candidate
    return None


def load_gene_table(path: Path):
    import pandas as pd

    columns = ["Chromosome", "Feature", "Start", "End", "Strand", "gene_name", "gene_id"]
    suffix = path.suffix.lower()
    if suffix == ".feather":
        table = pd.read_feather(path, columns=columns)
    elif suffix == ".parquet":
        table = pd.read_parquet(path, columns=columns)
    else:
        table = pd.read_csv(path, usecols=columns, sep=None, engine="python")
    return table[table["Feature"] == "gene"]


def default_gtf_candidates(output_dir: Path) -> List[Path]:
    candidates = [output_dir / "gtf_cache.feather"]
    try:
        from genomics.workspace import DEFAULT_DATASET_DIR

        candidates.append(Path(DEFAULT_DATASET_DIR) / "gtf_cache.feather")
    except Exception:
        pass
    return candidates


def resolve_targets(spec: Dict[str, Any], fasta_names: List[str], gtf_path: Optional[Path]) -> List[Target]:
    size = int(spec["window_size"])
    targets: List[Target] = []
    genes = [str(g).strip() for g in spec.get("genes") or [] if str(g).strip()]
    if genes:
        if gtf_path is None:
            raise ImportSpecError("Gene symbols need a GTF table (gtf_cache.feather); set 'gtf' or give explicit regions")
        table = load_gene_table(gtf_path)
        for gene in genes:
            column = "gene_id" if gene.upper().startswith("ENSG") else "gene_name"
            if column == "gene_id":
                rows = table[table["gene_id"].astype(str).str.split(".").str[0] == gene.split(".")[0]]
            else:
                rows = table[table["gene_name"] == gene]
            if rows.empty:
                raise ImportSpecError(f"Gene not found in {gtf_path.name}: {gene}")
            if len(rows) > 1:
                primary = rows[~rows["Chromosome"].astype(str).str.contains("_")]
                rows = primary if len(primary) == 1 else rows
            if len(rows) > 1:
                raise ImportSpecError(f"Gene {gene} has {len(rows)} entries in {gtf_path.name}; use its Ensembl id or an explicit region")
            row = rows.iloc[0]
            chrom = match_chrom(str(row["Chromosome"]), fasta_names)
            if chrom is None:
                raise ImportSpecError(f"Chromosome {row['Chromosome']} of {gene} is not in the reference FASTA")
            strand = str(row["Strand"])
            start, end = centred_window(int(row["Start"]), int(row["End"]), strand, size)
            targets.append(Target(gene, chrom, start, end, strand, int(row["Start"]), int(row["End"])))
    for region in spec.get("regions") or []:
        name = str(region.get("name") or "").strip()
        if not name:
            raise ImportSpecError("Every region needs a name")
        chrom = match_chrom(str(region["chrom"]), fasta_names)
        if chrom is None:
            raise ImportSpecError(f"Chromosome {region['chrom']} of region {name} is not in the reference FASTA")
        feature_start = int(region["start"]) - 1
        feature_end = int(region.get("end") or region["start"])
        start, end = centred_window(feature_start, feature_end, "+", size)
        targets.append(Target(name, chrom, start, end, "+", feature_start, feature_end, "region"))
    names = [t.name for t in targets]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ImportSpecError(f"Duplicate window names: {', '.join(duplicates)}")
    if not targets:
        raise ImportSpecError("No genes or regions to import")
    return targets


# ------------------------------------------------------------------------------------- building
class DatasetImporter:
    def __init__(self, spec: Dict[str, Any]):
        self.spec = dict(spec)
        self.out = Path(self.spec["output_dir"]).expanduser().resolve()
        self.vcf = str(self.spec["vcf"])
        self.vcf_overrides = {str(k): str(v) for k, v in (self.spec.get("vcf_overrides") or {}).items()}
        reference = str(self.spec["reference_fasta"])
        self.reference = reference if is_remote(reference) else str(Path(reference).expanduser().resolve())
        self.window_size = int(self.spec.get("window_size") or 524288)
        self.workers = max(1, int(self.spec.get("workers") or min(8, os.cpu_count() or 4)))
        self.haplotypes = [h for h in (self.spec.get("haplotypes") or ["H1", "H2"]) if h in ("H1", "H2")] or ["H1", "H2"]
        self.variant_filter = self.spec.get("variant_filter") or "all"
        self.keep_raw = bool(self.spec.get("keep_raw", True))
        self.extend = bool(self.spec.get("extend"))
        self.failures: List[Dict[str, str]] = []
        self._lock = threading.Lock()
        self._done = 0
        self._total = 1

    # -- inputs ---------------------------------------------------------------------------------
    def _metadata(self, samples: List[str]) -> Dict[str, Dict[str, Any]]:
        inline = self.spec.get("metadata")
        if isinstance(inline, dict):
            return {str(k): {str(f): v for f, v in (v or {}).items() if v is not None} for k, v in inline.items()}
        path = self.spec.get("metadata_file")
        if not path:
            return {}
        path = Path(path).expanduser()
        columns, rows = parse_table(path.read_text(encoding="utf-8"), path.name)
        id_column = self.spec.get("id_column") or guess_id_column(columns, rows, samples)
        if id_column not in columns:
            raise ImportSpecError(f"ID column {id_column!r} not in {path.name} (columns: {', '.join(columns)})")
        return records_by_sample(
            rows,
            id_column,
            family_column=self.spec.get("family_column") or guess_role(columns, FAMILY_COLUMN_NAMES),
            sex_column=self.spec.get("sex_column") or guess_role(columns, SEX_COLUMN_NAMES),
        )

    def _check_inputs(self) -> None:
        missing = missing_tools()
        if missing:
            raise ImportSpecError(f"Missing tools on PATH: {', '.join(missing)}")
        if not reference_available(self.reference):
            raise ImportSpecError(f"Reference FASTA not found: {self.reference}")
        if is_remote(self.reference):
            print(f"[INFO] Reading reference windows from {self.reference}", flush=True)
        elif not reference_fai(self.reference).exists():
            print(f"[INFO] Indexing reference {self.reference}", flush=True)
            run(["samtools", "faidx", self.reference])
        self.fai = reference_fai(self.reference)
        if self.window_size not in ALPHAGENOME_WINDOW_SIZES:
            print(f"[WARN] Window size {self.window_size} is not an AlphaGenome input length {ALPHAGENOME_WINDOW_SIZES}; predictions will not be possible", flush=True)

    def _vcf_for(self, chrom: str) -> Tuple[Any, str]:
        """(VCF file, contig name inside it) for a FASTA chromosome name."""
        names = (chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}")
        override = next((self.vcf_overrides[n] for n in names if n in self.vcf_overrides), None)
        if override:
            path = override if is_remote(override) else Path(override).expanduser()
        elif "{chrom}" in self.vcf:
            for candidate in names:
                path = self.vcf.replace("{chrom}", candidate)
                if path_exists(path):
                    path = path if is_remote(path) else Path(path).expanduser()
                    break
            else:
                raise ImportSpecError(f"No VCF for {chrom}: {self.vcf.replace('{chrom}', chrom)}")
        else:
            path = self.vcf if is_remote(self.vcf) else Path(self.vcf).expanduser()
        if not vcf_indexed(path):
            print(f"[INFO] Indexing {path}", flush=True)
            run(["bcftools", "index", "-t", path])
        contig = match_chrom(chrom, vcf_contigs(path))
        if contig is None:
            raise ImportSpecError(f"Chromosome {chrom} not found in {str(path).rsplit('/', 1)[-1]}")
        return path, contig

    # -- per window -------------------------------------------------------------------------
    def _reference_window(self, target: Target, fai_sizes: Dict[str, int]) -> Tuple[Path, Dict[str, Any]]:
        ref_dir = self.out / "references" / "windows" / target.name
        ref_dir.mkdir(parents=True, exist_ok=True)
        chrom_size = fai_sizes.get(target.chromosome)
        actual_start = max(0, target.start)
        actual_end = min(target.end, chrom_size) if chrom_size else target.end
        region = f"{target.chromosome}:{actual_start + 1}-{actual_end}"
        ref_path = ref_dir / "ref.window.fa"
        if not ref_path.exists():
            with open(ref_path.with_suffix(".fa.tmp"), "wb") as handle:
                run(["samtools", "faidx", self.reference, region], stdout=handle)
            os.replace(ref_path.with_suffix(".fa.tmp"), ref_path)
        meta = {
            "type": target.kind,
            "chromosome": target.chromosome,
            "start": actual_start + 1,
            "end": actual_end,
            "window_size": actual_end - actual_start,
            "requested_start": target.start,
            "requested_end": target.end,
            "requested_window_size": self.window_size,
            "strand": target.strand,
            "feature_start": target.feature_start + 1 if target.feature_start is not None else None,
            "feature_end": target.feature_end,
            "outputs": [],
            "ontologies": [],
            "added_at": now_iso(),
            # The window's cohort VCF (contigs named like the reference) is what the bcftools_chain
            # aligner streams variants from; the original VCF is kept for provenance.
            "raw_variant_source": {"chromosome": target.chromosome, "vcf_path": str(self.cohort_vcf_path(target)), "vcf_pattern": self.vcf},
        }
        meta_path = ref_dir / "window_metadata.json"
        if meta_path.exists():
            try:
                previous = json.loads(meta_path.read_text(encoding="utf-8"))
                meta["outputs"] = previous.get("outputs") or []
                meta["ontologies"] = previous.get("ontologies") or []
                meta["added_at"] = previous.get("added_at") or meta["added_at"]
            except ValueError:
                pass
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        return ref_path, {"region": region, "actual_start": actual_start, "actual_end": actual_end}

    def cohort_vcf_path(self, target: Target) -> Path:
        return self.out / "references" / "windows" / target.name / "cohort.window.vcf.gz"

    def _cohort_vcf(self, target: Target, samples: List[str]) -> Path:
        """The window's variants for every dataset sample, contigs named like the reference.

        Per-sample VCFs are split from it, and the bcftools_chain training aligner reads the
        cohort's indels from it (``raw_variant_source.vcf_path``), so the dataset does not depend
        on the original VCF staying where it was.
        """
        path = self.cohort_vcf_path(target)
        listing = path.with_name("cohort.window.samples.txt")
        wanted = "\n".join(sorted(samples)) + "\n"
        if Path(f"{path}.tbi").exists() and listing.exists() and listing.read_text(encoding="utf-8") == wanted:
            return path
        work = self.out / ".import_tmp" / target.name
        work.mkdir(parents=True, exist_ok=True)
        vcf_path, contig = self._vcf_for(target.chromosome)
        region = f"{contig}:{max(0, target.start) + 1}-{target.end}"
        samples_file = work / "samples.txt"
        samples_file.write_text("\n".join(samples) + "\n", encoding="utf-8")
        region_vcf = work / "region.vcf.gz"
        cmd = ["bcftools", "view", "-r", region, "-S", samples_file, "--force-samples"]
        if self.variant_filter == "snps":
            cmd += ["-v", "snps"]
        cmd += ["-Oz", "-o", region_vcf, vcf_path]
        run(cmd)
        if contig != target.chromosome:
            renamed = work / "region.renamed.vcf.gz"
            mapping = work / "chroms.txt"
            mapping.write_text(f"{contig}\t{target.chromosome}\n", encoding="utf-8")
            run(["bcftools", "annotate", "--rename-chrs", mapping, "-Oz", "-o", renamed, region_vcf])
            os.replace(renamed, region_vcf)
        run(["bcftools", "index", "-t", "-f", region_vcf])
        os.replace(region_vcf, path)
        os.replace(f"{region_vcf}.tbi", f"{path}.tbi")
        listing.write_text(wanted, encoding="utf-8")
        shutil.rmtree(work, ignore_errors=True)
        return path

    def _sample_window(self, sample: str, target: Target, ref_path: Path, ref_seq: str, bounds: Dict[str, Any], region_vcf: Path) -> None:
        from genomics.workflows.dataset_builders.non_longevous.build_window_and_predict import (
            adjust_to_target_size,
            read_fasta_seq_only,
            write_fasta,
        )

        case = self.out / "individuals" / sample / "windows" / target.name
        fixed = {h: case / f"{sample}.{h}.window.fixed.fa" for h in self.haplotypes}
        vcf_window = case / f"{sample}.window.vcf.gz"
        vcf_cons = case / f"{sample}.window.consensus_ready.vcf.gz"
        if all(p.exists() for p in fixed.values()) and Path(f"{vcf_cons}.tbi").exists():
            return
        case.mkdir(parents=True, exist_ok=True)
        if not Path(f"{vcf_window}.tbi").exists():
            run(["bcftools", "view", "-s", sample, "-Oz", "-o", vcf_window, region_vcf])
            run(["bcftools", "index", "-t", "-f", vcf_window])
        if not Path(f"{vcf_cons}.tbi").exists():
            run(["bcftools", "view", "-e", CONSENSUS_EXCLUDE, "-Oz", "-o", vcf_cons, vcf_window])
            run(["bcftools", "index", "-t", "-f", vcf_cons])
        for hap in self.haplotypes:
            if fixed[hap].exists():
                continue
            raw = case / f"{sample}.{hap}.window.raw.fa"
            proc = run(["bcftools", "consensus", "-H", hap[1], "-f", ref_path, vcf_cons])
            raw.write_bytes(proc.stdout)
            seq = adjust_to_target_size(
                read_fasta_seq_only(raw), ref_seq, self.window_size,
                target.start, target.end, bounds["actual_start"], bounds["actual_end"],
            )
            tmp = fixed[hap].with_suffix(".tmp")
            write_fasta(seq, tmp, header=f"{sample}_{hap}_{target.name}")
            os.replace(tmp, fixed[hap])
            if not self.keep_raw:
                raw.unlink(missing_ok=True)

    def _tick(self, message: str) -> None:
        with self._lock:
            self._done += 1
            if self._done % 5 == 0 or self._done == self._total:
                emit_progress(0.02 + 0.93 * self._done / self._total, message)

    def build_target(self, target: Target, samples: List[str], fai_sizes: Dict[str, int], cohort: Optional[List[str]] = None) -> List[str]:
        from genomics.workflows.dataset_builders.non_longevous.build_window_and_predict import adjust_to_target_size, read_fasta_seq_only

        ref_path, bounds = self._reference_window(target, fai_sizes)
        ref_raw = read_fasta_seq_only(ref_path)
        ref_seq = adjust_to_target_size(ref_raw, ref_raw, self.window_size, target.start, target.end, bounds["actual_start"], bounds["actual_end"], use_ref_only=True)
        todo = [s for s in samples if not all((self.out / "individuals" / s / "windows" / target.name / f"{s}.{h}.window.fixed.fa").exists() for h in self.haplotypes)]
        built = [s for s in samples if s not in set(todo)]
        print(f"[INFO] {target.name}: {bounds['region']} ({len(todo)} samples to build)", flush=True)
        region_vcf = self._cohort_vcf(target, sorted(set(cohort or []) | set(samples)))
        for _ in built:
            self._tick(f"{target.name}: {len(built)} sample windows already built")
        if not todo:
            return samples
        ok = list(built)
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            futures = {pool.submit(self._sample_window, s, target, ref_path, ref_seq, bounds, region_vcf): s for s in todo}
            for future in as_completed(futures):
                sample = futures[future]
                try:
                    future.result()
                    ok.append(sample)
                except Exception as exc:
                    with self._lock:
                        self.failures.append({"sample": sample, "window": target.name, "error": str(exc)[:500]})
                    print(f"[ERROR] {sample}/{target.name}: {exc}", flush=True)
                self._tick(f"{target.name}: {sample}")
        return ok

    # -- metadata -------------------------------------------------------------------------------
    def write_metadata(self, samples: List[str], metadata: Dict[str, Dict[str, Any]], targets: List[Target], built: Dict[str, List[str]], individuals_files: bool = True) -> None:
        meta_path = self.out / "dataset_metadata.json"
        previous: Dict[str, Any] = {}
        if meta_path.exists():
            try:
                previous = json.loads(meta_path.read_text(encoding="utf-8"))
            except ValueError:
                previous = {}
        windows_by_sample: Dict[str, List[str]] = {}
        for gene, members in built.items():
            for sample in members:
                windows_by_sample.setdefault(sample, []).append(gene)
        pedigree = dict(previous.get("individuals_pedigree") or {})
        for sample in samples:
            pedigree[sample] = {**pedigree.get(sample, {}), **metadata.get(sample, {})}
        individuals = list(dict.fromkeys([*(previous.get("individuals") or []), *[s for s in samples if s in windows_by_sample]]))
        for sample in individuals if individuals_files else []:
            ind_dir = self.out / "individuals" / sample
            ind_path = ind_dir / "individual_metadata.json"
            existing: Dict[str, Any] = {}
            if ind_path.exists():
                try:
                    existing = json.loads(ind_path.read_text(encoding="utf-8"))
                except ValueError:
                    existing = {}
            windows = sorted(set(existing.get("windows") or []) | set(windows_by_sample.get(sample, [])))
            if self.extend:
                if existing and sample in windows_by_sample:
                    ind_path.write_text(json.dumps({**existing, "windows": windows}, indent=2, default=str), encoding="utf-8")
                continue
            record = {**existing, "sample_id": sample, **pedigree.get(sample, {}), "windows": windows}
            record.setdefault("family_id", sample)
            ind_dir.mkdir(parents=True, exist_ok=True)
            ind_path.write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
        catalog = dict(previous.get("window_catalog") or {})
        strands = dict(previous.get("gene_strands") or {})
        for target in targets:
            if target.name not in built:
                continue
            window_meta = json.loads((self.out / "references" / "windows" / target.name / "window_metadata.json").read_text(encoding="utf-8"))
            catalog[target.name] = {**catalog.get(target.name, {}), **window_meta}
            strands[target.name] = target.strand
        if self.extend:
            payload = {
                **previous,
                "last_updated": now_iso(),
                "genes": sorted(set(str(g) for g in previous.get("genes") or []) | set(catalog)),
                "window_types": sorted({str(v.get("type") or "gene") for v in catalog.values()}),
                "gene_strands": strands,
                "window_catalog": catalog,
            }
            tmp = meta_path.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
            os.replace(tmp, meta_path)
            return
        records = {s: pedigree.get(s, {}) for s in individuals}
        fields = describe_fields(records)
        payload = {
            **previous,
            "dataset_name": self.spec.get("name") or previous.get("dataset_name") or self.out.name,
            "creation_date": previous.get("creation_date") or now_iso(),
            "last_updated": now_iso(),
            "total_individuals": len(individuals),
            "window_size": self.window_size,
            "individuals": individuals,
            "individuals_pedigree": {s: pedigree.get(s, {}) for s in individuals},
            "genes": sorted(catalog),
            "window_types": sorted({str(v.get("type") or "gene") for v in catalog.values()}),
            "alphagenome_outputs": previous.get("alphagenome_outputs") or [],
            "ontologies": previous.get("ontologies") or [],
            "ontology_details": previous.get("ontology_details") or {},
            "gene_strands": strands,
            "layout_version": 1,
            "window_catalog": catalog,
            "raw_variant_sources": {"vcf_pattern": self.vcf, "reference_fasta": str(self.reference), **({"vcf_overrides": self.vcf_overrides} if self.vcf_overrides else {})},
            "source": {"builder": "vcf_import", "spec": "import_spec.json"},
            "metadata_fields": [{"name": f["name"], "kind": f["kind"]} for f in fields],
        }
        for key in [k for k in payload if k.endswith("_distribution")]:
            payload.pop(key)
        for field in fields:
            if field["kind"] == "categorical" and field["distinct"] <= MAX_CATEGORIES:
                payload[f"{field['name']}_distribution"] = dict(field["counts"])
        tmp = meta_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
        os.replace(tmp, meta_path)
        layout = self.out / "layout_metadata.json"
        if not layout.exists():
            layout.write_text(json.dumps({"layout_version": 1, "source": "vcf_import", "created_at": now_iso()}, indent=2), encoding="utf-8")

    def _record_extension(self) -> None:
        """Log the added windows; an imported dataset's spec also lists them, for re-runs."""
        added = {"genes": list(self.spec.get("genes") or []), "regions": list(self.spec.get("regions") or [])}
        log_path = self.out / "window_extensions.json"
        try:
            log = json.loads(log_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            log = []
        log.append({"added_at": now_iso(), **added, "vcf": self.vcf, "reference_fasta": str(self.reference), "window_size": self.window_size})
        log_path.write_text(json.dumps(log, indent=2, default=str), encoding="utf-8")
        spec_path = self.out / "import_spec.json"
        if spec_path.exists():
            try:
                spec = json.loads(spec_path.read_text(encoding="utf-8"))
            except ValueError:
                return
            spec["genes"] = list(dict.fromkeys([*(spec.get("genes") or []), *added["genes"]]))
            names = {r.get("name") for r in spec.get("regions") or []}
            spec["regions"] = [*(spec.get("regions") or []), *[r for r in added["regions"] if r.get("name") not in names]]
            spec_path.write_text(json.dumps(spec, indent=2, default=str), encoding="utf-8")

    # -- main -----------------------------------------------------------------------------------
    def run(self) -> Dict[str, Any]:
        self._check_inputs()
        self.out.mkdir(parents=True, exist_ok=True)
        if self.extend:
            self._record_extension()
        else:
            (self.out / "import_spec.json").write_text(json.dumps({k: v for k, v in self.spec.items() if k != "metadata"}, indent=2, default=str), encoding="utf-8")
        emit_progress(0.0, "Reading VCF samples")
        files = resolve_vcf_files(self.vcf)
        if not files:
            raise ImportSpecError(f"No VCF found at {self.vcf}")
        available = vcf_samples(files[0])
        requested = [str(s) for s in self.spec.get("samples") or []]
        if requested and self.extend:
            # Windows for the dataset's own samples; ones the VCF lacks just get no new window.
            present = set(available)
            missing = [s for s in requested if s not in present]
            if missing:
                print(f"[WARN] {len(missing)} dataset samples are not in the VCF and get no new windows, e.g. {', '.join(missing[:5])}", flush=True)
            samples = [s for s in dict.fromkeys(requested) if s in present]
            if not samples:
                raise ImportSpecError("None of the dataset's samples are in the VCF")
        elif requested:
            unknown = [s for s in requested if s not in set(available)]
            if unknown:
                raise ImportSpecError(f"{len(unknown)} requested samples are not in the VCF, e.g. {', '.join(unknown[:5])}")
            samples = list(dict.fromkeys(requested))
        else:
            samples = available
        metadata = self._metadata(samples)
        fai = self.fai
        fasta_names = fasta_chrom_names(fai)
        fai_sizes = {}
        with open(fai, "r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.split("\t")
                if len(parts) > 1:
                    fai_sizes[parts[0]] = int(parts[1])
        gtf = None
        if self.spec.get("genes"):
            explicit = self.spec.get("gtf")
            candidates = [Path(explicit).expanduser()] if explicit else default_gtf_candidates(self.out)
            gtf = next((c for c in candidates if c.exists()), None)
        targets = resolve_targets(self.spec, fasta_names, gtf)
        print(f"[INFO] {len(samples)} samples x {len(targets)} windows of {self.window_size} bp -> {self.out}", flush=True)
        if gtf is not None and not (self.out / "gtf_cache.feather").exists() and gtf.suffix == ".feather":
            try:
                os.symlink(gtf.resolve(), self.out / "gtf_cache.feather")
            except OSError:
                pass
        self._total = max(1, len(samples) * len(targets))
        built: Dict[str, List[str]] = {}
        try:
            previous_individuals = [str(i) for i in json.loads((self.out / "dataset_metadata.json").read_text(encoding="utf-8")).get("individuals") or []]
        except (OSError, ValueError):
            previous_individuals = []
        for target in targets:
            members = self.build_target(target, samples, fai_sizes, cohort=previous_individuals)
            if members:
                built[target.name] = members
            # Dataset metadata after every window, so a long import can be opened while it runs.
            self.write_metadata(samples, metadata, targets, built, individuals_files=False)
        emit_progress(0.97, "Writing metadata")
        self.write_metadata(samples, metadata, targets, built)
        report = {
            "finished_at": now_iso(),
            "samples": len(samples),
            "windows": {name: len(members) for name, members in built.items()},
            "failures": self.failures[:1000],
            "failure_count": len(self.failures),
        }
        (self.out / "import_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        shutil.rmtree(self.out / ".import_tmp", ignore_errors=True)
        if not built:
            raise ImportSpecError("No window could be built; see the errors above")
        print(f"@@result {json.dumps({'dataset_dir': str(self.out), 'samples': len(samples), 'windows': len(built), 'failures': len(self.failures)})}", flush=True)
        emit_progress(1.0, f"Imported {len(samples)} samples x {len(built)} windows")
        return report


def load_spec(path: Path) -> Dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    if Path(path).suffix.lower() in (".yaml", ".yml"):
        import yaml

        return yaml.safe_load(text) or {}
    return json.loads(text)


def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="genomics dataset-builders vcf-import", description="Build a canonical dataset from a phased VCF and sample metadata")
    parser.add_argument("--spec", type=Path, required=True, help="Import spec (JSON or YAML): output_dir, vcf, reference_fasta, window_size, genes/regions, metadata_file or metadata")
    parser.add_argument("--workers", type=int, default=None, help="Parallel bcftools jobs (default: spec 'workers' or min(8, CPUs))")
    parser.add_argument("--inspect", action="store_true", help="Only print the VCF's samples, contigs and phasing")
    args = parser.parse_args(argv)
    spec = load_spec(args.spec)
    if args.workers:
        spec["workers"] = args.workers
    try:
        if args.inspect:
            print(json.dumps(inspect_vcf(str(spec["vcf"])), indent=2))
            return 0
        for key in ("output_dir", "vcf", "reference_fasta"):
            if not spec.get(key):
                raise ImportSpecError(f"The spec needs '{key}'")
        DatasetImporter(spec).run()
    except ImportSpecError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr, flush=True)
        return 2
    return 0
