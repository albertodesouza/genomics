"""Builds `notebooks/data/` (this notebook's self-contained data bundle) straight from
the public 1000 Genomes Project source data -- the same sources
`genomics.workflows.dataset_builders.non_longevous.build_non_longevous_dataset` and
`build_window_and_predict` use (metadata CSV with pedigree/sex/population, GENCODE GTF,
the reference FASTA, and the per-chromosome phased VCFs) -- rather than from the
already-built `1kG_high_coverage` dataset directory. That means this notebook has no
dependency on that pre-built dataset or on `genomics.*`.

By default the reference FASTA and per-chromosome VCFs are read directly off the public
1000 Genomes FTP (`DEFAULT_REF_FASTA`/`DEFAULT_VCF_PATTERN`, both plain HTTP URLs) --
`bcftools`/`samtools` stream just the small byte ranges each gene's window needs via
HTTP range requests (confirmed: ~5s per 512kb window, no multi-GB download of the whole
reference or chromosome VCF). So on a fresh clone with no local genomics data at all,
this still works end to end -- clone, install deps, run this notebook. If you happen to
have these files locally already (faster, no network), pass local paths instead.

Meant to be driven interactively from `notebooks/prepare_data.ipynb`; also runnable
standalone:

    python notebooks/utils/prepare_data.py

Needs `bcftools`/`samtools` on PATH, and either network access to the public 1000
Genomes FTP or local copies of the reference FASTA and per-chromosome VCFs.

Writes, under `notebooks/data/`:
  experiment.json          -- genes, ontology terms, pigmentation class_map, window size
  genes.json                -- {gene: {chrom, start, end}} (0-based half-open window,
                                centered on the gene's GTF "gene" feature, same as the
                                original dataset build)
  individuals.json          -- {sample_id: {population, superpopulation, sex, family_id}},
                                restricted to samples whose population is in class_map
  references/<gene>/ref.window.fa
  variants/<gene>.vcf.gz(+.tbi)  -- sliced to that gene's window, restricted to the
                                     same individuals, symbolic ALTs other than <DEL>
                                     filtered out (bcftools consensus can't use them)
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = NOTEBOOK_DIR.parent
DATA_DIR = NOTEBOOK_DIR / "data"

DEFAULT_METADATA_CSV = REPO_ROOT / "docs" / "historical" / "1000-genomes" / "1000_genomes_metadata.csv"

# Public 1000 Genomes FTP (EBI mirror) -- reachable over plain HTTP, supports byte-range
# requests, and each file has a remote .fai/.tbi index alongside it, so bcftools/samtools
# can randomly-access just the regions we need without downloading the whole file.
_ONE_KG_FTP = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
DEFAULT_REF_FASTA = f"{_ONE_KG_FTP}/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
DEFAULT_VCF_PATTERN = (
    f"{_ONE_KG_FTP}/data_collections/1000G_2504_high_coverage/working/"
    "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{chrom}."
    "filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)

# Opt-in only: this machine's own pre-built dataset, if present -- lets Step 2/5 reuse its
# exact window boundaries and precomputed AlphaGenome predictions (see import_precomputed.py)
# instead of talking to the FTP or the AlphaGenome API at all. Not needed for a fresh clone.
DEFAULT_PRECOMPUTED_DATASET_DIR = Path("/dados/GENOMICS_DATA/v1/1kG_high_coverage")

GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)

GENES = [
    "MC1R", "TYRP1", "TYR", "SLC45A2", "DDB1",
    "EDAR", "MFSD12", "OCA2", "HERC2", "SLC24A5", "TCHH",
]
ONTOLOGY_TERMS = ["CL:1000458", "CL:0000346", "CL:2000092"]
CLASS_MAP = {
    "strong pigmentation": ["YRI", "ESN", "LWK", "MSL", "GWD"],
    "weak pigmentation": ["FIN", "CEU", "GBR"],
}
WINDOW_CENTER_SIZE = 32768
WINDOW_SIZE = 524288  # alphagenome.models.dna_client.SEQUENCE_LENGTH_500KB

SEX_LABELS = {1: "male", 2: "female"}
METADATA_REQUIRED_COLUMNS = ["FamilyID", "SampleID", "FatherID", "MotherID", "Sex", "Population", "Superpopulation"]


def _run(cmd: list) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, proc.stdout, proc.stderr)
    return proc


def require_tools() -> None:
    missing = [tool for tool in ("bcftools", "samtools") if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(
            f"Missing required tool(s) on PATH: {', '.join(missing)}. "
            "Run `source scripts/env/start_genomics_universal.sh` or otherwise install them."
        )


def is_remote(path) -> bool:
    """True for http(s)/ftp URLs. Never wrap these in `pathlib.Path` -- `Path("http://a/b")`
    collapses the double slash (`http:/a/b`), silently breaking the URL."""
    return str(path).startswith(("http://", "https://", "ftp://"))


# ── Metadata (FamilyID,SampleID,FatherID,MotherID,Sex,Population,Superpopulation) ──

def load_metadata(csv_path: Path = DEFAULT_METADATA_CSV) -> pd.DataFrame:
    """Loads the 1000 Genomes pedigree/population CSV (same format and file
    `build_non_longevous_dataset.py` reads -- Sex: 1=male, 2=female)."""
    if not csv_path.exists():
        raise FileNotFoundError(f"1000 Genomes metadata CSV not found: {csv_path}")
    df = pd.read_csv(csv_path, dtype={"FamilyID": str, "SampleID": str, "FatherID": str, "MotherID": str})
    missing = [c for c in METADATA_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {csv_path}: {missing}")
    return df


def summarize_class_populations(df: pd.DataFrame, class_map: Dict[str, list]) -> pd.DataFrame:
    """Per-(class, population) individual/sex counts, for the populations in class_map."""
    pop_to_class = {pop: cls for cls, pops in class_map.items() for pop in pops}
    subset = df[df["Population"].isin(pop_to_class)].copy()
    subset["class"] = subset["Population"].map(pop_to_class)
    rows = []
    for (class_name, population), group in subset.groupby(["class", "Population"]):
        sex_counts = Counter(group["Sex"])
        rows.append({
            "class": class_name,
            "population": population,
            "n": len(group),
            "male": sex_counts.get(1, 0),
            "female": sex_counts.get(2, 0),
        })
    return pd.DataFrame(rows).sort_values(["class", "population"]).reset_index(drop=True)


# ── Reference genome (chromosome-prefix detection + window extraction) ──

def detect_chr_prefix(ref_fasta) -> str:
    """Reads just the first line of the FASTA's `.fai` index to detect `chr`-prefixed vs.
    bare chromosome naming. Works for both local paths and remote URLs -- for a URL this
    fetches only the small `.fai` file itself, never the FASTA."""
    fai_url_or_path = str(ref_fasta) + ".fai"
    if is_remote(ref_fasta):
        import requests

        resp = requests.get(fai_url_or_path, timeout=30)
        resp.raise_for_status()
        first = resp.text.splitlines()[0]
    else:
        with open(fai_url_or_path) as f:
            first = f.readline()
    first_col = first.strip().split("\t")[0]
    return "chr" if first_col.startswith("chr") else ""


def coerce_chromosome_name(chrom: str, desired_prefix: str) -> str:
    has_chr = chrom.startswith("chr")
    if desired_prefix == "chr" and not has_chr:
        return "chr" + chrom
    if desired_prefix == "" and has_chr:
        return chrom.replace("chr", "", 1)
    return chrom


def ensure_fasta_index(ref_fasta) -> None:
    """Local FASTA files need a `.fai` built once. Remote URLs need nothing here --
    `samtools` fetches the reference's existing remote `.fai` over HTTP transparently
    whenever a region is requested."""
    if is_remote(ref_fasta):
        return
    ref_fai = Path(str(ref_fasta) + ".fai")
    if not ref_fai.exists():
        print(f"[INFO] Indexing {ref_fasta} (one-time)...")
        _run(["samtools", "faidx", str(ref_fasta)])


def gene_window(gtf: pd.DataFrame, gene: str, chr_prefix: str, window_size: int = WINDOW_SIZE) -> Tuple[str, int, int]:
    """Returns (chrom, start, end) -- the gene's GTF "gene" feature, resized (centered)
    to `window_size`. Same algorithm `build_window_and_predict.py` uses to build the
    windows this repo's `1kG_high_coverage` dataset was originally constructed with."""
    from alphagenome.data import gene_annotation

    interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene).resize(window_size)
    chrom = coerce_chromosome_name(interval.chromosome, chr_prefix)
    return chrom, interval.start, interval.end


def extract_ref_window_fasta(ref_fasta: Path, chrom: str, start: int, end: int, out_path: Path) -> None:
    region = f"{chrom}:{start + 1}-{end}"
    proc = _run(["samtools", "faidx", str(ref_fasta), region])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(proc.stdout)


# ── Alternative: reuse the pre-built 1kG_high_coverage dataset's own window definition ──
#
# Only needed if you also want to import that dataset's precomputed AlphaGenome
# predictions (see `import_precomputed.py`) instead of calling the API -- those
# predictions and their consensus-ready VCFs are keyed to that dataset's own window
# boundaries, which can differ from a freshly-`resize()`d GTF interval by a
# rounding-induced base pair or two. Using this path instead of `gene_window()` above
# keeps everything bit-for-bit consistent with the predictions being imported.

def load_dataset_window(dataset_dir: Path, gene: str, window_size: int = WINDOW_SIZE) -> Tuple[str, int, int]:
    path = dataset_dir / "references" / "windows" / gene / "window_metadata.json"
    with open(path) as f:
        w = json.load(f)
    start = int(w["start"])
    return w["chromosome"], start, start + window_size


def copy_ref_window_fasta(dataset_dir: Path, gene: str, out_path: Path) -> None:
    src = dataset_dir / "references" / "windows" / gene / "ref.window.fa"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, out_path)


# ── Per-gene variant slicing ──

def slice_gene_variants(gene: str, chrom: str, start: int, end: int, vcf_path, sample_ids: list) -> None:
    out_dir = DATA_DIR / "variants"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_vcf = out_dir / f"{gene}.vcf.gz"
    if out_vcf.exists() and Path(str(out_vcf) + ".tbi").exists():
        print(f"[INFO] {gene}: variants VCF already exists, skipping slice")
        return
    if not is_remote(vcf_path) and not Path(vcf_path).exists():
        raise FileNotFoundError(f"VCF not found for {gene} ({chrom}): {vcf_path}")

    region = f"{chrom}:{start + 1}-{end}"
    samples_file = out_dir / f".{gene}.samples.txt"
    samples_file.write_text("\n".join(sample_ids))
    try:
        _run([
            "bcftools", "view",
            "-S", str(samples_file), "--force-samples",
            "-r", region,
            "-e", 'ALT~"<" && ALT!="<DEL>" && ALT!="<NON_REF>"',
            "-Oz", "-o", str(out_vcf),
            str(vcf_path),
        ])
        _run(["bcftools", "index", "-t", str(out_vcf)])
    finally:
        samples_file.unlink(missing_ok=True)
    print(f"[INFO] {gene}: sliced variants VCF ({out_vcf.stat().st_size / 1e6:.1f} MB)")


# ── notebooks/data/ writers ──

def write_experiment_json() -> None:
    path = DATA_DIR / "experiment.json"
    with open(path, "w") as f:
        json.dump(
            {
                "genes": GENES,
                "ontology_terms": ONTOLOGY_TERMS,
                "class_map": CLASS_MAP,
                "window_center_size": WINDOW_CENTER_SIZE,
                "window_size": WINDOW_SIZE,
            },
            f,
            indent=2,
        )
    print(f"[INFO] wrote {path}")


def write_individuals_json(df: pd.DataFrame) -> Dict[str, dict]:
    populations = {pop for pops in CLASS_MAP.values() for pop in pops}
    subset = df[df["Population"].isin(populations)]
    individuals = {
        row["SampleID"]: {
            "population": row["Population"],
            "superpopulation": row["Superpopulation"],
            "sex": SEX_LABELS.get(int(row["Sex"]), "unknown"),
            "family_id": row["FamilyID"],
        }
        for _, row in subset.iterrows()
    }
    path = DATA_DIR / "individuals.json"
    with open(path, "w") as f:
        json.dump(individuals, f, indent=2)
    print(f"[INFO] wrote {path} ({len(individuals)} individuals across {len(populations)} populations)")
    return individuals


def write_genes_json(genes_json: Dict[str, dict]) -> None:
    path = DATA_DIR / "genes.json"
    with open(path, "w") as f:
        json.dump(genes_json, f, indent=2)
    print(f"[INFO] wrote {path}")


def main(
    ref_fasta=DEFAULT_REF_FASTA,
    vcf_pattern: str = DEFAULT_VCF_PATTERN,
    metadata_csv: Path = DEFAULT_METADATA_CSV,
    precomputed_dataset_dir: Path | None = None,
) -> None:
    """`ref_fasta`/`vcf_pattern` may be local paths or http(s)/ftp URLs (the defaults are
    the public 1000 Genomes FTP -- no local genomics data required).

    If `precomputed_dataset_dir` is given (this machine's own pre-built
    `1kG_high_coverage` dataset), gene windows and reference FASTAs are copied from it
    directly instead of re-derived from the GTF -- required for exact alignment if
    you'll also import its precomputed predictions via `import_precomputed.py`."""
    require_tools()
    if not is_remote(ref_fasta) and not Path(ref_fasta).exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    df = load_metadata(metadata_csv)
    print(summarize_class_populations(df, CLASS_MAP).to_string(index=False))

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    write_experiment_json()
    individuals = write_individuals_json(df)
    sample_ids = sorted(individuals)

    if precomputed_dataset_dir is not None:
        gtf = None
        chr_prefix = None
    else:
        ensure_fasta_index(ref_fasta)
        chr_prefix = detect_chr_prefix(ref_fasta)
        print(f"[INFO] Loading GTF from: {GTF_URL}")
        gtf = pd.read_feather(GTF_URL)

    genes_json = {}
    for gene in GENES:
        if precomputed_dataset_dir is not None:
            chrom, start, end = load_dataset_window(precomputed_dataset_dir, gene)
            copy_ref_window_fasta(precomputed_dataset_dir, gene, DATA_DIR / "references" / gene / "ref.window.fa")
        else:
            chrom, start, end = gene_window(gtf, gene, chr_prefix)
            extract_ref_window_fasta(ref_fasta, chrom, start, end, DATA_DIR / "references" / gene / "ref.window.fa")
        print(f"[INFO] {gene}: window={chrom}:{start}-{end}")
        vcf_path = vcf_pattern.format(chrom=chrom)  # keep as str -- Path() mangles URLs
        slice_gene_variants(gene, chrom, start, end, vcf_path, sample_ids)
        genes_json[gene] = {"chrom": chrom, "start": start, "end": end}

    write_genes_json(genes_json)
    print("[DONE] notebooks/data is ready.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_REF_FASTA)
