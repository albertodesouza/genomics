"""plink2-based GWAS baseline over the pigmentation candidate-gene panel -- a classical-statistics
comparison point for the AlphaGenome-based classifiers evaluated in `05_evaluation.ipynb`: does a
standard variant-level association test, run with real QC (missingness, HWE, relatedness,
LD-pruning, multiple-testing correction), find the same genes/regions "significant"?

Important limitation, by design of the whole notebook series (not introduced here): this cohort's
phenotype (`notebooks/data/experiment.json`'s `class_map`) is defined by 1000 Genomes population
membership -- AFR populations = "strong pigmentation", EUR populations = "weak pigmentation". That
means genome-wide ancestry PCs would be almost perfectly collinear with the phenotype itself
(and we don't have genome-wide genotype data locally to compute them anyway -- only the 11
gene-window VCFs under `notebooks/data/variants/`). So this pipeline deliberately does NOT attempt
ancestry-PC adjustment; "significant" hits here reflect AFR/EUR allele-frequency differentiation
within these genes broadly, not fine-mapped individual causality -- the same limitation the
AlphaGenome side already has, from the same population-derived class labels. What this pipeline
*does* correct for, because these are independently justified regardless of the ancestry question:
sample relatedness (pedigree-based, one individual per 1000 Genomes family -- see
`unrelated_sample_ids`), variant/sample QC
(`--maf`/`--geno`/`--mind`/HWE-in-controls), LD structure (pruning for kinship estimation,
clumping of association results so one haplotype block -- e.g. the ~235kb OCA2/HERC2 window
overlap -- isn't reported as dozens of independent "hits"), and multiple testing (Bonferroni over
the actual number of tested variants, since this is a candidate-gene-scoped test, not a genome
scan).

Follows this package's established conventions: functions take an explicit `cache_dir: Path`
(`annotations.py`, `literature_variants.py`), shell out and read results back from disk rather
than parsing stdout (`literature_variants.py`'s `_run_gtex_cli`), and reuse the VCF-to-PLINK
conversion flags (`--double-id --snps-only just-acgt --max-alleles 2`) already established in
`src/genomics/workflows/genomes_analyzer/legacy.py`'s `_convert_vcf_to_bed`.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from . import data as data_utils

PLINK2_INSTALL_HINT = (
    "plink2 not found on PATH. No official upstream Linux ARM64 build exists (cog-genomics.org "
    "only ships Linux x86_64 and macOS ARM64), but Debian/Ubuntu package an aarch64 build. "
    "Install it with:\n"
    "  sudo apt-get update && sudo apt-get install -y plink2"
)


def require_plink2() -> str:
    plink2 = shutil.which("plink2")
    if plink2 is None:
        raise RuntimeError(PLINK2_INSTALL_HINT)
    return plink2


def _run(cmd: List[str]) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
        )
    return proc


def _count_data_lines(path: Path) -> int:
    """Counts non-header lines in a plink2 output file (.pvar/.psam/.prune.in/.snplist/.id --
    comment/header lines all start with `#`, whether or not the file has one)."""
    with open(path) as f:
        return sum(1 for line in f if not line.startswith("#"))


# ── 1. Merge the already-sliced per-gene VCFs into one candidate-panel VCF ──

def merge_gene_vcfs(genes: List[str], cache_dir: Path) -> Path:
    """bcftools-concatenates each gene's already-sliced VCF (`notebooks/data/variants/<gene>.vcf.gz`,
    same files the rest of this notebook series uses) into one VCF covering the whole candidate
    panel. Genes are ordered by (chrom, start) so that any pair of genes whose windows overlap
    (OCA2/HERC2, ~235kb shared) are adjacent inputs -- required for `-d snps` to correctly dedup
    the shared region's variants instead of double-counting them as if independent. A final
    `bcftools sort` guarantees a well-formed, tabix-indexable output regardless of input order."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    gene_rows = data_utils.load_genes()
    ordered_genes = sorted(genes, key=lambda g: (gene_rows[g]["chrom"], gene_rows[g]["start"]))

    vcf_paths = []
    for gene in ordered_genes:
        _, vcf_path = data_utils.gene_paths(gene)
        if not vcf_path.exists():
            raise FileNotFoundError(f"Missing sliced variants VCF for {gene}: {vcf_path}")
        vcf_paths.append(str(vcf_path))

    concat_vcf = cache_dir / "_merged_panel.unsorted.vcf.gz"
    out_vcf = cache_dir / "merged_panel.vcf.gz"
    _run(["bcftools", "concat", "-a", "-d", "snps", "-Oz", "-o", str(concat_vcf), *vcf_paths])
    _run(["bcftools", "sort", "-Oz", "-o", str(out_vcf), str(concat_vcf)])
    _run(["bcftools", "index", "-t", "-f", str(out_vcf)])
    concat_vcf.unlink(missing_ok=True)

    return out_vcf


# ── 2. Phenotype (case/control) and covariate (sex) files ──

def write_phenotype_covariate_files(
    individuals: Dict[str, dict],
    samples_by_class: Dict[str, List[str]],
    cache_dir: Path,
    case_class: str = "strong pigmentation",
    control_class: str = "weak pigmentation",
) -> Tuple[Path, Path]:
    """plink2-format phenotype (`PHENO1`: 2=case/`case_class`, 1=control/`control_class` -- the
    default plink 1/2 case-control coding) and covariate (`SEX`: 1=male, 2=female) files.
    FID==IID throughout, matching `--double-id` used when converting the merged VCF to pgen."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    sex_code = {"male": 1, "female": 2}

    case_ids = samples_by_class[case_class]
    control_ids = samples_by_class[control_class]

    pheno_path = cache_dir / "phenotype.tsv"
    with open(pheno_path, "w") as f:
        f.write("#FID\tIID\tPHENO1\n")
        for sample_id in case_ids:
            f.write(f"{sample_id}\t{sample_id}\t2\n")
        for sample_id in control_ids:
            f.write(f"{sample_id}\t{sample_id}\t1\n")

    covar_path = cache_dir / "covariates.tsv"
    with open(covar_path, "w") as f:
        f.write("#FID\tIID\tSEX\n")
        for sample_id in case_ids + control_ids:
            sex = sex_code.get(individuals[sample_id]["sex"])
            if sex is None:
                raise ValueError(f"{sample_id}: unrecognized sex {individuals[sample_id]['sex']!r}")
            f.write(f"{sample_id}\t{sample_id}\t{sex}\n")

    return pheno_path, covar_path


# ── 3. VCF -> pgen ──

def convert_to_pgen(merged_vcf: Path, cache_dir: Path) -> Path:
    """`--double-id` (FID=IID, matches the phenotype/covariate files above), `--snps-only
    just-acgt --max-alleles 2` (biallelic SNVs only -- same restriction score_variants and the
    rest of this notebook series already apply), `--set-all-var-ids` (guarantees unique,
    position-based IDs for --extract/--clump to key off), `--output-chr chrM` (keeps "chr"-prefixed
    contig names throughout, matching `windows_df`/`genes.json`'s convention)."""
    plink2 = require_plink2()
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = cache_dir / "panel"
    _run([
        plink2, "--vcf", str(merged_vcf), "--double-id",
        "--snps-only", "just-acgt", "--max-alleles", "2",
        "--set-all-var-ids", "@:#:$r:$a", "--output-chr", "chrM",
        "--make-pgen", "--out", str(out_prefix),
    ])
    return out_prefix


def _read_psam_ids(psam_path: Path) -> List[str]:
    with open(psam_path) as f:
        header = f.readline().lstrip("#").strip().split("\t")
        iid_col = header.index("IID")
        return [line.strip().split("\t")[iid_col] for line in f]


def unrelated_sample_ids(individuals: Dict[str, dict], sample_ids: List[str]) -> List[str]:
    """One representative sample per family, using 1000 Genomes' own published pedigree
    (`individuals[id]["family_id"]`) rather than a genotype-based kinship estimate.

    Genotype-based relatedness estimation (e.g. `--king-cutoff`) needs many independent loci
    spread across the genome to separate true pedigree relatedness from ordinary
    population-structure allele-sharing. The only genotype data available here is confined to 11
    highly population-differentiated candidate genes (that's the whole point of this panel --
    genes under strong recent selection, with the most extreme AFR/EUR frequency differences in
    the genome), so a kinship estimate computed from it would conflate "same population" with
    "related individual". Confirmed empirically during development: `--king-cutoff` run on
    LD-pruned markers from a single one of these genes flagged 69% of an unrelated 1000 Genomes
    cohort as "related" -- an implausible relatedness rate, and the signature of exactly this
    problem. Since 1000 Genomes' pedigree is exact ground truth (not an estimate), using it
    directly is both more correct and simpler for this dataset."""
    seen_families = set()
    kept = []
    for sample_id in sorted(sample_ids):
        family_id = individuals[sample_id].get("family_id") or sample_id
        if family_id in seen_families:
            continue
        seen_families.add(family_id)
        kept.append(sample_id)
    return kept


# ── 4. QC funnel: basic QC -> LD-prune -> relatedness -> HWE-in-controls -> final pgen ──

def run_qc_pipeline(
    raw_pgen_prefix: Path,
    pheno_path: Path,
    individuals: Dict[str, dict],
    cache_dir: Path,
    maf: float = 0.01,
    geno: float = 0.1,
    mind: float = 0.1,
    ld_window_kb: str = "50kb",
    ld_r2: float = 0.2,
    hwe_p: float = 1e-6,
) -> Dict[str, object]:
    """Standard GWAS QC funnel. Returns `{"pgen_prefix": <final pgen prefix>, "counts": {step:
    {"samples": n, "variants": n}}}` -- the counts make each filtering step's impact visible
    rather than silent.

    1. `--maf`/`--geno`/`--mind`: basic variant/sample missingness and allele-frequency QC.
    2. `--indep-pairwise`: LD-prune. Not used to restrict the association test itself (step 5) --
       kept as a real, reported QC step: its count gives a less-conservative alternative
       Bonferroni basis (effective independent test count) alongside the raw all-tested-variants
       one used by default (see `bonferroni_threshold`).
    3. `unrelated_sample_ids`: pedigree-based relatedness pruning, one sample per family --
       see that function's docstring for why this replaces genotype-based `--king-cutoff` here.
    4. `--hwe` restricted to controls: plink2's `--hwe` has no automatic case/control awareness
       (unlike some older tools) -- must explicitly `--keep-if PHENO1==control` first.
    5. Intersect unrelated samples x HWE-passing variants into the final analysis-ready pgen.
    """
    plink2 = require_plink2()
    cache_dir.mkdir(parents=True, exist_ok=True)
    counts: Dict[str, Dict[str, int]] = {
        "raw": {
            "samples": _count_data_lines(raw_pgen_prefix.with_suffix(".psam")),
            "variants": _count_data_lines(raw_pgen_prefix.with_suffix(".pvar")),
        }
    }

    qc1_prefix = cache_dir / "qc1_basic"
    _run([
        plink2, "--pfile", str(raw_pgen_prefix),
        "--maf", str(maf), "--geno", str(geno), "--mind", str(mind),
        "--make-pgen", "--out", str(qc1_prefix),
    ])
    counts["after_maf_geno_mind"] = {
        "samples": _count_data_lines(qc1_prefix.with_suffix(".psam")),
        "variants": _count_data_lines(qc1_prefix.with_suffix(".pvar")),
    }

    ld_prefix = cache_dir / "ld"
    _run([
        plink2, "--pfile", str(qc1_prefix),
        "--indep-pairwise", ld_window_kb, str(ld_r2),
        "--out", str(ld_prefix),
    ])
    prune_in = ld_prefix.with_suffix(".prune.in")
    counts["ld_pruned_variant_count"] = {"variants": _count_data_lines(prune_in)}

    unrelated = unrelated_sample_ids(individuals, _read_psam_ids(qc1_prefix.with_suffix(".psam")))
    unrelated_ids = cache_dir / "unrelated.id"
    with open(unrelated_ids, "w") as f:
        f.write("#FID\tIID\n")
        for sample_id in unrelated:
            f.write(f"{sample_id}\t{sample_id}\n")
    counts["unrelated"] = {"samples": len(unrelated)}

    hwe_prefix = cache_dir / "hwe"
    _run([
        plink2, "--pfile", str(qc1_prefix),
        "--pheno", str(pheno_path), "--keep", str(unrelated_ids),
        "--keep-if", "PHENO1==control", "--hwe", str(hwe_p),
        "--write-snplist", "--out", str(hwe_prefix),
    ])
    hwe_snplist = hwe_prefix.with_suffix(".snplist")
    counts["hwe_pass_in_controls"] = {"variants": _count_data_lines(hwe_snplist)}

    final_prefix = cache_dir / "final"
    _run([
        plink2, "--pfile", str(qc1_prefix),
        "--keep", str(unrelated_ids), "--extract", str(hwe_snplist),
        "--make-pgen", "--out", str(final_prefix),
    ])
    counts["final"] = {
        "samples": _count_data_lines(final_prefix.with_suffix(".psam")),
        "variants": _count_data_lines(final_prefix.with_suffix(".pvar")),
    }
    return {"pgen_prefix": final_prefix, "counts": counts}


# ── 5. Association test ──

def _find_glm_output(out_prefix: Path) -> Path:
    candidates = sorted(out_prefix.parent.glob(f"{out_prefix.name}.PHENO1.glm.*"))
    if not candidates:
        raise RuntimeError(f"No plink2 --glm output found for prefix {out_prefix}")
    return candidates[0]


def run_association(qc_pgen_prefix: Path, pheno_path: Path, covar_path: Path, cache_dir: Path) -> pd.DataFrame:
    """Logistic regression of phenotype on allele dosage, covarying for sex only -- no ancestry
    PCs, see this module's docstring for why. `firth-fallback`: standard logistic regression,
    falling back to Firth (bias-reduced) regression only where it fails to converge -- current
    plink2 default/recommended mode. `hide-covar`: keep only the genotype-term row per variant."""
    plink2 = require_plink2()
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = cache_dir / "assoc"
    _run([
        plink2, "--pfile", str(qc_pgen_prefix),
        "--pheno", str(pheno_path), "--covar", str(covar_path), "--covar-name", "SEX",
        "--glm", "firth-fallback", "hide-covar", "cols=+a1freq,+beta,+se",
        "--out", str(out_prefix),
    ])
    result_path = _find_glm_output(out_prefix)
    return pd.read_csv(result_path, sep="\t")


def orient_beta_to_alt(assoc_df: pd.DataFrame) -> pd.DataFrame:
    """Adds a `beta_alt` column: `BETA` re-signed to always represent the ALT allele's effect.

    plink2's `--glm` reports `BETA` for whichever allele it picked as the `A1` test allele
    (usually the minor allele, which is not always the panel's ALT) -- `REF`/`ALT`/`A1` are all
    present in `assoc_df` to check this. Every other variant-scoring method compared against this
    in `05_evaluation.ipynb` (score_variants' raw_score, the 3 trained models' delta-logits) is
    computed for a specific REF->ALT substitution, so an unoriented `BETA` would have an
    arbitrarily flipped sign relative to them for whichever SNPs plink2 happened to code A1=REF."""
    assoc_df = assoc_df.copy()
    assoc_df["beta_alt"] = np.where(assoc_df["A1"] == assoc_df["ALT"], assoc_df["BETA"], -assoc_df["BETA"])
    return assoc_df


def bonferroni_threshold(n_tested: int, alpha: float = 0.05) -> float:
    """Bonferroni significance threshold over the actual number of tested variants -- this is a
    candidate-gene-scoped test (11 genes), not a genome scan, so a genome-wide 5e-8 threshold
    would be the wrong (needlessly conservative) frame."""
    return alpha / n_tested


# ── 6. LD-based clumping (collapse correlated hits within a haplotype block to one lead SNP) ──

def clump_results(
    assoc_df: pd.DataFrame,
    qc_pgen_prefix: Path,
    cache_dir: Path,
    p1: float = 1e-4,
    r2: float = 0.5,
    kb: int = 250,
) -> pd.DataFrame:
    """LD-clumps `assoc_df` (plink2 --glm output) so that one haplotype block -- e.g. the
    ~235kb/45% OCA2/HERC2 window overlap documented in `01_overview_and_vep.ipynb` -- is reported
    as a single independent signal (one lead/"index" SNP plus its tagged/SP2 variants) rather than
    every correlated SNP in the block counted as if independent."""
    plink2 = require_plink2()
    cache_dir.mkdir(parents=True, exist_ok=True)
    assoc_path = cache_dir / "_assoc_for_clump.tsv"
    assoc_df.to_csv(assoc_path, sep="\t", index=False)

    out_prefix = cache_dir / "clump"
    _run([
        plink2, "--pfile", str(qc_pgen_prefix),
        "--clump", str(assoc_path),
        "--clump-p1", str(p1), "--clump-r2", str(r2), "--clump-kb", str(kb),
        "--clump-id-field", "ID", "--clump-p-field", "P",
        "--out", str(out_prefix),
    ])
    candidates = [
        p for p in cache_dir.glob(f"{out_prefix.name}.clump*")
        if p.suffix not in (".log",) and ".log" not in p.name
    ]
    if not candidates:
        raise RuntimeError(f"No plink2 --clump output found for prefix {out_prefix}")
    return pd.read_csv(candidates[0], sep="\t")


# ── 7. Attribute each tested variant to its candidate gene(s) ──

def _normalize_chrom(chrom) -> str:
    """`--output-chr` (used at VCF->pgen conversion time to get "chr"-prefixed contig names) only
    applies to the invocation that sets it -- every subsequent plink2 command (QC, --glm, --clump)
    re-derives its own output text from the internal chromosome code and reverts to plink2's
    default bare-numeric format ("19", not "chr19"). Normalizing here, at comparison time, is
    robust to that regardless of which of the two formats a given plink2 output actually used."""
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def assign_genes(assoc_df: pd.DataFrame, windows_df: pd.DataFrame) -> pd.DataFrame:
    """Adds a `genes` column: the list of candidate genes whose window contains this variant's
    position (a *list*, not a single gene -- OCA2/HERC2's windows overlap by ~235kb, so a variant
    in the shared region legitimately belongs to both, and forcing a single label would silently
    misattribute it)."""
    assoc_df = assoc_df.copy()

    def genes_at(chrom: str, pos: int) -> List[str]:
        chrom = _normalize_chrom(chrom)
        hits = windows_df[
            (windows_df["chrom"] == chrom) & (windows_df["start"] <= pos) & (pos < windows_df["end"])
        ]
        return hits["gene"].tolist()

    assoc_df["genes"] = [genes_at(chrom, pos) for chrom, pos in zip(assoc_df["#CHROM"], assoc_df["POS"])]
    return assoc_df
