"""Build a sample's diploid consensus FASTA for a genomic window at runtime.

Sibling of `genomics.predictors.genotype_based.analysis.individual_consensus`, adapted
to read from this notebook's own bundled data (`notebooks/data/variants/<gene>.vcf.gz`,
already sliced to the gene's window and restricted to the individuals this notebook
cares about) instead of a whole-chromosome VCF under the shared dataset mount -- so the
notebook has no dependency on `genomics.*` or any external dataset path at prediction
time. See `notebooks/utils/prepare_data.py` for how that per-gene VCF is produced.

Results are cached to a caller-supplied directory so repeat calls for the same
sample/gene/haplotype are cheap.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Tuple


def _run(cmd: list) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, proc.stdout, proc.stderr)
    return proc


def _read_fasta_seq_only(path: Path) -> str:
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


def _write_fasta(seq: str, path: Path, header: str) -> None:
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def _adjust_to_target_size(consensus_seq: str, ref_window_seq: str, target_size: int) -> str:
    length = len(consensus_seq)
    if length == target_size:
        return consensus_seq
    if length > target_size:
        return consensus_seq[:target_size]
    needed = target_size - length
    pad_from_ref = ref_window_seq[length:length + needed]
    if len(pad_from_ref) < needed:
        pad_from_ref += "N" * (needed - len(pad_from_ref))
    return consensus_seq + pad_from_ref


def build_haplotype_fasta(
    *,
    sample_id: str,
    gene: str,
    window_width: int,
    ref_window_fa: Path,
    gene_vcf_path: Path,
    haplotype: str,
    cache_dir: Path,
    force: bool = False,
) -> Tuple[str, Path]:
    """Builds (and caches) one haplotype's consensus FASTA for a sample+window.

    Args:
        sample_id: sample ID (must be a column in `gene_vcf_path`).
        gene: gene/target name, used only to namespace the cache directory.
        window_width: expected sequence length (the gene's `end - start`).
        ref_window_fa: path to the (shared, non-individual) reference FASTA for this
            exact window -- `notebooks/data/references/<gene>/ref.window.fa`.
        gene_vcf_path: `notebooks/data/variants/<gene>.vcf.gz` -- already sliced to
            this gene's window and restricted to the relevant individuals.
        haplotype: "H1" or "H2".
        cache_dir: root directory to cache intermediate/output files under.
        force: rebuild even if cached outputs already exist.

    Returns:
        (fixed_sequence, consensus_ready_vcf_path) -- the sequence is exactly
        `window_width` bases long, and the VCF path is suitable for
        `realignment.load_haplotype_variants`.
    """
    if shutil.which("bcftools") is None:
        raise RuntimeError(
            "bcftools not found on PATH. Run scripts/env/start_genomics_universal.sh "
            "or otherwise install bcftools."
        )
    if haplotype not in ("H1", "H2"):
        raise ValueError(f"haplotype must be 'H1' or 'H2', got {haplotype!r}")
    if not ref_window_fa.exists():
        raise FileNotFoundError(f"ref_window_fa not found: {ref_window_fa}")

    case_dir = cache_dir / "individuals" / sample_id / "windows" / gene
    case_dir.mkdir(parents=True, exist_ok=True)

    vcf_cons = case_dir / f"{sample_id}.window.consensus_ready.vcf.gz"
    if force or not vcf_cons.exists():
        vcf_cons.unlink(missing_ok=True)
        Path(str(vcf_cons) + ".tbi").unlink(missing_ok=True)
        _run([
            "bcftools", "view",
            "-s", sample_id,
            "-e", 'ALT~"<" && ALT!="<DEL>" && ALT!="<NON_REF>"',
            "-Oz", "-o", str(vcf_cons),
            str(gene_vcf_path),
        ])
        _run(["bcftools", "index", "-t", str(vcf_cons)])

    fixed_fa = case_dir / f"{sample_id}.{haplotype}.window.fixed.fa"
    if not force and fixed_fa.exists():
        return _read_fasta_seq_only(fixed_fa), vcf_cons

    hap_arg = "1" if haplotype == "H1" else "2"
    proc = _run(["bcftools", "consensus", "-H", hap_arg, "-f", str(ref_window_fa), str(vcf_cons)])
    raw_seq = "".join(line.strip() for line in proc.stdout.splitlines() if not line.startswith(">"))

    ref_seq = _read_fasta_seq_only(ref_window_fa)
    fixed_seq = _adjust_to_target_size(raw_seq, ref_seq, window_width)
    _write_fasta(fixed_seq, fixed_fa, header=f"{sample_id}_{haplotype}_{gene}")

    return fixed_seq, vcf_cons
