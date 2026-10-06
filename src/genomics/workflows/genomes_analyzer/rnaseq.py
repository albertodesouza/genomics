"""Optional RNA-seq alignment and transcript assembly steps."""

import re
from pathlib import Path
from typing import Optional

from . import legacy
from .fastq import stage_fastqs_from_local, stage_fastqs_from_sra


def rnaseq_pipeline(rna_samples, threads, assembly_name):
    if not rna_samples:
        return
    Path("rnaseq").mkdir(exist_ok=True)

    for sample in rna_samples:
        if sample["source"] == "sra":
            stage_fastqs_from_sra(sample["sra_ids"])
        else:
            stage_fastqs_from_local(sample["fastq1"], sample.get("fastq2"))

    def _find_r2_for_r1_file(r1: Path) -> Optional[Path]:
        name = r1.name
        for pat, rep in (
            (r"_1\.fastq\.gz$", "_2.fastq.gz"),
            (r"_1\.fq\.gz$", "_2.fq.gz"),
            (r"R1\.fastq\.gz$", "R2.fastq.gz"),
            (r"R1\.fq\.gz$", "R2.fq.gz"),
            (r"\.1\.fastq\.gz$", ".2.fastq.gz"),
            (r"\.1\.fq\.gz$", ".2.fq.gz"),
        ):
            if re.search(pat, name, flags=re.I):
                cand = r1.parent / re.sub(pat, rep, name, flags=re.I)
                if cand.exists():
                    return cand
        return None

    r1_list: list[Path] = []
    for sample in rna_samples:
        if sample["source"] == "sra":
            for acc in sample["sra_ids"]:
                path = Path("fastq") / f"{acc}_1.fastq.gz"
                if path.exists():
                    r1_list.append(path)
        else:
            fq1 = Path("fastq") / Path(sample["fastq1"]).name
            if fq1.exists():
                r1_list.append(fq1)

    for r1 in sorted(set(r1_list)):
        base = re.sub(r"_1\.fastq\.gz$", "", r1.name, flags=re.I)
        r2 = _find_r2_for_r1_file(r1)
        gtf_out = Path("rnaseq") / f"{base}.transcripts.gtf"
        bam_out = Path("rnaseq") / f"{base}.rnaseq.bam"

        need = not (gtf_out.exists() and legacy._is_newer(gtf_out, r1, r2 if r2 else r1, Path("refs/genes.gtf")))
        if not need:
            legacy.console.print(f"[RNA-seq:{base}] GTF → [bold]SKIP (atual)[/bold]")
            continue

        if r2 and r2.exists():
            cmd_list = [
                ["hisat2", "-p", str(threads), "-x", f"refs/{assembly_name}", "-1", str(r1), "-2", str(r2)],
                ["samtools", "sort", "-@", str(threads), "-o", str(bam_out)],
            ]
        else:
            cmd_list = [
                ["hisat2", "-p", str(threads), "-x", f"refs/{assembly_name}", "-U", str(r1)],
                ["samtools", "sort", "-@", str(threads), "-o", str(bam_out)],
            ]
        legacy.run_long_stream_pipeline(cmd_list, label=f"[RNA-seq:{base}] HISAT2 → sort")
        legacy.run(["samtools", "index", str(bam_out)])
        legacy.run(["stringtie", str(bam_out), "-G", "refs/genes.gtf", "-o", str(gtf_out), "-p", str(threads)])

    tgts = list(Path("rnaseq").glob("*.transcripts.gtf"))
    cmp_prefix = Path("rnaseq") / "cmp"
    need_cmp = tgts and (not (cmp_prefix.with_suffix(".stats").exists()))
    if tgts and need_cmp:
        legacy.run(["gffcompare", "-r", "refs/genes.gtf", "-o", str(cmp_prefix), *map(str, tgts)])

__all__ = ["rnaseq_pipeline"]
