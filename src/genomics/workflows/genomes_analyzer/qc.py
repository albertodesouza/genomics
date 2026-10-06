"""FASTQ quality control and trimming."""

import re
from pathlib import Path
from typing import Optional

from . import legacy


def _fastqc_base(fq: Path) -> str:
    name = fq.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fa.gz", ".fasta.gz", ".fa", ".fasta"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def fastqc_outputs_exist(fq: Path, outdir: Path = Path("qc")) -> bool:
    base = _fastqc_base(fq)
    zipf = outdir / f"{base}_fastqc.zip"
    html = outdir / f"{base}_fastqc.html"
    return zipf.exists() and html.exists()


def _find_r2_for_r1(r1: Path) -> Optional[Path]:
    name = r1.name
    candidates = [
        re.sub(r"_1\.ds\.fastq\.gz$", "_2.ds.fastq.gz", name, flags=re.I),
        re.sub(r"_1\.fastq\.gz$", "_2.fastq.gz", name, flags=re.I),
        re.sub(r"_1\.fq\.gz$", "_2.fq.gz", name, flags=re.I),
        re.sub(r"_1\.fastq$", "_2.fastq", name, flags=re.I),
        re.sub(r"_1\.fq$", "_2.fq", name, flags=re.I),
        re.sub(r"_1\.ds\.fq\.gz$", "_2.ds.fq.gz", name, flags=re.I),
    ]
    for cand in candidates:
        p = r1.with_name(cand)
        if p.exists():
            return p
    return None


def qc_and_trim(threads, adapters, read_type, use_ds: bool):
    fq_dir = Path("fastq_ds") if use_ds and any(Path("fastq_ds").glob("*.fastq.gz")) else Path("fastq")
    gz = sorted(list(fq_dir.glob("*.fastq.gz")) + list(fq_dir.glob("*.fq.gz")))
    Path("qc").mkdir(exist_ok=True)
    Path("trimmed").mkdir(exist_ok=True)

    to_run = [f for f in gz if not fastqc_outputs_exist(f)]
    if to_run:
        miss = ", ".join(_fastqc_base(f) for f in to_run)
        legacy.console.print(f"[yellow]FastQC faltando para:[/yellow] {miss}", style="dim")
        legacy.run(["fastqc", *map(str, to_run), "-o", "qc"])
    else:
        legacy.console.print("FastQC → [bold]SKIP (todos presentes)[/bold]")

    def _multiqc_uptodate(outdir: Path = Path("qc")) -> bool:
        report = outdir / "multiqc_report.html"
        if not report.exists():
            alt = sorted(outdir.glob("multiqc_report*.html"))
            if not alt:
                return False
            latest = max(alt, key=lambda p: p.stat().st_mtime)
            try:
                report.write_bytes(latest.read_bytes())
            except Exception:
                return False
        fq_outs = []
        for f in gz:
            base = _fastqc_base(f)
            fq_outs += [(outdir / f"{base}_fastqc.zip"), (outdir / f"{base}_fastqc.html")]
        fq_outs = [p for p in fq_outs if p.exists()]
        if not fq_outs:
            return False
        latest_fq = max(p.stat().st_mtime for p in fq_outs)
        return report.stat().st_mtime >= latest_fq

    if not _multiqc_uptodate():
        legacy.run(["multiqc", "qc", "-o", "qc", "--force"])
    else:
        legacy.console.print("MultiQC → [bold]SKIP (atual)[/bold]")

    if read_type != "short":
        legacy.console.print("Leituras longas: pulo trimming por padrão.", style="dim")
        return

    fwd, rev = adapters["fwd"], adapters["rev"]

    for r1 in sorted(
        list(fq_dir.glob("*_1*.fastq.gz"))
        + list(fq_dir.glob("*_1*.fq.gz"))
        + list(fq_dir.glob("*_1*.fastq"))
        + list(fq_dir.glob("*_1*.fq"))
    ):
        core = re.sub(r"_1(?:\.ds)?\.(?:fastq|fq)(?:\.gz)?$", "", r1.name, flags=re.IGNORECASE)
        out1 = Path("trimmed") / f"{core}_1.trim.fq.gz"

        r2 = _find_r2_for_r1(r1)
        if r2:
            out2 = Path("trimmed") / f"{core}_2.trim.fq.gz"
            if out1.exists() and out2.exists():
                legacy.console.print(f"{core}: trimming → [bold]SKIP (cache)[/bold]")
                continue
            legacy.run(
                [
                    "cutadapt",
                    "-j",
                    str(threads),
                    "-q",
                    "20,20",
                    "-m",
                    "30",
                    "-a",
                    fwd,
                    "-A",
                    rev,
                    "-o",
                    str(out1),
                    "-p",
                    str(out2),
                    str(r1),
                    str(r2),
                ]
            )
            legacy.print_meta(f"FASTQs pós-trimming ({core})", [out1, out2])
        else:
            if out1.exists():
                legacy.console.print(f"{core}: trimming (single) → [bold]SKIP (cache)[/bold]")
                continue
            legacy.run(["cutadapt", "-j", str(threads), "-q", "20", "-m", "30", "-a", fwd, "-o", str(out1), str(r1)])
            legacy.print_meta(f"FASTQ pós-trimming ({core})", [out1])

    trimmed = sorted(Path("trimmed").glob("*.fq.gz"))
    to_run_trim = [f for f in trimmed if not fastqc_outputs_exist(f)]
    if to_run_trim:
        legacy.run(["fastqc", *map(str, to_run_trim), "-o", "qc"])
        legacy.run(["multiqc", "qc", "-o", "qc", "--force"])


__all__ = ["fastqc_outputs_exist", "qc_and_trim"]
