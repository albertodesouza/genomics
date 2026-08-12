"""DNA alignment, duplicate marking, CRAM conversion, and coverage."""

import re
import shutil
import subprocess as sp
import zipfile
from pathlib import Path
from typing import Optional

from . import legacy


def align_one_sample(
    r1: Path,
    threads: int,
    read_type: str,
    cleanup,
    expected_r2_name: Optional[str] = None,
    sample_id: Optional[str] = None,
):
    def _bwa_index_ready(prefix: Path) -> bool:
        return all(Path(str(prefix) + ext).exists() for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"))

    def _mem2_index_present(prefix: Path) -> bool:
        return Path(str(prefix) + ".bwt.2bit.64").exists() or any(Path("refs").glob("reference.fa.*.0123"))

    def _fastqc_total_for_filename(target_name: str) -> Optional[int]:
        qc_dir = Path("qc")
        if not qc_dir.exists():
            return None
        for zpath in qc_dir.glob("*.zip"):
            try:
                with zipfile.ZipFile(zpath) as zf:
                    member = next((name for name in zf.namelist() if name.endswith("/fastqc_data.txt")), None)
                    if member is None:
                        continue
                    with zf.open(member) as fh:
                        fname_ok = False
                        total = None
                        for raw in fh:
                            line = raw.decode("utf-8", "ignore").rstrip("\n")
                            if line.startswith("Filename"):
                                parts = line.split("\t", 1)
                                if len(parts) == 2 and parts[1].strip() == target_name:
                                    fname_ok = True
                            elif line.startswith("Total Sequences"):
                                parts = line.split("\t", 1)
                                if len(parts) == 2:
                                    try:
                                        total = int(parts[1].replace(",", "").strip())
                                    except Exception:
                                        total = None
                            if fname_ok and total is not None:
                                return total
            except Exception:
                continue
        return None

    def _estimate_total_reads_for_alignment(r1_path: Path, r2_path: Optional[Path]) -> Optional[int]:
        n1 = _fastqc_total_for_filename(r1_path.name)
        n2 = _fastqc_total_for_filename(r2_path.name) if r2_path is not None else None
        if n1 is None and r2_path is not None:
            n1 = _fastqc_total_for_filename(r1_path.name)
            if n1 is not None:
                return n1 * 2
        if n1 is None:
            return None
        return n1 + (n2 or 0)

    sample = sample_id if sample_id else re.sub(r"(_R?1|\.1)(\.trim)?\.(fastq|fq)\.gz$", "", r1.name)

    out_sorted = Path("bam") / f"{sample}.sorted.bam"
    out_mkdup = Path("bam") / f"{sample}.mkdup.bam"
    out_cram = Path("bam") / f"{sample}.mkdup.cram"

    if out_mkdup.exists() or out_cram.exists():
        legacy.console.print(f"[{sample}] alinhado → [bold]SKIP (cache)[/bold]")
        return

    r2 = None
    candidates = []
    if expected_r2_name:
        candidates.append(r1.parent / expected_r2_name)
    name = r1.name
    for pat, rep in (
        (r"_1\.trim\.fq\.gz$", "_2.trim.fq.gz"),
        (r"_1\.fastq\.gz$", "_2.fastq.gz"),
        (r"_1\.fq\.gz$", "_2.fq.gz"),
        (r"R1\.fastq\.gz$", "R2.fastq.gz"),
        (r"R1\.fq\.gz$", "R2.fq.gz"),
        (r"\.1\.fastq\.gz$", ".2.fastq.gz"),
        (r"\.1\.fq\.gz$", ".2.fq.gz"),
    ):
        if re.search(pat, name):
            candidates.append(r1.parent / re.sub(pat, rep, name))
    for candidate in candidates:
        if candidate.exists():
            r2 = candidate
            break
    if r2 and r2.resolve() == r1.resolve():
        legacy.console.print(f"[yellow][{sample}] Atenção:[/yellow] R2==R1 detectado; prosseguindo como single-end.", style="italic")
        r2 = None

    g = legacy.cfg_global.get("general", {})
    aln_threads = max(1, int(g.get("aln_threads", threads)))
    sort_threads = int(g.get("sort_threads", max(1, aln_threads // 2)))
    sort_threads = max(1, min(sort_threads, aln_threads))
    sort_mem_mb = int(g.get("sort_mem_mb", 512))
    bwa_batch_k = str(int(g.get("bwa_batch_k", 20_000_000)))

    aligner_cfg = g.get("aligner", "bwa-mem2").lower()
    rg_field = f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    ref_prefix = Path("refs/reference.fa")

    aligner_effective = aligner_cfg
    if read_type == "short":
        if aligner_cfg in ("bwa-mem2", "bwa_mem2", "mem2"):
            if not _mem2_index_present(ref_prefix):
                legacy.console.print(
                    f"[orange3][{sample}] Índice do bwa-mem2 ausente → fallback automático para 'bwa mem'. "
                    "Use aligner: bwa + bwa_prebuilt_url no YAML OU gere o índice do mem2 se tiver RAM.[/orange3]"
                )
                aligner_effective = "bwa"
            else:
                aligner_effective = "bwa-mem2"
        else:
            aligner_effective = "bwa"

    legacy.console.print(
        f"[bold cyan][{sample}] Alinhador:[/bold cyan] {('BWA' if aligner_effective == 'bwa' else 'BWA-MEM2')}  "
        f"[bold cyan]aln_threads:[/bold cyan] {aln_threads}  "
        f"[bold cyan]-K:[/bold cyan] {bwa_batch_k}  "
        f"[bold cyan]sort -@:[/bold cyan] {sort_threads}  "
        f"[bold cyan]sort -m:[/bold cyan] {sort_mem_mb}M  "
        f"[bold cyan]R2:[/bold cyan] {('OK' if r2 else 'single-end')}"
    )

    if read_type == "short":
        if aligner_effective == "bwa" and not _bwa_index_ready(ref_prefix):
            raise RuntimeError(
                "Índice BWA ausente/incompleto em refs/reference.fa.{amb,ann,bwt,pac,sa}. "
                "Se estiver usando o pacote do GDC, garanta que os symlinks apontem para os arquivos extraídos."
            )

        if aligner_effective == "bwa":
            base_cmd = ["bwa", "mem", "-R", rg_field, "-t", str(aln_threads), "-K", bwa_batch_k, "refs/reference.fa"]
        else:
            base_cmd = ["bwa-mem2", "mem", "-R", rg_field, "-t", str(aln_threads), "-K", bwa_batch_k, "-Y", "refs/reference.fa"]

        if shutil.which("stdbuf"):
            base_cmd = ["stdbuf", "-oL", "-eL"] + base_cmd

        reads = [str(r1)] + ([str(r2)] if r2 else [])
        total_reads = _estimate_total_reads_for_alignment(r1, r2)
        cmd_list = [
            base_cmd + reads,
            ["samtools", "view", "-u", "-"],
            ["samtools", "sort", "-@", str(sort_threads), "-m", f"{sort_mem_mb}M", "-o", str(out_sorted)],
        ]

        try:
            legacy.run_long_stream_pipeline(
                cmd_list,
                label=f"[{sample}] {('BWA' if aligner_effective == 'bwa' else 'BWA-MEM2')} → sort",
                heartbeat_sec=60,
                progress_spec={"total_reads": total_reads, "print_every": 15},
            )
        except sp.CalledProcessError as exc:
            if aligner_effective == "bwa-mem2":
                legacy.console.print(
                    f"[orange3][{sample}] bwa-mem2 falhou (status {exc.returncode}). "
                    "Tentando fallback com 'bwa mem' e -K menor…[/orange3]"
                )
                small_k = str(max(5_000_000, int(bwa_batch_k) // 2))
                base_cmd_fb = ["bwa", "mem", "-R", rg_field, "-t", str(aln_threads), "-K", small_k, "refs/reference.fa"]
                if shutil.which("stdbuf"):
                    base_cmd_fb = ["stdbuf", "-oL", "-eL"] + base_cmd_fb
                cmd_list_fb = [
                    base_cmd_fb + reads,
                    ["samtools", "view", "-u", "-"],
                    ["samtools", "sort", "-@", str(sort_threads), "-m", f"{sort_mem_mb}M", "-o", str(out_sorted)],
                ]
                legacy.run_long_stream_pipeline(
                    cmd_list_fb,
                    label=f"[{sample}] BWA (fallback) → sort",
                    heartbeat_sec=60,
                    progress_spec={"total_reads": total_reads, "print_every": 15},
                )
            else:
                raise

        legacy.run(["samtools", "index", str(out_sorted)])
        legacy.run(
            [
                "picard",
                "MarkDuplicates",
                "-I",
                str(out_sorted),
                "-O",
                str(out_mkdup),
                "-M",
                str(out_sorted).replace(".sorted.bam", ".mkdup.metrics"),
                "--VALIDATION_STRINGENCY",
                "LENIENT",
            ]
        )
        legacy.run(["samtools", "index", str(out_mkdup)])

        if cleanup.get("remove_sorted_bam", True) and out_sorted.exists():
            out_sorted.unlink(missing_ok=True)
    else:
        cmd_list = [
            ["minimap2", "-ax", "map-ont", "-t", str(aln_threads), "refs/reference.fa", str(r1)],
            ["samtools", "sort", "-@", str(sort_threads), "-m", f"{sort_mem_mb}M", "-o", str(out_sorted)],
        ]
        legacy.run_long_stream_pipeline(cmd_list, label=f"[{sample}] minimap2 → sort", heartbeat_sec=60, progress_spec=None)
        legacy.run(["samtools", "index", str(out_sorted)])


def align_dna_for_all(dna_samples, threads, default_read_type, cleanup, use_ds):
    def _has_mkdup(sid: str) -> bool:
        bam_dir = Path("bam")
        return (bam_dir / f"{sid}.mkdup.bam").exists() or (bam_dir / f"{sid}.mkdup.cram").exists()

    sample_ids = [sample["id"] for sample in dna_samples if "id" in sample]
    unique_ids = sorted(set(sample_ids))
    if unique_ids and all(_has_mkdup(sid) for sid in unique_ids):
        legacy.console.print("Etapa 6 (alinhamento) → [bold]SKIP[/bold] (todas as amostras já têm mkdup)", style="dim")
        outs = sorted(Path("bam").glob("*.mkdup.bam")) + sorted(Path("bam").glob("*.mkdup.cram"))
        if outs:
            legacy.print_meta("Alinhados (mkdup)", outs)
        return

    def _pick_r1_for_run(acc: str) -> Optional[Path]:
        trimmed = sorted(Path("trimmed").glob(f"{acc}*_1.trim.fq.gz"))
        if trimmed:
            return trimmed[0]
        if use_ds:
            ds = sorted(Path("fastq_ds").glob(f"{acc}*_1.fastq.gz")) + sorted(Path("fastq_ds").glob(f"{acc}*_1.fq.gz"))
            if ds:
                return ds[0]
        fq = sorted(Path("fastq").glob(f"{acc}*_1.fastq.gz")) + sorted(Path("fastq").glob(f"{acc}*_1.fq.gz"))
        return fq[0] if fq else None

    for sample in dna_samples:
        read_type = sample.get("read_type", default_read_type)
        sid = sample["id"]
        if _has_mkdup(sid):
            legacy.console.print(f"[{sid}] alinhado → [bold]SKIP (cache mkdup)[/bold]")
            continue

        if sample.get("source") == "sra":
            for acc in sample.get("sra_ids", []):
                r1 = _pick_r1_for_run(acc)
                if not r1:
                    legacy.console.print(f"[yellow][{sid}] Nenhum R1 encontrado para {acc} (trimmed/fastq_ds/fastq).[/yellow]")
                    continue
                align_one_sample(r1, threads, read_type, cleanup, expected_r2_name=None, sample_id=sid)
        else:
            fastq1 = Path(sample.get("fastq1", ""))
            if fastq1:
                r1_trim = Path("trimmed") / fastq1.name.replace(".fastq.gz", ".trim.fq.gz")
                r1 = r1_trim if r1_trim.exists() else (Path("fastq_ds") / fastq1.name if use_ds else Path("fastq") / fastq1.name)
                expected_r2 = Path(sample.get("fastq2", ""))
                align_one_sample(r1, threads, read_type, cleanup, expected_r2_name=(expected_r2.name if expected_r2 else None), sample_id=sid)
            else:
                legacy.console.print(f"[yellow][{sid}] Sem 'fastq1' e sem 'sra_ids'; amostra ignorada.[/yellow]")

    outs = sorted(Path("bam").glob("*.mkdup.bam")) + sorted(Path("bam").glob("*.mkdup.cram"))
    if outs:
        legacy.print_meta("Alinhados (mkdup)", outs)


def to_cram_and_coverage(use_cram: bool, threads: int):
    bam_dir = Path("bam")
    ref_fa = Path("refs/reference.fa")
    if not bam_dir.exists():
        legacy.console.print("[yellow]Diretório 'bam' não existe — nada a fazer na etapa 7.[/yellow]")
        return

    mkdups = sorted(bam_dir.glob("*.mkdup.bam")) + sorted(bam_dir.glob("*.mkdup.cram"))
    if not mkdups:
        legacy.console.print("[yellow]Nenhum *.mkdup.(bam|cram) encontrado — nada a fazer.[/yellow]")
        return

    def _mosdepth_done(prefix: Path) -> bool:
        return prefix.with_suffix(".mosdepth.global.dist.txt").exists() and prefix.with_suffix(".mosdepth.summary.txt").exists()

    for mk in mkdups:
        sid = mk.name.split(".")[0]
        bam = bam_dir / f"{sid}.mkdup.bam"
        cram = bam_dir / f"{sid}.mkdup.cram"
        cov_prefix = bam_dir / sid

        if use_cram:
            if not cram.exists() and bam.exists():
                if not ref_fa.exists():
                    raise RuntimeError("refs/reference.fa ausente — necessário para converter BAM→CRAM.")
                legacy.console.print(f"[{sid}] Convertendo BAM→CRAM…")
                legacy.run(["samtools", "view", "-C", "-T", str(ref_fa), "-@", str(max(1, threads // 2)), "-o", str(cram), str(bam)])
            inp = cram if cram.exists() else (bam if bam.exists() else None)
        else:
            inp = bam if bam.exists() else (cram if cram.exists() else None)

        if not inp:
            legacy.console.print(f"[yellow][{sid}] Nenhum mkdup BAM/CRAM encontrado. Pulando.[/yellow]")
            continue

        if inp.suffix == ".cram":
            crai = Path(str(inp) + ".crai")
            if not crai.exists():
                legacy.run(["samtools", "index", "-@", str(max(1, threads // 2)), str(inp)])
        else:
            bai = Path(str(inp) + ".bai")
            if not bai.exists():
                legacy.run(["samtools", "index", "-@", str(max(1, threads // 2)), str(inp)])

        if _mosdepth_done(cov_prefix):
            legacy.console.print(f"[{sid}] mosdepth → [bold]SKIP[/bold] (cache)")
            continue

        cmd = ["mosdepth", "-t", str(threads)]
        if inp.suffix == ".cram":
            if not ref_fa.exists():
                raise RuntimeError("refs/reference.fa ausente — necessário para mosdepth com CRAM.")
            cmd += ["--fasta", str(ref_fa)]
        cmd += [str(cov_prefix), str(inp)]

        legacy.console.print(f"[{sid}] mosdepth em {inp.name}…")
        legacy.run(cmd)

    outs = sorted(Path("bam").glob("*.mosdepth.summary.txt"))
    if outs:
        legacy.print_meta("Cobertura (mosdepth)", outs)


__all__ = ["align_dna_for_all", "align_one_sample", "to_cram_and_coverage"]
