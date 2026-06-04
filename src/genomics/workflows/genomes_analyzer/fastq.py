"""FASTQ acquisition, SRA/ENA conversion, compression, and downsampling."""

import subprocess as sp
import time
import re
import shlex
import shutil
from pathlib import Path
from typing import List, Optional

from rich.panel import Panel
from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn, TimeRemainingColumn, TransferSpeedColumn

from . import legacy


def _sra_expected_size_bytes(acc: str) -> Optional[int]:
    try:
        proc = sp.run(["vdb-dump", acc, "--info"], capture_output=True, text=True, check=True)
        match = re.search(r"size\s*:\s*([0-9,]+)", proc.stdout, re.IGNORECASE)
        if match:
            return int(match.group(1).replace(",", ""))
    except Exception:
        pass
    return None


def prefetch_with_progress(acc: str, outdir: Path, interval_sec: float = 2.0, stall_warn_min: int = 10):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    Path("logs").mkdir(exist_ok=True)

    acc_dir = outdir / acc
    acc_file = outdir / f"{acc}.sra"

    def measure_bytes() -> int:
        if acc_dir.exists() and acc_dir.is_dir():
            return legacy._du_bytes(acc_dir)
        if acc_file.exists() and acc_file.is_file():
            return acc_file.stat().st_size
        matches = list(outdir.glob(f"{acc}*"))
        if matches:
            return sum(legacy._du_bytes(match) if match.is_dir() else match.stat().st_size for match in matches)
        return 0

    expected = _sra_expected_size_bytes(acc)
    expected_str = f"{expected:,} B" if expected else "desconhecido"
    legacy.console.print(f"[bold]prefetch[/bold] {acc} → {outdir}  (tamanho esperado: {expected_str})")

    cmd = ["prefetch", "--output-directory", str(outdir), "--max-size", "u", acc]
    legacy.console.print(f"[bold]>[/bold] {' '.join(cmd)}", style="dim")

    log_path = Path("logs") / f"prefetch_{acc}.log"
    with open(log_path, "w") as log_fh:
        proc = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, text=True)

        with Progress(
            TextColumn("[bold blue]{task.fields[acc]}[/]"),
            BarColumn(),
            TextColumn("{task.percentage:>5.1f}%"),
            TextColumn("• {task.description}"),
            TransferSpeedColumn(),
            TextColumn("• {task.completed:,.0f}/{task.total:,.0f} B"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=legacy.console,
            transient=False,
        ) as progress:
            total = expected if expected else 1
            task = progress.add_task("conectando...", total=total, acc=acc)
            last_bytes = 0
            last_change = time.time()

            while True:
                rc = proc.poll()
                bytes_now = measure_bytes()

                if not expected and bytes_now > total:
                    total = max(bytes_now * 2, total * 2)
                    progress.update(task, total=total)

                progress.update(task, completed=bytes_now, description="baixando")

                if bytes_now > last_bytes:
                    last_bytes = bytes_now
                    last_change = time.time()
                elif (time.time() - last_change) > stall_warn_min * 60 and rc is None:
                    legacy.console.print(
                        f"[yellow]Aviso:[/yellow] prefetch {acc} sem progresso por ~{stall_warn_min} min. "
                        f"Veja logs em {log_path}",
                        style="italic",
                    )
                    last_change = time.time()

                if rc is not None:
                    bytes_now = measure_bytes()
                    progress.update(task, completed=min(bytes_now, total), description="finalizando")
                    if rc != 0:
                        try:
                            tail = sp.run(
                                ["conda", "run", "-n", "genomics", "bash", "-lc", f"tail -n 50 {shlex.quote(str(log_path))}"],
                                capture_output=True,
                                text=True,
                                check=True,
                            ).stdout
                        except Exception:
                            tail = ""
                        raise sp.CalledProcessError(rc, cmd, output=f"(veja {log_path})\n{tail}")
                    break

                time.sleep(interval_sec)

    sra = None
    if acc_dir.exists():
        sra = next(acc_dir.glob(f"{acc}.sra"), None) or next(acc_dir.rglob("*.sra"), None)
    if not sra and acc_file.exists():
        sra = acc_file
    vdb = None
    if acc_dir.exists():
        vdb = next(acc_dir.glob(f"{acc}.vdbcache"), None) or next(acc_dir.rglob("*.vdbcache"), None)

    legacy.print_meta(f"Prefetch concluído ({acc})", [path for path in [sra, vdb] if path])
    legacy.console.rule(f"[green]Prefetch {acc} concluído[/green]")


def _find_local_sra(acc: str, raw_dir: Path) -> Optional[Path]:
    raw_dir = Path(raw_dir)
    candidates = [raw_dir / acc / f"{acc}.sra", raw_dir / f"{acc}.sra"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    hits = list(raw_dir.rglob(f"{acc}.sra"))
    return hits[0] if hits else None


def ena_get_fastq_urls(acc: str):
    fields = "fastq_http,fastq_ftp,fastq_md5"
    cmd = [
        "conda",
        "run",
        "-n",
        "genomics",
        "bash",
        "-lc",
        f"curl -fsSL 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields={fields}&format=tsv&download=false' | tail -n +2",
    ]
    result = sp.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        return [], []

    cols = result.stdout.strip().split("\t")
    http_s = cols[0] if len(cols) > 0 else ""
    ftp_s = cols[1] if len(cols) > 1 else ""
    md5_s = cols[2] if len(cols) > 2 else ""

    urls = []
    if http_s:
        urls.extend([url.strip() for url in http_s.split(";") if url.strip()])
    if not urls and ftp_s:
        for url in ftp_s.split(";"):
            url = url.strip()
            if not url:
                continue
            if url.startswith("ftp://"):
                url = "https://" + url[len("ftp://") :]
            urls.append(url)

    urls = [url for url in urls if url.endswith("_1.fastq.gz") or url.endswith("_2.fastq.gz")]
    md5s = [md5.strip() for md5 in md5_s.split(";")] if md5_s else []
    return urls, md5s


def ena_fetch_fastqs(acc: str, outdir: Path, threads: int = 8) -> bool:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    urls, md5s = ena_get_fastq_urls(acc)
    if not urls:
        legacy.console.print(f"[yellow]ENA não retornou URLs de FASTQ para {acc}.[/yellow]")
        return False

    ok_any = False
    for i, url in enumerate(urls):
        fname = url.strip().split("/")[-1]
        dst = outdir / fname

        if dst.exists():
            gz_ok = sp.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"gzip -t {shlex.quote(str(dst))} >/dev/null 2>&1"]).returncode == 0
            if gz_ok:
                legacy.console.print(f"{fname}: → [bold]SKIP (cache OK)[/bold]")
                ok_any = True
                continue
            legacy.console.print(f"{fname}: corrompido, refazendo…", style="yellow")
            dst.unlink(missing_ok=True)

        legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"wget -c -O {shlex.quote(str(dst))} {shlex.quote(url)}"])

        if sp.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"gzip -t {shlex.quote(str(dst))} >/dev/null 2>&1"]).returncode != 0:
            legacy.console.print(f"[red]{fname}: gzip inválido[/red]")
            continue

        if i < len(md5s) and md5s[i]:
            md5 = sp.run(
                ["conda", "run", "-n", "genomics", "bash", "-lc", f"md5sum {shlex.quote(str(dst))} | cut -d' ' -f1"],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
            if md5 != md5s[i]:
                legacy.console.print(f"[red]{fname}: MD5 não confere ({md5} != {md5s[i]})[/red]")
                continue

        ok_any = True

    outs = sorted(outdir.glob(f"{acc}_*.fastq.gz"))
    if outs:
        legacy.print_meta(f"FASTQs baixados do ENA ({acc})", outs)
    have_pair = any((outdir / f"{acc}_1.fastq.gz").exists() and (outdir / f"{acc}_2.fastq.gz").exists())
    return ok_any or have_pair


def fasterq_with_progress(
    source: str,
    acc: str,
    outdir: Path,
    threads: int = 8,
    tmp_dir: Path = Path("tmp"),
    interval_sec: float = 2.0,
    stall_warn_min: int = 10,
    stall_fail_min: int = 45,
) -> bool:
    outdir = Path(outdir).resolve()
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    Path("logs").mkdir(exist_ok=True)

    r1 = outdir / f"{acc}_1.fastq"
    r2 = outdir / f"{acc}_2.fastq"
    log_path = Path("logs") / f"fasterq_{acc}.log"

    cmd = ["fasterq-dump", "--split-files", "-e", str(threads), "-t", str(tmp_dir), str(source), "-O", str(outdir)]
    legacy.console.print(f"[bold]>[/bold] {' '.join(cmd)}", style="dim")

    cancel_on_convert_stall = bool(legacy.cfg_global["general"].get("cancel_on_convert_stall", False))
    min_delta = 32 * 1024 * 1024

    with open(log_path, "w") as log_fh:
        proc = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, text=True)
        io_last = legacy._proc_io(proc.pid)
        io_last_ts = time.time()

        with Progress(
            TextColumn("[bold blue]fasterq[/] " + acc),
            BarColumn(),
            TextColumn("{task.percentage:>5.1f}%"),
            TextColumn("• {task.description}"),
            TransferSpeedColumn(),
            TextColumn("• {task.completed:,.0f}/{task.total:,.0f} B"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=legacy.console,
            transient=False,
        ) as progress:
            total = 1
            phase = "preparando"
            task = progress.add_task(phase, total=total)
            last_done = 0
            last_tmp = 0
            last_change = time.time()
            next_warn_at = stall_warn_min

            while True:
                rc = proc.poll()
                have_r1 = r1.exists()
                b1 = r1.stat().st_size if have_r1 else 0
                have_r2 = r2.exists()
                b2 = r2.stat().st_size if have_r2 else 0
                done = b1 + b2
                tmp_bytes = legacy._du_bytes(tmp_dir)

                if have_r1 or have_r2:
                    if phase != "convertendo":
                        phase = "convertendo"
                        progress.update(task, description=phase)
                    if done > total:
                        total = int(done * 1.10)
                        progress.update(task, total=total)
                    progress.update(task, completed=done)
                else:
                    if tmp_bytes > total:
                        total = int(tmp_bytes * 1.10)
                        progress.update(task, total=total)
                    progress.update(task, completed=tmp_bytes, description="preparando (lendo .sra)")

                grew = False
                if done > last_done + min_delta:
                    last_done = done
                    grew = True
                if tmp_bytes > last_tmp + min_delta:
                    last_tmp = tmp_bytes
                    grew = True
                if grew:
                    last_change = time.time()

                io_now = legacy._proc_io(proc.pid)
                if io_now and io_last:
                    d_read = io_now[0] - io_last[0]
                    d_write = io_now[1] - io_last[1]
                    if d_read > 0 or d_write > 0:
                        last_change = time.time()
                if io_now:
                    if time.time() - io_last_ts >= 60:
                        legacy.console.print(
                            f"[dim]fasterq {acc} I/O heartbeat: +{legacy.sizeof_fmt(max(0, io_now[0] - (io_last[0] if io_last else 0)))} lidos, "
                            f"+{legacy.sizeof_fmt(max(0, io_now[1] - (io_last[1] if io_last else 0)))} escritos[/dim]"
                        )
                        io_last_ts = time.time()
                    io_last = io_now

                idle_min = (time.time() - last_change) / 60.0
                if idle_min >= next_warn_at and rc is None:
                    legacy.console.print(f"[yellow]Aviso:[/yellow] fasterq {acc} sem progresso há ~{int(idle_min)} min. Log: {log_path}", style="italic")
                    next_warn_at += 5

                if idle_min >= stall_fail_min and rc is None:
                    if phase == "preparando" or cancel_on_convert_stall:
                        legacy.console.print(f"[orange3]Sem progresso há ≥{stall_fail_min} min — cancelando fasterq ({acc}).[/orange3]")
                        proc.terminate()
                        try:
                            proc.wait(timeout=30)
                        except sp.TimeoutExpired:
                            proc.kill()
                        return False
                    next_warn_at += 5

                if rc is not None:
                    if (r1.exists() or r2.exists()) and rc == 0:
                        legacy.print_meta(f"FASTQs gerados ({acc})", [path for path in [r1, r2] if path.exists()])
                        legacy.console.rule(f"[green]fasterq-dump {acc} concluído[/green]")
                        return True
                    if rc != 0:
                        try:
                            tail = sp.run(
                                ["conda", "run", "-n", "genomics", "bash", "-lc", f"tail -n 50 {shlex.quote(str(log_path))}"],
                                capture_output=True,
                                text=True,
                                check=True,
                            ).stdout
                        except Exception:
                            tail = ""
                        legacy.console.print(f"[red]fasterq falhou ({acc})[/red]\n{tail}")
                        return False

                time.sleep(interval_sec)


def compress_fastqs_with_progress(acc: str, outdir: Path, threads: int = 8):
    outdir = Path(outdir).resolve()
    r1 = outdir / f"{acc}_1.fastq"
    r2 = outdir / f"{acc}_2.fastq"
    to_do = [path for path in [r1, r2] if path.exists() and not (Path(str(path) + ".gz").exists())]
    if not to_do:
        legacy.console.print(f"{acc}: compressão → [bold]SKIP (cache)[/bold]")
        return

    use_pigz = shutil.which("pigz") is not None
    comp = ["pigz", "-p", str(threads)] if use_pigz else ["gzip", "-9"]

    with Progress(
        TextColumn("[bold blue]compressão[/]"),
        BarColumn(),
        TextColumn("{task.percentage:>5.1f}%"),
        TextColumn("• {task.description}"),
        TransferSpeedColumn(),
        TextColumn("• {task.completed:,.0f}/{task.total:,.0f} B"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=legacy.console,
        transient=False,
    ) as progress:
        tasks = {}
        totals = {}
        for fq in to_do:
            totals[fq] = fq.stat().st_size
            tasks[fq] = progress.add_task(fq.name, total=max(1, totals[fq]), start=False)

        procs = {}
        for fq in to_do:
            procs[fq] = sp.Popen(comp + [str(fq)], stdout=sp.DEVNULL, stderr=sp.DEVNULL, text=False)
            progress.start_task(tasks[fq])

        while procs:
            finished = []
            for fq, proc in procs.items():
                gz = Path(str(fq) + ".gz")
                gz_bytes = gz.stat().st_size if gz.exists() else 0
                progress.update(tasks[fq], completed=min(gz_bytes, totals[fq]))
                if proc.poll() is not None:
                    finished.append(fq)
                    progress.update(tasks[fq], completed=totals[fq], description="ok")
            for fq in finished:
                procs.pop(fq, None)
            time.sleep(0.5)

    outs = []
    for fq in [r1, r2]:
        gz = Path(str(fq) + ".gz")
        if gz.exists():
            outs.append(gz)
    if outs:
        legacy.print_meta(f"FASTQs comprimidos ({acc})", outs)


def stage_fastqs_from_sra(sra_ids: List[str]):
    Path("raw").mkdir(exist_ok=True)
    Path("fastq").mkdir(exist_ok=True)

    threads = int(legacy.cfg_global.get("general", {}).get("threads", 8))
    tmp_root = Path(legacy.cfg_global.get("general", {}).get("temp_dir", "tmp")).expanduser().resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)

    for acc in sra_ids:
        tmp_dir = tmp_root / "sra" / acc
        tmp_dir.mkdir(parents=True, exist_ok=True)

        r1_gz = Path("fastq") / f"{acc}_1.fastq.gz"
        r2_gz = Path("fastq") / f"{acc}_2.fastq.gz"
        if r1_gz.exists() or r2_gz.exists():
            legacy.console.print(f"{acc}: FASTQ(.gz) → [bold]SKIP (cache)[/bold]")
            legacy.print_meta(f"FASTQs ({acc})", [r1_gz, r2_gz])
            continue

        tmp_root = Path(legacy.cfg_global.get("general", {}).get("temp_dir", "tmp")).expanduser().resolve()
        tmp_dir = tmp_root / "sra" / acc
        tmp_dir.mkdir(parents=True, exist_ok=True)

        threads = int(legacy.cfg_global.get("general", {}).get("threads", 8))

        if legacy.cfg_global["general"].get("prefer_ena_fastq", False):
            legacy.console.print(f"[dim]{acc}: prefer_ena_fastq=true → tentando ENA primeiro[/dim]")
            if ena_fetch_fastqs(acc, outdir=Path("fastq").resolve(), threads=threads):
                legacy.console.rule(f"[green]ENA HTTP ({acc}) concluído[/green]")
            else:
                legacy.console.print(f"[yellow]{acc}: ENA não disponível/sem URLs — caindo para SRA (prefetch+fasterq).[/yellow]")

        if not r1_gz.exists() and not r2_gz.exists():
            sra_target = next(Path("raw").glob(f"**/{acc}.sra"), None)
            if sra_target and sra_target.exists():
                legacy.console.print(f"{acc}: .sra → [bold]SKIP (cache)[/bold]")
            else:
                retries = int(legacy.cfg_global["general"].get("prefetch_retries", 2))
                backoff = 30
                last_err = None
                for attempt in range(retries + 1):
                    try:
                        prefetch_with_progress(acc, Path("raw").resolve())
                        last_err = None
                        break
                    except sp.CalledProcessError as exc:
                        last_err = exc
                        if attempt < retries:
                            wait = backoff * (attempt + 1)
                            legacy.console.print(
                                f"[orange3]{acc}: prefetch falhou (tentativa {attempt+1}/{retries}). "
                                f"Repetindo em {wait}s…[/orange3]"
                            )
                            time.sleep(wait)
                        else:
                            legacy.console.print(f"[red]{acc}: prefetch falhou após {retries+1} tentativas.[/red]")
                if last_err is not None and not legacy.cfg_global["general"].get("ena_fallback", False):
                    raise last_err

            if not r1_gz.exists() and not r2_gz.exists():
                sra_path = legacy._find_local_sra(acc, Path("raw"))
                if sra_path:
                    ok = fasterq_with_progress(
                        str(sra_path),
                        acc=acc,
                        outdir=Path("fastq").resolve(),
                        threads=threads,
                        tmp_dir=tmp_dir,
                        stall_warn_min=legacy.cfg_global["general"].get("stall_warn_min", 10),
                        stall_fail_min=legacy.cfg_global["general"].get("stall_fail_min", 45),
                    )
                else:
                    ok = fasterq_with_progress(
                        acc,
                        acc=acc,
                        outdir=Path("fastq").resolve(),
                        threads=threads,
                        tmp_dir=tmp_dir,
                        stall_warn_min=legacy.cfg_global["general"].get("stall_warn_min", 10),
                        stall_fail_min=legacy.cfg_global["general"].get("stall_fail_min", 45),
                    )

                if not ok:
                    if legacy.cfg_global["general"].get("ena_fallback", False):
                        if ena_fetch_fastqs(acc, outdir=Path("fastq").resolve(), threads=threads):
                            legacy.console.rule(f"[green]ENA HTTP ({acc}) concluído[/green]")
                        else:
                            raise RuntimeError(f"Não foi possível obter FASTQs para {acc} via fasterq nem via ENA.")
                    else:
                        raise RuntimeError(
                            f"fasterq {acc} não avançou e ena_fallback:false. "
                            f"Veja logs em logs/fasterq_{acc}.log e tente novamente."
                        )

        compress_fastqs_with_progress(acc, outdir=Path("fastq").resolve(), threads=threads)
        legacy.print_meta(f"FASTQs ({acc})", [r1_gz, r2_gz])


def stage_fastqs_from_local(fq1, fq2=None):
    Path("fastq").mkdir(exist_ok=True)
    dst1 = Path("fastq") / Path(fq1).name
    if not dst1.exists():
        legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"ln -s {Path(fq1).resolve()} {dst1} || cp {Path(fq1).resolve()} {dst1}"])
    else:
        legacy.console.print(f"{dst1.name}: → [bold]SKIP (cache)[/bold]")
    if fq2:
        dst2 = Path("fastq") / Path(fq2).name
        if not dst2.exists():
            legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"ln -s {Path(fq2).resolve()} {dst2} || cp {Path(fq2).resolve()} {dst2}"])
        else:
            legacy.console.print(f"{dst2.name}: → [bold]SKIP (cache)[/bold]")
    legacy.print_meta("FASTQs de entrada (locais)", [dst1] + ([dst2] if fq2 else []))


def downsample_fastqs(fraction: float, seed: int):
    if fraction <= 0 or fraction >= 1:
        return
    legacy.console.print(Panel.fit(f"Downsample FASTQ com seqtk (fração={fraction}, seed={seed})", border_style="cyan"))
    Path("fastq_ds").mkdir(exist_ok=True)

    paired_r1 = sorted(list(Path("fastq").glob("*_1.fastq.gz")))
    for r1 in paired_r1:
        r2 = Path(str(r1).replace("_1.fastq.gz", "_2.fastq.gz"))
        base = r1.name.replace("_1.fastq.gz", "")
        if r2.exists():
            out1 = Path("fastq_ds") / f"{base}_1.ds.fastq.gz"
            out2 = Path("fastq_ds") / f"{base}_2.ds.fastq.gz"
            if out1.exists() and out2.exists():
                legacy.console.print(f"{base}: downsample → [bold]SKIP (cache)[/bold]")
                continue
            legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"seqtk sample -s{seed} {r1} {fraction} | gzip > {out1}"])
            legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"seqtk sample -s{seed} {r2} {fraction} | gzip > {out2}"])
            legacy.print_meta(f"FASTQs downsample ({base})", [out1, out2])
        else:
            out1 = Path("fastq_ds") / f"{base}.ds.fastq.gz"
            if out1.exists():
                legacy.console.print(f"{base}: downsample (single) → [bold]SKIP (cache)[/bold]")
                continue
            legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"seqtk sample -s{seed} {r1} {fraction} | gzip > {out1}"])
            legacy.print_meta(f"FASTQ downsample ({base})", [out1])

__all__ = [
    "compress_fastqs_with_progress",
    "downsample_fastqs",
    "ena_fetch_fastqs",
    "ena_get_fastq_urls",
    "fasterq_with_progress",
    "prefetch_with_progress",
    "stage_fastqs_from_local",
    "stage_fastqs_from_sra",
]
