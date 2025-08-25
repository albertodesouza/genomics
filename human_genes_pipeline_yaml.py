#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, sys, time, yaml, shutil, subprocess as sp
from pathlib import Path
from datetime import datetime
from typing import List, Optional

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

import re, shlex
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn, TransferSpeedColumn

console = Console(highlight=False)

# cfg global para funções auxiliares
cfg_global = None

# ====================== Utilidades ======================

def sizeof_fmt(num, suffix="B"):
    for unit in ["","K","M","G","T","P"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}E{suffix}"

def file_meta(p: Path):
    if not p.exists():
        return {"exists": False}
    st = p.stat()
    return {"exists": True, "size_bytes": st.st_size,
            "size": sizeof_fmt(st.st_size),
            "mtime": datetime.fromtimestamp(st.st_mtime).isoformat(" ", "seconds"),
            "path": str(p.resolve())}

def print_meta(title: str, paths: List[Path]):
    tbl = Table(title=title, header_style="bold cyan")
    tbl.add_column("Arquivo"); tbl.add_column("Tamanho", justify="right")
    tbl.add_column("Modificado em"); tbl.add_column("Caminho")
    for p in paths:
        m = file_meta(p)
        if m["exists"]:
            tbl.add_row(p.name, m["size"], m["mtime"], m["path"])
        else:
            tbl.add_row(p.name, "—", "—", "(não existe)")
    console.print(tbl)

def run(cmd: List[str], cwd: Optional[str]=None, env=None, check=True):
    console.print(f"[bold]>[/bold] {' '.join(cmd)}", style="dim")
    sp.run(cmd, cwd=cwd, env=env, check=check)

def run_long_stream_pipeline(cmd_list: List[List[str]], label: str, heartbeat_sec: int = 60):
    """
    Executa pipeline encadeado (Popen) com heartbeats.
    Ex.: [["bwa-mem2","mem",...], ["samtools","view","-bS","-"], ["samtools","sort","-o","..."]]
    """
    console.print(Panel.fit(Text(label, style="bold yellow"), border_style="yellow"))
    console.print("[bold]>[/bold] " + "  |  ".join(" ".join(c) for c in cmd_list), style="dim")

    procs = []
    start = time.time()

    for i, cmd in enumerate(cmd_list):
        if i == 0:
            p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, text=True, bufsize=1, universal_newlines=True)
        else:
            p = sp.Popen(cmd, stdin=procs[-1].stdout, stdout=sp.PIPE if i < len(cmd_list)-1 else None)
        procs.append(p)

    last = start
    if procs and procs[0].stdout:
        for line in procs[0].stdout:
            if line.strip() and any(k in line for k in ["WARNING","ERROR","ERR","progress","reads","mapped","trimmed","Time","ETA","%","Mbp","Gbp"]):
                console.print(line.rstrip(), style="dim")
            now = time.time()
            if now - last > heartbeat_sec:
                elapsed = int(now - start)
                console.print(f"[heartbeat] {label}... {elapsed//60}m{elapsed%60:02d}s", style="magenta")
                last = now

    rc = 0
    for p in procs[::-1]:
        p.wait()
        rc = rc or p.returncode
    if rc != 0:
        raise sp.CalledProcessError(rc, cmd_list[0])
    elapsed = int(time.time() - start)
    console.print(f":check_mark_button: [bold green]{label} concluído[/bold green] em {elapsed//60}m{elapsed%60:02d}s.")

def ensure_dirs():
    for d in ["refs","raw","fastq","fastq_ds","qc","trimmed","bam","vcf","vep","genes","rnaseq","logs"]:
        Path(d).mkdir(parents=True, exist_ok=True)

def disk_space_report(base_dir: Path):
    total, used, free = shutil.disk_usage(base_dir)
    console.print(Panel.fit(
        f"[bold]Espaço em disco[/bold]\nTotal: {sizeof_fmt(total)}  Usado: {sizeof_fmt(used)}  Livre: {sizeof_fmt(free)}\nLocal: {base_dir}",
        border_style="blue"
    ))
    return total, used, free

def _du_bytes(path: Path) -> int:
    """Tamanho total em bytes (recursivo) de 'path' usando du -sb (rápido e robusto)."""
    try:
        out = sp.run(
            ["bash","-lc", f"du -sb {shlex.quote(str(path))} 2>/dev/null | cut -f1"],
            capture_output=True, text=True, check=True
        ).stdout.strip()
        return int(out) if out else 0
    except Exception:
        return 0

def _sra_expected_size_bytes(acc: str) -> Optional[int]:
    """Obtém o tamanho esperado do run via 'vdb-dump <acc> --info' (campo 'size : ...')."""
    try:
        p = sp.run(["vdb-dump", acc, "--info"], capture_output=True, text=True, check=True)
        m = re.search(r"size\s*:\s*([0-9,]+)", p.stdout, re.IGNORECASE)
        if m:
            return int(m.group(1).replace(",", ""))
    except Exception:
        pass
    return None

def prefetch_with_progress(acc: str, outdir: Path, interval_sec: float = 2.0, stall_warn_min: int = 10):
    """
    Baixa um accession SRA com 'prefetch' mostrando progresso.
    Evita PIPE (que pode bloquear) e grava logs em logs/prefetch_<acc>.log.
    """
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    Path("logs").mkdir(exist_ok=True)

    # Onde o prefetch pode escrever:
    #  - outdir/<acc>/... (padrão)
    #  - ou outdir/<acc>.sra diretamente (depende da config)
    acc_dir = outdir / acc
    acc_file = outdir / f"{acc}.sra"

    def measure_bytes() -> int:
        if acc_dir.exists() and acc_dir.is_dir():
            return _du_bytes(acc_dir)
        if acc_file.exists() and acc_file.is_file():
            return acc_file.stat().st_size
        # fallback: mede o que houver com prefixo do acc
        matches = list(outdir.glob(f"{acc}*"))
        if matches:
            return sum(_du_bytes(m) if m.is_dir() else m.stat().st_size for m in matches)
        return 0

    expected = _sra_expected_size_bytes(acc)  # pode ser None
    expected_str = f"{expected:,} B" if expected else "desconhecido"
    console.print(f"[bold]prefetch[/bold] {acc} → {outdir}  (tamanho esperado: {expected_str})")

    cmd = ["prefetch", "--output-directory", str(outdir), "--max-size", "u", acc]
    console.print(f"[bold]>[/bold] {' '.join(cmd)}", style="dim")

    # Redireciona stdout/err para arquivo (evita bloquear PIPE)
    log_path = Path("logs")/f"prefetch_{acc}.log"
    with open(log_path, "w") as log_fh:
        proc = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, text=True)

        from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn, TransferSpeedColumn
        with Progress(
            TextColumn("[bold blue]{task.fields[acc]}[/]"),
            BarColumn(),
            TextColumn("{task.percentage:>5.1f}%"),
            TextColumn("• {task.description}"),
            TransferSpeedColumn(),
            TextColumn("• {task.completed:,.0f}/{task.total:,.0f} B"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=False,
        ) as progress:
            total = expected if expected else 1
            task = progress.add_task("conectando...", total=total, acc=acc)
            last_bytes = 0
            last_change = time.time()

            while True:
                rc = proc.poll()
                bytes_now = measure_bytes()

                # barra "elástica" quando não se sabe o total
                if not expected and bytes_now > total:
                    total = max(bytes_now * 2, total * 2)
                    progress.update(task, total=total)

                # atualização
                progress.update(task, completed=bytes_now, description="baixando")

                # detector de estagnação (sem crescimento por X minutos)
                if bytes_now > last_bytes:
                    last_bytes = bytes_now
                    last_change = time.time()
                elif (time.time() - last_change) > stall_warn_min * 60 and rc is None:
                    console.print(f"[yellow]Aviso:[/yellow] prefetch {acc} sem progresso por ~{stall_warn_min} min. "
                                  f"Veja logs em {log_path}", style="italic")
                    last_change = time.time()  # evita spam do aviso

                if rc is not None:
                    # atualização final
                    bytes_now = measure_bytes()
                    progress.update(task, completed=min(bytes_now, total), description="finalizando")
                    if rc != 0:
                        # Mostra o fim do log pra ajudar no debug
                        try:
                            tail = sp.run(["bash","-lc", f"tail -n 50 {shlex.quote(str(log_path))}"],
                                          capture_output=True, text=True, check=True).stdout
                        except Exception:
                            tail = ""
                        raise sp.CalledProcessError(rc, cmd, output=f"(veja {log_path})\n{tail}")
                    break

                time.sleep(interval_sec)

    # Meta dos artefatos baixados
    sra = None
    if acc_dir.exists():
        sra = next(acc_dir.glob(f"{acc}.sra"), None) or next(acc_dir.rglob("*.sra"), None)
    if not sra and acc_file.exists():
        sra = acc_file
    vdb = None
    if acc_dir.exists():
        vdb = next(acc_dir.glob(f"{acc}.vdbcache"), None) or next(acc_dir.rglob("*.vdbcache"), None)

    print_meta(f"Prefetch concluído ({acc})", [p for p in [sra, vdb] if p])
    console.rule(f"[green]Prefetch {acc} concluído[/green]")

def _find_local_sra(acc: str, raw_dir: Path) -> Optional[Path]:
    raw_dir = Path(raw_dir)
    candidates = [
        raw_dir / acc / f"{acc}.sra",
        raw_dir / f"{acc}.sra",
    ]
    for c in candidates:
        if c.exists():
            return c
    # fallback: procura em subpastas
    hits = list(raw_dir.rglob(f"{acc}.sra"))
    return hits[0] if hits else None

def ena_get_fastq_urls(acc: str):
    """
    Consulta a API do ENA e retorna (urls, md5s).
    Tenta fastq_http; se vazio, usa fastq_ftp (convertendo ftp:// -> https://).
    Filtra apenas *_1.fastq.gz e *_2.fastq.gz.
    """
    fields = "fastq_http,fastq_ftp,fastq_md5"
    cmd = ["bash","-lc",
           f"curl -fsSL 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields={fields}&format=tsv&download=false' | tail -n +2"]
    r = sp.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        return [], []
    cols = r.stdout.strip().split("\t")
    # colunas esperadas: 0=http, 1=ftp, 2=md5
    http_s = cols[0] if len(cols) > 0 else ""
    ftp_s  = cols[1] if len(cols) > 1 else ""
    md5_s  = cols[2] if len(cols) > 2 else ""

    urls = []
    if http_s:
        urls.extend([u.strip() for u in http_s.split(";") if u.strip()])
    if not urls and ftp_s:
        for u in ftp_s.split(";"):
            u = u.strip()
            if not u:
                continue
            # usa HTTPS em vez de FTP (o host do ENA aceita)
            if u.startswith("ftp://"):
                u = "https://" + u[len("ftp://"):]
            urls.append(u)

    # só *_1.fastq.gz / *_2.fastq.gz
    urls = [u for u in urls if u.endswith("_1.fastq.gz") or u.endswith("_2.fastq.gz")]

    md5s = [m.strip() for m in md5_s.split(";")] if md5_s else []
    return urls, md5s

def ena_fetch_fastqs(acc: str, outdir: Path, threads: int = 8) -> bool:
    """
    Baixa *_1.fastq.gz e *_2.fastq.gz do ENA com resume (-c), valida gzip e MD5 (quando fornecido).
    Retorna True se ambos os arquivos válidos estiverem presentes ao final.
    """
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    urls, md5s = ena_get_fastq_urls(acc)
    if not urls:
        console.print(f"[yellow]ENA não retornou URLs de FASTQ para {acc}.[/yellow]")
        return False

    ok_any = False
    for i, url in enumerate(urls):
        fname = url.strip().split("/")[-1]
        dst = outdir / fname

        # Se já existe, valide antes de pular
        if dst.exists():
            gz_ok = sp.run(["bash","-lc", f"gzip -t {shlex.quote(str(dst))} >/dev/null 2>&1"]).returncode == 0
            if gz_ok:
                console.print(f"{fname}: → [bold]SKIP (cache OK)[/bold]")
                ok_any = True
                continue
            else:
                console.print(f"{fname}: corrompido, refazendo…", style="yellow")
                dst.unlink(missing_ok=True)

        run(["bash","-lc", f"wget -c -O {shlex.quote(str(dst))} {shlex.quote(url)}"])

        # Valida gzip
        if sp.run(["bash","-lc", f"gzip -t {shlex.quote(str(dst))} >/dev/null 2>&1"]).returncode != 0:
            console.print(f"[red]{fname}: gzip inválido[/red]")
            continue

        # Valida MD5 se fornecido pela API
        if i < len(md5s) and md5s[i]:
            md5 = sp.run(["bash","-lc", f"md5sum {shlex.quote(str(dst))} | cut -d' ' -f1"],
                         capture_output=True, text=True, check=True).stdout.strip()
            if md5 != md5s[i]:
                console.print(f"[red]{fname}: MD5 não confere ({md5} != {md5s[i]})[/red]")
                continue

        ok_any = True

    # imprime meta do que foi baixado
    outs = sorted(outdir.glob(f"{acc}_*.fastq.gz"))
    if outs:
        print_meta(f"FASTQs baixados do ENA ({acc})", outs)
    # Retorna True se temos pelo menos um par válido
    have_pair = any((outdir/f"{acc}_1.fastq.gz").exists() and (outdir/f"{acc}_2.fastq.gz").exists())
    return ok_any or have_pair

def _proc_io(pid: int):
    """Retorna (read_bytes, write_bytes) do /proc/<pid>/io (Linux)."""
    try:
        with open(f"/proc/{pid}/io") as fh:
            m = {}
            for line in fh:
                k, v = line.split(":")
                m[k.strip()] = int(v.strip())
        return m.get("read_bytes", 0), m.get("write_bytes", 0)
    except Exception:
        return None

def fasterq_with_progress(source: str, acc: str, outdir: Path, threads: int = 8, tmp_dir: Path = Path("tmp"),
                          interval_sec: float = 2.0, stall_warn_min: int = 10, stall_fail_min: int = 45) -> bool:
    """
    Executa fasterq-dump (source = accession ou caminho .sra) com barra em 2 fases.
    Retorna True se converteu; False se falhou/estagnou tempo demais (para permitir fallback).
    Detecção de progresso considera crescimento dos FASTQs e do tmp_dir.
    """
    outdir = Path(outdir).resolve()
    tmp_dir = Path(tmp_dir); tmp_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    Path("logs").mkdir(exist_ok=True)

    r1 = outdir / f"{acc}_1.fastq"
    r2 = outdir / f"{acc}_2.fastq"
    log_path = Path("logs")/f"fasterq_{acc}.log"

    cmd = ["fasterq-dump", "--split-files", "-e", str(threads), "-t", str(tmp_dir), str(source), "-O", str(outdir)]
    console.print(f"[bold]>[/bold] {' '.join(cmd)}", style="dim")

    cancel_on_convert_stall = bool(cfg_global["general"].get("cancel_on_convert_stall", False))
    min_delta = 32 * 1024 * 1024  # 32MB: variação mínima para considerar "progresso"

    with open(log_path, "w") as log_fh:
        proc = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, text=True)
        io_last = _proc_io(proc.pid)
        io_last_ts = time.time()

        with Progress(
            TextColumn("[bold blue]fasterq[/] "+acc),
            BarColumn(),
            TextColumn("{task.percentage:>5.1f}%"),
            TextColumn("• {task.description}"),
            TransferSpeedColumn(),
            TextColumn("• {task.completed:,.0f}/{task.total:,.0f} B"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=False,
        ) as progress:
            total = 1
            phase = "preparando"
            task = progress.add_task(phase, total=total)

            last_done = 0
            last_tmp  = 0
            last_change = time.time()
            next_warn_at = stall_warn_min  # minutos

            while True:
                rc = proc.poll()

                # bytes atuais dos FASTQs (se já existirem)
                have_r1 = r1.exists(); b1 = r1.stat().st_size if have_r1 else 0
                have_r2 = r2.exists(); b2 = r2.stat().st_size if have_r2 else 0
                done = b1 + b2

                # bytes atuais no tmp_dir
                tmp_bytes = _du_bytes(tmp_dir)

                # Seleciona descrição e barra
                if have_r1 or have_r2:
                    if phase != "convertendo":
                        phase = "convertendo"
                        progress.update(task, description=phase)
                    if done > total:
                        total = int(done * 1.10)  # elástico com +10%
                        progress.update(task, total=total)
                    progress.update(task, completed=done)
                else:
                    if tmp_bytes > total:
                        total = int(tmp_bytes * 1.10)
                        progress.update(task, total=total)
                    progress.update(task, completed=tmp_bytes, description="preparando (lendo .sra)")

                # Houve progresso real?
                grew = False
                if done > last_done + min_delta:
                    last_done = done; grew = True
                if tmp_bytes > last_tmp + min_delta:
                    last_tmp = tmp_bytes; grew = True
                if grew:
                    last_change = time.time()

                # Estagnação: sem crescimento de FASTQ nem tmp_dir
                # Heartbeat de I/O do processo: mesmo sem crescer arquivos, conta como progresso
                io_now = _proc_io(proc.pid)
                if io_now and io_last:
                    d_read  = io_now[0] - io_last[0]
                    d_write = io_now[1] - io_last[1]
                    if d_read > 0 or d_write > 0:
                        grew = True
                        last_change = time.time()
                if io_now:
                    # imprime um pulso de vida a cada ~60s
                    if time.time() - io_last_ts >= 60:
                        console.print(
                            f"[dim]fasterq {acc} I/O heartbeat: +{sizeof_fmt(max(0,io_now[0] - (io_last[0] if io_last else 0)))} lidos, "
                            f"+{sizeof_fmt(max(0,io_now[1] - (io_last[1] if io_last else 0)))} escritos[/dim]"
                        )
                        io_last_ts = time.time()
                    io_last = io_now

                idle_min = (time.time() - last_change)/60.0
                if idle_min >= next_warn_at and rc is None:
                    console.print(
                        f"[yellow]Aviso:[/yellow] fasterq {acc} sem progresso há ~{int(idle_min)} min. Log: {log_path}",
                        style="italic"
                    )
                    next_warn_at += 5  # evita spam

                if idle_min >= stall_fail_min and rc is None:
                    # Cancelar apenas se ainda estiver "preparando" OU se o usuário permitir cancelar na conversão
                    if phase == "preparando" or cancel_on_convert_stall:
                        console.print(f"[orange3]Sem progresso há ≥{stall_fail_min} min — cancelando fasterq ({acc}).[/orange3]")
                        proc.terminate()
                        try:
                            proc.wait(timeout=30)
                        except sp.TimeoutExpired:
                            proc.kill()
                        return False
                    else:
                        # Na fase convertendo, não cancelar: só estender o próximo aviso
                        next_warn_at += 5

                if rc is not None:
                    # fim do processo
                    if (r1.exists() or r2.exists()) and rc == 0:
                        print_meta(f"FASTQs gerados ({acc})", [p for p in [r1, r2] if p.exists()])
                        console.rule(f"[green]fasterq-dump {acc} concluído[/green]")
                        return True
                    if rc != 0:
                        try:
                            tail = sp.run(["bash","-lc", f"tail -n 50 {shlex.quote(str(log_path))}"],
                                          capture_output=True, text=True, check=True).stdout
                        except Exception:
                            tail = ""
                        console.print(f"[red]fasterq falhou ({acc})[/red]\n{tail}")
                        return False

                time.sleep(interval_sec)

def compress_fastqs_with_progress(acc: str, outdir: Path, threads: int = 8):
    """
    Comprime ERRxxxx_1.fastq/_2.fastq em .gz mostrando progresso (tamanho do .gz
    em relação ao .fastq original; não é exato, mas dá boa noção).
    Usa pigz se disponível, senão gzip.
    """
    outdir = Path(outdir).resolve()
    r1 = outdir / f"{acc}_1.fastq"
    r2 = outdir / f"{acc}_2.fastq"
    to_do = [p for p in [r1, r2] if p.exists() and not (Path(str(p) + ".gz").exists())]
    if not to_do:
        console.print(f"{acc}: compressão → [bold]SKIP (cache)[/bold]")
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
        console=console,
        transient=False,
    ) as progress:
        tasks = {}
        totals = {}
        for fq in to_do:
            totals[fq] = fq.stat().st_size
            tasks[fq] = progress.add_task(fq.name, total=max(1, totals[fq]), start=False)

        procs = {}
        for fq in to_do:
            gz = Path(str(fq) + ".gz")
            # inicia processo
            procs[fq] = sp.Popen(comp + [str(fq)], stdout=sp.DEVNULL, stderr=sp.DEVNULL, text=False)
            progress.start_task(tasks[fq])

        # monitora até todos acabarem
        while procs:
            finished = []
            for fq, p in procs.items():
                gz = Path(str(fq) + ".gz")
                gz_bytes = gz.stat().st_size if gz.exists() else 0
                # limita a barra ao tamanho do input (estimativa conservadora)
                progress.update(tasks[fq], completed=min(gz_bytes, totals[fq]))
                if p.poll() is not None:
                    finished.append(fq)
                    progress.update(tasks[fq], completed=totals[fq], description="ok")
            for fq in finished:
                procs.pop(fq, None)
            time.sleep(0.5)

    # imprime metadados finais
    outs = []
    for fq in [r1, r2]:
        gz = Path(str(fq) + ".gz")
        if gz.exists():
            outs.append(gz)
    if outs:
        print_meta(f"FASTQs comprimidos ({acc})", outs)


# ================== Normalização do YAML ==================

def normalize_config_schema(cfg_in: dict) -> dict:
    """
    Aceita:
      (A) antigo: {'general': {...}, 'dna_samples': [...], 'rna_samples': [...]}
      (B) novo  : {'project': {...}, 'storage': {...}, 'reference': {...},
                   'download': {...}, 'execution': {...}, 'size_control': {...},
                   'samples': [...], 'params': {...}}
    Retorna esquema (A).
    """
    if "general" in cfg_in:
        return cfg_in

    project   = cfg_in.get("project", {})
    storage   = cfg_in.get("storage", {})
    ref_top   = cfg_in.get("reference", {})
    ref_proj  = project.get("reference", {}) if isinstance(project, dict) else {}
    ref       = {**ref_proj, **ref_top}
    execv     = cfg_in.get("execution", {})
    download  = cfg_in.get("download", {})
    download_tool = (download.get("tool") or "sra_toolkit").lower()
    stall_warn_min = int(execv.get("stall_warn_min", 10))
    stall_fail_min = int(execv.get("stall_fail_min", 45))
    ena_fallback   = bool(execv.get("ena_fallback", False))
    sizec     = cfg_in.get("size_control", {})
    samples   = cfg_in.get("samples", [])
    params    = cfg_in.get("params", {})
    temp_dir = storage.get("temp_dir", "tmp")  # padrão: subpasta 'tmp' no projeto

    base_dir = storage.get("base_dir", ".")
    assembly = ref.get("name", "GRCh38")
    ref_fa_url = ref.get("fasta_url") or \
                 "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
    gtf_url    = ref.get("gtf_url") or \
                 "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"

    threads = int(execv.get("threads") or download.get("threads") or params.get("bwa_mem2_threads") or 16)
    mem_gb  = int(params.get("mem_gb", 64))

    ds_cfg = (sizec or {}).get("downsample", {})
    downsample_frac = float(ds_cfg.get("fraction", 0.0)) if ds_cfg.get("enabled", False) else 0.0
    downsample_seed = int(ds_cfg.get("seed", 123))

    keep_inter = bool(storage.get("keep_intermediates", False))
    use_cram = True
    cleanup = {
        "remove_sorted_bam": True,
        "remove_bam_after_cram": (not keep_inter)
    }

    adapters = {"fwd": "AGATCGGAAGAGC", "rev": "AGATCGGAAGAGC"}

    aligner = (params.get("aligner") or execv.get("aligner") or "bwa-mem2").lower()
    bwa_prebuilt = ref.get("bwa_index_url") or cfg_in.get("bwa_prebuilt_url") or ""
    limit_to_canonical = bool(cfg_in.get("limit_to_canonical", False))

    tool = (download.get("tool") or "sra_toolkit").lower()
    if tool != "sra_toolkit":
        console.print(f"[orange3]Aviso:[/orange3] download.tool='{tool}' ainda não é suportado; usando sra-tools.", style="italic")

    dna_samples = []
    for s in samples:
        runs = s.get("runs", [])
        if not runs:
            continue
        dna_samples.append({
            "id": s.get("sample_id") or s.get("id") or runs[0],
            "source": "sra",
            "sra_ids": runs,
            "read_type": "short"
        })

    general = {
        "base_dir": base_dir,
        "assembly_name": assembly,
        "ref_fa_url": ref_fa_url,
        "gtf_url": gtf_url,
        "threads": threads,
        "mem_gb": mem_gb,
        "default_read_type": "short",
        "adapters": adapters,
        "call_variants": True,
        "annotate_vars": True,
        "rnaseq": False,
        "force_refs": False,
        "force_indexes": False,
        "vep_cache_dir": "",
        "space_guard_gb_min": 100,
        "use_cram": use_cram,
        "cleanup": cleanup,
        "downsample_frac": downsample_frac,
        "downsample_seed": downsample_seed,
        "aligner": aligner,
        "bwa_prebuilt_url": bwa_prebuilt,
        "limit_to_canonical": limit_to_canonical,
        "download_tool": download_tool,
        "stall_warn_min": stall_warn_min,
        "stall_fail_min": stall_fail_min,
        "ena_fallback": ena_fallback,
        "temp_dir": temp_dir,
    }

    return {"general": general, "dna_samples": dna_samples, "rna_samples": []}

# ================ Estimativa de espaço =================

def estimate_outputs_and_warn(fastqs: List[Path], base_dir: Path, guard_gb: int):
    free_b = shutil.disk_usage(base_dir).free
    in_bytes = sum((file_meta(f)["size_bytes"] for f in fastqs if f.exists()))
    trimmed_est = int(in_bytes * 0.8)
    bam_est     = int(in_bytes * 1.2)
    vcf_est     = int(2e9)
    needed_bam_path = bam_est + trimmed_est + vcf_est
    msg = (f"Entrada FASTQ (comprimido): ~{sizeof_fmt(in_bytes)}\n"
           f"Pós-trimming: ~{sizeof_fmt(trimmed_est)}; BAM: ~{sizeof_fmt(bam_est)}; VCF: ~{sizeof_fmt(vcf_est)}\n"
           f"Necessário (sem CRAM/limpeza): ~{sizeof_fmt(needed_bam_path)} | Livre: {sizeof_fmt(free_b)}")
    style = "yellow" if (free_b < needed_bam_path or free_b < guard_gb*1024**3) else "green"
    console.print(Panel.fit(msg, title="Estimativa de espaço", border_style=style))
    if style == "yellow":
        console.print("[orange3]Sugestões:[/orange3] use [bold]use_cram: true[/bold], [bold]cleanup.remove_sorted_bam: true[/bold] e/ou [bold]downsample_frac[/bold].", style="italic")

# ================== Referências / índices ==================

def _download_and_place(url: str, dst_plain: Path):
    """
    Baixa e coloca o arquivo de referência.
    - Se for tar.gz (mesmo sem extensão): extrai e aponta dst_plain -> maior .fa/.fasta encontrado.
    - Se for gz simples (fa.gz): descompacta com gzip -dc.
    - Caso contrário: move para dst_plain.
    """
    tmp = dst_plain.with_suffix(dst_plain.suffix + ".tmp")
    run(["wget", "-O", str(tmp), url])

    # 1) Se for tar.gz (mesmo sem extensão), 'tar -tzf' retorna 0
    is_targz = sp.run(
        ["bash", "-lc", f"tar -tzf {tmp} >/dev/null 2>&1"],
        stdout=sp.DEVNULL, stderr=sp.DEVNULL
    ).returncode == 0

    if is_targz:
        run(["bash","-lc", f"mkdir -p refs/_extracted && tar -xzf {tmp} -C refs/_extracted"])
        fa_candidates = list(Path("refs/_extracted").rglob("*.fa")) + list(Path("refs/_extracted").rglob("*.fasta"))
        if not fa_candidates:
            raise RuntimeError("Nenhum FASTA (.fa/.fasta) encontrado dentro do tar.")
        fa = max(fa_candidates, key=lambda p: p.stat().st_size)
        if dst_plain.exists():
            dst_plain.unlink()
        dst_plain.symlink_to(fa.resolve())
        tmp.unlink(missing_ok=True)
        return

    # 2) Se for gzip simples (fa.gz), descompacta ignorando sufixo
    is_gzip = sp.run(
        ["bash", "-lc", f"file -b --mime-type {tmp} | grep -qi 'gzip'"],
        stdout=sp.DEVNULL, stderr=sp.DEVNULL
    ).returncode == 0
    if is_gzip:
        run(["bash","-lc", f"gzip -dc {tmp} > {dst_plain}"])
        tmp.unlink(missing_ok=True)
        return

    # 3) Caso contrário, mover como está (já é .fa descompactado)
    run(["bash","-lc", f"mv {tmp} {dst_plain}"])

def download_refs(ref_fa_url, gtf_url, force=False):
    Path("refs").mkdir(exist_ok=True)
    fa  = Path("refs/reference.fa")
    gtf = Path("refs/genes.gtf")

    if fa.exists() and not force:
        console.print(":floppy_disk: Referência já presente → [bold]SKIP[/bold]")
    else:
        _download_and_place(ref_fa_url, fa)

    if gtf.exists() and not force:
        console.print(":floppy_disk: Anotação (GTF) já presente → [bold]SKIP[/bold]")
    else:
        _download_and_place(gtf_url, gtf)

    if not Path("refs/reference.fa.fai").exists() or force:
        run(["samtools","faidx","refs/reference.fa"]) 
    if not Path("refs/reference.dict").exists() or force:
        run(["gatk","CreateSequenceDictionary","-R","refs/reference.fa","-O","refs/reference.dict"])

    print_meta("Referências", [fa, gtf])

def _bwa_index_ready(prefix: Path) -> bool:
    """Verifica se os 5 arquivos do índice BWA existem (e se symlinks apontam para arquivos reais)."""
    exts = [".amb",".ann",".bwt",".pac",".sa"]
    for ext in exts:
        p = Path(str(prefix) + ext)
        if not p.exists():
            return False
        if p.is_symlink() and not p.resolve().exists():
            return False
    return True

def _uuid_basename_from_url(url: str, default: str) -> str:
    """Se a URL terminar com um UUID do GDC, usa-o como nome de arquivo; senão usa 'default'."""
    m = re.search(r'([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$', url, re.I)
    return (m.group(1) + ".tar.gz") if m else default

def install_prebuilt_bwa_index(bwa_tar_url: str, expect_md5: str = "", force: bool = False):
    """
    Instala o índice BWA pré-construído do GDC com cache:
      - Se os symlinks refs/reference.fa.{amb,ann,bwt,pac,sa} já existirem → SKIP
      - Se o tar.gz já existir → SKIP download, apenas extrai (se necessário)
      - Caso contrário → baixa com 'wget -c' (resume), extrai e cria symlinks
    """
    if not bwa_tar_url:
        return

    prefix = Path("refs/reference.fa")
    if _bwa_index_ready(prefix) and not force:
        console.print("Índice BWA já presente → [bold]SKIP download[/bold]")
        return

    index_dir = Path("refs/_bwa"); index_dir.mkdir(parents=True, exist_ok=True)
    tar_name = _uuid_basename_from_url(bwa_tar_url, "index.tar.gz")
    tar_path = index_dir / tar_name

    # Se já temos o tar, não baixa de novo
    if tar_path.exists() and not force:
        console.print(f"{tar_path.name} já presente → [bold]SKIP download[/bold]")
    else:
        # -c = resume; -O define nome estável no cache
        run(["bash","-lc", f"wget -c -O {shlex.quote(str(tar_path))} {shlex.quote(bwa_tar_url)}"])

    # (Opcional) checagem de MD5 se você passar expect_md5
    if expect_md5:
        md5 = sp.run(["bash","-lc", f"md5sum {shlex.quote(str(tar_path))} | cut -d' ' -f1"],
                     capture_output=True, text=True, check=True).stdout.strip()
        if md5 != expect_md5:
            raise RuntimeError(f"MD5 não confere para {tar_path.name}: {md5} (esp. {expect_md5})")

    # Extrai apenas se ainda não houver arquivos de índice extraídos
    have_any = list(index_dir.rglob("*.amb")) or list(index_dir.rglob("*.bwt"))
    if not have_any or force:
        run(["bash","-lc", f"tar -xzf {shlex.quote(str(tar_path))} -C {shlex.quote(str(index_dir))}"])

    ambs = list(index_dir.rglob("*.amb"))
    if not ambs:
        raise RuntimeError("Índice BWA não encontrado após extração.")
    pref = ambs[0].with_suffix("")  # prefixo real do índice

    # Cria/atualiza symlinks com o prefixo padrão refs/reference.fa.*
    for ext in [".amb",".ann",".bwt",".pac",".sa"]:
        src = pref.with_suffix(ext)
        if src.exists():
            dst = Path(str(prefix) + ext)
            if dst.exists():
                dst.unlink()
            dst.symlink_to(src.resolve())

    if _bwa_index_ready(prefix):
        console.print("Índice BWA instalado e linkado em [bold]refs/reference.fa.*[/bold]")

def limit_reference_to_canonical_if_enabled():
    g = cfg_global["general"]
    if not g.get("limit_to_canonical", False):
        return
    if not Path("refs/reference.fa.fai").exists():
        run(["samtools","faidx","refs/reference.fa"])
    names = [l.split("\t")[0] for l in open("refs/reference.fa.fai")]
    if not names:
        console.print("[orange3]Aviso:[/orange3] .fai vazio, mantendo referência completa.")
        return
    chr_prefix = names[0].startswith("chr")
    canon = [f"{'chr' if chr_prefix else ''}{i}" for i in list(range(1,23))+["X","Y"]]
    canon += [f"{'chr' if chr_prefix else ''}M", f"{'chr' if chr_prefix else ''}MT"]
    canon = [c for c in canon if c in names]
    if not canon:
        console.print("[orange3]Aviso:[/orange3] Não foi possível detectar contigs canônicos.")
        return
    run(["bash","-lc", f"samtools faidx refs/reference.fa {' '.join(canon)} > refs/reference.canonical.fa"])
    Path("refs/reference.fa").unlink(missing_ok=True)
    Path("refs/reference.fa").symlink_to(Path("refs/reference.canonical.fa").resolve())
    run(["samtools","faidx","refs/reference.fa"]) 
    run(["gatk","CreateSequenceDictionary","-R","refs/reference.fa","-O","refs/reference.dict"])
    console.print("[green]Referência reduzida aos cromossomos canônicos.[/green]")


def build_indexes(default_read_type, assembly_name, need_rna_index, threads, force=False):
    g = cfg_global["general"]
    if g.get("bwa_prebuilt_url"):
        install_prebuilt_bwa_index(g["bwa_prebuilt_url"])
        console.print("Índice BWA pré-instalado → [bold]SKIP build[/bold]")
    elif default_read_type=="short" and g.get("aligner","bwa-mem2")=="bwa-mem2":
        try:
            if force or not Path("refs/reference.fa.bwt.2bit.64").exists():
                run(["bwa-mem2","index","refs/reference.fa"]) 
            else:
                console.print("Índice BWA-MEM2 → [bold]SKIP[/bold]", style="dim")
        except sp.CalledProcessError:
            console.print("[red]Falha ao indexar com bwa-mem2 (provável falta de RAM).[/red]")
            console.print("Sugestão: defina [bold]aligner: bwa[/bold] + [bold]bwa_prebuilt_url[/bold] no YAML ou ative [bold]limit_to_canonical: true[/bold].")
            raise
    else:
        console.print("Sem necessidade de indexar (usar BWA com índice pré-instalado ou long reads).", style="dim")

    if need_rna_index and (force or not Path(f"refs/{assembly_name}.1.ht2").exists()):
        run(["hisat2-build","-p",str(threads),"refs/reference.fa",f"refs/{assembly_name}"])
    elif need_rna_index:
        console.print("Índice HISAT2 → [bold]SKIP[/bold]", style="dim")

# ================== Entrada SRA / FASTQ ===================

def stage_fastqs_from_sra(sra_ids: List[str]):
    Path("raw").mkdir(exist_ok=True)
    Path("fastq").mkdir(exist_ok=True)

    # threads e pastas de temporários (em disco rápido, vindo do YAML)
    threads = int(cfg_global.get("general", {}).get("threads", 8))
    tmp_root = Path(cfg_global.get("general", {}).get("temp_dir", "tmp")).expanduser().resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)

    for acc in sra_ids:
        # temporários dedicados a este accession (evita misturar runs)
        tmp_dir = tmp_root / "sra" / acc
        tmp_dir.mkdir(parents=True, exist_ok=True)

        # Se já houver FASTQs comprimidos, não refaz
        r1_gz = Path("fastq")/f"{acc}_1.fastq.gz"
        r2_gz = Path("fastq")/f"{acc}_2.fastq.gz"
        if r1_gz.exists() or r2_gz.exists():
            console.print(f"{acc}: FASTQ(.gz) → [bold]SKIP (cache)[/bold]")
            print_meta(f"FASTQs ({acc})", [r1_gz, r2_gz])
            continue

        # Baixa o .sra (se ainda não existir)
        sra_target = next(Path("raw").glob(f"**/{acc}.sra"), None)
        if sra_target and sra_target.exists():
            console.print(f"{acc}: .sra → [bold]SKIP (cache)[/bold]")
        else:
            prefetch_with_progress(acc, Path("raw").resolve())
            # run(["prefetch", "--output-directory", str(Path("raw").resolve()), "--max-size", "u", acc])

        # CONVERSÃO com barra de progresso (já passa tmp_dir)
        sra_path = _find_local_sra(acc, Path("raw"))
        source = str(sra_path) if sra_path else acc
        ok = fasterq_with_progress(
            source, acc=acc, outdir=Path("fastq").resolve(),
            threads=threads, tmp_dir=tmp_dir,
            stall_warn_min=cfg_global["general"].get("stall_warn_min", 10),
            stall_fail_min=cfg_global["general"].get("stall_fail_min", 45),
        )

        # Se não converteu, decidir fallback ENA
        if not ok:
            if cfg_global["general"].get("ena_fallback", False):
                if ena_fetch_fastqs(acc, outdir=Path("fastq").resolve(), threads=threads):
                    console.rule(f"[green]ENA HTTP ({acc}) concluído[/green]")
                else:
                    raise RuntimeError(f"Não foi possível obter FASTQs para {acc} via fasterq nem via ENA.")
            else:
                raise RuntimeError(
                    f"fasterq {acc} não avançou e ena_fallback:false. "
                    f"Veja logs em logs/fasterq_{acc}.log e tente novamente."
                )

        # COMPACTAÇÃO com barra de progresso
        # (Se vier do ENA, já é .fastq.gz, então esta rotina apenas fará SKIP)
        compress_fastqs_with_progress(acc, outdir=Path("fastq").resolve(), threads=threads)

        print_meta(f"FASTQs ({acc})", [r1_gz, r2_gz])

def stage_fastqs_from_local(fq1, fq2=None):
    Path("fastq").mkdir(exist_ok=True)
    dst1 = Path("fastq")/Path(fq1).name
    if not dst1.exists():
        run(["bash","-lc",f"ln -s {Path(fq1).resolve()} {dst1} || cp {Path(fq1).resolve()} {dst1}"])
    else:
        console.print(f"{dst1.name}: → [bold]SKIP (cache)[/bold]")
    if fq2:
        dst2 = Path("fastq")/Path(fq2).name
        if not dst2.exists():
            run(["bash","-lc",f"ln -s {Path(fq2).resolve()} {dst2} || cp {Path(fq2).resolve()} {dst2}"])
        else:
            console.print(f"{dst2.name}: → [bold]SKIP (cache)[/bold]")
    print_meta("FASTQs de entrada (locais)", [dst1] + ([dst2] if fq2 else []))

# ======================= Downsample =======================

def downsample_fastqs(fraction: float, seed: int):
    if fraction <= 0 or fraction >= 1:
        return
    console.print(Panel.fit(f"Downsample FASTQ com seqtk (fração={fraction}, seed={seed})", border_style="cyan"))
    Path("fastq_ds").mkdir(exist_ok=True)

    paired_r1 = sorted(list(Path("fastq").glob("*_1.fastq.gz")))
    for r1 in paired_r1:
        r2 = Path(str(r1).replace("_1.fastq.gz","_2.fastq.gz"))
        base = r1.name.replace("_1.fastq.gz","")
        if r2.exists():
            out1 = Path("fastq_ds")/f"{base}_1.ds.fastq.gz"
            out2 = Path("fastq_ds")/f"{base}_2.ds.fastq.gz"
            if out1.exists() and out2.exists():
                console.print(f"{base}: downsample → [bold]SKIP (cache)[/bold]")
                continue
            run(["bash","-lc", f"seqtk sample -s{seed} {r1} {fraction} | gzip > {out1}"])
            run(["bash","-lc", f"seqtk sample -s{seed} {r2} {fraction} | gzip > {out2}"])
            print_meta(f"FASTQs downsample ({base})", [out1, out2])
        else:
            out1 = Path("fastq_ds")/f"{base}.ds.fastq.gz"
            if out1.exists():
                console.print(f"{base}: downsample (single) → [bold]SKIP (cache)[/bold]")
                continue
            run(["bash","-lc", f"seqtk sample -s{seed} {r1} {fraction} | gzip > {out1}"])
            print_meta(f"FASTQ downsample ({base})", [out1])

# ==================== QC / Trimming ====================

def qc_and_trim(threads, adapters, read_type, use_ds: bool):
    fq_dir = Path("fastq_ds") if use_ds and any(Path("fastq_ds").glob("*.fastq.gz")) else Path("fastq")
    gz = list(fq_dir.glob("*.fastq.gz")) + list(fq_dir.glob("*.fq.gz"))
    if gz:
        run(["fastqc", *map(str,gz), "-o","qc"])
        run(["multiqc","qc","-o","qc"])

    if read_type != "short":
        console.print("Leituras longas: pulo trimming por padrão.", style="dim")
        return

    fwd, rev = adapters["fwd"], adapters["rev"]
    for r1 in sorted(list(fq_dir.glob("*_1.fastq.gz")) + list(fq_dir.glob("*_1.fq.gz"))):
        base = r1.name.replace("_1.fastq.gz","").replace("_1.fq.gz","")
        r2 = Path(str(r1).replace("_1.fastq.gz","_2.fastq.gz").replace("_1.fq.gz","_2.fq.gz"))
        out1 = Path("trimmed")/f"{base}_1.trim.fq.gz"
        out2 = Path("trimmed")/f"{base}_2.trim.fq.gz" if r2.exists() else Path("trimmed")/f"{base}.trim.fq.gz"

        if out1.exists() and (not r2.exists() or out2.exists()):
            console.print(f"{base}: trimming → [bold]SKIP (cache)[/bold]")
            continue

        if r2.exists():
            run(["cutadapt","-j",str(threads),"-q","20,20","-m","30","-a",fwd,"-A",rev,"-o",str(out1),"-p",str(out2),str(r1),str(r2)])
            print_meta(f"FASTQs pós-trimming ({base})", [out1, out2])
        else:
            run(["cutadapt","-j",str(threads),"-q","20","-m","30","-a",fwd,"-o",str(out1),str(r1)])
            print_meta(f"FASTQ pós-trimming ({base})", [out1])

    trimmed = list(Path("trimmed").glob("*.fq.gz"))
    if trimmed:
        run(["fastqc", *map(str,trimmed), "-o","qc"])
        run(["multiqc","qc","-o","qc"])

# =================== Alinhamento DNA ===================

def align_dna_for_all(dna_samples, threads, default_read_type, cleanup, use_ds):
    for s in dna_samples:
        rtype = s.get("read_type", default_read_type)
        if s["source"] == "sra":
            for sra in s["sra_ids"]:
                candidates = list(Path("trimmed").glob(f"{sra}*_1.trim.fq.gz"))
                if not candidates:
                    src_dir = Path("fastq_ds") if use_ds else Path("fastq")
                    candidates = list(src_dir.glob(f"{sra}*_1.fastq.gz")) + list(src_dir.glob(f"{sra}*_1.fq.gz"))
                for r1 in candidates:
                    align_one_sample(r1, threads, rtype, cleanup)
        else:
            r1 = Path("trimmed")/Path(s["fastq1"]).name.replace(".fastq.gz",".trim.fq.gz")
            if not r1.exists():
                src_dir = Path("fastq_ds") if use_ds else Path("fastq")
                r1 = src_dir/Path(s["fastq1"]).name
            expected_r2 = Path(s.get("fastq2",""))
            align_one_sample(r1, threads, rtype, cleanup, expected_r2_name=expected_r2.name if expected_r2 else None)


def align_one_sample(r1: Path, threads: int, read_type: str, cleanup, expected_r2_name: Optional[str]=None):
    sample = r1.name.split("_")[0]
    out_sorted = Path("bam")/f"{sample}.sorted.bam"
    out_mkdup  = Path("bam")/f"{sample}.mkdup.bam"

    if out_mkdup.exists() or Path("bam")/f"{sample}.mkdup.cram".exists():
        console.print(f"[{sample}] alinhado → [bold]SKIP (cache)[/bold]")
        return

    # localizar R2
    if str(r1).endswith("_1.trim.fq.gz"):
        r2 = Path(str(r1).replace("_1.trim.fq.gz","_2.trim.fq.gz"))
    elif str(r1).endswith("_1.fastq.gz"):
        r2 = Path(str(r1).replace("_1.fastq.gz","_2.fastq.gz"))
    elif str(r1).endswith("_1.fq.gz"):
        r2 = Path(str(r1).replace("_1.fq.gz","_2.fq.gz"))
    else:
        r2 = None
    if expected_r2_name:
        cand = Path("trimmed")/expected_r2_name.replace(".fastq.gz",".trim.fq.gz").replace(".fq.gz",".trim.fq.gz")
        if cand.exists(): r2 = cand
        else:
            cand = Path("fastq")/expected_r2_name
            if cand.exists(): r2 = cand
            cand = Path("fastq_ds")/expected_r2_name.replace(".fastq.gz",".ds.fastq.gz")
            if cand.exists(): r2 = cand

    aligner = cfg_global["general"].get("aligner","bwa-mem2")
    if read_type == "short":
        base_cmd = ["bwa","mem","-t",str(threads),"refs/reference.fa"] if aligner=="bwa" else ["bwa-mem2","mem","-t",str(threads),"refs/reference.fa"]
        cmd_list = [ base_cmd + [str(r1)] + ([str(r2)] if r2 and r2.exists() else []),
                     ["samtools","view","-bS","-"],
                     ["samtools","sort","-@",str(threads),"-T", str(Path(cfg_global["general"].get("temp_dir","tmp")).expanduser().resolve() / f"{sample}.sorttmp"), "-o", str(out_sorted)] ]
        run_long_stream_pipeline(cmd_list, label=f"[{sample}] {aligner.upper()} → sort")
        run(["samtools","index",str(out_sorted)])
        run(["picard","MarkDuplicates","-I",str(out_sorted),"-O",str(out_mkdup),
             "-M",str(out_sorted).replace(".sorted.bam",".mkdup.metrics"),"--VALIDATION_STRINGENCY","LENIENT"])
        run(["samtools","index",str(out_mkdup)])
        if cleanup.get("remove_sorted_bam", True) and out_sorted.exists():
            out_sorted.unlink(missing_ok=True)
    else:
        cmd_list = [
            ["minimap2","-ax","map-ont","-t",str(threads),"refs/reference.fa",str(r1)],
            ["samtools","sort","-@",str(threads),"-T", str(Path(cfg_global["general"].get("temp_dir","tmp")).expanduser().resolve() / f"{sample}.sorttmp"), "-o", str(out_sorted)]
        ]
        run_long_stream_pipeline(cmd_list, label=f"[{sample}] minimap2 → sort")
        run(["samtools","index",str(out_sorted)])

# ============== CRAM + Cobertura (mosdepth) ==============

def to_cram_and_coverage(use_cram: bool, threads: int):
    for bam in sorted(Path("bam").glob("*.mkdup.bam")):
        sample = bam.name.replace(".mkdup.bam","")
        cram = Path("bam")/f"{sample}.mkdup.cram"
        if use_cram and not cram.exists():
            run(["samtools","view","-T","refs/reference.fa","-@",str(threads),"-C","-o",str(cram),str(bam)])
            run(["samtools","index",str(cram)])
        cov_prefix = Path("bam")/sample
        if not Path(str(cov_prefix)+".mosdepth.summary.txt").exists():
            run(["mosdepth","-t",str(threads),str(cov_prefix), str(cram if use_cram else bam)])
        summ = Path(str(cov_prefix)+".mosdepth.summary.txt")
        if summ.exists():
            with open(summ) as fh:
                lines = [l.strip().split("\t") for l in fh if l.strip() and not l.startswith("chrom")]
            mean_cov = [l for l in lines if l[0]=="total"]
            if mean_cov:
                console.print(f"[bold cyan][{sample}][/bold cyan] Cobertura média (mosdepth): [bold]{float(mean_cov[0][3]):.2f}×[/bold]")

# =================== Variantes / VEP ===================

def call_variants(samples, threads, mem_gb):
    declared_ids = {s["id"] for s in samples} if samples else set()
    aln = sorted(Path("bam").glob("*.mkdup.bam")) + sorted(Path("bam").glob("*.mkdup.cram"))
    for f in aln:
        sample = f.name.replace(".mkdup.bam","").replace(".mkdup.cram","")
        if declared_ids and sample not in declared_ids:
            continue
        gvcf = Path("vcf")/f"{sample}.g.vcf.gz"
        vcf  = Path("vcf")/f"{sample}.vcf.gz"
        if vcf.exists():
            console.print(f"[{sample}] VCF → [bold]SKIP[/bold]")
            print_meta(f"VCF ({sample})", [vcf])
            continue
        run(["gatk","--java-options",f"-Xmx{mem_gb}g","HaplotypeCaller",
             "-R","refs/reference.fa","-I",str(f),"-O",str(gvcf),"-ERC","GVCF"])
        run(["gatk","--java-options",f"-Xmx{mem_gb}g","GenotypeGVCFs",
             "-R","refs/reference.fa","-V",str(gvcf),"-O",str(vcf)])
        run(["bcftools","index","-t",str(vcf)])
        print_meta(f"VCF ({sample})", [vcf])


def annotate_variants(assembly_name, threads, vep_cache_dir=""):
    Path("vep").mkdir(exist_ok=True)
    env = os.environ.copy()
    if vep_cache_dir: env["VEP_CACHE_DIR"] = vep_cache_dir
    for vcf in sorted(Path("vcf").glob("*.vcf.gz")):
        sample = vcf.name.replace(".vcf.gz","")
        out = Path("vep")/f"{sample}.vep.vcf.gz"
        if out.exists():
            console.print(f"[{sample}] VEP → [bold]SKIP[/bold]")
            print_meta(f"VEP ({sample})", [out])
            continue
        cmd = ["vep","-i",str(vcf),"--cache","--offline","--assembly",assembly_name,
               "--format","vcf","--vcf","-o",str(out),
               "--fork",str(threads),"--everything","--no_stats"]
        try:
            run(cmd, env=env)
        except sp.CalledProcessError:
            console.print("[orange3]Aviso:[/orange3] VEP falhou (verifique cache).", style="bold")
        print_meta(f"VEP ({sample})", [out])

# =================== Lista de genes ===================

def gene_list_from_gtf():
    Path("genes").mkdir(exist_ok=True)
    out = Path("genes/gene_list.txt")
    if out.exists():
        console.print("Lista de genes → [bold]SKIP[/bold]")
        print_meta("Lista de genes", [out]); return
    cmd = r'''awk '$3=="gene"{print $0}' refs/genes.gtf | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort -u > genes/gene_list.txt'''
    run(["bash","-lc",cmd]); print_meta("Lista de genes", [out])

# =================== RNA-seq (opcional) ===================

def rnaseq_pipeline(rna_samples, threads, assembly_name):
    if not rna_samples: return
    Path("rnaseq").mkdir(exist_ok=True)
    for s in rna_samples:
        if s["source"]=="sra": stage_fastqs_from_sra(s["sra_ids"])
        else: stage_fastqs_from_local(s["fastq1"], s.get("fastq2"))

    for r1 in sorted(Path("fastq").glob("*_1.fastq.gz")):
        base = r1.name.replace("_1.fastq.gz","")
        gtf_out = Path("rnaseq")/f"{base}.transcripts.gtf"
        if gtf_out.exists():
            console.print(f"[RNA-seq:{base}] GTF → [bold]SKIP[/bold]"); continue
        r2 = Path(str(r1).replace("_1.fastq.gz","_2.fastq.gz"))
        if r2.exists():
            cmd_list = [["hisat2","-p",str(threads),"-x",f"refs/{assembly_name}","-1",str(r1),"-2",str(r2)],
                        ["samtools","sort","-@",str(threads),"-o",f"rnaseq/{base}.rnaseq.bam"]]
        else:
            cmd_list = [["hisat2","-p",str(threads),"-x",f"refs/{assembly_name}","-U",str(r1)],
                        ["samtools","sort","-@",str(threads),"-o",f"rnaseq/{base}.rnaseq.bam"]]
        run_long_stream_pipeline(cmd_list, label=f"[RNA-seq:{base}] HISAT2 → sort")
        run(["samtools","index",f"rnaseq/{base}.rnaseq.bam"])
        run(["stringtie",f"rnaseq/{base}.rnaseq.bam","-G","refs/genes.gtf","-o",str(gtf_out),"-p",str(threads)])
    tgts = list(Path("rnaseq").glob("*.transcripts.gtf"))
    if tgts and not Path("rnaseq/cmp.stats").exists():
        run(["gffcompare","-r","refs/genes.gtf","-o","rnaseq/cmp",*map(str,tgts)])

# ========================== Main ==========================

def main(cfg):
    g = cfg["general"]
    base_dir = Path(g.get("base_dir",".")).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(base_dir)
    # Configura TMPDIR global para todos os subprocessos (samtools, gatk, etc.)
    tmpdir = Path(g.get("temp_dir","tmp")).expanduser().resolve()
    tmpdir.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(tmpdir)
    console.print(f"Usando TMPDIR={tmpdir} (do YAML)", style="dim")

    console.rule(f"[bold]Pipeline Humano[/bold]  →  [cyan]{base_dir}[/cyan]")
    disk_space_report(base_dir)

    ensure_dirs()

    # Referências e índices
    download_refs(g["ref_fa_url"], g["gtf_url"], force=g.get("force_refs", False))
    limit_reference_to_canonical_if_enabled()
    build_indexes(g.get("default_read_type","short"), g["assembly_name"],
                  g.get("rnaseq",False) or bool(cfg.get("rna_samples",[])),
                  g["threads"], force=g.get("force_indexes", False))

    # Entrada DNA
    dna_samples = cfg.get("dna_samples", [])
    for s in dna_samples:
        if s["source"]=="sra": stage_fastqs_from_sra(s["sra_ids"])
        else: stage_fastqs_from_local(s["fastq1"], s.get("fastq2"))

    # Estimativa de espaço
    fq_list = list(Path("fastq").glob("*.fastq.gz"))
    estimate_outputs_and_warn(fq_list, base_dir, g.get("space_guard_gb_min", 100))

    # Downsample (opcional)
    downsample_fastqs(g.get("downsample_frac", 0.0), g.get("downsample_seed", 123))

    # QC + trimming
    qc_and_trim(g["threads"], g["adapters"], g.get("default_read_type","short"),
                use_ds=g.get("downsample_frac", 0.0) > 0)

    # Alinhamento DNA
    align_dna_for_all(dna_samples, g["threads"], g.get("default_read_type","short"),
                      cleanup=g.get("cleanup", {"remove_sorted_bam": True, "remove_bam_after_cram": True}),
                      use_ds=g.get("downsample_frac", 0.0) > 0)

    # CRAM + Cobertura
    to_cram_and_coverage(g.get("use_cram", True), g["threads"])

    # Remover BAM após CRAM
    if g.get("use_cram", True) and g.get("cleanup", {}).get("remove_bam_after_cram", True):
        for bam in Path("bam").glob("*.mkdup.bam"):
            console.print(f"Removendo {bam.name} (CRAM mantido)", style="dim")
            bam.unlink(missing_ok=True)
            idx = Path(str(bam)+".bai"); idx.unlink(missing_ok=True)

    # Variantes e VEP
    if g.get("call_variants", True):
        call_variants(dna_samples, g["threads"], g["mem_gb"])
    if g.get("annotate_vars", True):
        annotate_variants(g["assembly_name"], g["threads"], g.get("vep_cache_dir",""))

    # Lista de genes
    gene_list_from_gtf()

    # RNA-seq (opcional)
    if g.get("rnaseq", False):
        rnaseq_pipeline(cfg.get("rna_samples", []), g["threads"], g["assembly_name"])

    console.rule("[bold green]Pipeline concluído ✅[/bold green]")

if __name__=="__main__":
    ap = argparse.ArgumentParser(description="Pipeline humano — Low-Memory: índice BWA pré-pronto, verbose, cache, CRAM, cobertura e downsample")
    ap.add_argument("-c","--config",required=True,help="Arquivo YAML de configuração")
    args = ap.parse_args()
    with open(args.config, "r") as f:
        cfg_raw = yaml.safe_load(f)
    cfg = normalize_config_schema(cfg_raw)
    cfg_global = cfg
    console.print(Panel.fit(f"[bold]YAML normalizado[/bold]\nChaves topo: {list(cfg.keys())}\nGeneral → {', '.join(sorted(list(cfg['general'].keys()))[:12])}...", border_style="cyan"))
    main(cfg)

