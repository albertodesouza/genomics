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


console = Console(highlight=False, emoji=True)  # emoji=True por clareza

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

def run_long_stream_pipeline(
    cmd_list: List[List[str]],
    label: str,
    heartbeat_sec: int = 60,
    progress_spec: Optional[dict] = None,  # {"total_reads": int|None, "print_every": int}
):
    """
    Executa o pipeline com Popen, sem tocar no stdout do 1º comando (SAM),
    e monitora o *stderr* do 1º comando para inferir progresso (BWA/BWA-MEM2).

    progress_spec:
      - total_reads: total esperado de reads (int) ou None se desconhecido
      - print_every: segundos entre prints de progresso (padrão 15)
    """
    import threading, queue, re

    console.print(Panel.fit(Text(label, style="bold yellow"), border_style="yellow"))
    console.print("[bold]>[/bold] " + "  |  ".join(" ".join(map(str, c)) for c in cmd_list), style="dim")

    procs = []
    start = time.time()

    # 1) primeiro comando: stderr=PIPE (para monitor) / stdout=PIPE (vai para o próximo)
    p0 = sp.Popen(
        cmd_list[0],
        stdout=sp.PIPE,
        stderr=sp.PIPE,            # <— monitoramos só o stderr
        text=True,
        bufsize=1,
        universal_newlines=True,
    )
    procs.append(p0)

    # 2) demais comandos encadeados
    for i, cmd in enumerate(cmd_list[1:], start=1):
        p = sp.Popen(
            cmd,
            stdin=procs[-1].stdout,
            stdout=sp.PIPE if i < len(cmd_list)-1 else None,
            stderr=sp.STDOUT if i < len(cmd_list)-1 else None,
            text=False  # binário nos estágios seguintes
        )
        procs.append(p)

    # monitor de progresso (apenas se progress_spec for passado)
    stop_flag = {"stop": False}
    def _stderr_monitor():
        if not p0.stderr:
            return
        import math, time, re
        total = None
        done = 0
        last_print = time.time()
        last_done  = 0
        last_time  = last_print
        print_every = int(progress_spec.get("print_every", 15)) if progress_spec else 15
        if progress_spec and isinstance(progress_spec.get("total_reads", None), int):
            total = int(progress_spec["total_reads"])

        # BWA/BWA-MEM2 imprime lotes como: "[M::process] read 537168 sequences (80000209 bp)..."
        rx_batch = re.compile(r"\[M::process\]\s+read\s+(\d+)\s+sequences")

        while True:
            line = p0.stderr.readline()
            if not line:
                if p0.poll() is not None:
                    break
                time.sleep(0.05)
                continue

            m = rx_batch.search(line)
            if m:
                done += int(m.group(1))

            now = time.time()
            if now - last_print >= print_every:
                dt = max(1e-3, (now - last_time))
                delta = max(0, done - last_done)
                rate = delta / dt  # reads/s
                last_time = now
                last_done = done

                if total and total > 0:
                    # evita % > 100 e negativos
                    pct = min(100.0, max(0.0, 100.0 * min(done, total) / total))
                    rem = max(0, total - done)

                    # ETA só quando taxa é finita e > 0
                    if rate > 1e-6 and math.isfinite(rate):
                        eta_s = rem / rate
                        if math.isfinite(eta_s) and eta_s < 9e6:  # ~104 dias como limite de sanidade
                            eta_txt = f"{int(eta_s // 60)}m{int(eta_s % 60):02d}s"
                        else:
                            eta_txt = "--:--"
                    else:
                        eta_txt = "--:--"

                    console.print(
                        f"[magenta]{label}[/magenta] — {done:,}/{total:,} reads ({pct:4.1f}%) • "
                        f"{rate/1000:.1f} kreads/s • ETA {eta_txt}",
                        highlight=False
                    )
                else:
                    console.print(
                        f"[magenta]{label}[/magenta] — {done:,} reads • "
                        f"{rate/1000:.1f} kreads/s",
                        highlight=False
                    )
                last_print = now

    t_mon = None
    if progress_spec is not None:
        t_mon = threading.Thread(target=_stderr_monitor, daemon=True)
        t_mon.start()

    # heartbeats enquanto a *última* etapa roda
    last_hb = start
    while True:
        rc_last = procs[-1].poll()
        if rc_last is not None:
            break
        now = time.time()
        if now - last_hb >= heartbeat_sec:
            elapsed = int(now - start)
            console.print(f" {label}... {elapsed//60}m{elapsed%60:02d}s")
            last_hb = now
        time.sleep(0.5)

    # finalize/cheque return codes
    rc = 0
    for p in procs[::-1]:
        p.wait()
        rc = rc or p.returncode

    # encerra monitor
    if t_mon:
        try:
            t_mon.join(timeout=2)
        except Exception:
            pass

    if rc != 0:
        # opcional: coletar um tail do stderr do p0 para debug
        tail = ""
        try:
            if p0.stderr:
                # não dá para "tail" aqui, mas podemos informar que falhou
                tail = "BWA/BWA-MEM2 retornou erro — veja logs/ ou reexecute com '-v 3'."
        except Exception:
            pass
        raise sp.CalledProcessError(rc, cmd_list[0], output=tail)

    elapsed = int(time.time() - start)
    console.print(f"✅ [bold green]{label} concluído[/bold green] em {elapsed//60}m{elapsed%60:02d}s.")

# === Heartbeat p/ pipelines bcftools (Linux) =============================

def _descendant_pids(root_pid: int) -> set[int]:
    """Lista recursiva de PIDs descendentes via /proc/*/children (Linux)."""
    seen = set()
    stack = [root_pid]
    while stack:
        pid = stack.pop()
        try:
            with open(f"/proc/{pid}/task/{pid}/children") as fh:
                for tok in fh.read().strip().split():
                    try:
                        c = int(tok)
                    except Exception:
                        continue
                    if c not in seen:
                        seen.add(c)
                        stack.append(c)
        except Exception:
            pass
    return seen

def _sum_proc_io_tree(pid: int) -> tuple[int,int]:
    """Soma read_bytes/write_bytes do processo e seus descendentes."""
    total_r = total_w = 0
    for p in [pid, *list(_descendant_pids(pid))]:
        io = _proc_io(p)  # já existe no seu código
        if io:
            total_r += max(0, io[0])
            total_w += max(0, io[1])
    return total_r, total_w

def run_bcf_pipeline_with_heartbeat(cmd_str: str, label: str, out_watch: Path,
                                    heartbeat_sec: int = 30):
    """
    Executa 'bash -lc <cmd_str>' e, a cada heartbeat_sec, imprime:
      - tempo decorrido
      - tamanho atual de 'out_watch' (arquivo BGZF sendo escrito)
      - delta de I/O agregado (leituras/escritas) dos processos do pipeline
    """
    console.print(Panel.fit(Text(label, style="bold yellow"), border_style="yellow"))
    console.print("[bold]>[/bold] bash -lc " + cmd_str, style="dim")

    proc = sp.Popen(["bash","-lc", cmd_str], text=False)
    start = time.time()
    last = start
    last_r = last_w = 0

    while True:
        rc = proc.poll()
        now = time.time()

        if now - last >= heartbeat_sec:
            # tamanho atual do arquivo de saída (se já existir)
            out_bytes = out_watch.stat().st_size if Path(out_watch).exists() else 0
            # I/O agregado da árvore de processos
            r_tot, w_tot = _sum_proc_io_tree(proc.pid)
            d_r = max(0, r_tot - last_r); d_w = max(0, w_tot - last_w)
            last_r, last_w = r_tot, w_tot

            elapsed = int(now - start)
            console.print(
                f"{label}… {elapsed//60}m{elapsed%60:02d}s • "
                f"out {sizeof_fmt(out_bytes)} • "
                f"I/O +{sizeof_fmt(d_r)} read, +{sizeof_fmt(d_w)} write",
                highlight=False
            )
            last = now

        if rc is not None:
            break
        time.sleep(0.5)

    if proc.returncode != 0:
        raise sp.CalledProcessError(proc.returncode, ["bash","-lc", cmd_str])

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

def stage_banner(title: str, sub: str = ""):
    t = f"[bold white]{title}[/bold white]"
    if sub:
        t += f"\n[dim]{sub}[/dim]"
    console.rule()
    console.print(Panel.fit(t, border_style="bright_blue"))
    console.rule()

def step_done(msg: str):
    console.print(f"✅ [bold green]{msg}[/bold green]")

def _is_newer(target: Path, *sources: Path) -> bool:
    "True se target existe e é mais novo que TODAS as sources."
    if not target.exists(): return False
    t = target.stat().st_mtime
    return all((not s.exists()) or t >= s.stat().st_mtime for s in sources)

def _atomic_rename(tmp: Path, final: Path):
    tmp.replace(final)

def _ensure_tabix(vcf_gz: Path):
    tbi = Path(str(vcf_gz)+".tbi")
    if not tbi.exists():
        run(["tabix","-p","vcf",str(vcf_gz)])
    return tbi

def _read_fai(fai: Path):
    "Retorna lista de tuplas (contig, length). Gera .fai se faltar."
    if not fai.exists():
        run(["samtools","faidx","refs/reference.fa"])
    rows = []
    with fai.open() as fh:
        for line in fh:
            ctg, ln, *_ = line.strip().split("\t")
            rows.append((ctg, int(ln)))
    return rows

def _canonical_subset(contigs):
    "Seleciona apenas chr1-22,X,Y,M/MT ou sem 'chr' (1-22,X,Y,M/MT)."
    keep = []
    valid = {str(i) for i in range(1,23)} | {"X","Y","M","MT"}
    for ctg,_ in contigs:
        base = ctg[3:] if ctg.lower().startswith("chr") else ctg
        if base.upper() in valid:
            keep.append((ctg,_))
    return keep

def _filter_problematic_contigs(contigs):
    """
    Remove cromossomos/contigs problemáticos que podem causar falhas no BCFtools/VEP.
    Filtra cromossomos virais, vetores e contigs problemáticos conhecidos.
    """
    # Cromossomos problemáticos conhecidos (virais, vetores, etc.)
    problem_patterns = {
        # Vírus (mais abrangente)
        "EBV", "CMV", "HBV", "HCV", "HIV", "HPV", "HTLV", "KSHV", "SV40", "HSV",
        "HERPES", "POLYOMA", "ADENO", "PAPILLOMA", "RETRO", "LENTI",
        # Vírus específicos encontrados nos dados
        "HCV-1", "HCV-2", "HIV-1", "HIV-2", "HTLV-1", "HTLV-2",
        # Vetores e plasmídeos
        "pUC", "pBR", "pET", "pcDNA", "VECTOR", "PLASMID",
        # Contigs decoy problemáticos (causam reference mismatch)
        "_DECOY", "DECOY", "JTFH01", "KN707", "KI270", "GL000",
        # Prefixos de contigs não mapeados problemáticos
        "CHRUN_", "chrUn_",
        # Outros problemáticos
        "LAMBDA", "PHAGE", "SYNTHETIC", "ARTIFICIAL"
    }
    
    keep = []
    filtered_count = 0
    
    for ctg, length in contigs:
        ctg_upper = ctg.upper()
        is_problematic = any(pattern in ctg_upper for pattern in problem_patterns)
        
        # Verifica padrões específicos adicionais
        if not is_problematic:
            # Contigs decoy específicos
            if ctg.startswith("chrUn_") and ("_decoy" in ctg or "JTFH01" in ctg or "KN707" in ctg):
                is_problematic = True
            # Contigs KI270 problemáticos
            elif "KI270" in ctg_upper:
                is_problematic = True
            # Contigs GL000 problemáticos  
            elif "GL000" in ctg_upper:
                is_problematic = True
        
        if not is_problematic:
            keep.append((ctg, length))
        else:
            filtered_count += 1
    
    if filtered_count > 0:
        console.print(f"[yellow]⚠️  Filtrados {filtered_count} cromossomos problemáticos (virais/vetores)[/yellow]")
        console.print("[dim]Cromossomos filtrados podem causar falhas em BCFtools/VEP[/dim]")
    
    return keep

def _write_intervals_file(contigs, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Escreve BED: chrom  start(0)  end(len)
    with dest.open("w") as fo:
        for ctg, ln in contigs:
            fo.write(f"{ctg}\t0\t{ln}\n")
    return dest

def _human_bp(n: int) -> str:
    units = ["bp","kb","Mb","Gb","Tb"]; i = 0; x = float(n)
    while x >= 1000 and i < len(units)-1:
        x /= 1000.0; i += 1
    return f"{x:.1f} {units[i]}"

def _bp_in_bed(bed: Path) -> int:
    total = 0
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            a = line.split("\t")
            if len(a) >= 3:
                try: total += max(0, int(a[2]) - int(a[1]))
                except Exception: pass
    return total

def _split_balanced(contigs, parts: int):
    "Greedy balance por tamanho; retorna lista de listas de (ctg,len)."
    parts = max(1, int(parts))
    bins = [([], 0)] * parts
    bins = [ ([],0) for _ in range(parts) ]
    for ctg, ln in sorted(contigs, key=lambda x: x[1], reverse=True):
        i = min(range(parts), key=lambda k: bins[k][1])
        bins[i][0].append((ctg,ln))
        bins[i] = (bins[i][0], bins[i][1] + ln)
    return [b[0] for b in bins]

def _tmp_vcfgz_path(final: Path) -> Path:
    # garante que o sufixo final permaneça .vcf.gz (para GATK entender BGZF)
    suf2 = ''.join(final.suffixes[-2:])  # normalmente ".vcf.gz"
    if suf2 != ".vcf.gz":
        # fallback defensivo: ainda tenta manter .gz no fim
        if final.suffix == ".gz":
            base = final.name[:-len(".gz")]
            return final.with_name(base + ".tmp" + ".gz")
        # pior caso, sem .gz → semear .tmp.vcf.gz
        return final.with_name(final.name + ".tmp.vcf.gz")
    base = final.name[:-len(suf2)]
    return final.with_name(base + ".tmp" + suf2)

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
    # já normalizado?
    if "general" in cfg_in and ("dna_samples" in cfg_in or "rna_samples" in cfg_in):
        return cfg_in

    project   = cfg_in.get("project", {})
    storage   = cfg_in.get("storage", {})
    ref_top   = cfg_in.get("reference", {})
    ref_proj  = project.get("reference", {}) if isinstance(project, dict) else {}
    ref       = {**ref_proj, **ref_top}
    execv     = cfg_in.get("execution", {})
    download  = cfg_in.get("download", {})
    sizec     = cfg_in.get("size_control", {})
    samples   = cfg_in.get("samples", [])
    params    = cfg_in.get("params", {})

    # caminhos / referência
    base_dir = storage.get("base_dir", ".")
    temp_dir = storage.get("temp_dir", "tmp")
    assembly = ref.get("name", "GRCh38")
    ref_fa_url = ref.get("fasta_url") or \
                 "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
    gtf_url    = ref.get("gtf_url") or \
                 "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"
    bwa_prebuilt = ref.get("bwa_index_url") or cfg_in.get("bwa_prebuilt_url") or ""

    # opções gerais
    prefetch_retries = int(download.get("prefetch_retries", 2))
    prefer_ena_fastq = bool(execv.get("prefer_ena_fastq", False))
    cancel_on_convert_stall = bool(execv.get("cancel_on_convert_stall", False))
    stall_warn_min = int(execv.get("stall_warn_min", 10))
    stall_fail_min = int(execv.get("stall_fail_min", 45))
    ena_fallback   = bool(execv.get("ena_fallback", False))
    download_tool  = (download.get("tool") or "sra_toolkit").lower()
    if download_tool != "sra_toolkit":
        console.print(f"[orange3]Aviso:[/orange3] download.tool='{download_tool}' ainda não é suportado; usando sra-tools.", style="italic")

    # threads “globais” (download, fastqc, multiqc, mosdepth, etc.)
    # Mantemos a precedência antiga p/ compat: execution.threads > download.threads > (fallback 16)
    threads_global = int(execv.get("threads") or download.get("threads") or 16)

    # memória para GATK, etc.
    mem_gb  = int(params.get("mem_gb", 64))

    # downsample
    ds_cfg = (sizec or {}).get("downsample", {})
    downsample_frac = float(ds_cfg.get("fraction", 0.0)) if ds_cfg.get("enabled", False) else 0.0
    downsample_seed = int(ds_cfg.get("seed", 123))

    # limpeza / CRAM
    keep_inter = bool(storage.get("keep_intermediates", False))
    use_cram = bool(cfg_in.get("general", {}).get("use_cram", True))
    cleanup = dict(cfg_in.get("general", {}).get("cleanup", {})) or {
        "remove_sorted_bam": True,
        "remove_bam_after_cram": (not keep_inter)
    }

    # adaptadores default
    adapters = {"fwd": "AGATCGGAAGAGC", "rev": "AGATCGGAAGAGC"}

    # alinhador + threads específicas p/ alinhamento
    aligner = (params.get("aligner") or execv.get("aligner") or "bwa-mem2").lower()
    aligner = "bwa-mem2" if aligner in ("bwa-mem2", "bwa_mem2", "mem2") else "bwa"

    # >>>>> ÚNICA PARTE ALTERADA: precedência do aln_threads <<<<<
    # Prioridade: general.aln_threads > params.aln_threads >
    #             (params.bwa_mem2_threads|params.bwa_threads) >
    #             execution.aln_threads > execution.threads > download.threads > 16
    cand_aln_threads = [
        (cfg_in.get("general") or {}).get("aln_threads"),
        params.get("aln_threads"),
        params.get("bwa_mem2_threads") if aligner == "bwa-mem2" else params.get("bwa_threads"),
        execv.get("aln_threads"),
        execv.get("threads"),
        download.get("threads"),
        16,
    ]
    aln_threads = next(int(v) for v in cand_aln_threads if v is not None)
    aln_threads = max(1, aln_threads)
    # <<<<< FIM DA ALTERAÇÃO >>>>>

    # knobs de RAM do alinhamento / sort
    sort_mem_mb  = int(cfg_in.get("general", {}).get("sort_mem_mb", 384))
    bwa_batch_k  = int(cfg_in.get("general", {}).get("bwa_batch_k", 20000000))  # -K

    # canônicos
    limit_to_canonical = bool(cfg_in.get("limit_to_canonical", False))

    # amostras DNA
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

    # bloco general (saída)
    general = {
        "base_dir": base_dir,
        "assembly_name": assembly,
        "ref_fa_url": ref_fa_url,
        "gtf_url": gtf_url,

        # threads globais vs. threads de alinhamento
        "threads": threads_global,
        "aln_threads": aln_threads,

        "mem_gb": mem_gb,
        "default_read_type": "short",
        "adapters": adapters,

        "call_variants": True,
        "annotate_vars": True,
        "rnaseq": False,

        # flags de força
        "force_refs": bool(cfg_in.get("general", {}).get("force_refs", False)),
        "force_indexes": bool(cfg_in.get("general", {}).get("force_indexes", False)),

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
        "prefetch_retries": prefetch_retries,
        "prefer_ena_fastq": prefer_ena_fastq,
        "cancel_on_convert_stall": cancel_on_convert_stall,
        "ena_fallback": ena_fallback,

        "temp_dir": temp_dir,

        # novos knobs expostos para fases pesadas de RAM
        "sort_mem_mb": sort_mem_mb,
        "bwa_batch_k": bwa_batch_k,
    }

    # --- Pass-through adicional de chaves opcionais úteis vindas do YAML "novo" ---
    opt_keys = [
        # trio / filtros
        "trio_child_id","trio_parent_ids","trio_min_dp_child","trio_min_dp_parents",
        "trio_min_gq","trio_min_ab_het","trio_max_ab_het","trio_min_ab_hom","trio_max_parent_alt_frac",
        # presença de genes
        "gene_presence_min_mean_cov","gene_presence_min_breadth_1x",
        # robustez de download/conversão
        "stall_warn_min","stall_fail_min","cancel_on_convert_stall",
        "prefer_ena_fastq","ena_fallback","prefetch_retries",
        # tmp customizado e knobs de RAM
        "temp_dir","sort_mem_mb","bwa_batch_k",
        # força / limpeza
        "force_refs","force_indexes","use_cram","cleanup",
    ]
    for src in (execv, storage, cfg_in.get("general", {}), cfg_in):
        if isinstance(src, dict):
            for k in opt_keys:
                if k in src:
                    general[k] = src[k]

    return {
        "general": general,
        "dna_samples": dna_samples,
        "rna_samples": [],
        "params": dict(params)  # <-- preserva o bloco params
    }

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
    Instala o índice BWA pré-construído (GDC) de forma robusta e idempotente:
      - Fast-path: se os 5 links refs/reference.fa.{amb,ann,bwt,pac,sa} já existem
        e apontam para alvos válidos → SKIP imediato.
      - Baixa com cache (wget -c) só se precisar.
      - Extrai do tar SOMENTE se faltar arquivo (ou force=True).
      - Cria/reusa symlinks apontando para caminhos absolutos.
      - Só roda chmod -R quando houve extração.
    """
    import os, re, shlex
    from pathlib import Path
    import subprocess as sp

    exts = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    ref_prefix = Path("refs/reference.fa")
    index_dir = Path("refs/_bwa")
    index_dir.mkdir(parents=True, exist_ok=True)

    def _links_ok() -> bool:
        for e in exts:
            p = Path(str(ref_prefix) + e)
            if not p.exists(): 
                return False
            if p.is_symlink():
                try:
                    _ = p.resolve(strict=True)
                except FileNotFoundError:
                    return False
            else:
                # preferimos que sejam symlinks, mas se for arquivo regular válido, aceite
                if not p.is_file():
                    return False
        return True

    def _any_prefix_in_dir() -> Path | None:
        """Tenta achar rapidamente um prefixo *.fa sob refs/_bwa (usa só o .amb)."""
        amb = next(index_dir.rglob("*.fa.amb"), None)
        return amb.with_suffix("") if amb else None

    def _prefix_has_all(prefix: Path) -> bool:
        for e in exts:
            if not Path(str(prefix) + e).exists():
                return False
        return True

    def _uuid_basename_from_url(url: str, default: str) -> str:
        m = re.search(r'([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$', url, re.I)
        return (m.group(1) + ".tar.gz") if m else default

    # ---------- Fast path: tudo pronto? ----------
    if not force and _links_ok():
        console.print("Índice BWA já presente e válido → [bold]SKIP[/bold]")
        print_meta("Índice BWA linkado", [Path(str(ref_prefix)+e) for e in exts])
        return

    if not bwa_tar_url:
        # Sem URL e sem links prontos → erro claro
        if not _links_ok():
            raise RuntimeError("Sem bwa_tar_url e índice BWA ausente/incompleto.")
        return

    tar_name = _uuid_basename_from_url(bwa_tar_url, "index.tar.gz")
    tar_path = index_dir / tar_name

    # ---------- Download (cache) ----------
    if force or not tar_path.exists():
        run(["bash","-lc", f"wget -c -O {shlex.quote(str(tar_path))} {shlex.quote(bwa_tar_url)}"])
    else:
        console.print(f"{tar_name} já presente → [bold]SKIP download[/bold]")

    # ---------- MD5 opcional ----------
    if expect_md5:
        md5 = sp.run(
            ["bash","-lc", f"md5sum {shlex.quote(str(tar_path))} | cut -d' ' -f1"],
            capture_output=True, text=True, check=True
        ).stdout.strip()
        if md5 != expect_md5:
            raise RuntimeError(f"MD5 não confere para {tar_name}: {md5} (esp. {expect_md5})")

    # ---------- Precisamos extrair? ----------
    # Se já houver um prefixo completo em refs/_bwa e não for force, reaproveite.
    best_prefix: Path | None = None
    if not force:
        p = _any_prefix_in_dir()
        if p and _prefix_has_all(p):
            best_prefix = p

    extracted = False
    if best_prefix is None:
        # (1) listar rapidamente os caminhos que queremos (uma vez só)
        lst = sp.run(
            ["bash","-lc", f"tar -tzf {shlex.quote(str(tar_path))}"],
            capture_output=True, text=True, check=True
        ).stdout.splitlines()
        wanted = [x for x in lst if re.search(r"\.fa\.(amb|ann|bwt|pac|sa)$", x)]
        if not wanted:
            # diagnóstico amigável
            head = "\n".join(lst[:10])
            console.print(Panel.fit(
                "Tar não contém *.fa.{amb,ann,bwt,pac,sa}.\n"
                "Pré-visualização (head):\n" + head,
                border_style="red"
            ))
            raise RuntimeError("Tar sem arquivos de índice do BWA.")

        wanted = [re.sub(r"^\./", "", x) for x in wanted]
        names_quoted = " ".join(shlex.quote(x) for x in wanted)

        # (2) extração só dos 5 arquivos (overwrite se preciso)
        rc = sp.run(
            ["bash","-lc",
             f"tar -xzf {shlex.quote(str(tar_path))} -C {shlex.quote(str(index_dir))} --overwrite {names_quoted}"],
            stdout=sp.DEVNULL, stderr=sp.DEVNULL
        ).returncode
        if rc != 0:
            run(["bash","-lc",
                 f"tar -xzf {shlex.quote(str(tar_path))} -C {shlex.quote(str(index_dir))} {names_quoted}"])
        extracted = True

        # (3) deduz prefixo agora existente
        amb = next(index_dir.rglob("*.fa.amb"), None)
        if not amb:
            raise RuntimeError("Falha ao localizar *.fa.amb após extração.")
        best_prefix = amb.with_suffix("")
        if not _prefix_has_all(best_prefix):
            raise RuntimeError("Índice BWA incompleto após extração.")

    # ---------- (Re)link inteligente ----------
    from os.path import samefile
    for e in exts:
        target = Path(str(best_prefix) + e).resolve()
        if not target.exists():
            raise RuntimeError(f"Arquivo do índice ausente: {target}")
        dst = Path(str(ref_prefix) + e)
        if dst.is_symlink():
            try:
                if samefile(dst, target):
                    # link já correto → não mexe
                    pass
                else:
                    dst.unlink(missing_ok=True)
                    dst.symlink_to(target)
            except FileNotFoundError:
                # alvo sumiu → recria
                dst.unlink(missing_ok=True)
                dst.symlink_to(target)
        elif dst.exists():
            # arquivo regular (raro) → substitui por symlink
            dst.unlink(missing_ok=True)
            dst.symlink_to(target)
        else:
            dst.symlink_to(target)

    # permissões só quando houve extração (evita custo em cada run)
    if extracted:
        run(["bash","-lc", f"chmod -R a+r {shlex.quote(str(index_dir))} || true"])

    # ---------- Validação final + sumário ----------
    if not all(Path(str(ref_prefix)+e).exists() for e in exts):
        raise RuntimeError("Após relink, o índice BWA continua incompleto.")
    console.print("Índice BWA linkado → [bold]pronto[/bold]")
    print_meta("Índice BWA linkado", [Path(str(ref_prefix)+e) for e in exts])

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


def _build_bwa_index_optimized(ref_prefix: Path, label: str = "BWA"):
    """
    Constrói índice BWA otimizado baseado na RAM disponível e parâmetros YAML.
    - Detecta RAM total do sistema
    - Usa parâmetros configuráveis do YAML
    - Monitora progresso durante construção
    - Otimiza block size para máquinas monster
    """
    import subprocess as sp
    import time
    import psutil
    from pathlib import Path
    
    g = cfg_global.get("general", {})
    p = cfg_global.get("params", {})
    
    # Parâmetros configuráveis
    bwa_index_max_mem_gb = int(p.get("bwa_index_max_mem_gb", 100))  # Máximo de RAM a usar
    bwa_index_block_size = p.get("bwa_index_block_size", "auto")    # Block size ou "auto"
    bwa_index_algorithm = p.get("bwa_index_algorithm", "bwtsw")     # Algoritmo
    bwa_index_progress_sec = int(p.get("bwa_index_progress_sec", 60))  # Intervalo de progresso
    
    # Detecta RAM total
    total_ram_gb = psutil.virtual_memory().total / (1024**3)
    
    # Calcula block size otimizado
    if bwa_index_block_size == "auto":
        if total_ram_gb >= 200:
            # Máquinas monster (200GB+): usar até 100GB para indexação
            block_size = min(2000000000, bwa_index_max_mem_gb * 20000000)  # 20M por GB
        elif total_ram_gb >= 64:
            # Máquinas potentes (64GB+): usar proporcionalmente
            block_size = min(1000000000, int(total_ram_gb * 15000000))     # 15M por GB
        else:
            # Máquinas normais: usar padrão
            block_size = 10000000
    else:
        block_size = int(bwa_index_block_size)
    
    # Limita RAM usada ao configurado
    estimated_ram_gb = block_size / 50000000  # Aproximação: 50M block = 1GB RAM
    if estimated_ram_gb > bwa_index_max_mem_gb:
        block_size = bwa_index_max_mem_gb * 50000000
        estimated_ram_gb = bwa_index_max_mem_gb
    
    console.print(Panel.fit(
        f"[bold]Criando Índice BWA Otimizado ({label})[/bold]\n"
        f"• RAM total: {total_ram_gb:.0f}GB\n"
        f"• RAM para indexação: ~{estimated_ram_gb:.0f}GB\n"
        f"• Block size: {block_size:,}\n"
        f"• Algoritmo: {bwa_index_algorithm}\n"
        f"• Tempo estimado: {45 if total_ram_gb >= 200 else 90}-{90 if total_ram_gb >= 200 else 180}min",
        border_style="cyan"
    ))
    
    # Comando otimizado
    cmd = ["bwa", "index", "-a", bwa_index_algorithm, "-b", str(block_size), str(ref_prefix)]
    
    console.print(f"[bold]💻 Comando executado:[/bold]")
    console.print(f"[yellow]> {' '.join(cmd)}[/yellow]")
    
    # Executa com monitoramento de progresso
    start_time = time.time()
    last_progress = start_time
    
    console.print(f"[cyan]🚀 Iniciando criação do índice BWA...[/cyan]")
    
    # Executa em background para monitorar
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    
    while True:
        rc = proc.poll()
        current_time = time.time()
        
        # Progress heartbeat
        if current_time - last_progress >= bwa_index_progress_sec:
            elapsed = int(current_time - start_time)
            
            # Verifica arquivos BWA sendo criados e calcula total
            files_info = []
            total_size = 0
            
            for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
                file_path = Path(str(ref_prefix) + ext)
                if file_path.exists():
                    size = file_path.stat().st_size
                    total_size += size
                    files_info.append(f"{ext[1:]}({sizeof_fmt(size)})")
            
            # Status compacto com total
            if files_info:
                status = f"total {sizeof_fmt(total_size)} • {', '.join(files_info[:3])}"
                if len(files_info) > 3:
                    status += f" +{len(files_info)-3} mais"
            else:
                status = "iniciando..."
                
            console.print(
                f"[cyan]BWA index … {elapsed//60}m{elapsed%60:02d}s • {status}[/cyan]",
                highlight=False
            )
            last_progress = current_time
        
        if rc is not None:
            break
            
        time.sleep(5)
    
    # Verifica resultado
    if proc.returncode != 0:
        stdout, stderr = proc.communicate()
        console.print(f"[red]❌ BWA index falhou com código {proc.returncode}[/red]")
        if stderr:
            console.print(f"[red]Erro:[/red] {stderr}")
        raise sp.CalledProcessError(proc.returncode, cmd)
    
    # Relatório final
    elapsed = int(time.time() - start_time)
    total_size = sum(
        Path(str(ref_prefix) + ext).stat().st_size 
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]
        if Path(str(ref_prefix) + ext).exists()
    )
    
    console.print(
        f"[bold green]✅ Índice BWA criado em {elapsed//60}m{elapsed%60:02d}s • "
        f"tamanho total: {sizeof_fmt(total_size)}[/bold green]"
    )

def _build_bwa_mem2_index_optimized(ref_prefix: Path, label: str = "BWA-MEM2"):
    """
    Constrói índice BWA-MEM2 com monitoramento de progresso.
    NOTA: BWA-MEM2 não aceita parâmetros de otimização como -b (block size),
    mas usa automaticamente toda RAM disponível e é mais eficiente que BWA clássico.
    """
    import subprocess as sp
    import time
    import psutil
    from pathlib import Path
    
    g = cfg_global.get("general", {})
    p = cfg_global.get("params", {})
    
    # Parâmetros configuráveis
    bwa_mem2_progress_sec = int(p.get("bwa_index_progress_sec", 60))
    
    # Detecta RAM total
    total_ram_gb = psutil.virtual_memory().total / (1024**3)
    estimated_ram_gb = min(total_ram_gb * 0.8, 200)  # BWA-MEM2 usa automaticamente muita RAM
    
    console.print(Panel.fit(
        f"[bold]Criando Índice BWA-MEM2 ({label})[/bold]\n"
        f"• RAM total: {total_ram_gb:.0f}GB\n"
        f"• RAM que será usada: ~{estimated_ram_gb:.0f}GB (automático)\n"
        f"• Algoritmo: BWA-MEM2 (otimizado internamente)\n"
        f"• Tempo estimado: {20 if total_ram_gb >= 200 else 45}-{45 if total_ram_gb >= 200 else 90}min\n"
        f"• NOTA: BWA-MEM2 não aceita parâmetros de block size",
        border_style="green"
    ))
    
    # Comando BWA-MEM2 (sem parâmetros extras - ele otimiza internamente)
    cmd = ["bwa-mem2", "index", str(ref_prefix)]
    
    console.print(f"[bold]💻 Comando executado:[/bold]")
    console.print(f"[yellow]> {' '.join(cmd)}[/yellow]")
    
    # Executa com monitoramento de progresso
    start_time = time.time()
    last_progress = start_time
    
    console.print(f"[green]🚀 Iniciando criação do índice BWA-MEM2...[/green]")
    
    # Executa em background para monitorar
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    
    while True:
        rc = proc.poll()
        current_time = time.time()
        
        # Progress heartbeat
        if current_time - last_progress >= bwa_mem2_progress_sec:
            elapsed = int(current_time - start_time)
            
            # Verifica arquivos BWA-MEM2 sendo criados e calcula total
            files_info = []
            total_size = 0
            
            # BWA-MEM2 cria arquivos diferentes
            bwa_mem2_exts = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
            for ext in bwa_mem2_exts:
                file_path = Path(str(ref_prefix) + ext)
                if file_path.exists():
                    size = file_path.stat().st_size
                    total_size += size
                    # Simplifica nome da extensão para display
                    ext_name = ext.replace('.', '').replace('2bit64', '2bit')
                    files_info.append(f"{ext_name}({sizeof_fmt(size)})")
            
            # Verifica também padrões alternativos do BWA-MEM2
            for file_path in Path("refs").glob("reference.fa.*"):
                if file_path.is_file() and str(file_path) not in [str(ref_prefix) + ext for ext in bwa_mem2_exts]:
                    size = file_path.stat().st_size
                    total_size += size
                    ext_name = file_path.suffix.replace('.', '')
                    files_info.append(f"{ext_name}({sizeof_fmt(size)})")
            
            # Status compacto
            if files_info:
                status = f"total {sizeof_fmt(total_size)} • {', '.join(files_info[:3])}"
                if len(files_info) > 3:
                    status += f" +{len(files_info)-3} mais"
            else:
                status = "iniciando..."
                
            console.print(
                f"[green]BWA-MEM2 index … {elapsed//60}m{elapsed%60:02d}s • {status}[/green]",
                highlight=False
            )
            last_progress = current_time
        
        if rc is not None:
            break
            
        time.sleep(5)
    
    # Verifica resultado
    if proc.returncode != 0:
        stdout, stderr = proc.communicate()
        console.print(f"[red]❌ BWA-MEM2 index falhou com código {proc.returncode}[/red]")
        if stderr:
            console.print(f"[red]Erro:[/red] {stderr}")
        raise sp.CalledProcessError(proc.returncode, cmd)
    
    # Relatório final
    elapsed = int(time.time() - start_time)
    
    # Calcula tamanho total dos arquivos BWA-MEM2
    total_size = 0
    mem2_files = list(Path("refs").glob("reference.fa.*"))
    for file_path in mem2_files:
        if file_path.is_file():
            total_size += file_path.stat().st_size
    
    console.print(
        f"[bold green]✅ Índice BWA-MEM2 criado em {elapsed//60}m{elapsed%60:02d}s • "
        f"tamanho total: {sizeof_fmt(total_size)}[/bold green]"
    )

def build_indexes(default_read_type, assembly_name, need_rna_index, threads, force=False):
    """
    Prepara índices de referência para o(s) alinhador(es):
      - Se general.bwa_prebuilt_url estiver definido → instala/relinca índice BWA pré-pronto do GDC.
      - Se aligner == bwa-mem2:
          * usa índice do mem2 se já existir;
          * senão tenta "bwa-mem2 index"; se falhar por RAM, faz fallback p/ BWA clássico otimizado.
      - Se aligner == bwa (sem prebuilt): garante índice clássico otimizado.
    Também monta índice HISAT2 quando need_rna_index=True.

    Parâmetros YAML para otimização BWA index (params):
      bwa_index_max_mem_gb: int     # máximo de RAM a usar (default: 100GB)
      bwa_index_block_size: str|int # "auto" ou valor específico (default: "auto")
      bwa_index_algorithm: str      # "bwtsw" para genomas grandes (default: "bwtsw")
      bwa_index_progress_sec: int   # intervalo de updates (default: 60s)

    Observação importante:
      - NÃO use limit_to_canonical com índice prebuilt do GDC (mismatch de contigs).
      - Em máquinas monster (200GB+), pode usar até 120GB para indexação 3x mais rápida.
    """
    from pathlib import Path
    import subprocess as sp

    g = cfg_global["general"]

    def _bwa_index_ready(prefix: Path) -> bool:
        """5 arquivos do índice BWA, apontando para destinos válidos."""
        for ext in (".amb",".ann",".bwt",".pac",".sa"):
            p = Path(str(prefix) + ext)
            if not p.exists():
                return False
            if p.is_symlink() and not p.resolve().exists():
                return False
        return True

    def _mem2_index_present(prefix: Path) -> bool:
        """Sinais de índice do bwa-mem2 (qualquer variante gerada)."""
        if Path(str(prefix) + ".bwt.2bit.64").exists():
            return True
        # alguns builds geram sufixos *.0123 / *.amb.* etc.
        return any(Path("refs").glob("reference.fa.*.0123"))

    # nada a fazer para long reads (minimap2) se não for RNA index
    if default_read_type != "short" and not need_rna_index:
        console.print("Sem necessidade de indexar (long reads).", style="dim")
        return

    ref_prefix = Path("refs/reference.fa")
    aligner_cfg = (g.get("aligner") or "bwa-mem2").lower()
    prebuilt_url = g.get("bwa_prebuilt_url", "")

    # alerta de canônicos x prebuilt
    if prebuilt_url and g.get("limit_to_canonical", False):
        console.print(
            "[orange3]Aviso:[/orange3] você definiu [bold]limit_to_canonical:true[/bold] "
            "E também um [bold]bwa_prebuilt_url[/bold]. Isso causa mismatch entre FASTA e índice. "
            "Desative limit_to_canonical ao usar o índice prebuilt do GDC.", style="italic"
        )

    # ------------------ fluxo do índice BWA/BWA-MEM2 ------------------
    if default_read_type == "short":
        # Caso 1: índice prebuilt explícito — sempre priorize
        if prebuilt_url:
            install_prebuilt_bwa_index(prebuilt_url, force=force)
            if not _bwa_index_ready(ref_prefix):
                raise RuntimeError("Índice BWA não encontrado ou incompleto após instalar o prebuilt.")
            console.print("Índice BWA pré-instalado e linkado → [bold]pronto[/bold]")
        else:
            # Sem prebuilt: agir conforme o alinhador
            if aligner_cfg in ("bwa-mem2", "bwa_mem2", "mem2"):
                # se já existe índice do mem2 e não for force → SKIP
                if not force and _mem2_index_present(ref_prefix):
                    console.print("Índice BWA-MEM2 → [bold]SKIP[/bold]", style="dim")
                else:
                    # tentar construir índice do mem2 (pode falhar por RAM)
                    try:
                        _build_bwa_mem2_index_optimized(ref_prefix, "BWA-MEM2")
                        console.print("Índice BWA-MEM2 construído.", style="dim")
                    except sp.CalledProcessError:
                        console.print("[red]Falha ao indexar com bwa-mem2 (provável falta de RAM).[/red]")
                        console.print(
                            "Fallback: tentando preparar índice do [bold]BWA clássico[/bold]. "
                            "Sugestões: use [bold]bwa_prebuilt_url[/bold] no YAML, ou reduza a referência.",
                            style="italic"
                        )
                        # se houver prebuilt, usar agora
                        if g.get("bwa_prebuilt_url"):
                            install_prebuilt_bwa_index(g["bwa_prebuilt_url"], force=True)
                            if not _bwa_index_ready(ref_prefix):
                                raise RuntimeError("Prebuilt BWA apontado não ficou completo após fallback.")
                            console.print("Índice BWA pré-instalado e linkado (fallback) → [bold]pronto[/bold]")
                        else:
                            # tentar construir índice do BWA clássico localmente
                            try:
                                # Para genomas grandes, -a bwtsw é o recomendado p/ BWA clássico
                                _build_bwa_index_optimized(ref_prefix, "fallback BWA clássico")
                                if not _bwa_index_ready(ref_prefix):
                                    raise RuntimeError("Índice BWA clássico parece incompleto após 'bwa index'.")
                                console.print("Índice BWA clássico construído (fallback).", style="dim")
                            except sp.CalledProcessError:
                                raise RuntimeError(
                                    "Não foi possível preparar nenhum índice (mem2 falhou e 'bwa index' também). "
                                    "Defina 'project.reference.bwa_index_url' no YAML para usar o índice pré-pronto do GDC."
                                )
                # Observação: o alinhamento pode cair para BWA se o índice do mem2 não existir.
                # Sua align_one_sample já trata esse fallback em runtime.
            else:
                # aligner == bwa (clássico) — garantir índice clássico
                if force or not _bwa_index_ready(ref_prefix):
                    if not force and _bwa_index_ready(ref_prefix):
                        console.print("Índice BWA → [bold]SKIP[/bold]", style="dim")
                    else:
                        # construir localmente
                        _build_bwa_index_optimized(ref_prefix, "BWA clássico")
                        if not _bwa_index_ready(ref_prefix):
                            raise RuntimeError("Índice BWA clássico parece incompleto após 'bwa index'.")
                        console.print("Índice BWA clássico construído.", style="dim")
                else:
                    console.print("Índice BWA → [bold]SKIP[/bold]", style="dim")

    # ------------------ índice HISAT2 (RNA) ------------------
    if need_rna_index:
        ht2_base = Path(f"refs/{assembly_name}")
        if force or not Path(f"{ht2_base}.1.ht2").exists():
            run(["hisat2-build","-p",str(threads), str(ref_prefix), str(ht2_base)])
        else:
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

        # tmp dedicado (no disco rápido do YAML)
        tmp_root = Path(cfg_global.get("general", {}).get("temp_dir", "tmp")).expanduser().resolve()
        tmp_dir = tmp_root / "sra" / acc
        tmp_dir.mkdir(parents=True, exist_ok=True)

        threads = int(cfg_global.get("general", {}).get("threads", 8))

        # 0) Preferir ENA direto? (pula SRA por completo)
        if cfg_global["general"].get("prefer_ena_fastq", False):
            console.print(f"[dim]{acc}: prefer_ena_fastq=true → tentando ENA primeiro[/dim]")
            if ena_fetch_fastqs(acc, outdir=Path("fastq").resolve(), threads=threads):
                console.rule(f"[green]ENA HTTP ({acc}) concluído[/green]")
                # ainda roda compressão (vai dar SKIP) e segue o fluxo normal
            else:
                console.print(f"[yellow]{acc}: ENA não disponível/sem URLs — caindo para SRA (prefetch+fasterq).[/yellow]")

        # Se ainda não temos FASTQs, seguir com SRA
        if not r1_gz.exists() and not r2_gz.exists():
            # Já existe .sra local?
            sra_target = next(Path("raw").glob(f"**/{acc}.sra"), None)
            if sra_target and sra_target.exists():
                console.print(f"{acc}: .sra → [bold]SKIP (cache)[/bold]")
            else:
                # 1) PREFETCH com retries
                retries = int(cfg_global["general"].get("prefetch_retries", 2))
                backoff = 30  # segundos base
                last_err = None
                for attempt in range(retries + 1):
                    try:
                        prefetch_with_progress(acc, Path("raw").resolve())
                        last_err = None
                        break
                    except sp.CalledProcessError as e:
                        last_err = e
                        if attempt < retries:
                            wait = backoff * (attempt + 1)
                            console.print(f"[orange3]{acc}: prefetch falhou (tentativa {attempt+1}/{retries}). "
                                          f"Repetindo em {wait}s…[/orange3]")
                            time.sleep(wait)
                        else:
                            console.print(f"[red]{acc}: prefetch falhou após {retries+1} tentativas.[/red]")
                            # Se ena_fallback estiver ativo, tentaremos ENA logo abaixo
                if last_err is not None and not cfg_global["general"].get("ena_fallback", False):
                    raise last_err  # sem fallback, propaga o erro

            # 2) CONVERSÃO (se ainda precisamos)
            if not r1_gz.exists() and not r2_gz.exists():
                sra_path = _find_local_sra(acc, Path("raw"))
                if sra_path:
                    ok = fasterq_with_progress(
                        str(sra_path), acc=acc, outdir=Path("fastq").resolve(),
                        threads=threads, tmp_dir=tmp_dir,
                        stall_warn_min=cfg_global["general"].get("stall_warn_min", 10),
                        stall_fail_min=cfg_global["general"].get("stall_fail_min", 45),
                    )
                else:
                    ok = fasterq_with_progress(
                        acc, acc=acc, outdir=Path("fastq").resolve(),
                        threads=threads, tmp_dir=tmp_dir,
                        stall_warn_min=cfg_global["general"].get("stall_warn_min", 10),
                        stall_fail_min=cfg_global["general"].get("stall_fail_min", 45),
                    )

                if not ok:
                    # 3) Fallback ENA se habilitado
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

        # 4) COMPACTAÇÃO (se veio do ENA já é .gz e dá SKIP)
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

# === helpers p/ FastQC/pareado ===
def _fastqc_base(fq: Path) -> str:
    name = fq.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fa.gz", ".fasta.gz", ".fa", ".fasta"):
        if name.endswith(suf):
            return name[:-len(suf)]
    return name  # fallback

def fastqc_outputs_exist(fq: Path, outdir: Path = Path("qc")) -> bool:
    base = _fastqc_base(fq)
    zipf = outdir / f"{base}_fastqc.zip"
    html = outdir / f"{base}_fastqc.html"
    return zipf.exists() and html.exists()

def _find_r2_for_r1(r1: Path) -> Optional[Path]:
    name = r1.name
    candidates = [
        re.sub(r'_1\.ds\.fastq\.gz$', '_2.ds.fastq.gz', name, flags=re.I),
        re.sub(r'_1\.fastq\.gz$',    '_2.fastq.gz',    name, flags=re.I),
        re.sub(r'_1\.fq\.gz$',       '_2.fq.gz',       name, flags=re.I),
        re.sub(r'_1\.fastq$',        '_2.fastq',       name, flags=re.I),
        re.sub(r'_1\.fq$',           '_2.fq',          name, flags=re.I),
        re.sub(r'_1\.ds\.fq\.gz$',   '_2.ds.fq.gz',    name, flags=re.I),
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

    # 1) FASTQC (somente faltantes)
    to_run = [f for f in gz if not fastqc_outputs_exist(f)]
    if to_run:
        miss = ", ".join(_fastqc_base(f) for f in to_run)
        console.print(f"[yellow]FastQC faltando para:[/yellow] {miss}", style="dim")
        run(["fastqc", *map(str,to_run), "-o","qc"])
    else:
        console.print("FastQC → [bold]SKIP (todos presentes)[/bold]")

    # 2) MULTIQC (só se necessário)
    def _multiqc_uptodate(outdir: Path = Path("qc")) -> bool:
        report = outdir / "multiqc_report.html"
        if not report.exists():
            # aceita um report antigo com sufixo e cria um alias fixo
            alt = sorted(outdir.glob("multiqc_report*.html"))
            if not alt:
                return False
            latest = max(alt, key=lambda p: p.stat().st_mtime)
            try:
                report.write_bytes(latest.read_bytes())
            except Exception:
                return False
        # todos os outputs FastQC que deveriam existir
        fq_outs = []
        for f in gz:
            base = _fastqc_base(f)
            fq_outs += [(outdir/f"{base}_fastqc.zip"), (outdir/f"{base}_fastqc.html")]
        fq_outs = [p for p in fq_outs if p.exists()]
        if not fq_outs:
            return False
        latest_fq = max(p.stat().st_mtime for p in fq_outs)
        return report.stat().st_mtime >= latest_fq

    if not _multiqc_uptodate():
        run(["multiqc","qc","-o","qc","--force"])
    else:
        console.print("MultiQC → [bold]SKIP (atual)[/bold]")

    # 3) Trimming (apenas short reads)
    if read_type != "short":
        console.print("Leituras longas: pulo trimming por padrão.", style="dim")
        return

    fwd, rev = adapters["fwd"], adapters["rev"]

    for r1 in sorted(list(fq_dir.glob("*_1*.fastq.gz")) + list(fq_dir.glob("*_1*.fq.gz")) +
                     list(fq_dir.glob("*_1*.fastq"))    + list(fq_dir.glob("*_1*.fq"))):
        # prefixo até o _1 (mantém ERR..., NA..., etc.)
        core = re.sub(r'_1(?:\.ds)?\.(?:fastq|fq)(?:\.gz)?$', '', r1.name, flags=re.IGNORECASE)
        out1 = Path("trimmed")/f"{core}_1.trim.fq.gz"

        r2 = _find_r2_for_r1(r1)
        if r2:
            out2 = Path("trimmed")/f"{core}_2.trim.fq.gz"
            if out1.exists() and out2.exists():
                console.print(f"{core}: trimming → [bold]SKIP (cache)[/bold]")
                continue
            run(["cutadapt","-j",str(threads),"-q","20,20","-m","30",
                 "-a",fwd,"-A",rev,"-o",str(out1),"-p",str(out2),str(r1),str(r2)])
            print_meta(f"FASTQs pós-trimming ({core})", [out1, out2])
        else:
            if out1.exists():
                console.print(f"{core}: trimming (single) → [bold]SKIP (cache)[/bold]")
                continue
            run(["cutadapt","-j",str(threads),"-q","20","-m","30","-a",fwd,"-o",str(out1),str(r1)])
            print_meta(f"FASTQ pós-trimming ({core})", [out1])

    # FastQC extra nos trimmed faltantes
    trimmed = sorted(Path("trimmed").glob("*.fq.gz"))
    to_run_trim = [f for f in trimmed if not fastqc_outputs_exist(f)]
    if to_run_trim:
        run(["fastqc", *map(str,to_run_trim), "-o","qc"])
        run(["multiqc","qc","-o","qc","--force"])

# =================== Alinhamento DNA ===================

def align_one_sample(
    r1: Path,
    threads: int,
    read_type: str,
    cleanup,
    expected_r2_name: Optional[str] = None,
    sample_id: Optional[str] = None,
):
    """
    Alinha uma amostra short/long read com uso de RAM controlado e progresso (ETA) no BWA/BWA-MEM2.
    - Usa -K (bwa/bwa-mem2) e sort -m/-@ ajustáveis (general.sort_mem_mb, general.bwa_batch_k).
    - Respeita 'general.aln_threads' (se existir).
    - Detecta ausência de índice do bwa-mem2 e faz fallback para 'bwa mem'.
    - Mostra progresso/ETA lendo o stderr do BWA e, se possível, o Total Sequences do FastQC.
    """

    import re, shutil, subprocess as sp, zipfile

    # ---------- helpers internos ----------
    def _bwa_index_ready(prefix: Path) -> bool:
        for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            if not Path(str(prefix) + ext).exists():
                return False
        return True

    def _mem2_index_present(prefix: Path) -> bool:
        # sinais do índice do bwa-mem2
        return Path(str(prefix) + ".bwt.2bit.64").exists() or any(Path("refs").glob("reference.fa.*.0123"))

    def _fastqc_total_for_filename(target_name: str) -> Optional[int]:
        """
        Percorre qc/*.zip, abre *fastqc_data.txt* e procura por:
          Filename\t<target_name>
          Total Sequences\tN
        Retorna N se encontrar; caso contrário None.
        """
        qc_dir = Path("qc")
        if not qc_dir.exists():
            return None
        for z in qc_dir.glob("*.zip"):
            try:
                with zipfile.ZipFile(z) as zf:
                    # acha o arquivo fastqc_data.txt
                    member = None
                    for n in zf.namelist():
                        if n.endswith("/fastqc_data.txt"):
                            member = n
                            break
                    if member is None:
                        continue
                    with zf.open(member) as fh:
                        fname_ok = False
                        total = None
                        for raw in fh:
                            line = raw.decode("utf-8", "ignore").rstrip("\n")
                            if line.startswith("Filename"):
                                # "Filename\tERR3239334_1.trim.fq.gz" (ou nome original)
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
        """
        Tenta obter o total de reads via FastQC. Se existir R2, total = n1 + n2.
        Se não conseguir, retorna None (sem ETA, mas com taxa e contagem).
        """
        n1 = _fastqc_total_for_filename(r1_path.name)
        n2 = _fastqc_total_for_filename(r2_path.name) if r2_path is not None else None
        if n1 is None and r2_path is not None:
            # alguns setups só têm FastQC do R1; nesse caso, aproxima 2×
            n1 = _fastqc_total_for_filename(r1_path.name)
            if n1 is not None:
                return n1 * 2
        if n1 is None:
            return None
        return n1 + (n2 or 0)

    # ---------- nome canônico ----------
    sample = sample_id if sample_id else re.sub(
        r'(_R?1|\.1)(\.trim)?\.(fastq|fq)\.gz$', '', r1.name
    )

    out_sorted = Path("bam") / f"{sample}.sorted.bam"
    out_mkdup  = Path("bam") / f"{sample}.mkdup.bam"
    out_cram   = Path("bam") / f"{sample}.mkdup.cram"

    if out_mkdup.exists() or out_cram.exists():
        console.print(f"[{sample}] alinhado → [bold]SKIP (cache)[/bold]")
        return

    # ---------- localizar R2 ----------
    r2 = None
    cand = []
    if expected_r2_name:
        cand.append(r1.parent / expected_r2_name)
    name = r1.name
    for pat, rep in (
        (r"_1\.trim\.fq\.gz$", "_2.trim.fq.gz"),
        (r"_1\.fastq\.gz$",    "_2.fastq.gz"),
        (r"_1\.fq\.gz$",       "_2.fq.gz"),
        (r"R1\.fastq\.gz$",    "R2.fastq.gz"),
        (r"R1\.fq\.gz$",       "R2.fq.gz"),
        (r"\.1\.fastq\.gz$",   ".2.fastq.gz"),
        (r"\.1\.fq\.gz$",      ".2.fq.gz"),
    ):
        if re.search(pat, name):
            cand.append(r1.parent / re.sub(pat, rep, name))
    for c in cand:
        if c.exists():
            r2 = c
            break
    if r2 and r2.resolve() == r1.resolve():
        console.print(f"[yellow][{sample}] Atenção:[/yellow] R2==R1 detectado; prosseguindo como single-end.", style="italic")
        r2 = None

    # ---------- parâmetros ----------
    g = cfg_global.get("general", {})
    aln_threads = int(g.get("aln_threads", threads))
    aln_threads = max(1, aln_threads)

    # samtools sort: por padrão metade das threads do alinhamento (mín. 1)
    sort_threads = int(g.get("sort_threads", max(1, aln_threads // 2)))
    sort_threads = max(1, min(sort_threads, aln_threads))
    sort_mem_mb  = int(g.get("sort_mem_mb", 512))
    bwa_batch_k  = str(int(g.get("bwa_batch_k", 20_000_000)))

    aligner_cfg  = g.get("aligner", "bwa-mem2").lower()
    rg_field     = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    ref_prefix   = Path("refs/reference.fa")

    # escolher alinhador efetivo
    aligner_effective = aligner_cfg
    if read_type == "short":
        if aligner_cfg in ("bwa-mem2", "bwa_mem2", "mem2"):
            if not _mem2_index_present(ref_prefix):
                console.print(
                    f"[orange3][{sample}] Índice do bwa-mem2 ausente → fallback automático para 'bwa mem'. "
                    "Use aligner: bwa + bwa_prebuilt_url no YAML OU gere o índice do mem2 se tiver RAM.[/orange3]"
                )
                aligner_effective = "bwa"
            else:
                aligner_effective = "bwa-mem2"
        else:
            aligner_effective = "bwa"

    # sumário
    console.print(
        f"[bold cyan][{sample}] Alinhador:[/bold cyan] {('BWA' if aligner_effective=='bwa' else 'BWA-MEM2')}  "
        f"[bold cyan]aln_threads:[/bold cyan] {aln_threads}  "
        f"[bold cyan]-K:[/bold cyan] {bwa_batch_k}  "
        f"[bold cyan]sort -@:[/bold cyan] {sort_threads}  "
        f"[bold cyan]sort -m:[/bold cyan] {sort_mem_mb}M  "
        f"[bold cyan]R2:[/bold cyan] {('OK' if r2 else 'single-end')}"
    )

    # ---------- execução ----------
    if read_type == "short":
        # pré-cheque de índice BWA se for o caminho efetivo
        if aligner_effective == "bwa" and not _bwa_index_ready(ref_prefix):
            raise RuntimeError(
                "Índice BWA ausente/incompleto em refs/reference.fa.{amb,ann,bwt,pac,sa}. "
                "Se estiver usando o pacote do GDC, garanta que os symlinks apontem para os arquivos extraídos."
            )

        # comando base do alinhador
        if aligner_effective == "bwa":
            base_cmd = [
                "bwa", "mem",
                "-R", rg_field,
                "-t", str(aln_threads),
                "-K", bwa_batch_k,
                "refs/reference.fa",
            ]
        else:
            base_cmd = [
                "bwa-mem2", "mem",
                "-R", rg_field,
                "-t", str(aln_threads),
                "-K", bwa_batch_k,
                "-Y",  # CIGAR compatível com Picard (soft-clip em reads suplementares)
                "refs/reference.fa",
            ]

        # stdbuf para flush de stderr/stdout (progresso em tempo real)
        if shutil.which("stdbuf"):
            base_cmd = ["stdbuf", "-oL", "-eL"] + base_cmd

        reads = [str(r1)] + ([str(r2)] if r2 else [])
        total_reads = _estimate_total_reads_for_alignment(r1, r2)

        cmd_list = [
            base_cmd + reads,
            ["samtools", "view", "-u", "-"],  # BAM não comprimido no pipe → menos CPU/RAM
            ["samtools", "sort", "-@", str(sort_threads), "-m", f"{sort_mem_mb}M",
             "-o", str(out_sorted)],
        ]

        try:
            run_long_stream_pipeline(
                cmd_list,
                label=f"[{sample}] {('BWA' if aligner_effective=='bwa' else 'BWA-MEM2')} → sort",
                heartbeat_sec=60,
                progress_spec={"total_reads": total_reads, "print_every": 15},
            )
        except sp.CalledProcessError as e:
            # fallback extra se bwa-mem2 falhar em runtime
            if aligner_effective == "bwa-mem2":
                console.print(
                    f"[orange3][{sample}] bwa-mem2 falhou (status {e.returncode}). "
                    "Tentando fallback com 'bwa mem' e -K menor…[/orange3]"
                )
                small_k = str(max(5_000_000, int(bwa_batch_k) // 2))
                base_cmd_fb = [
                    "bwa", "mem",
                    "-R", rg_field,
                    "-t", str(aln_threads),
                    "-K", small_k,
                    "refs/reference.fa",
                ]
                if shutil.which("stdbuf"):
                    base_cmd_fb = ["stdbuf", "-oL", "-eL"] + base_cmd_fb
                cmd_list_fb = [
                    base_cmd_fb + reads,
                    ["samtools", "view", "-u", "-"],
                    ["samtools", "sort", "-@", str(sort_threads), "-m", f"{sort_mem_mb}M",
                     "-o", str(out_sorted)],
                ]
                run_long_stream_pipeline(
                    cmd_list_fb,
                    label=f"[{sample}] BWA (fallback) → sort",
                    heartbeat_sec=60,
                    progress_spec={"total_reads": total_reads, "print_every": 15},
                )
            else:
                raise

        run(["samtools", "index", str(out_sorted)])

        run([
            "picard", "MarkDuplicates",
            "-I", str(out_sorted),
            "-O", str(out_mkdup),
            "-M", str(out_sorted).replace(".sorted.bam", ".mkdup.metrics"),
            "--VALIDATION_STRINGENCY", "LENIENT",
        ])
        run(["samtools", "index", str(out_mkdup)])

        if cleanup.get("remove_sorted_bam", True) and out_sorted.exists():
            out_sorted.unlink(missing_ok=True)

    else:
        # long reads (ONT, por ex.) — sem monitor específico
        cmd_list = [
            ["minimap2", "-ax", "map-ont", "-t", str(aln_threads), "refs/reference.fa", str(r1)],
            ["samtools", "sort", "-@", str(sort_threads), "-m", f"{sort_mem_mb}M", "-o", str(out_sorted)],
        ]
        run_long_stream_pipeline(
            cmd_list,
            label=f"[{sample}] minimap2 → sort",
            heartbeat_sec=60,
            progress_spec=None,  # (não extraímos ETA do minimap2 aqui)
        )
        run(["samtools", "index", str(out_sorted)])

def align_dna_for_all(dna_samples, threads, default_read_type, cleanup, use_ds):
    """
    Alinha todas as amostras DNA de forma idempotente:
      - Pula a etapa inteira se todas as saídas mkdup já existem.
      - Pula cada amostra individualmente se mkdup já existe.
      - Escolhe um único R1 por run (trimmed > fastq_ds > fastq).
    """
    from pathlib import Path

    def _has_mkdup(sid: str) -> bool:
        b = Path("bam")
        return (b / f"{sid}.mkdup.bam").exists() or (b / f"{sid}.mkdup.cram").exists()

    # --- fast-exit da etapa 6: tudo pronto? ---
    sample_ids = [s["id"] for s in dna_samples if "id" in s]
    unique_ids = sorted(set(sample_ids))
    if unique_ids and all(_has_mkdup(sid) for sid in unique_ids):
        console.print("Etapa 6 (alinhamento) → [bold]SKIP[/bold] (todas as amostras já têm mkdup)", style="dim")
        outs = sorted(Path("bam").glob("*.mkdup.bam")) + sorted(Path("bam").glob("*.mkdup.cram"))
        if outs:
            print_meta("Alinhados (mkdup)", outs)
        return

    # --- helper para escolher UM R1 por run ---
    def _pick_r1_for_run(acc: str) -> Optional[Path]:
        # 1) trimmed
        trimmed = sorted(Path("trimmed").glob(f"{acc}*_1.trim.fq.gz"))
        if trimmed:
            return trimmed[0]
        # 2) fastq_ds (se habilitado)
        if use_ds:
            ds = sorted(Path("fastq_ds").glob(f"{acc}*_1.fastq.gz")) + \
                 sorted(Path("fastq_ds").glob(f"{acc}*_1.fq.gz"))
            if ds:
                return ds[0]
        # 3) fastq
        fq = sorted(Path("fastq").glob(f"{acc}*_1.fastq.gz")) + \
             sorted(Path("fastq").glob(f"{acc}*_1.fq.gz"))
        return fq[0] if fq else None

    # --- loop por amostra ---
    for s in dna_samples:
        rtype = s.get("read_type", default_read_type)
        sid = s["id"]

        # skip por amostra se mkdup já existe
        if _has_mkdup(sid):
            console.print(f"[{sid}] alinhado → [bold]SKIP (cache mkdup)[/bold]")
            continue

        if s.get("source") == "sra":
            for acc in s.get("sra_ids", []):
                r1 = _pick_r1_for_run(acc)
                if not r1:
                    console.print(f"[yellow][{sid}] Nenhum R1 encontrado para {acc} (trimmed/fastq_ds/fastq).[/yellow]")
                    continue
                # processa só UM R1 por run
                align_one_sample(r1, threads, rtype, cleanup, expected_r2_name=None, sample_id=sid)
        else:
            # modo "direto": tenta usar trimmed, senão o original
            fastq1 = Path(s.get("fastq1", ""))
            if fastq1:
                r1_trim = Path("trimmed") / fastq1.name.replace(".fastq.gz", ".trim.fq.gz")
                r1 = r1_trim if r1_trim.exists() else (Path("fastq_ds") / fastq1.name if use_ds else Path("fastq") / fastq1.name)
                expected_r2 = Path(s.get("fastq2", ""))
                align_one_sample(
                    r1, threads, rtype, cleanup,
                    expected_r2_name=(expected_r2.name if expected_r2 else None),
                    sample_id=sid
                )
            else:
                console.print(f"[yellow][{sid}] Sem 'fastq1' e sem 'sra_ids'; amostra ignorada.[/yellow]")

    # --- resumo final da etapa ---
    outs = sorted(Path("bam").glob("*.mkdup.bam")) + sorted(Path("bam").glob("*.mkdup.cram"))
    if outs:
        print_meta("Alinhados (mkdup)", outs)

# ============== CRAM + Cobertura (mosdepth) ==============

def to_cram_and_coverage(use_cram: bool, threads: int):
    """
    Converte BAM→CRAM (se configurado), garante índices e roda mosdepth.
    - Se a entrada do mosdepth for CRAM, passa --fasta refs/reference.fa.
    - Cria .crai/.bai se estiverem faltando.
    - Pula mosdepth se outputs já existem (cache).
    """
    from pathlib import Path

    bdir = Path("bam")
    ref_fa = Path("refs/reference.fa")
    if not bdir.exists():
        console.print("[yellow]Diretório 'bam' não existe — nada a fazer na etapa 7.[/yellow]")
        return

    # detecta amostras a partir dos mkdup
    mkdups = sorted(bdir.glob("*.mkdup.bam")) + sorted(bdir.glob("*.mkdup.cram"))
    if not mkdups:
        console.print("[yellow]Nenhum *.mkdup.(bam|cram) encontrado — nada a fazer.[/yellow]")
        return

    # helper: checa se mosdepth já foi feito
    def _mosdepth_done(prefix: Path) -> bool:
        # arquivos característicos do mosdepth
        gdist = prefix.with_suffix(".mosdepth.global.dist.txt")
        summ  = prefix.with_suffix(".mosdepth.summary.txt")
        # se você usa -n, não há per-base; estes 2 bastam p/ verificar cache
        return gdist.exists() and summ.exists()

    for mk in mkdups:
        sid = mk.name.split(".")[0]  # antes de .mkdup.*
        bam  = bdir / f"{sid}.mkdup.bam"
        cram = bdir / f"{sid}.mkdup.cram"
        cov_prefix = bdir / sid

        # escolhe entrada
        inp = None
        if use_cram:
            # se ainda não há CRAM mas há BAM, converte
            if not cram.exists() and bam.exists():
                if not ref_fa.exists():
                    raise RuntimeError("refs/reference.fa ausente — necessário para converter BAM→CRAM.")
                console.print(f"[{sid}] Convertendo BAM→CRAM…")
                run(["samtools","view","-C","-T",str(ref_fa),"-@",str(max(1, threads//2)),
                     "-o",str(cram), str(bam)])
            inp = cram if cram.exists() else (bam if bam.exists() else None)
        else:
            inp = bam if bam.exists() else (cram if cram.exists() else None)

        if not inp:
            console.print(f"[yellow][{sid}] Nenhum mkdup BAM/CRAM encontrado. Pulando.[/yellow]")
            continue

        # garante índice do arquivo de entrada
        if inp.suffix == ".cram":
            crai = Path(str(inp) + ".crai")
            if not crai.exists():
                run(["samtools","index","-@",str(max(1, threads//2)), str(inp)])
        else:
            bai = Path(str(inp) + ".bai")
            if not bai.exists():
                run(["samtools","index","-@",str(max(1, threads//2)), str(inp)])

        # mosdepth: usa --fasta apenas para CRAM
        if _mosdepth_done(cov_prefix):
            console.print(f"[{sid}] mosdepth → [bold]SKIP[/bold] (cache)")
            continue

        cmd = ["mosdepth","-t",str(threads)]
        if inp.suffix == ".cram":
            if not ref_fa.exists():
                raise RuntimeError("refs/reference.fa ausente — necessário para mosdepth com CRAM.")
            cmd += ["--fasta", str(ref_fa)]
        cmd += [str(cov_prefix), str(inp)]

        console.print(f"[{sid}] mosdepth em {inp.name}…")
        run(cmd)

    # (opcional) resumo do que foi gerado
    outs = sorted(Path("bam").glob("*.mosdepth.summary.txt"))
    if outs:
        print_meta("Cobertura (mosdepth)", outs)

# =================== Variantes / VEP ===================

def _haplotypecaller_one(
    bam_or_cram: Path, sample: str, out_gvcf: Path, intervals_file: Path|None,
    java_mem_gb: int, nat_threads: int, pcr_free: bool, extra_args: list[str]|None=None
):
    out_gvcf.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp_vcfgz_path(out_gvcf)
    if tmp.exists(): tmp.unlink()

    java_opts = f"-Xmx{int(java_mem_gb)}g -Dsamjdk.compression_level=2"
    cmd = [
        "gatk","--java-options",java_opts,"HaplotypeCaller",
        "-R","refs/reference.fa","-I",str(bam_or_cram),
        "-O",str(tmp), "-ERC","GVCF",
        "--native-pair-hmm-threads",str(max(1,int(nat_threads)))
    ]
    if pcr_free:
        cmd += ["--pcr-indel-model","NONE"]
    if intervals_file:
        cmd += ["-L", str(intervals_file)]
    if extra_args:
        cmd += list(map(str, extra_args))

    console.print(f"[{sample}] HC → mem={java_mem_gb}g, pair-hmm={nat_threads}, "
                  f"{'PCR-free' if pcr_free else 'PCR-default'} "
                  f"{'(intervals)' if intervals_file else '(genome inteiro)'}")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(max(1, int(nat_threads)))

    console.print(f"[{sample}] HC → mem={java_mem_gb}g, pair-hmm={nat_threads}, "
                f"{'PCR-free' if pcr_free else 'PCR-default'} "
                f"{'(intervals)' if intervals_file else '(genome inteiro)'}")

    run(cmd)
    _atomic_rename(tmp, out_gvcf)
    _ensure_tabix(out_gvcf)

def _gather_vcfs(inputs: list[Path], out: Path):
    # Gathers só se necessário
    if out.exists() and _is_newer(out, *inputs):
        return out
    tmp = _tmp_vcfgz_path(out)
    if tmp.exists(): tmp.unlink()
    args = []
    for f in sorted(inputs):
        args += ["-I", str(f)]
    run(["gatk","GatherVcfs", *args, "-O", str(tmp)])
    _atomic_rename(tmp, out)
    _ensure_tabix(out)
    return out

def _concat_vcfs(inputs: list[Path], out: Path):
    """
    Concatena shards VCF (ordenados por cromossomo/pos) com bcftools concat.
    Idempotente: só refaz se 'out' estiver faltando ou mais velho que os inputs.
    """
    if out.exists() and _is_newer(out, *inputs):
        _ensure_tabix(out)
        return out
    tmp = _tmp_vcfgz_path(out)
    if tmp.exists(): tmp.unlink()
    run(["bcftools", "concat", "-Oz", "-o", str(tmp), *map(str, inputs)])
    _atomic_rename(tmp, out)
    _ensure_tabix(out)
    return out


def _bcftools_call_one(
    bam_or_cram: Path, sample: str, out_vcf: Path, intervals_file: Path | None,
    min_baseq: int, min_mapq: int, threads: int,
    split_multiallelic: bool = True, emit_all_sites: bool = False,
    extra_mpileup_args: list[str] | None = None,
    extra_call_args: list[str] | None = None,
):
    """
    mpileup (C) -> call (C) -> norm (C) -> BGZF VCF + .tbi
    Mantém idempotência e compatibilidade de saídas com as fases 9–13.
    """
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp_vcfgz_path(out_vcf)
    if tmp.exists(): tmp.unlink()

    reg = f"-R {shlex.quote(str(intervals_file))}" if intervals_file else ""
    emit_flag = "-A" if emit_all_sites else "-v"
    split_flag = "-m +any" if split_multiallelic else ""
    xtra_mp = " ".join(map(shlex.quote, extra_mpileup_args or []))
    xtra_call = " ".join(map(shlex.quote, extra_call_args or []))

    console.print(
        f"[{sample}] BCFtools → mpileup/call (q={min_mapq}, Q={min_baseq}) "
        f"{'(intervals)' if intervals_file else '(genoma inteiro)'}"
    )

    # Obs.: --threads atua na compressão/IO; mpileup/call em si são single-threaded.
    script = f"""
set -eo pipefail
bcftools mpileup -Ou -f refs/reference.fa -q {int(min_mapq)} -Q {int(min_baseq)} \
  -a FORMAT/DP,FORMAT/AD,INFO/AD {reg} {xtra_mp} {shlex.quote(str(bam_or_cram))} \
| bcftools call -m {emit_flag} -Ou {xtra_call} \
| bcftools norm -f refs/reference.fa {split_flag} -Ou \
| bcftools view --threads {int(max(1,threads))} -Oz -o {shlex.quote(str(tmp))}
"""
    run(["bash","-lc", script])
    _atomic_rename(tmp, out_vcf)
    _ensure_tabix(out_vcf)

def _concat_vcfs_bcftools(inputs: list[Path], out: Path, io_threads: int = 2):
    """Concatena shards com bcftools concat (idempotente)."""
    if out.exists() and _is_newer(out, *inputs):
        _ensure_tabix(out); return out
    tmp = _tmp_vcfgz_path(out)
    if tmp.exists(): tmp.unlink()
    cmd = ["bcftools","concat","-Oz","-o",str(tmp), *map(str, sorted(inputs))]
    if io_threads > 1: cmd += ["--threads", str(io_threads)]
    run(cmd)
    _atomic_rename(tmp, out)
    _ensure_tabix(out)
    return out

def _bcftools_call_one_verbose(
    bam_or_cram: Path, sample: str, out_vcf: Path, intervals_file: Path | None,
    mapq: int, baseq: int, io_threads: int, emit_variants_only: bool, heartbeat_sec: int
):
    """Chama variantes com bcftools mpileup→call emitindo banners/heartbeats/throughput."""
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp_vcfgz_path(out_vcf)
    if tmp.exists(): tmp.unlink()

    # escopo alvo (para feedback)
    if intervals_file:
        target_bp = _bp_in_bed(intervals_file)
        scope = f"intervals ({_human_bp(target_bp)})"
    else:
        # soma do .fai
        target_bp = sum(ln for _, ln in _read_fai(Path("refs/reference.fa.fai")))
        scope = f"genoma inteiro ({_human_bp(target_bp)})"

    label = f"[{sample}] bcftools mpileup→call " + ("(intervals)" if intervals_file else "(genoma inteiro)")
    console.print(Panel.fit(
        f"{label}\nAlvo: {_human_bp(target_bp)}  •  MAPQ≥{mapq}  •  BQ≥{baseq}  •  IO-threads={io_threads}",
        border_style="magenta"
    ))

    # comandos
    mp = ["bcftools","mpileup", "-f","refs/reference.fa",
          "-q",str(mapq), "-Q",str(baseq),
          "-Ou", "-a","DP,AD", str(bam_or_cram)]
    if intervals_file:
        mp += ["-R", str(intervals_file)]
    if io_threads > 1:
        mp += ["--threads", str(io_threads)]

    call = ["bcftools","call","-m", "-Oz", "-o", str(tmp)]
    if emit_variants_only:
        call.insert(3, "-v")
    if io_threads > 1:
        call += ["--threads", str(io_threads)]

    import time as _time
    t0 = _time.time()
    run_long_stream_pipeline([mp, call], label=label, heartbeat_sec=max(5, int(heartbeat_sec)))
    dt = max(1, int(_time.time() - t0))
    speed_mbps = (target_bp / dt) / 1e6
    console.print(f"[bold cyan]{sample}[/bold cyan] terminou {scope} em {dt//60}m{dt%60:02d}s  (~{speed_mbps:.1f} Mbp/s)")
    _ensure_tabix(tmp)
    _atomic_rename(tmp, out_vcf)

def call_variants(samples, threads, mem_gb):
    """
    Gera variantes por amostra com GATK (gVCF->GenotypeGVCFs) ou BCFtools (mpileup->call->norm).
    - Se params.variant_caller == "gatk": usa HaplotypeCaller (com scatter opcional) + GenotypeGVCFs.
    - Se params.variant_caller == "bcftools": paraleliza shards e concatena em vcf/<sample>.vcf.gz.
    - Idempotente: pula quando outputs existem e estão atualizados (inclui .tbi).

    Parâmetros YAML (params):
      variant_caller: "gatk" | "bcftools"
      filter_problematic_contigs: bool # filtra cromossomos virais/problemáticos (default: true)
      # BCFtools
      bcf_scatter_parts: int          # nº de shards
      bcf_max_parallel: int           # nº de shards em paralelo
      bcf_threads_io: int             # --threads de htslib por processo (mpileup/norm)
      bcf_mapq: int                   # -q em mpileup (MAPQ)
      bcf_baseq: int                  # -Q em mpileup (BaseQ)
      bcf_max_depth: int              # --max-depth em mpileup
      bcf_heartbeat_sec: int          # intervalo do heartbeat (s) por shard (default 30)
      bcf_preview_max_contigs: int    # nº máx de contigs no preview dos BEDs (default 4)
      bcf_ignore_ref_mismatch: bool   # adiciona --ignore-RG para robustez (default: true)
      bcf_skip_indels_on_mismatch: bool # adiciona --skip-indels para robustez (default: true)
    """
    from pathlib import Path
    import os as _os, time, shutil, subprocess as sp, shlex

    Path("vcf").mkdir(exist_ok=True)
    Path("vcf/shards").mkdir(parents=True, exist_ok=True)
    Path("logs").mkdir(exist_ok=True)

    p = cfg_global.get("params", {})
    g = cfg_global.get("general", {})

    variant_caller = (p.get("variant_caller") or "gatk").lower().strip()
    if variant_caller not in ("gatk", "bcftools"):
        console.print(f"[yellow]variant_caller='{variant_caller}' desconhecido — usando GATK.[/yellow]")
        variant_caller = "gatk"

    # === Contigs/intervalos ===
    keep_alt_decoy = bool(p.get("keep_alt_decoy", True))
    filter_problematic = bool(p.get("filter_problematic_contigs", True))  # Novo parâmetro
    
    if g.get("limit_to_canonical", False):
        keep_alt_decoy = False
        
    contigs = _read_fai(Path("refs/reference.fa.fai"))
    
    # Aplica filtros em sequência
    if keep_alt_decoy:
        contigs_used = contigs
    else:
        contigs_used = _canonical_subset(contigs)
    
    # Filtra cromossomos problemáticos se habilitado
    if filter_problematic:
        contigs_used = _filter_problematic_contigs(contigs_used)
    intervals_dir = Path("intervals"); intervals_dir.mkdir(exist_ok=True)

    # === Entradas (mkdup: prefer BAM ou CRAM) ===
    prefer_bam_for_calling = bool(p.get("prefer_bam_for_calling", True))
    if prefer_bam_for_calling:
        aln = sorted(Path("bam").glob("*.mkdup.bam")) + sorted(Path("bam").glob("*.mkdup.cram"))
    else:
        aln = sorted(Path("bam").glob("*.mkdup.cram")) + sorted(Path("bam").glob("*.mkdup.bam"))

    if not aln:
        console.print("[yellow]Aviso:[/yellow] nenhum BAM/CRAM encontrado para chamada de variantes.", style="italic")
        return []

    # ======== Helpers ========
    def _mk_panel(msg, style="bright_cyan"):
        return Panel.fit(msg, border_style=style)

    def _fai_order_dict(fai_path="refs/reference.fa.fai"):
        od = {}
        with open(fai_path) as fh:
            for i, ln in enumerate(fh):
                od[ln.split("\t", 1)[0]] = i
        return od

    def _first_bed_contig(bed_path: Path) -> str:
        with open(bed_path) as fh:
            for ln in fh:
                if ln.strip():
                    return ln.split()[0]
        return ""

    def _read_bed_stats(bed: Path):
        """Retorna (contigs_em_ordem, num_linhas, total_bp)."""
        contigs_seen = []
        total_bp = 0
        n = 0
        with open(bed) as fh:
            for ln in fh:
                if not ln.strip():
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[0], parts[1], parts[2]
                try:
                    bp = max(0, int(end) - int(start))
                except Exception:
                    bp = 0
                total_bp += bp
                n += 1
                if (not contigs_seen) or contigs_seen[-1] != chrom:
                    contigs_seen.append(chrom)
        return contigs_seen, n, total_bp

    def _fmt_bp(n):
        # formato amigável para pares de bases
        units = [("Gbp", 10**9), ("Mbp", 10**6), ("kbp", 10**3)]
        for u, s in units:
            if n >= s:
                return f"{n/s:.2f} {u}"
        return f"{n} bp"

    # ---------- GATK (mesmo fluxo já usado) ----------
    def _gatk_call_one(
        bam_or_cram: Path, sample: str, out_gvcf: Path, intervals_file: Path|None,
        java_mem_gb: int, nat_threads: int, pcr_free: bool, extra_args: list[str]|None=None
    ):
        return _haplotypecaller_one(
            bam_or_cram, sample, out_gvcf, intervals_file,
            java_mem_gb=java_mem_gb, nat_threads=nat_threads, pcr_free=pcr_free, extra_args=extra_args
        )

    # ---------- bcftools: cmd + labels ----------
    def _bcf_make_cmd(sample: str, bam: Path, bed: Path, out_vcf: Path,
                      io_threads: int, mapq: int, baseq: int, max_depth: int) -> list[str]:
        tmp = out_vcf.with_suffix(".tmp.vcf.gz")
        
        # Parâmetros para lidar com reference mismatches
        p = cfg_global.get("params", {})
        skip_indels = bool(p.get("bcf_skip_indels_on_mismatch", True))
        ignore_ref_mismatch = bool(p.get("bcf_ignore_ref_mismatch", True))
        
        # Flags adicionais para robustez
        extra_mpileup_flags = ""
        if ignore_ref_mismatch:
            extra_mpileup_flags += " --ignore-RG"  # Ignora read groups problemáticos
        if skip_indels:
            extra_mpileup_flags += " --skip-indels"  # Pula indels em regiões problemáticas
            
        sh = (
            "set -euo pipefail; "
            f"bcftools mpileup -f refs/reference.fa -Ou --threads {io_threads} "
            f"-q {mapq} -Q {baseq} --max-depth {max_depth} "
            f"-a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD "
            f"{extra_mpileup_flags} "
            f"-R {shlex.quote(str(bed))} {shlex.quote(str(bam))} "
            f"| bcftools call -m -v -Ou "
            f"| bcftools norm -f refs/reference.fa --threads {io_threads} --multiallelics -both "
            f"-Oz -o {shlex.quote(str(tmp))}"
        )
        return ["bash","-lc", sh]

    def _bcf_log_command_and_bed(sample: str, bed: Path, cmd: list[str]):
        """Mostra comando bcftools e cromossomos do BED para debug."""
        console.print(f"[dim]💻 Comando BCFtools:[/dim]")
        console.print(f"[dim]> {' '.join(cmd)}[/dim]")
        
        # Mostra cromossomos no BED
        try:
            with open(bed) as f:
                chroms = set()
                for line in f:
                    if line.strip():
                        chrom = line.split('\t')[0]
                        chroms.add(chrom)
                if chroms:
                    console.print(f"[dim]📋 Cromossomos no shard: {', '.join(sorted(chroms))}[/dim]")
        except Exception:
            pass

    def _bcf_label(sample: str, io_threads: int, mapq: int, baseq: int, max_depth: int, shard_tag: str=""):
        return (f"[{sample}{('|' + shard_tag) if shard_tag else ''}] bcftools: mpileup→call→norm "
                f"(threads={io_threads}, MAPQ≥{mapq}, BaseQ≥{baseq}, max-depth={max_depth})")

    # ---------- helper serial (progresso por shard) ----------
    def _run_bcftools_shard_with_progress(sample: str, bam: Path, bed: Path, out_shard: Path,
                                          io_threads: int, mapq: int, baseq: int, max_depth: int,
                                          heartbeat_sec: int = 30):
        """Executa 1 shard (sem paralelismo externo) monitorando crescimento do .tmp e emitindo heartbeat."""
        tmp = out_shard.with_suffix(".tmp.vcf.gz")
        tmp.unlink(missing_ok=True)

        shard_tag = out_shard.name.replace(".vcf.gz","")
        console.print(_mk_panel(_bcf_label(sample, io_threads, mapq, baseq, max_depth, shard_tag), "white"))

        cmd = _bcf_make_cmd(sample, bam, bed, out_shard, io_threads, mapq, baseq, max_depth)
        _bcf_log_command_and_bed(sample, bed, cmd)  # 🔧 Mostra comando e cromossomos
        log_path = Path("logs")/f"bcftools_{sample}_{shard_tag}.log"
        with open(log_path, "w") as log_fh:
            if shutil.which("stdbuf"):
                cmd = ["stdbuf","-oL","-eL"] + cmd
            p = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, text=True)

            start = time.time()
            last_ts = start
            last_sz = 0
            while True:
                rc = p.poll()
                now = time.time()
                if now - last_ts >= max(3, int(heartbeat_sec)):
                    sz = tmp.stat().st_size if tmp.exists() else 0
                    dsz = max(0, sz - last_sz)
                    console.print(f"[{sample}|{shard_tag}] … {int((now-start)//60)}m{int((now-start)%60):02d}s "
                                  f"• out {sizeof_fmt(sz)} • Δ{sizeof_fmt(dsz)}  • log={log_path.name}", highlight=False)
                    last_ts = now
                    last_sz = sz
                if rc is not None:
                    break
                time.sleep(0.5)

        if p.returncode != 0:
            # mostra tail do log p/ debug
            try:
                tail = sp.run(["bash","-lc", f"tail -n 60 {shlex.quote(str(log_path))}"],
                              capture_output=True, text=True, check=True).stdout
            except Exception:
                tail = ""
            raise sp.CalledProcessError(p.returncode, p.args, output=tail)

        # finalize: renomeia tmp, indexa e reporta
        if tmp.exists():
            _atomic_rename(tmp, out_shard)
        run(["tabix","-p","vcf",str(out_shard)])
        size_mb = (out_shard.stat().st_size / (1024*1024)) if out_shard.exists() else 0.0
        dt = int(time.time() - start)
        # imprime resumo do log (linhas de estatística do bcftools, se houver)
        try:
            lines_line = sp.run(
                ["bash","-lc", f"grep -E '^[Ll]ines\\s' -m1 {shlex.quote(str(log_path))} || true"],
                capture_output=True, text=True, check=True
            ).stdout.strip()
        except Exception:
            lines_line = ""
        console.print(f"[{sample}] shard pronto → {out_shard.name}  •  {size_mb:.1f}MB  •  {dt//60}m{dt%60:02d}s"
                      + (f"\n{lines_line}" if lines_line else ""))

    # ======== Caminho GATK ========
    if variant_caller == "gatk":
        scatter_parts   = max(1, int(p.get("hc_scatter_parts", 1)))
        hc_mem_gb       = int(p.get("hc_java_mem_gb", 24))
        _hc_threads_cfg = p.get("hc_threads_native", g.get("hc_threads_native", None))
        try:
            hc_nat_threads = int(_hc_threads_cfg) if _hc_threads_cfg is not None else 4
        except Exception:
            hc_nat_threads = 4
        try:
            if _os.cpu_count():
                hc_nat_threads = max(1, min(hc_nat_threads, _os.cpu_count()))
        except Exception:
            pass
        hc_pcr_free     = bool(p.get("hc_pcr_free", True))
        hc_extra        = list(p.get("hc_extra_args", [])) if isinstance(p.get("hc_extra_args", []), list) else []

        gg_mem_gb       = int(p.get("gg_java_mem_gb", 12))
        gg_extra        = list(p.get("gg_extra_args", [])) if isinstance(p.get("gg_extra_args", []), list) else []
        force_recall    = bool(p.get("force_recall", False))

        vcfs_final = []

        for bam in aln:
            sample = bam.name.replace(".mkdup.cram","").replace(".mkdup.bam","")
            gvcf = Path("vcf")/f"{sample}.g.vcf.gz"
            vcf  = Path("vcf")/f"{sample}.vcf.gz"

            # ---- gVCF ----
            need_gvcf = force_recall or (not gvcf.exists()) or (not _is_newer(gvcf, bam)) or (not Path(str(gvcf)+".tbi").exists())

            if not need_gvcf and scatter_parts == 1:
                console.print(f"[{sample}] HaplotypeCaller → [bold]SKIP (cache ok)[/bold]")
            else:
                if scatter_parts == 1:
                    intervals_file = None
                    if not keep_alt_decoy:
                        intervals_file = intervals_dir / "canonical.bed"
                        if not intervals_file.exists():
                            _write_intervals_file(contigs_used, intervals_file)
                    _gatk_call_one(
                        bam, sample, gvcf, intervals_file,
                        java_mem_gb=hc_mem_gb, nat_threads=hc_nat_threads,
                        pcr_free=hc_pcr_free, extra_args=hc_extra
                    )
                else:
                    parts = _split_balanced(contigs_used, scatter_parts)
                    shard_dir = Path(f"vcf/shards/{sample}"); shard_dir.mkdir(parents=True, exist_ok=True)
                    shard_gvcfs = []
                    for i, part in enumerate(parts, start=1):
                        intr = shard_dir / f"part_{i:02d}.bed"
                        _write_intervals_file(part, intr)
                        out_shard = shard_dir / f"{sample}.part_{i:02d}.g.vcf.gz"
                        if not (out_shard.exists() and _is_newer(out_shard, bam) and Path(str(out_shard)+".tbi").exists()):
                            _gatk_call_one(
                                bam, sample, out_shard, intr,
                                java_mem_gb=max(24, hc_mem_gb//max(1,scatter_parts//2)),
                                nat_threads=max(1, min(16, hc_nat_threads)),
                                pcr_free=hc_pcr_free, extra_args=hc_extra
                            )
                        else:
                            _ensure_tabix(out_shard)
                        shard_gvcfs.append(out_shard)
                    _gather_vcfs(shard_gvcfs, gvcf)

            # ---- GenotypeGVCFs → VCF ----
            need_vcf = force_recall or (not vcf.exists()) or (not _is_newer(vcf, gvcf)) or (not Path(str(vcf)+".tbi").exists())
            if need_vcf:
                tmp = _tmp_vcfgz_path(vcf)
                if tmp.exists(): tmp.unlink()
                java_opts = f"-Xmx{gg_mem_gb}g -Dsamjdk.compression_level=2"
                cmd = [
                    "gatk","--java-options",java_opts,"GenotypeGVCFs",
                    "-R","refs/reference.fa","-V",str(gvcf),"-O",str(tmp)
                ] + gg_extra
                console.print(f"[{sample}] GenotypeGVCFs → mem={gg_mem_gb}g")
                run(cmd)
                _atomic_rename(tmp, vcf)
                _ensure_tabix(vcf)
            else:
                console.print(f"[{sample}] GenotypeGVCFs → [bold]SKIP (cache ok)[/bold]")

            vcfs_final.append(vcf)
            print_meta(f"VCFs ({sample})", [gvcf, vcf])

        if vcfs_final:
            print_meta("VCFs gerados (por amostra)", vcfs_final)
        else:
            console.print("[yellow]Aviso:[/yellow] Nenhum VCF novo foi gerado (talvez já existissem).", style="italic")
        return vcfs_final

    # ======== Caminho BCFtools (paralelo em shards) ========
    bcf_parts        = max(1, int(p.get("bcf_scatter_parts", p.get("hc_scatter_parts", 1))))
    bcf_max_parallel = max(1, int(p.get("bcf_max_parallel", 1)))
    # threads de I/O (htslib) por processo — prudente para não saturar disco
    guess_io = max(1, threads // max(1, min(bcf_max_parallel, threads)))
    bcf_io_threads   = max(1, int(p.get("bcf_threads_io", min(2, guess_io))))
    bcf_mapq         = int(p.get("bcf_mapq", 20))
    bcf_baseq        = int(p.get("bcf_baseq", 20))
    bcf_max_depth    = int(p.get("bcf_max_depth", 250))
    bcf_heartbeat    = int(p.get("bcf_heartbeat_sec", 30))
    bcf_prev_maxc    = int(p.get("bcf_preview_max_contigs", 4))
    force_recall     = bool(p.get("force_recall", False))

    vcfs_final = []
    fai_order = _fai_order_dict()

    for bam in aln:
        sample = bam.name.replace(".mkdup.cram","").replace(".mkdup.bam","")
        final_vcf = Path("vcf")/f"{sample}.vcf.gz"
        final_tbi = Path(str(final_vcf)+".tbi")

        # ---- gerar shards (listas de BED) ----
        parts = _split_balanced(contigs_used, bcf_parts)
        shard_dir = Path(f"vcf/shards/{sample}")
        shard_dir.mkdir(parents=True, exist_ok=True)

        shards = []
        for i, part in enumerate(parts, start=1):
            bed = shard_dir / f"part_{i:02d}.bed"
            _write_intervals_file(part, bed)
            out_shard = shard_dir / f"{sample}.part_{i:02d}.vcf.gz"
            shards.append((i, bed, out_shard))

        # ---- PREVIEW dos BEDs (antes de iniciar) ----
        tbl = Table(title=f"Shards {sample} — preview dos BEDs", header_style="bold cyan")
        tbl.add_column("part", justify="right")
        tbl.add_column("contigs (preview)")
        tbl.add_column("#regiões", justify="right")
        tbl.add_column("total", justify="right")
        tbl.add_column("BED")
        for (i, bed, _) in shards:
            contigs_seen, n_regions, total_bp = _read_bed_stats(bed)
            preview = ", ".join(contigs_seen[:bcf_prev_maxc]) + (f" +{len(contigs_seen)-bcf_prev_maxc}" if len(contigs_seen) > bcf_prev_maxc else "")
            tbl.add_row(f"{i:02d}", preview if preview else "—", str(n_regions), _fmt_bp(total_bp), bed.name)
        console.print(tbl)

        # Se não for force_recall, descarta da fila os shards que já estão OK
        todo = []
        for i, bed, out_shard in shards:
            ok = out_shard.exists() and Path(str(out_shard)+".tbi").exists() and _is_newer(out_shard, bam)
            if ok and not force_recall:
                console.print(f"[{sample}] shard pronto (cache) → {out_shard.name}")
            else:
                todo.append((i, bed, out_shard))

        # ---- exec paralela dos shards que faltam ----
        if todo:
            console.print(_mk_panel(
                f"[{sample}] Rodando {len(todo)}/{len(shards)} shard(s) (bcftools) "
                f"com paralelismo={bcf_max_parallel} e I/O threads por proc={bcf_io_threads}",
                "cyan"
            ))

            procs = {}   # i -> dict(p, bed, out, tmp, start, last_ts, last_sz, log)
            next_idx = 0

            def _launch(idx: int):
                i, bed, out_shard = todo[idx]
                shard_tag = f"part_{i:02d}"
                console.print(_mk_panel(_bcf_label(sample, bcf_io_threads, bcf_mapq, bcf_baseq, bcf_max_depth, shard_tag), "white"))
                cmd = _bcf_make_cmd(sample, bam, bed, out_shard, bcf_io_threads, bcf_mapq, bcf_baseq, bcf_max_depth)
                _bcf_log_command_and_bed(sample, bed, cmd)  # 🔧 Mostra comando e cromossomos
                if shutil.which("stdbuf"):
                    cmd = ["stdbuf","-oL","-eL"] + cmd
                log_path = Path("logs")/f"bcftools_{sample}_{shard_tag}.log"
                log_fh = open(log_path, "w")
                p = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, text=True)
                procs[i] = {
                    "p": p, "bed": bed, "out": out_shard, "tmp": out_shard.with_suffix(".tmp.vcf.gz"),
                    "start": time.time(), "last_ts": time.time(), "last_sz": 0, "log": log_path, "log_fh": log_fh
                }

            # inicial
            while next_idx < len(todo) and len(procs) < bcf_max_parallel:
                _launch(next_idx); next_idx += 1

            # loop
            while procs:
                finished = []
                now = time.time()
                for i, st in list(procs.items()):
                    p = st["p"]
                    rc = p.poll()
                    # heartbeat
                    if now - st["last_ts"] >= max(3, int(bcf_heartbeat)):
                        sz = st["tmp"].stat().st_size if st["tmp"].exists() else 0
                        dsz = max(0, sz - st["last_sz"])
                        elapsed = int(now - st["start"])
                        console.print(f"[{sample}|part_{i:02d}] … {elapsed//60}m{elapsed%60:02d}s "
                                      f"• out {sizeof_fmt(sz)} • Δ{sizeof_fmt(dsz)}  • log={st['log'].name}",
                                      highlight=False)
                        st["last_ts"] = now
                        st["last_sz"] = sz

                    if rc is None:
                        continue

                    # terminou
                    try:
                        st["log_fh"].flush(); st["log_fh"].close()
                    except Exception:
                        pass
                    finished.append(i)

                    if rc != 0:
                        # Informações detalhadas do erro
                        shard_info = f"part_{i:02d}"
                        console.print(f"[red]❌ BCFtools falhou no shard {shard_info} com código {rc}[/red]")
                        
                        # Mostra cromossomos do shard problemático
                        try:
                            with open(st["bed"]) as f:
                                chroms = set()
                                for line in f:
                                    if line.strip():
                                        chrom = line.split('\t')[0]
                                        chroms.add(chrom)
                                if chroms:
                                    console.print(f"[yellow]🧬 Cromossomos no shard {shard_info}: {', '.join(sorted(chroms))}[/yellow]")
                        except Exception:
                            pass
                        
                        # Mostra comando que falhou
                        console.print(f"[red]💻 Comando que falhou:[/red]")
                        console.print(f"[red]> {' '.join(p.args)}[/red]")
                        
                        # tail do log
                        try:
                            tail = sp.run(["bash","-lc", f"tail -n 60 {shlex.quote(str(st['log']))}"],
                                          capture_output=True, text=True, check=True).stdout
                            if tail.strip():
                                console.print(f"[red]📝 Log do erro (últimas 60 linhas):[/red]")
                                console.print(f"[dim]{tail}[/dim]")
                        except Exception:
                            tail = ""
                        
                        console.print(f"[yellow]💡 Log completo em: {st['log']}[/yellow]")
                        raise sp.CalledProcessError(rc, p.args, output=tail)

                    # rename tmp -> final e index
                    tmp = st["tmp"]
                    out_shard = st["out"]
                    if tmp.exists():
                        _atomic_rename(tmp, out_shard)
                    run(["tabix","-p","vcf",str(out_shard)])
                    dt = int(time.time() - st["start"])
                    size_mb = (out_shard.stat().st_size / (1024*1024)) if out_shard.exists() else 0.0
                    # extrai a linha "Lines ..." do log, se existir
                    try:
                        lines_line = sp.run(
                            ["bash","-lc", f"grep -E '^[Ll]ines\\s' -m1 {shlex.quote(str(st['log']))} || true"],
                            capture_output=True, text=True, check=True
                        ).stdout.strip()
                    except Exception:
                        lines_line = ""
                    console.print(f"[{sample}] shard pronto → {out_shard.name}  •  {size_mb:.1f}MB  •  {dt//60}m{dt%60:02d}s"
                                  + (f"\n{lines_line}" if lines_line else ""))

                    # lança próximo
                    if next_idx < len(todo):
                        _launch(next_idx); next_idx += 1

                    procs.pop(i, None)

                time.sleep(0.3)

        # ---- concatena shards → VCF final (ordem do .fai; -a) ----
        shard_vcfs = [out for (_, _, out) in shards if out.exists()]
        need_final = force_recall or (not final_vcf.exists()) or (not final_tbi.exists()) \
                     or any(out.stat().st_mtime > final_vcf.stat().st_mtime for out in shard_vcfs)

        if not shard_vcfs:
            console.print(f"[orange3][{sample}] Nenhum shard VCF encontrado — nada a concatenar.[/orange3]")
        elif need_final:
            tmp = _tmp_vcfgz_path(final_vcf)
            if tmp.exists(): tmp.unlink()

            # Ordena shards pela ordem do 1º contig do respectivo BED
            shards_for_concat = []
            for (i, bed, out) in shards:
                if out.exists():
                    c = _first_bed_contig(bed)
                    shards_for_concat.append((fai_order.get(c, 10**9), out))
            shards_for_concat.sort(key=lambda x: x[0])
            shard_vcfs_sorted = [str(v) for _, v in shards_for_concat]

            console.print(f"[{sample}] Concat shards → {final_vcf.name}  (n={len(shard_vcfs_sorted)})")
            # -a: permite blocos não contíguos; -D none: não tenta deduplicar registros
            run(["bcftools","concat","-a","-Oz","-o", str(tmp), *shard_vcfs_sorted])
            _atomic_rename(tmp, final_vcf)
            run(["tabix","-p","vcf",str(final_vcf)])
        else:
            console.print(f"[{sample}] VCF final → [bold]SKIP[/bold] (cache ok)")

        vcfs_final.append(final_vcf)
        print_meta(f"VCF ({sample})", [final_vcf])

    if vcfs_final:
        print_meta("VCFs gerados (por amostra)", vcfs_final)
    else:
        console.print("[yellow]Aviso:[/yellow] Nenhum VCF novo foi gerado (talvez já existissem).", style="italic")
    return vcfs_final

def annotate_with_vep(samples, threads: int):
    """
    Passo 9 — Anotação de variantes com Ensembl VEP (usa o VCF final de cada amostra).
    - Lê parâmetros de cfg_global["params"]:
        vep_species (str)       : espécie, ex. "homo_sapiens" (default)
        vep_assembly (str)      : "GRCh38" (default) ou "GRCh37"/"hg19" (normalizado)
        vep_dir_cache (str)     : diretório do cache do VEP (default: ~/.vep)
        vep_fork (int)          : processos internos do VEP (default: max(1, threads//2))
        vep_buffer_size (int)   : tamanho do buffer para processamento em lotes (default: 5000)
        vep_extra (list[str])   : flags extra para o VEP (opcional)
        vep_heartbeat_sec (int) : intervalo do heartbeat de progresso (default: 30)
    - Idempotente: se o .vep.vcf.gz existir e estiver mais novo que o VCF de entrada, faz SKIP.
    - Valida presença do cache (estrutura nova do VEP: ~/.vep/<species>/<version>_<assembly>/).
    - Verifica e corrige automaticamente o formato do FASTA (linhas ≤65536 chars para Bio::DB::Fasta).
    - Garante que o índice .fai existe para o arquivo FASTA de referência.
    - Remove arquivos temporários órfãos e usa --force_overwrite para robustez.
    - Verifica e filtra automaticamente cromossomos problemáticos (virais: chrEBV, chrCMV, etc.).
    - Recupera automaticamente trabalho anterior de arquivos temporários grandes (>100MB).
    - Monitora progresso durante execução com heartbeats informativos (tamanho de saída e tempo decorrido).
    """
    from pathlib import Path
    import os as _os
    import subprocess as sp

    p = cfg_global.get("params", {})

    # ---------------- Helpers ----------------
    def _vep_norm_assembly(a: str) -> str:
        a = (a or "").upper()
        if "GRCH38" in a: return "GRCh38"
        if "GRCH37" in a or "HG19" in a: return "GRCh37"
        return a or "GRCh38"

    def _mk_panel(msg, style="bright_cyan"):
        return Panel.fit(msg, border_style=style)
    
    def _check_and_filter_vcf_for_vep(vcf_in: Path, sample_id: str) -> Path:
        """
        Verifica se VCF contém cromossomos problemáticos para o VEP.
        Se encontrar, cria versão filtrada. Se não, retorna o original.
        Retorna o caminho do VCF a ser usado pelo VEP.
        """
        import subprocess as sp
        
        # Lista de cromossomos problemáticos conhecidos (virais, etc.)
        problem_chroms = {"chrEBV", "chrCMV", "chrHBV", "chrHCV", "chrHIV", "chrHPV", "chrSV40"}
        
        console.print(f"[cyan]🔍 Verificando cromossomos no VCF: {vcf_in.name}[/cyan]")
        
        try:
            # Verifica cromossomos presentes no VCF (rápido)
            cmd_check = ["bcftools", "index", "--stats", str(vcf_in)]
            result = sp.run(cmd_check, capture_output=True, text=True, check=False)
            
            if result.returncode != 0:
                # Fallback: usa zcat + grep se bcftools falhar
                cmd_check = f"zcat {vcf_in} | grep -v '^#' | cut -f1 | sort -u"
                result = sp.run(cmd_check, shell=True, capture_output=True, text=True, check=True)
                chroms_in_vcf = set(result.stdout.strip().split('\n')) if result.stdout.strip() else set()
            else:
                # Extrai cromossomos do stats do bcftools
                chroms_in_vcf = set()
                for line in result.stdout.split('\n'):
                    if line.startswith('CHR\t'):
                        chrom = line.split('\t')[1]
                        chroms_in_vcf.add(chrom)
            
            # Verifica se há cromossomos problemáticos
            problem_found = chroms_in_vcf & problem_chroms
            
            if not problem_found:
                console.print("[green]✅ Nenhum cromossomo problemático encontrado[/green]")
                return vcf_in  # Usa VCF original
                
            # Cromossomos problemáticos encontrados - precisa filtrar
            console.print(f"[yellow]⚠️  Cromossomos problemáticos: {', '.join(sorted(problem_found))}[/yellow]")
            console.print("[cyan]🔄 Criando VCF filtrado para VEP...[/cyan]")
            
            # Cria VCF filtrado
            vcf_filtered = vcf_in.parent / f"{sample_id}.filtered_for_vep.vcf.gz"
            
            # Remove arquivo filtrado anterior se existir
            if vcf_filtered.exists():
                vcf_filtered.unlink()
                
            # Constrói expressão de filtro (exclui cromossomos problemáticos)
            exclude_expr = " || ".join([f'CHROM=="{chrom}"' for chrom in problem_found])
            
            # Filtra VCF
            cmd_filter = [
                "bcftools", "view", 
                "-e", exclude_expr,  # Exclui (-e) cromossomos problemáticos
                "-O", "z", "-o", str(vcf_filtered),
                str(vcf_in)
            ]
            
            console.print(f"[dim]Executando: {' '.join(cmd_filter)}[/dim]")
            sp.run(cmd_filter, check=True)
            
            # Indexa VCF filtrado
            sp.run(["tabix", "-p", "vcf", str(vcf_filtered)], check=True)
            
            console.print(f"[green]✅ VCF filtrado criado: {vcf_filtered.name}[/green]")
            console.print(f"[dim]Excluídos: {', '.join(sorted(problem_found))}[/dim]")
            
            return vcf_filtered  # Usa VCF filtrado
            
        except Exception as e:
            console.print(f"[red]❌ Erro ao verificar/filtrar VCF: {e}[/red]")
            console.print("[yellow]⚠️  Continuando com VCF original (pode falhar em cromossomos problemáticos)[/yellow]")
            return vcf_in  # Fallback para VCF original
    
    def _run_vep_with_progress(cmd: list, sample_id: str, output_file: Path, species: str, assembly: str, fork: int, heartbeat_sec: int = 30):
        """
        Executa VEP monitorando progresso através do tamanho do arquivo de saída.
        """
        import time
        import subprocess as sp
        
        console.print(_mk_panel(f"[{sample_id}] VEP → {species} {assembly} (fork={fork})", "cyan"))
        
        # Mostra comando completo
        cmd_str = " ".join(cmd)
        console.print(f"[dim]> {cmd_str}[/dim]")
        
        start_time = time.time()
        last_heartbeat = start_time
        last_size = 0
        
        # Executa VEP em background
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        
        while True:
            # Verifica se o processo terminou
            rc = proc.poll()
            current_time = time.time()
            
            # Heartbeat periódico
            if current_time - last_heartbeat >= heartbeat_sec:
                elapsed = int(current_time - start_time)
                
                # Verifica tamanho do arquivo de saída
                current_size = output_file.stat().st_size if output_file.exists() else 0
                delta_size = current_size - last_size
                
                # Adiciona contexto sobre a fase do VEP
                phase_info = ""
                if current_size == 0:
                    phase_info = " (iniciando…)"
                elif delta_size == 0 and current_size > 0:
                    if elapsed < 300:  # 5 minutos
                        phase_info = " (carregando cache…)"
                    else:
                        phase_info = " (processando…)"
                elif delta_size > 0:
                    phase_info = " (escrevendo resultados)"
                
                console.print(
                    f"[{sample_id}] VEP … {elapsed//60}m{elapsed%60:02d}s • "
                    f"{output_file.name} {sizeof_fmt(current_size)} • Δ{sizeof_fmt(delta_size)}{phase_info}",
                    highlight=False
                )
                
                last_heartbeat = current_time
                last_size = current_size
            
            if rc is not None:
                break
                
            time.sleep(1)
        
        # Verifica se houve erro
        if proc.returncode != 0:
            stdout, stderr = proc.communicate()
            console.print(f"[red]❌ VEP falhou com código {proc.returncode}[/red]")
            if stderr and stderr.strip():
                console.print(f"[red]Erro VEP:[/red]\n{stderr.strip()}")
            if stdout and stdout.strip():
                console.print(f"[dim]Saída VEP:[/dim]\n{stdout.strip()}")
            raise sp.CalledProcessError(proc.returncode, cmd)
        
        # Relatório final
        elapsed = int(time.time() - start_time)
        final_size = output_file.stat().st_size if output_file.exists() else 0
        console.print(
            f"[bold cyan]{sample_id}[/bold cyan] VEP concluído em {elapsed//60}m{elapsed%60:02d}s • "
            f"{output_file.name} {sizeof_fmt(final_size)}"
        )
    
    def _ensure_fasta_format_for_vep(ref_fa: Path) -> bool:
        """
        Verifica se o arquivo FASTA tem linhas compatíveis com Bio::DB::Fasta (≤65536 chars).
        Se necessário, reformata o arquivo para ter linhas de 80 caracteres.
        Retorna True se alguma correção foi feita.
        """
        if not ref_fa.exists():
            return False
            
        console.print(f"[cyan]Verificando formato do FASTA: {ref_fa.name}[/cyan]")
        
        # Verifica se há linhas muito longas (usa awk para eficiência)
        needs_reformat = False
        try:
            import subprocess as sp
            result = sp.run(
                ["awk", "length($0) > 65536 && !/^>/ {print NR \":\" length($0); exit}", str(ref_fa)],
                capture_output=True, text=True, check=False
            )
            if result.stdout.strip():
                line_info = result.stdout.strip()
                needs_reformat = True
                console.print(f"[yellow]⚠️  FASTA linha {line_info} caracteres (>65536) — incompatível com VEP[/yellow]")
            else:
                console.print("[green]✅ FASTA já está no formato correto (linhas ≤65536 chars)[/green]")
        except Exception:
            return False
            
        if not needs_reformat:
            return False
            
        # Reformata o arquivo
        console.print(_mk_panel("[bold]Reformatando FASTA para compatibilidade com VEP[/bold]\n"
                               "• Criando backup do arquivo original\n"
                               "• Reformatando para linhas de 80 caracteres\n"
                               "• Recriando índice .fai", "yellow"))
        
        backup_path = ref_fa.with_suffix('.fa.backup')
        temp_path = ref_fa.with_suffix('.fa.tmp')
        
        # Backup do original (apenas se não existir)
        if not backup_path.exists():
            import shutil
            console.print(f"[cyan]📦 Criando backup: {backup_path.name}[/cyan]")
            shutil.copy2(ref_fa, backup_path)
        else:
            console.print(f"[dim]📦 Backup já existe: {backup_path.name}[/dim]")
            
        # Reformata usando seqtk
        console.print("[cyan]🔄 Reformatando FASTA (seqtk seq -l 80)...[/cyan]")
        cmd = ["seqtk", "seq", "-l", "80", str(backup_path)]
        try:
            with open(temp_path, 'w') as f:
                sp.run(cmd, stdout=f, check=True)
            
            # Substitui o original
            temp_path.replace(ref_fa)
            console.print(f"[green]✅ Arquivo reformatado: {ref_fa.name}[/green]")
            
            # Recria o índice
            fai_path = ref_fa.with_suffix('.fa.fai')
            if fai_path.exists():
                fai_path.unlink()
            console.print("[cyan]📇 Recriando índice .fai...[/cyan]")
            sp.run(["samtools", "faidx", str(ref_fa)], check=True)
            
            console.print("[green]✅ FASTA reformatado e indexado com sucesso![/green]")
            return True
            
        except Exception as e:
            console.print(f"[red]❌ Erro ao reformatar FASTA: {e}[/red]")
            # Restaura o backup se algo deu errado
            if backup_path.exists() and temp_path.exists():
                temp_path.unlink()
            return False

    # --------------- Parâmetros ---------------
    vep_species     = (p.get("vep_species") or "homo_sapiens").lower()
    vep_assembly    = _vep_norm_assembly(p.get("vep_assembly") or "GRCh38")
    vep_dir_cache   = p.get("vep_dir_cache") or str(Path.home()/".vep")
    vep_fork        = int(p.get("vep_fork", max(1, threads//2)))
    vep_extra       = list(p.get("vep_extra", []))
    vep_heartbeat   = int(p.get("vep_heartbeat_sec", 30))

    Path("vep").mkdir(exist_ok=True)

    # --------------- VEP disponível? ---------------
    try:
        rc = sp.run(["vep"], stdout=sp.PIPE, stderr=sp.STDOUT, text=True, check=False)
        vep_ver = (rc.stdout or "").strip()
        if rc.returncode != 0:
            console.print("[red]VEP não encontrado no PATH do ambiente atual.[/red] "
                          "Ative o env correto (ex.: 'conda activate genomics') ou reinstale.")
            raise SystemExit(1)
        else:
            console.print(_mk_panel(f"Ensembl VEP detectado:\n{vep_ver}", "green"))
    except FileNotFoundError:
        console.print("[red]VEP não encontrado (binário 'vep').[/red]")
        raise SystemExit(1)

    # --------------- Verificação do FASTA ---------------
    console.print(_mk_panel("[bold]Preparação do FASTA para VEP[/bold]\n"
                           "• Verificando compatibilidade com Bio::DB::Fasta\n"
                           "• Reformatando se necessário (linhas ≤65536 chars)", "cyan"))
    
    ref_fa = Path("refs/reference.fa")
    if ref_fa.exists():
        _ensure_fasta_format_for_vep(ref_fa)
    else:
        console.print(f"[yellow]⚠️  FASTA de referência não encontrado: {ref_fa}[/yellow]")

    # --------------- Checagem de cache ---------------
    cache_root = Path(vep_dir_cache)
    species_dir = cache_root / vep_species

    # Estrutura típica nova: ~/.vep/homo_sapiens/110_GRCh38
    cache_ok = False
    missing_msg = ""
    if species_dir.exists():
        # Procura subpastas que contenham o assembly no nome (ex.: "110_GRCh38")
        try:
            for child in species_dir.iterdir():
                if not child.is_dir():
                    continue
                name = child.name.upper()
                if vep_assembly.upper() in name:
                    cache_ok = True
                    break
            if not cache_ok:
                missing_msg = f"Cache de {vep_species} não tem subpasta para {vep_assembly} em {species_dir}"
        except Exception as _:
            pass
    else:
        missing_msg = f"Diretório da espécie não encontrado: {species_dir}"

    if not cache_ok:
        console.print(
            f"[red]VEP cache não encontrado para {vep_species} / {vep_assembly}.[/red]\n"
            f"{missing_msg}\n"
            f"Sugestão: instalar com [bold]vep_install[/bold], ex.:\n"
            f"  vep_install -a cf -s {vep_species} -y {vep_assembly} -c {vep_dir_cache} --NO_BIOPERL"
        )
        raise SystemExit(1)

    # --------------- Execução por amostra ---------------
    annotated = []
    for sample in samples:
        sample_id = sample["id"] if isinstance(sample, dict) else sample
        vcf_original = Path("vcf")/f"{sample_id}.vcf.gz"
        vcf_tbi = Path(str(vcf_original)+".tbi")
        if not vcf_original.exists():
            console.print(f"[yellow][{sample_id}] VCF de entrada não encontrado:[/yellow] {vcf_original}")
            continue

        # Verifica e filtra cromossomos problemáticos se necessário
        vcf_in = _check_and_filter_vcf_for_vep(vcf_original, sample_id)

        out_vep = Path("vep")/f"{sample_id}.vep.vcf.gz"
        tmp_vep = Path("vep")/f"{sample_id}.vep.vcf.gz.tmp"

        # Skip se já está atualizado
        if out_vep.exists() and _is_newer(out_vep, vcf_in):
            console.print(f"[{sample_id}] VEP → [bold]SKIP[/bold] (cache ok)")
            annotated.append(out_vep)
            continue

        # Gerencia arquivos temporários de forma inteligente
        if tmp_vep.exists():
            tmp_size = tmp_vep.stat().st_size
            # Se arquivo temporário é grande (>100MB), pode conter trabalho válido
            if tmp_size > 100 * 1024 * 1024:  # 100MB
                console.print(f"[yellow]⚠️  Arquivo temporário grande encontrado: {tmp_vep.name} ({sizeof_fmt(tmp_size)})[/yellow]")
                
                # Tenta recuperar como arquivo final se não existe
                if not out_vep.exists():
                    console.print(f"[cyan]🔄 Tentando recuperar trabalho anterior...[/cyan]")
                    try:
                        # Verifica se é VCF válido
                        result = sp.run(["head", "-50", str(tmp_vep)], capture_output=True, text=True)
                        if "##fileformat=VCF" in result.stdout:
                            # Conta variantes no arquivo
                            count_result = sp.run(["grep", "-v", "^#", str(tmp_vep)], 
                                                capture_output=True, text=True)
                            if count_result.returncode == 0:
                                variant_count = len(count_result.stdout.strip().split('\n')) if count_result.stdout.strip() else 0
                                console.print(f"[cyan]📊 Encontradas {variant_count:,} variantes anotadas[/cyan]")
                            
                            # Comprime e move para arquivo final
                            console.print(f"[cyan]📦 Comprimindo arquivo recuperado...[/cyan]")
                            with open(out_vep, 'wb') as f:
                                sp.run(["bgzip", "-c", str(tmp_vep)], stdout=f, check=True)
                            sp.run(["tabix", "-p", "vcf", str(out_vep)], check=True)
                            
                            # Remove temporário após sucesso
                            tmp_vep.unlink()
                            console.print(f"[green]✅ Trabalho anterior recuperado: {out_vep.name}[/green]")
                            annotated.append(out_vep)
                            continue
                        else:
                            console.print(f"[red]❌ Arquivo temporário não é VCF válido[/red]")
                    except Exception as e:
                        console.print(f"[red]❌ Erro ao recuperar arquivo: {e}[/red]")
                
                # Se não conseguiu recuperar ou arquivo final já existe, faz backup
                backup_path = tmp_vep.with_suffix('.tmp.backup')
                if not backup_path.exists():
                    console.print(f"[cyan]📦 Fazendo backup do arquivo temporário grande...[/cyan]")
                    import shutil
                    shutil.move(str(tmp_vep), str(backup_path))
                    console.print(f"[dim]Backup salvo como: {backup_path.name}[/dim]")
                else:
                    tmp_vep.unlink()
                    console.print(f"[cyan]🧹 Removido arquivo temporário (backup já existe)[/cyan]")
            else:
                # Arquivo pequeno (<100MB) - remove normalmente
                tmp_vep.unlink()
                console.print(f"[cyan]🧹 Removido arquivo temporário pequeno: {tmp_vep.name} ({sizeof_fmt(tmp_size)})[/cyan]")

        # Constrói comando
        ref_fa = Path("refs/reference.fa")
        
        # Garante que o índice .fai existe
        fai_path = ref_fa.with_suffix('.fa.fai')
        if ref_fa.exists() and not fai_path.exists():
            console.print(f"[cyan]📇 Criando índice .fai para {ref_fa.name}...[/cyan]")
            sp.run(["samtools", "faidx", str(ref_fa)], check=True)
            console.print(f"[green]✅ Índice criado: {fai_path.name}[/green]")
        
        # Parâmetros otimizados para performance
        buffer_size = p.get("vep_buffer_size", 5000)  # Default do VEP é 5000
        
        cmd = [
            "vep",
            "--cache", "--offline",
            "--species", vep_species,
            "--assembly", vep_assembly,
            "--dir_cache", vep_dir_cache,
            "--fork", str(vep_fork),
            "--format", "vcf", "--vcf",
            "-i", str(vcf_in),
            "-o", str(tmp_vep),
            "--fasta", str(ref_fa),
            "--force_overwrite",
            "--buffer_size", str(buffer_size)
        ]
        
        # Adiciona compressão se suportada (VEP mais novos)
        if p.get("vep_compress_output", False):  # Mudado para False por padrão
            cmd.extend(["--compress_output", "bgzip"])
            
        cmd += vep_extra

        # Executa VEP com monitoramento de progresso
        _run_vep_with_progress(cmd, sample_id, tmp_vep, vep_species, vep_assembly, vep_fork, heartbeat_sec=vep_heartbeat)

        # Atomic rename
        _atomic_rename(tmp_vep, out_vep)
        annotated.append(out_vep)

        # Mostra meta
        print_meta(f"VEP ({sample_id})", [out_vep])

    if annotated:
        print_meta("VCFs anotados (VEP) — por amostra", annotated)
    else:
        console.print("[yellow]Aviso:[/yellow] Nenhum VCF anotado (verifique entradas/erros).", style="italic")
    return annotated

# =================== Comparações par-a-par ===================

def _parse_gt(gt: str):
    # Ex.: "0/1", "1/1", "./.", "0|1"
    if not gt or gt == "." or gt.startswith("."):
        return None
    alleles = gt.replace("|","/").split("/")
    try:
        return tuple(sorted(int(a) for a in alleles if a != "."))
    except Exception:
        return None

def _share_allele(gt1, gt2) -> bool:
    if gt1 is None or gt2 is None:
        return False
    return bool(set(gt1) & set(gt2))

def _same_geno(gt1, gt2) -> bool:
    return gt1 is not None and gt2 is not None and tuple(gt1) == tuple(gt2)

def _is_het(gt) -> bool:
    return gt is not None and len(set(gt)) > 1

def pairwise_comparisons(dna_samples):
    Path("comparisons").mkdir(exist_ok=True)
    sample_ids = sorted(set(s["id"] for s in dna_samples))
    vcfs = {sid: Path("vcf")/f"{sid}.vcf.gz" for sid in sample_ids}
    for sid, v in vcfs.items():
        if not v.exists():
            console.print(f"[yellow]Aviso:[/yellow] VCF ausente para {sid}: {v}", style="italic")

    from itertools import combinations
    for a,b in combinations(sample_ids, 2):
        va, vb = vcfs[a], vcfs[b]
        if not (va.exists() and vb.exists()):
            console.print(f"[orange3]Pulo comparação {a} vs {b} (VCF faltando).[/orange3]")
            continue

        pair_name = f"{a}_vs_{b}"
        out_merge = Path("comparisons")/f"{pair_name}.merge.vcf.gz"
        need_merge = not (out_merge.exists() and _is_newer(out_merge, va, vb))
        if need_merge:
            run(["bcftools","merge","-m","all","-Oz","-o",str(out_merge), str(va), str(vb)])
            run(["bcftools","index","-t",str(out_merge)])

        # query e métricas (sempre recalculadas — idempotente por overwrite)
        q = sp.run(
            ["bcftools","query","-s",f"{a},{b}","-f","%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",str(out_merge)],
            capture_output=True, text=True, check=True
        ).stdout.splitlines()

        stats = {"pair": f"{a}_vs_{b}", "sites_total": 0, "sites_both_called": 0,
                 "geno_exact_match": 0, "het_concord": 0, "share_allele": 0, "ibs0": 0}

        for line in q:
            if not line or line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            gt_a = _parse_gt(parts[-2]); gt_b = _parse_gt(parts[-1])
            stats["sites_total"] += 1
            both = (gt_a is not None) and (gt_b is not None)
            if not both: continue
            stats["sites_both_called"] += 1
            if _same_geno(gt_a, gt_b): stats["geno_exact_match"] += 1
            if _is_het(gt_a) and _is_het(gt_b) and _same_geno(gt_a, gt_b): stats["het_concord"] += 1
            if _share_allele(gt_a, gt_b): stats["share_allele"] += 1
            else: stats["ibs0"] += 1

        bc = max(1, stats["sites_both_called"])
        frac_exact = stats["geno_exact_match"]/bc
        frac_share = stats["share_allele"]/bc
        frac_ibs0  = stats["ibs0"]/bc

        tsv = Path("comparisons")/f"{pair_name}.metrics.tsv"
        with open(tsv,"w") as fh:
            fh.write("metric\tvalue\n")
            for k,v in stats.items(): fh.write(f"{k}\t{v}\n")
            fh.write(f"geno_exact_match_frac\t{frac_exact:.6f}\n")
            fh.write(f"share_allele_frac\t{frac_share:.6f}\n")
            fh.write(f"ibs0_frac\t{frac_ibs0:.6f}\n")

        md = Path("comparisons")/f"{pair_name}.summary.md"
        with open(md,"w") as fh:
            fh.write(f"# {a} vs {b}\n\n")
            fh.write(f"- Sites com genótipo em ambos: **{stats['sites_both_called']:,}**\n")
            fh.write(f"- Concordância exata de genótipo: **{frac_exact:.3%}** "
                     f"({stats['geno_exact_match']:,}/{bc:,})\n")
            fh.write(f"- Compartilham ≥1 alelo (IBS≥1): **{frac_share:.3%}** "
                     f"({stats['share_allele']:,}/{bc:,})\n")
            fh.write(f"- **IBS0** (nenhum alelo em comum): **{frac_ibs0:.3%}** "
                     f"({stats['ibs0']:,}/{bc:,})\n\n")
            fh.write("> Regra prática: parentais-filhos apresentam IBS0 muito baixo e "
                     "alta fração de compartilhamento de alelos.\n")

        console.print(f"[bold cyan]Comparação {pair_name}:[/bold cyan] "
                      f"IBS0={frac_ibs0:.2%}, share≥1={frac_share:.2%}, exact={frac_exact:.2%}")

# =================== Consolida par-a-par ===================

def combine_pairwise_stats(dna_samples):
    """
    Lê comparisons/*_vs_*.metrics.tsv e escreve comparisons/summary_pairwise.csv
    com as métricas principais por par.
    """
    Path("comparisons").mkdir(exist_ok=True)
    rows = []
    for tsv in sorted(Path("comparisons").glob("*_vs_*.metrics.tsv")):
        pair = None
        d = {}
        with open(tsv) as fh:
            for line in fh:
                if not line.strip(): continue
                k, v = line.rstrip("\n").split("\t")
                if k == "pair":
                    pair = v
                else:
                    # números inteiros/floats quando possível
                    try:
                        if "." in v: d[k] = float(v)
                        else: d[k] = int(v)
                    except Exception:
                        d[k] = v
        if not pair: 
            continue
        a,b = pair.split("_vs_")
        rows.append({
            "sample_a": a, "sample_b": b,
            "sites_both_called": d.get("sites_both_called", 0),
            "geno_exact_match_frac": d.get("geno_exact_match_frac", 0.0),
            "share_allele_frac": d.get("share_allele_frac", 0.0),
            "ibs0_frac": d.get("ibs0_frac", 0.0),
        })
    # escreve CSV
    out_csv = Path("comparisons")/"summary_pairwise.csv"
    with open(out_csv,"w") as out:
        out.write("sample_a,sample_b,sites_both_called,geno_exact_match_frac,share_allele_frac,ibs0_frac\n")
        for r in rows:
            out.write(f"{r['sample_a']},{r['sample_b']},{r['sites_both_called']},"
                      f"{r['geno_exact_match_frac']:.6f},{r['share_allele_frac']:.6f},{r['ibs0_frac']:.6f}\n")
    console.print(f"[bold cyan]CSV par-a-par[/bold cyan] → {out_csv}")


# =================== Trio: potenciais de novo ===================

def _infer_trio_ids(dna_samples):
    """Define child, parents a partir do YAML (se houver) ou heurística."""
    g = cfg_global.get("general", {})
    sample_ids = [s["id"] for s in dna_samples]
    child = g.get("trio_child_id")
    parents = g.get("trio_parent_ids", [])

    if child and parents and len(parents)==2:
        return child, parents[0], parents[1]

    # heurística: se NA12878 existir, é filha
    if "NA12878" in sample_ids and len(sample_ids) >= 3:
        others = [x for x in sample_ids if x != "NA12878"]
        return "NA12878", others[0], others[1]

    # fallback: pega os 3 primeiros
    if len(sample_ids) >= 3:
        console.print("[yellow]Aviso:[/yellow] trio não especificado; assumindo primeiro como filho.", style="italic")
        return sample_ids[0], sample_ids[1], sample_ids[2]

    return None, None, None

def _parse_gt(gt: str):
    if not gt or gt == "." or gt.startswith("."):
        return None
    alleles = gt.replace("|","/").split("/")
    try:
        return tuple(sorted(int(a) for a in alleles if a != "."))
    except Exception:
        return None

def _is_hom_ref(gt): return gt is not None and len(gt)==1 and gt[0]==0
def _is_nonref(gt):   return gt is not None and not _is_hom_ref(gt)

def _parse_ad(ad_str):
    # AD = "ref,alt[,alt2...]"
    try:
        parts = [int(x) for x in ad_str.split(",") if x != "."]
        if not parts: return None, None, None
        ref = parts[0]; alt = sum(parts[1:]) if len(parts)>1 else 0
        tot = ref + alt
        ab = (alt/tot) if tot>0 else None
        return ref, alt, ab
    except Exception:
        return None, None, None

def trio_denovo_report(dna_samples):
    child, p1, p2 = _infer_trio_ids(dna_samples)
    if not child:
        console.print("[orange3]Trio não identificado — pulando relatório de de novo.[/orange3]")
        return

    vcfs = { s["id"]: Path("vcf")/f"{s['id']}.vcf.gz" for s in dna_samples }
    for sid in [child, p1, p2]:
        if not vcfs.get(sid, Path()).exists():
            console.print(f"[orange3]Sem VCF para {sid} — pulando trio.[/orange3]")
            return

    Path("trio").mkdir(exist_ok=True)
    merged = Path("trio")/"trio_merged.vcf.gz"
    need_merge = not (merged.exists() and _is_newer(merged, vcfs[child], vcfs[p1], vcfs[p2]))
    if need_merge:
        run(["bcftools","merge","-m","all","-Oz","-o",str(merged),
             str(vcfs[child]), str(vcfs[p1]), str(vcfs[p2])])
        run(["bcftools","index","-t",str(merged)])

    # thresholds (como antes)
    g = cfg_global.get("general", {})
    min_dp_child = int(g.get("trio_min_dp_child", 8))
    min_dp_par   = int(g.get("trio_min_dp_parents", 8))
    min_gq       = int(g.get("trio_min_gq", 20))
    min_ab_het   = float(g.get("trio_min_ab_het", 0.25))
    max_ab_het   = float(g.get("trio_max_ab_het", 0.75))
    min_ab_hom   = float(g.get("trio_min_ab_hom", 0.90))
    max_par_alt  = float(g.get("trio_max_parent_alt_frac", 0.02))

    q = sp.run(
        ["bash","-lc",
         f"bcftools view -f PASS {shlex.quote(str(merged))} | "
         f"bcftools query -s {child},{p1},{p2} -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP:%GQ:%AD]\\n'"],
        capture_output=True, text=True, check=True
    ).stdout.splitlines()

    out_tsv = Path("trio")/"trio_denovo_candidates.tsv"
    kept = 0; total = 0
    with open(out_tsv, "w") as out:
        out.write("chrom\tpos\tref\talt\tchild_GT\tchild_DP\tchild_GQ\tchild_AB\tp1_GT\tp1_DP\tp1_GQ\tp1_AB\tp2_GT\tp2_DP\tp2_GQ\tp2_AB\n")
        for line in q:
            if not line.strip(): continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:  # precisa ter 3 blocos de amostra
                continue
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            sa, sb, sc = parts[4], parts[5], parts[6]
            def split_block(b):
                xx = b.split(":")
                GT = xx[0] if len(xx)>0 else "."
                DP = int(xx[1]) if len(xx)>1 and xx[1].isdigit() else None
                GQ = int(xx[2]) if len(xx)>2 and xx[2].isdigit() else None
                AD = xx[3] if len(xx)>3 else ""
                return GT, DP, GQ, AD

            cGT,cDP,cGQ,cAD = split_block(sa)
            p1GT,p1DP,p1GQ,p1AD = split_block(sb)
            p2GT,p2DP,p2GQ,p2AD = split_block(sc)

            gt_c  = _parse_gt(cGT); gt_p1 = _parse_gt(p1GT); gt_p2 = _parse_gt(p2GT)
            _,_,ab_c  = _parse_ad(cAD);  _,_,ab_p1 = _parse_ad(p1AD);  _,_,ab_p2 = _parse_ad(p2AD)

            # pais
            if not _is_hom_ref(gt_p1) or not _is_hom_ref(gt_p2): continue
            if (p1DP is not None and p1DP < min_dp_par) or (p2DP is not None and p2DP < min_dp_par): continue
            if (p1GQ is not None and p1GQ < min_gq) or (p2GQ is not None and p2GQ < min_gq): continue
            if ab_p1 is not None and ab_p1 > max_par_alt: continue
            if ab_p2 is not None and ab_p2 > max_par_alt: continue

            # filho
            if not _is_nonref(gt_c): continue
            if (cDP is not None and cDP < min_dp_child): continue
            if (cGQ is not None and cGQ < min_gq): continue
            if ab_c is not None:
                if len(set(gt_c)) > 1:
                    if not (min_ab_het <= ab_c <= max_ab_het): continue
                else:
                    if ab_c < min_ab_hom: continue

            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{cGT}\t{cDP}\t{cGQ}\t{ab_c if ab_c is not None else ''}"
                      f"\t{p1GT}\t{p1DP}\t{p1GQ}\t{ab_p1 if ab_p1 is not None else ''}"
                      f"\t{p2GT}\t{p2DP}\t{p2GQ}\t{ab_p2 if ab_p2 is not None else ''}\n")
            kept += 1

    md = Path("trio")/"trio_denovo_summary.md"
    with open(md,"w") as fh:
        fh.write(f"# Trio de novo — {child} vs {p1},{p2}\n\n")
        fh.write(f"- Variantes elegíveis (PASS) avaliadas: **{total:,}**\n")
        fh.write(f"- Candidatos *de novo* após filtros: **{kept:,}**\n\n")
        fh.write(f"- TSV: `trio/{out_tsv.name}`\n- Merged VCF: `trio/{merged.name}`\n")

    console.print(f"[bold cyan]Trio de novo[/bold cyan] → {out_tsv}  (n={kept})")

# =================== Lista de genes ===================

def gene_list_from_gtf():
    Path("genes").mkdir(exist_ok=True)
    out = Path("genes/gene_list.txt")
    gtf = Path("refs/genes.gtf")
    # refaz se o out não existe OU é mais antigo que o GTF
    if out.exists() and _is_newer(out, gtf):
        console.print("Lista de genes → [bold]SKIP (atual)[/bold]")
        print_meta("Lista de genes", [out])
        return
    cmd = r'''awk '$3=="gene"{print $0}' refs/genes.gtf | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p' | sort -u > genes/gene_list.txt'''
    run(["bash","-lc",cmd])
    print_meta("Lista de genes", [out])

# ============ Presença de genes por amostra (cobertura) ============

def _build_genes_bed_from_gtf(gtf: Path, out_bed: Path):
    """Extrai regiões 'gene' do GTF como BED (chr start end name). Recria se o GTF for mais novo."""
    if out_bed.exists() and _is_newer(out_bed, gtf):
        return out_bed
    cmd = r"""awk '$3=="gene"{ 
        chr=$1; start=$4-1; end=$5; 
        match($0,/gene_id "([^"]+)"/,a); gid=a[1];
        match($0,/gene_name "([^"]+)"/,b); gname=b[1];
        match($0,/gene_type "([^"]+)"/,c); gtype=c[1];
        match($0,/gene_biotype "([^"]+)"/,d); if(d[1]=="") d[1]=c[1];
        match($0,/gene_description "([^"]+)"/,e); gdesc=e[1];
        name=gname!=""?gname:gid;
        print chr"\t"start"\t"end"\t"gid"\t"name"\t"gtype"\t"gdesc;
    }' """ + shlex.quote(str(gtf)) + " > " + shlex.quote(str(out_bed))
    run(["bash","-lc",cmd])
    return out_bed

def _parse_thresholds_bed(thr_bed_gz: Path, region_len: dict) -> dict:
    """
    Lê *.thresholds.bed.gz do mosdepth --by, retorna breadth_1x por gene_id.
    Formato: chrom start end region_name thresh1 thresh2 ... counts_above
    Usamos a 1a coluna de contagem como >=1x.
    """
    import gzip
    breadth = {}
    with gzip.open(thr_bed_gz, "rt") as fh:
        for line in fh:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5: continue
            gid = parts[3]
            length = max(1, int(parts[2]) - int(parts[1]))
            # coluna 4 é o nome; as demais são contagens por threshold
            # mosdepth escreve os thresholds depois das 4 primeiras colunas
            counts = [int(x) for x in parts[4:]] if len(parts) > 4 else [0]
            cov1 = counts[0] if counts else 0
            breadth[gid] = cov1 / length
            region_len[gid] = length
    return breadth

def gene_presence_reports(dna_samples, min_mean_cov: float = 5.0, min_breadth_1x: float = 0.8, threads: int = 8):
    Path("genes").mkdir(exist_ok=True)
    gtf = Path("refs/genes.gtf")
    genes_bed = _build_genes_bed_from_gtf(gtf, Path("genes/genes.bed"))

    # gene meta (id → (name,type,desc))
    meta = {}
    with open(genes_bed) as fh:
        for line in fh:
            if not line.strip(): continue
            chrom, start, end, gid, gname, gtype, gdesc = line.rstrip("\n").split("\t")
            meta[gid] = (gname, gtype if gtype else "", gdesc if gdesc else "")

    for bam in sorted(Path("bam").glob("*.mkdup.cram")) + sorted(Path("bam").glob("*.mkdup.bam")):
        sample = bam.name.split(".")[0]
        prefix = Path("genes")/f"{sample}.genes"
        regions_bed_gz = Path(str(prefix)+".regions.bed.gz")
        thr_bed_gz     = Path(str(prefix)+".thresholds.bed.gz")

        # (Re)executa mosdepth se outputs faltam ou estão mais antigos que BAM/CRAM ou genes.bed
        need_mosdepth = True
        if regions_bed_gz.exists() and thr_bed_gz.exists():
            newer_regions = _is_newer(regions_bed_gz, bam, genes_bed)
            newer_thr     = _is_newer(thr_bed_gz, bam, genes_bed)
            need_mosdepth = not (newer_regions and newer_thr)

        if need_mosdepth:
            cmd = ["mosdepth","-t",str(threads),"--by",str(genes_bed),
                   "--thresholds","1,5,10",str(prefix), str(bam)]
            if bam.suffix == ".cram":
                # por segurança, garante FASTA quando entrada é CRAM
                cmd = ["mosdepth","-t",str(threads),"--by",str(genes_bed),
                       "--thresholds","1,5,10","--fasta","refs/reference.fa",
                       str(prefix), str(bam)]
            run(cmd)

        # parse dos outputs (como antes)
        import gzip
        mean_cov = {}
        with gzip.open(regions_bed_gz, "rt") as fh:
            for line in fh:
                if not line.strip(): continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8: continue
                gid = parts[3]
                mean = float(parts[7])
                mean_cov[gid] = mean

        region_len = {}
        breadth1 = _parse_thresholds_bed(thr_bed_gz, region_len)

        out_tsv = Path("genes")/f"{sample}_gene_presence.tsv"
        with open(out_tsv,"w") as out:
            out.write("gene_id\tgene_name\tgene_type\tlength_bp\tmean_cov\tbreadth_1x\tpresent\tmini_summary\n")
            for gid,(gname,gtype,gdesc) in meta.items():
                length = region_len.get(gid, 0)
                mean   = mean_cov.get(gid, 0.0)
                br1x   = breadth1.get(gid, 0.0)
                present = (mean >= min_mean_cov) and (br1x >= min_breadth_1x)
                mini = gdesc if gdesc else (gtype if gtype else "")
                out.write(f"{gid}\t{gname}\t{gtype}\t{length}\t{mean:.2f}\t{br1x:.3f}\t{str(present)}\t{mini}\n")

        console.print(f"[bold cyan]Genes presentes ({sample})[/bold cyan] → {out_tsv}")

# =================== RNA-seq (opcional) ===================

def rnaseq_pipeline(rna_samples, threads, assembly_name):
    if not rna_samples: return
    Path("rnaseq").mkdir(exist_ok=True)

    # 1) garantir FASTQs só das amostras RNA
    for s in rna_samples:
        if s["source"]=="sra":
            stage_fastqs_from_sra(s["sra_ids"])
        else:
            stage_fastqs_from_local(s["fastq1"], s.get("fastq2"))

    # helper para achar R2 correspondente ao R1
    def _find_r2_for_r1_file(r1: Path) -> Optional[Path]:
        name = r1.name
        for pat, rep in (
            (r"_1\.fastq\.gz$",  "_2.fastq.gz"),
            (r"_1\.fq\.gz$",     "_2.fq.gz"),
            (r"R1\.fastq\.gz$",  "R2.fastq.gz"),
            (r"R1\.fq\.gz$",     "R2.fq.gz"),
            (r"\.1\.fastq\.gz$", ".2.fastq.gz"),
            (r"\.1\.fq\.gz$",    ".2.fq.gz"),
        ):
            if re.search(pat, name, flags=re.I):
                cand = r1.parent / re.sub(pat, rep, name, flags=re.I)
                if cand.exists(): return cand
        return None

    # 2) lista de R1s esperados só das amostras RNA
    r1_list: list[Path] = []
    for s in rna_samples:
        if s["source"]=="sra":
            for acc in s["sra_ids"]:
                p = Path("fastq")/f"{acc}_1.fastq.gz"
                if p.exists(): r1_list.append(p)
        else:
            fq1 = Path("fastq")/Path(s["fastq1"]).name
            if fq1.exists(): r1_list.append(fq1)

    # 3) processa cada R1
    for r1 in sorted(set(r1_list)):
        base = re.sub(r'_1\.fastq\.gz$', '', r1.name, flags=re.I)
        r2 = _find_r2_for_r1_file(r1)
        gtf_out = Path("rnaseq")/f"{base}.transcripts.gtf"
        bam_out = Path("rnaseq")/f"{base}.rnaseq.bam"

        # precisa reprocessar?
        need = not (gtf_out.exists() and _is_newer(gtf_out, r1, r2 if r2 else r1, Path("refs/genes.gtf")))
        if not need:
            console.print(f"[RNA-seq:{base}] GTF → [bold]SKIP (atual)[/bold]")
            continue

        # alinhamento HISAT2 → BAM
        if r2 and r2.exists():
            cmd_list = [["hisat2","-p",str(threads),"-x",f"refs/{assembly_name}","-1",str(r1),"-2",str(r2)],
                        ["samtools","sort","-@",str(threads),"-o",str(bam_out)]]
        else:
            cmd_list = [["hisat2","-p",str(threads),"-x",f"refs/{assembly_name}","-U",str(r1)],
                        ["samtools","sort","-@",str(threads),"-o",str(bam_out)]]
        run_long_stream_pipeline(cmd_list, label=f"[RNA-seq:{base}] HISAT2 → sort")
        run(["samtools","index",str(bam_out)])

        # StringTie (reconstrói se GTF faltando/antigo)
        run(["stringtie",str(bam_out),"-G","refs/genes.gtf","-o",str(gtf_out),"-p",str(threads)])

    # 4) gffcompare (idempotente: só refaz se outputs faltarem)
    tgts = list(Path("rnaseq").glob("*.transcripts.gtf"))
    cmp_prefix = Path("rnaseq")/"cmp"
    need_cmp = tgts and (not (cmp_prefix.with_suffix(".stats").exists()))
    if tgts and need_cmp:
        run(["gffcompare","-r","refs/genes.gtf","-o",str(cmp_prefix),*map(str,tgts)])

# ========================== Main ==========================

def main(cfg):
    g = cfg["general"]
    base_dir = Path(g.get("base_dir",".")).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(base_dir)

    # TMP global
    tmpdir = Path(g.get("temp_dir","tmp")).expanduser().resolve()
    tmpdir.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(tmpdir)
    console.print(f"Usando TMPDIR={tmpdir} (do YAML)", style="dim")

    console.rule(f"[bold]Pipeline Humano[/bold]  →  [cyan]{base_dir}[/cyan]")
    disk_space_report(base_dir)
    ensure_dirs()

    # Referências e índices
    stage_banner("1) Referências e Índices")
    download_refs(g["ref_fa_url"], g["gtf_url"], force=g.get("force_refs", False))
    limit_reference_to_canonical_if_enabled()
    build_indexes(g.get("default_read_type","short"), g["assembly_name"],
                  g.get("rnaseq",False) or bool(cfg.get("rna_samples",[])),
                  g["threads"], force=g.get("force_indexes", False))
    step_done("Referências/índices prontos")

    # Entrada DNA
    stage_banner("2) Entrada de Leituras (FASTQ)", "SRA→FASTQ, cache e compressão")
    dna_samples = cfg.get("dna_samples", [])
    for s in dna_samples:
        if s["source"]=="sra":
            stage_fastqs_from_sra(s["sra_ids"])
        else:
            stage_fastqs_from_local(s["fastq1"], s.get("fastq2"))
    step_done("FASTQs prontos")

    # Estimativa de espaço
    stage_banner("3) Estimativa de Espaço")
    fq_list = list(Path("fastq").glob("*.fastq.gz"))
    estimate_outputs_and_warn(fq_list, base_dir, g.get("space_guard_gb_min", 100))
    step_done("Estimativa impressa")

    # Downsample (opcional)
    if g.get("downsample_frac", 0.0) > 0:
        stage_banner("4) Downsample (opcional)")
        downsample_fastqs(g["downsample_frac"], g.get("downsample_seed", 123))
        step_done("Downsample concluído")

    # QC + trimming
    stage_banner("5) QC e Trimming")
    qc_and_trim(g["threads"], g["adapters"], g.get("default_read_type","short"),
                use_ds=g.get("downsample_frac", 0.0) > 0)
    step_done("QC + trimming concluídos")

    # Alinhamento DNA
    stage_banner("6) Alinhamento, Sort e Duplicatas")
    align_dna_for_all(dna_samples, g["threads"], g.get("default_read_type","short"),
                      cleanup=g.get("cleanup", {"remove_sorted_bam": True, "remove_bam_after_cram": True}),
                      use_ds=g.get("downsample_frac", 0.0) > 0)
    step_done("Alinhamento e marcação de duplicatas concluídos")

    # CRAM + Cobertura
    stage_banner("7) CRAM e Cobertura (mosdepth)")
    to_cram_and_coverage(g.get("use_cram", True), g["threads"])
    # Remover BAM após CRAM
    if g.get("use_cram", True) and g.get("cleanup", {}).get("remove_bam_after_cram", True):
        for bam in Path("bam").glob("*.mkdup.bam"):
            console.print(f"Removendo {bam.name} (CRAM mantido)", style="dim")
            bam.unlink(missing_ok=True)
            idx = Path(str(bam)+".bai"); idx.unlink(missing_ok=True)
    step_done("CRAM/index e cobertura prontos")

    # Variantes e VEP
    stage_banner("8) Chamadas de Variantes e Anotação")
    if g.get("call_variants", True):
        call_variants(dna_samples, g["threads"], g["mem_gb"])
    if g.get("annotate_vars", True):
        annotate_with_vep(dna_samples, g["threads"])
    step_done("VCF e VEP-VCF gerados")

    # Lista de genes (da referência)
    stage_banner("9) Lista de Genes (Referência)")
    gene_list_from_gtf()
    step_done("Lista de genes da referência pronta")

    # Relatórios de presença de genes por amostra (cobertura por gene)
    stage_banner("10) Presença de genes por amostra (cobertura por gene)")
    gene_presence_reports(dna_samples,
                          min_mean_cov=float(g.get("gene_presence_min_mean_cov", 5.0)),
                          min_breadth_1x=float(g.get("gene_presence_min_breadth_1x", 0.8)),
                          threads=g["threads"])
    print_meta("Genes BED (alvos p/ cobertura)", [Path("genes/genes.bed")])
    pres = sorted(Path("genes").glob("*_gene_presence.tsv"))
    if pres:
        print_meta("Presença de genes (por amostra)", pres)
    step_done("Relatórios de presença de genes gerados")

    # Comparações par-a-par
    stage_banner("11) Comparações par-a-par (vcf)")
    pairwise_comparisons(dna_samples)
    step_done("Relatórios de comparação gerados")

    # Comparações par-a-par (você já chama pairwise_comparisons antes)
    combine_pairwise_stats(dna_samples)

    # Trio: potenciais de novo (filho ≠ pais)
    stage_banner("12) Trio: candidatos de novo")
    trio_denovo_report(dna_samples)
    step_done("Relatório trio gerado")

    # RNA-seq (opcional)
    if g.get("rnaseq", False):
        stage_banner("13) RNA-seq (opcional)")
        rnaseq_pipeline(cfg.get("rna_samples", []), g["threads"], g["assembly_name"])
        step_done("RNA-seq concluído")

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

