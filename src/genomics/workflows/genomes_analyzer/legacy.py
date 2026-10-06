#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse, os, sys, time, yaml, shutil, subprocess as sp
from pathlib import Path
from datetime import datetime
from typing import List, Optional

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

import re, shlex, csv
import math
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn, TransferSpeedColumn

# Detecta se está rodando em terminal interativo ou background/nohup
import sys
is_interactive = sys.stdout.isatty() and sys.stderr.isatty()

# Função para configurar console (será reconfigurado após carregar YAML)
def _configure_console_for_mode(cfg_params=None):
    """Configura console baseado no modo de execução e parâmetros YAML"""
    global console
    
    if is_interactive:
        # Terminal interativo: usa largura automática
        console = Console(highlight=False, emoji=True)
    else:
        # Background/nohup: usa parâmetros configuráveis
        params = cfg_params or {}
        log_width = int(params.get("log_width_background", 180))
        log_emoji = bool(params.get("log_emoji_background", False))
        log_colors = bool(params.get("log_colors_background", False))
        
        console = Console(
            highlight=False, 
            emoji=log_emoji,
            force_terminal=log_colors,
            width=log_width,
            no_color=not log_colors,
            legacy_windows=False
        )

# Configuração inicial (será reconfigurada no main)
console = Console(highlight=False, emoji=True if is_interactive else False, 
                 width=None if is_interactive else 180)

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
    console.print(f"[bold cyan]>[/bold cyan] {' '.join(cmd)}")
    sp.run(cmd, cwd=cwd, env=env, check=check)

def print_cmd(cmd):
    """Função auxiliar para imprimir comandos antes da execução."""
    if isinstance(cmd, list):
        console.print(f"[bold cyan]>[/bold cyan] {' '.join(cmd)}")
    else:
        console.print(f"[bold cyan]>[/bold cyan] {cmd}")

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
    import threading
    import queue
    
    console.print(Panel.fit(Text(label, style="bold yellow"), border_style="yellow"))
    console.print("[bold]>[/bold] conda run -n genomics bash -lc " + cmd_str, style="dim")

    # Captura stderr para evitar bloqueio com WARNINGs
    proc = sp.Popen(
        ["conda", "run", "-n", "genomics", "bash", "-lc", cmd_str], 
        stdout=sp.PIPE, 
        stderr=sp.PIPE, 
        text=True,
        bufsize=1
    )
    
    start = time.time()
    last = start
    last_r = last_w = 0
    
    # Thread para consumir stderr e evitar bloqueio
    stderr_queue = queue.Queue()
    def stderr_reader():
        try:
            while True:
                line = proc.stderr.readline()
                if not line:
                    break
                stderr_queue.put(line)
        except Exception:
            pass
    
    stderr_thread = threading.Thread(target=stderr_reader, daemon=True)
    stderr_thread.start()

    while True:
        rc = proc.poll()
        now = time.time()

        # Processa mensagens de stderr sem bloqueio
        while not stderr_queue.empty():
            try:
                stderr_line = stderr_queue.get_nowait().strip()
                if stderr_line and ("WARNING" in stderr_line or "ERROR" in stderr_line):
                    # Mostra apenas WARNINGs/ERRORs importantes, truncados
                    console.print(f"[yellow]⚠️  {stderr_line[:120]}{'...' if len(stderr_line) > 120 else ''}[/yellow]")
            except queue.Empty:
                break

        if now - last >= heartbeat_sec:
            # tamanho atual do arquivo de saída (se já existir)
            out_bytes = out_watch.stat().st_size if Path(out_watch).exists() else 0
            # I/O agregado da árvore de processos
            try:
                r_tot, w_tot = _sum_proc_io_tree(proc.pid)
                d_r = max(0, r_tot - last_r); d_w = max(0, w_tot - last_w)
                last_r, last_w = r_tot, w_tot
            except Exception:
                # Se falhar ao ler I/O, continua sem mostrar delta
                d_r = d_w = 0

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

    # Aguarda thread de stderr terminar
    stderr_thread.join(timeout=1.0)
    
    # Processa stderr restante
    remaining_stderr = []
    while not stderr_queue.empty():
        try:
            remaining_stderr.append(stderr_queue.get_nowait().strip())
        except queue.Empty:
            break
    
    if proc.returncode != 0:
        # Mostra stderr completo em caso de erro
        if remaining_stderr:
            console.print(f"[red]❌ Stderr final:[/red]")
            for line in remaining_stderr[-10:]:  # últimas 10 linhas
                if line:
                    console.print(f"[dim]{line}[/dim]")
        raise sp.CalledProcessError(proc.returncode, ["conda", "run", "-n", "genomics", "bash","-lc", cmd_str])
    elif remaining_stderr:
        # Mostra resumo de WARNINGs em caso de sucesso
        warnings = [line for line in remaining_stderr if "WARNING" in line]
        if warnings:
            console.print(f"[yellow]⚠️  {len(warnings)} warning(s) gerados (comando executou com sucesso)[/yellow]")

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
            ["conda", "run", "-n", "genomics", "bash","-lc", f"du -sb {shlex.quote(str(path))} 2>/dev/null | cut -f1"],
            capture_output=True, text=True, check=True
        ).stdout.strip()
        return int(out) if out else 0
    except Exception:
        return 0

from .fastq import (
    _find_local_sra,
    _sra_expected_size_bytes,
    ena_fetch_fastqs,
    ena_get_fastq_urls,
    prefetch_with_progress,
)

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
from .fastq import compress_fastqs_with_progress, fasterq_with_progress

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

from .config import normalize_config_schema

# ================ Estimativa de espaço =================

from .space import estimate_outputs_and_warn

# ================== Referências / índices ==================

from .references import _download_and_place, download_refs

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

from .references import install_prebuilt_bwa_index

from .references import limit_reference_to_canonical_if_enabled


from .references import _build_bwa_index_optimized, _build_bwa_mem2_index_optimized, build_indexes
# ================== Entrada SRA / FASTQ ===================

from .fastq import downsample_fastqs, stage_fastqs_from_local, stage_fastqs_from_sra

# ==================== QC / Trimming ====================

from .qc import fastqc_outputs_exist, qc_and_trim

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

# Reexporta as implementações modularizadas; mantém o bloco legado acima apenas
# como referência temporária até a remoção física segura.
from .alignment import align_dna_for_all, align_one_sample, to_cram_and_coverage

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
    run(["conda", "run", "-n", "genomics", "bash","-lc", script])
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
        return ["conda", "run", "-n", "genomics", "bash","-lc", sh]

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
                tail = sp.run(["conda", "run", "-n", "genomics", "bash","-lc", f"tail -n 60 {shlex.quote(str(log_path))}"],
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
                ["conda", "run", "-n", "genomics", "bash","-lc", f"grep -E '^[Ll]ines\\s' -m1 {shlex.quote(str(log_path))} || true"],
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
                            tail = sp.run(["conda", "run", "-n", "genomics", "bash","-lc", f"tail -n 60 {shlex.quote(str(st['log']))}"],
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
                            ["conda", "run", "-n", "genomics", "bash","-lc", f"grep -E '^[Ll]ines\\s' -m1 {shlex.quote(str(st['log']))} || true"],
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
            
            print_cmd(cmd_filter)
            sp.run(cmd_filter, check=True)
            
            # Indexa VCF filtrado
            print_cmd(["tabix", "-p", "vcf", str(vcf_filtered)])
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
        import threading
        import queue
        
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, bufsize=1)
        
        # Thread para consumir stderr e evitar bloqueio
        stderr_queue = queue.Queue()
        stderr_lines = []
        
        def stderr_reader():
            try:
                while True:
                    line = proc.stderr.readline()
                    if not line:
                        break
                    stderr_queue.put(line)
                    stderr_lines.append(line)  # Guarda para debug se necessário
            except Exception:
                pass
        
        stderr_thread = threading.Thread(target=stderr_reader, daemon=True)
        stderr_thread.start()
        
        while True:
            # Verifica se o processo terminou
            rc = proc.poll()
            current_time = time.time()
            
            # Processa mensagens de stderr sem bloqueio
            while not stderr_queue.empty():
                try:
                    stderr_line = stderr_queue.get_nowait().strip()
                    if stderr_line and ("WARNING" in stderr_line or "ERROR" in stderr_line):
                        console.print(f"[yellow]⚠️  VEP: {stderr_line[:120]}{'...' if len(stderr_line) > 120 else ''}[/yellow]")
                except queue.Empty:
                    break
            
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
        
        # Aguarda thread de stderr terminar
        stderr_thread.join(timeout=2.0)
        
        # Verifica se houve erro
        if proc.returncode != 0:
            console.print(f"[red]❌ VEP falhou com código {proc.returncode}[/red]")
            # Mostra últimas linhas de stderr se disponível
            if stderr_lines:
                console.print(f"[red]Últimas mensagens de erro VEP:[/red]")
                for line in stderr_lines[-10:]:
                    if line.strip():
                        console.print(f"[dim]{line.strip()}[/dim]")
            raise sp.CalledProcessError(proc.returncode, cmd)
        elif stderr_lines:
            # Mostra resumo de warnings em caso de sucesso
            warnings = [line for line in stderr_lines if "WARNING" in line]
            if warnings:
                console.print(f"[yellow]⚠️  VEP concluído com {len(warnings)} warning(s)[/yellow]")
        
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
            print_cmd(["samtools", "faidx", str(ref_fa)])
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
                                print_cmd(["bgzip", "-c", str(tmp_vep)])
                                sp.run(["bgzip", "-c", str(tmp_vep)], stdout=f, check=True)
                            print_cmd(["tabix", "-p", "vcf", str(out_vep)])
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
            print_cmd(["samtools", "faidx", str(ref_fa)])
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
            console.print(f"[cyan]🔀 Fazendo merge de VCFs: {a} vs {b}[/cyan]")
            merge_cmd = ["bcftools","merge","-m","all","-Oz","-o",str(out_merge), str(va), str(vb)]
            console.print(f"[dim]💻 Comando merge:[/dim]")
            console.print(f"[dim]> {' '.join(merge_cmd)}[/dim]")
            run(merge_cmd)
            
            index_cmd = ["bcftools","index","-t",str(out_merge)]
            console.print(f"[dim]💻 Comando index:[/dim]")
            console.print(f"[dim]> {' '.join(index_cmd)}[/dim]")
            run(index_cmd)
            console.print(f"[green]✅ Merge concluído: {out_merge.name}[/green]")
        else:
            console.print(f"[{pair_name}] merge → [bold]SKIP[/bold] (cache ok)")

        # Verifica se arquivo merged é válido antes de query
        if not out_merge.exists():
            console.print(f"[red]❌ Arquivo merge não encontrado: {out_merge}[/red]")
            continue
            
        # Verifica e corrige sample names no arquivo merged
        console.print(f"[cyan]📊 Verificando amostras no arquivo merged...[/cyan]")
        samples_check = sp.run(
            ["bcftools","query","-l",str(out_merge)],
            capture_output=True, text=True, check=False
        )
        
        if samples_check.returncode != 0:
            console.print(f"[red]❌ Erro ao verificar amostras no arquivo merged[/red]")
            continue
            
        available_samples_raw = samples_check.stdout.strip().split('\n') if samples_check.stdout.strip() else []
        console.print(f"[dim]📋 Sample names no VCF: {', '.join(available_samples_raw)}[/dim]")
        
        # Função para mapear sample names problemáticos
        def _map_sample_name(expected_name, available_names):
            """Mapeia nome esperado para nome disponível no VCF"""
            # Busca exata
            if expected_name in available_names:
                return expected_name
            
            # Busca por substring (NA12878 em bam/NA12878.mkdup.bam)
            for avail_name in available_names:
                if expected_name in avail_name:
                    return avail_name
            
            # Busca por padrão de arquivo
            for avail_name in available_names:
                if avail_name.endswith(f'{expected_name}.mkdup.bam') or avail_name.endswith(f'{expected_name}.bam'):
                    return avail_name
            
            return None
        
        # Mapeia nomes das amostras
        mapped_a = _map_sample_name(a, available_samples_raw)
        mapped_b = _map_sample_name(b, available_samples_raw)
        
        if not mapped_a or not mapped_b:
            console.print(f"[yellow]⚠️  Não foi possível mapear amostras:[/yellow]")
            console.print(f"[yellow]   Esperado: {a}, {b}[/yellow]")
            console.print(f"[yellow]   Disponível: {', '.join(available_samples_raw)}[/yellow]")
            console.print(f"[yellow]Pulando comparação {a} vs {b}[/yellow]")
            continue
        
        console.print(f"[green]✅ Mapeamento: {a}→{mapped_a}, {b}→{mapped_b}[/green]")

        # query e métricas usando nomes mapeados
        console.print(f"[cyan]🔍 Extraindo genótipos para comparação...[/cyan]")
        query_cmd = ["bcftools","query","-s",f"{mapped_a},{mapped_b}","-f","%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",str(out_merge)]
        console.print(f"[dim]💻 Comando query:[/dim]")
        console.print(f"[dim]> {' '.join(query_cmd)}[/dim]")
        
        try:
            q = sp.run(query_cmd, capture_output=True, text=True, check=True).stdout.splitlines()
            console.print(f"[green]✅ Query concluído: {len(q):,} variantes extraídas[/green]")
        except sp.CalledProcessError as e:
            console.print(f"[red]❌ bcftools query falhou com código {e.returncode}[/red]")
            console.print(f"[red]Comando: {' '.join(query_cmd)}[/red]")
            if hasattr(e, 'stderr') and e.stderr:
                console.print(f"[red]Erro: {e.stderr}[/red]")
            continue

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

        # Interpretação dos resultados
        def _interpret_relationship(ibs0_pct, share_pct, exact_pct, sample_a, sample_b):
            """Interpreta métricas de comparação genética"""
            if ibs0_pct < 0.001:  # <0.1%
                if exact_pct > 0.7:  # >70%
                    return "👨‍👩‍👧 Relação parental forte"
                elif exact_pct > 0.6:  # >60%
                    return "👥 Relação familiar provável"
                else:
                    return "🤝 Mesma população"
            elif ibs0_pct < 0.01:  # <1%
                return "👨‍👩‍👧 Família ou população próxima"
            elif ibs0_pct < 0.05:  # <5%
                return "🌍 Mesma ancestralidade"
            else:
                return "🌐 Populações diferentes"
        
        relationship = _interpret_relationship(frac_ibs0, frac_share, frac_exact, a, b)
        
        console.print(Panel.fit(
            f"[bold]Comparação Genética: {a} vs {b}[/bold]\n"
            f"📊 Variantes analisadas: {stats['sites_both_called']:,} (ambas amostras)\n"
            f"🧬 Genótipos idênticos: {frac_exact:.2%} ({stats['geno_exact_match']:,} sites)\n"
            f"🤝 Compartilham ≥1 alelo: {frac_share:.2%} ({stats['share_allele']:,} sites)\n"
            f"🚫 IBS0 (zero alelos): {frac_ibs0:.2%} ({stats['ibs0']:,} sites)\n"
            f"🎯 Interpretação: {relationship}",
            border_style="green"
        ))
        
        # Resumo compacto para logs
        console.print(f"[bold cyan]Resumo {pair_name}:[/bold cyan] "
                      f"IBS0={frac_ibs0:.2%}, share≥1={frac_share:.2%}, exact={frac_exact:.2%} • {relationship}")

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

def _print_trio_interpretation_report(total, kept, gt_debug, rejected_counts, 
                                     min_dp_child, min_dp_par, max_par_alt, 
                                     min_ab_het, max_ab_het, min_ab_hom):
    # Imprime relatório interpretativo detalhado dos resultados da análise trio de novo
    console.print("\n" + "="*100)
    console.print(f"[bold blue]🧬 RELATÓRIO INTERPRETATIVO - ANÁLISE TRIO DE NOVO[/bold blue]")
    console.print("="*100)
    
    # 1. Resumo dos resultados
    console.print(f"\n[bold green]📊 RESUMO DOS RESULTADOS:[/bold green]")
    console.print(f"[green]✅ Candidatos de novo identificados: {kept:,}[/green]")
    console.print(f"[dim]   De um total de {total:,} variantes processadas[/dim]")
    console.print(f"[dim]   Taxa de candidatos: {(kept/total*100):.2f}%[/dim]")
    
    # 2. Interpretação biológica
    console.print(f"\n[bold yellow]🔬 INTERPRETAÇÃO BIOLÓGICA:[/bold yellow]")
    if kept > 10000:
        console.print(f"[yellow]⚠️  ALTO número de candidatos ({kept:,}) - provavelmente muitos falsos positivos[/yellow]")
        console.print(f"[dim]   • Genomas humanos típicos: ~50-100 mutações de novo verdadeiras[/dim]")
        console.print(f"[dim]   • Maioria dos candidatos resulta de dados missing nos pais[/dim]")
        console.print(f"[dim]   • Recomenda-se filtros mais rigorosos para análise downstream[/dim]")
    elif kept > 1000:
        console.print(f"[yellow]⚠️  MODERADO número de candidatos ({kept:,}) - filtros adicionais recomendados[/yellow]")
    elif kept > 100:
        console.print(f"[green]✅ Número razoável de candidatos ({kept:,}) - compatível com dados reais[/green]")
    else:
        console.print(f"[blue]ℹ️  Poucos candidatos ({kept:,}) - filtros podem estar muito restritivos[/blue]")
    
    # 3. Análise da qualidade dos dados
    console.print(f"\n[bold cyan]📈 QUALIDADE DOS DADOS:[/bold cyan]")
    
    # Análise do filho
    child_total = sum(gt_debug[k] for k in gt_debug if k.startswith('child_'))
    if child_total > 0:
        child_missing_pct = (gt_debug['child_missing'] / child_total) * 100
        child_het_pct = (gt_debug['child_het'] / child_total) * 100
        child_hom_alt_pct = (gt_debug['child_hom_alt'] / child_total) * 100
        
        console.print(f"[cyan]👶 Filho (NA12878):[/cyan]")
        console.print(f"[dim]   • Missing: {gt_debug['child_missing']:,} ({child_missing_pct:.1f}%) - {'Alto' if child_missing_pct > 40 else 'Normal' if child_missing_pct > 20 else 'Baixo'}[/dim]")
        console.print(f"[dim]   • Heterozigotos: {gt_debug['child_het']:,} ({child_het_pct:.1f}%) - {'Normal' if 35 <= child_het_pct <= 45 else 'Fora do esperado'}[/dim]")
        console.print(f"[dim]   • Homozigotos alt: {gt_debug['child_hom_alt']:,} ({child_hom_alt_pct:.1f}%)[/dim]")
    
    # Análise dos pais
    parents_total = sum(gt_debug[k] for k in gt_debug if k.startswith('parents_'))
    if parents_total > 0:
        parents_missing_pct = (gt_debug['parents_both_missing'] / parents_total) * 100
        parents_mixed_pct = (gt_debug['parents_mixed'] / parents_total) * 100
        
        console.print(f"[cyan]👨‍👩‍👧‍👦 Pais (NA12891 + NA12892):[/cyan]")
        console.print(f"[dim]   • Ambos missing: {gt_debug['parents_both_missing']:,} ({parents_missing_pct:.1f}%)[/dim]")
        console.print(f"[dim]   • Mistos (≥1 com variante): {gt_debug['parents_mixed']:,} ({parents_mixed_pct:.1f}%)[/dim]")
        console.print(f"[dim]   • Ambos hom_ref: {gt_debug['parents_both_hom_ref']:,} (0.0%) - {'⚠️  Problemático' if gt_debug['parents_both_hom_ref'] == 0 else 'OK'}[/dim]")
    
    # 4. Análise dos filtros aplicados
    console.print(f"\n[bold magenta]🔧 EFICÁCIA DOS FILTROS:[/bold magenta]")
    total_rejected = sum(rejected_counts.values())
    
    for reason, count in rejected_counts.items():
        if count > 0:
            pct = (count / total) * 100
            reason_pt = {
                "parents_not_hom_ref": "Pais não são hom_ref/missing",
                "parents_low_dp": "Baixa cobertura nos pais", 
                "parents_high_alt_freq": "Alta freq. alternativa nos pais",
                "child_not_nonref": "Filho não tem variante",
                "child_low_dp": "Baixa cobertura no filho",
                "child_bad_ab": "Frequências alélicas inadequadas",
                "parsing_error": "Erros de parsing"
            }.get(reason, reason)
            
            console.print(f"[dim]   • {reason_pt}: {count:,} ({pct:.1f}%)[/dim]")
    
    # 5. Recomendações
    console.print(f"\n[bold green]💡 RECOMENDAÇÕES:[/bold green]")
    
    if kept > 10000:
        console.print(f"[green]🔧 Para reduzir falsos positivos:[/green]")
        console.print(f"[dim]   • Aumentar min_dp_child para 15-20 (atual: {min_dp_child})[/dim]")
        console.print(f"[dim]   • Aumentar min_dp_parents para 10-15 (atual: {min_dp_par})[/dim]")
        console.print(f"[dim]   • Exigir pelo menos um pai com genótipo 0/0 explícito[/dim]")
    
    console.print(f"[green]🧬 Para análise downstream:[/green]")
    console.print(f"[dim]   • Priorizar variantes em genes codificantes[/dim]")
    console.print(f"[dim]   • Filtrar por impacto funcional (VEP annotations)[/dim]")
    console.print(f"[dim]   • Verificar contra bancos de dados (ClinVar, gnomAD)[/dim]")
    console.print(f"[dim]   • Validar candidatos top por Sanger sequencing[/dim]")
    
    console.print(f"[green]🔍 Para validação:[/green]")
    console.print(f"[dim]   • Inspecionar no IGV: variantes com alta cobertura[/dim]")
    console.print(f"[dim]   • Verificar padrões de herança mendeliana[/dim]")
    console.print(f"[dim]   • Considerar CNVs e rearranjos estruturais[/dim]")
    
    # 6. Arquivos gerados
    console.print(f"\n[bold blue]📁 ARQUIVOS GERADOS:[/bold blue]")
    console.print(f"[blue]   • trio/trio_denovo_candidates.tsv - {kept:,} candidatos[/blue]")
    console.print(f"[blue]   • trio/trio_denovo_summary.md - Resumo em markdown[/blue]")
    console.print(f"[blue]   • trio/trio_merged.vcf.gz - VCF trio original[/blue]")
    
    console.print("="*100 + "\n")
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
        console.print(f"[cyan]🧬 Fazendo merge trio: {child} (filho) + {p1}, {p2} (pais)[/cyan]")
        merge_cmd = ["bcftools","merge","-m","all","-Oz","-o",str(merged),
                     str(vcfs[child]), str(vcfs[p1]), str(vcfs[p2])]
        console.print(f"[dim]💻 Comando merge trio:[/dim]")
        console.print(f"[dim]> {' '.join(merge_cmd)}[/dim]")
        run(merge_cmd)
        
        index_cmd = ["bcftools","index","-t",str(merged)]
        console.print(f"[dim]💻 Comando index trio:[/dim]")
        console.print(f"[dim]> {' '.join(index_cmd)}[/dim]")
        run(index_cmd)
        console.print(f"[green]✅ Merge trio concluído: {merged.name}[/green]")
    else:
        console.print(f"Trio merge → [bold]SKIP[/bold] (cache ok)")

    # Verifica e mapeia sample names no trio merged
    console.print(f"[cyan]📊 Verificando sample names no trio merged...[/cyan]")
    samples_check = sp.run(
        ["bcftools","query","-l",str(merged)],
        capture_output=True, text=True, check=False
    )
    
    if samples_check.returncode != 0:
        console.print(f"[red]❌ Erro ao verificar amostras no trio merged[/red]")
        return
        
    available_samples = samples_check.stdout.strip().split('\n') if samples_check.stdout.strip() else []
    console.print(f"[dim]📋 Sample names no trio VCF: {', '.join(available_samples)}[/dim]")
    
    # Função para mapear sample names (reutilizada da etapa 11)
    def _map_trio_sample_name(expected_name, available_names):
        """Mapeia nome esperado para nome disponível no VCF trio"""
        if expected_name in available_names:
            return expected_name
        for avail_name in available_names:
            if expected_name in avail_name:
                return avail_name
        return None
    
    # Mapeia nomes do trio
    mapped_child = _map_trio_sample_name(child, available_samples)
    mapped_p1 = _map_trio_sample_name(p1, available_samples)
    mapped_p2 = _map_trio_sample_name(p2, available_samples)
    
    if not all([mapped_child, mapped_p1, mapped_p2]):
        console.print(f"[red]❌ Não foi possível mapear todas as amostras do trio:[/red]")
        console.print(f"[red]   Esperado: {child}, {p1}, {p2}[/red]")
        console.print(f"[red]   Disponível: {', '.join(available_samples)}[/red]")
        return
    
    console.print(f"[green]✅ Mapeamento trio: {child}→{mapped_child}, {p1}→{mapped_p1}, {p2}→{mapped_p2}[/green]")

    # thresholds (como antes) - valores mais permissivos para dados reais
    g = cfg_global.get("general", {})
    min_dp_child = int(g.get("trio_min_dp_child", 5))    # Reduzido de 15 para 5
    min_dp_par   = int(g.get("trio_min_dp_parents", 5))  # Reduzido de 15 para 5  
    min_gq       = int(g.get("trio_min_gq", 20))  # Não usado (GQ não disponível no VCF)
    min_ab_het   = float(g.get("trio_min_ab_het", 0.20))  # Mais permissivo: 0.20-0.80
    max_ab_het   = float(g.get("trio_max_ab_het", 0.80))
    min_ab_hom   = float(g.get("trio_min_ab_hom", 0.85))  # Mais permissivo: 0.85
    max_par_alt  = float(g.get("trio_max_parent_alt_frac", 0.05))  # Mais permissivo: 0.05
    
    console.print(f"[yellow]⚠️  GQ (Genotype Quality) não disponível no VCF - filtros de qualidade baseados apenas em DP[/yellow]")
    console.print(f"[cyan]🔧 Filtros ajustados para dados reais: pais podem ser 0/0 ou missing, DP mínimo reduzido[/cyan]")

    # Diagnóstico do arquivo trio antes de query
    console.print(f"[cyan]🔍 Diagnosticando arquivo trio merged...[/cyan]")
    
    # Verifica total de variantes no arquivo
    total_variants = sp.run(
        ["bcftools", "view", "-H", str(merged)],
        capture_output=True, text=True, check=False
    )
    
    if total_variants.returncode == 0:
        total_count = len(total_variants.stdout.strip().split('\n')) if total_variants.stdout.strip() else 0
        console.print(f"[dim]📊 Total de variantes no trio: {total_count:,}[/dim]")
    else:
        console.print(f"[yellow]⚠️  Não foi possível contar variantes totais[/yellow]")
        total_count = 0
    
    # Verifica variantes PASS
    pass_variants = sp.run(
        ["bcftools", "view", "-f", "PASS", "-H", str(merged)],
        capture_output=True, text=True, check=False
    )
    
    if pass_variants.returncode == 0:
        pass_count = len(pass_variants.stdout.strip().split('\n')) if pass_variants.stdout.strip() else 0
        console.print(f"[dim]📊 Variantes PASS no trio: {pass_count:,}[/dim]")
    else:
        console.print(f"[yellow]⚠️  Não foi possível contar variantes PASS[/yellow]")
        pass_count = 0
    
    # Verifica filtros disponíveis no VCF
    filter_info = sp.run(
        ["bcftools", "view", "-h", str(merged)],
        capture_output=True, text=True, check=False
    )
    
    if filter_info.returncode == 0:
        filters_found = []
        for line in filter_info.stdout.split('\n'):
            if line.startswith('##FILTER='):
                filter_match = re.search(r'ID=([^,>]+)', line)
                if filter_match:
                    filters_found.append(filter_match.group(1))
        console.print(f"[dim]📋 Filtros disponíveis: {', '.join(filters_found) if filters_found else 'nenhum'}[/dim]")
    
    if pass_count == 0:
        console.print(f"[yellow]⚠️  Nenhuma variante com FILTER=PASS - usando todas as variantes...[/yellow]")
        # Query sem filtro PASS (mais permissivo) - removendo GQ que não existe
        trio_query_cmd = f"bcftools query -s {mapped_child},{mapped_p1},{mapped_p2} -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP:%AD]\\n' {shlex.quote(str(merged))}"
    else:
        console.print(f"[green]✅ {pass_count:,} variantes passaram no filtro PASS[/green]")
        # Query com filtro PASS - removendo GQ que não existe
        trio_query_cmd = (
            f"bcftools view -f PASS {shlex.quote(str(merged))} | "
            f"bcftools query -s {mapped_child},{mapped_p1},{mapped_p2} -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP:%AD]\\n'"
        )
    
    console.print(f"[cyan]🔍 Extraindo genótipos do trio para análise de novo...[/cyan]")
    console.print(f"[dim]💻 Comando trio query:[/dim]")
    console.print(f"[dim]> {trio_query_cmd}[/dim]")
    
    try:
        q = sp.run(
            ["conda", "run", "-n", "genomics", "bash","-lc", trio_query_cmd],
            capture_output=True, text=True, check=True
        ).stdout.splitlines()
        console.print(f"[green]✅ Trio query concluído: {len(q):,} variantes extraídas[/green]")
        
        # Mostra exemplo das primeiras variantes para debug
        if len(q) > 0:
            console.print(f"[dim]🔍 Exemplo das primeiras 3 variantes:[/dim]")
            for i, line in enumerate(q[:3]):
                if line.strip():
                    console.print(f"[dim]   {i+1}: {line[:100]}{'...' if len(line) > 100 else ''}[/dim]")
        
        if len(q) == 0:
            console.print(f"[yellow]⚠️  Query retornou 0 variantes - verificando possíveis causas:[/yellow]")
            console.print(f"[yellow]   • Total no arquivo: {total_count:,}[/yellow]")
            console.print(f"[yellow]   • Variantes PASS: {pass_count:,}[/yellow]")
            console.print(f"[yellow]   • Sample names: {mapped_child}, {mapped_p1}, {mapped_p2}[/yellow]")
            
            # Tenta query simples para debug
            debug_cmd = f"bcftools query -f '%CHROM\\t%POS[\\t%SAMPLE=%GT]\\n' {shlex.quote(str(merged))} | head -5"
            console.print(f"[dim]🔍 Debug query:[/dim]")
            console.print(f"[dim]> {debug_cmd}[/dim]")
            
            debug_result = sp.run(["conda", "run", "-n", "genomics", "bash", "-lc", debug_cmd], capture_output=True, text=True, check=False)
            if debug_result.stdout.strip():
                console.print(f"[dim]📋 Exemplo de dados no arquivo:[/dim]")
                for line in debug_result.stdout.strip().split('\n')[:3]:
                    console.print(f"[dim]   {line}[/dim]")
            
    except sp.CalledProcessError as e:
        console.print(f"[red]❌ bcftools trio query falhou com código {e.returncode}[/red]")
        console.print(f"[red]Comando: {trio_query_cmd}[/red]")
        if hasattr(e, 'stderr') and e.stderr:
            console.print(f"[red]Erro: {e.stderr}[/red]")
        return

    out_tsv = Path("trio")/"trio_denovo_candidates.tsv"
    kept = 0; total = 0
    
    # Contadores de diagnóstico para entender onde as variantes são filtradas
    rejected_counts = {
        "parents_not_hom_ref": 0,
        "parents_low_dp": 0, 
        "parents_high_alt_freq": 0,
        "child_not_nonref": 0,
        "child_low_dp": 0,
        "child_bad_ab": 0,
        "parsing_error": 0
    }
    
    # Contadores para debug de genótipos
    gt_debug = {
        "child_missing": 0,
        "child_hom_ref": 0,
        "child_het": 0,
        "child_hom_alt": 0,
        "parents_both_missing": 0,
        "parents_both_hom_ref": 0,
        "parents_mixed": 0
    }
    
    with open(out_tsv, "w") as out:
        out.write("chrom\tpos\tref\talt\tchild_GT\tchild_DP\tchild_GQ\tchild_AB\tp1_GT\tp1_DP\tp1_GQ\tp1_AB\tp2_GT\tp2_DP\tp2_GQ\tp2_AB\n")
        for line in q:
            if not line.strip(): continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:  # precisa ter 3 blocos de amostra
                rejected_counts["parsing_error"] += 1
                continue
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            sa, sb, sc = parts[4], parts[5], parts[6]
            def split_block(b):
                xx = b.split(":")
                GT = xx[0] if len(xx)>0 else "."
                DP = int(xx[1]) if len(xx)>1 and xx[1].isdigit() else None
                AD = xx[2] if len(xx)>2 else ""
                # GQ não disponível no VCF - será None para todos
                GQ = None
                return GT, DP, GQ, AD

            cGT,cDP,cGQ,cAD = split_block(sa)
            p1GT,p1DP,p1GQ,p1AD = split_block(sb)
            p2GT,p2DP,p2GQ,p2AD = split_block(sc)

            gt_c  = _parse_gt(cGT); gt_p1 = _parse_gt(p1GT); gt_p2 = _parse_gt(p2GT)
            _,_,ab_c  = _parse_ad(cAD);  _,_,ab_p1 = _parse_ad(p1AD);  _,_,ab_p2 = _parse_ad(p2AD)
           
            # Debug de genótipos
            if gt_c is None:
                gt_debug["child_missing"] += 1
            elif _is_hom_ref(gt_c):
                gt_debug["child_hom_ref"] += 1
            elif len(set(gt_c)) > 1:
                gt_debug["child_het"] += 1
            else:
                gt_debug["child_hom_alt"] += 1
                
            if gt_p1 is None and gt_p2 is None:
                gt_debug["parents_both_missing"] += 1
            elif _is_hom_ref(gt_p1) and _is_hom_ref(gt_p2):
                gt_debug["parents_both_hom_ref"] += 1
            else:
                gt_debug["parents_mixed"] += 1

            # pais - aceitar hom_ref (0/0) OU missing (None)
            p1_ok = _is_hom_ref(gt_p1) or gt_p1 is None
            p2_ok = _is_hom_ref(gt_p2) or gt_p2 is None
            
            if not (p1_ok and p2_ok):
                rejected_counts["parents_not_hom_ref"] += 1
                continue
            # Filtro DP para pais - mais flexível para dados missing
            p1_dp_ok = p1DP is None or p1DP >= min_dp_par
            p2_dp_ok = p2DP is None or p2DP >= min_dp_par
            if not (p1_dp_ok and p2_dp_ok):
                rejected_counts["parents_low_dp"] += 1
                continue
            # Removido filtro GQ para pais (não disponível no VCF)
            if ab_p1 is not None and ab_p1 > max_par_alt: 
                rejected_counts["parents_high_alt_freq"] += 1
                continue
            if ab_p2 is not None and ab_p2 > max_par_alt: 
                rejected_counts["parents_high_alt_freq"] += 1
                continue

            # filho
            if not _is_nonref(gt_c): 
                rejected_counts["child_not_nonref"] += 1
                continue
            # Filtro DP para filho - aceitar se DP >= min ou se missing mas com variante clara
            if cDP is not None and cDP < min_dp_child:
                rejected_counts["child_low_dp"] += 1
                continue
            # Removido filtro GQ para filho (não disponível no VCF)
            if ab_c is not None:
                if len(set(gt_c)) > 1:
                    if not (min_ab_het <= ab_c <= max_ab_het): 
                        rejected_counts["child_bad_ab"] += 1
                        continue
                else:
                    if ab_c < min_ab_hom: 
                        rejected_counts["child_bad_ab"] += 1
                        continue

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

    # Relatório de diagnóstico
    console.print(f"[cyan]📊 Diagnóstico de genótipos processados:[/cyan]")
    console.print(f"[dim]   Filho - Missing: {gt_debug['child_missing']:,}, Hom_ref: {gt_debug['child_hom_ref']:,}, Het: {gt_debug['child_het']:,}, Hom_alt: {gt_debug['child_hom_alt']:,}[/dim]")
    console.print(f"[dim]   Pais - Ambos missing: {gt_debug['parents_both_missing']:,}, Ambos hom_ref: {gt_debug['parents_both_hom_ref']:,}, Mistos: {gt_debug['parents_mixed']:,}[/dim]")
    
    console.print(f"[cyan]📊 Diagnóstico de filtros trio de novo:[/cyan]")
    console.print(f"[dim]   Total processadas: {total:,}[/dim]")
    console.print(f"[dim]   Candidatos finais: {kept:,}[/dim]")
    console.print(f"[dim]   Rejeitadas por:[/dim]")
    for reason, count in rejected_counts.items():
        if count > 0:
            console.print(f"[dim]     • {reason}: {count:,}[/dim]")
    
    # Mostra configuração dos filtros
    console.print(f"[dim]📋 Configuração dos filtros:[/dim]")
    console.print(f"[dim]   • min_dp_child: {min_dp_child}[/dim]")
    console.print(f"[dim]   • min_dp_parents: {min_dp_par}[/dim]")
    console.print(f"[dim]   • max_parent_alt_frac: {max_par_alt}[/dim]")
    console.print(f"[dim]   • AB het range: {min_ab_het}-{max_ab_het}[/dim]")
    console.print(f"[dim]   • AB hom min: {min_ab_hom}[/dim]")
    
    console.print(f"[bold cyan]Trio de novo[/bold cyan] → {out_tsv}  (n={kept})")
    
    # Relatório interpretativo detalhado
    _print_trio_interpretation_report(total, kept, gt_debug, rejected_counts, 
                                    min_dp_child, min_dp_par, max_par_alt, 
                                    min_ab_het, max_ab_het, min_ab_hom)

# =================== Paternidade (SNPs, trio) ===================

def _gt_to_alt_count(gt_tuple) -> int | None:
    # Converte tuple de alelos (0/1) para contagem de ALT (0..2)
    if gt_tuple is None:
        return None
    return sum(1 for a in gt_tuple if a == 1)

def _passes_qc_for_paternity(gt: str, dp: int|None, ad: tuple[int,int]|None, is_child: bool, g: dict) -> bool:
    # Reaproveita thresholds do trio; AD=(ref,alt)
    if not gt or gt == "." or gt.startswith("."):
        return False
    min_dp = int(g.get("trio_min_dp_child",5) if is_child else g.get("trio_min_dp_parents",5))
    if dp is not None and dp < min_dp:
        return False
    # Se não houver AD, aceita — confiando no GT
    if ad is None or (ad[0]+ad[1]) == 0:
        return True
    ab = ad[1] / max(1, ad[0]+ad[1])
    min_ab_het = float(g.get("trio_min_ab_het",0.20))
    max_ab_het = float(g.get("trio_max_ab_het",0.80))
    min_ab_hom = float(g.get("trio_min_ab_hom",0.85))
    parsed = _parse_gt(gt)
    altc = _gt_to_alt_count(parsed)
    if altc == 1:
        return (min_ab_het <= ab <= max_ab_het)
    if altc == 0:
        return (ab <= (1.0 - min_ab_hom))
    if altc == 2:
        return (ab >= min_ab_hom)
    return True

def _clip_freq(q: float) -> float:
    # Evita 0/1 exatos
    q = max(1e-4, min(1.0-1e-4, q))
    return q

def _load_vep_csq_index(vep_vcf: Path) -> dict | None:
    """Lê o header do VEP-VCF e retorna um mapa {field_name: index} da tag CSQ.
    Procura por gnomADg_AF, gnomAD_AF, AF nessa ordem.
    """
    cmd = f"bcftools view -h {shlex.quote(str(vep_vcf))} | grep '^##INFO=<ID=CSQ' | tail -n1"
    res = sp.run(["conda","run","-n","genomics","bash","-lc", cmd], stdout=sp.PIPE, stderr=sp.STDOUT, text=True, check=False)
    if res.returncode != 0 or not res.stdout.strip():
        return None
    line = res.stdout.strip()
    i = line.find("Format:")
    if i == -1:
        return None
    fmt = line[i+7:].strip().strip('"').strip()
    fields = [f.strip() for f in fmt.split('|')]
    return {name: idx for idx, name in enumerate(fields)}

def _build_af_map_from_vep(dna_samples) -> dict:
    """Constrói um dicionário (chrom,pos,ref,alt) -> AF (gnomAD) a partir dos VEP-VCFs existentes.
    Usa gnomADg_AF ou gnomAD_AF; se ambos ausentes, ignora.
    """
    af_map: dict[tuple[str,str,str,str], float] = {}
    for s in dna_samples:
        sid = s["id"]
        # Tenta primeiro o caminho relativo, depois absoluto
        vep_vcf = Path("vep")/f"{sid}.vep.vcf.gz"
        if not vep_vcf.exists():
            # Tenta caminho absoluto baseado no diretório de trabalho do config
            base_dir = cfg_global.get("storage", {}).get("base_dir", ".")
            vep_vcf = Path(base_dir) / "vep" / f"{sid}.vep.vcf.gz"
        if not vep_vcf.exists():
            continue
        idx = _load_vep_csq_index(vep_vcf)
        if not idx:
            continue
        key_candidates = [k for k in ("gnomADg_AF","gnomAD_AF","AF") if k in idx]
        if not key_candidates:
            continue
        af_key = key_candidates[0]
        # extrai CHROM,POS,REF,ALT,CSQ
        cmd = ["conda","run","-n","genomics","bash","-lc",
               f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' {shlex.quote(str(vep_vcf))}"]
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt, csq = parts[0], parts[1], parts[2], parts[3], parts[4]
            if not csq:
                continue
            best_af = None
            for entry in csq.split(','):
                cols = entry.split('|')
                if len(cols) <= max(idx.values()):
                    continue
                allele = cols[0]
                if allele != alt:
                    continue
                val = cols[idx[af_key]] if idx[af_key] < len(cols) else ""
                try:
                    if val and val != '.':
                        af_val = float(val)
                        best_af = max(best_af, af_val) if best_af is not None else af_val
                except Exception:
                    continue
            if best_af is not None:
                af_map[(chrom, pos, ref, alt)] = max(1e-4, min(1-1e-4, best_af))
        proc.wait()
    return af_map

def _load_samples_from_vcf(vcf_path: Path) -> list[str]:
    cmd = f"bcftools query -l {shlex.quote(str(vcf_path))}"
    res = sp.run(["conda","run","-n","genomics","bash","-lc", cmd], stdout=sp.PIPE, stderr=sp.STDOUT, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"bcftools query -l falhou: {res.stdout}")
    return [x.strip() for x in res.stdout.strip().split("\n") if x.strip()]

def _map_name(expected: str, avail: list[str]) -> str | None:
    if expected in avail:
        return expected
    cand = [x for x in avail if expected.lower() in x.lower()]
    if len(cand) == 1:
        return cand[0]
    return None

def _parse_sample_field(tok: str) -> tuple[str|None, int|None, tuple[int,int]|None]:
    # Espera "%GT:%DP:%AD"
    gt, dp, ad = None, None, None
    parts = tok.split(":")
    if len(parts) >= 1 and parts[0]:
        gt = parts[0]
    if len(parts) >= 2 and parts[1] and parts[1] != ".":
        try: dp = int(parts[1])
        except: dp = None
    if len(parts) >= 3 and parts[2] and parts[2] != ".":
        try:
            ad_parts = parts[2].split(",")
            if len(ad_parts) >= 2:
                ad = (int(ad_parts[0]), int(ad_parts[1]))
        except:
            ad = None
    return gt, dp, ad

def _estimate_af_from_trio(child_alt: int|None, mom_alt: int|None, ap_alt: int|None) -> float:
    counts = [x for x in (child_alt, mom_alt, ap_alt) if x is not None]
    if not counts:
        return 0.5
    alt = sum(counts)
    an = 2 * len(counts)
    return alt / max(1, an)

def _mendel_child_prob(m_alt: int, f_alt: int, c_alt: int) -> float:
    # Probabilidade do filho (ALT count) dado mãe/pai (ALT counts) sob Mendel simples
    def dist(g):
        if g == 0: return {0:1.0}
        if g == 1: return {0:0.5, 1:0.5}
        return {1:1.0}
    dm, df = dist(m_alt), dist(f_alt)
    p = 0.0
    for tm, pm in dm.items():
        for tf, pf in df.items():
            if tm + tf == c_alt:
                p += pm * pf
    return p
def _paternity_for_pair(merged: Path, child_name:str, mom_name:str, ap_name:str, g:dict, prior:float, eps:float, require_pass:bool) -> dict:
    # Formato por-site com bloco por amostra; sem TAB extra ao final
    # Não requer INFO/AF (nem sempre presente nos VCFs do pipeline)
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT:%DP:%AD]\\n"

    def _run_and_collect(pass_only: bool):
        if pass_only:
            cmd_str = (
                f"bcftools view -f PASS {shlex.quote(str(merged))} | "
                f"bcftools query --force-samples -s {child_name},{mom_name},{ap_name} -f '{fmt}'"
            )
        else:
            cmd_str = (
                f"bcftools query --force-samples -s {child_name},{mom_name},{ap_name} -f '{fmt}' {shlex.quote(str(merged))}"
            )
        proc = sp.Popen(["conda","run","-n","genomics","bash","-lc", cmd_str], stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        return proc

    # 1ª tentativa: com PASS (se configurado)
    proc = _run_and_collect(require_pass)
    log10_sum = 0.0
    n_used = n_excluded = 0
    rows = []
    read_any = False
    while True:
        line = proc.stdout.readline()
        if not line:
            break
        read_any = True
        line = line.rstrip("\n")
        cols = line.split("\t")
        if len(cols) < 4:
            continue
        chrom, pos, ref, alt = cols[0], cols[1], cols[2], cols[3]
        # Apenas SNP bialélico
        if len(ref) != 1 or len(alt) != 1:
            n_excluded += 1
            continue
        samp_toks = cols[4:]
        if len(samp_toks) < 3:
            continue
        gt_c, dp_c, ad_c = _parse_sample_field(samp_toks[0].strip())
        gt_m, dp_m, ad_m = _parse_sample_field(samp_toks[1].strip())
        gt_f, dp_f, ad_f = _parse_sample_field(samp_toks[2].strip())
        # QC
        if (not _passes_qc_for_paternity(gt_c, dp_c, ad_c, True, g)) or (not _passes_qc_for_paternity(gt_m, dp_m, ad_m, False, g)) or (not _passes_qc_for_paternity(gt_f, dp_f, ad_f, False, g)):
            n_excluded += 1
            continue
        c_alt = _gt_to_alt_count(_parse_gt(gt_c))
        m_alt = _gt_to_alt_count(_parse_gt(gt_m))
        f_alt = _gt_to_alt_count(_parse_gt(gt_f))
        if c_alt is None or m_alt is None or f_alt is None:
            n_excluded += 1
            continue
        # Filtro leve: ignorar trios todos heterozigotos (pouco informativos)
        if bool(g.get("paternity_skip_all_hets", False)) and (c_alt == 1 and m_alt == 1 and f_alt == 1):
            n_excluded += 1
            continue
        # AF: tenta VEP (gnomAD), senão estima a partir do trio
        q = None
        try:
            key = (chrom, pos, ref, alt)
            if "_paternity_af_map" in globals() and key in globals()["_paternity_af_map"]:
                q = globals()["_paternity_af_map"][key]
        except Exception:
            q = None
        if q is None:
            q = _estimate_af_from_trio(c_alt, m_alt, f_alt)
        q = _clip_freq(q); p = 1.0 - q
        # P(C|H1) e P(C|H2)
        p_h1 = _mendel_child_prob(m_alt, f_alt, c_alt)
        if p_h1 == 0.0:
            p_h1 = eps
        geno_probs = {0: p*p, 1: 2*p*q, 2: q*q}
        p_h2 = 0.0
        for fa_alt, pg in geno_probs.items():
            p_h2 += pg * _mendel_child_prob(m_alt, fa_alt, c_alt)
        p_h2 = max(eps, p_h2)
        pi = p_h1 / p_h2
        log10_pi = math.log10(pi) if pi > 0 else math.log10(eps)
        log10_sum += log10_pi
        n_used += 1
        rows.append((chrom,pos,ref,alt,gt_c,gt_m,gt_f,dp_c,dp_m,dp_f,
                     ad_c[0] if ad_c else ".", ad_c[1] if ad_c else ".",
                     ad_m[0] if ad_m else ".", ad_m[1] if ad_m else ".",
                     ad_f[0] if ad_f else ".", ad_f[1] if ad_f else ".",
                     f"{q:.6f}", f"{p_h1:.6g}", f"{p_h2:.6g}", f"{pi:.6g}", f"{log10_pi:.6g}"))
    stderr = proc.stderr.read() if proc.stderr else ""
    retcode = proc.wait()

    # Fallback: se nada lido com PASS, tenta sem PASS
    if require_pass and not read_any:
        console.print("[yellow]⚠️  0 linhas com FILTER=PASS; tentando sem PASS...[/yellow]")
        proc2 = _run_and_collect(False)
        read_any = False
        while True:
            line = proc2.stdout.readline()
            if not line:
                break
            read_any = True
            line = line.rstrip("\n")
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            chrom, pos, ref, alt = cols[0], cols[1], cols[2], cols[3]
            if len(ref) != 1 or len(alt) != 1:
                n_excluded += 1
                continue
            samp_toks = cols[4:]
            if len(samp_toks) < 3:
                continue
            gt_c, dp_c, ad_c = _parse_sample_field(samp_toks[0].strip())
            gt_m, dp_m, ad_m = _parse_sample_field(samp_toks[1].strip())
            gt_f, dp_f, ad_f = _parse_sample_field(samp_toks[2].strip())
            if (not _passes_qc_for_paternity(gt_c, dp_c, ad_c, True, g)) or (not _passes_qc_for_paternity(gt_m, dp_m, ad_m, False, g)) or (not _passes_qc_for_paternity(gt_f, dp_f, ad_f, False, g)):
                n_excluded += 1
                continue
            c_alt = _gt_to_alt_count(_parse_gt(gt_c))
            m_alt = _gt_to_alt_count(_parse_gt(gt_m))
            f_alt = _gt_to_alt_count(_parse_gt(gt_f))
            if c_alt is None or m_alt is None or f_alt is None:
                n_excluded += 1
                continue
            if bool(g.get("paternity_skip_all_hets", False)) and (c_alt == 1 and m_alt == 1 and f_alt == 1):
                n_excluded += 1
                continue
            # AF: tenta VEP (gnomAD), senão estima a partir do trio
            q = None
            try:
                key = (chrom, pos, ref, alt)
                if "_paternity_af_map" in globals() and key in globals()["_paternity_af_map"]:
                    q = globals()["_paternity_af_map"][key]
            except Exception:
                q = None
            if q is None:
                q = _estimate_af_from_trio(c_alt, m_alt, f_alt)
            q = _clip_freq(q); p = 1.0 - q
            p_h1 = _mendel_child_prob(m_alt, f_alt, c_alt)
            if p_h1 == 0.0:
                p_h1 = eps
            geno_probs = {0: p*p, 1: 2*p*q, 2: q*q}
            p_h2 = 0.0
            for fa_alt, pg in geno_probs.items():
                p_h2 += pg * _mendel_child_prob(m_alt, fa_alt, c_alt)
            p_h2 = max(eps, p_h2)
            pi = p_h1 / p_h2
            log10_pi = math.log10(pi) if pi > 0 else math.log10(eps)
            log10_sum += log10_pi
            n_used += 1
            rows.append((chrom,pos,ref,alt,gt_c,gt_m,gt_f,dp_c,dp_m,dp_f,
                         ad_c[0] if ad_c else ".", ad_c[1] if ad_c else ".",
                         ad_m[0] if ad_m else ".", ad_m[1] if ad_m else ".",
                         ad_f[0] if ad_f else ".", ad_f[1] if ad_f else ".",
                         f"{q:.6f}", f"{p_h1:.6g}", f"{p_h2:.6g}", f"{pi:.6g}", f"{log10_pi:.6g}"))
        stderr2 = proc2.stderr.read() if proc2.stderr else ""
        proc2.wait()
        if not read_any:
            console.print("[yellow]⚠️  0 linhas mesmo sem PASS. Verifique trio_merged.vcf.gz.[/yellow]")
            if stderr2:
                console.print(f"[dim]{stderr2[:400]}[/dim]")
    else:
        if not read_any:
            console.print("[yellow]⚠️  0 linhas lidas do bcftools.")
            if stderr:
                console.print(f"[dim]{stderr[:400]}[/dim]")

    lr_log10 = log10_sum
    
    # Cálculo numericamente estável do posterior para evitar overflow/underflow
    # Fórmula original: posterior = (LR * prior) / (LR * prior + (1-prior))
    # Versão estável: posterior = 1 / (1 + (1-prior)/(LR * prior))
    #                           = 1 / (1 + (1-prior) * 10^(-lr_log10) / prior)
    
    if lr_log10 > 100:  # LR muito grande (evidência forte para paternidade)
        posterior = 1.0  # Praticamente 100%
        lr = float('inf')  # Representa LR infinito
    elif lr_log10 < -100:  # LR muito pequeno (evidência forte contra paternidade)
        posterior = 0.0  # Praticamente 0%
        lr = 0.0
    else:
        # Caso intermediário: usa cálculo logarítmico estável
        # log(posterior) = log(LR * prior) - log(LR * prior + (1-prior))
        # Para evitar overflow, reescrevemos como:
        # posterior = 1 / (1 + (1-prior)/(LR * prior))
        #           = 1 / (1 + (1-prior)/prior * 10^(-lr_log10))
        
        log10_term = -lr_log10 + math.log10((1.0 - prior) / prior)
        
        if log10_term < -15:  # Termo desprezível
            posterior = 1.0
        elif log10_term > 15:  # Termo dominante
            posterior = 0.0
        else:
            term = 10**log10_term
            posterior = 1.0 / (1.0 + term)
        
        # Calcula LR apenas se necessário e seguro
        if abs(lr_log10) < 300:  # Limite seguro para double precision
            lr = 10**lr_log10
        else:
            lr = float('inf') if lr_log10 > 0 else 0.0
    return {
        "n_used": n_used, "n_excluded": n_excluded,
        "lr_log10": lr_log10, "lr": lr, "posterior": posterior,
        "rows": rows
    }

def paternity_analysis(dna_samples, prior: float|None=None, eps: float|None=None, require_pass: bool|None=None):
    g = cfg_global.get("general", {})
    prior = float(g.get("paternity_prior", 0.5) if prior is None else prior)
    eps = float(g.get("paternity_epsilon", 1e-3) if eps is None else eps)
    require_pass = bool(g.get("paternity_require_pass", True) if require_pass is None else require_pass)

    # Opcional: mapa de AFs populacionais via VEP (pode ser custoso). Controlado por YAML.
    af_map = {}
    use_vep_af = bool(g.get("paternity_use_vep_af", False))
    if use_vep_af:
        console.print("[cyan]Carregando AFs populacionais do VEP (pode consumir RAM/tempo)...[/cyan]")
        af_map = _build_af_map_from_vep(dna_samples)
        if len(af_map) > 0:
            console.print(f"[green]✅ Carregadas {len(af_map):,} frequências do VEP/gnomAD[/green]")
        else:
            console.print("[yellow]⚠️ Nenhuma frequência VEP carregada! Usando estimativas do trio.[/yellow]")
    else:
        console.print("[dim]VEP AF desabilitado, usando estimativas do trio.[/dim]")

    # Garante trio merged
    vcfs = { s["id"]: Path("vcf")/f"{s['id']}.vcf.gz" for s in dna_samples }
    ids = [s["id"] for s in dna_samples]
    if len(ids) < 3:
        console.print("[yellow]Menos de 3 amostras — paternidade requer trio.[/yellow]")
        return
    Path("trio").mkdir(exist_ok=True)
    merged = Path("trio")/"trio_merged.vcf.gz"
    need_merge = (not merged.exists()) or any(not _is_newer(merged, v) for v in vcfs.values())
    if need_merge:
        run(["bcftools","merge","-m","all","-Oz","-o",str(merged), *(str(vcfs[i]) for i in ids)])
        run(["tabix","-p","vcf",str(merged)])
    # Mapear nomes
    avail = _load_samples_from_vcf(merged)
    child_id, p1_id, p2_id = _infer_trio_ids(dna_samples)
    mapped = {}
    for wanted in (child_id, p1_id, p2_id):
        got = _map_name(wanted, avail)
        if not got:
            console.print(f"[red]Não foi possível mapear {wanted} no trio merged.[/red]")
            return
        mapped[wanted] = got

    # Todas as direções AP→Child (o terceiro é a mãe)
    sids = [child_id, p1_id, p2_id]
    triples = []
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            k = 3 - i - j
            ap, child, mom = sids[i], sids[j], sids[k]
            triples.append((ap, child, mom))

    # Checa rapidamente se há variantes PASS; respeita paternity_force_pass
    require_pass_effective = require_pass
    if require_pass and not bool(g.get("paternity_force_pass", False)):
        check_cmd = ["conda","run","-n","genomics","bash","-lc",
                     f"bcftools view -H -f PASS {shlex.quote(str(merged))} | head -1"]
        chk = sp.run(check_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, text=True, check=False)
        if chk.returncode != 0 or (not chk.stdout.strip()):
            require_pass_effective = False
            console.print("[dim]Nenhuma variante com FILTER=PASS no trio; usando todas as variantes.[/dim]")

    Path("paternity").mkdir(exist_ok=True)
    summary_lines = []
    results = []  # coleta resultados para sumarização
    for ap, child, mom in triples:
        ap_m = mapped[ap]; ch_m = mapped[child]; mo_m = mapped[mom]
        console.print(f"[cyan]👨‍👩‍👧 Paternidade: AP={ap} → Child={child} (Mãe={mom})[/cyan]")
        # Monkey-patch: injeta af_map via global temporário
        globals()["_paternity_af_map"] = af_map
        res = _paternity_for_pair(merged, ch_m, mo_m, ap_m, g, prior, eps, require_pass_effective)
        # salvar TSV
        tsv = Path("paternity")/f"{ap}_to_{child}.paternity.tsv"
        with open(tsv, "w") as fh:
            fh.write("CHROM\tPOS\tREF\tALT\tGT_child\tGT_mom\tGT_ap\tDP_child\tDP_mom\tDP_ap\tADc_ref\tADc_alt\tADm_ref\tADm_alt\tADf_ref\tADf_alt\tAF\tP_H1\tP_H2\tPI\tlog10_PI\n")
            for r in res["rows"]:
                fh.write("\t".join(map(str,r))+"\n")
        # imprime resultado imediato deste trio
        console.print(f"   → n={res['n_used']} | log10(LR)={res['lr_log10']:.3f} | posterior(prior={prior:.2f})={res['posterior']*100:.6f}% | excluídas={res['n_excluded']}\n")
        summary_lines.append(f"- {ap} → {child} | n={res['n_used']}, excl={res['n_excluded']}, log10(LR)={res['lr_log10']:.3f}, Posterior(prior={prior:.2f})={res['posterior']*100:.6f}%")
        results.append({
            "ap": ap, "child": child, "mom": mom,
            "n_used": res['n_used'], "n_excluded": res['n_excluded'],
            "lr_log10": res['lr_log10'], "posterior": res['posterior']
        })

    md = Path("paternity")/"paternity_summary.md"
    with open(md,"w") as fh:
        fh.write("# Paternidade (SNP, trio, HWE)\n")
        fh.write(f"- prior={prior}, epsilon={eps}, require_pass={require_pass_effective}\n")
        for ln in summary_lines:
            fh.write(ln+"\n")
    console.print(f"[green]✅ Paternidade concluída. Resumo: paternity/paternity_summary.md[/green]")

    # Impressão amigável no console
    console.print("\n[bold cyan]Resumo paternidade (SNP/HWE)[/bold cyan]")
    for ln in summary_lines:
        console.print(ln)

    # # Destaques por filho: melhor AP por LR
    # by_child: dict[str, list[dict]] = {}
    # for r in results:
    #     by_child.setdefault(r["child"], []).append(r)
    # for child, lst in by_child.items():
    #     if not lst:
    #         continue
    #     best = max(lst, key=lambda x: x["lr_log10"])  # maior LR favorece paternidade
    #     best_pct = best["posterior"]*100.0
    #     console.print(f"[bold]Filho {child}[/bold]: melhor AP = {best['ap']} (mãe={best['mom']}) | log10(LR)={best['lr_log10']:.3f} | posterior={best_pct:.6f}% | n={best['n_used']}")

# =================== Ancestralidade (ADMIXTURE supervisionado) ===================

def _ensure_ancestry_refs(anc_cfg):
    """Baixa e prepara o painel HGDP+1KG em PLINK (refs/ancestry)."""
    ref_dir = Path("refs/ancestry"); ref_dir.mkdir(parents=True, exist_ok=True)
    tar_url = anc_cfg["reference"]["plink_tar_url"]
    tsv_url = anc_cfg["reference"]["sample_info_url"]
    tar_path = ref_dir/"HGDP+1KG_SNPData.tar.gz"
    tsv_path = ref_dir/"hgdp_1kg_sample_info.tsv"
    
    # Download do tar se não existir
    if not tar_path.exists():
        console.print(f"[yellow]Baixando dados de referência HGDP+1KG...[/yellow]")
        run(["conda", "run", "-n", "genomics", "bash", "-lc", f"curl -L --retry 5 -o {tar_path} '{tar_url}'"])
    
    # Download do TSV se não existir
    if not tsv_path.exists():
        console.print(f"[yellow]Baixando informações das amostras...[/yellow]")
        run(["conda", "run", "-n", "genomics", "bash", "-lc", f"curl -L --retry 5 -o {tsv_path} '{tsv_url}'"])
    
    out_dir = ref_dir/"hgdp1kg"
    
    # Verificar se já foi extraído corretamente procurando qualquer arquivo .fam
    fam_files = list(out_dir.glob("*.fam")) if out_dir.exists() else []
    
    if not fam_files:
        console.print(f"[yellow]Extraindo dados de referência...[/yellow]")
        # Limpar diretório se existir e estiver vazio ou corrompido
        if out_dir.exists():
            import shutil
            shutil.rmtree(out_dir)
        out_dir.mkdir(exist_ok=True)
        
        # Extrair usando conda run para garantir ambiente correto
        run(["conda", "run", "-n", "genomics", "bash", "-lc", 
             f"tar -xzf {tar_path} -C {out_dir}"])
        
        # Verificar se a extração foi bem-sucedida
        fam_files = list(out_dir.glob("*.fam"))
        if not fam_files:
            raise RuntimeError(f"Falha na extração do tar. Nenhum arquivo .fam encontrado em {out_dir}.")
        console.print(f"[green]Extração concluída. Encontrados {len(fam_files)} arquivos .fam[/green]")
    
    return out_dir, tsv_path

def _make_ref_subset_and_popfile(ref_dir: Path, sample_info_tsv: Path, anc_cfg):
    """
    Cria subpainel PLINK com as populações desejadas e gera:
      - ref_sub.{bed,bim,fam}
      - ref.pop (categoria por indivíduo de referência na mesma ordem)
    """
    plink = (anc_cfg.get("tools") or {}).get("plink","plink")
    categories = list(anc_cfg["categories"].keys())
    cat_to_pops = anc_cfg["categories"]
    pop_to_cat = {p.lower():cat for cat, pops in cat_to_pops.items() for p in pops}

    # Verificar se o subpainel já foi criado (idempotência)
    ref_sub_bed = ref_dir/"ref_sub.bed"
    ref_pop_file = ref_dir/"ref.pop"
    
    if ref_sub_bed.exists() and ref_pop_file.exists():
        console.print("[bold]Subpainel de referência → SKIP (já existe)[/bold]")
        return ref_dir/"ref_sub", ref_pop_file, categories

    # === NOVO: detectar automaticamente o prefixo dos arquivos PLINK extraídos ===
    fam_list = sorted(ref_dir.glob("*.fam"))
    if not fam_list:
        raise RuntimeError(f"Nenhum .fam encontrado em {ref_dir}. Verifique a extração do tar.")
    # preferir nomes que contenham 'hgdp' ou '1kg', senão pega o primeiro
    preferred = [p for p in fam_list if ("hgdp" in p.stem.lower() or "1kg" in p.stem.lower())]
    fam_path = preferred[0] if preferred else fam_list[0]
    ref_prefix = fam_path.with_suffix("")  # caminho sem extensão para --bfile

    # Carregar IDs presentes no .fam (Family ID e Individual ID)
    fam_data = []
    with open(fam_path) as f:
        for ln in f:
            parts = ln.split()
            if len(parts) >= 2:
                fam_id, ind_id = parts[0], parts[1]
                fam_data.append((fam_id, ind_id))
    
    # Criar um mapeamento de Individual ID para Family ID
    ind_to_fam = {ind_id: fam_id for fam_id, ind_id in fam_data}
    fam_ids = [ind_id for fam_id, ind_id in fam_data]
    possible_id_cols = ["s","sample","sample_id","Sample","ID","iid","IID"]
    possible_pop_cols= ["pop","population","Population","research_population","SuperPop"]

    rows=[]
    with open(sample_info_tsv,"r") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        id_idx = next((i for i,h in enumerate(header) if h in possible_id_cols), None)
        pop_idx= next((i for i,h in enumerate(header) if h in possible_pop_cols), None)
        if id_idx is None or pop_idx is None:
            raise RuntimeError("Não encontrei colunas de ID/pop em sample_info TSV.")
        for r in rd:
            rows.append((r[id_idx], r[pop_idx].lower()))

    keep=[]
    fam_set=set(fam_ids)
    for iid,pop in rows:
        if iid in fam_set and pop in pop_to_cat:
            keep.append((iid, pop_to_cat[pop]))
    if not keep:
        raise RuntimeError("Nenhuma amostra do painel corresponde às populações pedidas.")

    keep_path = ref_dir/"keep.ids"
    pop_path  = ref_dir/"ref.pop"
    with open(keep_path,"w") as fh, open(pop_path,"w") as ph:
        for iid,cat in keep:
            # Usar Family ID e Individual ID do arquivo .fam
            fam_id = ind_to_fam.get(iid, iid)  # fallback para o próprio ID se não encontrar
            fh.write(f"{fam_id} {iid}\n")
            ph.write(f"{cat}\n")

    # usar o prefixo detectado e conda run
    console.print(f"[yellow]Criando subpainel de referência com {len(keep)} amostras...[/yellow]")
    run(["conda", "run", "-n", "genomics", plink, "--bfile", str(ref_prefix), "--keep", str(keep_path),
                 "--make-bed", "--out", str(ref_dir/"ref_sub")])

    return ref_dir/"ref_sub", pop_path, categories

def _convert_vcf_to_bed(vcf: Path, out_prefix: Path, anc_cfg):
    """Converte VCF para formato PLINK BED com idempotência."""
    bed_file = out_prefix.with_suffix(".bed")
    
    if bed_file.exists():
        console.print(f"[bold]Conversão VCF→BED {out_prefix.name} → SKIP (já existe)[/bold]")
        return
        
    plink2 = anc_cfg["tools"].get("plink2","plink2")
    console.print(f"[yellow]Convertendo {vcf.name} para formato PLINK...[/yellow]")
    
    # Primeira conversão: apenas filtros básicos
    temp_prefix = out_prefix.with_name(out_prefix.stem + ".temp")
    run(["conda", "run", "-n", "genomics", plink2, "--vcf", str(vcf), "--double-id", "--snps-only", "just-acgt",
                 "--max-alleles", "2", "--chr", "1-22", "--make-bed", "--out", str(temp_prefix)])
    
    # Segunda etapa: normalizar IDs usando PLINK1 com contador sequencial
    plink = anc_cfg["tools"].get("plink","plink")
    console.print(f"[yellow]Normalizando IDs das variantes...[/yellow]")
    run(["conda", "run", "-n", "genomics", plink, "--bfile", str(temp_prefix), 
                 "--set-missing-var-ids", "@:#_\\$1_\\$2", "--make-bed", "--out", str(out_prefix)])
    
    # Limpar apenas arquivos temporários (com .temp no nome)
    for suf in [".bed", ".bim", ".fam", ".log", ".nosex"]:
        temp_file = temp_prefix.with_suffix(suf)
        if temp_file.exists() and ".temp" in temp_file.name:
            temp_file.unlink()

def _merge_and_prune(ref_prefix: Path, samp_prefix: Path, anc_cfg):
    """Merge dados de referência com amostra e aplica QC com idempotência."""
    final_bed = samp_prefix.with_name(samp_prefix.stem + ".pruned.bed")
    
    if final_bed.exists():
        console.print(f"[bold]Merge e QC {samp_prefix.name} → SKIP (já existe)[/bold]")
        return samp_prefix.with_name(samp_prefix.stem + ".pruned")
    
    plink = anc_cfg["tools"].get("plink","plink")
    qc = anc_cfg["qc"]; maf = str(qc.get("maf",0.01)); geno = str(qc.get("geno_missing",0.05))
    mind = str(qc.get("mind",0.99))
    w,s,r2 = map(str, qc.get("indep_pairwise",[200,50,0.2]))
    
    console.print(f"[yellow]Fazendo merge e QC para {samp_prefix.name}...[/yellow]")
    
    # Merge (com tratamento de inconsistências de strand)
    merge_cmd = ["conda", "run", "-n", "genomics", plink, "--bfile", str(ref_prefix), "--bmerge", f"{samp_prefix}.bed", f"{samp_prefix}.bim", f"{samp_prefix}.fam",
                 "--make-bed", "--out", f"{samp_prefix}.merged"]
    
    try:
        run(merge_cmd)
    except sp.CalledProcessError:
        # Se o merge falhar, aplicar estratégia de correção
        console.print(f"[yellow]Merge falhou. Aplicando correções de compatibilidade...[/yellow]")
        flip_file = f"{samp_prefix}.merged-merge.missnp"
        
        if Path(flip_file).exists():
            # Tentar flip nos dados originais
            try:
                run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{samp_prefix}", "--flip", flip_file,
                             "--make-bed", "--out", f"{samp_prefix}.flipped"])
                # Merge com dados flipped
                run(["conda", "run", "-n", "genomics", plink, "--bfile", str(ref_prefix), "--bmerge", f"{samp_prefix}.flipped.bed", f"{samp_prefix}.flipped.bim", f"{samp_prefix}.flipped.fam",
                             "--make-bed", "--out", f"{samp_prefix}.merged"])
            except sp.CalledProcessError:
                # Fallback: excluir SNPs problemáticos
                console.print(f"[yellow]Flip falhou. Excluindo SNPs problemáticos de ambos os painéis...[/yellow]")
                run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{samp_prefix}", "--exclude", flip_file,
                             "--make-bed", "--out", f"{samp_prefix}.clean"])
                run(["conda", "run", "-n", "genomics", plink, "--bfile", str(ref_prefix), "--exclude", flip_file,
                             "--make-bed", "--out", f"{ref_prefix}.clean"])
                run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{ref_prefix}.clean", "--bmerge", f"{samp_prefix}.clean.bed", f"{samp_prefix}.clean.bim", f"{samp_prefix}.clean.fam",
                             "--make-bed", "--out", f"{samp_prefix}.merged"])
        else:
            raise RuntimeError(f"Merge falhou e arquivo de flip {flip_file} não foi gerado")
    
    # Filtros de QC
    run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{samp_prefix}.merged", "--maf", maf, "--geno", geno,
                 "--make-bed", "--out", f"{samp_prefix}.flt"])
    
    # Pruning
    run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{samp_prefix}.flt", "--indep-pairwise", w, s, r2, "--out", f"{samp_prefix}.flt"])
    run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{samp_prefix}.flt", "--extract", f"{samp_prefix}.flt.prune.in",
                 "--make-bed", "--out", f"{samp_prefix}.pruned"])
    
    # Filtro adicional: remover amostras com muitos genótipos faltantes
    console.print(f"[yellow]Aplicando filtro de qualidade das amostras (mind={mind})...[/yellow]")
    run(["conda", "run", "-n", "genomics", plink, "--bfile", f"{samp_prefix}.pruned", "--mind", mind,
                 "--make-bed", "--out", f"{samp_prefix}.pruned"])
    
    # Copiar arquivo .pop para corresponder aos dados pruned
    base_pop = f"{samp_prefix}.pop"  # arquivo .pop original (antes do merge)
    pruned_pop = f"{samp_prefix}.pruned.pop"
    if Path(base_pop).exists():
        console.print(f"[yellow]Copiando arquivo .pop para dados pruned...[/yellow]")
        import shutil
        shutil.copy2(base_pop, pruned_pop)
    
    return samp_prefix.with_name(samp_prefix.stem + ".pruned")

def _map_columns_to_categories(q_lines, popfile_path: Path, categories):
    """Mapeia colunas do Q para categorias usando médias nas amostras de referência."""
    pop_labels=[ln.strip() for ln in open(popfile_path) if ln.strip()]
    n_ref=len(pop_labels)
    Qref=[[float(x) for x in q_lines[i].split()] for i in range(n_ref)]
    K=len(Qref[0])
    means={cat:[0.0]*K for cat in categories}
    counts={cat:0 for cat in categories}
    for row,cat in zip(Qref, pop_labels):
        if cat in means:
            counts[cat]+=1
            for j,v in enumerate(row): means[cat][j]+=v
    for cat in categories:
        if counts[cat]>0: means[cat]=[v/counts[cat] for v in means[cat]]
    # atribuição gulosa coluna↔categoria
    used=set(); mapping={}
    prefs=[]
    for cat in categories:
        order=sorted(range(K), key=lambda j: means[cat][j], reverse=True)
        prefs.extend((cat,j,means[cat][j]) for j in order)
    prefs.sort(key=lambda t:t[2], reverse=True)
    for cat,j,_ in prefs:
        if cat in mapping or j in used: continue
        mapping[cat]=j; used.add(j)
        if len(mapping)==K: break
    for cat in categories:
        if cat not in mapping:
            for j in range(K):
                if j not in used: mapping[cat]=j; used.add(j); break
    return [mapping[cat] for cat in categories]

def ancestry_admixture_step(cfg):
    if not cfg["steps"].get("ancestry", False): return
    anc = cfg.get("ancestry", {})
    threads = int(anc.get("threads", 8))
    admixture = anc["tools"].get("admixture","admixture")
    categories = list(anc.get("categories", {}).keys())
    K = int(anc.get("k", max(2, len(categories))))  # K = populações de referência apenas
    
    Path("ancestry").mkdir(exist_ok=True)
    
    # Verificar se os resultados finais já existem (idempotência completa)
    out_k = Path("ancestry")/f"ancestry_summary_K{K}.tsv"
    collapse = anc.get("collapse", None)
    out_col = Path("ancestry")/"ancestry_summary_collapsed.tsv" if collapse else None
    
    all_outputs_exist = out_k.exists() and (out_col is None or out_col.exists())
    
    if all_outputs_exist:
        console.print("[bold]Ancestralidade → SKIP (resultados já existem)[/bold]")
        outs = [out_k]
        if out_col: outs.append(out_col)
        print_meta(f"Ancestralidade (ADMIXTURE K={K})", outs)
        return
    
    ref_dir, sample_info_tsv = _ensure_ancestry_refs(anc)
    ref_sub, ref_pop, categories = _make_ref_subset_and_popfile(ref_dir, sample_info_tsv, anc)
    
    samples = [s["id"] for s in cfg["dna_samples"]]
    header_k = ["sample"] + categories
    lines_k = ["\t".join(header_k)]
    lines_col = ["\t".join(["sample"] + list(collapse.keys()))] if isinstance(collapse, dict) else None
    
    for sid in samples:
        vcf = Path("vcf")/f"{sid}.vcf.gz"
        if not vcf.exists():
            console.print(f"[yellow]VCF da amostra {sid} não encontrado (pulando).[/yellow]")
            continue
            
        # Verificar se os resultados já existem para esta amostra
        outpfx = Path("ancestry")/sid
        qfile_expected = outpfx.with_name(outpfx.stem + f".pruned.{K}.Q")
        
        if qfile_expected.exists():
            console.print(f"[bold]ADMIXTURE {sid} → SKIP (resultado já existe)[/bold]")
        else:
            _convert_vcf_to_bed(vcf, outpfx, anc)
            
            # Criar arquivo .pop para ADMIXTURE ANTES do merge/prune
            temp_popfile = outpfx.with_suffix(".pop")
            with open(ref_pop,"r") as rp, open(temp_popfile,"w") as ph:
                for ln in rp: ph.write(ln)
                # NÃO adicionar "0" - ADMIXTURE supervisionado não precisa
            
            merged = _merge_and_prune(ref_sub, outpfx, anc)
            popfile = merged.with_suffix(".pop")
            
            console.print(f"[yellow]Executando ADMIXTURE supervisionado para {sid}...[/yellow]")
            # ADMIXTURE deve executar no diretório ancestry/ para gerar arquivos no local correto
            current_dir = os.getcwd()
            ancestry_dir = os.path.join(current_dir, "ancestry")
            run(["conda", "run", "-n", "genomics", admixture, "--supervised", f"-j{threads}", 
                 f"{merged.name}.bed", str(K)], cwd=ancestry_dir)
        
        # Processar resultados
        qfile = qfile_expected
        if qfile.exists():
            q_lines=[l.strip() for l in open(qfile) if l.strip()]
            col_order = _map_columns_to_categories(q_lines, qfile.parent/(qfile.stem.replace(f".{K}", "") + ".pop"), categories)
            sample_vec=[float(x) for x in q_lines[-1].split()]
            vec_ord=[sample_vec[j] for j in col_order]
            lines_k.append("\t".join([sid] + [f"{100.0*v:.2f}" for v in vec_ord]))
            
            if lines_col is not None:
                cat_val={cat:val for cat,val in zip(categories, vec_ord)}
                collapse_vals=[]
                for grp, cats in collapse.items():
                    s=sum(cat_val.get(c,0.0) for c in cats)
                    collapse_vals.append(f"{100.0*s:.2f}")
                lines_col.append("\t".join([sid]+collapse_vals))
            
            console.print(f"[green]→ {sid}: " + ", ".join([f"{cat} {100.0*v:.1f}%" for cat,v in zip(categories, vec_ord)]) + "[/green]")
        else:
            console.print(f"[red]Erro: arquivo de resultados {qfile} não foi gerado para {sid}[/red]")
    
    # Salvar resultados finais
    with open(out_k,"w") as fh: fh.write("\n".join(lines_k)+"\n")
    outs=[out_k]
    
    if lines_col is not None:
        with open(out_col,"w") as fh: fh.write("\n".join(lines_col)+"\n")
        outs.append(out_col)
    
    print_meta(f"Ancestralidade (ADMIXTURE K={K})", outs)

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
    run(["conda", "run", "-n", "genomics", "bash","-lc",cmd])
    print_meta("Lista de genes", [out])
# ============ Presença de genes por amostra (cobertura) ============

def _build_genes_bed_from_gtf(gtf: Path, out_bed: Path):
    """Extrai regiões 'gene' do GTF como BED (chr start end name). Recria se o GTF for mais novo."""
    if out_bed.exists() and _is_newer(out_bed, gtf):
        return out_bed
    # Comando AWK robusto que valida dados antes de escrever
    cmd = r"""awk '$3=="gene" && NF>=9 && $4~/^[0-9]+$/ && $5~/^[0-9]+$/ { 
        chr=$1; start=$4-1; end=$5; 
        
        # Valida coordenadas
        if(start<0 || end<=start) next;
        
        # Extrai atributos com validação
        match($0,/gene_id "([^"]+)"/,a); gid=a[1];
        match($0,/gene_name "([^"]+)"/,b); gname=b[1];
        match($0,/gene_type "([^"]+)"/,c); gtype=c[1];
        match($0,/gene_biotype "([^"]+)"/,d); if(d[1]=="") d[1]=c[1];
        match($0,/gene_description "([^"]+)"/,e); gdesc=e[1];
        
        # Só escreve se gene_id for válido
        if(gid!="" && chr!="") {
            name=gname!=""?gname:gid;
            print chr"\t"start"\t"end"\t"gid"\t"name"\t"gtype"\t"gdesc;
        }
    }' """ + shlex.quote(str(gtf)) + " > " + shlex.quote(str(out_bed))
    
    console.print(f"[cyan]🧬 Extraindo genes do GTF...[/cyan]")
    console.print(f"[dim]💻 Comando GTF parsing:[/dim]")
    console.print(f"[dim]> {cmd}[/dim]")
    run(["conda", "run", "-n", "genomics", "bash","-lc",cmd])
    
    # Verifica se arquivo BED foi criado corretamente
    if not out_bed.exists():
        raise RuntimeError(f"Falha ao criar arquivo BED: {out_bed}")
    
    # Conta genes extraídos e valida formato
    gene_count = 0
    with open(out_bed) as f:
        for line_num, line in enumerate(f, 1):
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                try:
                    # Valida que coordenadas são números
                    int(parts[1])  # start
                    int(parts[2])  # end
                    gene_count += 1
                except ValueError:
                    console.print(f"[red]❌ Linha {line_num} inválida no BED: {line.strip()[:80]}...[/red]")
                    raise RuntimeError(f"Arquivo BED contém dados inválidos na linha {line_num}")
    
    console.print(f"[green]✅ Arquivo BED criado: {gene_count:,} genes extraídos[/green]")
    return out_bed

def _parse_thresholds_bed(thr_bed_gz: Path, region_len: dict) -> dict:
    """
    Lê *.thresholds.bed.gz do mosdepth --by, retorna breadth_1x por gene_id.
    Formato: chrom start end region_name thresh1 thresh2 ... counts_above
    Usamos a 1a coluna de contagem como >=1x.
    """
    import gzip
    breadth = {}
    line_num = 0
    valid_lines = 0
    
    with gzip.open(thr_bed_gz, "rt") as fh:
        for line in fh:
            line_num += 1
            line = line.strip()
            
            if not line:
                continue
                
            # Pula headers e comentários
            if line.startswith("#") or line.startswith("chrom") or "start" in line[:20]:
                continue
                
            parts = line.split("\t")
            if len(parts) < 5:
                continue
                
            try:
                # Valida que as coordenadas são números válidos
                start_coord = int(parts[1])
                end_coord = int(parts[2])
                
                # Valida coordenadas lógicas
                if start_coord < 0 or end_coord <= start_coord:
                    continue
                    
                gid = parts[3]
                if not gid or gid == "":
                    continue
                    
                length = max(1, end_coord - start_coord)
                
                # Processa contagens de threshold
                counts = []
                for count_str in parts[4:]:
                    try:
                        counts.append(int(count_str))
                    except ValueError:
                        counts.append(0)
                
                cov1 = counts[0] if counts else 0
                breadth[gid] = cov1 / length
                region_len[gid] = length
                valid_lines += 1
                
            except (ValueError, IndexError) as e:
                # Log apenas as primeiras linhas problemáticas para não spam
                if line_num <= 10:
                    console.print(f"[yellow]⚠️  Linha {line_num} inválida em {thr_bed_gz.name}: {str(e)[:50]}...[/yellow]")
                    console.print(f"[dim]Conteúdo: {line[:80]}...[/dim]")
                continue
                
    console.print(f"[cyan]📊 {thr_bed_gz.name}: {valid_lines:,} genes válidos de {line_num:,} linhas[/cyan]")
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
            console.print(f"[cyan]🧬 Calculando cobertura por gene para {sample}...[/cyan]")
            
            cmd = ["mosdepth","-t",str(threads),"--by",str(genes_bed),
                   "--thresholds","1,5,10",str(prefix), str(bam)]
            if bam.suffix == ".cram":
                # por segurança, garante FASTA quando entrada é CRAM
                cmd = ["mosdepth","-t",str(threads),"--by",str(genes_bed),
                       "--thresholds","1,5,10","--fasta","refs/reference.fa",
                       str(prefix), str(bam)]
            
            # Mostra comando completo que será executado
            console.print(f"[dim]💻 Comando mosdepth:[/dim]")
            console.print(f"[dim]> {' '.join(cmd)}[/dim]")
            
            # Executa e aguarda conclusão
            console.print(f"[cyan]⏳ Executando mosdepth (pode demorar 10-30min com {threads} threads)...[/cyan]")
            run(cmd)
            console.print(f"[green]✅ mosdepth concluído para {sample}[/green]")
        else:
            console.print(f"[{sample}] mosdepth → [bold]SKIP[/bold] (cache)")
        
        # Verifica se arquivos foram criados corretamente antes de processar
        if not regions_bed_gz.exists():
            raise RuntimeError(f"mosdepth não criou arquivo regions: {regions_bed_gz}")
        if not thr_bed_gz.exists():
            raise RuntimeError(f"mosdepth não criou arquivo thresholds: {thr_bed_gz}")
        
        console.print(f"[cyan]📊 Processando resultados mosdepth para {sample}...[/cyan]")

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

from .rnaseq import rnaseq_pipeline
