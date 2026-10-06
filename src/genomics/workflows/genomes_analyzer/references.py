"""Reference FASTA/GTF download and aligner index preparation."""

import subprocess as sp
import re
import shlex
import time
from pathlib import Path

from . import legacy


def _bwa_index_ready(prefix: Path) -> bool:
    for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        path = Path(str(prefix) + ext)
        if not path.exists():
            return False
        if path.is_symlink() and not path.resolve().exists():
            return False
    return True


def _mem2_index_present(prefix: Path) -> bool:
    return Path(str(prefix) + ".bwt.2bit.64").exists() or any(Path("refs").glob("reference.fa.*.0123"))


def _build_bwa_index_optimized(ref_prefix: Path, label: str = "BWA"):
    import queue
    import threading

    import psutil

    params = legacy.cfg_global.get("params", {})
    max_mem_gb = int(params.get("bwa_index_max_mem_gb", 100))
    block_size_cfg = params.get("bwa_index_block_size", "auto")
    algorithm = params.get("bwa_index_algorithm", "bwtsw")
    progress_sec = int(params.get("bwa_index_progress_sec", 60))

    total_ram_gb = psutil.virtual_memory().total / (1024**3)
    if block_size_cfg == "auto":
        if total_ram_gb >= 200:
            block_size = min(2_000_000_000, max_mem_gb * 20_000_000)
        elif total_ram_gb >= 64:
            block_size = min(1_000_000_000, int(total_ram_gb * 15_000_000))
        else:
            block_size = 10_000_000
    else:
        block_size = int(block_size_cfg)

    estimated_ram_gb = block_size / 50_000_000
    if estimated_ram_gb > max_mem_gb:
        block_size = max_mem_gb * 50_000_000
        estimated_ram_gb = max_mem_gb

    from rich.panel import Panel

    legacy.console.print(
        Panel.fit(
            f"[bold]Criando Índice BWA Otimizado ({label})[/bold]\n"
            f"• RAM total: {total_ram_gb:.0f}GB\n"
            f"• RAM para indexação: ~{estimated_ram_gb:.0f}GB\n"
            f"• Block size: {block_size:,}\n"
            f"• Algoritmo: {algorithm}\n"
            f"• Tempo estimado: {45 if total_ram_gb >= 200 else 90}-{90 if total_ram_gb >= 200 else 180}min",
            border_style="cyan",
        )
    )

    cmd = ["bwa", "index", "-a", algorithm, "-b", str(block_size), str(ref_prefix)]
    legacy.console.print("[bold]💻 Comando executado:[/bold]")
    legacy.console.print(f"[yellow]> {' '.join(cmd)}[/yellow]")
    legacy.console.print("[cyan]🚀 Iniciando criação do índice BWA...[/cyan]")

    start_time = time.time()
    last_progress = start_time
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, bufsize=1)
    stderr_queue = queue.Queue()
    stderr_lines = []

    def stderr_reader():
        try:
            while True:
                line = proc.stderr.readline()
                if not line:
                    break
                stderr_queue.put(line)
                stderr_lines.append(line)
        except Exception:
            pass

    stderr_thread = threading.Thread(target=stderr_reader, daemon=True)
    stderr_thread.start()

    while True:
        rc = proc.poll()
        now = time.time()

        while not stderr_queue.empty():
            try:
                stderr_line = stderr_queue.get_nowait().strip()
                if stderr_line and ("WARNING" in stderr_line or "ERROR" in stderr_line):
                    legacy.console.print(f"[yellow]⚠️  BWA: {stderr_line[:100]}{'...' if len(stderr_line) > 100 else ''}[/yellow]")
            except queue.Empty:
                break

        if now - last_progress >= progress_sec:
            elapsed = int(now - start_time)
            files_info = []
            total_size = 0
            for ext, display_name in ((".amb", "amb"), (".ann", "ann"), (".bwt", "bwt"), (".pac", "pac"), (".sa", "sa")):
                path = Path(str(ref_prefix) + ext)
                if path.exists() and path.stat().st_size > 0:
                    total_size += path.stat().st_size
                    files_info.append(f"{display_name}({legacy.sizeof_fmt(path.stat().st_size)})")
            status = f"total {legacy.sizeof_fmt(total_size)} • {', '.join(files_info[:3])}" if files_info else "iniciando..."
            legacy.console.print(f"[cyan]BWA index … {elapsed//60}m{elapsed%60:02d}s • {status}[/cyan]", highlight=False)
            last_progress = now

        if rc is not None:
            break
        time.sleep(5)

    stderr_thread.join(timeout=2.0)
    if proc.returncode != 0:
        legacy.console.print(f"[red]❌ BWA index falhou com código {proc.returncode}[/red]")
        for line in stderr_lines[-5:]:
            if line.strip():
                legacy.console.print(f"[dim]{line.strip()}[/dim]")
        raise sp.CalledProcessError(proc.returncode, cmd)

    warnings = [line for line in stderr_lines if "WARNING" in line]
    if warnings:
        legacy.console.print(f"[yellow]⚠️  BWA index concluído com {len(warnings)} warning(s)[/yellow]")

    elapsed = int(time.time() - start_time)
    total_size = sum(Path(str(ref_prefix) + ext).stat().st_size for ext in (".amb", ".ann", ".bwt", ".pac", ".sa") if Path(str(ref_prefix) + ext).exists())
    legacy.console.print(f"[bold green]✅ Índice BWA criado em {elapsed//60}m{elapsed%60:02d}s • tamanho total: {legacy.sizeof_fmt(total_size)}[/bold green]")


def _build_bwa_mem2_index_optimized(ref_prefix: Path, label: str = "BWA-MEM2"):
    import queue
    import threading

    import psutil
    from rich.panel import Panel

    params = legacy.cfg_global.get("params", {})
    progress_sec = int(params.get("bwa_index_progress_sec", 60))
    total_ram_gb = psutil.virtual_memory().total / (1024**3)
    estimated_ram_gb = min(total_ram_gb * 0.8, 200)

    legacy.console.print(
        Panel.fit(
            f"[bold]Criando Índice BWA-MEM2 ({label})[/bold]\n"
            f"• RAM total: {total_ram_gb:.0f}GB\n"
            f"• RAM que será usada: ~{estimated_ram_gb:.0f}GB (automático)\n"
            f"• Algoritmo: BWA-MEM2 (otimizado internamente)\n"
            f"• Tempo estimado: {20 if total_ram_gb >= 200 else 45}-{45 if total_ram_gb >= 200 else 90}min\n"
            f"• NOTA: BWA-MEM2 não aceita parâmetros de block size",
            border_style="green",
        )
    )

    cmd = ["bwa-mem2", "index", str(ref_prefix)]
    legacy.console.print("[bold]💻 Comando executado:[/bold]")
    legacy.console.print(f"[yellow]> {' '.join(cmd)}[/yellow]")
    legacy.console.print("[green]🚀 Iniciando criação do índice BWA-MEM2...[/green]")

    start_time = time.time()
    last_progress = start_time
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, bufsize=1)
    stderr_queue = queue.Queue()
    stderr_lines = []

    def stderr_reader():
        try:
            while True:
                line = proc.stderr.readline()
                if not line:
                    break
                stderr_queue.put(line)
                stderr_lines.append(line)
        except Exception:
            pass

    stderr_thread = threading.Thread(target=stderr_reader, daemon=True)
    stderr_thread.start()

    while True:
        rc = proc.poll()
        now = time.time()

        while not stderr_queue.empty():
            try:
                stderr_line = stderr_queue.get_nowait().strip()
                if stderr_line and ("WARNING" in stderr_line or "ERROR" in stderr_line):
                    legacy.console.print(f"[yellow]⚠️  BWA-MEM2: {stderr_line[:100]}{'...' if len(stderr_line) > 100 else ''}[/yellow]")
            except queue.Empty:
                break

        if now - last_progress >= progress_sec:
            elapsed = int(now - start_time)
            files_info = []
            total_size = 0
            for ext, display_name in ((".0123", "idx"), (".amb", "amb"), (".ann", "ann"), (".bwt.2bit.64", "bwt"), (".pac", "pac")):
                path = Path(str(ref_prefix) + ext)
                if path.exists() and path.stat().st_size > 0:
                    total_size += path.stat().st_size
                    files_info.append(f"{display_name}({legacy.sizeof_fmt(path.stat().st_size)})")
            status = f"total {legacy.sizeof_fmt(total_size)} • {', '.join(files_info[:3])}" if files_info else "iniciando..."
            legacy.console.print(f"[green]BWA-MEM2 index … {elapsed//60}m{elapsed%60:02d}s • {status}[/green]", highlight=False)
            last_progress = now

        if rc is not None:
            break
        time.sleep(5)

    stderr_thread.join(timeout=2.0)
    if proc.returncode != 0:
        legacy.console.print(f"[red]❌ BWA-MEM2 index falhou com código {proc.returncode}[/red]")
        for line in stderr_lines[-5:]:
            if line.strip():
                legacy.console.print(f"[dim]{line.strip()}[/dim]")
        raise sp.CalledProcessError(proc.returncode, cmd)

    warnings = [line for line in stderr_lines if "WARNING" in line]
    if warnings:
        legacy.console.print(f"[yellow]⚠️  BWA-MEM2 index concluído com {len(warnings)} warning(s)[/yellow]")

    elapsed = int(time.time() - start_time)
    total_size = sum(path.stat().st_size for path in Path("refs").glob("reference.fa.*") if path.is_file())
    legacy.console.print(f"[bold green]✅ Índice BWA-MEM2 criado em {elapsed//60}m{elapsed%60:02d}s • tamanho total: {legacy.sizeof_fmt(total_size)}[/bold green]")


def build_indexes(default_read_type, assembly_name, need_rna_index, threads, force=False):
    g = legacy.cfg_global["general"]

    if default_read_type != "short" and not need_rna_index:
        legacy.console.print("Sem necessidade de indexar (long reads).", style="dim")
        return

    ref_prefix = Path("refs/reference.fa")
    aligner_cfg = (g.get("aligner") or "bwa-mem2").lower()
    prebuilt_url = g.get("bwa_prebuilt_url", "")

    if prebuilt_url and g.get("limit_to_canonical", False):
        legacy.console.print(
            "[orange3]Aviso:[/orange3] você definiu [bold]limit_to_canonical:true[/bold] "
            "E também um [bold]bwa_prebuilt_url[/bold]. Isso causa mismatch entre FASTA e índice. "
            "Desative limit_to_canonical ao usar o índice prebuilt do GDC.",
            style="italic",
        )

    if default_read_type == "short":
        if prebuilt_url:
            install_prebuilt_bwa_index(prebuilt_url, force=force)
            if not _bwa_index_ready(ref_prefix):
                raise RuntimeError("Índice BWA não encontrado ou incompleto após instalar o prebuilt.")
            legacy.console.print("Índice BWA pré-instalado e linkado → [bold]pronto[/bold]")
        else:
            if aligner_cfg in ("bwa-mem2", "bwa_mem2", "mem2"):
                if not force and _mem2_index_present(ref_prefix):
                    legacy.console.print("Índice BWA-MEM2 → [bold]SKIP[/bold]", style="dim")
                else:
                    try:
                        _build_bwa_mem2_index_optimized(ref_prefix, "BWA-MEM2")
                        legacy.console.print("Índice BWA-MEM2 construído.", style="dim")
                    except sp.CalledProcessError:
                        legacy.console.print("[red]Falha ao indexar com bwa-mem2 (provável falta de RAM).[/red]")
                        legacy.console.print(
                            "Fallback: tentando preparar índice do [bold]BWA clássico[/bold]. "
                            "Sugestões: use [bold]bwa_prebuilt_url[/bold] no YAML, ou reduza a referência.",
                            style="italic",
                        )
                        if g.get("bwa_prebuilt_url"):
                            install_prebuilt_bwa_index(g["bwa_prebuilt_url"], force=True)
                            if not _bwa_index_ready(ref_prefix):
                                raise RuntimeError("Prebuilt BWA apontado não ficou completo após fallback.")
                            legacy.console.print("Índice BWA pré-instalado e linkado (fallback) → [bold]pronto[/bold]")
                        else:
                            try:
                                _build_bwa_index_optimized(ref_prefix, "fallback BWA clássico")
                                if not _bwa_index_ready(ref_prefix):
                                    raise RuntimeError("Índice BWA clássico parece incompleto após 'bwa index'.")
                                legacy.console.print("Índice BWA clássico construído (fallback).", style="dim")
                            except sp.CalledProcessError:
                                raise RuntimeError(
                                    "Não foi possível preparar nenhum índice (mem2 falhou e 'bwa index' também). "
                                    "Defina 'project.reference.bwa_index_url' no YAML para usar o índice pré-pronto do GDC."
                                )
            else:
                if force or not _bwa_index_ready(ref_prefix):
                    if not force and _bwa_index_ready(ref_prefix):
                        legacy.console.print("Índice BWA → [bold]SKIP[/bold]", style="dim")
                    else:
                        _build_bwa_index_optimized(ref_prefix, "BWA clássico")
                        if not _bwa_index_ready(ref_prefix):
                            raise RuntimeError("Índice BWA clássico parece incompleto após 'bwa index'.")
                        legacy.console.print("Índice BWA clássico construído.", style="dim")
                else:
                    legacy.console.print("Índice BWA → [bold]SKIP[/bold]", style="dim")

    if need_rna_index:
        ht2_base = Path(f"refs/{assembly_name}")
        if force or not Path(f"{ht2_base}.1.ht2").exists():
            legacy.run(["hisat2-build", "-p", str(threads), str(ref_prefix), str(ht2_base)])
        else:
            legacy.console.print("Índice HISAT2 → [bold]SKIP[/bold]", style="dim")


def install_prebuilt_bwa_index(bwa_tar_url: str, expect_md5: str = "", force: bool = False):
    exts = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    ref_prefix = Path("refs/reference.fa")
    index_dir = Path("refs/_bwa")
    index_dir.mkdir(parents=True, exist_ok=True)

    def _links_ok() -> bool:
        for ext in exts:
            path = Path(str(ref_prefix) + ext)
            if not path.exists():
                return False
            if path.is_symlink():
                try:
                    path.resolve(strict=True)
                except FileNotFoundError:
                    return False
            elif not path.is_file():
                return False
        return True

    def _any_prefix_in_dir() -> Path | None:
        amb = next(index_dir.rglob("*.fa.amb"), None)
        return amb.with_suffix("") if amb else None

    def _prefix_has_all(prefix: Path) -> bool:
        return all(Path(str(prefix) + ext).exists() for ext in exts)

    def _uuid_basename_from_url(url: str, default: str) -> str:
        match = re.search(r"([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$", url, re.I)
        return (match.group(1) + ".tar.gz") if match else default

    if not force and _links_ok():
        legacy.console.print("Índice BWA já presente e válido → [bold]SKIP[/bold]")
        legacy.print_meta("Índice BWA linkado", [Path(str(ref_prefix) + ext) for ext in exts])
        return

    if not bwa_tar_url:
        if not _links_ok():
            raise RuntimeError("Sem bwa_tar_url e índice BWA ausente/incompleto.")
        return

    tar_name = _uuid_basename_from_url(bwa_tar_url, "index.tar.gz")
    tar_path = index_dir / tar_name

    if force or not tar_path.exists():
        legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"wget -c -O {shlex.quote(str(tar_path))} {shlex.quote(bwa_tar_url)}"])
    else:
        legacy.console.print(f"{tar_name} já presente → [bold]SKIP download[/bold]")

    if expect_md5:
        md5 = sp.run(
            ["conda", "run", "-n", "genomics", "bash", "-lc", f"md5sum {shlex.quote(str(tar_path))} | cut -d' ' -f1"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        if md5 != expect_md5:
            raise RuntimeError(f"MD5 não confere para {tar_name}: {md5} (esp. {expect_md5})")

    best_prefix = None
    if not force:
        prefix = _any_prefix_in_dir()
        if prefix and _prefix_has_all(prefix):
            best_prefix = prefix

    extracted = False
    if best_prefix is None:
        listing = sp.run(
            ["conda", "run", "-n", "genomics", "bash", "-lc", f"tar -tzf {shlex.quote(str(tar_path))}"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.splitlines()
        wanted = [path for path in listing if re.search(r"\.fa\.(amb|ann|bwt|pac|sa)$", path)]
        if not wanted:
            from rich.panel import Panel

            head = "\n".join(listing[:10])
            legacy.console.print(
                Panel.fit(
                    "Tar não contém *.fa.{amb,ann,bwt,pac,sa}.\n" + "Pré-visualização (head):\n" + head,
                    border_style="red",
                )
            )
            raise RuntimeError("Tar sem arquivos de índice do BWA.")

        wanted = [re.sub(r"^\./", "", path) for path in wanted]
        names_quoted = " ".join(shlex.quote(path) for path in wanted)

        rc = sp.run(
            [
                "conda",
                "run",
                "-n",
                "genomics",
                "bash",
                "-lc",
                f"tar -xzf {shlex.quote(str(tar_path))} -C {shlex.quote(str(index_dir))} --overwrite {names_quoted}",
            ],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
        ).returncode
        if rc != 0:
            legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"tar -xzf {shlex.quote(str(tar_path))} -C {shlex.quote(str(index_dir))} {names_quoted}"])
        extracted = True

        amb = next(index_dir.rglob("*.fa.amb"), None)
        if not amb:
            raise RuntimeError("Falha ao localizar *.fa.amb após extração.")
        best_prefix = amb.with_suffix("")
        if not _prefix_has_all(best_prefix):
            raise RuntimeError("Índice BWA incompleto após extração.")

    from os.path import samefile

    for ext in exts:
        target = Path(str(best_prefix) + ext).resolve()
        if not target.exists():
            raise RuntimeError(f"Arquivo do índice ausente: {target}")
        dst = Path(str(ref_prefix) + ext)
        if dst.is_symlink():
            try:
                if samefile(dst, target):
                    pass
                else:
                    dst.unlink(missing_ok=True)
                    dst.symlink_to(target)
            except FileNotFoundError:
                dst.unlink(missing_ok=True)
                dst.symlink_to(target)
        elif dst.exists():
            dst.unlink(missing_ok=True)
            dst.symlink_to(target)
        else:
            dst.symlink_to(target)

    if extracted:
        legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"chmod -R a+r {shlex.quote(str(index_dir))} || true"])

    if not all(Path(str(ref_prefix) + ext).exists() for ext in exts):
        raise RuntimeError("Após relink, o índice BWA continua incompleto.")
    legacy.console.print("Índice BWA linkado → [bold]pronto[/bold]")
    legacy.print_meta("Índice BWA linkado", [Path(str(ref_prefix) + ext) for ext in exts])


def _download_and_place(url: str, dst_plain: Path):
    tmp = dst_plain.with_suffix(dst_plain.suffix + ".tmp")
    legacy.run(["wget", "-O", str(tmp), url])

    is_targz = (
        sp.run(
            ["conda", "run", "-n", "genomics", "bash", "-lc", f"tar -tzf {tmp} >/dev/null 2>&1"],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
        ).returncode
        == 0
    )

    if is_targz:
        legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"mkdir -p refs/_extracted && tar -xzf {tmp} -C refs/_extracted"])
        fa_candidates = list(Path("refs/_extracted").rglob("*.fa")) + list(Path("refs/_extracted").rglob("*.fasta"))
        if not fa_candidates:
            raise RuntimeError("Nenhum FASTA (.fa/.fasta) encontrado dentro do tar.")
        fa = max(fa_candidates, key=lambda p: p.stat().st_size)
        if dst_plain.exists():
            dst_plain.unlink()
        dst_plain.symlink_to(fa.resolve())
        tmp.unlink(missing_ok=True)
        return

    is_gzip = (
        sp.run(
            ["conda", "run", "-n", "genomics", "bash", "-lc", f"file -b --mime-type {tmp} | grep -qi 'gzip'"],
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
        ).returncode
        == 0
    )
    if is_gzip:
        legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"gzip -dc {tmp} > {dst_plain}"])
        tmp.unlink(missing_ok=True)
        return

    legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"mv {tmp} {dst_plain}"])


def download_refs(ref_fa_url, gtf_url, force=False):
    Path("refs").mkdir(exist_ok=True)
    fa = Path("refs/reference.fa")
    gtf = Path("refs/genes.gtf")

    if fa.exists() and not force:
        legacy.console.print(":floppy_disk: Referência já presente → [bold]SKIP[/bold]")
    else:
        _download_and_place(ref_fa_url, fa)

    if gtf.exists() and not force:
        legacy.console.print(":floppy_disk: Anotação (GTF) já presente → [bold]SKIP[/bold]")
    else:
        _download_and_place(gtf_url, gtf)

    if not Path("refs/reference.fa.fai").exists() or force:
        legacy.run(["samtools", "faidx", "refs/reference.fa"])
    if not Path("refs/reference.dict").exists() or force:
        legacy.run(["gatk", "CreateSequenceDictionary", "-R", "refs/reference.fa", "-O", "refs/reference.dict"])

    legacy.print_meta("Referências", [fa, gtf])


def limit_reference_to_canonical_if_enabled():
    g = legacy.cfg_global["general"]
    if not g.get("limit_to_canonical", False):
        return
    if not Path("refs/reference.fa.fai").exists():
        legacy.run(["samtools", "faidx", "refs/reference.fa"])
    names = [line.split("\t")[0] for line in open("refs/reference.fa.fai")]
    if not names:
        legacy.console.print("[orange3]Aviso:[/orange3] .fai vazio, mantendo referência completa.")
        return
    chr_prefix = names[0].startswith("chr")
    canon = [f"{'chr' if chr_prefix else ''}{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    canon += [f"{'chr' if chr_prefix else ''}M", f"{'chr' if chr_prefix else ''}MT"]
    canon = [contig for contig in canon if contig in names]
    if not canon:
        legacy.console.print("[orange3]Aviso:[/orange3] Não foi possível detectar contigs canônicos.")
        return
    legacy.run(["conda", "run", "-n", "genomics", "bash", "-lc", f"samtools faidx refs/reference.fa {' '.join(canon)} > refs/reference.canonical.fa"])
    Path("refs/reference.fa").unlink(missing_ok=True)
    Path("refs/reference.fa").symlink_to(Path("refs/reference.canonical.fa").resolve())
    legacy.run(["samtools", "faidx", "refs/reference.fa"])
    legacy.run(["gatk", "CreateSequenceDictionary", "-R", "refs/reference.fa", "-O", "refs/reference.dict"])
    legacy.console.print("[green]Referência reduzida aos cromossomos canônicos.[/green]")

__all__ = [
    "build_indexes",
    "download_refs",
    "install_prebuilt_bwa_index",
    "limit_reference_to_canonical_if_enabled",
]
