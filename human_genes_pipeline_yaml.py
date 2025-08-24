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

console = Console(highlight=False)

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
    Executa pipeline de processos encadeados (Popen) com heartbeats.
    Ex.: [["bwa-mem2","mem",...], ["samtools","view","-bS","-"], ["samtools","sort","-o","..."]]
    """
    console.print(Panel.fit(Text(label, style="bold yellow"), border_style="yellow"))
    console.print("[bold]>[/bold] " + "  |  ".join(" ".join(c) for c in cmd_list), style="dim")

    procs = []
    start = time.time()

    # encadeia stdout->stdin
    for i, cmd in enumerate(cmd_list):
        if i == 0:
            p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, text=True, bufsize=1, universal_newlines=True)
        else:
            p = sp.Popen(cmd, stdin=procs[-1].stdout, stdout=sp.PIPE if i < len(cmd_list)-1 else None)
        procs.append(p)

    last = start
    # Lê a saída do primeiro processo (que costuma imprimir progresso)
    if procs and procs[0].stdout:
        for line in procs[0].stdout:
            if line.strip() and any(k in line for k in ["WARNING","ERROR","ERR","progress","reads","mapped","trimmed","Time","ETA","%","Mbp","Gbp"]):
                console.print(line.rstrip(), style="dim")
            now = time.time()
            if now - last > heartbeat_sec:
                elapsed = int(now - start)
                console.print(f"[heartbeat] {label}... {elapsed//60}m{elapsed%60:02d}s", style="magenta")
                last = now

    # aguarda todos
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

def bytes_free(p: Path) -> int:
    return shutil.disk_usage(p).free

# ================== Normalização do YAML ==================

def normalize_config_schema(cfg_in: dict) -> dict:
    """
    Aceita:
      (A) antigo: {'general': {...}, 'dna_samples': [...], 'rna_samples': [...]}
      (B) novo  : {'project': {...}, 'storage': {...}, 'download': {...},
                   'execution': {...}, 'size_control': {...}, 'samples': [...], 'params': {...}}
      No (B), a referência pode aparecer como 'reference' no topo OU dentro de 'project.reference'.
    Retorna esquema (A).
    """
    if "general" in cfg_in:
        return cfg_in

    project   = cfg_in.get("project", {})
    storage   = cfg_in.get("storage", {})
    ref_top   = cfg_in.get("reference", {})
    ref_proj  = project.get("reference", {}) if isinstance(project, dict) else {}
    ref       = {**ref_proj, **ref_top}  # proj > topo
    download  = cfg_in.get("download", {})
    execv     = cfg_in.get("execution", {})
    sizec     = cfg_in.get("size_control", {})
    samples   = cfg_in.get("samples", [])
    params    = cfg_in.get("params", {})

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
    }

    return {"general": general, "dna_samples": dna_samples, "rna_samples": []}

# ================ Estimativa de espaço =================

def estimate_outputs_and_warn(fastqs: List[Path], base_dir: Path, guard_gb: int):
    free_b = bytes_free(base_dir)
    in_bytes = sum((file_meta(f)["size_bytes"] for f in fastqs if f.exists()))
    trimmed_est = int(in_bytes * 0.8)
    bam_est     = int(in_bytes * 1.2)
    vcf_est     = int(2e9)  # ~2 GB por amostra (chute)
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
    tmp = dst_plain.with_suffix(dst_plain.suffix + ".tmp")
    run(["wget", "-O", str(tmp), url])
    # Se URL termina com .gz, descompacta; senão, tenta detectar gzip via `file`
    if str(url).endswith(".gz"):
        gz = tmp if tmp.suffix == ".gz" else Path(str(tmp) + ".gz")
        if gz != tmp: tmp.rename(gz)
        run(["gunzip", "-f", str(gz)])
    else:
        try:
            run(["bash", "-lc", f"file -b {tmp} | grep -qi gzip && gunzip -f {tmp} || mv {tmp} {dst_plain}"])
        except sp.CalledProcessError:
            tmp.rename(dst_plain)

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

def build_indexes(default_read_type, assembly_name, need_rna_index, threads, force=False):
    if default_read_type=="short" and (force or not Path("refs/reference.fa.bwt.2bit.64").exists()):
        run(["bwa-mem2","index","refs/reference.fa"])
    else:
        console.print("Índice BWA → [bold]SKIP[/bold]", style="dim")
    if need_rna_index and (force or not Path(f"refs/{assembly_name}.1.ht2").exists()):
        run(["hisat2-build","-p",str(threads),"refs/reference.fa",f"refs/{assembly_name}"])
    elif need_rna_index:
        console.print("Índice HISAT2 → [bold]SKIP[/bold]", style="dim")

# ================== Entrada SRA / FASTQ ===================

def stage_fastqs_from_sra(sra_ids: List[str]):
    Path("raw").mkdir(exist_ok=True); Path("fastq").mkdir(exist_ok=True)
    for acc in sra_ids:
        sra_target = next(Path("raw").glob(f"**/{acc}.sra"), None)
        if sra_target and sra_target.exists():
            console.print(f"{acc}: .sra → [bold]SKIP (cache)[/bold]")
        else:
            run(["prefetch","--output-directory",str(Path("raw").resolve()), acc])
        r1 = Path("fastq")/f"{acc}_1.fastq.gz"
        r2 = Path("fastq")/f"{acc}_2.fastq.gz"
        if r1.exists() or r2.exists():
            console.print(f"{acc}: FASTQ → [bold]SKIP (cache)[/bold]")
        else:
            # conversão pode ser lenta; usa heartbeat interno
            console.print(Panel.fit(Text(f"Convertendo {acc} (fasterq-dump)", style="bold yellow"), border_style="yellow"))
            run(["fasterq-dump","--split-files","--gzip",acc,"-O",str(Path("fastq").resolve())])
        print_meta(f"FASTQs ({acc})", [r1, r2])

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
            run_long_stream_pipeline(
                [
                    ["cutadapt","-j",str(threads),"-q","20,20","-m","30","-a",fwd,"-A",rev,str(r1),str(r2)],
                    ["tee","/dev/stderr"],  # para ver logs de progresso do cutadapt
                    ["bash","-lc", f"cat > /dev/null"]  # placeholder para manter pipeline visível
                ],
                label=f"Cutadapt (paired) {base}"
            )
            # salvamos via parâmetros -o/-p; então só exibimos metadados:
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
                # preferir trimmed; senão usar fastq_ds (se houver) ou fastq
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
            expected_r2 = Path(s.get("fastq2","")).name if s.get("fastq2") else None
            align_one_sample(r1, threads, rtype, cleanup, expected_r2_name=expected_r2)

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

    if read_type == "short":
        if r2 and r2.exists():
            cmd_list = [
                ["bwa-mem2","mem","-t",str(threads),"refs/reference.fa",str(r1),str(r2)],
                ["samtools","view","-bS","-"],
                ["samtools","sort","-@",str(threads),"-o",str(out_sorted)]
            ]
            run_long_stream_pipeline(cmd_list, label=f"[{sample}] BWA-MEM2 (paired) → sort")
        else:
            cmd_list = [
                ["bwa-mem2","mem","-t",str(threads),"refs/reference.fa",str(r1)],
                ["samtools","view","-bS","-"],
                ["samtools","sort","-@",str(threads),"-o",str(out_sorted)]
            ]
            run_long_stream_pipeline(cmd_list, label=f"[{sample}] BWA-MEM2 (single) → sort")

        run(["samtools","index",str(out_sorted)])
        run(["picard","MarkDuplicates","-I",str(out_sorted),"-O",str(out_mkdup),
             "-M",str(out_sorted).replace(".sorted.bam",".mkdup.metrics"),"--VALIDATION_STRINGENCY","LENIENT"])
        run(["samtools","index",str(out_mkdup)])
        if cleanup.get("remove_sorted_bam", True) and out_sorted.exists():
            out_sorted.unlink(missing_ok=True)

    else:
        # long reads (ONT/PacBio) → minimap2 | sort
        cmd_list = [
            ["minimap2","-ax","map-ont","-t",str(threads),"refs/reference.fa",str(r1)],
            ["samtools","sort","-@",str(threads),"-o",str(out_sorted)]
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
            cmd_list = [
                ["hisat2","-p",str(threads),"-x",f"refs/{assembly_name}","-1",str(r1),"-2",str(r2)],
                ["samtools","sort","-@",str(threads),"-o",f"rnaseq/{base}.rnaseq.bam"]
            ]
        else:
            cmd_list = [
                ["hisat2","-p",str(threads),"-x",f"refs/{assembly_name}","-U",str(r1)],
                ["samtools","sort","-@",str(threads),"-o",f"rnaseq/{base}.rnaseq.bam"]
            ]
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

    console.rule(f"[bold]Pipeline Humano[/bold]  →  [cyan]{base_dir}[/cyan]")
    disk_space_report(base_dir)

    ensure_dirs()

    # Referências e índices
    download_refs(g["ref_fa_url"], g["gtf_url"], force=g.get("force_refs", False))
    build_indexes(g.get("default_read_type","short"), g["assembly_name"],
                  g.get("rnaseq",False) or bool(cfg.get("rna_samples",[])),
                  g["threads"], force=g.get("force_indexes", False))

    # Entrada DNA
    dna_samples = cfg.get("dna_samples", [])
    for s in dna_samples:
        if s["source"]=="sra": stage_fastqs_from_sra(s["sra_ids"])
        else: stage_fastqs_from_local(s["fastq1"], s.get("fastq2"))

    # Estimativa de espaço (com base no que já baixou)
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

    # Remover BAM após CRAM (economia de espaço)
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
    ap = argparse.ArgumentParser(description="Pipeline humano (YAML) — Verbose, cores, cache, CRAM, cobertura e downsample")
    ap.add_argument("-c","--config",required=True,help="Arquivo YAML de configuração")
    args = ap.parse_args()
    with open(args.config, "r") as f:
        cfg_raw = yaml.safe_load(f)
    cfg = normalize_config_schema(cfg_raw)
    console.print(Panel.fit(f"[bold]YAML normalizado[/bold]\nChaves topo: {list(cfg.keys())}\nGeneral → {', '.join(sorted(list(cfg['general'].keys()))[:10])}...", border_style="cyan"))
    main(cfg)

