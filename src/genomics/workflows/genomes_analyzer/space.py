"""Disk-space estimation helpers."""

import shutil

from rich.panel import Panel


def estimate_outputs_and_warn(fastqs, base_dir, guard_gb: int):
    from . import legacy

    free_b = shutil.disk_usage(base_dir).free
    in_bytes = sum((legacy.file_meta(f)["size_bytes"] for f in fastqs if f.exists()))
    trimmed_est = int(in_bytes * 0.8)
    bam_est = int(in_bytes * 1.2)
    vcf_est = int(2e9)
    needed_bam_path = bam_est + trimmed_est + vcf_est
    msg = (
        f"Entrada FASTQ (comprimido): ~{legacy.sizeof_fmt(in_bytes)}\n"
        f"Pós-trimming: ~{legacy.sizeof_fmt(trimmed_est)}; BAM: ~{legacy.sizeof_fmt(bam_est)}; "
        f"VCF: ~{legacy.sizeof_fmt(vcf_est)}\n"
        f"Necessário (sem CRAM/limpeza): ~{legacy.sizeof_fmt(needed_bam_path)} | "
        f"Livre: {legacy.sizeof_fmt(free_b)}"
    )
    style = "yellow" if (free_b < needed_bam_path or free_b < guard_gb * 1024**3) else "green"
    legacy.console.print(Panel.fit(msg, title="Estimativa de espaço", border_style=style))
    if style == "yellow":
        legacy.console.print(
            "[orange3]Sugestões:[/orange3] use [bold]use_cram: true[/bold], "
            "[bold]cleanup.remove_sorted_bam: true[/bold] e/ou [bold]downsample_frac[/bold].",
            style="italic",
        )

__all__ = ["estimate_outputs_and_warn"]
