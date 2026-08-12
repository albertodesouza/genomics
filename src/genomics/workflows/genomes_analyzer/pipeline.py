"""Top-level orchestration for the human genomics pipeline."""

import os
from pathlib import Path

from . import legacy
from .alignment import align_dna_for_all, to_cram_and_coverage
from .ancestry import ancestry_admixture_step
from .config import normalize_config_schema
from .fastq import downsample_fastqs, stage_fastqs_from_local, stage_fastqs_from_sra
from .qc import qc_and_trim
from .references import build_indexes, download_refs, limit_reference_to_canonical_if_enabled
from .reports import (
    combine_pairwise_stats,
    gene_list_from_gtf,
    gene_presence_reports,
    pairwise_comparisons,
    paternity_analysis,
    trio_denovo_report,
)
from .rnaseq import rnaseq_pipeline
from .space import estimate_outputs_and_warn
from .state import set_config
from .utils import disk_space_report, ensure_dirs, print_meta, stage_banner, step_done
from .variants import annotate_with_vep, call_variants


def _setup_runtime(cfg):
    set_config(cfg)
    g = cfg["general"]

    legacy._configure_console_for_mode(cfg.get("params", {}))

    base_dir = Path(g.get("base_dir", ".")).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(base_dir)

    tmpdir = Path(g.get("temp_dir", "tmp")).expanduser().resolve()
    tmpdir.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(tmpdir)
    legacy.console.print(f"Usando TMPDIR={tmpdir} (do YAML)", style="dim")

    legacy.console.rule(f"[bold]Pipeline Humano[/bold]  →  [cyan]{base_dir}[/cyan]")
    disk_space_report(base_dir)
    ensure_dirs()
    return base_dir


def _run_refs_and_indexes(cfg, g):
    stage_banner("1) Referências e Índices")
    download_refs(g["ref_fa_url"], g["gtf_url"], force=g.get("force_refs", False))
    limit_reference_to_canonical_if_enabled()
    build_indexes(
        g.get("default_read_type", "short"),
        g["assembly_name"],
        g.get("rnaseq", False) or bool(cfg.get("rna_samples", [])),
        g["threads"],
        force=g.get("force_indexes", False),
    )
    step_done("Referências/índices prontos")


def _run_fetch_fastqs(dna_samples):
    stage_banner("2) Entrada de Leituras (FASTQ)", "SRA→FASTQ, cache e compressão")
    for sample in dna_samples:
        if sample["source"] == "sra":
            stage_fastqs_from_sra(sample["sra_ids"])
        else:
            stage_fastqs_from_local(sample["fastq1"], sample.get("fastq2"))
    step_done("FASTQs prontos")


def _run_estimate_space(base_dir, g):
    stage_banner("3) Estimativa de Espaço")
    fq_list = list(Path("fastq").glob("*.fastq.gz"))
    estimate_outputs_and_warn(fq_list, base_dir, g.get("space_guard_gb_min", 100))
    step_done("Estimativa impressa")


def _run_downsample(g):
    stage_banner("4) Downsample (opcional)")
    downsample_fastqs(g["downsample_frac"], g.get("downsample_seed", 123))
    step_done("Downsample concluído")


def _run_qc_and_trimming(g):
    stage_banner("5) QC e Trimming")
    qc_and_trim(
        g["threads"],
        g["adapters"],
        g.get("default_read_type", "short"),
        use_ds=g.get("downsample_frac", 0.0) > 0,
    )
    step_done("QC + trimming concluídos")


def _run_align_and_sort(dna_samples, g):
    stage_banner("6) Alinhamento, Sort e Duplicatas")
    align_dna_for_all(
        dna_samples,
        g["threads"],
        g.get("default_read_type", "short"),
        cleanup=g.get("cleanup", {"remove_sorted_bam": True, "remove_bam_after_cram": True}),
        use_ds=g.get("downsample_frac", 0.0) > 0,
    )
    step_done("Alinhamento e marcação de duplicatas concluídos")


def _run_cram_and_coverage(g):
    stage_banner("7) CRAM e Cobertura (mosdepth)")
    to_cram_and_coverage(g.get("use_cram", True), g["threads"])
    if g.get("use_cram", True) and g.get("cleanup", {}).get("remove_bam_after_cram", True):
        for bam in Path("bam").glob("*.mkdup.bam"):
            legacy.console.print(f"Removendo {bam.name} (CRAM mantido)", style="dim")
            bam.unlink(missing_ok=True)
            idx = Path(str(bam) + ".bai")
            idx.unlink(missing_ok=True)
    step_done("CRAM/index e cobertura prontos")


def _run_variants_and_vep(dna_samples, g):
    stage_banner("8) Chamadas de Variantes e Anotação")
    if g.get("call_variants", True):
        call_variants(dna_samples, g["threads"], g["mem_gb"])
    if g.get("annotate_vars", True):
        annotate_with_vep(dna_samples, g["threads"])
    step_done("VCF e VEP-VCF gerados")


def _run_gene_list():
    stage_banner("9) Lista de Genes (Referência)")
    gene_list_from_gtf()
    step_done("Lista de genes da referência pronta")


def _run_gene_presence(dna_samples, g):
    stage_banner("10) Presença de genes por amostra (cobertura por gene)")
    gene_presence_reports(
        dna_samples,
        min_mean_cov=float(g.get("gene_presence_min_mean_cov", 5.0)),
        min_breadth_1x=float(g.get("gene_presence_min_breadth_1x", 0.8)),
        threads=g["threads"],
    )
    print_meta("Genes BED (alvos p/ cobertura)", [Path("genes/genes.bed")])
    pres = sorted(Path("genes").glob("*_gene_presence.tsv"))
    if pres:
        print_meta("Presença de genes (por amostra)", pres)
    step_done("Relatórios de presença de genes gerados")


def _run_pairwise(dna_samples):
    stage_banner("11) Comparações par-a-par (vcf)")
    pairwise_comparisons(dna_samples)
    step_done("Relatórios de comparação gerados")
    combine_pairwise_stats(dna_samples)


def _run_trio_denovo(dna_samples):
    stage_banner("12) Trio: candidatos de novo")
    trio_denovo_report(dna_samples)
    step_done("Relatório trio gerado")


def _run_paternity(dna_samples):
    stage_banner("13) Paternidade (trio, SNPs)")
    paternity_analysis(dna_samples)
    step_done("Paternidade concluída")


def _run_ancestry(cfg):
    stage_banner("14) Ancestralidade (ADMIXTURE supervisionado)")
    ancestry_admixture_step(cfg)
    step_done("Ancestralidade concluída")


def _run_rnaseq(cfg, g):
    stage_banner("15) RNA-seq (opcional)")
    rnaseq_pipeline(cfg.get("rna_samples", []), g["threads"], g["assembly_name"])
    step_done("RNA-seq concluído")


def main(cfg):
    g = cfg["general"]
    steps = cfg.get("steps", {})
    base_dir = _setup_runtime(cfg)

    if steps.get("refs_and_indexes", True):
        _run_refs_and_indexes(cfg, g)

    dna_samples = cfg.get("dna_samples", [])
    if steps.get("fetch_fastqs", True):
        _run_fetch_fastqs(dna_samples)

    if steps.get("estimate_space", True):
        _run_estimate_space(base_dir, g)

    if steps.get("downsample", False) and (g.get("downsample_frac", 0.0) > 0):
        _run_downsample(g)

    if steps.get("qc_and_trimming", True):
        _run_qc_and_trimming(g)

    if steps.get("align_and_sort", True):
        _run_align_and_sort(dna_samples, g)

    if steps.get("cram_and_coverage", True):
        _run_cram_and_coverage(g)

    if steps.get("variants_and_vep", True):
        _run_variants_and_vep(dna_samples, g)

    if steps.get("gene_list", True):
        _run_gene_list()

    if steps.get("gene_presence", True):
        _run_gene_presence(dna_samples, g)

    if steps.get("pairwise", True):
        _run_pairwise(dna_samples)

    if steps.get("trio_denovo", True):
        _run_trio_denovo(dna_samples)

    if steps.get("paternity", True):
        _run_paternity(dna_samples)

    if steps.get("ancestry", False):
        _run_ancestry(cfg)

    if steps.get("rnaseq", False) and g.get("rnaseq", False):
        _run_rnaseq(cfg, g)

    legacy.console.rule("[bold green]Pipeline concluído ✅[/bold green]")


__all__ = ["main", "normalize_config_schema"]
