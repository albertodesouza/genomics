"""Configuration loading and schema normalization."""


def normalize_config_schema(cfg_in: dict) -> dict:
    """
    Accepts the legacy schema and the newer project/storage/reference schema.

    Returns the normalized legacy shape used by the pipeline internals:
    ``general``, ``dna_samples``, ``rna_samples``, ``params``, ``steps`` and
    ``ancestry``.
    """
    if "general" in cfg_in and ("dna_samples" in cfg_in or "rna_samples" in cfg_in):
        return cfg_in

    project = cfg_in.get("project", {})
    storage = cfg_in.get("storage", {})
    ref_top = cfg_in.get("reference", {})
    ref_proj = project.get("reference", {}) if isinstance(project, dict) else {}
    ref = {**ref_proj, **ref_top}
    execv = cfg_in.get("execution", {})
    download = cfg_in.get("download", {})
    sizec = cfg_in.get("size_control", {})
    samples = cfg_in.get("samples", [])
    params = cfg_in.get("params", {})

    base_dir = storage.get("base_dir", ".")
    temp_dir = storage.get("temp_dir", "tmp")
    assembly = ref.get("name", "GRCh38")
    ref_fa_url = (
        ref.get("fasta_url")
        or "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
    )
    gtf_url = (
        ref.get("gtf_url")
        or "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz"
    )
    bwa_prebuilt = ref.get("bwa_index_url") or cfg_in.get("bwa_prebuilt_url") or ""

    prefetch_retries = int(download.get("prefetch_retries", 2))
    prefer_ena_fastq = bool(execv.get("prefer_ena_fastq", False))
    cancel_on_convert_stall = bool(execv.get("cancel_on_convert_stall", False))
    stall_warn_min = int(execv.get("stall_warn_min", 10))
    stall_fail_min = int(execv.get("stall_fail_min", 45))
    ena_fallback = bool(execv.get("ena_fallback", False))
    download_tool = (download.get("tool") or "sra_toolkit").lower()
    if download_tool != "sra_toolkit":
        from . import legacy

        legacy.console.print(
            f"[orange3]Aviso:[/orange3] download.tool='{download_tool}' ainda não é suportado; usando sra-tools.",
            style="italic",
        )

    threads_global = int(execv.get("threads") or download.get("threads") or 16)
    mem_gb = int(params.get("mem_gb", 64))

    ds_cfg = (sizec or {}).get("downsample", {})
    downsample_frac = float(ds_cfg.get("fraction", 0.0)) if ds_cfg.get("enabled", False) else 0.0
    downsample_seed = int(ds_cfg.get("seed", 123))

    keep_inter = bool(storage.get("keep_intermediates", False))
    use_cram = bool(cfg_in.get("general", {}).get("use_cram", True))
    cleanup = dict(cfg_in.get("general", {}).get("cleanup", {})) or {
        "remove_sorted_bam": True,
        "remove_bam_after_cram": (not keep_inter),
    }

    adapters = {"fwd": "AGATCGGAAGAGC", "rev": "AGATCGGAAGAGC"}

    aligner = (params.get("aligner") or execv.get("aligner") or "bwa-mem2").lower()
    aligner = "bwa-mem2" if aligner in ("bwa-mem2", "bwa_mem2", "mem2") else "bwa"

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

    sort_mem_mb = int(cfg_in.get("general", {}).get("sort_mem_mb", 384))
    bwa_batch_k = int(cfg_in.get("general", {}).get("bwa_batch_k", 20000000))
    limit_to_canonical = bool(cfg_in.get("limit_to_canonical", False))

    dna_samples = []
    for sample in samples:
        runs = sample.get("runs", [])
        if not runs:
            continue
        dna_samples.append(
            {
                "id": sample.get("sample_id") or sample.get("id") or runs[0],
                "source": "sra",
                "sra_ids": runs,
                "read_type": "short",
            }
        )

    general = {
        "base_dir": base_dir,
        "assembly_name": assembly,
        "ref_fa_url": ref_fa_url,
        "gtf_url": gtf_url,
        "threads": threads_global,
        "aln_threads": aln_threads,
        "mem_gb": mem_gb,
        "default_read_type": "short",
        "adapters": adapters,
        "call_variants": True,
        "annotate_vars": True,
        "rnaseq": False,
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
        "sort_mem_mb": sort_mem_mb,
        "bwa_batch_k": bwa_batch_k,
    }

    opt_keys = [
        "trio_child_id",
        "trio_parent_ids",
        "trio_min_dp_child",
        "trio_min_dp_parents",
        "trio_min_gq",
        "trio_min_ab_het",
        "trio_max_ab_het",
        "trio_min_ab_hom",
        "trio_max_parent_alt_frac",
        "gene_presence_min_mean_cov",
        "gene_presence_min_breadth_1x",
        "stall_warn_min",
        "stall_fail_min",
        "cancel_on_convert_stall",
        "prefer_ena_fastq",
        "ena_fallback",
        "prefetch_retries",
        "temp_dir",
        "sort_mem_mb",
        "bwa_batch_k",
        "force_refs",
        "force_indexes",
        "use_cram",
        "cleanup",
        "paternity_prior",
        "paternity_epsilon",
        "paternity_require_pass",
        "paternity_force_pass",
        "paternity_use_vep_af",
        "paternity_skip_all_hets",
        "aligner",
        "bwa_prebuilt_url",
        "limit_to_canonical",
        "sort_threads",
        "hc_threads_native",
        "space_guard_gb_min",
        "downsample_frac",
        "downsample_seed",
        "call_variants",
        "annotate_vars",
        "rnaseq",
    ]
    for src in (execv, storage, cfg_in.get("general", {}), cfg_in):
        if isinstance(src, dict):
            for key in opt_keys:
                if key in src:
                    general[key] = src[key]

    steps_cfg = dict(cfg_in.get("steps", {})) if isinstance(cfg_in.get("steps", {}), dict) else {}

    return {
        "general": general,
        "dna_samples": dna_samples,
        "rna_samples": [],
        "params": dict(params),
        "steps": steps_cfg,
        "ancestry": dict(cfg_in.get("ancestry", {})),
    }

__all__ = ["normalize_config_schema"]
