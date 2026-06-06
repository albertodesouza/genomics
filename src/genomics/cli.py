from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

import yaml

from genomics.core.data_registry import classify_dataset_path, registered_dataset_ids, resolve_dataset
from genomics.workspace import (
    DEFAULT_CONSENSUS_DATASET_DIR,
    DEFAULT_DATASET_DIR,
    DEFAULT_GENOTYPE_RUNS_ROOT,
    DEFAULT_VARIANT_TRANSFORMER_DATASET_DIR,
    results_path,
)


PathLike = Union[str, Path]


def _run_module(module: str, args: Iterable[PathLike]) -> int:
    command = [sys.executable, "-m", module, *[str(arg) for arg in args]]
    try:
        return subprocess.call(command)
    finally:
        _restore_terminal()


def _restore_terminal() -> None:
    if not sys.stdout.isatty():
        return
    try:
        sys.stdout.write("\033[0m\033[?25h")
        sys.stdout.flush()
    except Exception:
        pass
    try:
        with open("/dev/tty") as tty:
            subprocess.run(["stty", "sane"], stdin=tty, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass


def _load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Config YAML invalido: {path}")
    return payload


def _write_temp_yaml(payload: Dict[str, Any], prefix: str) -> Path:
    temp_dir = Path(tempfile.mkdtemp(prefix=prefix))
    path = temp_dir / "config.yaml"
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)
    return path


def _set_if_present(mapping: Dict[str, Any], key: str, value: Optional[PathLike]) -> None:
    if value is not None:
        mapping[key] = str(Path(value).resolve())


def _resolve_dataset_id_field(mapping: Dict[str, Any], dir_key: str = "dataset_dir") -> None:
    dataset_id = mapping.get("dataset_id") or mapping.get("source_dataset_id")
    if dataset_id and not mapping.get(dir_key):
        mapping[dir_key] = str(resolve_dataset(str(dataset_id)).path)


def _absolutize_existing_path(payload: Dict[str, Any], section: str, key: str, base_dir: Path) -> None:
    mapping = payload.get(section)
    if not isinstance(mapping, dict):
        return
    value = mapping.get(key)
    if not value:
        return
    path = Path(value)
    if not path.is_absolute():
        mapping[key] = str((base_dir / path).resolve())


def _genotype_config_with_overrides(args: argparse.Namespace) -> Path:
    source_config = Path(args.config).resolve()
    payload = _load_yaml(source_config)
    _absolutize_existing_path(payload, "dataset_input", "view_path", source_config.parent)
    dataset_input = payload.setdefault("dataset_input", {})
    _resolve_dataset_id_field(dataset_input, "dataset_dir")
    _set_if_present(dataset_input, "dataset_dir", getattr(args, "dataset_dir", None))
    if getattr(args, "dataset_id", None):
        dataset_input["dataset_id"] = args.dataset_id
        if not getattr(args, "dataset_dir", None):
            dataset_input.pop("dataset_dir", None)
            _resolve_dataset_id_field(dataset_input, "dataset_dir")
    _set_if_present(dataset_input, "view_path", getattr(args, "view_path", None))
    _set_if_present(dataset_input, "processed_cache_dir", getattr(args, "processed_cache_dir", None))
    _set_if_present(dataset_input, "results_dir", getattr(args, "results_dir", None))
    _set_if_present(dataset_input, "consensus_dataset_dir", getattr(args, "consensus_dataset_dir", None))
    return _write_temp_yaml(payload, "genomics_genotype_config_")


def _variant_config_with_overrides(args: argparse.Namespace) -> Path:
    payload = _load_yaml(Path(args.config).resolve())
    dataset = payload.setdefault("dataset", {})
    _resolve_dataset_id_field(dataset, "source_dataset_dir")
    _set_if_present(dataset, "processed_dir", getattr(args, "processed_dir", None))
    if getattr(args, "source_dataset_id", None):
        dataset["source_dataset_id"] = args.source_dataset_id
        dataset.pop("source_dataset_dir", None)
        _resolve_dataset_id_field(dataset, "source_dataset_dir")
    _set_if_present(dataset, "results_dir", getattr(args, "results_dir", None))
    return _write_temp_yaml(payload, "genomics_variant_config_")


def _neural_config_with_overrides(args: argparse.Namespace) -> Path:
    payload = _load_yaml(Path(args.config).resolve())
    dataset_input = payload.setdefault("dataset_input", {})
    _resolve_dataset_id_field(dataset_input, "dataset_dir")
    _set_if_present(dataset_input, "dataset_dir", getattr(args, "dataset_dir", None))
    if getattr(args, "dataset_id", None):
        dataset_input["dataset_id"] = args.dataset_id
        if not getattr(args, "dataset_dir", None):
            dataset_input.pop("dataset_dir", None)
            _resolve_dataset_id_field(dataset_input, "dataset_dir")
    _set_if_present(dataset_input, "processed_cache_dir", getattr(args, "processed_cache_dir", None))
    return _write_temp_yaml(payload, "genomics_neural_config_")


def _add_genotype_config_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config", type=Path)
    parser.add_argument("--dataset-dir", type=Path, default=None)
    parser.add_argument("--dataset-id", default=None)
    parser.add_argument("--view-path", type=Path, default=None)
    parser.add_argument("--processed-cache-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--consensus-dataset-dir", type=Path, default=None)


def _add_variant_config_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config", type=Path)
    parser.add_argument("--processed-dir", type=Path, default=None)
    parser.add_argument("--source-dataset-id", default=None)
    parser.add_argument("--results-dir", type=Path, default=None)


def _add_neural_config_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config", type=Path)
    parser.add_argument("--dataset-dir", type=Path, default=None)
    parser.add_argument("--dataset-id", default=None)
    parser.add_argument("--processed-cache-dir", type=Path, default=None)


def cmd_genotype_prepare_cache(args: argparse.Namespace) -> int:
    from pathlib import Path

    from genomics.predictors.genotype_based.config import get_dataset_cache_dir, load_config
    from genomics.predictors.genotype_based.data.pipeline import prepare_cache_only

    config_path = _genotype_config_with_overrides(args)
    config = load_config(config_path)
    prepare_dir = Path(config.dataset_input.processed_cache_dir).resolve() / "prepare_cache"
    prepare_dir.mkdir(parents=True, exist_ok=True)
    prepare_cache_only(config, prepare_dir)
    print(f"Dataset cache: {get_dataset_cache_dir(config)}")
    return 0


def cmd_genotype_split(args: argparse.Namespace) -> int:
    return cmd_genotype_prepare_cache(args)


def cmd_genotype_train(args: argparse.Namespace) -> int:
    return _run_module("genomics.predictors.genotype_based.experiments.train", [_genotype_config_with_overrides(args)])


def cmd_genotype_evaluate(args: argparse.Namespace) -> int:
    config_path = _genotype_config_with_overrides(args)
    command_args: list[PathLike] = [config_path, "--checkpoint", args.checkpoint, "--split", args.split]
    if args.experiment_dir:
        command_args.extend(["--experiment-dir", args.experiment_dir])
    if args.output_name:
        command_args.extend(["--output-name", args.output_name])
    return _run_module("genomics.predictors.genotype_based.experiments.evaluate_checkpoint", command_args)


def cmd_genotype_test(args: argparse.Namespace) -> int:
    config_path = _genotype_config_with_overrides(args)
    command_args: list[PathLike] = [config_path, "--checkpoint", args.checkpoint, "--split", "test"]
    if args.experiment_dir:
        command_args.extend(["--experiment-dir", args.experiment_dir])
    if args.output_name:
        command_args.extend(["--output-name", args.output_name])
    else:
        command_args.extend(["--output-name", "test_best_accuracy"])
    return _run_module("genomics.predictors.genotype_based.experiments.evaluate_checkpoint", command_args)


def cmd_genotype_search(args: argparse.Namespace) -> int:
    return _run_module("genomics.predictors.genotype_based.experiments.search_sklearn", [_genotype_config_with_overrides(args)])


def cmd_genotype_pca_variance(args: argparse.Namespace) -> int:
    config_path = _genotype_config_with_overrides(args)
    command_args: list[PathLike] = [
        config_path,
        "--output",
        args.output,
        "--json-output",
        args.json_output,
    ]
    if args.max_components is not None:
        command_args.extend(["--max-components", str(args.max_components)])
    if args.force:
        command_args.append("--force")
    return _run_module("genomics.predictors.genotype_based.analysis.plot_sklearn_pca_variance", command_args)


def cmd_genotype_workbench(args: argparse.Namespace) -> int:
    command_args: list[PathLike] = [
        "--dataset-dir",
        args.dataset_dir,
        "--runs-root",
        args.runs_root,
        "--consensus-dataset-dir",
        args.consensus_dataset_dir,
        "--host",
        args.host,
        "--port",
        args.port,
    ]
    if args.aligned_tsv_root:
        command_args.extend(["--aligned-tsv-root", args.aligned_tsv_root])
    return _run_module("genomics.predictors.genotype_based.apps.genomics_workbench", command_args)


def cmd_genotype_sync_bcftools(args: argparse.Namespace) -> int:
    command_args: list[PathLike] = [
        "--source-dir",
        args.source_dir,
        "--target-dir",
        args.target_dir,
        "--link-mode",
        args.link_mode,
    ]
    if args.apply:
        command_args.append("--apply")
    if args.sample_limit is not None:
        command_args.extend(["--sample-limit", str(args.sample_limit)])
    if args.genes:
        command_args.append("--genes")
        command_args.extend(args.genes)
    if args.json:
        command_args.append("--json")
    return _run_module("genomics.predictors.genotype_based.tools.sync_bcftools_chain_artifacts", command_args)


def cmd_genotype_single_gene_screen(args: argparse.Namespace) -> int:
    config_path = _genotype_config_with_overrides(args)
    command_args: list[PathLike] = [config_path]
    for name in ["dataset_dir", "genes_file", "output_root", "run_name_prefix"]:
        value = getattr(args, name)
        if value is not None:
            command_args.extend([f"--{name.replace('_', '-')}", value])
    if args.genes:
        command_args.append("--genes")
        command_args.extend(args.genes)
    if args.split_by_ontology:
        command_args.append("--split-by-ontology")
    if args.ontology:
        for ontology in args.ontology:
            command_args.extend(["--ontology", ontology])
    if args.force:
        command_args.append("--force")
    if args.dry_run:
        command_args.append("--dry-run")
    return _run_module("genomics.predictors.genotype_based.experiments.run_single_gene_screen", command_args)


def cmd_variant_materialize(args: argparse.Namespace) -> int:
    output_dir = args.output_dir or DEFAULT_VARIANT_TRANSFORMER_DATASET_DIR
    if args.dataset_id and not args.dataset_dir:
        args.dataset_dir = resolve_dataset(args.dataset_id).path
    command_args: list[PathLike] = ["--output-dir", output_dir]
    for name in [
        "dataset_dir",
        "regions_bed",
        "samples_metadata",
        "vcf_pattern",
        "vcf_root_dir",
        "target",
        "l_max",
        "max_indel_size",
        "sample_batch_size",
        "train_split",
        "val_split",
        "test_split",
        "random_seed",
        "family_split_mode",
        "max_samples",
        "central_window_size",
    ]:
        value = getattr(args, name)
        if value is not None:
            command_args.extend([f"--{name.replace('_', '-')}", str(value)])
    if args.classes:
        command_args.append("--classes")
        command_args.extend(args.classes)
    if args.genes:
        command_args.append("--genes")
        command_args.extend(args.genes)
    return _run_module("genomics.predictors.variant_transformer.materialize_dataset", command_args)


def cmd_variant_train(args: argparse.Namespace) -> int:
    return _run_module("genomics.predictors.variant_transformer.train", [_variant_config_with_overrides(args)])


def cmd_variant_evaluate(args: argparse.Namespace) -> int:
    config_path = _variant_config_with_overrides(args)
    command_args: list[PathLike] = [config_path, "--checkpoint", args.checkpoint, "--split", args.split]
    if args.experiment_dir:
        command_args.extend(["--experiment-dir", args.experiment_dir])
    return _run_module("genomics.predictors.variant_transformer.evaluate_checkpoint", command_args)


def cmd_variant_analyze_counts(args: argparse.Namespace) -> int:
    command_args: list[PathLike] = [args.processed_dir, "--central-window-size", str(args.central_window_size)]
    if args.output_dir:
        command_args.extend(["--output-dir", args.output_dir])
    return _run_module("genomics.predictors.variant_transformer.analyze_variant_counts", command_args)


def cmd_convert_vcf_to_23andme(args: argparse.Namespace) -> int:
    return _run_module("genomics.converters.vcf_to_23andme", ["--config", args.config])


def cmd_snp_ancestry_run(args: argparse.Namespace) -> int:
    command_args: list[PathLike] = ["--config", args.config]
    if args.individual:
        command_args.extend(["--individual", args.individual])
    return _run_module("genomics.predictors.snp_ancestry", command_args)


def cmd_genomes_analyzer_run(args: argparse.Namespace) -> int:
    return _run_module("genomics.workflows.genomes_analyzer", ["--config", args.config])


def cmd_non_longevous_build(args: argparse.Namespace) -> int:
    return _run_module("genomics.workflows.dataset_builders.non_longevous", ["--config", args.config])


def cmd_non_longevous_build_window(args: argparse.Namespace) -> int:
    forwarded_args = list(args.forwarded_args)
    if forwarded_args[:1] == ["--"]:
        forwarded_args = forwarded_args[1:]
    return _run_module("genomics.workflows.dataset_builders.non_longevous.build_window_and_predict", forwarded_args)


def cmd_non_longevous_visualize(args: argparse.Namespace) -> int:
    return _run_module("genomics.workflows.dataset_builders.non_longevous.alphagenome_output_visualization", [args.config])


def cmd_alphagenome_analyze(args: argparse.Namespace) -> int:
    forwarded_args = list(args.forwarded_args)
    if forwarded_args[:1] == ["--"]:
        forwarded_args = forwarded_args[1:]
    return _run_module("genomics.workflows.alphagenome", forwarded_args)


def cmd_alphagenome_integrate(args: argparse.Namespace) -> int:
    forwarded_args = list(args.forwarded_args)
    if forwarded_args[:1] == ["--"]:
        forwarded_args = forwarded_args[1:]
    return _run_module("genomics.workflows.alphagenome.neural_integration", forwarded_args)


def cmd_alphagenome_tracks(args: argparse.Namespace) -> int:
    command_args: list[PathLike] = []
    if args.api_key:
        command_args.extend(["--api-key", args.api_key])
    if args.output:
        command_args.extend(["--output", args.output])
    return _run_module("genomics.workflows.alphagenome.alphagenome_output_tracks", command_args)


def cmd_neural_train(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    return _run_module("legacy.neural_ancestry_predictor_deprecated.neural_ancestry_predictor_deprecated", ["--config", config_path, "--mode", "train"])


def cmd_neural_test(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    return _run_module("legacy.neural_ancestry_predictor_deprecated.neural_ancestry_predictor_deprecated", ["--config", config_path, "--mode", "test"])


def cmd_neural_summarize(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    return _run_module(
        "legacy.neural_ancestry_predictor_deprecated.neural_ancestry_predictor_deprecated",
        ["--config", config_path, "--summarize_results", "--sort_by", args.sort_by],
    )


def cmd_neural_pca_cache(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    command_args: list[PathLike] = ["--config", config_path]
    if args.force:
        command_args.append("--force")
    return _run_module("legacy.neural_ancestry_predictor_deprecated.sklearn_pca_cache", command_args)


def cmd_completion_bash(args: argparse.Namespace) -> int:
    print(resources.read_text("genomics.completions", "genomics.bash"))
    return 0


def cmd_audit_configs(args: argparse.Namespace) -> int:
    roots = [
        Path("configs/predictors/genotype_based"),
        Path("configs/predictors/variant_transformer"),
    ]
    rows = []
    for root in roots:
        for path in sorted(root.glob("*.yaml")):
            text = path.read_text(encoding="utf-8")
            try:
                payload = _load_yaml(path)
            except Exception:
                payload = {}
            dataset_section = payload.get("dataset_input") if isinstance(payload.get("dataset_input"), dict) else payload.get("dataset")
            dataset_section = dataset_section if isinstance(dataset_section, dict) else {}
            metadata_section = payload.get("metadata") if isinstance(payload.get("metadata"), dict) else {}
            active = metadata_section.get("active", True) is not False
            data_paths = []
            for key in ("dataset_dir", "source_dataset_dir", "processed_dir", "consensus_dataset_dir", "reference_predictions_dataset_dir"):
                if dataset_section.get(key):
                    data_paths.append((key, str(dataset_section[key]), classify_dataset_path(dataset_section[key])))
            rows.append(
                {
                    "path": str(path),
                    "dataset_id": dataset_section.get("dataset_id") or dataset_section.get("source_dataset_id"),
                    "active": active,
                    "data_paths": data_paths,
                    "uses_dados": "/dados" in text,
                    "uses_repo_results": "results/" in text or "results_dir" in text,
                    "has_results_dir": "results_dir" in text,
                    "has_cache_dir": "cache_dir" in text,
                    "has_processed_cache_dir": "processed_cache_dir" in text,
                    "has_processed_dir": "processed_dir" in text,
                    "needs_results_migration": (
                        ("results_dir" not in text and ("processed_cache_dir" in text or "processed_dir" in text))
                        or "non_longevous_results_runs" in text
                    ),
                }
            )
    legacy_rows = [
        row for row in rows
        if any(classification == "legacy-top3" for _key, _value, classification in row["data_paths"])
    ]
    inactive_tokens = ("verify_", "random", "snps_only")
    active_legacy_rows = [
        row for row in legacy_rows
        if row.get("active", True) and not any(token in Path(row["path"]).name for token in inactive_tokens)
    ]
    if args.fail_on_active_legacy:
        output_rows = active_legacy_rows
    elif args.fail_on_legacy or args.legacy_only:
        output_rows = legacy_rows
    else:
        output_rows = rows
    if args.json:
        print(json.dumps(output_rows, indent=2, ensure_ascii=False))
    else:
        for row in output_rows:
            flags = []
            if not row.get("active", True):
                flags.append("inactive")
            if row["dataset_id"]:
                flags.append(f"dataset_id={row['dataset_id']}")
            for key, _value, classification in row["data_paths"]:
                flags.append(f"{key}:{classification}")
            if row["uses_dados"]:
                flags.append("dados-ok")
            if row["needs_results_migration"]:
                flags.append("migrar-resultados")
            if row["has_results_dir"]:
                flags.append("results_dir")
            if row["has_cache_dir"]:
                flags.append("cache_dir")
            if row["has_processed_cache_dir"]:
                flags.append("processed_cache_dir")
            if row["has_processed_dir"]:
                flags.append("processed_dir")
            print(f"{row['path']}: {', '.join(flags) if flags else 'sem paths explicitos'}")
    if args.fail_on_legacy and legacy_rows:
        print(f"ERRO: {len(legacy_rows)} config(s) ainda referenciam paths legacy-top3", file=sys.stderr)
        return 2
    if args.fail_on_active_legacy and active_legacy_rows:
        print(f"ERRO: {len(active_legacy_rows)} config(s) ativos ainda referenciam paths legacy-top3", file=sys.stderr)
        return 2
    return 0


BCFTOOLS_CHAIN_REQUIRED_TEMPLATES = (
    "{sample}.H1.window.raw.fa",
    "{sample}.H2.window.raw.fa",
    "{sample}.window.consensus_ready.vcf.gz",
    "{sample}.window.consensus_ready.vcf.gz.tbi",
    "{sample}.window.vcf.gz",
    "{sample}.window.vcf.gz.tbi",
)


def _load_json_file(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    return payload if isinstance(payload, dict) else {}


def _audit_dataset(dataset_id: str, *, check_bcftools_chain: bool, sample_limit: int, genes: Optional[list[str]]) -> Dict[str, Any]:
    ref = resolve_dataset(dataset_id)
    path = ref.path
    row: Dict[str, Any] = {
        "dataset_id": dataset_id,
        "role": ref.role,
        "path": str(path),
        "exists": path.exists(),
        "metadata_exists": ref.metadata_path.exists(),
        "layout_metadata_exists": (path / "layout_metadata.json").exists(),
        "individuals_dir_exists": (path / "individuals").exists(),
        "references_dir_exists": (path / "references").exists(),
        "reference_windows_dir_exists": (path / "references" / "windows").exists(),
        "status": "ok",
        "errors": [],
        "warnings": [],
    }
    if not row["exists"]:
        row["status"] = "error"
        row["errors"].append("dataset directory missing")
        return row
    if not row["metadata_exists"]:
        row["status"] = "error"
        row["errors"].append("dataset_metadata.json missing")
        return row
    if not row["individuals_dir_exists"]:
        row["status"] = "error"
        row["errors"].append("individuals directory missing")
    if not row["reference_windows_dir_exists"]:
        row["status"] = "error"
        row["errors"].append("references/windows directory missing")

    metadata = _load_json_file(ref.metadata_path)
    individuals = [str(x) for x in metadata.get("individuals", [])]
    metadata_genes = [str(x) for x in metadata.get("genes", [])]
    catalog = metadata.get("window_catalog", {}) or {}
    if not catalog:
        ref_windows = path / "references" / "windows"
        catalog = {p.parent.name: {} for p in sorted(ref_windows.glob("*/window_metadata.json"))}
    row.update(
        {
            "individual_count": len(individuals),
            "gene_count": len(metadata_genes),
            "window_catalog_count": len(catalog),
        }
    )
    if not individuals:
        row["warnings"].append("metadata has no individuals")
    if not metadata_genes:
        row["warnings"].append("metadata has no genes")
    if not catalog:
        row["warnings"].append("no window_catalog and no reference window metadata found")
    else:
        missing_ref_windows = []
        for window_name in sorted(catalog.keys()):
            ref_window_dir = path / "references" / "windows" / str(window_name)
            if not (ref_window_dir / "ref.window.fa").exists() or not (ref_window_dir / "window_metadata.json").exists():
                missing_ref_windows.append(str(window_name))
        if missing_ref_windows:
            row["status"] = "error"
            row["errors"].append("reference window artifacts missing")
            row["missing_reference_windows"] = missing_ref_windows[:50]

    if check_bcftools_chain:
        row["bcftools_chain"] = _audit_bcftools_chain_artifacts(path, individuals, metadata_genes, catalog, sample_limit, genes)
        if row["bcftools_chain"]["missing"]:
            row["status"] = "error"
            row["errors"].append("bcftools_chain artifacts missing")
    if row["status"] == "ok" and row["warnings"]:
        row["status"] = "warning"
    return row


def _audit_bcftools_chain_artifacts(
    dataset_dir: Path,
    individuals: list[str],
    metadata_genes: list[str],
    catalog: Dict[str, Any],
    sample_limit: int,
    genes: Optional[list[str]],
) -> Dict[str, Any]:
    selected_samples = individuals[:sample_limit] if sample_limit > 0 else individuals
    selected_genes = genes or metadata_genes or sorted(catalog)
    stats: Dict[str, Any] = {
        "sample_limit": sample_limit,
        "samples_checked": len(selected_samples),
        "genes_checked": len(selected_genes),
        "required": 0,
        "present": 0,
        "missing": 0,
        "examples": [],
    }
    for sample_id in selected_samples:
        for gene in selected_genes:
            window_dir = dataset_dir / "individuals" / sample_id / "windows" / gene
            for template in BCFTOOLS_CHAIN_REQUIRED_TEMPLATES:
                stats["required"] += 1
                path = window_dir / template.format(sample=sample_id)
                if path.exists():
                    stats["present"] += 1
                else:
                    stats["missing"] += 1
                    if len(stats["examples"]) < 20:
                        stats["examples"].append(str(path))
    return stats


def cmd_audit_data(args: argparse.Namespace) -> int:
    dataset_ids = args.dataset_id or list(registered_dataset_ids())
    rows = [
        _audit_dataset(
            dataset_id,
            check_bcftools_chain=args.check_bcftools_chain,
            sample_limit=args.sample_limit,
            genes=args.genes,
        )
        for dataset_id in dataset_ids
    ]
    if args.json:
        print(json.dumps(rows, indent=2, ensure_ascii=False))
    else:
        for row in rows:
            flags = [row["status"], row["role"]]
            flags.append(f"exists={row['exists']}")
            flags.append(f"metadata={row['metadata_exists']}")
            if "individual_count" in row:
                flags.append(f"individuals={row['individual_count']}")
                flags.append(f"genes={row['gene_count']}")
                flags.append(f"windows={row['window_catalog_count']}")
            if row.get("bcftools_chain"):
                bc = row["bcftools_chain"]
                flags.append(f"bcftools_missing={bc['missing']}/{bc['required']}")
            print(f"{row['dataset_id']}: {', '.join(flags)}")
            for message in row.get("errors", []):
                print(f"  ERROR: {message}")
            for message in row.get("warnings", []):
                print(f"  WARN: {message}")
            if row.get("bcftools_chain", {}).get("examples"):
                print("  Missing examples:")
                for example in row["bcftools_chain"]["examples"][:5]:
                    print(f"    {example}")
    if args.fail_on_missing and any(row["status"] == "error" for row in rows):
        return 2
    return 0


def _resolve_config_kind(kind: Optional[str], config_path: Optional[Path] = None) -> str:
    from genomics.core.config_schema import infer_schema_kind, schema_kinds

    if kind:
        return kind
    if config_path is not None:
        inferred = infer_schema_kind(config_path)
        if inferred:
            return inferred
    raise ValueError(f"Config kind is required. Choices: {', '.join(schema_kinds())}")


def cmd_config_describe(args: argparse.Namespace) -> int:
    from genomics.core.config_schema import describe_config_schema, get_schema_spec

    spec = get_schema_spec(args.kind)
    rows = describe_config_schema(args.kind)
    if args.json:
        print(json.dumps({"kind": spec.kind, "title": spec.title, "fields": rows}, indent=2, ensure_ascii=False))
        return 0
    print(f"{spec.title} config schema ({spec.kind})")
    for row in rows:
        required = "required" if row["required"] else f"default={row['default']}"
        print(f"{row['path']}")
        print(f"  type: {row['type']}")
        print(f"  {required}")
        if row["description"]:
            print(f"  {row['description']}")
    return 0


def cmd_config_schema(args: argparse.Namespace) -> int:
    from genomics.core.config_schema import json_schema_for_kind

    payload = json_schema_for_kind(args.kind)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


def cmd_config_validate(args: argparse.Namespace) -> int:
    from pydantic import ValidationError
    from genomics.core.config_schema import load_typed_config

    try:
        kind = _resolve_config_kind(args.kind, args.config)
        load_typed_config(kind, args.config)
    except ValidationError as exc:
        if args.json:
            print(json.dumps({"status": "error", "errors": exc.errors()}, indent=2, ensure_ascii=False, default=str))
        else:
            print(f"ERROR: config validation failed for {args.config}", file=sys.stderr)
            print(str(exc), file=sys.stderr)
        return 2
    except Exception as exc:
        if args.json:
            print(json.dumps({"status": "error", "error": str(exc)}, indent=2, ensure_ascii=False))
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps({"status": "ok", "kind": kind, "path": str(args.config)}, indent=2, ensure_ascii=False))
    else:
        print(f"ok: {args.config} ({kind})")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CLI comum para pipelines de pesquisa genômica")
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit = subparsers.add_parser("audit-configs", help="Lista configs e status de migracao de resultados")
    audit.add_argument("--json", action="store_true")
    audit.add_argument("--legacy-only", action="store_true")
    audit.add_argument("--fail-on-legacy", action="store_true")
    audit.add_argument("--fail-on-active-legacy", action="store_true")
    audit.set_defaults(func=cmd_audit_configs)

    audit_data = subparsers.add_parser("audit-data", help="Valida existencia fisica dos datasets registrados")
    audit_data.add_argument("--dataset-id", action="append", choices=registered_dataset_ids(), default=None)
    audit_data.add_argument("--check-bcftools-chain", action="store_true")
    audit_data.add_argument("--sample-limit", type=int, default=3)
    audit_data.add_argument("--genes", nargs="*", default=None)
    audit_data.add_argument("--json", action="store_true")
    audit_data.add_argument("--fail-on-missing", action="store_true")
    audit_data.set_defaults(func=cmd_audit_data)

    completion = subparsers.add_parser("completion", help="Gera scripts de auto-completion")
    completion_sub = completion.add_subparsers(dest="completion_shell", required=True)
    completion_bash = completion_sub.add_parser("bash")
    completion_bash.set_defaults(func=cmd_completion_bash)

    config = subparsers.add_parser("config", help="Describe and validate typed config files")
    config_sub = config.add_subparsers(dest="config_command", required=True)
    config_describe = config_sub.add_parser("describe", help="Describe fields for a config kind")
    config_describe.add_argument("kind", choices=["genotype", "variant"])
    config_describe.add_argument("--json", action="store_true")
    config_describe.set_defaults(func=cmd_config_describe)
    config_schema = config_sub.add_parser("schema", help="Print JSON schema for a config kind")
    config_schema.add_argument("kind", choices=["genotype", "variant"])
    config_schema.set_defaults(func=cmd_config_schema)
    config_validate = config_sub.add_parser("validate", help="Validate a typed YAML config")
    config_validate.add_argument("config", type=Path)
    config_validate.add_argument("--kind", choices=["genotype", "variant"], default=None)
    config_validate.add_argument("--json", action="store_true")
    config_validate.set_defaults(func=cmd_config_validate)

    convert = subparsers.add_parser("convert", help="Conversores de formatos genomicos")
    convert_sub = convert.add_subparsers(dest="convert_command", required=True)
    c23 = convert_sub.add_parser("vcf-to-23andme")
    c23.add_argument("--config", type=Path, required=True)
    c23.set_defaults(func=cmd_convert_vcf_to_23andme)

    snp = subparsers.add_parser("snp-ancestry", help="Pipeline snp_ancestry_predictor")
    snp_sub = snp.add_subparsers(dest="snp_ancestry_command", required=True)
    snp_run = snp_sub.add_parser("run")
    snp_run.add_argument("--config", type=Path, required=True)
    snp_run.add_argument("--individual", type=Path, default=None)
    snp_run.set_defaults(func=cmd_snp_ancestry_run)

    genomes_analyzer = subparsers.add_parser("genomes-analyzer", help="Pipeline operacional FASTQ/BAM/CRAM/VCF")
    genomes_analyzer_sub = genomes_analyzer.add_subparsers(dest="genomes_analyzer_command", required=True)
    ga_run = genomes_analyzer_sub.add_parser("run")
    ga_run.add_argument("--config", "-c", type=Path, required=True)
    ga_run.set_defaults(func=cmd_genomes_analyzer_run)

    builders = subparsers.add_parser("dataset-builders", help="Workflows de construcao de datasets")
    builders_sub = builders.add_subparsers(dest="dataset_builder_command", required=True)
    non_longevous = builders_sub.add_parser("non-longevous", help="Builder 1000G/AlphaGenome non-longevous")
    non_longevous_sub = non_longevous.add_subparsers(dest="non_longevous_command", required=True)
    nlb = non_longevous_sub.add_parser("build")
    nlb.add_argument("--config", type=Path, required=True)
    nlb.set_defaults(func=cmd_non_longevous_build)
    nlw = non_longevous_sub.add_parser("build-window")
    nlw.add_argument("forwarded_args", nargs=argparse.REMAINDER)
    nlw.set_defaults(func=cmd_non_longevous_build_window)
    nlv = non_longevous_sub.add_parser("visualize")
    nlv.add_argument("config", type=Path)
    nlv.set_defaults(func=cmd_non_longevous_visualize)

    alphagenome = subparsers.add_parser("alphagenome", help="Workflow AlphaGenome/neural_module")
    alphagenome_sub = alphagenome.add_subparsers(dest="alphagenome_command", required=True)
    aga = alphagenome_sub.add_parser("analyze")
    aga.add_argument("forwarded_args", nargs=argparse.REMAINDER)
    aga.set_defaults(func=cmd_alphagenome_analyze)
    agi = alphagenome_sub.add_parser("integrate")
    agi.add_argument("forwarded_args", nargs=argparse.REMAINDER)
    agi.set_defaults(func=cmd_alphagenome_integrate)
    agt = alphagenome_sub.add_parser("tracks")
    agt.add_argument("--api-key", default=None)
    agt.add_argument("--output", type=Path, default=None)
    agt.set_defaults(func=cmd_alphagenome_tracks)

    genotype = subparsers.add_parser("genotype", help="Pipeline genotype_based_predictor")
    genotype_sub = genotype.add_subparsers(dest="genotype_command", required=True)
    gp_cache = genotype_sub.add_parser("prepare-cache")
    _add_genotype_config_args(gp_cache)
    gp_cache.set_defaults(func=cmd_genotype_prepare_cache)
    gp_split = genotype_sub.add_parser("split", help="Materializa split/cache e relatórios sem treinar")
    _add_genotype_config_args(gp_split)
    gp_split.set_defaults(func=cmd_genotype_split)
    gp_train = genotype_sub.add_parser("train")
    _add_genotype_config_args(gp_train)
    gp_train.set_defaults(func=cmd_genotype_train)
    gp_eval = genotype_sub.add_parser("evaluate")
    _add_genotype_config_args(gp_eval)
    gp_eval.add_argument("--checkpoint", default="best_accuracy")
    gp_eval.add_argument("--split", choices=["train", "val", "test"], default="test")
    gp_eval.add_argument("--experiment-dir", type=Path, default=None)
    gp_eval.add_argument("--output-name", default=None)
    gp_eval.set_defaults(func=cmd_genotype_evaluate)
    gp_test = genotype_sub.add_parser("test")
    _add_genotype_config_args(gp_test)
    gp_test.add_argument("--checkpoint", default="best_accuracy")
    gp_test.add_argument("--experiment-dir", type=Path, default=None)
    gp_test.add_argument("--output-name", default=None)
    gp_test.set_defaults(func=cmd_genotype_test)
    gp_search = genotype_sub.add_parser("search")
    _add_genotype_config_args(gp_search)
    gp_search.set_defaults(func=cmd_genotype_search)
    gp_pca = genotype_sub.add_parser("pca-variance")
    _add_genotype_config_args(gp_pca)
    gp_pca.add_argument("--output", type=Path, required=True)
    gp_pca.add_argument("--json-output", type=Path, required=True)
    gp_pca.add_argument("--max-components", type=int, default=None)
    gp_pca.add_argument("--force", action="store_true")
    gp_pca.set_defaults(func=cmd_genotype_pca_variance)
    gp_workbench = genotype_sub.add_parser("workbench")
    gp_workbench.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    gp_workbench.add_argument("--runs-root", type=Path, default=DEFAULT_GENOTYPE_RUNS_ROOT)
    gp_workbench.add_argument("--consensus-dataset-dir", type=Path, default=DEFAULT_CONSENSUS_DATASET_DIR)
    gp_workbench.add_argument("--aligned-tsv-root", type=Path, default=None)
    gp_workbench.add_argument("--host", default="127.0.0.1")
    gp_workbench.add_argument("--port", type=int, default=8780)
    gp_workbench.set_defaults(func=cmd_genotype_workbench)
    gp_sync = genotype_sub.add_parser("sync-bcftools-artifacts")
    gp_sync.add_argument("--source-dir", type=Path, default=DEFAULT_CONSENSUS_DATASET_DIR)
    gp_sync.add_argument("--target-dir", type=Path, default=DEFAULT_DATASET_DIR)
    gp_sync.add_argument("--link-mode", choices=["hardlink", "symlink", "copy"], default="hardlink")
    gp_sync.add_argument("--apply", action="store_true")
    gp_sync.add_argument("--sample-limit", type=int, default=None)
    gp_sync.add_argument("--genes", nargs="*", default=None)
    gp_sync.add_argument("--json", action="store_true")
    gp_sync.set_defaults(func=cmd_genotype_sync_bcftools)
    gp_single = genotype_sub.add_parser("single-gene-screen")
    _add_genotype_config_args(gp_single)
    gp_single.add_argument("--genes-file", type=Path, default=None)
    gp_single.add_argument("--genes", nargs="*", default=None)
    gp_single.add_argument("--output-root", type=Path, default=results_path("genotype_based_predictor", "single_gene_screen"))
    gp_single.add_argument("--run-name-prefix", default="single_gene")
    gp_single.add_argument("--split-by-ontology", action="store_true")
    gp_single.add_argument("--ontology", action="append", default=None)
    gp_single.add_argument("--force", action="store_true")
    gp_single.add_argument("--dry-run", action="store_true")
    gp_single.set_defaults(func=cmd_genotype_single_gene_screen)

    variant = subparsers.add_parser("variant", help="Pipeline variant_transformer_predictor")
    variant_sub = variant.add_subparsers(dest="variant_command", required=True)
    vt_mat = variant_sub.add_parser("materialize")
    vt_mat.add_argument("--output-dir", type=Path, default=None)
    vt_mat.add_argument("--dataset-dir", type=Path, default=None)
    vt_mat.add_argument("--dataset-id", default=None)
    vt_mat.add_argument("--regions-bed", type=Path, default=None)
    vt_mat.add_argument("--genes", nargs="*", default=None)
    vt_mat.add_argument("--samples-metadata", type=Path, default=None)
    vt_mat.add_argument("--vcf-pattern", default=None)
    vt_mat.add_argument("--vcf-root-dir", default=None)
    vt_mat.add_argument("--target", default="superpopulation")
    vt_mat.add_argument("--classes", nargs="+", default=None)
    vt_mat.add_argument("--l-max", type=int, default=16)
    vt_mat.add_argument("--max-indel-size", type=int, default=50)
    vt_mat.add_argument("--sample-batch-size", type=int, default=128)
    vt_mat.add_argument("--train-split", type=float, default=0.7)
    vt_mat.add_argument("--val-split", type=float, default=0.15)
    vt_mat.add_argument("--test-split", type=float, default=0.15)
    vt_mat.add_argument("--random-seed", type=int, default=13)
    vt_mat.add_argument("--family-split-mode", choices=["family_aware", "ignore"], default="family_aware")
    vt_mat.add_argument("--max-samples", type=int, default=None)
    vt_mat.add_argument("--central-window-size", type=int, default=None)
    vt_mat.set_defaults(func=cmd_variant_materialize)
    vt_train = variant_sub.add_parser("train")
    _add_variant_config_args(vt_train)
    vt_train.set_defaults(func=cmd_variant_train)
    vt_eval = variant_sub.add_parser("evaluate")
    _add_variant_config_args(vt_eval)
    vt_eval.add_argument("--checkpoint", default="best_accuracy")
    vt_eval.add_argument("--split", choices=["train", "val", "test"], default="test")
    vt_eval.add_argument("--experiment-dir", type=Path, default=None)
    vt_eval.set_defaults(func=cmd_variant_evaluate)
    vt_counts = variant_sub.add_parser("analyze-counts")
    vt_counts.add_argument("processed_dir", type=Path)
    vt_counts.add_argument("--central-window-size", type=int, default=32768)
    vt_counts.add_argument("--output-dir", type=Path, default=None)
    vt_counts.set_defaults(func=cmd_variant_analyze_counts)

    neural = subparsers.add_parser("neural", help="Pipeline neural_ancestry_predictor_deprecated legado")
    neural_sub = neural.add_subparsers(dest="neural_command", required=True)
    na_train = neural_sub.add_parser("train")
    _add_neural_config_args(na_train)
    na_train.set_defaults(func=cmd_neural_train)
    na_test = neural_sub.add_parser("test")
    _add_neural_config_args(na_test)
    na_test.set_defaults(func=cmd_neural_test)
    na_sum = neural_sub.add_parser("summarize")
    _add_neural_config_args(na_sum)
    na_sum.add_argument("--sort-by", default="test_acc")
    na_sum.set_defaults(func=cmd_neural_summarize)
    na_pca = neural_sub.add_parser("pca-cache")
    _add_neural_config_args(na_pca)
    na_pca.add_argument("--force", action="store_true")
    na_pca.set_defaults(func=cmd_neural_pca_cache)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


def cli_main() -> int:
    return main()


if __name__ == "__main__":
    raise SystemExit(main())
