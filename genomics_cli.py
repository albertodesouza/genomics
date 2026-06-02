from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

import yaml

from genomics_pipeline.data_registry import classify_dataset_path, resolve_dataset
from genomics_workspace import (
    DEFAULT_CONSENSUS_DATASET_DIR,
    DEFAULT_DATASET_DIR,
    DEFAULT_GENOTYPE_RUNS_ROOT,
    DEFAULT_VARIANT_TRANSFORMER_DATASET_DIR,
    results_path,
)


PathLike = Union[str, Path]


def _run_module(module: str, args: Iterable[PathLike]) -> int:
    command = [sys.executable, "-m", module, *[str(arg) for arg in args]]
    return subprocess.call(command)


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
    from genotype_based_predictor.config import get_dataset_cache_dir, get_experiment_runs_dir, load_config
    from genotype_based_predictor.data_pipeline import prepare_data

    config_path = _genotype_config_with_overrides(args)
    config = load_config(config_path)
    run_dir = get_experiment_runs_dir(config).resolve() / "_prepare_cache"
    run_dir.mkdir(parents=True, exist_ok=True)
    prepare_data(config, run_dir)
    print(f"Dataset cache: {get_dataset_cache_dir(config)}")
    return 0


def cmd_genotype_train(args: argparse.Namespace) -> int:
    return _run_module("genotype_based_predictor.train", [_genotype_config_with_overrides(args)])


def cmd_genotype_evaluate(args: argparse.Namespace) -> int:
    config_path = _genotype_config_with_overrides(args)
    command_args: list[PathLike] = [config_path, "--checkpoint", args.checkpoint, "--split", args.split]
    if args.experiment_dir:
        command_args.extend(["--experiment-dir", args.experiment_dir])
    if args.output_name:
        command_args.extend(["--output-name", args.output_name])
    return _run_module("genotype_based_predictor.evaluate_checkpoint", command_args)


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
    return _run_module("genotype_based_predictor.plot_sklearn_pca_variance", command_args)


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
    return _run_module("genotype_based_predictor.genomics_workbench", command_args)


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
    return _run_module("genotype_based_predictor.sync_bcftools_chain_artifacts", command_args)


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
    return _run_module("genotype_based_predictor.run_single_gene_screen", command_args)


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
    return _run_module("variant_transformer_predictor.materialize_dataset", command_args)


def cmd_variant_train(args: argparse.Namespace) -> int:
    return _run_module("variant_transformer_predictor.train", [_variant_config_with_overrides(args)])


def cmd_variant_evaluate(args: argparse.Namespace) -> int:
    config_path = _variant_config_with_overrides(args)
    command_args: list[PathLike] = [config_path, "--checkpoint", args.checkpoint, "--split", args.split]
    if args.experiment_dir:
        command_args.extend(["--experiment-dir", args.experiment_dir])
    return _run_module("variant_transformer_predictor.evaluate_checkpoint", command_args)


def cmd_variant_analyze_counts(args: argparse.Namespace) -> int:
    command_args: list[PathLike] = [args.processed_dir, "--central-window-size", str(args.central_window_size)]
    if args.output_dir:
        command_args.extend(["--output-dir", args.output_dir])
    return _run_module("variant_transformer_predictor.analyze_variant_counts", command_args)


def cmd_neural_train(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    return _run_module("neural_ancestry_predictor_deprecated.neural_ancestry_predictor_deprecated", ["--config", config_path, "--mode", "train"])


def cmd_neural_test(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    return _run_module("neural_ancestry_predictor_deprecated.neural_ancestry_predictor_deprecated", ["--config", config_path, "--mode", "test"])


def cmd_neural_summarize(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    return _run_module(
        "neural_ancestry_predictor_deprecated.neural_ancestry_predictor_deprecated",
        ["--config", config_path, "--summarize_results", "--sort_by", args.sort_by],
    )


def cmd_neural_pca_cache(args: argparse.Namespace) -> int:
    config_path = _neural_config_with_overrides(args)
    command_args: list[PathLike] = ["--config", config_path]
    if args.force:
        command_args.append("--force")
    return _run_module("neural_ancestry_predictor_deprecated.sklearn_pca_cache", command_args)


def cmd_audit_configs(args: argparse.Namespace) -> int:
    roots = [
        Path("genotype_based_predictor/configs"),
        Path("neural_ancestry_predictor_deprecated/configs"),
        Path("variant_transformer_predictor/configs"),
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CLI comum para pipelines de pesquisa genômica")
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit = subparsers.add_parser("audit-configs", help="Lista configs e status de migracao de resultados")
    audit.add_argument("--json", action="store_true")
    audit.add_argument("--legacy-only", action="store_true")
    audit.add_argument("--fail-on-legacy", action="store_true")
    audit.add_argument("--fail-on-active-legacy", action="store_true")
    audit.set_defaults(func=cmd_audit_configs)

    genotype = subparsers.add_parser("genotype", help="Pipeline genotype_based_predictor")
    genotype_sub = genotype.add_subparsers(dest="genotype_command", required=True)
    gp_cache = genotype_sub.add_parser("prepare-cache")
    _add_genotype_config_args(gp_cache)
    gp_cache.set_defaults(func=cmd_genotype_prepare_cache)
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


if __name__ == "__main__":
    raise SystemExit(main())
