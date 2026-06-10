from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any, Dict, List

import torch
from rich.console import Console
from rich.table import Table

from genomics.core import update_manifest
from genomics.core.metrics import save_results_json
from genomics.core.run_utils import select_device, training_manifest_fields
from genomics.core.wandb_utils import finish_wandb, init_wandb_if_enabled
from genomics.predictors.genotype_based.config import PipelineConfig, generate_dataset_name, load_config
from genomics.predictors.genotype_based.config import get_experiment_runs_dir, save_config
from genomics.predictors.genotype_based.experiments.evaluation import run_test_and_save
from genomics.predictors.genotype_based.experiments.experiment import setup_experiment_dir
from genomics.predictors.genotype_based.experiments.train import _build_model
from genomics.predictors.genotype_based.experiments.training import Trainer
from genomics.predictors.genotype_based.utils import set_random_seeds


console = Console()


def _search_dir(config: PipelineConfig) -> Path:
    search = config.hyperparameter_search
    if search.output_dir:
        return Path(search.output_dir)
    run_name = search.run_name or f"search_{generate_dataset_name(config)}"
    return get_experiment_runs_dir(config) / run_name


def _set_path(payload: Dict[str, Any], dotted_path: str, value: Any) -> None:
    parts = dotted_path.split(".")
    if not parts or any(not part for part in parts):
        raise ValueError(f"Override invalido: {dotted_path!r}")
    current: Dict[str, Any] = payload
    for part in parts[:-1]:
        next_value = current.get(part)
        if not isinstance(next_value, dict):
            raise ValueError(f"Override '{dotted_path}' aponta para campo inexistente ou nao-dicionario em '{part}'")
        current = next_value
    current[parts[-1]] = value


def _candidate_config(config: PipelineConfig, overrides: Dict[str, Any], candidate_results_dir: Path, symbol: str) -> PipelineConfig:
    payload = config.model_dump(mode="python")
    for dotted_path, value in overrides.items():
        _set_path(payload, dotted_path, value)
    payload["hyperparameter_search"]["enabled"] = False
    payload["dataset_input"]["results_dir"] = str(candidate_results_dir)
    payload["wandb"]["run_name"] = f"{payload['wandb'].get('run_name') or 'genotype-pytorch-search'}-{symbol}"
    return PipelineConfig.model_validate(payload)


def _safe_symbol(symbol: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in symbol.strip())
    return cleaned or "candidate"


def _load_best_checkpoint(model: torch.nn.Module, candidate_dir: Path, checkpoint_name: str, device: torch.device) -> Path:
    checkpoint_path = candidate_dir / "models" / f"{checkpoint_name}.pt"
    if not checkpoint_path.exists():
        checkpoint_path = candidate_dir / "models" / "final.pt"
    if checkpoint_path.exists():
        checkpoint = torch.load(checkpoint_path, map_location=device)
        state_dict = checkpoint.get("model_state_dict", checkpoint)
        model.load_state_dict(state_dict)
    return checkpoint_path


def _train_candidate(
    *,
    base_config: PipelineConfig,
    candidate: Any,
    search_dir: Path,
    rank_index: int,
    config_path: Path,
    device: torch.device,
) -> Dict[str, Any]:
    from genomics.predictors.genotype_based.data.pipeline import prepare_data

    symbol = _safe_symbol(candidate.symbol)
    candidate_root = search_dir / "candidates"
    candidate_root.mkdir(parents=True, exist_ok=True)
    candidate_config = _candidate_config(base_config, candidate.overrides, candidate_root / f"{rank_index:03d}_{symbol}", symbol)
    if candidate_config.data_split.random_seed is not None and candidate_config.data_split.random_seed != -1:
        set_random_seeds(candidate_config.data_split.random_seed, candidate_config.data_split.strict_determinism)

    temp_config_path = search_dir / "candidate_configs" / f"{rank_index:03d}_{symbol}.json"
    save_config(candidate_config, temp_config_path)

    experiment_dir = setup_experiment_dir(candidate_config, str(temp_config_path))
    full_ds, train_loader, val_loader, _test_loader = prepare_data(candidate_config, experiment_dir)

    wandb_run = init_wandb_if_enabled(candidate_config.wandb, candidate_config, experiment_dir, config_path, console)
    try:
        model = _build_model(candidate_config, full_ds).to(device)
        parameter_count = model.count_parameters() if hasattr(model, "count_parameters") else sum(p.numel() for p in model.parameters() if p.requires_grad)
        trainer = Trainer(
            model=model,
            train_loader=train_loader,
            val_loader=val_loader,
            config=candidate_config,
            device=device,
            experiment_dir=experiment_dir,
            wandb_run=wandb_run,
        )
        history = trainer.train()
        update_manifest(
            experiment_dir,
            status="interrupted" if history.get("interrupted") else "completed",
            **training_manifest_fields(history),
            search_symbol=symbol,
            search_description=candidate.description,
        )

        checkpoint_name = base_config.hyperparameter_search.pytorch.copy_best_checkpoint
        checkpoint_path = _load_best_checkpoint(model, experiment_dir, checkpoint_name, device)
        results = run_test_and_save(model, val_loader, full_ds, candidate_config, device, "val_best_accuracy", experiment_dir, wandb_run)
    finally:
        finish_wandb(wandb_run)

    row = {
        "candidate": f"{rank_index:03d}_{symbol}",
        "symbol": symbol,
        "description": candidate.description,
        "baseline": bool(candidate.baseline),
        "model_type": candidate_config.model.type,
        "params": candidate.overrides,
        "parameters": int(parameter_count),
        "val_accuracy": float(results.get("accuracy", 0.0)),
        "val_precision": float(results.get("precision", 0.0)),
        "val_recall": float(results.get("recall", 0.0)),
        "val_f1": float(results.get("f1", 0.0)),
        "artifact_path": str(checkpoint_path),
        "candidate_dir": str(experiment_dir),
    }
    return row


def _save_search_outputs(rows: List[Dict[str, Any]], best: Dict[str, Any], search_dir: Path, metric: str, paper_rows: List[Dict[str, Any]]) -> None:
    search_dir.mkdir(parents=True, exist_ok=True)
    with open(search_dir / "search_results.json", "w", encoding="utf-8") as f:
        json.dump({"selection_metric": metric, "best": best, "results": rows}, f, indent=2, ensure_ascii=False)
    with open(search_dir / "search_results.csv", "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "rank",
            "symbol",
            "description",
            "model_type",
            "val_accuracy",
            "val_precision",
            "val_recall",
            "val_f1",
            "parameters",
            "params",
            "artifact_path",
            "candidate_dir",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            csv_row = {key: row.get(key) for key in fieldnames}
            csv_row["rank"] = rank
            csv_row["params"] = json.dumps(row["params"], sort_keys=True)
            writer.writerow(csv_row)
    with open(search_dir / "paper_table.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Symbol", "Experiment Description", "Accuracy", "Parameters"])
        writer.writeheader()
        for row in paper_rows:
            writer.writerow(
                {
                    "Symbol": row["symbol"],
                    "Experiment Description": row["description"],
                    "Accuracy": f"{row['val_accuracy']:.4f}",
                    "Parameters": row["parameters"],
                }
            )
    save_results_json(best, search_dir / "best_summary.json", console)
    _save_search_plots(rows, search_dir / "plots")


def _save_search_plots(rows: List[Dict[str, Any]], plots_dir: Path) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError:
        console.print("[yellow]matplotlib indisponivel; plots da busca nao foram gerados.[/yellow]")
        return

    labels = [row["symbol"] for row in rows]
    accuracies = [row["val_accuracy"] for row in rows]
    parameters = [row["parameters"] for row in rows]

    fig, ax = plt.subplots(figsize=(max(10, len(rows) * 0.45), 5))
    bars = ax.bar(labels, accuracies)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Validation accuracy")
    ax.set_title("PyTorch hyperparameter search validation accuracy")
    ax.tick_params(axis="x", rotation=90)
    for bar, value in zip(bars, accuracies):
        ax.annotate(f"{value:.3f}", xy=(bar.get_x() + bar.get_width() / 2, value), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    fig.savefig(plots_dir / "validation_accuracy_by_candidate.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(10, len(rows) * 0.45), 5))
    ax.bar(labels, parameters)
    ax.set_ylabel("Trainable parameters")
    ax.set_title("PyTorch candidate parameter count")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(plots_dir / "parameters_by_candidate.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(parameters, accuracies)
    for row in rows:
        ax.annotate(row["symbol"], (row["parameters"], row["val_accuracy"]), fontsize=8)
    ax.set_xlabel("Trainable parameters")
    ax.set_ylabel("Validation accuracy")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Validation accuracy vs parameters")
    fig.tight_layout()
    fig.savefig(plots_dir / "accuracy_vs_parameters.png", dpi=160)
    plt.close(fig)


def _copy_best_artifact(best: Dict[str, Any], search_dir: Path) -> Path:
    best_dir = search_dir / "best"
    models_dir = best_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = Path(best["artifact_path"])
    if artifact_path.exists():
        shutil.copy2(artifact_path, models_dir / artifact_path.name)
    candidate_dir = Path(best["candidate_dir"])
    plots_dir = candidate_dir / "plots"
    if plots_dir.exists():
        shutil.copytree(plots_dir, best_dir / "plots", dirs_exist_ok=True)
    for name in ("val_best_accuracy_results.json", "config.yaml", "resolved_config.json"):
        src = candidate_dir / name
        if src.exists():
            shutil.copy2(src, best_dir / name)
    return best_dir


def run_search(config_path: Path) -> Path:
    config = load_config(config_path)
    if not config.hyperparameter_search.enabled:
        raise ValueError("hyperparameter_search.enabled deve ser true para genotype search")
    if not config.hyperparameter_search.pytorch.enabled:
        raise ValueError("hyperparameter_search.pytorch.enabled deve ser true para search_pytorch")
    candidates = list(config.hyperparameter_search.pytorch.candidates)
    if not candidates:
        raise ValueError("Nenhum candidato em hyperparameter_search.pytorch.candidates")

    search_dir = _search_dir(config).resolve()
    search_dir.mkdir(parents=True, exist_ok=True)
    device = select_device()
    console.print(f"[green]Device:[/green] {device}")
    console.print(f"[green]Busca PyTorch:[/green] {search_dir}")

    rows = []
    for idx, candidate in enumerate(candidates, start=1):
        console.print(f"[cyan]Treinando candidato {idx}/{len(candidates)}:[/cyan] {candidate.symbol} {candidate.overrides}")
        rows.append(
            _train_candidate(
                base_config=config,
                candidate=candidate,
                search_dir=search_dir,
                rank_index=idx,
                config_path=config_path,
                device=device,
            )
        )

    metric_key = f"val_{config.hyperparameter_search.selection_metric}"
    ranked_rows = sorted(rows, key=lambda row: row.get(metric_key, 0.0), reverse=True)
    best = ranked_rows[0]
    best_dir = _copy_best_artifact(best, search_dir)
    _save_search_outputs(ranked_rows, {**best, "best_dir": str(best_dir)}, search_dir, metric_key, rows)
    update_manifest(search_dir, status="completed", best_dir=str(best_dir), selection_metric=metric_key)

    table = Table(title="PyTorch hyperparameter search results", show_header=True)
    table.add_column("Rank", justify="right")
    table.add_column("Symbol")
    table.add_column("Val Acc", justify="right")
    table.add_column("Val F1", justify="right")
    table.add_column("Parameters", justify="right")
    for rank, row in enumerate(ranked_rows, start=1):
        table.add_row(str(rank), row["symbol"], f"{row['val_accuracy']:.4f}", f"{row['val_f1']:.4f}", f"{row['parameters']:,}")
    console.print(table)
    console.print(f"[bold green]Melhor modelo:[/bold green] {best['symbol']} -> {best_dir}")
    return search_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Named hyperparameter search/ablation for genotype PyTorch models")
    parser.add_argument("config_path", type=Path)
    args = parser.parse_args()
    run_search(args.config_path.resolve())


if __name__ == "__main__":
    main()
