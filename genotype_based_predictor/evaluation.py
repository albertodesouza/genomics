# -*- coding: utf-8 -*-
"""
evaluation.py — Tester, run_test_and_save, summarize_experiments.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch
import torch.nn as nn
from rich.console import Console
from rich.table import Table
from torch.utils.data import DataLoader

from genotype_based_predictor.config import PipelineConfig
from genomics_pipeline.metrics import classification_metrics, print_classification_metrics, save_results_json

console = Console()


class Tester:
    """
    Avalia um modelo PyTorch em um DataLoader e calcula métricas padrão.

    Parameters
    ----------
    model : nn.Module
    loader : DataLoader
    full_dataset : CachedProcessedDataset
        Dataset com ``idx_to_target``.
    config : PipelineConfig
    device : torch.device
    wandb_run : Any, optional
    description : str
        Descrição exibida nos logs (ex: ``"Test (best_accuracy)"``).
    """

    def __init__(
        self,
        model: nn.Module,
        loader: DataLoader,
        full_dataset: Any,
        config: PipelineConfig,
        device: torch.device,
        wandb_run: Optional[Any] = None,
        description: str = "Test",
    ):
        self.model = model
        self.loader = loader
        self.full_dataset = full_dataset
        self.config = config
        self.device = device
        self.wandb_run = wandb_run
        self.description = description
        self.is_classification = config.output.prediction_target != "frog_likelihood"

    def test(self) -> Dict[str, Any]:
        """
        Executa avaliação completa.

        Returns
        -------
        Dict com: accuracy, precision, recall, f1, confusion_matrix,
        classification_report, num_samples.
        """
        self.model.eval()
        all_preds: List[int] = []
        all_targets: List[int] = []

        with torch.no_grad():
            for batch in self.loader:
                if len(batch) == 3:
                    features, targets, _idx = batch
                else:
                    features, targets = batch
                features = features.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)
                outputs = self.model(features)

                if self.is_classification:
                    probs = torch.softmax(outputs, dim=1)
                    preds = probs.argmax(dim=1)
                    valid = targets >= 0
                    all_preds.extend(preds[valid].cpu().tolist())
                    all_targets.extend(targets[valid].cpu().tolist())
                else:
                    all_preds.extend(outputs.cpu().tolist())
                    all_targets.extend(targets.cpu().tolist())

        if not all_targets:
            console.print("[red]Nenhuma amostra válida para avaliação![/red]")
            return {}

        if self.is_classification:
            num_classes = self.full_dataset.get_num_classes()
            labels = list(range(num_classes))
            target_names = [self.full_dataset.idx_to_target[i] for i in labels]

            results = classification_metrics(all_targets, all_preds, target_names)
            print_classification_metrics(results, f"📊 {self.description}", console)

            if self.wandb_run:
                prefix = (
                    self.description.lower()
                    .replace(" ", "_").replace("(", "").replace(")", "")
                )
                self.wandb_run.log({
                    f"{prefix}/accuracy": results["accuracy"],
                    f"{prefix}/precision": results["precision"],
                    f"{prefix}/recall": results["recall"],
                    f"{prefix}/f1": results["f1"],
                })

            return results

        else:
            preds_arr = np.array(all_preds)
            targets_arr = np.array(all_targets)
            mse = float(np.mean((preds_arr - targets_arr) ** 2))
            mae = float(np.mean(np.abs(preds_arr - targets_arr)))
            console.print(f"[{self.description}] MSE={mse:.4f} | MAE={mae:.4f}")
            return {"mse": mse, "mae": mae, "num_samples": len(all_targets)}


def run_test_and_save(
    model: nn.Module,
    loader: DataLoader,
    full_dataset: Any,
    config: PipelineConfig,
    device: torch.device,
    split_name: str,
    experiment_dir: Path,
    wandb_run: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Avalia o modelo e salva resultados JSON em ``experiment_dir``.

    Parameters
    ----------
    split_name : str
        Nome do split (ex: ``'test'``, ``'val_best_loss'``).

    Returns
    -------
    Dict com os resultados da avaliação.
    """
    tester = Tester(model, loader, full_dataset, config, device, wandb_run, description=split_name)
    results = tester.test()

    if results:
        out_path = experiment_dir / f"{split_name}_results.json"
        save_results_json(results, out_path, console)

    return results


def summarize_experiments(config: PipelineConfig, sort_by: str = "test_acc") -> None:
    """
    Lê resultados de todos os experimentos e exibe tabela comparativa.

    Parameters
    ----------
    config : PipelineConfig
        Configuração usada para descobrir o ``processed_cache_dir``.
    sort_by : str
        Métrica(s) para ordenação, simples (``'test_acc'``) ou
        composta separada por vírgula (``'val_acc,test_acc'``).
    """
    base_dir = Path(config.dataset_input.processed_cache_dir)
    if not base_dir.exists():
        console.print(f"[red]Diretório não encontrado: {base_dir}[/red]")
        return

    rows: List[Dict[str, Any]] = []
    for exp_dir in sorted(base_dir.iterdir()):
        if not exp_dir.is_dir():
            continue
        row: Dict[str, Any] = {"experiment": exp_dir.name}
        for prefix in ("test", "val", "train"):
            for checkpoint in ("", "_best_loss"):
                fname = exp_dir / f"{prefix}{checkpoint}_results.json"
                if fname.exists():
                    try:
                        with open(fname) as f:
                            data = json.load(f)
                        row[f"{prefix}{checkpoint}_acc"] = data.get("accuracy")
                        row[f"{prefix}{checkpoint}_f1"] = data.get("f1")
                    except Exception:
                        pass
        if len(row) > 1:
            rows.append(row)

    if not rows:
        console.print("[yellow]Nenhum resultado encontrado.[/yellow]")
        return

    sort_keys = [s.strip() for s in sort_by.split(",")]

    def sort_key(r: Dict[str, Any]):
        return tuple(-(r.get(k) or 0.0) for k in sort_keys)

    rows.sort(key=sort_key)

    table = Table(title="📊 Comparação de Experimentos", show_header=True)
    table.add_column("Experimento", style="cyan")
    for col in ("test_acc", "val_acc", "train_acc", "test_f1", "val_f1"):
        table.add_column(col, justify="right")

    for r in rows:
        table.add_row(
            r["experiment"],
            f"{r.get('test_acc', 0.0):.4f}" if r.get("test_acc") is not None else "—",
            f"{r.get('val_acc', 0.0):.4f}" if r.get("val_acc") is not None else "—",
            f"{r.get('train_acc', 0.0):.4f}" if r.get("train_acc") is not None else "—",
            f"{r.get('test_f1', 0.0):.4f}" if r.get("test_f1") is not None else "—",
            f"{r.get('val_f1', 0.0):.4f}" if r.get("val_f1") is not None else "—",
        )

    console.print(table)
    console.print(f"[green]Total: {len(rows)} experimentos[/green]")
