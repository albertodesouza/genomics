#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reporting helpers for neural ancestry experiments."""

from typing import Dict
import json

import numpy as np
from rich.console import Console
from rich.table import Table

from neural_ancestry_predictor_deprecated.config import get_results_dir


console = Console()


def summarize_experiments(config: Dict, sort_by: str = 'test_acc'):
    """
    Sumariza resultados de todos os experimentos e cria gráfico comparativo.

    Args:
        config: Configuração (para obter results_dir)
        sort_by: Métrica(s) para ordenação, separada por vírgula para ordenação composta
                 Ex: 'test_acc' ou 'val_acc,test_acc'
    """
    base_cache_dir = get_results_dir(config)

    if not base_cache_dir.exists():
        console.print(f"[red]Erro: Diretório base não encontrado: {base_cache_dir}[/red]")
        return

    experiments = []
    for exp_dir in base_cache_dir.iterdir():
        if exp_dir.is_dir():
            train_json = exp_dir / 'train_results.json'
            val_json = exp_dir / 'val_results.json'
            test_json = exp_dir / 'test_results.json'

            if train_json.exists() and val_json.exists() and test_json.exists():
                try:
                    with open(train_json) as f:
                        train_results = json.load(f)
                    with open(val_json) as f:
                        val_results = json.load(f)
                    with open(test_json) as f:
                        test_results = json.load(f)

                    experiments.append({
                        'name': exp_dir.name,
                        'train_accuracy': train_results.get('accuracy', 0),
                        'val_accuracy': val_results.get('accuracy', 0),
                        'test_accuracy': test_results.get('accuracy', 0)
                    })
                except Exception as e:
                    console.print(f"[yellow]⚠ Erro ao ler resultados de {exp_dir.name}: {e}[/yellow]")

    if not experiments:
        console.print("[yellow]Nenhum experimento completo encontrado![/yellow]")
        return

    console.print(f"[green]Encontrados {len(experiments)} experimentos completos[/green]")

    sort_metrics = [s.strip() for s in sort_by.split(',')]
    sort_key_map = {
        'train_acc': 'train_accuracy',
        'val_acc': 'val_accuracy',
        'test_acc': 'test_accuracy'
    }

    sort_keys = []
    for metric in sort_metrics:
        if metric not in sort_key_map:
            console.print(f"[red]Erro: Métrica inválida '{metric}'. Use: train_acc, val_acc ou test_acc[/red]")
            return
        sort_keys.append(sort_key_map[metric])

    def sort_key_func(exp):
        return tuple(-exp[key] for key in sort_keys)

    experiments = sorted(experiments, key=sort_key_func)

    if len(sort_metrics) == 1:
        console.print(f"[cyan]Ordenando por: {sort_metrics[0]} (maior para menor)[/cyan]")
    else:
        sort_desc = " → ".join(sort_metrics)
        console.print(f"[cyan]Ordenando por: {sort_desc} (prioridade da esquerda para direita, maior para menor)[/cyan]")

    import matplotlib.pyplot as plt

    exp_names = [exp['name'] for exp in experiments]
    train_accs = [exp['train_accuracy'] for exp in experiments]
    val_accs = [exp['val_accuracy'] for exp in experiments]
    test_accs = [exp['test_accuracy'] for exp in experiments]

    x = np.arange(len(exp_names))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(12, len(exp_names) * 2), 8))

    bars1 = ax.bar(x - width, train_accs, width, label='Train', color='#1f77b4')
    bars2 = ax.bar(x, val_accs, width, label='Validation', color='#ff7f0e')
    bars3 = ax.bar(x + width, test_accs, width, label='Test', color='#2ca02c')

    ax.set_xlabel('Experimento', fontsize=12, fontweight='bold')
    ax.set_ylabel('Accuracy', fontsize=12, fontweight='bold')
    ax.set_title('Comparação de Experimentos - Accuracy', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(exp_names, rotation=45, ha='right', fontsize=8)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.0)

    def add_value_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}',
                    ha='center', va='bottom', fontsize=7)

    add_value_labels(bars1)
    add_value_labels(bars2)
    add_value_labels(bars3)

    plt.tight_layout()

    output_path = base_cache_dir / 'experiments_summary.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    console.print(f"[green]✓ Gráfico salvo em: {output_path}[/green]")

    plt.show()

    table = Table(title="Resumo dos Experimentos")
    table.add_column("Experimento", style="cyan")
    table.add_column("Train Acc", justify="right", style="blue")
    table.add_column("Val Acc", justify="right", style="yellow")
    table.add_column("Test Acc", justify="right", style="green")

    for exp in experiments:
        table.add_row(
            exp['name'],
            f"{exp['train_accuracy']:.4f}",
            f"{exp['val_accuracy']:.4f}",
            f"{exp['test_accuracy']:.4f}"
        )

    console.print(table)
