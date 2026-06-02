#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
neural_ancestry_predictor_deprecated.py

Rede Neural para Predição de Ancestralidade a partir de Dados AlphaGenome
==========================================================================

Este módulo implementa uma rede neural configurável via YAML que prediz
ancestralidade (superpopulation, population ou FROG likelihood) a partir
de predições AlphaGenome armazenadas em um dataset PyTorch.

Uso:
    # Treino
    python3 neural_ancestry_predictor_deprecated.py --config configs/default.yaml

    # Teste
    python3 neural_ancestry_predictor_deprecated.py --config configs/default.yaml --mode test

Author: ChatGPT (for Alberto)
Created: 2025-11-14
"""

import argparse
import json
import os
import sys
import signal
from pathlib import Path

import torch
from rich.console import Console
from rich.panel import Panel
import matplotlib

_mpl = os.environ.get("MPLBACKEND")
if _mpl:
    matplotlib.use(_mpl)
else:
    try:
        matplotlib.use("TkAgg")
    except Exception:
            matplotlib.use("Agg")
from genomics_pipeline import set_random_seeds, update_manifest
from neural_ancestry_predictor_deprecated.config import (
    generate_experiment_name,
    get_results_dir,
    load_config,
)
from neural_ancestry_predictor_deprecated.data_pipeline import (
    prepare_data,
    validate_cache,
)
from neural_ancestry_predictor_deprecated.models import (
    CNN2AncestryPredictor,
    CNNAncestryPredictor,
    NNAncestryPredictor,
)
from neural_ancestry_predictor_deprecated.interpretability import DeepLIFT, GradCAM
from neural_ancestry_predictor_deprecated.experiment import setup_experiment_dir
from neural_ancestry_predictor_deprecated.interrupts import interrupt_state
from neural_ancestry_predictor_deprecated.evaluation import (
    Tester,
    extract_dna_sequence,
    run_test_and_save,
)
from neural_ancestry_predictor_deprecated.reporting import summarize_experiments
from neural_ancestry_predictor_deprecated.sklearn_baselines import (
    SKLEARN_ARTIFACT_FILENAME,
    SKLEARN_BASELINE_TYPES,
    _ensure_sklearn_classification_target,
    _sklearn_effective_random_seed,
    build_sklearn_classifier,
    load_sklearn_baseline_artifact,
    run_sklearn_eval_and_save,
    run_sklearn_test_mode,
    sklearn_metrics_dict,
    sklearn_predict_labels,
    train_sklearn_baseline,
)
from neural_ancestry_predictor_deprecated.training import Trainer

console = Console()


def main():
    """Função principal."""
    console.print(
        "[yellow]DEPRECATED:[/yellow] neural_ancestry_predictor_deprecated is kept for historical reproducibility. "
        "Use `python3 -m genomics_cli genotype train genotype_based_predictor/configs/neural_legacy/genes_1000_all.yaml` "
        "or another genotype_based_predictor config for new runs."
    )
    parser = argparse.ArgumentParser(
        description="Neural Ancestry Predictor - Predição de ancestralidade a partir de dados AlphaGenome"
    )
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Caminho para arquivo de configuração YAML'
    )
    parser.add_argument(
        '--mode',
        type=str,
        choices=['train', 'test'],
        help='Modo de operação (sobrescreve config)'
    )
    parser.add_argument(
        '--summarize_results',
        action='store_true',
        help='Sumariza resultados de todos os experimentos e gera gráfico comparativo'
    )
    parser.add_argument(
        '--sort_by',
        type=str,
        default='test_acc',
        help='Métrica(s) para ordenar experimentos (padrão: test_acc). '
             'Aceita ordenação simples (ex: val_acc) ou composta separada por vírgula '
             '(ex: val_acc,test_acc ordena por val_acc, depois test_acc como desempate)'
    )
    
    args = parser.parse_args()
    
    # Carregar configuração
    config = load_config(Path(args.config))
    
    # Limpeza agressiva de estado CUDA ANTES de qualquer coisa
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
        try:
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.reset_accumulated_memory_stats()
        except:
            pass  # Versões antigas do PyTorch podem não ter essas funções
        console.print(f"[green]✓ Estado CUDA limpo (pré-inicialização)[/green]")
    
    # Configurar semente randômica para reprodutibilidade (ANTES de qualquer operação)
    random_seed = config['data_split']['random_seed']
    if random_seed is not None and random_seed != -1:
        strict_determinism = config['data_split'].get('strict_determinism', True)
        set_random_seeds(random_seed, strict_determinism)
    elif random_seed == -1:
        console.print("[yellow]⚠ MODO DEBUG: random_seed=-1 (dados não serão embaralhados)[/yellow]")
        # Ainda configurar determinismo para outras operações
        strict_determinism = config['data_split'].get('strict_determinism', True)
        set_random_seeds(0, strict_determinism)  # Usar seed 0 para outras operações
    else:
        console.print("[yellow]⚠ Semente randômica não configurada - resultados não serão reprodutíveis[/yellow]")
    
    # Limpeza adicional APÓS configurar seeds
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
        console.print(f"[green]✓ Estado CUDA re-sincronizado (pós-seeds)[/green]")
    
    # Se flag --summarize_results, executar e sair
    if args.summarize_results:
        summarize_experiments(config, sort_by=args.sort_by)
        return
    
    # Sobrescrever modo se fornecido
    if args.mode:
        config['mode'] = args.mode
    
    # Configurar signal handler para CTRL+C (diferente para train vs test)
    if config['mode'] == 'train':
        def signal_handler(sig, frame):
            """Handler para capturar CTRL+C durante treinamento."""
            console.print("\n[yellow]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/yellow]")
            console.print("[yellow]⚠ CTRL+C detectado - Finalizando treino graciosamente...[/yellow]")
            console.print("[yellow]━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━[/yellow]")
            interrupt_state.interrupted = True
    else:
        def signal_handler(sig, frame):
            """Handler para capturar CTRL+C durante teste - interrompe imediatamente."""
            console.print("\n[yellow]⚠ CTRL+C detectado - Interrompendo...[/yellow]")
            sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    
    # Banner
    banner_text = (
        "[bold cyan]Neural Ancestry Predictor[/bold cyan]\n"
        f"Modo: {config['mode']}\n"
    )
    
    # Adicionar info do conjunto de teste se modo for test
    if config['mode'] == 'test':
        test_dataset = config.get('test_dataset', 'test')
        banner_text += f"Conjunto: {test_dataset}\n"
    
    banner_text += (
        f"Target: {config['output']['prediction_target']}\n"
        f"Config: {args.config}"
    )
    
    console.print(Panel.fit(banner_text, title="🧬 Genomics"))
    
    # Configurar device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    console.print(f"[green]Device: {device}[/green]")
    
    # Setup do diretório do experimento (apenas para treino)
    if config['mode'] == 'train':
        experiment_dir = setup_experiment_dir(config, args.config)
    else:
        # Para teste, precisa reconstruir o experiment_dir a partir dos parâmetros
        experiment_name = generate_experiment_name(config)
        experiment_dir = get_results_dir(config) / experiment_name
        
        if not experiment_dir.exists():
            console.print(f"[red]Erro: Experimento não encontrado: {experiment_dir}[/red]")
            console.print("[yellow]Execute o treinamento primeiro![/yellow]")
            sys.exit(1)
    
    # Limpeza CUDA antes de preparar dados
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
    
    # Preparar dados
    full_dataset, train_loader, val_loader, test_loader = prepare_data(config, experiment_dir)
    
    # Calcular dimensões de entrada
    input_shape = full_dataset.get_input_shape()
    num_classes = full_dataset.get_num_classes() if config['output']['prediction_target'] != 'frog_likelihood' else 150
    
    console.print(f"[green]Input shape: {input_shape[0]} x {input_shape[1]} (2D)[/green]")
    console.print(f"[green]Number of classes: {num_classes}[/green]")
    
    # Limpeza CUDA antes de criar modelo
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
        console.print(f"[green]✓ Estado CUDA limpo (pré-modelo)[/green]")
    
    # Criar modelo baseado no tipo configurado
    model_type = config['model'].get('type', 'NN').upper()
    use_sklearn_baseline = model_type in SKLEARN_BASELINE_TYPES
    model = None
    
    if use_sklearn_baseline:
        console.print(
            f"[cyan]Modo baseline sklearn: {model_type} "
            f"(StandardScaler + IncrementalPCA + classificador)[/cyan]"
        )
    elif model_type == 'NN':
        console.print(f"[cyan]Criando modelo: Neural Network (NN) totalmente conectada[/cyan]")
        model = NNAncestryPredictor(config, input_shape, num_classes).to(device)
    elif model_type == 'CNN':
        console.print(f"[cyan]Criando modelo: Convolutional Neural Network (CNN)[/cyan]")
        model = CNNAncestryPredictor(config, input_shape, num_classes).to(device)
    elif model_type == 'CNN2':
        console.print(f"[cyan]Criando modelo: CNN2 (Multi-layer with Global Pooling)[/cyan]")
        model = CNN2AncestryPredictor(config, input_shape, num_classes).to(device)
    else:
        raise ValueError(
            f"Tipo de modelo não suportado: {model_type}. "
            "Use 'NN', 'CNN', 'CNN2', 'SVM', 'RF' ou 'XGBOOST'."
        )

    if model is not None:
        console.print(f"[green]Modelo alocado em: {next(model.parameters()).device}[/green]")
    
    # Inicializar W&B
    wandb_run = None
    if config['wandb']['use_wandb']:
        try:
            import wandb
            
            # Se run_name não for especificado, usar o nome do experimento
            run_name = config['wandb'].get('run_name')
            if run_name is None:
                run_name = experiment_dir.name  # Nome do diretório do experimento
            
            wandb_run = wandb.init(
                project=config['wandb']['project_name'],
                name=run_name,
                config=config
            )
            console.print("[green]✓ Weights & Biases inicializado[/green]")
            console.print(f"  • Run name: {run_name}")
            if getattr(wandb_run, 'url', None):
                console.print(f"  • URL: {wandb_run.url}")
        except ImportError:
            console.print("[yellow]⚠ Weights & Biases não disponível. Instale com: pip install wandb[/yellow]")
        except Exception as e:
            console.print(f"[yellow]⚠ Erro ao inicializar W&B: {e}[/yellow]")
    
    # Modo de operação
    if config['mode'] == 'train':
        if use_sklearn_baseline:
            try:
                history = train_sklearn_baseline(
                    config,
                    model_type,
                    train_loader,
                    val_loader,
                    test_loader,
                    full_dataset,
                    experiment_dir,
                    wandb_run,
                )
            except (ValueError, ImportError) as e:
                console.print(f"[red]{e}[/red]")
                sys.exit(1)
            history_path = experiment_dir / 'models' / 'training_history.json'
            with open(history_path, 'w') as f:
                json.dump({'sklearn_baseline': True, **history}, f, indent=2)
            update_manifest(
                experiment_dir,
                status="completed",
                model_type=model_type,
                training_history="models/training_history.json",
                sklearn_baseline=True,
            )
            console.print(f"[green]✓ Histórico salvo em {history_path}[/green]")
            console.print(
                "\n[bold green]✓ Treino baseline sklearn concluído "
                "(train/val/test_results.json já gerados).[/bold green]"
            )
        else:
            # Carregar checkpoint se fornecido
            if config['checkpointing'].get('load_checkpoint'):
                checkpoint_path = Path(config['checkpointing']['load_checkpoint'])
                if checkpoint_path.exists():
                    console.print(f"[yellow]Carregando checkpoint: {checkpoint_path}[/yellow]")
                    
                    # Limpar memória GPU e carregar na CPU primeiro
                    if device.type == 'cuda':
                        torch.cuda.empty_cache()
                    
                    checkpoint = torch.load(checkpoint_path, map_location='cpu')
                    model.load_state_dict(checkpoint['model_state_dict'])
            
            # Limpeza CUDA final antes de treinar
            if device.type == 'cuda':
                torch.cuda.empty_cache()
                torch.cuda.synchronize()
                console.print(f"[green]✓ Estado CUDA limpo (pré-treino)[/green]")
            
            # Treinar
            trainer = Trainer(model, train_loader, val_loader, config, device, experiment_dir, wandb_run)
            history = trainer.train()
            
            # Salvar histórico
            history_path = experiment_dir / 'models' / 'training_history.json'
            with open(history_path, 'w') as f:
                json.dump(history, f, indent=2)
            update_manifest(
                experiment_dir,
                status="interrupted" if history.get('interrupted') else "completed",
                model_type=model_type,
                training_history="models/training_history.json",
                last_epoch=history.get('last_epoch'),
                best_val_accuracy=max(history.get('val_accuracy', [0.0]) or [0.0]),
                best_val_loss=min(history.get('val_loss', [float('inf')]) or [float('inf')]),
            )
            console.print(f"[green]✓ Histórico salvo em {history_path}[/green]")
            
            # Executar testes automáticos após o treino (executar MESMO se CTRL+C)
            console.print("\n[bold cyan]═══════════════════════════════════════════════[/bold cyan]")
            if interrupt_state.interrupted:
                console.print("[bold cyan]Executando Testes Após CTRL+C[/bold cyan]")
            else:
                console.print("[bold cyan]Executando Testes Automáticos Após Treinamento[/bold cyan]")
            console.print("[bold cyan]═══════════════════════════════════════════════[/bold cyan]\n")
            
            models_dir = experiment_dir / 'models'
            
            # Lista de checkpoints para testar (em ordem de prioridade)
            checkpoints_to_test = []
            if (models_dir / 'best_accuracy.pt').exists():
                checkpoints_to_test.append(('best_accuracy', models_dir / 'best_accuracy.pt'))
            if (models_dir / 'best_loss.pt').exists():
                checkpoints_to_test.append(('best_loss', models_dir / 'best_loss.pt'))
            if not checkpoints_to_test and (models_dir / 'final.pt').exists():
                checkpoints_to_test.append(('final', models_dir / 'final.pt'))
            
            if not checkpoints_to_test:
                console.print("[yellow]⚠ Nenhum checkpoint encontrado para teste automático[/yellow]")
            
            for checkpoint_name, checkpoint_path in checkpoints_to_test:
                console.print(f"\n[bold magenta]{'═' * 50}[/bold magenta]")
                console.print(f"[bold magenta]Testando com checkpoint: {checkpoint_name}.pt[/bold magenta]")
                console.print(f"[bold magenta]{'═' * 50}[/bold magenta]")
                
                # Limpar memória GPU antes de carregar checkpoint
                if device.type == 'cuda':
                    torch.cuda.empty_cache()
                
                # Carregar checkpoint na CPU, depois aplicar ao modelo (já na GPU)
                checkpoint = torch.load(checkpoint_path, map_location='cpu')
                model.load_state_dict(checkpoint['model_state_dict'])
                
                # Sufixo para arquivos de resultado
                suffix = f'_{checkpoint_name}' if checkpoint_name != 'best_accuracy' else ''
                
                # Teste no conjunto de treino
                console.print("\n[cyan]━━━ Testando no conjunto de TREINO ━━━[/cyan]")
                train_results = run_test_and_save(model, train_loader, full_dataset, config, device, f'train{suffix}', experiment_dir)
                
                # Teste no conjunto de validação
                console.print("\n[cyan]━━━ Testando no conjunto de VALIDAÇÃO ━━━[/cyan]")
                val_results = run_test_and_save(model, val_loader, full_dataset, config, device, f'val{suffix}', experiment_dir)
                
                # Teste no conjunto de teste
                console.print("\n[cyan]━━━ Testando no conjunto de TESTE ━━━[/cyan]")
                test_results = run_test_and_save(model, test_loader, full_dataset, config, device, f'test{suffix}', experiment_dir)
            
            console.print("\n[bold green]✓ Testes automáticos concluídos![/bold green]")
        
    elif config['mode'] == 'test':
        sklearn_artifact_path = experiment_dir / 'models' / SKLEARN_ARTIFACT_FILENAME
        if sklearn_artifact_path.exists():
            try:
                run_sklearn_test_mode(
                    config,
                    train_loader,
                    val_loader,
                    test_loader,
                    full_dataset,
                    experiment_dir,
                    wandb_run,
                )
            except (ValueError, FileNotFoundError) as e:
                console.print(f"[red]{e}[/red]")
                sys.exit(1)
        elif model_type in SKLEARN_BASELINE_TYPES:
            console.print(
                f"[red]ERRO: Esperado artefato sklearn em {sklearn_artifact_path} "
                f"(treine com model.type={model_type} primeiro).[/red]"
            )
            sys.exit(1)
        else:
            models_dir = experiment_dir / 'models'
            
            # Lista de checkpoints para testar
            checkpoints_to_test = []
            if (models_dir / 'best_accuracy.pt').exists():
                checkpoints_to_test.append(('best_accuracy', models_dir / 'best_accuracy.pt'))
            if (models_dir / 'best_loss.pt').exists():
                checkpoints_to_test.append(('best_loss', models_dir / 'best_loss.pt'))
            if not checkpoints_to_test and (models_dir / 'final.pt').exists():
                checkpoints_to_test.append(('final', models_dir / 'final.pt'))
            
            if not checkpoints_to_test:
                console.print(f"[red]ERRO: Nenhum checkpoint encontrado em: {models_dir}[/red]")
                sys.exit(1)
            
            # Selecionar conjunto de dados para teste
            test_dataset_choice = config.get('test_dataset', 'test').lower()
            
            if test_dataset_choice == 'train':
                selected_loader = train_loader
                dataset_name = "Train"
            elif test_dataset_choice == 'val':
                selected_loader = val_loader
                dataset_name = "Validation"
            else:  # 'test' is the default
                selected_loader = test_loader
                dataset_name = "Test"
            
            for checkpoint_name, checkpoint_path in checkpoints_to_test:
                console.print(f"\n[bold magenta]{'═' * 50}[/bold magenta]")
                console.print(f"[bold magenta]Testando com checkpoint: {checkpoint_name}.pt[/bold magenta]")
                console.print(f"[bold magenta]{'═' * 50}[/bold magenta]")
                console.print(f"[yellow]Carregando checkpoint: {checkpoint_path}[/yellow]")
                
                # Limpar memória GPU antes de carregar checkpoint
                if device.type == 'cuda':
                    torch.cuda.empty_cache()
                
                # Carregar checkpoint na CPU, depois aplicar ao modelo (já na GPU)
                checkpoint = torch.load(checkpoint_path, map_location='cpu')
                model.load_state_dict(checkpoint['model_state_dict'])
                
                console.print(f"[cyan]Testando no conjunto de: {dataset_name}[/cyan]")
                
                # Limpeza CUDA final antes de testar
                if device.type == 'cuda':
                    torch.cuda.empty_cache()
                    torch.cuda.synchronize()
                    console.print(f"[green]✓ Estado CUDA limpo (pré-teste)[/green]")
                
                # Testar
                tester = Tester(model, selected_loader, full_dataset, config, device, wandb_run, f"{dataset_name} ({checkpoint_name})")
                results = tester.test()
    
    # Finalizar W&B
    update_manifest(
        experiment_dir,
        status="completed" if not interrupt_state.interrupted else "interrupted",
        mode=config['mode'],
        wandb_enabled=wandb_run is not None,
    )

    if wandb_run:
        wandb_run.finish()
    
    console.print("\n[bold green]✓ Execução concluída![/bold green]")


if __name__ == '__main__':
    main()
