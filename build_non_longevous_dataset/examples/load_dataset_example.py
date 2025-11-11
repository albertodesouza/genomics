#!/usr/bin/env python3
"""
load_dataset_example.py

Exemplo completo de como carregar e usar GenomicLongevityDataset com PyTorch DataLoader.

Usage:
    python3 examples/load_dataset_example.py --dataset-dir non_longevous_results

Author: ChatGPT (for Alberto)
Created: 2025-11-11
"""

import argparse
import sys
from pathlib import Path

# Adicionar diretório pai ao path para imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import torch
    from torch.utils.data import DataLoader
    TORCH_AVAILABLE = True
except ImportError:
    print("[ERROR] PyTorch não está instalado.")
    print("Instale com: pip install torch")
    sys.exit(1)

from genomic_dataset import GenomicLongevityDataset, collate_fn_variable_windows
import numpy as np


def print_section(title: str):
    """Helper para imprimir seções formatadas."""
    print("\n" + "="*80)
    print(title)
    print("="*80)


def example_basic_loading(dataset_dir: Path):
    """
    Exemplo 1: Carregamento básico do dataset.
    """
    print_section("EXEMPLO 1: Carregamento Básico do Dataset")
    
    # Carregar dataset
    dataset = GenomicLongevityDataset(
        dataset_dir=dataset_dir,
        load_predictions=True,
        load_sequences=False,  # Não carregar sequências para este exemplo
        cache_sequences=False
    )
    
    # Imprimir informações do dataset
    dataset.print_summary()
    
    print(f"Total de samples: {len(dataset)}")
    
    return dataset


def example_single_sample(dataset: GenomicLongevityDataset):
    """
    Exemplo 2: Acessar um único sample.
    """
    print_section("EXEMPLO 2: Acessar um Único Sample")
    
    if len(dataset) == 0:
        print("[WARN] Dataset vazio, pulando exemplo.")
        return
    
    # Acessar primeiro sample
    input_data, output_data = dataset[0]
    
    print(f"\n📋 Informações do Indivíduo:")
    print(f"  • Sample ID: {output_data['sample_id']}")
    print(f"  • Sexo: {'Masculino' if output_data['sex'] == 1 else 'Feminino'}")
    print(f"  • População: {output_data['population']}")
    print(f"  • Superpopulação: {output_data['superpopulation']}")
    print(f"  • Longevity: {output_data['longevity']}")
    
    if output_data['frog_likelihood'] is not None:
        print(f"  • FROG likelihood shape: {output_data['frog_likelihood'].shape}")
        print(f"  • Top 3 populações FROG:")
        
        likelihood = output_data['frog_likelihood']
        pop_names = output_data['frog_population_names']
        
        # Ordenar por likelihood (decrescente)
        top_indices = np.argsort(likelihood)[::-1][:3]
        for i, idx in enumerate(top_indices, 1):
            print(f"      {i}. {pop_names[idx]}: {likelihood[idx]:.2e}")
    
    print(f"\n🧬 Janelas Genômicas:")
    windows = input_data['windows']
    print(f"  • Total de janelas: {len(windows)}")
    
    for window_name, window_data in windows.items():
        print(f"\n  📍 Janela: {window_name}")
        
        # Predições H1
        if 'predictions_h1' in window_data:
            pred_h1 = window_data['predictions_h1']
            print(f"    Predições H1:")
            for output_type, array in pred_h1.items():
                print(f"      • {output_type}: shape={array.shape}, "
                      f"mean={array.mean():.6f}, std={array.std():.6f}")
        
        # Predições H2
        if 'predictions_h2' in window_data:
            pred_h2 = window_data['predictions_h2']
            print(f"    Predições H2:")
            for output_type, array in pred_h2.items():
                print(f"      • {output_type}: shape={array.shape}, "
                      f"mean={array.mean():.6f}, std={array.std():.6f}")


def example_dataloader(dataset: GenomicLongevityDataset):
    """
    Exemplo 3: Usar DataLoader para iteração em batches.
    """
    print_section("EXEMPLO 3: Usar DataLoader para Batching")
    
    if len(dataset) == 0:
        print("[WARN] Dataset vazio, pulando exemplo.")
        return
    
    # Criar DataLoader
    dataloader = DataLoader(
        dataset,
        batch_size=2,  # Pequeno batch para exemplo
        shuffle=False,
        num_workers=0,  # 0 para evitar problemas de multiprocessing no exemplo
        collate_fn=collate_fn_variable_windows
    )
    
    print(f"\n📦 DataLoader configurado:")
    print(f"  • Batch size: 2")
    print(f"  • Total de batches: {len(dataloader)}")
    
    # Iterar sobre primeiro batch
    print(f"\n🔄 Iterando sobre primeiro batch:")
    
    for batch_idx, (batch_input, batch_output) in enumerate(dataloader):
        if batch_idx >= 1:  # Apenas primeiro batch para exemplo
            break
        
        print(f"\n  Batch {batch_idx + 1}:")
        print(f"  • Sample IDs: {batch_output['sample_id']}")
        print(f"  • Populações: {batch_output['population']}")
        print(f"  • Sexos: {batch_output['sex']}")
        print(f"  • Longevity tensor: {batch_output['longevity']}")
        
        # Examinar input de cada sample no batch
        for i, input_data in enumerate(batch_input):
            sample_id = batch_output['sample_id'][i]
            n_windows = len(input_data['windows'])
            print(f"  • {sample_id}: {n_windows} janelas")


def example_sample_by_id(dataset: GenomicLongevityDataset):
    """
    Exemplo 4: Acessar sample por ID.
    """
    print_section("EXEMPLO 4: Acessar Sample por ID")
    
    if len(dataset) == 0:
        print("[WARN] Dataset vazio, pulando exemplo.")
        return
    
    # Obter lista de samples disponíveis
    summary = dataset.get_summary()
    available_ids = summary['dataset_name']
    
    # Pegar primeiro ID disponível
    sample_id = dataset.individuals[0]
    
    print(f"\n🔍 Buscando sample: {sample_id}")
    
    try:
        input_data, output_data = dataset.get_sample_by_id(sample_id)
        
        print(f"✓ Sample encontrado!")
        print(f"  • População: {output_data['population']}")
        print(f"  • Janelas: {len(input_data['windows'])}")
        
    except ValueError as e:
        print(f"✗ Erro: {e}")


def example_prediction_analysis(dataset: GenomicLongevityDataset):
    """
    Exemplo 5: Análise de predições AlphaGenome.
    """
    print_section("EXEMPLO 5: Análise de Predições AlphaGenome")
    
    if len(dataset) == 0:
        print("[WARN] Dataset vazio, pulando exemplo.")
        return
    
    # Acessar primeiro sample
    input_data, output_data = dataset[0]
    
    windows = input_data['windows']
    if len(windows) == 0:
        print("[WARN] Nenhuma janela disponível.")
        return
    
    # Pegar primeira janela
    window_name = list(windows.keys())[0]
    window_data = windows[window_name]
    
    print(f"\n🧬 Analisando janela: {window_name}")
    print(f"  • Sample: {output_data['sample_id']}")
    
    if 'predictions_h1' in window_data and 'predictions_h2' in window_data:
        pred_h1 = window_data['predictions_h1']
        pred_h2 = window_data['predictions_h2']
        
        # Comparar H1 vs H2 para cada output type
        print(f"\n  📊 Comparação H1 vs H2:")
        
        for output_type in pred_h1.keys():
            if output_type in pred_h2:
                array_h1 = pred_h1[output_type]
                array_h2 = pred_h2[output_type]
                
                # Calcular diferença
                diff = np.abs(array_h1 - array_h2)
                
                print(f"\n    {output_type}:")
                print(f"      H1 mean: {array_h1.mean():.6f}")
                print(f"      H2 mean: {array_h2.mean():.6f}")
                print(f"      |H1-H2| mean: {diff.mean():.6f}")
                print(f"      |H1-H2| max: {diff.max():.6f}")
                print(f"      Positions with |diff| > 0.1: {(diff > 0.1).sum()}")


def example_custom_transform(dataset_dir: Path):
    """
    Exemplo 6: Dataset com transformação customizada.
    """
    print_section("EXEMPLO 6: Dataset com Transformação Customizada")
    
    def my_transform(input_data, output_data):
        """
        Exemplo de transformação customizada.
        Normaliza predições para média=0, std=1.
        """
        for window_name, window_data in input_data['windows'].items():
            # Normalizar predições H1
            if 'predictions_h1' in window_data:
                for output_type, array in window_data['predictions_h1'].items():
                    mean = array.mean()
                    std = array.std()
                    if std > 0:
                        window_data['predictions_h1'][output_type] = (array - mean) / std
            
            # Normalizar predições H2
            if 'predictions_h2' in window_data:
                for output_type, array in window_data['predictions_h2'].items():
                    mean = array.mean()
                    std = array.std()
                    if std > 0:
                        window_data['predictions_h2'][output_type] = (array - mean) / std
        
        return input_data, output_data
    
    # Carregar dataset com transformação
    dataset_transformed = GenomicLongevityDataset(
        dataset_dir=dataset_dir,
        load_predictions=True,
        load_sequences=False,
        transform=my_transform
    )
    
    print(f"\n✨ Dataset com transformação customizada carregado")
    print(f"  • Total de samples: {len(dataset_transformed)}")
    
    if len(dataset_transformed) > 0:
        input_data, output_data = dataset_transformed[0]
        
        # Verificar que predições foram normalizadas
        windows = input_data['windows']
        if len(windows) > 0:
            window_name = list(windows.keys())[0]
            window_data = windows[window_name]
            
            if 'predictions_h1' in window_data:
                pred_h1 = window_data['predictions_h1']
                for output_type, array in pred_h1.items():
                    print(f"\n  {output_type} (normalizado):")
                    print(f"    mean: {array.mean():.6f} (esperado ~0)")
                    print(f"    std: {array.std():.6f} (esperado ~1)")


def main():
    parser = argparse.ArgumentParser(
        description="Exemplos de uso do GenomicLongevityDataset"
    )
    parser.add_argument(
        "--dataset-dir",
        type=str,
        default="non_longevous_results",
        help="Diretório do dataset (default: non_longevous_results)"
    )
    parser.add_argument(
        "--example",
        type=int,
        choices=[1, 2, 3, 4, 5, 6],
        help="Executar apenas um exemplo específico (1-6)"
    )
    
    args = parser.parse_args()
    
    dataset_dir = Path(args.dataset_dir)
    
    # Verificar se dataset existe
    if not dataset_dir.exists():
        print(f"[ERROR] Dataset não encontrado: {dataset_dir}")
        print("[INFO] Execute build_non_longevous_dataset.py primeiro para criar o dataset.")
        sys.exit(1)
    
    print("="*80)
    print("EXEMPLOS DE USO: GenomicLongevityDataset")
    print("="*80)
    print(f"Dataset: {dataset_dir}")
    
    try:
        if args.example is None:
            # Executar todos os exemplos
            dataset = example_basic_loading(dataset_dir)
            example_single_sample(dataset)
            example_dataloader(dataset)
            example_sample_by_id(dataset)
            example_prediction_analysis(dataset)
            example_custom_transform(dataset_dir)
        else:
            # Executar exemplo específico
            if args.example == 1:
                dataset = example_basic_loading(dataset_dir)
            elif args.example == 2:
                dataset = example_basic_loading(dataset_dir)
                example_single_sample(dataset)
            elif args.example == 3:
                dataset = example_basic_loading(dataset_dir)
                example_dataloader(dataset)
            elif args.example == 4:
                dataset = example_basic_loading(dataset_dir)
                example_sample_by_id(dataset)
            elif args.example == 5:
                dataset = example_basic_loading(dataset_dir)
                example_prediction_analysis(dataset)
            elif args.example == 6:
                example_custom_transform(dataset_dir)
        
        print_section("CONCLUÍDO")
        print("✓ Todos os exemplos executados com sucesso!")
        print("\nPara mais informações, consulte:")
        print("  • docs/PYTORCH_DATASET.md")
        print("  • README.md")
        
    except Exception as e:
        print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

