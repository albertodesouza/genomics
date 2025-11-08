#!/usr/bin/env python3
"""
read_alphagenome_predictions.py

Script para ler e analisar as predições do AlphaGenome salvas por 
build_non_longevous_dataset/build_window_and_predict.py

Exemplo de uso:
  python3 read_alphagenome_predictions.py \
    alphagenome/HG00096__CYP2B6/predictions_H1/atac.npz

Author: ChatGPT (for Alberto)
Date: 2025-11-04
"""

import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_and_analyze(npz_file: Path):
    """Carrega e analisa um arquivo .npz de predições."""
    
    print(f"\n{'='*70}")
    print(f"Analisando: {npz_file}")
    print(f"{'='*70}\n")
    
    # Carregar arrays
    data = np.load(npz_file)
    
    print(f"Arrays disponíveis: {data.files}\n")
    
    # Analisar cada array
    for array_name in data.files:
        array_data = data[array_name]
        
        print(f"Array: {array_name}")
        print(f"  Shape: {array_data.shape}")
        print(f"  Dtype: {array_data.dtype}")
        print(f"  Mean:   {array_data.mean():.6f}")
        print(f"  Std:    {array_data.std():.6f}")
        print(f"  Min:    {array_data.min():.6f}")
        print(f"  Max:    {array_data.max():.6f}")
        print(f"  Median: {np.median(array_data):.6f}")
        print()
    
    # Verificar metadados se existirem
    metadata_file = npz_file.parent / (npz_file.stem + "_metadata.json")
    if metadata_file.exists():
        print(f"Metadados encontrados em: {metadata_file}")
        with open(metadata_file) as f:
            metadata = json.load(f)
        print(f"Metadados:")
        for key, value in metadata.items():
            print(f"  {key}: {value}")
        print()
    
    return data


def plot_track(npz_file: Path, track_name: str = "values", 
               start: int = 0, end: int = 10000, output_file: Path = None):
    """Plota uma região específica de um track."""
    
    data = np.load(npz_file)
    
    if track_name not in data.files:
        print(f"Erro: Track '{track_name}' não encontrado.")
        print(f"Tracks disponíveis: {data.files}")
        return
    
    track_data = data[track_name]
    
    # Validar região
    if end > len(track_data):
        end = len(track_data)
    
    region_data = track_data[start:end]
    positions = np.arange(start, end)
    
    # Criar plot
    plt.figure(figsize=(15, 4))
    plt.plot(positions, region_data, linewidth=0.5, alpha=0.7)
    plt.xlabel('Posição no genoma (bp)')
    plt.ylabel('Valor de predição')
    plt.title(f'{npz_file.stem} - {track_name} (região {start:,}-{end:,})')
    plt.grid(True, alpha=0.3)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot salvo em: {output_file}")
    else:
        plt.show()


def compare_haplotypes(h1_npz: Path, h2_npz: Path, track_name: str = "values"):
    """Compara predições entre haplótipos."""
    
    print(f"\n{'='*70}")
    print(f"Comparando haplótipos")
    print(f"{'='*70}\n")
    
    h1_data = np.load(h1_npz)
    h2_data = np.load(h2_npz)
    
    if track_name not in h1_data.files or track_name not in h2_data.files:
        print(f"Erro: Track '{track_name}' não encontrado em ambos os arquivos.")
        return
    
    h1_track = h1_data[track_name]
    h2_track = h2_data[track_name]
    
    # Diferença absoluta
    diff = np.abs(h1_track - h2_track)
    
    print(f"H1 - Média: {h1_track.mean():.6f}, Std: {h1_track.std():.6f}")
    print(f"H2 - Média: {h2_track.mean():.6f}, Std: {h2_track.std():.6f}")
    print(f"\nDiferença absoluta:")
    print(f"  Média: {diff.mean():.6f}")
    print(f"  Max:   {diff.max():.6f}")
    print(f"  Posições com |diff| > 0.1: {(diff > 0.1).sum()} ({100*(diff > 0.1).sum()/len(diff):.2f}%)")
    print()
    
    # Plot comparação
    fig, axes = plt.subplots(3, 1, figsize=(15, 10))
    
    # Primeiros 10kb
    region = slice(0, 10000)
    positions = np.arange(10000)
    
    axes[0].plot(positions, h1_track[region], label='H1', alpha=0.7)
    axes[0].plot(positions, h2_track[region], label='H2', alpha=0.7)
    axes[0].set_ylabel('Predição')
    axes[0].set_title('Comparação H1 vs H2 (primeiros 10kb)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(positions, diff[region], color='red', alpha=0.7)
    axes[1].set_ylabel('|H1 - H2|')
    axes[1].set_title('Diferença absoluta')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].hist(diff, bins=100, alpha=0.7, edgecolor='black')
    axes[2].set_xlabel('|H1 - H2|')
    axes[2].set_ylabel('Frequência')
    axes[2].set_title('Distribuição das diferenças')
    axes[2].set_yscale('log')
    
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Ler e analisar predições do AlphaGenome"
    )
    parser.add_argument("npz_file", type=Path, help="Arquivo .npz com predições")
    parser.add_argument("--track", default="track_0", help="Nome do track a analisar")
    parser.add_argument("--plot", action="store_true", help="Gerar plot da região")
    parser.add_argument("--start", type=int, default=0, help="Posição inicial para plot")
    parser.add_argument("--end", type=int, default=10000, help="Posição final para plot")
    parser.add_argument("--output", type=Path, help="Arquivo de saída para plot (PNG)")
    parser.add_argument("--compare", type=Path, help="Arquivo H2 .npz para comparar com H1")
    
    args = parser.parse_args()
    
    if not args.npz_file.exists():
        print(f"Erro: Arquivo não encontrado: {args.npz_file}")
        return 1
    
    # Análise básica
    data = load_and_analyze(args.npz_file)
    
    # Plot se solicitado
    if args.plot:
        plot_track(args.npz_file, args.track, args.start, args.end, args.output)
    
    # Comparação se solicitada
    if args.compare:
        if not args.compare.exists():
            print(f"Erro: Arquivo de comparação não encontrado: {args.compare}")
            return 1
        compare_haplotypes(args.npz_file, args.compare, args.track)
    
    return 0


if __name__ == "__main__":
    exit(main())

