#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para gerar grafico de barras duplas com Accuracy e Parameter Count
para diferentes experimentos de CNN.
"""

import matplotlib.pyplot as plt
import numpy as np

# Dados da tabela de experimentos
data = {
    'Symbol': [
        'K11',
        'K12',
        'K13',
        'K21',
        'K22',
        'K23',
        'FC128',
        'FC256',
        'FC512',
        'W8',
        'W32',
        'W128',
        'W512',
        'AG',
        '9G',
        '7G',
    ],
    'Accuracy': [
        0.78,  # K11
        0.80,  # K12
        0.83,  # K13
        0.81,  # K21
        0.80,  # K22
        0.81,  # K23
        0.78,  # FC128
        0.80,  # FC256
        0.80,  # FC512
        0.79,  # W8
        0.80,  # W32
        0.80,  # W128
        0.85,  # W512
        0.80,  # AG
        0.80,  # 9G
        0.80,  # 7G
    ],
    'Parameter Count': [
        376437,  # K11
        197749,  # K12
        200821,  # K13
        192629,  # K21
        197749,  # K22
        213109,  # K23
        106869,  # FC128
        197749,  # FC256
        379509,  # FC512
        197749,  # W8
        197749,  # W32
        197749,  # W128
        197749,  # W512
        197749,  # AG
        197749,  # 9G
        197749,  # 7G
    ],
}

def create_dual_bar_chart():
    """Cria grafico de barras duplas com dois eixos Y."""
    
    # Configuracao do grafico
    fig, ax1 = plt.subplots(figsize=(14, 7))
    
    # Posicoes das barras no eixo X
    x = np.arange(len(data['Symbol']))
    width = 0.35  # Largura das barras
    
    # Eixo Y esquerdo - Accuracy
    color_accuracy = '#2E86AB'  # Azul
    ax1.set_xlabel('Symbol', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Accuracy', color=color_accuracy, fontsize=12, fontweight='bold')
    bars1 = ax1.bar(x - width/2, data['Accuracy'], width, 
                    label='Accuracy', color=color_accuracy, alpha=0.8)
    ax1.tick_params(axis='y', labelcolor=color_accuracy)
    ax1.set_ylim(0.0, 1.0)
    ax1.set_yticks(np.arange(0.0, 1.1, 0.1))
    
    # Eixo Y direito - Parameter Count
    ax2 = ax1.twinx()
    color_params = '#E94F37'  # Vermelho
    ax2.set_ylabel('Parameter Count', color=color_params, fontsize=12, fontweight='bold')
    bars2 = ax2.bar(x + width/2, data['Parameter Count'], width,
                    label='Parameter Count', color=color_params, alpha=0.8)
    ax2.tick_params(axis='y', labelcolor=color_params)
    ax2.set_ylim(0, 400000)
    ax2.set_yticks(np.arange(0, 450000, 50000))
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
    
    # Configuracao do eixo X
    ax1.set_xticks(x)
    ax1.set_xticklabels(data['Symbol'], rotation=45, ha='right', fontsize=10)
    
    # Titulo
    plt.title('Comparacao de Accuracy e Parameter Count por Experimento', 
              fontsize=14, fontweight='bold', pad=20)
    
    # Legenda combinada
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
    
    # Grid apenas no eixo Y esquerdo (accuracy)
    ax1.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax1.set_axisbelow(True)
    
    # Ajuste de layout
    plt.tight_layout()
    
    # Salvar figura
    output_path = '/home/lume2/genomics/neural_ancestry_predictor/experiment_comparison.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f'Grafico salvo em: {output_path}')
    
    plt.show()
    
    return output_path

if __name__ == '__main__':
    create_dual_bar_chart()
