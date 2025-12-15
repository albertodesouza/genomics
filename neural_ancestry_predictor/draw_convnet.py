# -*- coding: utf-8 -*-
"""
CNN2AncestryPredictor Architecture Visualization

Desenha a arquitetura da rede neural CNN2AncestryPredictor para predição de ancestralidade.

Arquitetura:
  • Input shape: 66 x 32768 (1 canal) - 11 genes, 6 linhas por gene
  • Stage 1: 16 filters, kernel=(6, 32), stride=(6, 32) → Output: 16 x 11 x 1024
  • Stage 2: 32 filters, kernel=(1, 5), stride=(1, 2) → Output: 32 x 11 x 512
  • Stage 3: 64 filters, kernel=(1, 5), stride=(1, 2) → Output: 64 x 11 x 256
  • Global MaxPool: kernel=(1, 256) → 64 x 11 x 1
  • Flatten: 704
  • FC: 704 → 256 → 5

Genes: MC1R, TYRP1, TYR, SLC45A2, DDB1, EDAR, MFSD12, OCA2, HERC2, SLC24A5, TCHH

Based on original code by Gavin Weiguang Ding (BSD-3-Clause License)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.patches import Circle
import matplotlib.patches as mpatches

# Configurações visuais
NumDots = 4
NumConvMax = 8
NumFcMax = 20
White = 1.
Light = 0.7
Medium = 0.5
Dark = 0.3
Darker = 0.15
Black = 0.

# Cores personalizadas para a visualização
COLOR_INPUT = '#4ECDC4'      # Turquesa para input
COLOR_CONV = '#FF6B6B'       # Vermelho coral para todas as convoluções
COLOR_POOL = '#A8D8EA'       # Azul claro para pooling
COLOR_FC = '#DDA0DD'         # Roxo claro para FC layers
COLOR_OUTPUT = '#F7DC6F'     # Amarelo dourado para output


def hex_to_rgb(hex_color):
    """Converte cor hexadecimal para RGB normalizado (0-1)."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def draw_layer_3d(ax, x_pos, height, width, depth, color, label_text=None, 
                  label_pos='top', alpha=0.8):
    """
    Desenha uma camada 3D representando feature maps.
    
    Args:
        ax: Eixo matplotlib
        x_pos: Posição x do centro da camada
        height: Altura da camada (dimensão espacial H)
        width: Largura visual da camada
        depth: Profundidade (número de canais/filtros)
        color: Cor da camada
        label_text: Texto do label
        label_pos: Posição do label ('top' ou 'bottom')
        alpha: Transparência
    """
    # Ajustar profundidade visual
    depth_visual = min(depth, 15)  # Limitar profundidade visual
    depth_offset = 1.5
    
    rgb = hex_to_rgb(color)
    
    # Desenhar as "fatias" de trás para frente
    for i in range(depth_visual - 1, -1, -1):
        offset = i * depth_offset
        shade = 0.7 + 0.3 * (i / max(depth_visual - 1, 1))  # Gradiente de cor
        
        rect = FancyBboxPatch(
            (x_pos + offset, -height/2 + offset * 0.3),
            width, height,
            boxstyle="round,pad=0.02,rounding_size=0.5",
            facecolor=tuple(c * shade for c in rgb),
            edgecolor='black',
            linewidth=0.5,
            alpha=alpha
        )
        ax.add_patch(rect)
    
    # Adicionar label
    if label_text:
        if label_pos == 'top':
            ax.text(x_pos + width/2 + (depth_visual-1) * depth_offset / 2, 
                   height/2 + (depth_visual-1) * depth_offset * 0.3 + 5,
                   label_text, ha='center', va='bottom', fontsize=12,
                   fontweight='bold', wrap=True)
        else:
            ax.text(x_pos + width/2 + (depth_visual-1) * depth_offset / 2, 
                   -height/2 - 5,
                   label_text, ha='center', va='top', fontsize=12,
                   fontweight='bold', wrap=True)
    
    return x_pos + width + depth_visual * depth_offset


def draw_arrow(ax, x_start, x_end, y=0, label_text=None, color='gray'):
    """Desenha uma seta entre camadas com label opcional."""
    ax.annotate('', xy=(x_end, y), xytext=(x_start, y),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
    
    if label_text:
        mid_x = (x_start + x_end) / 2
        ax.text(mid_x, y + 8, label_text, ha='center', va='bottom', 
                fontsize=11, style='italic', color='darkgray')


def draw_fc_layer(ax, x_pos, num_neurons, max_show=10, color='#DDA0DD', 
                  label_text=None, neuron_size=3):
    """
    Desenha uma camada fully connected como círculos verticais.
    """
    rgb = hex_to_rgb(color)
    
    num_show = min(num_neurons, max_show)
    spacing = neuron_size * 1.5
    total_height = num_show * spacing
    
    # Se há mais neurônios do que podemos mostrar, indicar com ...
    show_dots = num_neurons > max_show
    
    if show_dots:
        # Mostrar alguns no topo, pontos no meio, alguns embaixo
        top_neurons = num_show // 2
        bottom_neurons = num_show - top_neurons - 2  # -2 para os pontos
        
        # Neurônios do topo
        for i in range(top_neurons):
            y = total_height/2 - i * spacing
            circle = Circle((x_pos, y), neuron_size/2, 
                           facecolor=rgb, edgecolor='black', linewidth=0.5)
            ax.add_patch(circle)
        
        # Pontos no meio
        for i in range(3):
            y = total_height/2 - (top_neurons + 0.5) * spacing - i * spacing * 0.4
            dot = Circle((x_pos, y), 0.8, facecolor='black')
            ax.add_patch(dot)
        
        # Neurônios de baixo
        for i in range(bottom_neurons):
            y = -total_height/2 + (bottom_neurons - 1 - i) * spacing
            circle = Circle((x_pos, y), neuron_size/2,
                           facecolor=rgb, edgecolor='black', linewidth=0.5)
            ax.add_patch(circle)
    else:
        for i in range(num_show):
            y = total_height/2 - i * spacing - spacing/2
            circle = Circle((x_pos, y), neuron_size/2,
                           facecolor=rgb, edgecolor='black', linewidth=0.5)
            ax.add_patch(circle)
    
    # Label
    if label_text:
        ax.text(x_pos, total_height/2 + 8, label_text, 
                ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    return x_pos + spacing


def draw_cnn2_architecture():
    """
    Desenha a arquitetura completa do CNN2AncestryPredictor.
    """
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    
    # Título
    ax.set_title('CNN2AncestryPredictor Architecture\n'
                 '11 Genes → 5 Ancestry Classes', 
                 fontsize=18, fontweight='bold', pad=20)
    
    x_pos = 0
    spacing = 25
    
    # ========== INPUT ==========
    # 66 x 32768 (representado como retângulo horizontal - largo e baixo)
    input_height = 8   # Altura menor para ficar horizontal
    input_width = 20   # Largura maior para representar as 32768 colunas
    
    end_x = draw_layer_3d(ax, x_pos, input_height, input_width, 1, 
                          COLOR_INPUT, 'Input\n66×32768\n(1 canal)')
    
    # Seta para Conv1
    x_pos = end_x + 5
    draw_arrow(ax, end_x, end_x + spacing - 5, 0, 
               'Conv2d\nk=(6,32)\ns=(6,32)', 'darkblue')
    x_pos = end_x + spacing
    
    # ========== CONV1 + ReLU ==========
    # Output: 16 x 11 x 1024
    conv1_height = 18
    conv1_width = 6
    
    end_x = draw_layer_3d(ax, x_pos, conv1_height, conv1_width, 16, 
                          COLOR_CONV, 'Stage 1\n16×11×1024\n+ ReLU')
    
    # Seta para Conv2
    draw_arrow(ax, end_x, end_x + spacing - 5, 0,
               'Conv2d\nk=(1,5)\ns=(1,2)', 'darkblue')
    x_pos = end_x + spacing
    
    # ========== CONV2 + ReLU ==========
    # Output: 32 x 11 x 512
    conv2_height = 16
    conv2_width = 5
    
    end_x = draw_layer_3d(ax, x_pos, conv2_height, conv2_width, 32, 
                          COLOR_CONV, 'Stage 2\n32×11×512\n+ ReLU')
    
    # Seta para Conv3
    draw_arrow(ax, end_x, end_x + spacing - 5, 0,
               'Conv2d\nk=(1,5)\ns=(1,2)', 'darkblue')
    x_pos = end_x + spacing
    
    # ========== CONV3 + ReLU ==========
    # Output: 64 x 11 x 256
    conv3_height = 14
    conv3_width = 4
    
    end_x = draw_layer_3d(ax, x_pos, conv3_height, conv3_width, 64, 
                          COLOR_CONV, 'Stage 3\n64×11×256\n+ ReLU')
    
    # Seta para Global Pool
    draw_arrow(ax, end_x, end_x + spacing - 5, 0,
               'MaxPool2d\nk=(1,256)', 'darkgreen')
    x_pos = end_x + spacing
    
    # ========== GLOBAL POOL ==========
    # Output: 64 x 11 x 1
    pool_height = 12
    pool_width = 3
    
    end_x = draw_layer_3d(ax, x_pos, pool_height, pool_width, 64, 
                          COLOR_POOL, 'Global Pool\n64×11×1')
    
    # Seta para Flatten
    draw_arrow(ax, end_x, end_x + spacing - 5, 0, 'Flatten', 'purple')
    x_pos = end_x + spacing
    
    # ========== FLATTEN / FC LAYERS ==========
    # 704 neurons
    fc_spacing = 20
    
    end_x = draw_fc_layer(ax, x_pos, 704, max_show=7, color=COLOR_FC,
                          label_text='Flatten\n704', neuron_size=2.5)
    
    # Seta para FC
    draw_arrow(ax, x_pos + 5, x_pos + fc_spacing, 0, 'Linear', 'darkred')
    x_pos = x_pos + fc_spacing + 5
    
    # ========== FC ==========
    # 256 neurons
    end_x = draw_fc_layer(ax, x_pos, 256, max_show=6, color=COLOR_FC,
                          label_text='FC\n256\n+ ReLU\n+ Dropout', neuron_size=2.5)
    
    # Seta para Output
    draw_arrow(ax, x_pos + 5, x_pos + fc_spacing, 0, 'Linear', 'darkred')
    x_pos = x_pos + fc_spacing + 5
    
    # ========== OUTPUT ==========
    # 5 classes
    output_neuron_size = 3.5
    end_x = draw_fc_layer(ax, x_pos, 5, max_show=5, color=COLOR_OUTPUT,
                          label_text='Output\n5 classes', neuron_size=output_neuron_size)
    
    # Adicionar labels das classes de ancestralidade
    ancestry_classes = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    spacing_output = output_neuron_size * 1.5  # neuron_size * 1.5 para output
    for i, cls in enumerate(ancestry_classes):
        y = (5 * spacing_output)/2 - i * spacing_output - spacing_output/2
        ax.text(x_pos + 8, y, cls, fontsize=12, va='center', fontweight='bold',
                color='darkblue')
    
    # ========== LEGENDA ==========
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_INPUT, edgecolor='black', label='Input'),
        mpatches.Patch(facecolor=COLOR_CONV, edgecolor='black', label='Conv + ReLU'),
        mpatches.Patch(facecolor=COLOR_POOL, edgecolor='black', label='MaxPool'),
        mpatches.Patch(facecolor=COLOR_FC, edgecolor='black', label='FC + ReLU + Dropout'),
        mpatches.Patch(facecolor=COLOR_OUTPUT, edgecolor='black', label='Output'),
    ]
    ax.legend(handles=legend_elements, loc='lower center', fontsize=11,
              ncol=5, framealpha=0.9)
    
    # Configurações do plot
    ax.set_xlim(-10, x_pos + 35)
    ax.set_ylim(-35, 50)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    
    return fig


if __name__ == '__main__':
    # Desenhar a arquitetura
    fig = draw_cnn2_architecture()
    
    # Salvar figura
    fig_dir = os.path.dirname(os.path.abspath(__file__))
    fig_name = 'cnn2_ancestry_predictor_architecture.png'
    fig_path = os.path.join(fig_dir, fig_name)
    
    fig.savefig(fig_path, dpi=150, bbox_inches='tight', 
                pad_inches=0.2, facecolor='white')
    print(f"Figura salva em: {fig_path}")
    
    # Também salvar em PDF para alta qualidade
    fig_path_pdf = os.path.join(fig_dir, 'cnn2_ancestry_predictor_architecture.pdf')
    fig.savefig(fig_path_pdf, bbox_inches='tight', pad_inches=0.2)
    print(f"PDF salvo em: {fig_path_pdf}")
    
    plt.show()
