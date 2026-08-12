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
  • FC: 704 → 256 → N (N=5 ancestralidade ou N=2 pigmentação binária)

Genes: MC1R, TYRP1, TYR, SLC45A2, DDB1, EDAR, MFSD12, OCA2, HERC2, SLC24A5, TCHH

Uso:
  python draw_convnet.py                # 5 classes de ancestralidade (default)
  python draw_convnet.py --binary       # 2 classes: Strong P. / Weak P.

Based on original code by Gavin Weiguang Ding (BSD-3-Clause License)
"""

import argparse
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
import yaml

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


DEFAULT_ARCH = {
    'title_suffix': '11 Genes → 5 Ancestry Classes',
    'input_rows': 66,
    'input_width': 32768,
    'stage1_filters': 16,
    'stage1_kernel': [6, 32],
    'stage1_stride': [6, 32],
    'stage2_filters': 32,
    'stage3_filters': 64,
    'stages23_kernel': [1, 5],
    'stages23_stride': [1, 2],
    'stages23_padding': [0, 2],
    'pool_type': 'max',
    'fc_hidden_size': 256,
    'genes': [
        'MC1R', 'TYRP1', 'TYR', 'SLC45A2', 'DDB1', 'EDAR',
        'MFSD12', 'OCA2', 'HERC2', 'SLC24A5', 'TCHH',
    ],
    'output_classes': ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
}


def _conv_out(size, kernel, stride, padding=0):
    return (size + 2 * padding - kernel) // stride + 1


def _pair(value, default):
    if value is None:
        return list(default)
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return [int(value[0]), int(value[1])]
    raise ValueError(f'Esperado par [altura, largura], recebido: {value!r}')


def _load_arch_from_config(config_path, binary=False, rows_per_gene=6):
    with open(config_path, 'r', encoding='utf-8') as handle:
        config = yaml.safe_load(handle)

    dataset_input = config.get('dataset_input', {})
    model = config.get('model', {})
    cnn2 = model.get('cnn2', {})
    output = config.get('output', {})

    genes = dataset_input.get('genes_to_use') or []
    if not genes:
        genes = DEFAULT_ARCH['genes']
    input_rows = len(genes) * int(rows_per_gene)
    input_width = int(dataset_input.get('window_center_size', DEFAULT_ARCH['input_width']))
    downsample_factor = int(dataset_input.get('downsample_factor', 1) or 1)
    input_width = input_width // downsample_factor

    known_classes = output.get('known_classes')
    if binary:
        output_classes = ['Strong P.', 'Weak P.']
    elif known_classes:
        output_classes = [str(item) for item in known_classes]
    elif output.get('prediction_target') == 'superpopulation':
        output_classes = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    else:
        output_classes = DEFAULT_ARCH['output_classes']

    title_target = 'Binary Pigmentation' if binary else f'{len(output_classes)} Classes'
    return {
        'title_suffix': f"{len(genes)} Gene{'s' if len(genes) != 1 else ''} → {title_target}",
        'input_rows': input_rows,
        'input_width': input_width,
        'stage1_filters': int(cnn2.get('num_filters_stage1', DEFAULT_ARCH['stage1_filters'])),
        'stage1_kernel': _pair(cnn2.get('kernel_stage1'), DEFAULT_ARCH['stage1_kernel']),
        'stage1_stride': _pair(cnn2.get('stride_stage1'), DEFAULT_ARCH['stage1_stride']),
        'stage2_filters': int(cnn2.get('num_filters_stage2', DEFAULT_ARCH['stage2_filters'])),
        'stage3_filters': int(cnn2.get('num_filters_stage3', DEFAULT_ARCH['stage3_filters'])),
        'stages23_kernel': _pair(cnn2.get('kernel_stages23'), DEFAULT_ARCH['stages23_kernel']),
        'stages23_stride': _pair(cnn2.get('stride_stages23'), DEFAULT_ARCH['stages23_stride']),
        'stages23_padding': _pair(cnn2.get('padding_stages23'), DEFAULT_ARCH['stages23_padding']),
        'pool_type': str(cnn2.get('global_pool_type', DEFAULT_ARCH['pool_type'])),
        'fc_hidden_size': int(cnn2.get('fc_hidden_size', DEFAULT_ARCH['fc_hidden_size'])),
        'genes': genes,
        'output_classes': output_classes,
    }


def _with_shapes(arch):
    h0 = arch['input_rows']
    w0 = arch['input_width']
    k1h, k1w = arch['stage1_kernel']
    s1h, s1w = arch['stage1_stride']
    k23h, k23w = arch['stages23_kernel']
    s23h, s23w = arch['stages23_stride']
    p23h, p23w = arch['stages23_padding']

    h1 = _conv_out(h0, k1h, s1h)
    w1 = _conv_out(w0, k1w, s1w)
    h2 = _conv_out(h1, k23h, s23h, p23h)
    w2 = _conv_out(w1, k23w, s23w, p23w)
    h3 = _conv_out(h2, k23h, s23h, p23h)
    w3 = _conv_out(w2, k23w, s23w, p23w)
    if min(h1, w1, h2, w2, h3, w3) < 1:
        raise ValueError(f'Dimensao convolucional invalida: {(h1, w1, h2, w2, h3, w3)}')

    arch = dict(arch)
    arch.update({
        'stage1_shape': (arch['stage1_filters'], h1, w1),
        'stage2_shape': (arch['stage2_filters'], h2, w2),
        'stage3_shape': (arch['stage3_filters'], h3, w3),
        'pool_shape': (arch['stage3_filters'], h3, 1),
        'flatten_size': arch['stage3_filters'] * h3,
    })
    return arch


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
    depth_offset = 3.0  # Aumentado para maior visibilidade
    
    rgb = hex_to_rgb(color)
    
    # Desenhar as "fatias" de trás para frente
    for i in range(depth_visual - 1, -1, -1):
        offset = i * depth_offset
        shade = 0.7 + 0.3 * (i / max(depth_visual - 1, 1))  # Gradiente de cor
        
        rect = FancyBboxPatch(
            (x_pos + offset, -height/2 + offset * 0.3),
            width, height,
            boxstyle="round,pad=0.02,rounding_size=1.0",
            facecolor=tuple(c * shade for c in rgb),
            edgecolor='black',
            linewidth=1.0,
            alpha=alpha
        )
        ax.add_patch(rect)
    
    # Adicionar label
    if label_text:
        if label_pos == 'top':
            ax.text(x_pos + width/2 + (depth_visual-1) * depth_offset / 2, 
                   height/2 + (depth_visual-1) * depth_offset * 0.3 + 10,
                   label_text, ha='center', va='bottom', fontsize=24,
                   fontweight='bold', wrap=True)
        else:
            ax.text(x_pos + width/2 + (depth_visual-1) * depth_offset / 2, 
                   -height/2 - 10,
                   label_text, ha='center', va='top', fontsize=24,
                   fontweight='bold', wrap=True)
    
    return x_pos + width + depth_visual * depth_offset


def draw_arrow(ax, x_start, x_end, y=0, label_text=None, color='gray'):
    """Desenha uma seta entre camadas com label opcional embaixo."""
    ax.annotate('', xy=(x_end, y), xytext=(x_start, y),
                arrowprops=dict(arrowstyle='->', color=color, lw=3.0))
    
    if label_text:
        mid_x = (x_start + x_end) / 2
        ax.text(mid_x, y - 18, label_text, ha='center', va='top', 
                fontsize=18, style='italic', color='dimgray')


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
                           facecolor=rgb, edgecolor='black', linewidth=1.0)
            ax.add_patch(circle)
        
        # Pontos no meio
        for i in range(3):
            y = total_height/2 - (top_neurons + 0.5) * spacing - i * spacing * 0.4
            dot = Circle((x_pos, y), 1.5, facecolor='black')
            ax.add_patch(dot)
        
        # Neurônios de baixo
        for i in range(bottom_neurons):
            y = -total_height/2 + (bottom_neurons - 1 - i) * spacing
            circle = Circle((x_pos, y), neuron_size/2,
                           facecolor=rgb, edgecolor='black', linewidth=1.0)
            ax.add_patch(circle)
    else:
        for i in range(num_show):
            y = total_height/2 - i * spacing - spacing/2
            circle = Circle((x_pos, y), neuron_size/2,
                           facecolor=rgb, edgecolor='black', linewidth=1.0)
            ax.add_patch(circle)
    
    # Label
    if label_text:
        ax.text(x_pos, total_height/2 + 16, label_text, 
                ha='center', va='bottom', fontsize=24, fontweight='bold')
    
    return x_pos + spacing


def draw_cnn2_architecture(binary=False, arch=None):
    """
    Desenha a arquitetura completa do CNN2AncestryPredictor.
    Otimizado para colunas de artigos IEEE (fontes grandes).

    Args:
        binary: Se True, desenha a versão com 2 saídas (Strong P. / Weak P.)
                para classificação binária de pigmentação. Se False (default),
                desenha a versão com 5 classes de ancestralidade.
    """
    fig, ax = plt.subplots(1, 1, figsize=(20, 12))

    if arch is None:
        arch = dict(DEFAULT_ARCH)
        if binary:
            arch['output_classes'] = ['Strong P.', 'Weak P.']
            arch['title_suffix'] = '11 Genes → Binary Pigmentation'
    arch = _with_shapes(arch)

    output_classes = arch['output_classes']
    output_label = f"Output\n{len(output_classes)} classes"
    title_suffix = arch['title_suffix']
    num_output = len(output_classes)

    # Título
    ax.set_title('CNN2AncestryPredictor Architecture\n'
                 f'{title_suffix}',
                 fontsize=36, fontweight='bold', pad=30)
    
    x_pos = 0
    spacing = 55  # Espaçamento maior entre camadas
    
    # ========== INPUT ==========
    # 66 x 32768 (representado como retângulo horizontal - largo e baixo)
    input_height = 16   # Altura menor para ficar horizontal
    input_width = 40    # Largura maior para representar as 32768 colunas
    
    end_x = draw_layer_3d(ax, x_pos, input_height, input_width, 1,
                          COLOR_INPUT, f"Input\n{arch['input_rows']}×{arch['input_width']:,}")
    
    # Seta para Conv1
    x_pos = end_x + 10
    draw_arrow(ax, end_x, end_x + spacing - 10, 0,
               f"Conv2d\nk={tuple(arch['stage1_kernel'])}\ns={tuple(arch['stage1_stride'])}", 'darkblue')
    x_pos = end_x + spacing
    
    # ========== CONV1 + ReLU ==========
    # Output: 16 x 11 x 1024
    conv1_height = 36
    conv1_width = 12
    
    s1c, s1h, s1w = arch['stage1_shape']
    end_x = draw_layer_3d(ax, x_pos, conv1_height, conv1_width, s1c,
                          COLOR_CONV, f'Stage 1\n{s1c}×{s1h}×{s1w}\n+ ReLU')
    
    # Seta para Conv2
    draw_arrow(ax, end_x, end_x + spacing - 10, 0,
               f"Conv2d\nk={tuple(arch['stages23_kernel'])}\ns={tuple(arch['stages23_stride'])}", 'darkblue')
    x_pos = end_x + spacing
    
    # ========== CONV2 + ReLU ==========
    # Output: 32 x 11 x 512
    conv2_height = 32
    conv2_width = 10
    
    s2c, s2h, s2w = arch['stage2_shape']
    end_x = draw_layer_3d(ax, x_pos, conv2_height, conv2_width, s2c,
                          COLOR_CONV, f'Stage 2\n{s2c}×{s2h}×{s2w}\n+ ReLU')
    
    # Seta para Conv3
    draw_arrow(ax, end_x, end_x + spacing - 10, 0,
               f"Conv2d\nk={tuple(arch['stages23_kernel'])}\ns={tuple(arch['stages23_stride'])}", 'darkblue')
    x_pos = end_x + spacing
    
    # ========== CONV3 + ReLU ==========
    # Output: 64 x 11 x 256
    conv3_height = 28
    conv3_width = 8
    
    s3c, s3h, s3w = arch['stage3_shape']
    end_x = draw_layer_3d(ax, x_pos, conv3_height, conv3_width, s3c,
                          COLOR_CONV, f'Stage 3\n{s3c}×{s3h}×{s3w}\n+ ReLU')
    
    # Seta para Global Pool
    pool_name = 'AvgPool2d' if arch['pool_type'].lower() == 'avg' else 'MaxPool2d'
    draw_arrow(ax, end_x, end_x + spacing - 10, 0,
               f'{pool_name}\nglobal width', 'darkgreen')
    x_pos = end_x + spacing
    
    # ========== GLOBAL POOL ==========
    # Output: 64 x 11 x 1
    pool_height = 24
    pool_width = 6
    
    pc, ph, pw = arch['pool_shape']
    end_x = draw_layer_3d(ax, x_pos, pool_height, pool_width, pc,
                          COLOR_POOL, f'Global {arch["pool_type"].upper()} Pool\n{pc}×{ph}×{pw}')
    
    # Seta para Flatten
    draw_arrow(ax, end_x, end_x + spacing - 10, 0, 'Flatten', 'purple')
    x_pos = end_x + spacing
    
    # ========== FLATTEN / FC LAYERS ==========
    # 704 neurons
    fc_spacing = 45
    
    end_x = draw_fc_layer(ax, x_pos, arch['flatten_size'], max_show=7, color=COLOR_FC,
                          label_text=f"Flatten\n{arch['flatten_size']}", neuron_size=5)
    
    # Seta para FC
    draw_arrow(ax, x_pos + 10, x_pos + fc_spacing, 0, 'Linear', 'darkred')
    x_pos = x_pos + fc_spacing + 10
    
    # ========== FC ==========
    # 256 neurons
    end_x = draw_fc_layer(ax, x_pos, arch['fc_hidden_size'], max_show=6, color=COLOR_FC,
                          label_text=f"F.C.\n{arch['fc_hidden_size']}", neuron_size=5)
    
    # Seta para Output
    draw_arrow(ax, x_pos + 10, x_pos + fc_spacing, 0, 'Linear', 'darkred')
    x_pos = x_pos + fc_spacing + 10
    
    # ========== OUTPUT ==========
    output_neuron_size = 9
    end_x = draw_fc_layer(ax, x_pos, num_output, max_show=num_output,
                          color=COLOR_OUTPUT,
                          label_text=output_label,
                          neuron_size=output_neuron_size)

    # Adicionar labels das classes de saída (rótulos mais longos no modo binário
    # precisam de um pequeno deslocamento extra para não sobrepor os neurônios)
    spacing_output = output_neuron_size * 1.5  # neuron_size * 1.5 para output
    label_offset = 18 if binary else 16
    for i, cls in enumerate(output_classes):
        y = (num_output * spacing_output)/2 - i * spacing_output - spacing_output/2
        ax.text(x_pos + label_offset, y, cls, fontsize=24, va='center',
                fontweight='bold', color='darkblue')
    
    # ========== LEGENDA ==========
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_INPUT, edgecolor='black', label='Input'),
        mpatches.Patch(facecolor=COLOR_CONV, edgecolor='black', label='Conv + ReLU'),
        mpatches.Patch(facecolor=COLOR_POOL, edgecolor='black', label='MaxPool'),
        mpatches.Patch(facecolor=COLOR_FC, edgecolor='black', label='F.C.'),
        mpatches.Patch(facecolor=COLOR_OUTPUT, edgecolor='black', label='Output'),
    ]
    ax.legend(handles=legend_elements, loc='upper center', fontsize=28,
              ncol=5, framealpha=0.9, bbox_to_anchor=(0.5, 0.15))
    
    # Configurações do plot
    ax.set_xlim(-20, x_pos + 70)
    ax.set_ylim(-90, 100)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    
    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description='Desenha a arquitetura do CNN2AncestryPredictor.'
    )
    parser.add_argument(
        '--binary', '--pigmentation',
        dest='binary',
        action='store_true',
        help='Gera a versão da CNN com 2 saídas (Strong P. / Weak P.) '
             'para classificação binária de pigmentação. Sem este flag, '
              'gera a versão padrão com 5 classes de ancestralidade.'
    )
    parser.add_argument(
        '--config',
        type=Path,
        help='YAML de genotype_based_predictor usado para calcular a arquitetura CNN2.'
    )
    parser.add_argument(
        '--rows-per-gene',
        type=int,
        default=6,
        help='Numero de linhas/tracks por gene para inferir o input a partir da config (default: 6).'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Diretorio de saida das figuras. Default: diretorio deste script.'
    )
    parser.add_argument(
        '--no-show',
        action='store_true',
        help='Nao abre janela matplotlib; util em execucao remota/headless.'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    arch = None
    if args.config:
        arch = _load_arch_from_config(args.config, binary=args.binary, rows_per_gene=args.rows_per_gene)

    # Desenhar a arquitetura (5 classes ou 2 classes conforme o modo)
    fig = draw_cnn2_architecture(binary=args.binary, arch=arch)

    # Nome base dos arquivos: usa sufixo "_binary" no modo binário
    fig_dir = str(args.output_dir.resolve()) if args.output_dir else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(fig_dir, exist_ok=True)
    base_name = 'cnn2_ancestry_predictor_architecture'
    if args.config:
        base_name = f'{args.config.stem}_architecture'
    if args.binary:
        base_name += '_binary'

    fig_path = os.path.join(fig_dir, base_name + '.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight',
                pad_inches=0.2, facecolor='white')
    print(f"Figura salva em: {fig_path}")

    # Também salvar em PDF para alta qualidade
    fig_path_pdf = os.path.join(fig_dir, base_name + '.pdf')
    fig.savefig(fig_path_pdf, bbox_inches='tight', pad_inches=0.2)
    print(f"PDF salvo em: {fig_path_pdf}")

    if not args.no_show:
        plt.show()
