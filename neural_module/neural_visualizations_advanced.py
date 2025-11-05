#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualizações Avançadas para AlphaGenome
Baseado na documentação oficial: https://www.alphagenomedocs.com/colabs/quick_start.html
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional
from rich.console import Console

console = Console()

# ====================== Visualizações Melhoradas ======================

def create_enhanced_track_visualization(
    output_data,
    seq_id: str,
    output_name: str,
    output_dir: Path,
    config: Dict,
    show_metadata: bool = True
):
    """
    Cria visualização melhorada de track com metadados.
    
    Baseado em: https://www.alphagenomedocs.com/colabs/quick_start.html
    """
    try:
        if not hasattr(output_data, 'values'):
            console.print(f"[yellow]⚠ {output_name} não tem dados plotáveis[/yellow]")
            return
        
        data_array = output_data.values
        
        # Criar figura com tamanho maior para melhor visualização
        fig, axes = plt.subplots(
            nrows=min(data_array.shape[1], 4),  # Máximo 4 subplots
            ncols=1,
            figsize=(config['plot_width'], min(config['plot_height'], 3 * data_array.shape[1])),
            squeeze=False
        )
        
        # Se tem metadados, mostrar informações de cada track
        metadata = None
        if hasattr(output_data, 'metadata'):
            metadata = output_data.metadata
        
        # Plotar cada track
        num_tracks = min(data_array.shape[1], 4)
        for i in range(num_tracks):
            ax = axes[i, 0]
            
            track_data = data_array[:, i] if len(data_array.shape) > 1 else data_array
            
            # Plot da track
            ax.fill_between(
                range(len(track_data)),
                track_data,
                alpha=0.6,
                color=plt.cm.tab10(i % 10),
                linewidth=0
            )
            ax.plot(
                track_data,
                linewidth=0.5,
                color=plt.cm.tab10(i % 10),
                alpha=0.8
            )
            
            # Título com metadados se disponível
            if metadata is not None and i < len(metadata):
                track_info = metadata.iloc[i]
                title_parts = []
                
                if 'biosample_name' in track_info:
                    title_parts.append(f"Tecido: {track_info['biosample_name']}")
                if 'strand' in track_info and track_info['strand'] != '.':
                    title_parts.append(f"Strand: {track_info['strand']}")
                if 'Assay title' in track_info:
                    title_parts.append(f"Assay: {track_info['Assay title']}")
                
                ax.set_title(' | '.join(title_parts) if title_parts else f'Track {i+1}',
                           fontsize=10, pad=5)
            else:
                ax.set_title(f'Track {i+1}', fontsize=10)
            
            ax.set_ylabel('Sinal', fontsize=9)
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # Apenas último plot tem xlabel
            if i == num_tracks - 1:
                ax.set_xlabel('Posição (bp)', fontsize=10)
            else:
                ax.set_xticklabels([])
        
        # Título geral
        fig.suptitle(
            f'{seq_id} - {output_name}',
            fontsize=14,
            fontweight='bold',
            y=0.995
        )
        
        plt.tight_layout()
        
        # Salvar
        for fmt in config['output_formats']:
            output_file = output_dir / f"{seq_id}_{output_name}_enhanced.{fmt}"
            plt.savefig(output_file, dpi=config['plot_resolution'], bbox_inches='tight')
            console.print(f"[green]  ✓ Salvo (enhanced): {output_file.name}[/green]")
        
        plt.close(fig)
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro em visualização enhanced: {e}[/yellow]")


def create_multi_output_comparison(
    outputs: Dict,
    seq_id: str,
    output_names: List[str],
    output_dir: Path,
    config: Dict
):
    """
    Cria visualização comparativa de múltiplos outputs.
    
    Inspirado em: https://www.alphagenomedocs.com/colabs/quick_start.html
    """
    try:
        # Filtrar apenas outputs válidos
        valid_outputs = {}
        for name in output_names:
            data = getattr(outputs, name.lower(), None)
            if data is not None and hasattr(data, 'values'):
                # Pegar primeira track de cada output
                values = data.values
                if len(values.shape) > 1:
                    values = values[:, 0]
                valid_outputs[name] = values
        
        if not valid_outputs:
            console.print("[yellow]⚠ Nenhum output válido para comparação[/yellow]")
            return
        
        # Criar figura com subplots empilhados
        n_outputs = len(valid_outputs)
        fig, axes = plt.subplots(
            nrows=n_outputs,
            ncols=1,
            figsize=(config['plot_width'], 2.5 * n_outputs),
            sharex=True
        )
        
        if n_outputs == 1:
            axes = [axes]
        
        colors = plt.cm.Set2(np.linspace(0, 1, n_outputs))
        
        for idx, (name, data) in enumerate(valid_outputs.items()):
            ax = axes[idx]
            
            # Plot com preenchimento
            ax.fill_between(
                range(len(data)),
                data,
                alpha=0.4,
                color=colors[idx],
                linewidth=0
            )
            ax.plot(data, linewidth=0.7, color=colors[idx], alpha=0.9)
            
            ax.set_ylabel(name, fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.2, linestyle='--')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # Mostrar estatísticas
            stats_text = f'μ={np.mean(data):.3f} σ={np.std(data):.3f}'
            ax.text(0.02, 0.95, stats_text, transform=ax.transAxes,
                   fontsize=8, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        axes[-1].set_xlabel('Posição (bp)', fontsize=11)
        
        fig.suptitle(
            f'{seq_id} - Comparação de Outputs',
            fontsize=15,
            fontweight='bold'
        )
        
        plt.tight_layout()
        
        # Salvar
        for fmt in config['output_formats']:
            output_file = output_dir / f"{seq_id}_comparison.{fmt}"
            plt.savefig(output_file, dpi=config['plot_resolution'], bbox_inches='tight')
            console.print(f"[green]  ✓ Salvo (comparison): {output_file.name}[/green]")
        
        plt.close(fig)
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro em comparação: {e}[/yellow]")


def create_heatmap_visualization(
    output_data,
    seq_id: str,
    output_name: str,
    output_dir: Path,
    config: Dict
):
    """
    Cria visualização em heatmap para múltiplas tracks.
    
    Útil para comparar múltiplos tecidos/condições simultaneamente.
    """
    try:
        if not hasattr(output_data, 'values'):
            return
        
        data_array = output_data.values
        
        # Apenas criar heatmap se tem múltiplas tracks
        if len(data_array.shape) < 2 or data_array.shape[1] < 2:
            return
        
        # Limitar número de tracks para visualização
        max_tracks = 20
        if data_array.shape[1] > max_tracks:
            data_array = data_array[:, :max_tracks]
        
        # Transpor para ter tracks nas linhas
        data_transposed = data_array.T
        
        # Criar figura
        fig, ax = plt.subplots(figsize=(config['plot_width'], max(6, data_transposed.shape[0] * 0.3)))
        
        # Heatmap
        im = ax.imshow(
            data_transposed,
            aspect='auto',
            cmap='viridis',
            interpolation='nearest'
        )
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Sinal', rotation=270, labelpad=15, fontsize=10)
        
        # Labels
        ax.set_xlabel('Posição (bp)', fontsize=11)
        ax.set_ylabel('Track', fontsize=11)
        ax.set_title(f'{seq_id} - {output_name} Heatmap', fontsize=14, fontweight='bold', pad=10)
        
        # Se tem metadados, usar como labels
        if hasattr(output_data, 'metadata'):
            metadata = output_data.metadata
            if len(metadata) == data_transposed.shape[0]:
                labels = []
                for _, row in metadata.iterrows():
                    if 'biosample_name' in row:
                        labels.append(row['biosample_name'][:20])  # Limitar tamanho
                    else:
                        labels.append(f"Track {len(labels)+1}")
                
                ax.set_yticks(range(len(labels)))
                ax.set_yticklabels(labels, fontsize=8)
        
        plt.tight_layout()
        
        # Salvar
        for fmt in config['output_formats']:
            output_file = output_dir / f"{seq_id}_{output_name}_heatmap.{fmt}"
            plt.savefig(output_file, dpi=config['plot_resolution'], bbox_inches='tight')
            console.print(f"[green]  ✓ Salvo (heatmap): {output_file.name}[/green]")
        
        plt.close(fig)
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro em heatmap: {e}[/yellow]")


def create_summary_dashboard(
    results: Dict,
    output_dir: Path,
    config: Dict
):
    """
    Cria dashboard resumindo todas as análises.
    """
    try:
        outputs = results['outputs']
        seq_id = results['sequence_id']
        requested_outputs = results['requested_outputs']
        
        # Coletar estatísticas
        stats = []
        for output_name in requested_outputs:
            data = getattr(outputs, output_name.lower(), None)
            if data is not None and hasattr(data, 'values'):
                values = data.values
                if len(values.shape) > 1:
                    values = values[:, 0]  # Primeira track
                
                stats.append({
                    'output': output_name,
                    'mean': np.mean(values),
                    'std': np.std(values),
                    'max': np.max(values),
                    'min': np.min(values),
                    'num_tracks': data.values.shape[1] if len(data.values.shape) > 1 else 1
                })
        
        if not stats:
            return
        
        # Criar dashboard
        fig = plt.figure(figsize=(config['plot_width'], 8))
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
        
        # 1. Comparação de médias
        ax1 = fig.add_subplot(gs[0, :])
        output_names = [s['output'] for s in stats]
        means = [s['mean'] for s in stats]
        stds = [s['std'] for s in stats]
        
        x = range(len(output_names))
        bars = ax1.bar(x, means, yerr=stds, capsize=5, alpha=0.7, color=plt.cm.Set3(range(len(output_names))))
        ax1.set_xticks(x)
        ax1.set_xticklabels(output_names, rotation=45, ha='right')
        ax1.set_ylabel('Sinal Médio', fontsize=10)
        ax1.set_title('Comparação de Sinais Médios por Output', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # 2. Número de tracks
        ax2 = fig.add_subplot(gs[1, 0])
        num_tracks = [s['num_tracks'] for s in stats]
        ax2.barh(output_names, num_tracks, alpha=0.7, color='steelblue')
        ax2.set_xlabel('Número de Tracks', fontsize=10)
        ax2.set_title('Tracks por Output', fontsize=11, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
        
        # 3. Range de valores
        ax3 = fig.add_subplot(gs[1, 1])
        ranges = [s['max'] - s['min'] for s in stats]
        ax3.barh(output_names, ranges, alpha=0.7, color='coral')
        ax3.set_xlabel('Range (max - min)', fontsize=10)
        ax3.set_title('Range de Valores', fontsize=11, fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='x')
        
        # 4. Tabela de estatísticas
        ax4 = fig.add_subplot(gs[2, :])
        ax4.axis('tight')
        ax4.axis('off')
        
        table_data = []
        for s in stats:
            table_data.append([
                s['output'],
                f"{s['mean']:.4f}",
                f"{s['std']:.4f}",
                f"{s['min']:.4f}",
                f"{s['max']:.4f}",
                s['num_tracks']
            ])
        
        table = ax4.table(
            cellText=table_data,
            colLabels=['Output', 'Média', 'Desvio Padrão', 'Mínimo', 'Máximo', 'N Tracks'],
            cellLoc='center',
            loc='center',
            colWidths=[0.2, 0.15, 0.15, 0.15, 0.15, 0.1]
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.5)
        
        # Estilizar header
        for i in range(6):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        fig.suptitle(
            f'{seq_id} - Dashboard de Análise',
            fontsize=16,
            fontweight='bold',
            y=0.98
        )
        
        # Salvar
        for fmt in config['output_formats']:
            output_file = output_dir / f"{seq_id}_dashboard.{fmt}"
            plt.savefig(output_file, dpi=config['plot_resolution'], bbox_inches='tight')
            console.print(f"[green]  ✓ Salvo (dashboard): {output_file.name}[/green]")
        
        plt.close(fig)
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro em dashboard: {e}[/yellow]")


def create_variant_comparison_enhanced(
    results: Dict,
    output_dir: Path,
    config: Dict
):
    """
    Visualização melhorada de comparação de variantes REF vs ALT.
    
    Baseado em: https://www.alphagenomedocs.com/colabs/quick_start.html
    """
    try:
        seq_id = results['sequence_id']
        outputs = results['outputs']
        variant = results['variant']
        
        ref_data = outputs.reference.rna_seq
        alt_data = outputs.alternate.rna_seq
        
        if not (hasattr(ref_data, 'values') and hasattr(alt_data, 'values')):
            return
        
        ref_values = ref_data.values
        alt_values = alt_data.values
        
        if len(ref_values.shape) > 1:
            ref_values = ref_values[:, 0]
            alt_values = alt_values[:, 0]
        
        # Calcular diferença
        diff = alt_values - ref_values
        
        # Criar figura com 3 subplots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(config['plot_width'], 12), sharex=True)
        
        # 1. Sobreposição REF vs ALT
        ax1.plot(ref_values, linewidth=1.2, color='dimgrey', label='Referência', alpha=0.8)
        ax1.plot(alt_values, linewidth=1.2, color='red', label='Alternativa', alpha=0.8)
        ax1.fill_between(range(len(ref_values)), ref_values, alpha=0.2, color='dimgrey')
        ax1.fill_between(range(len(alt_values)), alt_values, alpha=0.2, color='red')
        
        # Marcar variante
        if variant.position >= ref_data.interval.start and variant.position <= ref_data.interval.end:
            var_pos = variant.position - ref_data.interval.start
            ax1.axvline(x=var_pos, color='orange', linestyle='--', linewidth=2.5, alpha=0.7, label='Variante')
            ax1.axvspan(var_pos-50, var_pos+50, alpha=0.1, color='orange')
        
        ax1.set_ylabel('Sinal RNA-seq', fontsize=11, fontweight='bold')
        ax1.legend(loc='upper right', framealpha=0.9)
        ax1.grid(True, alpha=0.3, linestyle='--')
        ax1.set_title('Comparação Referência vs Alternativa', fontsize=12, fontweight='bold')
        
        # 2. Diferença (ALT - REF)
        ax2.fill_between(range(len(diff)), diff, where=(diff>=0), interpolate=True, 
                        alpha=0.6, color='green', label='Aumento')
        ax2.fill_between(range(len(diff)), diff, where=(diff<0), interpolate=True,
                        alpha=0.6, color='red', label='Diminuição')
        ax2.plot(diff, linewidth=0.8, color='black', alpha=0.5)
        ax2.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3)
        
        if variant.position >= ref_data.interval.start and variant.position <= ref_data.interval.end:
            ax2.axvline(x=var_pos, color='orange', linestyle='--', linewidth=2.5, alpha=0.7)
        
        ax2.set_ylabel('Δ Sinal (ALT - REF)', fontsize=11, fontweight='bold')
        ax2.legend(loc='upper right', framealpha=0.9)
        ax2.grid(True, alpha=0.3, linestyle='--')
        ax2.set_title('Diferença de Sinal', fontsize=12, fontweight='bold')
        
        # 3. Zoom na região da variante
        if variant.position >= ref_data.interval.start and variant.position <= ref_data.interval.end:
            zoom_start = max(0, var_pos - 200)
            zoom_end = min(len(ref_values), var_pos + 200)
            
            ax3.plot(range(zoom_start, zoom_end), ref_values[zoom_start:zoom_end],
                    linewidth=2, color='dimgrey', label='Referência', alpha=0.8, marker='o', markersize=2)
            ax3.plot(range(zoom_start, zoom_end), alt_values[zoom_start:zoom_end],
                    linewidth=2, color='red', label='Alternativa', alpha=0.8, marker='s', markersize=2)
            ax3.axvline(x=var_pos, color='orange', linestyle='--', linewidth=3, alpha=0.8, label='Variante')
            ax3.axvspan(var_pos-10, var_pos+10, alpha=0.15, color='orange')
            
            ax3.set_xlim(zoom_start, zoom_end)
            ax3.set_ylabel('Sinal (Zoom)', fontsize=11, fontweight='bold')
            ax3.legend(loc='upper right', framealpha=0.9)
            ax3.grid(True, alpha=0.3)
            ax3.set_title(f'Zoom: ±200bp ao redor da variante', fontsize=12, fontweight='bold')
        
        ax3.set_xlabel('Posição (bp)', fontsize=11)
        
        fig.suptitle(
            f'{seq_id} - Efeito de Variante\n{variant.reference_bases}>{variant.alternate_bases} @ chr{variant.chromosome}:{variant.position}',
            fontsize=15,
            fontweight='bold',
            y=0.995
        )
        
        plt.tight_layout()
        
        # Salvar
        for fmt in config['output_formats']:
            output_file = output_dir / f"{seq_id}_variant_enhanced.{fmt}"
            plt.savefig(output_file, dpi=config['plot_resolution'], bbox_inches='tight')
            console.print(f"[green]  ✓ Salvo (variant enhanced): {output_file.name}[/green]")
        
        plt.close(fig)
        
    except Exception as e:
        console.print(f"[yellow]⚠ Erro em variante enhanced: {e}[/yellow]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")


if __name__ == '__main__':
    console.print("[bold cyan]Módulo de Visualizações Avançadas para AlphaGenome[/bold cyan]")
    console.print("Baseado em: https://www.alphagenomedocs.com/colabs/quick_start.html")
    console.print("\nFunções disponíveis:")
    console.print("  • create_enhanced_track_visualization()")
    console.print("  • create_multi_output_comparison()")
    console.print("  • create_heatmap_visualization()")
    console.print("  • create_summary_dashboard()")
    console.print("  • create_variant_comparison_enhanced()")

