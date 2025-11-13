#!/usr/bin/env python3
"""
alphagenome_output_visualization.py

Visualizador interativo de predições do AlphaGenome (.npz files).

Este programa permite navegar pelos resultados do AlphaGenome com controles de teclado:
  - Space / → : Próximo arquivo
  - ← : Arquivo anterior
  - a : Ativar modo overlay (sobrepor próximo arquivo)
  - A : Desativar modo overlay
  - d : Avançar tipo de output (atac → rna_seq → ...)
  - D : Retroceder tipo de output
  - g : Ativar/alternar modo grupo (atac + rna_seq juntos)
  - G : Desativar modo grupo
  - q / ESC : Sair

Uso:
    python3 alphagenome_output_visualization.py configs/small.yaml

Author: ChatGPT (for Alberto)
Date: 2025-11-13
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


# ═══════════════════════════════════════════════════════════════════
# Funções Auxiliares
# ═══════════════════════════════════════════════════════════════════

def load_config(config_path: Path) -> dict:
    """
    Carrega arquivo de configuração YAML.
    
    Caminhos relativos são resolvidos em relação ao diretório do config.
    """
    config_dir = config_path.parent
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Expandir variáveis de ambiente e resolver caminhos relativos
    def expand_env_vars(obj, is_path_field=False):
        if isinstance(obj, dict):
            return {k: expand_env_vars(v, k in ['metadata_csv', 'fasta', 'vcf_pattern', 'snp_list_file', 'gene_list_file', 'gtf_feather', 'output_dir']) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [expand_env_vars(item) for item in obj]
        elif isinstance(obj, str):
            # Expandir variáveis de ambiente
            if obj.startswith("${") and obj.endswith("}"):
                env_var = obj[2:-1]
                obj = os.environ.get(env_var, obj)
            
            # Resolver caminhos relativos (apenas para campos de path conhecidos)
            if is_path_field and not obj.startswith("/"):
                obj = str((config_dir / obj).resolve())
            
            return obj
        return obj
    
    return expand_env_vars(config)


def extract_file_info(path: Path) -> Dict[str, str]:
    """
    Extrai informações do arquivo a partir do path.
    
    Estrutura esperada:
    non_longevous_results/individuals/{SAMPLE}/windows/{SNP}/predictions_{HAPLOTYPE}/{OUTPUT}.npz
    
    Returns:
        Dict com keys: sample_id, snp, haplotype, output_type
    """
    parts = path.parts
    
    info = {
        'sample_id': 'Unknown',
        'snp': 'Unknown',
        'haplotype': 'Unknown',
        'output_type': path.stem  # Nome do arquivo sem extensão
    }
    
    try:
        # Encontrar índices das partes relevantes
        if 'individuals' in parts:
            idx = parts.index('individuals')
            info['sample_id'] = parts[idx + 1]
        
        if 'windows' in parts:
            idx = parts.index('windows')
            info['snp'] = parts[idx + 1]
        
        # Haplótipo está no nome do diretório predictions_H1 ou predictions_H2
        for part in parts:
            if part.startswith('predictions_'):
                info['haplotype'] = part.replace('predictions_', '')
                break
    
    except (IndexError, ValueError):
        pass
    
    return info


def plot_tracks(data: np.ndarray, ax: Axes, color, label: str, alpha: float = 0.7, 
                linewidth: float = 0.5):
    """
    Plota tracks de dados do AlphaGenome em um eixo.
    
    Args:
        data: Array numpy com shape (n_positions,) ou (n_positions, n_tracks)
        ax: Eixo matplotlib
        color: Cor para o plot
        label: Label para legenda
        alpha: Transparência
        linewidth: Largura da linha
    """
    if len(data.shape) == 1:
        # Apenas uma track
        ax.fill_between(range(len(data)), data, alpha=alpha*0.6, color=color, linewidth=0)
        ax.plot(data, linewidth=linewidth, color=color, alpha=alpha, label=label)
    else:
        # Múltiplas tracks - usar apenas a primeira
        track_data = data[:, 0]
        ax.fill_between(range(len(track_data)), track_data, alpha=alpha*0.6, color=color, linewidth=0)
        ax.plot(track_data, linewidth=linewidth, color=color, alpha=alpha, label=label)


# ═══════════════════════════════════════════════════════════════════
# Classe Principal - AlphaGenomeViewer
# ═══════════════════════════════════════════════════════════════════

class AlphaGenomeViewer:
    """
    Visualizador interativo de arquivos .npz do AlphaGenome.
    
    Permite navegação por teclado, overlay de múltiplos arquivos,
    e agrupamento de diferentes tipos de output.
    """
    
    def __init__(self, config_path: Path):
        """
        Inicializa o visualizador.
        
        Args:
            config_path: Caminho para arquivo YAML de configuração
        """
        print(f"[INFO] Carregando configuração: {config_path}")
        self.config = load_config(config_path)
        
        # Extrair output_dir do config
        self.output_dir = Path(self.config['project']['output_dir'])
        
        if not self.output_dir.exists():
            raise FileNotFoundError(f"Diretório de saída não encontrado: {self.output_dir}")
        
        print(f"[INFO] Diretório de saída: {self.output_dir}")
        
        # Encontrar todos os arquivos .npz
        self.all_files = self._find_npz_files()
        
        if not self.all_files:
            raise ValueError(f"Nenhum arquivo .npz encontrado em {self.output_dir}")
        
        print(f"[INFO] Total de arquivos .npz encontrados: {len(self.all_files)}")
        
        # Detectar tipos de output disponíveis
        self.output_types = self._detect_output_types()
        print(f"[INFO] Tipos de output detectados: {self.output_types}")
        
        # Estado inicial
        self.current_output_type_idx = 0
        self.current_files = self._filter_files(self.output_types[0])
        self.current_index = 0
        
        # Modos de visualização
        self.overlay_mode = False
        self.overlay_data = []  # Lista de (file_path, data, info)
        self.group_mode = False
        self.group_display_mode = 'stacked'  # 'stacked' ou 'overlaid'
        
        # Matplotlib
        self.fig = None
        self.current_data_cache = {}  # Cache para dados carregados
        
        print(f"[INFO] Arquivos do tipo '{self.output_types[0]}': {len(self.current_files)}")
        print("\n" + "="*70)
        print("CONTROLES:")
        print("  Space / → : Próximo arquivo")
        print("  ← : Arquivo anterior")
        print("  a : Ativar overlay (sobrepor próximo)")
        print("  A : Desativar overlay")
        print("  d : Próximo tipo de output")
        print("  D : Tipo de output anterior")
        print("  g : Ativar/alternar modo grupo (atac+rna_seq)")
        print("  G : Desativar modo grupo")
        print("  q / ESC : Sair")
        print("="*70 + "\n")
    
    def _find_npz_files(self) -> List[Path]:
        """
        Encontra todos os arquivos .npz no diretório de saída.
        
        Returns:
            Lista de Paths para arquivos .npz, ordenada
        """
        npz_files = sorted(self.output_dir.rglob("*.npz"))
        return npz_files
    
    def _detect_output_types(self) -> List[str]:
        """
        Detecta os tipos de output disponíveis (atac, rna_seq, etc).
        
        Returns:
            Lista de tipos de output únicos
        """
        types = set()
        for file_path in self.all_files:
            types.add(file_path.stem)  # Nome do arquivo sem extensão
        
        return sorted(list(types))
    
    def _filter_files(self, output_type: str) -> List[Path]:
        """
        Filtra arquivos por tipo de output.
        
        Args:
            output_type: Tipo de output (ex: 'atac', 'rna_seq')
        
        Returns:
            Lista filtrada de Paths
        """
        filtered = [f for f in self.all_files if f.stem == output_type]
        return filtered
    
    def _load_npz(self, file_path: Path) -> Optional[np.ndarray]:
        """
        Carrega arquivo .npz e retorna os dados.
        
        Args:
            file_path: Caminho para arquivo .npz
        
        Returns:
            Array numpy com os dados, ou None se erro
        """
        try:
            # Verificar cache
            if file_path in self.current_data_cache:
                return self.current_data_cache[file_path]
            
            data_npz = np.load(file_path)
            
            # Procurar por arrays de tracks
            # Prioridade: 'values', depois 'track_0', depois primeiro array disponível
            if 'values' in data_npz:
                data = data_npz['values']
            elif 'track_0' in data_npz:
                # Se tem múltiplas tracks, concatenar as primeiras 4
                tracks = []
                for i in range(min(4, len(data_npz.files))):
                    track_name = f'track_{i}'
                    if track_name in data_npz:
                        tracks.append(data_npz[track_name])
                
                if tracks:
                    # Empilhar tracks como colunas
                    data = np.column_stack(tracks) if len(tracks) > 1 else tracks[0]
                else:
                    data = None
            else:
                # Pegar primeiro array disponível
                keys = [k for k in data_npz.files if not k.endswith('_metadata')]
                if keys:
                    data = data_npz[keys[0]]
                else:
                    data = None
            
            # Adicionar ao cache
            if data is not None:
                self.current_data_cache[file_path] = data
            
            return data
        
        except Exception as e:
            print(f"[ERRO] Falha ao carregar {file_path}: {e}")
            return None
    
    def _plot_single(self):
        """
        Plota um único arquivo atual.
        """
        if self.current_index >= len(self.current_files):
            print("[AVISO] Índice fora do range")
            return
        
        file_path = self.current_files[self.current_index]
        print(f"\n[PLOT] {file_path}")
        
        # Carregar dados
        data = self._load_npz(file_path)
        if data is None:
            print("[ERRO] Não foi possível carregar dados")
            return
        
        # Extrair informações
        info = extract_file_info(file_path)
        
        # Criar figura
        if self.fig is not None:
            plt.close(self.fig)
        
        # Determinar número de subplots
        if len(data.shape) == 1:
            n_tracks = 1
            data_to_plot = [data]
        else:
            n_tracks = min(4, data.shape[1])  # Limitar a 4 tracks
            data_to_plot = [data[:, i] for i in range(n_tracks)]
        
        self.fig, axes = plt.subplots(n_tracks, 1, figsize=(14, 3*n_tracks), squeeze=False)
        
        # Plotar cada track
        for i, track_data in enumerate(data_to_plot):
            ax = axes[i, 0]
            
            color = plt.cm.tab10(i % 10)
            ax.fill_between(range(len(track_data)), track_data, alpha=0.4, color=color, linewidth=0)
            ax.plot(track_data, linewidth=0.7, color=color, alpha=0.8)
            
            ax.set_ylabel(f'Track {i}', fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # Estatísticas
            stats_text = f'μ={np.mean(track_data):.4f} σ={np.std(track_data):.4f}'
            ax.text(0.02, 0.95, stats_text, transform=ax.transAxes,
                   fontsize=8, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Xlabel apenas no último
        axes[-1, 0].set_xlabel('Posição (bp)', fontsize=10)
        
        # Título
        title = f"{info['sample_id']} | {info['snp']} | {info['haplotype']} | {info['output_type']}"
        self.fig.suptitle(title, fontsize=14, fontweight='bold')
        
        # Info no console
        print(f"  Sample: {info['sample_id']}")
        print(f"  SNP: {info['snp']}")
        print(f"  Haplótipo: {info['haplotype']}")
        print(f"  Output: {info['output_type']}")
        print(f"  Shape: {data.shape}")
        print(f"  [{self.current_index + 1}/{len(self.current_files)}]")
        
        plt.tight_layout()
    
    def _plot_overlay(self):
        """
        Plota múltiplos arquivos sobrepostos.
        """
        if not self.overlay_data:
            print("[AVISO] Nenhum dado para overlay")
            return
        
        print(f"\n[OVERLAY] {len(self.overlay_data)} arquivos")
        
        # Criar figura
        if self.fig is not None:
            plt.close(self.fig)
        
        # Determinar número máximo de tracks
        max_tracks = 1
        for _, data, _ in self.overlay_data:
            if len(data.shape) > 1:
                max_tracks = max(max_tracks, min(4, data.shape[1]))
        
        self.fig, axes = plt.subplots(max_tracks, 1, figsize=(14, 3*max_tracks), squeeze=False)
        
        # Cores contrastantes para overlay (azul, vermelho, verde, laranja, roxo, marrom, rosa, cinza, etc.)
        color_list = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#9467bd', 
                      '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        colors = [color_list[i % len(color_list)] for i in range(len(self.overlay_data))]
        
        # Plotar cada arquivo
        for file_idx, (file_path, data, info) in enumerate(self.overlay_data):
            print(f"  [{file_idx+1}] {file_path}")
            
            # Preparar dados
            if len(data.shape) == 1:
                data_to_plot = [data]
            else:
                n_tracks = min(max_tracks, data.shape[1])
                data_to_plot = [data[:, i] for i in range(n_tracks)]
            
            # Plotar em cada subplot
            for track_idx, track_data in enumerate(data_to_plot):
                ax = axes[track_idx, 0]
                
                color = colors[file_idx]
                label = f"{info['sample_id']}_{info['snp']}_{info['haplotype']}"
                
                ax.fill_between(range(len(track_data)), track_data, 
                               alpha=0.2, color=color, linewidth=0)
                ax.plot(track_data, linewidth=0.8, color=color, 
                       alpha=0.7, label=label)
                
                if file_idx == 0:  # Configurar apenas na primeira iteração
                    ax.set_ylabel(f'Track {track_idx}', fontsize=10, fontweight='bold')
                    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                
                # Legenda
                ax.legend(loc='upper right', fontsize=8, framealpha=0.8)
        
        # Xlabel
        axes[-1, 0].set_xlabel('Posição (bp)', fontsize=10)
        
        # Título
        output_type = self.overlay_data[0][2]['output_type']
        title = f"Overlay: {len(self.overlay_data)} arquivos | {output_type}"
        self.fig.suptitle(title, fontsize=14, fontweight='bold')
        
        plt.tight_layout()
    
    def _plot_grouped(self):
        """
        Plota diferentes tipos de output juntos (atac + rna_seq).
        """
        if len(self.output_types) < 2:
            print("[AVISO] Necessário pelo menos 2 tipos de output para agrupar")
            self.group_mode = False
            return
        
        # Pegar arquivo atual de cada tipo
        file_groups = {}
        
        for output_type in self.output_types:
            filtered = self._filter_files(output_type)
            if self.current_index < len(filtered):
                file_groups[output_type] = filtered[self.current_index]
        
        if len(file_groups) < 2:
            print("[AVISO] Não há arquivos correspondentes para agrupar")
            self.group_mode = False
            return
        
        print(f"\n[GRUPO - {self.group_display_mode.upper()}]")
        
        # Criar figura
        if self.fig is not None:
            plt.close(self.fig)
        
        if self.group_display_mode == 'stacked':
            # Subplots empilhados
            n_types = len(file_groups)
            self.fig, axes = plt.subplots(n_types, 1, figsize=(14, 4*n_types), squeeze=False)
            
            for idx, (output_type, file_path) in enumerate(file_groups.items()):
                print(f"  [{output_type}] {file_path}")
                
                data = self._load_npz(file_path)
                if data is None:
                    continue
                
                info = extract_file_info(file_path)
                ax = axes[idx, 0]
                
                # Pegar primeira track
                if len(data.shape) == 1:
                    track_data = data
                else:
                    track_data = data[:, 0]
                
                color = plt.cm.Set2(idx)
                ax.fill_between(range(len(track_data)), track_data, 
                               alpha=0.5, color=color, linewidth=0)
                ax.plot(track_data, linewidth=0.8, color=color, alpha=0.8)
                
                ax.set_ylabel(output_type.upper(), fontsize=11, fontweight='bold')
                ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                
                # Estatísticas
                stats_text = f'μ={np.mean(track_data):.4f} σ={np.std(track_data):.4f}'
                ax.text(0.02, 0.95, stats_text, transform=ax.transAxes,
                       fontsize=8, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            axes[-1, 0].set_xlabel('Posição (bp)', fontsize=10)
            
            # Usar info do primeiro arquivo para título
            first_info = extract_file_info(list(file_groups.values())[0])
            title = f"{first_info['sample_id']} | {first_info['snp']} | {first_info['haplotype']} | GRUPO"
            
        else:  # overlaid
            # Mesmo subplot, eixos Y diferentes
            self.fig, ax1 = plt.subplots(1, 1, figsize=(14, 6))
            
            axes_list = [ax1]
            colors = plt.cm.Set2(np.linspace(0, 1, len(file_groups)))
            
            for idx, (output_type, file_path) in enumerate(file_groups.items()):
                print(f"  [{output_type}] {file_path}")
                
                data = self._load_npz(file_path)
                if data is None:
                    continue
                
                # Pegar primeira track
                if len(data.shape) == 1:
                    track_data = data
                else:
                    track_data = data[:, 0]
                
                if idx == 0:
                    ax = ax1
                    color = colors[idx]
                    ax.plot(track_data, linewidth=1.0, color=color, alpha=0.8, label=output_type.upper())
                    ax.set_ylabel(output_type.upper(), fontsize=11, fontweight='bold', color=color)
                    ax.tick_params(axis='y', labelcolor=color)
                    ax.spines['top'].set_visible(False)
                else:
                    ax = ax1.twinx()
                    axes_list.append(ax)
                    color = colors[idx]
                    ax.plot(track_data, linewidth=1.0, color=color, alpha=0.8, label=output_type.upper())
                    ax.set_ylabel(output_type.upper(), fontsize=11, fontweight='bold', color=color)
                    ax.tick_params(axis='y', labelcolor=color)
                    ax.spines['top'].set_visible(False)
                    
                    # Posicionar eixos adicionais
                    if idx > 1:
                        ax.spines['right'].set_position(('outward', 60 * (idx - 1)))
            
            ax1.set_xlabel('Posição (bp)', fontsize=10)
            ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            
            # Legenda combinada
            lines, labels = ax1.get_legend_handles_labels()
            for ax in axes_list[1:]:
                l, lab = ax.get_legend_handles_labels()
                lines.extend(l)
                labels.extend(lab)
            ax1.legend(lines, labels, loc='upper right', fontsize=10, framealpha=0.9)
            
            # Usar info do primeiro arquivo para título
            first_info = extract_file_info(list(file_groups.values())[0])
            title = f"{first_info['sample_id']} | {first_info['snp']} | {first_info['haplotype']} | GRUPO (Overlaid)"
        
        self.fig.suptitle(title, fontsize=14, fontweight='bold')
        plt.tight_layout()
    
    def _on_key_press(self, event):
        """
        Handler de eventos de teclado.
        
        Args:
            event: Evento de teclado matplotlib
        """
        key = event.key
        
        # Sair
        if key in ['q', 'escape']:
            print("\n[INFO] Encerrando visualizador...")
            plt.close('all')
            sys.exit(0)
        
        # Navegação
        elif key in [' ', 'right']:  # Space ou seta direita
            # Avançar primeiro
            self.current_index = min(self.current_index + 1, len(self.current_files) - 1)
            
            if self.overlay_mode:
                # Adicionar o NOVO arquivo ao overlay (após avançar)
                file_path = self.current_files[self.current_index]
                data = self._load_npz(file_path)
                if data is not None:
                    info = extract_file_info(file_path)
                    self.overlay_data.append((file_path, data, info))
                    print(f"[OVERLAY] Adicionado: {file_path.name}")
            
            self._update_plot()
        
        elif key == 'left':  # Seta esquerda
            self.current_index = max(self.current_index - 1, 0)
            self._update_plot()
        
        # Overlay
        elif key == 'a':
            self.overlay_mode = True
            # Adicionar arquivo atual (sem avançar)
            file_path = self.current_files[self.current_index]
            data = self._load_npz(file_path)
            if data is not None:
                info = extract_file_info(file_path)
                self.overlay_data.append((file_path, data, info))
                print(f"[OVERLAY ATIVADO] Arquivo base: {file_path.name}")
                print(f"[INFO] Pressione Space para adicionar mais arquivos ao overlay")
            self._update_plot()
        
        elif key == 'A':
            self.overlay_mode = False
            self.overlay_data.clear()
            print("[OVERLAY DESATIVADO]")
            self._update_plot()
        
        # Trocar tipo de output
        elif key == 'd':
            self.overlay_mode = False
            self.overlay_data.clear()
            self.group_mode = False
            
            self.current_output_type_idx = (self.current_output_type_idx + 1) % len(self.output_types)
            new_type = self.output_types[self.current_output_type_idx]
            self.current_files = self._filter_files(new_type)
            self.current_index = min(self.current_index, len(self.current_files) - 1)
            
            print(f"\n[TIPO] {new_type} ({len(self.current_files)} arquivos)")
            self._update_plot()
        
        elif key == 'D':
            self.overlay_mode = False
            self.overlay_data.clear()
            self.group_mode = False
            
            self.current_output_type_idx = (self.current_output_type_idx - 1) % len(self.output_types)
            new_type = self.output_types[self.current_output_type_idx]
            self.current_files = self._filter_files(new_type)
            self.current_index = min(self.current_index, len(self.current_files) - 1)
            
            print(f"\n[TIPO] {new_type} ({len(self.current_files)} arquivos)")
            self._update_plot()
        
        # Modo grupo
        elif key == 'g':
            self.overlay_mode = False
            self.overlay_data.clear()
            
            if self.group_mode:
                # Alternar modo de display
                self.group_display_mode = 'overlaid' if self.group_display_mode == 'stacked' else 'stacked'
                print(f"[GRUPO] Modo alterado para: {self.group_display_mode}")
            else:
                self.group_mode = True
                print("[GRUPO ATIVADO] Modo: stacked")
            
            self._update_plot()
        
        elif key == 'G':
            self.group_mode = False
            print("[GRUPO DESATIVADO]")
            self._update_plot()
        
        else:
            # Tecla não reconhecida
            pass
    
    def _update_plot(self):
        """
        Atualiza o plot baseado no estado atual.
        """
        if self.overlay_mode and self.overlay_data:
            self._plot_overlay()
        elif self.group_mode:
            self._plot_grouped()
        else:
            self._plot_single()
        
        # Reconectar eventos de teclado (necessário após criar nova figura)
        if self.fig is not None:
            self.fig.canvas.mpl_connect('key_press_event', self._on_key_press)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
            plt.show(block=False)
            plt.pause(0.001)  # Pequena pausa para processar eventos
    
    def show(self):
        """
        Inicia a visualização interativa.
        """
        # Ativar modo interativo
        plt.ion()
        
        # Plot inicial
        self._update_plot()
        
        # Manter o programa rodando
        print("\n[INFO] Visualizador ativo. Use as teclas para navegar.")
        print("[INFO] Pressione 'q' ou ESC para sair.\n")
        
        # Bloquear até que todas as janelas sejam fechadas ou usuário saia
        try:
            plt.show(block=True)
        except KeyboardInterrupt:
            print("\n[INFO] Interrompido pelo usuário")
            plt.close('all')


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Visualizador interativo de predições do AlphaGenome",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Controles de teclado:
  Space / → : Próximo arquivo
  ← : Arquivo anterior
  a : Ativar overlay (sobrepor próximo arquivo)
  A : Desativar overlay
  d : Próximo tipo de output (atac → rna_seq → ...)
  D : Tipo de output anterior
  g : Ativar/alternar modo grupo (atac + rna_seq juntos)
  G : Desativar modo grupo
  q / ESC : Sair

Exemplos:
  python3 alphagenome_output_visualization.py configs/small.yaml
        """
    )
    
    parser.add_argument(
        'config',
        type=Path,
        help='Arquivo de configuração YAML (ex: configs/small.yaml)'
    )
    
    args = parser.parse_args()
    
    if not args.config.exists():
        print(f"[ERRO] Arquivo de configuração não encontrado: {args.config}")
        return 1
    
    try:
        viewer = AlphaGenomeViewer(args.config)
        viewer.show()
    except KeyboardInterrupt:
        print("\n[INFO] Interrompido pelo usuário")
        return 0
    except Exception as e:
        print(f"[ERRO] {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

