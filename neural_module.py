#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Neural Module - Análise de DNA usando AlphaGenome
==================================================

Este módulo utiliza a API do AlphaGenome (Google DeepMind) para análise
avançada de sequências de DNA, incluindo:
- Predição de expressão gênica
- Padrões de splicing
- Características de cromatina
- Mapas de contato
- Efeitos de variantes

Uso:
    python neural_module.py -i input.fasta -o output_dir -k API_KEY [opções]
"""

import argparse
import os
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Tuple
import json

from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.table import Table
from rich import box

console = Console()

# ====================== Configuração ======================

DEFAULT_CONFIG = {
    'supported_lengths': [2048, 16384, 131072, 524288, 1048576],  # Tamanhos suportados pelo AlphaGenome
    'max_sequence_length': 1048576,  # AlphaGenome suporta até 1M bp
    'min_sequence_length': 2048,  # Tamanho mínimo: 2kb
    'output_formats': ['png', 'pdf', 'svg'],
    'default_outputs': [
        'RNA_SEQ',
        'CAGE',
        'ATAC',
        'CHIP_HISTONE',  # Marcadores de histonas (H3K27AC, H3K4ME3, etc.)
        'CHIP_TF',       # Fatores de transcrição (CTCF, etc.)
    ],
    'plot_resolution': 300,  # DPI
    'plot_width': 15,  # polegadas
    'plot_height': 10,  # polegadas
    'use_advanced_viz': True,  # Visualizações avançadas habilitadas por padrão
    'show_ontology_info': True,  # Mostrar informações de ontologia por padrão
    'save_metadata': True,  # Salvar metadados de ontologia por padrão
}

# ====================== Parse de FASTA ======================

def parse_fasta(fasta_path: Path) -> List[Dict[str, str]]:
    """
    Faz parse de um arquivo FASTA e retorna lista de sequências.
    
    Args:
        fasta_path: Caminho para o arquivo FASTA
        
    Returns:
        Lista de dicionários com 'id', 'description' e 'sequence'
    """
    sequences = []
    current_id = None
    current_desc = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Salvar sequência anterior se existir
                if current_id is not None:
                    sequences.append({
                        'id': current_id,
                        'description': current_desc,
                        'sequence': ''.join(current_seq),
                        'length': len(''.join(current_seq))
                    })
                
                # Parse do header
                header = line[1:].split(None, 1)
                current_id = header[0]
                current_desc = header[1] if len(header) > 1 else ''
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        # Adicionar última sequência
        if current_id is not None:
            sequences.append({
                'id': current_id,
                'description': current_desc,
                'sequence': ''.join(current_seq),
                'length': len(''.join(current_seq))
            })
    
    return sequences


def validate_sequence(seq: str) -> Tuple[bool, Optional[str]]:
    """
    Valida uma sequência de DNA.
    
    Args:
        seq: Sequência de DNA
        
    Returns:
        Tupla (válida, mensagem_erro)
    """
    # Verifica caracteres válidos (ACGT e ambiguidades)
    valid_chars = set('ACGTNRYSWKMBDHV-')
    invalid_chars = set(seq) - valid_chars
    
    if invalid_chars:
        return False, f"Caracteres inválidas encontrados: {', '.join(sorted(invalid_chars))}"
    
    # Verifica se o tamanho é um dos suportados
    supported = DEFAULT_CONFIG['supported_lengths']
    if len(seq) not in supported:
        supported_str = ', '.join([f"{x:,} bp" for x in supported])
        return False, f"Tamanho {len(seq):,} bp não suportado. Tamanhos válidos: {supported_str}"
    
    return True, None


# ====================== Integração com AlphaGenome ======================

class AlphaGenomeAnalyzer:
    """
    Classe para integração com a API do AlphaGenome.
    """
    
    def __init__(self, api_key: str, config: Dict = None):
        """
        Inicializa o analisador AlphaGenome.
        
        Args:
            api_key: Chave da API do AlphaGenome
            config: Configurações adicionais
        """
        self.api_key = api_key
        self.config = config or DEFAULT_CONFIG
        self.model = None
        self._initialized = False
        
    def initialize(self) -> bool:
        """
        Inicializa a conexão com a API do AlphaGenome.
        
        Returns:
            True se inicialização bem-sucedida
        """
        try:
            from alphagenome.models import dna_client
            from alphagenome.data import genome
            
            console.print("[cyan]Inicializando conexão com AlphaGenome...[/cyan]")
            self.model = dna_client.create(self.api_key)
            self.dna_client = dna_client
            self.genome = genome
            self._initialized = True
            console.print("[green]✓ Conexão estabelecida com sucesso![/green]")
            return True
            
        except ImportError as e:
            console.print(f"[red]✗ Erro: AlphaGenome não está instalado.[/red]")
            console.print(f"[yellow]Execute: git clone https://github.com/google-deepmind/alphagenome.git && pip install ./alphagenome[/yellow]")
            return False
        except Exception as e:
            console.print(f"[red]✗ Erro ao inicializar AlphaGenome: {e}[/red]")
            return False
    
    def predict_sequence(self, 
                        sequence: str, 
                        seq_id: str,
                        chromosome: str = 'chr1',
                        start: int = 1000000,
                        requested_outputs: List[str] = None,
                        ontology_terms: List[str] = None) -> Optional[Dict]:
        """
        Faz predições para uma sequência de DNA.
        
        Args:
            sequence: Sequência de DNA
            seq_id: Identificador da sequência
            chromosome: Cromossomo (para contexto)
            start: Posição inicial (para contexto)
            requested_outputs: Lista de tipos de saída desejados
            
        Returns:
            Dicionário com resultados das predições
        """
        if not self._initialized:
            console.print("[red]Erro: AlphaGenome não está inicializado[/red]")
            return None
        
        try:
            # Preparar outputs solicitados
            if requested_outputs is None:
                requested_outputs = self.config['default_outputs']
            
            # Tentar converter outputs para tipos do AlphaGenome
            output_types = []
            for out in requested_outputs:
                try:
                    output_type = getattr(self.dna_client.OutputType, out)
                    output_types.append(output_type)
                except AttributeError:
                    console.print(f"[yellow]⚠ Output '{out}' não disponível, pulando...[/yellow]")
                    continue
            
            if not output_types:
                console.print(f"[red]✗ Nenhum output válido especificado[/red]")
                return None
            
            # Criar intervalo genômico
            end = start + len(sequence)
            interval = self.genome.Interval(
                chromosome=chromosome,
                start=start,
                end=end
            )
            
            console.print(f"[cyan]Fazendo predições para {seq_id} ({len(sequence)} bp)...[/cyan]")
            console.print(f"[dim]Outputs solicitados: {', '.join([str(ot).split('.')[-1] for ot in output_types])}[/dim]")
            
            # Ontology terms (tecidos/tipos celulares)
            if ontology_terms is None:
                # Usar termos padrão: cérebro, fígado, coração
                ontology_terms = ['UBERON:0000955', 'UBERON:0002107', 'UBERON:0000948']
            
            # Fazer predição usando predict_interval
            outputs = self.model.predict_interval(
                interval=interval,
                ontology_terms=ontology_terms,
                requested_outputs=output_types
            )
            
            return {
                'sequence_id': seq_id,
                'interval': interval,
                'outputs': outputs,
                'requested_outputs': [str(ot).split('.')[-1] for ot in output_types],
                'sequence_length': len(sequence)
            }
            
        except Exception as e:
            console.print(f"[red]✗ Erro ao processar sequência {seq_id}: {e}[/red]")
            import traceback
            console.print(f"[dim]{traceback.format_exc()}[/dim]")
            return None
    
    def predict_variant(self,
                       sequence: str,
                       seq_id: str,
                       variant_position: int,
                       ref_base: str,
                       alt_base: str,
                       chromosome: str = 'chr1',
                       start: int = 1000000,
                       ontology_terms: List[str] = None) -> Optional[Dict]:
        """
        Prediz o efeito de uma variante.
        
        Args:
            sequence: Sequência de DNA de referência
            seq_id: Identificador da sequência
            variant_position: Posição da variante (relativa ao início da sequência)
            ref_base: Base de referência
            alt_base: Base alternativa
            chromosome: Cromossomo
            start: Posição inicial
            ontology_terms: Lista de termos de ontologia UBERON (opcional)
            
        Returns:
            Dicionário com predições para referência e alternativa
        """
        if not self._initialized:
            console.print("[red]Erro: AlphaGenome não está inicializado[/red]")
            return None
        
        try:
            end = start + len(sequence)
            interval = self.genome.Interval(
                chromosome=chromosome,
                start=start,
                end=end
            )
            
            # Posição absoluta da variante
            variant_abs_position = start + variant_position
            
            variant = self.genome.Variant(
                chromosome=chromosome,
                position=variant_abs_position,
                reference_bases=ref_base,
                alternate_bases=alt_base,
            )
            
            console.print(f"[cyan]Analisando variante {ref_base}>{alt_base} na posição {variant_position}...[/cyan]")
            
            # Ontology terms (tecidos/tipos celulares)
            if ontology_terms is None:
                # Usar termos padrão: cérebro, fígado, coração
                ontology_terms = ['UBERON:0000955', 'UBERON:0002107', 'UBERON:0000948']
            
            outputs = self.model.predict_variant(
                interval=interval,
                variant=variant,
                ontology_terms=ontology_terms,
                requested_outputs=[self.dna_client.OutputType.RNA_SEQ]
            )
            
            return {
                'sequence_id': seq_id,
                'variant': variant,
                'interval': interval,
                'outputs': outputs
            }
            
        except Exception as e:
            console.print(f"[red]✗ Erro ao analisar variante: {e}[/red]")
            return None


# ====================== Informações de Ontologia ======================

def display_ontology_info(ontology_terms: List[str] = None):
    """
    Exibe informações sobre os termos de ontologia usados na análise.
    
    Args:
        ontology_terms: Lista de termos UBERON
    """
    # Dicionário de termos comuns
    ontology_dict = {
        'UBERON:0000955': 'Brain (Cérebro)',
        'UBERON:0002107': 'Liver (Fígado)',
        'UBERON:0000948': 'Heart (Coração)',
        'UBERON:0002048': 'Lung (Pulmão)',
        'UBERON:0000178': 'Blood (Sangue)',
        'UBERON:0002113': 'Kidney (Rim)',
        'UBERON:0002106': 'Spleen (Baço)',
        'UBERON:0000970': 'Eye (Olho)',
        'UBERON:0001723': 'Tongue (Língua)',
        'UBERON:0001255': 'Urinary bladder (Bexiga)',
        'UBERON:0002367': 'Prostate gland (Próstata)',
        'UBERON:0000473': 'Testis (Testículo)',
        'UBERON:0000992': 'Ovary (Ovário)',
        'UBERON:0001264': 'Pancreas (Pâncreas)',
        'UBERON:0002097': 'Skin of body (Pele)',
        'UBERON:0002371': 'Bone marrow (Medula óssea)',
        'UBERON:0002240': 'Spinal cord (Medula espinhal)',
        'UBERON:0000945': 'Stomach (Estômago)',
        'UBERON:0002110': 'Gall bladder (Vesícula biliar)',
    }
    
    if ontology_terms is None:
        ontology_terms = ['UBERON:0000955', 'UBERON:0002107', 'UBERON:0000948']
    
    table = Table(title="🧬 Contextos Biológicos (Ontologia UBERON)", 
                  box=box.ROUNDED, 
                  header_style="bold cyan")
    table.add_column("Termo UBERON", style="cyan")
    table.add_column("Tecido/Órgão", style="green")
    table.add_column("Descrição", style="white")
    
    descriptions = {
        'UBERON:0000955': 'Tecido neural - expressão em neurônios e glia',
        'UBERON:0002107': 'Tecido hepático - metabolismo e detoxificação',
        'UBERON:0000948': 'Tecido cardíaco - função cardiovascular'
    }
    
    for term in ontology_terms:
        tissue_name = ontology_dict.get(term, 'Desconhecido')
        desc = descriptions.get(term, 'Contexto tecido-específico')
        table.add_row(term, tissue_name, desc)
    
    console.print("\n")
    console.print(table)
    console.print("\n[dim]ℹ️  AlphaGenome usa esses contextos para predizer padrões tecido-específicos[/dim]")
    console.print("[dim]   As visualizações mostram como cada output varia entre esses tecidos[/dim]\n")


def save_metadata_to_file(results: Dict, output_dir: Path):
    """
    Extrai e salva metadados de ontologia dos outputs.
    
    Args:
        results: Resultados das predições com outputs
        output_dir: Diretório para salvar metadados
    """
    try:
        import pandas as pd
        
        seq_id = results['sequence_id']
        outputs = results['outputs']
        
        console.print(f"\n[cyan]Extraindo metadados de ontologia para {seq_id}...[/cyan]")
        
        metadata_found = False
        
        # Para cada output, extrair metadata
        for output_name in results['requested_outputs']:
            output_data = getattr(outputs, output_name.lower(), None)
            
            if output_data is None:
                continue
            
            # Verificar se tem metadata
            if hasattr(output_data, 'metadata') and output_data.metadata is not None:
                metadata_found = True
                metadata_df = output_data.metadata
                
                # Salvar como CSV
                csv_file = output_dir / f"{seq_id}_{output_name}_metadata.csv"
                metadata_df.to_csv(csv_file, index=False)
                console.print(f"[green]  ✓ Metadados salvos: {csv_file.name}[/green]")
                
                # Salvar como JSON (mais legível)
                json_file = output_dir / f"{seq_id}_{output_name}_metadata.json"
                metadata_df.to_json(json_file, orient='records', indent=2)
                console.print(f"[green]  ✓ Metadados JSON salvos: {json_file.name}[/green]")
                
                # Exibir resumo no terminal
                console.print(f"\n[bold cyan]Metadados de {output_name}:[/bold cyan]")
                
                # Criar tabela com principais informações
                metadata_table = Table(show_header=True, header_style="bold cyan", box=box.ROUNDED)
                
                # Adicionar colunas relevantes
                relevant_cols = ['biosample_name', 'strand', 'Assay title', 'File accession']
                available_cols = [col for col in relevant_cols if col in metadata_df.columns]
                
                if not available_cols:
                    available_cols = metadata_df.columns[:5].tolist()  # Primeiras 5 colunas
                
                for col in available_cols:
                    metadata_table.add_column(col, overflow='fold')
                
                # Adicionar linhas (máximo 10 tracks)
                for idx, row in metadata_df.head(10).iterrows():
                    values = [str(row[col])[:30] for col in available_cols]  # Limitar tamanho
                    metadata_table.add_row(*values)
                
                if len(metadata_df) > 10:
                    console.print(f"[dim]  (Mostrando 10 de {len(metadata_df)} tracks)[/dim]")
                
                console.print(metadata_table)
                console.print("")
            else:
                console.print(f"[yellow]  ⚠ {output_name}: Sem metadados disponíveis[/yellow]")
        
        if not metadata_found:
            console.print("[yellow]  ⚠ Nenhum metadado de ontologia encontrado nos outputs[/yellow]")
    
    except ImportError as ie:
        console.print(f"[yellow]⚠ Biblioteca pandas necessária para salvar metadados: {ie}[/yellow]")
    except Exception as e:
        console.print(f"[yellow]⚠ Erro ao salvar metadados: {e}[/yellow]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")


# ====================== Visualização ======================

def create_visualizations(results: Dict, output_dir: Path, config: Dict):
    """
    Cria visualizações para os resultados das predições.
    
    Args:
        results: Resultados das predições
        output_dir: Diretório de saída
        config: Configurações de visualização
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        seq_id = results['sequence_id']
        outputs = results['outputs']
        use_advanced = config.get('use_advanced_viz', False)
        
        console.print(f"[cyan]Gerando visualizações para {seq_id}...[/cyan]")
        
        # Se modo avançado, usar visualizações melhoradas
        if use_advanced:
            try:
                from neural_visualizations_advanced import (
                    create_enhanced_track_visualization,
                    create_heatmap_visualization,
                    create_multi_output_comparison,
                    create_summary_dashboard
                )
                
                # Visualizações avançadas para cada output
                for output_name in results['requested_outputs']:
                    output_data = getattr(outputs, output_name.lower(), None)
                    if output_data is not None:
                        create_enhanced_track_visualization(
                            output_data, seq_id, output_name, output_dir, config
                        )
                        create_heatmap_visualization(
                            output_data, seq_id, output_name, output_dir, config
                        )
                
                # Comparação multi-output
                create_multi_output_comparison(
                    outputs, seq_id, results['requested_outputs'], output_dir, config
                )
                
                # Dashboard resumo
                create_summary_dashboard(results, output_dir, config)
                
                console.print("[green]  ✓ Visualizações avançadas criadas![/green]")
                
            except ImportError:
                console.print("[yellow]  ⚠ Modo avançado não disponível, usando visualizações básicas[/yellow]")
                use_advanced = False
        
        # Visualizações básicas (sempre criar)
        for output_name in results['requested_outputs']:
            try:
                output_data = getattr(outputs, output_name.lower(), None)
                
                if output_data is None:
                    console.print(f"[yellow]  ⚠ Output {output_name} não encontrado[/yellow]")
                    continue
                
                # Extrair valores do TrackData
                if hasattr(output_data, 'values'):
                    data_array = output_data.values
                    
                    # Criar plot
                    fig, ax = plt.subplots(figsize=(config['plot_width'], config['plot_height']))
                    
                    # Se tem múltiplas tracks (2D), plotar cada uma
                    if len(data_array.shape) > 1 and data_array.shape[1] > 1:
                        for i in range(min(data_array.shape[1], 5)):  # Máximo 5 tracks
                            ax.plot(data_array[:, i], linewidth=0.5, alpha=0.7, label=f'Track {i+1}')
                        if data_array.shape[1] > 5:
                            ax.plot([], [], ' ', label=f'... e mais {data_array.shape[1]-5} tracks')
                        ax.legend(loc='best', fontsize=8)
                    else:
                        # Uma única track
                        if len(data_array.shape) > 1:
                            data_array = data_array[:, 0]
                        ax.plot(data_array, linewidth=0.5, color='#2E86AB')
                    
                    ax.set_title(f'{seq_id} - {output_name}', fontsize=16, fontweight='bold')
                    ax.set_xlabel('Posição (bp)', fontsize=12)
                    ax.set_ylabel('Sinal', fontsize=12)
                    ax.grid(True, alpha=0.3, linestyle='--')
                    
                    plt.tight_layout()
                    
                    # Salvar em múltiplos formatos
                    for fmt in config['output_formats']:
                        suffix = '' if use_advanced else ''
                        output_file = output_dir / f"{seq_id}_{output_name}{suffix}.{fmt}"
                        plt.savefig(
                            output_file,
                            dpi=config['plot_resolution'],
                            bbox_inches='tight'
                        )
                        console.print(f"[green]  ✓ Salvo: {output_file.name}[/green]")
                    
                    plt.close(fig)
                else:
                    console.print(f"[yellow]  ⚠ {output_name} não tem dados plotáveis[/yellow]")
                
            except Exception as e:
                console.print(f"[yellow]  ⚠ Aviso: Não foi possível criar plot para {output_name}: {e}[/yellow]")
        
    except ImportError as ie:
        console.print(f"[yellow]⚠ Biblioteca necessária não disponível: {ie}[/yellow]")
    except Exception as e:
        console.print(f"[red]✗ Erro ao criar visualizações: {e}[/red]")


def create_variant_visualization(results: Dict, output_dir: Path, config: Dict):
    """
    Cria visualização para análise de variante.
    
    Args:
        results: Resultados da predição de variante
        output_dir: Diretório de saída
        config: Configurações
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        seq_id = results['sequence_id']
        outputs = results['outputs']
        variant = results['variant']
        use_advanced = config.get('use_advanced_viz', False)
        
        console.print(f"[cyan]Gerando visualização de variante para {seq_id}...[/cyan]")
        
        # Se modo avançado, usar visualização melhorada
        if use_advanced:
            try:
                from neural_visualizations_advanced import create_variant_comparison_enhanced
                create_variant_comparison_enhanced(results, output_dir, config)
                console.print("[green]  ✓ Visualização de variante avançada criada![/green]")
            except ImportError:
                console.print("[yellow]  ⚠ Modo avançado não disponível para variantes[/yellow]")
                use_advanced = False
        
        # Visualização básica (sempre criar)
        ref_data = outputs.reference.rna_seq
        alt_data = outputs.alternate.rna_seq
        
        if hasattr(ref_data, 'values') and hasattr(alt_data, 'values'):
            ref_values = ref_data.values
            alt_values = alt_data.values
            
            # Se multidimensional, pegar primeira dimensão
            if len(ref_values.shape) > 1:
                ref_values = ref_values[:, 0]
                alt_values = alt_values[:, 0]
            
            # Criar plot
            fig, ax = plt.subplots(figsize=(config['plot_width'], config['plot_height']))
            
            ax.plot(ref_values, linewidth=0.8, color='dimgrey', label='Referência', alpha=0.7)
            ax.plot(alt_values, linewidth=0.8, color='red', label='Alternativa', alpha=0.7)
            
            # Marcar posição da variante (se dentro do plot)
            if variant.position >= ref_data.interval.start and variant.position <= ref_data.interval.end:
                var_pos = variant.position - ref_data.interval.start
                ax.axvline(x=var_pos, color='orange', linestyle='--', linewidth=2, alpha=0.5, label='Variante')
            
            ax.set_title(
                f'{seq_id} - Efeito de Variante\n{variant.reference_bases}>{variant.alternate_bases} @ pos {variant.position}',
                fontsize=16, 
                fontweight='bold'
            )
            ax.set_xlabel('Posição (bp)', fontsize=12)
            ax.set_ylabel('Sinal RNA-seq', fontsize=12)
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3, linestyle='--')
            
            plt.tight_layout()
            
            # Salvar
            for fmt in config['output_formats']:
                output_file = output_dir / f"{seq_id}_variant.{fmt}"
                plt.savefig(output_file, dpi=config['plot_resolution'], bbox_inches='tight')
                console.print(f"[green]  ✓ Salvo: {output_file.name}[/green]")
            
            plt.close(fig)
        else:
            console.print(f"[yellow]⚠ Dados de variante não têm formato plotável[/yellow]")
        
    except Exception as e:
        console.print(f"[red]✗ Erro ao criar visualização de variante: {e}[/red]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")


# ====================== Relatório ======================

def generate_report(sequences: List[Dict], results: List[Dict], output_dir: Path):
    """
    Gera relatório JSON com resumo das análises.
    
    Args:
        sequences: Lista de sequências processadas
        results: Lista de resultados
        output_dir: Diretório de saída
    """
    report = {
        'timestamp': datetime.now().isoformat(),
        'total_sequences': len(sequences),
        'successful_analyses': len(results),
        'sequences': []
    }
    
    for seq, res in zip(sequences, results):
        if res is None:
            report['sequences'].append({
                'id': seq['id'],
                'length': seq['length'],
                'status': 'failed'
            })
        else:
            report['sequences'].append({
                'id': seq['id'],
                'length': seq['length'],
                'status': 'success',
                'outputs': res.get('requested_outputs', [])
            })
    
    report_file = output_dir / 'analysis_report.json'
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    console.print(f"[green]✓ Relatório salvo: {report_file}[/green]")


def print_summary(sequences: List[Dict], results: List[Dict]):
    """
    Imprime resumo das análises.
    
    Args:
        sequences: Lista de sequências
        results: Lista de resultados
    """
    table = Table(title="Resumo das Análises", box=box.ROUNDED)
    table.add_column("ID da Sequência", style="cyan")
    table.add_column("Tamanho (bp)", justify="right")
    table.add_column("Status", justify="center")
    
    for seq, res in zip(sequences, results):
        status = "[green]✓ Sucesso[/green]" if res else "[red]✗ Falhou[/red]"
        table.add_row(seq['id'], f"{seq['length']:,}", status)
    
    console.print(table)


# ====================== Main ======================

def main():
    parser = argparse.ArgumentParser(
        description='Neural Module - Análise de DNA usando AlphaGenome',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  %(prog)s -i sequence.fasta -k YOUR_API_KEY -o results/
  %(prog)s -i sequence.fasta -k YOUR_API_KEY -o results/ --outputs RNA_SEQ ATAC H3K27AC
  %(prog)s -i sequence.fasta -k YOUR_API_KEY -o results/ --variant 1000 A C

Para obter uma API key: https://www.alphagenomedocs.com/
        """
    )
    
    # Argumentos obrigatórios
    parser.add_argument('-i', '--input', required=True, type=Path,
                       help='Arquivo FASTA de entrada')
    parser.add_argument('-k', '--api-key', required=True,
                       help='Chave da API do AlphaGenome')
    parser.add_argument('-o', '--output', required=True, type=Path,
                       help='Diretório de saída')
    
    # Argumentos opcionais
    parser.add_argument('--outputs', nargs='+',
                       choices=['RNA_SEQ', 'CAGE', 'ATAC', 'CHIP_HISTONE', 'CHIP_TF',
                               'DNASE', 'PROCAP', 'CONTACT_MAPS', 
                               'SPLICE_JUNCTIONS', 'SPLICE_SITES', 'SPLICE_SITE_USAGE'],
                       help='Tipos de output desejados (padrão: RNA_SEQ CAGE ATAC CHIP_HISTONE CHIP_TF)')
    parser.add_argument('--chromosome', default='chr1',
                       help='Cromossomo de referência (padrão: chr1)')
    parser.add_argument('--start', type=int, default=1000000,
                       help='Posição inicial de referência (padrão: 1000000)')
    parser.add_argument('--variant', nargs=3, metavar=('POS', 'REF', 'ALT'),
                       help='Analisar variante na posição POS (relativa) com bases REF>ALT')
    parser.add_argument('--formats', nargs='+', choices=['png', 'pdf', 'svg'],
                       default=['png'],
                       help='Formatos de saída para gráficos (padrão: png)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Resolução dos gráficos (padrão: 300)')
    parser.add_argument('--no-plots', action='store_true',
                       help='Não gerar gráficos (apenas análise)')
    parser.add_argument('--save-metadata', action='store_true',
                       help='Salvar metadados de ontologia em CSV/JSON (padrão: ativado)')
    parser.add_argument('--no-metadata', action='store_true',
                       help='Não salvar metadados de ontologia')
    
    args = parser.parse_args()
    
    # Banner
    console.print(Panel.fit(
        "[bold cyan]Neural Module - Análise de DNA com AlphaGenome[/bold cyan]\n"
        "[dim]Powered by Google DeepMind AlphaGenome[/dim]",
        border_style="cyan"
    ))
    
    # Validar entrada
    if not args.input.exists():
        console.print(f"[red]✗ Erro: Arquivo não encontrado: {args.input}[/red]")
        sys.exit(1)
    
    # Criar diretório de saída
    args.output.mkdir(parents=True, exist_ok=True)
    console.print(f"[green]✓ Diretório de saída: {args.output}[/green]")
    
    # Configuração
    config = DEFAULT_CONFIG.copy()
    config['output_formats'] = args.formats
    config['plot_resolution'] = args.dpi
    # use_advanced_viz e show_ontology_info já estão definidos em DEFAULT_CONFIG como True
    # Controlar save_metadata
    if args.no_metadata:
        config['save_metadata'] = False
    elif args.save_metadata:
        config['save_metadata'] = True
    # else: usa o valor de DEFAULT_CONFIG (True)
    
    if args.outputs:
        config['default_outputs'] = args.outputs
    
    # Parse FASTA
    console.print(f"\n[cyan]Lendo arquivo FASTA: {args.input}[/cyan]")
    try:
        sequences = parse_fasta(args.input)
        console.print(f"[green]✓ {len(sequences)} sequência(s) encontrada(s)[/green]")
    except Exception as e:
        console.print(f"[red]✗ Erro ao ler FASTA: {e}[/red]")
        sys.exit(1)
    
    # Validar sequências
    console.print("\n[cyan]Validando sequências...[/cyan]")
    valid_sequences = []
    for seq in sequences:
        valid, error = validate_sequence(seq['sequence'])
        if valid:
            console.print(f"[green]  ✓ {seq['id']}: {seq['length']:,} bp[/green]")
            valid_sequences.append(seq)
        else:
            console.print(f"[yellow]  ⚠ {seq['id']}: {error}[/yellow]")
    
    if not valid_sequences:
        console.print("[red]✗ Nenhuma sequência válida encontrada[/red]")
        sys.exit(1)
    
    # Inicializar AlphaGenome
    console.print("\n[cyan]Inicializando AlphaGenome...[/cyan]")
    analyzer = AlphaGenomeAnalyzer(args.api_key, config)
    if not analyzer.initialize():
        sys.exit(1)
    
    # Mostrar informações de ontologia se configurado
    if config.get('show_ontology_info', True):
        display_ontology_info()
    
    # Processar sequências
    console.print("\n[bold cyan]Processando sequências...[/bold cyan]")
    results = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console
    ) as progress:
        
        task = progress.add_task("Analisando...", total=len(valid_sequences))
        
        for seq in valid_sequences:
            # Análise de variante ou predição normal
            if args.variant:
                pos, ref, alt = args.variant
                result = analyzer.predict_variant(
                    sequence=seq['sequence'],
                    seq_id=seq['id'],
                    variant_position=int(pos),
                    ref_base=ref,
                    alt_base=alt,
                    chromosome=args.chromosome,
                    start=args.start
                )
                
                if result and not args.no_plots:
                    create_variant_visualization(result, args.output, config)
                
                # Salvar metadados se configurado
                if result and config.get('save_metadata', True):
                    save_metadata_to_file(result, args.output)
            else:
                result = analyzer.predict_sequence(
                    sequence=seq['sequence'],
                    seq_id=seq['id'],
                    chromosome=args.chromosome,
                    start=args.start,
                    requested_outputs=config['default_outputs']
                )
                
                if result and not args.no_plots:
                    create_visualizations(result, args.output, config)
                
                # Salvar metadados se configurado
                if result and config.get('save_metadata', True):
                    save_metadata_to_file(result, args.output)
            
            results.append(result)
            progress.update(task, advance=1)
    
    # Gerar relatório
    console.print("\n[cyan]Gerando relatório...[/cyan]")
    generate_report(valid_sequences, results, args.output)
    
    # Resumo
    console.print("\n")
    print_summary(valid_sequences, results)
    
    # Mensagem final
    successful = sum(1 for r in results if r is not None)
    console.print(f"\n[bold green]✓ Análise concluída![/bold green]")
    console.print(f"[green]  {successful}/{len(valid_sequences)} sequências processadas com sucesso[/green]")
    console.print(f"[green]  Resultados salvos em: {args.output}[/green]")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n[yellow]⚠ Análise interrompida pelo usuário[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[red]✗ Erro fatal: {e}[/red]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")
        sys.exit(1)

