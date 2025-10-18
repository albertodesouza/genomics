#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para verificar estrutura dos outputs do AlphaGenome
"""

import sys

try:
    from alphagenome.models import dna_client
    from alphagenome.data import genome
    from rich.console import Console
    from rich.table import Table
    import pprint
    
    console = Console()
    
    API_KEY = sys.argv[1] if len(sys.argv) > 1 else "test"
    
    console.print("\n[bold cyan]Verificando Estrutura dos Outputs[/bold cyan]\n")
    
    # Criar cliente
    client = dna_client.create(API_KEY)
    
    # Criar uma sequência pequena de teste (2048 bp)
    seq = "ATCG" * 512  # 2048 bp
    
    # Criar intervalo
    interval = genome.Interval(chromosome='chr1', start=1000000, end=1002048)
    
    # Ontology terms
    ontology_terms = ['UBERON:0000955']  # cérebro
    
    # Fazer predição
    console.print("[cyan]Fazendo predição de teste...[/cyan]")
    outputs = client.predict_interval(
        interval=interval,
        ontology_terms=ontology_terms,
        requested_outputs=[dna_client.OutputType.RNA_SEQ]
    )
    
    console.print("[green]✓ Predição concluída[/green]\n")
    
    # Analisar estrutura
    console.print("[bold]Tipo do objeto outputs:[/bold]")
    console.print(f"  {type(outputs)}\n")
    
    console.print("[bold]Atributos disponíveis:[/bold]")
    for attr in dir(outputs):
        if not attr.startswith('_'):
            console.print(f"  • {attr}")
    
    console.print("\n[bold]Tentando acessar RNA_SEQ:[/bold]")
    try:
        rna_seq_data = outputs.rna_seq
        console.print(f"  Tipo: {type(rna_seq_data)}")
        console.print(f"  Atributos: {[a for a in dir(rna_seq_data) if not a.startswith('_')][:10]}")
        
        # Tentar acessar interval
        if hasattr(rna_seq_data, 'interval'):
            console.print(f"  Tem intervalo: {rna_seq_data.interval}")
        
        # Tentar ver dados
        if hasattr(rna_seq_data, 'data'):
            console.print(f"  Tem data: {type(rna_seq_data.data)}")
        
        # Tentar ver shape
        if hasattr(rna_seq_data, 'shape'):
            console.print(f"  Shape: {rna_seq_data.shape}")
            
    except Exception as e:
        console.print(f"  [red]Erro: {e}[/red]")
    
except Exception as e:
    print(f"Erro: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

