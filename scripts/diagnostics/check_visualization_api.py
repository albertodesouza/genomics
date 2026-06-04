#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para verificar API de visualização do AlphaGenome
"""

import sys

try:
    from alphagenome.visualization import plot_components
    from rich.console import Console
    from rich.table import Table
    
    console = Console()
    
    console.print("\n[bold cyan]Verificando API de Visualização[/bold cyan]\n")
    
    # Listar todos os atributos disponíveis
    table = Table(title="Componentes Disponíveis em plot_components")
    table.add_column("Nome", style="cyan")
    table.add_column("Tipo", style="green")
    
    components_found = []
    
    for attr in dir(plot_components):
        if not attr.startswith('_'):
            try:
                value = getattr(plot_components, attr)
                components_found.append(attr)
                attr_type = type(value).__name__
                table.add_row(attr, attr_type)
            except:
                pass
    
    console.print(table)
    console.print(f"\n[green]✓ Total de {len(components_found)} componentes encontrados[/green]\n")
    
    # Tentar ver funções específicas de plot
    console.print("[bold]Funções de plot disponíveis:[/bold]")
    for attr in components_found:
        obj = getattr(plot_components, attr)
        if callable(obj) and 'plot' in attr.lower():
            console.print(f"  • {attr}")
    
except ImportError as e:
    print(f"Erro: Módulo de visualização não está disponível: {e}")
    sys.exit(1)
except Exception as e:
    print(f"Erro: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

