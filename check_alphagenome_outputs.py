#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para verificar quais outputs estão disponíveis no AlphaGenome
"""

import sys

try:
    from alphagenome.models import dna_client
    from rich.console import Console
    from rich.table import Table
    
    console = Console()
    
    console.print("\n[bold cyan]Verificando Outputs Disponíveis no AlphaGenome[/bold cyan]\n")
    
    # Listar todos os outputs disponíveis
    table = Table(title="Outputs Disponíveis")
    table.add_column("Nome", style="cyan")
    table.add_column("Valor", style="green")
    
    outputs_found = []
    
    # Tentar obter todos os atributos de OutputType
    for attr in dir(dna_client.OutputType):
        if not attr.startswith('_'):
            try:
                value = getattr(dna_client.OutputType, attr)
                outputs_found.append(attr)
                table.add_row(attr, str(value))
            except:
                pass
    
    console.print(table)
    console.print(f"\n[green]✓ Total de {len(outputs_found)} outputs encontrados[/green]\n")
    
    # Sugestão de uso
    console.print("[bold]Exemplos de uso:[/bold]")
    console.print(f"  python neural_module.py -i input.fasta -k API_KEY -o results/ \\")
    console.print(f"      --outputs {' '.join(outputs_found[:3])}")
    
except ImportError as e:
    print(f"Erro: AlphaGenome não está instalado.")
    print(f"Execute: bash install_alphagenome.sh")
    sys.exit(1)
except Exception as e:
    print(f"Erro: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

