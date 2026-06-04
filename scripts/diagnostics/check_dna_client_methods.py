#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para verificar métodos disponíveis no DnaClient
"""

import sys

try:
    from alphagenome.models import dna_client
    from rich.console import Console
    from rich.table import Table
    
    console = Console()
    
    console.print("\n[bold cyan]Verificando Métodos do DnaClient[/bold cyan]\n")
    
    # Criar cliente de exemplo
    try:
        client = dna_client.create("test_key")
        
        # Listar métodos públicos
        table = Table(title="Métodos Públicos do DnaClient")
        table.add_column("Método", style="cyan")
        table.add_column("Tipo", style="green")
        
        methods_found = []
        
        for attr in dir(client):
            if not attr.startswith('_') and callable(getattr(client, attr)):
                methods_found.append(attr)
                attr_type = type(getattr(client, attr)).__name__
                table.add_row(attr, attr_type)
        
        console.print(table)
        console.print(f"\n[green]✓ Total de {len(methods_found)} métodos públicos encontrados[/green]\n")
        
    except Exception as e:
        console.print(f"[yellow]Aviso: Não foi possível criar cliente: {e}[/yellow]")
        console.print("[cyan]Listando métodos da classe mesmo assim...[/cyan]\n")
        
        # Listar métodos da classe
        table = Table(title="Métodos da Classe DnaClient")
        table.add_column("Método", style="cyan")
        
        for attr in dir(dna_client.DnaClient):
            if not attr.startswith('_') and attr not in ['create']:
                table.add_row(attr)
        
        console.print(table)
    
except ImportError as e:
    print(f"Erro: AlphaGenome não está instalado.")
    print(f"Execute: bash install_alphagenome.sh")
    sys.exit(1)
except Exception as e:
    print(f"Erro: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

