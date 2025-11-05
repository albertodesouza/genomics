#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exemplo de Uso Programático do Neural Module

Este script demonstra como usar o neural_module.py como uma biblioteca
Python em seus próprios projetos.
"""

import sys
from pathlib import Path

# Importar componentes do neural_module
from neural_module import (
    AlphaGenomeAnalyzer,
    parse_fasta,
    validate_sequence,
    create_visualizations,
    DEFAULT_CONFIG,
    console
)


def exemplo_analise_basica(api_key: str):
    """
    Exemplo 1: Análise básica de uma sequência
    """
    console.print("\n[bold cyan]═══ Exemplo 1: Análise Básica ═══[/bold cyan]\n")
    
    # Sequência de exemplo (500 bp)
    sequence = "ATCG" * 125  # 500 bp
    
    # Inicializar analisador
    analyzer = AlphaGenomeAnalyzer(api_key)
    if not analyzer.initialize():
        return False
    
    # Fazer predição
    resultado = analyzer.predict_sequence(
        sequence=sequence,
        seq_id="exemplo_1",
        requested_outputs=["RNA_SEQ", "ATAC"]
    )
    
    if resultado:
        console.print("[green]✓ Análise concluída com sucesso![/green]")
        console.print(f"Outputs gerados: {resultado['requested_outputs']}")
        return True
    return False


def exemplo_analise_fasta(api_key: str, fasta_path: str):
    """
    Exemplo 2: Análise de arquivo FASTA
    """
    console.print("\n[bold cyan]═══ Exemplo 2: Análise de FASTA ═══[/bold cyan]\n")
    
    # Parse do FASTA
    sequences = parse_fasta(Path(fasta_path))
    console.print(f"Sequências encontradas: {len(sequences)}")
    
    # Validar sequências
    valid_seqs = []
    for seq in sequences:
        valid, error = validate_sequence(seq['sequence'])
        if valid:
            console.print(f"[green]✓ {seq['id']}: {seq['length']} bp[/green]")
            valid_seqs.append(seq)
        else:
            console.print(f"[yellow]⚠ {seq['id']}: {error}[/yellow]")
    
    if not valid_seqs:
        console.print("[red]Nenhuma sequência válida[/red]")
        return False
    
    # Analisar primeira sequência
    analyzer = AlphaGenomeAnalyzer(api_key)
    if not analyzer.initialize():
        return False
    
    seq = valid_seqs[0]
    resultado = analyzer.predict_sequence(
        sequence=seq['sequence'],
        seq_id=seq['id'],
        requested_outputs=["RNA_SEQ"]
    )
    
    if resultado:
        console.print(f"[green]✓ Análise de '{seq['id']}' concluída![/green]")
        return True
    return False


def exemplo_analise_variante(api_key: str):
    """
    Exemplo 3: Análise de variante
    """
    console.print("\n[bold cyan]═══ Exemplo 3: Análise de Variante ═══[/bold cyan]\n")
    
    # Sequência com variante
    sequence = "ATCG" * 250  # 1000 bp
    
    # Inicializar
    analyzer = AlphaGenomeAnalyzer(api_key)
    if not analyzer.initialize():
        return False
    
    # Analisar variante A>C na posição 500
    console.print("Analisando variante A→C na posição 500...")
    resultado = analyzer.predict_variant(
        sequence=sequence,
        seq_id="exemplo_variante",
        variant_position=500,
        ref_base="A",
        alt_base="C"
    )
    
    if resultado:
        console.print("[green]✓ Análise de variante concluída![/green]")
        console.print(f"Variante: {resultado['variant']}")
        return True
    return False


def exemplo_analise_multipla(api_key: str):
    """
    Exemplo 4: Análise de múltiplas sequências
    """
    console.print("\n[bold cyan]═══ Exemplo 4: Análise Múltipla ═══[/bold cyan]\n")
    
    # Várias sequências
    sequences = [
        {"id": "seq1", "sequence": "ATCG" * 125},
        {"id": "seq2", "sequence": "GCTA" * 125},
        {"id": "seq3", "sequence": "TACG" * 125},
    ]
    
    # Inicializar
    analyzer = AlphaGenomeAnalyzer(api_key)
    if not analyzer.initialize():
        return False
    
    # Analisar todas
    resultados = []
    for seq in sequences:
        console.print(f"\nProcessando {seq['id']}...")
        resultado = analyzer.predict_sequence(
            sequence=seq['sequence'],
            seq_id=seq['id'],
            requested_outputs=["RNA_SEQ"]
        )
        resultados.append(resultado)
    
    # Resumo
    sucesso = sum(1 for r in resultados if r is not None)
    console.print(f"\n[green]✓ {sucesso}/{len(sequences)} sequências analisadas com sucesso[/green]")
    return True


def exemplo_configuracao_customizada(api_key: str):
    """
    Exemplo 5: Configuração customizada
    """
    console.print("\n[bold cyan]═══ Exemplo 5: Configuração Customizada ═══[/bold cyan]\n")
    
    # Configuração customizada
    config = DEFAULT_CONFIG.copy()
    config['default_outputs'] = ['RNA_SEQ', 'CAGE', 'ATAC', 'H3K27AC']
    config['plot_resolution'] = 600  # Alta resolução
    config['output_formats'] = ['png', 'pdf', 'svg']
    
    console.print(f"Outputs configurados: {config['default_outputs']}")
    console.print(f"Resolução: {config['plot_resolution']} DPI")
    console.print(f"Formatos: {config['output_formats']}")
    
    # Inicializar com configuração customizada
    analyzer = AlphaGenomeAnalyzer(api_key, config)
    if not analyzer.initialize():
        return False
    
    # Análise
    sequence = "ATCG" * 125
    resultado = analyzer.predict_sequence(
        sequence=sequence,
        seq_id="exemplo_custom",
        requested_outputs=config['default_outputs']
    )
    
    if resultado:
        console.print("[green]✓ Análise com configuração customizada concluída![/green]")
        return True
    return False


def exemplo_integracao_genomes_analyzer():
    """
    Exemplo 6: Integração com genomes_analyzer.py
    """
    console.print("\n[bold cyan]═══ Exemplo 6: Integração com genomes_analyzer.py ═══[/bold cyan]\n")
    
    console.print("""
Este exemplo demonstra como integrar neural_module com genomes_analyzer:

1. Executar pipeline de análise básica:
   python genomes_analyzer.py config.yaml

2. Extrair sequências de interesse do VCF/BAM

3. Converter para FASTA

4. Analisar com AlphaGenome:
   python neural_module.py -i extracted.fasta -k API_KEY -o neural_results/

5. Correlacionar resultados do neural_module com variantes do pipeline

Fluxo de trabalho recomendado:
┌─────────────────┐
│ genomes_analyzer│
│   Pipeline      │
└────────┬────────┘
         │ VCF + BAM
         ▼
┌─────────────────┐
│ Extrair Regiões │
│   de Interesse  │
└────────┬────────┘
         │ FASTA
         ▼
┌─────────────────┐
│ neural_module   │
│  AlphaGenome    │
└────────┬────────┘
         │ Predições
         ▼
┌─────────────────┐
│   Análise       │
│   Integrada     │
└─────────────────┘
    """)
    
    return True


def exemplo_batch_processing(api_key: str, fasta_dir: str, output_dir: str):
    """
    Exemplo 7: Processamento em lote
    """
    console.print("\n[bold cyan]═══ Exemplo 7: Processamento em Lote ═══[/bold cyan]\n")
    
    from pathlib import Path
    
    # Encontrar todos os arquivos FASTA
    fasta_files = list(Path(fasta_dir).glob("*.fasta")) + list(Path(fasta_dir).glob("*.fa"))
    
    if not fasta_files:
        console.print(f"[yellow]Nenhum arquivo FASTA encontrado em {fasta_dir}[/yellow]")
        return False
    
    console.print(f"Encontrados {len(fasta_files)} arquivos FASTA")
    
    # Inicializar
    analyzer = AlphaGenomeAnalyzer(api_key)
    if not analyzer.initialize():
        return False
    
    # Processar cada arquivo
    for fasta_file in fasta_files:
        console.print(f"\n[cyan]Processando: {fasta_file.name}[/cyan]")
        
        sequences = parse_fasta(fasta_file)
        
        for seq in sequences:
            valid, error = validate_sequence(seq['sequence'])
            if not valid:
                console.print(f"[yellow]  ⚠ {seq['id']}: {error}[/yellow]")
                continue
            
            resultado = analyzer.predict_sequence(
                sequence=seq['sequence'],
                seq_id=f"{fasta_file.stem}_{seq['id']}",
                requested_outputs=["RNA_SEQ"]
            )
            
            if resultado:
                # Salvar visualização
                output_path = Path(output_dir) / fasta_file.stem
                output_path.mkdir(parents=True, exist_ok=True)
                console.print(f"[green]  ✓ {seq['id']} processado[/green]")
    
    console.print(f"\n[green]✓ Processamento em lote concluído![/green]")
    return True


def main():
    """
    Função principal - executa todos os exemplos
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Exemplos de uso programático do Neural Module"
    )
    parser.add_argument('-k', '--api-key', required=True,
                       help='Chave da API do AlphaGenome')
    parser.add_argument('-e', '--exemplo', type=int, choices=range(1, 8),
                       help='Executar exemplo específico (1-7)')
    parser.add_argument('-f', '--fasta',
                       help='Arquivo FASTA para exemplos 2 e 7')
    parser.add_argument('-d', '--directory',
                       help='Diretório com FASTAs para exemplo 7')
    parser.add_argument('-o', '--output', default='output_examples',
                       help='Diretório de saída')
    
    args = parser.parse_args()
    
    console.print("""
[bold cyan]═══════════════════════════════════════════════════════════
    Neural Module - Exemplos de Uso Programático
═══════════════════════════════════════════════════════════[/bold cyan]
    """)
    
    exemplos = {
        1: lambda: exemplo_analise_basica(args.api_key),
        2: lambda: exemplo_analise_fasta(args.api_key, args.fasta or 'example_sequence.fasta'),
        3: lambda: exemplo_analise_variante(args.api_key),
        4: lambda: exemplo_analise_multipla(args.api_key),
        5: lambda: exemplo_configuracao_customizada(args.api_key),
        6: lambda: exemplo_integracao_genomes_analyzer(),
        7: lambda: exemplo_batch_processing(args.api_key, args.directory or '.', args.output),
    }
    
    if args.exemplo:
        # Executar exemplo específico
        exemplos[args.exemplo]()
    else:
        # Executar exemplos 1-6 (não 7 que precisa de diretório)
        for i in range(1, 7):
            try:
                if i == 2 and not Path(args.fasta or 'example_sequence.fasta').exists():
                    console.print(f"[yellow]⚠ Pulando exemplo {i} - FASTA não encontrado[/yellow]")
                    continue
                exemplos[i]()
            except Exception as e:
                console.print(f"[red]✗ Erro no exemplo {i}: {e}[/red]")
    
    console.print("""
[bold green]
═══════════════════════════════════════════════════════════
    Exemplos concluídos!
═══════════════════════════════════════════════════════════[/bold green]

Para mais informações, consulte:
• NEURAL_MODULE_README.md
• python neural_module.py --help
    """)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n[yellow]⚠ Interrompido pelo usuário[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[red]✗ Erro: {e}[/red]")
        sys.exit(1)

