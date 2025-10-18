#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Neural Integration - Ponte entre genomes_analyzer.py e neural_module.py

Este módulo facilita a integração entre o pipeline de análise genômica
(genomes_analyzer.py) e a análise neural (neural_module.py usando AlphaGenome).

Funcionalidades:
- Extração de regiões de interesse do VCF
- Conversão para FASTA
- Análise com AlphaGenome
- Correlação de resultados
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import subprocess as sp
import json

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
from rich.table import Table
from rich.panel import Panel

console = Console()


# ====================== Extração de Sequências ======================

def extract_sequences_from_vcf(vcf_path: Path, 
                              ref_fasta: Path,
                              region: Optional[str] = None,
                              output_fasta: Path = None) -> Path:
    """
    Extrai sequências de referência para regiões de interesse do VCF.
    
    Args:
        vcf_path: Caminho para arquivo VCF
        ref_fasta: Referência genômica FASTA
        region: Região específica (chr:start-end) ou None para todas
        output_fasta: Arquivo de saída (auto se None)
        
    Returns:
        Path para arquivo FASTA gerado
    """
    console.print(f"[cyan]Extraindo sequências do VCF: {vcf_path.name}[/cyan]")
    
    if output_fasta is None:
        output_fasta = vcf_path.parent / f"{vcf_path.stem}_sequences.fasta"
    
    try:
        # Usar bcftools para extrair regiões
        cmd = ['bcftools', 'query', '-f', '%CHROM\t%POS\n', str(vcf_path)]
        
        if region:
            cmd.extend(['-r', region])
        
        result = sp.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse posições
        positions = []
        for line in result.stdout.strip().split('\n'):
            if line:
                chrom, pos = line.split('\t')
                positions.append((chrom, int(pos)))
        
        console.print(f"[green]✓ {len(positions)} posições encontradas[/green]")
        
        # Extrair sequências com samtools faidx
        with open(output_fasta, 'w') as fout:
            for i, (chrom, pos) in enumerate(positions):
                # Região de ±5kb ao redor da variante
                start = max(1, pos - 5000)
                end = pos + 5000
                region_str = f"{chrom}:{start}-{end}"
                
                # Extrair com samtools
                cmd = ['samtools', 'faidx', str(ref_fasta), region_str]
                result = sp.run(cmd, capture_output=True, text=True, check=True)
                
                # Escrever no FASTA
                header = f">variant_{i+1}_{chrom}_{pos}\n"
                fout.write(header)
                # Pular header do samtools output
                seq_lines = result.stdout.split('\n')[1:]
                fout.write(''.join(seq_lines) + '\n')
        
        console.print(f"[green]✓ Sequências salvas em: {output_fasta}[/green]")
        return output_fasta
        
    except Exception as e:
        console.print(f"[red]✗ Erro ao extrair sequências: {e}[/red]")
        return None


def extract_sequences_from_bed(bed_path: Path,
                               ref_fasta: Path,
                               output_fasta: Path = None) -> Path:
    """
    Extrai sequências para regiões especificadas em arquivo BED.
    
    Args:
        bed_path: Arquivo BED com regiões
        ref_fasta: Referência genômica FASTA
        output_fasta: Arquivo de saída
        
    Returns:
        Path para arquivo FASTA gerado
    """
    console.print(f"[cyan]Extraindo sequências do BED: {bed_path.name}[/cyan]")
    
    if output_fasta is None:
        output_fasta = bed_path.parent / f"{bed_path.stem}_sequences.fasta"
    
    try:
        # Usar bedtools getfasta
        cmd = [
            'bedtools', 'getfasta',
            '-fi', str(ref_fasta),
            '-bed', str(bed_path),
            '-fo', str(output_fasta),
            '-name'
        ]
        
        sp.run(cmd, check=True, capture_output=True)
        
        console.print(f"[green]✓ Sequências extraídas: {output_fasta}[/green]")
        return output_fasta
        
    except Exception as e:
        console.print(f"[red]✗ Erro ao extrair sequências: {e}[/red]")
        return None


def extract_gene_sequences(gene_list: List[str],
                          gtf_path: Path,
                          ref_fasta: Path,
                          output_fasta: Path,
                          flank: int = 10000) -> Path:
    """
    Extrai sequências de genes específicos.
    
    Args:
        gene_list: Lista de nomes de genes
        gtf_path: Arquivo GTF com anotações
        ref_fasta: Referência genômica
        output_fasta: Arquivo de saída
        flank: Bases flanqueadoras (padrão: 10kb)
        
    Returns:
        Path para arquivo FASTA gerado
    """
    console.print(f"[cyan]Extraindo sequências de {len(gene_list)} gene(s)[/cyan]")
    
    try:
        # Parse GTF para encontrar genes
        import gzip
        
        gene_regions = {}
        
        open_func = gzip.open if gtf_path.suffix == '.gz' else open
        
        with open_func(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                if fields[2] != 'gene':
                    continue
                
                # Parse atributos
                attrs = dict(item.strip().split(' ', 1) 
                           for item in fields[8].split(';') 
                           if item.strip())
                
                gene_name = attrs.get('gene_name', '').strip('"')
                
                if gene_name in gene_list:
                    chrom = fields[0]
                    start = max(1, int(fields[3]) - flank)
                    end = int(fields[4]) + flank
                    
                    gene_regions[gene_name] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end
                    }
        
        console.print(f"[green]✓ {len(gene_regions)} gene(s) encontrado(s) no GTF[/green]")
        
        # Extrair sequências
        with open(output_fasta, 'w') as fout:
            for gene_name, region in gene_regions.items():
                region_str = f"{region['chrom']}:{region['start']}-{region['end']}"
                
                cmd = ['samtools', 'faidx', str(ref_fasta), region_str]
                result = sp.run(cmd, capture_output=True, text=True, check=True)
                
                # Escrever com nome do gene
                header = f">{gene_name}_{region['chrom']}_{region['start']}_{region['end']}\n"
                fout.write(header)
                seq_lines = result.stdout.split('\n')[1:]
                fout.write(''.join(seq_lines) + '\n')
        
        console.print(f"[green]✓ Sequências salvas em: {output_fasta}[/green]")
        return output_fasta
        
    except Exception as e:
        console.print(f"[red]✗ Erro ao extrair genes: {e}[/red]")
        return None


# ====================== Análise Integrada ======================

def run_integrated_analysis(vcf_path: Path,
                           ref_fasta: Path,
                           api_key: str,
                           output_dir: Path,
                           outputs: List[str] = None):
    """
    Executa análise integrada: VCF → FASTA → AlphaGenome.
    
    Args:
        vcf_path: Arquivo VCF de entrada
        ref_fasta: Referência genômica
        api_key: API key do AlphaGenome
        output_dir: Diretório de saída
        outputs: Tipos de output do AlphaGenome
    """
    console.print(Panel.fit(
        "[bold cyan]Análise Integrada: Pipeline + Neural[/bold cyan]",
        border_style="cyan"
    ))
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Passo 1: Extrair sequências
    console.print("\n[bold cyan]Passo 1: Extração de Sequências[/bold cyan]")
    fasta_path = extract_sequences_from_vcf(vcf_path, ref_fasta)
    
    if not fasta_path or not fasta_path.exists():
        console.print("[red]✗ Falha na extração de sequências[/red]")
        return False
    
    # Passo 2: Análise com AlphaGenome
    console.print("\n[bold cyan]Passo 2: Análise Neural com AlphaGenome[/bold cyan]")
    
    neural_output = output_dir / "neural_results"
    neural_output.mkdir(exist_ok=True)
    
    cmd = [
        'python', 'neural_module.py',
        '-i', str(fasta_path),
        '-k', api_key,
        '-o', str(neural_output)
    ]
    
    if outputs:
        cmd.extend(['--outputs'] + outputs)
    
    try:
        sp.run(cmd, check=True)
        console.print("[green]✓ Análise neural concluída[/green]")
    except Exception as e:
        console.print(f"[red]✗ Erro na análise neural: {e}[/red]")
        return False
    
    # Passo 3: Correlacionar resultados
    console.print("\n[bold cyan]Passo 3: Correlação de Resultados[/bold cyan]")
    correlate_results(vcf_path, neural_output, output_dir)
    
    console.print(f"\n[bold green]✓ Análise integrada concluída![/bold green]")
    console.print(f"[green]Resultados em: {output_dir}[/green]")
    
    return True


def correlate_results(vcf_path: Path, neural_dir: Path, output_dir: Path):
    """
    Correlaciona resultados do VCF com predições do AlphaGenome.
    
    Args:
        vcf_path: Arquivo VCF original
        neural_dir: Diretório com resultados do neural_module
        output_dir: Diretório para salvar correlação
    """
    console.print("[cyan]Correlacionando resultados...[/cyan]")
    
    # Carregar relatório neural
    report_path = neural_dir / "analysis_report.json"
    
    if not report_path.exists():
        console.print("[yellow]⚠ Relatório neural não encontrado[/yellow]")
        return
    
    with open(report_path) as f:
        neural_report = json.load(f)
    
    # Criar tabela de correlação
    correlation_table = Table(title="Correlação: Variantes × Predições Neurais")
    correlation_table.add_column("Sequência", style="cyan")
    correlation_table.add_column("Status", style="green")
    correlation_table.add_column("Outputs", style="yellow")
    
    for seq in neural_report['sequences']:
        status = "✓" if seq['status'] == 'success' else "✗"
        outputs = ", ".join(seq.get('outputs', []))
        correlation_table.add_row(seq['id'], status, outputs)
    
    console.print(correlation_table)
    
    # Salvar correlação
    correlation_file = output_dir / "correlation_report.json"
    with open(correlation_file, 'w') as f:
        json.dump({
            'vcf_source': str(vcf_path),
            'neural_results': neural_report,
            'summary': {
                'total_sequences': neural_report['total_sequences'],
                'successful_predictions': neural_report['successful_analyses']
            }
        }, f, indent=2)
    
    console.print(f"[green]✓ Correlação salva em: {correlation_file}[/green]")


# ====================== Main ======================

def main():
    parser = argparse.ArgumentParser(
        description='Integração entre genomes_analyzer e neural_module',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:

1. Análise integrada completa:
   %(prog)s --vcf variants.vcf --ref genome.fa --api-key KEY --output results/

2. Extrair sequências de genes:
   %(prog)s --genes BRCA1 TP53 --gtf genes.gtf --ref genome.fa --output genes.fasta

3. Extrair de arquivo BED:
   %(prog)s --bed regions.bed --ref genome.fa --output regions.fasta
        """
    )
    
    # Modos de operação
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--integrated', action='store_true',
                     help='Análise integrada completa')
    mode.add_argument('--extract-vcf', action='store_true',
                     help='Apenas extrair sequências do VCF')
    mode.add_argument('--extract-bed', action='store_true',
                     help='Apenas extrair sequências do BED')
    mode.add_argument('--extract-genes', action='store_true',
                     help='Apenas extrair sequências de genes')
    
    # Arquivos de entrada
    parser.add_argument('--vcf', type=Path,
                       help='Arquivo VCF')
    parser.add_argument('--bed', type=Path,
                       help='Arquivo BED')
    parser.add_argument('--ref', type=Path, required=True,
                       help='Referência genômica FASTA')
    parser.add_argument('--gtf', type=Path,
                       help='Arquivo GTF (para extração de genes)')
    parser.add_argument('--genes', nargs='+',
                       help='Lista de genes para extrair')
    
    # Parâmetros
    parser.add_argument('--api-key',
                       help='API key do AlphaGenome (para análise integrada)')
    parser.add_argument('--output', type=Path, required=True,
                       help='Arquivo/diretório de saída')
    parser.add_argument('--outputs', nargs='+',
                       help='Outputs do AlphaGenome')
    parser.add_argument('--flank', type=int, default=10000,
                       help='Bases flanqueadoras para genes (padrão: 10000)')
    
    args = parser.parse_args()
    
    # Banner
    console.print(Panel.fit(
        "[bold cyan]Neural Integration[/bold cyan]\n"
        "[dim]Ponte entre genomes_analyzer.py e neural_module.py[/dim]",
        border_style="cyan"
    ))
    
    # Executar modo apropriado
    if args.integrated:
        if not args.vcf or not args.api_key:
            console.print("[red]✗ Modo integrado requer --vcf e --api-key[/red]")
            sys.exit(1)
        
        run_integrated_analysis(
            args.vcf,
            args.ref,
            args.api_key,
            args.output,
            args.outputs
        )
    
    elif args.extract_vcf:
        if not args.vcf:
            console.print("[red]✗ Modo extract-vcf requer --vcf[/red]")
            sys.exit(1)
        
        extract_sequences_from_vcf(args.vcf, args.ref, output_fasta=args.output)
    
    elif args.extract_bed:
        if not args.bed:
            console.print("[red]✗ Modo extract-bed requer --bed[/red]")
            sys.exit(1)
        
        extract_sequences_from_bed(args.bed, args.ref, output_fasta=args.output)
    
    elif args.extract_genes:
        if not args.genes or not args.gtf:
            console.print("[red]✗ Modo extract-genes requer --genes e --gtf[/red]")
            sys.exit(1)
        
        extract_gene_sequences(
            args.genes,
            args.gtf,
            args.ref,
            args.output,
            args.flank
        )


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n[yellow]⚠ Interrompido pelo usuário[/yellow]")
        sys.exit(130)
    except Exception as e:
        console.print(f"\n[red]✗ Erro: {e}[/red]")
        import traceback
        console.print(f"[dim]{traceback.format_exc()}[/dim]")
        sys.exit(1)

