#!/usr/bin/env python3
"""
Script para selecionar 11 genes aleatórios do genoma humano,
excluindo cromossomos X e Y.

Uso:
    python select_random_genes.py [--seed SEED] [--n N_GENES] [--output FILE]

Autor: ChatGPT (para Alberto)
Data: 2025-12-29
"""

import pandas as pd
import random
import argparse
from pathlib import Path

# Configurações padrão
DEFAULT_GTF_CACHE = Path("/dados/GENOMICS_DATA/top3/non_longevous_results_genes/gtf_cache.feather")
DEFAULT_OUTPUT = Path("/home/lume2/genomics/build_non_longevous_dataset/docs/random_11_genes.txt")
DEFAULT_N_GENES = 11
DEFAULT_SEED = 42


def main():
    parser = argparse.ArgumentParser(
        description="Seleciona genes aleatórios do genoma humano (excluindo chrX e chrY)"
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"Seed para reprodutibilidade (default: {DEFAULT_SEED})")
    parser.add_argument("--n", type=int, default=DEFAULT_N_GENES,
                        help=f"Número de genes a selecionar (default: {DEFAULT_N_GENES})")
    parser.add_argument("--output", type=str, default=str(DEFAULT_OUTPUT),
                        help=f"Arquivo de saída (default: {DEFAULT_OUTPUT})")
    parser.add_argument("--gtf", type=str, default=str(DEFAULT_GTF_CACHE),
                        help=f"Arquivo GTF cache (default: {DEFAULT_GTF_CACHE})")
    parser.add_argument("--all-types", action="store_true",
                        help="Incluir todos os tipos de genes (não apenas protein_coding)")
    
    args = parser.parse_args()
    
    gtf_path = Path(args.gtf)
    output_path = Path(args.output)
    n_genes = args.n
    seed = args.seed
    
    # Carregar GTF cache
    print(f"[INFO] Carregando GTF de: {gtf_path}")
    
    if not gtf_path.exists():
        print(f"[ERROR] Arquivo GTF não encontrado: {gtf_path}")
        print("[INFO] Tentando baixar do AlphaGenome...")
        gtf = pd.read_feather(
            'https://storage.googleapis.com/alphagenome/reference/gencode/'
            'hg38/gencode.v46.annotation.gtf.gz.feather'
        )
        # Salvar cache para uso futuro
        gtf_path.parent.mkdir(parents=True, exist_ok=True)
        gtf.to_feather(gtf_path)
        print(f"[INFO] Cache salvo em: {gtf_path}")
    else:
        gtf = pd.read_feather(gtf_path)
    
    print(f"[INFO] Colunas disponíveis: {gtf.columns.tolist()}")
    print(f"[INFO] Total de linhas no GTF: {len(gtf)}")
    
    # Detectar nomes das colunas (podem variar entre versões do GTF)
    feature_col = 'Feature' if 'Feature' in gtf.columns else 'feature'
    seqname_col = 'Chromosome' if 'Chromosome' in gtf.columns else 'seqname'
    start_col = 'Start' if 'Start' in gtf.columns else 'start'
    end_col = 'End' if 'End' in gtf.columns else 'end'
    
    # Filtrar apenas linhas do tipo 'gene'
    genes_df = gtf[gtf[feature_col] == 'gene'].copy()
    print(f"[INFO] Total de genes: {len(genes_df)}")
    
    # Excluir cromossomos X e Y
    genes_df = genes_df[~genes_df[seqname_col].isin(['chrX', 'chrY', 'X', 'Y'])]
    print(f"[INFO] Genes após excluir X e Y: {len(genes_df)}")
    
    # Filtrar apenas genes protein_coding (a menos que --all-types)
    if not args.all_types:
        if 'gene_type' in genes_df.columns:
            genes_df = genes_df[genes_df['gene_type'] == 'protein_coding']
            print(f"[INFO] Genes protein_coding: {len(genes_df)}")
        else:
            print("[WARN] Coluna 'gene_type' não encontrada, usando todos os genes")
    
    # Garantir que temos gene_name
    if 'gene_name' not in genes_df.columns:
        print(f"[ERROR] Coluna 'gene_name' não encontrada!")
        print(f"[ERROR] Colunas disponíveis: {genes_df.columns.tolist()}")
        return 1
    
    # Remover genes sem nome
    genes_df = genes_df[genes_df['gene_name'].notna() & (genes_df['gene_name'] != '')]
    print(f"[INFO] Genes com nome válido: {len(genes_df)}")
    
    # Selecionar genes únicos
    unique_genes = genes_df['gene_name'].unique()
    print(f"[INFO] Genes únicos: {len(unique_genes)}")
    
    if len(unique_genes) < n_genes:
        print(f"[ERROR] Não há genes suficientes! Disponíveis: {len(unique_genes)}, Solicitados: {n_genes}")
        return 1
    
    # Selecionar N genes aleatórios
    random.seed(seed)
    selected_genes = random.sample(list(unique_genes), n_genes)
    
    print(f"\n{'='*70}")
    print(f"GENES SELECIONADOS ({n_genes}) - Seed: {seed}")
    print(f"{'='*70}")
    
    # Obter informações adicionais para cada gene
    gene_details = []
    for gene in selected_genes:
        gene_info = genes_df[genes_df['gene_name'] == gene].iloc[0]
        chrom = gene_info[seqname_col]
        start = int(gene_info.get(start_col, 0))
        end = int(gene_info.get(end_col, 0))
        length = end - start
        gene_type = gene_info.get('gene_type', 'unknown')
        
        gene_details.append({
            'name': gene,
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': length,
            'type': gene_type
        })
        
        print(f"  {gene:15s} - {chrom}:{start:,}-{end:,} ({length:,} bp)")
    
    # Salvar arquivo
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        for g in gene_details:
            f.write(f"{g['name']} # {g['chrom']}:{g['start']}-{g['end']}, {g['length']} bp, {g['type']}\n")
    
    print(f"\n{'='*70}")
    print(f"✓ Lista salva em: {output_path}")
    print(f"{'='*70}")
    print(f"\nPara usar no config, edite gene_1000.yaml:")
    print(f"  gene_list_file: ../docs/random_11_genes.txt")
    
    return 0


if __name__ == "__main__":
    exit(main())

