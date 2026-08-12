#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, csv, time
from collections import defaultdict

# ---------- Config ----------
IN_PAI   = "fasta/NA12891.genes.consensus.fa"
IN_MAE   = "fasta/NA12892.genes.consensus.fa"
IN_FILHA = "fasta/NA12878.genes.consensus.fa"
OUT_CSV  = "fasta/family_pairwise_differences.csv"
MAX_ALIGN_LEN = 200_000  # acima disso, usa caminho rápido (aprox.)
GAP_PENALTY = -1
MATCH_SCORE = 1
MISMATCH_SCORE = -1

# ---------- IUPAC ----------
IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'}, 'W': {'A','T'},
    'K': {'G','T'}, 'M': {'A','C'},
    'B': {'C','G','T'}, 'D': {'A','G','T'}, 'H': {'A','C','T'}, 'V': {'A','C','G'},
    'N': {'A','C','G','T'}
}
def compat(a,b):
    a=a.upper(); b=b.upper()
    if a == '-' or b == '-':
        return False
    if a not in IUPAC or b not in IUPAC:
        return a == b
    return len(IUPAC[a] & IUPAC[b]) > 0

# ---------- FASTA parser (retorna {gene: seq}) ----------
def read_fasta_genes(path):
    d = {}
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    d[name] = ''.join(buf).upper()
                hdr = line[1:].strip()
                # Seus headers estão no formato: Sample|GENE|chr:start-end(strand)
                # então pegamos a 2a parte (índice 1).
                parts = hdr.split('|')
                gene = parts[1] if len(parts) >= 2 else parts[0]
                name = gene
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            d[name] = ''.join(buf).upper()
    return d

# ---------- Alinhamento (Hirschberg: memória linear, tempo O(n*m)) ----------
def score(a,b):
    if a == '-' or b == '-':
        return GAP_PENALTY
    return MATCH_SCORE if compat(a,b) else MISMATCH_SCORE

def nw_last_row(a, b):
    # DP de uma linha (para Hirschberg)
    prev = [j * GAP_PENALTY for j in range(len(b)+1)]
    for i in range(1, len(a)+1):
        curr = [i * GAP_PENALTY] + [0]*len(b)
        ai = a[i-1]
        for j in range(1, len(b)+1):
            bj = b[j-1]
            curr[j] = max(
                prev[j] + GAP_PENALTY,           # gap em a
                curr[j-1] + GAP_PENALTY,         # gap em b
                prev[j-1] + (MATCH_SCORE if compat(ai,bj) else MISMATCH_SCORE)
            )
        prev = curr
    return prev

def hirschberg(a, b):
    # retorna alinhamento (a_aln, b_aln)
    if len(a) == 0:
        return ('-'*len(b), b)
    if len(b) == 0:
        return (a, '-'*len(a))
    if len(a) == 1 or len(b) == 1:
        # Needleman-Wunsch completo em pequeno
        n, m = len(a), len(b)
        dp = [[0]*(m+1) for _ in range(n+1)]
        for i in range(n+1): dp[i][0] = i*GAP_PENALTY
        for j in range(m+1): dp[0][j] = j*GAP_PENALTY
        for i in range(1,n+1):
            for j in range(1,m+1):
                dp[i][j] = max(
                    dp[i-1][j] + GAP_PENALTY,
                    dp[i][j-1] + GAP_PENALTY,
                    dp[i-1][j-1] + (MATCH_SCORE if compat(a[i-1], b[j-1]) else MISMATCH_SCORE)
                )
        # traceback
        i, j = n, m
        aa, bb = [], []
        while i>0 or j>0:
            if i>0 and dp[i][j] == dp[i-1][j] + GAP_PENALTY:
                aa.append(a[i-1]); bb.append('-'); i-=1
            elif j>0 and dp[i][j] == dp[i][j-1] + GAP_PENALTY:
                aa.append('-'); bb.append(b[j-1]); j-=1
            else:
                aa.append(a[i-1] if i>0 else '-')
                bb.append(b[j-1] if j>0 else '-')
                i-=1; j-=1
        return (''.join(reversed(aa)), ''.join(reversed(bb)))
    # caso geral
    mid = len(a)//2
    scoreL = nw_last_row(a[:mid], b)
    scoreR = nw_last_row(a[mid:][::-1], b[::-1])
    # escolhe a coluna de corte melhor
    best, cut = None, 0
    for j in range(len(b)+1):
        val = scoreL[j] + scoreR[len(b)-j]
        if best is None or val > best:
            best, cut = val, j
    left = hirschberg(a[:mid], b[:cut])
    right = hirschberg(a[mid:], b[cut:])
    return (left[0]+right[0], left[1]+right[1])

# ---------- Contador a partir de alinhamento ----------
def count_diffs_from_alignment(a_aln, b_aln):
    assert len(a_aln) == len(b_aln)
    match = diff_sub = diff_gap = ambig_N = 0
    for x, y in zip(a_aln, b_aln):
        if x == '-' or y == '-':
            diff_gap += 1
            continue
        if x == 'N' or y == 'N':
            ambig_N += 1
            continue
        if compat(x,y):
            match += 1
        else:
            diff_sub += 1
    total_cols = len(a_aln)
    informative = total_cols - ambig_N
    return {
        "aligned_cols": total_cols,
        "matches": match,
        "subs": diff_sub,
        "indels": diff_gap,
        "ambiguous_N": ambig_N,
        "differences_total": diff_sub + diff_gap,
        "informative_cols": informative
    }

# ---------- Caminho rápido (aprox, sem alinhar) ----------
def fast_diff_estimate(a, b):
    n = min(len(a), len(b))
    match = diff_sub = ambig_N = 0
    for i in range(n):
        x, y = a[i], b[i]
        if x == 'N' or y == 'N':
            ambig_N += 1
            continue
        if compat(x,y):
            match += 1
        else:
            diff_sub += 1
    # diferença de tamanho ~ indels
    diff_gap = abs(len(a) - len(b))
    total_cols = n + diff_gap
    return {
        "aligned_cols": total_cols,
        "matches": match,
        "subs": diff_sub,
        "indels": diff_gap,
        "ambiguous_N": ambig_N,
        "differences_total": diff_sub + diff_gap,
        "informative_cols": n - ambig_N
    }

# ---------- Pairwise por gene ----------
def compare_sequences(a, b):
    if len(a) == len(b) and len(a) <= MAX_ALIGN_LEN:
        a_aln, b_aln = hirschberg(a, b)
        return count_diffs_from_alignment(a_aln, b_aln)
    # se curtos mas com tamanhos diferentes, ainda alinha (se não muito longos)
    if max(len(a), len(b)) <= MAX_ALIGN_LEN:
        a_aln, b_aln = hirschberg(a, b)
        return count_diffs_from_alignment(a_aln, b_aln)
    # fallback rápido
    return fast_diff_estimate(a, b)

def main():
    print("🧬 Iniciando análise de diferenças genéticas...")
    print(f"📁 Carregando arquivos FASTA...")
    
    start_time = time.time()
    pai = read_fasta_genes(IN_PAI)
    print(f"   ✓ Pai: {len(pai)} genes carregados de {IN_PAI}")
    
    mae = read_fasta_genes(IN_MAE)
    print(f"   ✓ Mãe: {len(mae)} genes carregados de {IN_MAE}")
    
    filha = read_fasta_genes(IN_FILHA)
    print(f"   ✓ Filha: {len(filha)} genes carregados de {IN_FILHA}")

    genes_comuns = sorted(set(pai) & set(mae) & set(filha))
    if not genes_comuns:
        print("[ERRO] Nenhum gene em comum entre os três FASTAs.")
        sys.exit(2)
    
    print(f"\n📊 Encontrados {len(genes_comuns)} genes em comum para análise")
    print(f"🔬 Cada gene será comparado em 3 pares: pai×mãe, pai×filha, mãe×filha")
    print(f"⏱️  Total de comparações: {len(genes_comuns) * 3}")
    print("─" * 60)

    totals = {
        ("pai","mae"): defaultdict(int),
        ("pai","filha"): defaultdict(int),
        ("mae","filha"): defaultdict(int),
    }

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    total_comparisons = len(genes_comuns) * 3
    comparison_count = 0
    
    with open(OUT_CSV, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["gene","pair","aligned_cols","informative_cols","matches","subs","indels","ambiguous_N","differences_total"])
        
        for gene_idx, g in enumerate(genes_comuns, 1):
            print(f"\n🧬 Gene {gene_idx}/{len(genes_comuns)}: {g}")
            p = pai[g]; m = mae[g]; f = filha[g]
            
            for (lab, A, B) in [
                (("pai","mae"), p, m),
                (("pai","filha"), p, f),
                (("mae","filha"), m, f),
            ]:
                comparison_count += 1
                pair_name = f"{lab[0]} × {lab[1]}"
                
                # Informação sobre o tamanho das sequências
                seq_info = f"({len(A)} vs {len(B)} bases)"
                method = "alinhamento completo" if max(len(A), len(B)) <= MAX_ALIGN_LEN else "estimativa rápida"
                
                print(f"   🔬 {pair_name} {seq_info} - usando {method}")
                
                # Calcular estatísticas
                comp_start = time.time()
                stats = compare_sequences(A, B)
                comp_time = time.time() - comp_start
                
                # Mostrar resultado da comparação
                diff_total = stats["differences_total"]
                diff_percent = (diff_total / stats["informative_cols"] * 100) if stats["informative_cols"] > 0 else 0
                print(f"      ➤ {diff_total:,} diferenças ({diff_percent:.2f}%) em {comp_time:.3f}s")
                
                # Progresso geral
                progress_percent = (comparison_count / total_comparisons) * 100
                elapsed_time = time.time() - start_time
                avg_time_per_comp = elapsed_time / comparison_count
                remaining_comps = total_comparisons - comparison_count
                eta_seconds = remaining_comps * avg_time_per_comp
                
                print(f"   📊 Progresso: {comparison_count}/{total_comparisons} ({progress_percent:.1f}%) - ETA: {eta_seconds:.0f}s")
                
                w.writerow([g, f"{lab[0]}_vs_{lab[1]}",
                            stats["aligned_cols"], stats["informative_cols"],
                            stats["matches"], stats["subs"], stats["indels"],
                            stats["ambiguous_N"], stats["differences_total"]])
                for k, v in stats.items():
                    totals[lab][k] += v

    # Calcular tempo total
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"✅ ANÁLISE CONCLUÍDA!")
    print(f"⏱️  Tempo total: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"📊 {len(genes_comuns)} genes analisados em {total_comparisons} comparações")
    print(f"⚡ Velocidade média: {total_comparisons/total_time:.1f} comparações/segundo")
    print("="*60)
    
    # imprime sumário bonito
    print("\n🧬 SUMÁRIO FINAL (tudo somado sobre os genes em comum)")
    print("─" * 60)
    for lab in [("pai","mae"), ("pai","filha"), ("mae","filha")]:
        t = totals[lab]
        diff_percent = (t['differences_total'] / t['informative_cols'] * 100) if t['informative_cols'] > 0 else 0
        
        print(f"\n👥 {lab[0].upper()} × {lab[1].upper()}:")
        print(f"   📏 Bases alinhadas (colunas):  {t['aligned_cols']:,}")
        print(f"   🔍 Informativas (sem 'N'):     {t['informative_cols']:,}")
        print(f"   ✅ Iguais (compatíveis):       {t['matches']:,}")
        print(f"   🔄 Diferentes (substituições): {t['subs']:,}")
        print(f"   📝 Diferentes (indels):        {t['indels']:,}")
        print(f"   ❓ Ambíguas ignoradas (N):     {t['ambiguous_N']:,}")
        print(f"   🎯 TOTAL diferenças:           {t['differences_total']:,} ({diff_percent:.2f}%)")

    print(f"\n📄 Resultado detalhado salvo em: {OUT_CSV}")
    print("🎉 Análise de diferenças genéticas finalizada com sucesso!")

if __name__ == "__main__":
    main()
