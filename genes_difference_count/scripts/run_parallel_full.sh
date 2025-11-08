#!/bin/bash
# Script para executar a análise completa com paralelização

echo "🚀 ANÁLISE GENÉTICA PARALELA - Dataset Completo"
echo "=============================================================="

# Verificar se estamos no diretório correto
if [ ! -d "/dados/GENOMICS_DATA/top3/fasta" ]; then
    echo "❌ Erro: Diretório de dados não encontrado!"
    echo "   Execute este script de: /dados/GENOMICS_DATA/top3"
    exit 1
fi

cd /dados/GENOMICS_DATA/top3

echo ""
echo "🔧 Configurações do sistema:"
echo "   CPU cores disponíveis: $(nproc)"
echo "   Memória RAM: $(free -h | grep '^Mem:' | awk '{print $2}')"
echo ""

echo "📋 Verificando arquivos FASTA..."
for file in fasta/NA1289*.genes.consensus.fa; do
    if [ -f "$file" ]; then
        genes=$(grep -c "^>" "$file")
        size=$(du -h "$file" | cut -f1)
        echo "   ✓ $(basename "$file"): $genes genes, $size"
    else
        echo "   ❌ $(basename "$file"): não encontrado"
        exit 1
    fi
done

echo ""
echo "🚀 Iniciando análise paralela com dataset completo..."
echo "   Usando todos os $(($(nproc))) cores disponíveis"
echo "   Estimativa: 15-30x mais rápido que versão Python"
echo "   Estimativa: 8-12x mais rápido que versão C++ serial"
echo ""
echo "=============================================================="

# Executar com monitoramento de recursos
echo "⏱️  Iniciando cronômetro..."
START_TIME=$(date +%s)

# Executar o programa
~/genomics/conta_diferencas_genes

# Calcular tempo total
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "=============================================================="
echo "✅ ANÁLISE PARALELA CONCLUÍDA!"
echo ""
echo "⏱️  Tempo total: ${HOURS}h ${MINUTES}m ${SECONDS}s"

# Verificar resultado
if [ -f "fasta/family_pairwise_differences.csv" ]; then
    lines=$(wc -l < "fasta/family_pairwise_differences.csv")
    size=$(du -h "fasta/family_pairwise_differences.csv" | cut -f1)
    echo "📄 Arquivo CSV: $lines linhas, $size"
    
    echo ""
    echo "📊 Estatísticas do arquivo de saída:"
    echo "   Total de linhas: $lines"
    echo "   Genes únicos: $((($lines - 1) / 3))"
    echo "   Comparações por gene: 3 (pai×mãe, pai×filha, mãe×filha)"
    
    echo ""
    echo "🎯 Primeiras comparações:"
    head -10 "fasta/family_pairwise_differences.csv" | column -t -s ','
    
else
    echo "❌ Erro: Arquivo de saída não foi gerado!"
    exit 1
fi

echo ""
echo "🎉 Análise genética paralela finalizada com sucesso!"
echo "   Arquivo de resultados: fasta/family_pairwise_differences.csv"
