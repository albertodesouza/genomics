#!/bin/bash
# Script para executar a análise de diferenças genéticas em C++

echo "🚀 Iniciando análise de diferenças genéticas (versão C++ ultra-otimizada)"
echo "=================================================================="
echo ""

# Verificar se o executável existe
if [ ! -f "./conta_diferencas_genes" ]; then
    echo "❌ Executável não encontrado. Compilando..."
    make
    if [ $? -ne 0 ]; then
        echo "❌ Falha na compilação!"
        exit 1
    fi
fi

# Verificar se os arquivos FASTA existem
echo "📋 Verificando arquivos FASTA..."
if [ ! -f "fasta/NA12891.genes.consensus.fa" ]; then
    echo "❌ Arquivo fasta/NA12891.genes.consensus.fa não encontrado!"
    exit 1
fi
if [ ! -f "fasta/NA12892.genes.consensus.fa" ]; then
    echo "❌ Arquivo fasta/NA12892.genes.consensus.fa não encontrado!"
    exit 1
fi
if [ ! -f "fasta/NA12878.genes.consensus.fa" ]; then
    echo "❌ Arquivo fasta/NA12878.genes.consensus.fa não encontrado!"
    exit 1
fi
echo "✅ Todos os arquivos FASTA encontrados!"
echo ""

# Criar diretório de saída se não existir
mkdir -p fasta

# Executar análise
echo "⚡ Executando análise C++ ultra-otimizada..."
echo "=================================================================="
time ./conta_diferencas_genes
EXIT_CODE=$?

echo ""
echo "=================================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Análise concluída com sucesso!"
    echo "📄 Resultado salvo em: fasta/family_pairwise_differences.csv"
    
    # Mostrar tamanho do arquivo de saída
    if [ -f "fasta/family_pairwise_differences.csv" ]; then
        LINES=$(wc -l < "fasta/family_pairwise_differences.csv")
        SIZE=$(du -h "fasta/family_pairwise_differences.csv" | cut -f1)
        echo "📊 Arquivo CSV: $LINES linhas, $SIZE"
    fi
else
    echo "❌ Análise falhou com código de saída: $EXIT_CODE"
fi

echo ""
echo "💡 Para comparar com a versão Python:"
echo "   time python3 conta_diferencas_genes.py"
echo ""
echo "🏁 Script finalizado!"
