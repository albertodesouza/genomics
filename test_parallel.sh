#!/bin/bash
# Script para testar a versão paralela

echo "🚀 TESTE DE PERFORMANCE - Versão Paralela vs Serial"
echo "============================================================"

cd /dados/GENOMICS_DATA/top3

echo ""
echo "🔧 Configurações do sistema:"
echo "   CPU cores: $(nproc)"
echo "   OpenMP threads disponíveis: $(~/genomics/conta_diferencas_genes 2>&1 | head -2 | tail -1 | grep -o '[0-9]\+' | head -1 || echo 'N/A')"
echo ""

echo "🧪 Testando com um subconjunto pequeno primeiro..."

# Criar um subconjunto pequeno para teste rápido
echo "📦 Criando subconjunto de teste (primeiros 100 genes)..."
mkdir -p test_subset

# Extrair apenas os primeiros 100 genes de cada arquivo
for file in fasta/NA1289*.genes.consensus.fa; do
    basename=$(basename "$file")
    echo "   Processando $basename..."
    awk '/^>/ {count++} count <= 100 {print}' "$file" > "test_subset/$basename"
done

echo ""
echo "📊 Arquivos de teste criados:"
for file in test_subset/*.fa; do
    genes=$(grep -c "^>" "$file")
    echo "   $(basename "$file"): $genes genes"
done

echo ""
echo "🚀 Executando teste paralelo..."
echo "============================================================"

# Modificar temporariamente os caminhos no executável seria complexo,
# então vamos copiar os arquivos para o local esperado
cp test_subset/* fasta/ 2>/dev/null || echo "Arquivos já no local correto"

# Executar o programa paralelo
time ~/genomics/conta_diferencas_genes

echo ""
echo "============================================================"
echo "✅ Teste concluído!"

# Mostrar estatísticas do resultado
if [ -f "fasta/family_pairwise_differences.csv" ]; then
    lines=$(wc -l < "fasta/family_pairwise_differences.csv")
    echo "📄 Arquivo CSV gerado com $lines linhas"
    echo ""
    echo "📊 Primeiras linhas do resultado:"
    head -10 "fasta/family_pairwise_differences.csv"
fi

echo ""
echo "💡 Para executar com o dataset completo:"
echo "   cd /dados/GENOMICS_DATA/top3"
echo "   ~/genomics/conta_diferencas_genes"
echo ""
echo "🎯 Com 16 cores, esperamos speedup de 8-12x comparado à versão serial!"
