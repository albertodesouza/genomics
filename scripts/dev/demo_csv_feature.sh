#!/bin/bash
# Script para demonstrar a funcionalidade de leitura de CSV

echo "🧪 DEMONSTRAÇÃO - Funcionalidade de Leitura de CSV"
echo "=================================================================="
echo ""

cd /dados/GENOMICS_DATA/top3

echo "📋 Verificando arquivos..."
if [ -f "fasta/family_pairwise_differences.csv" ]; then
    csv_size=$(du -h "fasta/family_pairwise_differences.csv" | cut -f1)
    csv_lines=$(wc -l < "fasta/family_pairwise_differences.csv")
    echo "   ✓ CSV existente: $csv_lines linhas, $csv_size"
else
    echo "   ❌ CSV não encontrado"
    exit 1
fi

echo ""
echo "🚀 Teste 1: Executando com CSV existente (modo rápido)"
echo "=================================================================="
time ~/genomics/native/genes_difference_count/genes_difference_count

echo ""
echo "=================================================================="
echo "💡 EXPLICAÇÃO DA FUNCIONALIDADE:"
echo ""
echo "🔍 O programa agora verifica automaticamente se o arquivo CSV já existe:"
echo "   • Se EXISTE: Lê os dados do CSV e apresenta o resumo (< 1 segundo)"
echo "   • Se NÃO EXISTE: Executa o processamento completo (vários minutos)"
echo ""
echo "✅ VANTAGENS:"
echo "   • Evita reprocessamento desnecessário"
echo "   • Apresenta resultados corretos a partir do CSV"
echo "   • Permite verificação rápida de resultados"
echo "   • Mantém compatibilidade com processamento completo"
echo ""
echo "🎯 COMO USAR:"
echo "   • Para ver resultados: ~/genomics/native/genes_difference_count/genes_difference_count"
echo "   • Para reprocessar: rm fasta/family_pairwise_differences.csv && ~/genomics/native/genes_difference_count/genes_difference_count"
echo ""
echo "🎉 Funcionalidade implementada com sucesso!"
