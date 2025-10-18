#!/usr/bin/env bash
# Script de teste para o neural_module.py

set -e

echo "=========================================="
echo "Teste do Neural Module"
echo "=========================================="
echo ""

# Verificar se neural_module.py existe
if [ ! -f "neural_module.py" ]; then
    echo "❌ Erro: neural_module.py não encontrado"
    exit 1
fi

# Verificar se API key foi fornecida
if [ -z "$1" ]; then
    echo "❌ Erro: API key não fornecida"
    echo ""
    echo "Uso: $0 YOUR_API_KEY"
    echo ""
    echo "Para obter uma API key, acesse:"
    echo "https://www.alphagenomedocs.com/"
    exit 1
fi

API_KEY="$1"

# Criar diretório de testes
TEST_DIR="neural_test_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo "📁 Diretório de teste: $TEST_DIR"
echo ""

# Teste 1: Análise básica
echo "🧪 Teste 1: Análise básica com RNA_SEQ"
python neural_module.py \
    -i example_sequence.fasta \
    -k "$API_KEY" \
    -o "$TEST_DIR/test1_basic" \
    --outputs RNA_SEQ

echo ""
echo "✅ Teste 1 concluído"
echo ""

# Teste 2: Múltiplos outputs
echo "🧪 Teste 2: Análise com múltiplos outputs"
python neural_module.py \
    -i example_sequence.fasta \
    -k "$API_KEY" \
    -o "$TEST_DIR/test2_multi" \
    --outputs RNA_SEQ ATAC H3K27AC

echo ""
echo "✅ Teste 2 concluído"
echo ""

# Teste 3: Múltiplos formatos
echo "🧪 Teste 3: Análise com múltiplos formatos de saída"
python neural_module.py \
    -i example_sequence.fasta \
    -k "$API_KEY" \
    -o "$TEST_DIR/test3_formats" \
    --outputs RNA_SEQ \
    --formats png pdf

echo ""
echo "✅ Teste 3 concluído"
echo ""

# Teste 4: Sem gráficos
echo "🧪 Teste 4: Análise sem geração de gráficos"
python neural_module.py \
    -i example_sequence.fasta \
    -k "$API_KEY" \
    -o "$TEST_DIR/test4_noplot" \
    --outputs RNA_SEQ \
    --no-plots

echo ""
echo "✅ Teste 4 concluído"
echo ""

# Resumo
echo "=========================================="
echo "✅ Todos os testes concluídos com sucesso!"
echo "=========================================="
echo ""
echo "Resultados salvos em: $TEST_DIR/"
echo ""
echo "Estrutura de saída:"
tree "$TEST_DIR/" 2>/dev/null || find "$TEST_DIR/" -type f
echo ""

