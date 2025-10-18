#!/usr/bin/env bash
# Script para instalação do AlphaGenome

set -e

echo "=========================================="
echo "Instalação do AlphaGenome"
echo "=========================================="
echo ""

# Verificar se conda está disponível
if ! command -v conda &> /dev/null; then
    echo "❌ Erro: conda não encontrado. Por favor, instale o Miniconda ou Anaconda primeiro."
    exit 1
fi

# Verificar se ambiente genomics existe
if conda env list | grep -q "^genomics "; then
    echo "✓ Ambiente conda 'genomics' encontrado"
    ACTIVATE_CMD="conda activate genomics"
else
    echo "⚠ Ambiente 'genomics' não encontrado. Usando ambiente base."
    ACTIVATE_CMD=""
fi

# Diretório temporário
TEMP_DIR=$(mktemp -d)
cd "$TEMP_DIR"

echo ""
echo "📦 Clonando repositório do AlphaGenome..."
git clone https://github.com/google-deepmind/alphagenome.git

echo ""
echo "📦 Instalando AlphaGenome..."
if [ -n "$ACTIVATE_CMD" ]; then
    eval "$ACTIVATE_CMD"
fi

pip install ./alphagenome

echo ""
echo "🧹 Limpando arquivos temporários..."
cd -
rm -rf "$TEMP_DIR"

echo ""
echo "=========================================="
echo "✅ Instalação concluída com sucesso!"
echo "=========================================="
echo ""
echo "Próximos passos:"
echo "1. Obtenha sua API key em: https://www.alphagenomedocs.com/"
echo "2. Execute o neural_module.py:"
echo "   python neural_module.py -i seu_arquivo.fasta -k SUA_API_KEY -o resultados/"
echo ""

