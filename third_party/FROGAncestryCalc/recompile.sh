#!/bin/bash

# Script para recompilar FROGAncestryCalc a partir dos fontes
# Uso: ./recompile.sh

cd "$(dirname "$0")"

echo "=========================================="
echo "Recompilando FROGAncestryCalc"
echo "=========================================="
echo ""

# Remover binários antigos
echo "→ Removendo binários antigos..."
rm -rf bin
mkdir -p bin

# Compilar todos os arquivos Java
echo "→ Compilando código fonte..."
javac -d bin -sourcepath src $(find src -name "*.java") 2>&1 | grep -v "^Note:" | grep -v "warning: \[removal\]"

if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "✓ Compilação concluída com sucesso!"
    
    # Copiar arquivos de dados necessários
    echo "→ Copiando arquivos de dados..."
    cp -r src/read/data bin/read/
    
    echo ""
    echo "=========================================="
    echo "✓ Recompilação completa!"
    echo "=========================================="
    echo ""
    echo "Para executar: ./run.sh"
else
    echo ""
    echo "✗ Erro na compilação!"
    exit 1
fi

