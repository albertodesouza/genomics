#!/bin/bash

# Script para executar FROGAncestryCalc a partir dos fontes compilados
# Uso: ./run.sh

cd "$(dirname "$0")"

echo "FROG-kb Batch Likelihood Computation Tool (compilado de fontes)"
echo "================================================================"

# Limpar arquivos temporários anteriores
rm -f input/ind/* 
rm -f input/indGenotype/*
rm -f output/*.txt

# Executar com locale inglês para evitar problemas de formatação
LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 java -cp bin main.ComputeBatchAnalysis

# Verificar resultados
if [ -f output/*_likelihood.txt ]; then
    echo ""
    echo "✓ Processamento concluído com sucesso!"
    echo "Arquivos gerados em output/:"
    ls -lh output/*.txt
else
    echo ""
    echo "⚠ Verifique o arquivo log/workingLog.txt para detalhes"
fi

