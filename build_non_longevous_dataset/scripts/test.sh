#!/bin/bash
# Script de teste para build_non_longevous_dataset.py
# 
# Este script demonstra como usar o programa em modo de teste
# Execute a partir do diretório scripts/

set -e  # Exit on error

# Mover para o diretório do módulo (pai de scripts/)
cd "$(dirname "$0")/.."

echo "========================================================================"
echo "TESTE: Non-Longevous Dataset Builder"
echo "========================================================================"
echo ""
echo "[INFO] Diretório de trabalho: $(pwd)"
echo ""

# Verificar se o arquivo de configuração existe
if [ ! -f "configs/default.yaml" ]; then
    echo "[ERROR] Arquivo configs/default.yaml não encontrado!"
    exit 1
fi

echo "[INFO] Arquivos encontrados. Iniciando teste..."
echo ""

# PASSO 1: Análise de Metadados
echo "========================================================================"
echo "PASSO 1: Analisando metadados do CSV..."
echo "========================================================================"
echo ""

python3 build_non_longevous_dataset.py --config configs/default.yaml

echo ""
echo "========================================================================"
echo "TESTE CONCLUÍDO!"
echo "========================================================================"
echo ""
echo "Próximos passos:"
echo "  1. Revise as estatísticas acima"
echo "  2. Edite configs/default.yaml para configurar:"
echo "     - Critérios de seleção de amostras"
echo "     - Parâmetros do gene a analisar"
echo "     - Caminhos para VCF e referência"
echo "  3. Habilite os passos adicionais em pipeline.steps:"
echo "     - select_samples: true"
echo "     - run_predictions: true"
echo "  4. Execute novamente:"
echo "     python3 build_non_longevous_dataset.py --config configs/default.yaml"
echo ""
echo "Para mais informações, veja: README.md"
echo ""

